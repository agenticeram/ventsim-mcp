"""Engineering checks on a model: fan selection from the model's fan library, fan energy
and where it is spent, and limit checks (velocity, diesel air, wet bulb, escape distance)."""
from __future__ import annotations

import math

import numpy as np

from . import escape, vsm


def fan_library(model: vsm.Model) -> dict[str, dict]:
    """Fan curves stored in the file: {type id: {name, diameter_m, ref_density, points:
    [(flow m3/s, static Pa, total Pa, static efficiency %, power kW)]}}."""
    names = vsm.fan_types(model)
    out: dict[str, dict] = {}
    cur = None
    for r in model.section("FANS"):
        if len(r) > 3 and r[1] and not vsm._is_number(r[1]):
            try:
                d = float(r[2])
            except ValueError:
                d = 0.0
            try:
                ref = float(r[7]) or 1.2
            except (ValueError, IndexError):
                ref = 1.2
            cur = out.setdefault(r[0], {"name": names.get(r[0], r[1]), "diameter_m": d,
                                        "ref_density": ref, "points": []})
        elif cur is not None and len(r) > 5:
            try:
                q, ps, pt, eff, pw = (float(r[3]), float(r[1]), float(r[2]), float(r[4]),
                                      float(r[5]))
            except ValueError:
                continue
            if q > 0 or ps > 0 or pt > 0:
                cur["points"].append((q, ps, pt, eff, pw))
    for f in out.values():
        f["points"].sort()
    return {k: f for k, f in out.items() if len(f["points"]) >= 2}


def _interp(xs, ys, x):
    return float(np.interp(x, xs, ys))


def fan_select(model: vsm.Model, airflow: float, pressure: float, fans_in_parallel: int = 1,
               density: float = 1.2, max_speed_pct: float = 100.0, pressure_type: str = "static",
               limit: int = 10) -> dict:
    """Fans from the model's library that deliver ``airflow`` (m3/s, total for the set) at
    ``pressure`` (Pa) with ``fans_in_parallel`` units. For each fan: the speed (affinity
    laws) at which its curve passes through the duty point, efficiency and power there, and
    the margin to the stall peak. Ranked by power."""
    if airflow <= 0 or pressure <= 0:
        raise ValueError("airflow and pressure must be positive")
    col = 1 if pressure_type == "static" else 2
    q_fan = airflow / max(int(fans_in_parallel), 1)
    ok, other = [], []
    for k, f in fan_library(model).items():
        pts = f["points"]
        qs = [p[0] for p in pts]
        ps = [p[col] for p in pts]
        if max(ps) <= 0 and col == 1:
            ps = [p[2] for p in pts]  # curve given as total pressure only
        if max(ps) <= 0:
            continue
        kd = density / f["ref_density"]
        # pressure of the curve at speed s for flow q_fan: s^2 * kd * P(q_fan / s)
        def at(s):
            x = q_fan / s
            if x < qs[0] or x > qs[-1]:
                return None
            return s * s * kd * _interp(qs, ps, x)
        lo, hi = 0.2, 2.5
        grid = [lo + i * (hi - lo) / 200 for i in range(201)]
        speed = None
        for a, b in zip(grid, grid[1:]):
            pa, pb = at(a), at(b)
            if pa is None or pb is None:
                continue
            if (pa - pressure) * (pb - pressure) <= 0:
                for _ in range(50):
                    mid = (a + b) / 2
                    pm = at(mid)
                    if pm is None:
                        break
                    if (pa - pressure) * (pm - pressure) <= 0:
                        b = mid
                    else:
                        a, pa = mid, pm
                speed = (a + b) / 2
                break
        p100 = at(1.0)
        row = {"fan_type": k, "name": f["name"], "diameter_m": f["diameter_m"],
               "pressure_at_100pct_Pa": None if p100 is None else round(p100, 1)}
        if speed is None:
            row["note"] = "duty point outside this curve at any speed from 20 to 250 %"
            other.append(row)
            continue
        x = q_fan / speed
        effs = [p[3] for p in pts]
        eff = _interp(qs, effs, x) if max(effs) > 0 else None
        peak_i = int(np.argmax(ps))
        stall_margin = (ps[peak_i] * speed * speed * kd) / pressure - 1
        row.update({"speed_pct": round(speed * 100, 1),
                    "static_efficiency_pct" if col == 1 else "efficiency_pct":
                        None if eff is None else round(eff, 1),
                    "power_per_fan_kW": None if not eff else round(q_fan * pressure / (eff / 100) / 1000, 1),
                    "stall_margin_pct": round(100 * stall_margin, 1),
                    "left_of_stall_peak": x < qs[peak_i]})
        (ok if speed * 100 <= max_speed_pct + 1e-9 and not row["left_of_stall_peak"] else other).append(row)
    key = lambda r: (r.get("power_per_fan_kW") is None, r.get("power_per_fan_kW") or 0, -r["stall_margin_pct"])
    ok.sort(key=key)
    other.sort(key=lambda r: r.get("speed_pct") or 1e9)
    return {"duty": {"airflow_m3s": airflow, "pressure_Pa": pressure, "fans_in_parallel": fans_in_parallel,
                     "airflow_per_fan_m3s": round(q_fan, 3), "density": density,
                     "pressure_type": pressure_type},
            "suitable": ok[:limit],
            "not_suitable": other[:limit],
            "library_size": len(fan_library(model))}


def energy(model: vsm.Model, price_per_kwh: float, hours_per_year: float = 8760.0,
           limit: int = 20) -> dict:
    """Fan power from the saved results, its yearly energy and cost, and the airways where
    the air power is spent (frictional pressure drop x airflow)."""
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    fans, losses = [], []
    for r in rows:
        uid = model.value(r, "Unique No")
        if model.value(r, "Fan Type") not in ("", "0") and model.value(r, "Fan Off") != "True":
            kw = model.num(r, "Fan Power")
            fans.append({"Unique No": uid, "fan": vsm.airway(model, r).get("Fan Name", ""),
                         "power_kW": round(kw, 2),
                         "energy_MWh_per_year": round(kw * hours_per_year / 1000, 1),
                         "cost_per_year": round(kw * hours_per_year * price_per_kwh, 0)})
        else:
            q = model.num(r, "Quantity")
            dp = abs(model.num(r, "Resistance Total") * q * q)
            if dp * abs(q) > 0:
                losses.append({"Unique No": uid, "Branch Name": model.value(r, "Branch Name"),
                               "air_power_kW": dp * abs(q) / 1000, "pressure_drop_Pa": round(dp, 1),
                               "airflow_m3s": round(abs(q), 2),
                               "resistance_type": model.value(r, "Resistance Type")})
    total_fan = sum(f["power_kW"] for f in fans)
    total_loss = sum(x["air_power_kW"] for x in losses) or 1.0
    losses.sort(key=lambda x: -x["air_power_kW"])
    # doors, regulators, stoppings: airways using one of the resistance presets
    reg = sum(x["air_power_kW"] for x in losses
              if x["resistance_type"].isdigit() and int(x["resistance_type"]) > 1)
    for x in losses:
        x["share_pct"] = round(100 * x["air_power_kW"] / total_loss, 2)
        x["air_power_kW"] = round(x["air_power_kW"], 2)
    return {"fan_power_kW": round(total_fan, 1),
            "energy_MWh_per_year": round(total_fan * hours_per_year / 1000, 1),
            "cost_per_year": round(total_fan * hours_per_year * price_per_kwh, 0),
            "price_per_kwh": price_per_kwh, "hours_per_year": hours_per_year,
            "fan_count": len(fans),
            "fans": sorted(fans, key=lambda f: -f["power_kW"])[:limit],
            "air_power_lost_in_airways_kW": round(total_loss, 1),
            "of_which_in_doors_regulators_stoppings_kW": round(reg, 1),
            "largest_losses": losses[:limit]}


DEFAULT_LIMITS = {
    "min_velocity_ms": 0.25,
    "max_velocity_ms": 6.0,
    "max_velocity_shaft_ms": 20.0,
    "air_per_diesel_kw_m3s": 0.06,
    "max_wet_bulb_C": 28.0,
    "max_escape_distance_m": None,
}


def compliance(model: vsm.Model, limits: dict | None = None, wet_bulb: dict[str, float] | None = None,
               limit: int = 30) -> dict:
    """Airways outside the given limits. Defaults are common practice values, not any
    particular regulation: pass your own. Shafts (steeper than 45 degrees) use
    max_velocity_shaft_ms. Closed and impassable airways are skipped."""
    L = dict(DEFAULT_LIMITS)
    for k, v in (limits or {}).items():
        if k not in L:
            raise ValueError(f"unknown limit {k!r}; known: {sorted(L)}")
        L[k] = v
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    imp = escape.impassable_types(model) | {"1"}
    issues: dict[str, list] = {"low_velocity": [], "high_velocity": [], "diesel_air": [],
                               "wet_bulb": [], "escape_distance": []}
    dist = escape.distances(model) if L["max_escape_distance_m"] else {}
    for r in rows:
        uid = model.value(r, "Unique No")
        if model.value(r, "Resistance Type") in imp:
            continue
        name = model.value(r, "Branch Name")
        q, a = abs(model.num(r, "Quantity")), model.num(r, "Area")
        v = q / a if a > 0 else 0.0
        dz = abs(model.num(r, "Z2") - model.num(r, "Z1"))
        horiz = math.hypot(model.num(r, "X2") - model.num(r, "X1"), model.num(r, "Y2") - model.num(r, "Y1"))
        shaft = dz > horiz
        vmax = L["max_velocity_shaft_ms"] if shaft else L["max_velocity_ms"]
        if L["min_velocity_ms"] is not None and v < L["min_velocity_ms"]:
            issues["low_velocity"].append({"Unique No": uid, "Branch Name": name, "velocity_ms": round(v, 3)})
        if vmax is not None and v > vmax:
            issues["high_velocity"].append({"Unique No": uid, "Branch Name": name, "velocity_ms": round(v, 2),
                                            "limit": vmax, "shaft": shaft})
        kw = model.num(r, "dieselkw")
        if kw > 0 and L["air_per_diesel_kw_m3s"]:
            need = kw * L["air_per_diesel_kw_m3s"]
            if q < need:
                issues["diesel_air"].append({"Unique No": uid, "Branch Name": name, "diesel_kW": kw,
                                             "airflow_m3s": round(q, 2), "required_m3s": round(need, 2)})
        if wet_bulb is not None and L["max_wet_bulb_C"] is not None and uid in wet_bulb:
            if wet_bulb[uid] > L["max_wet_bulb_C"]:
                issues["wet_bulb"].append({"Unique No": uid, "Branch Name": name,
                                           "wet_bulb_C": round(wet_bulb[uid], 2)})
        if dist:
            d = dist.get(uid, {}).get("to_refuge")
            if d is None or d > L["max_escape_distance_m"]:
                issues["escape_distance"].append({"Unique No": uid, "Branch Name": name,
                                                  "to_refuge_or_surface_m": None if d is None else round(d, 1)})
    issues["low_velocity"].sort(key=lambda x: x["velocity_ms"])
    issues["high_velocity"].sort(key=lambda x: -x["velocity_ms"])
    issues["wet_bulb"].sort(key=lambda x: -x["wet_bulb_C"])
    issues["escape_distance"].sort(key=lambda x: -(x["to_refuge_or_surface_m"] or 1e9))
    return {"limits": L, "airways_checked": len(rows),
            "counts": {k: len(v) for k, v in issues.items()},
            "wet_bulb_checked": wet_bulb is not None,
            **{k: v[:limit] for k, v in issues.items()}}
