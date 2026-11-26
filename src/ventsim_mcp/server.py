"""MCP server exposing Ventsim .vsm model files to an AI assistant."""
from __future__ import annotations

import copy
import csv
import io
import os
from collections import defaultdict

from mcp.server.fastmcp import FastMCP

from . import calib, design, edit, escape, render, solver, thermo, vsm

mcp = FastMCP("ventsim")

_cache: dict[str, tuple[float, vsm.Model]] = {}

# SI value * factor = imperial value (US mining practice)
_FT, _INWG = 3.28084, 0.00401463
IMPERIAL = {
    "Quantity": (2118.88, "cfm"), "Velocity m/s": (196.85, "fpm"),
    **{c + s: (_FT, "ft") for c in "XYZ" for s in "12"},
    "Width": (_FT, "ft"), "Height": (_FT, "ft"), "Length": (_FT, "ft"),
    "Perimeter": (_FT, "ft"), "Area": (10.7639, "ft2"),
    "Pressure Branch": (_INWG, "in. w.g."), "Pressure Fan": (_INWG, "in. w.g."),
    "Static Pressure": (_INWG, "in. w.g."), "Fan Power": (1.34102, "hp"),
    "Resistance Total": (_INWG / 2118.88 ** 2 * 1e10, "in. w.g. per (100 kcfm)^2"),
}


def _to_units(w: dict, units: str) -> dict:
    u = (units or "SI").lower()
    if u in ("si", "metric"):
        return w
    if u != "imperial":
        raise ValueError("units must be 'SI' or 'imperial'")
    out = dict(w)
    for k, (f, _u) in IMPERIAL.items():
        if isinstance(out.get(k), (int, float)) and not isinstance(out.get(k), bool):
            out[k] = round(out[k] * f, 4)
    if "Velocity m/s" in out:
        out["Velocity fpm"] = out.pop("Velocity m/s")
    out["units"] = "imperial: Quantity cfm, Velocity fpm, lengths ft, Area ft2, pressures in. w.g., Fan Power hp"
    return out


def _load(path: str, stage: int | None = None) -> vsm.Model:
    path = os.path.abspath(os.path.expanduser(path))
    mtime = os.path.getmtime(path)
    hit = _cache.get(path)
    if hit and hit[0] == mtime:
        model = hit[1]
    else:
        model = vsm.parse(path)
        _cache[path] = (mtime, model)
    if stage is None:
        return model
    present = {k for b in vsm._stage_bits(model) for k in range(31) if b >> k & 1}
    if int(stage) not in present:
        raise ValueError(f"stage {stage} not in the model (stages: {sorted(present)}; see list_stages)")
    view = copy.copy(model)  # shares the lines, own stage choice
    view.__dict__ = dict(model.__dict__)
    view.__dict__["_active_stage"] = int(stage)
    return view


@mcp.tool()
def model_summary(path: str) -> dict:
    """Summarise a Ventsim .vsm file: title, version, airway/fan counts, total length."""
    return vsm.summary(_load(path))


@mcp.tool()
def list_airways(
    path: str,
    name_contains: str = "",
    min_velocity: float | None = None,
    sort_by: str = "Unique No",
    descending: bool = False,
    limit: int = 50,
    units: str = "SI",
) -> list[dict]:
    """List airways. Filter by name/notes text or minimum velocity (m/s); sort by any
    returned field (e.g. 'Quantity', 'Length', 'Velocity m/s', 'Pressure Branch').
    units='imperial' returns cfm, fpm, ft, ft2, in. w.g. and hp (filters stay in SI)."""
    ws = vsm.airways(_load(path))
    if name_contains:
        s = name_contains.lower()
        ws = [w for w in ws if s in (w["Branch Name"] + " " + w["Notes"]).lower()]
    if min_velocity is not None:
        ws = [w for w in ws if (w["Velocity m/s"] or 0) >= min_velocity]
    if ws and sort_by not in ws[0]:
        raise ValueError(f"sort_by must be one of {list(ws[0])}")

    def key(w):
        v = w[sort_by]
        return (v is None, v if v is not None else 0)

    if sort_by != "Unique No" or descending:
        ws.sort(key=key, reverse=descending)
    return [_to_units(w, units) for w in ws[: max(1, min(limit, 500))]]


@mcp.tool()
def get_airway(path: str, unique_no: str) -> dict:
    """Every column of one airway (by 'Unique No'), including heat and gas fields."""
    m = _load(path)
    for r in m.rows():
        if m.value(r, "Unique No") == str(unique_no):
            return {c: r[i] for i, c in enumerate(m.columns) if c and i < len(r) and r[i] != ""}
    raise ValueError(f"no airway with Unique No {unique_no}")


@mcp.tool()
def list_fans(path: str) -> list[dict]:
    """Airways that carry a fan, with duty airflow, fan pressure and power."""
    return vsm.fans(_load(path))


@mcp.tool()
def check_flow_balance(path: str, tolerance: float = 0.5, limit: int = 30) -> dict:
    """Mass-balance check of simulated airflow at every junction (mass flow expressed as
    m3/s of standard 1.2 kg/m3 air when the file has it, else volume flow).

    Only airways of the simulated stage are used. Junctions are found by matching airway
    end coordinates. Nodes touched by a single
    airway (portals, dead ends) are skipped. A large imbalance usually means the model
    needs to be re-simulated, contains recirculation fans, or has disconnected nodes."""
    m = _load(path)
    net: dict[int, float] = defaultdict(float)
    deg: dict[int, int] = defaultdict(int)
    where: dict[int, tuple] = {}
    # balance mass flow where Ventsim saved it (heated air expands), else volume
    use_mass = any(m.num(r, "Dry Massflow") for r in m.rows())
    for r, (a, b), active in zip(m.rows(), vsm.node_ids(m), vsm.active_mask(m)):
        if not active:
            continue
        q = m.num(r, "Quantity")
        if use_mass:
            q = m.num(r, "Dry Massflow") / 1.2 * (1 if q >= 0 else -1)
        net[a] -= q  # leaves node a
        net[b] += q  # arrives at node b
        deg[a] += 1
        deg[b] += 1
        where.setdefault(a, tuple(round(m.num(r, c), 2) for c in ("X1", "Y1", "Z1")))
        where.setdefault(b, tuple(round(m.num(r, c), 2) for c in ("X2", "Y2", "Z2")))
    bad = [
        {"x": where[k][0], "y": where[k][1], "z": where[k][2], "imbalance_m3s": round(v, 3),
         "airways": deg[k]}
        for k, v in net.items()
        if deg[k] > 1 and abs(v) > tolerance
    ]
    bad.sort(key=lambda d: -abs(d["imbalance_m3s"]))
    return {
        "junctions_checked": sum(1 for k in net if deg[k] > 1),
        "unbalanced": len(bad),
        "worst": bad[:limit],
    }


@mcp.tool()
def velocity_check(path: str, max_velocity: float = 20.0, min_velocity: float = 0.0,
                   limit: int = 50) -> dict:
    """Find airways whose simulated air velocity (m/s) is above max_velocity or, if
    min_velocity > 0, below it while still carrying air (e.g. dust/gas settling limits).
    Airways with zero airflow are ignored."""
    over, under = [], []
    for w in vsm.airways(_load(path)):
        v = w["Velocity m/s"]
        if not v:
            continue
        item = {k: w[k] for k in ("Unique No", "Branch Name", "Quantity", "Area", "Velocity m/s")}
        if v > max_velocity:
            over.append(item)
        elif min_velocity and v < min_velocity:
            under.append(item)
    over.sort(key=lambda d: -d["Velocity m/s"])
    under.sort(key=lambda d: d["Velocity m/s"])
    return {"over_max": len(over), "under_min": len(under),
            "highest": over[:limit], "lowest": under[:limit]}


@mcp.tool()
def fan_summary(path: str) -> dict:
    """Totals for all fans: count by model, airflow, total fan power (kW)."""
    fs = vsm.fans(_load(path))
    by_model: dict[str, dict] = {}
    for f in fs:
        name = f["Fan Name"] or f["Fan Type"]
        n = f["Fan Numbers"] or 1
        e = by_model.setdefault(name, {"airways": 0, "fans": 0, "power_kw": 0.0})
        e["airways"] += 1
        e["fans"] += n
        e["power_kw"] = round(e["power_kw"] + (f["Fan Power"] or 0), 2)
    return {"total_power_kw": round(sum(e["power_kw"] for e in by_model.values()), 2),
            "by_model": by_model}


@mcp.tool()
def fan_duty_check(path: str, only_problems: bool = False) -> dict:
    """Where every fan of the simulated stage runs on its curve: airflow per fan, position
    between the first and last curve point (0..1), margin to the stall peak, static
    efficiency as saved by Ventsim and a status: ok, near_stall (within 10 % of the peak
    flow), stall_region (left of the pressure peak), beyond_curve, below_curve, no_flow,
    fixed_flow (airflow fixed, the curve is not used), off or no_curve. Uses the results
    saved in the file."""
    fans = solver.fan_duty(_load(path))
    counts: dict[str, int] = defaultdict(int)
    for f in fans:
        counts[f["status"]] += 1
    if only_problems:
        fans = [f for f in fans if f["status"] not in ("ok", "off", "fixed_flow")]
    return {"fans": len(fans), "by_status": dict(counts), "items": fans}


@mcp.tool()
def compare_models(path_a: str, path_b: str, min_change: float = 0.1, limit: int = 50) -> dict:
    """Compare two saved versions of a model by Unique No: added/removed airways and
    those whose airflow (m3/s) or branch pressure (Pa) changed by more than min_change."""
    a = {w["Unique No"]: w for w in vsm.airways(_load(path_a))}
    b = {w["Unique No"]: w for w in vsm.airways(_load(path_b))}
    changed = []
    for k in a.keys() & b.keys():
        dq = (b[k]["Quantity"] or 0) - (a[k]["Quantity"] or 0)
        dp = (b[k]["Pressure Branch"] or 0) - (a[k]["Pressure Branch"] or 0)
        if abs(dq) > min_change or abs(dp) > min_change:
            changed.append({"Unique No": k, "Branch Name": b[k]["Branch Name"],
                            "airflow_a": a[k]["Quantity"], "airflow_b": b[k]["Quantity"],
                            "delta_airflow": round(dq, 3), "delta_pressure": round(dp, 3)})
    changed.sort(key=lambda d: -abs(d["delta_airflow"]))
    return {"added": sorted(b.keys() - a.keys(), key=int), "removed": sorted(a.keys() - b.keys(), key=int),
            "changed": len(changed), "largest_changes": changed[:limit]}


@mcp.tool()
def what_if(path: str, resistance_scale: dict[str, float] | None = None,
            fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
            limit: int = 30) -> dict:
    """EXPERIMENTAL what-if airflow estimate. Re-solves the network with the stored fan
    curves and resistances, scaling resistance of given airways (Unique No -> factor, e.g.
    {"12": 50} for a new door), switching fans off (list of Unique No) and/or changing fan
    speed (Unique No -> speed in %, e.g. {"6": 90}; affinity laws), and returns
    the airflow change per airway relative to an unmodified solve of the same model.
    'baseline_fit' tells you how well the solver reproduces the flows saved by Ventsim
    (mean_rel_error = total absolute flow error / total flow). If 'reliable' is false, or
    junctions are unbalanced in the saved flows (ducting or leakage links the solver cannot
    see), do not trust the numbers. Always confirm in Ventsim."""
    m = _load(path)
    base = solver.solve(m)
    mod = solver.solve(m, resistance_scale or {}, set(fans_off or []),
                       fan_speed={str(k): float(v) for k, v in (fan_speed or {}).items()})
    act = [r for r, k in zip(m.rows(), vsm.active_mask(m)) if k]
    stored = [m.num(r, "Quantity") for r in act]
    names = [m.value(r, "Branch Name") for r in act]
    total = sum(abs(s) for s in stored)
    rel = [sum(abs(abs(b) - abs(s)) for b, s in zip(base["Q"], stored)) / total] if total else []
    out = []
    for uid, name, s, b, q in zip(base["uid"], names, stored, base["Q"], mod["Q"]):
        d = float(q - b)
        if abs(d) > 0.01:
            out.append({"Unique No": uid, "Branch Name": name, "stored_airflow": round(s, 3),
                        "estimated_airflow": round(s + d, 3), "delta": round(d, 3)})
    out.sort(key=lambda x: -abs(x["delta"]))
    return {"baseline_fit": {"mean_rel_error": round(float(rel[0]), 4) if rel else None,
                             "max_abs_error_m3s": round(float(max((abs(abs(b) - abs(s)) for b, s in zip(base["Q"], stored)), default=0)), 3),
                             "residual_m3s": round(mod["residual"], 4),
                             "converged": mod["residual"] < 0.05,
                             "junctions_unbalanced_in_saved_flows": check_flow_balance(path)["unbalanced"],
                             "reliable": bool(mod["residual"] < 0.05 and rel and rel[0] < 0.05)},
            "airways_changed": len(out), "largest_changes": out[:limit]}


@mcp.tool()
def size_regulator(path: str, unique_no: str, target_airflow: float, limit: int = 20) -> dict:
    """EXPERIMENTAL. Resistance to add in one airway (a regulator, door or damper) so that it
    carries target_airflow (m3/s), solved with the what_if network model. Returns the added
    resistance (Ns2/m8), the pressure it drops, an approximate regulator opening area
    (McPherson: A_o = A / (1 + 0.759 A sqrt(R_add)), standard air) and the airflow changes
    elsewhere. A target above the airway's present airflow cannot be reached by a
    regulator; the result then says so."""
    m = _load(path)
    uid = str(unique_no)
    base = solver.solve(m)
    if uid not in base["uid"]:
        raise ValueError(f"no airway with Unique No {uid} in the simulated stage")
    k = base["uid"].index(uid)
    q0 = abs(float(base["Q"][k]))
    target = abs(float(target_airflow))
    R0 = base["branches"][k][2]
    out = {"Unique No": uid, "present_airflow": round(q0, 3), "target_airflow": round(target, 3)}
    if target >= q0:
        out["note"] = "target is not below the present airflow; a regulator can only reduce it"
        return out
    # bisection on log(scale): airflow falls monotonically as resistance rises
    lo, hi, sol = 0.0, 1.0, None
    while True:
        r = solver.solve(m, {uid: 10 ** hi})
        if abs(float(r["Q"][k])) <= target or hi > 9:
            break
        lo, hi = hi, hi + 1.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        sol = solver.solve(m, {uid: 10 ** mid})
        qm = abs(float(sol["Q"][k]))
        if abs(qm - target) <= max(1e-4 * target, 1e-4):
            break
        lo, hi = (mid, hi) if qm > target else (lo, mid)
    scale = 10 ** mid
    r_add = R0 * (scale - 1.0)
    q = abs(float(sol["Q"][k]))
    row = next(r for r in m.rows() if m.value(r, "Unique No") == uid)
    area = m.num(row, "Area")
    out.update({
        "achieved_airflow": round(q, 3),
        "added_resistance_Ns2m8": round(r_add, 5),
        "regulator_pressure_drop_pa": round(r_add * q * q, 1),
        "regulator_open_area_m2": round(area / (1 + 0.759 * area * r_add ** 0.5), 3) if area and r_add > 0 else None,
        "airway_area_m2": area or None,
        "converged": sol["residual"] < 0.05,
    })
    changes = []
    for u, b, q2 in zip(base["uid"], base["Q"], sol["Q"]):
        d = float(q2 - b)
        if abs(d) > 0.01 and u != uid:
            changes.append({"Unique No": u, "delta": round(d, 3)})
    changes.sort(key=lambda x: -abs(x["delta"]))
    out["other_airways_changed"] = len(changes)
    out["largest_other_changes"] = changes[:limit]
    return out


@mcp.tool()
def trace_airflow(path: str, unique_no: str, max_steps: int = 200) -> dict:
    """Follow the main airflow path through one airway: upstream to where the air enters
    the mine and downstream to where it leaves, taking at every junction the airway that
    carries the most air (saved results, simulated stage). Returns the path in flow order
    with length, airflow and the frictional pressure drop (Pa) of each airway and in total."""
    m = _load(path)
    act = [(r, e) for r, e, a in zip(m.rows(), vsm.node_ids(m), vsm.active_mask(m)) if a]
    ways = []
    for r, (a, b) in act:
        q = m.num(r, "Quantity")
        src, dst = (a, b) if q >= 0 else (b, a)
        ways.append({"uid": m.value(r, "Unique No"), "name": m.value(r, "Branch Name"),
                     "from": src, "to": dst, "q": abs(q), "len": m.num(r, "Length"),
                     "dp": abs(m.num(r, "Pressure Branch")), "fan": m.value(r, "Fan Type") not in ("", "0")})
    into: dict[int, list] = defaultdict(list)
    outof: dict[int, list] = defaultdict(list)
    for w in ways:
        if w["q"] > 1e-6:
            into[w["to"]].append(w)
            outof[w["from"]].append(w)
    start = next((w for w in ways if w["uid"] == str(unique_no)), None)
    if start is None:
        raise ValueError(f"no airway with Unique No {unique_no} in the simulated stage")
    seen = {start["uid"]}
    up, node = [], start["from"]
    for _ in range(max_steps):
        cand = [w for w in into.get(node, []) if w["uid"] not in seen]
        if not cand:
            break
        w = max(cand, key=lambda w: w["q"])
        seen.add(w["uid"]); up.append(w); node = w["from"]
    down, node = [], start["to"]
    for _ in range(max_steps):
        cand = [w for w in outof.get(node, []) if w["uid"] not in seen]
        if not cand:
            break
        w = max(cand, key=lambda w: w["q"])
        seen.add(w["uid"]); down.append(w); node = w["to"]
    path_ = list(reversed(up)) + [start] + down
    return {
        "airways": len(path_),
        "total_length_m": round(sum(w["len"] for w in path_), 1),
        "total_friction_pressure_pa": round(sum(w["dp"] for w in path_), 1),
        "fans_on_path": [w["uid"] for w in path_ if w["fan"]],
        "path": [{"Unique No": w["uid"], "Branch Name": w["name"], "airflow": round(w["q"], 3),
                  "length_m": round(w["len"], 1), "pressure_drop_pa": round(w["dp"], 1),
                  "fan": w["fan"]} for w in path_],
    }


@mcp.tool()
def simulate_airflow(path: str, resistance_scale: dict[str, float] | None = None,
                     fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                     use_saved_heat: bool = False, stage: int | None = None,
                     limit: int = 30) -> dict:
    """EXPERIMENTAL. Solve the airflow from the model inputs alone, without the results
    Ventsim saved: airway resistance, fan curves, fan speed and count, fixed-flow and
    fixed-pressure airways, closed airways, surface connections and dead ends. Works on a
    model that was never simulated. Optional changes as in what_if. Natural ventilation and
    heat are not modelled; use_saved_heat=True takes them from the saved results instead.
    If the file has saved airflow, 'vs_saved' compares the two (error = total absolute
    airflow difference / total airflow). stage: simulate another stage (see list_stages)
    instead of the one whose results are saved."""
    m = _load(path, stage)
    s = solver.simulate(m, resistance_scale={str(k): float(v) for k, v in (resistance_scale or {}).items()},
                        fans_off={str(u) for u in (fans_off or [])},
                        fan_speed={str(k): float(v) for k, v in (fan_speed or {}).items()},
                        nvp_from_file=use_saved_heat, density_from_file=use_saved_heat)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    items = []
    for u, q in zip(s["uid"], s["Q"]):
        r = rows[u]
        items.append({"Unique No": u, "Branch Name": m.value(r, "Branch Name"),
                      "airflow": round(float(q), 3), "saved_airflow": round(m.num(r, "Quantity"), 3)})
    out = {"airways": len(items), "converged": s["residual"] < 0.05, "residual_m3s": round(s["residual"], 4),
           "open_ends": s["open_ends"], "dead_ends": s["dead_ends"]}
    tot = sum(abs(i["saved_airflow"]) for i in items)
    if tot:
        diff = [abs(abs(i["airflow"]) - abs(i["saved_airflow"])) for i in items]
        out["vs_saved"] = {"error": round(sum(diff) / tot, 4), "max_abs_diff_m3s": round(max(diff), 3)}
        items.sort(key=lambda i: -abs(abs(i["airflow"]) - abs(i["saved_airflow"])))
    else:
        items.sort(key=lambda i: -abs(i["airflow"]))
    out["items"] = items[:limit]
    return out


def _flows(m: vsm.Model, resistance_scale, fans_off, fan_speed, airflow: str):
    """Airflow per Unique No for the heat and gas tools: None (= saved airflow), re-solved
    with the what_if solver when changes are given, or simulated from the inputs."""
    rs = {str(k): float(v) for k, v in (resistance_scale or {}).items()}
    off = {str(u) for u in (fans_off or [])}
    sp = {str(k): float(v) for k, v in (fan_speed or {}).items()}
    if airflow not in ("saved", "simulated"):
        raise ValueError("airflow must be 'saved' or 'simulated'")
    if airflow == "simulated":
        s = solver.simulate(m, resistance_scale=rs, fans_off=off, fan_speed=sp)
    elif rs or off or sp:
        s = solver.solve(m, resistance_scale=rs, fans_off=off, fan_speed=sp)
    else:
        return None
    return {u: float(q) for u, q in zip(s["uid"], s["Q"])}


def _flow_label(flow, airflow: str) -> str:
    return "saved" if flow is None else ("simulated" if airflow == "simulated" else "re-solved")


@mcp.tool()
def simulate_heat(path: str, extra_heat_kw: dict[str, float] | None = None,
                  resistance_scale: dict[str, float] | None = None,
                  fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                  airflow: str = "saved", surface_temp: float | None = None,
                  extra_latent_kw: dict[str, float] | None = None, stage: int | None = None,
                  limit: int = 30) -> dict:
    """EXPERIMENTAL. Dry- and wet-bulb air temperature and relative humidity through the
    network, mixed at junctions.
    Heat per airway is what Ventsim saved (rock, diesel, fans, other sources) plus
    'extra_heat_kw' ({Unique No: kW}: a new diesel fleet, a fire, a cooler as a negative
    value); auto-compression follows the air. With resistance_scale / fans_off / fan_speed
    the airflow is re-solved first (as in what_if); airflow='simulated' takes it from
    simulate_airflow. A stopped fan's power is taken out of its airway's heat, a speed
    change scales it with speed cubed. Moisture: the water Ventsim saved as evaporating in
    each airway plus 'extra_latent_kw' ({Unique No: kW of evaporation}). The change of
    rock heat with air temperature is not modelled. Returns the hottest airways (by wet
    bulb) and the largest changes against the unchanged model (baseline_fit says how well
    the baseline matches Ventsim)."""
    m = _load(path, stage)
    flow = _flows(m, resistance_scale, fans_off, fan_speed, airflow)
    extra = {str(k): float(v) for k, v in (extra_heat_kw or {}).items()}
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    # a fan's power ends up as heat in its airway: take it out when the fan stops, scale it
    # with speed cubed (affinity laws) when the speed changes
    for u in fans_off or []:
        r = rows.get(str(u))
        if r is not None:
            extra[str(u)] = extra.get(str(u), 0.0) - m.num(r, "Fan Power")
    for u, pct in (fan_speed or {}).items():
        r = rows.get(str(u))
        old = m.num(r, "Fan Speed", 100.0) if r is not None else 0.0
        if r is not None and old > 0 and str(u) not in {str(x) for x in fans_off or []}:
            extra[str(u)] = extra.get(str(u), 0.0) + m.num(r, "Fan Power") * ((float(pct) / old) ** 3 - 1)
    base = thermo.temperatures(m, None, None, surface_temp)
    t = thermo.temperatures(m, flow, extra, surface_temp,
                            {str(k): float(v) for k, v in (extra_latent_kw or {}).items()})
    items, fit_e, fit_w, fit_wb = [], 0.0, 0.0, 0.0
    for u, d in t.items():
        r, b = rows[u], base[u]
        saved = m.num(r, "db air out")
        if b["massflow"] > 0.5:
            fit_e += abs(b["t_out"] - saved) * b["massflow"]
            fit_wb += abs(b["wb_out"] - m.num(r, "wb air out")) * b["massflow"]
            fit_w += b["massflow"]
        items.append({"Unique No": u, "Branch Name": m.value(r, "Branch Name"),
                      "temp_in_C": round(d["t_in"], 2), "temp_out_C": round(d["t_out"], 2),
                      "wet_bulb_out_C": round(d["wb_out"], 2), "rh_out_pct": round(d["rh_out"], 1),
                      "baseline_temp_out_C": round(b["t_out"], 2),
                      "change_C": round(d["t_out"] - b["t_out"], 2),
                      "wet_bulb_change_C": round(d["wb_out"] - b["wb_out"], 2),
                      "baseline_misfit_C": round(b["t_out"] - saved, 2),
                      "massflow_kgs": round(d["massflow"], 2), "heat_kw": round(d["heat_kw"], 2),
                      "latent_kw": round(d["latent_kw"], 2)})
    moving = [i for i in items if i["massflow_kgs"] > 0.5]
    heat_run = any(m.num(r, "Calc Sensible") for r in rows.values())
    return {"airways": len(items), "airflow": _flow_label(flow, airflow),
            "baseline_fit": {
                "heat_simulation_in_file": heat_run,
                "weighted_mean_abs_error_C": round(fit_e / fit_w, 3) if fit_w and heat_run else None,
                "wet_bulb_weighted_mean_abs_error_C": round(fit_wb / fit_w, 3) if fit_w and heat_run else None,
                "note": "baseline = this model with the saved airflow and no changes, compared "
                        "with the temperatures Ventsim saved; changes are reported against "
                        "the baseline. Airways with a large baseline_misfit_C are unreliable."},
            "max_temp_C": max((i["temp_out_C"] for i in moving), default=None),
            "max_wet_bulb_C": max((i["wet_bulb_out_C"] for i in moving), default=None),
            "hottest": sorted(moving, key=lambda i: -i["wet_bulb_out_C"])[:limit],
            "largest_change": sorted(moving, key=lambda i: -max(abs(i["change_C"]),
                                                                abs(i["wet_bulb_change_C"])))[:limit]}


@mcp.tool()
def contaminant_spread(path: str, sources: dict[str, float] | None = None,
                       resistance_scale: dict[str, float] | None = None,
                       fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                       airflow: str = "saved", include_file_sources: bool = False,
                       threshold_ppm: float = 1.0, stage: int | None = None,
                       limit: int = 50) -> dict:
    """EXPERIMENTAL. Steady-state spread of a gas or smoke released in airways:
    'sources' = {Unique No: release rate in m3/s of pure gas} (methane make, CO from a fire,
    averaged blasting fumes...). Fully mixed at junctions, recirculation included.
    include_file_sources=True also reads the model's own 'Base Injection' (m3/s) x 'Base
    Concentration' (ppm) inputs (this reading of Ventsim's fields is not verified). Airflow
    changes as in simulate_heat. Returns the airways above threshold_ppm with their
    concentration (ppm and %), most contaminated first."""
    m = _load(path, stage)
    flow = _flows(m, resistance_scale, fans_off, fan_speed, airflow)
    c = thermo.concentrations(m, {str(k): float(v) for k, v in (sources or {}).items()}, flow,
                              include_file_sources)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    items = []
    for u, d in c.items():
        ppm = d["c_out"] * 1e6
        if ppm >= threshold_ppm and d["airflow"] > 1e-6:
            r = rows[u]
            items.append({"Unique No": u, "Branch Name": m.value(r, "Branch Name"),
                          "ppm": round(ppm, 1), "percent": round(ppm / 1e4, 4),
                          "airflow_m3s": round(d["airflow"], 3)})
    items.sort(key=lambda i: -i["ppm"])
    return {"airways_above_threshold": len(items), "threshold_ppm": threshold_ppm,
            "airflow": _flow_label(flow, airflow), "items": items[:limit]}


@mcp.tool()
def render_network(path: str, out_path: str, view: str = "plan", color_by: str = "Velocity m/s",
                   highlight: list[str] | None = None, width_px: int = 1600) -> str:
    """Draw the airways of the simulated stage to an SVG file (opens in any browser).
    view: 'plan' (X-Y), 'section_x' (X-Z) or 'section_y' (Y-Z). color_by: any numeric
    airway field, e.g. 'Velocity m/s', 'Quantity', 'Pressure Branch'. Line width follows
    airflow, fans are circles, airways in 'highlight' (Unique No, e.g. from trace_airflow)
    are drawn on top in magenta. Hover an airway for its values. Never overwrites a file."""
    dst = os.path.abspath(os.path.expanduser(out_path))
    if not dst.lower().endswith(".svg"):
        raise ValueError("out_path must end with .svg")
    if os.path.exists(dst):
        raise ValueError("out_path already exists; choose a new file name")
    text, n = render.svg(_load(path), view, color_by, {str(u) for u in (highlight or [])}, width_px)
    with open(dst, "w", encoding="utf-8") as f:
        f.write(text)
    return f"{n} airways drawn to {dst}"


@mcp.tool()
def export_csv(path: str, out_path: str, columns: list[str] | None = None) -> str:
    """Write airways to a CSV file. 'columns' may name any column of the .vsm table;
    by default the common fields are exported."""
    m = _load(path)
    out_path = os.path.abspath(os.path.expanduser(out_path))
    cols = columns or list(vsm.BASIC_FIELDS)
    for c in cols:
        m.col(c)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for r in m.rows():
        w.writerow([m.value(r, c) for c in cols])
    with open(out_path, "w", encoding="utf-8-sig", newline="") as f:
        f.write(buf.getvalue())
    return f"{len(m.row_indexes)} airways written to {out_path}"


@mcp.tool()
def edit_airways(path: str, out_path: str, edits: list[dict]) -> str:
    """EXPERIMENTAL. Copy the model to 'out_path' (never overwrites the source) after
    changing fields. Each edit: {"unique_no": "12", "field": "Notes", "value": "..."}.
    Allowed fields: Notes, Group, Branch Name (text), Fan Speed (% of curve speed, 1-200,
    fan airways only), Fan Off (true/false, fan airways only) and the geometry fields Width,
    Height, Length (m) and kFactor (kg/m3): area, perimeter and resistance are rescaled with
    them (Atkinson, R ~ k*P*L/A^3). "Closed": true closes an airway (no air, as a sealed
    stopping; cannot be undone here). Saved results are not recomputed: run
    simulate_airflow on the copy for the new airflow, then confirm in Ventsim."""
    allowed = {"Notes", "Group", "Branch Name", "Fan Speed", "Fan Off", "Closed", *edit.GEOMETRY_FIELDS}
    src, dst = _out_path(path, out_path)
    m = vsm.parse(src)  # fresh copy, not the cached one
    index = {m.value(m.row(n), "Unique No"): n for n in range(len(m.row_indexes))}
    for e in edits:
        if e["field"] not in allowed:
            raise ValueError(f"field not editable: {e['field']} (allowed: {sorted(allowed)})")
        value = str(e["value"])
        if "\t" in value or "\r" in value or "\n" in value:
            raise ValueError("value must not contain tabs or line breaks")
        n = index.get(str(e["unique_no"]))
        if n is None:
            raise ValueError(f"no airway with Unique No {e['unique_no']}")
        cells = m.row(n)
        if e["field"] == "Closed":
            if value.lower() != "true":
                raise ValueError("Closed only accepts true (reopening needs the original resistance)")
            cells[m.col("Resistance Type")] = "1"
            m.lines[m.row_indexes[n]] = "	".join(cells)
            continue
        if e["field"] in edit.GEOMETRY_FIELDS:
            try:
                v = float(value)
            except ValueError:
                raise ValueError(f"{e['field']} must be a number") from None
            edit.set_geometry(m, cells, **{e["field"]: v})
            m.lines[m.row_indexes[n]] = "\t".join(cells)
            continue
        if e["field"] in ("Fan Speed", "Fan Off"):
            if m.value(cells, "Fan Type") in ("", "0"):
                raise ValueError(f"airway {e['unique_no']} has no fan")
            if e["field"] == "Fan Speed":
                try:
                    v = float(value)
                except ValueError:
                    raise ValueError("Fan Speed must be a number (percent)") from None
                if not 1 <= v <= 200:
                    raise ValueError("Fan Speed must be between 1 and 200 %")
                value = f"{v:g}"
            elif value.lower() in ("true", "false"):
                value = "True" if value.lower() == "true" else "False"
            else:
                raise ValueError("Fan Off must be true or false")
        cells[m.col(e["field"])] = value
        m.lines[m.row_indexes[n]] = "\t".join(cells)
    vsm.dump(m, dst)
    return f"{len(edits)} edit(s) written to {dst}"


def _out_path(path: str, out_path: str) -> tuple[str, str]:
    src = os.path.abspath(os.path.expanduser(path))
    dst = os.path.abspath(os.path.expanduser(out_path))
    if src == dst:
        raise ValueError("out_path must differ from the source file")
    if os.path.exists(dst):
        raise ValueError("out_path already exists; choose a new file name")
    return src, dst


@mcp.tool()
def add_airways(path: str, out_path: str, airways: list[dict]) -> dict:
    """EXPERIMENTAL. Copy the model to 'out_path' with new airways added. Each airway:
    {"like": "12", "start": [x, y, z], "end": [x, y, z], "name": "new drive",
    "width": 5.0, "height": 5.0, "open_ends": ["end"]} (name, width, height, open_ends
    optional). It copies size, shape, friction factor and stage from airway 'like' (not its
    fan, fixed flow or pressure, heat or gas sources); resistance follows the new length and
    size. An end within 5 cm of an existing airway end joins the network there (not at a
    surface connection); an end that joins nothing is a dead end unless listed in
    'open_ends' ("start" / "end"), then it is open to the surface. Saved results are not
    recomputed: run simulate_airflow on the copy, then check the copy in Ventsim."""
    src, dst = _out_path(path, out_path)
    m = vsm.parse(src)
    new = [edit.add_airway(m, str(a["like"]), a["start"], a["end"], str(a.get("name", "")),
                           a.get("width"), a.get("height"), tuple(a.get("open_ends") or ()))
           for a in airways]
    vsm.dump(m, dst)
    check = vsm.parse(dst)
    ends = dict(zip((check.value(r, "Unique No") for r in check.rows()), vsm.node_ids(check)))
    deg: dict[int, int] = defaultdict(int)
    for a_, b_ in ends.values():
        deg[a_] += 1
        deg[b_] += 1
    return {"written": dst, "new_unique_no": new,
            "ends_joined": {u: [deg[k] > 1 for k in ends[u]] for u in new}}


@mcp.tool()
def list_stages(path: str) -> dict:
    """Stages (mine plan periods) in the model, their airway counts, and the stage whose
    results are saved (used by default). Other stages can be studied with the 'stage'
    argument of simulate_airflow, simulate_heat, contaminant_spread, escape_routes,
    recirculation_check and compliance_check."""
    m = _load(path)
    bits = vsm._stage_bits(m)
    present = sorted({k for b in bits for k in range(31) if b >> k & 1})
    k = vsm.active_stage(m)
    return {"saved_results_stage": k, "saved_results_stage_name": vsm.stage_name(m, k),
            "stages": [{"stage": s, "name": vsm.stage_name(m, s),
                        "airways": sum(1 for b in bits if b >> s & 1)} for s in present]}


@mcp.tool()
def escape_routes(path: str, unique_no: str | None = None, blocked: list[str] | None = None,
                  stage: int | None = None, limit: int = 30) -> dict:
    """Walking distance from each airway to the surface and to the nearest refuge (or the
    surface, if nearer), as Ventsim computes it: shortest path along airways from the
    middle of the airway; stoppings and other impassable resistances cannot be walked.
    'blocked': extra airways to treat as impassable (a fire, a fall). With unique_no, the
    distances for that airway; otherwise the farthest airways first."""
    m = _load(path, stage)
    d = escape.distances(m, {str(u) for u in (blocked or [])})
    rows = {m.value(r, "Unique No"): r for r in m.rows()}

    def item(u):
        x = d[u]
        return {"Unique No": u, "Branch Name": m.value(rows[u], "Branch Name"),
                "to_surface_m": None if x["to_surface"] is None else round(x["to_surface"], 1),
                "to_refuge_or_surface_m": None if x["to_refuge"] is None else round(x["to_refuge"], 1),
                "is_refuge": x["refuge"], "impassable": x["impassable"]}

    if unique_no is not None:
        if str(unique_no) not in d:
            raise ValueError(f"airway {unique_no} not in the simulated stage")
        return item(str(unique_no))
    walk = [u for u in d if not d[u]["impassable"]]
    cut = [u for u in walk if d[u]["to_surface"] is None]
    far = sorted((u for u in walk if d[u]["to_refuge"] is not None),
                 key=lambda u: -d[u]["to_refuge"])
    return {"airways": len(d), "refuges": sum(1 for x in d.values() if x["refuge"]),
            "cut_off_from_surface": len(cut), "cut_off_examples": [item(u) for u in cut[:limit]],
            "farthest": [item(u) for u in far[:limit]]}


@mcp.tool()
def recirculation_check(path: str, resistance_scale: dict[str, float] | None = None,
                        fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                        airflow: str = "saved", min_flow: float = 0.1, stage: int | None = None) -> dict:
    """Closed loops in the airflow: air that comes back to a point it already passed (fan
    recirculation and other circulating loops). Each group lists its airways, the fans in
    it and the smallest airflow around the loop. Changes and airflow as in simulate_heat."""
    m = _load(path, stage)
    flow = _flows(m, resistance_scale, fans_off, fan_speed, airflow)
    g = escape.recirculation(m, flow, min_flow)
    return {"loops": len(g), "loops_with_fans": sum(1 for x in g if x["fans"]),
            "airflow": _flow_label(flow, airflow), "groups": g}


@mcp.tool()
def contaminant_arrival(path: str, sources: list[str], resistance_scale: dict[str, float] | None = None,
                        fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                        airflow: str = "saved", stage: int | None = None, limit: int = 50) -> dict:
    """Time for smoke or gas released at the start of the 'sources' airways (Unique Nos) to
    reach each airway, carried at the mean air velocity (Ventsim's 'Spread Time' rule; no
    dispersion). Release at t = 1 s as in Ventsim. Returns the reached airways, earliest
    first, and the number never reached."""
    m = _load(path, stage)
    flow = _flows(m, resistance_scale, fans_off, fan_speed, airflow)
    t = thermo.arrival_times(m, [str(s) for s in sources], flow)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    items = sorted(({"Unique No": u, "Branch Name": m.value(rows[u], "Branch Name"),
                     "arrives_s": round(d["arrive_s"], 1), "leaves_s": round(d["leave_s"], 1)}
                    for u, d in t.items()), key=lambda i: i["arrives_s"])
    return {"reached": len(items), "not_reached": len(rows) - len(items),
            "last_arrival_s": items[-1]["leaves_s"] if items else None,
            "airflow": _flow_label(flow, airflow), "items": items[:limit]}


@mcp.tool()
def diesel_exhaust(path: str, emission_g_per_kwh: float | None = None, load_factor: float = 1.0,
                   air_per_kw_m3s: float = 0.06, stage: int | None = None, limit: int = 30) -> dict:
    """EXPERIMENTAL, not compared with Ventsim. Diesel particulate (DPM) concentration from
    the diesel engines in the model ('dieselkw' x load_factor x emission factor, g/kWh;
    the file's own 'dieselDPM' factor per airway when emission_g_per_kwh is not given),
    mixed along the saved airflow; and the airflow against air_per_kw_m3s for each engine."""
    m = _load(path, stage)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    src, engines = {}, []
    for u, r in rows.items():
        kw = m.num(r, "dieselkw")
        if kw <= 0:
            continue
        ef = float(emission_g_per_kwh) if emission_g_per_kwh is not None else m.num(r, "dieselDPM")
        src[u] = kw * load_factor * ef / 3600.0  # g/s
        q = abs(m.num(r, "Quantity"))
        engines.append({"Unique No": u, "Branch Name": m.value(r, "Branch Name"), "diesel_kW": kw,
                        "emission_g_per_kwh": ef, "airflow_m3s": round(q, 2),
                        "required_m3s": round(kw * air_per_kw_m3s, 2), "enough_air": q >= kw * air_per_kw_m3s})
    c = thermo.concentrations(m, src, None, False)
    items = sorted(({"Unique No": u, "Branch Name": m.value(rows[u], "Branch Name"),
                     "dpm_mg_m3": round(d["c_out"] * 1000, 4)} for u, d in c.items()
                    if d["c_out"] > 0 and d["airflow"] > 1e-6), key=lambda i: -i["dpm_mg_m3"])
    return {"engines": engines, "airways_with_dpm": len(items), "items": items[:limit]}


@mcp.tool()
def compliance_check(path: str, limits: dict | None = None, check_wet_bulb: bool = True,
                     stage: int | None = None, limit: int = 30) -> dict:
    """Airways outside ventilation limits: low and high air velocity (shafts have their own
    maximum), airflow per kW of diesel, wet-bulb temperature (from simulate_heat's model,
    when the file has a heat simulation) and walking distance to a refuge or the surface.
    The defaults are common practice values, not a particular regulation; pass your own in
    'limits': min_velocity_ms, max_velocity_ms, max_velocity_shaft_ms,
    air_per_diesel_kw_m3s, max_wet_bulb_C, max_escape_distance_m (null switches a check
    off). Closed and impassable airways are skipped."""
    m = _load(path, stage)
    wb = None
    if check_wet_bulb and any(m.num(r, "Calc Sensible") for r in m.rows()):
        wb = {u: d["wb_out"] for u, d in thermo.temperatures(m).items()}
    return design.compliance(m, limits, wb, limit)


@mcp.tool()
def fan_select(path: str, airflow: float, pressure: float, fans_in_parallel: int = 1,
               density: float = 1.2, max_speed_pct: float = 100.0, pressure_type: str = "static",
               limit: int = 10) -> dict:
    """Pick fans from the fan library stored in the model for a duty: total 'airflow' (m3/s)
    at 'pressure' (Pa, static or total) with 'fans_in_parallel' units. For each fan: the
    speed at which it meets the duty (affinity laws), efficiency, power per fan and margin
    to the stall peak; suitable fans (speed <= max_speed_pct, right of the stall peak)
    ranked by power."""
    if pressure_type not in ("static", "total"):
        raise ValueError("pressure_type must be 'static' or 'total'")
    return design.fan_select(_load(path), airflow, pressure, fans_in_parallel, density,
                             max_speed_pct, pressure_type, limit)


@mcp.tool()
def energy_report(path: str, price_per_kwh: float, hours_per_year: float = 8760.0,
                  limit: int = 20) -> dict:
    """Fan power (as saved), yearly energy and cost at 'price_per_kwh', and where the air
    power goes: the airways with the largest frictional pressure drop x airflow, and the
    part lost in doors, regulators and stoppings."""
    return design.energy(_load(path), price_per_kwh, hours_per_year, limit)


@mcp.tool()
def calibrate_airflow(path: str, measured: dict[str, float], adjustable: list[str] | None = None,
                      max_factor: float = 20.0, from_inputs: bool = False) -> dict:
    """EXPERIMENTAL. Fit airway resistances to survey measurements: 'measured' = {Unique
    No: airflow m3/s}; 'adjustable' = airways whose resistance may change (default: the
    measured ones). Returns the resistance factors, the airflow before and after at the
    measured airways, and the remaining misfit. Feed the factors to what_if /
    simulate_airflow as resistance_scale, or write them with edit_airways (kFactor).
    Measurements that cannot tell airways apart (airways in series) get the smallest
    change; check 'at_limit'."""
    r = calib.calibrate(_load(path), {str(k): float(v) for k, v in measured.items()},
                        [str(u) for u in adjustable] if adjustable else None, max_factor,
                        from_inputs)
    r.pop("flows")
    r["factors"] = {u: round(v, 4) for u, v in r["factors"].items()}
    for k in ("before", "after"):
        r[k] = {u: round(v, 3) for u, v in r[k].items()}
    r["rms_rel_error_before"] = round(r["rms_rel_error_before"], 5)
    r["rms_rel_error_after"] = round(r["rms_rel_error_after"], 5)
    return r


@mcp.tool()
def write_results(path: str, out_path: str, resistance_scale: dict[str, float] | None = None,
                  fans_off: list[str] | None = None, fan_speed: dict[str, float] | None = None,
                  airflow: str = "simulated", with_heat: bool = True, stage: int | None = None) -> str:
    """EXPERIMENTAL. Copy the model to 'out_path' with this server's results written into
    it: airflow ('Quantity', 'Dry Massflow') from simulate_airflow (airflow='simulated',
    default) or the what_if solver (airflow='saved' plus changes), and with_heat the dry-
    and wet-bulb temperatures from simulate_heat. Other saved results are left as they
    were; re-simulate in Ventsim before relying on the copy. Never overwrites a file."""
    src, dst = _out_path(path, out_path)
    m = vsm.parse(src)
    if stage is not None:
        m.__dict__["_active_stage"] = int(_load(path, stage).__dict__["_active_stage"])
    flow = _flows(m, resistance_scale, fans_off, fan_speed, airflow)
    if flow is None:
        flow = {m.value(r, "Unique No"): m.num(r, "Quantity") for r in m.rows()}
    temps = thermo.temperatures(m, flow) if with_heat and any(m.num(r, "Calc Sensible") for r in m.rows()) else {}
    index = {m.value(m.row(n), "Unique No"): n for n in range(len(m.row_indexes))}
    n_written = 0
    for u, q in flow.items():
        n = index[u]
        cells = m.row(n)
        rho = thermo._density(m, cells)
        for col, val in (("Quantity", q), ("Dry Massflow", rho * q)):
            if col in m._col:
                cells[m.col(col)] = f"{val:.9g}"
        if u in temps:
            t = temps[u]
            for col, key in (("db air In", "t_in"), ("db air out", "t_out"), ("wb air In", "wb_in"),
                             ("wb air out", "wb_out")):
                if col in m._col:
                    cells[m.col(col)] = f"{t[key]:.7g}"
        m.lines[m.row_indexes[n]] = "\t".join(cells)
        n_written += 1
    vsm.dump(m, dst)
    return (f"{n_written} airways written to {dst} (airflow: {_flow_label(flow, airflow)}"
            f"{', temperatures' if temps else ''})")


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
