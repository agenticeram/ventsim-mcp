"""Fan speed search: the lowest estimated fan power at which every airway stays within the
velocity limits. Airflow comes from the what_if solver (incompressible, saved state);
power is scaled from the fan power saved in the file (power ~ flow x speed^2 at the same
efficiency). Both are approximations: check the chosen copy in Ventsim before relying on it.
"""
from __future__ import annotations

import math
from itertools import product

from . import design, solver, vsm

DEFAULT_LIMITS = {"min_velocity_ms": 0.25, "max_velocity_ms": 6.0, "max_velocity_shaft_ms": 20.0}


def _areas_and_shaft(model: vsm.Model):
    mask = vsm.active_mask(model)
    out = {}
    for r, k in zip(model.rows(), mask):
        if not k:
            continue
        uid = model.value(r, "Unique No")
        dz = abs(model.num(r, "Z2") - model.num(r, "Z1"))
        horiz = math.hypot(model.num(r, "X2") - model.num(r, "X1"), model.num(r, "Y2") - model.num(r, "Y1"))
        out[uid] = (model.num(r, "Area"), dz > horiz, model.value(r, "Resistance Type"))
    return out


def _violations(model, flows, limits, areas):
    bad = []
    for uid, q in flows.items():
        if uid not in areas:
            continue
        a, shaft, rtype = areas[uid]
        if rtype == "1" or a <= 0:
            continue  # closed airways carry no air
        v = abs(q) / a
        vmax = limits["max_velocity_shaft_ms"] if shaft else limits["max_velocity_ms"]
        if v > vmax or v < limits["min_velocity_ms"]:
            bad.append({"Unique No": uid, "velocity_ms": round(v, 2)})
    return bad


def _has_fixed_flow(model: vsm.Model, uid: str) -> bool:
    for r in model.rows():
        if model.value(r, "Unique No") == uid:
            qf = model.num(r, "QuantityFix")
            flags = int(float(model.value(r, "Branch Flags") or 0))
            return bool(qf) and bool(flags & solver.QFIX_MASK)
    return False


def _power(model: vsm.Model, fans: dict[str, float], base: dict[str, float], flows: dict[str, float]) -> float:
    """Estimated shaft power (kW) of the tuned fans: saved power x flow ratio x speed ratio^2."""
    total = 0.0
    for r in model.rows():
        uid = model.value(r, "Unique No")
        if uid not in fans or model.value(r, "Fan Type") in ("", "0"):
            continue
        p0 = model.num(r, "Fan Power")
        q0 = abs(model.num(r, "Quantity")) or 1e-9
        q = abs(flows.get(uid, q0))
        s0 = max(base[uid], 1e-6)
        total += p0 * (q / q0) * (fans[uid] / s0) ** 2
    return total


def optimize_fans(model: vsm.Model, fans: list[str] | None = None, limits: dict | None = None,
                  min_speed: float = 50.0, max_speed: float = 100.0, step: float = 5.0,
                  max_evaluations: int = 400, strict: bool = False,
                  min_flow: dict[str, float] | None = None) -> dict:
    """Search the speeds of the given fans (default: all running fans). One fan: every
    speed on the grid. Several fans: each gets its own speed on the grid, up to
    max_evaluations combinations; beyond that a common speed (% of each fan's saved speed)
    is searched instead."""
    L = dict(DEFAULT_LIMITS)
    L.update(limits or {})
    if not set(L) <= set(DEFAULT_LIMITS):
        raise ValueError(f"unknown limit {sorted(set(L) - set(DEFAULT_LIMITS))}")
    tune = [str(u) for u in fans] if fans else [
        model.value(r, "Unique No") for r in model.rows()
        if model.value(r, "Fan Type") not in ("", "0") and model.value(r, "Fan Off") != "True"]
    if not tune:
        raise ValueError("no running fans to tune")
    fixed = [u for u in tune if _has_fixed_flow(model, u)]
    if fixed:
        raise ValueError(f"fan airway(s) {fixed} carry a fixed flow: their airflow does not follow "
                         "the fan speed, so they cannot be tuned here")
    if not min_flow and not strict:
        raise ValueError("give min_flow (airflow each named airway must keep, m3/s) or strict=True: "
                         "without a requirement the lowest speed always wins")
    base = {model.value(r, "Unique No"): model.num(r, "Fan Speed", 100.0) or 100.0 for r in model.rows()
            if model.value(r, "Unique No") in tune}
    grid = [min_speed + i * step for i in range(int(round((max_speed - min_speed) / step)) + 1)]
    areas = _areas_and_shaft(model)
    # guard: a fan whose speed does not change the solved airflow cannot be tuned by this search
    lo_s = solver.solve(model, fan_speed={u: base[u] * min_speed / 100.0 for u in tune})
    hi_s = solver.solve(model, fan_speed={u: base[u] * max_speed / 100.0 for u in tune})
    dq = max(abs(a - b) for a, b in zip(map(float, lo_s["Q"]), map(float, hi_s["Q"])))
    if dq < 1e-3 * max(1.0, max(abs(float(q)) for q in hi_s["Q"])):
        raise ValueError("changing these fans' speed does not change the solved airflow in this "
                         "model (check for fixed flows or pressures on the fan airways)")
    baseline = solver.solve(model)
    base_flows = dict(zip(baseline["uid"], map(float, baseline["Q"])))
    base_power = _power(model, {u: base[u] for u in tune}, base, base_flows)
    # Many example models already break the velocity limits at their saved speeds. Unless
    # strict, a speed set is feasible when it creates no violation the saved state did not have.
    already = {v["Unique No"] for v in _violations(model, base_flows, L, areas)}

    def evaluate(speeds: dict[str, float]):
        s = solver.solve(model, fan_speed=speeds)
        flows = dict(zip(s["uid"], map(float, s["Q"])))
        bad = _violations(model, flows, L, areas)
        new = [v for v in bad if v["Unique No"] not in already or strict]
        short = [{"Unique No": u, "airflow_m3s": round(abs(flows.get(u, 0.0)), 2), "required_m3s": f}
                 for u, f in (min_flow or {}).items() if abs(flows.get(u, 0.0)) < f]
        new += [{"Unique No": x["Unique No"], "airflow_short": True} for x in short]
        return {"speeds_pct": {u: round(v, 1) for u, v in speeds.items()},
                "power_kW": round(_power(model, speeds, base, flows), 1),
                "violations": bad[:10], "violation_count": len(bad),
                "new_violations": new[:10], "ok": not new}

    if len(tune) == 1:
        candidates = [{tune[0]: g} for g in grid]
        mode = "grid"
    elif len(grid) ** len(tune) <= max_evaluations:
        candidates = [dict(zip(tune, combo)) for combo in product(grid, repeat=len(tune))]
        mode = "grid"
    else:
        factors = [min_speed + i * step for i in range(int(round((max_speed - min_speed) / step)) + 1)]
        candidates = [{u: base[u] * f / 100.0 for u in tune} for f in factors]
        mode = "common_factor"
    tried = []
    for speeds in candidates:
        if len(tried) >= max_evaluations:
            break
        tried.append(evaluate(speeds))
    ok = [t for t in tried if t["ok"]]
    best = min(ok, key=lambda t: t["power_kW"]) if ok else None
    return {"mode": mode, "fans": tune, "limits": L, "strict": strict, "min_flow": min_flow or {},
            "baseline_power_kW": round(base_power, 1),
            "violations_at_saved_state": len(already),
            "evaluated": len(tried), "best": best,
            "feasible_count": len(ok),
            "note": ("Power is estimated from the saved fan power, not from the fan curves; "
                     "airflow is incompressible. Verify the chosen speeds in Ventsim.")}
