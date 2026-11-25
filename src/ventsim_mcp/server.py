"""MCP server exposing Ventsim .vsm model files to an AI assistant."""
from __future__ import annotations

import csv
import io
import os
from collections import defaultdict

from mcp.server.fastmcp import FastMCP

from . import solver, vsm

mcp = FastMCP("ventsim")

_cache: dict[str, tuple[float, vsm.Model]] = {}


def _load(path: str) -> vsm.Model:
    path = os.path.abspath(os.path.expanduser(path))
    mtime = os.path.getmtime(path)
    hit = _cache.get(path)
    if hit and hit[0] == mtime:
        return hit[1]
    model = vsm.parse(path)
    _cache[path] = (mtime, model)
    return model


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
) -> list[dict]:
    """List airways. Filter by name/notes text or minimum velocity (m/s); sort by any
    returned field (e.g. 'Quantity', 'Length', 'Velocity m/s', 'Pressure Branch')."""
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
    return ws[: max(1, min(limit, 500))]


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
    changing text fields. Each edit: {"unique_no": "12", "field": "Notes", "value": "..."}.
    Allowed fields: Notes, Group, Branch Name. Geometry and results are not touched, so
    open the copy in Ventsim and re-simulate."""
    allowed = {"Notes", "Group", "Branch Name"}
    src = os.path.abspath(os.path.expanduser(path))
    dst = os.path.abspath(os.path.expanduser(out_path))
    if src == dst:
        raise ValueError("out_path must differ from the source file")
    if os.path.exists(dst):
        raise ValueError("out_path already exists; choose a new file name")
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
        cells[m.col(e["field"])] = value
        m.lines[m.row_indexes[n]] = "\t".join(cells)
    vsm.dump(m, dst)
    return f"{len(edits)} edit(s) written to {dst}"


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
