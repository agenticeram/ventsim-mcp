"""Model checks that reproduce the Ventsim warnings and errors which can be found from the
file alone (Ventsim manual, Appendix E), and fixes for the ones with a safe, mechanical
repair.

Found here: loose ends without a surface or closed-end setting (Ventsim: "no entry",
"no exit", "no entry or exit"), duplicate airways (same two junctions), zero-length
airways, running fans in closed airways ("operating fan blocked"), fans in series ("fan
interfering with another fan"), stoppings in series ("airway stopping redundant"), and
fixed flows that conflict or cannot run. Warnings that need a simulation (stalled fans,
recirculation, density variation, heat) are not covered here: fan_duty_check and
recirculation_check report those from saved results.
"""
from __future__ import annotations

import math
from collections import defaultdict

from . import escape, solver, vsm

SURFACE_FLAGS = (1 << 12, 1 << 13)  # entry end, exit end open to the surface
DEAD_END_FLAG = 1 << 18             # 'Closed End' in the Ventsim edit box


def _int(model: vsm.Model, r: list[str], col: str) -> int:
    try:
        return int(float(model.value(r, col) or 0))
    except ValueError:
        return 0


def _closed_types(model: vsm.Model) -> set[str]:
    return escape.impassable_types(model) | {"1"}


def _fan_running(model: vsm.Model, r: list[str]) -> bool:
    return model.value(r, "Fan Type") not in ("", "0") and model.value(r, "Fan Off") != "True"


def _fix(model: vsm.Model, r: list[str]) -> float | None:
    q = model.num(r, "QuantityFix")
    if not q:
        return None
    mask = solver.QFIX_MASK
    if mask is not None and not _int(model, r, "Branch Flags") & mask:
        return None
    return q


def _network(model: vsm.Model):
    mask = vsm.active_mask(model)
    idx = [n for n, k in enumerate(mask) if k]
    rows = [model.row(n) for n in idx]
    ends = [e for e, k in zip(vsm.node_ids(model), mask) if k]
    deg: dict[int, int] = defaultdict(int)
    for a, b in ends:
        deg[a] += 1
        deg[b] += 1
    return idx, rows, ends, deg


def _chains(ends, deg) -> list[list[int]]:
    """Airways joined end to end through junctions of exactly two airways (one airflow
    path with no branches), as lists of positions in ``ends``."""
    at: dict[int, list[int]] = defaultdict(list)
    for i, (a, b) in enumerate(ends):
        at[a].append(i)
        at[b].append(i)
    seen: set[int] = set()
    out = []
    for i in range(len(ends)):
        if i in seen:
            continue
        chain, stack = [], [i]
        seen.add(i)
        while stack:
            j = stack.pop()
            chain.append(j)
            for node in ends[j]:
                if deg[node] == 2:
                    for k in at[node]:
                        if k not in seen:
                            seen.add(k)
                            stack.append(k)
        out.append(chain)
    return out


def _stuck(n_nodes_hint, ends, passable, open_nodes) -> set[int]:
    """Airways (positions) that cannot carry air: not on any closed path through the
    network and the atmosphere, once the closed airways are taken out. Those are the
    bridges of the graph where every open end joins one 'atmosphere' node."""
    ATM = -1
    adj: dict[int, list[tuple[int, int]]] = defaultdict(list)
    edges = []
    for i, (a, b) in enumerate(ends):
        if not passable[i] or a == b:
            continue
        edges.append(i)
        adj[a].append((b, i))
        adj[b].append((a, i))
    for k in open_nodes:
        e = ("atm", k)
        adj[k].append((ATM, e))
        adj[ATM].append((k, e))
    disc, low, bridges = {}, {}, set()
    t = 0
    for root in list(adj):
        if root in disc:
            continue
        disc[root] = low[root] = t; t += 1
        stack = [(root, None, iter(adj[root]))]
        while stack:
            v, via, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                stack.pop()
                if stack:
                    u = stack[-1][0]
                    low[u] = min(low[u], low[v])
                    if low[v] > disc[u] and isinstance(via, int):
                        bridges.add(via)
                continue
            w, e = nxt
            if e == via:
                continue
            if w in disc:
                low[v] = min(low[v], disc[w])
            else:
                disc[w] = low[w] = t; t += 1
                stack.append((w, e, iter(adj[w])))
    return bridges


def check(model: vsm.Model, limit: int = 50) -> dict:
    """Ventsim-style warnings and errors found in the file. Each item names the airway(s),
    what Ventsim would report, and the fix (``fix_model`` applies the ones marked
    automatic)."""
    idx, rows, ends, deg = _network(model)
    uid = [model.value(r, "Unique No") for r in rows]
    name = [model.value(r, "Branch Name") for r in rows]
    closed = _closed_types(model)
    issues: dict[str, list] = defaultdict(list)

    # loose ends: Ventsim's 'no entry' / 'no exit' / 'no entry or exit'. A loose end with
    # no flag at all that carried air in the saved results was treated as open by Ventsim
    # (duct outlets in HVAC models), so only ends without airflow are reported.
    for i, (r, (a, b)) in enumerate(zip(rows, ends)):
        f = _int(model, r, "Branch Flags")
        if f & DEAD_END_FLAG or abs(model.num(r, "Quantity")) > 1e-3:
            continue
        loose = [deg[k] == 1 and not f & SURFACE_FLAGS[s] for s, k in enumerate((a, b))]
        if not any(loose):
            continue
        kind = ("no_entry_or_exit" if all(loose) else "no_entry" if loose[0] else "no_exit")
        issues[kind].append({"Unique No": uid[i], "Branch Name": name[i],
                             "fix": "closed_end (automatic) or connect / open to surface"})

    # duplicates: two airways between the same two junctions
    by_pair: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i, (a, b) in enumerate(ends):
        if a != b:
            by_pair[(min(a, b), max(a, b))].append(i)
    for group in by_pair.values():
        if len(group) < 2:
            continue
        auto = not any(_fan_running(model, rows[i]) or _fix(model, rows[i]) is not None
                       or model.value(rows[i], "Resistance Type") in closed for i in group)
        issues["duplicate"].append({
            "Unique No": [uid[i] for i in group], "Branch Name": [name[i] for i in group],
            "reversed": [ends[i] != ends[group[0]] for i in group],
            "fix": ("merge into one airway with Pathways = number of airways (automatic)" if auto
                    else "has a fan, fixed flow or stopping: merge or delete by hand")})

    # zero length: both ends at the same junction
    for i, (a, b) in enumerate(ends):
        if a == b:
            issues["zero_length"].append({"Unique No": uid[i], "Branch Name": name[i],
                                          "fix": "move one end or delete the airway"})

    # a running fan or a fixed flow on an airway that air cannot pass through once the
    # closed airways are taken out: Ventsim errors 'operating fan blocked' and
    # 'fix stopped, check connections are joined/open'
    passable = [model.value(r, "Resistance Type") not in closed for r in rows]
    stuck = _stuck(None, ends, passable, escape.open_nodes(model, rows, ends))
    for i, r in enumerate(rows):
        if not passable[i] or i in stuck:
            if _fan_running(model, r):
                issues["fan_blocked"].append({
                    "Unique No": uid[i], "Branch Name": name[i],
                    "fix": "turn the fan off (automatic) or open a path for its air"})
            if _fix(model, r) is not None and _fix(model, r) != 0:
                issues["fixed_flow_blocked"].append({
                    "Unique No": uid[i], "Branch Name": name[i],
                    "fix": "remove the fixed flow or open a path for its air"})

    # one airflow path (no branches): fans in series, stoppings in series, fixed flows
    for chain in _chains(ends, deg):
        fans = [i for i in chain if _fan_running(model, rows[i])]
        stops = [i for i in chain if model.value(rows[i], "Resistance Type") in closed]
        fixes = [(i, _fix(model, rows[i])) for i in chain if _fix(model, rows[i]) is not None]
        if len(fans) > 1:
            issues["fans_in_series"].append({
                "Unique No": [uid[i] for i in fans],
                "fix": "remove a fan or give the air another path between them"})
        if len(stops) > 1:
            issues["stoppings_in_series"].append({
                "Unique No": [uid[i] for i in stops], "fix": "none needed; one stopping is enough"})
        qs = [abs(q) for _i, q in fixes]
        if qs and max(qs) - min(qs) > 1e-3 * max(qs):
            issues["conflicting_fixed_flow"].append({
                "Unique No": [uid[i] for i, _q in fixes],
                "fixed_flows_m3s": [q for _i, q in fixes],
                "fix": "keep one fixed flow on this path"})
    ERRORS = {"fan_blocked", "conflicting_fixed_flow", "fixed_flow_blocked", "zero_length"}
    counts = {k: len(v) for k, v in issues.items()}
    return {"airways_checked": len(rows),
            "errors": {k: n for k, n in counts.items() if k in ERRORS},
            "warnings": {k: n for k, n in counts.items() if k not in ERRORS},
            **{k: v[:limit] for k, v in issues.items()}}


def _set_int(model: vsm.Model, cells: list[str], col: str, value: int) -> None:
    cells[model.col(col)] = str(value)


def delete_airways(model: vsm.Model, uids: set[str]) -> int:
    """Remove airway rows (in place) and renumber the row index column."""
    keep_lines = []
    gone = 0
    first, last = model.row_indexes[0], model.row_indexes[-1]
    for n, i in enumerate(model.row_indexes):
        cells = model.lines[i].split("\t")
        if model.value(cells, "Unique No") in uids:
            gone += 1
            continue
        keep_lines.append(cells)
    for k, cells in enumerate(keep_lines, start=1):
        cells[0] = str(k)
    model.lines[first:last + 1] = ["\t".join(c) for c in keep_lines]
    model.row_indexes = list(range(first, first + len(keep_lines)))
    return gone


def fix(model: vsm.Model, kinds: set[str]) -> dict:
    """Apply automatic fixes in place. kinds: 'loose_ends' (mark as closed end),
    'duplicates' (merge into one airway with Pathways), 'blocked_fans' (turn off)."""
    report = check(model, limit=10 ** 9)
    done: dict[str, list] = defaultdict(list)
    index = {model.value(model.row(n), "Unique No"): n for n in range(len(model.row_indexes))}

    def put(n, cells):
        model.lines[model.row_indexes[n]] = "\t".join(cells)

    if "loose_ends" in kinds:
        for k in ("no_entry", "no_exit", "no_entry_or_exit"):
            for it in report.get(k, []):
                n = index[it["Unique No"]]
                cells = model.row(n)
                _set_int(model, cells, "Branch Flags", _int(model, cells, "Branch Flags") | DEAD_END_FLAG)
                put(n, cells)
                done["closed_end"].append(it["Unique No"])
    if "blocked_fans" in kinds:
        for it in report.get("fan_blocked", []):
            n = index[it["Unique No"]]
            cells = model.row(n)
            cells[model.col("Fan Off")] = "True"
            put(n, cells)
            done["fan_off"].append(it["Unique No"])
    if "duplicates" in kinds:
        drop: set[str] = set()
        for it in report.get("duplicate", []):
            if "automatic" not in it["fix"]:
                continue
            ids = it["Unique No"]
            rs = [model.row(index[u]) for u in ids]
            inv = sum(1 / math.sqrt(max(model.num(r, "Resistance Total"), 1e-12)) for r in rs)
            r_eq = 1 / inv ** 2
            n0 = index[ids[0]]
            cells = model.row(n0)
            ratio = r_eq / max(model.num(cells, "Resistance Total"), 1e-12)
            for col in ("Resistance", "Resistance Total"):
                if col in model._col and model.num(cells, col) > 0:
                    cells[model.col(col)] = f"{model.num(cells, col) * ratio:.9g}"
            paths = sum(max(int(model.num(r, "Pathways", 1)), 1) for r in rs)
            if "Pathways" in model._col:
                cells[model.col("Pathways")] = str(paths)
            put(n0, cells)
            drop.update(ids[1:])
            done["merged"].append({"kept": ids[0], "removed": ids[1:], "Pathways": paths})
        if drop:
            delete_airways(model, drop)
    return dict(done)
