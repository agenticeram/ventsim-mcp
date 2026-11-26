"""Travel distance to the surface and to refuges, and recirculation (air that returns to a
point it already passed).

Distances follow Ventsim's own rule (reproduces its 'distancetosurface' / 'distancetorefuge'
columns): shortest path along airways, measured from the middle of an airway, from any
surface connection or from the middle of a refuge airway; airways whose resistance type is
marked impassable (stoppings, roof falls, ...) cannot be walked and get no distance.
"""
from __future__ import annotations

import heapq
from collections import defaultdict

from . import vsm

SURFACE_FLAGS = (1 << 12, 1 << 13)
DEAD_END_FLAG = 1 << 18
REFUGE_FLAG2 = 1 << 5  # 'Branch Flags2' bit of a refuge station airway


def _int(model: vsm.Model, r: list[str], col: str) -> int:
    try:
        return int(float(model.value(r, col) or 0))
    except ValueError:
        return 0


def impassable_types(model: vsm.Model) -> set[str]:
    """Resistance preset numbers flagged as impassable in the RESISTANCES table."""
    return {x[0] for x in model.section("RESISTANCES") if len(x) > 1 and x[-1] == "True"}


def _active(model: vsm.Model):
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    ends = [e for e, k in zip(vsm.node_ids(model), mask) if k]
    return rows, ends


def open_nodes(model: vsm.Model, rows, ends) -> set[int]:
    deg: dict[int, int] = defaultdict(int)
    for a, b in ends:
        deg[a] += 1
        deg[b] += 1
    out = set()
    for r, (a, b) in zip(rows, ends):
        f = _int(model, r, "Branch Flags")
        for s, k in enumerate((a, b)):
            if deg[k] == 1 and (f & SURFACE_FLAGS[s] or not f & DEAD_END_FLAG):
                out.add(k)
    return out


def _dijkstra(adj, seeds: dict[int, float]) -> dict[int, float]:
    d = dict(seeds)
    h = [(x, n) for n, x in seeds.items()]
    heapq.heapify(h)
    while h:
        x, n = heapq.heappop(h)
        if x > d.get(n, float("inf")):
            continue
        for v, L in adj[n]:
            if x + L < d.get(v, float("inf")):
                d[v] = x + L
                heapq.heappush(h, (x + L, v))
    return d


def distances(model: vsm.Model, blocked: set[str] | None = None) -> dict[str, dict]:
    """{Unique No: {"to_surface": m or None, "to_refuge": m or None, "refuge": bool,
    "impassable": bool}}. ``to_refuge`` is the nearer of a refuge and the surface, as in
    Ventsim. ``blocked``: extra airways to treat as impassable (e.g. a fire)."""
    blocked = blocked or set()
    rows, ends = _active(model)
    imp = impassable_types(model)
    adj = defaultdict(list)
    info = []
    for r, (a, b) in zip(rows, ends):
        uid = model.value(r, "Unique No")
        L = model.num(r, "Length")
        no_pass = model.value(r, "Resistance Type") in imp or uid in blocked
        refuge = bool(_int(model, r, "Branch Flags2") & REFUGE_FLAG2)
        if not no_pass:
            adj[a].append((b, L))
            adj[b].append((a, L))
        info.append((uid, a, b, L, no_pass, refuge))
    surf = {n: 0.0 for n in open_nodes(model, rows, ends)}
    ds = _dijkstra(adj, surf)
    seeds = dict(surf)
    for uid, a, b, L, no_pass, refuge in info:
        if refuge and not no_pass:
            for n in (a, b):
                seeds[n] = min(seeds.get(n, float("inf")), L / 2)
    dr = _dijkstra(adj, seeds)
    out = {}
    inf = float("inf")
    for uid, a, b, L, no_pass, refuge in info:
        if no_pass:
            out[uid] = {"to_surface": None, "to_refuge": None, "refuge": refuge, "impassable": True}
            continue
        s = min(ds.get(a, inf), ds.get(b, inf)) + L / 2
        q = 0.0 if refuge else min(dr.get(a, inf), dr.get(b, inf)) + L / 2
        out[uid] = {"to_surface": None if s == inf else s, "to_refuge": None if q == inf else q,
                    "refuge": refuge, "impassable": False}
    return out


def recirculation(model: vsm.Model, flow: dict[str, float] | None = None,
                  min_flow: float = 0.1) -> list[dict]:
    """Groups of airways where the air goes round in a closed loop (strongly connected parts
    of the directed airflow graph), with the airflow that circulates and the fans in them."""
    rows, ends = _active(model)
    edges = []
    for r, (a, b) in zip(rows, ends):
        uid = model.value(r, "Unique No")
        q = flow[uid] if flow is not None and uid in flow else model.num(r, "Quantity")
        if abs(q) < min_flow:
            continue
        u, v = (a, b) if q > 0 else (b, a)
        edges.append((u, v, uid, abs(q), model.value(r, "Fan Type") not in ("", "0")))
    succ = defaultdict(list)
    for u, v, *_ in edges:
        succ[u].append(v)
    # Tarjan, iterative
    index, low, on, stack, comp = {}, {}, set(), [], {}
    counter = [0]
    nodes = {n for e in edges for n in e[:2]}
    for root in nodes:
        if root in index:
            continue
        work = [(root, iter(succ[root]))]
        index[root] = low[root] = counter[0]; counter[0] += 1
        stack.append(root); on.add(root)
        while work:
            n, it = work[-1]
            nxt = next(it, None)
            if nxt is not None:
                if nxt not in index:
                    index[nxt] = low[nxt] = counter[0]; counter[0] += 1
                    stack.append(nxt); on.add(nxt)
                    work.append((nxt, iter(succ[nxt])))
                elif nxt in on:
                    low[n] = min(low[n], index[nxt])
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[n])
            if low[n] == index[n]:
                while True:
                    x = stack.pop(); on.discard(x); comp[x] = n
                    if x == n:
                        break
    groups = defaultdict(list)
    for u, v, uid, q, fan in edges:
        if comp[u] == comp[v]:
            groups[comp[u]].append((uid, q, fan))
    out = []
    for g in groups.values():
        out.append({"airways": [u for u, _q, _f in g],
                    "fans": [u for u, _q, f in g if f],
                    "min_loop_airflow_m3s": round(min(q for _u, q, _f in g), 3)})
    out.sort(key=lambda g: -g["min_loop_airflow_m3s"])
    return out
