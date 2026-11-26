"""Analysis tools from Ventsim's analysis menu: shortest path, resistance sensitivity,
goal seek and steady-state radon. Experimental; each states what it approximates."""
from __future__ import annotations

import heapq
import math
import random

import numpy as np

from . import solver, thermo, vsm

RADON_DECAY = 2.1e-6  # 1/s, Rn-222 (3.82 day half-life); decay along an airway neglected


def shortest_path(model: vsm.Model, src: str, dst: str, weight: str = "distance") -> dict:
    """Shortest route between two airways (through their nearest ends) over the simulated
    stage. weight: 'distance' (Length), 'resistance' (Resistance Total) or 'hops'."""
    if weight not in ("distance", "resistance", "hops"):
        raise ValueError("weight must be distance, resistance or hops")
    net, _N = thermo._network(model, None)
    info = {uid: (r, a, b) for r, uid, a, b, _q in net}
    if src not in info or dst not in info:
        raise ValueError("source or destination airway is not in the simulated stage")
    adj: dict[int, list] = {}
    for r, uid, a, b, _q in net:
        if weight == "distance":
            w = model.num(r, "Length")
        elif weight == "resistance":
            w = model.num(r, "Resistance Total")
        else:
            w = 1.0
        w = max(w, 1e-9)
        adj.setdefault(a, []).append((b, w, uid))
        adj.setdefault(b, []).append((a, w, uid))
    _rs, sa, sb = info[src]
    _rd, da, db = info[dst]
    dist: dict[int, float] = {}
    prev: dict[int, tuple[int, str]] = {}
    heap: list[tuple[float, int]] = []
    for n in (sa, sb):
        dist[n] = 0.0
        heapq.heappush(heap, (0.0, n))
    targets = {da, db}
    end = None
    while heap:
        d, n = heapq.heappop(heap)
        if d > dist.get(n, math.inf):
            continue
        if n in targets:
            end = n
            break
        for m, w, uid in adj.get(n, []):
            nd = d + w
            if nd < dist.get(m, math.inf):
                dist[m] = nd
                prev[m] = (n, uid)
                heapq.heappush(heap, (nd, m))
    if end is None:
        raise ValueError("no route between the two airways")
    path: list[str] = []
    n = end
    while n in prev:
        p, uid = prev[n]
        path.append(uid)
        n = p
    path.reverse()
    seq = [src] + [u for u in path if u not in (src, dst)] + [dst]
    seq = list(dict.fromkeys(seq))
    length = sum(model.num(info[u][0], "Length") for u in seq)
    resistance = sum(model.num(info[u][0], "Resistance Total") for u in seq)
    return {"weight": weight, "airways": seq, "count": len(seq),
            "length_m": round(length, 2), "resistance_sum": round(resistance, 6),
            "note": "shortest by the chosen weight over airway ends; Ventsim's menu rules "
                    "(horizontal-only, duct-free, follow airflow) are not all reproduced"}


def sensitivity(model: vsm.Model, targets: list[str] | None = None, runs: int = 100,
                spread_pct: float = 10.0, seed: int = 1) -> dict:
    """Monte Carlo on resistances: each airway's resistance is multiplied by a uniform
    factor in 1 +/- spread_pct %, the network re-solved, and the flow statistics reported.
    Confidence is the share of runs whose flow lies within 5 % of the mean."""
    if runs < 2 or runs > 1000:
        raise ValueError("runs must be between 2 and 1000")
    rng = random.Random(seed)
    base = solver.solve(model)
    uids = list(base["uid"])
    tgt = targets or uids[:]
    for t in tgt:
        if t not in uids:
            raise ValueError(f"airway {t} is not in the simulated stage")
    samples: dict[str, list[float]] = {t: [] for t in tgt}
    for _ in range(runs):
        scale = {u: 1.0 + rng.uniform(-spread_pct, spread_pct) / 100.0 for u in uids}
        s = solver.solve(model, resistance_scale=scale)
        q = dict(zip(s["uid"], map(float, s["Q"])))
        for t in tgt:
            samples[t].append(abs(q[t]))
    out = {}
    for t in tgt:
        arr = np.array(samples[t])
        mean = float(arr.mean())
        std = float(arr.std(ddof=1))
        within = float(np.mean(np.abs(arr - mean) <= 0.05 * abs(mean))) if mean else 0.0
        out[t] = {"mean_m3s": round(mean, 4), "std_m3s": round(std, 4),
                  "deviation_pct": round(100 * std / mean, 2) if mean else None,
                  "min_m3s": round(float(arr.min()), 4), "max_m3s": round(float(arr.max()), 4),
                  "confidence_within_5pct": round(100 * within, 1)}
    return {"runs": runs, "spread_pct": spread_pct, "seed": seed, "airways": out,
            "note": "uniform resistance perturbation; Ventsim's own distribution and its "
                    "'Direction Confidence' are not reproduced"}


def goal_seek(model: vsm.Model, fan_uid: str, target_uid: str, target_m3s: float,
              speed_min: float = 20.0, speed_max: float = 200.0, tol: float = 0.01,
              max_iter: int = 30) -> dict:
    """Find the fan speed (% of curve speed) that gives ``target_m3s`` in airway
    ``target_uid``, by bisection. Flow is assumed to change monotonically with speed over
    the bracket."""
    def flow(speed: float) -> float:
        s = solver.solve(model, fan_speed={fan_uid: speed})
        return abs(dict(zip(s["uid"], map(float, s["Q"])))[target_uid])

    lo, hi = speed_min, speed_max
    f_lo, f_hi = flow(lo), flow(hi)
    if not (min(f_lo, f_hi) <= target_m3s <= max(f_lo, f_hi)):
        raise ValueError(f"target {target_m3s} m3/s is outside the flow range "
                         f"{min(f_lo, f_hi):.3f}-{max(f_lo, f_hi):.3f} m3/s over speeds "
                         f"{speed_min}-{speed_max}%")
    rising = f_hi > f_lo
    it, mid, f_mid = 0, lo, f_lo
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        f_mid = flow(mid)
        if abs(f_mid - target_m3s) <= tol:
            break
        if (f_mid < target_m3s) == rising:
            lo = mid
        else:
            hi = mid
    return {"fan": fan_uid, "target_airway": target_uid, "target_m3s": target_m3s,
            "speed_pct": round(mid, 3), "achieved_m3s": round(f_mid, 4), "iterations": it,
            "converged": abs(f_mid - target_m3s) <= tol,
            "note": "result is not saved to the model; apply it with edit_airways (Fan Speed)"}


def radon_steady(model: vsm.Model, surface_bq_s_m2: dict[str, float] | None = None,
                 point_bq_s: dict[str, float] | None = None) -> dict:
    """Steady-state radon in every airway: surface exhalation (Bq/s/m2 over the airway's
    wall area Perimeter x Length) plus point sources (Bq/s), carried by the airflow and
    mixed at junctions. Concentration in Bq/m3 and pCi/l (1 pCi/l = 37 Bq/m3). Radon decay
    along an airway is neglected (transit times are minutes)."""
    net, N = thermo._network(model, None)
    src_bq: dict[str, float] = {}
    for r, uid, _a, _b, _q in net:
        s = 0.0
        if surface_bq_s_m2 and uid in surface_bq_s_m2:
            s += surface_bq_s_m2[uid] * model.num(r, "Perimeter") * model.num(r, "Length")
        if point_bq_s and uid in point_bq_s:
            s += point_bq_s[uid]
        if s:
            src_bq[uid] = s
    links, info = [], []
    for r, uid, a, b, q in net:
        u, v = (a, b) if q >= 0 else (b, a)
        aq = abs(q)
        s = src_bq.get(uid, 0.0)
        inc = s / aq if aq > thermo.MIN_MASSFLOW and s else 0.0
        links.append((u, v, aq, inc))
        info.append((uid, u, aq, inc))
    C = thermo._mix(N, links, {}, 0.0)
    out = {}
    for uid, u, aq, inc in info:
        c = float(C[u]) + inc
        out[uid] = {"bq_m3": round(c, 3), "pci_l": round(c / 37.0, 4), "airflow_m3s": round(aq, 4)}
    return {"airways": out, "max_bq_m3": max((d["bq_m3"] for d in out.values()), default=0.0),
            "note": "steady state, well mixed at junctions; progeny and occupancy dose not modelled"}
