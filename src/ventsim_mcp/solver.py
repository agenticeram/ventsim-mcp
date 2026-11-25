"""Experimental steady-state airflow solver for what-if checks on a saved model.

Model: every airway obeys  P_from - P_to = R*Q*|Q| - Ps(Q/n)  where R is the 'Resistance
Total' stored in the file (it already includes shock losses) and Ps is the static pressure
curve of the fan on that airway (n fans in parallel). Node pressures are solved with a
Newton iteration. Boundary nodes (an end touched by one airway, i.e. portals and
dead-ends) are held at the pressures stored in the file.

NOT modelled: natural ventilation pressure, density and heat changes, fixed-flow or
fixed-pressure airways, recirculation, regulators set by flow. Use it to see how a change
in resistance or fans redistributes air, then confirm in Ventsim.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from . import vsm


def _curve(model: vsm.Model) -> dict[str, list[tuple[float, float]]]:
    """Fan type id -> [(flow m3/s, static pressure Pa)] sorted by flow."""
    out: dict[str, list[tuple[float, float]]] = {}
    cur = None
    for r in model.section("FANS"):
        if len(r) > 3 and r[1] and not vsm._is_number(r[1]):
            cur = out.setdefault(r[0], [])
        elif cur is not None and len(r) > 3:
            try:
                cur.append((float(r[3]), float(r[1])))
            except ValueError:
                pass
    return {k: sorted(v) for k, v in out.items() if len(v) > 1}


def _fan_p(curve, q):
    """Static pressure and its slope at flow q (clamped to the tabulated range)."""
    qs = [p[0] for p in curve]
    ps = [p[1] for p in curve]
    if q <= qs[0]:
        return ps[0], 0.0
    if q >= qs[-1]:
        return 0.0 if q > qs[-1] else ps[-1], 0.0
    i = int(np.searchsorted(qs, q)) - 1
    slope = (ps[i + 1] - ps[i]) / (qs[i + 1] - qs[i])
    return ps[i] + slope * (q - qs[i]), slope


def _flow(dp, R, curve, n):
    """Solve  dp = R*q*|q| - Ps(q/n)  for q; return (q, dq/d(dp))."""
    # f(q) = R q|q| - Ps(q/n) - dp is monotone increasing, so bisect.
    lo, hi = -1e4, 1e4
    for _ in range(80):
        q = 0.5 * (lo + hi)
        ps = _fan_p(curve, q / n)[0] if curve else 0.0
        f = R * q * abs(q) - ps - dp
        if f > 0:
            hi = q
        else:
            lo = q
    q = 0.5 * (lo + hi)
    slope = _fan_p(curve, q / n)[1] / n if curve else 0.0
    d = 2 * R * abs(q) - slope
    return q, 1.0 / max(d, 1e-9)


def solve(model: vsm.Model, resistance_scale: dict[str, float] | None = None,
          fans_off: set[str] | None = None, iterations: int = 60, tol: float = 1e-4):
    resistance_scale = resistance_scale or {}
    fans_off = fans_off or set()
    curves = _curve(model)
    rows = list(model.rows())
    key = lambda r, s: tuple(round(model.num(r, f"{c}{s}"), 2) for c in "XYZ")
    ids, nodes = {}, defaultdict(int)
    br = []
    for r in rows:
        a, b = key(r, 1), key(r, 2)
        for k in (a, b):
            ids.setdefault(k, len(ids))
            nodes[k] += 1
        uid = model.value(r, "Unique No")
        R = max(model.num(r, "Resistance Total"), 1e-7) * resistance_scale.get(uid, 1.0)
        ft = model.value(r, "Fan Type")
        curve = curves.get(ft) if ft not in ("", "0") and uid not in fans_off else None
        n = max(model.num(r, "Fan Numbers", 1.0), 1.0)
        br.append((ids[a], ids[b], R, curve, n, uid, a, b, r))
    N = len(ids)
    P = np.zeros(N)
    fixed = np.zeros(N, bool)
    for i, j, R, c, n, uid, a, b, r in br:
        for node, k, col in ((i, a, "Pressure From"), (j, b, "Pressure To")):
            if nodes[k] == 1:
                P[node] = model.num(r, col)
                fixed[node] = True
    free = np.where(~fixed)[0]
    pos = {n: k for k, n in enumerate(free)}
    for it in range(iterations):
        Q = np.zeros(len(br))
        G = np.zeros(len(br))
        for k, (i, j, R, c, n, *_rest) in enumerate(br):
            Q[k], G[k] = _flow(P[i] - P[j], R, c, n)
        imb = np.zeros(N)
        np.add.at(imb, [b[0] for b in br], -Q)
        np.add.at(imb, [b[1] for b in br], Q)
        # imb[i] = inflow - outflow at i; Newton step solves  L dP = imb  on free nodes
        if len(free) == 0 or np.max(np.abs(imb[free])) < tol:
            break
        A = np.zeros((len(free), len(free)))
        for k, (i, j, *_r) in enumerate(br):
            for u, v in ((i, j), (j, i)):
                if u in pos:
                    A[pos[u], pos[u]] += G[k]
                    if v in pos:
                        A[pos[u], pos[v]] -= G[k]
        P[free] += np.linalg.solve(A + 1e-12 * np.eye(len(free)), imb[free])
    return {"uid": [b[5] for b in br], "Q": Q, "iterations": it + 1,
            "residual": float(np.max(np.abs(imb[free]))) if len(free) else 0.0}
