"""Experimental steady-state airflow solver for what-if checks on a saved model.

Model: every airway obeys  P_from - P_to = R*Q*|Q| - Ps(Q/n)  where R is the 'Resistance
Total' stored in the file (it already includes shock losses) and Ps is the static pressure
curve of the fan on that airway (n fans in parallel). Node pressures are solved with a
Newton iteration. Boundary nodes (an end touched by one airway, i.e. portals and
dead-ends) are held at the pressures stored in the file.

Natural ventilation pressure is taken as stored per airway (frozen, so heat or density
changes are not re-evaluated). NOT modelled: heat/density/fire changes, recirculation. Use it to see how a change
in resistance or fans redistributes air, then confirm in Ventsim.
"""
from __future__ import annotations

from collections import defaultdict

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from . import vsm


def _curve(model: vsm.Model) -> dict[str, dict]:
    """Fan type id -> [(flow m3/s, static pressure Pa)] sorted by flow.

    Static pressure is column 1 of each curve row; fans whose column 1 is zero throughout
    keep it in column 2."""
    raw: dict[str, list[tuple[float, float, float]]] = {}
    meta: dict[str, tuple[float, float]] = {}
    cur = None
    for r in model.section("FANS"):
        if len(r) > 3 and r[1] and not vsm._is_number(r[1]):
            cur = raw.setdefault(r[0], [])
            try:
                d = float(r[2])
            except ValueError:
                d = 0.0
            try:
                ref = float(r[7])
            except (ValueError, IndexError):
                ref = 1.2
            meta[r[0]] = (math.pi / 4 * d * d, ref or 1.2)
        elif cur is not None and len(r) > 3:
            try:
                cur.append((float(r[3]), float(r[1]), float(r[2])))
            except ValueError:
                pass
    out = {}
    for k, pts in raw.items():
        if len(pts) < 2:
            continue
        # column 1 is the pressure unless the fan leaves it unset (0 or -999) throughout
        use1 = all(p[1] >= 0 for p in pts) and any(p[1] > 0 for p in pts)
        curve = sorted((q, p1 if use1 else p2) for q, p1, p2 in pts)
        # a fan curve must not rise with flow (stall hump), or the airway equation
        # stops being monotone and the solver wanders
        fixed, top = [], float("inf")
        for q, pr in curve:
            top = min(top, pr)
            fixed.append((q, top))
        out[k] = {"points": fixed, "total": not use1, "area": meta[k][0], "ref": meta[k][1]}
    return out


def _fan_p(curve, q):
    """Static pressure and slope at flow q: flat below the first point, linear
    interpolation inside, and the last slope continued (down to zero) beyond the end."""
    qs = [p[0] for p in curve]
    ps = [p[1] for p in curve]
    if q <= qs[0]:
        return ps[0], 0.0
    if q >= qs[-1]:
        slope = (ps[-1] - ps[-2]) / (qs[-1] - qs[-2]) if qs[-1] > qs[-2] else 0.0
        if slope >= 0:
            return 0.0, 0.0
        p = ps[-1] + slope * (q - qs[-1])
        return (p, slope) if p > 0 else (0.0, 0.0)
    i = int(np.searchsorted(qs, q)) - 1
    slope = (ps[i + 1] - ps[i]) / (qs[i + 1] - qs[i])
    return ps[i] + slope * (q - qs[i]), slope


def _scaled(fan, speed, const, rho=1.2, gain=1.0, eq_total=False):
    """Fan static-pressure function Ps(q) with affinity-law speed scaling and density
    correction, or a constant pressure rise (fixed-pressure airway). Curves given as total
    pressure lose the velocity pressure at the fan outlet."""
    if const is not None:
        return lambda q: (const, 0.0)
    if fan is None:
        return None
    s = max(speed, 1e-3)
    pts = fan["points"]
    k = rho / fan["ref"]
    area = fan["area"]
    # a total-pressure curve loses the outlet velocity pressure unless the airway
    # equation itself uses total pressure
    total = fan["total"] and area > 0 and not eq_total

    def f(q):
        p, d = _fan_p(pts, q / s)
        p, d = k * s * s * p, k * s * d
        if total:
            p -= 0.5 * rho * (q / area) ** 2
            d -= rho * q / (area * area)
        p, d = p * gain, d * gain  # total-pressure fans: see solve()
        return p, d

    return f


def _flow(dp, R, ps_fn, n, qfix=None):
    """Solve  dp = R*q*|q| - Ps(q/n)  for q; return (q, dq/d(dp))."""
    if qfix is not None:
        return qfix, 0.0
    curve = ps_fn
    # f(q) = R q|q| - Ps(q/n) - dp is monotone increasing, so bisect.
    lo, hi = -1e4, 1e4
    for _ in range(80):
        q = 0.5 * (lo + hi)
        ps = curve(q / n)[0] if curve else 0.0
        f = R * q * abs(q) - ps - dp
        if f > 0:
            hi = q
        else:
            lo = q
    q = 0.5 * (lo + hi)
    slope = curve(q / n)[1] / n if curve else 0.0
    d = 2 * R * abs(q) - slope
    return q, 1.0 / max(d, 1e-9)


def solve(model: vsm.Model, resistance_scale: dict[str, float] | None = None,
          fans_off: set[str] | None = None, use_nvp: bool = True, iterations: int = 100, tol: float = 1e-3):
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
        raw = curves.get(ft) if ft not in ("", "0") and uid not in fans_off else None
        # an airway with a pressure source but no fan curve (fixed-pressure airway): keep
        # the stored pressure rise constant
        pfan = model.num(r, "Pressure Fan")
        const = pfan if raw is None and pfan and uid not in fans_off else None
        speed = model.num(r, "Fan Speed", 100.0) / 100.0
        T = model.num(r, "db air In") + 273.15
        pa = model.num(r, "p air In")
        rho = pa / (287.05 * T) if T > 200 and pa > 5e4 else 1.2
        n = max(model.num(r, "Fan Numbers", 1.0), 1.0)
        gain, eq_total = 1.0, False
        if raw is not None:
            # Some fan types apply total pressure in the airway equation, others static.
            # Decide from the saved results of this airway; if total, scale the static
            # curve by the saved total/static ratio (smooth and monotone, unlike adding
            # the velocity pressure explicitly).
            q_s, pf, st = (model.num(r, c) for c in ("Quantity", "Pressure Fan", "Static Pressure"))
            dp_s = model.num(r, "Pressure From") - model.num(r, "Pressure To")
            rq = model.num(r, "Resistance Total") * q_s * abs(q_s)
            if q_s and pf > st and abs(dp_s - (rq - pf)) < abs(dp_s - (rq - st)):
                eq_total = True
                gain = 1.0 if raw["total"] or st <= 1 else pf / st
        curve = _scaled(raw, speed, const, rho, gain, eq_total)
        qf = model.num(r, "QuantityFix")
        qfix = qf if qf and abs(qf - model.num(r, "Quantity")) <= 0.001 * abs(qf) else None
        nvp = model.num(r, "Calc NVP Pressure") if use_nvp else 0.0
        br.append((ids[a], ids[b], R, curve, n, uid, a, b, r, qfix, nvp))
    N = len(ids)
    M = len(br)
    P = np.zeros(N)
    fixed = np.zeros(N, bool)
    for i, j, R, c, n, uid, a, b, r, _qf, _nv in br:
        for node, k, col in ((i, a, "Pressure From"), (j, b, "Pressure To")):
            P[node] = model.num(r, col)  # start from the pressures saved by Ventsim
            if nodes[k] == 1:
                fixed[node] = True
    free = np.where(~fixed)[0]
    pos = -np.ones(N, int)
    pos[free] = np.arange(len(free))
    I = np.array([x[0] for x in br], dtype=int)
    J = np.array([x[1] for x in br], dtype=int)
    Rv = np.array([x[2] for x in br])
    NV = np.array([x[10] for x in br])
    plain = np.array([x[3] is None and x[9] is None for x in br])
    special = [k for k in range(M) if not plain[k]]

    def evaluate(Pv):
        dp = Pv[I] - Pv[J] + NV
        # airways without fan or fixed flow: closed form q = sign(dp) sqrt(|dp|/R)
        Q = np.sign(dp) * np.sqrt(np.abs(dp) / Rv)
        G = 1.0 / np.maximum(2.0 * np.sqrt(np.abs(dp) * Rv), 1e-9)
        for k in special:
            _i, _j, R, c, n, _u, _a, _b, _r, qf, _nv = br[k]
            Q[k], G[k] = _flow(dp[k], R, c, n, qf)
        imb = np.bincount(J, weights=Q, minlength=N) - np.bincount(I, weights=Q, minlength=N)
        return Q, G, imb

    def jacobian(G):
        rows, cols, vals = [], [], []
        for u, v in ((I, J), (J, I)):
            m_u = pos[u] >= 0
            rows.append(pos[u][m_u]); cols.append(pos[u][m_u]); vals.append(G[m_u])
            m_uv = m_u & (pos[v] >= 0)
            rows.append(pos[u][m_uv]); cols.append(pos[v][m_uv]); vals.append(-G[m_uv])
        A = sparse.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(len(free), len(free))).tocsc()
        return A + 1e-9 * sparse.identity(len(free), format="csc")

    Q, G, imb = evaluate(P)
    it = 0
    for it in range(iterations):
        # imb[i] = inflow - outflow at i; Newton step solves  L dP = imb  on free nodes
        if len(free) == 0 or np.max(np.abs(imb[free])) < tol:
            break
        A = jacobian(G)
        step = spsolve(A, imb[free])
        if not np.all(np.isfinite(step)):  # disconnected sub-networks
            step = np.linalg.lstsq(A.toarray(), imb[free], rcond=None)[0]
        norm0 = np.linalg.norm(imb[free])
        alpha = 1.0
        for _ in range(12):  # backtracking line search
            Pn = P.copy()
            Pn[free] += alpha * step
            Qn, Gn, imbn = evaluate(Pn)
            if np.linalg.norm(imbn[free]) < norm0:
                break
            alpha *= 0.5
        P, Q, G, imb = Pn, Qn, Gn, imbn
    return {"uid": [b[5] for b in br], "Q": Q, "iterations": it + 1,
            "residual": float(np.max(np.abs(imb[free]))) if len(free) else 0.0}
