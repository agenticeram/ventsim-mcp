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
        out[k] = {"points": fixed, "raw": curve, "total": not use1, "area": meta[k][0],
                  "ref": meta[k][1]}
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
          fans_off: set[str] | None = None, use_nvp: bool = True, iterations: int = 100, tol: float = 5e-3,
          fan_speed: dict[str, float] | None = None):
    """fan_speed: Unique No -> new speed in % of the curve's speed (affinity laws)."""
    resistance_scale = resistance_scale or {}
    fan_speed = fan_speed or {}
    fans_off = fans_off or set()
    curves = _curve(model)
    # only the airways of the simulated stage form the network
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    ends = [e for e, k in zip(vsm.node_ids(model), mask) if k]
    nodes: dict[int, int] = defaultdict(int)
    br = []
    dens: list[float] = []
    base: list[tuple] = []  # unmodified resistance and fan curve of each airway
    for r, (a, b) in zip(rows, ends):
        for k in (a, b):
            nodes[k] += 1
        uid = model.value(r, "Unique No")
        R = max(model.num(r, "Resistance Total"), 1e-7)
        q_saved = model.num(r, "Quantity")
        pb = model.num(r, "Pressure Branch")
        if abs(q_saved) > 1e-3 and pb and -pb / (q_saved * abs(q_saved)) > 0:
            # effective resistance including fittings, dampers and orifices: Ventsim saves
            # the frictional/shock pressure drop as 'Pressure Branch' (= -R_eff*Q|Q|)
            R = max(-pb / (q_saved * abs(q_saved)), 1e-7)
        R *= resistance_scale.get(uid, 1.0)
        ft = model.value(r, "Fan Type")
        # A fan that was stopped when Ventsim saved the results is a plain airway. The
        # 'Fan Off' flag is not reliable on its own (it can belong to another stage), so
        # trust the saved state: the pressure drop is explained by resistance alone.
        q_sv = model.num(r, "Quantity")
        dp_sv = model.num(r, "Pressure From") - model.num(r, "Pressure To") + (
            model.num(r, "Calc NVP Pressure") if use_nvp else 0.0)
        r_sv = model.num(r, "Resistance Total") * q_sv * abs(q_sv)
        off_in_file = (model.value(r, "Fan Off") == "True"
                       and abs(dp_sv - r_sv) < 0.05 * max(model.num(r, "Pressure Fan"), 1.0))
        raw_all = curves.get(ft) if ft not in ("", "0") and not off_in_file else None
        raw = raw_all if uid not in fans_off else None
        # an airway with a pressure source but no fan curve (fixed-pressure airway): keep
        # the stored pressure rise constant
        pfan = model.num(r, "Pressure Fan")
        const = pfan if raw is None and pfan and uid not in fans_off and not off_in_file else None
        speed = model.num(r, "Fan Speed", 100.0) / 100.0
        T = model.num(r, "db air In") + 273.15
        pa = model.num(r, "p air In")
        rho = pa / (287.05 * T) if T > 200 and pa > 5e4 else 1.2
        n = max(model.num(r, "Fan Numbers", 1.0), 1.0)
        gain, eq_total = 1.0, False
        if raw_all is not None:
            # Some fan types apply total pressure in the airway equation, others static.
            # Decide from the saved results of this airway; if total, scale the static
            # curve by the saved total/static ratio (smooth and monotone, unlike adding
            # the velocity pressure explicitly).
            q_s, pf, st = (model.num(r, c) for c in ("Quantity", "Pressure Fan", "Static Pressure"))
            dp_s = model.num(r, "Pressure From") - model.num(r, "Pressure To")
            rq = model.num(r, "Resistance Total") * q_s * abs(q_s)
            if q_s and pf > st and abs(dp_s - (rq - pf)) < abs(dp_s - (rq - st)):
                eq_total = True
                gain = 1.0 if raw_all["total"] or st <= 1 else pf / st
        if raw_all is not None and q_saved:
            # calibrate the curve so it passes through the saved duty point (covers
            # speed control, blade settings and pressure definitions the file does not
            # spell out); the curve shape is kept for what-if changes
            f0 = _scaled(raw_all, speed, None, rho, gain, eq_total)
            need = (R / resistance_scale.get(uid, 1.0)) * q_saved * abs(q_saved) - (
                model.num(r, "Pressure From") - model.num(r, "Pressure To")
                + (model.num(r, "Calc NVP Pressure") if use_nvp else 0.0))
            have = f0(q_saved / n)[0]
            if have > 1 and need > 0:
                gain *= min(max(need / have, 0.2), 5.0)
        new_speed = fan_speed.get(uid, speed * 100.0) / 100.0
        if const is not None and uid in fan_speed and speed > 0:
            const *= (new_speed / speed) ** 2  # fixed-pressure source scales like a fan
        curve = _scaled(raw, new_speed, const, rho, gain, eq_total)
        pfan_c = pfan if raw_all is None and pfan and not off_in_file else None
        base.append((R / resistance_scale.get(uid, 1.0), _scaled(raw_all, speed, pfan_c, rho, gain, eq_total)))
        qf = model.num(r, "QuantityFix")
        qfix = qf if qf and abs(qf - model.num(r, "Quantity")) <= 0.001 * abs(qf) else None
        dp_saved = model.num(r, "Pressure From") - model.num(r, "Pressure To")
        closed = model.value(r, "Resistance Type") == "1" or (abs(q_saved) < 1e-6 and abs(dp_saved) > 0.01)
        if closed and uid not in resistance_scale:
            qfix = 0.0  # closed airway (stopping/sealed): carries no air
        nvp = model.num(r, "Calc NVP Pressure") if use_nvp else 0.0
        if qfix is None and abs(q_saved) > 1e-6:
            # Pressure terms the file does not itemise (velocity-pressure changes, Ventsim's
            # own convergence residual): whatever the saved state needs beyond resistance,
            # fan and natural ventilation is kept as a fixed pressure, like the NVP.
            R0 = R / resistance_scale.get(uid, 1.0)
            c0 = _scaled(raw_all, speed, pfan if raw_all is None and pfan and not off_in_file else None,
                         rho, gain, eq_total)
            ps0 = c0(q_saved / n)[0] if c0 else 0.0
            nvp -= (dp_saved + nvp) - (R0 * q_saved * abs(q_saved) - ps0)
        # Ventsim conserves mass, not volume: air heated or compressed along the way changes
        # volume. Keep each airway's saved density and balance mass flow at the nodes.
        q_s, mf = model.num(r, "Quantity"), model.num(r, "Dry Massflow")
        dens.append(mf / q_s if q_s and mf and 0.3 < mf / q_s < 3 else rho)
        br.append((a, b, R, curve, n, uid, a, b, r, qfix, nvp))
    N = (max(max(e) for e in ends) + 1) if ends else 0
    P = np.zeros(N)
    fixed = np.zeros(N, bool)
    for i, j, R, c, n, uid, a, b, r, _qf, _nv in br:
        for node, k, col in ((i, a, "Pressure From"), (j, b, "Pressure To")):
            P[node] = model.num(r, col)  # start from the pressures saved by Ventsim
            if nodes[k] == 1:
                fixed[node] = True
    # A boundary node (portal, dead end) is held at the pressure that makes its airway
    # carry the saved flow: saved boundary pressures are sometimes reported on another
    # datum (elevation, surface reference) than the interior ones.
    for (i, j, R, c, n, uid, a, b, r, qf, nv), (R0, c0) in zip(br, base):
        ea, eb = nodes[a] == 1, nodes[b] == 1
        if ea == eb or qf is not None:
            continue
        q = model.num(r, "Quantity")
        drop = R0 * q * abs(q) - (c0(q / n)[0] if c0 else 0.0) - nv  # = P_i - P_j
        if ea:
            P[i] = P[j] + drop
        else:
            P[j] = P[i] - drop
    return _network(br, dens, P, fixed, N, iterations, tol)



def _network(br, dens, P, fixed, N, iterations, tol, dp_floor=0.0):
    """Newton iteration on node pressures for the airways in ``br`` (mass balance with the
    densities ``dens``); nodes marked in ``fixed`` keep their pressure."""
    M = len(br)
    free = np.where(~fixed)[0]
    pos = -np.ones(N, int)
    pos[free] = np.arange(len(free))
    I = np.array([x[0] for x in br], dtype=int)
    J = np.array([x[1] for x in br], dtype=int)
    Rv = np.array([x[2] for x in br])
    NV = np.array([x[10] for x in br])
    RHO = np.array(dens) if dens else np.zeros(0)
    plain = np.array([x[3] is None and x[9] is None for x in br])
    special = [k for k in range(M) if not plain[k]]

    def evaluate(Pv):
        dp = Pv[I] - Pv[J] + NV
        # airways without fan or fixed flow: closed form q = sign(dp) sqrt(|dp|/R)
        Q = np.sign(dp) * np.sqrt(np.abs(dp) / Rv)
        G = 1.0 / np.maximum(2.0 * np.sqrt(np.maximum(np.abs(dp), dp_floor) * Rv), 1e-9)
        for k in special:
            _i, _j, R, c, n, _u, _a, _b, _r, qf, _nv = br[k]
            Q[k], G[k] = _flow(dp[k], R, c, n, qf)
        Mf = Q * RHO
        imb = np.bincount(J, weights=Mf, minlength=N) - np.bincount(I, weights=Mf, minlength=N)
        return Q, G * RHO, imb

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
        if len(free) == 0 or np.max(np.abs(imb[free])) < tol * 1.2:
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
    return {"uid": [b[5] for b in br], "Q": Q, "P": P, "branches": br, "iterations": it + 1,
            "residual": float(np.max(np.abs(imb[free]))) / 1.2 if len(free) else 0.0}

DEAD_END_FLAG = 1 << 18  # 'Branch Flags' bit Ventsim sets on airways ending in a dead end
SURFACE_FLAGS = (1 << 12, 1 << 13)  # 'Branch Flags' bits: entry / exit end open to surface
PRESSURE_FIX_FLAG = 1 << 19  # 'Pressure Fix' is only applied when this bit is set
TOTAL_AS_IS = True  # total-pressure curves act as total pressure (no outlet loss)
QFIX_MASK = (1 << 2) | (1 << 4)  # Branch Flags bits that switch a QuantityFix on


def simulate(model: vsm.Model, iterations: int = 200, tol: float = 5e-3,
             resistance_scale: dict[str, float] | None = None, fans_off: set[str] | None = None,
             fan_speed: dict[str, float] | None = None, nvp_from_file: bool = False,
             density_from_file: bool = False, nvp: dict[str, float] | None = None,
             density: dict[str, float] | None = None) -> dict:
    """Solve the airflow from the model inputs only, without any saved results: airway
    resistance ('Resistance Total'), fan curves at their speed and count, fixed-flow and
    fixed-pressure airways, closed airways ('Resistance Type' 1) and fans switched off.
    Loose ends are at atmospheric pressure unless Ventsim flagged them as dead ends; air
    is incompressible at 1.2 kg/m3 and natural ventilation is ignored (no heat model)."""
    resistance_scale = resistance_scale or {}
    fans_off = fans_off or set()
    fan_speed = fan_speed or {}
    curves = _curve(model)
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    ends = [list(e) for e, k in zip(vsm.node_ids(model), mask) if k]

    def flags_of(r):
        try:
            return int(float(model.value(r, "Branch Flags") or 0))
        except ValueError:
            return 0

    def degrees():
        d: dict[int, int] = defaultdict(int)
        for a, b in ends:
            d[a] += 1
            d[b] += 1
        return d

    deg = degrees()
    br, dead, open_ = [], set(), set()
    for r, (a, b) in zip(rows, ends):
        uid = model.value(r, "Unique No")
        R = max(model.num(r, "Resistance Total"), 1e-7) * resistance_scale.get(uid, 1.0)
        ft = model.value(r, "Fan Type")
        off = model.value(r, "Fan Off") == "True" or uid in fans_off
        raw = curves.get(ft) if ft not in ("", "0") and not off else None
        speed = fan_speed.get(uid, model.num(r, "Fan Speed", 100.0) or 100.0) / 100.0
        flags = flags_of(r)
        for s_, k in enumerate((a, b)):
            # a loose end is open to the atmosphere when flagged as a surface connection, or
            # when it carries no dead-end flag at all (e.g. duct outlets in HVAC models)
            if deg[k] == 1 and (flags & SURFACE_FLAGS[s_] or not flags & DEAD_END_FLAG):
                open_.add(k)
        pfix = model.num(r, "Pressure Fix") if flags & PRESSURE_FIX_FLAG else 0.0
        const = pfix if raw is None and pfix and not off else None
        n = max(model.num(r, "Fan Numbers", 1.0), 1.0)
        curve = _scaled(raw, speed, const, 1.2, eq_total=TOTAL_AS_IS)
        qf = model.num(r, "QuantityFix")
        qfix = qf if qf and (QFIX_MASK is None or flags & QFIX_MASK) else None
        if model.value(r, "Resistance Type") == "1" and uid not in resistance_scale:
            qfix = 0.0
        nv = model.num(r, "Calc NVP Pressure") if nvp_from_file else 0.0
        if nvp is not None:
            nv = nvp.get(uid, 0.0)
        br.append((a, b, R, curve, n, uid, a, b, r, qfix, nv))
    N = (max(max(e) for e in ends) + 1) if ends else 0
    P = np.zeros(N)
    fixed = np.zeros(N, bool)
    for k in open_:
        fixed[k] = True  # open to the atmosphere; other loose ends are dead ends
    dens = [1.2] * len(br)
    if density_from_file:
        dens = []
        for b_ in br:
            q_s, mf = model.num(b_[8], "Quantity"), model.num(b_[8], "Dry Massflow")
            dens.append(mf / q_s if q_s and mf and 0.3 < mf / q_s < 3 else 1.2)
    if density is not None:
        dens = [density.get(b_[5], 1.2) for b_ in br]
    out = _network(br, dens, P, fixed, N, iterations, tol, dp_floor=1e-3)
    out["open_ends"] = int(fixed.sum())
    out["dead_ends"] = sum(1 for k, d in deg.items() if d == 1 and k not in open_)
    return out


def fan_duty(model: vsm.Model) -> list[dict]:
    """Where each fan of the simulated stage runs on its curve.

    The duty flow per fan (airway flow / number of fans) is brought back to the curve's
    reference speed by the affinity laws and placed on the curve as published (before the
    stall hump is flattened for the solver). Left of the pressure peak the curve rises with
    flow: that is the stall region, where parallel fans can share air unevenly and
    airflow can surge. Beyond the last point the fan runs off its published curve."""
    curves = _curve(model)
    names = vsm.fan_types(model)
    out = []
    for r, active in zip(model.rows(), vsm.active_mask(model)):
        ft = model.value(r, "Fan Type")
        if not active or ft in ("", "0"):
            continue
        q = model.num(r, "Quantity")
        n = max(model.num(r, "Fan Numbers", 1.0), 1.0)
        speed = model.num(r, "Fan Speed", 100.0) or 100.0
        pfan = model.num(r, "Pressure Fan")
        power = model.num(r, "Fan Power")
        d = {"Unique No": model.value(r, "Unique No"), "Branch Name": model.value(r, "Branch Name"),
             "Fan Name": names.get(ft, ft), "fans": int(n), "speed_pct": round(speed, 1),
             "airflow_per_fan_m3s": round(abs(q) / n, 3), "fan_pressure_pa": round(pfan, 1),
             "power_kw": round(power, 2),
             "static_efficiency_pct": round(model.num(r, "Pressure Eff"), 1) or None}
        c = curves.get(ft)
        qf = model.num(r, "QuantityFix")
        dp = model.num(r, "Pressure From") - model.num(r, "Pressure To") + model.num(r, "Calc NVP Pressure")
        rq = model.num(r, "Resistance Total") * q * abs(q)
        if qf and abs(qf - q) <= 0.001 * abs(qf):
            d["status"] = "fixed_flow"  # pressure is whatever forces the fixed airflow
            d["static_efficiency_pct"] = None
        elif model.value(r, "Fan Off") == "True" and abs(dp - rq) < 0.05 * max(abs(pfan), 1.0):
            d["status"] = "off"
        elif c is None:
            d["status"] = "no_curve"
        else:
            raw = c["raw"]
            qr = abs(q) / n / (speed / 100.0)
            qmin, qmax = raw[0][0], raw[-1][0]
            peak = max(raw, key=lambda t: t[1])[0]
            d["curve_flow_range_m3s"] = [round(qmin, 3), round(qmax, 3)]
            d["position_on_curve"] = round((qr - qmin) / (qmax - qmin), 3) if qmax > qmin else None
            if peak > qmin:
                d["stall_peak_flow_m3s"] = round(peak, 3)
                d["stall_margin"] = round((qr - peak) / peak, 3)
            if abs(q) < 1e-6:
                d["status"] = "no_flow"
            elif peak > qmin and qr < peak:
                d["status"] = "stall_region"
            elif qr > qmax * 1.001:
                d["status"] = "beyond_curve"
            elif qr < qmin * 0.999:
                d["status"] = "below_curve"
            elif peak > qmin and qr < peak * 1.1:
                d["status"] = "near_stall"
            else:
                d["status"] = "ok"
        out.append(d)
    return out
