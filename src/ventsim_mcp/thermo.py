"""Steady-state transport along the airflow: air temperature and gas/contaminant
concentration, fully mixed at junctions (recirculation included, solved as one linear system).

Temperature: each airway adds its sensible heat (kW) and the auto-compression of air moving
down (g*dz/cp). The heat per airway is the sensible heat Ventsim saved ('Calc Sensible':
rock, diesel, fans, point and linear sources together, auto-compression taken out; below
1 kg/s the heat that explains the saved temperature rise) plus any
extra source given, so after an airflow change the heat loads stay as they were and the
temperatures follow the new flow. Moisture (latent heat) is not modelled.

Concentration: a source of gas (m3/s) in an airway adds G/Q to the concentration of the air
passing through it.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from . import vsm

CP = 1005.0  # J/kg/K, dry air
G = 9.81
LOW_FLOW = 1.0  # kg/s
MAX_HEAT_RISE = 30.0  # K per airway
MIN_MASSFLOW = 1e-3  # kg/s; below this an airway carries no heat or gas of its own


def _mix(N: int, links: list[tuple[int, int, float, float]], fixed: dict[int, float],
         default: float) -> np.ndarray:
    """Node values for links (u -> v, weight w, increment b): at every node with inflow,
    value = sum(w * (x_u + b)) / sum(w). Nodes without inflow take ``fixed`` or ``default``."""
    W = np.zeros(N)
    for u, v, w, b in links:
        W[v] += w
    rows, cols, vals = [], [], []
    rhs = np.zeros(N)
    for k in range(N):
        if W[k] <= 0 or k in fixed:
            rows.append(k); cols.append(k); vals.append(1.0)
            rhs[k] = fixed.get(k, default)
        else:
            rows.append(k); cols.append(k); vals.append(W[k])
    for u, v, w, b in links:
        if W[v] > 0 and v not in fixed:
            rows.append(v); cols.append(u); vals.append(-w)
            rhs[v] += w * b
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()
    return spsolve(A, rhs)


def _network(model: vsm.Model, flow: dict[str, float] | None):
    """Rows of the simulated stage with node ids and signed airflow (end 1 -> end 2)."""
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    ends = [e for e, k in zip(vsm.node_ids(model), mask) if k]
    out = []
    for r, (a, b) in zip(rows, ends):
        uid = model.value(r, "Unique No")
        q = flow[uid] if flow is not None and uid in flow else model.num(r, "Quantity")
        out.append((r, uid, a, b, q))
    N = max((max(a, b) for _r, _u, a, b, _q in out), default=-1) + 1
    return out, N


def _density(model: vsm.Model, r: list[str]) -> float:
    q, mf = model.num(r, "Quantity"), model.num(r, "Dry Massflow")
    rho = mf / q if q and mf else 0.0
    return rho if 0.3 < rho < 3 else 1.2


HFG = 2450.0  # kJ/kg, latent heat of evaporation used for moisture pick-up


def p_sat(t: float) -> float:
    """Saturation vapour pressure over water, Pa (Magnus form)."""
    return 610.94 * math.exp(17.625 * t / (t + 243.04))


def moisture(db: float, wb: float, p: float) -> float:
    """Moisture content, kg water / kg dry air, from dry and wet bulb (ASHRAE)."""
    ps = p_sat(wb)
    xs = 0.622 * ps / max(p - ps, 1.0)
    return ((2501 - 2.326 * wb) * xs - 1.006 * (db - wb)) / (2501 + 1.86 * db - 4.186 * wb)


def wet_bulb(db: float, x: float, p: float) -> float:
    """Wet-bulb temperature for dry bulb, moisture content and pressure (bisection);
    capped at saturation."""
    db = min(max(db, -40.0), 120.0)
    lo, hi = -40.0, db
    if moisture(db, hi, p) <= x:
        return db
    for _ in range(60):
        mid = (lo + hi) / 2
        if moisture(db, mid, p) < x:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def rel_humidity(db: float, x: float, p: float) -> float:
    pv = x * p / (0.622 + x)
    return min(100.0 * pv / p_sat(db), 100.0)


def temperatures(model: vsm.Model, flow: dict[str, float] | None = None,
                 extra_heat_kw: dict[str, float] | None = None,
                 surface_temp: float | None = None,
                 extra_latent_kw: dict[str, float] | None = None) -> dict[str, dict]:
    """Dry- and wet-bulb temperature in and out of every airway of the simulated stage.

    ``flow``: airflow per Unique No (m3/s, sign along end 1 -> end 2); saved airflow if
    omitted. ``surface_temp``: intake air temperature; by default the temperature Ventsim
    saved at each intake. Moisture: each airway adds the latent heat Ventsim saved
    ('Calc Latent', water evaporated = kW / 2450) plus ``extra_latent_kw``; moisture mixes
    by mass at junctions; wet bulb follows from dry bulb, moisture and pressure."""
    extra = extra_heat_kw or {}
    extra_l = extra_latent_kw or {}
    xlinks, xintake = [], {}
    net, N = _network(model, flow)
    links, info = [], []
    intake: dict[int, list[float]] = {}
    # In a heat simulation Ventsim's 'Calc Sensible' already contains the auto-compression
    # of the air (m*g*dz); split it off so it follows the new mass flow.
    heat_run = any(model.num(r, "Calc Sensible") for r, *_ in net)
    for r, uid, a, b, q in net:
        u, v = (a, b) if q >= 0 else (b, a)
        dz = model.num(r, "Z2") - model.num(r, "Z1")
        q_s = model.num(r, "Quantity")
        dz_saved = -dz if q_s < 0 else dz
        if q < 0:
            dz = -dz
        rho = _density(model, r)
        m = rho * abs(q)
        heat = 0.0
        if heat_run:
            m_s = rho * abs(q_s)
            if m_s >= LOW_FLOW:
                # saved sensible heat includes the auto-compression: take it out
                heat = model.num(r, "Calc Sensible") + m_s * G * dz_saved / 1000.0
            else:
                # at low flow Ventsim's saved heat and temperatures disagree; use the heat
                # that explains the saved temperature rise
                dT_s = model.num(r, "db air out") - model.num(r, "db air In")
                heat = m_s * (CP * dT_s + G * dz_saved) / 1000.0
        heat += extra.get(uid, 0.0)
        dT_heat = heat * 1000.0 / m / CP if m > MIN_MASSFLOW else 0.0
        # an airway whose airflow dropped far below the saved one would get its saved heat
        # in a trickle of air; cap the rise (the rock cannot heat air beyond its own
        # temperature anyway)
        dT_heat = max(min(dT_heat, MAX_HEAT_RISE), -MAX_HEAT_RISE)
        dT = dT_heat - G * dz / CP
        # moisture
        pin, pout = model.num(r, "p air In"), model.num(r, "p air out")
        if not 5e4 < pin < 2e5:
            pin = pout = 101325.0
        if q_s < 0:
            pin, pout = pout, pin
        if q * q_s < 0:  # flow reversed against the saved direction
            pin, pout = pout, pin
        xi_s = moisture(model.num(r, "db air In"), model.num(r, "wb air In"), pin)
        latent = 0.0
        if heat_run:
            m_s = rho * abs(q_s)
            if m_s >= LOW_FLOW:
                latent = model.num(r, "Calc Latent")
            else:
                xo_s = moisture(model.num(r, "db air out"), model.num(r, "wb air out"), pout)
                latent = m_s * (xo_s - xi_s) * HFG
        latent += extra_l.get(uid, 0.0)
        dx = latent / (m * HFG) if m > MIN_MASSFLOW else 0.0
        links.append((u, v, m, dT))
        xlinks.append((u, v, m, dx))
        info.append((uid, u, v, m, dT, heat, dx, latent, pin, pout))
        intake.setdefault(u, []).append(model.num(r, "db air In"))
        xintake.setdefault(u, []).append(xi_s)
    has_in = np.zeros(N, bool)
    for _u, v, m, _d in links:
        if m > 0:
            has_in[v] = True
    fixed = {}
    for k, ts in intake.items():
        if not has_in[k]:
            if surface_temp is not None:
                fixed[k] = float(surface_temp)
            else:
                ts = [t for t in ts if -60 < t < 80]
                fixed[k] = float(np.median(ts)) if ts else 20.0
    T = _mix(N, links, fixed, 20.0 if surface_temp is None else float(surface_temp))
    xfixed = {k: float(np.median(xs)) for k, xs in xintake.items() if not has_in[k] and xs}
    X = _mix(N, xlinks, xfixed, 0.01)
    out = {}
    for uid, u, v, m, dT, heat, dx, latent, pin, pout in info:
        ti, to = float(T[u]), float(T[u] + dT)
        xi, xo = float(X[u]), float(X[u] + dx)
        out[uid] = {"t_in": ti, "t_out": to, "heat_kw": heat, "latent_kw": latent,
                    "massflow": m, "x_out": xo,
                    "wb_in": wet_bulb(ti, xi, pin), "wb_out": wet_bulb(to, xo, pout),
                    "rh_out": rel_humidity(to, xo, pout)}
    return out


def concentrations(model: vsm.Model, sources_m3s: dict[str, float] | None = None,
                   flow: dict[str, float] | None = None, include_file_sources: bool = False
                   ) -> dict[str, dict]:
    """Steady-state volume fraction of a gas released in airways (m3/s of pure gas per
    Unique No), mixed at junctions: c_out = c_in + G/Q. Ventsim's own source input ('Base
    Injection' m3/s of a stream at 'Base Concentration' ppm) is added when
    ``include_file_sources``."""
    src = dict(sources_m3s or {})
    net, N = _network(model, flow)
    if include_file_sources:
        for r, uid, *_ in net:
            g = model.num(r, "Base Injection") * model.num(r, "Base Concentration") / 1e6
            if g > 0:
                src[uid] = src.get(uid, 0.0) + g
    links, info = [], []
    for r, uid, a, b, q in net:
        u, v = (a, b) if q >= 0 else (b, a)
        g, aq = src.get(uid, 0.0), abs(q)
        inc = g / max(aq, g) if g > 0 else 0.0
        links.append((u, v, aq, inc))
        info.append((uid, u, aq, inc))
    C = _mix(N, links, {}, 0.0)
    return {uid: {"c_in": float(C[u]), "c_out": min(float(C[u]) + inc, 1.0), "airflow": aq}
            for uid, u, aq, inc in info}


def arrival_times(model: vsm.Model, sources: list[str], flow: dict[str, float] | None = None,
                  start_s: float = 1.0) -> dict[str, dict]:
    """First arrival of a contaminant released at the start of the ``sources`` airways
    (time ``start_s``), carried at the mean air velocity (plug flow, no dispersion): the
    rule behind Ventsim's 'Spread Time' (arrival at the airway inlet) and 'End Time'
    (at its outlet). Airways the air never reaches from a source are left out."""
    import heapq

    net, N = _network(model, flow)
    succ: dict[int, list] = {}
    info = {}
    for r, uid, a, b, q in net:
        if abs(q) < 1e-9:
            continue
        u, v = (a, b) if q > 0 else (b, a)
        t = model.num(r, "Length") * model.num(r, "Area") / abs(q)
        succ.setdefault(u, []).append((v, t))
        info[uid] = (u, v, t)
    T: dict[int, float] = {}
    h = []
    for s in sources:
        if s not in info:
            raise ValueError(f"airway {s} carries no air (or is not in the simulated stage)")
        u = info[s][0]
        if start_s < T.get(u, float("inf")):
            T[u] = start_s
            h.append((start_s, u))
    heapq.heapify(h)
    while h:
        x, n = heapq.heappop(h)
        if x > T.get(n, float("inf")):
            continue
        for v, t in succ.get(n, []):
            if x + t < T.get(v, float("inf")):
                T[v] = x + t
                heapq.heappush(h, (x + t, v))
    return {uid: {"arrive_s": T[u], "leave_s": T[u] + t} for uid, (u, v, t) in info.items()
            if u in T}
