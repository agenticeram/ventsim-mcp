"""Fit airway resistances to measured airflows (ventilation survey calibration)."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from . import solver, vsm

REGULARISATION = 1e-3


def calibrate(model: vsm.Model, measured: dict[str, float], adjustable: list[str] | None = None,
              max_factor: float = 20.0, use_simulate: bool = False) -> dict:
    """Resistance factors for the ``adjustable`` airways (default: the measured ones) that
    bring the solved airflow closest to ``measured`` ({Unique No: m3/s}), in relative
    terms. Factors are limited to 1/max_factor .. max_factor. Solves with the what_if solver
    (saved state) or, with ``use_simulate``, from inputs only."""
    if not measured:
        raise ValueError("give at least one measured airflow")
    adj = [str(u) for u in (adjustable or list(measured))]
    meas = {str(k): float(v) for k, v in measured.items()}
    run = solver.simulate if use_simulate else solver.solve

    def flows(logs):
        s = run(model, resistance_scale={u: float(np.exp(x)) for u, x in zip(adj, logs)})
        return dict(zip(s["uid"], s["Q"]))

    base = flows(np.zeros(len(adj)))
    missing = [u for u in list(meas) + adj if u not in base]
    if missing:
        raise ValueError(f"airways not in the simulated stage: {missing[:10]}")
    keys = list(meas)
    scale = np.array([max(abs(meas[u]), 0.05 * max(abs(v) for v in meas.values()), 1e-3)
                      for u in keys])

    def resid(logs):
        q = flows(logs)
        r = np.array([abs(q[u]) - abs(meas[u]) for u in keys]) / scale
        # a light pull towards the model's own resistances picks the smallest change when
        # the measurements cannot tell airways apart (e.g. airways in series)
        return np.concatenate([r, REGULARISATION * np.asarray(logs)])

    lim = float(np.log(max_factor))
    r0 = resid(np.zeros(len(adj)))[:len(keys)]
    fit = least_squares(resid, np.zeros(len(adj)), bounds=(-lim, lim), diff_step=1e-4,
                        xtol=1e-6, ftol=1e-10, max_nfev=60)
    q = flows(fit.x)
    return {
        "factors": {u: float(np.exp(x)) for u, x in zip(adj, fit.x)},
        "before": {u: float(base[u]) for u in keys},
        "after": {u: float(q[u]) for u in keys},
        "rms_rel_error_before": float(np.sqrt(np.mean(r0 ** 2))),
        "rms_rel_error_after": float(np.sqrt(np.mean(fit.fun[:len(keys)] ** 2))),
        "at_limit": [u for u, x in zip(adj, fit.x) if abs(abs(x) - lim) < 1e-3],
        "flows": q,
    }
