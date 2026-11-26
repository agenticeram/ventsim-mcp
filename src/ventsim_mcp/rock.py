"""Rock-to-air heat exchange that responds to a change in airflow (experimental).

Ventsim's saved heat for an airway already contains the heat from the rock at the saved
airflow and air temperature. When the airflow or the air temperature changes, that rock
heat changes too. This module gives the change only:

    Q_rock = U' * L * (T_rock - T_air)          [kW, U' = U * Perimeter, W/(m K)]
    dQ     = Q_rock(new) - Q_rock(saved)

so an unchanged model reproduces the saved heat exactly. The rock temperature is the mean
of the airway's 'Rock VRT1/VRT2' if set, else the ``rock_temp`` given by the caller.

U per unit wall area combines the air film (Dittus-Boelter, Re from the airway's own
velocity and hydraulic diameter 4A/P, laminar floor Nu = 3.66) with the rock conduction
film k / sqrt(pi * alpha * t), where alpha is the rock's diffusivity and t the rock age
(years, at least 1 year: Ventsim's saved age is often 0). Latent (wet-rock) heat is not
modelled. These are standard correlations, not Ventsim's own rock model; confirm changes
in Ventsim.
"""
from __future__ import annotations

import math

from . import vsm

AIR_K = 0.026  # W/(m K)
AIR_MU = 1.8e-5  # Pa s
AIR_PR = 0.71
SECONDS_PER_YEAR = 3.156e7
MIN_AGE_YEARS = 1.0


def _rock_alpha(model: vsm.Model, r: list[str]) -> float:
    """Rock thermal diffusivity, m2/s: 'rock diffuse' (1e-6 m2/s) if saved, else k/(rho cp)."""
    diffuse = model.num(r, "rock diffuse")
    if diffuse > 0:
        return diffuse * 1e-6
    k = model.num(r, "rock conductivity", 2.0)
    rho = model.num(r, "rock density", 2700.0) or 2700.0
    cp = model.num(r, "rock spec heat", 790.0) or 790.0
    return k / (rho * cp)


def rock_temperature(model: vsm.Model, r: list[str], default: float | None) -> float | None:
    """Mean of the saved virgin rock temperatures, or ``default`` if the airway has none."""
    v1, v2 = model.num(r, "Rock VRT1"), model.num(r, "Rock VRT2")
    if v1 > 0 or v2 > 0:
        vals = [v for v in (v1, v2) if v > 0]
        return sum(vals) / len(vals)
    return default


def wall_conductance(model: vsm.Model, r: list[str], flow_m3s: float, rho: float) -> float:
    """U' = U * Perimeter, W/(m K), for an airway carrying ``flow_m3s``."""
    A = model.num(r, "Area")
    P = model.num(r, "Perimeter")
    if A <= 0 or P <= 0 or rho <= 0:
        return 0.0
    v = abs(flow_m3s) / A
    dh = 4.0 * A / P
    mu = AIR_MU
    re = rho * v * dh / mu
    nu = max(3.66, 0.023 * re ** 0.8 * AIR_PR ** 0.4) if re > 0 else 3.66
    h_air = nu * AIR_K / dh
    k = model.num(r, "rock conductivity", 2.0)
    t = max(model.num(r, "rock age"), MIN_AGE_YEARS) * SECONDS_PER_YEAR
    h_rock = k / math.sqrt(math.pi * _rock_alpha(model, r) * t)
    U = 1.0 / (1.0 / h_air + 1.0 / h_rock)
    return U * P


def rock_heat_change_kw(model: vsm.Model, r: list[str], flow_new: float, rho_new: float,
                        t_air_new: float, flow_saved: float, rho_saved: float,
                        t_air_saved: float, t_rock: float) -> float:
    """Change of the rock-to-air heat (kW) of one airway, against the saved state."""
    L = model.num(r, "Length")
    if L <= 0:
        return 0.0
    up_new = wall_conductance(model, r, flow_new, rho_new)
    up_saved = wall_conductance(model, r, flow_saved, rho_saved)
    q_new = up_new * L * (t_rock - t_air_new) / 1000.0
    q_saved = up_saved * L * (t_rock - t_air_saved) / 1000.0
    return q_new - q_saved
