"""Fan upgrade economics: yearly energy saving and net present value (experimental).

A plain discounted-cash-flow calculator for a change of fan power, not Ventsim's costing
model (which also prices shafts, excavation and fixed costs)."""
from __future__ import annotations


def annuity_factor(rate: float, years: int) -> float:
    """Present value of 1 per year for ``years`` at discount ``rate`` (0.1 = 10 %)."""
    if years <= 0:
        return 0.0
    if rate == 0:
        return float(years)
    return (1.0 - (1.0 + rate) ** -years) / rate


def fan_upgrade_npv(power_now_kw: float, power_new_kw: float, hours_per_year: float,
                    price_per_kwh: float, capex: float, life_years: int,
                    discount_rate: float) -> dict:
    """Net present value of replacing a fan (or changing its speed) that moves from
    ``power_now_kw`` to ``power_new_kw``: yearly energy cost saving over ``life_years``,
    discounted at ``discount_rate``, minus the up-front ``capex``. Payback is simple
    (capex / yearly saving)."""
    if min(power_now_kw, power_new_kw, hours_per_year, price_per_kwh, capex) < 0:
        raise ValueError("power, hours, price and capex must not be negative")
    if life_years < 1:
        raise ValueError("life_years must be at least 1")
    saving_kwh = (power_now_kw - power_new_kw) * hours_per_year
    saving = saving_kwh * price_per_kwh
    pv = saving * annuity_factor(discount_rate, life_years)
    npv = pv - capex
    payback = capex / saving if saving > 0 else None
    return {"energy_saving_kwh_per_year": round(saving_kwh, 1),
            "cost_saving_per_year": round(saving, 2),
            "present_value_of_savings": round(pv, 2),
            "capex": capex, "npv": round(npv, 2),
            "simple_payback_years": round(payback, 2) if payback else None,
            "note": "constant power, price and rate; no fan maintenance, escalation or tax; "
                    "not Ventsim's costing model"}
