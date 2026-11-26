"""Fan upgrade NPV and DXF export."""
import pytest

from ventsim_mcp import dxf, finance, vsm

from test_vsm import model_path  # noqa: F401  (fixture)


def test_annuity_factor_matches_closed_form():
    assert finance.annuity_factor(0.0, 5) == 5.0
    assert abs(finance.annuity_factor(0.1, 10) - 6.1446) < 1e-4


def test_npv_positive_when_savings_beat_capex():
    r = finance.fan_upgrade_npv(100, 60, 8000, 0.1, 50000, 10, 0.1)
    # saving 40 kW * 8000 h = 320,000 kWh = 32,000 per year
    assert r["energy_saving_kwh_per_year"] == 320000.0
    assert r["npv"] > 0
    assert r["simple_payback_years"] == round(50000 / 32000, 2)


def test_npv_negative_when_no_saving():
    r = finance.fan_upgrade_npv(50, 80, 8000, 0.1, 1000, 5, 0.1)
    assert r["npv"] < 0 and r["simple_payback_years"] is None


def test_dxf_export_writes_lines(model_path, tmp_path):
    pytest.importorskip("ezdxf")
    m = vsm.parse(model_path)
    out = tmp_path / "a.dxf"
    n = dxf.export_airways_dxf(m, str(out))
    import ezdxf
    lines = [e for e in ezdxf.readfile(str(out)).modelspace() if e.dxftype() == "LINE"]
    assert n == len(lines) > 0
