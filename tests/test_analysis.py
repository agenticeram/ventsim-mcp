"""Analysis tools: shortest path, sensitivity, goal seek, radon. Use the synthetic model."""
import pytest

from ventsim_mcp import analysis, vsm

from test_vsm import model_path  # noqa: F401  (fixture)


def test_shortest_path_connects_the_ends(model_path):
    m = vsm.parse(model_path)
    p = analysis.shortest_path(m, "1", "3")
    assert p["airways"][0] == "1" and p["airways"][-1] == "3"


def test_sensitivity_reports_each_target(model_path):
    m = vsm.parse(model_path)
    s = analysis.sensitivity(m, targets=["1", "2"], runs=5)
    assert set(s["airways"]) == {"1", "2"}
    assert s["airways"]["1"]["std_m3s"] >= 0


def test_goal_seek_rejects_unreachable_target(model_path):
    m = vsm.parse(model_path)
    with pytest.raises(ValueError):
        analysis.goal_seek(m, "1", "2", 1e6)


def test_radon_zero_without_sources(model_path):
    m = vsm.parse(model_path)
    r = analysis.radon_steady(m)
    assert r["max_bq_m3"] == 0.0


def test_radon_rises_with_surface_source(model_path):
    m = vsm.parse(model_path)
    r = analysis.radon_steady(m, surface_bq_s_m2={"1": 1.0})
    assert r["airways"]["1"]["bq_m3"] > 0
