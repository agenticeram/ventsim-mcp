"""Rock-heat response: unchanged model is untouched; less airflow in a hot airway gives the
air more time with the rock and so cools it. Uses the synthetic model from test_vsm."""
from ventsim_mcp import rock, thermo, vsm
from test_vsm import model_path  # noqa: F401  (fixture)


def test_wall_conductance_grows_with_velocity(model_path):
    m = vsm.parse(model_path)
    r = list(m.rows())[0]
    slow = rock.wall_conductance(m, r, 0.5, 1.2)
    fast = rock.wall_conductance(m, r, 5.0, 1.2)
    assert 0 < slow < fast


def test_unchanged_model_identical_with_rock(model_path):
    m = vsm.parse(model_path)
    base = thermo.temperatures(m)
    with_rock = thermo.temperatures_with_rock(m, rock_temp=30.0)
    for uid in base:
        assert abs(base[uid]["t_out"] - with_rock[uid]["t_out"]) < 1e-9


def test_less_flow_means_less_rock_cooling(model_path):
    # saved state: 5 m3/s cooled by a rock at -10 C; new state: 1 m3/s. The rock then cools
    # less than it did, so the heat change against the saved state is positive (warmer).
    m = vsm.parse(model_path)
    r = list(m.rows())[0]
    d = rock.rock_heat_change_kw(m, r, 1.0, 1.2, 0.0, 5.0, 1.2, 0.0, -10.0)
    assert d > 0
    # and more flow cools more: the change is negative
    assert rock.rock_heat_change_kw(m, r, 5.0, 1.2, 0.0, 1.0, 1.2, 0.0, -10.0) < 0
