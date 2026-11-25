"""Tests run on a small synthetic .vsm built in-memory (no Ventsim files shipped)."""
import gzip

import pytest

from ventsim_mcp import server, vsm

COLS = ["MAIN", "Unique No", "Warning Key", "GUID", "Branch Name", "Entry Node", "Exit Node",
        "Notes", "Group", "X1", "Y1", "Z1", "X2", "Y2", "Z2", "Width", "Height", "Area",
        "Perimeter", "Length", "kFactor", "Quantity", "Fan Type", "Fan Numbers", "Stage",
        "Resistance Total", "Pressure From", "Pressure To", "QuantityFix"]


def _row(i, name, p1, p2, q, fan="0"):
    d = dict.fromkeys(COLS, "")
    d.update({"MAIN": str(i), "Unique No": str(i), "Branch Name": name, "Notes": "", "Group": "0",
              "X1": p1[0], "Y1": p1[1], "Z1": p1[2], "X2": p2[0], "Y2": p2[1], "Z2": p2[2],
              "Width": "5", "Height": "5", "Area": "25", "Perimeter": "20", "Length": "100",
              "kFactor": "0.012", "Quantity": q, "Fan Type": fan, "Fan Numbers": "1", "Stage": "1"})
    return "\t".join(d[c] for c in COLS)


@pytest.fixture()
def model_path(tmp_path):
    lines = ["4.0\tTest\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t 6.0.0.0", "PRESETS", "END\tPRESETS",
             "\t".join(COLS),
             _row(1, "Decline", ("0", "0", "0"), ("100", "0", "0"), "50", fan="1"),
             _row(2, "Left", ("100", "0", "0"), ("200", "0", "0"), "30"),
             _row(3, "Right", ("100", "0", "0"), ("200", "50", "0"), "30"),  # 60 out, 50 in
             "FANS", "1\tTest Fan\t2.8", "0\t2100\t0", "END\tFANS",
             "ROLLBACK\tx", "ENDFILE"]
    p = tmp_path / "t.vsm"
    p.write_bytes(gzip.compress("\r\n".join(lines).encode("latin-1")))
    return str(p)


def test_summary(model_path):
    s = vsm.summary(vsm.parse(model_path))
    assert s["airways"] == 3 and s["fans"] == 1 and s["total_length_m"] == 300.0
    assert s["ventsim_version"] == "6.0.0.0"


def test_fan_name(model_path):
    assert vsm.fans(vsm.parse(model_path))[0]["Fan Name"] == "Test Fan"


def test_roundtrip(model_path, tmp_path):
    out = tmp_path / "o.vsm"
    vsm.dump(vsm.parse(model_path), str(out))
    assert gzip.open(out).read() == gzip.open(model_path).read()


def test_balance(model_path):
    r = server.check_flow_balance(model_path)
    assert r["unbalanced"] == 1 and r["worst"][0]["imbalance_m3s"] == -10.0


def test_edit_never_overwrites(model_path, tmp_path):
    out = str(tmp_path / "e.vsm")
    server.edit_airways(model_path, out, [{"unique_no": "2", "field": "Notes", "value": "hi"}])
    assert server.get_airway(out, "2")["Notes"] == "hi"
    with pytest.raises(ValueError):
        server.edit_airways(model_path, out, [])
    with pytest.raises(ValueError):
        server.edit_airways(model_path, str(tmp_path / "x.vsm"),
                            [{"unique_no": "2", "field": "Width", "value": "9"}])


def test_velocity_and_fans_and_compare(model_path, tmp_path):
    r = server.velocity_check(model_path, max_velocity=1.5)
    assert r["over_max"] == 1 and r["highest"][0]["Unique No"] == "1"  # 50/25 = 2 m/s
    assert server.fan_summary(model_path)["by_model"]["Test Fan"]["fans"] == 1
    other = str(tmp_path / "b.vsm")
    server.edit_airways(model_path, other, [{"unique_no": "2", "field": "Notes", "value": "x"}])
    assert server.compare_models(model_path, other)["changed"] == 0


def test_solver_series_what_if(tmp_path):
    lines = ["4.0\tT\t1 6.0.0.0", "\t".join(COLS)]
    specs = [(0, 1, "1000", "0"), (1, 2, "0", "0")]  # node 0 held at 1000 Pa, node 2 at 0
    for i, (a, b, pf, pt) in enumerate(specs, start=1):
        d = dict.fromkeys(COLS, "")
        d.update({"MAIN": str(i), "Unique No": str(i), "X1": str(a), "Y1": "0", "Z1": "0",
                  "X2": str(b), "Y2": "0", "Z2": "0", "Quantity": "70.7",
                  "Resistance Total": "0.1", "Pressure From": pf, "Pressure To": pt,
                  "Fan Type": "0", "Fan Numbers": "1"})
        lines.append("\t".join(d[c] for c in COLS))
    p = tmp_path / "s.vsm"
    p.write_bytes(gzip.compress("\r\n".join(lines + ["ENDFILE"]).encode("latin-1")))
    r = server.what_if(str(p), resistance_scale={"2": 2.0})
    assert r["baseline_fit"]["mean_rel_error"] < 0.01
    # series R 0.2 -> sqrt(1000/0.2)=70.7; after scaling R = 0.3 -> 57.7
    assert abs(r["largest_changes"][0]["delta"] + 13.0) < 0.1


def test_solver_fixed_flow_stays_fixed(tmp_path):
    lines = ["4.0\tT\t1 6.0.0.0", "\t".join(COLS)]
    for i, (a, b, pf, pt) in enumerate([(0, 1, "1000", "0"), (1, 2, "0", "0")], start=1):
        d = dict.fromkeys(COLS, "")
        d.update({"MAIN": str(i), "Unique No": str(i), "X1": str(a), "Y1": "0", "Z1": "0",
                  "X2": str(b), "Y2": "0", "Z2": "0", "Quantity": "50", "QuantityFix": "50",
                  "Resistance Total": "0.1", "Pressure From": pf, "Pressure To": pt,
                  "Fan Type": "0", "Fan Numbers": "1"})
        lines.append("\t".join(d[c] for c in COLS))
    p = tmp_path / "f.vsm"
    p.write_bytes(gzip.compress("\r\n".join(lines + ["ENDFILE"]).encode("latin-1")))
    r = server.what_if(str(p), resistance_scale={"2": 5.0})
    assert r["baseline_fit"]["mean_rel_error"] == 0 and r["airways_changed"] == 0


def _stage_model(tmp_path):
    """Y-junction where airway 3 exists twice: stage bit 0 (flow 10) and bit 1 (flow 30).
    Saved flows balance only for stage 1 (bit 1): 50 in = 20 + 30 out."""
    lines = ["4.0\tT\t1 6.0.0.0", "\t".join(COLS)]
    spec = [("1", (0, 0), (1, 0), "50", "3"), ("2", (1, 0), (2, 0), "20", "3"),
            ("3", (1, 0), (2, 1), "10", "1"), ("4", (1, 0), (2, 1), "30", "2")]
    for uid, p1, p2, q, stage in spec:
        d = dict.fromkeys(COLS, "")
        d.update({"MAIN": uid, "Unique No": uid, "X1": str(p1[0]), "Y1": str(p1[1]), "Z1": "0",
                  "X2": str(p2[0]), "Y2": str(p2[1]), "Z2": "0", "Quantity": q, "Stage": stage,
                  "Fan Type": "0", "Fan Numbers": "1", "Resistance Total": "0.01"})
        lines.append("\t".join(d[c] for c in COLS))
    p = tmp_path / "st.vsm"
    p.write_bytes(gzip.compress("\r\n".join(lines + ["ENDFILE"]).encode("latin-1")))
    return str(p)


def test_active_stage_picks_balanced_stage(tmp_path):
    m = vsm.parse(_stage_model(tmp_path))
    assert vsm.active_stage(m) == 1
    assert vsm.active_mask(m) == [True, True, False, True]
    assert server.check_flow_balance(str(m.path))["unbalanced"] == 0


def test_close_points_are_one_node(tmp_path):
    lines = ["4.0\tT\t1 6.0.0.0", "\t".join(COLS)]
    for uid, x1, x2 in (("1", "0", "10"), ("2", "10.007", "20")):
        d = dict.fromkeys(COLS, "")
        d.update({"MAIN": uid, "Unique No": uid, "X1": x1, "Y1": "0", "Z1": "0", "X2": x2,
                  "Y2": "0", "Z2": "0", "Quantity": "5", "Fan Type": "0"})
        lines.append("\t".join(d[c] for c in COLS))
    p = tmp_path / "n.vsm"
    p.write_bytes(gzip.compress("\r\n".join(lines + ["ENDFILE"]).encode("latin-1")))
    ends = vsm.node_ids(vsm.parse(str(p)))
    assert ends[0][1] == ends[1][0]
