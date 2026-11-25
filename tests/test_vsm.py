"""Tests run on a small synthetic .vsm built in-memory (no Ventsim files shipped)."""
import gzip

import pytest

from ventsim_mcp import server, vsm

COLS = ["MAIN", "Unique No", "Warning Key", "GUID", "Branch Name", "Entry Node", "Exit Node",
        "Notes", "Group", "X1", "Y1", "Z1", "X2", "Y2", "Z2", "Width", "Height", "Area",
        "Perimeter", "Length", "kFactor", "Quantity", "Fan Type", "Fan Numbers", "Stage"]


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
