"""Export airways as DXF lines (experimental). Needs the optional 'ezdxf' package."""
from __future__ import annotations

from . import vsm


def export_airways_dxf(model: vsm.Model, out_path: str, layer: str = "VENTSIM_AIRWAYS") -> int:
    """Write every airway of the simulated stage as a LINE from its end 1 to end 2
    (X1,Y1,Z1 -> X2,Y2,Z2, metres). Returns the number of lines written."""
    try:
        import ezdxf
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("DXF export needs the optional package ezdxf (pip install ezdxf)") from e
    doc = ezdxf.new("R2010")
    doc.layers.add(layer)
    msp = doc.modelspace()
    n = 0
    mask = vsm.active_mask(model)
    for r, k in zip(model.rows(), mask):
        if not k:
            continue
        p1 = (model.num(r, "X1"), model.num(r, "Y1"), model.num(r, "Z1"))
        p2 = (model.num(r, "X2"), model.num(r, "Y2"), model.num(r, "Z2"))
        msp.add_line(p1, p2, dxfattribs={"layer": layer})
        n += 1
    doc.saveas(out_path)
    return n
