"""Geometry edits on an in-memory model: airway size, length, friction factor, new airways.

Resistance is rescaled with the Atkinson relation R ~ k * P * L / A**3 relative to the stored
value, so whatever Ventsim folded into it (air density, shock losses entered as a factor) is
kept in proportion. Saved results of edited airways become stale; simulate_airflow does not
use them.
"""
from __future__ import annotations

import math
import uuid

from . import vsm

GEOMETRY_FIELDS = ("Width", "Height", "Length", "kFactor")
# Branch Flags bits a new airway must not copy: fixed flow (2, 4), surface ends (12, 13),
# dead end (18), fixed pressure (19)
CLEARED_FLAGS = (1 << 2) | (1 << 4) | (1 << 12) | (1 << 13) | (1 << 18) | (1 << 19)
SURFACE_FLAGS = (1 << 12, 1 << 13)
DEAD_END_FLAG = 1 << 18
JOIN_DISTANCE = 0.05  # m, as vsm.node_ids
# template inputs a new airway does not take over
NOT_COPIED = ("Fan Type", "QuantityFix", "Pressure Fix", "Heat Type", "heatpointkw",
              "heatlinearkw", "dieselkw", "latentheatkw", "Base Injection", "Base Concentration",
              "Gas Source", "Dust Source")
# saved results cleared on a new airway (it has not been simulated)
RESULT_FIELDS = (
    "Quantity", "Pressure Branch", "Pressure From", "Pressure To", "PressureRes",
    "Calc Sim Pressure", "Calc NVP Pressure", "Dry Massflow", "Calc Sensible", "Calc Latent",
    "Calc Condensate", "Pressure Fan", "Static Pressure", "Fan Power", "Sim Concentration",
    "Sim Concentration Out",
)


def _fmt(v: float) -> str:
    return f"{v:.7g}"


def _shape(model: vsm.Model, cells: list[str], w: float, h: float) -> tuple[float, float]:
    """Area and perimeter for width w and height h, keeping the airway's shape."""
    shape = model.value(cells, "Shape Type")
    if shape == "0":  # rectangle
        return w * h, 2 * (w + h)
    if shape == "1":  # round: width is the diameter
        return math.pi * w * w / 4, math.pi * w
    # other profiles (arched, custom): scale the stored area and perimeter
    w0, h0 = model.num(cells, "Width"), model.num(cells, "Height")
    a0, p0 = model.num(cells, "Area"), model.num(cells, "Perimeter")
    if w0 > 0 and h0 > 0 and a0 > 0 and p0 > 0:
        return a0 * (w * h) / (w0 * h0), p0 * (w + h) / (w0 + h0)
    return w * h, 2 * (w + h)


def _set(model: vsm.Model, cells: list[str], name: str, value: str) -> None:
    i = model._col.get(name)
    if i is not None and i < len(cells):
        cells[i] = value


def set_geometry(model: vsm.Model, cells: list[str], **new: float) -> dict:
    """Change Width / Height / Length / kFactor of one airway row (in place) and rescale its
    area, perimeter and resistance. Returns the old and new resistance."""
    for k, v in new.items():
        if k not in GEOMETRY_FIELDS:
            raise ValueError(f"not a geometry field: {k}")
        if not (isinstance(v, (int, float)) and v > 0 and math.isfinite(v)):
            raise ValueError(f"{k} must be a positive number")
    w0, h0 = model.num(cells, "Width"), model.num(cells, "Height")
    a0, p0 = model.num(cells, "Area"), model.num(cells, "Perimeter")
    l0, k0 = model.num(cells, "Length"), model.num(cells, "kFactor")
    w, h = new.get("Width", w0), new.get("Height", h0)
    L, k = new.get("Length", l0), new.get("kFactor", k0)
    if "Width" in new or "Height" in new:
        a, p = _shape(model, cells, w, h)
    else:
        a, p = a0, p0
    r0 = model.num(cells, "Resistance Total")
    r = r0
    closed = model.value(cells, "Resistance Type") == "1"
    if not closed and min(a0, p0, l0, k0) > 0:
        ratio = (k * p * L / a ** 3) / (k0 * p0 * l0 / a0 ** 3)
        r = r0 * ratio
        _set(model, cells, "Resistance Total", _fmt(r))
        if model.num(cells, "Resistance") > 0:
            _set(model, cells, "Resistance", _fmt(model.num(cells, "Resistance") * ratio))
    for name, v in (("Width", w), ("Height", h), ("Area", a), ("Perimeter", p),
                    ("Length", L), ("kFactor", k)):
        if name in new or name in ("Area", "Perimeter") and ("Width" in new or "Height" in new):
            _set(model, cells, name, _fmt(v))
    # Ventsim recomputes an AUTO resistance (type 0) from size, friction preset and the
    # length between the end coordinates when it simulates. A friction factor is kept only
    # as a custom value (kFactor type -1), and a length that differs from the coordinates
    # only through a custom resistance (type -1). Checked in Ventsim 6.0.5.
    if "kFactor" in new:
        _set(model, cells, "kFactor type", "-1")
    if "Length" in new and not closed:
        set_custom_resistance(model, cells, r)
    return {"resistance_before": r0, "resistance_after": r, "closed": closed}


def set_custom_resistance(model: vsm.Model, cells: list[str], r: float) -> None:
    """Store ``r`` as the airway's custom resistance (Ventsim keeps a custom value when it
    simulates; an AUTO one it recomputes from the airway size and friction)."""
    _set(model, cells, "Resistance Type", "-1")
    _set(model, cells, "Resistance", _fmt(r))
    _set(model, cells, "Resistance Total", _fmt(r))


def _existing_ends(model: vsm.Model) -> list[tuple[tuple[float, float, float], bool, str]]:
    """Ends of the simulated-stage airways: (point, is an open loose end, Unique No)."""
    mask = vsm.active_mask(model)
    rows = [r for r, k in zip(model.rows(), mask) if k]
    ids = [e for e, k in zip(vsm.node_ids(model), mask) if k]
    deg: dict[int, int] = {}
    for a, b in ids:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    out = []
    for r, e in zip(rows, ids):
        try:
            flags = int(float(model.value(r, "Branch Flags") or 0))
        except ValueError:
            flags = 0
        for s_, c in enumerate("12"):
            pt = tuple(model.num(r, x + c) for x in "XYZ")
            loose = deg[e[s_]] == 1
            is_open = loose and bool(flags & SURFACE_FLAGS[s_] or not flags & DEAD_END_FLAG)
            out.append((pt, is_open, model.value(r, "Unique No")))
    return out


def add_airway(model: vsm.Model, like: str, start: list[float], end: list[float],
               name: str = "", width: float | None = None, height: float | None = None,
               open_ends: tuple[str, ...] = ()) -> str:
    """Append a new airway copied from airway ``like`` (size, shape, friction, stage) between
    two points. An end within 5 cm of an existing airway end joins it there; joining an end
    that is open to the surface is refused (it would close that portal). An end that joins
    nothing is a dead end unless named in ``open_ends`` ("start" / "end"), which makes it
    open to the surface. Fans, fixed flows and pressures, heat and gas sources of the
    template are not copied. Returns the new Unique No."""
    if any(e not in ("start", "end") for e in open_ends):
        raise ValueError("open_ends may only contain 'start' and 'end'")
    if len(start) != 3 or len(end) != 3:
        raise ValueError("start and end must be [x, y, z]")
    pts = [float(v) for v in (*start, *end)]
    if not all(math.isfinite(v) for v in pts):
        raise ValueError("coordinates must be finite numbers")
    L = math.dist(pts[:3], pts[3:])
    if L < 0.1:
        raise ValueError("airway must be at least 0.1 m long")
    index = {model.value(model.row(n), "Unique No"): n for n in range(len(model.row_indexes))}
    if str(like) not in index:
        raise ValueError(f"no airway with Unique No {like}")
    cells = model.row(index[str(like)])
    uids = [int(u) for u in index if u.lstrip("-").isdigit()]
    uid = str(max(uids, default=0) + 1)
    idx = max((int(model.row(n)[0]) for n in range(len(model.row_indexes))), default=0) + 1
    ends = _existing_ends(model)
    joined = []
    for which, p in (("start", pts[:3]), ("end", pts[3:])):
        near = [(o, u) for q, o, u in ends if math.dist(p, q) <= JOIN_DISTANCE]
        if near and which in open_ends:
            raise ValueError(f"{which} joins existing airways; it cannot also be open to surface")
        if any(o for o, _u in near):
            raise ValueError(f"{which} point is the surface connection of airway "
                             f"{[u for o, u in near if o][0]}; joining there would close it")
        joined.append(bool(near))
    cells = list(cells)
    cells[0] = str(idx)
    _set(model, cells, "Unique No", uid)
    _set(model, cells, "GUID", str(uuid.uuid4()))
    _set(model, cells, "Branch Name", name.replace("\t", " "))
    for c, v in zip(("X1", "Y1", "Z1", "X2", "Y2", "Z2"), pts):
        _set(model, cells, c, _fmt(v))
    for c in ("Entry Node", "Exit Node"):
        _set(model, cells, c, "")
    # do not inherit end conditions or fixed flow / pressure from the template
    try:
        flags = int(float(model.value(cells, "Branch Flags") or 0)) & ~CLEARED_FLAGS
    except ValueError:
        flags = 0
    for s_, which in enumerate(("start", "end")):
        if which in open_ends:
            flags |= SURFACE_FLAGS[s_]
    if not all(j or w in open_ends for j, w in zip(joined, ("start", "end"))):
        flags |= DEAD_END_FLAG
    _set(model, cells, "Branch Flags", str(flags))
    for c in NOT_COPIED:
        if c in model._col and model.value(cells, c) not in ("", "0"):
            _set(model, cells, c, "0")
    for c in RESULT_FIELDS:
        if c in model._col and model.value(cells, c) != "":
            _set(model, cells, c, "0")
    geo = {"Length": L}
    if width:
        geo["Width"] = float(width)
    if height:
        geo["Height"] = float(height)
    set_geometry(model, cells, **geo)
    at = model.row_indexes[-1] + 1
    model.lines.insert(at, "\t".join(cells))
    model.row_indexes.append(at)
    model.__dict__.pop("_fan_names", None)
    return uid
