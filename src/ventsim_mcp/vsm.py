"""Reader/writer for Ventsim ``.vsm`` model files.

A ``.vsm`` file is a gzip-compressed, CRLF separated, tab delimited text file.
The airway table starts at a header line whose first cell is ``MAIN``; the
header names every column, so the parser is driven by those names instead of
hard-coded positions.  Values are kept as the original strings, which makes
``parse -> dump`` byte-exact when nothing was edited.

Ventsim stores everything in SI units (m, m3/s, Pa, kW) regardless of the
display units chosen in the program.
"""
from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass, field
from typing import Iterable

EOL = "\r\n"
ENCODING = "latin-1"  # byte-transparent, so unknown text survives a round trip


@dataclass
class Model:
    path: str
    lines: list[str]
    header_index: int
    columns: list[str]
    row_indexes: list[int]  # line numbers of the airway rows
    _col: dict[str, int] = field(default_factory=dict, repr=False)

    # ---- basic access -------------------------------------------------
    def col(self, name: str) -> int:
        try:
            return self._col[name]
        except KeyError:
            raise KeyError(f"column not in file: {name!r}") from None

    def row(self, n: int) -> list[str]:
        return self.lines[self.row_indexes[n]].split("\t")

    def rows(self) -> Iterable[list[str]]:
        for i in self.row_indexes:
            yield self.lines[i].split("\t")

    def value(self, row: list[str], name: str) -> str:
        if name == "Unique No" and "Unique No" not in self._col:
            return row[0]  # older files: the running row number identifies the airway
        i = self._col.get(name)
        return row[i] if i is not None and i < len(row) else ""

    def num(self, row: list[str], name: str, default: float = 0.0) -> float:
        try:
            return float(self.value(row, name))
        except ValueError:
            return default

    # ---- sections -----------------------------------------------------
    def section(self, name: str) -> list[list[str]]:
        """Rows of a simple ``NAME ... END<TAB>NAME`` block, split on tabs."""
        out: list[list[str]] = []
        inside = False
        for line in self.lines:
            if not inside:
                if line == name:
                    inside = True
                continue
            if line.startswith("END\t") and line.split("\t")[1] == name:
                break
            out.append(line.split("\t"))
        return out

    @property
    def version(self) -> str:
        for cell in self.lines[0].split("	"):
            if re.fullmatch(r"\s*\d+\.\d+\.\d+\.\d+\s*", cell):
                return cell.strip()
        return ""

    @property
    def title(self) -> str:
        head = self.lines[0].split("\t")
        return head[1] if len(head) > 1 else ""


def read_lines(path: str) -> list[str]:
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode(ENCODING).split(EOL)


def parse(path: str) -> Model:
    path = os.path.abspath(path)
    lines = read_lines(path)
    header = next(
        (i for i, l in enumerate(lines) if l.startswith("MAIN\t") and l.count("\t") > 10), None
    )
    if header is None:
        raise ValueError("not a Ventsim model: airway table header not found")
    columns = lines[header].split("\t")
    width = len(columns)
    rows: list[int] = []
    for i in range(header + 1, len(lines)):
        line = lines[i]
        cells = line.split("\t")
        # airway rows are numeric-indexed and as wide as the header (±1 for
        # a trailing empty cell); the table ends at the first line that is not.
        if len(cells) >= width - 2 and cells[0].lstrip("-").isdigit():
            rows.append(i)
        else:
            break
    m = Model(path, lines, header, columns, rows)
    m._col = {name: i for i, name in enumerate(columns) if name}
    return m


def dump(model: Model, path: str) -> None:
    data = EOL.join(model.lines).encode(ENCODING)
    with open(path, "wb") as f:
        f.write(gzip.compress(data, mtime=0))


# ---- airway view ---------------------------------------------------------
BASIC_FIELDS = (
    "Unique No", "Branch Name", "Notes", "Group", "Stage",
    "X1", "Y1", "Z1", "X2", "Y2", "Z2",
    "Width", "Height", "Area", "Perimeter", "Length", "kFactor",
    "Quantity", "Pressure Branch", "Resistance Total",
    "Fan Type", "Fan Numbers", "Pressure Fan", "Static Pressure", "Fan Power",
)


def _is_number(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def fan_types(model: Model) -> dict[str, str]:
    """Fan type id -> name, from the FANS section (header rows carry a text name)."""
    out: dict[str, str] = {}
    for r in model.section("FANS"):
        if len(r) > 2 and r[1] and not _is_number(r[1]):
            out[r[0]] = r[1]
    return out


def air_types(model: Model) -> dict[str, str]:
    return {r[0]: r[1] for r in model.section("AIRTYPES") if len(r) > 1}


def airway(model: Model, row: list[str]) -> dict:
    d: dict = {}
    for name in BASIC_FIELDS:
        v = model.value(row, name)
        if name in ("Unique No", "Branch Name", "Notes", "Group", "Stage", "Fan Type"):
            d[name] = v
        else:
            try:
                d[name] = round(float(v), 4)
            except ValueError:
                d[name] = None
    area = d.get("Area") or 0
    q = d.get("Quantity") or 0
    names = model.__dict__.get("_fan_names")
    if names is None:  # scanning the FANS section per airway made large models slow
        names = model.__dict__["_fan_names"] = fan_types(model)
    d["Fan Name"] = names.get(d["Fan Type"], "") if d["Fan Type"] not in ("", "0") else ""
    d["Velocity m/s"] = round(abs(q) / area, 3) if area else None
    return d


def airways(model: Model) -> list[dict]:
    out = []
    for r, active in zip(model.rows(), active_mask(model)):
        d = airway(model, r)
        d["In Simulated Stage"] = active
        out.append(d)
    return out


def fans(model: Model) -> list[dict]:
    out = []
    for r in model.rows():
        ft = model.value(r, "Fan Type")
        if ft not in ("", "0"):
            out.append(airway(model, r))
    return out


def summary(model: Model) -> dict:
    ws = airways(model)
    note = None
    if not ws:
        ref = sum(1 for l in model.lines if l[:1].isdigit()) if "REFERENCE" in model.lines else 0
        note = ("the airway table is empty" + (
            "; the file holds only reference graphics (e.g. imported CAD lines) that Ventsim "
            "has not converted to airways" if ref else ""))
    total_len = sum(w["Length"] or 0 for w in ws)
    return {
        "file": model.path,
        "title": model.title,
        "ventsim_version": model.version,
        "airways": len(ws),
        "fans": len(fans(model)),
        "total_length_m": round(total_len, 1),
        "z_range_m": [
            round(min(min(w["Z1"] or 0, w["Z2"] or 0) for w in ws), 1),
            round(max(max(w["Z1"] or 0, w["Z2"] or 0) for w in ws), 1),
        ] if ws else None,
        "has_simulation_results": any((w["Quantity"] or 0) for w in ws),
        "simulated_stage": stage_name(model, active_stage(model)) or None,
        "airways_in_simulated_stage": sum(active_mask(model)),
        "units": "SI (m, m3/s, Pa, kW) as stored in file",
        **({"note": note} if note else {}),
    }


def node_ids(model: Model, tolerance: float = 0.05) -> list[tuple[int, int]]:
    """(entry node, exit node) index per airway. Ends closer than ``tolerance`` metres are
    the same junction (Ventsim stores coordinates with small rounding differences)."""
    pts = []
    for r in model.rows():
        for s in "12":
            pts.append(tuple(model.num(r, f"{c}{s}") for c in "XYZ"))
    parent = list(range(len(pts)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    cell: dict[tuple, list[int]] = {}
    for k, p in enumerate(pts):
        g = tuple(int(v // tolerance) for v in p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for o in cell.get((g[0] + dx, g[1] + dy, g[2] + dz), ()):
                        q = pts[o]
                        if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 <= tolerance ** 2:
                            parent[find(k)] = find(o)
        cell.setdefault(g, []).append(k)
    roots: dict[int, int] = {}
    ids = [roots.setdefault(find(k), len(roots)) for k in range(len(pts))]
    return [(ids[2 * n], ids[2 * n + 1]) for n in range(len(pts) // 2)]


def _stage_bits(model: Model) -> list[int]:
    out = []
    for r in model.rows():
        try:
            out.append(int(float(model.value(r, "Stage") or 1)))
        except ValueError:
            out.append(1)
    return out


def active_stage(model: Model) -> int | None:
    """Index of the stage whose results are saved, or None for a single-stage model.

    'Stage' is a bit mask: an airway exists in stage k when bit k is set, and an airway
    that changes between stages is stored once per version, on top of each other. The file
    does not mark the stage that was simulated last, so pick the stage whose airways give a
    mass balance of the saved airflow (ties: the stage with more airways)."""
    cached = model.__dict__.get("_active_stage", "unset")
    if cached != "unset":
        return cached
    bits = _stage_bits(model)
    present = sorted({k for b in bits for k in range(31) if b >> k & 1})
    result = None
    if len(present) > 1 and len(set(bits)) > 1:
        rows = list(model.rows())
        ends = node_ids(model)
        best = None
        for k in present:
            net: dict[int, float] = {}
            deg: dict[int, int] = {}
            count = 0
            for r, (a, b), m in zip(rows, ends, bits):
                if not m >> k & 1:
                    continue
                count += 1
                q = model.num(r, "Dry Massflow") or model.num(r, "Quantity")
                if model.num(r, "Quantity") < 0:
                    q = -abs(q)
                net[a] = net.get(a, 0.0) - q
                net[b] = net.get(b, 0.0) + q
                deg[a] = deg.get(a, 0) + 1
                deg[b] = deg.get(b, 0) + 1
            bad = sum(abs(v) for n, v in net.items() if deg[n] > 1)
            score = (round(bad, 3), -count)
            if best is None or score < best[0]:
                best = (score, k)
        result = best[1]
    model.__dict__["_active_stage"] = result
    return result


def active_mask(model: Model) -> list[bool]:
    """True for airways that belong to the simulated stage."""
    k = active_stage(model)
    bits = _stage_bits(model)
    return [True] * len(bits) if k is None else [bool(b >> k & 1) for b in bits]


def stage_name(model: Model, k: int | None) -> str:
    if k is None:
        return ""
    for r in model.section("STAGES"):
        if r and r[0] == str(k) and len(r) > 1:
            return r[1]
    return str(k)
