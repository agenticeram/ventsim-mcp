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
        (i for i, l in enumerate(lines) if l.startswith("MAIN\tUnique No\t")), None
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
    d["Velocity m/s"] = round(abs(q) / area, 3) if area else None
    return d


def airways(model: Model) -> list[dict]:
    return [airway(model, r) for r in model.rows()]


def fans(model: Model) -> list[dict]:
    out = []
    for r in model.rows():
        ft = model.value(r, "Fan Type")
        if ft not in ("", "0"):
            out.append(airway(model, r))
    return out


def summary(model: Model) -> dict:
    ws = airways(model)
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
        "units": "SI (m, m3/s, Pa, kW) as stored in file",
    }
