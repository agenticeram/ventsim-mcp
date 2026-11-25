"""MCP server exposing Ventsim .vsm model files to an AI assistant."""
from __future__ import annotations

import csv
import io
import os
from collections import defaultdict

from mcp.server.fastmcp import FastMCP

from . import vsm

mcp = FastMCP("ventsim")

_cache: dict[str, tuple[float, vsm.Model]] = {}


def _load(path: str) -> vsm.Model:
    path = os.path.abspath(os.path.expanduser(path))
    mtime = os.path.getmtime(path)
    hit = _cache.get(path)
    if hit and hit[0] == mtime:
        return hit[1]
    model = vsm.parse(path)
    _cache[path] = (mtime, model)
    return model


@mcp.tool()
def model_summary(path: str) -> dict:
    """Summarise a Ventsim .vsm file: title, version, airway/fan counts, total length."""
    return vsm.summary(_load(path))


@mcp.tool()
def list_airways(
    path: str,
    name_contains: str = "",
    min_velocity: float | None = None,
    sort_by: str = "Unique No",
    descending: bool = False,
    limit: int = 50,
) -> list[dict]:
    """List airways. Filter by name/notes text or minimum velocity (m/s); sort by any
    returned field (e.g. 'Quantity', 'Length', 'Velocity m/s', 'Pressure Branch')."""
    ws = vsm.airways(_load(path))
    if name_contains:
        s = name_contains.lower()
        ws = [w for w in ws if s in (w["Branch Name"] + " " + w["Notes"]).lower()]
    if min_velocity is not None:
        ws = [w for w in ws if (w["Velocity m/s"] or 0) >= min_velocity]
    if ws and sort_by not in ws[0]:
        raise ValueError(f"sort_by must be one of {list(ws[0])}")

    def key(w):
        v = w[sort_by]
        return (v is None, v if v is not None else 0)

    if sort_by != "Unique No" or descending:
        ws.sort(key=key, reverse=descending)
    return ws[: max(1, min(limit, 500))]


@mcp.tool()
def get_airway(path: str, unique_no: str) -> dict:
    """Every column of one airway (by 'Unique No'), including heat and gas fields."""
    m = _load(path)
    for r in m.rows():
        if m.value(r, "Unique No") == str(unique_no):
            return {c: r[i] for i, c in enumerate(m.columns) if c and i < len(r) and r[i] != ""}
    raise ValueError(f"no airway with Unique No {unique_no}")


@mcp.tool()
def list_fans(path: str) -> list[dict]:
    """Airways that carry a fan, with duty airflow, fan pressure and power."""
    return vsm.fans(_load(path))


@mcp.tool()
def check_flow_balance(path: str, tolerance: float = 0.5, limit: int = 30) -> dict:
    """Mass-balance check of simulated airflow (m3/s) at every junction.

    Junctions are found by matching airway end coordinates. Nodes touched by a single
    airway (portals, dead ends) are skipped. A large imbalance usually means the model
    needs to be re-simulated, contains recirculation fans, or has disconnected nodes."""
    m = _load(path)
    net: dict[tuple, float] = defaultdict(float)
    deg: dict[tuple, int] = defaultdict(int)
    for r in m.rows():
        a = tuple(round(m.num(r, c), 2) for c in ("X1", "Y1", "Z1"))
        b = tuple(round(m.num(r, c), 2) for c in ("X2", "Y2", "Z2"))
        q = m.num(r, "Quantity")
        net[a] -= q  # leaves node a
        net[b] += q  # arrives at node b
        deg[a] += 1
        deg[b] += 1
    bad = [
        {"x": k[0], "y": k[1], "z": k[2], "imbalance_m3s": round(v, 3), "airways": deg[k]}
        for k, v in net.items()
        if deg[k] > 1 and abs(v) > tolerance
    ]
    bad.sort(key=lambda d: -abs(d["imbalance_m3s"]))
    return {
        "junctions_checked": sum(1 for k in net if deg[k] > 1),
        "unbalanced": len(bad),
        "worst": bad[:limit],
    }


@mcp.tool()
def export_csv(path: str, out_path: str, columns: list[str] | None = None) -> str:
    """Write airways to a CSV file. 'columns' may name any column of the .vsm table;
    by default the common fields are exported."""
    m = _load(path)
    out_path = os.path.abspath(os.path.expanduser(out_path))
    cols = columns or list(vsm.BASIC_FIELDS)
    for c in cols:
        m.col(c)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for r in m.rows():
        w.writerow([m.value(r, c) for c in cols])
    with open(out_path, "w", encoding="utf-8-sig", newline="") as f:
        f.write(buf.getvalue())
    return f"{len(m.row_indexes)} airways written to {out_path}"


@mcp.tool()
def edit_airways(path: str, out_path: str, edits: list[dict]) -> str:
    """EXPERIMENTAL. Copy the model to 'out_path' (never overwrites the source) after
    changing text fields. Each edit: {"unique_no": "12", "field": "Notes", "value": "..."}.
    Allowed fields: Notes, Group, Branch Name. Geometry and results are not touched, so
    open the copy in Ventsim and re-simulate."""
    allowed = {"Notes", "Group", "Branch Name"}
    src = os.path.abspath(os.path.expanduser(path))
    dst = os.path.abspath(os.path.expanduser(out_path))
    if src == dst:
        raise ValueError("out_path must differ from the source file")
    if os.path.exists(dst):
        raise ValueError("out_path already exists; choose a new file name")
    m = vsm.parse(src)  # fresh copy, not the cached one
    index = {m.value(m.row(n), "Unique No"): n for n in range(len(m.row_indexes))}
    for e in edits:
        if e["field"] not in allowed:
            raise ValueError(f"field not editable: {e['field']} (allowed: {sorted(allowed)})")
        value = str(e["value"])
        if "\t" in value or "\r" in value or "\n" in value:
            raise ValueError("value must not contain tabs or line breaks")
        n = index.get(str(e["unique_no"]))
        if n is None:
            raise ValueError(f"no airway with Unique No {e['unique_no']}")
        cells = m.row(n)
        cells[m.col(e["field"])] = value
        m.lines[m.row_indexes[n]] = "\t".join(cells)
    vsm.dump(m, dst)
    return f"{len(edits)} edit(s) written to {dst}"


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
