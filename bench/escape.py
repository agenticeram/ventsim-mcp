"""escape.distances() vs. Ventsim's saved 'distancetosurface' / 'distancetorefuge' (share of
airways within 0.1 % / 1 m), and recirculation groups found, for every .vsm in a folder."""
import glob, os, sys
from ventsim_mcp import escape, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes or "distancetosurface" not in m._col:
        continue
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    if not any(m.num(r, "distancetosurface") for r in rows.values()):
        continue
    d = escape.distances(m)
    res = []
    for key, col in (("to_surface", "distancetosurface"), ("to_refuge", "distancetorefuge")):
        ok = n = 0
        for u, x in d.items():
            sv = m.num(rows[u], col)
            mine = 9999.0 if x[key] is None else x[key]
            n += 1
            ok += abs(mine - sv) <= max(1e-3 * sv, 1.0) or (sv >= 9999 and x[key] is None)
        res.append(f"{key} {ok}/{n}")
    print(f"{os.path.basename(f):24s} " + "  ".join(res))
