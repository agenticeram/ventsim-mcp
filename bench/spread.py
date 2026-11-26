"""thermo.arrival_times() vs. Ventsim's saved 'Spread Time' for every .vsm in a folder that
has spread results (sources = airways whose spread time is the start time, 1 s)."""
import glob, os, sys
from ventsim_mcp import thermo, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes or "Spread Time" not in m._col:
        continue
    st = {m.value(r, "Unique No"): m.num(r, "Spread Time") for r in m.rows()}
    src = [u for u, t in st.items() if t == 1]
    if not src:
        continue
    try:
        t = thermo.arrival_times(m, src)
    except ValueError:
        continue
    ok = sum(abs(d["arrive_s"] - st[u]) <= 1e-3 * st[u] for u, d in t.items() if 0 < st[u] < 99999)
    n = sum(1 for u in t if 0 < st[u] < 99999)
    print(f"{os.path.basename(f):24s} sources {len(src)}  within 0.1 %: {ok}/{n}")
