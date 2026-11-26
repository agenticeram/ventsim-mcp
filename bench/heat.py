"""thermo.temperatures() with the saved airflow vs. the dry-bulb temperatures Ventsim saved,
for every .vsm in a folder: mass-flow weighted mean and 95th percentile of |error| in K."""
import glob, os, sys
import numpy as np
from ventsim_mcp import thermo, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes:
        continue
    if not any(m.num(r, "Calc Sensible") for r in m.rows()):
        print(f"{os.path.basename(f):32s} no heat simulation saved"); continue
    t = thermo.temperatures(m)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    e, w = [], []
    for uid, d in t.items():
        if d["massflow"] <= 0.5:
            continue
        e.append(abs(d["t_out"] - m.num(rows[uid], "db air out")))
        w.append(d["massflow"])
    if not e:
        continue
    e, w = np.array(e), np.array(w)
    print(f"{os.path.basename(f):32s} n={len(e):5d}  weighted={np.average(e, weights=w):6.3f} K"
          f"  p95={np.percentile(e, 95):6.3f} K  max={e.max():6.2f} K")
