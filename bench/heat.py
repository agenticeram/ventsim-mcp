"""thermo.temperatures() with the saved airflow vs. the dry- and wet-bulb temperatures
Ventsim saved, for every .vsm in a folder: mass-flow weighted mean, 95th percentile and
worst |error| in K, airways above 0.5 kg/s. A reproduction check (the heat and moisture
loads come from the same saved results), not an accuracy figure for changes."""
import glob, os, sys
import numpy as np
from ventsim_mcp import thermo, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes:
        continue
    if not any(m.num(r, "Calc Sensible") for r in m.rows()):
        print(f"{os.path.basename(f):32s} no heat simulation saved")
        continue
    t = thermo.temperatures(m)
    rows = {m.value(r, "Unique No"): r for r in m.rows()}
    line = f"{os.path.basename(f):32s}"
    for key, col in (("t_out", "db air out"), ("wb_out", "wb air out")):
        e, w = [], []
        for uid, d in t.items():
            if d["massflow"] > 0.5:
                e.append(abs(d[key] - m.num(rows[uid], col)))
                w.append(d["massflow"])
        if not e:
            break
        e, w = np.array(e), np.array(w)
        line += (f"  {col.split()[0]}: {np.average(e, weights=w):6.3f} / "
                 f"{np.percentile(e, 95):6.3f} / {e.max():6.2f} K")
    print(line + f"  n={len(e)}" if len(e) else line + "  no airways above 0.5 kg/s")
