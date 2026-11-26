"""Check fan_select() against the fans Ventsim actually runs: for every running fan with a
library curve, ask for its saved duty (per-fan airflow, static pressure, saved air density)
and see whether its own fan type is offered, at what speed vs. its 'Fan Speed'."""
import glob, os, sys
import numpy as np
from ventsim_mcp import design, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes:
        continue
    lib = design.fan_library(m)
    mask = vsm.active_mask(m)
    found = n = 0
    dev = []
    for r, k in zip(m.rows(), mask):
        ft = m.value(r, "Fan Type")
        if not k or ft not in lib or m.value(r, "Fan Off") == "True":
            continue
        q, ps = abs(m.num(r, "Quantity")), m.num(r, "Static Pressure")
        nf = max(int(m.num(r, "Fan Numbers", 1) or 1), 1)
        if q < 0.1 or ps < 1:
            continue
        rho = m.num(r, "Dry Massflow") / m.num(r, "Quantity") if m.num(r, "Quantity") else 1.2
        n += 1
        res = design.fan_select(m, q, ps, nf, rho if 0.3 < rho < 3 else 1.2, 250, limit=1000)
        own = [x for x in res["suitable"] + res["not_suitable"] if x["fan_type"] == ft and "speed_pct" in x]
        if own:
            found += 1
            dev.append(own[0]["speed_pct"] - (m.num(r, "Fan Speed", 100) or 100))
    if n:
        d = np.abs(dev) if dev else np.array([np.nan])
        print(f"{os.path.basename(f):28s} fans {n:3d}  own curve meets its duty: {found:3d}  "
              f"speed difference vs saved: median {np.median(d):.2f} pct-points, max {np.max(d):.2f}")
