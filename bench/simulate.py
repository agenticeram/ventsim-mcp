"""simulate() from model inputs only vs. the airflow Ventsim saved, for every .vsm in a folder.
Columns: error without heat, error with the saved natural ventilation and density."""
import glob, os, sys
import numpy as np
from ventsim_mcp import solver, vsm

for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.vsm"))):
    m = vsm.parse(f)
    if not m.row_indexes:
        continue
    row = []
    for heat in (False, True):
        s = solver.simulate(m, nvp_from_file=heat, density_from_file=heat)
        saved = {m.value(r, "Unique No"): m.num(r, "Quantity") for r in m.rows()}
        q0 = np.array([saved[u] for u in s["uid"]])
        tot = np.abs(q0).sum()
        row.append(float(np.abs(np.abs(s["Q"]) - np.abs(q0)).sum() / tot) if tot else 0.0)
    print(f"{os.path.basename(f):32s} n={len(s['uid']):5d}  cold={row[0]:7.2%}  with saved heat={row[1]:7.2%}")
