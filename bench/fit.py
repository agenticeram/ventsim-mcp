"""Solver fit vs. flows saved by Ventsim, for every .vsm in a folder."""
import glob, os, sys, time
import numpy as np
from ventsim_mcp import solver, vsm

folder = sys.argv[1]
for f in sorted(glob.glob(os.path.join(folder, "*.vsm"))):
    m = vsm.parse(f)
    if len(m.row_indexes) > 1500:
        continue
    t = time.time()
    s = solver.solve(m)
    saved = {m.value(r, "Unique No"): m.num(r, "Quantity") for r in m.rows()}
    q0 = np.array([saved[u] for u in s["uid"]])
    err = np.abs(abs(s["Q"]) - abs(q0))
    big = q0 != 0
    rel = float(err.sum() / np.abs(q0).sum()) if big.any() else 0.0
    print(f"{os.path.basename(f):32s} n={len(q0):5d} it={s['iterations']:3d} res={s['residual']:.4f} mean_rel={rel:7.4f} max={(err.max() if len(err) else 0):8.3f} {time.time()-t:.1f}s")
