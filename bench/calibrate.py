"""Independent check for calibrate(): treat model B's saved airflow in the airways whose
resistance differs from model A as survey measurements, fit those resistances on A, and
compare the recovered factors and the whole airflow with B."""
import sys, time
import numpy as np
from ventsim_mcp import calib, vsm

a, b = sys.argv[1], sys.argv[2]
ma, mb = vsm.parse(a), vsm.parse(b)
A = {ma.value(r, "Unique No"): r for r in ma.rows()}
B = {mb.value(r, "Unique No"): r for r in mb.rows()}
true = {}
for u in A.keys() & B.keys():
    ra, rb = ma.num(A[u], "Resistance Total"), mb.num(B[u], "Resistance Total")
    if ra > 0 and abs(ra - rb) > 1e-6 * ra:
        true[u] = rb / ra
meas = {u: mb.num(B[u], "Quantity") for u in true}
t = time.time()
r = calib.calibrate(ma, meas)
fe = [abs(r["factors"][u] / true[u] - 1) for u in true]
q = r["flows"]
tot = sum(abs(mb.num(B[u], "Quantity")) for u in q)
err = sum(abs(abs(q[u]) - abs(mb.num(B[u], "Quantity"))) for u in q)
print(f"{len(true)} airways changed; measured airflow misfit rms {r['rms_rel_error_before']:.2%} -> "
      f"{r['rms_rel_error_after']:.3%}; recovered factor error median {np.median(fe):.2%}, max "
      f"{max(fe):.1%}; whole-network airflow error vs B {err / tot:.3%}; {time.time() - t:.0f} s")
for u in sorted(true, key=lambda u: -abs(true[u] - 1))[:3]:
    print(f"  airway {u}: true factor {true[u]:.3f}, fitted {r['factors'][u]:.3f}")
