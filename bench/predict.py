"""Independent check: predict model B's saved airflow from model A plus the resistance
changes between them (same network, same Unique Nos)."""
import sys
from ventsim_mcp import server, vsm

a, b = sys.argv[1], sys.argv[2]
ma, mb = vsm.parse(a), vsm.parse(b)
A = {ma.value(r, "Unique No"): r for r in ma.rows()}
B = {mb.value(r, "Unique No"): r for r in mb.rows()}
scale = {}
for u in A.keys() & B.keys():
    ra, rb = ma.num(A[u], "Resistance Total"), mb.num(B[u], "Resistance Total")
    if ra > 0 and abs(ra - rb) > 1e-6 * ra:
        scale[u] = rb / ra
print("resistance changes:", {u: round(v, 3) for u, v in scale.items()})
r = server.what_if(a, resistance_scale=scale, limit=100000)
pred = {u: ma.num(A[u], "Quantity") for u in A}
for c in r["largest_changes"]:
    pred[c["Unique No"]] = c["estimated_airflow"]
common = [u for u in A.keys() & B.keys()]
tot = sum(abs(mb.num(B[u], "Quantity")) for u in common)
err = sum(abs(pred[u] - mb.num(B[u], "Quantity")) for u in common)
naive = sum(abs(ma.num(A[u], "Quantity") - mb.num(B[u], "Quantity")) for u in common)
print(f"airflow change Ventsim A->B: {naive:.2f} m3/s summed; prediction error: {err:.2f} m3/s "
      f"= {err / tot:.4%} of total flow ({1 - err / naive:.1%} of the change explained)")
