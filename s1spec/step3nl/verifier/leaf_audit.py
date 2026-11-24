"""What the leaves certify: for production-type cells report U, lb, number of leaves accepted by tolerance vs minimal width, and the gap U-lb vs max(tol,rtol|lb|).
Usage: ADD=abcd python3 leaf_audit.py"""
import bb, time
cells = [(0.0, 1e-12, 1e-12), (1e-3, 1.4e-13, 1.5e-13), (5e-4, 3e-13, 3.2e-13), (7e-5, 2.4e-13, 2.45e-13), (2.5e-4, 5e-13, 5.3e-13)]
for a, L0, L1 in cells:
    t = time.time(); r = bb.bound_cell(0.0, a, L0, L1, tol=5e-4, rtol=0.005)
    gap = r["U"] - r["lb"]; allowed = max(5e-4, 0.005 * abs(r["lb"]))
    print(f"alpha<={a:g} Ly in [{L0:g},{L1:g}]: U={r['U']:+.5f} lb={r['lb']:+.5f} U-lb={gap:.2e} (allowed {allowed:.2e})  leaves={r['nleaf']} minwidth-leaves={r['nminw']}  evals={r['nev']}  {time.time()-t:.1f}s", flush=True)
