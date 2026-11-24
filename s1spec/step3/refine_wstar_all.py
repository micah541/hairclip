"""Float Newton for the even adjoint eigenfunctions at lambda = 1, 1/4 (and lamt), in the same frame as
refine_wstar.py: initial data alpha*wa - wb at sigma = 32 (decaying WKB pair), even conditions at the tip.
The symmetry eigenvalues must come out as exactly 1 and 1/4 (consistency check of the adjoint machinery)."""
import numpy as np, json
import refine_wstar as R          # reuses lt_frame/tip (its own Newton runs once on import, harmlessly)
out = {}
for lam0 in [1.0, 0.25, -0.25626]:
    best = None
    for al0 in np.linspace(-3, 3, 25):
        al, lt = al0, lam0
        try:
            for it in range(12):
                y = R.tip(al, lt); r = np.array([y[1], y[3]])
                J = np.column_stack([(np.array([R.tip(al + 1e-7, lt)[1], R.tip(al + 1e-7, lt)[3]]) - r) / 1e-7,
                                     (np.array([R.tip(al, lt + 1e-7)[1], R.tip(al, lt + 1e-7)[3]]) - r) / 1e-7])
                d = np.linalg.solve(J, -r); al, lt = al + d[0], lt + d[1]
                if np.abs(d).max() < 1e-13: break
            res = np.abs(np.array([R.tip(al, lt)[1], R.tip(al, lt)[3]])).max() / np.abs(R.tip(al, lt)).max()
            if abs(lt - lam0) < 0.05 and (best is None or res < best[2]): best = (al, lt, res)
        except Exception:
            pass
    out[str(lam0)] = {"alpha": best[0], "lam": best[1], "rel_res": best[2]}
    print(f"lambda ~ {lam0}: alpha = {best[0]:+.15f}, lambda = {best[1]:.15f}, relative tip residual {best[2]:.1e}")
json.dump(out, open('wstar_all.json', 'w'), indent=1)
