"""T7 float: projected-form coercivity (C1) along the T1 base.  c(theta) = -top eig of B + theta ||u''||^2 on {<u,phi_j> = 0}.
usage: python3 c1float.py [h] [end]   ->  c1float_<end>_<h>.json"""
import sys, json, time, numpy as np
import bgeo, fem
h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1; end = sys.argv[2] if len(sys.argv) > 2 else 'clamped'
X = 34.25 if end == 'clamped' else 40.0
B = bgeo.base()
KS = [(20.4, 0), (40, None), (200, None), (2e3, None), (2e4, None)]
st = fem.Setup(X, h, end); out = []
for kk, i in KS:
    if i is None: tau, i = B.tau_of_k(kk); th = (tau - B.tau0[i]) / B.dt[i]
    else: th = 0.5
    t0 = time.time(); g = bgeo.BaseGeom(i, th); op = fem.operator(st, g); Bm, Phiv = fem.projected(st, op)
    row = dict(k=g.k, slab=i, theta=th, cond_M3=float(np.linalg.cond(op['M3'] / np.sqrt(np.outer(np.diag(op['M3'] @ op['M3'].T), np.ones(3))))))
    for thc in (0.0, 0.02, 0.03, 0.05):
        w, _ = fem.top_eigs(st, Bm, Phiv, thc, nev=3); row[f'top_{thc}'] = [float(v) for v in w]
    # unprojected sym form on the same complement (FEAS2-type statement) and with base modes as constraints
    w0, _ = fem.top_eigs(st, op['K'], Phiv, 0.0, nev=2); row['unproj'] = [float(v) for v in w0]
    row['time'] = time.time() - t0; out.append(row)
    print(f"k={g.k:9.2f} slab {i:3d}: c(0)={-row['top_0.0'][0]:.4f} c(.02)={-row['top_0.02'][0]:.4f} c(.03)={-row['top_0.03'][0]:.4f} "
          f"c(.05)={-row['top_0.05'][0]:.4f} | 2nd {row['top_0.0'][1]:+.4f} | unprojected top {row['unproj'][0]:+.4f} | {row['time']:.1f}s", flush=True)
json.dump(out, open(f'c1float_{end}_{h}.json', 'w'), indent=1)
