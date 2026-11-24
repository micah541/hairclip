"""Residual feasibility: arb-polished chains ending at representative k; T5 norms per slab (max over theta samples and nodes).
usage: python3 resfeas.py TAG NCHAIN k1 [k2 ...]   (k <= first slab k: chain starts at the seed)   -> resfeas_<TAG>.json"""
import sys, time, json, numpy as np
from base2 import Base
from polish_chain import polish_chain, slab_norms, KEYS

tag, nc = sys.argv[1], int(sys.argv[2]); ks = [float(v) for v in sys.argv[3:]]
B = Base(tag); out = []
k0 = np.array([B.k(t) for t in B.tau0])
for kk in ks:
    iend = int(np.argmin(np.abs(k0 - kk))); i0 = max(0, iend - nc + 1)
    t0 = time.time(); SPs = polish_chain(B, i0, iend - i0 + 1); tpol = (time.time() - t0) / (iend - i0 + 1)
    for j, SP in enumerate(SPs):
        t1 = time.time(); mx, rows, th = slab_norms(SP, 6); tres = time.time() - t1
        rec = dict(slab=i0 + j, k=float(k0[i0 + j]), n=SP.R.n, N=SP.R.N, t_polish=tpol, t_resid=tres, **mx)
        out.append(rec)
        print(f"slab {i0+j:3d} k={k0[i0+j]:9.1f} n={SP.R.n}: " + ' '.join(f'{k}={mx[k]:.1e}' for k in KEYS) + f"  (polish {tpol:.0f}s, resid {tres:.0f}s)", flush=True)
    json.dump(out, open(f'resfeas_{tag}.json', 'w'), indent=1)
