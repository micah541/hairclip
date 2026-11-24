import sys, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
from f3geom import *; from f3fem import *
sp = s1_splines(Smax=60, dx=0.01)
for Sigma, h in [(30, 0.1), (30, 0.05), (40, 0.05), (30, 0.025)]:
    t = time.time(); m = Mesh(Sigma, h); M = m.mass(); K = m.form(coeffs(m, sp, use_tau=False))
    w, _ = eig_nonsym(K, M); ws, _ = eig_sym(K, M)
    print(f"Sigma={Sigma} h={h} ndof={len(M)} nonsym top: {np.round(w[:8].real,5)} imag max {np.max(np.abs(w[:8].imag)):.1e} | sym top: {np.round(ws[:8],4)}  #sym>-0.05: {np.sum(ws>-0.05)} #>-0.1: {np.sum(ws>-0.1)} [{time.time()-t:.0f}s]", flush=True)
