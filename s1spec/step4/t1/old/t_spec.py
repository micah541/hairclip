import os, sys, time, numpy as np
os.environ.setdefault('OMP_NUM_THREADS', '1')
from rep import *; from ops import *; from seed import HLSHalf; from projtest import project
S = HLSHalf('../feas2/snaps/p2_k20.npz')
for P, JT, NO in [(9, 40, 80), (11, 40, 80)]:
    R = Rep(P, Mesh(JT=JT, NO=NO)); D = Dyn(R)
    y = project(R, S.kap, S.ell)
    Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
    t0 = time.time(); F0, p = D.F(Y, True); t1 = time.time()
    Mm = p['Mm'].toarray(); print(f"P={P} n={R.n}: F eval {1e3*(t1-t0):.1f} ms, cond(M) = {np.linalg.cond(Mm):.2e}, a = {p['a']:.3e}, b = {p['b']:.3e}, elld = {p['elld']:.4f} (ell/4 = {S.ell/4:.4f})")
    t0 = time.time(); J = D.jac(Y, F0); print('  jac', time.time() - t0, 's')
    w = np.linalg.eigvals(J[:R.n + 1, :R.n + 1])
    o = np.argsort(-w.real)
    print('  top eigenvalues:', np.round(w[o][:8], 4))
    print('  most negative:', w.real.min(), ' max |imag|/|real| for |w|>100:', np.max(np.abs(w.imag[np.abs(w) > 100]) / np.abs(w.real[np.abs(w) > 100])))
    # compare time derivative with HLS: dkappa/dtau at fixed x from two HLS snapshots is not available; print kdot at tip region
    yd = F0[:R.n]
    print('  |ydot| max', np.abs(yd).max())
