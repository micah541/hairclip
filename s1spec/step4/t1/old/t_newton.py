import os, numpy as np
from rep import *; from ops import *; from radau import Radau; from seed import HLSHalf; from projtest import project
S = HLSHalf('../feas2/snaps/p2_k20.npz')
z = np.load('../feas2/snaps/p2_k20.npz')
from numpy.fft import irfft
N = int(z['N']); th = 2*np.pi*np.arange(N)/N + irfft(z['th_hat'], N); mean = np.mean(np.exp(1j*th))
print('HLS closure |mean|*L =', abs(mean)*float(z['L']), ' rescaled:', abs(mean)*float(z['L'])/S.lam)
R = Rep(9, Mesh()); D = Dyn(R); RD = Radau(8)
y = project(R, S.kap, S.ell); Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1.0)
F0 = D.F(Y)
for eps, mode in [(1e-7, 'rel'), (1e-7, 'abs'), (1e-6, 'abs')]:
    if mode == 'rel': J = D.jac(Y, F0, eps)
    else:
        n = D.ntot; J = np.empty((n, n))
        for j in range(n):
            Yp = Y.copy(); dj = eps * max(abs(Y[j]), 1.0); Yp[j] += dj; J[:, j] = (D.F(Yp) - F0) / dj
    RD.factor(J, 0.15)
    Ys, hist = RD.solve(D.F, Y, 0.15, scale=sc, tol=1e-15, maxit=25)
    print(eps, mode, ' '.join(f'{h:.1e}' for h in hist))
