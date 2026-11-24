"""compare the float base with HLS snapshots (even part) at the same k"""
import sys, numpy as np
from base2 import Base; from seed import HLSHalf
B = Base(sys.argv[1])
for f in sys.argv[2:]:
    S = HLSHalf(f); tau, i = B.tau_of_k(S.kmax); Y = B.Y(tau, i)
    ell = Y[-4]
    xs = np.r_[np.linspace(0, 10, 201), np.linspace(10, 40, 121)[1:], np.geomspace(40, min(ell, S.ell) * 0.999, 200)]
    kb = B.kappa(tau, xs, 0, i); kh = S.kap(xs)
    kb2 = B.kappa(tau, xs, 2, i); kh2 = S.kap(xs, 2)
    e = np.abs(kb - kh); e2 = np.abs(kb2 - kh2)
    print(f"{f.split('/')[-1]}: k_HLS={S.kmax:.4f} tau={tau:.4f} ell base/HLS={ell:.6f}/{S.ell:.6f}  |dk| x<=10 {e[xs<=10].max():.2e}  10-40 {e[(xs>10)&(xs<=40)].max():.2e}  >40 {e[xs>40].max():.2e} | |dk''| x<=10 {e2[xs<=10].max():.2e} >40 {e2[xs>40].max():.2e}")
