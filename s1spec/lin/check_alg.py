"""Numerical sanity checks of the block-basis algebra used in lintail.py:
 (1) U^{-1} A_zeta U = blockdiag(Bs, z2, z3);  (2) T^{-1} T' = Mt/(3 sigma)  (finite differences)."""
import numpy as np
from block2 import G0
def groots(eps):
    z = G0.copy()
    for t in np.linspace(0, 1, 801)[1:]:
        for _ in range(50): z = z - (z**4 + z - t*eps)/(4*z**3 + 1)
    return z
def U_of(eps, w=0.5):
    z2, z3 = groots(eps); sp_ = z2 + z3; pp = z2*z3; s = -sp_; p = sp_**2 - pp
    v = lambda z: np.array([1, z, z*z, z**3])
    ep = np.array([1, s/2, (s*s - 2*p)/2, (s**3 - 3*s*p)/2]); em = np.array([0, 1, s, s*s - p])
    U = np.column_stack([ep, em, v(z2), v(z3)])
    return U, s, p, z2, z3
for c in [0.375, -1.75, -0.75+0.4j, 0.375+3.9j]:
    for sig in [8.0, 20.0]:
        mu = (sig/4)**(1/3); eps = c/mu**4
        U, s, p, z2, z3 = U_of(eps)
        Az = np.zeros((4, 4), complex); Az[0,1] = Az[1,2] = Az[2,3] = 1; Az[3,0] = eps; Az[3,1] = -1
        blk = np.linalg.inv(U) @ Az @ U
        ref = np.zeros((4, 4), complex); ref[:2,:2] = [[s/2, 1], [(s*s-4*p)/4, s/2]]; ref[2,2] = z2; ref[3,3] = z3
        # finite-difference check of T^{-1}T'
        def T(sg):
            m = (sg/4)**(1/3); U_, *_ = U_of(c/m**4); return np.diag([1, m, m*m, m**3]) @ U_
        h = 1e-5
        TT = np.linalg.inv(T(sig)) @ (T(sig + h) - T(sig - h))/(2*h)
        # analytic Mt
        z2_, z3_ = z2, z3
        dz = 1/(4*np.array([z2_, z3_])**3 + 1); sp_ = z2_ + z3_
        dsp = dz[0] + dz[1]; dpp = dz[0]*z3_ + z2_*dz[1]; ds = -dsp; dp = 2*sp_*dsp - dpp
        dep = np.array([0, ds/2, (2*s*ds - 2*dp)/2, (3*s*s*ds - 3*ds*p - 3*s*dp)/2]); dem = np.array([0, 0, ds, 2*s*ds - dp])
        dv = lambda z, d: np.array([0, d, 2*z*d, 3*z*z*d])
        dU = np.column_stack([dep, dem, dv(z2_, dz[0]), dv(z3_, dz[1])])
        Ui = np.linalg.inv(U)
        Mt = Ui @ np.diag([0, 1, 2, 3]) @ U - 4*eps*Ui @ dU
        print(f"c={c} sig={sig}: |blk-ref|={np.abs(blk-ref).max():.2e}   |T^-1T' - Mt/(3s)|={np.abs(TT - Mt/(3*sig)).max():.2e}  |Mt|={np.abs(Mt).max():.3f}")
