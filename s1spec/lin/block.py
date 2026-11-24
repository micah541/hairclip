"""Block (S|G) basis for the frozen tail eigen-ODE.  S = {slow, decaying} pair (may collide), G = growing pair.
 U(eps) = [e+, e-, v(z2), v(z3)],  e+ = (v(z0)+v(z1))/2, e- = (v(z0)-v(z1))/(z0-z1), v(z)=(1,z,z^2,z^3).
 With s = z0+z1, p = z0 z1:  e+ = (1, s/2, (s^2-2p)/2, (s^3-3sp)/2),  e- = (0,1,s,s^2-p).
 Frozen matrix in this basis (zeta units): blockdiag(Bs, diag(z2,z3)), Bs = [[s/2,1],[(s^2-4p)/4, s/2]].
 Z' = [mu*blockdiag - Mt/(3 sigma)] Z,  Mt = U^{-1} diag(0,1,2,3) U - 4 eps U^{-1} dU/deps."""
import numpy as np
mu = lambda s: (s/4)**(1/3)
G0 = np.array([np.exp(1j*np.pi/3), np.exp(-1j*np.pi/3)])
def groots(eps, n=200):
    z = G0.copy()
    for t in np.linspace(0, 1, n+1)[1:]:
        e = t*eps
        for _ in range(40):
            z = z - (z**4 + z - e)/(4*z**3 + 1)
    return z
def basis(eps):
    z2, z3 = groots(eps)
    sp_ = z2 + z3; pp = z2*z3
    s = -sp_; p = sp_**2 - pp
    # derivatives
    dz = 1/(4*np.array([z2, z3])**3 + 1)
    dsp = dz[0] + dz[1]; dpp = dz[0]*z3 + z2*dz[1]
    ds = -dsp; dp = 2*sp_*dsp - dpp
    v = lambda z: np.array([1, z, z*z, z**3])
    dv = lambda z, d: np.array([0, d, 2*z*d, 3*z*z*d])
    ep = np.array([1, s/2, (s*s - 2*p)/2, (s**3 - 3*s*p)/2])
    em = np.array([0, 1, s, s*s - p])
    dep = np.array([0, ds/2, (2*s*ds - 2*dp)/2, (3*s*s*ds - 3*ds*p - 3*s*dp)/2])
    dem = np.array([0, 0, ds, 2*s*ds - dp])
    U = np.column_stack([ep, em, v(z2), v(z3)])
    dU = np.column_stack([dep, dem, dv(z2, dz[0]), dv(z3, dz[1])])
    Ui = np.linalg.inv(U)
    Mt = Ui @ np.diag([0, 1, 2, 3]) @ U - 4*eps*Ui @ dU
    Bs = np.array([[s/2, 1], [(s*s - 4*p)/4, s/2]])
    return U, Mt, Bs, (z2, z3), (s, p)
def Bfull(sig, c, w=1.0):
    """Z-system matrix; the e- coordinate is rescaled by w (Zs = diag(1,1/w) ...) to tame the S-block."""
    m = mu(sig); eps = c/m**4
    U, Mt, Bs, zg, sp = basis(eps)
    Bm = -Mt/(3*sig)
    Bm[:2, :2] += m*Bs
    Bm[2, 2] += m*zg[0]; Bm[3, 3] += m*zg[1]
    Wm = np.diag([1, w, 1, 1]); Wi = np.diag([1, 1/w, 1, 1])
    return Wi @ Bm @ Wm, U @ Wm
def lognorm(A):
    return np.linalg.eigvalsh((A + A.conj().T)/2).max()
def lp_bound(c, S=8.0, w=1.0, sig_max=2000, n=3000):
    sigs = np.exp(np.linspace(np.log(S), np.log(sig_max), n))
    ratios = []; gaps = []
    for s_ in sigs:
        B, _ = Bfull(s_, c, w)
        lS = lognorm(B[:2, :2]); mG = -lognorm(-B[2:, 2:])
        g = mG - lS
        ratios.append((np.linalg.norm(B[2:, :2], 2), np.linalg.norm(B[:2, 2:], 2), g))
    # eta from sup_{tau>=S} (bGS + eta^2 bSG)/g  (fixed point)
    eta = 0.0
    for _ in range(100):
        eta = max((a + eta**2*b)/g for a, b, g in ratios)
    Lip = max(2*eta*b/g for a, b, g in ratios)
    return eta, Lip, min(g for a, b, g in ratios)
if __name__ == "__main__":
    for lam in [1.0, 0.25, -0.125, 2.0, 1.8+0.4j, -0.125+4j, 1+2j, 0.5+4.5j, -0.125+0.0001j]:
        c = 0.25 - lam
        best = min((lp_bound(c, w=w) + (w,) for w in [0.5, 1.0, 2.0, 3.0]), key=lambda t: t[0])
        print(f"lam={lam}: eta={best[0]:.4f} Lip={best[1]:.3f} min gap={best[2]:.3f} (w={best[3]})", flush=True)
