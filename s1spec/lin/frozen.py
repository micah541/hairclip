"""Frozen-root diagonalisation of the tail eigen-ODE (model part, k = 0, q = sigma):
     u'''' = a1 u' + a0 u,   a1 = -sigma/4,  a0 = c = 1/4 - lam.
 Y = (u,u',u'',u''') = T Z,  T = D_mu V(zeta),  nu = mu zeta,  zeta^4 + zeta = eps = c/mu^4,
 Z' = [mu diag(zeta) - M(eps)/(3 sigma)] Z,  M = V^{-1} diag(0,1,2,3) V - 4 eps V^{-1} dV/deps.
 Modes: 0 slow (zeta->0), 1 decaying (zeta->-1), 2,3 growing (zeta->e^{+-i pi/3})."""
import numpy as np
Z0 = np.array([0, -1, np.exp(1j*np.pi/3), np.exp(-1j*np.pi/3)], dtype=complex)
mu = lambda s: (s/4)**(1/3)
def zeta(eps, n=200):
    """roots tracked continuously along t*eps, t in [0,1] (Newton steps)"""
    z = Z0.copy()
    for t in np.linspace(0, 1, n+1)[1:]:
        e = t*eps
        for _ in range(30):
            z = z - (z**4 + z - e)/(4*z**3 + 1)
    return z
def VM(eps, z=None):
    if z is None: z = zeta(eps)
    V = np.array([z**i for i in range(4)])
    dz = 1/(4*z**3 + 1)
    dV = np.array([i*z**(i-1)*dz if i > 0 else 0*z for i in range(4)])
    Vi = np.linalg.inv(V)
    M = Vi @ np.diag([0, 1, 2, 3]) @ V - 4*eps*Vi @ dV
    return V, M, z
def B(sig, c):
    m = mu(sig); eps = c/m**4
    V, M, z = VM(eps)
    return m*np.diag(z) - M/(3*sig), V, z
if __name__ == "__main__":
    for lam in [1.0, 0.25, -0.125, 2.0+0j, -0.125+4j, 1+2j]:
        c = 0.25 - lam
        for sig in [8, 16, 40]:
            Bm, V, z = B(sig, c)
            m = mu(sig)
            gap = min(Bm[2,2].real, Bm[3,3].real) - max(Bm[0,0].real, Bm[1,1].real)
            nGS = np.linalg.norm(Bm[2:, :2], 2); nSG = np.linalg.norm(Bm[:2, 2:], 2)
            print(f"lam={lam} sig={sig}: zeta={np.round(z,3)}  gap={gap:.3f}  |B_GS|={nGS:.4f} |B_SG|={nSG:.4f} eta~{nGS/gap:.4f}")
