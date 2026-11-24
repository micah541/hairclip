"""Numerical (non-rigorous) numerical range of L in L^2(rho), rho = exp(-beta (1+x^2)^(2/3)).
Gives the best |Im lam| bound obtainable from numerical-range arguments:
   |Im lam| <= min_th (Lam(th) + delta cos th)/sin th,  Lam(th) = sup spec Re(e^{i th} L)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, sys
from evans import PROF, STIP
def coeffs(x):
    s = np.abs(x); sig = s + STIP
    P = np.zeros((4, len(x)))
    inside = sig <= 22
    for i in np.where(inside)[0]:
        P[:, i] = PROF(sig[i])
    P[1, ~inside] = sig[~inside]
    p, q, k, m = P
    sg = np.sign(x); sg[sg == 0] = 1
    q = q * sg; m = m * sg                 # q, m odd; p, k even
    return p, q, k, m
def build(beta, X=80.0, h=0.04):
    x = np.arange(-X, X + h/2, h); n = len(x)
    p, q, k, m = coeffs(x)
    kpp = p / 4
    a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - 4*k*kpp
    e = np.ones(n)
    D1 = sp.diags([-e[:-1], e[:-1]], [-1, 1]) / (2*h)
    D2 = sp.diags([e[:-1], -2*e, e[:-1]], [-1, 0, 1]) / h**2
    D4 = sp.diags([e[:-2], -4*e[:-1], 6*e, -4*e[:-1], e[:-2]], [-2, -1, 0, 1, 2]) / h**4
    L = -D4 + sp.diags(a2) @ D2 + sp.diags(a1) @ D1 + sp.diags(a0)
    rho = np.exp(-beta*(1 + x*x)**(2/3))
    r = np.sqrt(rho)
    A = sp.diags(r) @ L @ sp.diags(1/r)
    return x, A.tocsr()
def Lam(A, th):
    H = (np.exp(1j*th)*A + np.exp(-1j*th)*A.conj().T)/2
    Hc = H.tocsc()
    for sig in [60.0, 400.0, 3000.0, 3e4]:
        v = sla.eigsh(Hc, k=1, sigma=sig, which='LM', tol=1e-10)[0][0]
        if v < sig:
            # make sure it is the top: re-solve with sigma just above v
            v2 = sla.eigsh(Hc, k=1, sigma=v + 1.0, which='LM', tol=1e-10)[0][0]
            return max(v, v2)
    return np.inf
if __name__ == "__main__":
    delta = 0.125
    for beta in [0.1, 0.2, 0.3, 0.4]:
        x, A = build(beta)
        L0 = Lam(A, 0.0)
        best = []
        for th in np.linspace(0.02, 1.2, 30):
            best.append(((Lam(A, th) + delta*np.cos(th))/np.sin(th), th))
        b = min(best)
        print(f"beta={beta}: sup Re W(L) = {L0:.4f};  |Im| bound = {b[0]:.3f} at th={b[1]:.3f}", flush=True)
