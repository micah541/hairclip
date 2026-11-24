"""Numerical prototype of the rigorous scheme (m = 3, plain L^2(rho), beta = 1/4).

Constraint functions f1 = chi*p, f2 = chi*(-sin theta), f3 = chi*psi3 (chi = smooth cutoff,
1 on [0,S0-1], 0 beyond S0). With A = -H_sym - omega (H_sym = symmetric part of L in L^2_even(rho)):
   n(A)  = # eigenvalues of H_sym above -omega        (expect 3)
   G     = [<A^{-1} f_i, f_j>]                           (expect 3 negative eigenvalues)
   top   = sup Re<Lv,v>/|v|^2 over v perp f1,f2,f3       (expect about -0.082 <= -omega)
"""
import sys, numpy as np, scipy.linalg as sl
from compress import build

def cutoff(x, S0):
    t = np.clip((np.abs(x) - (S0 - 1.0)), 0, 1)
    s = 1 - t**4*(35 - 84*t + 70*t**2 - 20*t**3)      # C^3 smoothstep, 1 -> 0
    return np.where(np.abs(x) <= S0 - 1.0, 1.0, s)

def run(beta=0.25, X=30.0, h=0.04, S0=26.0, omega=0.05):
    x, B, A, L, r, psi1, psi2 = build(beta, X, h)
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2           # conjugated, even, standard inner product
    ev, vr = sl.eig(Ae)
    i3 = int(np.argmin(abs(ev + 0.2563)))
    psi3c = np.real(B @ vr[:, i3])                  # conjugated eigenvector on full grid
    psi3 = psi3c/r
    chi = cutoff(x, S0)
    F = np.column_stack([chi*psi1, chi*psi2, chi*psi3])
    Fc = B.T @ (r[:, None]*F)                       # conjugated, even basis
    Aop = -He - omega*np.eye(He.shape[0])
    w = sl.eigh(Aop, eigvals_only=True)
    nA = int(np.sum(w < 0))
    G = Fc.T @ np.linalg.solve(Aop, Fc)
    G = (G + G.T)/2
    gev = np.linalg.eigvalsh(G)
    Qm, _ = np.linalg.qr(Fc, mode='complete'); Z = Qm[:, 3:]
    top = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-1, Z.shape[1]-1])[0]
    # normalized Gram (scale-free check)
    D = np.sqrt(np.abs(np.diag(G)))
    return dict(nA=nA, Aeigs=np.sort(w)[:5], Geig=gev, Gnorm=G/np.outer(D, D), top=top)

if __name__ == "__main__":
    for S0 in [20.0, 26.0]:
        for om in [0.03, 0.05, 0.07]:
            R = run(S0=S0, omega=om)
            print(f"S0={S0} omega={om}: n(A)={R['nA']}  A lowest={np.round(R['Aeigs'][:4],4)}  "
                  f"eig(G)={R['Geig']}  top(perp f)={R['top']:+.4f}", flush=True)
    R = run()
    print("normalized G (S0=26, omega=0.05):\n", np.round(R['Gnorm'], 4))
