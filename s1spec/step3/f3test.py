import numpy as np, scipy.linalg as sl
from compress import build, profile_on
from proto import cutoff
beta, X, h, omega = 0.25, 30.0, 0.04, 0.05
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
p, q, k, m = profile_on(x)
Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
Aop = -He - omega*np.eye(He.shape[0])
nA = int(np.sum(np.linalg.eigvalsh(Aop) < 0))
ev, vr = sl.eig(Ae); i3 = int(np.argmin(abs(ev + 0.2563))); psi3 = np.real(B @ vr[:, i3])/r
cands = {'psi3': psi3, 'k': k, 'x^2': x**2, 'x^4': x**4, 'q*m': q*m, 'p^2': p*p, 'const': np.ones_like(x),
         'k*x^2': k*x*x, 'cos(th)^2': None, 'q^2': q*q, 'm^2': m*m}
from scipy.integrate import cumulative_trapezoid
n = (len(x)-1)//2
th = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0); th -= th[n] + np.pi/2
cands['cos(th)^2'] = np.cos(th)**2
cands['p*x^2'] = p*x*x
print("n(A) =", nA)
for S0 in [30.0]:
    chi = cutoff(x, S0)
    for name, f3 in cands.items():
        F = np.column_stack([chi*psi1, chi*psi2, chi*f3])
        Fc = B.T @ (r[:, None]*F)
        G = Fc.T @ np.linalg.solve(Aop, Fc); G = (G + G.T)/2
        g = np.linalg.eigvalsh(G)
        D = np.sqrt(abs(np.diag(G))); gn = np.linalg.eigvalsh(G/np.outer(D, D))
        Qm, _ = np.linalg.qr(Fc, mode='complete'); Z = Qm[:, 3:]
        top = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-1, Z.shape[1]-1])[0]
        print(f"f3={name:10s} eig(G)={np.round(g,3)}  normalized min |eig|={min(abs(gn)):.3f}  top(perp)={top:+.4f}")
