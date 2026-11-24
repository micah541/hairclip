"""Float: constrained abscissa of Re<Lv,v> + int rho theta(x) v''^2 on {f1,f2,f3}^perp (original cutoff at 30),
theta(x) = th1 + (th2 - th1) S((|x| - X0)/(X1 - X0))."""
import numpy as np, scipy.linalg as sl
from compress import build
from proto import cutoff
beta, X, h = 0.25, 34.0, 0.04
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
N = len(x); I = np.eye(N)
D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
D2r = (r[:, None]*D2)/r[None, :]
Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
ev, vl, vr = sl.eig(Ae, left=True, right=True)
i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3])
chi = cutoff(x, 32.0)
F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s])
Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
def S(u): u = np.clip(u, 0, 1); return u**4*(35 - 84*u + 70*u*u - 20*u**3)
for th1, th2, X0, X1 in [(1/9, 1/9, 10, 20), (1/9, 0.25, 12, 22), (1/9, 0.25, 15, 25), (1/9, 0.35, 15, 25), (1/9, 0.25, 8, 18), (1/9, 0.5, 15, 25), (1/9, 0.5, 20, 30)]:
    th = th1 + (th2 - th1)*S((np.abs(x) - X0)/(X1 - X0))
    S2 = B.T @ (D2r.T @ np.diag(th) @ D2r) @ B
    H = Z.T @ (He + S2) @ Z
    top = sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0]
    print(f"th1={th1:.3f} th2={th2:.3f} transition [{X0},{X1}]: abscissa {top:+.4f}", flush=True)
