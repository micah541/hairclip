"""Float: compressed numerical abscissa of Re<Lv,v>_rho + theta ||v''||_rho^2 on {chi psi1, chi psi2, psi3*}^perp."""
import numpy as np, scipy.linalg as sl
from compress import build
from proto import cutoff
beta, X, h = 0.25, 30.0, 0.04
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
N = len(x); I = np.eye(N)
D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
D2r = (r[:, None]*D2)/r[None, :]
S2 = B.T @ (D2r.T @ D2r) @ B
Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
ev, vl, vr = sl.eig(Ae, left=True, right=True)
i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3])
chi = cutoff(x, 30.0)
F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s])
Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
for th in [0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]:
    H = Z.T @ (He + th*S2) @ Z
    top = sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0]
    print(f"theta={th:5.2f}: sup (Re<Lv,v> + theta||v''||^2)/||v||^2 = {top:+.4f}", flush=True)
