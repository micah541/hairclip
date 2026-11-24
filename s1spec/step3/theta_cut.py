"""Float: constrained abscissa of Re<Lv,v> + theta ||v''||^2 on {chi_R psi1, chi_R psi2, psi3*}^perp vs R
(chi_R = C^3 smoothstep from 1 at |x| = R to 0 at |x| = R + 3)."""
import numpy as np, scipy.linalg as sl
from compress import build
beta, X, h = 0.25, 30.0, 0.04
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
N = len(x); I = np.eye(N)
D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
D2r = (r[:, None]*D2)/r[None, :]
S2 = B.T @ (D2r.T @ D2r) @ B
Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
ev, vl, vr = sl.eig(Ae, left=True, right=True)
i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3])
def chi(R):
    u = np.clip((R + 3 - np.abs(x)) / 3, 0, 1); return u**4 * (35 - 84*u + 70*u*u - 20*u**3)
for R in [6, 8, 10, 12, 15, 20, 26]:
    F = np.column_stack([B.T @ (r*chi(R)*psi1), B.T @ (r*chi(R)*psi2), psi3s])
    Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
    out = []
    for th in [0, 1/9, 0.2]:
        H = Z.T @ (He + th*S2) @ Z
        out.append(sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0])
    print(f"R={R:3d}: theta=0 {out[0]:+.4f}  theta=1/9 {out[1]:+.4f}  theta=0.2 {out[2]:+.4f}", flush=True)
