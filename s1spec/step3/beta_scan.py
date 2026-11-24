"""Float: constrained abscissa of Re<Lv,v>_rho_b + theta||v''||^2 on {chi psi1, chi psi2, psi3*_b}^perp for
rho_b = exp(-b (1+x^2)^(2/3)), and the Sobolev loss factor rho_b(34)^(-1/2)."""
import numpy as np, scipy.linalg as sl
from compress import build
from proto import cutoff
X, h = 40.0, 0.04
for b in [0.06, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25]:
    x, B, A, L, r, psi1, psi2 = build(b, X, h)
    N = len(x); I = np.eye(N)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D2r = (r[:, None]*D2)/r[None, :]
    S2 = B.T @ (D2r.T @ D2r) @ B
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    ev, vl, vr = sl.eig(Ae, left=True, right=True)
    rev = np.sort(ev.real[abs(ev.imag) < 1e-8])[::-1][:5]
    i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3])
    chi = cutoff(x, 32.0)
    F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s])
    Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
    out = []
    for th in [0, 1/9]:
        H = Z.T @ (He + th*S2) @ Z
        out.append(sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0])
    loss = np.exp(b * (1 + 34**2) ** (2/3) / 2)
    print(f"beta={b:.3f}: eigs {np.round(rev, 4)}  abscissa theta=0 {out[0]:+.4f}  theta=1/9 {out[1]:+.4f}  loss(34) {loss:.2e}", flush=True)
