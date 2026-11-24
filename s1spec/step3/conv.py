import sys, numpy as np, scipy.linalg as sl
from compress import build, top_compressed
beta = 0.25
for X, h in [(30, 0.05), (30, 0.04), (30, 0.03), (40, 0.04)]:
    x, B, A, L, r, psi1, psi2 = build(beta, X, h)
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    ev, vr = sl.eig(Ae)
    order = [i for i in np.argsort(-ev.real) if abs(ev[i].imag) < 1e-8][:4]
    R = np.real(vr[:, order])
    C2 = B.T @ np.column_stack([r*psi1, r*psi2])
    w0 = sl.eigh(He, eigvals_only=True, subset_by_index=[He.shape[0]-1, He.shape[0]-1])[0]
    print(f"X={X} h={h}: lam={np.round(ev[order].real,4)} w0={w0:.4f} "
          f"orth2={top_compressed(He, C2)[0]:+.4f} orth3={top_compressed(He, R[:, :3])[0]:+.4f} "
          f"orth3'={top_compressed(He, np.column_stack([C2, R[:, 2]]))[0]:+.4f} orth4={top_compressed(He, R)[0]:+.4f}", flush=True)
