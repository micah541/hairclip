import sys, numpy as np, scipy.linalg as sl
from compress import build, top_compressed
X = float(sys.argv[1]); h = float(sys.argv[2])
for beta in [0.15, 0.2, 0.25, 0.3]:
    x, B, A, L, r, psi1, psi2 = build(beta, X, h)
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    ev, vl, vr = sl.eig(Ae, left=True, right=True)
    order = np.argsort(-ev.real)
    reals = [i for i in order if abs(ev[i].imag) < 1e-8][:5]
    lam = ev[reals].real
    R = np.real(vr[:, reals]); Lf = np.real(vl[:, reals])
    out = []
    for m in [2, 3, 4]:
        wo = top_compressed(He, R[:, :m])[0]       # orthogonal to eigenvectors
        ws = top_compressed(He, Lf[:, :m])[0]      # orthogonal to adjoint eigenvectors
        out.append(f"m={m}: orth {wo:+.4f} spec {ws:+.4f}")
    print(f"beta={beta}  lam={np.round(lam,4)}  " + " | ".join(out), flush=True)
