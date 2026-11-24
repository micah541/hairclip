"""Modulated operator Lm = (I - Pi_ob) L on W = {f1,f2}^perp (Pi_ob: oblique projection onto span{psi1,psi2}
along W). Decompose w = c phi + v, phi = unit vector of (I - Pi_F2) f3, v perp f1,f2,f3.
E_kappa(w) = kappa |c|^2 + |v|^2. Report: a = <Lm phi,phi>, |g| = |P_v Lm phi|, |l| = |P_phi Lm^* restricted|,
omega' = sup <Lm v,v>/|v|^2 on v perp f, and the exact best numerical abscissa of Lm in E_kappa over kappa."""
import numpy as np, scipy.linalg as sl
from compress import build, profile_on
from proto import cutoff
beta, X, h = 0.25, 30.0, 0.04
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
p, q, k, m = profile_on(x)
Ae = B.T @ A @ B                      # conjugated L on even subspace (standard inner product)
ev, vr = sl.eig(Ae); i3 = int(np.argmin(abs(ev + 0.2563)))
psi3c = np.real(vr[:, i3])
chi = cutoff(x, 30.0)
P1c = B.T @ (r*psi1); P2c = B.T @ (r*psi2)              # exact symmetry modes (conjugated)
def analyse(f3c, label):
    F2 = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2)])
    Psi = np.column_stack([P1c, P2c])
    # oblique projection onto span Psi along F2^perp:  Pi u = Psi (F2^T Psi)^{-1} F2^T u
    Piob = Psi @ np.linalg.solve(F2.T @ Psi, F2.T)
    Lm = (np.eye(Ae.shape[0]) - Piob) @ Ae
    Q2, _ = np.linalg.qr(F2)
    phi = f3c - Q2 @ (Q2.T @ f3c); phi /= np.linalg.norm(phi)
    F3 = np.column_stack([F2, f3c]); Q3, _ = np.linalg.qr(F3, mode='complete')
    Zv = Q3[:, 3:]                                       # basis of v-space
    a = phi @ Lm @ phi
    g = Zv.T @ (Lm @ phi)
    lvec = Zv.T @ (Lm.T @ phi)                           # l(v) = <Lm v, phi> = v . Lm^T phi
    Hv = Zv.T @ Lm @ Zv; Hv = (Hv + Hv.T)/2
    wv = sl.eigh(Hv, eigvals_only=True, subset_by_index=[Hv.shape[0]-1, Hv.shape[0]-1])[0]
    # exact: basis of W = [phi, Zv]; E_kappa Gram = diag(kappa, I)
    Wb = np.column_stack([phi, Zv])
    Lw = Wb.T @ Lm @ Wb
    best = (np.inf, None)
    for kap in np.exp(np.linspace(np.log(1e-3), np.log(1e3), 61)):
        S = np.ones(Wb.shape[1]); S[0] = np.sqrt(kap)
        Mk = (S[:, None]*Lw)/S[None, :]; Hk = (Mk + Mk.T)/2
        t = sl.eigh(Hk, eigvals_only=True, subset_by_index=[Hk.shape[0]-1, Hk.shape[0]-1])[0]
        if t < best[0]: best = (t, kap)
    print(f"{label:6s}: a={a:+.4f} |g|={np.linalg.norm(g):.4f} |l|={np.linalg.norm(lvec):.4f}  "
          f"omega'={wv:+.4f}  crude test |l||g| < |a| |omega'|: {np.linalg.norm(lvec)*np.linalg.norm(g):.4f} vs {abs(a)*abs(wv):.4f}"
          f"   exact best abscissa {best[0]:+.4f} at kappa={best[1]:.3g}", flush=True)
analyse(B.T @ (r*chi*k), "k")
analyse(psi3c, "psi3")
analyse(B.T @ (r*chi*p*p), "p^2")
