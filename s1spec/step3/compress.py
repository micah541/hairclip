"""Numerical (non-rigorous) study for step 3: energy decay of the modulated linear flow.

Even class, H = L^2(rho), rho = exp(-beta (1+x^2)^(2/3)).
Unstable eigenfunctions (even): psi1 = p (eigenvalue 1), psi2 = <e1,nu> = -sin(theta) (eigenvalue 1/4).
We compute
  w0   = sup Re <Lu,u>/<u,u>                    (numerical abscissa, even class)
  wc   = sup over u orthogonal to psi1, psi2    (orthogonal modulation)
  ws   = sup over u orthogonal to adjoint eigenfunctions (spectral complement)
plus the leading even eigenvalues as a sanity check.
"""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP
from scipy.integrate import cumulative_trapezoid

def profile_on(x):
    s = np.abs(x) + STIP
    P = np.zeros((4, len(x)))
    ins = s <= 22
    P[:, ins] = np.array([PROF(si) for si in s[ins]]).T
    P[1, ~ins] = s[~ins]
    p, q, k, m = P
    sg = np.where(x < 0, -1.0, 1.0)
    return p, q*sg, k, m*sg

def build(beta, X=40.0, h=0.04):
    n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x)
    p, q, k, m = profile_on(x)
    a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*p     # 4 k k'' = k p
    I = np.eye(N)
    D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
    L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
    phi = beta*(1 + x*x)**(2/3)
    r = np.exp(-phi/2)
    A = (r[:, None]*L)/r[None, :]          # r L r^{-1}: numerical range of L in L^2(rho)
    # even basis
    B = np.zeros((N, n+1)); B[n, 0] = 1.0
    for j in range(1, n+1):
        B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
    theta = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0) - cumulative_trapezoid(k, x, initial=0.0)[n]
    psi1 = p; psi2 = -np.sin(theta)
    return x, B, A, L, r, psi1, psi2

def top_compressed(H, C):
    """max of v^T H v / v^T v over v orthogonal to columns of C (H symmetric)."""
    Q, _ = np.linalg.qr(C, mode='complete')
    Z = Q[:, C.shape[1]:]
    return sl.eigh(Z.T @ H @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-3, Z.shape[1]-1])[::-1]

if __name__ == "__main__":
    X = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
    h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
    for beta in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]:
        x, B, A, L, r, psi1, psi2 = build(beta, X, h)
        Ae = B.T @ A @ B
        He = (Ae + Ae.T)/2
        w0 = sl.eigh(He, eigvals_only=True, subset_by_index=[He.shape[0]-3, He.shape[0]-1])[::-1]
        C = B.T @ np.column_stack([r*psi1, r*psi2])
        wc = top_compressed(He, C)
        # spectral complement: orthogonal to adjoint eigenvectors of Ae (eigs 1, 1/4)
        ev, vl, vr = sl.eig(Ae, left=True, right=True)
        idx = [int(np.argmin(abs(ev - 1.0))), int(np.argmin(abs(ev - 0.25)))]
        Cs = np.real(vl[:, idx])
        ws = top_compressed(He, Cs)
        evs = np.sort(ev.real[abs(ev.imag) < 1e-8])[::-1][:5]
        print(f"beta={beta:.2f}  eig={np.round(evs,4)}  w0={np.round(w0,4)}  wc={np.round(wc,4)}  ws={np.round(ws,4)}", flush=True)
