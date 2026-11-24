"""Float: constrained abscissa in a hybrid weight: rho for |x| <= x1 (C^1 match), polynomial (1+x^2)^(-p) beyond,
optionally switching to a steep tail exp(-gamma (|x| - x2)^2) beyond x2."""
import numpy as np, scipy.linalg as sl, sys
sys.path.insert(0, '/home/claude/s1spec')
from compress import profile_on
from proto import cutoff
from scipy.integrate import cumulative_trapezoid
def run(p, X=60.0, h=0.05, x2=None, gam=0.5):
    n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x)
    pp_, q, k, m = profile_on(x)
    a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*pp_
    I = np.eye(N)
    D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
    L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
    ax = np.abs(x)
    # C^1 match point: (1+x^2)^(2/3) = 6p   (phi'_rho = 2p x/(1+x^2))
    x1 = np.sqrt((6*p) ** 1.5 - 1)
    phi_r = 0.25*(1 + ax**2)**(2/3)
    phi_p = 0.25*(1 + x1**2)**(2/3) + p*np.log((1 + ax**2)/(1 + x1**2))
    phi = np.where(ax <= x1, phi_r, phi_p)
    if x2 is not None: phi = phi + np.where(ax > x2, gam*(ax - x2)**2, 0.0)
    r = np.exp(-phi/2)
    A = (r[:, None]*L)/r[None, :]
    B = np.zeros((N, n+1)); B[n, 0] = 1.0
    for j in range(1, n+1): B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
    D2r = (r[:, None]*D2)/r[None, :]
    S2 = B.T @ (D2r.T @ D2r) @ B
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    ev, vl, vr = sl.eig(Ae, left=True, right=True)
    rev = np.sort(ev.real[abs(ev.imag) < 1e-8])[::-1][:4]
    i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3])
    theta = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0); theta -= theta[n] + np.pi/2
    psi1, psi2 = pp_, -np.sin(theta)
    chi = cutoff(x, 32.0)
    F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s])
    Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
    out = []
    for th in [0, 1/9]:
        H = Z.T @ (He + th*S2) @ Z
        out.append(sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0])
    loss = np.exp(np.interp(34, ax[n:], phi[n:]) / 2)
    return x1, rev, out, loss
for p, x2 in [(2.5, None), (3.0, None), (3.5, None), (3.0, 40.0), (4.0, None)]:
    x1, rev, out, loss = run(p, x2=x2)
    print(f"p={p} (x1={x1:.2f}) tail x2={x2}: eigs {np.round(rev, 4)}  abscissa theta=0 {out[0]:+.4f} theta=1/9 {out[1]:+.4f}  loss(34) {loss:.2e}", flush=True)
