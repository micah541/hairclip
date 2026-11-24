"""Float: hybrid weight rho -> polynomial (p) -> rho-like tail, smooth phi' (C^3 blend); constrained abscissa with the
constraint cutoff at x_c; also slow-mode geometry (cos angle between psi3 and psi3*) and the Sobolev loss at 34."""
import numpy as np, scipy.linalg as sl, sys
sys.path.insert(0, '/home/claude/s1spec')
from compress import profile_on
from scipy.integrate import cumulative_trapezoid
def S(u): u = np.clip(u, 0, 1); return u**4*(35 - 84*u + 70*u*u - 20*u**3)
def phip(ax, p, a1, b1, a2, b2):
    """phi' = blend: rho' (|x|<a1) -> 2p x/(1+x^2) on [b1, a2] -> rho' (|x|>b2)"""
    fr = (1/3)*ax*(1 + ax**2)**(-1/3); fp = 2*p*ax/(1 + ax**2)
    s1 = S((ax - a1)/(b1 - a1)); s2 = S((ax - a2)/(b2 - a2))
    return (1 - s1)*fr + s1*((1 - s2)*fp + s2*fr)
def run(p=3.5, a1=8.0, b1=11.0, a2=38.0, b2=46.0, xc=58.0, X=70.0, h=0.05):
    n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x); ax = np.abs(x)
    pp_, q, k, m = profile_on(x)
    a2c = -k*k; a1c = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*pp_
    I = np.eye(N)
    D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
    L = -D4 + a2c[:, None]*D2 + a1c[:, None]*D1 + np.diag(a0)
    dph = phip(ax, p, a1, b1, a2, b2)
    phi = cumulative_trapezoid(dph[n:], ax[n:], initial=0.0); phi = np.concatenate([phi[:0:-1], phi])
    r = np.exp(-phi/2)
    A = (r[:, None]*L)/r[None, :]
    B = np.zeros((N, n+1)); B[n, 0] = 1.0
    for j in range(1, n+1): B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
    D2r = (r[:, None]*D2)/r[None, :]; S2 = B.T @ (D2r.T @ D2r) @ B
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    ev, vl, vr = sl.eig(Ae, left=True, right=True)
    rev = np.sort(ev.real[abs(ev.imag) < 1e-8])[::-1][:4]
    i3 = int(np.argmin(abs(ev + 0.2563))); psi3s = np.real(vl[:, i3]); psi3 = np.real(vr[:, i3])
    cosang = abs(psi3 @ psi3s)/np.linalg.norm(psi3)/np.linalg.norm(psi3s)
    theta = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0); theta -= theta[n] + np.pi/2
    psi1, psi2 = pp_, -np.sin(theta)
    chi = 1 - S((ax - (xc - 3))/3)
    F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s*(B.T @ chi)])
    Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
    out = []
    for th in [0, 1/9]:
        H = Z.T @ (He + th*S2) @ Z
        out.append(sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-1, H.shape[0]-1])[0])
    loss34 = np.exp(np.interp(34, ax[n:], phi[n:])/2); loss58 = np.exp(np.interp(58, ax[n:], phi[n:])/2)
    return rev, out, cosang, loss34, loss58
for args in [dict(), dict(p=3.0), dict(p=4.0), dict(a2=30, b2=38, xc=50), dict(p=3.5, a1=6, b1=10)]:
    rev, out, ca, l34, l58 = run(**args)
    print(f"{args}: eigs {np.round(rev, 4)} abscissa th=0 {out[0]:+.4f} th=1/9 {out[1]:+.4f}  cos(psi3,psi3*) {ca:.3f}  loss34 {l34:.2e} loss58 {l58:.2e}", flush=True)
