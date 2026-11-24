import numpy as np, scipy.linalg as sl, sys
import wopt as W
x, h, B, L, psi, N, n = W.x, W.h, W.B, W.L, W.psi, W.N, W.n
Df = (np.eye(N, k=1) - np.eye(N))/h; Df[-1, :] = 0
xm = x + h/2
def omega2_G(a, eps_fun, nt=2, ret=False):
    phi = W.beta*(1 + x*x)**(2/3) + W.G @ a
    rho = np.exp(-phi)
    phim = W.beta*(1 + xm*xm)**(2/3) + np.array([np.exp(-(xm-c)**2/(2*W.S**2)) + np.exp(-(xm+c)**2/(2*W.S**2)) for c in W.XK]).T @ a
    G = np.diag(rho*h) + Df.T @ np.diag(np.exp(-phim)*h*eps_fun(xm)) @ Df
    Ge = B.T @ G @ B; Le = B.T @ L @ B
    R = np.linalg.cholesky(Ge).T
    M = R @ Le @ np.linalg.inv(R)
    He = (M + M.T)/2
    C = R @ (B.T @ psi)
    Qm, _ = np.linalg.qr(C, mode='complete'); Z = Qm[:, 2:]
    w = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-nt, Z.shape[1]-1])[::-1]
    return w
if __name__ == "__main__":
    a0 = np.zeros(len(W.XK)); abest = np.load('wopt_best.npy')
    for e0 in [0.0, 0.3, 1.0, 2.0, 4.0, 8.0, 16.0]:
        f = lambda z, e0=e0: e0*np.ones_like(z)
        print(f"eps={e0:5.2f}  plain: {np.round(omega2_G(a0, f),4)}   optw: {np.round(omega2_G(abest, f),4)}", flush=True)
