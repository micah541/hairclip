import numpy as np, scipy.linalg as sl
import wopt as W
x, h, B, L, psi, N, n = W.x, W.h, W.B, W.L, W.psi, W.N, W.n
D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
def omega2_G(a, eps_fun, nt=1):
    phi = W.beta*(1 + x*x)**(2/3) + W.G @ a
    rho = np.exp(-phi)
    G = np.diag(rho*h) + D1.T @ np.diag(rho*h*eps_fun(x)) @ D1
    Ge = B.T @ G @ B; Le = B.T @ L @ B
    R = np.linalg.cholesky(Ge).T            # Ge = R^T R
    M = R @ np.linalg.solve(R, Le.T).T if False else R @ Le @ np.linalg.inv(R)
    He = (M + M.T)/2
    C = R @ (B.T @ psi)
    Qm, _ = np.linalg.qr(C, mode='complete'); Z = Qm[:, 2:]
    return sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-nt, Z.shape[1]-1])[::-1]
a0 = np.zeros(len(W.XK)); abest = np.load('wopt_best.npy')
for e0 in [0.0, 0.01, 0.03, 0.1, 0.3, 1.0]:
    f = lambda z, e0=e0: e0*np.ones_like(z)
    print(f"eps={e0:5.2f}  plain weight: {omega2_G(a0, f, 2)}   opt weight: {omega2_G(abest, f, 2)}", flush=True)
