"""Search for an even weight rho = exp(-phi), phi = beta(1+x^2)^(2/3) + sum_k a_k G_k(x),
G_k = exp(-(x-x_k)^2/(2 s^2)) + exp(-(x+x_k)^2/(2 s^2)), that makes the plain weighted energy
dissipative on the rho-orthogonal complement of the two symmetry modes (omega_2 < 0).
Numerical only."""
import sys, numpy as np, scipy.linalg as sl
from scipy.optimize import minimize
from compress import profile_on
from scipy.integrate import cumulative_trapezoid

X, h, beta = 30.0, 0.05, 0.25
n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x)
p, q, k, m = profile_on(x)
a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*p
I = np.eye(N)
D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
B = np.zeros((N, n+1)); B[n, 0] = 1.0
for j in range(1, n+1):
    B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
theta = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0); theta -= theta[n] + np.pi/2
psi = np.column_stack([p, -np.sin(theta)])
XK = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.5, 7.0]); S = 0.8
G = np.array([np.exp(-(x-c)**2/(2*S*S)) + np.exp(-(x+c)**2/(2*S*S)) for c in XK]).T

def omega2(a, nt=1, ret_all=False):
    phi = beta*(1 + x*x)**(2/3) + G @ a
    r = np.exp(-phi/2)
    A = (r[:, None]*L)/r[None, :]
    Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
    C = B.T @ (r[:, None]*psi)
    Qm, _ = np.linalg.qr(C, mode='complete'); Z = Qm[:, 2:]
    w = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-nt, Z.shape[1]-1])[::-1]
    return w if ret_all else w[0]

if __name__ == "__main__":
    a0v = np.zeros(len(XK))
    print("start", omega2(a0v), flush=True)
    best = (omega2(a0v), a0v)
    for trial, start in enumerate([a0v, -0.5*np.ones(len(XK)), 0.5*np.ones(len(XK))]):
        res = minimize(omega2, start, method='Nelder-Mead',
                       options=dict(maxiter=1500, xatol=1e-3, fatol=1e-5, adaptive=True))
        print(trial, res.fun, np.round(res.x, 3), flush=True)
        if res.fun < best[0]:
            best = (res.fun, res.x)
    print("best", best[0], np.round(best[1], 4), "top3:", omega2(best[1], nt=3, ret_all=True))
    np.save('wopt_best.npy', best[1])
