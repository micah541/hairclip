import numpy as np, scipy.linalg as sl, sys
from compress import profile_on
from scipy.integrate import cumulative_trapezoid
def setup(X, h, beta=0.25):
    n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x)
    p, q, k, m = profile_on(x)
    a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*p
    I = np.eye(N)
    D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
    L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
    B = np.zeros((N, n+1)); B[n, 0] = 1.0
    for j in range(1, n+1): B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
    th = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0); th -= th[n] + np.pi/2
    psi = np.column_stack([p, -np.sin(th)])
    Df = (np.eye(N, k=1) - np.eye(N))/h; Df[-1, :] = 0
    return dict(x=x, h=h, N=N, n=n, L=L, B=B, psi=psi, Df=Df, beta=beta)
def omega2(S, eps, nt=3, vec=False):
    x, h, B, L, Df = S['x'], S['h'], S['B'], S['L'], S['Df']
    xm = x + h/2
    rho = np.exp(-S['beta']*(1 + x*x)**(2/3)); rhom = np.exp(-S['beta']*(1 + xm*xm)**(2/3))
    G = np.diag(rho*h) + Df.T @ np.diag(rhom*h*eps(xm)) @ Df
    Ge = B.T @ G @ B; Le = B.T @ L @ B
    R = np.linalg.cholesky(Ge).T
    M = R @ Le @ np.linalg.inv(R); He = (M + M.T)/2
    C = R @ (B.T @ S['psi'])
    Qm, _ = np.linalg.qr(C, mode='complete'); Z = Qm[:, 2:]
    w, V = sl.eigh(Z.T @ He @ Z, subset_by_index=[Z.shape[1]-nt, Z.shape[1]-1])
    w = w[::-1]; V = V[:, ::-1]
    if vec:
        U = np.linalg.solve(R, Z @ V)          # back to nodal values on x >= 0
        return w, U, x[S['n']:]
    return w
if __name__ == "__main__":
    S = setup(30.0, 0.05)
    w, U, xh = omega2(S, lambda z: 4.0*np.ones_like(z), vec=True)
    for j in range(3):
        u = U[:, j]; i = np.argmax(abs(u))
        print(f"eig {w[j]:+.4e}  peak at x={xh[i]:.2f}")
    for X, h in [(30, 0.05), (30, 0.04), (40, 0.04), (30, 0.03)]:
        S = setup(X, h)
        for e0 in [3.0, 4.0, 6.0]:
            w = omega2(S, lambda z, e0=e0: e0*np.ones_like(z))
            print(f"X={X} h={h} eps={e0}: {np.round(w,4)}", flush=True)
