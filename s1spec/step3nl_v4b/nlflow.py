"""Numerical (non-rigorous) study of the modulated rescaled CDF near S1, even (reflection-symmetric) class.
Gamma = Gamma1 + u nu1 (normal graph over S1), x = sigma - sigma_tip in [0, X].
  <Gamma_tau, N> = -k_ss + (1/4 - a)<Gamma, N> - b <e1, N>,   gamma(t) = x(tau) + lam(tau) Gamma(tau),
  dt/dtau = lam^4, (log lam)_tau = -1/4 + a, x_tau = lam b e1.
u_tau = (J/alpha) * RHS, with (a,b) fixed by <u, f1>_rho = <u, f2>_rho = 0 (f1 = chi p, f2 = chi <tau,tau_tip>).
Semi-implicit: the leading term -J^{-4} u'''' is implicit with its coefficient frozen at the old step.
BCs: even at the tip; free end (u'' = u''' = 0) at x = X (outflow for the drift)."""
import sys, numpy as np, sympy as sp, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP

# ---- symbolic right-hand side in terms of the jet (u0..u4) and the profile (p,q,k,m), plus cos/sin theta
u = sp.symbols('u0:5'); p, q, k, m, ct, st = sp.symbols('p q k m ct st')


def D(f):
    r = sum(sp.diff(f, u[i]) * u[i + 1] for i in range(4))
    return r + sp.diff(f, p) * (-k * q) + sp.diff(f, q) * (1 + k * p) + sp.diff(f, k) * m + sp.diff(f, m) * p / 4 \
        + sp.diff(f, ct) * (-k * st) + sp.diff(f, st) * (k * ct)


al = 1 - k * u[0]
J = sp.sqrt(al**2 + u[1]**2)
kap = (k * J**2 + al * u[2] - u[1] * D(al)) / J**3
ks = D(kap) / J
kss = D(ks) / J
GN = (al * (p + u[0]) - u[1] * q) / J            # <Gamma, N>
E1N = (al * (-st) - u[1] * ct) / J               # <e1, N>;  nu1 = (-sin, cos), tau1 = (cos, sin)
R0 = (J / al) * (-kss + GN / 4)
Ra = -(J / al) * GN
Rb = -(J / al) * E1N
lead = 1 / J**4
args = list(u) + [p, q, k, m, ct, st]
fR0, fRa, fRb, flead, fkap = [sp.lambdify(args, e, 'numpy') for e in (R0, Ra, Rb, lead, kap)]


def profile(x):
    s = x + STIP
    P = np.zeros((4, len(x)))
    ins = s <= 22
    P[:, ins] = PROF(s[ins])
    P[1, ~ins] = s[~ins]
    return P


class Flow:
    def __init__(self, X=34.0, h=0.05, beta=0.25):
        self.n = n = int(round(X / h)); self.h = h
        self.x = x = h * np.arange(n + 1)                 # unknowns u_0..u_n
        self.P = profile(x)
        pp, qq, kk, mm = self.P
        th = -np.pi / 2 + np.concatenate([[0], np.cumsum((kk[1:] + kk[:-1]) / 2) * h])
        self.ct, self.st = np.cos(th), np.sin(th)
        self.rho = np.exp(-beta * (1 + x * x) ** (2 / 3))
        w = np.full(n + 1, h); w[0] = h / 2; w[-1] = h / 2
        self.wq = 2 * w * self.rho                       # full-line L^2(rho) product of even functions
        sig = x + STIP; uu = np.clip((32 - sig) / 3, 0, 1)
        chi = uu**4 * (35 - 84 * uu + 70 * uu**2 - 20 * uu**3)
        self.F = np.array([chi * pp, chi * (-self.st)])   # f1 = chi p, f2 = chi <tau,tau_tip> = -chi sin(theta)

    def ext(self, v):
        n = self.n
        left = v[[2, 1]]                                 # even reflection at the tip
        g1 = 2 * v[n] - v[n - 1]                         # free end at X: u'' = u''' = 0
        g2 = 4 * v[n] - 4 * v[n - 1] + v[n - 2]
        return np.concatenate([left, v, [g1, g2]])

    def jet(self, v):
        e = self.ext(v); h = self.h
        d1 = (e[3:-1] - e[1:-3]) / (2 * h)
        d2 = (e[3:-1] - 2 * e[2:-2] + e[1:-3]) / h**2
        d3 = (e[4:] - 2 * e[3:-1] + 2 * e[1:-3] - e[:-4]) / (2 * h**3)
        d4 = (e[4:] - 4 * e[3:-1] + 6 * e[2:-2] - 4 * e[1:-3] + e[:-4]) / h**4
        return e[2:-2], d1, d2, d3, d4

    def D4mat(self):
        N = self.n + 1
        M = np.zeros((N, N))
        for i in range(N):
            e = np.zeros(N); e[i] = 1.0
            M[:, i] = self.jet(e)[4]
        return M

    def ip(self, f, g):
        return np.sum(self.wq * f * g)

    def args(self, v):
        J0, J1, J2, J3, J4 = self.jet(v)
        pp, qq, kk, mm = self.P
        return (J0, J1, J2, J3, J4, pp, qq, kk, mm, self.ct, self.st), J4

    def rhs_parts(self, v):
        A, J4 = self.args(v)
        return fR0(*A), fRa(*A), fRb(*A), flead(*A), J4


def run(fl, u0, dt=1e-3, T=20.0, every=1.0, verbose=True, keep=False):
    v = u0.copy(); D4 = fl.D4mat(); n = fl.n + 1
    F = fl.F
    loglam, xc, t = 0.0, 0.0, 0.0
    hist, snaps = [], []
    nsteps = int(round(T / dt)); out_every = int(round(every / dt))
    for it in range(nsteps + 1):
        R0v, Rav, Rbv, c, d4 = fl.rhs_parts(v)
        if it % out_every == 0 or it == nsteps:
            nrm = np.sqrt(fl.ip(v, v)); sup = np.max(np.abs(v[fl.x < 10]))
            hist.append((it * dt, nrm, sup, loglam, xc, t))
            if keep: snaps.append(v.copy())
            if verbose:
                print(f"tau={it*dt:7.3f} |u|_rho={nrm:.3e} sup(x<10)={sup:.3e} log lam={loglam:+.4f} x0={xc:+.3e}", flush=True)
            if not np.isfinite(nrm) or nrm > 5: break
        if it == nsteps: break
        M = np.eye(n) + dt * c[:, None] * D4
        ab = np.zeros((5, n))
        for o in range(-2, 3):
            dg = np.diagonal(M, o)
            if o >= 0: ab[2 - o, o:] = dg
            else: ab[2 - o, :n + o] = dg
        sol = sl.solve_banded((2, 2), ab, np.column_stack([v + dt * (R0v + c * d4), dt * Rav, dt * Rbv]))
        base, ya, yb = sol[:, 0], sol[:, 1], sol[:, 2]
        G = np.array([[fl.ip(F[i], ya), fl.ip(F[i], yb)] for i in range(2)])
        rhs = -np.array([fl.ip(F[i], base) for i in range(2)])
        a, b = np.linalg.solve(G, rhs)
        v = base + a * ya + b * yb
        t += dt * np.exp(4 * loglam); xc += dt * np.exp(loglam) * b; loglam += dt * (-0.25 + a)
    return (np.array(hist), v, snaps) if keep else (np.array(hist), v)


if __name__ == "__main__":
    fl = Flow()
    R0v, Rav, Rbv, c, d4 = fl.rhs_parts(np.zeros(fl.n + 1))
    print("max|R0(0)| =", np.max(np.abs(R0v)))
