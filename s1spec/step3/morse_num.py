"""Non-rigorous check of the Morse count for A = -H_sym - omega (even class):
   n(A_even) = #{zeros of W01 of the admissible plane on (0, x0)} + [W13(0) W01(0) < 0]."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP
from scipy.integrate import solve_ivp
beta = 0.25
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
IDX = {pq: n for n, pq in enumerate(PAIRS)}

import sympy as _sp
_X = _sp.symbols('x')
_PHI = _sp.Rational(1, 4)*(1 + _X**2)**_sp.Rational(2, 3)
_PD = [_sp.lambdify(_X, _sp.diff(_PHI, _X, j), 'numpy') for j in range(1, 5)]
def phis(x):
    return tuple(float(f(x)) for f in _PD)

def prof(x):
    sig = x + STIP
    if sig <= 22:
        return PROF(sig)
    return np.array([0.0, sig, 0.0, 0.0])

def coef(x, lam):
    p, q, k, m = prof(x)
    P1, P2, P3, P4 = phis(x)
    A = 2*(P1**2 - P2) + k*k
    V = (-P1**4/2 + 3*P1**2*P2 - P1**2*k*k/2 - 2*P1*P3 - P1*k*m/2 - P1*q/8 - 1.5*P2**2 + P2*k*k/2
         + P4/2 - k*p/2 - 1.5*m*m + 3/8)
    M = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, -A, P1, -1], [-(V - lam), 0, 0, P1]], float)
    return M

def compound(M):
    C = np.zeros((6, 6))
    for n, (i, j) in enumerate(PAIRS):
        for kk in range(4):
            # d/dx W_ij = sum_k M_ik W_kj + M_jk W_ik
            for (a, b), sgn in (((kk, j), 1), ((i, kk), 1)):
                pass
        for kk in range(4):
            if kk != j:
                a, b = (kk, j) if kk < j else (j, kk); s = 1 if kk < j else -1
                C[n, IDX[(a, b)]] += M[i, kk]*s
            if kk != i:
                a, b = (i, kk) if i < kk else (kk, i); s = 1 if i < kk else -1
                C[n, IDX[(a, b)]] += M[j, kk]*s
    return C

def stable_plane(x, lam):
    M = coef(x, lam)
    ev, V = np.linalg.eig(M)
    order = np.argsort(ev.real)
    vs = V[:, order[0]]      # one of the stable conjugate pair
    Y = np.column_stack([vs.real, vs.imag])
    return Y, ev[order]

def run(lam=-0.05, x0=34.0):
    Y, ev = stable_plane(x0, lam)
    W0 = np.array([Y[i, 0]*Y[j, 1] - Y[j, 0]*Y[i, 1] for (i, j) in PAIRS])
    f = lambda x, W: compound(coef(x, lam)) @ W
    sol = solve_ivp(f, [x0, 0.0], W0/np.max(abs(W0)), method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True,
                    max_step=0.02)
    xs = np.linspace(x0, 0.0, 200001)
    Ws = sol.sol(xs)
    w01 = Ws[0]
    zeros = xs[:-1][np.sign(w01[:-1]) != np.sign(w01[1:])]
    Wt = sol.y[:, -1]
    return ev, zeros, Wt

if __name__ == "__main__":
    for lam in [-0.03, -0.05, -0.07]:
        ev, zeros, Wt = run(lam)
        w01, w13 = Wt[IDX[(0, 1)]], Wt[IDX[(1, 3)]]
        delta = 1 if w01*w13 < 0 else 0
        print(f"lam={lam}: frozen eigs at x0 {np.round(ev,3)}; zeros of W01 at x = {np.round(zeros,4)};"
              f" W01(0)={w01:+.4e} W13(0)={w13:+.4e} -> count = {len(zeros)} + {delta} = {len(zeros)+delta}")
