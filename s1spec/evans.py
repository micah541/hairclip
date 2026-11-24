"""
Numerical (non-rigorous) Evans function for the linearization of the rescaled CDF at S1:
    L u = -u'''' - k^2 u'' - (5 k k' + q/4) u' + (1/4 - 3 k'^2 - 4 k k'') u .
Eigenvalue problem (L - lam) u = 0 with u admissible at the conical end (power-law mode
sigma^(1-4 lam) plus the decaying Airy mode; no growing Airy modes).  The admissible plane is
carried from sigma_far back to the tip with the second compound system (2-forms W_ij), rescaled by
the lambda-independent positive factor e^{psi}, psi' = mu(max(sigma,1)).
Even class (reflection-symmetric perturbations): D_even = W_13(tip)  (u' = u''' = 0)
Odd class:                                        D_odd  = W_02(tip)  (u = u'' = 0)
Tail normalisation: q - sigma -> 0; tip at sigma ~ -2.2491.
"""
import numpy as np
from scipy.integrate import solve_ivp

A_STAR = 4.7072774791559315e-14      # decaying-mode amplitude of S1 at sigma=22 (inward shooting)
SFAR = 22.0
w3 = np.exp(2j * np.pi * np.arange(3) / 3)
mu = lambda s: (s / 4) ** (1 / 3)


def _profile():
    u = mu(SFAR); lam = -u * w3
    A = A_STAR
    y0 = [4 * (lam[0] ** 2 * A).real, SFAR, A, (lam[0] * A).real]
    f = lambda s, y: [-y[2] * y[1], 1 + y[2] * y[0], y[3], y[0] / 4]
    ev = lambda s, y: y[1]
    ev.direction = 1
    sol = solve_ivp(f, [SFAR, -3.0], y0, method='DOP853', rtol=1e-13, atol=1e-30, dense_output=True, events=ev)
    return sol.sol, sol.t_events[0][0]


PROF, STIP = _profile()


def psi_prime(s):
    return mu(max(s, 1.0))


def rhs(s, W, lam):
    p, q, k, m = PROF(s)
    W01, W02, W03, W12, W13, W23 = W
    a0 = 0.25 - 3 * m * m - k * p - lam          # 4 k k'' = k p
    a1 = -(5 * k * m + q / 4)
    a2 = -k * k
    d = np.array([W02,
                  W03 + W12,
                  a1 * W01 + a2 * W02 + W13,
                  W13,
                  -a0 * W01 + a2 * W12 + W23,
                  -a0 * W02 - a1 * W12], dtype=complex)
    return d + psi_prime(s) * W        # Wtilde = e^{psi} W


def admissible_plane(s, lam, nterms=4, scaled=False):
    """scaled=True divides the power-law mode by s**alpha (analytic, nonvanishing in lam: winding unchanged)."""
    al = 1 - 4 * lam
    def ff(a, j):
        r = 1.0 + 0j
        for i in range(j): r *= (a - i)
        return r
    U = np.zeros(4, dtype=complex)
    c = 1.0 + 0j
    for n in range(nterms):
        e = al - 4 * n
        for d in range(4):
            U[d] += c * ff(e, d) * (s ** (-4 * n - d) if scaled else s ** (e - d))
        c = c * ff(al - 4 * n, 4) / (n + 1)
    a0 = 0.25 - lam; a1 = -s / 4
    roots = np.roots([1, 0, 0, -a1, -a0])
    nu = roots[np.argmin(abs(roots + mu(s)))]
    V = np.array([1, nu, nu ** 2, nu ** 3], dtype=complex)
    return np.array([U[i] * V[j] - U[j] * V[i] for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]])


def evans(lam, sfar=SFAR, rtol=1e-11):
    W0 = admissible_plane(sfar, lam)
    sol = solve_ivp(rhs, [sfar, STIP], W0, args=(lam,), method='DOP853', rtol=rtol, atol=1e-14)
    if sol.status != 0:
        raise RuntimeError(sol.message)
    Wt = sol.y[:, -1]
    return Wt[4], Wt[1], Wt


if __name__ == "__main__":
    print("tip sigma =", STIP)
    for lam in [1.0, 0.25, 0.0, 0.5, -0.5, 2.0, 0.3 + 0.2j]:
        De, Do, Wt = evans(lam)
        n = np.linalg.norm(Wt)
        print(f"lam={lam!s:>10}: |W|={n:.3e}  D_even/|W|={abs(De)/n:.3e}  D_odd/|W|={abs(Do)/n:.3e}")
