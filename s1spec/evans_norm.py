import numpy as np
import evans as E
from scipy.integrate import solve_ivp

_P = E.PROF
E.PROF = lambda s: _P(s) if s <= 22.0 else np.array([0.0, s, 0.0, 0.0])


def evans_n(lam, seg=0.5):
    """Evans function up to a positive real factor (renormalised every 'seg'); its argument is exact."""
    sfar = max(22.0, (40 * abs(lam)) ** 0.75)
    W = E.admissible_plane(sfar, lam, scaled=True); s = sfar
    W = W / np.linalg.norm(W)
    while s > E.STIP:
        s1 = max(E.STIP, s - seg)
        sol = solve_ivp(lambda t, w: E.rhs(t, w, lam) - E.psi_prime(t) * w, [s, s1], W,
                        method='DOP853', rtol=1e-10, atol=1e-13)
        W = sol.y[:, -1]; W = W / np.linalg.norm(W); s = s1
    return W[4], W[1]


if __name__ == "__main__":
    import sys
    from winding import winding
    which = sys.argv[1]; idx = 0 if which == 'even' else 1
    f = lambda lam: evans_n(lam)[idx]
    boxes = [(3.0, 25.0, -25.0, 25.0), (-0.2, 3.0, -25.0, -4.0)]
    for (x0, x1, y0, y1) in boxes:
        w, n = winding(f, [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)], n0=40)
        print(f"{which}: Re[{x0},{x1}] x Im[{y0},{y1}] winding = {w:.4f} ({n} evals)", flush=True)
