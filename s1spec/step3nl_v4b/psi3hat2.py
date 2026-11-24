"""psi3 by inward integration of the admissible plane (polynomial + decaying branch) from sigma = 22 to the tip,
then the even conditions u'(tip) = u'''(tip) = 0 (2x2 null vector).  Untrusted part of the certificate."""
import sys, json, numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/home/claude/s1spec')
import evans as E
LT = -0.2561532224577463
def rhs(s, U):
    p, q, k, m = E.PROF(s); u, u1, u2, u3 = U
    return [u1, u2, u3, -k * k * u2 - (5 * k * m + q / 4) * u1 + (0.25 - 3 * m * m - k * p - LT) * u]
S0 = 22.0; A = 1 - 4 * LT; z = -(S0 / 4) ** (1 / 3)
P = np.array([[S0**A, A * S0**(A - 1), A * (A - 1) * S0**(A - 2), A * (A - 1) * (A - 2) * S0**(A - 3)],
              [1.0, z, z * z, z**3]])
sols = [solve_ivp(rhs, [22.0, E.STIP], P[i], method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True) for i in range(2)]
Mt = np.array([[sols[i].y[1, -1], sols[i].y[3, -1]] for i in range(2)]).T
_, sv, vt = np.linalg.svd(Mt); c = vt[-1]
u0 = sum(c[i] * sols[i].y[0, -1] for i in range(2)); c = c / u0
psi = lambda s: sum(c[i] * sols[i].sol(s)[0] for i in range(2))
print("singular values of the tip matrix:", sv, " (second ~ 0 means lamt is an eigenvalue)")
print("psi3 at x = 2, 5, 10:", [round(float(psi(E.STIP + x)), 4) for x in (2, 5, 10)], " (FD: -0.452, 5.0, 40.79)")
np.save('psi3hat_coef.npy', np.array(c)); json.dump({"P": P.tolist(), "c": np.array(c).tolist()}, open('psi3hat.json', 'w'))
