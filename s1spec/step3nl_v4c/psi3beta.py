"""Float beta for psi3^ = U(J1 + beta J2 at sigma = 32): u'(tip) = 0; report u'''(tip) (should be ~0)."""
import sys, json, numpy as np, sympy as sp
from scipy.integrate import solve_ivp
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
LT = -0.2561532224577463; A = 1 - 4 * LT; S0 = 32.0
sg = sp.symbols('s'); Sp = -(sg / 4) ** sp.Rational(1, 3)
d1, d2, d3 = Sp, sp.diff(Sp, sg) + Sp**2, sp.diff(Sp, sg, 2) + 3 * Sp * sp.diff(Sp, sg) + Sp**3
aj = [1.0]
for j in range(1, 4):
    aj.append(aj[-1] * (A - 4*j + 4) * (A - 4*j + 3) * (A - 4*j + 2) * (A - 4*j + 1) / j)
ser = sum(sp.Float(aj[j], 30) * sg ** (sp.Float(A, 30) - 4 * j) for j in range(4))   # polynomial branch, 4 terms
J1 = np.array([float(sp.diff(ser, sg, n).subs(sg, S0)) for n in range(4)])
print("asymptotic coefficients a_j:", aj)
J2 = np.array([1.0] + [float(e.subs(sg, S0)) for e in (d1, d2, d3)])
def prof(s): return E.PROF(s) if s <= 22 else np.array([0.0, s, 0.0, 0.0])
def rhs(s, U):
    p, q, k, m = prof(s); u, u1, u2, u3 = U
    return [u1, u2, u3, -k * k * u2 - (5 * k * m + q / 4) * u1 + (0.25 - 3 * m * m - k * p - LT) * u]
def tip(beta):
    return solve_ivp(rhs, [S0, E.STIP], J1 + beta * J2, method='DOP853', rtol=1e-13, atol=1e-30).y[:, -1]
y1, y2 = tip(0.0), tip(1.0) - tip(0.0)                 # linear in beta
beta = -y1[1] / y2[1]; yt = y1 + beta * y2
print(f"beta = {beta:.15e};  tip jet / u: u'={yt[1]/yt[0]:.2e}, u'''={yt[3]/yt[0]:.2e};  u(tip) = {yt[0]:.6e}")
json.dump({"beta": beta, "J1": J1.tolist(), "J2": J2.tolist(), "lamt": LT, "utip": yt[0]}, open('psi3beta.json', 'w'))
