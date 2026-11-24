"""Float: sign structure of the inner cutoff commutator on the strip sigma in [29,32] after integrating by parts down to
w^2 and w'^2:  -<[L,chi]w, chi w>_rho = int A w^2 + B w'^2  (k, m ~ 0 on the strip).  A/rho, B/rho printed."""
import sympy as sp, numpy as np
s, X = sp.symbols('s X')          # s = sigma, X = x = sigma - sigma_tip
STIP = -2.2490776510
u = (32 - s) / 3
S = u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)
chi = S
c1, c2, c3, c4 = [sp.diff(chi, s, j) for j in (1, 2, 3, 4)]
xx = s - STIP
phi = (1 + xx**2) ** sp.Rational(2, 3) / 4
f1, f2, f3 = [sp.diff(phi, s, j) for j in (1, 2, 3)]
def d(e, n=1): return sp.diff(e, s, n)
g1, g2, g3, g4 = chi * c1, chi * c2, chi * c3, chi * c4
R3 = lambda g: d(g, 3) - 3 * f1 * d(g, 2) + 3 * (f1**2 - f2) * d(g) + (3 * f1 * f2 - f1**3 - f3) * g
R2 = lambda g: d(g, 2) - 2 * f1 * d(g) + (f1**2 - f2) * g
R1 = lambda g: d(g) - f1 * g
q = s
A = -2 * R3(g1) + 3 * R2(g2) - 2 * R1(g3) + g4 + q / 4 * g1
B = 6 * (R1(g1) - g2)
Af, Bf = sp.lambdify(s, A, 'numpy'), sp.lambdify(s, B, 'numpy')
ss = np.linspace(29, 32, 3001)
a, b = Af(ss), Bf(ss)
print("A/rho: max %.3f at sigma %.3f, min %.3f;  part without drift: max %.3f" % (a.max(), ss[a.argmax()], a.min(),
      (a - sp.lambdify(s, q / 4 * g1, 'numpy')(ss)).max()))
print("B/rho: max %.3f at sigma %.3f, min %.3f" % (b.max(), ss[b.argmax()], b.min()))
for sv in (29.0, 29.5, 30.0, 30.5, 31.0, 31.5, 32.0):
    print(f"  sigma {sv}: A/rho {Af(sv):+.3f}  B/rho {Bf(sv):+.3f}  chi {float(chi.subs(s, sv)):.3f}")
