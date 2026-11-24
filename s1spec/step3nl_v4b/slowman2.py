"""Third order of the slow manifold u = c psi3 + c^2 h2 + c^3 h3 + ... (numerical; norms are true L^2(rho) norms)."""
import numpy as np, scipy.linalg as sl
import slowman as S          # reruns the second-order computation and exposes its objects
h = S.h; sq = np.sqrt(h)
nr = lambda e: sq * np.sqrt(e @ e)                       # true L^2(rho) norm of an even r-scaled vector
eps = 2e-3
def Q(f):                                                 # ~ N2(f, f)
    return (S.R0(eps * f) + S.R0(-eps * f)) / (2 * eps**2)
def N2b(f, g): return (Q(f + g) - Q(f - g)) / 4
p3, h2 = S.psi[2], S.h2
cut = S.cut
F3 = S.from_fun((2 * N2b(p3, h2) + S.N3) * cut)
F3 = F3 - S.R @ (S.Lf.T @ F3) - 2 * S.g2 * S.h2e           # P[...] - 2 g2 h2
g3full = S.Lf[:, 2] @ S.from_fun((2 * N2b(p3, h2) + S.N3) * cut)
h3e = np.linalg.solve(S.Ae - 3 * S.lams[2] * np.eye(len(S.Ae)), -F3)
print(f"true norms:  ||psi3|| = {nr(S.R[:,2]):.2f}   ||P N2|| = {nr(S.PN2):.2f}   ||h2|| = {nr(S.h2e):.2f}")
print(f"cubic:  c^3 coefficient of c' = {g3full:+.4f};  ||F3|| = {nr(F3):.2f};  ||h3|| = {nr(h3e):.2f}")
h3 = S.to_fun(h3e)
for xx in [0, 2, 4, 6, 8, 10]:
    j = int(round(xx / h)); print(f"   x={xx:4.1f}: psi3={p3[j]:+8.3f}  h2={h2[j]:+9.3f}  h3={h3[j]:+10.3f}")
for c in [0.01, 0.02, 0.04]:
    print(f"c={c}: |c|psi3 = {c*nr(S.R[:,2]):.3e}, c^2 h2 = {c*c*nr(S.h2e):.3e}, c^3 h3 = {c**3*nr(h3e):.3e}; "
          f"rate {S.lams[2] + S.g2*c:+.4f}; a = {S.a2*c*c:+.2e}, b = {S.b2*c*c:+.2e}")
