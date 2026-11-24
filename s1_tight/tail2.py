"""Tightened conical-end enclosure at sigma = S (replaces Prop. 2's first-order state by the EXACT linear solution).
Full tail solution:  Y = (k, k', k'') = eps * Yhat(sigma) + Y_nl(sigma),
  Yhat = exact decaying solution of k''' + (sigma/4)k = 0 (Laplace integral, laplace.py), normalised so that its
         decaying Z-coordinate at S is 1 (same normalisation of eps as in step 1: Z_0(S) = eps);
  Y_nl = P Z_nl, Z_nl the Lyapunov-Perron fixed point of the Z-system with forcing  -(C/3 sigma) Z_nl + (w/3mu^2) N(Y),
         Z_nl,0(S) = 0,  Z_nl,1/2 -> 0 at oo.
Bounds (same weighted norms and constants c, kappa, nu0, nu1 as step 1; tail.py):
  ||Zhat_0 - Phi_0||_0 <= r0L, ||Zhat_1||_1 <= r1L  (linear LP, N = 0, per unit eps);  Rhat = 1 + r0L + 2 r1L;
  R = eps_max Rhat + rn0 + 2 rn1;  rn0 >= (c/2) rn1 + nu0 R^3,  rn1 >= kappa (rn0 + rn1) + nu1 R^3;  L < 1.
The linear part of the step-1 bound (kappa * eps) is gone: the remainder is cubic."""
import sys
sys.path.insert(0, __import__('os').environ.get('S1DIR', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 's1_cap')))
from flint import arb, acb, ctx
import tail as T
from laplace import K
ctx.prec = 200
def linear_radii(S, inflate=arb("1.0001")):
    Kc = T.constants(S)
    r1 = Kc["kappa"] / (1 - Kc["kappa"]); r0 = T.C_ABS / 2 * r1
    for _ in range(60):
        r1n = Kc["kappa"] * (1 + r0 + r1) * inflate; r0n = T.C_ABS / 2 * r1n * inflate
        r0, r1 = arb(r0n.upper()), arb(r1n.upper())
    ok = bool(Kc["kappa"] * (1 + r0 + r1) <= r1 and T.C_ABS / 2 * r1 <= r0)
    return r0, r1, ok
def nl_radii(S, eps_max, inflate=arb("1.0001")):
    Kc = T.constants(S); r0L, r1L, okL = linear_radii(S); assert okL
    Rhat = 1 + r0L + 2 * r1L; e = arb(eps_max)
    rn0, rn1 = arb(0), arb(0)
    for _ in range(80):
        R = e * Rhat + rn0 + 2 * rn1
        rn1n = (Kc["kappa"] * (rn0 + rn1) + Kc["nu1"] * R ** 3) * inflate
        rn0n = (T.C_ABS / 2 * rn1n + Kc["nu0"] * R ** 3) * inflate
        rn0, rn1 = arb(rn0n.upper()), arb(rn1n.upper())
    R = e * Rhat + rn0 + 2 * rn1
    ok0 = (T.C_ABS / 2 * rn1 + Kc["nu0"] * R ** 3) <= rn0
    ok1 = (Kc["kappa"] * (rn0 + rn1) + Kc["nu1"] * R ** 3) <= rn1
    Lc = T.C_ABS / 2 + Kc["kappa"] + 6 * (Kc["nu0"] + Kc["nu1"]) * R ** 2
    return dict(r0L=r0L, r1L=r1L, Rhat=Rhat, rn0=rn0, rn1=rn1, R=R, L=Lc, verified=bool(ok0 and ok1 and Lc < 1), **Kc)
def Pmat(S):
    """P(S) = diag(1, mu, mu^2) V,  V_ij = (-w_j)^i"""
    m = T.mu(S); ws = [acb(1), T.W1, T.W2]
    return [[(acb(-1) * ws[j]) ** i * m ** i for j in range(3)] for i in range(3)]
def solve3(Pm, y):
    """solve P z = y (3x3 complex) by Cramer"""
    def det(M): return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    D = det(Pm); out = []
    for c in range(3):
        M = [[(y[i] if j == c else Pm[i][j]) for j in range(3)] for i in range(3)]
        out.append(det(M) / D)
    return out
def yhat(S):
    """exact linear solution at S normalised by Zhat_0(S) = 1: returns (k, k', k''), Zhat(S)"""
    Ks = [acb(K(S, j)) for j in range(3)]
    Z = solve3(Pmat(S), Ks)
    cn = 1 / Z[0]
    return [Ks[j] * cn for j in range(3)], [z * cn for z in Z]
if __name__ == "__main__":
    S = 8
    Y, Z = yhat(S)
    print("Yhat(S) (k,k',k''):", [y.real for y in Y], " imag parts:", [y.imag for y in Y])
    print("Zhat(S):", Z, "   compare eps*gamma1 (step 1):", T.gamma1(S))
    r = nl_radii(S, 1.62e-4 * (1 + 1e-12))
    print({k: float(v.upper()) if hasattr(v, 'upper') else v for k, v in r.items() if k in ('r0L', 'r1L', 'Rhat', 'rn0', 'rn1', 'R', 'L', 'verified', 'kappa')})
    old = T.radii(arb(8), arb(1.62e-4 * (1 + 1e-12)))
    print("step-1 delta1 =", float(old['delta1'].upper()), "   new state uncertainty rn1 =", float(r['rn1'].upper()))
