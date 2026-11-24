"""Float estimate of the residual r3 = (L* - lamt) f3 for f3 = chi F - a eta1 - b eta3 (L^2(rho), half line).
L*psi = rho^{-1} L^T(rho psi),  L^T w = -w'''' - k^2 w'' + (k m + q/4) w' + w/2."""
import numpy as np, sympy as sp, json
from scipy.integrate import solve_ivp
import gram_float as GF
from gram_float import coef, MLT, SIG0
import gen_hs3 as G3
STIP = GF.STIP; T0 = G3.T0; c_ = (8 - STIP) / T0; LAMT = G3.LAMT
# --- F = w/rho on the half line from the scaled-frame solution
def f(tau, Y):
    c = coef(SIG0 - tau, -GF.OMEGA); return -(MLT(c) @ Y) - c['P1'] * Y
sol = solve_ivp(f, [0, SIG0 - STIP], np.array(GF.wstar_init()), method='DOP853', rtol=1e-12, atol=1e-16, dense_output=True)
Xg = np.linspace(0, SIG0 - STIP, 40001)
W = np.array([sol.sol(SIG0 - (x + STIP)) for x in Xg])          # ws_i = w^(i)/rho
rho = np.exp(-(1 + Xg**2) ** (2 / 3) / 4)
F = W[:, 0]
# --- symbolic L* applied to an explicit psi(X): returns expression in psi derivs, k, m, q, X
X = sp.symbols('X'); k, m, q, lam = sp.symbols('k m q lam')
psi = sp.Function('psi')(X)
r = sp.exp(-(1 + X**2) ** sp.Rational(2, 3) / 4)
wv = r * psi
LT = -sp.diff(wv, X, 4) - k**2 * sp.diff(wv, X, 2) + (k * m + q / 4) * sp.diff(wv, X) + wv / 2
Lstar = sp.simplify(LT / r) - lam * psi
def apply(expr_fun, Xv, prof):
    """evaluate (L* - lam) on psi = expr_fun (sympy expr in X)"""
    e = Lstar.subs(psi, expr_fun).doit()
    fn = sp.lambdify((X, k, m, q, lam), e, 'numpy')
    return fn(Xv, prof[:, 2], prof[:, 3], prof[:, 1], LAMT)
prof = np.array([GF.prof(x + STIP) for x in Xg])
s = sp.symbols('s')
def bumps(zeta_expr, Lb):
    sv = X / (c_ * Lb)
    z = zeta_expr.subs(s, sv)
    return sv * z * Lb, sv**3 * z * Lb**3 / 6        # scaled so eta1'(0)=1/c, eta3'''(0)=1/c^3
def norm(v, mask=None):
    v = np.where(mask, v, 0.0) if mask is not None else v
    return np.sqrt(np.trapezoid(rho * v**2, Xg))
nf3 = norm(F)
ws1, ws3 = W[0, 1], W[0, 3]
print(f"||F||_rho (half line) = {nf3:.4e};  true tip data ws1 = {ws1:.3e}, ws3 = {ws3:.3e}")
for name, zeta, Lb in [("(1-s^2)^4, L=1", (1 - s**2)**4, 1), ("(1-s^4)^4, L=1", (1 - s**4)**4, 1),
                       ("(1-s^2)^4, L=2", (1 - s**2)**4, 2), ("(1-s^4)^4, L=2", (1 - s**4)**4, 2)]:
    e1, e3 = bumps(zeta, Lb)
    mask = Xg <= c_ * Lb
    n1 = norm(apply(e1, Xg, prof), mask); n3 = norm(apply(e3, Xg, prof), mask)
    d3_e1 = float(sp.diff(e1, X, 3).subs(X, 0))
    print(f"  bump {name}: ||(L*-lam)eta1|| = {n1:.3e}, ||(L*-lam)eta3|| = {n3:.3e}, eta1'''(0) = {d3_e1:+.3f}")
# --- far-field cutoff residual [L*, chi] F on sigma in [29, 32]
sig = Xg + STIP
u = np.clip((SIG0 - sig) / 3.0, 0, 1)
chiS = sp.Function('chi')
# commutator via finite differences of chi*F is fragile; use exact formula: (L*-lam)(chi F) - chi (L*-lam) F, with L* in terms of derivatives
uu = sp.symbols('uu'); Sm = uu**4 * (35 - 84*uu + 70*uu**2 - 20*uu**3)
# derivatives of chi in sigma: chi(sigma) = S((32 - sigma)/3)
dchi = [sp.lambdify(uu, sp.diff(Sm, uu, j), 'numpy') for j in range(5)]
inzone = (u > 0) & (u < 1)
chi_d = [dchi[0](u)] + [np.where(inzone, dchi[j](u) * (-1/3)**j, 0.0) for j in range(1, 5)]
# F derivatives: F = w e^phi; use scaled frame: F^(n) = sum_j binom(n,j) w^(j) (e^phi)^(n-j); (e^phi)^(i) = e^phi * B_i(phi', ...)
Xs = sp.symbols('Xs'); ph = (1 + Xs**2) ** sp.Rational(2, 3) / 4
Bpoly = [sp.lambdify(Xs, sp.simplify(sp.diff(sp.exp(ph), Xs, i) / sp.exp(ph)), 'numpy') for i in range(4)]
Bv = [Bpoly[i](Xg) * np.ones_like(Xg) for i in range(4)]
from math import comb
Fd = [sum(comb(n, j) * W[:, j] * Bv[n - j] for j in range(n + 1)) for n in range(4)]   # F^(n), n<=3, ws_j = w^(j) e^phi
# L* (in psi) has the form -psi'''' + a3 psi''' + a2 psi'' + a1 psi' + a0 psi; commutator only involves psi^(n<=3)
P = sp.Function('P')(X)
coeffs_L = sp.Poly(sp.expand(Lstar.subs(psi, P).doit()).subs({sp.Derivative(P, (X, 4)): sp.Symbol('d4'), sp.Derivative(P, (X, 3)): sp.Symbol('d3'),
                    sp.Derivative(P, (X, 2)): sp.Symbol('d2'), sp.Derivative(P, X): sp.Symbol('d1'), P: sp.Symbol('d0')}),
                    *[sp.Symbol(f'd{i}') for i in range(5)])
aco = {i: sp.lambdify((X, k, m, q, lam), coeffs_L.coeff_monomial(sp.Symbol(f'd{i}')), 'numpy') for i in range(5)}
A = {i: aco[i](Xg, prof[:, 2], prof[:, 3], prof[:, 1], LAMT) * np.ones_like(Xg) for i in range(5)}
# commutator [L*,chi]F = sum_n a_n sum_{j>=1} binom(n,j) chi^(j) F^(n-j)
comm = np.zeros_like(Xg)
for n in range(1, 5):
    for j in range(1, n + 1):
        comm += A[n] * comb(n, j) * chi_d[j] * Fd[n - j]
print(f"far-field cutoff residual ||[L*,chi]F||_rho = {norm(comm):.3e}   (relative {norm(comm)/nf3:.2e})")
json.dump({"normF": nf3}, open("resid_float.json", "w"))
i = np.argmax(np.abs(comm) * np.sqrt(rho))
print("argmax sigma", sig[i], "comm", comm[i], "rho", rho[i], "chi_d", [cd[i] for cd in chi_d], "F^(n)", [fd[i] for fd in Fd], "A", [A[n][i] for n in range(5)])
mask = (sig > 28.9)
print("max|comm| on cutoff zone", np.abs(comm[mask]).max(), " ; max |comm| elsewhere", np.abs(comm[~mask]).max())
