"""Float: ||L* f_j|| for f_j = chi psi_j, and the residual after subtracting the best combination of f1, f2, f3."""
import numpy as np, sympy as sp
from math import comb
from scipy.integrate import solve_ivp, cumulative_trapezoid
import gram_float as GF
from gram_float import coef, SIG0
import gen_hs3 as G3
STIP = GF.STIP; LAMT = G3.LAMT
Xg = np.linspace(0, SIG0 - STIP, 60001); sig = Xg + STIP
prof = np.array([GF.prof(s) for s in sig]); p_, q_, k_, m_ = prof.T
rho = np.exp(-(1 + Xg**2) ** (2 / 3) / 4)
u = np.clip((SIG0 - sig) / 3, 0, 1); chi = np.where(sig <= 29, 1.0, u**4 * (35 - 84*u + 70*u*u - 20*u**3))
# sympy: L* in terms of psi derivatives, and derivatives of p and of cos(theta - theta_tip) via the profile ODE
X, k, m, q, p = sp.symbols('X k m q p'); psi = sp.Function('psi')(X)
r = sp.exp(-(1 + X**2) ** sp.Rational(2, 3) / 4)
LT = -sp.diff(r*psi, X, 4) - k**2*sp.diff(r*psi, X, 2) + (k*m + q/4)*sp.diff(r*psi, X) + r*psi/2
Ls = sp.expand(sp.simplify(LT / r))
ds = [sp.Derivative(psi, (X, n)) if n else psi for n in range(5)]
cS = [sp.lambdify((X, k, m, q, p), Ls.coeff(ds[n]) if n else Ls.subs({ds[n]: 0 for n in range(1, 5)}).coeff(psi), 'numpy') for n in range(5)]
Cs, Ss = sp.symbols('Cs Ss'); rules = {p: -k*q, q: 1 + k*p, k: m, m: p/4, Cs: k*Ss, Ss: -k*Cs}
D = lambda e: sp.expand(sum(sp.diff(e, v) * rv for v, rv in rules.items()))
def ders(e):
    o = [e]
    for _ in range(4): o.append(D(o[-1]))
    return [sp.lambdify((k, m, q, p, Cs, Ss), x, 'numpy') for x in o]
th = cumulative_trapezoid(k_, sig, initial=0.0)      # theta - theta_tip
Cv, Sv = np.cos(th), np.sin(th)                       # psi2 = cos(theta - theta_tip) = C_tip C + S_tip S in the run's variables
dp, dc = ders(p), ders(Cs)
ps1 = [f(k_, m_, q_, p_, Cv, Sv) * np.ones_like(Xg) for f in dp]
ps2 = [f(k_, m_, q_, p_, Cv, -Sv) * np.ones_like(Xg) for f in dc]   # cos(th): d/dsigma cos = -k sin  => use (Cs,Ss)=(cos th, -sin th)
uu = sp.symbols('u'); Sm = uu**4 * (35 - 84*uu + 70*uu**2 - 20*uu**3)
inz = (sig > 29) & (sig < SIG0)
chid = [chi] + [np.where(inz, sp.lambdify(uu, sp.diff(Sm, uu, j))(u) * (-1/3)**j, 0.0) for j in range(1, 5)]
coefs = [cS[n](Xg, k_, m_, q_, p_) * np.ones_like(Xg) for n in range(5)]
def Lstar_chi(ps):
    tot = np.zeros_like(Xg)
    for n in range(5):
        for j in range(n + 1):
            tot += coefs[n] * comb(n, j) * chid[j] * ps[n - j]
    return tot
ip = lambda a, b: np.trapezoid(rho * a * b, Xg)
f1, f2 = chi * ps1[0], chi * ps2[0]
# f3 (float) from the adjoint solution
def fw(tau, Y):
    c = coef(SIG0 - tau, -GF.OMEGA); kk, mm, qq, P1 = c['k'], c['m'], c['q'], c['P1']
    M = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [(0.5 - LAMT), (kk*mm + qq/4), -kk*kk, 0]])
    return -(M @ Y) - P1 * Y
sol = solve_ivp(fw, [0, SIG0 - STIP], np.array(G3.wstar_init()), method='DOP853', rtol=1e-11, atol=1e-18, dense_output=True)
f3 = chi * np.array([sol.sol(SIG0 - s)[0] for s in sig])
for name, ps, f in (("f1", ps1, f1), ("f2", ps2, f2)):
    Lf = Lstar_chi(ps)
    B = np.vstack([f1, f2, f3]); Gm = np.array([[ip(a, b) for b in B] for a in B]); bv = np.array([ip(a, Lf) for a in B])
    mu = np.linalg.solve(Gm, bv); Res = Lf - B.T @ mu
    mu2 = np.linalg.solve(Gm[:2, :2], bv[:2]); Res2 = Lf - B[:2].T @ mu2
    print(f"{name}: ||f|| = {np.sqrt(ip(f, f)):.4f}  ||L* f|| = {np.sqrt(ip(Lf, Lf)):.4f}   minus best f1,f2: {np.sqrt(ip(Res2, Res2)):.4f} (mu={mu2})   minus best f1,f2,f3: {np.sqrt(ip(Res, Res)):.4f}")
print("check: <f1, f1> =", ip(f1, f1), " <psi1, f1> =", ip(ps1[0], f1), "  <psi2,f2>=", ip(ps2[0], f2), " <psi1,f3>=", ip(ps1[0], f3), " <psi2,f3>=", ip(ps2[0], f3))
