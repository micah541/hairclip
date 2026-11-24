"""Float comparison of tip-bump choices: ||(L*-lam) eta||, ||(L-lam) eta||, ||eta|| in L^2(rho) on the half line,
and the resulting bound |a| n1 + |b| n3 with the certified-size tip data |ws1|, |ws3| <= 7.2e-7."""
import numpy as np, sympy as sp
import gram_float as GF
import gen_hs3 as G3
STIP = GF.STIP; T0 = G3.T0; c_ = (8 - STIP) / T0; LAMT = G3.LAMT
X, k, m, q, p, lam = sp.symbols('X k m q p lam')
psi = sp.Function('psi')(X)
r = sp.exp(-(1 + X**2) ** sp.Rational(2, 3) / 4)
wv = r * psi
LT = -sp.diff(wv, X, 4) - k**2 * sp.diff(wv, X, 2) + (k * m + q / 4) * sp.diff(wv, X) + wv / 2
Lstar = sp.simplify(LT / r) - lam * psi
Lop = -sp.diff(psi, X, 4) - k**2 * sp.diff(psi, X, 2) - (5 * k * m + q / 4) * sp.diff(psi, X) + (sp.Rational(1, 4) - 3 * m**2 - k * p) * psi - lam * psi
s = sp.symbols('s')
def run(zeta, Lb, W1=7.2e-7, W3=6.8e-7):
    Xg = np.linspace(0, c_ * Lb, 6001)
    prof = np.array([GF.prof(x + STIP) for x in Xg])
    rho = np.exp(-(1 + Xg**2) ** (2 / 3) / 4)
    sv = X / (c_ * Lb); z = zeta.subs(s, sv)
    e1 = Lb * sv * z; e3 = Lb**3 * sv**3 * z / 6
    out = []
    for e in (e1, e3):
        vals = []
        for op in (Lstar, Lop):
            fn = sp.lambdify((X, k, m, q, p, lam), op.subs(psi, e).doit(), 'numpy')
            v = fn(Xg, prof[:, 2], prof[:, 3], prof[:, 1], prof[:, 0], LAMT) * np.ones_like(Xg)
            vals.append(np.sqrt(np.trapezoid(rho * v**2, Xg)))
        ef = sp.lambdify(X, e, 'numpy')(Xg)
        vals.append(np.sqrt(np.trapezoid(rho * ef**2, Xg)))
        out.append(vals)
    d3 = float(sp.diff(e1, X, 3).subs(X, 0))
    a = c_ * W1; b = c_**3 * (W3 + W1) + abs(d3) * c_**3 * W1
    return out, d3, a * out[0][0] + b * out[1][0]
for name, zeta, Lb in [("(1-s^2)^4", (1 - s**2)**4, 2), ("(1-s^2)^4", (1 - s**2)**4, 3), ("(1-s^2)^4", (1 - s**2)**4, 4),
                       ("(1-s^2)^5", (1 - s**2)**5, 3), ("(1-s^2)^6", (1 - s**2)**6, 4)]:
    (o1, o3), d3, bound = run(zeta, Lb)
    print(f"{name} L={Lb}: eta1 [L*:{o1[0]:.3g} L:{o1[1]:.3g} |.|:{o1[2]:.3g}]  eta3 [L*:{o3[0]:.3g} L:{o3[1]:.3g} |.|:{o3[2]:.3g}]  eta1'''(0)={d3:+.3f}  r3 bump bound {bound:.2e}")
