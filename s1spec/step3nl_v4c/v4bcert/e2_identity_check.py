"""Checks of the E2 identity used by the E2 couplings (E2_AUDIT.md section 1, [hand] there):
   (1/2) d/dtau ||w~''||^2_rho = <d_tau w~, D_rho w~>_rho,   D_rho g := rho^-1 (rho g'')'' = g'''' - 2 phi' g''' + (phi'^2 - phi'') g'',
   and the symmetry <f, D_rho g>_rho = <f'', g''>_rho = <D_rho f, g>_rho used for -a<chi psi1, D w~> = -a<D f1, w~>  (no boundary term at x = 0 for even f, g;
   rho = exp(-phi), phi = (1+x^2)^(2/3)/4).
 C1 [sym] D_rho expansion: rho^-1 (rho g'')'' - (g'''' - 2 phi' g''' + (phi'^2 - phi'') g'') == 0 for generic phi, g.
 C2/C3 [sym] rho f'' g'' - rho f D_rho g = d/dx B exactly, B = rho f' g'' - f (rho g'')' (generic phi, f, g), and B(0) = 0 for even f, g, phi
       (even polynomials of degree 10, 10, 6 with symbolic coefficients).
 C4 [num] float quadrature on [0, 40] with the true phi, random even f, g = P(x^2) exp(-c x^2) (6 seeds): |<f,D g> - <f'',g''>| and |<f,D g> - <D f,g>| relative.
 C5 [num] ||D_rho g|| <= ||g''''|| + 2.2 ||g'''|| + 1.2 ||g''|| on the random g, with g supported in x <= 34.25 (cutoff), as used in Dq.
Output out/e2_identity_check.json."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump
import sympy as sp, numpy as np
x = sp.symbols('x'); ph = sp.Function('phi')(x); f = sp.Function('f')(x); g = sp.Function('g')(x)
rho = sp.exp(-ph)
Dg = sp.diff(rho * sp.diff(g, x, 2), x, 2) / rho
exp_ = sp.diff(g, x, 4) - 2 * sp.diff(ph, x) * sp.diff(g, x, 3) + (sp.diff(ph, x) ** 2 - sp.diff(ph, x, 2)) * sp.diff(g, x, 2)
res = {'C1_sym': bool(sp.simplify(Dg - exp_) == 0)}
Lag = rho * sp.diff(f, x, 2) * sp.diff(g, x, 2) - rho * f * Dg
# C3: Lag = d/dx B with B = rho f' g'' - f (rho g'')'  (two integrations by parts)
Bt = rho * sp.diff(f, x) * sp.diff(g, x, 2) - f * sp.diff(rho * sp.diff(g, x, 2), x)
res['C3_total_derivative_form'] = bool(sp.simplify(sp.diff(Bt, x) - Lag) == 0)
# even f, g, phi: f', g''' , phi' odd -> B(0) = rho(0) [f'(0) g''(0) - f(0)(g'''(0) - phi'(0) g''(0))] = 0
a = sp.symbols('a0:6'); b = sp.symbols('b0:6'); cc = sp.symbols('c0:4')
fe = sum(a[i] * x ** (2 * i) for i in range(6)); ge = sum(b[i] * x ** (2 * i) for i in range(6)); pe = sum(cc[i] * x ** (2 * i) for i in range(4))
B0 = Bt.subs({f: fe, g: ge, ph: pe}).doit().subs(x, 0)
res['C3_sym_boundary_at_0'] = bool(sp.simplify(B0) == 0)
# C4, C5 numerics
X = np.linspace(0.0, 40.0, 400001); hh = X[1] - X[0]
phi = (1 + X ** 2) ** (2 / 3) / 4; r = np.exp(-phi)
p1 = X * (1 + X ** 2) ** (-1 / 3) / 3; p2 = (1 + X ** 2) ** (-1 / 3) / 3 - 2 * X ** 2 * (1 + X ** 2) ** (-4 / 3) / 9
def simpson(y): return hh / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())
xs = sp.symbols('xs')
def fun(seed, cut):
    rng = np.random.default_rng(seed); co = rng.normal(size=4); c = rng.uniform(0.02, 0.2)
    e = sum(co[i] * xs ** (2 * i) for i in range(4)) * sp.exp(-c * xs ** 2)
    if cut: e = e * sp.exp(-(xs / 30) ** 12)
    return [sp.lambdify(xs, sp.diff(e, xs, j), 'numpy')(X) for j in range(5)]
worst = dict(sym_FDg_vs_f2g2=0.0, sym_FDg_vs_DfG=0.0, C5_ratio=0.0)
for seed in range(6):
    F_ = fun(seed, True); G_ = fun(100 + seed, True)
    D = lambda J: J[4] - 2 * p1 * J[3] + (p1 ** 2 - p2) * J[2]
    I1 = simpson(r * F_[0] * D(G_)); I2 = simpson(r * F_[2] * G_[2]); I3 = simpson(r * D(F_) * G_[0])
    sc = max(abs(I1), abs(I2), 1e-300)
    worst['sym_FDg_vs_f2g2'] = max(worst['sym_FDg_vs_f2g2'], abs(I1 - I2) / sc); worst['sym_FDg_vs_DfG'] = max(worst['sym_FDg_vs_DfG'], abs(I1 - I3) / sc)
    m = X <= 34.25; nrm = lambda y: np.sqrt(simpson(np.where(m, r * y * y, 0.0)))
    lhs = nrm(D(G_)); rhs = nrm(G_[4]) + 2.2 * nrm(G_[3]) + 1.2 * nrm(G_[2])
    worst['C5_ratio'] = max(worst['C5_ratio'], lhs / rhs)
res['C4_num'] = worst; res['C4_ok'] = bool(worst['sym_FDg_vs_f2g2'] < 1e-8 and worst['sym_FDg_vs_DfG'] < 1e-8); res['C5_ok'] = bool(worst['C5_ratio'] <= 1.0)
res['all_ok'] = all(v for k, v in res.items() if k.endswith('_sym') or k.endswith('_ok') or k == 'C3_total_derivative_form')
res['args'] = args_block([])
print(json.dumps({k: v for k, v in res.items() if k != 'args'}, indent=1)); assert res['all_ok']
print('sha256', dump(os.path.join(OUT, 'e2_identity_check.json'), res))
