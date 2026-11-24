"""Offline check of p2/t9chart.py (not part of the certified chain; sympy is used here only).
(1) Linear part: the first bracket of (*) plus tan(phi) dQ^b_j/du eta_0 (second bracket at first order) must equal
    D^j(sec^2 phi e' + sec^2 phi tan phi k e) (review B5's first variation), D = (y/cos phi) d/dx, derived here with sympy.
(2) Remainder: for random jets in the boxes, |Q^t_j(eta) - Q^b_j - sum l_i eta_i| <= sum r_i |eta_i| (plain arb evaluation of the same
    Taylor pipeline at numeric jets).
usage: python3 -I t9chart_check.py"""
import sys, random
import sympy as sp
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/p2')
import t9chart as TC

# a point with S1-ray-like data (values of slab 180 near x = 200, mid-points)
y0 = arb('197.75'); phi0 = arb('0.0712'); k = [arb('-1.83e-4'), arb('9.1e-7'), arb('-9.0e-9'), arb('1.3e-10')]; ratio = arb('0.99731')
ph, kh, S, C, Yb = TC.base_series(y0, phi0, k, ratio)
Qb = TC.base_Q(S, C, Yb)

# (1) sympy reference, in x (not u)
y, f = sp.symbols('y phi'); K = sp.symbols('k0:5'); W = sp.symbols('W0:5'); E0 = sp.symbols('E0')
def ddx(z):
    out = sp.diff(z, y) * sp.cos(f) + sp.diff(z, f) * K[0] + sum(sp.diff(z, K[i]) * K[i + 1] for i in range(4))
    return out + sum(sp.diff(z, W[i]) * W[i + 1] for i in range(4)) + sp.diff(z, E0) * W[0]
z = sp.sec(f) ** 2 * W[0] + sp.sec(f) ** 2 * sp.tan(f) * K[0] * E0
F_ = lambda a: sp.Float(a.mid().str(35, radius=False), 40)
num = {y: F_(ratio * y0), f: F_(phi0), **{K[i]: F_(k[i]) for i in range(4)}}
TC.BOX[:] = [arb('1e-6'), arb('1.5e-3'), arb('3e-4'), arb('1e-2'), None]
Qt, Mt, Yt = TC.true_Q(kh, S, C, Yb)
ok = True
for j in range(4):
    zz = sp.expand(z)
    ref = {m: float(zz.coeff(W[m]).subs(num)) for m in range(j + 1)}; refE = float(zz.coeff(E0).subs(num))
    lin = Qt[j].l
    for m in range(j + 1):
        mine = float((lin.get(m + 1, arb(0)) * y0.mid() ** m).mid())
        rel = abs(mine - ref[m]) / max(abs(ref[m]), 1e-300)
        ok &= rel < 1e-12; print(f'j={j} w_{m}: chain {ref[m]: .15e}  taylor {mine: .15e}  rel {rel:.1e}')
    mineE = float(((lin.get(0, arb(0)) + phi0.mid().tan() * Qb[j][1]) / y0.mid()).mid())
    rel = abs(mineE - refE) / abs(refE); ok &= rel < 1e-12
    print(f'j={j} e  : chain {refE: .15e}  taylor {mineE: .15e}  rel {rel:.1e};  Q^b_j {float(Qb[j][0].mid()): .6e} vs E.c {float(Qt[j].c.mid()): .6e}')
    z = (y / sp.cos(f)) * ddx(z)

# (2) remainder bounds at random jets (eta_4 has no box: drawn of size 1)
random.seed(1); worst = 0.0
for trial in range(200):
    eta = [arb(random.uniform(-1, 1)) * TC.BOX[i] for i in range(4)] + [arb(random.uniform(-1, 1))]
    vals, _, _ = TC.true_Q(kh, S, C, Yb, eta)
    for j in range(4):
        d = vals[j].c - Qt[j].c - sum((Qt[j].l.get(i, arb(0)) * eta[i] for i in range(5)), arb(0))
        bnd = sum((Qt[j].r.get(i, arb(0)) * abs(eta[i]) for i in range(5)), arb(0))
        ratio_ = float(abs(d).upper()) / max(float(bnd.lower()), 1e-300)
        worst = max(worst, ratio_)
        ok &= float(abs(d).lower()) <= float(bnd.upper())
print(f'remainder: max |NL| / bound over 200 draws = {worst:.3f} (must be <= 1)')
print('r_i / |l_i| at j = 3:', {i: f'{float((Qt[3].r.get(i, arb(0)) / abs(Qt[3].l[i])).mid()):.2e}' for i in Qt[3].l if Qt[3].l[i] != 0})
print('ALL OK' if ok else 'MISMATCH')
