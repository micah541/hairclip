"""T9 (part): certified base cone jets (H3 K4, H5 cone cylinders) and the chart factor KCH, in arb, from the corrected T6 jet
hulls (p2/t6/jets; each hull is a sup over its sub-cell and over the whole slab), slabs 177-180 (the good-time window).
Cone cylinders (Lemma O, H5): centre radius r_c in [200, 380], L = r_c/4, X = 3 (half-width 0.75 r_c), graph over the
  centre tangent: h'' = k sec^3, h''' = sec^4 (k_s + 3 k^2 t), h'''' = sec^5 (k_ss + 10 k k_s t + 12 k^3 t^2 + 3 k^3 sec^2);
  |psi| <= sup|k| * (0.75 r_c + 2) * 1.01; K_j = L^(j-1) sup|h^(j)|, N3 = L^1.5 ||h'''||_L2(|x|<=XL) <= sqrt(6) K3.
  r = |Gamma~| = sqrt(p^2 + q^2) per cell; any frame with |x_c| <= 0.5, |lambda - 1| <= 1e-3 is covered by the slack 1.
  Error parts at Ybar (C4 at the cylinder's inner end y = 0.25 r_c - X0, constants as in boxes.py), with D2(tau*) from camp.
KCH (LEMMA_F C3): the error slope over the S1 ray is g e' + (nonlinear), g = sec^2(phi), phi = theta~ - Theta; with
  D = y d/dy = r d/dx, r = y_ray/cos(phi) (y_ray' = cos phi), D^j(g e') = sum_m c_jm w_m (sympy, exact), and
  sqrt(E_j[D^j(g e')]) <= sum_m sup|c_jm| y_T8^-m (y_T8/y_ray) sqrt(cos phi) sqrt(F_m)  (T_m >= 28 y_T8^(2m-2) off the flattening).
  KCH = max over j, m >= 1 of that factor / S(j,m); the m = 0 terms (S(j,0) = 0 for j >= 1) are returned separately (c0_j).
  y_ray/y_T8 from q - x monotone on a cell (q' = 1 + k p): q - x in [q_lo - a, q_hi - b].
Inputs used as stated: Theta = THINF = 1.8617364190 +/- 1e-3 (CAPD gives int_tip^8 k_S1; the tail int_8^oo k_S1 = 1.13e-4 is the
  float value, covered by the allowance), X0 = 2.2490776515 (T8's y = x - X0).
usage: python3 -I t9.py AMP_JSON"""
import sys, os, json, math
import numpy as np
import sympy as sp
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
sys.path.insert(0, '/home/claude/s1spec/step4/p2'); import curvjets as CJ
B0 = json.load(open(os.environ.get('N3JSON', '/home/claude/s1spec/step4/p2/t8se/constants/n3_fin.json')))['args']['b0']   # |e'| box (n3boxes.py)
JD = '/home/claude/s1spec/step4/p2/t6/jets'
SLABS = range(177, 181)
X0 = arb('2.2490776515'); TH = arb('1.8617364190179404', '1e-3')
amp = json.load(open(sys.argv[1])); hb = amp['handoff']['at_Ybar']
sF = [arb(repr(v)) for v in hb['sqrtF']]; F = [v * v for v in sF]; D2 = arb(repr(hb['D2_goodtime']))
U = lambda v: math.nextafter(float(v.upper()), math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))
def absb(a): return arb(U(abs(a)))

# ---- cells over the window slabs
cells = []
for s in SLABS:
    J = RD.jets(s, d=JD)
    for e in range(J['x'].shape[0]):
        for c in range(J['x'].shape[1]):
            g = lambda nm: iv(*J[nm][e, c])
            cells.append(dict(slab=s, a=float(J['x'][e, c, 0]), b=float(J['x'][e, c, 1]), k0=g('k0'), k1=g('k1'), k2=g('k2'), k3=g('k3'),
                              th=g('theta'), p=g('p'), q=g('q')))
for c in cells: c['r'] = (c['p'] ** 2 + c['q'] ** 2).sqrt()

# ---- error pointwise bounds at y (C4; T_0 >= 28/(1.5y)^2, T_m >= 28 (0.8828 y/1.5)^(2m-2))
rr = arb('1.5') / arb('0.8828'); c2 = rr * rr; c3 = c2 * c2
def err(y):
    y = arb(y)
    e2 = ((arb(3) / 28 * F[1] + 2 * rr / 28 * (F[1] * F[2]).sqrt()) / y).sqrt()
    e3 = ((3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt()) / y ** 3).sqrt()
    e4 = (3 * c3 * F[3] / (28 * y ** 5) + 2 * (c2 * c3).sqrt() / 28 * (F[3] * D2).sqrt() / y ** 3).sqrt()
    return e2, e3, e4

# ---- cone jets
out = dict(cone=[])
worst = dict(K2=0.0, K3=0.0, K4=0.0, N3=0.0, slope=0.0)
for rlo in range(200, 380, 10):
    rhi = rlo + 10; L = arb(rhi) / 4
    lo, hi = 0.25 * rlo - 1, 1.75 * rhi + 1
    sel = [c for c in cells if float(c['r'].upper()) >= lo and float(c['r'].lower()) <= hi]
    k0 = max((absb(c['k0']) for c in sel), key=U); k1 = max((absb(c['k1']) for c in sel), key=U); k2 = max((absb(c['k2']) for c in sel), key=U)
    e2, e3, e4 = err(max(0.25 * rlo - 1 - 2.25, 40.0))
    k3 = max((absb(c['k3']) for c in sel), key=U)
    yhi = arb(1.75 * rhi + 1) ; E0 = arb('1e-5') + F[0].sqrt() * (yhi ** 3 / 84).sqrt()   # |e| <= T5 tip box + core integral (LEMMA_F C4)
    K0, K1, K2_ = (arb(v) for v in CJ.kmajor([arb(0, U(k0)), arb(0, U(k1)), arb(0, U(k2)), arb(0, U(k3))], [U(E0), B0, U(e2), U(e3), U(e4)], 3))
    psi = K0 * (arb('0.75') * rhi + 2) * arb('1.01'); t = psi.tan(); sec = 1 / psi.cos()
    assert U(sec) <= 1.01                                  # true half-length <= 1.01 (0.75 r + 2), so psi <= K0 * that
    h2 = K0 * sec ** 3; h3 = sec ** 4 * (K1 + 3 * K0 ** 2 * t)
    h4 = sec ** 5 * (K2_ + 10 * K0 * K1 * t + 12 * K0 ** 3 * t * t + 3 * K0 ** 3 * sec * sec)
    row = dict(r=[rlo, rhi], L=U(L), slope=U(t), K2=U(L * h2), K3=U(L ** 2 * h3), K4=U(L ** 3 * h4), N3=U(arb(6).sqrt() * L ** 2 * h3),
               base_k=[U(k0), U(k1), U(k2)], err=[U(e2), U(e3), U(e4)])
    out['cone'].append(row)
    for k in worst: worst[k] = max(worst[k], row[k])
out['cone_max'] = worst
h3k4 = out['cone'][0]                                       # r_c in [200, 210]: y in [50, 350] (H3's K4)
print(f"cone (r_c 200-380, base + error at Ybar): K2 <= {worst['K2']:.3e}, K3 <= {worst['K3']:.3e}, K4 <= {worst['K4']:.3e}, "
      f"N3 <= {worst['N3']:.3e}, slope <= {worst['slope']:.4f}  (H5: 0.12, 0.12, 0.12, 0.14);  H3 K4 (r_c = 200): {h3k4['K4']:.3e} (<= 0.06)")

# ---- KCH: exact chain/Leibniz coefficients
x = sp.symbols('x'); yv, ph = sp.symbols('y phi'); kk = sp.symbols('k0:4'); W = sp.symbols('W0:5')
def ddx(expr):
    out_ = sp.diff(expr, yv) * sp.cos(ph) + sp.diff(expr, ph) * kk[0]
    out_ += sum(sp.diff(expr, kk[i]) * kk[i + 1] for i in range(3)) + sum(sp.diff(expr, W[i]) * W[i + 1] for i in range(4))
    return out_
r_ = yv / sp.cos(ph); g_ = 1 / sp.cos(ph) ** 2
expr = g_ * W[0]; coef = {}
for j in range(4):
    poly = sp.expand(expr)
    coef[j] = {m: sp.simplify(poly.coeff(W[m])) for m in range(j + 1)}
    expr = r_ * ddx(expr)
Yv, rh = sp.symbols('Y rho', positive=True)                # y_ray = Y rho, Y = y_T8: divide by Y^m symbolically (no wrapping)
fns = {(j, m): sp.lambdify((Yv, rh, ph, kk[0], kk[1], kk[2]), sp.expand(coef[j][m].subs(yv, Yv * rh) / Yv ** m),
       modules=[{'cos': lambda a: a.cos(), 'sin': lambda a: a.sin(), 'tan': lambda a: a.tan(), 'sec': lambda a: 1 / a.cos()}, 'math'])
       for j in range(4) for m in range(j + 1)}
S = {(1, 1): 1, (2, 1): 1, (2, 2): 1, (3, 1): 1, (3, 2): 3, (3, 3): 1, (0, 0): 1}
kch = 0.0; c0 = {}; worstw = None
ylo_w, yhi_w = 50.0, 1.5 * 253
for c in cells:
    if c['b'] < ylo_w + 2.25 - 2 or c['a'] > yhi_w + 2.25 + 2: continue
    a_, b_ = arb(c['a']), arb(c['b'])
    phi = c['th'] - TH
    xa = iv(c['a'], c['b'])
    qmx = iv(float((c['q'].lower() - a_).lower()), float((c['q'].upper() - b_).upper()))   # q - x (monotone in x)
    yT8 = xa - X0
    ratio = phi.cos() * (1 + (qmx + X0) / yT8) - c['p'] * phi.sin() / yT8                # y_ray / y_T8
    yr = yT8 * ratio
    for j in range(4):
        for m in range(j + 1):
            v = fns[(j, m)](yT8, ratio, phi, c['k0'], c['k1'], c['k2'])
            # m >= 1: T_m >= 28 y_T8^(2m-2) -> factor sup|c_jm| sqrt(cos)/(ratio y_T8^m);  m = 0: T_0 >= 28/y_T8^2 -> sup|c_j0| sqrt(cos)/ratio
            b = absb(v) * phi.cos().sqrt() / ratio
            if (j, m) in S:
                fac = U(b / S[(j, m)])
                if fac > kch: kch, worstw = fac, (c['slab'], c['a'], j, m)
            else:
                c0[j] = max(c0.get(j, 0.0), U(b))
nl = arb('1.0001')                                          # nonlinear + reparametrisation corrections (|k e|, e'^2 <= 2.5e-3; see report)
out['KCH'] = U(arb(kch) * nl); out['KCH_where'] = worstw; out['c0_chain'] = dict(c0)
out['c0'] = {0: 0.002, **{j: math.nextafter(v + 0.01, math.inf) for j, v in c0.items()}}   # + reparametrisation allowance (as used since the v2 run)
out['phi_max'] = max(U(abs(c['th'] - TH)) for c in cells if ylo_w + 2.25 - 2 <= c['b'] and c['a'] <= yhi_w + 2.25 + 2)
print(f"KCH (linear chain/Leibniz factor x 1.0001) <= {out['KCH']:.5f} at (slab, x, j, m) = {worstw}; sup|phi| = {out['phi_max']:.4f}; "
      f"m = 0 terms c0_j = " + ', '.join(f'{j}: {v:.3e}' for j, v in sorted(c0.items())))
out['args'] = dict(amp=sys.argv[1], window=amp['handoff']['window_slabs'], SLABS=list(SLABS), b0=B0, env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'H1TAG', 'H1NSUB', 'H1REF', 'TREF', 'PDOT', 'T9TC', 'T9IFT', 'WINOUT', 'SUBW', 'T9JSON', 'TOPLOSS', 'N3JSON') if k in os.environ})
json.dump(out, open(os.environ.get('T9OUT', '/home/claude/s1spec/step4/p2/out/t9_cone_kch.json'), 'w'), indent=1)
