"""T9 item 5 (H5 far cylinders), certified composition: X = 2, unit L = 113 (lambda_1 units), centres r >= 376, so x >= 150.
Base: corrected T6 jet hulls of slab 180 (p2/t6/jets; sup over sub-cell and slab, contains W'), cells reaching x >= 148.
Error (Ybar, hand-off good time in W'): C4 pointwise bounds at y = 147 - X0 (decreasing in y) from camp's at_Ybar (sqrt F_j, D2),
|e'| (far core box) and |e| (far end) as in boxes.py.  Truth curvature jets (lambda in [1 - 1e-5, 1]):
  K0 = k0 + e2 + k0^2|e|,  K1 = k1 + e3 + 2 k0 k1 |e| + k0^2 e1,  K2 = k2 + e4 + (2 k1^2 + 2 k0 k2)|e| + 4 k0 k1 e1 + k0^2 e2,
  times (1.003)^(j+1) (arclength / normal-graph factors, |k e| <= 2.4e-3).  Graph over the centre tangent, |psi| <= K0 (X L + 2) 1.01 + e1:
  h'' = K0 sec^3, h''' = sec^4 (K1 + 3 K0^2 t), h'''' = sec^5 (K2 + 10 K0 K1 t + 12 K0^3 t^2 + 3 K0^3 sec^2), ||h'''||_L2(|x|<=2L) <= h''' sqrt(4L).
All arithmetic in arb; floats leave only through nextafter (outward).  usage: python3 -I t9_far.py AMP_JSON CONSTS_JSON"""
import sys, json, math
import numpy as np
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
amp, cons = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
hb = amp['handoff']['at_Ybar']; sF = [arb(repr(v)) for v in hb['sqrtF']]; F = [v * v for v in sF]; D2 = arb(repr(hb['D2_goodtime']))
X0 = arb('2.2490776515'); r = arb('1.5') / arb('0.8828'); c2 = r * r; c3 = c2 * c2
J = RD.jets(180, d='/home/claude/s1spec/step4/p2/t6/jets')
m = J['x'][..., 1] >= 148.0
k0, k1, k2 = (arb(math.nextafter(float(np.abs(J[n][m]).max()), math.inf)) for n in ('k0', 'k1', 'k2'))
y = arb(147) - X0
e2 = ((arb(3) / 28 * F[1] + 2 * r / 28 * (F[1] * F[2]).sqrt()) / y).sqrt()
e3 = ((3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt()) / (y * y * y)).sqrt()
e4 = (3 * c3 * F[3] / (28 * y ** 5) + 2 * (c2 * c3).sqrt() / 28 * (F[3] * D2).sqrt() / (y * y * y)).sqrt()
# |e'| and |e| on the far core at slab 180 (as boxes.py): y_l = ell - X0
i = amp['slabs'].index(180); ps = amp['per_slab']; Y = arb('1e-12')
Fs = [(arb(repr(ps[str(j)]['aF'][i])) * Y + arb(repr(ps[str(j)]['bF'][i]))) ** 2 for j in range(2)]
yl = arb(repr([p for p in cons['per_slab'] if p['slab'] == 180][0]['scalars']['ell'][1])) - X0
e1 = (arb('6.75') / 28 * yl * Fs[0] + arb(3) / 28 * yl * (Fs[0] * Fs[1]).sqrt()).sqrt()
e0 = arb('1.06e-8') + Fs[0].sqrt() * yl ** arb('1.5') / arb(84).sqrt()
K0 = (k0 + e2 + k0 * k0 * e0) * arb('1.003')
K1 = (k1 + e3 + 2 * k0 * k1 * e0 + k0 * k0 * e1) * arb('1.003') ** 2
K2 = (k2 + e4 + (2 * k1 * k1 + 2 * k0 * k2) * e0 + 4 * k0 * k1 * e1 + k0 * k0 * e2) * arb('1.003') ** 3
L = arb(113); XL = 2 * L
psi = K0 * (XL + 2) * arb('1.01') + e1; t = psi.tan(); sec = 1 / psi.cos()
h2 = K0 * sec ** 3; h3 = sec ** 4 * (K1 + 3 * K0 * K0 * t); h4 = sec ** 5 * (K2 + 10 * K0 * K1 * t + 12 * K0 ** 3 * t * t + 3 * K0 ** 3 * sec * sec)
n3 = h3 * (4 * L).sqrt(); thr = dict(h2=8.8e-3, h3=7.0e-5, h4=6.2e-7, slope=0.21, L2=U(L ** arb(-1.5)))
out = dict(base=dict(k0=U(k0), k1=U(k1), k2=U(k2)), err=dict(e0=U(e0), e1=U(e1), e2=U(e2), e3=U(e3), e4=U(e4)),
           h2=U(h2), h3=U(h3), h4=U(h4), slope=U(t), L2=U(n3), thr=thr)
out['margins'] = {k: thr[k] / out[k] for k in ('h2', 'h3', 'h4', 'slope', 'L2')}
json.dump(out, open('/home/claude/s1spec/step4/p2/out/t9_far.json', 'w'), indent=1)
print('item 5 (far cylinders, x >= 150):', {k: f'{out[k]:.3e}' for k in ('h2', 'h3', 'h4', 'slope', 'L2')}, 'margins', {k: f'{v:.1f}' for k, v in out['margins'].items()})
