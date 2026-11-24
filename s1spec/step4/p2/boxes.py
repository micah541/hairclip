"""Nonlinear bootstrap boxes along the certified run (T5 2.5 box B; far box on the core y >= 28), at Ybar.
Far core (LEMMA_F C4 on I of length >= y/3 inside [y/1.5, 1.5y], T_0 >= 28/(1.5y)^2, T_m >= 28 (0.8828 y/1.5)^(2m-2), x <= 1.0804 y):
  |e'|^2      <= (6.75/28) y F0 + (3/28) y sqrt(F0 F1)
  x^2|e''|^2  <= 1.0804^2 y [(3/28) F1 + (2*1.5/0.8828/28) sqrt(F1 F2)]
  x^4|e'''|^2 <= 1.0804^4 y [(3 c2/28) F2 + (2 sqrt(c2 c3)/28) sqrt(F2 F3)],  c2 = (1.5/0.8828)^2, c3 = c2^2
with sqrt F_j <= aF_j Ybar + bF_j per slab (camp), y up to y_l = ell - X0 (sup over the slab).
usage: python3 -I boxes.py AMP_JSON CONSTS_JSON [YBAR]"""
import sys, json
from flint import arb, ctx
ctx.prec = 128
amp, cons = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
Y = arb(sys.argv[3] if len(sys.argv) > 3 else '1e-12')
X0 = arb('2.2490776515'); r = arb('1.5') / arb('0.8828'); c2 = r * r; c3 = c2 * c2; xf = 1 + X0 / 28
ell = {p['slab']: arb(repr(p['scalars']['ell'][1])) for p in cons['per_slab']}
ps = amp['per_slab']; slabs = amp['slabs']
import math
U = lambda v: math.nextafter(float(v.upper()), math.inf)
worst = [(0.0, None)] * 3
for i, s in enumerate(slabs):
    F = [(arb(repr(ps[str(j)]['aF'][i])) * Y + arb(repr(ps[str(j)]['bF'][i]))) ** 2 for j in range(4)]
    y = ell[s] - X0
    b1 = (arb('6.75') / 28 * y * F[0] + arb(3) / 28 * y * (F[0] * F[1]).sqrt()).sqrt()
    b2 = (xf ** 2 * y * (arb(3) / 28 * F[1] + 2 * r / 28 * (F[1] * F[2]).sqrt())).sqrt()
    b3 = (xf ** 4 * y * (3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt())).sqrt()
    for k, b in enumerate((b1, b2, b3)):
        if U(b) > worst[k][0]: worst[k] = (U(b), s)
# review B2: the boxes are those Lemma N3' was certified on (n3far args b0, b1, b2; 4th argument), asserted here
if len(sys.argv) > 4:
    na = json.load(open(sys.argv[4]))['args']; boxn = [na['b0'], na['b1'], na['b2']]
    assert float(Y.upper()) <= float(na['Y']), 'Ybar above the Y of the n3 run'
else:
    boxn = [1e-2, 1e-2, 1e-2]
for k in range(3): assert worst[k][0] <= boxn[k], ('far core box fails', k, worst[k][0], boxn[k])
mg = [math.nextafter(boxn[k] / worst[k][0], 0.0) for k in range(3)]
print(f'Ybar = {Y.str(3)}: far core box vs the N3\' boxes {boxn}:  sup|e\'| <= {worst[0][0]:.4g} (slab {worst[0][1]}), '
      f'sup x|e\'\'| <= {worst[1][0]:.4g} (slab {worst[1][1]}), sup x^2|e\'\'\'| <= {worst[2][0]:.4g} (slab {worst[2][1]}); '
      f'margins {mg[0]:.3f}, {mg[1]:.3f}, {mg[2]:.3f}')
eps, mue, epsI = arb('1e-5'), arb('0.036'), arb('3e-5')
print(f'B <= Ybar/sqrt(eps) = {U(Y / eps.sqrt()):.3e} (eta_B = 1e-4);  tip (T5 2.5: 2e3 Ybar/sqrt(mu_e)) = {U(2000 * Y / mue.sqrt()):.3e} '
      f'(eta_1 = 1e-5);  |e| <= K_e Ybar = {U(arb("890.6") * Y):.3e} on x <= 16.25')
# certified |e| on the core: |e(y)| <= |e(28)| + sqrt(F0) (int_28^y Y_l^2/28)^(1/2) <= 1.06e-8 + sqrt(F0) y^(3/2)/sqrt(84);
# |k~ e| = the top-coefficient perturbation scale (A4 - A4(0) ~ 4 k~ e), with k~ from the corrected T6 jet hulls per cell
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD, numpy as np
wk = (0.0, None, None); we = (0.0, None)
for i, s in enumerate(slabs):
    sF0 = arb(repr(ps['0']['aF'][i])) * Y + arb(repr(ps['0']['bF'][i]))
    J = RD.jets(s, d='/home/claude/s1spec/step4/p2/t6/jets')
    xh = J['x'][..., 1].ravel(); kh = np.abs(J['k0']).max(axis=-1).ravel()
    yh = np.maximum(xh - 2.2490776515, 28.0)
    eb = 1.06e-8 + U(sF0) * yh ** 1.5 / 84 ** 0.5 * (1 + 1e-12)
    ke = kh * eb; m = int(np.argmax(ke))
    if ke[m] > wk[0]: wk = (float(ke[m]), s, float(xh[m]))
    if eb.max() > we[0]: we = (float(eb.max()), s)
print(f'far |e| <= {we[0]:.3e} (slab {we[1]}); sup |k~ e| <= {wk[0]:.3e} (slab {wk[1]}, x ~ {wk[2]:.4g}); '
      f'4|k~e| = {4 * wk[0]:.3e} vs theta_3 = 0.238 (T8 cconst theta_min, j = 3)')
