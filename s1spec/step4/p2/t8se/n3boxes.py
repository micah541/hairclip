"""Review B2: check the bootstrap boxes that n3far.py assumed against the certified far-core box with N3' included.

The far-core bounds are those of p2/boxes.py (LEMMA_F C4 on intervals of length >= y/3 inside [y/1.5, 1.5y]):
  |e'|^2      <= (6.75/28) y F0 + (3/28) y sqrt(F0 F1)
  x^2|e''|^2  <= 1.0804^2 y [(3/28) F1 + (2*1.5/0.8828/28) sqrt(F1 F2)]
  x^4|e'''|^2 <= 1.0804^4 y [(3 c2/28) F2 + (2 sqrt(c2 c3)/28) sqrt(F2 F3)],  c2 = (1.5/0.8828)^2, c3 = c2^2,
with sqrt F_j <= aF_j Ybar + bF_j per slab from the camp run WITH --n3, and y up to y_l (sup over the slab).  Asserted at the Ybar of
the n3 file (the boxes are monotone in Ybar, so every smaller Ybar is covered):
  sup|e'| <= b0,  sup x|e''| <= b1,  sup x^2|e'''| <= b2,  and per slab sqrt F0 (with N3') <= fac * sqrt F0 (the one in the |e| box).
The tip box |e|, |e'| <= eta1 on x <= 34.25 is checked against T7 Lemma 4b (t7/constants.json meta.tip_box, same Ybar).
usage: python3 -I n3boxes.py N3_JSON AMP_WITH_N3_JSON CONSTS_JSON [YBAR_REPORT]"""
import sys, json, math
from flint import arb, ctx
ctx.prec = 128
n3, amp, cons = (json.load(open(p)) for p in sys.argv[1:4])
YR = arb(sys.argv[4] if len(sys.argv) > 4 else '1e-12')
args = n3['args']; Yn = arb(args['Y']); fac = arb(repr(args['fac']))
b0, b1, b2 = arb(repr(args['b0'])), arb(repr(args['b1'])), arb(repr(args['b2']))
used = json.load(open(args['amp']))
X0 = arb('2.2490776515'); r = arb('1.5') / arb('0.8828'); c2 = r * r; c3 = c2 * c2; xf = 1 + X0 / 28
U = lambda v: math.nextafter(float(v.upper()), math.inf)
ell = {p['slab']: arb(repr(p['scalars']['ell'][1])) for p in cons['per_slab']}
ps = amp['per_slab']; slabs = amp['slabs']; psu = used['per_slab']; slu = used['slabs']


def boxes(Y):
    worst = [0.0, 0.0, 0.0]; ratio0 = 0.0
    for i, s in enumerate(slabs):
        sF = [arb(repr(ps[str(j)]['aF'][i])) * Y + arb(repr(ps[str(j)]['bF'][i])) for j in range(4)]
        F = [v * v for v in sF]
        y = ell[s] - X0
        b = [(arb('6.75') / 28 * y * F[0] + arb(3) / 28 * y * (F[0] * F[1]).sqrt()).sqrt(),
             (xf ** 2 * y * (arb(3) / 28 * F[1] + 2 * r / 28 * (F[1] * F[2]).sqrt())).sqrt(),
             (xf ** 4 * y * (3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt())).sqrt()]
        for k in range(3): worst[k] = max(worst[k], U(b[k]))
        k_ = slu.index(s)
        s0 = arb(repr(psu['0']['aF'][k_])) * Y + arb(repr(psu['0']['bF'][k_]))
        ratio0 = max(ratio0, U(sF[0] / s0))
    return worst, ratio0


wn, rn = boxes(Yn)
ok = [wn[0] <= float(b0.lower()), wn[1] <= float(b1.lower()), wn[2] <= float(b2.lower()), rn <= float(fac.lower())]
print(f"at the n3 Ybar = {float(Yn.mid()):.3g}: |e'| {wn[0]:.3e} <= {float(b0.mid()):.3g}, x|e''| {wn[1]:.3e} <= {float(b1.mid()):.3g}, "
      f"x^2|e'''| {wn[2]:.3e} <= {float(b2.mid()):.3g}, sqrtF0 ratio {rn:.5f} <= {float(fac.mid()):.3g}: {'all hold' if all(ok) else 'FAIL'}")
assert all(ok), ok
wr, rr = boxes(YR)
m = [float(b0.lower()) / wr[0], float(b1.lower()) / wr[1], float(b2.lower()) / wr[2], float(fac.lower()) / rr]
import os
tb = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 't7', 'constants.json')))['meta'].get('tip_box', {})
tip_ok = bool(tb) and abs(tb.get('Ybox', -1) - float(Yn.mid())) <= 1e-6 * float(Yn.mid()) and tb['eta1_at_Ybox'] <= args['eta1']
print(f"margins at Ybar = {float(YR.mid()):.3g}: |e'| x{m[0]:.2f}, x|e''| x{m[1]:.2f}, x^2|e'''| x{m[2]:.2f}, |e| (sqrtF0) x{m[3]:.2f}; "
      f"thinnest x{min(m):.2f}.  Tip box: |e|, |e'| <= {args['eta1']:.1e} assumed; T7 Lemma 4b gives "
      f"{tb.get('eta1_at_Ybox', float('nan')):.3e} at Ybar = {tb.get('Ybox', float('nan')):.3g}: {'holds' if tip_ok else 'NOT CHECKED'}")
assert tip_ok, 'tip box: T7 meta.tip_box missing, at another Ybox, or larger than eta1'
