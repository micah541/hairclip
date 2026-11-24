"""Lemma O H6 (a-priori bounds close) for the step-3 families, from the hand-off data (t9.py, t9_far.py) and lemmaO_v3_consts.json.
The bound recovered at a point p from the lemma centred at p is M' + C2 sqrt(t) in p-units (M' = the certified data bound there);
in the units of a cylinder c containing p it is (L_c/L_p) (M' + C2 sqrt(t)), to be compared with K2 (Lemma O hypothesis (A) on [-X, X]).
cone: L = r/4, X = 3, half-width 0.75 r, so L_c/L_p <= 4 (inner end); needed time 1.6e-7 (C2 sqrt t scales as sqrt t from the table's 1.2e-7).
far: L = 113, X = 2 (ratio 1). Cone cylinders beyond the last cone centre get their bound from the far family: ratio L_cone/113."""
import json, math
from flint import arb
U = lambda v: math.nextafter(float(v.upper()), math.inf)
O = '/home/claude/s1spec/step4/p2/out'; LC = json.load(open('/home/claude/s1spec/writeup/lemmaO/lemmaO_v3_consts.json'))
cone = [c for c in LC['cases'] if c['name'] == 'S1 cone, X=3'][0]; far = LC['far_margin2']['case']
C2c = arb(repr(cone['deviation_upper_at_needed_time']['C2'])) * (arb('1.6e-7') / arb(repr(cone['needed']))).sqrt()
C2f = arb(repr(far['deviation_upper_at_needed_time']['C2']))
K2c, K2f = arb('0.18'), arb('1.5')
Mc = arb(repr(json.load(open(f'{O}/t9_cone_kch.json'))['cone_max']['K2']))         # cone data |h0''| (own units)
Mf = arb(113) * arb(repr(json.load(open(f'{O}/t9_far.json'))['h2']))               # far data |h0''| (far units)
rc = dict(inner_ratio=4, recovered_p_units=U(Mc + C2c), recovered_c_units=U(4 * (Mc + C2c)), K2=0.18, margin=float((K2c / (4 * (Mc + C2c))).lower()),
          max_ratio_ok=float((K2c / (Mc + C2c)).lower()))
rf = dict(ratio=1, recovered=U(Mf + C2f), K=1.5, margin=float((K2f / (Mf + C2f)).lower()))
outer = dict(L_cone_last=95, far_recovered_in_cone_units=U(arb(95) / 113 * (Mf + C2f)), K2=0.18, margin=float((K2c / (arb(95) / 113 * (Mf + C2f))).lower()),
             far_C2sqrt_needed=float((K2c * 113 / 95 - Mf).lower()))
fin = dict(L_p_min=37.5, cone_recovered_in_far_units=U(arb(113) / arb('37.5') * (Mc + C2c)), K=1.5, margin=float((K2f / (arb(113) / arb('37.5') * (Mc + C2c))).lower()))
out = dict(args=dict(lemmaO='writeup/lemmaO/lemmaO_v3_consts.json', t9=f'{O}/t9_cone_kch.json', t9_far=f'{O}/t9_far.json', C2c_at_1p6e7=U(C2c), C2f=U(C2f)),
           cone_by_overlap=rc, far_by_overlap=rf, cone_outer_from_far=outer, far_inner_from_cone=fin)
json.dump(out, open(f'{O}/t9_h6.json', 'w'), indent=1); print(json.dumps(out))
