"""H6 margins for the proposed step-3 amendment (analysis only): cone K2 = K3 = 0.25 (X = 3, M = 0.12), far (M, K) = (0.1, 0.2) (X = 2, L = 113).
Covering: Lemma O's conclusion holds on |x| <= 1 of each member (lemmaO_v2 l.19). For p in cylinder c (r_p >= r_c/4) the member with the
largest unit whose |x| <= 1 contains p has r' = 4 r_p/3, so L_c/L' = 3 r_c/(4 r_p) <= 3. Bounds recovered at p: |h''| <= data + C2 sqrt t (C2 at M).
Inner strip (r_p < 150, no cone member covers it): step 3's time-uniform core bound M1 = sup y|f'| <= 0.1844 gives L_c |k| <= (y + 2.26)/y M1 (L_c <= r_p).
Outer strip (r_p > 475): the far member centred at p, in cone units L_c/113 <= 95/113."""
import json, math
from flint import arb
U = lambda v: math.nextafter(float(arb(v).upper()), math.inf); Lw = lambda v: math.nextafter(float(arb(v).lower()), -math.inf)
C = json.load(open('/home/claude/s1spec/step4/p2/h6/h6_cases.json')); O = '/home/claude/s1spec/step4/p2/out'
cone = C['cone X=3 M=0.12 K=0.25']['deviation_upper_at_needed_time']; farc = C['far X=2 M=0.1 K=0.2']['deviation_upper_at_needed_time']
t9 = json.load(open(f'{O}/t9_cone_kch.json'))['cone_max']; tf = json.load(open(f'{O}/t9_far.json')); h4 = json.load(open(f'{O}/h4loc_A1_w180.json'))
Kc, Kf, Mf_new = arb('0.25'), arb('0.2'), arb('0.1')
c2c, c3c = arb(repr(cone['C2'])), arb(repr(cone['sqrtC3t'])); c2f, c3f = arb(repr(farc['C2'])), arb(repr(farc['sqrtC3t']))
Mc = arb(repr(t9['K2'])); Mf = arb(113) * arb(repr(tf['h2']))
rec_c = Mc + c2c; rec_f = Mf + c2f
m = {}
m['cone_K2_by_cover (ratio 3)'] = Lw(Kc / (3 * rec_c))
m['cone_K3_L2 (pieces sum ratio^3 <= 37)'] = Lw(Kc / (arb(37).sqrt() * c3c + arb(repr(t9['N3']))))
m['cone_inner_strip (M1core 0.1844, y >= 50)'] = Lw(Kc / (arb('0.1844') * (arb(50) + arb('2.26')) / 50 * arb('1.006')))
m['cone_outer_strip_from_far'] = Lw(Kc / (arb(95) / 113 * rec_f))
m['far_self_cover (ratio 1)'] = Lw(Kf / rec_f)
m['far_inner_strip_from_cone (113/50)'] = Lw(Kf / (arb(113) / 50 * rec_c))
m['far_data_slope_vs_M'] = Lw(Mf_new / arb(repr(tf['slope']))); m['far_data_j2_vs_M'] = Lw(Mf_new / Mf)
m['far_data_j4_vs_M'] = Lw(Mf_new / (arb(113) ** 3 * arb(repr(tf['h4']))))
m['far_H5_L2 (K3 - sqrt(2 C3 t) vs data)'] = Lw((Kf - arb(2).sqrt() * c3f) / (arb(113) ** arb('1.5') * arb(repr(tf['L2']))))
m['cone_H5_L2 (K3 - sqrt(3 C3 t) vs data)'] = Lw((Kc - arb(3).sqrt() * c3c) / arb(repr(t9['N3'])))
m['far_H4 (0.0052 M2_old=1 -> M2_new=0.01)'] = Lw(arb('0.01') / arb('0.005169'))
m['tstar_cone (need 1.6e-7)'] = C['cone X=3 M=0.12 K=0.25']['margin_factor_lower']; m['tstar_far (need 6.13e-9)'] = C['far X=2 M=0.1 K=0.2']['margin_factor_lower']
out = dict(args=dict(cases='h6_cases.json', t9=f'{O}/t9_cone_kch.json', t9_far=f'{O}/t9_far.json'), recovered=dict(cone=U(rec_c), far=U(rec_f)), margins=m)
json.dump(out, open('/home/claude/s1spec/step4/p2/h6/h6_margins.json', 'w'), indent=1)
for k, v in m.items(): print(f'{k}: x{v:.3f}')
