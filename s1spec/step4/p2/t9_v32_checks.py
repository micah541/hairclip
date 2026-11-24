"""v3.2 cone H6-data checks (DESIGN sec. 4/7): whole-cylinder data thresholds in own cylinder units, from the values t9.py already computed
(out/t9_cone_kch.json, 'cone_max' = max over the cone rows r_c in [200, 380] of slope, L sup|h''| (K2), L^1.5 ||h'''||_L2(|x|<=3L) (N3)).
Thresholds: slope <= 0.08, sup|h0''| <= 0.03, ||h0'''||_L2 <= 0.004.  Certified arb comparison (margin = threshold / bound, rounded down).
Also re-derives each max from the per-row values and asserts it equals cone_max.  usage: python3 -I t9_v32_checks.py"""
import json, math
from flint import arb, ctx
ctx.prec = 128
O = '/home/claude/s1spec/step4/p2/out'
d = json.load(open(f'{O}/t9_cone_kch.json')); cm = d['cone_max']
for k in ('slope', 'K2', 'N3'): assert max(r[k] for r in d['cone']) == cm[k], k
TH = dict(slope='0.08', K2='0.03', N3='0.004')
Lw = lambda v: math.nextafter(float(arb(v).lower()), -math.inf)
res = {}
for k, name in (('slope', 'cone_slope'), ('K2', 'cone_h2'), ('N3', 'cone_N3')):
    v, t = arb(repr(cm[k])), arb(TH[k])
    res[name] = dict(value=cm[k], threshold=float(t.mid()), margin=Lw(t / v), ok=bool(v < t))
out = dict(checks=res, all_ok=all(r['ok'] for r in res.values()),
           args=dict(source=f'{O}/t9_cone_kch.json', thresholds=TH, units='own cylinder units, whole cylinder (cone r_c 200-380)'))
json.dump(out, open(f'{O}/t9_v32_checks.json', 'w'), indent=1); print(json.dumps(out))
