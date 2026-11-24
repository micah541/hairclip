"""v3.3 task 6: Lemma O H7 (time ratio) for the fixed-L cone family, certified in arb.
  (T - t_1)/L^4 <= lambda_1^4/(1 - 4 abar) / (50 lambda_1)^4 = 50^-4/(1 - 4 abar)      (step3_proof_v32.tex l.346-348; every member has L = 50)
abar = sup|(a, b)| of the v3.3 verifier rerun (step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json, alpha1 = 1.5e-4: sup_ab;
unchanged from v3.2, 2.19103e-7).  Asserts ratio <= 1.6001e-7 (the needed ratio of lemmaO_v32_consts.json, case 'S1 cone, X=3') and
reports the margin t*/ratio against that case's certified t*_lower (3.68334e-7), and t*/1.6001e-7.  Writes p2/out_v33/ratio_v33.json.
usage: python3 -I ratio_v33.py"""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
MS = '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'
LO = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json'
abar = arb(repr(json.load(open(MS))['1.5e-4']['sup_ab']))
case = next(c for c in json.load(open(LO))['cases'] if c['name'] == 'S1 cone, X=3')
need, tstar = arb(repr(case['needed'])), arb(repr(case['t_star_lower']))
ratio = arb(50) ** -4 / (1 - 4 * abar)
out = dict(item='Lemma O H7, fixed-L cone family (L = 50)', label='certified', abar=float(abar.mid()), ratio_upper=U(ratio), needed=case['needed'],
           ratio_le_needed=bool(ratio < need), t_star_lower=case['t_star_lower'], margin_tstar_over_ratio=Lw(tstar / ratio), margin_tstar_over_needed=Lw(tstar / need),
           C0_sqrt_t_at_needed=case['deviation_upper_at_needed_time']['C0'],
           args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), MS, LO)}))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v33/ratio_v33.json', 'w'), indent=1); print(json.dumps(out))
assert out['ratio_le_needed']
