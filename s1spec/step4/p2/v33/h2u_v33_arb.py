"""v3.3 review M9 (arb copy of h2u_v33.py; the only change: E_k[w'] from h2t_v33_arb.json, i.e. on the arb shell sums of m9_balls.py).
Original docstring:
v3.3 fix B0, stage 2: the tightened hand-off threshold (H2) = (1e-8, 1e-6, 1e-3) for the v3.3 deviation f - U_1 on the shells [14, 60],
true curve at t_1 (Ybar = 1e-12), uniformly over slab 180.  At t_1, U_1(r) = 1.3056 alpha(tau_1) (r/lam_1)^p = U(alpha(tau_1), y), so
    f - U_1 = w' + d,   d = alpha (psi3^'(x) - 1.3056 y^p) + alpha^2 h2^'(x)      (b0_dprof.py),
and by Minkowski in L^2((28/y) ds) on [14, 60], level by level (D^k linear),
    sqrt G_k[f - U_1] <= sqrt E_k[w']  +  sqrt E_k[d],
  E_k[w']: the certified true-curve shell energies of step 4 (p2/out_v33/h2t_v33.json 'with_error': v3.2 base shells + the shells' own
           error transfer at Ybar = 1e-12; the shell cells carry the level-0 chart allowance 2e-7 as a ball containing 0, so they bound
           the energy of f - alpha psi3^' - alpha^2 h2^' with f = tan(phi) the ray slope, the f of v3.3);
  E_k[d]:  b0_dprof.json 'at_alpha_t1' E_14_60 (|alpha| <= max over the widened balls, which contain alpha(tau_1) of the true curve).
Margin = threshold / bound (rounded down).  Requirement recorded: margin >= 2 at every level.
Writes p2/out_v33/h2u_v33.json.   usage: python3 -I -B h2u_v33.py"""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
OV = '/home/claude/s1spec/step4/p2/out_v33'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
H2T, DP = f'{OV}/h2t_v33_arb.json', f'{OV}/b0_dprof.json'
h, d = json.load(open(H2T)), json.load(open(DP))
assert h['all_ok'] and h['label'] == 'certified'
THR = [arb('1e-8'), arb('1e-6'), arb('1e-3')]
rows = []
for j in range(3):
    ew = arb(repr(h['rows'][j]['with_error'])); ed = arb(repr(d['at_alpha_t1']['E_14_60'][j]))
    tot = (ew.sqrt() + ed.sqrt()) ** 2
    rows.append(dict(j=j, threshold=float(THR[j].mid()), E_wprime_with_error=h['rows'][j]['with_error'], E_d=d['at_alpha_t1']['E_14_60'][j],
                     bound=U(tot), margin=Lw(THR[j] / tot), ok=bool(tot < THR[j]), margin_ge_2=bool(Lw(THR[j] / tot) >= 2)))
out = dict(item='tightened (H2) for f - U_1 on [14, 60], v3.3 thresholds (1e-8, 1e-6, 1e-3), true curve at t_1, slab 180', label='certified',
           rule='sqrt G_k[f - U_1] <= sqrt E_k[w\'] (h2t_v33_arb with_error) + sqrt E_k[d] (b0_dprof at_alpha_t1)', alpha_t1=d['alpha_t1'], rows=rows,
           all_ok=all(r['ok'] for r in rows), all_margin_ge_2=all(r['margin_ge_2'] for r in rows), min_margin=min(r['margin'] for r in rows),
           Ybar=h['Ybar'], args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), H2T, DP)}))
tmp = f'{OV}/h2u_v33_arb.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, f'{OV}/h2u_v33_arb.json')
for r in rows: print(f"(H2) for f - U_1, E_{r['j']}: w' {r['E_wprime_with_error']:.4e} + d {r['E_d']:.4e} -> {r['bound']:.4e} vs {r['threshold']:.0e}: x{r['margin']:.2f}")
assert out['all_ok'] and out['all_margin_ge_2']
