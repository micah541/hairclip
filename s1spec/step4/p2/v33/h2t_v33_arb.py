"""v3.3 review M9 (arb copy of h2t_v33.py; the only change: the base shells are the arb sums of m9_balls.py,
p2/out_v33/m9/win_arb_summary.json 'combined' (max over the 32 t_c balls and both parameter-ball sets), not win.py's float sums).
Original docstring:
v3.3 task 6 (CONTRACTS_r2, task 6: 'the tightened (H2), as planned'): the shells [14, 60] of step 4 against the tightened hand-off
threshold (H2) = (1e-8, 1e-6, 1e-3) of v3.3, true curve at Ybar = 1e-12, uniformly over slab 180.  Pure re-assembly of pinned v3.2
outputs (nothing recomputed): base shells from p2/out/h1slab/win_summary.json, the shells' own error constants (t9.py block 'shell',
Lambda_e route: KCH = 0, c0 = 0, c0tip per unit Ybar) from p2/out/t9_cone_kch.json, and camp.py's item-2 formula
  sqrt E_j[true] <= sqrt E_j[base] + KCH sum S(j,m) sqrt F_m + c0_j sqrt F_0 + c0tip_j Ybar,  sqrt F_m <= Ybar A_m + B_m (amp_A1_w180.json).
Writes p2/out_v33/h2t_v33.json (margin = threshold / bound, rounded down; Ybar limit as camp)."""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
O = '/home/claude/s1spec/step4/p2/out'; WS, T9, AMP = '/home/claude/s1spec/step4/p2/out_v33/m9/win_arb_summary.json', f'{O}/t9_cone_kch.json', '/home/claude/s1spec/step4/p2/t8se/constants/amp_A1_w180.json'
w, t, h = json.load(open(WS))['combined'], json.load(open(T9))['shell'], json.load(open(AMP))['handoff']
YB = arb(repr(h['at_Ybar']['Ybar'])); THR = [arb('1e-8'), arb('1e-6'), arb('1e-3')]
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}}
A = {j: arb(repr(h['A'][str(j)])) for j in range(3)}; B = {j: arb(repr(h['B'][str(j)])) for j in range(3)}
rows = []
for j in range(3):
    ks, c0s, cts = arb(repr(t['KCH'])), arb(repr(t['c0'][str(j)])), arb(repr(t['c0tip'][str(j)]))
    amp_ = ks * sum((s * A[m] for m, s in STIR[j].items()), arb(0)) + c0s * A[0] + cts
    bres = ks * sum((s * B[m] for m, s in STIR[j].items()), arb(0)) + c0s * B[0]
    base = arb(repr(math.nextafter(float(w['shells'][j]), math.inf)))
    tot = (base.sqrt() + amp_ * YB + bres) ** 2
    rows.append(dict(j=j, threshold=float(THR[j].mid()), base=U(base), with_error=U(tot), margin=Lw(THR[j] / tot), Ybar_limit=Lw((THR[j].sqrt() - base.sqrt() - bres) / amp_), ok=bool(tot < THR[j])))
out = dict(item='tightened (H2) shells [14, 60], v3.3 thresholds (1e-8, 1e-6, 1e-3)', label='certified', base_source='m9 arb cells (win_arb_summary combined)', rows=rows, all_ok=all(r['ok'] for r in rows),
           min_margin=min(r['margin'] for r in rows), Ybar=U(YB), args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), WS, T9, AMP)}))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v33/h2t_v33_arb.json', 'w'), indent=1)
for r in rows: print(f"(H2) tightened E_{r['j']}: base {r['base']:.3e}, with error {r['with_error']:.3e} vs {r['threshold']:.0e}: x{r['margin']:.1f}")
assert out['all_ok']
