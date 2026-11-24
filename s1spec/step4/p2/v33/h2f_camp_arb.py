"""v3.3 review M9 (arb copy of h2f_camp.py; the only change: G_j[base] is the arb cell sum of m9_balls.py (win_arb.py: arb
half-width, exact cell breakpoints, arb sums), p2/out_v33/m9/win_arb_summary.json 'combined' = max over the 32 t_c balls and both
parameter-ball sets (the widened (H2f) balls and v3.2's balls); parts from the same file).
Original docstring:
v3.3 task 6, (H2f) stage 4: the hand-off item (H2f) for the TRUE curve at t_1, uniformly over slab 180:
    G_j = int_{ln 60}^{ln 700} (28/y) (d_s^j (f - U_1))^2 ds <= H2f_j = (1e-6, 5e-6, 1e-5),  j = 0, 1, 2  (v3.3, beta_A = 2).
Error transfer exactly as camp.py's item 3 does for the (H3) windows (LEMMA_F C3), with the constants of t9_h2f.py on [60, 700]:
    sqrt G_j[true] <= sqrt G_j[base] + KCH_f sum_m S(j,m) sqrt F_m + c0_f[j] sqrt F_0 + c0tip_f[j] Ybar,   sqrt F_m <= Ybar A_m + B_m,
A_m, B_m: the window-180 camp run (t8se/constants/amp_A1_w180.json 'handoff', pinned v3.2), Ybar = 1e-12 (ycomp: Y(tau_a) <= Ybar).
G_j[base]: h2f_balls.py (max over the 32 t_c balls, on parameter balls widened by the true-curve frame shift and C_alpha YBOX, so the
true curve's frame and alpha are covered).  All boxes used (n3 boxes, c_e, frame shift) hold for Ybar <= YBOX = 1.7e-12.
Reported: base, with error at Ybar = 1e-12, margin = threshold / bound (rounded down), and the Ybar limit
(sqrt(thr) - sqrt(base) - bres)/amp as in camp.py (the result is certified only for Ybar <= YBOX, recorded).  Also the (H3)-style
informational parts on [60, 379.5] and [379.5, 700] (base).  Stop trigger (CONTRACTS_r2): margin < 1.2.
usage: python3 -I h2f_camp.py      (writes p2/out_v33/h2f_result.json)"""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
BASE, T9F, AMP, YC = f'{OV}/m9/win_arb_summary.json', f'{OV}/t9_h2f.json', f'{P2}/t8se/constants/amp_A1_w180.json', f'{P2}/out/ycomp_final.json'
b, t, a = json.load(open(BASE))['combined'], json.load(open(T9F)), json.load(open(AMP))
h = a['handoff']; assert h['window_slabs'] == [180, 180]
YB = arb(repr(h['at_Ybar']['Ybar'])); YBOX = arb(t['boxes']['Y'])
assert YB <= YBOX
THR = [arb('1e-6'), arb('5e-6'), arb('1e-5')]
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}}
A = {j: arb(repr(h['A'][str(j)])) for j in range(4)}; B = {j: arb(repr(h['B'][str(j)])) for j in range(4)}
KCH = arb(repr(t['KCH'])); C0 = {int(k): arb(repr(v)) for k, v in t['c0'].items()}; CT = {int(k): arb(repr(v)) for k, v in t['c0tip'].items()}
rows = []
for j in range(3):
    gb = arb(repr(b['G_60_700'][j]))
    amp_ = KCH * sum((s * A[m] for m, s in STIR[j].items()), arb(0)) + C0[j] * A[0] + CT[j]
    bres = KCH * sum((s * B[m] for m, s in STIR[j].items()), arb(0)) + C0[j] * B[0]
    err = amp_ * YB + bres
    gt = (gb.sqrt() + err) ** 2
    lim = (THR[j].sqrt() - gb.sqrt() - bres) / amp_
    rows.append(dict(j=j, threshold=float(THR[j].mid()), base=U(gb), sqrt_err_at_Ybar=U(err), with_error=U(gt), margin=Lw(THR[j] / gt),
                     margin_base_only=Lw(THR[j] / gb), Ybar_limit=Lw(lim), Ybar_limit_capped_at_YBOX=min(Lw(lim), float(YBOX.mid())),
                     ok=bool(gt < THR[j]), stop_trigger_margin_lt_1p2=bool(Lw(THR[j] / gt) < 1.2)))
yc = json.load(open(YC))
out = dict(item='(H2f) v3.3: G_j on [60, 700] at t_1, uniformly over slab 180, true curve', label='certified', base_source='m9 arb cells (win_arb_summary combined)',
           thresholds=[1e-6, 5e-6, 1e-5], rows=rows, all_ok=all(r['ok'] for r in rows), min_margin=min(r['margin'] for r in rows),
           Ybar=U(YB), YBOX=str(t['boxes']['Y']), constants=dict(KCH=t['KCH'], c0=t['c0'], c0tip=t['c0tip']), A={j: U(A[j]) for j in range(3)},
           B={j: U(B[j]) for j in range(3)}, base_parts=dict(G_60_379p5=b['G_60_379p5'][:3], G_379p5_700=b['G_379p5_700'][:3]),
           probe_for_comparison=[3.126224097280483e-07, 1.0030690031212077e-06, 2.4301664752850864e-06],
           ycomp=dict(file=YC, note='ycomp asserts Ybar = 1e-12 >= Y(tau_a)'),
           args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), BASE, T9F, AMP, YC)}))
json.dump(out, open(f'{OV}/h2f_result_arb.json.tmp', 'w'), indent=1); os.replace(f'{OV}/h2f_result_arb.json.tmp', f'{OV}/h2f_result_arb.json')
for r in rows:
    print(f"(H2f) G_{r['j']}: base {r['base']:.4e}, with error {r['with_error']:.4e} vs {r['threshold']:.0e}: margin x{r['margin']:.3f} "
          f"(Ybar limit {r['Ybar_limit']:.3e})")
assert out['all_ok'], 'H2f FAILS'
