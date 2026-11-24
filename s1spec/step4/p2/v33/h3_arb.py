"""v3.3 review M9, item 1: the v3.2 hand-off items (H3) (windows [y, 1.5 y], y in [50, 253], levels 0-3) and (H2) (shells [14, 60],
levels 0-2, v3.2 thresholds) re-assembled on the arb base sums of m9_balls.py (p2/out_v33/m9/win_arb_summary.json 'combined'), with
camp.py's formulas (p2/t8se/camp.py handoff(), items 3 and 2; closed v3.2, not imported) and its pinned inputs:
  sqrt E_j[true] <= sqrt E_j[base] + amp_j Ybar + bres_j,
  amp_j = KCH sum_m S(j,m) A_m + c0_j A_0 + c0tip_j,   bres_j = KCH sum_m S(j,m) B_m + c0_j B_0,
  A_m, B_m: amp_A1_w180.json 'handoff' (window slab 180); KCH, c0, c0tip: t9_cone_kch.json (windows: top level; shells: block 'shell'),
  thresholds (camp THR): windows E = (2e-7, 1e-5, 1e-5, 1e-5), shells (2e-7, 1e-4, 0.2);
  Ybar limit (camp 'lim'): (sqrt(thr) - sqrt(base) - bres) / amp;  the hand-off value Ybar = 1e-12 (amp 'at_Ybar').
Reported for each item: arb base vs the float base camp read (win_summary.json, nextafter up), E with error at Ybar = 1e-12, margin
thr / E (rounded down), Ybar limit and its ratio to Ybar; and, as a fidelity check, the Ybar limits recomputed here from the float base
against the pinned amp_A1_w180.json 'lim' (camp ran at 64 bits; agreement to ~1e-15 relative expected).
Writes p2/out_v33/h3_arb.json.   usage: python3 -I -B h3_arb.py"""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
UP = lambda v: math.nextafter(float(v.upper()), math.inf); LO = lambda v: math.nextafter(float(v.lower()), -math.inf)
WA, WS, T9, AMP = f'{OV}/m9/win_arb_summary.json', f'{P2}/out/h1slab/win_summary.json', f'{P2}/out/t9_cone_kch.json', f'{P2}/t8se/constants/amp_A1_w180.json'
wa, ws, t9, h = json.load(open(WA)), json.load(open(WS)), json.load(open(T9)), json.load(open(AMP))['handoff']
assert h['window_slabs'] == [180, 180]
YB = arb(repr(h['at_Ybar']['Ybar']))
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
A = {j: arb(repr(h['A'][str(j)])) for j in range(4)}; B = {j: arb(repr(h['B'][str(j)])) for j in range(4)}
THR = dict(win=[arb('2e-7'), arb('1e-5'), arb('1e-5'), arb('1e-5')], shell=[arb('2e-7'), arb('1e-4'), arb('0.2')])
def consts(kind, j):
    if kind == 'win': return arb(repr(t9['KCH'])), arb(repr(t9['c0'].get(str(j), 0.0))), arb(repr(t9['c0tip'].get(str(j), 0.0)))
    s = t9['shell']; return arb(repr(s['KCH'])), arb(repr(s['c0'][str(j)])), arb(repr(s['c0tip'][str(j)]))
rows = {}; fid = {}
for kind, key, n, item in (('win', 'windows_max', 4, 'item3'), ('shell', 'shells', 3, 'item2')):
    rows[kind] = []
    for j in range(n):
        ks, c0, ct = consts(kind, j)
        amp_ = ks * sum((s * A[m] for m, s in STIR[j].items()), arb(0)) + c0 * A[0] + ct
        bres = ks * sum((s * B[m] for m, s in STIR[j].items()), arb(0)) + c0 * B[0]
        base = arb(repr(wa['combined'][key][j])); fbase = arb(repr(math.nextafter(float(ws[key][j]), math.inf)))
        tot = (base.sqrt() + amp_ * YB + bres) ** 2; ftot = (fbase.sqrt() + amp_ * YB + bres) ** 2
        lim = (THR[kind][j].sqrt() - base.sqrt() - bres) / amp_; flim = (THR[kind][j].sqrt() - fbase.sqrt() - bres) / amp_
        m_a, m_f = LO(THR[kind][j] / tot), LO(THR[kind][j] / ftot)
        rows[kind].append(dict(j=j, threshold=float(THR[kind][j].mid()), base_arb=UP(base), base_float=UP(fbase), base_rel=float((base / fbase - 1).mid()),
                               with_error_arb=UP(tot), with_error_float=UP(ftot), margin_arb=m_a, margin_float=m_f, margin_rel_change=m_a / m_f - 1,
                               Ybar_limit_arb=LO(lim), Ybar_limit_float=LO(flim), Ybar_margin_arb=LO(lim / YB), ok=bool(tot < THR[kind][j] and lim > YB)))
        pin = h['lim'][f'{item}_E{j}']; fid[f'{item}_E{j}'] = dict(recomputed_float=LO(flim), pinned=pin, rel=LO(flim) / pin - 1)
allr = rows['win'] + rows['shell']
out = dict(item='M9: v3.2 (H3) windows and (H2) shells (v3.2 thresholds) on the arb base sums', label='certified', rows=rows,
           all_ok=all(r['ok'] for r in allr), binding_Ybar_limit_arb=min(r['Ybar_limit_arb'] for r in allr),
           binding_Ybar_limit_float=min(r['Ybar_limit_float'] for r in allr), max_abs_margin_rel_change=max(abs(r['margin_rel_change']) for r in allr),
           fidelity_vs_pinned_camp_lim=fid, max_fidelity_rel=max(abs(v['rel']) for v in fid.values()), Ybar=UP(YB),
           args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), WA, WS, T9, AMP)}))
fo = f'{OV}/h3_arb.json'; json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
for kind in rows:
    for r in rows[kind]:
        print(f"{kind} E{r['j']}: base arb {r['base_arb']:.6e} (float {r['base_float']:.6e}, rel {r['base_rel']:+.2e}); with error {r['with_error_arb']:.4e} "
              f"vs {r['threshold']:.0e}: x{r['margin_arb']:.3f} (float x{r['margin_float']:.3f}); Ybar limit {r['Ybar_limit_arb']:.5e} (x{r['Ybar_margin_arb']:.2f})")
print('binding Ybar limit: arb %.6e, float %.6e; fidelity vs pinned camp lim: max rel %.2e' % (out['binding_Ybar_limit_arb'], out['binding_Ybar_limit_float'], out['max_fidelity_rel']))
assert out['all_ok']
