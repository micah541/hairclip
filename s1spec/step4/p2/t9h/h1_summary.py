"""H1 on slab 180, summary (Lemma H1 (iii) of the step-4 document): from the 16 exact points t_k = -1 + (2k+1)/16 (h1.py, Krawczyk +
Lambda; Lambda.sqrtLam_total includes the tip-offset (rotation-mode) term) and the 32 derivative balls [t_k - 1/16, t_k], [t_k, t_k + 1/16]
(h1dt.py, reference t_k; lip_term = dist sup||w_t||_Lambda):
  sqrt Lam(t) <= sqrtLam_total(t_k) + max(lip_term of the two balls of t_k) on [t_k - 1/16, t_k + 1/16]   (rows: t_k, the two terms, sum rounded up, sum <= 5e-8);
  P_hull: union of the balls' parameter boxes P_B (they contain the Krawczyk boxes);  alpha~(t) in alpha(t_k) + [-dist, dist] alpha_dot over each
  ball, so |alpha~| <= alpha_abs_max (floats with nextafter, outward).
Written to replace the earlier unsaved summary step; on the A1 inputs it reproduces out/h1slab/summary.json byte for byte.
usage: python3 -I h1_summary.py   (reads p2/out/t9_h1_p<t_k>_coh_s180.json, p2/out/h1slab/dt_<tc>_0.03125.json; writes p2/out/h1slab/summary.json)"""
import json, math, glob
O = '/home/claude/s1spec/step4/p2/out'
up = lambda v: math.nextafter(v, math.inf)
H = {}
for k in range(16):
    tk = -1 + (2 * k + 1) / 16; H[tk] = json.load(open(f'{O}/t9_h1_p{tk}_coh_s180.json'))
    assert H[tk]['Lambda']['krawczyk']['ok']
D = [json.load(open(f)) for f in sorted(glob.glob(f'{O}/h1slab/dt_*_0.03125.json'))]
assert len(D) == 32 and sorted(d['TREF'] for d in D) == sorted(t for t in H for _ in (0, 1))
rows = []
for tk in sorted(H):
    sl = H[tk]['Lambda']['sqrtLam_total']; lip = max(d['lip_term'] for d in D if d['TREF'] == tk); tot = up(sl + lip)
    rows.append([tk, sl, lip, tot, tot <= 5e-8])
for d in D:                                         # each ball started from its point's box (h1dt: P_B contains it): the points are the ones the balls read
    for key in ('xc', 'lam', 'alpha'): assert d['P_B'][key][0] <= H[d['TREF']][key][0] and H[d['TREF']][key][1] <= d['P_B'][key][1], (d['TREF'], key)
P = {key: [min(d['P_B'][key][0] for d in D), max(d['P_B'][key][1] for d in D)] for key in ('xc', 'lam', 'alpha')}
for tk in H:
    for key in ('xc', 'lam', 'alpha'): assert P[key][0] <= H[tk][key][0] and H[tk][key][1] <= P[key][1]
am = 0.0
for d in D:
    a = H[d['TREF']]['alpha']; pa = d['pdot'][2]; dist = d['dist']
    lo = min(a[0] + dist * pa[0], a[0] - dist * pa[1], a[0]); hi = max(a[1] + dist * pa[1], a[1] - dist * pa[0], a[1])
    am = max(am, abs(lo), abs(hi))
out = dict(rows=rows, sup_sqrtLam_base=max(r[3] for r in rows), P_hull=P, alpha_abs_max=up(am),
           args=dict(points=16, balls=32, slab=180, u4_cell='exact (review B5)'))
json.dump(out, open(f'{O}/h1slab/summary.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'rows'}))
assert all(r[4] for r in rows)
