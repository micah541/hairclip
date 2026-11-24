"""hand-off table (v3 thresholds) from handoff_<TAG>.json, plus the M1 comparison"""
import sys, json, glob
H = json.load(open(f'handoff_{sys.argv[1]}.json'))
M1 = {}
for f in glob.glob('/home/claude/s1spec/step4/m1/res_*.json'):
    d = json.load(open(f)); M1.setdefault(round(d['kmax']), []).append((f.split('/')[-1], d))
def ok(v, t): return f"{v:.2e} {'OK' if v <= t else 'FAIL'} x{t/v:.1f}" if v <= t else f"{v:.2e} FAIL /{v/t:.1f}"
for h in H:
    k = h['k']; print(f"\n=== k = {k:.0f} (tau = {h['tau']:.3f}), alpha*k = {h['alpha_k']:.4f}, chart ext {h['chart_ext']:.0f}, data to y = {h['ymax']:.0f}")
    print('|alpha|  (1.5e-4):', ok(abs(h['alpha']), 1.5e-4), '   sqrt(Lambda) (5e-8):', ok(h['sqrtLam'], 5e-8), f"  A = {h['A']:.3e}  A/|a|^3 = {h['A']/abs(h['alpha'])**3:.0f}")
    for lbl in ('v3', 'M1'):
        s = h[f'shell14_60_{lbl}']
        print(f'shells [14,60] U_{lbl}: E0 {ok(s[0], 2e-7)} | E1 {ok(s[1], 1e-4)} | E2 {ok(s[2], 0.2)}')
        sl = h[f'slide_{lbl}']; at = h[f'slide_at_{lbl}']
        print(f'windows [y,1.5y], y<=253, U_{lbl}: E0 {ok(sl[0], 2e-7)} | E1 {ok(sl[1], 1e-5)} | E2 {ok(sl[2], 1e-5)} | E3 {ok(sl[3], 1e-5)}   (max at y = {at[0]:.0f},{at[3]:.0f}; E4 = {sl[4]:.2e})')
    c = h['cone']; print(f"cone r={c['r']:.0f} L={c['L']:.1f} X=3: K4 {ok(c['K4'], 0.06)} | K2 {c['K2']:.3e} K3 {c['K3']:.3e} slope {c['slope']:.3f} N3 {c['N3']:.2e}")
    print('  -v/(y/k)^3 (U_M1):', {y: round(v, 2) for y, v in h['v_over_cubic'].items()})
    for nm, d in M1.get(round(k), []):
        print(f"  M1 {nm}: alpha {d['alpha']:.5e} sqrtLam {d['inner']['sqrtLam']:.3e} A {d['inner']['A']:.3e} | [14,60] {['%.3e' % v for v in d['shell']['14-60'][:3]]} | slide {['%.3e' % d['slide'][str(j)][0] for j in range(4)]} | cone K4 {d['cone_sup']['200']['K4']:.3e}")
        print(f"  base       : alpha {h['alpha']:.5e} sqrtLam {h['sqrtLam']:.3e} A {h['A']:.3e} | [14,60] {['%.3e' % v for v in h['shell14_60_M1'][:3]]} | slide {['%.3e' % v for v in h['slide_M1'][:4]]} | cone K4 {c['K4']:.3e}")
