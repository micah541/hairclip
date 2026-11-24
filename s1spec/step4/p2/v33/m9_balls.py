"""v3.3 review M9, item 1, stage 1: win_arb.py on the 32 t_c balls of slab 180 for both parameter-ball sets (v32: v3.2's IFT balls
p2/out/h1slab/ift_M.json, as win_cmds.txt; wide: the (H2f) balls p2/out_v33/h2f/ift_M.json, widened by the true-curve frame shift and
C_alpha YBOX), one subprocess at a time (python3 -I -B, nice inherited).  Resumable: a run is skipped if its output exists and records
the current sha256 of win_arb.py and of its ball file.  Then the slab-wide values (max over balls; each ball's sum already rounded up),
per set and combined (max over both sets: the value passed on), and the comparison with the float values they replace:
  v32 set:  windows_max, shells  vs  p2/out/h1slab/win_M.json (win.py, pinned v3.2) and win_summary.json;
  wide set: G_60_700, parts      vs  p2/out_v33/h2f/win_M.json 'G_*' (h2f_win.py, float half-width) and h2f_base.json.
Writes p2/out_v33/m9/win_arb_summary.json.   usage: python3 -I -B m9_balls.py"""
import sys, os, json, math, hashlib, subprocess
P2 = '/home/claude/s1spec/step4/p2'; OM = f'{P2}/out_v33/m9'; HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
WA = f'{HERE}/win_arb.py'; os.makedirs(OM, exist_ok=True)
IFT = {'v32': lambda m: f'{P2}/out/h1slab/ift_{m}.json', 'wide': lambda m: f'{P2}/out_v33/h2f/ift_{m}.json'}
reused = []
for S in ('v32', 'wide'):
    for m in range(32):
        o = f'{OM}/win_arb_{S}_{m}.json'
        if os.path.exists(o):
            a = json.load(open(o))['args']['inputs']
            if a.get(WA) == sha(WA) and a.get(IFT[S](m)) == sha(IFT[S](m)): reused.append(f'{S}_{m}'); continue
        r = subprocess.run([sys.executable, '-I', '-B', WA, S, str(m)], cwd=HERE, stdout=open(f'{OM}/win_arb_{S}_{m}.log', 'w'), stderr=subprocess.STDOUT)
        assert r.returncode == 0, f'{S} ball {m} failed'
        print(open(f'{OM}/win_arb_{S}_{m}.log').read().strip().splitlines()[-1], flush=True)
B = {S: [json.load(open(f'{OM}/win_arb_{S}_{m}.json')) for m in range(32)] for S in ('v32', 'wide')}
for S in B:
    assert [b['ball'] for b in B[S]] == list(range(32)) and all(b['set'] == S for b in B[S])
    cov = sorted(b['tc_ball'] for b in B[S]); assert cov[0][0] <= -1 and cov[-1][1] >= 1 and all(cov[i][1] >= cov[i + 1][0] for i in range(31))
KEYS = {'windows_max': 4, 'shells': 3, 'G_60_700': 4, 'G_60_379p5': 4, 'G_379p5_700': 4}
mx = lambda L, k: [max(b[k][j] for b in L) for j in range(KEYS[k])]
per = {S: {k: mx(B[S], k) for k in KEYS} for S in B}
for S in B: per[S]['worst_ball'] = {k: [max(B[S], key=lambda b: b[k][j])['ball'] for j in range(KEYS[k])] for k in KEYS}
for S in B: per[S]['max_rel_rad_hh'] = max(b['max_rel_rad_hh'] for b in B[S]); per[S]['min_cos_phi'] = min(b['coverage']['h2f']['min_cos_phi'] for b in B[S])
comb = {k: [max(per['v32'][k][j], per['wide'][k][j]) for j in range(KEYS[k])] for k in KEYS}
# float values replaced
FV = [json.load(open(f'{P2}/out/h1slab/win_{m}.json')) for m in range(32)]; WS = json.load(open(f'{P2}/out/h1slab/win_summary.json'))
FH = [json.load(open(f'{P2}/out_v33/h2f/win_{m}.json')) for m in range(32)]; HB = json.load(open(f'{P2}/out_v33/h2f/h2f_base.json'))
rel = lambda a, b: a / b - 1
cmp = dict(
    v32_windows_per_ball_max_rel=max(abs(rel(a['windows_max'][j], f['windows_max'][j])) for a, f in zip(B['v32'], FV) for j in range(4)),
    v32_shells_per_ball_max_rel=max(abs(rel(a['shells'][j], f['shells'][j])) for a, f in zip(B['v32'], FV) for j in range(3)),
    wide_G_per_ball_max_rel=max(abs(rel(a[k][j], f[k][j])) for a, f in zip(B['wide'], FH) for k in ('G_60_700', 'G_60_379p5', 'G_379p5_700') for j in range(4)),
    windows_at_same=all(a['windows_at'] == f['windows_at'] for a, f in zip(B['v32'], FV)),
    slab=dict(windows_max=dict(arb_v32=per['v32']['windows_max'], arb_combined=comb['windows_max'], float_win_summary=WS['windows_max'],
                               rel_combined=[rel(a, b) for a, b in zip(comb['windows_max'], WS['windows_max'])]),
              shells=dict(arb_v32=per['v32']['shells'], arb_combined=comb['shells'], float_win_summary=WS['shells'],
                          rel_combined=[rel(a, b) for a, b in zip(comb['shells'], WS['shells'])]),
              G_60_700=dict(arb_wide=per['wide']['G_60_700'], arb_combined=comb['G_60_700'], float_h2f_base=HB['G_60_700'],
                            rel_combined=[rel(a, b) for a, b in zip(comb['G_60_700'], HB['G_60_700'])])))
ins = {WA: sha(WA), f'{HERE}/m9_balls.py': sha(f'{HERE}/m9_balls.py'), f'{P2}/out/h1slab/win_summary.json': sha(f'{P2}/out/h1slab/win_summary.json'),
       f'{P2}/out_v33/h2f/h2f_base.json': sha(f'{P2}/out_v33/h2f/h2f_base.json')}
for S in B:
    for m in range(32): ins[f'{OM}/win_arb_{S}_{m}.json'] = sha(f'{OM}/win_arb_{S}_{m}.json')
out = dict(item='M9 item 1: arb window / shell / (H2f) base sums, slab 180 (32 t_c balls x 2 parameter-ball sets)', label='certified (arb, outward)',
           passed_on='combined = max over both sets', combined=comb, per_set=per, compare_float=cmp, reused_runs=reused,
           args=dict(argv=sys.argv[1:], env={}, inputs=ins))
fo = f'{OM}/win_arb_summary.json'; json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
print('combined windows', ['%.6e' % v for v in comb['windows_max']], 'shells', ['%.6e' % v for v in comb['shells']], 'G', ['%.6e' % v for v in comb['G_60_700']])
print('max rel vs float: windows %.2e shells %.2e G %.2e; rel_hh %.2e' % (cmp['v32_windows_per_ball_max_rel'], cmp['v32_shells_per_ball_max_rel'],
      cmp['wide_G_per_ball_max_rel'], max(per[S]['max_rel_rad_hh'] for S in B)))
