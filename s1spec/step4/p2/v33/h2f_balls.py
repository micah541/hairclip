"""v3.3 task 6, (H2f) stage 2: the 32 t_c balls of slab 180, one subprocess per ball (h2f_win.py M, sequential, nice inherited), then
the slab-wide base values G_j = max over the balls (each ball's sum is rounded up).  Resumable: a ball is skipped if its win_M.json
exists and records the current sha256 of h2f_win.py and of its widened ball ift_M.json.  Writes p2/out_v33/h2f/h2f_base.json.
usage: python3 -I h2f_balls.py"""
import sys, os, json, math, hashlib, subprocess, time
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33/h2f'; HERE = os.path.dirname(os.path.abspath(__file__))
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
WIN = f'{HERE}/h2f_win.py'; done = []
for m in range(32):
    o, ift = f'{OV}/win_{m}.json', f'{OV}/ift_{m}.json'
    if os.path.exists(o):
        a = json.load(open(o))['args']['inputs']
        if a.get(WIN) == sha(WIN) and a.get(ift) == sha(ift): done.append(m); print(f'ball {m}: reused', flush=True); continue
    r = subprocess.run([sys.executable, '-I', '-B', WIN, str(m)], cwd=HERE, stdout=open(f'{OV}/win_{m}.log', 'w'), stderr=subprocess.STDOUT)
    assert r.returncode == 0, f'ball {m} failed'
    print(open(f'{OV}/win_{m}.log').read().strip().splitlines()[-1], flush=True)
B = [json.load(open(f'{OV}/win_{m}.json')) for m in range(32)]
assert sorted(b['ball'] for b in B) == list(range(32)) and all(b['fidelity_vs_win_energies_45_390']['identical'] for b in B)
cov = sorted((b['tc_ball'] for b in B)); assert cov[0][0] <= -1 and cov[-1][1] >= 1 and all(cov[i][1] >= cov[i + 1][0] for i in range(31))
mx = lambda key: [max(b[key][j] for b in B) for j in range(4)]
out = dict(label='certified (arb cells, outward), base curve on the widened parameter balls; error part not included',
           G_60_700=mx('G_60_700'), G_60_379p5=mx('G_60_379p5'), G_379p5_700=mx('G_379p5_700'),
           worst_ball=[max(B, key=lambda b: b['G_60_700'][j])['ball'] for j in range(4)],
           info_379p5_700=dict(sup_f=max(b['info_379p5_700']['sup_f'] for b in B), sup_yfy=max(b['info_379p5_700']['sup_yfy'] for b in B),
                               sup_z=[max(b['info_379p5_700']['sup_z'][j] for b in B) for j in range(4)]),
           min_cos_phi=min(b['min_cos_phi'] for b in B), balls=32, tc_cover=[cov[0][0], cov[-1][1]], reused_balls=done,
           args=dict(argv=sys.argv[1:], env={}, inputs={**{f'{OV}/win_{m}.json': sha(f'{OV}/win_{m}.json') for m in range(32)}, WIN: sha(WIN)}))
json.dump(out, open(f'{OV}/h2f_base.json', 'w'), indent=1)
print('G_60_700 (base, max over 32 balls):', ['%.4e' % v for v in out['G_60_700']])
