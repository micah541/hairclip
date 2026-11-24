"""H1 on slab 180: the 16 exact points t_k = -1 + (2k+1)/16 (h1.py p<t_k>: Krawczyk IFT, Lambda, rotation-mode term), as run_h1slab.sh does
for them (same H1TAG = _s180 and the same starting points H1REF, a linear interpolation between two reference parameter triples), with H1JOBS
runs in parallel.  Each run writes p2/out/t9_h1_p<t_k>_coh_s180.json and its log p2/out/h1slab/h1_p<t_k>.log; any failure is an error.
The derivative balls (h1dt.py) are not rerun here: they read only the points' parameter boxes (xc, lam, alpha), which h1_summary.py checks.
Resumable: the marker p2/out/h1slab/h1_points.pass exists while a pass of this stage is in progress (created at its start, removed when all 16
points are done).  A point whose output is newer than the marker and complete (Krawczyk ok, rotation-mode term, this H1REF) was computed in
this pass and is not rerun; so a killed pass continues from its finished points (listed in h1slab/h1_points_resume.log, not in the stage log).
env H1ONLY = 'k,k,...' runs only those points (chunks, e.g. in the foreground) and keeps the marker; the stage itself runs without it.
usage: python3 -I h1_points.py   (env COHERENT, S1PROF, S1TT as for h1.py; H1JOBS, default 2; H1ONLY)"""
import os, sys, json, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); O = os.path.normpath(os.path.join(HERE, '..', 'out'))
A = (1.04655427e-4, 0.99999491111, -1.17386340e-4); B = (1.00712248e-4, 0.99999510274, -1.12963590e-4)   # as run_h1slab.sh
JOBS = int(os.environ.get('H1JOBS', '2')); ONLY = os.environ.get('H1ONLY')
KS = [int(v) for v in ONLY.split(',')] if ONLY else list(range(16))
MARK = f'{O}/h1slab/h1_points.pass'
if not os.path.exists(MARK): open(MARK, 'w').write('pass started ' + time.strftime('%Y-%m-%d %H:%M:%S') + '\n')
T_PASS = os.path.getmtime(MARK)
def done(tk, ref):
    f = f'{O}/t9_h1_p{tk}_coh_s180.json'
    if not os.path.exists(f) or os.path.getmtime(f) <= T_PASS: return False
    try:
        d = json.load(open(f)); L = d['Lambda']
        return L['krawczyk']['ok'] and 'norms_EI' in L and d['args']['env'].get('H1REF') == ref and d['args']['env'].get('H1TAG') == '_s180'
    except Exception:
        return False
todo, rcs, resumed = [], {}, []
for k in KS:
    tk = -1 + (2 * k + 1) / 16; f = (tk + 1) / 1.9990234375
    ref = ','.join(repr(x + f * (y - x)) for x, y in zip(A, B))
    if done(tk, ref): rcs[tk] = 0; resumed.append(tk)
    else: todo.append((tk, dict(os.environ, H1TAG='_s180', H1REF=ref)))
if resumed:
    with open(f'{O}/h1slab/h1_points_resume.log', 'a') as rl: rl.write(time.strftime('%Y-%m-%d %H:%M:%S') + f' resumed (done in this pass): {resumed}\n')
run = []
while todo or run:
    while todo and len(run) < JOBS:
        tk, env = todo.pop(0); lf = open(f'{O}/h1slab/h1_p{tk}.log', 'w')
        run.append((tk, lf, subprocess.Popen([sys.executable, '-I', 'h1.py', f'p{tk}'], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)))
    time.sleep(1)
    for item in [r for r in run if r[2].poll() is not None]:
        tk, lf, pr = item; lf.close(); run.remove(item); rcs[tk] = pr.returncode
for tk in sorted(rcs): print(f'h1 p{tk}: rc {rcs[tk]}')                          # in t_k order (the runs finish in any order)
fails = [tk for tk in sorted(rcs) if rcs[tk] != 0]
print(f'{len(KS)} points, {len(fails)} failures: {fails}')
assert not fails
if not ONLY: os.remove(MARK)                                                     # the pass of this stage is complete
