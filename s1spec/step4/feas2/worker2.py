"""F2b job queue.  usage: python3 worker2.py   (lock = mkdir; run several for parallelism)"""
import os, sys, json, time, shutil, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from run2 import *
import perturb2
D = '/home/claude/s1spec/step4/feas2'
LEV = [100, 100.3, 400, 401.2, 1600, 1604.8]
TJ = ['4.5', '3.0', '1.5', '0.3', '0']
def nref(tj):
    z = np.load(f'{OUT}/state_t{tj}.npz'); L = float(z['L']); N = int(z['N'])
    n = N
    while L / n > 0.05 / 6: n *= 2
    return n
JOBS = [(tj, 'ctl', 0, 0.0) for tj in TJ]
for rnd, dl in ((1, 1e-6), (2, 2e-6)):
    for sh in (('peak', 'fixed', 'far') if rnd == 1 else ('peak', 'far', 'fixed')):
        for tj in TJ: JOBS.append((tj, sh, rnd, dl))
def ref_state(tj):
    p = f'{OUT}/ref_t{tj}.npz'
    if not os.path.exists(p):
        S = load_state(f'{OUT}/state_t{tj}.npz')
        while S.N < nref(tj): S.refine()
        save(S, p + '.tmp.npz'); os.replace(p + '.tmp.npz', p)
    return p
for tj, sh, rnd, dl in JOBS:
    tag = f't{tj}_{sh}' + (f'_{rnd}' if sh != 'ctl' else '')
    try: os.mkdir(f'{D}/logs/lock_{tag}')
    except FileExistsError: continue
    t0 = time.time(); ref = ref_state(tj); start = f'{OUT}/{tag}_start.npz'
    if sh == 'ctl': shutil.copy(ref, start)
    else:
        info = perturb2.make(ref, start, sh, dl); json.dump(info, open(f'{D}/logs/pert_{tag}.json', 'w'), indent=1)
    log = open(f'{D}/logs/run_{tag}.log', 'a'); S = load_state(start); log.write(f"start {tag} k={S.kmax():.4f} N={S.N} t={S.t:.8f}\n")
    for lv in LEV:
        S.run_to(kstop=lv); save(S, f"{OUT}/{tag}_{lv:g}.npz")
        log.write(f"snap {lv} kmax={S.kmax():.4f} N={S.N} L={S.L:.8f} t={S.t:.12f} wall={time.time()-t0:.0f}s\n"); log.flush()
    log.write(f"DONE {time.time()-t0:.0f}s\n"); log.close(); os.remove(start)
