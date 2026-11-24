"""Job queue worker for F2 perturbation runs (lock = mkdir).  usage: python3 worker.py"""
import os, sys, json, time
import numpy as np
from run_f2 import *
import perturb
OUTD = '/home/claude/s1spec/step4/feas'
JOBS = [('ctl400', 400, None, 0.0)]
for kj in (30, 100):
    for sh in ('tip', 'chamber'):
        JOBS.append((f'p{kj}_{sh}_1', kj, sh, 1e-6))
for kj in (30, 100):
    for sh in ('tip', 'chamber'):
        JOBS.append((f'p{kj}_{sh}_2', kj, sh, 2e-6))
for kj in (10, 400):
    for sh in ('tip', 'chamber'):
        for a, d in ((1, 1e-6), (2, 2e-6)):
            JOBS.append((f'p{kj}_{sh}_{a}', kj, sh, d))
for tag, kj, sh, d in JOBS:
    try: os.mkdir(f'{OUTD}/logs/lock_{tag}')
    except FileExistsError: continue
    t0 = time.time()
    start = f'{OUT}/{tag}_{kj}.npz'
    if sh is None:
        import shutil; shutil.copy(f'{OUT}/base_{kj}.npz', start)
    else:
        info = perturb.make(f'{OUT}/base_{kj}.npz', start, sh, d)
        json.dump({k: float(v) if not isinstance(v, str) else v for k, v in info.items()}, open(f'{OUTD}/logs/pert_{tag}.json', 'w'), indent=1)
    log = open(f'{OUTD}/logs/run_{tag}.log', 'a')
    S = load_state(start); log.write(f"start {tag} k={S.kmax():.4f} N={S.N}\n")
    go(S, tag, kj, log); log.write(f"DONE {time.time()-t0:.0f}s\n"); log.close()
