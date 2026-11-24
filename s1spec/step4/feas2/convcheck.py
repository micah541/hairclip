"""time-step floor check: from-scratch base with cdt (arg) and MACRO clippings; snapshots at k = 100, 400, 1600"""
import sys, time, numpy as np; sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
import run2
from run2 import *
cdt = float(sys.argv[1]); tag = sys.argv[2]; tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
S = build_base(); S.cdt = cdt; S.tail_tol = tol; log = open(f'logs/conv_{tag}.log', 'w'); t0 = time.time()
for te in [0.3, 0.5, 1.5, 2.0, 3.0, 4.0, 4.5, 4.59]:
    S.run_to(tend=te)
for k in [100, 400, 1600]:
    S.run_to(kstop=k); save(S, f'{OUT}/conv_{tag}_{k}.npz'); log.write(f"k={k} kmax={S.kmax():.4f} t={S.t:.12f} N={S.N} wall={time.time()-t0:.0f}s\n"); log.flush()
