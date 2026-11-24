import sys, time, numpy as np; sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from run2 import *
TS = []
KS = [2560, 3620, 5120, 6600]
log = open('logs/gen_states2.log', 'a'); t0 = time.time()
S = load_state(f"{OUT}/p2_k1810.npz")
for te in TS:
    S.run_to(tend=te); save(S, f'{OUT}/p2_t{te}.npz'); log.write(f"t={S.t:.8f} kmax={S.kmax():.4f} N={S.N} wall={time.time()-t0:.0f}s\n"); log.flush()
for k in KS:
    S.run_to(kstop=k); save(S, f'{OUT}/p2_k{k}.npz'); log.write(f"k={k} kmax={S.kmax():.4f} t={S.t:.12f} N={S.N} wall={time.time()-t0:.0f}s\n"); log.flush()
log.write('DONE\n')
