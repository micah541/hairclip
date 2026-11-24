"""Continue the clip 24:0.4:0.05:0.85 from a snapshot (resume), track k (not t).
usage: python3 cont.py TAG SNAP CDT TOL NMAX LEVELS(comma)"""
import sys, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4')
from hls import HLS
tag, snap, cdt, tol, nmax, levels = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), [float(x) for x in sys.argv[6].split(',')]
z = np.load(snap)
S = HLS(np.zeros(int(z['N'])), float(z['L']), z['X0'], omega=1, cdt=cdt, tail_tol=tol, Nmax=nmax)
S.th_hat = z['th_hat']; S.t = float(z['t'])
log = open(f"run_{tag}.log", "a"); t0 = time.time(); nxt = 0; steps = 0
log.write(f"start {tag} from {snap} k={S.kmax():.2f} N={S.N} L={S.L:.6f} cdt={cdt} tol={tol} Nmax={nmax}\n"); log.flush()
while nxt < len(levels):
    km = S.kmax()
    if km >= levels[nxt]:
        np.savez(f"snap_{tag}_{int(levels[nxt])}.npz", th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=km, steps=steps)
        log.write(f"snap {levels[nxt]:.0f} k={km:.3f} N={S.N} L={S.L:.8f} tail={S.tail():.2e} steps={steps} wall={time.time()-t0:.0f}s\n"); log.flush(); nxt += 1; continue
    S.run_to(maxsteps=50); steps += 50
    if steps % 1000 == 0: np.savez(f'ckpt_{tag}.npz', th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=km, steps=steps)
    if steps % 2000 == 0: log.write(f"  k={km:.2f} N={S.N} tail={S.tail():.2e} steps={steps} wall={time.time()-t0:.0f}s\n"); log.flush()
log.write("DONE\n"); log.close()
