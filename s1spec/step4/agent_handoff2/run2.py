"""Keyhole CDF run (same initial curve as keyhole_run.py) with adjustable cdt / tail_tol / N0 / snapshot levels; can resume from a snapshot.
usage: python3 run2.py TAG CHAMBER cdt N0 tail_tol kfinal [resume_snap]   (levels: env LEVELS, comma list)"""
import os, sys, time, subprocess, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec')
tag, R, cdt, N0, tol, kfinal = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
resume = sys.argv[7] if len(sys.argv) > 7 else None
os.environ['CHAMBER'] = R
import hairpin
from hls import HLS
from run_hairpin import axis_diag
levels = [float(x) for x in os.environ.get('LEVELS', '12,15,20,25,35,50,70,100,140,200,280,400').split(',')]
if resume:
    z = np.load(resume); S = HLS(np.zeros(int(z['N'])), float(z['L']), z['X0'], omega=1, cdt=cdt, tail_tol=tol, Nmax=1 << 21)
    S.th_hat = z['th_hat']; S.t = float(z['t'])
else:
    def defect(d):
        os.environ['CHDW'] = repr(float(d))
        s, k, L = hairpin.curvature_fn(1 << 20, 24, 0.025, 1.0, 0.05)
        ds = L / len(s); th = -np.pi / 2 + np.concatenate([[0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
        return float((np.sum(np.exp(1j * th)) * ds).imag)
    d0, d1 = 0.0, 0.003; f0, f1 = defect(d0), defect(d1)
    for _ in range(6):
        d2 = d1 - f1 * (d1 - d0) / (f1 - f0); d0, f0, d1, f1 = d1, f1, d2, defect(d2)
    os.environ.update(CHDW=repr(float(d1)))
    th, L, X0, tot, clo = hairpin.hairpin_theta(N0, arm=24.0, a_moll=0.025, W=1.0, w=0.05)
    S = HLS(th, L, X0, omega=1, cdt=cdt, tail_tol=tol, Nmax=1 << 21)
t0 = time.time(); nxt = 0
while nxt < len(levels) and levels[nxt] <= S.kmax() and resume: nxt += 1
log = open(f"run2_{tag}.log", "a")
while nxt < len(levels):
    km = S.kmax()
    if km >= levels[nxt]:
        try: d = axis_diag(S); c1, c2, AX, AL = d['c1'], d['c2'], d['AX'], d['AL']
        except Exception: c1 = c2 = AX = AL = np.nan
        np.savez(f"snap_{tag}_{int(levels[nxt])}.npz", th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=km, c1=c1, c2=c2, AX=AX, AL=AL)
        log.write(f"snap {levels[nxt]} kmax={km:.3f} N={S.N} t={S.t:.12f} wall={time.time()-t0:.0f}s\n"); log.flush()
        nxt += 1; continue
    if S.t > 50: break
    S.run_to(maxsteps=20)
log.write("DONE\n"); log.close()
