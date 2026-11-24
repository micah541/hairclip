"""usage: python3 run3.py TAG KIND KFINAL  with KIND = clip:ARM:CHAMBER:SLOT[:W]  or  disk:Rd:SLOT:Rc:phi_deg[:Lw[:l0]]
Snapshots saved (after the double points exist) at LEVELS (env, default 100,200,400,800); coarse settings cdt=1e-2, tail_tol=1e-10, N0=16384."""
import os, sys, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec'); sys.path.insert(0, '/home/claude/s1spec/step4/design')
tag, kind, kfinal = sys.argv[1], sys.argv[2].split(':'), float(sys.argv[3])
cdt = float(os.environ.get('CDT', '1e-2')); tol = float(os.environ.get('TOL', '1e-10')); N0 = int(os.environ.get('N0', '16384'))
from hls import HLS
from run_hairpin import axis_diag
from scipy.spatial import cKDTree
if kind[0] == 'clip':
    arm, ch, slot = float(kind[1]), kind[2], float(kind[3]); W = float(kind[4]) if len(kind) > 4 else 1.0
    os.environ['CHAMBER'] = ch
    import hairpin
    def defect(d):
        os.environ['CHDW'] = repr(float(d))
        s, k, L = hairpin.curvature_fn(1 << 20, arm, 0.025, W, slot)
        ds = L / len(s); th = -np.pi / 2 + np.concatenate([[0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
        return float((np.sum(np.exp(1j * th)) * ds).imag)
    d0, d1 = 0.0, 0.003; f0, f1 = defect(d0), defect(d1)
    for _ in range(8):
        if abs(f1 - f0) < 1e-18: break
        d2 = d1 - f1 * (d1 - d0) / (f1 - f0); d0, f0, d1, f1 = d1, f1, d2, defect(d2)
    os.environ['CHDW'] = repr(float(d1))
    th, L, X0, tot, clo = hairpin.hairpin_theta(N0, arm=arm, a_moll=0.025, W=W, w=slot)
else:
    import curves
    Rd, slot, Rc, phi = float(kind[1]), float(kind[2]), float(kind[3]), np.radians(float(kind[4]))
    Lw = float(kind[5]) if len(kind) > 5 else None; l0 = float(kind[6]) if len(kind) > 6 else 0.3
    th, L, X0, beta, err = curves.disk_theta(N0, Rd, slot, Rc, phi, Lw, l0)
S = HLS(th, L, X0, omega=1, cdt=cdt, tail_tol=tol, Nmax=1 << 21)
al, t_, Z, clo = S.curve(2); sub = Z[::4]; tr = cKDTree(np.c_[sub.real, sub.imag]); n = len(sub)
dd, ii = tr.query(np.c_[sub.real, sub.imag], k=40)
sep = np.abs(np.arange(n)[:, None] - ii); sep = np.minimum(sep, n - sep) * (S.L / n)
mind = np.min(np.where(sep > 1.0, dd, np.inf))
levels = [float(x) for x in os.environ.get('LEVELS', '100,200,400,800').split(',')]
log = open(f"run3_{tag}.log", "w"); log.write(f"{tag} {sys.argv[2]} L0={S.L:.4f} N0={S.N} kmax0={S.kmax():.2f} closure={clo:.2e} min-gap(non-adjacent)={mind:.4f} area={S.area():.4f}\n"); log.flush()
t0 = time.time(); nxt = 0; AX0 = None; steps = 0
while nxt < len(levels):
    km = S.kmax()
    if km >= levels[nxt] and steps > 0:
        try:
            d = axis_diag(S); ok = 'c2' in d
        except Exception: ok = False
        if ok and AX0 is not None and abs(d['AX']) < 0.3 * AX0:
            np.savez(f"snap_{tag}_{int(levels[nxt])}.npz", th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=km, c1=d['c1'], c2=d['c2'], AX=d['AX'], AL=d['AL'])
            log.write(f"snap {levels[nxt]:.0f} kmax={km:.2f} N={S.N} L={S.L:.4f} t={S.t:.10f} wall={time.time()-t0:.0f}s |X|={abs(d['AX']):.3e} |Lam|={abs(d['AL']):.3e}\n"); log.flush(); nxt += 1; continue
    if steps % 400 == 0:
        try:
            d = axis_diag(S)
            if 'c2' in d and AX0 is None: AX0 = abs(d['AX'])
            if steps % 4000 == 0: log.write(f"  t={S.t:.8f} kmax={km:.2f} N={S.N} ncross={len(d['cross'])} " + (f"|X|={abs(d['AX']):.3e}" if 'c2' in d else '') + f" wall={time.time()-t0:.0f}s\n"); log.flush()
        except Exception: pass
    if S.t > 60 or time.time() - t0 > 2400: log.write("ABORT (time/t)\n"); break
    S.run_to(maxsteps=20); steps += 20
log.write("DONE\n"); log.close()
