"""F2 runner. usage:
  python3 run_f2.py base                       -> base clip run from scratch, snapshots at LEVELS (default 10,30,100,400,1600)
  python3 run_f2.py resume TAG SNAPFILE KFROM  -> resume from SNAPFILE (any state incl. perturbed), run through remaining levels > KFROM
Same HLS settings as design/run3.py (cdt 1e-2, tail_tol 1e-10, N0 32768, Nmax 2^21).  Each segment is one run_to(kstop=level) call,
so a resumed run that starts exactly at a snapshot reproduces the base run bit-for-bit (hprev is stored in the snapshot)."""
import os, sys, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec'); sys.path.insert(0, '/home/claude/s1spec/step4/design')
from hls import HLS
LEVELS = [float(x) for x in os.environ.get('LEVELS', '10,30,100,400,1600').split(',')]
OUT = '/home/claude/s1spec/step4/feas/snaps'
CDT, TOL = 1e-2, 1e-10

def save(S, path, extra=None):
    np.savez(path, th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=S.kmax(), hprev=getattr(S, 'hprev', np.nan), **(extra or {}))

def load_state(path):
    z = np.load(path)
    S = HLS(np.zeros(int(z['N'])), float(z['L']), z['X0'], omega=1, cdt=CDT, tail_tol=TOL, Nmax=1 << 21)
    S.th_hat = z['th_hat']; S.t = float(z['t'])
    hp = float(z['hprev']) if 'hprev' in z.files else np.nan
    if np.isfinite(hp): S.hprev = hp
    return S

def build_base():
    os.environ['CHAMBER'] = '0.4'
    import hairpin
    arm, slot, W = 24.0, 0.05, 0.85
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
    th, L, X0, tot, clo = hairpin.hairpin_theta(32768, arm=arm, a_moll=0.025, W=W, w=slot)
    return HLS(th, L, X0, omega=1, cdt=CDT, tail_tol=TOL, Nmax=1 << 21)

def go(S, tag, kfrom, log):
    t0 = time.time()
    for lv in LEVELS:
        if lv <= kfrom: continue
        S.run_to(kstop=lv)
        p = f"{OUT}/{tag}_{int(lv)}.npz"; save(S, p)
        log.write(f"snap {lv:.0f} kmax={S.kmax():.4f} N={S.N} L={S.L:.8f} t={S.t:.12f} wall={time.time()-t0:.0f}s\n"); log.flush()

if __name__ == '__main__':
    if sys.argv[1] == 'base':
        log = open(f"{os.path.dirname(OUT)}/logs/run_base.log", 'a')
        S = build_base(); log.write(f"base L0={S.L:.6f} kmax0={S.kmax():.3f} N={S.N}\n")
        S.run_to(tend=4.0)
        log.write(f"phase A done t={S.t:.6f} kmax={S.kmax():.3f}\n"); log.flush()
        save(S, f"{OUT}/base_A.npz")
        go(S, 'base', 0, log); log.write("DONE\n")
    else:
        tag, snap, kfrom = sys.argv[2], sys.argv[3], float(sys.argv[4])
        log = open(f"{os.path.dirname(OUT)}/logs/run_{tag}.log", 'a')
        S = load_state(snap); log.write(f"resume {tag} from {snap} k={S.kmax():.4f} N={S.N}\n")
        go(S, tag, kfrom, log); log.write("DONE\n")
