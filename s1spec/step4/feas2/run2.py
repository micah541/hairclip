"""feas2 runner (copy of feas/run_f2.py with OUT changed).  usage:
  python3 run2.py macro            -> base from scratch, saves states at MACRO times
"""
import os, sys, time, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec'); sys.path.insert(0, '/home/claude/s1spec/step4/design')
from hls import HLS
OUT = '/home/claude/s1spec/step4/feas2/snaps'
CDT, TOL = 1e-2, 1e-10
MACRO = [0.3, 0.5, 1.5, 2.0, 3.0, 4.0, 4.5, 4.59]

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

if __name__ == '__main__':
    log = open('/home/claude/s1spec/step4/feas2/logs/run_macro.log', 'a')
    S = build_base(); save(S, f'{OUT}/state_t0.npz'); log.write(f"t=0 kmax={S.kmax():.4f} N={S.N} L={S.L:.6f}\n"); log.flush()
    t0 = time.time()
    for te in MACRO:
        S.run_to(tend=te); save(S, f'{OUT}/state_t{te}.npz')
        log.write(f"t={S.t:.8f} kmax={S.kmax():.4f} N={S.N} L={S.L:.6f} wall={time.time()-t0:.0f}s\n"); log.flush()
    log.write("DONE\n")
