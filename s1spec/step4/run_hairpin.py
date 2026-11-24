"""Run the hairpin, track the double points c1, c2 (on the axis), pocket and lens areas; snapshots near blow-up."""
import sys, time, numpy as np
from hls import HLS
from hairpin import hairpin_theta

def axis_diag(S, up=2):
    al, th, Z, clo = S.curve(up)
    N = len(al); half = slice(1, N // 2)
    y = Z.imag
    idx = np.where(np.sign(y[1:N // 2]) != np.sign(y[2:N // 2 + 1]))[0] + 1   # crossings in (0, pi)
    cr = []
    for i in idx:
        if i + 1 >= N // 2 - 2 or i <= 2: continue          # skip the self-symmetric points alpha = 0, pi
        f = y[i] / (y[i] - y[i + 1])
        cr.append((al[i] + f * (al[i + 1] - al[i]), Z.real[i] + f * (Z.real[i + 1] - Z.real[i])))
    out = {'clo': clo, 'cross': cr, 'tip_x': Z.real[N // 2]}
    if len(cr) >= 2:
        cr = sorted(cr)                                  # by alpha; c2 = closest to pi
        a2, x2 = cr[-1]; a1, x1 = cr[-2]
        def loop_area(a_lo, a_hi):
            m = (al > a_lo) & (al < a_hi)
            z = np.concatenate([[x2 if a_lo == a2 else x1], Z[m]])
            return 0.5 * np.sum((z * np.roll(np.conj(z), 1)).imag)
        m = (al > a2) & (al < 2 * np.pi - a2)
        z = np.concatenate([[x2 + 0j], Z[m], [x2 + 0j]])
        AX = 0.5 * np.sum(np.imag(np.conj(z[:-1]) * z[1:]))
        m1 = (al > a1) & (al < a2); m2 = (al > 2 * np.pi - a2) & (al < 2 * np.pi - a1)
        z = np.concatenate([[x1 + 0j], Z[m1], [x2 + 0j], Z[m2], [x1 + 0j]])
        AL = 0.5 * np.sum(np.imag(np.conj(z[:-1]) * z[1:]))
        out.update(c1=x1, c2=x2, a1=a1, a2=a2, AX=AX, AL=AL)
    return out

if __name__ == "__main__":
    tag = sys.argv[1]; cdt = float(sys.argv[2]); N0 = int(sys.argv[3]); kfinal = float(sys.argv[4])
    import os
    th, L, X0, tot, clo = hairpin_theta(N0, arm=float(os.environ.get("ARM", "12")), a_moll=float(os.environ.get("AMOLL", "0.025")), W=float(os.environ.get("WID", "1.0")), w=float(os.environ.get("SLOT", "0.1")))
    S = HLS(th, L, X0, omega=1, cdt=cdt, tail_tol=1e-10, Nmax=1 << 18)
    A0 = S.area()
    log = open(f"hp_{tag}.log", "w")
    snaps_k = [12, 15, 20, 25, 35, 50, 70, 100, 140, 200, 280, 400, 560, 800]
    nsnap = 0; t0 = time.time(); last = 0.0; AX0 = None
    def report(S, force=False):
        global nsnap, last
        km = S.kmax()
        d = axis_diag(S)
        line = (f"t={S.t:.10f} N={S.N} L={S.L:.6f} dA={S.area()/A0-1:+.1e} kmax={km:9.3f} "
                f"ncross={len(d['cross'])} tipx={d['tip_x']:.5f}")
        if 'c2' in d:
            line += f" c1={d['c1']:.5f} c2={d['c2']:.5f} |X|={abs(d['AX']):.6e} |Lam|={abs(d['AL']):.6e}"
        log.write(line + f" wall={time.time()-t0:.0f}s\n"); log.flush()
        global AX0
        if 'c2' in d and AX0 is None: AX0 = abs(d['AX'])
        late = ('c2' in d) and (AX0 is not None) and abs(d['AX']) < 0.3 * AX0
        if late and nsnap < len(snaps_k) and km >= snaps_k[nsnap]:
            np.savez(f"snap_{tag}_{snaps_k[nsnap]}.npz", th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t, kmax=km,
                     c1=d['c1'], c2=d['c2'], AX=d['AX'], AL=d['AL'])
            nsnap += 1
    # phase 1: dense reporting
    step = 0
    while True:
        km = S.kmax()
        if km > kfinal and nsnap > 0: break
        if S.t > 50.0: break
        S.run_to(maxsteps=20)
        step += 20
        report(S)
        if step % 2000 == 0:
            np.savez(f"ckpt_{tag}.npz", th_hat=S.th_hat, N=S.N, L=S.L, X0=S.X0, t=S.t)
    report(S); log.write("DONE\n"); log.close()
