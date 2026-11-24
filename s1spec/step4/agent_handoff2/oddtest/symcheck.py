"""asymmetry of theta-spectrum: symmetric <=> mean(th) = -pi/2 and Re th_hat[k>=1] = 0 (th periodic part, theta = alpha + th)."""
import sys, os, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4'); sys.path.insert(0, '/home/claude/s1spec')
def asym(th_hat, N):
    mean = th_hat[0].real / N
    c = th_hat[1:].real * 2 / N; s = -th_hat[1:].imag * 2 / N     # th = mean + sum c cos + s sin ... (rfft sign convention)
    return mean + np.pi / 2, np.max(np.abs(c)), np.max(np.abs(s))
if __name__ == "__main__":
    os.environ['CHAMBER'] = '0.3'
    import hairpin
    from hls import HLS
    def defect(d):
        os.environ['CHDW'] = repr(float(d))
        s, k, L = hairpin.curvature_fn(1 << 20, 24, 0.025, 1.0, 0.05)
        ds = L / len(s); th = -np.pi / 2 + np.concatenate([[0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
        return float((np.sum(np.exp(1j * th)) * ds).imag)
    d0, d1 = 0.0, 0.003; f0, f1 = defect(d0), defect(d1)
    for _ in range(6):
        d2 = d1 - f1 * (d1 - d0) / (f1 - f0); d0, f0, d1, f1 = d1, f1, d2, defect(d2)
    os.environ['CHDW'] = repr(float(d1))
    th, L, X0, tot, clo = hairpin.hairpin_theta(16384, arm=24.0, a_moll=0.025, W=1.0, w=0.05)
    S = HLS(th, L, X0, omega=1)
    print("t=0: mean(th)+pi/2 = %.3e, max|cos coeff| = %.3e, max|sin coeff| = %.3e" % asym(S.th_hat, S.N)); print("X0 =", X0)
    al, thh, Z, clo2 = S.curve(1); N = len(al)
    m = np.conj(Z[(-np.arange(N)) % N]) ; y0 = Z[0].imag
    print("t=0 curve: max|Z - mirror(Z)|/2 =", np.max(np.abs(Z - (np.conj(Z[(-np.arange(N)) % N] - 1j*y0) + 1j*y0)))/2, " closure", clo2)
    for f in sys.argv[1:]:
        z = np.load(f); print(os.path.basename(f), "kmax %.0f" % float(z['kmax']), "mean+pi/2=%.3e cos=%.3e sin=%.3e  X0y=%.3e" % (*asym(z['th_hat'], int(z['N'])), float(z['X0'][1])))
