"""even-class remainder: fit (u_up(x)+u_low(x))/2 (axis through the tip) -> removes any rotation / perpendicular translation (odd) artefact."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2'); sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2/oddtest')
from handoff2 import *
from odd import Mirror
from scipy.optimize import least_squares
def even_r2(f):
    S, z = load(f); Zs, _ = strand(S); tip = Zs[0]; al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2]); lam0 = K0 / kt; xcA = tip.real + D_S1 * lam0 + 1j * tip.imag
    n = int(round(27.7 / HX)) + 1; sw = np.sqrt(wh[:n]); B = np.column_stack([PSI3[:n], H2[:n]]) * sw[:, None]; Mi = Mirror(S, tip.imag)
    def res(p, ret=False):
        xc = xcA + p[0] * lam0; lam = lam0 * np.exp(p[1])
        a = hermite_graph(S, xc, lam, n); b = hermite_graph(Mi, xc, lam, n)
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): return np.full(n, 1e3)
        u = 0.5 * (a + b); y = u * sw; c, *_ = np.linalg.lstsq(B, y, rcond=None); r = y - B @ c
        return (r, c, 0.5 * (a - b)) if ret else r
    sol = least_squares(res, [0.0, 0.0], x_scale=[1e-3, 1e-3], diff_step=1e-6, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    r, c, od = res(sol.x, True)
    return float(z['kmax']), c[0], np.linalg.norm(r), float(np.sqrt(np.sum(wh[:n] * od**2)))
if __name__ == "__main__":
    for f in sys.argv[1:]:
        k, a, r2, odn = even_r2(f)
        print(f"{f:34s} k={k:7.1f} alpha={a:+.3e} a*k={a*k:+.3f} EVEN r2={r2:.3e} r2/|a|^3={r2/abs(a)**3:7.1f}  ||odd part||={odn:.2e}", flush=True)
