import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2')
from handoff2 import *
from scipy.optimize import least_squares

class Mirror:
    """curve reflected about the line y = y0 (strand 'lower' becomes 'upper')."""
    def __init__(self, S, y0): self.S, self.y0, self.L = S, y0, S.L
    def curve(self, up=1):
        al, th, Z, clo = self.S.curve(up); N = len(al); j = (-np.arange(N)) % N
        return al, np.pi - th[j], np.conj(Z[j] - 1j * self.y0) + 1j * self.y0, clo

def fitc(S, xc0, lam0, xe, yfree):
    n = int(round(xe / HX)) + 1; sw = np.sqrt(wh[:n]); B = np.column_stack([PSI3[:n], H2[:n]]) * sw[:, None]
    def res(p, ret=False):
        xc = xc0 + (p[0] + 1j * (p[2] if yfree else 0.0)) * lam0
        u = hermite_graph(S, xc, lam0 * np.exp(p[1]), n)
        if not np.all(np.isfinite(u)): return np.full(n, 1e3)
        y = u * sw; c, *_ = np.linalg.lstsq(B, y, rcond=None); r = y - B @ c
        return (r, c, u) if ret else r
    p0 = [0.0, 0.0, 0.0] if yfree else [0.0, 0.0]
    sol = least_squares(res, p0, x_scale=[1e-3, 1e-3, 1e-2][:len(p0)], diff_step=1e-6, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    r, c, u = res(sol.x, True)
    xc = xc0 + (sol.x[0] + 1j * (sol.x[2] if yfree else 0.0)) * lam0
    return dict(xc=xc, lam=lam0 * np.exp(sol.x[1]), alpha=c[0], beta2=c[1], r2=float(np.linalg.norm(r)), u=u)

def run(f, fo):
    S, z = load(f); Zs, _ = strand(S); tip = Zs[0]
    al, th = S.theta(8); N = len(al); kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    lam0 = K0 / kt; xc0 = tip.real + D_S1 * lam0
    # chart end at the start modulation (axis through the tip)
    xcA = xc0 + 1j * tip.imag
    ug = hermite_graph(S, xcA, lam0, int(round(XT / HX)) + 1); fin = np.isfinite(ug); xend = xs[:len(ug)][fin][-1]; xe = min(XT, xend - 0.3)
    # local asymmetry: curve vs its mirror about y = tip.imag, over the window |Z - tip| < 30 lam
    al1, th1_, Z, _ = S.curve(1); Nn = len(al1); j = (-np.arange(Nn)) % Nn
    Zm = np.conj(Z[j] - 1j * tip.imag) + 1j * tip.imag
    win = np.abs(Z - tip) < 30 * lam0
    asym_loc = np.max(np.abs(Z - Zm)[win]) / 2 / lam0
    asym_all = np.max(np.abs(Z - Zm)) / 2
    # baseline (y-axis y=0), upper & lower strands in the same frame, odd/even decomposition
    fB = fitc(S, xc0 + 0j, lam0, xe, False)          # real xc: axis y = 0 (as in handoff2.py)
    n = len(fB['u']); u_up = fB['u']
    u_low = hermite_graph(Mirror(S, 0.0), fB['xc'], fB['lam'], n)
    ok = np.isfinite(u_low); m = min(n, ok.sum())
    od = 0.5 * (u_up - u_low)[:m]; ev = 0.5 * (u_up + u_low)[:m]
    cth = tau1.real[:m]; swm = wh[:m]
    cc = np.sum(swm * od * cth) / np.sum(swm * cth**2); rel = np.sqrt(np.sum(swm * (od - cc * cth)**2) / np.sum(swm * od**2))
    # frame fits: B) axis through the tip (dy = 0), C) free dy
    fT = fitc(S, xcA, lam0, xe, False); fC = fitc(S, xcA, lam0, xe, True)
    # lower strand with the tip-axis frame
    uL = hermite_graph(Mirror(S, tip.imag), fT['xc'], fT['lam'], n)
    B = np.column_stack([PSI3[:n], H2[:n]]); okL = np.isfinite(uL)
    cL, *_ = np.linalg.lstsq(B[okL] * np.sqrt(wh[:n][okL])[:, None], (uL * np.sqrt(wh[:n]))[okL], rcond=None)
    rL = np.linalg.norm((uL - B @ cL)[okL] * np.sqrt(wh[:n][okL]))
    line = (f"{float(z['kmax']):8.1f} lam={fT['lam']:.3e} tip.imag={tip.imag:+.3e} (={tip.imag/lam0:+.3e} lam) asym_loc={asym_loc:.2e} lam asym_max={asym_all:.2e} | "
            f"odd part (u_up-u_low)/2 vs cos(th1): coef={cc:+.3e} rel.resid={rel:.2e} ||odd||={np.sqrt(np.sum(swm*od**2)):.2e} ||even-fit w||={fB['r2']:.2e} | "
            f"r2: y=0 frame {fB['r2']:.3e}; axis-through-tip {fT['r2']:.3e}; free dy {fC['r2']:.3e} (dy={fC['xc'].imag/fT['lam']-tip.imag/fT['lam']:+.2e}, alpha={fC['alpha']:+.3e}, "
            f"b2/a2={fC['beta2']/fC['alpha']**2:+.3f}, r2/|a|^3={fC['r2']/abs(fC['alpha'])**3:.2f}); lower strand r2 {rL:.3e}")
    print(line, flush=True); fo.write(line + '\n'); fo.flush()
if __name__ == "__main__":
    with open(sys.argv[1], 'w') as fo:
        for f in sys.argv[2:]: run(f, fo)
