"""T5 float check 3: (1) theta-form of the projected form with FIXED S1 constraints Phi = {chi p, chi n1, psi3}_S1, directions G = base modes;
(2) far-zone pointwise rate of the level-0 W-energy of v = e' over the base (W=(1+x^2)^-2 and x^-2), whole half curve; (3) base size in the far zone."""
import sys, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
from f3geom import *; from f3fem import *
Sigma, h = 40.0, 0.1
mesh = Mesh(Sigma, h); M = mesh.mass(); K2 = mesh.stiff2(); xn = np.arange(mesh.ne + 1) * mesh.h
def cut(x, a=22.0, b=30.0):
    t = np.clip((x - a) / (b - a), 0, 1); return 1 - t**3 * (10 - 15 * t + 6 * t * t)
def fem_vec(x, f):
    sp = CubicSpline(x, f * cut(x)); d = np.zeros(mesh.ndof); d[0::2] = sp(xn); d[1::2] = sp(xn, 1); return d[mesh.free]
def gs(C):
    Q = []
    for c in C.T:
        for q in Q: c = c - (q @ M @ c) * q
        Q.append(c / np.sqrt(c @ M @ c))
    return np.array(Q).T
def proj_top(K, Phi, G, th0):
    Z = sla.null_space(Phi.T @ M)
    B = K - M @ G @ np.linalg.solve(Phi.T @ M @ G, Phi.T @ K)
    Bz = Z.T @ (B + th0 * K2) @ Z; w = sla.eigh(0.5 * (Bz + Bz.T), Z.T @ M @ Z, eigvals_only=True)
    return w[-1]
sp1 = s1_splines(Smax=Sigma + 5, dx=0.01); K1 = mesh.form(coeffs(mesh, sp1, use_tau=False))
x1 = np.arange(0, Sigma + 5, 0.01); g1 = s1_geometry(x1); th = CubicSpline(x1, g1['kap']).antiderivative()(x1)
w, vr = eig_nonsym(K1, M); psi3 = vr[:, 2].real
PhiS1 = gs(np.c_[fem_vec(x1, g1['p']), fem_vec(x1, -np.cos(th)), psi3])
ths = (0.0, 0.05, 1 / 9, 0.2)
print('theta0      :', ths)
print('S1          :', [round(proj_top(K1, PhiS1, PhiS1, t), 4) for t in ths])
for nm in sys.argv[1:]:
    S = Snap(f'/home/claude/s1spec/step4/feas2/snaps/p2_{nm}.npz'); lam0 = None
    g = numeric_geometry(f'/home/claude/s1spec/step4/feas2/snaps/p2_{nm}.npz', Smax=Sigma + 5, dsig=0.01)
    K = mesh.form(coeffs(mesh, as_splines(g), use_tau=True))
    G = gs(np.c_[fem_vec(g['x'], g['p']), fem_vec(g['x'], (np.conj(1j * g['T'])).real), fem_vec(g['x'], g['kap2'])])
    print(f"{nm:6s} proj  :", [round(proj_top(K, PhiS1, G, t), 4) for t in ths], flush=True)
    # ---- far zone, whole half curve
    lam = g['lam']; xend = 0.5 * S.L / lam
    gf = numeric_geometry(f'/home/claude/s1spec/step4/feas2/snaps/p2_{nm}.npz', Smax=min(xend, 4000.0) - 0.5, dsig=0.02)
    x = gf['x']; k, k1, k2, q, tau = gf['kap'], gf['kap1'], gf['kap2'], gf['q'], gf['tau']
    b2 = -k**2; b1 = -5 * k * k1 - q / 4 + tau; b0 = 0.25 - 3 * k1**2 - 4 * k * k2
    sb2, sb1 = CubicSpline(x, b2), CubicSpline(x, b1)
    for name, m in (('x^-4', 4), ('x^-2', 2)):
        W = (1 + x**2) ** (-m / 2); sW = CubicSpline(x, W)
        # multiplier of v^2 from the first-order and zeroth-order terms and from the weight:  m0 = 1/2 (W b2)''/W - 1/2 (W (b2'+b1))'/W + b1' + b0 ; -int W v v'''' = -int W v''^2 + 2 int W'' v'^2 - 1/2 int W'''' v^2
        Wb2 = CubicSpline(x, W * b2); Wc = CubicSpline(x, W * (sb2(x, 1) + b1))
        m0 = 0.5 * Wb2(x, 2) / W - 0.5 * Wc(x, 1) / W + sb1(x, 1) + b0 - 0.5 * sW(x, 4) / W
        c1 = 2 * sW(x, 2) / W - b2            # coefficient of v'^2 (bad if >0):  2 W''/W + kappa^2
        # worst case with v'^2 <= |v||v''|:  c1 v'^2 <= v''^2/2 + c1^2/2 v^2  (pointwise heuristic)
        rate = m0 + 0.5 * np.maximum(c1, 0) ** 2
        for lo, hi in ((16.25, 60), (60, 1e9)):
            mk = (x >= lo) & (x <= hi) & (x <= x[-1] - 1.0)
            if mk.any(): print(f"   {name} zone [{lo},{min(hi, x[-1]):.0f}]: sup m0 = {m0[mk].max():+.4f}  sup rate(heur) = {rate[mk].max():+.4f}  (asymptotic {(1-m)/8:+.4f})   sup|b1+x/4|/x = {np.max(np.abs(b1[mk] + x[mk] / 4) / x[mk]):.3e}")
    for lo, hi in ((10.25, 16.25), (16.25, 34.25), (34.25, 60), (60, 1e9)):
        mk = (x >= lo) & (x <= hi)
        if mk.any(): print(f"   base on [{lo},{min(hi, x[-1]):.0f}]: sup|kap|={np.abs(k[mk]).max():.3e} sup|kap'|={np.abs(k1[mk]).max():.3e} sup|kap''|={np.abs(k2[mk]).max():.3e} sup x|kap|={np.max(x[mk]*np.abs(k[mk])):.3f} sup|b0'| x={np.max(np.abs(CubicSpline(x,b0)(x,1))[mk]*x[mk]):.3e}")
    print(f"   lam={lam:.4e} xend={xend:.1f}")
