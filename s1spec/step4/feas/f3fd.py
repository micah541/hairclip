"""Finite-difference validation of  delta V = L e  (without transport) against the nonlinear V of the displaced curve  x = Gamma + eps e N  (same sigma labels).
Nonlinear V computed from the displaced parametrization with exact formulas  a = 1 - eps e kappa, b = eps e',
   x_s = aT + bN, x_ss = (a' - b kappa)T + (a kappa + b')N, kappa_eps = (a(a kappa + b') - b(a' - b kappa))/v^3, v^2 = a^2 + b^2,
   kappa_ss = (1/v) d/ds[(1/v) d kappa_eps/ds]-in-sigma-derivatives (8th-order FD in sigma),  Y.N_eps with N_eps = i x_s/v,  V = -kappa_ss + (Y + eps e N).N_eps/4.
base: S1 (grid HX=0.01 of measure.py) or numeric curve geometry."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
from f3geom import *

def D1(f, h):   # 8th order central first derivative, interior only (edges set nan)
    c = np.array([1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280])
    out = np.full_like(f, np.nan)
    out[4:-4] = sum(c[j] * f[j:len(f) - 8 + j] for j in range(9)) / h
    return out

def run_fd(T, N, kap, kap1, kap2, q, p, Y, h, label, eps_list=(1e-3, 5e-4)):
    x = np.arange(len(kap)) * h
    e0 = np.exp(-x**2 / (2 * 2.0**2)) * (1 + 0.3 * np.cos(0.7 * x))       # even test function, smooth
    e1 = D1(e0, h); e2 = D1(e1, h); e3 = D1(e2, h); e4 = D1(e3, h)
    # analytic derivatives of e0 to avoid FD noise: use sympy-free formulas via FD on analytic expression with complex step? use FD (8th order) with careful h
    Le = -e4 + (-kap**2) * e2 + (-5 * kap * kap1 - q / 4) * e1 + (0.25 - 3 * kap1**2 - 4 * kap * kap2) * e0
    def Vnl(eps):
        e = eps * e0; ep = eps * e1
        a = 1 - e * kap; b = ep
        a1 = D1(a, h); b1 = D1(b, h)
        v = np.sqrt(a**2 + b**2)
        kk = (a * (a * kap + b1) - b * (a1 - b * kap)) / v**3
        k_s = D1(kk, h) / v; k_ss = D1(k_s, h) / v
        Ne = 1j * (a * T + b * N * 1.0) / v           # N_eps = i x_s / v ; x_s = a T + b N,  N = iT  -> use complex arithmetic below
        return k_ss, Ne, v
    out = []
    for eps in eps_list:
        res = []
        for sgn in (+1, -1):
            e = sgn * eps * e0; ep = sgn * eps * e1
            a = 1 - e * kap; b = ep
            xs_ = a * T + b * (1j * T)                # x_sigma (complex)
            v = np.abs(xs_)
            a1 = D1(a, h); b1 = D1(b, h)
            kk = (a * (a * kap + b1) - b * (a1 - b * kap)) / v**3
            k_ss = D1(D1(kk, h) / v, h) / v
            Ne = 1j * xs_ / v
            Ye = Y + e * (1j * T)
            res.append(-k_ss + 0.25 * (Ye * np.conj(Ne)).real)
        dV = (res[0] - res[1]) / (2 * eps)
        sl = slice(20, len(x) - 20)
        err = np.max(np.abs(dV - Le)[sl]); sc = np.max(np.abs(Le)[sl])
        out.append((eps, err, sc))
        print(f"[{label}] eps={eps:g}: max|FD dV - L e| = {err:.3e}  (max|L e| = {sc:.3e}, rel {err/sc:.2e})", flush=True)
    return out

if __name__ == '__main__':
    # S1 base
    h = 0.01; x = np.arange(0, 20 + h / 2, h); g = s1_geometry(x)
    # S1 tangent: theta from curvature spline antiderivative, as in measure.py handoff
    th = -np.pi / 2 + CubicSpline(x, g['kap']).antiderivative()(x); T = np.exp(1j * th); Nn = 1j * T
    Y = g['q'] * T + g['p'] * Nn
    run_fd(T, Nn, g['kap'], g['kap1'], g['kap2'], g['q'], g['p'], Y, h, 'S1')
    # numeric base k=100 (axis-through-tip frame, tip-anchored)
    for k in (100, 1600):
        gn = numeric_geometry(f'/home/claude/s1spec/step4/feas/snaps/base_{k}.npz', Smax=20, dsig=0.01)
        hh = gn['x'][1] - gn['x'][0]
        run_fd(gn['T'], 1j * gn['T'], gn['kap'], gn['kap1'], gn['kap2'], gn['q'], gn['p'], gn['Y'], hh, f'num k={k}')
