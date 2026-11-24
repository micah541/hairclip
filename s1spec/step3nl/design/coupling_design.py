"""FLOAT design evaluator for the inner/intermediate coupling (v = zeta w' energy in W = (1+x^2)^-P on the overlap).
Strip (inner cutoff chi_32, sigma in [29,32]) and inner ramp of zeta on [XA, XB] are written, after integrating by
parts, as  int A w^2 + B w'^2  (strip, rho-weighted)  and  int A_z g^2 + B_z g'^2  (ramp, W-weighted, g = w').
  strip:  int B w'^2 <= sB ||v||^2_strip,  sB = sup B/W;   int A+ w^2 <= 2 int A+ dx K_rho(xh)(..E..) + 2 (int A+ Phi) ||v||^2
          Phi(x) = int_xh^x W^-1  (Hardy from xh = XB).   mu_min = (sB + sA) / (thI * mu1)
  ramp:   mu (a_f ||w~'||^2 + b_f ||w~''||^2) <= thE (om0 ||w~||^2 + 1/9 ||w~''||^2), ||w~'||^2 <= ||w~|| ||w~''|| + c1 ||w~||^2
          mu_max = thE min_t (om0 + t^2/9) / (a_f (t + c1) + b_f t^2)."""
import sympy as sp, numpy as np, itertools, sys
X = sp.symbols('X'); STIP = -2.2490776510
def smooth_down(z, z0, z1):            # 1 for z <= z0, 0 for z >= z1 (C^3 smoothstep)
    u = (z1 - z) / (z1 - z0); return u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)
def ibp_coeffs(wt_expr, cut, drift):
    """returns A, B with  -<[L,cut] f, cut f>_wt = int A f^2 + B f'^2  (weight wt = e^{-psi}, k ~ 0)"""
    c1, c2, c3, c4 = [sp.diff(cut, X, j) for j in (1, 2, 3, 4)]
    g1, g2, g3, g4 = wt_expr * cut * c1, wt_expr * cut * c2, wt_expr * cut * c3, wt_expr * cut * c4
    A = -2 * sp.diff(g1, X, 3) + 3 * sp.diff(g2, X, 2) - 2 * sp.diff(g3, X) + g4 + drift / 4 * g1
    B = 6 * sp.diff(g1, X) - 6 * g2
    return A, B
rho = sp.exp(-(1 + X**2) ** sp.Rational(2, 3) / 4)
q = X + STIP
OM0, C1 = 1 / 32 - 0.01607, 0.5804
def evaluate(P, XA, XB, mu1_min, mu1_strip, thE=1/3, thI=1/3, xh=None):
    W = (1 + X**2) ** (-sp.Rational(P).limit_denominator(100))
    # strip: chi_32 in sigma -> x = sigma - STIP
    chi = smooth_down(X + STIP, 29, 32)
    A, B = ibp_coeffs(rho, chi, q)
    xs = np.linspace(29 - STIP, 32 - STIP, 3001)
    Af, Bf = sp.lambdify(X, A, 'numpy')(xs), sp.lambdify(X, B, 'numpy')(xs)
    Wn = sp.lambdify(X, W, 'numpy'); rhon = sp.lambdify(X, rho, 'numpy')
    sB = np.max(np.maximum(Bf, 0) / Wn(xs))
    xh = XB if xh is None else xh
    Phi = lambda x: np.array([np.trapezoid(1 / Wn(np.linspace(xh, xi, 400)), np.linspace(xh, xi, 400)) for xi in np.atleast_1d(x)])
    Ap = np.maximum(Af, 0)
    sA = 2 * np.trapezoid(Ap * Phi(xs), xs)
    Kxh = 1 / rhon(xh)
    sA0 = 2 * np.trapezoid(Ap, xs) * Kxh                      # times (||w~||^2 + 2||w~|| ||w~'||): goes to E's budget
    mu_min = sB / (thI * mu1_strip) + sA / (thI * mu1_min)
    # ramp: zeta up from XA to XB
    zeta = 1 - smooth_down(X, XA, XB)
    Az, Bz = ibp_coeffs(W, zeta, q)
    xr = np.linspace(XA, XB, 3001)
    Azf, Bzf = sp.lambdify(X, Az, 'numpy')(xr), sp.lambdify(X, Bz, 'numpy')(xr)
    af = np.max(np.maximum(Azf, 0) / rhon(xr)); bf = np.max(np.maximum(Bzf, 0) / rhon(xr))
    ts = np.logspace(-4, 3, 2000)
    mu_max = thE * np.min((OM0 + ts**2 / 9) / (af * (ts + C1) + bf * ts**2))
    return dict(P=P, XA=XA, XB=XB, sB=sB, sA=sA, sA0=sA0, mu_min=mu_min, af=af, bf=bf, mu_max=mu_max, ratio=mu_max / mu_min)
if __name__ == "__main__":
    # certified mu1 (agent_lint): (P, XA) -> (poly min, blend min); strip x in [31.25, 34.25] lies in the blend [30.25, 38.25]
    MU = {(3, 12): (0.2435, 0.562), (2.5, 12): (0.2553, 0.450), (2.5, 10): (0.1011, 0.450), (2, 12): (0.2210, 0.336),
          (2, 10): (0.1415, 0.336), (1.5, 12): (0.1555, 0.221), (1.5, 10): (0.1178, 0.221), (1.5, 8): (0.0338, 0.221)}
    for (P, XA), (m1, m1s) in MU.items():
        for wdt in (2, 4, 6, 8):
            r = evaluate(P, XA, XA + wdt, m1, m1s)
            print(f"P={P} ramp=[{XA},{XA+wdt}]  sB={r['sB']:.2e} sA={r['sA']:.2e} sA0={r['sA0']:.1e}  mu_min={r['mu_min']:.2e}  "
                  f"a_f={r['af']:.2e} b_f={r['bf']:.2e} mu_max={r['mu_max']:.2e}  ratio={r['ratio']:.2f}", flush=True)
