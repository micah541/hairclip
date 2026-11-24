"""Same as coupling_design.py but with the ACTUAL blended weight W = exp(-phi_I),
phi_I' = (1-s) 2P x/(1+x^2) + s (1/3) x (1+x^2)^(-1/3), s = S((x-A2)/(B2-A2)), phi_I(0) = 0 (agent_lint convention)."""
import numpy as np, sympy as sp, sys
from scipy.integrate import cumulative_trapezoid as ctz
sys.path.insert(0, '.'); import coupling_design as CD
STIP = CD.STIP
def S(t): t = np.clip(t, 0, 1); return t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3)
def weight(P, A2, B2, xmax=70, h=1e-3):
    x = np.arange(0, xmax + h / 2, h); s = S((x - A2) / (B2 - A2))
    dphi = (1 - s) * 2 * P * x / (1 + x**2) + s * x * (1 + x**2) ** (-1 / 3) / 3
    phi = ctz(dphi, x, initial=0); return x, np.exp(-phi), dphi
def evaluate(P, XA, XB, mu1_min, mu1_strip, A2=30.25, B2=38.25, thE=1/3, thI=1/3):
    x, W, dphi = weight(P, A2, B2)
    rho = np.exp(-(1 + x**2) ** (2 / 3) / 4)
    X = CD.X
    chi = CD.smooth_down(X + STIP, 29, 32)
    A, B = CD.ibp_coeffs(CD.rho, chi, CD.q)
    st = (x >= 29 - STIP) & (x <= 32 - STIP)
    Af, Bf = sp.lambdify(X, A, 'numpy')(x[st]), sp.lambdify(X, B, 'numpy')(x[st])
    sB = np.max(np.maximum(Bf, 0) / W[st])
    i0 = np.searchsorted(x, XB); Phi = np.zeros_like(x); Phi[i0:] = ctz(1 / W[i0:], x[i0:], initial=0)
    sA = 2 * np.trapezoid(np.maximum(Af, 0) * Phi[st], x[st])
    # ramp: needs W derivatives -> use the sympy formula with W replaced by the local polynomial (blend starts at A2 >> XB)
    Wp = (1 + X**2) ** (-sp.Rational(P).limit_denominator(100))
    zeta = 1 - CD.smooth_down(X, XA, XB)
    Az, Bz = CD.ibp_coeffs(Wp, zeta, CD.q)
    xr = np.linspace(XA, XB, 3001); rr = np.exp(-(1 + xr**2) ** (2 / 3) / 4)
    af = np.max(np.maximum(sp.lambdify(X, Az, 'numpy')(xr), 0) / rr); bf = np.max(np.maximum(sp.lambdify(X, Bz, 'numpy')(xr), 0) / rr)
    ts = np.logspace(-4, 3, 2000)
    mu_max = thE * np.min((CD.OM0 + ts**2 / 9) / (af * (ts + CD.C1) + bf * ts**2))
    mu_min = sB / (thI * mu1_strip) + sA / (thI * mu1_min)
    iR = lambda R: W[np.searchsorted(x, R)]
    return dict(sB=sB, sA=sA, mu_min=mu_min, af=af, bf=bf, mu_max=mu_max, ratio=mu_max / mu_min,
                rhoW_strip=np.max(rho[st] / W[st]), W40=iR(40), W60=iR(60), W70=iR(70), KW34=1 / iR(34.25))
if __name__ == "__main__":
    for (P, XA, XB, m1, m1s, A2, B2) in [(1.5, 10, 16, 0.1178, 0.221, 30.25, 38.25), (1.5, 10.25, 16.25, 0.1178, 0.221, 30.25, 38.25),
                                         (2, 10, 16, 0.1415, 0.336, 30.25, 38.25), (2, 10.25, 16.25, 0.1415, 0.336, 30.25, 38.25),
                                         (2.5, 10, 16, 0.1011, 0.450, 30.25, 38.25), (1.5, 10, 16, 0.1178, 0.221, 34.25, 42.25)]:
        r = evaluate(P, XA, XB, m1, m1s, A2, B2)
        print(f"P={P} ramp=[{XA},{XB}] blend=[{A2},{B2}]: rho/W strip {r['rhoW_strip']:.1e} sB {r['sB']:.1e} sA {r['sA']:.1e} "
              f"mu in [{r['mu_min']:.2e}, {r['mu_max']:.2e}] ratio {r['ratio']:.1f};  W(40) {r['W40']:.1e} W(60) {r['W60']:.1e} "
              f"W(70) {r['W70']:.1e} K_W(34.25) {r['KW34']:.1e}", flush=True)
