"""J-cascade (LEMMA_F.md, Lemma J): time-weighted ramp data of the E_I norms from the tip closure, per unit Ybar.
For an end time tau and a weight omega(s) = exp(-lam (tau - s)), J_m := int_{tau0}^tau omega V_m^2 ds, Ybar := sup Ybar = 1.
Inputs are the T5/T7 constants (here: S1 reference values / T5 defaults, flagged in T8_REPORT); all enter monotonically.

(J-Lambda) from 1/2 Lambda' <= -rho Lambda - kD Dsurplus + a_c Y + Xi'^2 and Ybar' = -rho Ybar + Phi:
   int omega Dsurplus <= B := (1/2 + (rho - lam/2)^+ / lam)/kD,   Dsurplus = th0/4 B^2 + eps/8 D4^2 + muE/8 V2^2 + muE epsI/8 V4^2.
(I0) 1/2 (V0^2)' <= -mu1 V0^2 - 1/2 V2^2 + af W1^2 + bf B^2 + |x| gam0 V0   (|x| <= Cx A).
(I2) 1/2 (V2^2)' <= -1/2 V4^2 - mu2 V2^2 + X (V4 + 2 Psi1 V3 + 2 c1W V2),  X = c3 D4 + c2 W3 + c1 B + c0 W1 + |x| gam2.
Interpolation: W1^2 <= A B + c1 A^2, W3^2 <= B D4 + c1 B^2, V1^2 <= V0 V2 + c1W V0^2, V3^2 <= V2 V4 + c1W V2^2."""
import numpy as np

T7 = dict(muE=0.02, eps=1e-5, epsI=3e-5, th0=0.02, rho=0.03, c1=0.5804, af=0.0696, bf=0.0394, mu1=0.1552, mu2=0.6555,
          cI=(0.425, 0.424, 0.526, 0.442), c1W=0.8477, Psi1=0.39, Cx=4.0, gam0=0.05, gam2=0.01, Ke=1.4e4, kD=0.5)

def cascade(lam, dt_from_start=np.inf, P=None, D4_pointwise=None, verbose=False):
    """returns dict of J_m (m=0..4) and the tip integrals, per unit Ybar^2.
    dt_from_start: tau - tau0 (initial terms decay like exp(-lam dt)).
    D4_pointwise: if given (fallback 1), D4 <= D4_pointwise * Ybar pointwise in time (time-differentiated energies)."""
    p = dict(T7); p.update(P or {})
    w0 = np.exp(-lam * dt_from_start) if np.isfinite(dt_from_start) else 0.0
    B = (0.5 + max(p['rho'] - lam / 2, 0.0) / lam) / p['kD']   # reserved surplus fraction kD (Xi^2 -> Xi^2/(1-kD) in T5 (Y))
    JA = 1.0 / lam                                       # A <= sqrt(Lambda) <= 1
    JB = min(4 * B / p['th0'], 1.0 / (lam * p['eps']))
    JD4 = 8 * B / p['eps']
    if D4_pointwise is not None: JD4 = min(JD4, D4_pointwise**2 / lam)
    JW1 = np.sqrt(JA * JB) + p['c1'] * JA
    JW3 = np.sqrt(JB * JD4) + p['c1'] * JB
    # V0 (I0 damping), V2 (I0 dissipation), with the Lambda-routes as alternatives
    K = 0.5 * w0 / p['muE'] + p['af'] * JW1 + p['bf'] * JB
    mg = p['Cx'] * p['gam0'] * np.sqrt(JA)
    a = p['mu1'] - lam / 2
    JV0_L = 1.0 / (lam * p['muE'])
    if a > 0:
        s = (mg + np.sqrt(mg**2 + 4 * a * K)) / (2 * a); JV0 = min(JV0_L, s**2)
    else:
        JV0 = JV0_L
    JV2_I0 = 2 * K + 2 * mg * np.sqrt(JV0) + 2 * max(lam / 2 - p['mu1'], 0.0) * JV0
    JV2 = min(JV2_I0, 8 * B / p['muE'])
    # V4 via (I2)
    sqX = p['cI'][3] * np.sqrt(JD4) + p['cI'][2] * np.sqrt(JW3) + p['cI'][1] * np.sqrt(JB) + p['cI'][0] * np.sqrt(JW1) \
        + p['Cx'] * p['gam2'] * np.sqrt(JA)
    best = np.inf
    for kap in np.exp(np.linspace(-3, 3, 121)):
        al = 1 + p['Psi1'] / kap
        bh = p['Psi1'] * kap + 2 * p['Psi1'] * np.sqrt(p['c1W']) + 2 * p['c1W']
        den = p['mu2'] - lam / 2
        if den <= 0: continue
        best = min(best, 4 * (al**2 + bh**2 / (4 * den)))
    JV4_I2 = 2 * w0 / (p['muE'] * p['epsI']) + best * sqX**2
    JV4 = min(JV4_I2, 8 * B / (p['muE'] * p['epsI']))
    J = {0: JV0, 2: JV2, 4: JV4}
    J[1] = np.sqrt(J[0] * J[2]) + p['c1W'] * J[0]
    J[3] = np.sqrt(J[2] * J[4]) + p['c1W'] * J[2]
    out = dict(J=J, JA=JA, JB=JB, JD4=JD4, JW1=JW1, JW3=JW3, B=B, route=dict(V2='I0' if JV2_I0 < 8 * B / p['muE'] else 'Lambda',
               V4='I2' if JV4_I2 < 8 * B / (p['muE'] * p['epsI']) else 'Lambda'),
               JV4_shares=dict(D4=(p['cI'][3] * np.sqrt(JD4) / sqX)**2 * best * sqX**2 if sqX > 0 else 0, total=JV4))
    if verbose: print(out)
    return out

if __name__ == '__main__':
    for lam in (0.22, 0.25):
        o = cascade(lam)
        print(f"lam={lam}: J0..J4 = " + ', '.join(f'{o["J"][m]:.3e}' for m in range(5)), ' routes', o['route'],
              f" JB={o['JB']:.3g} JD4={o['JD4']:.3g} JW1={o['JW1']:.3g} JW3={o['JW3']:.3g} B={o['B']:.3f}")
        o = cascade(lam, D4_pointwise=200.0)
        print(f"   fallback-1 (D4 <= 200 Ybar pointwise): J0..J4 = " + ', '.join(f'{o["J"][m]:.3e}' for m in range(5)))
