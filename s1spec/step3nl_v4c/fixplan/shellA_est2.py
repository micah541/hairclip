#!/usr/bin/env python3
"""FLOAT estimate (not rigorous) of Lemma S(A), version 2: deviation shells v = w', levels 0-2, core y in [28,60], ramp y in [YA,28].
Difference to shellA_est.py: ramp commutators use ORDER-j ramp data only.  For Theta = omega eta^2 e^{-4s} and w = d_s^j v:
   -int Theta w w'''' = -int Theta w''^2 - 2 int Theta'' w w'' - (1/2) int Theta'''' w^2          (exact; no Theta' w' w'' term)
   |2 int Theta'' w w''| <= delta int Theta w''^2 + (1/delta) int (Theta''^2/Theta) w^2           (Theta''^2/Theta is bounded for eta = smoothstep)
so every ramp term is  (bounded coefficient) x int_ramp (d_s^j v)^2, bounded by sup(coef y^(2m-1)/W) J_m, m <= j.  No J3, J4.
Nonlinear terms on the ramp are treated as data too: eps*G_j(M1_ramp)*Theta*(d_s^j v)^2 with G_j supplied (the 'eps-split').
usage: shellA_est2.py J0 J1 J2 [YA] [G0 G1 G2] [E0h E1h E2h]"""
import sys, numpy as np
X0 = 2.2490776515
def run(J0, J1, J2, YA=14.0, YC=28.0, YO=60.0, gam=1.0, cbar=0.2501, delta=0.25, G=(4.2e3, 6.5e4, 2.5e6), Eh=(0, 0, 0), slop=2.0, verbose=True):
    J = [J0, J1, J2]
    L = np.log(YC/YA); y = np.linspace(YA, YC, 400001)[1:-1]; t = np.log(y/YA)/L
    S = t**3*(10 - 15*t + 6*t*t); d1 = 30*t*t*(1 - t)**2/L; d2 = 60*t*(1 - t)*(1 - 2*t)/L**2
    d3 = 60*(1 - 6*t + 6*t*t)/L**3; d4 = 60*(-6 + 12*t)/L**4
    om = (y/YC)**-gam; W = (1 + (y + X0)**2)**-2; eps = y**-4
    # F = eta^2 and its s-derivatives
    F0 = S*S; F1 = 2*S*d1; F2 = 2*d1*d1 + 2*S*d2; F3 = 6*d1*d2 + 2*S*d3; F4 = 6*d2*d2 + 8*d1*d3 + 2*S*d4
    k = gam + 4.0                                            # Theta = F e^{-k s} (up to a constant):  Theta^(n) e^{ks} = (D-k)^n F
    T0 = F0; T1 = F1 - k*F0; T2 = F2 - 2*k*F1 + k*k*F0
    T4 = F4 - 4*k*F3 + 6*k*k*F2 - 4*k**3*F1 + k**4*F0
    th_s = om*(F1 - gam*F0)
    ST = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}}
    def ramp_int(coef, j): n = len(ST[j]); return n*sum(c*c*np.max(coef*y**(2*m - 1)/W)*J[m] for m, c in ST[j].items())
    out = []
    for j in range(3):
        infl = cbar*ramp_int(np.maximum(th_s, 0), j)
        # 4th-order part (a <= 1) + the lower-order parts of P4 = D^4 - 6D^3 + 11D^2 - 6D, all reduced to order j with one Young each
        c4 = om*eps*((T2*T2/np.maximum(T0, 1e-300))/delta + 0.5*np.abs(T4))
        c3 = om*eps*6*(1.5*(T1*T1/np.maximum(T0, 1e-300))/delta + 0.5*np.abs(F3 - 3*k*F2 + 3*k*k*F1 - k**3*F0))
        c2 = om*eps*11*(np.abs(T2))
        comm = slop*ramp_int(c4 + c3 + c2, j)
        nl = ramp_int(om*eps*G[j]*F0, j)
        out.append((infl, comm, nl))
    conv = (YO/YC)**gam; ell = np.log(YO/YC)
    Eu = [conv*(sum(o) + h) for o, h in zip(out, Eh)]
    S0v = np.sqrt(Eu[0]/ell + 2*np.sqrt(Eu[0]*Eu[1])); S1v = np.sqrt(Eu[1]/ell + 2*np.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    if verbose:
        for j in range(3): print(f"  level {j}: inflow {out[j][0]:.3e}  commutators(order j, x{slop:g}) {out[j][1]:.3e}  ramp nonlinear (G={G[j]:.1e}) {out[j][2]:.3e}  hand-off {Eh[j]:.1e} -> sup E_{j} <= {Eu[j]:.3e}")
        print(f"  S0_v <= {S0v:.4f}   S1_v <= {S1v:.3f}   sup|w''| on core <= {S1v/YC:.4f}")
    return S0v, S1v, Eu
if __name__ == "__main__":
    a = [float(x) for x in sys.argv[1:]]
    kw = {}
    if len(a) > 3: kw['YA'] = a[3]
    if len(a) > 6: kw['G'] = tuple(a[4:7])
    if len(a) > 9: kw['Eh'] = tuple(a[7:10])
    run(*a[:3], **kw)
