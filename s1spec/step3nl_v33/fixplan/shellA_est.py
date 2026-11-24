#!/usr/bin/env python3
"""FLOAT estimate (not rigorous) of the redesigned Lemma S(A): deviation shells v = w' on y in [28,60], inner ramp y in [YA,28]
(x = y + 2.249 in [16.25,30.25]: zeta = 1 and W = (1+x^2)^-2 exactly), eta = quintic smoothstep in t = ln(y/YA)/ln(28/YA).
Energies E_j = int omega eta^2 (d_s^j v)^2 ds, s = ln y, omega = (y/28)^-gam.  Ramp terms by Cauchy-Schwarz against the
time integrals J_k = int V_k^2 dtau (k <= j+1), NO Young against the shell dissipation.
usage: shellA_est.py J0 J1 J2 J4 [YA]"""
import sys, numpy as np
X0 = 2.2490776515
def run(J0, J1, J2, J4, YA=14.0, YC=28.0, YO=60.0, gam=1.0, cbar=0.2501, C1W=0.8477, verbose=True):
    J3 = np.sqrt(J2*J4) + C1W*J2
    J = [J0, J1, J2, J3]
    L = np.log(YC/YA); y = np.linspace(YA, YC, 200001); t = np.log(y/YA)/L
    S = t**3*(10 - 15*t + 6*t*t); S1 = 30*t*t*(1 - t)**2/L; S2 = 60*t*(1 - t)*(1 - 2*t)/L**2
    om = (y/YC)**-gam; W = (1 + (y + X0)**2)**-2
    th = om*S*S                                   # weight on the ramp
    th_s = om*(2*S*S1 - gam*S*S)                  # d/ds
    # d^2/ds^2 of Theta = th*e^{-4s} (divided by e^{-4s}); crude sup of |Theta''|/e^{-4s}
    e2 = np.abs(om*(2*S1*S1 + 2*S*S2) - 2*(gam + 4)*om*2*S*S1 + (gam + 4)**2*th)
    e1 = np.abs(th_s - 4*th)                      # |Theta'|/e^{-4s}
    eps = y**-4
    # d_s^j v in terms of y-derivatives: j=0: v; 1: y v'; 2: y^2 v'' + y v'; 3: y^3 v''' + 3 y^2 v'' + y v'
    ST = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
    def ramp_int(coef, k):                        # int int coef (d_s^k v)^2 ds dtau  <=  n_k sum_m c_m^2 sup(coef y^(2m-1)/W) J_m
        n = len(ST[k]); return n*sum(c*c*np.max(coef*y**(2*m - 1)/W)*J[m] for m, c in ST[k].items())
    E = []
    for j in range(3):
        infl = cbar*ramp_int(np.maximum(th_s, 0), j)
        # linear model  -eps a P4(d_s):  2 int Theta'' w'^2 + (1/2) int |Theta''''| w^2 + 9 int |Theta'| w'^2 + ... ; keep the two leading ramp terms, x3 slop
        comm = 3.0*(2*ramp_int(eps*e2, j + 1) + 9*ramp_int(eps*e1, j + 1) + 50*ramp_int(eps*e2, j))
        E.append((infl, comm))
    conv = (YO/YC)**gam; ell = np.log(YO/YC)
    Eu = [conv*(a + b) for a, b in E]
    S0v = np.sqrt(Eu[0]/ell + 2*np.sqrt(Eu[0]*Eu[1])); S1v = np.sqrt(Eu[1]/ell + 2*np.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    if verbose:
        print(f"ramp [{YA},{YC}] gam={gam}: J3={J3:.2e}")
        for j in range(3): print(f"  level {j}: inflow {E[j][0]:.3e}  ramp commutators {E[j][1]:.3e}  -> sup E_{j} (unweighted core) <= {Eu[j]:.3e}  (+ hand-off E_{j}(tau1) x {conv:.2f})")
        print(f"  S0_v = sup|w'| on y in [28,60] <= {S0v:.4f};  S1_v = sup y|w''| <= {S1v:.3f}  (sup|w''| <= {S1v/YC:.4f})")
    return S0v, S1v, Eu
if __name__ == "__main__":
    a = [float(x) for x in sys.argv[1:5]]; YA = float(sys.argv[5]) if len(sys.argv) > 5 else 14.0
    run(*a, YA=YA)
