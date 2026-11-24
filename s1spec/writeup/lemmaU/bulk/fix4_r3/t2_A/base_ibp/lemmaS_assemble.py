#!/usr/bin/env python3
"""Assemble Lemma S (A): energies, S0_v, S1_v for given inputs.  Uses prop1_bounds.constants (machine-verified forms + arb sups)."""
import math, json, sys
sys.path.insert(0, '.')
from flint import arb, ctx
ctx.prec = 100
from lemmaS_interval import sup_arb
import lemmaS_prop1_bounds as PB
X0 = 2.2490776515
Sq = lambda t: t**3*(10 - 15*t + 6*t*t); Sp = lambda t: 30*t*t*(1 - t)**2
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}}
BELL = {0: 1, 1: 1, 2: 2}
def Wf(y): return (1 + (y + arb(str(X0)))**2)**-2
def ramp_K(La):
    DS = math.log(24.0/La)
    def theta_y(y):       # d vartheta/dy, vartheta = S(t)^4, t = ln(y/La)/DS
        t = (y/La).log()/DS
        return 4*Sq(t)**3*Sp(t)/(DS*y)
    Kth = [float(sup_arb(lambda y, m=m: theta_y(y)*y**(2*m)/Wf(y), La, 24).upper()) for m in range(3)]
    Kc = [float(sup_arb(lambda y, m=m: y**(2*m - 1)/Wf(y), La, 24).upper()) for m in range(3)]
    return Kth, Kc
def run(J, sqLam, gam, La=19.0, alpha1=1e-3, Rout=250.0, E1init=(0, 0, 0), verbose=False, with_ramp=True, cache={}):
    key = (gam, La)
    if key not in cache:
        out, G = PB.constants(gam=gam, La=La); cache[key] = (out, G, ramp_K(La))
    out, G, (Kth, Kc) = cache[key]
    Abar = sqLam; abar = 7.57*Abar; bbar = abar
    cbar = 0.25 + abar + bbar/La; cmin = 0.25 - abar - bbar/La
    eps_in = 24.0**-4; eps_a = La**-4; w_a = (24.0/La)**gam
    E = []
    p = 1.0246
    for j in range(3):
        GC, GR = G[j]
        gp = gam*cmin - (abs(1 - 2*j) + 1)*bbar/La - eps_in*GC[j]
        infl = cbar*w_a*BELL[j]*sum(STIR[j][m]*Kth[m]*J[m] for m in STIR[j])
        Rk = lambda k: BELL[k]*sum(STIR[k][m]*Kc[m]*J[m] for m in STIR[k])
        rampc = eps_a*w_a*sum(GR[k]*Rk(k) for k in range(j + 1)) if with_ramp else 0.0
        coup = sum(GC[k]*eps_in*E[k] for k in range(j))
        beta = (0.00244 + 63.5*alpha1)*Abar
        r0 = 1.4*(beta + 1.02*alpha1*abar)
        # F_j^2 = int omega (p^j r0 y^p)^2 ds,  omega <= (y/24)^-gam, y in [La, Rout]
        e = 2*p - gam
        Iint = 24.0**gam*((Rout**e - La**e)/e if abs(e) > 1e-9 else math.log(Rout/La))
        F2 = (r0*p**j)**2*Iint
        Ej = E1init[j] + infl + rampc + 4*F2/gp**2 + 2*coup/gp
        E.append(Ej)
        if verbose: print(j, 'gamma_prime', gp, 'inflow', infl, 'ramp_comm', rampc, 'forcing', 4*F2/gp**2, 'coupling', 2*coup/gp)
    conv = 2.5**gam; ell = math.log(2.5)
    Eu = [conv*e_ for e_ in E]
    S0 = math.sqrt(Eu[0]/ell + 2*math.sqrt(Eu[0]*Eu[1]))
    S1 = math.sqrt(Eu[1]/ell + 2*math.sqrt(Eu[1]*Eu[2]))
    return dict(E=E, S0=S0, S1=S1, S1_over_y24=S1/24)
if __name__ == '__main__':
    table = {"2.5e-7": (4.19e-11, 5.89e-11, 1.30e-11), "3.5e-7": (7.21e-11, 1.01e-10, 2.24e-11),
             "5e-7": (1.32e-10, 1.86e-10, 4.11e-11), "1e-6": (4.69e-10, 6.58e-10, 1.45e-10)}
    res = {}
    for La in (19.0, 18.0):
        for gam in (1.0, 2.0):
            for lv, J in table.items():
                sq = float(lv); Apk = 6.39e-7*sq/3.5e-7
                r1 = run(J, Apk, gam, La, with_ramp=False); r2 = run(J, Apk, gam, La, with_ramp=True)
                res[f'La{La}_gam{gam}_{lv}'] = dict(no_ramp_commutators=r1, with_ramp_commutators=r2)
                print(f'La={La} gam={gam} sqrtLam1={lv}: NO-COMM E=({r1["E"][0]:.1e},{r1["E"][1]:.1e},{r1["E"][2]:.1e}) S0={r1["S0"]:.4f} S1={r1["S1"]:.3f} | '
                      f'WITH-COMM E=({r2["E"][0]:.1e},{r2["E"][1]:.1e},{r2["E"][2]:.1e}) S0={r2["S0"]:.3f} S1={r2["S1"]:.2f}', flush=True)
    json.dump(res, open('lemmaS_assembled.json', 'w'), indent=1)
