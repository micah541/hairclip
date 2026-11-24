#!/usr/bin/env python3
"""Nested level-2 shell: ramp [45,50], core [50,60] (+ continuing into (B)); levels 0-1 as in F1 (core [28,60]).
Ramp data for level 2 (time-integrated): Jz_k = int int_{ramp} omega z_k^2 ds dtau
   k=0,1: from damping  (int gamma'_k E_k <= bud_k)      k=2: level-0 dissipation  Jz_2 <= bud_0/(2 amin eps(50))
   k=3  : level-1 dissipation  Jz_3 <= bud_1/(2 amin eps(50)).
Ramp commutators: sum over pairs (a,b<=3) P_ab sqrt(Jz_a Jz_b) (Cauchy-Schwarz), P_ab = sup eps*sum|co|Ab|rem_m|/omega.
Core rate of level 2 uses eps(45)=45^-4.  INVALID if gamma'_2<=0."""
import sys, math, json
sys.path.insert(0, '.')
import lemmaS_F1 as F, run_F1 as R
from flint import arb
gam = 1.0; y0, y1 = 45., 50.; amin_f = lambda M0: (1 + M0**2)**-2
def run(M0, M1, Mh, Tn=1e3, sqL=2.49e-7, al=5e-4, ho=(5e-7, 3e-4, 20.), J=(4.1e-11, 5.8e-11, 1.3e-11), J4=4.3e-7, tables=None, lowE=None):
    T = tables or F.build_tables()
    bx = R.box(Mh, Mh, Mh, M0=M0, M1=M1)
    c, G = F.constants(gam, bx, 'quasi', tables=T)
    Ab = F.A_bound_set('quasi', bx); amin = amin_f(M0)
    J3 = math.sqrt(J[2]*J4) + 0.8477*J[2]
    r = F.assemble(G, J + (J3,), sqL, gam, al, handoff=ho)           # levels 0,1 valid (F1 core [28,60])
    E = [float(x) for x in r['E']]; parts = r['parts']
    p = 1.0246; Rout = 350.
    cF = 1.4*(0.00244 + 63.5*al + 1.02*al*7.57)*p*math.sqrt(28.**gam*(((Rout**(2*p-gam)) - 14.**(2*p-gam))/(2*p-gam)))
    intA = 2*sqL/2.6e-3
    bud = []
    for j in (0, 1):
        pj = parts[j]
        bud.append(float(pj['handoff']) + float(pj['inflow']) + float(pj['ramp_comm']) + 2*cF*p**j*math.sqrt(E[j])*intA)
    gp = [float(parts[j]['gamma_prime']) for j in (0, 1)]
    eps50 = 50.**-4; eps45 = 45.**-4
    Jz = {0: bud[0]/gp[0], 1: bud[1]/gp[1], 2: bud[0]/(2*amin*eps50), 3: bud[1]/(2*amin*eps50)}
    # level 2 ramp sups
    Rm = F.Ramp(y0, y1, gam, True)
    n = 4; tab = T[2]; YC = G[2]['YC']; Ghat = G[2]['Ghat']
    P = {}
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        def f(t, mons=mons):
            tot = arb(0)
            for (m, i, rr), co in mons.items():
                if m == 0: continue
                tot += arb(abs(float(co))*Ab[(i, rr)])*Rm.rem(t, m)
            return tot*Rm.eps(t)/Rm.om(t)
        P[(a, b)] = F.sup_grid(f, 600)
    rho = F.sup_grid(lambda t: Rm.rho(t)/Rm.om(t), 600)
    P[(2, 2)] = P.get((2, 2), 0) + float(YC[3])*rho
    comm = sum(v*math.sqrt(Jz[a]*Jz[b]) for (a, b), v in P.items())
    # inflow
    cbar = 0.25 + 7.57*sqL + 7.57*sqL/y0
    infl_sup = F.sup_grid(lambda t: abs(Rm.e2(t, 1)), 600)
    infl = cbar*infl_sup*Jz[2]
    # damping & forcing
    cmin = 0.25 - 7.57*sqL - 7.57*sqL/y0
    Gjj = float(Ghat[2])
    gp2 = gam*cmin - 3*7.57*sqL/y0 - eps45*Gjj
    e = 2*p - gam; Iint = y0**gam*((Rout**e - y0**e)/e)
    F2 = (cF*p*p*sqL/math.sqrt(28.**gam*(((Rout**e) - 14.**e)/e))*0 + 0)  # placeholder
    r0 = 1.4*((0.00244 + 63.5*al)*sqL + 1.02*al*7.57*sqL)
    F2 = (r0*p**2)**2*Iint
    coup = sum(float(Ghat.get(k, 0))*eps45*E[k] for k in (0, 1))
    E2 = ho[2] + infl + comm + (4*F2/gp2**2 + 2*coup/gp2 if gp2 > 0 else float('nan'))
    # sup z_1 on [45,60]
    ell = math.log(60/45.); conv = (60/28.)**gam
    Eu1 = conv*E[1]; Eu2 = conv*E2 if gp2 > 0 else float('nan')
    S1 = math.sqrt(Eu1/ell + 2*math.sqrt(Eu1*(Eu1 + Eu2)))
    out = dict(M=(M0, M1, Mh), G22=Gjj, gp2=gp2, valid=gp2 > 0, E2=E2, infl=infl, comm=comm, P={str(k): v for k, v in P.items()}, Jz=Jz, S1_on_45_60=S1, E01=E[:2], bud=bud)
    return out
if __name__ == '__main__':
    T = F.build_tables(); res = {}
    for nm, (M0, M1, Mh) in {'honest': (0.07, 1.4, 3), 'GNideal_M1_1.4': (0.07, 1.4, 0.0)}.items():
        o = run(M0, M1, Mh, tables=T); res[nm] = o
        print(nm, 'G22 %.3g gp2 %.4f' % (o['G22'], o['gp2']), '' if o['valid'] else 'INVALID', 'E2 %.4g (infl %.3g comm %.3g)  S1[45,60] %.3f' % (o['E2'], o['infl'], o['comm'], o['S1_on_45_60']))
        print('   P', {k: '%.3g' % v for k, v in o['P'].items()}, ' Jz', {k: '%.3g' % v for k, v in o['Jz'].items()}, flush=True)
    json.dump(res, open('lemmaS_nested.json', 'w'), indent=1, default=str)
