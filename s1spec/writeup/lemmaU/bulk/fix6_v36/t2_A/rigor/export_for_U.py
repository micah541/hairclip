#!/usr/bin/env python3
"""Interface file for agent B (U_ramp/U_model2.py --inputs): out of out/assembly_rep.json (tightened (H2)), certified (arb upper ends).
D_k = bound on int_{tau_rho}^{tau_1+T} int_{[28,56]} eps*omega z_{k+2}^2 ds dtau for ANY tau_rho >= tau_1 (the integral over [tau_1, tau_1+T] dominates):
   D_k <= D0_k + D1_k * T,   D0_k = (H2_k + Dsrc_k)/(r_k a_min^c),   D1_k = (2F_k^2/gamma'_k + coup_k + outer_k + rho_J Dsrc_k)/(r_k a_min^c)
(integrate dE_k + gamma'_k E_k + r_k a_min^c int_core eps*omega z_{k+2}^2 <= data sources + 2 sqrt(E_k) F_k + coup_k + outer_k from tau_1, E_k(tau_1) <= H2_k).
r = (1249/1000, 1249/1000, 999/1000000) a_min^c (level 2: only the eps4 reserve; its D is honest but useless).  [28,56] is inside the core [28,350].
"D" holds D0 (horizon-free part); "D_slope" holds D1 per unit tau (NOT zero: forcing, coupling, outer rate at the H-O placeholders)."""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, Q_, arb
R = L.load(L.OUT + '/assembly_rep.json')['res']['tight']; OR = L.load(L.OUT + '/outer_rep.json')
amin = (1 + A(Q_(11, 100))**2)**-2; r = [A(Q_(1249, 1000)), A(Q_(1249, 1000)), A(Q_(999, 1000000))]; H2 = [A(Q_(1, 10**8)), A(Q_(1, 10**6)), A(Q_(1, 10**3))]
rhoJ = R['export']['0']['rhoJ']
D0, D1, D0a = [], [], []
for k in range(3):
    inf = R['info'][k]; Ds = inf['D']; gp = R['gp'][k]; out = OR['res'][str(k)]['rate_at_placeholder']
    D0.append((H2[k] + Ds)/(r[k]*amin)); D0a.append((R['E'][k] + Ds)/(r[k]*amin))
    D1.append((2*inf['F2']/gp + inf['coup_rate'] + out + rhoJ*Ds)/(r[k]*amin))
res = dict(EA=[L.fup(x) for x in R['E']], gp=[L.flo(x) for x in R['gp']], D=[L.fup(x) for x in D0], D_slope=[L.fup(x) for x in D1],
           D_from_any_tau_a_with_E_tau_a_le_EA=[L.fup(x) for x in D0a], amin=L.flo(amin), kappa=[1.249, 1.249, 0.000999],
           cD_constant_part=[L.fup((H2[k] + R['info'][k]['D'])/R['E'][k]) for k in range(3)],
           note='D_k <= D[k] + D_slope[k]*T, T = horizon measured from tau_1 (any tau_rho >= tau_1). Horizon-free D is NOT certified: the slope '
                '(forcing 2F^2/gamma\', coupling, outer rate R5 at placeholder H-O values) is a persistent per-unit-tau source. Level 2 exports only the eps4 reserve. '
                'Labels: [cert] given v3.1 hypotheses + H-O1/H-O2 placeholders (outer rate only).')
L.dump(res, L.ROOT + '/export_for_U.json', [L.OUT + '/assembly_rep.json', L.OUT + '/outer_rep.json', os.path.abspath(__file__)])
print(json.dumps({k: v for k, v in res.items() if k != 'note'}, indent=0))
