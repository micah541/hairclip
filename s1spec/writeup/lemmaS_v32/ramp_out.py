"""SCRATCH (step3_amend_v32, Risk 2): commutator coefficients of an OUTER cutoff ramp for the (B) energies, levels j = 0, 1, 2.
Same divergence-form tables and remainder formulas as lemmaS_F1.constants (inner ramp), with the decreasing cutoff
eta = 1 - S(t), y = y0 (y1/y0)^t, omega = (y/28)^-1, box (M0, M1) = (U(y1) + 0.15, 1.05 U(y1) + 0.15) as in lemmaS_Bbox (quasi mode, no L^4 extras).
Output c[j][k] >= sup_t of the coefficient of int_ramp z_k^2 ds in dE_j/dtau (k <= j+1), at entry coordinates; without the
'ramp share of the core rate' (that term multiplies eta^2 omega z_k^2, i.e. the energy itself, and is in the Gronwall factor).
At later tau the coefficient is multiplied by (lambda/lambda_1)^4 (eps = y^-4; frozen weight), so int dtau gives the factor 1/(1 - 4 abar).
usage: python3 -B ramp_out.py y0 y1"""
import sys, json, math
from flint import arb
import lemmaS_F1 as F, run_F1 as R, lemmaS_prop1_bounds as PB
y0, y1 = float(sys.argv[1]), float(sys.argv[2]); gam = 1.0; al = 1.5e-4; p = 1.0246
U = lambda y: 1.3056*al*y**p
M0 = U(y1) + 0.15; M1 = 1.05*U(y1) + 0.15
box = R.box(0, 0, 0, M0=M0, M1=M1); Ab = F.A_bound_set('quasi', box); amin = (1 + M0**2)**-2; kap = gam + 4
Rm = F.Ramp(y0, y1, gam, increasing=False); T = F.build_tables(); out = {}
for j in (0, 1, 2):
    n = j + 2; tab = T[j]; YC = {k: 0.0 for k in range(n + 1)}; pairs = []
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        gC = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
        pairs.append((a, b, mons))
    row = {}
    for k in range(j + 2):
        def Yrem(t, k=k):
            tot = arb(0)
            for (a, b, mons) in pairs:
                w = 0.0
                if a == k: w += 1.0 if a == b else 0.5
                elif b == k: w += 0.5
                if w == 0.0: continue
                for (m, i, r), co in mons.items():
                    if m == 0: continue
                    tot += arb(w*co*Ab[(i, r)])*Rm.rem(t, m)
            tot = tot*Rm.eps(t)
            if k == j: tot += arb(YC[j + 1])*Rm.rho(t)
            return tot
        row[str(k)] = F.sup_grid(Yrem)
    out[str(j)] = row; print(j, {k: '%.3e' % v for k, v in row.items()}, flush=True)
json.dump(dict(y0=y0, y1=y1, box=[M0, M1], c=out, eps_edge=y0**-4), open(f'ramp_out_{y0:g}_{y1:g}.json', 'w'), indent=1)
