"""FEASIBILITY ONLY (float/arb-sup mix, not a certificate): constants for a far-field energy of the FULL slope f with a co-moving power weight
Theta_2 = (r/Y1)^mu (i.e. gamma = -mu in the tables of lemmaS_F1), inner cutoff ramp fixed in rescaled coordinates on [28, 56].
Prints, for mu in (3, 4): core rate constants Ghat_jj (box M0 = 0.11, M1 = 0.52, no L^4 extras), eps(28) Ghat_22, and the ramp commutator
sups divided by eps*Theta at the ramp (dimensionless K_jk).   run from writeup/lemmaS_v32:  python3 -B ../step3_gap_v32/feas_P.py"""
import sys, json, math
sys.path.insert(0, '.')
from flint import arb
import lemmaS_F1 as F, run_F1 as R
T = F.build_tables(); out = {}
for mu in (3.0, 4.0):
    gam = -mu; box = R.box(0, 0, 0, Umax=1e-12, M0=0.11, M1=0.52)   # U = 0: energies of the full slope f
    o, G = F.constants(gam, box, 'quasi', y0=28.0, y1=56.0, tables=T)
    Gh = [G[j]['Ghat'][j] for j in range(3)]
    Rm = F.Ramp(28.0, 56.0, gam, increasing=True)
    Ab = F.A_bound_set('quasi', box); K = {}
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
                if k == j: tot += arb(YC[j + 1])*Rm.rho(t)/Rm.eps(t)
                return tot/Rm.om(t)                       # divided by eps*Theta: dimensionless
            row[k] = F.sup_grid(Yrem, n=300)
        K[j] = row
    out[mu] = dict(Ghat=Gh, epsG=[g*28.0**-4 for g in Gh], K=K)
    print('mu', mu, 'Ghat_jj', ['%.3g' % g for g in Gh], 'eps(28)*Ghat', ['%.3g' % (g*28.0**-4) for g in Gh], 'K_2k', {k: '%.3g' % v for k, v in K[2].items()}, 'K_1k', {k: '%.3g' % v for k, v in K[1].items()}, flush=True)
json.dump({str(k): v for k, v in out.items()}, open('../step3_gap_v32/feas_P.json', 'w'), indent=1, default=str)
