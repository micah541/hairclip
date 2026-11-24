#!/usr/bin/env python3
"""(B) windows with the bent-cone profile: Ghat_jj at box (M0,M1) (|f|<=M0, y|f'|<=M1, full f=U+v), gamma'_j at y_E, deviation energies propagated, derived sups and Lemma O jets K2,K3,K4 (U-part + deviation).
usage: lemmaS_Bbox.py e [dev_allow]  (per-window boxes M0_w = U(1.5y_w)+dev, M1_w = 1.05U(1.5y_w)+dev)   (hand-off deviation energies E1=E2=E3=E4=e, E0=min(e,2e-7))"""
import sys, math, json, pickle
sys.path.insert(0, '.')
import lemmaS_F1 as F, lemmaS_L4terms as L4, run_F1 as R, lemmaS_prop1_bounds as PB, lemmaS_B3 as B
e0 = float(sys.argv[1]); dv0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
import os
gam = 1.0; p = 1.0246; al = float(os.environ.get('AL', 5e-4)); sq = float(os.environ.get('SQ', 8.19e-8))
ABAR = os.environ.get('ABAR'); BETA = os.environ.get('BETA')
U = lambda y: 1.3056*al*y**p
def sup(a_, b_, ell): return math.sqrt(a_/ell + 2*math.sqrt(a_*(a_ + b_)))
T = F.build_tables(); f2 = 0.5; L = 50.; res = []
Q = pickle.load(open('prop1_forms3.pkl', 'rb'))[3][0]; tab3 = PB.table(B.div_form3(Q, 5))
for yE in (50., 75., 112.5, 168.75, 253.125):
    yb = 1.5*yE; M0 = U(yb) + dv0; M1 = 1.05*U(yb) + dv0; amin = (1 + M0**2)**-2
    bc = R.box(0, 0, 0, M0=M0, M1=M1); cc, G = F.constants(gam, bc, 'quasi', tables=T); Gk = {}
    for j in range(3):
        add, log, dis, spare = L4.l4_extra_YC(T, j, M0, M1, gam, s33='ibp'); assert dis <= spare*(1 + 1e-9)
        YC = {int(k): float(v) for k, v in G[j]['YC'].items()}; YC = {k: YC[k] + add.get(k, 0) for k in YC}; Gk[j] = L4.ghat_from_YC(YC, j, amin)
    Ab = F.A_bound_set('quasi', bc); YC = {k: 0.0 for k in range(6)}
    for (a, b), mons in tab3.items():
        if (a, b) == (5, 5): continue
        g = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
        if a == b: YC[a] += g
        else: YC[a] += g/2; YC[b] += g/2
    add, log, dis, spare = L4.l4_extra_YC({3: tab3}, 3, M0, M1, gam, fallback={2: f2, 3: 0.0}, s33='ibp')
    YC3 = {k: YC[k] + add.get(k, 0) for k in YC}; Gk[3] = L4.ghat_from_YC(YC3, 3, amin)
    abar = 7.57*sq if ABAR is None else float(ABAR); cmin = 0.25 - abar - abar/50; eps = yE**-4
    gp = [gam*cmin - (abs(1 - 2*j) + 1)*abar/yE - eps*Gk[j][j] for j in range(4)]
    Rout = yb; ee = 2*p - gam; IB = yE**gam*((Rout**ee - yE**ee)/ee); r0 = 1.4*(((0.00244 + 63.5*al)*sq if BETA is None else float(BETA)) + 1.02*al*abar)
    Es = [min(e0, 2e-7), e0, e0, e0]; Cexp = max(Gk[j][j] for j in range(4))*eps; EB = []
    for j in range(4):
        coup = sum(Gk[j].get(k, 0)*eps*EB[k] for k in range(j)) if j < 3 else sum(YC3[k]*eps*EB[k] for k in range(3))
        EB.append(Es[j]*math.exp(Cexp)*math.exp(8*abar/yE) + 4*(r0*p**j)**2*IB/gp[j]**2 + 2*coup/gp[j])
    cv = yb/28.; ell = math.log(1.5); Eu = [cv*x for x in EB]
    S = [sup(Eu[0], Eu[1], ell), sup(Eu[1], Eu[2], ell), sup(Eu[2], Eu[3], ell), sup(Eu[3], Eu[3], ell)]
    der = (U(yb) + S[0], 1.05*U(yb) + S[1], p*p*U(yb) + S[2])
    y = yE
    K2 = (L/y)*(1.05*U(y) + S[1]); K3 = (L/y)**2*(S[2] + S[1] + U(y)*abs(p*p - p)); K4 = (L/y)**3*(S[3] + 3*S[2] + 2*S[1] + U(y)*abs(p*(p - 1)*(p - 2)))
    K2 = max(K2, (L/yb)*(1.05*U(yb) + S[1]))
    res.append(dict(yE=yE, box=[M0, M1], Ghat=[float(Gk[j][j]) for j in range(4)], gp=[float(x) for x in gp], EB=[float(x) for x in EB], S=[float(x) for x in S], derived=[float(x) for x in der], K=[float(K2), float(K3), float(K4)]))
    print('window [%.0f,%.0f] box (%.3f,%.3f) Ghat %s gamma %s %s' % (yE, yb, M0, M1, [round(float(Gk[j][j])) for j in range(4)], [round(float(x), 3) for x in gp], 'VALID' if min(gp) > 0 else 'INVALID'))
    print('   E_B %s  derived |f|<=%.3f y|f\'|<=%.3f |z2|<=%.3f ; K2 %.3f K3 %.3f K4 %.3f' % (['%.2g' % x for x in EB], *der, K2, K3, K4), flush=True)
mk = [max(r['K'][i] for r in res) for i in range(3)]
print('max K2 %.3f (x%.2f vs M=0.12), K3 %.3f (x%.2f), K4 %.3f (x%.2f)' % (mk[0], .12/mk[0], mk[1], .12/mk[1], mk[2], .12/mk[2]))
json.dump(res, open(f'lemmaS_Bbox_{e0:g}' + ('' if al == 5e-4 else f'_a{al:g}_dv{dv0:g}') + '.json', 'w'), indent=1)
