#!/usr/bin/env python3
"""F4 absolute CFLUX constants + (B) entry/level-2 + S0_f,S1_f + H6 K2,K3, all from a lemmaS_<tag>.json (F2 run).  usage: lemmaS_F4B.py TAG"""
import json, math, sys
import os
tag = sys.argv[1] if len(sys.argv) > 1 else 'F2'
RT = os.environ.get('RT', '../../step3nl/verifier/rate_table_unc9_m25_L5e-8.json'); ABAR = os.environ.get('ABAR'); BETA = os.environ.get('BETA')
d = json.load(open(f'lemmaS_{tag}.json')); c = json.load(open('lemmaS_cflux.json'))
a = d['args']; gam = a['gam']; sq = a['sqLam']; al = a['alpha1']
E = [float(x) for x in d['result']['E']]; parts = d['result']['parts']; Gc = [float(d['Ghat_core'][f'Ghat_{j}'][str(j)]) for j in range(3)]
Gk = {j: {int(k): float(v) for k, v in d['Ghat_core'][f'Ghat_{j}'].items()} for j in range(3)}
M0c = a['core'][0]; amin = (1 + M0c**2)**-2; p = 1.0246; Rout = 350.
# ---- F4 ----
W = c['W57_upper']; k0 = c['kappa_I0']; k2 = c['kappa_I2']; conv = (60/28.)**gam
Eu = [conv*x for x in E]; rho = a.get('rho', 2.6e-3); intA = 2*sq/rho
e_ = 2*p - gam; Iint = 28.**gam*((Rout**e_ - 14.**e_)/e_)
cF = 1.4*(((0.00244 + 63.5*al)*sq if BETA is None else float(BETA)) + 1.02*al*(7.57*sq if ABAR is None else float(ABAR)))/sq*p*math.sqrt(Iint)
out = {}
# uniform-in-T (affine) form: no horizon.  Phi(T) = int_{tau1}^T sqrt(Lambda) <= PHI_C + PHI_S (T-tau1) from the verifier tail lemma:
#   Lambda <= L* e^{-rho(tau-tau_t)} + L_F  =>  sqrt <= sqrt(L*) e^{-rho(.)/2} + sqrt(L_F);  steps k<K: sum sqrt(Lup_k) D_k (rate_table_fix5_L1e-7)
tb = json.load(open(RT)); tl = tb['tail']
PHI_C = sum(max(math.sqrt(x['Ly_hi']), x['sqrtL_next'])*x['D'] for x in tb['steps']) + 2*math.sqrt(tl['L_star'])/tl['rho']; PHI_S = tl['sqrt_LF']
pr = [(float(pa['inflow']) + float(pa['ramp_comm']), float(pa['gamma_prime']), float(pa['forcing']) + float(pa['coupling']), float(pa['handoff'])) for pa in parts]
Jv = [float(x) for x in a['J'][:3]]; Sv = [float(x) for x in a['slopes'][:3]]
Es = []; Ec = []; ssl = []
for j in range(3):
    ir, g, fc, ho_ = pr[j]; r = max(Sv[m]/(g*Jv[m]) for m in range(3)); fr = r/(1 + r)
    Es.append(fr*ir + fc); Ec.append(E[j] - Es[-1]); ssl.append(g*fr*ir)      # Es: slope-driven floor part of the pointwise bound; ssl: rate of growth of the impulse total
IE_c = [Ec[j]/pr[j][1] for j in range(3)]; IE_s = Es                          # int E_j <= IE_c + IE_s*(T-tau1)
e_ = 2*p - gam; Iint = 28.**gam*((Rout**e_ - 14.**e_)/e_)
cF = 1.4*(((0.00244 + 63.5*al)*sq if BETA is None else float(BETA)) + 1.02*al*(7.57*sq if ABAR is None else float(ABAR)))/sq*p*math.sqrt(Iint)
pa = parts[1]; bud_c = float(pa['handoff']) + float(pa['inflow']) + float(pa['ramp_comm']) + 2*cF*math.sqrt(E[1])*PHI_C
bud_s = ssl[1] + 2*cF*math.sqrt(E[1])*PHI_S          # dropped: -gamma'_1 int E1 <= 0 (the eps*Ghat_11 term is absorbed in gamma'_1)
yw = 57.75; fac = 1/((amin*yw**-4*(28./yw))*2)*(yw/28.)**gam
Eu3c, Eu3s = bud_c*fac, bud_s*fac
N0 = W*57.75*Eu[0]; N1 = W*Eu[1]/54.75; N2 = W*2*(Eu[2] + Eu[1])/54.75**3
N3c = W*3*(Eu3c + conv*(9*IE_c[2] + 4*IE_c[1]))/54.75**5; N3s = W*3*(Eu3s + conv*(9*IE_s[2] + 4*IE_s[1]))/54.75**5
K0 = k0[0]*N0 + k0[1]*N1 + k0[2]*N2; K2 = k2[0]*N0 + k2[1]*N1 + k2[2]*N2
out['F4_affine'] = dict(K0_rate=K0, K2_rate=K2, I3_const=k2[3]*N3c, I3_slope=k2[3]*N3s, PHI_C=PHI_C, PHI_S=PHI_S, IE_c=IE_c, IE_s=IE_s)
print(f'F4 uniform: K0={K0:.2e}/tau K2={K2:.2e}/tau  I3 <= {k2[3]*N3c:.2e} + {k2[3]*N3s:.2e}(T-tau1)   [Phi={PHI_C:.3e}+{PHI_S:.2e}(T-tau1)]')
# ---- (B) ----
def sup1(Eu0, Eu1, Eu2, ell):
    return math.sqrt(Eu0/ell + 2*math.sqrt(Eu0*Eu1)), math.sqrt(Eu1/ell + 2*math.sqrt(Eu1*(Eu1 + Eu2)))
abar = 7.57*sq if ABAR is None else float(ABAR); cmin = 0.25 - abar - abar/50; r0 = 1.4*(((0.00244 + 63.5*al)*sq if BETA is None else float(BETA)) + 1.02*al*abar)
for yE in (48., 50.):
    eps = yE**-4; e = 2*p - gam; IB = yE**gam*((Rout**e - yE**e)/e); Cexp = max(Gc)*eps
    EB = []; gp = []
    for j in range(3):
        g = gam*cmin - (abs(1 - 2*j) + 1)*abar/yE - eps*Gc[j]
        coup = sum(Gk[j].get(k, 0)*eps*EB[k] for k in range(j))
        EB.append(E[j]*math.exp(Cexp)*math.exp(4*2*abar/yE) + 4*(r0*p**j)**2*IB/g**2 + 2*coup/g); gp.append(g)
    cv = (Rout/28.)**gam; S0, S1 = sup1(cv*EB[0], cv*EB[1], cv*EB[2], math.log(Rout/yE))
    out[f'B_y{yE:g}'] = dict(gamma_prime=gp, E_B=EB, S0_B=S0, S1_B=S1, valid=min(gp) > 0)
    print(f'B yE={yE:g}: gamma_prime={[round(x,4) for x in gp]} E_B={[f"{x:.3g}" for x in EB]} S0_B(v)={S0:.4f} S1_B(v)={S1:.3f}', '' if min(gp) > 0 else 'INVALID')
# ---- (A) core S0_f,S1_f and H6 strip (hand-off data x e^C, windows present at tau1, y in [50,150]) ----
Um = lambda y: 1.3056*al*y**p
S0A, S1A = d['result']['S0'], d['result']['S1']
print(f'A core [28,60]: S0_v={float(S0A):.4f} S1_v={float(S1A):.3f}; S0_f=|U|+S0_v={Um(60)+float(S0A):.4f}; S1_f=y|f\'\'|<=1.05U+S1_v={1.05*Um(60)+float(S1A):.3f}')
ho = a['handoff']; eC = math.exp(max(Gc)*50.**-4); cv = (150/28.)**gam; ell = math.log(150/50.)
S0h, S1h = sup1(cv*ho[0]*eC, cv*ho[1]*eC, cv*ho[2]*eC, ell)
K = {y: 50./y*(1.05*Um(y)) + 50./y*S1h for y in (50., 100., 150.)}
print(f'H6 strip [50,150] from hand-off E(tau1)={ho}: S0_h={S0h:.4f} S1_h={S1h:.3f}; curvature*L (L=50) = {K}; limit 1.5M=0.18')
out['H6'] = dict(S0h=S0h, S1h=S1h, curvL=K)
json.dump(out, open(f'lemmaS_F4B_{tag}.json', 'w'), indent=1)
