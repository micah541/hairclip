"""(B) proposition numbers: windows fixed in physical space y in [yE,350], levels 0-2, no ramp inflow.
Entry data: F1 core energies (E_j^{F1}, including hand-off) at y>=yE; Ghat from run_F1 'scan' (M0=.07,M1=1.4,M2..4=3) quasilinear."""
import json, math
F1 = json.load(open('lemmaS_F1.json'))['lin_gam1.0']['result']
S = json.load(open('lemmaS_F1_scan.json'))['gam1.0']
G = [float(S['consts'][f'Ghat_{j}'][str(j)]) for j in range(3)]
Gc = {j: {int(k): float(v) for k, v in S['consts'][f'Ghat_{j}'].items()} for j in range(3)}
sq = 2.49e-7; al = 5e-4; gam = 1.0; p = 1.0246; Rout = 350.
abar = 7.57*sq; bbar = abar; cmin = 0.25 - abar - bbar/50
res = {}
for yE in (48., 50.):
    eps = yE**-4; e = 2*p - gam; Iint = yE**gam*((Rout**e - yE**e)/e)
    Cexp = max(G)*eps*1.0     # e^{C yE^-4}, C = max Ghat_jj
    E = []; gp = []
    for j in range(3):
        g = gam*cmin - (abs(1 - 2*j) + 1)*bbar/yE - eps*G[j]
        coup = sum(Gc[j].get(k, 0)*eps*E[k] for k in range(j))
        r0 = 1.4*((0.00244 + 63.5*al)*sq + 1.02*al*abar)
        F2 = (r0*p**j)**2*Iint
        Ein = float(F1['E'][j])*math.exp(Cexp)*math.exp(4*(abar + bbar)/yE)
        E.append(Ein + 4*F2/g**2 + 2*coup/g); gp.append(g)
    res[yE] = dict(eps=eps, gamma_prime=gp, E=E, Cexp=Cexp, Ghat=G)
    print(f'yE={yE}: eps={eps:.2e} gamma\'={[round(x,4) for x in gp]} E_B={[f"{x:.3g}" for x in E]} e^C={math.exp(Cexp):.3f}')
json.dump({str(k): v for k, v in res.items()}, open('lemmaS_B.json', 'w'), indent=1)
