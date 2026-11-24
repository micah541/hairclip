"""GO/NO-GO test T1 (FLOAT shell model): does the crude uniformly local family close?
For any centre rho the weighted energy obeys  E^rho_j(tau) <= data_j + inflow_j + ramp_j [+ foot_j + outer_j]  (non-vanishing weight: no neighbour terms; rate <= 0 by the damping factor).
  weight   W(sigma) = 1/(e^{-mu(sigma+a)} + 1 + e^{nu(sigma-ln1.5-a)}), sigma = ln(r/rho);  >= w_min on the window [rho, 1.5 rho];
  IW       = int W(sigma(tau)) dtau = 4 int W dsigma  (the rescaled ramp [28,56] sweeps the weight at speed 1/4)
  Gw       = 1/(1 - Ghat_22/28^4): upper bound of the damping factor (1 + g y^-4)
  inflow_j = 1.70 IW Gw (sqrt E^A_j + sqrt dU_j)^2          (c max(2 S S')/ln2 * sup(y/28) = 0.25*3.40*2; Minkowski for v_rho = v_A + dU on the ramp)
  ramp_j   = Gw IW 2 KS [ sum_{k<=j} K_jk eps(28) (sqrt E^A_k + sqrt dU_k)^2 + K_{j,j+1} X_{j+1} ],  X_1 = eps(28) E^A_1,  X_k = gamma'_{k-2} E^A_{k-2}/(2 a_min) (level-(k-2) dissipation of (A), per unit tau)
  data_j   = Gw (60/28) H2_j               (hand-off shell thresholds; windows are smaller)
Then Eu_j = B_j / w_min on the window,  S0^2 = Eu0/l + 2 sqrt(Eu0 Eu1),  S1^2 = Eu1/l + 2 sqrt(Eu1 (Eu1 + Eu2)),  l = ln 1.5,
and the box needs  U + S0 <= 0.11,  1.05 U + S1 <= 0.52  with U = 1.3056 alpha_1 84^p = 0.0184 (window [56, 84] at launch).
K_jk: dimensionless ramp constants of feas_P.py (gamma = -4 tables, no L^4 extras) times a safety factor KS.
CLEAN-GO case (Nov 12, this agent; replaces the previous 'FINAL' row, which used typed values):
  E^A, gamma'  : CERTIFIED (arb) from lemmaS_A_rerun/rigor/assembly_arb.json (Lemma S(A) rerun at (H2) = (1e-8, 1e-6, 1e-3))
  Ghat_22      : t8_L4_cost.json (general-weight L^4, honest Gamma_m, kappa, kappa2 incl. the damping factor, fixed point in g)
  dU           : profile mismatch U(alpha(tau)) - U_rho on the ramp, from Prop. B' (step3_proof_v32 l.236): d/dtau[alpha lambda^-p] = (beta - p a alpha) lambda^-p, so
                 |dU| <= 1.3056 y^p sup|beta - p a alpha| / (p (1/4 - abar))  (certified sup|beta| <= 3.33811e-9, abar = 2.19103e-7; |alpha| <= alpha_1 = 1.5e-4)
  foot_j       : inner-ramp L^4 foot terms (t8 estimate, per weight)
  outer_j      : outer ramp (Lemma R) contribution, T6's unweighted Delta_j (t5_t6_fine.out) times W <= 1 (conservative; covers every rho <= 380 lambda_1)"""
import math, json, itertools
import numpy as np
K = {0: {0: 1.0e3, 1: 1.0e3}, 1: {0: 2e-12, 1: 1270., 2: 725.}, 2: {0: 2.3e-11, 1: 582., 2: 3100., 3: 704.}}
Kf = json.load(open('../feas_P.json'))['4.0']['K']
for j in ('0', '1', '2'):
    for k, v in Kf[j].items(): K[int(j)][int(k)] = float(v)
amin = (16/17)**2; eps28 = 28.**-4; GP0 = [0.246, 0.246, 0.143]; l = math.log(1.5); p = 1.0246; U = 1.3056*1.5e-4*84**p
def Wfun(mu, nu, a, x): return 1/(np.exp(-mu*(x + a)) + 1 + np.exp(nu*(x - l - a)))
def run(EA, H2, mu, nu, a, KS, Ghat=1.79e5, gp=GP0, dU=(0., 0., 0.), foot=(0., 0., 0.), outer=(0., 0., 0.)):
    xs = np.linspace(-15, 15, 300001); W = Wfun(mu, nu, a, xs)
    IW = 4*float(np.sum(W)*(xs[1] - xs[0])); wmin = float(np.min(W[(xs >= 0) & (xs <= l)])); Gw = 1/(1 - Ghat/28.**4)
    Ek = [(math.sqrt(EA[k]) + math.sqrt(dU[k]))**2 for k in range(3)]
    X = {1: eps28*Ek[1], 2: gp[0]*EA[0]/(2*amin) + eps28*p**4*dU[0], 3: gp[1]*EA[1]/(2*amin) + eps28*p**6*dU[0]}
    B = []; parts = []
    for j in range(3):
        infl = 1.70*IW*Gw*Ek[j]
        ramp = Gw*IW*2*KS*(sum(K[j][k]*eps28*Ek[k] for k in range(j + 1)) + K[j][j + 1]*X[j + 1])
        data = Gw*(60/28)*H2[j]
        B.append(data + infl + ramp + foot[j] + outer[j]); parts.append((data, infl, ramp, foot[j], outer[j]))
    Eu = [b/wmin for b in B]
    S0 = math.sqrt(Eu[0]/l + 2*math.sqrt(Eu[0]*Eu[1])); S1 = math.sqrt(Eu[1]/l + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(mu=mu, nu=nu, a=a, KS=KS, IW=IW, wmin=wmin, Gw=Gw, Ghat=Ghat, B=B, parts=parts, S0=S0, S1=S1, f=U + S0, yfp=1.05*U + S1, margin_M0=0.11/(U + S0), margin_M1=0.52/(1.05*U + S1))
def show(r):
    pp = r['parts'][2]
    return (f"  mu={r['mu']} nu={r['nu']} a={r['a']} KS={r['KS']}: IW={r['IW']:.2f} wmin={r['wmin']:.2f} Gw={r['Gw']:.3f} B=({r['B'][0]:.2e}, {r['B'][1]:.2e}, {r['B'][2]:.3f}) "
            f"[level-2 data/inflow/ramp/foot/outer = {pp[0]:.4f}/{pp[1]:.3f}/{pp[2]:.3f}/{pp[3]:.4f}/{pp[4]:.4f}; level-1 ramp share {r['parts'][1][2]/r['B'][1]:.2f}]  "
            f"|f|<={r['f']:.3f} (x{r['margin_M0']:.2f})  y|f'|<={r['yfp']:.3f} (x{r['margin_M1']:.2f})")
EA31 = (3.11e-7, 1.81e-4, 0.247); H31 = (2e-7, 1e-4, 0.2)
# tightened (H2): thresholds (1e-8, 1e-6, 1e-3); (A) floors from lemmaS_F10p parts (inflow + ramp + foot)  [previous agent's estimate, kept for continuity]
Ht = (1e-8, 1e-6, 1e-3); EAt = (1.1e-7 + Ht[0], 8.1e-5 + Ht[1], 0.047 + Ht[2])
out = []
for name, EA, H2 in (('v3.1 thresholds (H2) = (2e-7, 1e-4, 0.2)', EA31, H31), ('tightened (H2) = (1e-8, 1e-6, 1e-3), estimated floor', EAt, Ht)):
    print(name)
    for mu, nu, a in ((4, 2, 0.5), (4, 4, 0.3), (6, 4, 0.3), (4, 4, 0.0)):
        for KS in (1, 2, 4):
            r = run(EA, H2, mu, nu, a, KS); r['case'] = name; out.append(r); print(show(r))
# ---------------- CLEAN-GO case: certified E^A, honest Ghat, dU, foot, outer -----------------------------------------------------
AS = json.load(open('lemmaS_A_rerun/rigor/assembly_arb.json')); EAc = tuple(AS['E']); gpc = tuple(AS['gp'])
T8 = json.load(open('t8_L4_cost.json')); T7 = json.load(open('t7_step4_H2_tight.json'))
assert tuple(float(x) for x in T7['tight']['THR']) == Ht and T7['tight']['still_not_binding']
bbar, abar, alph = 3.33811e-9, 2.19103e-7, 1.5e-4
dUmax = 1.3056*56**p*(bbar + p*abar*alph)/(p*(0.25 - abar)); dU = tuple((p**j*dUmax)**2*math.log(2) for j in range(3))
DT6 = (9.86e-08, 1.04e-06, 0.0161)                                    # t5_t6_fine.out, T6 outer-ramp contribution (unweighted)
footw = {1: 1.5707250797145307, 2: 2862.4867820255}                   # rigor/weights_exact (certified, power weight, (A) ramp)
name = 'CLEAN-GO: certified E^A (arb), Ghat from t8 (general-weight L^4), dU (Prop. B\'), L^4 ramp foot, outer ramp'
print(name); print(f'  E^A = {EAc}, gamma\' = {gpc}; |dU| <= {dUmax:.4g} on [28,56] -> dU_j = {["%.2e" % x for x in dU]}')
flips = {}
for mu, nu, a in ((4, 2, 0.5), (4, 4, 0.3), (4, 4, 0.0)):
    c = T8['cases'].get(f'crude|mu={mu},nu={nu},a={a}')
    if not c or not c['closed']: print(f'  mu={mu} nu={nu} a={a}: rate criterion fails (t8)'); continue
    G22 = c['Ghat_L4'][2]; Gw = 1/(1 - G22/28.**4); rr = max(c['L4_extras_ratio_general_over_power'], 1.)
    xs = np.linspace(-15, 15, 300001); IW = 4*float(np.sum(Wfun(mu, nu, a, xs))*(xs[1] - xs[0]))
    F2 = eps28*Gw*IW*2*EAc[2]; JG = 2*(gpc[1]*EAc[1]/(2*amin))*IW*Gw
    s33 = 1.17e4*rr*math.sqrt(eps28*Gw*2*EAc[2])*(F2/math.sqrt(math.log(2)) + math.sqrt(2)*F2**0.75*JG**0.25)
    foot = (0., footw[1]*rr*F2, footw[2]*rr*F2 + s33)
    for KS in (1, 2, 4):
        r = run(EAc, Ht, mu, nu, a, KS, Ghat=G22, gp=gpc, dU=dU, foot=foot, outer=DT6); r['case'] = name; out.append(r); print(show(r))
    # what flips it (margin_M1 = 1 or margin_M0 = 1), KS = 4 (pessimistic) and KS = 1
    def bis(fun, lo, hi):
        for _ in range(80):
            mid = math.sqrt(lo*hi); (lo, hi) = (mid, hi) if fun(mid) > 1 else (lo, mid)
        return lo
    for KS in (1, 4):
        mm = lambda r: min(r['margin_M0'], r['margin_M1'])
        e2 = bis(lambda e: mm(run((EAc[0], EAc[1], e), Ht, mu, nu, a, KS, G22, gpc, dU, foot, DT6)), EAc[2], 100.)
        gh = bis(lambda g: mm(run(EAc, Ht, mu, nu, a, KS, g, gpc, dU, foot, DT6)), G22, 28.**4*(1 - 1e-12))
        ks = bis(lambda k: mm(run(EAc, Ht, mu, nu, a, k, G22, gpc, dU, foot, DT6)), KS, 1e4)
        flips[f'mu={mu},nu={nu},a={a},KS={KS}'] = dict(E2A_flip=e2, Ghat22_flip=gh, KS_flip=ks)
        print(f'    flips (KS={KS}): E2^A -> {e2:.4f} (x{e2/EAc[2]:.2f});  Ghat_22 -> {gh:.4g} (x{gh/G22:.2f}; Gw {1/(1 - gh/28.**4):.2f});  ramp constants -> KS = {ks:.1f}')
json.dump(dict(rows=out, flips=flips, dU=dict(max=dUmax, dU=dU)), open('t1_crude_shell_model.json', 'w'), indent=1)
