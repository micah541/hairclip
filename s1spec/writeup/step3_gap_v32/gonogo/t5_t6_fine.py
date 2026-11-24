"""GO/NO-GO tests T5 + T6 (FLOAT): the fine regime on [150, 250) lambda_1 with the fixed-L family (L = 50, X = 3, centres >= 300).
Centres rho = 150 and 225 (windows [150,225], [225,337.5]); weight W <= e^{mu(sigma+a)} on the inner tail, mu = 4, a = 0.3, >= wmin on the window.
  data    : hand-off shells (H2) on r <= 56 with the tail weight, plus the (H3) windows (unweighted, weight <= 1)
  inflow  : 1.70 IWf Gw * 2 (E^A_j + dU_j),  IWf = (4/mu) e^{mu(ln(56/rho)+a)},  dU_j = (p^j 0.0245)^2 ln 2   (z_j of U(alpha(tau)) - U_1 on the ramp)
  ramp_in : Gw IWf 2 KS [ sum_{k<=j} K_jk eps(28) 2 (E^A_k + dU_k) + K_{j,j+1} X_{j+1} ]                      (as in T1)
  ramp_out: T6, ramp [250, 379.5], members of unit 50 (pieces [250,330] member 300, [330,379.5] member 380), Lemma R R4(L=50) = 0.366,
            coefficients ramp_out_250_379.5.json (omega = 28/y normalisation) converted to the unweighted normalisation by y/28 <= 13.6
Outputs: Eu_j on the windows, S0, S1, and the collar contributions to the Lemma O a-priori bounds in units L = 50:
  K2 = 1.03 (50/y)(1.05 U_1(1.5 y) + S1) at y = 150;   K3_strip = sum_w (50/y_w)^1.5 (sqrt Eu2 + sqrt Eu1 + 0.0252 U_1 sqrt(ln 1.5));   K3_total = sqrt(K3_strip^2 + 2*0.03^2)."""
import json, math
K = {0: {}, 1: {}, 2: {}}
for j, d in json.load(open('../feas_P.json'))['4.0']['K'].items():
    for k, v in d.items(): K[int(j)][int(k)] = float(v)
amin = (16/17)**2; eps28 = 28.**-4; gp = [0.246, 0.246, 0.143]; l = math.log(1.5); p = 1.0246; al = 1.5e-4
U1 = lambda y: 1.3056*al*y**p; mu, a, wmin, Gw = 4, 0.3, 0.73, 1.41
# ---- T6: outer ramp with unit-50 members
co = json.load(open('t6_ramp_out_250_379.5.json'))['c']; C1, C2, G3, R4 = 0.00747719, 0.0395462, 0.0239123, 0.366035
s0, m2, n3, M = 0.08, 0.03, 0.004, 0.12; N = [0.0]*4
for (a_, b_) in ((250., 330.), (330., 379.5)):
    L = 50.; k2 = m2 + C2; k3 = n3 + G3; k4 = math.sqrt(2)*M + math.sqrt(R4); pm = math.tan(math.atan(s0) + C1); Fb = math.tan(math.atan(0.0928) + C1)
    sp_, sf = 1 + pm*pm, 1 + Fb*Fb; ln = math.sqrt((b_ - a_)/L); ell = math.log(b_/a_); Ub = U1(b_); Ks = k3 + 3*pm*k2**2*ln; D1 = sf**1.5*(b_/L)*k2
    N2f = (b_/L)**1.5*sf**2*(sp_**0.25*Ks + 3*Fb*k2**2*ln); N3f = (b_/L)**2.5*sf**2.5*sp_**0.25*(k4 + 10*(pm + Fb)*k2*Ks + (6 + 15*(pm*pm + Fb*Fb))*k2**3*ln)
    n = [(Fb + Ub)*math.sqrt(ell), (D1 + p*Ub)*math.sqrt(ell), N2f + (D1 + p*p*Ub)*math.sqrt(ell), N3f + 3*N2f + (D1 + p**3*Ub)*math.sqrt(ell)]
    N = [x + y*y for x, y in zip(N, n)]
conv = 379.5/28; Dout = [conv*sum(co[str(j)][str(k)]*N[k] for k in range(j + 2)) for j in range(3)]
print('T6: N_k^2 =', ['%.3g' % x for x in N], ' outer-ramp contribution (unweighted) Delta_j =', ['%.3g' % x for x in Dout])
out = dict(T6=dict(N_sq=N, Delta_unweighted=Dout))
# CLEAN-GO row (this agent): E^A, gamma' certified from lemmaS_A_rerun/rigor/assembly_arb.json; Gw from t8 (crude Ghat_22 binds at y = 28 for every centre);
# profile mismatch from Prop. B' (|dU| <= 1.3056 y^p sup|beta - p a alpha|/(p(1/4 - abar)) = 1.06e-6 at y = 56, instead of 0.0245) and Minkowski (sqrt E + sqrt dU)^2
_AS = json.load(open('lemmaS_A_rerun/rigor/assembly_arb.json')); _T8 = json.load(open('t8_L4_cost.json'))['cases']['crude|mu=4,nu=4,a=0.3']
_dUcg = 1.3056*56**p*(3.33811e-9 + p*2.19103e-7*1.5e-4)/(p*(0.25 - 2.19103e-7))
CASES = (('v3.1 thresholds (H2) = (2e-7, 1e-4, 0.2)', (3.11e-7, 1.81e-4, 0.247), (2e-7, 1e-4, 0.2)), ('tightened (H2) = (1e-8, 1e-6, 1e-3)', (1.2e-7, 8.2e-5, 0.048), (1e-8, 1e-6, 1e-3)),
         ('FINAL: certified E^A, Gw = 1.64', (1.2110055e-7, 8.1741747e-5, 0.048044818), (1e-8, 1e-6, 1e-3)),
         ('CLEAN-GO: certified E^A (json), Gw (t8), dU (Prop. B\')', tuple(_AS['E']), (1e-8, 1e-6, 1e-3)))
for name, EA, H2 in CASES:
    CG = name.startswith('CLEAN-GO'); Gw = _T8['Gw'] if CG else (1.64 if name.startswith('FINAL') else 1.41); gp = list(_AS['gp']) if CG else [0.246, 0.246, 0.143]
    comb = (lambda e, d: (math.sqrt(e) + math.sqrt(d))**2) if CG else (lambda e, d: 2*(e + d))
    for KS in (1, 4):
        Eus = {}
        for rho in (150., 225.):
            tail = math.exp(mu*(math.log(56/rho) + a)); IWf = 4/mu*tail; dU = [(p**j*(_dUcg if CG else 0.0245))**2*math.log(2) for j in range(3)]
            X = {1: eps28*comb(EA[1], dU[1]), 2: gp[0]*EA[0]/(2*amin), 3: gp[1]*EA[1]/(2*amin)}
            H3 = [7e-6, 3.5e-4, 3.5e-4]; B = []
            for j in range(3):
                data = Gw*(tail*(60/28)*H2[j] + H3[j])
                infl = 1.70*IWf*Gw*comb(EA[j], dU[j])
                rin = Gw*IWf*2*KS*(sum(K[j][k]*eps28*comb(EA[k], dU[k]) for k in range(j + 1)) + K[j][j + 1]*X[j + 1])
                B.append((data + infl + rin + Dout[j], data, infl, rin, Dout[j]))
            Eus[rho] = [b[0]/wmin for b in B]; parts = B
            if rho == 150.: p150 = B
        Eu = Eus[150.]; S0 = math.sqrt(Eu[0]/l + 2*math.sqrt(Eu[0]*Eu[1])); S1 = math.sqrt(Eu[1]/l + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
        K2 = 1.03*(50/150)*(1.05*U1(225) + S1)
        K3s = sum((50/y)**1.5*(math.sqrt(Eus[y][2]) + math.sqrt(Eus[y][1]) + 0.0252*U1(1.5*y)*math.sqrt(l)) for y in (150., 225.)); K3 = math.sqrt(K3s**2 + 2*0.03**2)
        r = dict(case=name, KS=KS, Eu_150=Eu, Eu_225=Eus[225.], parts_150_level2=dict(total=p150[2][0], data=p150[2][1], inflow=p150[2][2], ramp_in=p150[2][3], ramp_out=p150[2][4]),
                 S0=S0, S1=S1, K2=K2, K2_margin=0.18/K2, K3_strip=K3s, K3_total=K3, K3_margin=0.18/K3, slope_f=U1(225) + S0)
        out[f'{name}|KS={KS}'] = r
        print(f"{name}, KS={KS}: Eu(150) = ({Eu[0]:.2e}, {Eu[1]:.2e}, {Eu[2]:.3g}) [level 2: data {p150[2][1]:.2e}, inflow {p150[2][2]:.2e}, inner ramp {p150[2][3]:.2e}, outer ramp {p150[2][4]:.2e}]"
              f"  S0 = {S0:.3f}, S1 = {S1:.3f};  K2 = {K2:.3f} (x{0.18/K2:.2f}),  K3 = {K3:.3f} (x{0.18/K3:.2f}) [strip {K3s:.3f}]")
json.dump(out, open('t5_t6_fine.json', 'w'), indent=1)
