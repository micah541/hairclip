"""Phase 1, task P3: launch model for the crude uniformly local family (Lemma U).  FLOAT, every number below is [est].
This is the 'plan's way' accounting of GONOGO error 4 (launch at y_L = 56), to be replaced by the certified assembly (Phase 2).
Notation and formulas: launch_lemma.tex (this directory), sections 2-4.  Weight W(sigma) = 1/(e^{-mu(sigma+a)} + 1 + e^{nu(sigma - ln1.5 - a)}),
sigma = ln(r/rho); Theta = W (1 + g y^-4) <= Gw W, Gw = 1/(1 - Ghat_22/28^4); inner ramp [28,56] fixed in y; outer ramp [190,379.5] lambda_1 (Lemma R).

Centres.
  LATE   (rho < 56 lambda_1): launched at tau_rho with rho/lambda(tau_rho) = 56.  At launch v_rho = v_A (same alpha, same y), so
         D_j = Gw [ C_A E^A_j + ( sqrt(q B_j) + sqrt(Dd_j) )^2 ],
         C_A = sup_{28<=y<=350} W(ln(y/56)) y/28                (E^A is omega-normalised, omega = 28/y; (A) uncut on [14,350])
         q   = sum_k sup_{window k} W / w_min, windows [350*1.5^k, 350*1.5^(k+1)] (neighbour centres, all launched), Dd_j: profile
               mismatch U_rho - U_rho' between neighbours, |.| <= 1.3056 y^p sup|beta - p a alpha|/(p(1/4 - abar)) (Prop. B'), weighted.
  PRESENT(rho = y_c lambda_1, y_c in [56, 253]): launched at tau_1 with the hand-off data
         D_j = Gw [ sup_{28<=y<=60} W(ln(y/y_c)) y/28 * H2_j + sum_w sup_{[y_w,1.5y_w]} W(ln(y/y_c)) (1.5 y_w/28) H3_j ],  y_w = 50*1.5^k, k=0..4.
After launch (both): the weight's INNER tail only sweeps the ramp [28,56]:  IW_in(y_c) = 4 int_{-inf}^{ln(56/y_c)} W dsigma  (y_c = 56 for late)
         inflow_j = 1.70 IW_in Gw Ek_j,  ramp_j = Gw IW_in 2 KS [sum_k K_jk eps(28) Ek_k + K_{j,j+1} X_{j+1}],  foot_j (t8 formula with IW_in),
         outer_j = T6 Delta_j (Lemma R ramp, W <= 1),  Ek_j = (sqrt E^A_j + sqrt dU_j)^2  -- exactly as t1_crude_shell_model.py, with IW -> IW_in.
         B_j = D_j + inflow_j + ramp_j + foot_j + outer_j  (fixed point in B for the LATE centres: neighbour factor 1/(1 - Gw q)).
Box on the centre's own window (Sobolev, l = ln 1.5):  Eu = B/w_min, S0^2 = Eu0/l + 2 sqrt(Eu0 Eu1), S1^2 = Eu1/l + 2 sqrt(Eu1(Eu1 + Eu2)),
         |f| <= U_win + S0 <= 0.11,  y|f'| <= 1.05 U_win + S1 <= 0.52;  U_win = 1.3056 alpha 84^p (late), 1.3056 alpha (1.5 y_c)^p (present).
Also: T1's sweep model (whole weight, hand-off data) reproduced, with and without (H3) (GONOGO error 3), and the T3 Young-route addition
to E_2^A recomputed with the launch-model window energies.  Run: python3 -B launch_model.py  (reads only; writes launch_model.{json,out} here)."""
import sys, os, json, math, hashlib
sys.dont_write_bytecode = True
import numpy as np
GO = '/home/claude/s1spec/writeup/step3_gap_v32/gonogo'
INP = {'assembly_arb': f'{GO}/lemmaS_A_rerun/rigor/assembly_arb.json', 't8': f'{GO}/t8_L4_cost.json', 'feas_P': f'{GO}/../feas_P.json',
       't1': f'{GO}/t1_crude_shell_model.json', 't3': f'{GO}/t3_A_outer.json', 't5_t6_out': f'{GO}/t5_t6_fine.out'}
AS = json.load(open(INP['assembly_arb'])); T8 = json.load(open(INP['t8'])); T1 = json.load(open(INP['t1'])); T3 = json.load(open(INP['t3']))
EA = tuple(AS['E']); GP = tuple(AS['gp'])                                           # [cert] (conditional on the box, GONOGO error 6)
K = {0: {0: 1.0e3, 1: 1.0e3}, 1: {0: 2e-12, 1: 1270., 2: 725.}, 2: {0: 2.3e-11, 1: 582., 2: 3100., 3: 704.}}
for j, row in json.load(open(INP['feas_P']))['4.0']['K'].items():                   # [est] power-weight ramp constants (GONOGO error 2)
    for k, v in row.items(): K[int(j)][int(k)] = float(v)
amin = (16/17)**2; eps28 = 28.**-4; l = math.log(1.5); p = 1.0246; A1 = 1.5e-4; A1_step4 = 1.1739e-4
H2 = (1e-8, 1e-6, 1e-3); H3 = (2e-7, 1e-5, 1e-5)                                    # hand-off thresholds (H2), (H3) [paper: theorem hypotheses]
bbar, abar = 3.33811e-9, 2.19103e-7                                                 # [cert] verifier sups
DT6 = (9.86e-08, 1.04e-06, 0.0161)                                                  # [est] t5_t6_fine.out, T6 outer ramp (Lemma R), unweighted
footw = {1: 1.5707250797145307, 2: 2862.4867820255}                                 # [cert] power weight; use here [est] (error 6)
c_dU = 1.3056*(bbar + p*abar*A1)/(p*(0.25 - abar))                                  # |dU| <= c_dU y^p   (Prop. B', launch_lemma.tex Lemma 3)
dU = tuple((p**j*c_dU*56**p)**2*math.log(2) for j in range(3))                      # on the inner ramp [28,56] (as t1)
xs = np.linspace(-20, 20, 400001); dx = xs[1] - xs[0]
def Wf(mu, nu, a, x): return 1/(np.exp(-mu*(x + a)) + 1 + np.exp(nu*(x - l - a)))
def supW(mu, nu, a, s0, s1, fac=None, n=4001):
    x = np.linspace(s0, s1, n); w = Wf(mu, nu, a, x); w = w*fac(x) if fac else w; return float(w.max())*(1 + 1e-3)   # grid sup, +0.1% [est]
def weight(mu, nu, a):
    W = Wf(mu, nu, a, xs); return dict(IW=4*float(W.sum()*dx), wmin=float(W[(xs >= 0) & (xs <= l)].min()))
def IW_in(mu, nu, a, yc): m = xs <= math.log(56/yc); return 4*float(Wf(mu, nu, a, xs[m]).sum()*dx)
def foot_terms(IW, Gw, rr):
    F2 = eps28*Gw*IW*2*EA[2]; JG = 2*(GP[1]*EA[1]/(2*amin))*IW*Gw
    s33 = 1.17e4*rr*math.sqrt(eps28*Gw*2*EA[2])*(F2/math.sqrt(math.log(2)) + math.sqrt(2)*F2**0.75*JG**0.25)
    return (0., footw[1]*rr*F2, footw[2]*rr*F2 + s33)
def evol(IW, Gw, KS, rr):
    Ek = [(math.sqrt(EA[k]) + math.sqrt(dU[k]))**2 for k in range(3)]
    X = {1: eps28*Ek[1], 2: GP[0]*EA[0]/(2*amin) + eps28*p**4*dU[0], 3: GP[1]*EA[1]/(2*amin) + eps28*p**6*dU[0]}
    ft = foot_terms(IW, Gw, rr); out = []
    for j in range(3):
        infl = 1.70*IW*Gw*Ek[j]; ramp = Gw*IW*2*KS*(sum(K[j][k]*eps28*Ek[k] for k in range(j + 1)) + K[j][j + 1]*X[j + 1])
        out.append(dict(inflow=infl, ramp=ramp, foot=ft[j], outer=DT6[j]))
    return out
def box(B, wmin, Uwin):
    Eu = [b/wmin for b in B]
    S0 = math.sqrt(Eu[0]/l + 2*math.sqrt(Eu[0]*Eu[1])); S1 = math.sqrt(Eu[1]/l + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(S0=S0, S1=S1, f=Uwin + S0, yfp=1.05*Uwin + S1, M0=0.11/(Uwin + S0), M1=0.52/(1.05*Uwin + S1))
def neighbour(mu, nu, a, wmin):
    s0 = math.log(350/56); q = 0.0; Dd = [0.0, 0.0, 0.0]
    for k in range(60):
        q += supW(mu, nu, a, s0 + k*l, s0 + (k + 1)*l)/wmin
    # weighted energy of the profile mismatch beyond y = 350 (centre at 56):  int W(s - ln56) (p^j c_dU e^{ps})^2 ds
    s = np.linspace(math.log(350), math.log(350) + 30, 300001)
    base = float(np.sum(Wf(mu, nu, a, s - math.log(56))*(c_dU*np.exp(p*s))**2)*(s[1] - s[0]))
    return q, [p**(2*j)*base for j in range(3)]
U84 = 1.3056*A1*84**p
res = dict(inputs={}, late={}, present={}, t1_check={}, t3_young={}, dU=dict(c_dU=c_dU, dUmax56=c_dU*56**p, dU=dU))
for k, f in INP.items(): res['inputs'][f] = hashlib.sha256(open(f, 'rb').read()).hexdigest()
res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version)
print(f'E^A = {EA} [cert, conditional], gamma\' = {GP}; |dU| <= {c_dU:.4g} y^p ({c_dU*56**p:.3g} at y = 56) [paper; inputs cert]')
for (mu, nu, a) in ((4, 4, 0.3), (4, 2, 0.5), (4, 4, 0.0)):
    c = T8['cases'][f'crude|mu={mu},nu={nu},a={a}']; G22 = c['Ghat_L4'][2]; Gw = 1/(1 - G22/28.**4); rr = max(c['L4_extras_ratio_general_over_power'], 1.)
    wt = weight(mu, nu, a); wmin = wt['wmin']; IWi = IW_in(mu, nu, a, 56.)
    CA = supW(mu, nu, a, math.log(28/56), math.log(350/56), fac=lambda x: 2*np.exp(x))
    q, Dd = neighbour(mu, nu, a, wmin)
    key = f'mu={mu},nu={nu},a={a}'
    print(f'\n== weight ({mu},{nu},{a}): IW = {wt["IW"]:.3f}, IW_in(56) = {IWi:.3f}, w_min = {wmin:.3f}, Gw = {Gw:.3f}, C_A = {CA:.3f}, '
          f'q = {q:.4g}, neighbour factor 1/(1 - Gw q) = {1/(1 - Gw*q):.4f}, Dd_2 = {Dd[2]:.2e}  [est]')
    res['late'][key] = dict(IW=wt['IW'], IW_in=IWi, wmin=wmin, Gw=Gw, CA=CA, q=q, nfac=1/(1 - Gw*q), Dd=Dd, rows={})
    for KS in (1, 2, 4):
        ev = evol(IWi, Gw, KS, rr); B = [0.0, 0.0, 0.0]
        for it in range(200):                                                        # fixed point (neighbour term)
            D = [Gw*(CA*EA[j] + (math.sqrt(q*B[j]) + math.sqrt(Dd[j]))**2) for j in range(3)]
            Bn = [D[j] + sum(ev[j].values()) for j in range(3)]
            if max(abs(Bn[j] - B[j])/max(Bn[j], 1e-300) for j in range(3)) < 1e-14: B = Bn; break
            B = Bn
        bx = box(B, wmin, U84)
        res['late'][key]['rows'][KS] = dict(B=B, D=D, evol=ev, **bx)
        # reviewer variant [rev] (GONOGO error 4): same launch data and inflow, but ramp and foot terms with the FULL sweep IW (conservative)
        evr = evol(wt['IW'], Gw, KS, rr); evv = [dict(ev[j], ramp=evr[j]['ramp'], foot=evr[j]['foot']) for j in range(3)]; Br = [0.0]*3
        for it in range(200):
            Dr = [Gw*(CA*EA[j] + (math.sqrt(q*Br[j]) + math.sqrt(Dd[j]))**2) for j in range(3)]; Br = [Dr[j] + sum(evv[j].values()) for j in range(3)]
        bxr = box(Br, wmin, U84); res['late'][key]['rows'][KS]['rev_variant'] = dict(B=Br, **bxr)
        e = ev[2]
        print(f'  LATE   KS={KS}: B = ({B[0]:.2e}, {B[1]:.2e}, {B[2]:.3f})  [level 2: launch data {D[2]:.3f} / inflow {e["inflow"]:.3f} / ramp {e["ramp"]:.3f} / '
              f'foot {e["foot"]:.4f} / outer {e["outer"]:.4f}]  |f| <= {bx["f"]:.3f} (x{bx["M0"]:.2f})  y|f\'| <= {bx["yfp"]:.3f} (x{bx["M1"]:.2f})  [est]')
        print(f'         reviewer variant (ramp + foot with full IW): B_2 = {Br[2]:.3f}, M0 x{bxr["M0"]:.2f}, M1 x{bxr["M1"]:.2f}  [est]')
    # present centres: radius dependence; alpha_1 and step 4's certified |alpha|
    res['present'][key] = {}
    for KS in (1, 4):
        rows = []
        for yc in (56., 60., 70., 84., 100., 126.6, 150., 190/1.5, 200., 253.):
            IWy = IW_in(mu, nu, a, yc); ev = evol(IWy, Gw, KS, rr)
            dH2 = supW(mu, nu, a, math.log(28/yc), math.log(60/yc), fac=lambda x, yc=yc: yc*np.exp(x)/28)
            dH3 = sum(supW(mu, nu, a, math.log(yw/yc), math.log(1.5*yw/yc))*(1.5*yw/28) for yw in (50*1.5**k for k in range(5)))
            D = [Gw*(dH2*H2[j] + dH3*H3[j]) for j in range(3)]; B = [D[j] + sum(ev[j].values()) for j in range(3)]
            Uw = 1.3056*A1*(1.5*yc)**p; Uw4 = 1.3056*A1_step4*(1.5*yc)**p
            bx = box(B, wmin, Uw); bx4 = box(B, wmin, Uw4)
            rows.append(dict(yc=yc, IW_in=IWy, D=D, B=B, Uwin=Uw, **bx, M0_step4alpha=bx4['M0'], M1_step4alpha=bx4['M1']))
        res['present'][key][KS] = rows
        if (mu, nu, a) == (4, 4, 0.3):
            print(f'  PRESENT centres (launched at tau_1, hand-off (H2)+(H3) data), KS={KS}:  y_c, IW_in, B_2, U_win, |f| (xM0 vs 0.11; with step-4 |alpha|), y|f\'| (xM1)  [est]')
            for r in rows:
                print(f'    y_c = {r["yc"]:6.1f}: IW_in = {r["IW_in"]:.2e}, B = ({r["B"][0]:.1e}, {r["B"][1]:.1e}, {r["B"][2]:.2e}), U_win = {r["Uwin"]:.4f}, '
                      f'|f| <= {r["f"]:.4f} (x{r["M0"]:.2f}; x{r["M0_step4alpha"]:.2f}), y|f\'| <= {r["yfp"]:.4f} (x{r["M1"]:.2f})')
    # T1 sweep model reproduced (check), and with (H3) added to the data (GONOGO error 3)
    for KS in (1, 2, 4):
        ev = evol(wt['IW'], Gw, KS, rr)
        for withH3 in (False, True):
            dH3 = sum(supW(mu, nu, a, math.log(yw/56), math.log(1.5*yw/56))*(1.5*yw/28) for yw in (50*1.5**k for k in range(5))) if withH3 else 0.
            D = [Gw*((60/28)*H2[j] + dH3*H3[j]) for j in range(3)]; B = [D[j] + sum(ev[j].values()) for j in range(3)]; bx = box(B, wmin, U84)
            res['t1_check'].setdefault(key, {})[f'KS={KS},H3={withH3}'] = dict(B=B, **bx)
        r0 = res['t1_check'][key][f'KS={KS},H3=False']; r1 = res['t1_check'][key][f'KS={KS},H3=True']
        t1row = [x for x in T1['rows'] if x['case'].startswith('CLEAN-GO') and (x['mu'], x['nu'], x['a'], x['KS']) == (mu, nu, a, KS)][0]
        print(f'  T1 sweep model KS={KS}: M1 x{r0["M1"]:.3f} (t1 on disk x{t1row["margin_M1"]:.3f}), with (H3) x{r1["M1"]:.3f}; M0 x{r0["M0"]:.2f} -> x{r1["M0"]:.2f}  [est]')
# T3 Young route with the launch-model window energies (late centre, (4,4,0.3)), at the bootstrap level 2*B (hypothesis level)
c3, K3 = T3['c'], T3['K']
for KS in (2, 4):
    r = res['late']['mu=4,nu=4,a=0.3']; B = r['rows'][KS]['B']
    for lev in (1, 2):
        N = [1.71*lev*b/r['wmin'] for b in B]
        base = sum(c3['2'][str(k)]*N[k] for k in range(3))
        yg = (base + K3['2']['3']**2/(4*0.1*amin)*350.**-4*(28/350)*N[2])/GP[2]
        res['t3_young'][f'KS={KS},level={lev}xB'] = yg
        print(f'T3 Young route with launch-model N_k (KS={KS}, {lev}xB): adds {yg:.3g} to E_2^A = {EA[2]:.4g} (x{EA[2]/yg:.0f} below)  [est]')
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'launch_model.json'), 'w'), indent=1, default=float)
