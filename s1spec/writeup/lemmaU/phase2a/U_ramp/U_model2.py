#!/usr/bin/env python3
"""Lemma U, Phase 2a (agent B): float model U_model2 -- late-centre and present-centre M1/M0 margins with
  (i)  honest time-integrated ramp data (U_ramp.tex Prop. 1): X-terms from (A)'s exported dissipation integrals D_k, no rate model;
  (ii) the Ghat design of U_ramp.tex Sec. 4 ('data56': rate criterion from y = 56, ramp share of the rate as data; box (0.15, 0.52));
  (iii) pointwise damping factor D(y) = 1 + g y^-4 and joint (y, tau) sweeps instead of products of sups;
  (iv) E^A and D_k read from a JSON path (default: certified IBP E^A times the reviewer's repair factors (1.41, 1.13, 1.25));
  (v)  certified K, footw, C_eta^Xi, S33 coefficient from phase1/ramp/ramp_K_cert.json; Ghat from ghat_design_arb.json (this dir) / phase1.

  PYTHONHASHSEED=0 nice -n 19 python3 -B U_model2.py [--inputs agentA.json] [--tag name]

  --inputs JSON (all keys optional):  {"EA": [E0, E1, E2], "gp": [g0, g1, g2], "D": [D0, D1, D2]  (bounds on the time integrals
      iint_{tau_rho}^inf int_{[28,56]} eps*omega z_{k+2}^2 ds dtau exported by (A)), or instead "kappa": [k0, k1, k2] and "cD": c
      (then D_k = c * E^A_k / (kappa_k * a_min)), "amin": a_min, "repair": [r0, r1, r2] (multiplies EA; default (1,1,1) when EA is given),
      "args": {...}}.  The file's sha256 goes into the args block.
Labels: everything printed is [est] (float model) on [cert] inputs, except the typed items listed in the .md (H2, H3, DT6, abar, bbar, alpha_1).
SELF-CHECK (asserted): with every new ingredient off (cfg 'CHECKPOINT'), the late-centre launch row of phase1/checkpoint/t1_checkpoint.json
  (worst Ghat22, with the (2,2,3) placeholder, KS = 1 and 2) and its present-centre rows are reproduced to 1e-9 relative.
Outputs: U_model2[_tag].json (args block with input sha256), stdout -> U_model2[_tag].out."""
import sys, os, json, math, hashlib, time, argparse, platform
sys.dont_write_bytecode = True
import numpy as np
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
R = '/home/claude/s1spec/writeup'; GO = f'{R}/step3_gap_v32/gonogo'; PH = f'{R}/lemmaU/phase1'
ap = argparse.ArgumentParser(); ap.add_argument('--inputs', default=None); ap.add_argument('--tag', default='')
cli = ap.parse_args()
INP = dict(assembly_ibp=f'{GO}/lemmaS_A_ibp/rigor/assembly_arb.json', ghat=f'{PH}/ghat/ghat_cert.json', ramp=f'{PH}/ramp/ramp_K_cert.json',
           t1chk=f'{PH}/checkpoint/t1_checkpoint.json', design=f'{HERE}/ghat_design_arb.json', toppair=f'{HERE}/top_pair_est.json')
if cli.inputs: INP['agentA_inputs'] = os.path.abspath(cli.inputs)
J = {k: json.load(open(v)) for k, v in INP.items()}

# ---------------- typed constants (as t1_checkpoint / launch_model) ----------------
amin_t1 = (16/17)**2; eps28 = 28.**-4; l = math.log(1.5); p = 1.0246; A1 = 1.5e-4
H2 = (1e-8, 1e-6, 1e-3); H3 = (2e-7, 1e-5, 1e-5)                   # [paper] hand-off hypotheses
bbar, abar = 3.33811e-9, 2.19103e-7                                # [cert] verifier sups (typed); bbar here is sup|beta| (as t1)
bsup = 2.19103e-7                                                  # [cert] sup|b| (D1 of v3.2, typed): translation parameter, for the b-term
DT6 = (9.86e-08, 1.04e-06, 0.0161)                                 # [est] Lemma R outer ramp (t5_t6_fine.out), typed
c_dU = 1.3056*(bbar + p*abar*A1)/(p*(0.25 - abar))
dU = tuple((p**j*c_dU*56**p)**2*math.log(2) for j in range(3))
xs = np.linspace(-20, 20, 400001); dx = xs[1] - xs[0]
cm, cp = 0.25 - abar, 0.25 + abar + bsup/28                         # transport speed bounds (U_ramp.tex Lemma T): c- for sweeps, c+ for the inflow
KB = 1e-3                                                          # b-weight phi = exp(KB e^{-s}) (U_ramp.tex Lemma T(iii)); phi <= e^{KB/28}
PHI = math.exp(KB/28)
RRAT = math.exp(math.log(2)*3.5212668029853376)                    # [cert] sup over the ramp of (eps Theta)max/(eps Theta)min <= e^{ln2 sup|L|} (P2 geom|L)
LRAMP = 3.5212668029853376                                         # [cert] sup |L| on the ramp (P2)

def Wf(mu, nu, a, x): return 1/(np.exp(-mu*(x + a)) + 1 + np.exp(nu*(x - l - a)))
def supW(mu, nu, a, s0, s1, fac=None, n=4001):
    x = np.linspace(s0, s1, n); w = Wf(mu, nu, a, x); w = w*fac(x) if fac else w; return float(w.max())*(1 + 1e-3)
def IW_in(mu, nu, a, yc): m = xs <= math.log(56/yc); return 4*float(Wf(mu, nu, a, xs[m]).sum()*dx)
_wc = {}
def weight(mu, nu, a):
    k = (mu, nu, a)
    if k not in _wc:
        W = Wf(mu, nu, a, xs); wmin = float(W[(xs >= 0) & (xs <= l)].min())
        CA = supW(mu, nu, a, math.log(28/56), math.log(350/56), fac=lambda x: 2*np.exp(x))
        q = sum(supW(mu, nu, a, math.log(350/56) + i*l, math.log(350/56) + (i + 1)*l)/wmin for i in range(60))
        s = np.linspace(math.log(350), math.log(350) + 30, 300001)
        base = float(np.sum(Wf(mu, nu, a, s - math.log(56))*(c_dU*np.exp(p*s))**2)*(s[1] - s[0]))
        _wc[k] = dict(IW=4*float(W.sum()*dx), IWin56=IW_in(mu, nu, a, 56.), wmin=wmin, CA=CA, q=q, Dd=[p**(2*j)*base for j in range(3)])
    return _wc[k]

# ---------------- ramp geometry for the joint sweeps (s = ln y on [ln28, ln56], eta = S(t), t = (s - ln28)/ln2) ----------------
NY = 601; sR = np.linspace(math.log(28), math.log(56), NY); yR = np.exp(sR); tR = (sR - math.log(28))/math.log(2)
S = tR**3*(10 - 15*tR + 6*tR**2); Sp = 30*tR**2*(1 - tR)**2
eta2 = S**2; deta2 = 2*S*Sp/math.log(2)                            # (eta^2)' in s
XG = np.linspace(0., 16., 3201)                                    # x = c-(tau - tau_rho) shift of sigma
def sweep(f, yc, mu=4, nu=4, a=0.3):
    """(1/c-) int_0^inf sup_{y in ramp} f(y) W(ln(y/yc) - x) dx  (W increasing on sigma <= l/2, so W(sigma(tau)) <= W(ln(y/yc) - x)); grid +0.5%"""
    sig = np.log(yR/yc)[None, :] - XG[:, None]
    sup = (f[None, :]*Wf(mu, nu, a, sig)).max(axis=1)
    return float(np.sum(0.5*(sup[1:] + sup[:-1]))*(XG[1] - XG[0]))/cm*1.005
def supR(f, yc, mu=4, nu=4, a=0.3): return float((f*Wf(mu, nu, a, np.log(yR/yc))).max())*1.001

# ---------------- certified inputs ----------------
AS = J['assembly_ibp']; RK = J['ramp']; DS = J['design']
K_cert = {int(j): {int(k): float(v) for k, v in row.items()} for j, row in RK['4.0']['K'].items()}
CETA = RK['C_eta_Xi']; COEF_S33 = RK['foot']['2']['Kf_unnorm_S33'] + RK['cubic_remaining']['2|(2, 2, 2)']
W223 = RK['cubic_remaining']['2|(2, 2, 3)']
FOOTW = {1: RK['foot']['1']['footw'], 2: RK['foot']['2']['footw']}
G22_W = max(J['ghat']['sensitivity_Umax_0.09']['Ghat_jj'][2][1], RK['foot']['2']['Ghat_jj_with_foot_adds'])
# (2,2,3) in S33 form (U_ramp.tex Lemma C): coefficient (1/3) sup|(c eps Theta)'/(eps Theta)| <= (1/3) w223 (C_d + sup|L|), C_d = sup|c'|/sup|c| [est]
_PE = np.polynomial.Polynomial([0, 0, 0, 10, -15, 6])**2; _tt = np.linspace(0, 1, 200001)
_de = [np.abs(_PE.deriv(k)(_tt)).max()/math.log(2)**k for k in range(5)]
CD223 = max(_de[k + 1]/_de[k] for k in (1, 2, 3)) + 4.0                                  # [est] max_k sup|(eta^2)^(k+1)|/sup|(eta^2)^(k)| (= 17.23) + Gamma_1
COEF_223 = W223*(CD223 + LRAMP)/3
COEF_223_PLACEHOLDER = W223*CETA/3

# ---------------- E^A, D_k ----------------
EA_cert = tuple(AS['E']); GP_cert = tuple(AS['gp'])
REPAIR = (1.41, 1.13, 1.25)                                         # [rev] reviewer's repair factors (CHECKPOINT finding 1)
def inputs_default(kappa=(1.25, 1.25, 0.0), cD=2.0, amin=amin_t1, repair=REPAIR):
    EA = tuple(e*r for e, r in zip(EA_cert, repair))
    D = tuple(cD*EA[k]/(kappa[k]*amin) if kappa[k] > 0 else math.inf for k in range(3))
    return dict(EA=EA, gp=GP_cert, D=D, src=f'default: certified IBP E^A x {repair}; D_k = {cD} E^A_k/(kappa_k a_min), kappa = {kappa}, a_min = {amin:.4f}')
def inputs_json(d):
    rep = d.get('repair', (1., 1., 1.)); EA = tuple(e*r for e, r in zip(d.get('EA', EA_cert), rep)); gp = tuple(d.get('gp', GP_cert)); amin = d.get('amin', amin_t1)
    if 'D' in d: D = tuple(d['D'])
    else:
        kap = d.get('kappa', (1.25, 1.25, 0.0)); cD = d.get('cD', 2.0); D = tuple(cD*EA[k]/(kap[k]*amin) if kap[k] > 0 else math.inf for k in range(3))
    return dict(EA=EA, gp=gp, D=D, src=f'JSON {INP["agentA_inputs"]}')

# ---------------- configurations ----------------
def cfg_checkpoint(G22=G22_W, KS=1):
    """t1_checkpoint 'PHASE1 G22=WORST +(2,2,3)' row: every new ingredient off"""
    return dict(name='CHECKPOINT', EA=EA_cert, gp=GP_cert, D=None, K=K_cert, G22=G22, design='global28', KS=KS,
                pointwise=False, joint=False, honestX=False, energy12=False, coupling=False, bterm=False, honestS33=False,
                coef=(COEF_S33 + COEF_223_PLACEHOLDER), Ghat=None)
def cfg_new(inp, design='data56', box='widened (0.15,0.52), U .09, M0z .16', KS=1, **over):
    c = DS['cases'][box]; G22 = c['Ghat_jj'][2][1]
    Ghat = {j: {k: c['Ghat_full'][str(j)][str(k)] for k in range(j + 1)} for j in range(3)}
    cfg = dict(name=f'NEW {design} {box}', EA=inp['EA'], gp=inp['gp'], D=inp['D'], K=K_cert, G22=G22, design=design, KS=KS,
               pointwise=True, joint=True, honestX=True, energy12=True, coupling=True, bterm=True, honestS33=True,
               coef=(COEF_S33 + COEF_223), Ghat=Ghat, box=box)
    cfg.update(over); return cfg

def gD(cfg):
    G = cfg['G22']
    if cfg['design'] == 'global28':
        g = G/(1 - G/28.**4) if not cfg['pointwise'] else G/(1 - 4*abar - G/28.**4)
    else: g = G/(1 - 4*abar - G/56.**4)
    return g
def Dy(cfg, y): return 1 + gD(cfg)/y**4
def Gwof(cfg): return 1/(1 - cfg['G22']/28.**4) if not cfg['pointwise'] else Dy(cfg, 28.)*PHI

def foot_terms(cfg, Fj, JG, Ek2, Gw):
    """foot_1 = footw_1 F2, foot_2 = footw_2 F2 + S33 (sup-free form, U_ramp.tex eq. (S33))"""
    if cfg['honestS33']: pre = cfg['coef']*math.sqrt(2*RRAT*PHI*Dy(cfg, 28.)*Ek2)           # coef (eps Theta)_min^{-1/2} sqrt(F_max) <= coef sqrt(2 R D(28) E)
    else: pre = cfg['coef']*28.**2*math.sqrt(eps28*Gw*2*cfg['EA'][2])
    s33 = pre*(Fj/math.sqrt(math.log(2)) + math.sqrt(2)*Fj**0.75*JG**0.25)
    return (0., FOOTW[1]*Fj, FOOTW[2]*Fj + s33), s33

def evol(cfg, yc, mu, nu, a, IWin, B_prev=None):
    """level-j right-hand-side terms (time-integrated) for a centre at y_c (late: y_c = 56 at tau_rho; present: y_c at tau_1)"""
    EA, gp, K, KS = cfg['EA'], cfg['gp'], cfg['K'], cfg['KS']
    Ek = [(math.sqrt(EA[k]) + math.sqrt(dU[k]))**2 for k in range(3)]
    Gw = Gwof(cfg)
    if not cfg['joint']:
        X = {1: eps28*Ek[1], 2: gp[0]*EA[0]/(2*amin_t1) + eps28*p**4*dU[0], 3: gp[1]*EA[1]/(2*amin_t1) + eps28*p**6*dU[0]}
        F2 = eps28*Gw*IWin*2*EA[2]; JG = 2*(gp[1]*EA[1]/(2*amin_t1))*IWin*Gw
        ft, s33 = foot_terms(cfg, F2, JG, Ek[2], Gw); out = []
        for j in range(3):
            br = sum(K[j][k]*(eps28*Ek[k] if k <= j else X[k]) for k in K[j])
            out.append(dict(inflow=1.70*IWin*Gw*Ek[j], ramp=Gw*IWin*2*KS*br, foot=ft[j], outer=DT6[j]))
        return out, dict(Gw=Gw, s33=s33)
    # ---- joint sweeps with pointwise D (U_ramp.tex Prop. 1) ----
    Dr = (1 + gD(cfg)/yR**4)*PHI if cfg['pointwise'] else Gw*np.ones_like(yR)
    S_in = sweep(cp*deta2*Dr*yR/28, yc, mu, nu, a)                    # inflow sweep
    S_E = sweep((28/yR)**3*Dr, yc, mu, nu, a)                          # eps-weighted energy data sweep (per eps28)
    S_X = supR(Dr*yR/28, yc, mu, nu, a)                                # sup_{tau, ramp} Theta/omega (dissipation data)
    D = cfg['D']
    # dissipation-data integrals on the ramp (Minkowski with the delta-U part, swept)
    def Xint(k):        # bound on iint_ramp eps Theta z_k^2 for k = 3 (from D_1) or k = 2 (from D_0)
        Dk = D[k - 2]
        return (math.sqrt(S_X*Dk) + math.sqrt(eps28*S_E*p**(2*k)*dU[0]))**2
    Xi = {1: eps28*S_E*Ek[1], 2: Xint(2), 3: Xint(3)}
    if cfg['energy12']: Xi[2] = min(Xi[2], eps28*S_E*Ek[2])            # z_2 on the ramp: (A)'s level-2 ENERGY (sup in time), swept
    if not cfg['honestX']:                                             # rate model (t1), but joint sweeps
        Xi = {1: eps28*S_E*Ek[1], 2: (gp[0]*EA[0]/(2*amin_t1) + eps28*p**4*dU[0])*IWin*2*Gw, 3: (gp[1]*EA[1]/(2*amin_t1) + eps28*p**6*dU[0])*IWin*2*Gw}
    F2 = eps28*S_E*Ek[2]                                              # J_F = iint_ramp eps Theta z_2^2
    JG = 2*(Xi[3] + (LRAMP/2)**2*F2)                                   # J_G = iint |((eps Theta)^{1/2} z_2)'|^2 / ... (U_ramp.tex eq. (S33))
    ft, s33 = foot_terms(cfg, F2, JG, Ek[2], Gw)
    out = []
    g = gD(cfg); G = cfg['G22']
    for j in range(3):
        brE = sum(K[j][k]*eps28*S_E*Ek[k] for k in K[j] if k <= j)
        brX = sum(K[j][k]*Xi[k] for k in K[j] if k > j)
        if cfg.get('ibp2') and j == 2:
            # U_ramp.tex Lemma Y: mixed pairs (a,3), a < 3, by exact IBP to energy data; top pair (3,3) by X (D_1) or by Young (delta), whichever is smaller [est]
            TP = J['toppair']['levels']['2']; ftop = TP['f_top']
            xtop = K[2][3]*ftop*Xi[3]
            ytop = eps28*TP[f"delta={cfg['ibp2']}"]['sweep_Y']*Ek[2] if cfg['ibp2'] in (0.1, 0.25, 0.5) else math.inf
            brX = min(xtop, ytop) + TP['mixed_ibp_coef']*eps28*S_E*Ek[2]
            d_ibp = dict(top_X=xtop, top_Y=ytop, mixed=TP['mixed_ibp_coef']*eps28*S_E*Ek[2])
        d = dict(inflow=S_in*Ek[j], ramp=KS*brE, rampX=KS*brX, foot=ft[j], outer=DT6[j]*(1 + g/190.**4)*PHI if cfg['pointwise'] else DT6[j])
        # ramp share of the rate as data (design data56): excess eps[Ghat_jj D - (1-4a)g]_+ on y < 56, carried by eta^2 W phi (U_ramp.tex eq. (rs))
        if cfg['design'] == 'data56':
            Gjj = cfg['Ghat'][j][j] if cfg['Ghat'] else G
            exc = np.maximum(Gjj*(1 + g/yR**4) - (1 - 4*abar)*g, 0.)/yR**4
            d['rs'] = sweep(exc*eta2*PHI*yR/28, yc, mu, nu, a)*Ek[j]
        if cfg['coupling'] and cfg['Ghat'] and j >= 1 and B_prev is not None:
            # Ghat_jk, k < j: ramp part as data, bulk part through the level-k damping surplus (U_ramp.tex Lemma D)
            c = 0.
            for k in range(j):
                Gjk = cfg['Ghat'][j][k]; Gkk = cfg['Ghat'][k][k]; Dm = Dy(cfg, 28.)*PHI
                surplus = (1 - 4*abar)*g - Gkk*Dm
                assert surplus > 0
                c += Gjk*(eps28*S_E*Ek[k] + Dm*B_prev[k]/surplus)
            d['cpl'] = c
        if cfg['bterm'] and j == 2 and B_prev is not None:
            # b-term lower-order coupling at level 2: bsup * iint e^{-s} Psi z_1^2 <= bsup B_1/(c- KB - bsup C_1)  (U_ramp.tex Lemma T(iii))
            d['b'] = bsup*B_prev[1]/(cm*KB - bsup*5.0)
        out.append(d)
    return out, dict(Gw=Gw, S_in=S_in, S_E=S_E, S_X=S_X, Xi=Xi, F2=F2, JG=JG, s33=s33)

def box(B, wmin, Uwin):
    Eu = [b/wmin for b in B]
    S0 = math.sqrt(Eu[0]/l + 2*math.sqrt(Eu[0]*Eu[1])); S1 = math.sqrt(Eu[1]/l + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(f=Uwin + S0, yfp=1.05*Uwin + S1, M0=0.11/(Uwin + S0), M1=0.52/(1.05*Uwin + S1),
                M0w=(0.15 if cfg_box_is_wide[0] else 0.11)/(Uwin + S0))
cfg_box_is_wide = [False]

def late(cfg, mu=4, nu=4, a=0.3):
    w = weight(mu, nu, a); Gw = Gwof(cfg); EA = cfg['EA']; q, Dd, CA = w['q'], w['Dd'], w['CA']
    if cfg['pointwise']:
        ylin = np.exp(np.linspace(math.log(28), math.log(350), 4001))
        CAD = float((Wf(mu, nu, a, np.log(ylin/56))*(1 + gD(cfg)/ylin**4)*ylin/28).max())*1.001*PHI
        D350 = (1 + gD(cfg)/350.**4)*PHI
    B = [0.]*3; ev = None; info = None
    for it in range(300):
        if cfg['pointwise']: Dl = [CAD*EA[j] + D350*(math.sqrt(q*B[j]) + math.sqrt(Dd[j]))**2 for j in range(3)]
        else: Dl = [Gw*(CA*EA[j] + (math.sqrt(q*B[j]) + math.sqrt(Dd[j]))**2) for j in range(3)]
        ev, info = evol(cfg, 56., mu, nu, a, w['IWin56'], B_prev=B if it > 0 else None)
        Bn = [Dl[j] + sum(ev[j].values()) for j in range(3)]
        if it > 2 and max(abs(Bn[j] - B[j])/max(Bn[j], 1e-300) for j in range(3)) < 1e-14: B = Bn; break
        B = Bn
    return dict(B=B, D=Dl, ev=ev, info=info, Gw=Gw, **box(B, w['wmin'], 1.3056*A1*84**p))

YC = (56., 60., 70., 84., 100., 126.6, 150., 200., 253.)
def present(cfg, mu=4, nu=4, a=0.3):
    w = weight(mu, nu, a); Gw = Gwof(cfg); rows = []
    for yc in YC:
        if cfg['pointwise']:
            ylin = np.exp(np.linspace(math.log(28), math.log(60), 2001))
            dH2 = float((Wf(mu, nu, a, np.log(ylin/yc))*(1 + gD(cfg)/ylin**4)*ylin/28).max())*1.001*PHI
            dH3 = sum(supW(mu, nu, a, math.log(yw/yc), math.log(1.5*yw/yc))*(1.5*yw/28)*(1 + gD(cfg)/yw**4)*PHI for yw in (50*1.5**k for k in range(5)))
            pre = [dH2*H2[j] + dH3*H3[j] for j in range(3)]
        else:
            dH2 = supW(mu, nu, a, math.log(28/yc), math.log(60/yc), fac=lambda x, yc=yc: yc*np.exp(x)/28)
            dH3 = sum(supW(mu, nu, a, math.log(yw/yc), math.log(1.5*yw/yc))*(1.5*yw/28) for yw in (50*1.5**k for k in range(5)))
            pre = [Gw*(dH2*H2[j] + dH3*H3[j]) for j in range(3)]
        B = [0.]*3
        for it in range(100):
            ev, info = evol(cfg, yc, mu, nu, a, IW_in(mu, nu, a, yc), B_prev=B if it > 0 else None)
            Bn = [pre[j] + sum(ev[j].values()) for j in range(3)]
            if it > 2 and max(abs(Bn[j] - B[j])/max(Bn[j], 1e-300) for j in range(3)) < 1e-14: B = Bn; break
            B = Bn
        rows.append(dict(yc=yc, B=B, ev=ev, **box(B, w['wmin'], 1.3056*A1*(1.5*yc)**p)))
    return rows
def present_worst(rows):
    r = [x for x in rows if x['yc'] <= 126.7]
    return dict(M0=min(x['M0'] for x in r), M1=min(x['M1'] for x in r), M0w=min(x['M0w'] for x in r),
                yc_M0=min(r, key=lambda x: x['M0'])['yc'], yc_M1=min(r, key=lambda x: x['M1'])['yc'])

# ================= SELF-CHECK against t1_checkpoint.json =================
res = dict(selfcheck={}, rows=[], args={})
res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), numpy=np.__version__,
                   env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                   inputs={k: dict(path=v, sha256=hashlib.sha256(open(v, 'rb').read()).hexdigest()) for k, v in INP.items()},
                   script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
print('SELF-CHECK: all new ingredients off -> t1_checkpoint.json rows (weight (4,4,0.3), G22 = WORST, with (2,2,3) placeholder)')
ok = True
for KS in (1, 2):
    ref = [r for r in J['t1chk']['rows'] if r['weight'] == [4, 4, 0.3] and r['G22_name'] == 'WORST' and r['KS'] == KS and r['cubic223']][0]
    lt = late(cfg_checkpoint(KS=KS)); pw = present_worst(present(cfg_checkpoint(KS=KS)))
    devs = dict(late_M1=abs(lt['M1'] - ref['late_launch']['M1'])/ref['late_launch']['M1'], late_M0=abs(lt['M0'] - ref['late_launch']['M0'])/ref['late_launch']['M0'],
                late_B2=abs(lt['B'][2] - ref['late_launch']['B'][2])/ref['late_launch']['B'][2],
                pres_M1=abs(pw['M1'] - ref['present_launch']['M1'])/ref['present_launch']['M1'], pres_M0=abs(pw['M0'] - ref['present_launch']['M0'])/ref['present_launch']['M0'])
    ok &= max(devs.values()) < 1e-9
    print(f"  KS={KS}: late M1 x{lt['M1']:.6f} (ref x{ref['late_launch']['M1']:.6f}), M0 x{lt['M0']:.6f} (ref x{ref['late_launch']['M0']:.6f}), B2 {lt['B'][2]:.6f};"
          f" present M1 x{pw['M1']:.6f} (ref x{ref['present_launch']['M1']:.6f}), M0 x{pw['M0']:.6f}; max rel dev {max(devs.values()):.1e}")
    res['selfcheck'][f'KS={KS}'] = dict(late_M1=lt['M1'], late_M0=lt['M0'], ref_late_M1=ref['late_launch']['M1'], devs=devs)
res['selfcheck']['passed'] = bool(ok); print('  SELF-CHECK', 'PASSED' if ok else 'FAILED'); assert ok

# ================= runs =================
inp0 = inputs_json(J['agentA_inputs']) if cli.inputs else inputs_default()
print(f"\nINPUTS: {inp0['src']}\n  E^A = {['%.4g' % x for x in inp0['EA']]}, gamma' = {['%.4g' % x for x in inp0['gp']]}, D_k = {['%.4g' % x for x in inp0['D']]}")
print(f"  K = certified (ramp_K_cert['4.0']['K']); footw = {FOOTW}; S33 coef = {COEF_S33:.3f} + (2,2,3) S33-form {COEF_223:.3f} [est C_d = {CD223:.2f}] (placeholder was {COEF_223_PLACEHOLDER:.3f})")
w = weight(4, 4, 0.3)
def describe(tag, cfg):
    cfg_box_is_wide[0] = cfg.get('box', '').startswith('widened (0.15')
    lt = late(cfg); pw = present_worst(present(cfg)); e2 = lt['ev'][2]; e1 = lt['ev'][1]
    parts2 = dict(launch=lt['D'][2], **{k: v for k, v in e2.items()}); parts1 = dict(launch=lt['D'][1], **{k: v for k, v in e1.items()})
    row = dict(tag=tag, name=cfg['name'], design=cfg['design'], G22=cfg['G22'], g=gD(cfg), Gw=lt['Gw'], KS=cfg['KS'],
               late=dict(M0=lt['M0'], M0w=lt['M0w'], M1=lt['M1'], B=lt['B'], parts2=parts2, parts1=parts1, info={k: v for k, v in lt['info'].items()}),
               present=pw, flags={k: cfg[k] for k in ('pointwise', 'joint', 'honestX', 'energy12', 'coupling', 'bterm', 'honestS33')}, D=cfg['D'], EA=cfg['EA'])
    res['rows'].append(row)
    M0lab = 'M0(0.15)' if cfg_box_is_wide[0] else 'M0(0.11)'
    print(f"  {tag:58s} g {gD(cfg):.3g}  LATE: M1 x{lt['M1']:.3f}  {M0lab} x{lt['M0w']:.2f}  B=({lt['B'][0]:.2e},{lt['B'][1]:.2e},{lt['B'][2]:.4f}) | PRESENT(yc<=126.7): M1 x{pw['M1']:.2f}  {M0lab} x{pw['M0w']:.2f} (yc {pw['yc_M0']})")
    print('      level-2 parts: ' + ', '.join(f'{k} {v:.4g}' for k, v in parts2.items()))
    print('      level-1 parts: ' + ', '.join(f'{k} {v:.3g}' for k, v in parts1.items()))
    return lt, pw

print('\n== ladder (weight (4,4,0.3), KS=1): switch the new ingredients on one at a time  [est]')
base = cfg_checkpoint(); describe('0 CHECKPOINT (t1_checkpoint late row)', base)
c = dict(base, EA=inp0['EA'], D=inp0['D']); describe('1 + E^A repair factors', c)
c = dict(c, joint=True); describe('2 + joint (y,tau) sweeps, rate-model X', c)
c = dict(c, honestX=True); describe('3 + honest X (D_k), z_2 by dissipation', c)
c = dict(c, energy12=True); describe('4 + z_2 on ramp by (A) energy where smaller', c)
c = dict(c, honestS33=True, coef=COEF_S33 + COEF_223); describe('5 + honest S33 normalisation, (2,2,3) in S33 form', c)
c = dict(c, pointwise=True, coupling=True, bterm=True, Ghat={j: {k: (J['design']['cases']['crude (0.11,0.52), U .09, M0z .16']['Ghat_full'][str(j)][str(k)]) for k in range(j + 1)} for j in range(3)}); describe('6 + pointwise D, couplings, b-term (global28, worst G22)', c)
for box_ in ('crude (0.11,0.52), U .09, M0z .16', 'widened (0.15,0.52), U .09, M0z .16'):
    describe(f'7 data56, {box_}', cfg_new(inp0, 'data56', box_))
print('\n== headline: design data56, widened box (0.15,0.52), KS = 1 / 2, default inputs  [est]')
HL = {}
for KS in (1, 2):
    lt, pw = describe(f'HEADLINE KS={KS}', cfg_new(inp0, 'data56', KS=KS)); HL[KS] = dict(late_M1=lt['M1'], late_M0w=lt['M0w'], pres_M1=pw['M1'], pres_M0w=pw['M0w'])
print('\n== variant (U_ramp.tex Lemma Y, [est] coefficients from top_pair_est.py): level-2 mixed pairs by exact IBP, top pair by min(X, Young)')
for dl in (0.25, 0.5):
    # Young against Lemma U's own level-2 dissipation: B_c 1.5 -> 1.5 - delta; Ghat_22 rescaled conservatively by 1.5/(1.5 - delta) (whole Ghat_22)
    cv = cfg_new(inp0, 'data56', ibp2=dl); cv['G22'] = cv['G22']*1.5/(1.5 - dl); cv['Ghat'] = {jj: dict(v) for jj, v in cv['Ghat'].items()}; cv['Ghat'][2][2] = cv['G22']
    describe(f'VARIANT ibp2 delta={dl} (Ghat22 x{1.5/(1.5 - dl):.2f})', cv)
    cv['KS'] = 2; describe(f'VARIANT ibp2 delta={dl}, KS=2', cv)
print('\n== sensitivity of the headline (KS=1) to the D_k export  [est]')
for lab, kap, cD in (('(A) as is: kappa = (0.3, 0.2), cD = 2', (0.3, 0.2, 0.), 2.), ('kappa = 1.25, cD = 2 (default)', (1.25, 1.25, 0.), 2.),
                     ('kappa = 1.25, cD = 1 (additive bound, U_ramp.tex Lemma D(i))', (1.25, 1.25, 0.), 1.)):
    describe(f'SENS {lab}', cfg_new(inputs_default(kappa=kap, cD=cD) if not cli.inputs else dict(inp0), 'data56'))
print('\n== sensitivity to E^A_2 (repair factor on level 2), headline config  [est]')
for f2 in (1.0, 1.25, 1.65, 2.5):
    describe(f'SENS E^A x(1.41,1.13,{f2})', cfg_new(inputs_default(repair=(1.41, 1.13, f2)), 'data56'))
# honest X-term factor (Fable's normalisation: model X = gamma' E/(2 a_min) swept by 2 IW_in; honest = sup(W y/28) D_{j-1})
SXnoD = supR(yR/28, 56.); f_X = {k: SXnoD*inp0['D'][k - 2]/(2*w['IWin56']*inp0['gp'][k - 2]*inp0['EA'][k - 2]/(2*amin_t1)) for k in (2, 3)}
print(f"\nHONEST X-TERM FACTOR (vs t1 rate model, no D, no Gw): level-1 data z_2: x{f_X[2]:.2f}; level-2 data z_3: x{f_X[3]:.2f}  (= {SXnoD:.4f} cD/(kappa IW_in(56) gamma'), sup W y/28 = {SXnoD:.4f}, IW_in = {w['IWin56']:.4f})")
res['honest_X_factor'] = dict(f_X={str(k): v for k, v in f_X.items()}, supWy28=SXnoD, IWin56=w['IWin56'])
res['headline'] = HL
res['wall_s'] = time.time() - T0
fn = os.path.join(HERE, f'U_model2{"_" + cli.tag if cli.tag else ""}.json')
json.dump(res, open(fn + '.tmp', 'w'), indent=1, default=float); os.replace(fn + '.tmp', fn)
print(f'\nwall {res["wall_s"]:.1f} s -> {fn}')
