#!/usr/bin/env python3
"""Task 4, U6: certified (arb) assembly of Lemma U (Prop. U with Prop. 8 data term, Remark 2 launch data, Lemma N, Prop. 16 (K2), (K3), loop gain).

  PYTHONHASHSEED=0 nice -n 19 python3 -B assemble_U.py [--A PATH] [--setting design] [--zone3 PATH] [--case primary] [--h2f-setting NAME] [--out NAME]
      -> out/U_cert.json (primary) / out/U_cert_<case>.json, and out/U_cert_for_zone3.json (the U-fine Bhat in the key zone3.py --ucert reads)

Inputs (paths on the command line; sha256 in the args block):
  --A      t2_A/out/A_interface.json (task 2).  If absent: a PROVISIONAL interface built below from phase2a's certified export and the r2 expectations
           (labels est / hyp); the output is then marked provisional: true.
  --zone3  t5_zone3/out/zone3_interface.json (task 5): DeltaR, DeltaR_b, n_III, box_III_380.b1.
  u1_ghat.json (Ghat, g, design), u2_K[_case].json (K, footw, c33, c223, L_max, R), u34_<case>.json (sweeps, launch constants, present cells).
Formulas (STATEMENTS_r2.tex; U_ramp.tex Prop. U; all in arb, outward):
  E^k_j = (sqrt(EA_j) + sqrt(dU_j))^2, dU_j = (p^j c_dU 56^p)^2 ln 2;  D_k = S_X C0_k + C1_k Sigma_X (Prop. 8, k = 0, 1);  J = eps_28 Sigma_E;
  X_1 = J E^k_1;  X_2 = min(J E^k_2, (sqrt D_0 + sqrt(J p^4 dU_0))^2);  X_3 = (sqrt D_1 + sqrt(J p^6 dU_0))^2;
  B_j = ( sqrt( L_j + I_j + R^E_j + R^X_j + R^rs_j + F_j + C_j + Bb_j + DR_j ) + Phi_j )^2 with
    L_j  late: (sqrt(C_A^D EA_j) + sqrt(Delta_in_j))^2 + q^D beta Bhat_j;  present: c_2 H2_j + c_3 H3_j
    I_j  = c^+ Sigma_in E^k_j ;  R^E_j = sum_{k<=j} K_jk J E^k_k ;  R^X_j = sum_{k>j} K_jk X_k ;  R^rs_2 = Sigma_rs2 E^k_2
    F_1  = footw_1 J_F;  F_2 = footw_2 J_F + S33 + T223,  J_F = J E^k_2,  J_G = 2 (X_3 + L_max^2 J_F/4),  P2 = 2 R D(28) phi(28) E^k_2
         S33 = c33 sqrt(P2) (J_F/sqrt(ln2) + sqrt2 J_F^{3/4} J_G^{1/4});  T223 = c223 sqrt(P2 (J_F/ln2 + 2 sqrt(J_F J_G))) sqrt(X_3)  (Cauchy-Schwarz)
    C_j  = sum_{k<j} Ghat_jk [J E^k_k + D(28) phi(28) B_k/((1 - 4abar) g - Ghat_kk D(28) phi(28))]
    Bb_2 = bbar B_1/(c^- K_b - bbar (5 + K_b/28));  DR_j = DeltaR_j D(253) phi(253) + DeltaR_b_j ;  Phi_j from u34.
  Lemma N: N2_k(beta) = 2 (sqrt(beta Bhat_k/w_min) + p^k c_dU 700^p sqrt(ln1.5))^2 + n_III_k;  EA_j(beta) = EA0_j + sum_k kappa_o_jk N2_k/gp_j;
           C0_k = (EA_k + Dsrc_k)/(r_k a_c);  C1_k(beta) = C1_0_k + sum_i kappa_o_ki N2_i/(r_k a_c).
  Closure (Prop. 16): Bhat := (1 + 1e-6) x fixed point of Bhat = sup_centres U(2; Bhat), EA := A(2; Bhat) (rounded up); then asserted in arb:
    (K1, A part) A_j(2) <= EA_j;  (K2) U_j(2) <= Bhat_j for the late centre and every present cell y_c in [56, 253];
    (K3) window Sobolev: Ubar + S_0(B) < M0, p Ubar + S_1(B) < M1, S_0(B) < M0z for the late centre and every present cell y_c in [56, 168.67],
         S_0(B)^2 = B_0/(w_min l) + 2 sqrt(B_0 B_1)/w_min, S_1(B)^2 = B_1/(w_min l) + 2 sqrt(B_1 B_2)/w_min (B the centre's own bound).
    Also the hypothesis-level |z_1| bound used by u1 (M1z = 0.534): S_1(2 Bhat) <= M1z (zones I/II), b_1 + p U_1(379.5) <= M1z (zone III).
  Loop gain (Def. 17): G_j = (U_j(2) - U_j(0))/Bhat_j, split: neighbour (beta in the launch term only) and through (A) (beta in N2 only)."""
import sys, os, json, math, time, argparse
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
from flint import arb, ctx
from t4lib import A, Q_, fup, flo, num
ctx.prec = 128
ap = argparse.ArgumentParser()
ap.add_argument('--A', default=L.T2 + '/out/A_interface.json'); ap.add_argument('--setting', default='tight')
ap.add_argument('--zone3', default=L.T5 + '/out/zone3_interface.json'); ap.add_argument('--case', default='primary')
ap.add_argument('--u1', default=os.path.join(L.OUTD, 'u1_ghat.json')); ap.add_argument('--out', default=None)
ap.add_argument('--stressC1', default='1'); ap.add_argument('--EA2', default=None)
args = ap.parse_args()
T0 = time.time()
tag = '' if args.case == 'primary' else '_' + args.case
U2P = os.path.join(L.OUTD, f'u2_K{tag}.json'); U34P = os.path.join(L.OUTD, f'u34_{args.case}.json')
U1 = json.load(open(args.u1)); CASE = U1['cases'][args.case]; U2 = json.load(open(U2P)); U34 = json.load(open(U34P)); Z3 = json.load(open(args.zone3))
assert U2['case'] == args.case and U34['case'] == args.case and U34['g'] == CASE['design']['g']
G_LO, G_HI = (Q_(x) for x in U2['g_range']); gq = Q_(CASE['design']['g']); assert G_LO <= gq <= G_HI, 'K not certified at this g'
def a_(x):
    if isinstance(x, dict): return arb(x['double'])
    if isinstance(x, str): return A(Q_(x))
    return arb(x)
lo_ = a_

# ---------------------------------------------------------------- constants
VS = L.verifier_sups(); ABAR, BBAR = VS['abar'], VS['bbar']
CM = arb(1)/4 - ABAR; CP = arb(1)/4 + ABAR + BBAR/28
P = A(Q_(10246, 10000)); ALPHA1 = A(Q_(15, 100000)); KB = A(Q_(1, 1000)); LN2 = arb(2).log(); LL = (arb(3)/2).log()
EPS28 = arb(28)**-4; g = A(gq)
def Dy(y): return 1 + g/arb(y)**4
def phy(y): return (KB/arb(y)).exp()
D28p = Dy(28)*phy(28)
U1f = lambda y: A(Q_(13056, 10000))*ALPHA1*(arb(y).log()*P).exp()
BOX = CASE['box']; M0B, M1B, M0ZB, M1ZB = (A(Q_(BOX[k])) for k in ('M0', 'M1', 'M0z', 'M1z'))
GH = [[a_(CASE['Ghat'][j][k]) for k in range(j + 1)] for j in range(3)]
KA = {int(j): {int(k): arb(v) for k, v in d.items()} for j, d in U2['K_aug'].items()}
FOOTW = [a_(U2['foot'][str(j)]['footw']) for j in range(3)]
C33, C223, LMAX, RR = a_(U2['c33']), a_(U2['c223']), a_(U2['L_max']), a_(U2['R_ratio'])
LA = U34['launch']; CAD, QD, WMIN, CDU = a_(LA['CA_D']), a_(LA['qD']), lo_(LA['wmin']), a_(LA['c_dU'])
DIN = [a_(x) for x in LA['Delta_in']]; PHI_LATE = [a_(x) for x in LA['Phi_late']]
SW = {k: {kk: (a_(v) if isinstance(v, dict) else v) for kk, v in d.items()} for k, d in U34['sweeps'].items()}
DR = [a_(Z3['DeltaR'][j])*Dy(253)*phy(253) + a_(Z3['DeltaR_b'][j]) for j in range(3)]
NIII = [a_(Z3['n_III'][k]) for k in range(3)]
H2 = [A(Q_('1e-8')), A(Q_('1e-6')), A(Q_('1e-3'))]; H3 = [A(Q_('2e-7')), A(Q_('1e-5')), A(Q_('1e-5'))]
dU = [(P**j*CDU*(arb(56).log()*P).exp())**2*LN2 for j in range(3)]
SURP = [(1 - 4*ABAR)*g - GH[k][k]*D28p for k in range(3)]
assert SURP[0] > 0 and SURP[1] > 0
BDEN = CM*KB - BBAR*(5 + KB/28); assert (CM*KB - 2*BBAR*(8 + KB/28)) > 0 and BDEN > 0

# ---------------------------------------------------------------- (A) interface
PROVISIONAL = not os.path.exists(args.A)
if not PROVISIONAL:
    AI = json.load(open(args.A)); res = AI['res'] if 'res' in AI else AI
    S = res.get(args.setting) or res.get('tight'); assert S is not None, ('setting not found', args.setting, list(res))
    EA0 = [a_(x) for x in S['EA0']]; KO = [[a_(x) for x in row] for row in S['kappa_o']]; GP = [lo_(x) for x in S['gp']]
    DSRC = [a_(x) for x in S['Dsrc']]; RX = [lo_(x) for x in S['r_export']]; AC = lo_(S['amin_c']); C10 = [a_(x) for x in S['C1_0']]
    BETA_A = a_(AI.get('beta_A', S.get('beta_A', '5/4')))
    A_FLUX = {k: (v['double'] if isinstance(v, dict) and 'double' in v else v) for k, v in S.get('core', {}).get('margins', {}).items()}
    N2REF_A = [A(Q_(x)) for x in S['N2_ref']]
    A_SRC = dict(path=args.A, setting=args.setting, label='cert (task 2)')
else:
    # provisional interface [est / hyp]: phase2a certified export (tightened H2) + r2 expected E (R2_DECISIONS.md); outer source at the placeholder rates
    EXP = json.load(open(L.T2 + '/export_for_U.json'))
    GP = [arb(x) for x in EXP['gp']]; RX = [A(Q_('1.249'))]*2; AC = (1 + A(Q_(11, 100))**2)**-2
    EAexp = [A(Q_('7.0e-7')), A(Q_('1.0e-4')), A(Q_('0.0726'))]                    # r2 expected Ehat [est]
    N2REF = [A(Q_('0.1')), A(Q_('0.2')), A(Q_('5.0'))]; RATE1 = [A(Q_('6.165e-9')), A(Q_('7.120e-9')), A(Q_('1.2614e-7'))]   # phase2a outer rates [hyp]
    KO = [[RATE1[j]/N2REF[j] if k == j else arb(0) for k in range(3)] for j in range(3)]
    EA0 = [EAexp[j] - RATE1[j]/GP[j] for j in range(3)]
    C0x = [arb(x) for x in EXP['D_from_any_tau_a_with_E_tau_a_le_EA'][:2]]; Ex = [arb(x) for x in EXP['EA'][:2]]
    DSRC = [C0x[k]*RX[k]*AC - Ex[k] for k in range(2)]
    C10 = [arb(EXP['D_slope'][k])*A(Q_(args.stressC1)) - RATE1[k]/(RX[k]*AC) for k in range(2)]
    BETA_A = A(Q_(5, 4)); A_FLUX = None; N2REF_A = N2REF
    A_SRC = dict(path=L.T2 + '/export_for_U.json', label='PROVISIONAL [est/hyp]: EA0 = r2 expected Ehat (7.0e-7, 1.0e-4, 0.0726) minus placeholder outer rate/gp; '
                 'kappa_o = diag(phase2a placeholder rates / N2_ref); Dsrc and C1_0 from phase2a export (C0 r a_c - E, D_slope - rate/(r a_c))')
if args.EA2: EA0[2] = A(Q_(args.EA2)) - sum((KO[2][k]*A(Q_('5.0')) for k in range(3)), arb(0))/GP[2]

def N2(beta, Bh):
    return [2*((arb(beta)*Bh[k]/WMIN).sqrt() + P**k*CDU*(arb(700).log()*P).exp()*LL.sqrt())**2 + NIII[k] for k in range(3)]
def Aside(beta, Bh):
    n2 = N2(beta, Bh)
    EA = [EA0[j] + sum((KO[j][k]*n2[k] for k in range(3)), arb(0))/GP[j] for j in range(3)]
    C0 = [(EA[k] + DSRC[k])/(RX[k]*AC) for k in range(2)]
    C1 = [C10[k] + sum((KO[k][i]*n2[i] for i in range(3)), arb(0))/(RX[k]*AC) for k in range(2)]
    return dict(N2=n2, EA=EA, C0=C0, C1=C1)

def S01(B):
    S0 = (B[0]/(WMIN*LL) + 2*(B[0]*B[1]).sqrt()/WMIN).sqrt(); S1 = (B[1]/(WMIN*LL) + 2*(B[1]*B[2]).sqrt()/WMIN).sqrt(); return S0, S1

def centre(kind, sw, launch, phi, beta, Bh, Ad, parts=False):
    """kind 'late' or 'present'; sw: sweeps at y_c; launch: callable(EA, beta, Bh) -> [3]; Ad: Aside(...) dict"""
    EA, C0, C1 = Ad['EA'], Ad['C0'], Ad['C1']
    J = EPS28*sw['Sigma_E']
    Ek = [(EA[j].sqrt() + dU[j].sqrt())**2 for j in range(3)]
    Dk = [sw['S_X']*C0[k] + C1[k]*sw['Sigma_X'] for k in range(2)]
    X = {1: J*Ek[1]}
    x2a = J*Ek[2]; x2b = (Dk[0].sqrt() + (J*P**4*dU[0]).sqrt())**2
    X[2] = x2a if x2a.upper() <= x2b.upper() else x2b
    X[3] = (Dk[1].sqrt() + (J*P**6*dU[0]).sqrt())**2
    JF = J*Ek[2]; JG = 2*(X[3] + LMAX**2*JF/4); P2 = 2*RR*D28p*Ek[2]
    S33 = C33*P2.sqrt()*(JF/LN2.sqrt() + arb(2).sqrt()*JF**arb(0.75)*JG**arb(0.25))
    T223 = C223*(P2*(JF/LN2 + 2*(JF*JG).sqrt())).sqrt()*X[3].sqrt()
    Ln = launch(EA, beta, Bh); B = []; PT = []
    for j in range(3):
        d = dict(launch=Ln[j], inflow=CP*sw['Sigma_in']*Ek[j],
                 rampE=sum((KA[j][k]*J*Ek[k] for k in KA[j] if k <= j), arb(0)),
                 rampX=sum((KA[j][k]*X[k] for k in KA[j] if k > j), arb(0)),
                 rs=(sw['Sigma_rs2']*Ek[2] if j == 2 else arb(0)),
                 foot=(FOOTW[j]*JF if j >= 1 else arb(0)), S33=(S33 if j == 2 else arb(0)), T223=(T223 if j == 2 else arb(0)),
                 coupling=sum((GH[j][k]*(J*Ek[k] + D28p*B[k]/SURP[k]) for k in range(j)), arb(0)),
                 b=(BBAR*B[1]/BDEN if j == 2 else arb(0)), outer=DR[j])
        tot = sum(d.values(), arb(0))
        Bj = (tot.sqrt() + phi[j])**2
        B.append(arb(fup(Bj)))                       # exact upper bound, keeps radii from growing
        d['forcing_Phi'] = phi[j]; PT.append(d)
    return (B, PT, dict(X=X, D=Dk, Ek=Ek, JF=JF, JG=JG)) if parts else B

def late_launch(EA, beta, Bh): return [(CAD*EA[j]).sqrt().__add__(DIN[j].sqrt())**2 + QD*arb(beta)*Bh[j] for j in range(3)]
SW56 = SW['56']
PCELLS = U34['present_cells']
def pres_launch_fn(c):
    c2, c3 = a_(c['c2']), a_(c['c3'])
    return lambda EA, beta, Bh: [c2*H2[j] + c3*H3[j] for j in range(3)]
EAF = None
def all_centres(beta_nb, beta_A, Bh, parts=False, freeze=False):
    """returns dict name -> B (or (B, parts, info)); beta_nb multiplies the neighbour launch term, beta_A enters Lemma N;
    freeze: EA frozen at the certified EAF (>= A(2; Bhat)), C0 from it"""
    Ad = Aside(beta_A, Bh); out = {}
    if freeze: Ad['EA'] = EAF; Ad['C0'] = [(EAF[k] + DSRC[k])/(RX[k]*AC) for k in range(2)]
    out['late'] = centre('late', SW56, late_launch, PHI_LATE, beta_nb, Bh, Ad, parts)
    for c in PCELLS:
        out['present ' + '-'.join(c['yc'])] = centre('present', SW[c['sweeps_at']], pres_launch_fn(c), [a_(x) for x in c['Phi']], beta_nb, Bh, Ad, parts)
    return out, Ad
def supB(d): return [max((v[j] if isinstance(v, list) else v[0][j]) for v in d.values()) for j in range(3)]

# ---------------------------------------------------------------- fixed point
Bh = [arb(0)]*3
for it in range(60):
    d, Ad = all_centres(2, 2, Bh); Bn = supB(d)
    if it > 3 and max(fup(abs(Bn[j] - Bh[j])/Bn[j]) for j in range(3)) < 1e-12: Bh = Bn; break
    Bh = Bn
NIT = it
BH = [arb(fup(Bh[j]*(1 + arb(10)**-6))) for j in range(3)]
AdF = Aside(2, BH); EAF = [arb(fup(x)) for x in AdF['EA']]
def run_all(beta_nb, beta_A, Bh, parts=False, freezeEA=True):
    return all_centres(beta_nb, beta_A, Bh, parts, freeze=(freezeEA and beta_A == 2))
DF, AdK = run_all(2, 2, BH, parts=True)
K1_ok = all(fup(AdF['EA'][j]) <= fup(EAF[j]) for j in range(3))
K2_rows = {n: [fup(v[0][j]/BH[j]) for j in range(3)] for n, v in DF.items()}
K2_ok = all(v[0][j].upper() <= BH[j] for v in DF.values() for j in range(3))
assert K1_ok and K2_ok, (K1_ok, K2_rows)
# loop gain
D0, _ = run_all(0, 0, BH, freezeEA=False); Dnb, _ = run_all(2, 0, BH, freezeEA=False); DA, _ = run_all(0, 2, BH, freezeEA=False)
Uof = lambda d: supB(d)
U2v, U0v, Unb, UAv = supB({k: v[0] for k, v in DF.items()}), Uof(D0), Uof(Dnb), Uof(DA)
GAIN = [dict(total=fup((U2v[j] - U0v[j])/BH[j]), neighbour=fup((Unb[j] - U0v[j])/BH[j]), through_A=fup((UAv[j] - U0v[j])/BH[j])) for j in range(3)]

# ---------------------------------------------------------------- (K3) window Sobolev, margins
Blate = DF['late'][0]; S0l, S1l = S01(Blate); Ul = U1f(84)
late_M1 = M1B/(P*Ul + S1l); late_M0 = M0B/(Ul + S0l); late_M0z = M0ZB/S0l
pres = []
for c in PCELLS:
    a, b = Q_(c['yc'][0]), Q_(c['yc'][1])
    Bc = DF['present ' + '-'.join(c['yc'])][0]; S0, S1 = S01(Bc); Ub = a_(c['Ubar'])
    row = dict(yc=c['yc'], B=[fup(x) for x in Bc], S0=fup(S0), S1=fup(S1), Ubar=fup(Ub), M0=flo(M0B/(Ub + S0)), M1=flo(M1B/(P*Ub + S1)), M0z=flo(M0ZB/S0),
               in_zone_I_II=bool(b <= Q_(16867, 100)), ufine=bool(a >= 150 and b <= Q_(16867, 100)))
    pres.append(row)
PZ = [r for r in pres if r['in_zone_I_II']]
pmin_M0 = min(PZ, key=lambda r: r['M0']); pmin_M1 = min(PZ, key=lambda r: r['M1']); pmin_M0z = min(r['M0z'] for r in PZ)
K3_ok = late_M0.lower() > 1 and late_M1.lower() > 1 and late_M0z.lower() > 1 and all(r['M0'] > 1 and r['M1'] > 1 and r['M0z'] > 1 for r in PZ)
S0h, S1h = S01([2*x for x in BH])
zoneIII_z1 = a_(Z3['box_III_380']['b1']) + P*U1f(arb('379.5'))
M1z_ok = dict(zones_I_II=dict(S1_2Bhat=fup(S1h), M1z=float(M1ZB), ok=bool(S1h.upper() < M1ZB)),
              zone_III=dict(b1_plus_pU1=fup(zoneIII_z1), ok=bool(zoneIII_z1.upper() < M1ZB)), core_28_60='u1 assertion (0.52 + p U_1(60) + p c_dU 60^p < M1z)')
N2_ratio = [fup(AdF['N2'][k]/(4*N2REF_A[k])) for k in range(3)]
N2_ratio_ref = [fup(AdF['N2'][k]/N2REF_A[k]) for k in range(3)]
UF = [r for r in pres if r['ufine']]
BH_UFINE = [max(r['B'][j] for r in UF) for j in range(3)]

# ---------------------------------------------------------------- output
pl = DF['late'][1]; info = DF['late'][2]
parts2 = {k: fup(v) for k, v in pl[2].items()}
parts2_grouped = dict(launch=parts2['launch'], inflow=parts2['inflow'], rampE=parts2['rampE'], X_terms=parts2['rampX'], rs=parts2['rs'],
                      foot_S33=parts2['foot'] + parts2['S33'] + parts2['T223'], couplings=parts2['coupling'], b=parts2['b'], outer=parts2['outer'], forcing=parts2['forcing_Phi'])
slope_share = [fup(AdK['C1'][k]*SW56['Sigma_X']/info['D'][k]) for k in range(2)]
box_ag = dict(M0=float(M0B), M1=float(M1B), M0z=float(M0ZB))
LBL = 'est' if PROVISIONAL else 'cert'
def sw_key(k): return dict(yc56=num(SW56[k], LBL), cover_sup=num(SW56[k], LBL, rule='nonincreasing in y_c on [56, 253] (u34 doc): sup = value at 56'))
out = dict(provisional=PROVISIONAL, case=args.case, A_interface=A_SRC,
           box=dict(M0=BOX['M0'], M1=BOX['M1'], Umax=BOX['Umax'], M0z=BOX['M0z'], M1z=BOX['M1z'], design='data56'),
           Gamma=[U1['weight']['Theta_phi']['Gamma'][m] for m in range(3)], kappa=U1['weight']['Theta_phi']['kappa'], kappa2=U1['weight']['Theta_phi']['kappa2'],
           Ghat=CASE['Ghat'], g=CASE['design']['g'], D56=CASE['design']['D56'], D28=CASE['design']['D28'], theta28=CASE['design']['theta28'],
           rate_margin=CASE['design']['rate_margin'], R_lt2_margin=CASE['design']['R_lt2_margin'],
           K=U2['K_aug'], K_raw=U2['K'], footw=[U2['foot'][str(j)]['footw'] for j in range(3)], c33=U2['c33'], c223=U2['c223'], g_range=U2['g_range'],
           L_max=U2['L_max'], R_ratio=U2['R_ratio'], C_eta=U2['C_eta'],
           sweeps=dict(Sigma_in=sw_key('Sigma_in'), Sigma_E=sw_key('Sigma_E'), Sigma_rs=[dict(yc56=num(arb(0), 'cert', rule='(R<2)'), cover_sup=num(arb(0), 'cert')),
                       dict(yc56=num(arb(0), 'cert', rule='(R<2)'), cover_sup=num(arb(0), 'cert')), sw_key('Sigma_rs2')], S_X=sw_key('S_X'), Sigma_X=sw_key('Sigma_X'),
                       per_yc={k: {kk: (v['double'] if isinstance(v, dict) else v) for kk, v in d.items()} for k, d in U34['sweeps'].items()}),
           launch=dict(CA_D=LA['CA_D'], qD=LA['qD'], q=LA['q'], wmin=LA['wmin'], c2=[c['c2'] for c in PCELLS], c3=[c['c3'] for c in PCELLS], c23_cells=[c['yc'] for c in PCELLS],
                       Delta_in=LA['Delta_in'], c_dU=LA['c_dU'], Phi=LA['Phi_late'], sobolev=LA['sobolev']),
           supersolution=dict(EA=[num(x, LBL) for x in EAF], Bhat=[fup(x) for x in BH], Bhat_num=[num(x, LBL) for x in BH], N2=[num(x, LBL) for x in AdF['N2']],
                              C0=[num(x, LBL) for x in AdK['C0']], C1=[num(x, LBL) for x in AdK['C1']], beta=2, beta_A=num(BETA_A, 'cert' if not PROVISIONAL else 'hyp'),
                              K1_ok=K1_ok, K1_note='A_j(2) <= EA_j checked here; core-strictness and fluxes at beta_A EA are task 2 rows: ' + (json.dumps(A_FLUX) if A_FLUX else 'not available (provisional)'),
                              N2_over_4N2ref=N2_ratio, N2_over_N2ref=N2_ratio_ref, N2_le_4N2ref=bool(max(N2_ratio) <= 1),
                              K2_ok=K2_ok, K3_ok=bool(K3_ok), slack_K2=K2_rows, iterations=NIT, M1z_hyp_check=M1z_ok),
           loop_gain=GAIN,
           margins=dict(late_M1=num(late_M1, LBL, 'lo', against=box_ag['M1']), late_M0=num(late_M0, LBL, 'lo', against=box_ag['M0']), late_M0z=num(late_M0z, LBL, 'lo', against=box_ag['M0z']),
                        present_M1_min=dict(value=pmin_M1['M1'], at_yc=pmin_M1['yc'], against=box_ag['M1'], label=LBL),
                        present_M0_min=dict(value=pmin_M0['M0'], at_yc=pmin_M0['yc'], against=box_ag['M0'], label=LBL),
                        present_M0z_min=dict(value=pmin_M0z, against=box_ag['M0z'], label=LBL), present_rows=pres,
                        late_window=dict(B=[fup(x) for x in Blate], S0=fup(S0l), S1=fup(S1l), Ubar=fup(Ul))),
           parts_B2_late=parts2_grouped, parts_B2_late_detail=parts2, parts_B1_late={k: fup(v) for k, v in pl[1].items()}, parts_B0_late={k: fup(v) for k, v in pl[0].items()},
           data_terms_late=dict(X={str(k): fup(v) for k, v in info['X'].items()}, D=[fup(x) for x in info['D']], JF=fup(info['JF']), JG=fup(info['JG']),
                                slope_share_C1SigmaX_over_Dk=slope_share),
           S33_foot_over_Bhat2=fup((pl[2]['foot'] + pl[2]['S33'] + pl[2]['T223'])/BH[2]),
           ufine=dict(Bhat=BH_UFINE, rule='max over present cells y_c in [150, 168.67] of U_j(2) (conclusion level), cells ' + str([r['yc'] for r in UF]), label=LBL,
                      requirement_task5='Bhat_2 <= 0.106 with Bhat_0 <= 5e-6, Bhat_1 <= 1e-3', ok=bool(BH_UFINE[2] <= 0.106 and BH_UFINE[0] <= 5e-6 and BH_UFINE[1] <= 1e-3)),
           inputs=dict(H2=['1e-8', '1e-6', '1e-3'], H3=['2e-7', '1e-5', '1e-5'], alpha1='1.5e-4', abar=VS['strings']['sup_ab'], betabar=VS['strings']['sup_beta'],
                       DeltaR=[Z3['DeltaR'][j]['double'] for j in range(3)], n_III=[Z3['n_III'][k]['double'] for k in range(3)], zone3_provisional=Z3.get('provisional')))
ins = [args.u1, U2P, U34P, args.zone3, L.VERIF_ALPHA, os.path.abspath(__file__), os.path.join(HERE, 't4lib.py')] + ([args.A] if not PROVISIONAL else [L.T2 + '/export_for_U.json'])
out['args'] = L.args_block(ins, dict(setting=args.setting, stressC1=args.stressC1, EA2=args.EA2))
name = args.out or ('U_cert.json' if args.case == 'primary' else f'U_cert_{args.case}.json')
L.dump(out, os.path.join(L.OUTD, name))
if args.case == 'primary' and args.out is None:
    L.dump(dict(supersolution=dict(Bhat=BH_UFINE), note='U-fine Bhat (present centres rho in [150, 168.67] lambda_1, conclusion level) in the key zone3.py --ucert reads '
                '(supersolution.Bhat). The crude global Bhat is U_cert.json supersolution.Bhat.', provisional=PROVISIONAL, label=LBL, args=out['args']),
           os.path.join(L.OUTD, 'U_cert_for_zone3.json'))
print(f'[{args.case}] provisional={PROVISIONAL}  iterations {NIT}  Bhat {[fup(x) for x in BH]}  EA {[fup(x) for x in EAF]}')
print(f'  late: B {[fup(x) for x in Blate]}  M1 x{flo(late_M1):.4f} (vs {float(M1B)})  M0 x{flo(late_M0):.4f} (vs {float(M0B)})  M0z x{flo(late_M0z):.3f}')
print(f'  present (zones I/II): M1 min x{pmin_M1["M1"]:.4f} at {pmin_M1["yc"]}, M0 min x{pmin_M0["M0"]:.4f} at {pmin_M0["yc"]}, M0z min x{pmin_M0z:.3f}')
print(f'  loop gain {GAIN}')
print(f'  N2 (beta=2) {[fup(x) for x in AdF["N2"]]}  /N2_ref {N2_ratio_ref}  /(4 N2_ref) {N2_ratio}')
print(f'  U-fine Bhat {BH_UFINE}   M1z check {M1z_ok}')
print('  late B2 parts:', {k: '%.4g' % v for k, v in parts2_grouped.items()})
print('  late B1 parts:', {k: '%.3g' % fup(v) for k, v in pl[1].items()})
print(f'  slope share C1 Sigma_X / D_k: {slope_share};  S33+foot / Bhat2 {out["S33_foot_over_Bhat2"]:.3g};  K1 {K1_ok} K2 {K2_ok} K3 {K3_ok}   wall {time.time() - T0:.1f} s')
