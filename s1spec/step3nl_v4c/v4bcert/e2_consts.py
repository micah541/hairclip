"""v4b promotion of fix5/e2_audit/e2_consts.py (paths -> step3nl_v4b, output -> v4bcert/out). E2 audit: certified constants and comparisons (arb, 200 bits, outward).  Reads only; writes e2_consts.json in the cwd.
Sections: (a) LAMT against the lamt of the residual norms; (b) the D_rho constants 2.2, 1.2; (c) rounding of the typed Gaarding and
nonlinear constants against their own certified JSONs; (d) modulation constants of the E2 level; (e) strip weight ratio.
usage: python3 -B e2_consts.py"""
import sys, os, json, hashlib
from flint import arb, ctx
ctx.prec = 200
S3 = '/home/claude/s1spec/step3nl_v4c/'            # fin13v4 tree (read only); E2 files are byte-identical to step3nl / step3nl_v33
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod/ibp'); sys.path.insert(0, '/home/claude/s1spec/step3')
up = lambda z: float(arb(z).upper()); lo = lambda z: float(arb(z).lower())
mag = lambda z: arb(abs(z).upper())
out = {}; ins = []
def J(p):
    ins.append(p); return json.load(open(p))

# ---------------------------------------------------------------- (a) LAMT
import gen_hs3 as G3
ins.append('/home/claude/s1spec/step3/gen_hs3.py')
lam_res = arb(G3.LAMT)                           # the binary double used by psi3fun.LAM, h2fun.LAM (2x), ibp_boxes.LAMT, lyap_consts2.LAMT, c4fun
lam_ver = arb('-0.2561532')                      # verifier/constants.py:14, converted by closure.py:18 as arb(str(K.LAMT))
d = lam_ver - lam_res
assert d.lower() > 0
import ast, re
src = open(S3 + 'verifier/constants.py').read(); ins.append(S3 + 'verifier/constants.py')
assert re.search(r'^LAMT = -0\.2561532; PI3', src, re.M)
out['a_LAMT'] = dict(lamt_residual_double=G3.LAMT, lamt_residual_repr=repr(G3.LAMT), LAMT_verifier='-0.2561532',
                     verifier_minus_residual=[lo(d), up(d)], relative=up(d / abs(lam_res)),
                     direction="LAMT_verifier > lamt_residual (less negative): |alpha|' <= lamt_res|alpha| + |beta| <= LAMT_ver|alpha| + |beta|; conservative",
                     forcing_equivalent_of_mismatch_per_alpha=up(d * arb('7.6970')))
# ---------------------------------------------------------------- (b) ||D_rho g|| <= ||g''''|| + sup 2phi' ||g'''|| + sup|phi'^2 - phi''| ||g''||
def phi1(x): return x * (1 + x * x) ** (arb(-1) / 3) / 3
def phi2(x): return (1 + x * x) ** (arb(-1) / 3) / 3 - arb(2) / 9 * x * x * (1 + x * x) ** (arb(-4) / 3)
s1 = s2 = 0.0; N = 34250; XM = arb('34.25')
for i in range(N):
    X = (XM * i / N).union(XM * (i + 1) / N)
    s1 = max(s1, up(2 * abs(phi1(X)))); s2 = max(s2, up(abs(phi1(X) ** 2 - phi2(X))))
out['b_Drho'] = dict(sup_2phi1=s1, typed_t22=2.2, sup_abs_phi1sq_minus_phi2=s2, typed_t12=1.2, ok=bool(s1 <= 2.2 and s2 <= 1.2))
# ---------------------------------------------------------------- (c) rounding of typed constants
gc = J(S3 + 'garding_cert.json')['sup']
g = {k: max(v[k] for v in gc.values()) for k in ('-mu2', '|e1|', '|e0|')}
typed = dict(r143=14.34, r163=16.3, r107=10.7)
out['c_garding'] = {n: dict(typed=t, raw=g[k], ok=bool(t >= g[k]), typed_minus_raw=t - g[k]) for (n, t), k in zip(typed.items(), ('-mu2', '|e1|', '|e0|'))}
nl = J(S3 + 'nlin_E2_0.001_0.0001_0.0001.json')['sup']; T, Fr = nl['x<=12'], nl['12<x<=34.25']
tab = [('nd1', 6.9e-3, T['(0, 0, 1)']), ('nd2', 1.2e-2, T['(0, 1, 0)']), ('nd3', 4.2e-2, T['(1, 0, 0)']), ('nd4', 0.54, T['A0']), ('nd5', 0.62, T['A1']),
       ('nd6', 4.75, T['(2, 0, 0)']), ('nd7', 0.136, T['(1, 1, 0)']), ('nd8', 3.05, T['(3, 0, 0)']),
       ('nf1', 4.3e-3, Fr['(0, 0, 1)']), ('nf2', 6.9e-4, Fr['(0, 1, 0)']), ('nf3', 2.9e-5, Fr['(1, 0, 0)']), ('nf4', 0.014, Fr['(2, 0, 0)']),
       ('nf5', 0.46, Fr['(1, 1, 0)']), ('nf6', 3.06, Fr['(3, 0, 0)']), ('nf7(A0)', 8.6e-6, Fr['A0']), ('nf7(A1)', 8.6e-6, Fr['A1'])]
out['c_nlin_typed_vs_own_json'] = {n: dict(typed=t, raw=r, ok=bool(t >= r), rel=(t - r) / r) for n, t, r in tab}
mf = J(S3 + 'agent_mod/mod_forms.json')
out['c_modE2'] = dict(typed=43.33, raw_L2max_upper=mf['L2max_upper'], ok=bool(43.33 >= mf['L2max_upper']), slack=43.33 - mf['L2max_upper'], L2min=mf['L2min_enclosure'])
# ---------------------------------------------------------------- (d) modulation constants at the E2 level
import ibp_consts as IK                      # read-only import: Gram, H = (G^-1)^T, Minv, eps1, c_j (as the modulation lemma)
ins += ['/home/claude/s1spec/step3nl/agent_mod/ibp/ibp_consts.py', '/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json', '/home/claude/s1spec/step3nl/agent_mod/lemma/lem_integrals.json',
        '/home/claude/s1spec/step3/lyap_consts2.json', '/home/claude/s1spec/step3nl/agent_mod/lemma/lyap_extra.json', '/home/claude/s1spec/step3nl/psi3ip.json', '/home/claude/s1spec/step3nl/h2res.json']
H, Mi, GR, eps1, cj = IK.H, IK.Mi, IK.GR, IK.eps1, IK.c_j
def lam_max2(S):
    t = S[0][0] + S[1][1]; dd = S[0][0] * S[1][1] - S[0][1] * S[1][0]; disc = arb(max((t * t - 4 * dd).upper(), 0))
    return arb((t.upper() + disc.sqrt().upper()) / 2)
def proj(a, b): return GR[a][b] - sum((GR[a][6 + i] * Mi[i][j] * GR[6 + j][b] for i in range(2) for j in range(2)), arb(0))
GD = [[proj(a, b) for b in (4, 5)] for a in (4, 5)]                # Pi_{12}-projected Gram of (D f1, D f2)
NDF = lam_max2(GD).sqrt()                                           # |a <w,Df1> + b <w,Df2>| <= |(a,b)| NDF ||w||, w in {f1,f2}^perp
NDF_un = lam_max2([[GR[a][b] for b in (4, 5)] for a in (4, 5)]).sqrt()
NF12 = lam_max2([[GR[a][b] for b in (6, 7)] for a in (6, 7)]).sqrt()   # |a <f1,g> + b <f2,g>| <= |(a,b)| NF12 ||g||
NGAB = lam_max2([[GR[a][b] for b in (2, 3)] for a in (2, 3)]).sqrt()   # same for (gA, gB)
# lead used by Q2:  y_i = sum_j Mi[j][i] l_j ;  lemma lead: sum_{j<=3} H[i][j] R_j.   difference:
dH = [[mag(H[i][j] - Mi[j][i]) for j in range(2)] for i in range(2)]
dl = [sum((dH[i][j] * cj[j] for j in range(2)), arb(0)) + mag(H[i][2]) * eps1 for i in range(2)]
dab = (dl[0] ** 2 + dl[1] ** 2).sqrt()
dab_H3 = (mag(H[0][2]) ** 2 + mag(H[1][2]) ** 2).sqrt() * eps1
g7 = J(S3 + 'agent_mod/bloop/gram7_g3.json'); assert g7['names'][4] == 'g0'
n3 = arb(g7['G'][4][4][1]).sqrt(); nh2 = IK.h2n
c4 = J(S3 + 'agent_int/c4switch/c4switch_S1.json')
n3c4 = n3 + arb(c4['psi3']['rho'][0]); nh2c4 = nh2 + arb(c4['h2']['rho'][0])     # C^4 versions (rho-norm of the difference, half line; chi <= 1)
F5 = J(S3 + 'agent_mod/ibp/final_5e-4.json')['final']['core']
CAB0, DLab, MAB, QAB, PHIAB = (arb(F5['ab'][k]) for k in ('lead', 'dlead', 'm_total', 'q_total', 'cPhi_total'))
B0, DLb, B1, QB, BF = (arb(F5['beta'][k]) for k in ('lead', 'dlead', 'm_total', 'q_total', 'cPhi_total'))
# v4b certified w-part / alpha^2 constants (read from the running fin13v4b work if present; else the audit's own crude values)
V4B = '/home/claude/s1spec/step3nl_v4c/v4bcert/out/e0v4b_constants.json'
RA2 = RB2 = None
try:
    vb = json.load(open(V4B)); out['d_v4b_constants_seen'] = {k: (v.get('value') if isinstance(v, dict) else v) for k, v in (vb.get('constants', vb)).items() if k in ('RA2', 'RB2', 'CA1', 'CA0', 'CB1', 'CB0', 'AG_MIN', 'NGAB', 'IR', 'IRW', 'RW')}
except Exception as e:
    out['d_v4b_constants_seen'] = 'not readable: ' + repr(e)
eps = arb('1e-5'); Q2 = arb(mf['L2max_upper'])
D = dict(NDF_proj12=up(NDF), NDF_unprojected=up(NDF_un), NF12=up(NF12), NGAB=up(NGAB), n3=up(n3), n3_c4=up(n3c4), nh2=up(nh2), nh2_c4=up(nh2c4),
         dH_2x2=[[up(x) for x in r] for r in dH], H_i3=[up(mag(H[0][2])), up(mag(H[1][2]))], eps1=up(eps1), c_j=[up(x) for x in cj],
         dab_lead_minus_Q2lead=up(dab), dab_part_H3=up(dab_H3), CAB0=up(CAB0), dlead_ab=up(DLab), MAB=up(MAB), QAB=up(QAB), PHIAB=up(PHIAB),
         BETA0=up(B0), dlead_beta=up(DLb), BETA1=up(B1), QB=up(QB), BETAF=up(BF), Q2_cert=up(Q2))
D['lead_A2_needed_IBPform'] = up(Q2 + NDF * (dab + DLab)); D['lead_A2_shortfall_vs_43.33'] = up(Q2 + NDF * (dab + DLab) - arb('43.33'))
D['forcing_factor_E2_needed_CS'] = up(1 + BF * (n3c4 + 2 * arb('2.5e-4') * nh2c4) + PHIAB * NF12 + PHIAB * arb('2.5e-4') * (NGAB + arb('2.5e-4') * 33))
rows = {}
for a_ in ('3.6846e-5', '1.6e-4', '2.5e-4'):
    a = arb(a_)
    rab = dab + DLab + MAB * a                         # |(a,b) - Minv l| / A      (x-part, linear in A)
    bt = B0 + DLb + B1 * a                             # |beta| / A
    ab = CAB0 + DLab + MAB * a                         # |(a,b)| / A
    cs_ab = rab * NF12; cs_bt = bt * (n3c4 + 2 * a * nh2c4); cs_u0 = ab * a * (NGAB + a * 33)       # 33 >= RA2 + RB2 (v4b: 26.902 + 5.4997)
    cRM = cs_ab + cs_bt + cs_u0
    ibp_ab = NDF * rab
    rows[a_] = dict(rab=up(rab), beta_over_A=up(bt), ab_over_A=up(ab), CS_ab_rest=up(cs_ab), CS_beta=up(cs_bt), CS_A1B1_u0=up(cs_u0), c_RM_coefficient_of_A_Dq=up(cRM),
                    relative_to_coded_nd4_plus_nf7=up(cRM / arb('0.5400086')), IBP_alt_ab_rest_A2=up(ibp_ab),
                    A2_equiv_Young_est=up(eps * ((arb('0.5400086') + cRM) ** 2 - arb('0.5400086') ** 2) / 2), A2_equiv_IBP_alt=up(eps * ibp_ab))
D['per_alpha'] = rows
D['cubic_ab_w_part_coeff_of_A_times_W1_Dq'] = up((CAB0 + DLab + MAB * arb('2.5e-4')) * arb('33.04'))      # (CA1 + CB1) = 32.000 + 1.03219 (v4b)
D['quadratic_Q_coefficient_times_Dq'] = up(QAB * NF12 + QB * (n3c4 + 2 * arb('2.5e-4') * nh2c4))
out['d_modulation'] = D
# ---------------------------------------------------------------- (e) strip: sup over sigma in [29,32] of (rho/W)^(1/2) and of rho
# x >= 31.24 on the strip (x = sigma - sigma_tip, sigma_tip = -2.24908, strip enclosure [31.249, 34.249]);  psi - phi is decreasing there
# (psi' - phi' = (1-s)(4x/(1+x^2) - phi') < 0), so the sup is at the left end.  psi(x) <= 2 log(1+x^2) + s(x) * int_{30.25}^{x} phi'.
xl = arb('31.24'); t = (xl - arb('30.25')) / 8; S7 = t ** 4 * (35 - 84 * t + 70 * t * t - 20 * t ** 3)
phi = lambda x: (1 + x * x) ** (arb(2) / 3) / 4
psi_up = 2 * (1 + xl * xl).log() + S7 * (phi(xl) - phi(arb('30.25')))
SRW = ((psi_up - phi(xl)) / 2).exp()
co = J(S3 + 'agent_coup/coupling_constants.json')
out['e_strip'] = dict(SRW_sup_sqrt_rho_over_W=up(SRW), rho_at_31p24=up((-phi(xl)).exp()), IR_crude_3rho=up(3 * (-phi(xl)).exp()),
                      check_d3_over_4supchi1=co['P']['2']['rigorous']['d3'] / (4 * 35 / 48), note="d3 = sup 4|chi'|(rho/W)^(1/2) <= 4 (35/48) SRW must hold: " + str(bool(co['P']['2']['rigorous']['d3'] <= 4 * 35 / 48 * up(SRW))))
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, prec=ctx.prec,
                   inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in sorted(set(ins + [os.path.abspath(__file__)]))})
json.dump(out, open('/home/claude/s1spec/step3nl_v4c/v4bcert/out/e2_consts.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'args'}, indent=1))
