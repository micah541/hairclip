"""fix11 review, part A: independent arb recomputation of the fin13v4b E0/E2 lemma constants and of the R-3 whole-tail alpha-condition.
Read-only on the proof tree; writes a_recompute.json next to this script.  [cert] = arb (200 bit, outward) on the certified JSON inputs taken as given.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B a_recompute.py"""
import sys, os, json, hashlib, math
sys.dont_write_bytecode = True
from flint import arb, ctx
ctx.prec = 200
V = '/home/claude/s1spec/step3nl_v4b'
HERE = os.path.dirname(os.path.abspath(__file__))
INP = dict(gram7=V + '/agent_mod/bloop/gram7_g3.json', pair_f3=V + '/agent_mod/bloop/pair_f3.json', bloop=V + '/agent_mod/bloop/bloop.json',
           final=V + '/agent_mod/ibp/final_5e-4.json', lyap2='/home/claude/s1spec/step3/lyap_consts2.json', lyap_extra=V + '/agent_mod/lemma/lyap_extra.json',
           psi3ip='/home/claude/s1spec/step3nl/psi3ip.json', lem_int=V + '/agent_mod/lemma/lem_integrals.json', gram8='/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json',
           c4=V + '/v4cert/out/constants_v4.json', h2res=V + '/h2res.json', e0v4b=V + '/v4bcert/out/e0v4b_constants.json', e2v4b=V + '/v4bcert/out/e2v4b_constants.json',
           e0lc=V + '/v4bcert/out/e0_lemma_consts.json', sd=V + '/v4bcert/out/e0_strip_drift.json', r2=V + '/v4bcert/out/r2_strip_ibp.json',
           nl=V + '/v4bcert/out/nl_regions_a2.5e-4.json', gard=V + '/garding_cert.json', alin=V + '/alin_E0_spec_0.001.json',
           rt15=V + '/runs/P/rate_table_fin13v4b_a0.00015.json', rt125=V + '/runs/P/rate_table_fin13v4b_a0.000125.json')
J = {k: json.load(open(p)) for k, p in INP.items()}
up = lambda z: float(arb(z).upper()); lo = lambda z: float(arb(z).lower())
mag = lambda z: arb(abs(arb(z)).upper())
def ball(l, h):
    l, h = arb(repr(l)) if not isinstance(l, arb) else l, arb(repr(h)) if not isinstance(h, arb) else h
    return l.union(h)
R = lambda v: arb(repr(float(v)))
out = {}; checks = []
def chk(name, ok, **kw): checks.append(dict(name=name, ok=bool(ok), **kw)); print(('PASS ' if ok else 'FAIL ') + name, kw)
# ---------------------------------------------------------------- 1. projected norm of g0 = chi psi3^ (own formula: Gram determinant ratio)
G = [[ball(*J['gram7']['G'][i][j]) for j in range(7)] for i in range(7)]      # Lf1 Lf2 gA gB g0 f1 f2
g00 = arb(G[4][4].upper()); n3 = g00.sqrt()
M = [[G[5][5], G[5][6]], [G[6][5], G[6][6]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; assert det.lower() > 0
b = [G[4][5], G[4][6]]
P12 = (b[0] * b[0] * M[1][1] - 2 * b[0] * b[1] * M[0][1] + b[1] * b[1] * M[0][0]) / det       # |P_{f1,f2} g0|^2
n12sq = G[4][4] - P12
# f3: <g0,f3> >= lemma lower bound (psi3ip - SP y_F); <f_i,f3> = +-(psif3 + SP); ||f3|| <= nf3_hi
SP = J['lem_int']['SP']; P3 = J['psi3ip']; LY = J['lyap2']; LX = J['lyap_extra']
g33_lo = arb(repr(P3['psi3_f3'][0])) - arb(repr(SP['y_F'])); g33_hi = arb(repr(P3['psi3_f3'][1])) + arb(repr(SP['y_F']))
c13 = [arb(0, LY['psif3'][0] + SP['1_F']), arb(0, LY['psif3'][1] + SP['2_F'])]
bMc = (b[0] * (M[1][1] * c13[0] - M[0][1] * c13[1]) + b[1] * (-M[1][0] * c13[0] + M[0][0] * c13[1])) / det
num = arb(g33_lo.lower()) - mag(bMc); assert num.lower() > 0
nf3 = arb(repr(LX['nf3_hi']))
n123 = arb((n12sq - arb((num * num / (nf3 * nf3)).lower())).upper()).sqrt()
E0 = {k: v['value'] for k, v in J['e0v4b']['constants'].items()}; E2 = {k: v['value'] for k, v in J['e2v4b']['constants'].items()}
C4 = {k: v['value'] for k, v in J['c4']['constants'].items()}
out['n3_true'] = up(n3); out['n3_123'] = up(n123)
chk('N3P coded >= my n3_123', E0['N3P'] >= up(n123), coded=E0['N3P'], mine=up(n123))
chk('N3E2 coded >= ||g0||', E2['N3E2'] >= up(n3), coded=E2['N3E2'], mine=up(n3))
chk('7.52 >= n3_123 (qb752 literal)', 7.52 >= up(n123))
# ---------------------------------------------------------------- 2. H row 3 and the beta lead (own block inverse)
GR = [[ball(*x) for x in row] for row in J['gram8']['G']] if isinstance(J['gram8'], dict) and 'G' in J['gram8'] else None
out['gram8_keys'] = list(J['gram8'].keys()) if isinstance(J['gram8'], dict) else 'list'
fb = J['final']['final']['core']
H33_up = 1 / arb((arb(g33_lo.lower()) - mag(arb(0, 1) * 0)).lower())       # crude: 1/g33_lo without the (r M^-1 c) correction; compared below
lc = J['e0lc']
H33_l = arb(repr(lc['beta_lead_split']['lemma']['H33'])); off_l = arb(repr(lc['beta_lead_split']['lemma']['offdiag']))
out['H33_crude_1_over_g33lo'] = up(H33_up); out['H33_lemma'] = up(H33_l)
chk('H33 lemma value consistent with 1/g33_lo (rel < 1e-3)', abs(up(H33_up) / up(H33_l) - 1) < 1e-3, crude=up(H33_up), lemma=up(H33_l))
eps1 = arb(repr(LY['eps1'])); r3_far = arb(repr(LY['r3_far'])); r3_bump = arb(repr(LY['r3_bump'])); r3 = arb(repr(LY['r3']))
chk('r3 = r3_far + r3_bump (support split of u3) and r3 <= eps1', abs(up(r3_far + r3_bump) / up(r3) - 1) < 1e-9 and up(r3) <= lo(eps1), r3=up(r3), far=up(r3_far), bump=up(r3_bump), eps1=up(eps1))
# X >= |<u3, g0>|: strip part r3_far * ||g0||_{strip}, bump part r3_bump * ||g0||;  ||g0||_strip <= sup|psi3^| sqrt(len * rho(x_lo)), x in [31.24, 34.26]
sup_psi = arb(repr(lc['Pi3']['strip_sup_psi_jets'][0]))                    # [read] 719.62 (psi3fun ext_jet); not recomputed
xl = arb('31.24'); rho_s = (-(1 + xl * xl) ** (arb(2) / 3) / 4).exp(); g0s = sup_psi * (arb('3.02') * rho_s).sqrt()
X = r3_far * g0s + r3_bump * n3
dleadb = R(C4['BETA0_dlead'])
Pi3_true = n123 * (off_l + dleadb) + H33_l * (eps1 * n3 + X) / 2
out['X'] = up(X); out['Pi3_true'] = up(Pi3_true)
chk('Pi3_true <= 0.01607', up(Pi3_true) <= 0.01607, mine=up(Pi3_true), reported=0.0103854)
# ---------------------------------------------------------------- 3. A^2 budget at alpha_E = 2.5e-4, with the two pieces the lemma omits
aE = arb('2.5e-4')
Qfull = R(J['bloop']['Qfull_f123perp']['abs_upper']); m_new = R(fb['beta']['m_new'])
MAB = R(C4['MAB']); BETA0 = R(C4['BETA0']) + dleadb; BETAF = R(C4['BETAF']); PHIAB = R(C4['PHIAB'])
lead_a = R(fb['a']['lead']); lead_b = R(fb['b']['lead'])
r1 = R(J['bloop']['alpha2_remainders']['r1']); r2 = R(J['bloop']['alpha2_remainders']['r2'])
e13 = mag(ball(*J['bloop']['e13'])); e23 = mag(ball(*J['bloop']['e23'])); E33a = R(J['bloop']['E33_over_alpha=2<chi h2,f3>'][1])
Gg = [[G[2][2], G[2][3]], [G[3][2], G[3][3]]]; t = Gg[0][0] + Gg[1][1]; d = Gg[0][0] * Gg[1][1] - Gg[0][1] * Gg[1][0]
gAB = arb(((t + arb(max(up(t * t - 4 * d), 0)).sqrt()) / 2).upper()).sqrt()
h2n = arb(repr(J['h2res']['h2_norm2_half'][1])).sqrt()
en = (e13 * e13 + e23 * e23).sqrt(); rn = (r1 * r1 + r2 * r2).sqrt()
CDR = R(J['sd']['C_DR_CORE'])
parts = dict(Qfull=Qfull, beta_Lin=n123 * m_new, loop_2nd=H33_l * (r1 * lead_a + r2 * lead_b) * n123 * aE, ab_corr_Q1=MAB * gAB * aE,
             ab_corr_loop=H33_l * en * MAB * n123 * aE, h2_direct=2 * h2n * BETA0, E33_fb=n123 * H33_l * E33a * BETA0, drift=CDR)
need = sum(parts.values(), arb(0))
DAB = R(E2['DAB']); dleadab = R(C4['CAB0_dlead'])
LOOPC = n123 * H33_l * (en + aE * rn)
extra = dict(lead_mismatch_dab=(DAB + dleadab) * (gAB + LOOPC),                      # (a,b)_lemma-lead - Minv l, coupled through alpha (gA,gB) and the loop: NOT in the lemma's budget
             loop_r_times_MAB=H33_l * aE * rn * MAB * n123 * aE)                    # alpha^2 |r| * MAB alpha A: third order, NOT in the lemma's budget
need2 = need + sum(extra.values(), arb(0))
slack = arb('0.01607') - Pi3_true
exc = (need - 25) * aE; exc2 = (need2 - 25) * aE
out['A2_budget'] = dict(parts={k: up(v) for k, v in parts.items()}, need_reported_form=up(need), omitted={k: up(v) for k, v in extra.items()}, need_with_omitted=up(need2),
                        slack=lo(slack), excess=up(exc), excess_with_omitted=up(exc2), margin=lo(slack / exc), margin_with_omitted=lo(slack / exc2))
chk('A^2 budget (reported form) margin x1.679 reproduced', abs(lo(slack / exc) - 1.6792) < 2e-3, mine=lo(slack / exc))
chk('A^2 budget still holds with the omitted pieces', lo(slack) > up(exc2), margin=lo(slack / exc2))
chk('LOOPC coded >= mine', E0['LOOPC'] >= up(LOOPC) * (1 - 1e-12), coded=E0['LOOPC'], mine=up(LOOPC))
chk('GABN coded >= mine', E0['GABN'] >= up(gAB) * (1 - 1e-12), coded=E0['GABN'], mine=up(gAB))
HOB = 2 * h2n + n123 * H33_l * E33a
chk('HOB coded >= mine', E0['HOB'] >= up(HOB) * (1 - 1e-12), coded=E0['HOB'], mine=up(HOB))
# ---------------------------------------------------------------- 4. forcing budgets
RA2 = R(E0['RA2']); RB2 = R(E0['RB2']); H12AB = R(E0['H12AB']); H12B = R(E0['H12B'])
f0 = dict(T2=1 + n123 * BETAF, abQ1=PHIAB * gAB * aE, h2=2 * aE * h2n * BETAF, loop=aE * LOOPC * PHIAB, E33=aE * n123 * H33_l * E33a * BETAF,
          H12=aE * (H12AB * PHIAB + H12B * BETAF), A1rem=aE * aE * (RA2 + RB2) * PHIAB)
ftot = sum(f0.values(), arb(0)); out['E0_forcing_budget'] = dict(parts={k: up(v) for k, v in f0.items()}, total=up(ftot), margin=lo(arb('5.72') / ftot))
chk('E0 forcing budget <= 5.72', up(ftot) <= 5.72, total=up(ftot))
FE2 = 1 + BETAF * (R(E2['N3E2']) + 2 * aE * R(E2['NH2R'])) + PHIAB * (R(E2['NF12']) + aE * (R(E0['GABN']) + aE * (RA2 + RB2)))
out['FE2_at_2.5e-4'] = up(FE2); chk('FE2 = 13.98 at alpha_E', abs(up(FE2) - 13.984) < 2e-3, mine=up(FE2))
# ---------------------------------------------------------------- 5. rounding directions of the read constants
AL = J['alin']; raw = {nm: max(AL['sup'][r][key] for r in AL['sup']) for nm, key in (('ALIN_sq', 'spec01'), ('ALIN_AB', '|a02+a20|'), ('ALIN_W1B', '|a12+a21|'), ('ALIN_BB', '|a22|'))}
chk('ALIN constants >= raw region maxima', all(E0[k] >= v for k, v in raw.items()), raw=raw)
g = J['gard']['sup']; chk('R143 >= max(-mu2)', E2['R143'] >= max(v['-mu2'] for v in g.values()), raw=max(v['-mu2'] for v in g.values()))
chk('typed r163 16.3, r107 10.7 >= raw', 16.3 >= max(v['|e1|'] for v in g.values()) and 10.7 >= max(v['|e0|'] for v in g.values()), e1=max(v['|e1|'] for v in g.values()), e0=max(v['|e0|'] for v in g.values()))
S = J['nl']['sup']; col = {1: '(0, 0, 1)', 2: '(0, 1, 0)', 3: '(1, 0, 0)'}
ok = True
for i, key in col.items(): ok &= E2[f'ND{i}'] >= max(S['R1'][key], S['R2a'][key]) and E2[f'NF{i}'] >= max(S['R2b'][key], S['R3'][key])
ok &= E2['ND4'] >= max(S['R1']['A0'], S['R2a']['A0']) and E2['ND5'] >= max(S['R1']['A1'], S['R2a']['A1'])
for i, key in ((6, '(2, 0, 0)'), (7, '(1, 1, 0)'), (8, '(3, 0, 0)')): ok &= E2[f'ND{i}'] >= max(S['R1'][key], S['R2a'][key])
for i, key in ((4, '(2, 0, 0)'), (5, '(1, 1, 0)'), (6, '(3, 0, 0)')): ok &= E2[f'NF{i}'] >= max(S['R2b'][key], S['R3'][key])
ok &= E2['NF7'] >= max(S['R2b']['A0'], S['R3']['A0'], S['R2b']['A1'], S['R3']['A1'])
chk('nd/nf >= region maxima of nl_regions_a2.5e-4.json (column map as rhs_iv ND)', ok)
# the L2-type columns: the coded coefficient of each L2 norm is nd_k + nf_k >= max over ALL regions; A0 multiplies A with nd4 + nf7, A1 multiplies W1 with nd5 + nf7
chk('A-coefficient nd4+nf7 >= max A0 all regions; W1-coefficient nd5+nf7 >= max A1', E2['ND4'] + E2['NF7'] >= max(S[r]['A0'] for r in ('R1', 'R2a', 'R2b', 'R3')) and E2['ND5'] + E2['NF7'] >= max(S[r]['A1'] for r in ('R1', 'R2a', 'R2b', 'R3')))
# ---------------------------------------------------------------- 6. R-3: whole-tail alpha-condition, both alpha points
LAMT = arb('-0.2561532'); BETA1 = R(C4['BETA1']); QB = R(C4['QB']); FN = [R(C4['FN0']), R(C4['FN1']), R(C4['C3'])]
SB = R(J['r2']['STRIP_BETA_R2'])
out['R3'] = {}
for key in ('rt15', 'rt125'):
    T = J[key]['tail']; a_s = R(T['a_star']); Ls = R(T['L_star']); rho = R(T['rho']); yF = R(T['sqrt_LF']); kap = a_s / Ls.sqrt()
    # QBL: read back from the stored anchor rhs (rhs_at_anchor = BETA0 + BETA1 a* + (BETAF Fn(a*) + STRIP)/sq + QBL sq)
    ac = T['alpha_condition_whole_tail']; sq = Ls.sqrt()
    Fn = lambda a: FN[0] * a + FN[1] * a * a + FN[2] * a * a * a
    base = BETA0 + BETA1 * a_s + BETAF * Fn(a_s) / sq
    QBLsq = R(ac['rhs_at_anchor']) - base - SB / sq                        # >= true QBL*sq (stored value is rounded up)
    lhs = (-LAMT - rho / 2) * kap
    def rhs(l):                                                             # beta(kappa sqrt(l), l)/sqrt(l), l = L* * s^2
        s = (l / Ls).sqrt(); y = l.sqrt()
        return BETA0 + BETA1 * kap * y + BETAF * Fn(kap * y) / y + QBLsq * s + SB / y
    worst = None
    n = 4000
    for i in range(n + 1):
        y = yF * (sq / yF) ** (arb(i) / n); v = rhs(y * y)
        worst = v if worst is None or v.upper() > worst.upper() else worst
    rhs_all = base + mag(QBLsq) + SB / yF
    out['R3'][key] = dict(lhs=lo(lhs), rhs_all_monotone=up(rhs_all), stored_rhs_whole_tail=ac['rhs_whole_tail'], grid_max=up(worst), difference=lo(lhs - rhs_all), ratio=lo(lhs / rhs_all),
                          strip_at_LF=up(SB / yF), QBL_sq_backed_out=up(QBLsq))
    chk(f'R-3 whole-tail alpha-condition {key}', lo(lhs) > up(rhs_all) and up(rhs_all) <= ac['rhs_whole_tail'] * (1 + 1e-9), lhs=lo(lhs), rhs=up(rhs_all), ratio=lo(lhs / rhs_all))
    chk(f'R-3 box/hypothesis ratios {key} all <= 1', all(r <= 1 for row in J[key]['box_hypotheses'] for r in (row['ratio'].values() if isinstance(row, dict) and 'ratio' in row else [0])))
# R-2: new constants below the old ones and S3 coefficient zero
chk('R-2: S3 coefficients are 0 and constants consumed by the verifier are the R2 ones', J['r2']['beta']['const_coef_S0_S3'][3] <= 1e-300 and J['r2']['ab']['const_coef_S0_S3'][3] <= 1e-300 and J['r2']['STRIP_BETA_R2'] == __import__('json').load(open(V + '/v4bcert/out_A/r2_strip_ibp.json'))['STRIP_BETA_R2'], STRIP_BETA_R2=J['r2']['STRIP_BETA_R2'], STRIP_AB_R2=J['r2']['STRIP_AB_R2'])
out['checks'] = checks; out['n_fail'] = sum(1 for c in checks if not c['ok'])
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, prec=ctx.prec,
                   inputs={k: dict(path=p, sha256=hashlib.sha256(open(p, 'rb').read()).hexdigest()) for k, p in INP.items()},
                   script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
json.dump(out, open(os.path.join(HERE, 'a_recompute.json'), 'w'), indent=1, sort_keys=True)
print('n_fail', out['n_fail']); print(json.dumps({k: out[k] for k in ('n3_true', 'n3_123', 'X', 'Pi3_true', 'A2_budget', 'E0_forcing_budget', 'FE2_at_2.5e-4', 'R3')}, indent=1))
