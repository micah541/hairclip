"""v4b (fin13v4b) promotion of the E0 audit's e0_cert.py (fix5/e0_audit, Nov 19) to a certification stage. Changes against the audit script:
  inputs read from the step3nl_v4b copies (byte-identical to step3nl_v4; sha256 pinned), the [est] drift allowance replaced by the
  certified C_DR_CORE of e0_strip_drift.py, the absorption check made at the asserted level AMAX_E0 = 2.5e-4 (and 1.6e-4 for reference),
  and sections (8)-(12) added: the higher-order modulation pieces coded as explicit E0 terms, the forcing budget, and the verifier
  constants (raw certified value, outward-rounded value, source) written to out/e0v4b_constants.json.
Original docstring:
E0 audit, certified arithmetic part (python-flint arb, outward).  Read-only on the proof tree; writes e0_cert.json in its own directory.
Run:  cd <this dir>; PYTHONHASHSEED=0 nice -n 19 python3 -B e0_cert.py

Computes, from the certified Gram / pairing data of the modulation lemma:
  (1) n3_true = ||chi psi3^||_rho                (the typed n3 = 7.52 of ASSEMBLY.md sec. 3 is NOT this norm)
  (2) n3_12, n3_123 = norm of chi psi3^ projected off {f1,f2}, {f1,f2,f3}   (w~ is orthogonal to f1,f2,f3)
  (3) the split of the lemma's beta lead  BETA0 = offdiag(j=1,2) + |H33| eps1
  (4) the E0 constant Pi3 under four treatments  (T0 probe, T1 true norm, T2 projected norm, T3 projected + rank-one symmetrisation)
  (5) the E0 forcing factor under the same treatments (oblique projection identity)
  (6) the alpha-linear E0 modulation coefficient that is actually needed, against the coded MODA = 25
  (7) the no-rerun criterion:  Pi3_T3 + (MOD_need - 25) alpha_max <= 0.01607  and  f_T2 <= 5.72
Labels: every number printed is [cert] (arb, from certified JSON inputs taken as given) unless marked [hand] or [est]."""
import sys, os, json, hashlib, time, math
from flint import arb, ctx
ctx.prec = 200
V4 = '/home/claude/s1spec/step3nl_v4b'
sys.dont_write_bytecode = True
sys.path.insert(0, V4 + '/agent_mod/ibp')
import ibp_consts as K          # the lemma's own G, H, Gp, eps1, cF (reads step3nl/..., step3/lyap_consts2.json; read only)
t0 = time.time()
up = lambda z: float(arb(z).upper()); lo = lambda z: float(arb(z).lower())
mag = lambda z: arb(abs(z).upper())
def ball(l, h): return arb((arb(l) + arb(h)) / 2, 0).union(arb(l)).union(arb(h))
INP = {'gram7': V4 + '/agent_mod/bloop/gram7_g3.json', 'pair_f3': V4 + '/agent_mod/bloop/pair_f3.json', 'bloop': V4 + '/agent_mod/bloop/bloop.json',
       'final': V4 + '/agent_mod/ibp/final_5e-4.json', 'lyap2': '/home/claude/s1spec/step3/lyap_consts2.json',
       'lyap_extra': V4 + '/agent_mod/lemma/lyap_extra.json', 'lem_int': V4 + '/agent_mod/lemma/lem_integrals.json',
       'constants_v4': V4 + '/v4cert/out/constants_v4.json', 'h2res': V4 + '/h2res.json', 'q2_ip': V4 + '/q2_ip.json',
       'ibp_consts.py': V4 + '/agent_mod/ibp/ibp_consts.py'}
sha = {k: hashlib.sha256(open(p, 'rb').read()).hexdigest() for k, p in INP.items()}
J = {k: json.load(open(p)) for k, p in INP.items() if p.endswith('.json')}
G7 = [[ball(*J['gram7']['G'][i][j]) for j in range(7)] for i in range(7)]     # Lf1 Lf2 gA gB g0 f1 f2
C4 = {k: v['value'] for k, v in J['constants_v4']['constants'].items()}
A_ = lambda v: arb(repr(float(v))) if not isinstance(v, str) else arb(v)
out = {}
# ---------------------------------------------------------------- (1) true norm
n3t2 = G7[4][4]; n3_true = arb(n3t2.upper()).sqrt()
out['n3_true'] = dict(lo=lo(arb(n3t2.lower()).sqrt()), hi=up(n3_true), typed=7.52, note='||chi psi3^||_rho from gram7_g3.json G[g0][g0]')
# ---------------------------------------------------------------- (2) projected norms
M = [[G7[5][5], G7[5][6]], [G7[6][5], G7[6][6]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
b = [G7[4][5], G7[4][6]]
P12 = sum((b[i] * Mi[i][j] * b[j] for i in range(2) for j in range(2)), arb(0))
n12sq = n3t2 - P12; n3_12 = arb(n12sq.upper()).sqrt()
# f3 data: lemma's (psi3ip + strip) and the independent bloop pairing; both certified enclosures of <chi psi3^, f3> -> either may be used
g33_lem = K.g33; g33_bl = ball(*J['pair_f3']['g0'])
g33_int_lo = max(lo(g33_lem), lo(g33_bl)); g33_int_hi = min(up(g33_lem), up(g33_bl)); assert g33_int_lo < g33_int_hi
c13 = K.cvec           # [arb(0, g13), arb(0, g23)]: <f1,f3>, <f2,f3> (lemma)
nf3_hi = arb(J['lyap_extra']['nf3_hi'])
def n123(g33lo):
    bMc = sum((b[i] * Mi[i][j] * c13[j] for i in range(2) for j in range(2)), arb(0))
    num = arb(g33lo) - mag(bMc)                     # <Pi12perp g0, f3> >= g33_lo - |b M^-1 c|
    assert num.lower() > 0
    sub = (num * num) / (nf3_hi * nf3_hi)           # ||Pi12perp f3|| <= ||f3|| <= nf3_hi
    return arb((n12sq - arb(sub.lower())).upper()).sqrt()
n3_123_lem = n123(lo(g33_lem)); n3_123_int = n123(g33_int_lo)
out['n3_proj'] = dict(P12_lo=lo(P12), n3_12=up(n3_12), n3_123_lemmaG33=up(n3_123_lem), n3_123_intersectedG33=up(n3_123_int),
                      G33_lemma=[lo(g33_lem), up(g33_lem)], G33_bloop=[lo(g33_bl), up(g33_bl)], G33_intersection=[g33_int_lo, g33_int_hi])
# ---------------------------------------------------------------- (3) beta lead split (same formulas as ibp_consts.py:24-42, as a function of the G33 enclosure)
def lead_beta(g33):
    Mc = [sum((K.Mi_[a][bb] * K.cvec[bb] for bb in range(2)), arb(0)) for a in range(2)]
    s_ = g33 - sum((K.rvec[a] * Mc[a] for a in range(2)), arb(0)); assert s_.lower() > 0
    H31, H32, H33 = -Mc[0] / s_, -Mc[1] / s_, 1 / s_
    off = K.qform([H31, H32], K.Gp); return off, mag(H33), (H31, H32, H33)
off_l, H33_l, _ = lead_beta(g33_lem)
off_i, H33_i, _ = lead_beta(ball(g33_int_lo, g33_int_hi))
eps1 = K.eps1; r3 = arb(repr(J['lyap2']['r3']))
BETA0 = A_(C4['BETA0']) + A_(C4['BETA0_dlead']); BETAF = A_(C4['BETAF']); dlead = A_(C4['BETA0_dlead'])
chk = off_l + H33_l * eps1
assert abs(up(chk) - J['final']['final']['core']['beta']['lead']) < 1e-12, (up(chk), J['final']['final']['core']['beta']['lead'])
out['beta_lead_split'] = dict(lemma=dict(offdiag=up(off_l), H33=up(H33_l), diag_eps1=up(H33_l * eps1), total=up(chk)),
                              intersectedG33=dict(offdiag=up(off_i), H33=up(H33_i), diag_eps1=up(H33_i * eps1), total=up(off_i + H33_i * eps1)),
                              eps1=up(eps1), r3_read=up(r3), note='lead = qform([H31,H32],Gp) + |H33| eps1 ; reproduces final_5e-4.json core.beta.lead')
# ---------------------------------------------------------------- (4) Pi3 under the four treatments
PI3_coded = arb('0.01607')
# cross term X >= |<u, chi psi3^>|, u = (L* - lamt) f3 = u_far (supported on the strip sigma in [29,32], norm r3_far) + u_bump (norm r3_bump)   (lyap_consts2.py:143-163, 223)
#   X_direct [cert]: Cauchy-Schwarz piecewise:  r3_far * ||chi psi3^||_{L2(rho), strip} + r3_bump * ||chi psi3^||
#   X_dual   [hand]: <u, g0> = <f3, (L - lamt)(chi psi3~)>  =>  ||f3|| FN0 + sum_i kappa_3,i sup_strip|psi3^(i)|   (cross-check only)
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as PF
xs_lo = arb(29) - PF.STIP - arb('0.01'); xs_hi = arb(32) - PF.STIP + arb('0.01')        # strip in x = sigma - sigma_tip, padded
strip_sup_psi = [arb(0)] * 4; nb = 400
for i in range(nb):
    Xb = (xs_lo + (xs_hi - xs_lo) * i / nb).union(xs_lo + (xs_hi - xs_lo) * (i + 1) / nb)
    assert Xb.lower() > PF.XC.upper()
    jt = PF.ext_jet(Xb, 4)          # derivatives 0..3
    for j in range(4):
        v = arb(abs(jt[j]).upper()); strip_sup_psi[j] = v if v.upper() > strip_sup_psi[j].upper() else strip_sup_psi[j]
rho_s = (-(1 + xs_lo * xs_lo) ** (arb(2) / 3) / 4).exp()                              # rho decreasing: sup on the strip
g0_strip = strip_sup_psi[0] * ((xs_hi - xs_lo) * rho_s).sqrt()                          # ||chi psi3^||_{L2(rho)} on the strip (chi <= 1)
r3_far = arb(repr(J['lyap2']['r3_far'])); r3_bump = arb(repr(J['lyap2']['r3_bump']))
FN0 = A_(C4['FN0']); kap3 = [A_(x) for x in J['lem_int']['kappa_3']]
X_dual = nf3_hi * FN0 + sum((kap3[i] * strip_sup_psi[i] for i in range(4)), arb(0))
X = r3_far * g0_strip + r3_bump * n3_true
def T3(n123v, off, H33, nrm_u, Xv): return n123v * (off + dlead) + H33 * (nrm_u * n3_true + Xv) / 2
T = {}
T['T0_probe_7.52*BETA0'] = arb('7.52') * BETA0
T['T1_true_norm*BETA0'] = n3_true * BETA0
T['T2_proj123*BETA0_lemmaG33'] = n3_123_lem * BETA0
T['T2_proj123*BETA0_intersG33'] = n3_123_int * (off_i + H33_i * eps1 + dlead)
T['T3_sym_lemmaG33_eps1'] = T3(n3_123_lem, off_l, H33_l, eps1, X)
T['T3_sym_intersG33_eps1'] = T3(n3_123_int, off_i, H33_i, eps1, X)
T['T3_sym_intersG33_r3'] = T3(n3_123_int, off_i, H33_i, r3, X)
# robustness of T3: largest cross term X for which T3 (lemma G33, eps1) still equals the coded Pi3
Xmax = (PI3_coded - n3_123_lem * (off_l + dlead)) * 2 / H33_l - eps1 * n3_true
out['Pi3'] = dict(coded=0.01607, X_direct_cert=up(X), X_direct_parts=dict(far=up(r3_far * g0_strip), bump=up(r3_bump * n3_true)), X_dual_hand=up(X_dual),
                  strip_sup_psi_jets=[up(v) for v in strip_sup_psi], g0_strip_norm=up(g0_strip), eps1_n3_true=up(eps1 * n3_true),
                  treatments={k: up(v) for k, v in T.items()}, covered_by_coded={k: bool(v.upper() <= PI3_coded.lower()) for k, v in T.items()},
                  X_max_for_T3_le_coded=lo(Xmax), X_max_over_eps1_n3=lo(Xmax / (eps1 * n3_true)))
# ---------------------------------------------------------------- (5) forcing factor  (1 + n3 * c_Phi) Fn A
cF3 = K.cF[2]                                   # ||h||, h = sum_j H3j f_j  (0.62614)
Ff = {}
Ff['T0_probe_1+7.52*BETAF'] = 1 + arb('7.52') * BETAF
Ff['T1_1+n3true*BETAF'] = 1 + n3_true * BETAF
Ff['T2_1+n3_123*BETAF'] = 1 + n3_123_lem * BETAF
Ff['T3_oblique'] = n3_true * cF3 + n3_123_lem * (BETAF - cF3)        # ||I - h (x) g0|| = ||h|| ||g0|| since <g0,h> = 1 exactly; feedback part with the projected norm
# extra L-group pieces that E0 does not code: (a,b) forcing part * alpha ||(gA,gB)||, and -2 alpha beta_F <chi h2, w~>
Gg = [[G7[2][2], G7[2][3]], [G7[3][2], G7[3][3]]]; gABn = K.lam_max_sym2(Gg).sqrt()      # operator norm of w -> (<w,gA>,<w,gB>)
h2n = arb(J['h2res']['h2_norm2_half'][1]).sqrt()
AMAXS = ('1.6e-4', '2.5e-4')        # the E0 lemma is stated for |alpha| <= AMAX_E0 = 2.5e-4 (asserted in closure.py); 1.6e-4 for reference
PHIAB = A_(C4['PHIAB'])
extraF = {a: PHIAB * gABn * arb(a) + 2 * arb(a) * h2n * BETAF for a in AMAXS}
out['forcing'] = dict(coded=5.72, cPhi_norm_h=up(cF3), BETAF=up(BETAF), treatments={k: up(v) for k, v in Ff.items()},
                      extra_ab_h2_at_alpha={a: up(v) for a, v in extraF.items()},
                      T2_plus_extra={a: up(Ff['T2_1+n3_123*BETAF'] + v) for a, v in extraF.items()},
                      covered_by_coded={a: bool((Ff['T2_1+n3_123*BETAF'] + v).upper() <= 5.72) for a, v in extraF.items()})
# ---------------------------------------------------------------- (6) alpha-linear E0 modulation coefficient actually needed
fb = J['final']['final']['core']
Qfull = arb(repr(J['bloop']['Qfull_f123perp']['abs_upper']))            # Q1 + (E13,E23)-loop on {f1,f2,f3}^perp
m_new_b = arb(repr(fb['beta']['m_new']))                                # beta's alpha-linear part from <chi Lin(w), f_j> (NOT in Qfull)
m3 = arb(repr(J['final']['per_pairing']['m_j'][2]))
MAB = A_(C4['MAB']); lead_a = arb(repr(fb['a']['lead'])); lead_b = arb(repr(fb['b']['lead']))
rem = J['bloop']['alpha2_remainders']; r1, r2 = arb(repr(rem['r1'])), arb(repr(rem['r2']))
e13m = mag(ball(*J['bloop']['e13'])); e23m = mag(ball(*J['bloop']['e23'])); E33a = arb(repr(J['bloop']['E33_over_alpha=2<chi h2,f3>'][1]))
n3p = n3_123_lem
need = {}
for a in AMAXS:
    al = arb(a)
    parts = dict(
        Qfull=Qfull,
        beta_Lin_crude=n3p * m_new_b,                                            # n3_123 * 4.48
        E_loop_2nd_order=H33_l * (r1 * lead_a + r2 * lead_b) * n3p * al,          # (E13 - alpha e13, E23 - alpha e23) remainders
        ab_alpha_corr_in_Q1=MAB * gABn * al,                                      # (a,b) - leading, times alpha ||(gA,gB)||
        ab_alpha_corr_in_loop=H33_l * (e13m * e13m + e23m * e23m).sqrt() * MAB * n3p * al,
        h2_direct=2 * h2n * BETA0,                                                # -2 alpha beta <chi h2, w~>, beta lead part
        E33_feedback=n3p * H33_l * E33a * BETA0)
    tot = sum(parts.values(), arb(0))
    need[a] = dict(parts={k: up(v) for k, v in parts.items()}, total_without_drift=up(tot))
# symmetrised variant of the beta_Lin part: (1/2) H33 (||m3|| n3_true + |<m3, g0>|) ; |<m3,g0>| ~ 2 |<chi Q2(psi3), f3>| <= 0.0192 [est; q2_ip.json F hull, identity not certified]
q2F = max(abs(J['q2_ip']['F'][0]), abs(J['q2_ip']['F'][1]))
sym_lin = H33_l * (m3 * n3_true + 2 * arb(repr(q2F))) / 2
out['alpha_form'] = dict(coded_MODA=25.0, need=need, beta_Lin_symmetrised_est=up(sym_lin), gAB_opnorm=up(gABn), h2_norm=up(h2n),
                         note='drift alpha-form [est] is added in E0_AUDIT.md from e0_drift_est.py')
# ---------------------------------------------------------------- (7) no-rerun criterion (drift allowance DRIFT [est] read from e0_drift_est.json if present)
SD_PATH = os.path.join(V4, 'v4bcert', 'out', 'e0_strip_drift.json'); SD = json.load(open(SD_PATH))
assert SD['params']['AMAX'] == 2.5e-4
DR = SD['C_DR_CORE']            # [cert] drift alpha-form, core (sigma <= 29); the strip part is coded in state form
crit = {}
for a in AMAXS:
    al = arb(a)
    exc = (arb(repr(need[a]['total_without_drift'])) + arb(repr(DR)) - 25) * al
    for key in ('T3_sym_lemmaG33_eps1', 'T2_proj123*BETA0_lemmaG33', 'T2_proj123*BETA0_intersG33'):
        slack = PI3_coded - T[key]
        crit[f'{key}@alpha<={a}'] = dict(Pi3=up(T[key]), slack=lo(slack), excess_alpha_terms=up(exc), ok=bool(slack.lower() >= exc.upper()),
                                         margin=(lo(slack / exc) if exc.lower() > 0 else None))
out['no_rerun_criterion'] = dict(drift_alpha_form_cert=DR, cases=crit)
# ================================================================ v4b additions
AMAX_E0 = arb('2.5e-4'); WB = '0.02'
up_ = lambda z: math.nextafter(float(arb(z).upper()), math.inf)
offL, H33L, (H31, H32, H33s) = lead_beta(g33_lem)
# (8) higher-order modulation pieces, CODED as explicit E0 terms (state form; see E0_LEMMA_v4b.tex, Lemma part (iii))
#   Q1 (a,b) beyond lead:   |(a,b) - lead| alpha ||(gA,gB)|| A,    |(a,b) - lead| <= MAB|alpha| A + QAB Q + ABs + PHIAB Fn  (MAB part: E0-4 / absorbed)
#   loop (a,b) beyond lead: n3_123 H33 |(E13,E23)| |(a,b) - lead| A,  |(E13,E23)| <= |alpha| (|e| + AMAX_E0 |r|)
#   beta beyond lead in -2 alpha beta <chi h2,w~> and in the E33 feedback: alpha (2||h2|| + n3_123 H33 E33/alpha) |beta - beta_lead| A,
#        |beta - beta_lead| <= BETA1 |alpha| A + QB Q + BTs + BETAF Fn
#   H31, H32 part of beta_E:  n3_123 sum_{j=1,2} |H3j| (|(E1j,E2j)| |(a,b)| + |E3j| |beta|) A,  |E_kj| <= |alpha| Ea_kj  (lemma's E^alpha matrix)
#   alpha^2 remainders of A1(u0), B1(u0) in Q1:  |(a,b)| alpha^2 (RA2 + RB2) A
e_norm = (e13m * e13m + e23m * e23m).sqrt(); r_norm = (r1 * r1 + r2 * r2).sqrt()
LOOPC = n3p * H33L * (e_norm + AMAX_E0 * r_norm)
HOB = 2 * h2n + n3p * H33L * E33a
Ea = K.Ea
H12AB = n3p * (mag(H31) * (Ea[0][0] ** 2 + Ea[1][0] ** 2).sqrt() + mag(H32) * (Ea[0][1] ** 2 + Ea[1][1] ** 2).sqrt())
H12B = n3p * (mag(H31) * Ea[2][0] + mag(H32) * Ea[2][1])
RA2 = arb(repr(SD['RA2'])); RB2 = arb(repr(SD['RB2']))
# (9) forcing budget: every Fn-proportional piece of the E0 coupling against the coded factor 5.72, at |alpha| = AMAX_E0 (all pieces nondecreasing in alpha)
al = AMAX_E0
fparts = dict(T2=Ff['T2_1+n3_123*BETAF'], ab_forcing_Q1=PHIAB * gABn * al, h2_betaF=2 * al * h2n * BETAF, loop_ab_forcing=al * LOOPC * PHIAB,
              E33_betaF=al * n3p * H33L * E33a * BETAF, H12_forcing=al * (H12AB * PHIAB + H12B * BETAF), A1rem_forcing=al * al * (RA2 + RB2) * PHIAB)
ftot = sum(fparts.values(), arb(0))
assert ftot.upper() <= arb('5.72').lower(), ftot
out['v4b_forcing_budget'] = dict(alpha_max=0.00025, parts={k_: up(v) for k_, v in fparts.items()}, total=up(ftot), coded=5.72, margin=lo(arb('5.72') / ftot), ok=True)
# (10) A^2 budget (absorption by the Pi3 slack): Pi3_T3 + (need(alpha) + C_DR_CORE - 25) alpha <= 0.01607 for |alpha| <= AMAX_E0
#      need(alpha) = c0 + c1 alpha with c1 >= 0 (section (6)), so (need - 25) alpha is nondecreasing once need > 25: check at AMAX_E0.
ab_ = {}
for a_ in AMAXS:
    nd = arb(repr(need[a_]['total_without_drift'])) + arb(repr(DR))
    exc = (nd - 25) * arb(a_); sl = PI3_coded - T['T3_sym_lemmaG33_eps1']
    assert nd.lower() > 25 and sl.lower() > exc.upper()
    ab_[a_] = dict(need_total=up(nd), excess=up(exc), slack=lo(sl), margin=lo(sl / exc))
out['v4b_A2_budget'] = dict(cases=ab_, note='need = Qfull + n3_123 m_new + E0-4 pieces (section 6) + C_DR_CORE [cert]; coded MODA = 25, PI3 = 0.01607')
# (11) the coupling constant of the E0 beta-strip term (rhs_iv.py:184): -beta_S <g0, w~> = -beta_S <Pi g0, w~> since w~ is orthogonal to f1,f2,f3
N3P = n3_123_lem
# (12) verifier constants: raw certified value, outward value, source
AL = json.load(open(V4 + '/alin_E0_spec_0.001.json')); assert AL['alpha_max'] == 0.001
alin_raw = {nm: max(AL['sup'][r][key] for r in AL['sup']) for nm, key in (('ALIN_sq', 'spec01'), ('ALIN_AB', '|a02+a20|'), ('ALIN_W1B', '|a12+a21|'), ('ALIN_BB', '|a22|'))}
CONST = {}
def put(name, raw, src, d='up'):
    r = float(arb(raw).upper()) if isinstance(raw, arb) else float(raw)
    v = math.nextafter(r, math.inf) if d == 'up' else math.nextafter(r, -math.inf)
    assert (v >= r) if d == 'up' else (v <= r)
    CONST[name] = dict(raw=r, value=v, dir=d, source=src)
for nm, r in alin_raw.items(): put(nm, r, 'alin_E0_spec_0.001.json: max over regions (sum over regions <= max * global norms; Cauchy-Schwarz)')
put('N3P', N3P, 'n3_123 = ||Pi_{f1,f2,f3}^perp chi psi3^|| (lemma G33), this script (2)')
put('PHX', A_(C4['PHX']), 'constants_v4.json PHX (Phi1 = phi\'(34.25), SOB_NOTES (I2))')
for nm in ('M_STRIP_00', 'M_STRIP_01', 'M_STRIP_02', 'M_STRIP_11', 'M_STRIP_12', 'M_STRIP_22', 'IR', 'IRW', 'RW', 'XW_COEF', 'CD1', 'CD0', 'CA1', 'CA0', 'CB1', 'CB0', 'RA2', 'RB2', 'C_DR_CORE'):
    put(nm, SD[nm], 'e0_strip_drift.json ' + nm)
put('GABN', gABn, '||(gA,gB)|| operator norm, gram7 (section 5)'); put('LOOPC', LOOPC, 'n3_123 H33 (|(e13,e23)| + AMAX_E0 |(r1,r2)|), bloop.json')
put('HOB', HOB, '2||chi h2|| + n3_123 H33 E33/alpha'); put('H12AB', H12AB, 'n3_123 sum_{j=1,2}|H3j| |(Ea_1j,Ea_2j)| (ibp_consts Ea)'); put('H12B', H12B, 'n3_123 sum_{j=1,2}|H3j| Ea_3j')
CONST['AMAX_E0'] = dict(raw=2.5e-4, value=2.5e-4, dir='exact', source='hypothesis of the E0 lemma, asserted in closure.box_assert')
CONST['WB'] = dict(raw=0.02, value=0.02, dir='exact', source='|w| <= WB on [0,34.25]: hypothesis of e0_strip_drift (alpha_g lower bound), asserted in closure.box_assert')
out['v4b_constants'] = CONST
out['v4b_H3j'] = dict(H31=[lo(H31), up(H31)], H32=[lo(H32), up(H32)], H33=[lo(H33s), up(H33s)])
INP['e0_strip_drift'] = SD_PATH; INP['alin_spec'] = V4 + '/alin_E0_spec_0.001.json'
sha = {k: hashlib.sha256(open(p, 'rb').read()).hexdigest() for k, p in INP.items()}
# the audit's inputs are byte-identical copies of step3nl_v4 (pinned against the audit's e0_cert.json)
AUD = json.load(open('/home/claude/s1spec/writeup/lemmaU/bulk/fix5/e0_audit/e0_cert.json'))['args']['inputs']
for k_, d_ in AUD.items(): assert sha[k_] == d_['sha256'], k_
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, prec=ctx.prec,
                   inputs={k: dict(path=INP[k], sha256=sha[k]) for k in INP}, script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
OUTD = os.path.join(V4, 'v4bcert', 'out')
json.dump(out, open(os.path.join(OUTD, 'e0_lemma_consts.json'), 'w'), indent=1, sort_keys=True)
json.dump(dict(constants=CONST, sources=dict(e0_lemma_consts=os.path.join(OUTD, 'e0_lemma_consts.json'))), open(os.path.join(OUTD, 'e0v4b_constants.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps({k: out[k] for k in ('v4b_forcing_budget', 'v4b_A2_budget', 'v4b_H3j')}, indent=1))
print(json.dumps({k: (v['raw'], v['value']) for k, v in CONST.items()}, indent=0))
