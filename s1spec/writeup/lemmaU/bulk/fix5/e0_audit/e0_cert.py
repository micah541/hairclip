"""E0 audit, certified arithmetic part (python-flint arb, outward).  Read-only on the proof tree; writes e0_cert.json in its own directory.
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
import sys, os, json, hashlib, time
from flint import arb, ctx
ctx.prec = 200
V4 = '/home/claude/s1spec/step3nl_v4'
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
AMAXS = ('1.6e-4', '2.5e-4')
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
DR = 0.0
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e0_drift_est.json')
if os.path.exists(p): DR = json.load(open(p))['alpha_form_allowance']
crit = {}
for a in AMAXS:
    al = arb(a)
    exc = (arb(repr(need[a]['total_without_drift'])) + arb(repr(DR)) - 25) * al
    for key in ('T3_sym_lemmaG33_eps1', 'T2_proj123*BETA0_lemmaG33', 'T2_proj123*BETA0_intersG33'):
        slack = PI3_coded - T[key]
        crit[f'{key}@alpha<={a}'] = dict(Pi3=up(T[key]), slack=lo(slack), excess_alpha_terms=up(exc), ok=bool(slack.lower() >= exc.upper()),
                                         margin=(lo(slack / exc) if exc.lower() > 0 else None))
out['no_rerun_criterion'] = dict(drift_allowance_est=DR, cases=crit)
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, prec=ctx.prec, inputs={k: dict(path=INP[k], sha256=sha[k]) for k in INP},
                   runtime_s=round(time.time() - t0, 2))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e0_cert.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'args'}, indent=1))
