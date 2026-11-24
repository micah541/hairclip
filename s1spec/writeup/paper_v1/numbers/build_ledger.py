#!/usr/bin/env python3
"""Builds ledger.json and numbers_paper.tex (unit U0).  Reads JSON sources only; no heavy computation.
usage: nice -n 19 python3 -I build_ledger.py"""
import json, hashlib, os, re, sys
from fractions import Fraction as Fr
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_numbers as C

B = '/home/claude/s1spec/writeup/lemmaU/bulk'
P4 = '/home/claude/s1spec/step4'
W = '/home/claude/s1spec/writeup'
SRC = {
 # constants_needed_v34.json sources (v35 chain)
 'A':  B + '/fix2/n12/t2_A_v35/out/A_interface.json',
 'A4': B + '/fix2/n12/t2_A_v35/out/asm_design_N2ref4.json',
 'U':  B + '/t4_cert/out/U_cert_v35.json',
 'U1': B + '/t4_cert/out/u1_ghat.json',
 'Z':  B + '/t5_zone3/out/zone3_interface.json',
 'R':  B + '/t5_zone3/out/ramp_U.json',
 'S4': P4 + '/p2/out_v33/FINAL_manifest_v33.json',
 'F12': P4 + '/p2/out_v33/f12/f12_sliver.json',
 # step 4
 'S4C': P4 + '/p2/out/FINAL_manifest.json',
 'T1': P4 + '/p2/out/t1clock.json',
 'H1S': P4 + '/p2/out/h1slab/summary.json',
 'CLAM': P4 + '/p2/out/t9_h1_clam_s180.json',
 'SW': P4 + '/p2/out/switch_s2.json',
 'G0': P4 + '/t2/gamma0_final.json',
 'P1S': P4 + '/t3/p1/status.json',
 'AMP': P4 + '/p2/t8se/constants/amp_A1_w180.json',
 'H4L': P4 + '/p2/out/h4loc_A1_w180.json',
 'H4F': P4 + '/p2/t8se/constants/h4loc_v32.json',
 'T9F': P4 + '/p2/out/t9_far_v32.json',
 'T9C': P4 + '/p2/out/t9_cone_kch.json',
 'T7C': P4 + '/t7/constants.json',
 # step 3 verifier and Lemma O
 'VR': '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json',
 'MS': '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json',
 'LO': W + '/lemmaO_v32/lemmaO_v32_consts.json',
 # step 1, step 2
 'BR2': '/home/claude/s1_tight/bracket2.json',
 'EXTT': '/home/claude/s1_tight/ext_tip/ext_tip_check.json',
 'E813': '/home/claude/s1_tight/ext813/ext813_check.json',
 'CHK': '/home/claude/s1_tight/consistency/check_outA5.json',
 'CONT': '/home/claude/s1spec/cap/contour_result.json',
 'INT': '/home/claude/s1spec/cap/interior_result.json',
 'ETA': '/home/claude/s1spec/cap/eta_cache.json',
 'CE0': '/home/claude/s1spec/apr/cert_b0.25_th0.0_tau0.03.json',
 'CE1': '/home/claude/s1spec/apr/cert_b0.25_th0.48_tau0.03.json',
}
DOCS = {
 'MAIN3': B + '/t7_assembly/main_v3.tex',
 'V33I': B + '/t7_assembly/V33_INDEX.md',
 'S4RES': P4 + '/STEP4_RESULT.md',
 'S4ADD': B + '/t7_assembly/STEP4_RESULT_v33_addendum.md',
 'S1TIGHT': '/home/claude/s1_tight/README.md',
 'S1TIP': '/home/claude/s1_tight/ext_tip/README.md',
 'S3V32': W + '/step3_proof_v32/step3_proof_v32.tex',
 'N3LOG': P4 + '/p2/out/n3boxes_final_B.log',
 'CLAUDE': '/home/claude/CLAUDE.md',
 'S2TEX': W + '/step2_proof.tex',
}
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()

ENTRIES = []; NAMES = set()
def fmt_of(sig=None, dec_=None, style=None):
    if sig: return {'mode': 'sig', 'n': sig, 'style': style or 'sci'}
    return {'mode': 'dec', 'n': dec_, 'style': 'plain'}
def E(name, direction, fmt, expr='a', inputs=None, kind='cert', note='', src=None, key=None, flag=False, auto=False):
    """inputs: list of (var, srctag, key); or src/key for the single-input default."""
    if inputs is None: inputs = [('a', src, key)]
    name = 'N' + name
    assert name not in NAMES, name; NAMES.add(name)
    ENTRIES.append(dict(name=name, kind=kind, direction=direction, fmt=fmt, expr=expr,
        inputs=[dict(var=v, source=s, key=k) for v, s, k in inputs], note=note, flag=flag, auto=auto))
def DOC(name, direction, value, doc, snippet, cls, note, fmt=None):
    name = 'N' + name
    assert name not in NAMES, name; NAMES.add(name)
    ENTRIES.append(dict(name=name, kind='doc', direction=direction, doc_value=value, doc_file=DOCS[doc], doc_snippet=snippet,
        doc_class=cls, note=note, flag=True, fmt=fmt, inputs=[]))
def HAND(name, value, note, fmt=None):
    """hypothesis constant of the paper (exact decimal literal)"""
    x = Fr(Decimal(value)); sig = len(Decimal(value).normalize().as_tuple().digits)
    style = 'sci' if 'e' in value.lower() else 'plain'
    if style == 'plain':
        nd = max(0, -Decimal(value).as_tuple().exponent); f = {'mode': 'dec', 'n': nd, 'style': 'plain'}
    else:
        f = {'mode': 'sig', 'n': sig, 'style': 'sci'}
    E(name, 'exact', fmt or f, expr=value.lower(), inputs=[], kind='hand', note=note)

# ---------------------------------------------------------------- 1. quotes of constants_needed_v34.json (the Lemma S/U/zone-III chain)
CN = json.load(open(B + '/fix1/lemmaU_rev/constants_needed_v34.json'))
DIGIT = ['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine']
def words(s):
    toks = re.findall(r'[A-Z]+(?![a-z])|[A-Z]?[a-z]+|\d', s); out = ''
    for t in toks: out += DIGIT[int(t)] if t.isdigit() else t[0].upper() + t[1:]
    return out
TAGNAME = {'A': 'SA', 'A4': 'SAfour', 'U': 'U', 'U1': 'Uone', 'Z': 'Z', 'R': 'R', 'S4': 'Sfour'}
CURATED = {
 ('U', 'supersolution.EA.0'): 'UEhatZero', ('U', 'supersolution.EA.1'): 'UEhatOne', ('U', 'supersolution.EA.2'): 'UEhatTwo',
 ('U', 'supersolution.Bhat.0'): 'UBhatZero', ('U', 'supersolution.Bhat.1'): 'UBhatOne', ('U', 'supersolution.Bhat.2'): 'UBhatTwo',
 ('U', 'ufine.Bhat.0'): 'UfineBhatZero', ('U', 'ufine.Bhat.1'): 'UfineBhatOne', ('U', 'ufine.Bhat.2'): 'UfineBhatTwo',
 ('U', 'margins.late_M1'): 'UMarginLateMOne', ('U', 'margins.late_M0'): 'UMarginLateMZero',
 ('U', 'margins.present_M1_min.value'): 'UMarginPresentMOne', ('U', 'margins.present_M0_min.value'): 'UMarginPresentMZero',
 ('U', 'loop_gain.2.total'): 'ULoopGain', ('U', 'rate_margin'): 'URateMargin', ('U', 'launch.c_dU'): 'UcdeltaU',
 ('A4', 'res.tight.margins.M1core'): 'SAMarginCoreMOne', ('A4', 'res.tight.margins.M0core'): 'SAMarginCoreMZero',
 ('A4', 'res.tight.margins.S0v_vs_0.0153'): 'SAMarginCoreSZeroV', ('A4', 'res.tight.margins.S1v_vs_0.6'): 'SAMarginCoreSOneV',
 ('A4', 'res.tight.margins.Ghat22'): 'SAMarginGhat', ('A4', 'res.tight.margins.ramp M1'): 'SAMarginRampMOne',
 ('A4', 'res.tight.flux.beta.2.margin_FLUX0'): 'SAMarginFluxZero', ('A4', 'res.tight.flux.beta.2.margin_FLUX2'): 'SAMarginFluxTwo',
 ('Z', 'box_III_380.b0'): 'ZoneIIIBoxZero', ('Z', 'box_III_380.b1'): 'ZoneIIIBoxOne', ('Z', 'box_III_380.bz'): 'ZoneIIIBoxZ',
 ('Z', 'box_III_380.margins.b1'): 'ZoneIIIMarginThinnest', ('Z', 'cone_family.H6.at_Bhat_ref.min_margin'): 'RecoveryMinMargin',
 ('Z', 'acceptance.tstar_margin'): 'ConeTstarMargin',
}
seen = {}
for q, f, k, d in CN['numbers']:
    seen.setdefault((f, k), []).append((q, d))
def parse_quote(q):
    style = 'sci' if 'e' in q.lower() else 'plain'
    m = q.lower().split('e')[0]
    digs = m.replace('.', '').lstrip('0')
    return len(digs), style
for (f, k), lst in seen.items():
    q, d = max(lst, key=lambda t: parse_quote(t[0])[0])
    sig, style = parse_quote(q)
    nm = CURATED.get((f, k))
    if nm is None: nm = TAGNAME[f] + words(k.replace('res.tight.', ''))
    base = nm; i = 0
    while 'N' + nm in NAMES: i += 1; nm = base + 'V' + DIGIT[i]
    direction = 'upper' if d == 'up' else ('margin' if re.search(r'margin', k) else 'lower')
    E(nm, direction, fmt_of(sig, style=style), src=f, key=k, auto=True,
      note='quoted in constants_needed_v34.json as %s%s' % (', '.join(sorted(set(x[0] for x in lst))), '' if True else ''))

# ---------------------------------------------------------------- 2. curated: hypotheses (hand), margins against them
HAND_LIST = [
 ('HOneAlphaThr', '1.5e-4', '(H1) |alpha| threshold alpha_1'),
 ('HOneSqrtLamThr', '5e-8', '(H1) sqrt(Lambda) threshold'),
 ('YbarStep', '1e-12', 'step-4 error majorant Ybar'),
 ('HTwoThrZero', '1e-8', 'tightened (H2) level 0'), ('HTwoThrOne', '1e-6', 'tightened (H2) level 1'), ('HTwoThrTwo', '1e-3', 'tightened (H2) level 2'),
 ('HThreeThrZero', '1e-6', '(H3) (archive H2f) level 0, beta_A = 2'), ('HThreeThrOne', '5e-6', '(H3) level 1'), ('HThreeThrTwo', '1e-5', '(H3) level 2'),
 ('HFourThrZero', '2e-7', '(H4) (archive H3) window threshold level 0'), ('HFourThrOne', '1e-5', '(H4) level 1'),
 ('HFourThrTwo', '1e-5', '(H4) level 2'), ('HFourThrThree', '1e-5', '(H4) level 3'), ('HFourKfourThr', '0.06', '(H4) K4 threshold'),
 ('HSixConeSlopeThr', '0.08', '(H6) fixed-L cone slope datum threshold'), ('HSixConeKtwoThr', '0.03', '(H6) cone sup|h0\'\'| threshold'),
 ('HSixConeNthreeThr', '0.004', '(H6) cone L2 norm of h0\'\'\' threshold'), ('HSixConeJetThr', '0.12', '(H6) cone jets K2,K3,K4 vs M = 0.12'),
 ('HSixConeHFiveThr', '0.0144', '(H6) cone Lemma O H4 threshold M^2 with M = 0.12'),
 ('HSixFarHFiveThr', '0.0225', '(H6) far Lemma O H4 threshold M^2 with M = 0.15'),
 ('LemmaOConeM', '0.12', 'Lemma O cone family M_O'), ('LemmaOConeK', '0.18', 'Lemma O cone family K_O'),
 ('LemmaOFarM', '0.15', 'Lemma O far family M_O'), ('LemmaOFarK', '0.3', 'Lemma O far family K_O'),
 ('CoercivityConst', '0.046875', 'Re<Lv,v>_rho <= -(3/64)||v||^2 ; 3/64'),
 ('MOnez', '0.534', 'M_{1z}: hypothesis constant of the bootstrap (lemmaU_v34 Sec. 9.2)'),
 ('BetaA', '2', 'beta_A of the repaired Lemma S(A)'),
]
for nm, v, note in HAND_LIST: HAND(nm, v, note)
# the dependency of the 0.046875 literal: 3/64 exact
# ---------------------------------------------------------------- 3. step 4 (Theorem D, T3)
def Eh(name, direction, sig, src, key, style='sci', **kw): E(name, direction, fmt_of(sig, style=style), src=src, key=key, **kw)
E('TOneLower', 'lower', fmt_of(dec_=10), src='T1', key='t1.0', note='t1 interval, lower end (t1clock, slab 180)')
E('TOneUpper', 'upper', fmt_of(dec_=10), src='T1', key='t1.1', note='t1 interval, upper end')
E('LambdaOneLower', 'lower', fmt_of(3, style='sci'), expr='exp(a)', src='T1', key='lnlam.0', note='lambda_1 = exp(ln lambda), lower end; ln(lambda) in [lnlam]; exp evaluated to 50 digits')
E('LambdaOneUpper', 'upper', fmt_of(3, style='sci'), expr='exp(a)', src='T1', key='lnlam.1', note='lambda_1 upper end')
E('TMinusTOneUpper', 'upper', fmt_of(2, style='sci'), expr='exp(4*a)/(1-4*b)',
  inputs=[('a', 'T1', 'lnlam.1'), ('b', 'MS', '1.5e-4.sup_ab')], note='T - t1 <= lambda_1^4/(1-4 abar)')
E('AbarUpper', 'upper', fmt_of(6, style='sci'), src='MS', key='1.5e-4.sup_ab', note='sup |a|,|b| of the verifier, rounded up to 6 digits in the file (modsups_fin12f4)')
E('BetabarUpper', 'upper', fmt_of(6, style='sci'), src='MS', key='1.5e-4.sup_beta', note='sup |beta| (modsups_fin12f4)')
E('BaseKmaxLower', 'lower', fmt_of(4, style='sci'), src='T1', key='base_kmax.0', note='base curve max curvature at t1, lower end')
E('BaseKmaxUpper', 'upper', fmt_of(4, style='sci'), src='T1', key='base_kmax.1', note='base curve max curvature at t1, upper end')
E('HOneAlphaBase', 'upper', fmt_of(5, style='sci'), src='H1S', key='alpha_abs_max', note='(H1) |alpha| base, slab 180')
E('HOneAlphaErr', 'upper', fmt_of(3, style='sci'), expr='a*b', inputs=[('a', 'CLAM', 'C_alpha'), ('b', 'AMP', 'inputs.YSTAR')],
  note='C_alpha * Ybar (C_alpha=0.7207 up); main_v3 Table 3 quotes 7.3e-13 (looser)')
E('HOneAlphaMargin', 'margin', fmt_of(3, style='plain'), expr='1.5e-4/(a+b*c)', inputs=[('a', 'H1S', 'alpha_abs_max'), ('b', 'CLAM', 'C_alpha'), ('c', 'AMP', 'inputs.YSTAR')],
  note='1.5e-4 / (|alpha| + C_alpha Ybar), Ybar = YSTAR = 1e-12')
E('HOneSqrtLamBase', 'upper', fmt_of(4, style='sci'), src='CLAM', key='sqrtLam_base_Wp', note='base sqrt(Lambda), slab 180')
E('HOneSqrtLamTotal', 'upper', fmt_of(3, style='sci'), src='CLAM', key='sqrtLam_total_Wp', note='sqrt(Lambda) with error at Ybar = 1e-12')
E('HOneSqrtLamMargin', 'margin', fmt_of(3, style='plain'), expr='5e-8/a', src='CLAM', key='sqrtLam_total_Wp', note='5e-8 / sqrtLam_total')
E('HOneYbarLimit', 'lower', fmt_of(3, style='sci'), src='CLAM', key='Ybar_limit_H1', note='H1 holds for Ybar <= this (was misquoted 1.58e-9)')
E('SwitchYUpper', 'upper', fmt_of(3, style='sci'), src='SW', key='modulated.Y_raw', note='Y(tau_a) at the switch, base slab 2')
E('SwitchYMargin', 'margin', fmt_of(2, style='plain'), expr='1e-12/a', src='SW', key='modulated.Y_raw', note='Ybar / Y(tau_a); STEP4_RESULT quotes x1.6')
E('MinCI', 'lower', fmt_of(4, style='plain'), src='T7C', key='meta.summary.min_c_I', note='min over the 26 certificates; docs quote 0.05284 (nearest), floor is 0.05283')
E('NumCertificates', 'exact', {'mode': 'dec', 'n': 0, 'style': 'plain'}, src='T7C', key='meta.summary.n_cert', note='number of coercivity certificates')
E('MinRhoC', 'lower', fmt_of(4, style='plain'), src='T7C', key='meta.closure.rho_c_min', note='min varrho_c; the file stores the already-rounded 0.04711')
E('HTwoBaseZero', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2_tightened.rows.0.base', note='tightened (H2), w\', base')
E('HTwoBaseOne', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2_tightened.rows.1.base', note='')
E('HTwoBaseTwo', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2_tightened.rows.2.base', note='')
for i, nm in enumerate(['Zero', 'One', 'Two']):
    E('HTwoErrW' + nm, 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2_tightened.rows.%d.with_error' % i, note="tightened (H2) for w', with error at Ybar = 1e-12")
    E('HTwoMarginW' + nm, 'margin', fmt_of(5, style='plain'), src='S4', key='results.M9_H2_tightened.rows.%d.margin' % i, note="margin for w'")
E('HThreeBaseZero', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2f.rows.0.base', note='(H3) (archive H2f) base')
E('HThreeBaseOne', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2f.rows.1.base', note='')
E('HThreeBaseTwo', 'upper', fmt_of(4, style='sci'), src='S4', key='results.M9_H2f.rows.2.base', note='')
for i, nm in enumerate(['Zero', 'One', 'Two']):
    E('HThreeErr' + nm, 'upper', fmt_of(5, style='sci'), src='S4', key='results.M9_H2f.rows.%d.with_error' % i, note='(H3) with error at Ybar = 1e-12 (stage M9_H2f, arb)')
    E('HThreeMargin' + nm, 'margin', fmt_of(5, style='plain'), src='S4', key='results.M9_H2f.rows.%d.margin' % i, note='(H3) margin (arb)')
for j, nm in enumerate(['Zero', 'One', 'Two', 'Three']):
    E('HFourBase' + nm, 'upper', fmt_of(3, style='sci'), src='S4', key='results.M9_base_sums.combined.windows_max.%d' % j, note='(H4) (archive H3) window energy, base, arb')
    E('HFourErr' + nm, 'upper', fmt_of(5, style='sci'), src='S4', key='results.M9_H3_H2v32.rows.win.%d.with_error_arb' % j, note='(H4) window energy with error at Ybar = 1e-12 (arb)')
    E('HFourMargin' + nm, 'margin', fmt_of(4, style='plain'), src='S4', key='results.M9_H3_H2v32.rows.win.%d.margin_arb' % j, note='(H4) margin = threshold / with-error value (arb). R10: the settled values')
E('HFourYbarLimit', 'lower', fmt_of(4, style='sci'), src='S4', key='results.M9_H3_H2v32.binding_Ybar_limit', note='binding (E_3) Ybar limit of (H4)')
E('HFourYbarMargin', 'margin', fmt_of(4, style='plain'), expr='a/1e-12', src='S4', key='results.M9_H3_H2v32.binding_Ybar_limit', note='Ybar limit / 1e-12 = 15.75 (main_v3 prints x15.7)')
E('HFourKfourBase', 'upper', fmt_of(5, style='sci'), src='AMP', key='inputs.BASE_T1.K4', note='(H4) K4 base value at t1 (the step-4 base polynomials)')
# (H6) cone / far data
E('HSixConeSlope', 'upper', fmt_of(4, style='plain'), src='S4', key='results.cone_fixedL.checks.slope.value', note='fixed-L cone slope datum')
E('HSixConeSlopeMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.cone_fixedL.checks.slope.margin', note='main_v3 prints x2.64')
E('HSixConeKtwo', 'upper', fmt_of(4, style='plain'), src='S4', key='results.cone_fixedL.checks.K2.value', note='fixed-L cone sup|h0\'\'|')
E('HSixConeKtwoMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.cone_fixedL.checks.K2.margin')
E('HSixConeNthree', 'upper', fmt_of(4, style='sci'), src='S4', key='results.cone_fixedL.checks.N3.value')
E('HSixConeNthreeMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.cone_fixedL.checks.N3.margin')
E('HSixConeHFive', 'upper', fmt_of(5, style='sci'), src='S4', key='results.H4_fixedL.total', note='cone Lemma O H4 integral, worst cylinder')
E('HSixConeHFiveMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.H4_fixedL.margin')
E('HSixConeTstarNeeded', 'upper', fmt_of(8, style='sci'), src='S4', key='results.H7_fixedL.ratio_upper', note='needed time ratio (R10: the doc quotes 1.6000014e-7, rounded down)')
E('HSixConeTstarLower', 'lower', fmt_of(6, style='sci'), src='S4', key='results.H7_fixedL.t_star_lower')
E('HSixConeTstarMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.H7_fixedL.margin_tstar_over_ratio')
E('HSixConeHFiveRatioPeak', 'upper', fmt_of(4, style='plain'), src='H4L', key='windows.constants/amp_A1_w180.json.worst.ratio', note='v3.2 unit-L cone H4 / M^2 (superseded geometry; kept for D15 part K)')
E('HSixConeHOneSup', 'upper', fmt_of(4, style='plain'), src='S4', key='results.H1_fixedL.sup_h0_max', note='Lemma O H1 sup|h0|, fixed-L cone')
E('HSixConeHOneH', 'upper', fmt_of(4, style='plain'), src='S4', key='results.H1_fixedL.H_max')
E('HSixFarHFiveRatio', 'upper', fmt_of(4, style='plain'), src='H4F', key='windows.constants/amp_A1_w180.json.cylinders.18.ratio', note='v3.2 far H4/M^2 at M=0.15 (centre 376)')
E('HSixFarHFiveMarginOld', 'margin', fmt_of(3, style='plain'), expr='1/a', src='H4F', key='windows.constants/amp_A1_w180.json.cylinders.18.ratio', note='x5.68 of STEP4_RESULT (v3.2, centre 376)')
E('HSixFarHFiveMargin', 'margin', fmt_of(4, style='plain'), src='S4', key='results.M9_H4_far.min_margin', note='far-family Lemma O H4, every centre >= 490 (37 rows). R10: x15.58, reports print x15.59')
E('HSixConeHFiveMarginOld', 'margin', fmt_of(3, style='plain'), expr='1/a', src='H4L', key='windows.constants/amp_A1_w180.json.worst.ratio', note='x4.18 (v3.2 unit-L cone, closed manifest)')
for key, nm, sg in [('h2', 'Ktwo', 4), ('h3', 'Kthree', 4), ('h4', 'Kfour', 4), ('slope', 'Slope', 4), ('L2', 'Nthree', 4)]:
    E('HSixFar' + nm, 'upper', fmt_of(sg, style='sci' if key != 'slope' else 'plain'), src='T9F', key=key, note='far-family datum (lambda_1 units), slab 180 hulls, arb')
    E('HSixFar' + nm + 'Thr', 'exact', fmt_of(3, style='sci' if key != 'slope' else 'plain'), src='T9F', key='thr.' + key, kind='hand', note='far-family threshold (zone3 handoff_thresholds, v3.2 amendment)')
    E('HSixFar' + nm + 'Margin', 'margin', fmt_of(4, style='plain'), src='T9F', key='margins.' + key, note='margin = threshold / bound')
E('FarTstarMargin', 'margin', fmt_of(3, style='plain'), src='LO', key='far_v32.case.margin_factor_lower', note='far family t* margin (R10: x15.5)')
E('FarTstarLower', 'lower', fmt_of(5, style='sci'), src='LO', key='far_v32.case.t_star_lower')
E('FarTstarNeeded', 'upper', fmt_of(4, style='sci'), src='LO', key='far_v32.case.needed')
# step 4 Gamma0 and Phase I
E('GammaModes', 'exact', fmt_of(dec_=0), src='G0', key='M', note='sine modes of Gamma_0')
E('GammaHalfwidth', 'upper', fmt_of(4, style='sci'), src='G0', key='b_j_halfwidth', note='half-width of the b_2 Krawczyk interval (R10: docs print 8.137e-24, rounded down; up is 8.138e-24)')
E('GammaContraction', 'upper', fmt_of(3, style='sci'), src='G0', key='krawczyk_contraction')
E('GammaArcs', 'exact', fmt_of(dec_=0), src='G0', key='embedded.arcs')
E('GammaMaxTurn', 'upper', fmt_of(4, style='plain'), src='G0', key='embedded.maxturn', note='total turning of the half-curve < pi/2')
E('GammaMaxPair', 'upper', fmt_of(4, style='plain'), src='G0', key='embedded.maxpair', note='adjacent pairs turn by <= this < pi')
E('GammaMinSep', 'lower', fmt_of(3, style='sci'), src='G0', key='embedded.minsep', note='minimal separation of non-adjacent arc boxes')
E('GammaHalfGap', 'lower', fmt_of(dec_=7), src='G0', key='embedded.gap_cert', note='slot half-gap')
E('PhaseISlabs', 'exact', fmt_of(dec_=0), src='P1S', key='slab', note='last Phase-I slab (slabs 0..96)')
E('PhaseIModes', 'exact', fmt_of(dec_=0), src='P1S', key='last_slab.M', note='Fourier modes on the last Phase-I slab')
E('PhaseIKmaxLower', 'lower', fmt_of(dec_=9), src='P1S', key='kmax_enclosure.0', note='k_max at slab 96, lower')
E('PhaseIKmaxUpper', 'upper', fmt_of(dec_=9), src='P1S', key='kmax_enclosure.1')
E('PhaseITLower', 'lower', fmt_of(dec_=15), src='P1S', key='t', note='t_96 = 4.59785053948583083 +- 7.77e-19')
E('PhaseITUpper', 'upper', fmt_of(dec_=15), src='P1S', key='t')
# ---------------------------------------------------------------- 4. step 3 tip zone (verifier fin12f4)
E('PeakSqrtLam', 'upper', fmt_of(4, style='sci'), src='VR', key='sqrtL_max', note='peak of sqrt(Lambda) along the comparison trajectory (alpha_1 = 1.5e-4)')
E('TailRate', 'lower', fmt_of(4, style='sci'), src='VR', key='tail.rho', note='tail decay rate of the comparison system (lower bound)')
E('FloorSqrtLam', 'upper', fmt_of(4, style='sci'), src='VR', key='tail.sqrt_LF', note='floor of sqrt(Lambda)')
E('AlphaFloor', 'upper', fmt_of(4, style='sci'), src='VR', key='tail.alpha_F', note='alpha_F, the floor of |alpha|')
E('SupBetaVerifier', 'upper', fmt_of(6, style='sci'), src='MS', key='1.5e-4.sup_beta')
E('MOneRamp', 'upper', fmt_of(4, style='plain'), src='MS', key='1.5e-4.M1', note='M1 on the ramp (verifier)')
E('SITwoPrime', 'upper', fmt_of(5, style='plain'), src='S4', key='results.F12_sliver.sI2_prime', note="sup|w''| <= sI2' on the sliver (stage f12_sliver)")
E('LoopGainCoarse', 'upper', fmt_of(4, style='plain'), src='U', key='loop_gain.2.total', note='0.03904 as quoted in the outline')
E('RateMarginCoarse', 'margin', fmt_of(3, style='plain'), src='U', key='rate_margin', note='x28.2')
# ---------------------------------------------------------------- 5. step 1 and step 2
E('SOneEpsLower', 'lower', fmt_of(14, style='sci'), src='BR2', key='eps_a', note='eps* bracket (tight)')
E('SOneEpsUpper', 'upper', fmt_of(14, style='sci'), src='BR2', key='eps_b')
E('SOneRowsTip', 'exact', fmt_of(dec_=0), src='EXTT', key='rows', note='51,246 CAPD boxes')
E('SOneTipWidth', 'upper', fmt_of(2, style='sci'), src='EXTT', key='max_pointwise_width.th', note='tangent-angle pointwise width near the tip (conditioned): the 3.2e-12 of main_v3')
E('SOneExtWidth', 'upper', fmt_of(2, style='sci'), src='E813', key='max_pointwise_width.p', note='width of the sigma in [8,13] continuation (p is the widest)')
E('SOneRowsMC', 'exact', fmt_of(dec_=0), src='CHK', key='prof_tip_vs_MC.rows_compared', note='stage-C rows compared; the paper says 1,025')
E('SOneRowsMB', 'exact', fmt_of(dec_=0), src='CHK', key='ext813_vs_MB.rows_compared', note='stage-B rows compared; the paper says 501')
E('SOneMinMargin', 'margin', fmt_of(2, style='sci'), expr='a',
  inputs=[('a', 'CHK', 'min:*.min_relative_margin.*')],
  note='smallest relative margin over all three comparisons (min over keys); main_v3 prints 1.3e-9 (nearest); floor is 1.2e-9. R7: check_outA5.json was rewritten by a kit run, hash differs from the closed manifest pin')
E('TwoSegments', 'exact', fmt_of(dec_=0), expr='a', inputs=[('a', 'CONT', 'len:segments')], note='Evans contour segments')
E('TwoInteriorBoxes', 'exact', fmt_of(dec_=0), expr='a', inputs=[('a', 'INT', 'len:')], note='analyticity boxes')
E('TwoContourBoxes', 'exact', fmt_of(dec_=0), expr='a', inputs=[('a', 'ETA', 'len:')], note='boxes along the contour for the eta bound')
E('TwoEtaMax', 'upper', fmt_of(4, style='plain'), expr='a', inputs=[('a', 'ETA', 'max:*.eta')], note='eta <= (docs 0.1163); float values in the cache')
E('TwoLipMax', 'upper', fmt_of(3, style='plain'), expr='a', inputs=[('a', 'ETA', 'max:*.Lip')], note='contraction (docs 0.272); float values in the cache')
E('TwoLambdaZero', 'upper', fmt_of(7, style='plain'), src='CE0', key='Lam', kind='est', note='SDP value (float) of the numerical-range certificate theta=0; verify.py log gives the certified check')
E('TwoLambdaOne', 'upper', fmt_of(7, style='plain'), src='CE1', key='Lam', kind='est', note='SDP value (float), theta = 0.48; status optimal_inaccurate')

# ---------------------------------------------------------------- 6. doc-only entries (flagged)
DOC('SOneDTightLower', 'lower', '2.50080627654', 'S1TIGHT', 'D ∈ [2.50080627654, 2.50080627668]', 'cert', 'README only; no JSON of the enclosure of D')
DOC('SOneDTightUpper', 'upper', '2.50080627668', 'S1TIGHT', 'D ∈ [2.50080627654, 2.50080627668]', 'cert', 'README only')
DOC('SOneKTightLower', 'lower', '1.53982063762', 'S1TIGHT', 'k(0) ∈ [1.53982063762, 1.53982063822]', 'cert', 'README only')
DOC('SOneKTightUpper', 'upper', '1.53982063822', 'S1TIGHT', 'k(0) ∈ [1.53982063762, 1.53982063822]', 'cert', 'README only')
DOC('SOneDCondLower', 'lower', '2.500806276619695', 'S1TIP', 'D = −p(tip) ∈ [2.500806276619695, 2.500806276620049]', 'cert', 'README only (ext_tip)')
DOC('SOneDCondUpper', 'upper', '2.500806276620049', 'S1TIP', 'D = −p(tip) ∈ [2.500806276619695, 2.500806276620049]', 'cert', 'README only (ext_tip)')
DOC('SOneKCondLower', 'lower', '1.539820637871836', 'S1TIP', 'k(0) ∈ [1.539820637871836, 1.5398206378721015]', 'cert', 'README only (ext_tip)')
DOC('SOneKCondUpper', 'upper', '1.5398206378721015', 'S1TIP', 'k(0) ∈ [1.539820637871836, 1.5398206378721015]', 'cert', 'README only (ext_tip)')
DOC('SOneTauTipLower', 'lower', '10.249077651120778', 'S1TIP', 'τ_tip ∈ [10.249077651120778, 10.249077651121361]', 'cert', 'README only (conditioned)')
DOC('SOneTauTipUpper', 'upper', '10.249077651121361', 'S1TIP', 'τ_tip ∈ [10.249077651120778, 10.249077651121361]', 'cert', 'README only (conditioned)')
DOC('SOneStateUncert', 'upper', '4.4e-13', 'S1TIGHT', 'State uncertainty at σ = 8: 4.4e−13', 'cert', 'README (with a unicode minus); the outline quotes |d| <= 4.37e-13 which has no source found')
DOC('FarBoxEprime', 'margin', '6.4287', 'S4RES', '(6.4287, 4.6917, 1.9940, 1.69999995)', 'cert', 'n3boxes prints only a log; unrounded 6.4287 floors to 6.42 (log prints x6.43 nearest)', fmt={'mode': 'sig', 'n': 3, 'style': 'plain'})
DOC('FarBoxXEsecond', 'margin', '4.6917', 'S4RES', '(6.4287, 4.6917, 1.9940, 1.69999995)', 'cert', 'unrounded 4.6917', fmt={'mode': 'sig', 'n': 3, 'style': 'plain'})
DOC('FarBoxXsqEthird', 'margin', '1.9940', 'S4RES', '(6.4287, 4.6917, 1.9940, 1.69999995)', 'cert', 'unrounded 1.9940', fmt={'mode': 'sig', 'n': 3, 'style': 'plain'})
DOC('FarBoxE', 'margin', '1.69999995', 'S4RES', '(6.4287, 4.6917, 1.9940, 1.69999995)', 'cert', 'unrounded 1.69999995 (log prints 1.70)', fmt={'mode': 'sig', 'n': 3, 'style': 'plain'})
DOC('HFourKfourErr', 'upper', '5.60e-3', 'S4RES', '| 5.60e-3 against 0.06 | ×10.7 |', 'cert', 'with-error K4; the amp JSON gives only the base 3.1633e-4 and K4 at Ybar (2.17e-3); the 5.60e-3 matches cone.0.K4 of t9_cone_kch.json (0.0055898) by value only')
DOC('HFourKfourMargin', 'margin', '10.7', 'S4RES', '| 5.60e-3 against 0.06 | ×10.7 |', 'cert', 'see HFourKfourErr (0.06/0.0055898 = 10.73)')
DOC('CoercivityCoefficient', 'upper', '14.8', 'S3V32', r'\norm{e^{\tau\mathcal L}u}_\rho\le14.8\,e^{-0.0213\tau}\norm{u}_\rho', 'cert', 'constant of the modulated semigroup bound; no JSON found (step3_linear_certificate)')
DOC('CoercivityRate', 'lower', '0.0213', 'S3V32', r'\norm{e^{\tau\mathcal L}u}_\rho\le14.8\,e^{-0.0213\tau}\norm{u}_\rho', 'cert', 'decay rate of the semigroup bound; no JSON found')

# ---------------------------------------------------------------- finalise
# regions (definition constants of Table T2)
for nm, v, note in [('RegionRampStart', '14', 'tip/shell interface y = 14 (rescaled)'), ('RegionInterface', '28', 'interface y = 28 of Lemma S(A)'),
                    ('RegionShellEnd', '60', 'shells [14,60]'), ('RegionSegEnd', '350', 'segments (60,350]'), ('RegionAnnEnd', '700', 'outer annulus [350,700]'),
                    ('ZoneIIILow', '253', 'zone III lower end, units of lambda_1'), ('ZoneIIIHigh', '379.5', 'zone III upper end (box_III_380)'),
                    ('ConeCentreMin', '300', 'fixed-L cone family centres >= 300 lambda_1'), ('ConeCentreMax', '380', 'fixed-L cone family centres <= 380'),
                    ('FarCentreMin', '490', 'far family centres >= 490 (v3.3)'), ('ConeL', '50', 'fixed-L cone unit L = 50'), ('ConeX', '3', 'cone X = 3'), ('FarL', '113', 'far unit L = 113'), ('FarX', '2', 'far X = 2')]:
    HAND(nm, v, note + ' (definition, not a computed result)')
def _dedupe():
    cur = {(i['source'], i['key']) for e in ENTRIES if not e.get('auto') and e.get('expr') == 'a' and len(e['inputs']) == 1 for i in e['inputs']}
    keep = []
    for e in ENTRIES:
        if e.get('auto') and (e['inputs'][0]['source'], e['inputs'][0]['key']) in cur:
            NAMES.discard(e['name']); continue
        keep.append(e)
    ENTRIES[:] = keep
_dedupe()
def finalise():
    L = {'note': 'Number ledger of paper_v1 (unit U0). Printed values are exact outward roundings of the raw JSON values; see check_numbers.py.',
         'sources': {}, 'entries': []}
    used = set()
    for e in ENTRIES:
        for i in e['inputs']: used.add(i['source'])
    pins = C.manifest_pins()
    for t in sorted(used): L['sources'][t] = {'path': SRC[t], 'sha256': sha(SRC[t]), 'pin_status': C.pin_status(SRC[t], sha(SRC[t]), pins)}
    D = {t: json.load(open(SRC[t])) for t in used if SRC[t].endswith('.json')}
    for e in ENTRIES:
        e = dict(e)
        if e['kind'] == 'doc':
            e['doc_sha256'] = sha(e['doc_file']); e['sources'] = [{'path': e['doc_file']}]; e['raw'] = e['doc_value']
            e['raw_ledger'] = None
            d = e['direction']
            if e['fmt'] is None: e['fmt'] = {'mode': 'sig', 'n': len(Decimal(e['doc_value']).normalize().as_tuple().digits), 'style': 'sci' if 'e' in e['doc_value'].lower() else 'plain'}; e['printed'] = e['doc_value']
            else: e['printed'] = C.round_fmt(Fr(Decimal(e['doc_value'])), e['fmt'], 'up' if d == 'upper' else 'down')
        else:
            bnd, (lo, hi), raws = C.entry_bound(D, e)
            d = e['direction']
            e['printed'] = C.round_fmt(hi if d == 'upper' else lo, e['fmt'], 'up' if d == 'upper' else ('exact' if d == 'exact' else 'down'))
            e['sources'] = [{'path': SRC[i['source']], 'sha256': L['sources'][i['source']]['sha256'], 'key': i['key']} for i in e['inputs']]
            for i, r in zip(e['inputs'], raws): i['raw'] = r if not isinstance(r, (dict, list)) else json.dumps(r)[:200]
            e['raw'] = e['inputs'][0]['raw'] if e['inputs'] else e['printed']
            e['raw_ledger'] = [json.dumps(r, sort_keys=True) for r in raws]
            e['bound_lo'] = '%.17g' % float(lo); e['bound_hi'] = '%.17g' % float(hi)
        L['entries'].append(e)
    # kind of dependence: entries whose any input is est stay what they were set to
    json.dump(L, open(os.path.join(HERE, 'ledger.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'numbers_paper.tex'), 'w').write(C.macro_text(L))
    print('entries', len(L['entries']), 'sources', len(L['sources']))
finalise()
