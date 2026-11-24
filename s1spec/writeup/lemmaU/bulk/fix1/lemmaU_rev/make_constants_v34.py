#!/usr/bin/env python3
"""Writes constants_needed_v34.json (same format as t3_proof/constants_needed.json, plus a 'numbers' list of quoted values
with their outward-rounding direction).  Run: nice -n 19 python3 -I make_constants_v34.py ; then check_keys_v34.py."""
import json, hashlib, os
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
FILES = {
 'A':  B + 'fix2/n12/t2_A_v35/out/A_interface.json',
 'A4': B + 'fix2/n12/t2_A_v35/out/asm_design_N2ref4.json',
 'U':  B + 't4_cert/out/U_cert_v35.json',
 'U1': B + 't4_cert/out/u1_ghat.json',
 'Z':  B + 't5_zone3/out/zone3_interface.json',
 'R':  B + 't5_zone3/out/ramp_U.json',
 'S4': '/home/claude/s1spec/step4/p2/out_v33/FINAL_manifest_v33.json',
}
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
C = []
def row(sym, f, keys, task, used, label='cert'):
    C.append(dict(symbol=sym, file=f, key=' / '.join(keys), producer_task=task, used_in=used, label=label))
T = 'res.tight.'
# ---- A (repaired Lemma S(A), v34)
row('kappa_J', 'A', ['HJW.kappa_J'], '2 (t2_A_v35)', 'Lemma HJ, Lemma D (D2)')
row('E^{A,0}', 'A', [T+'EA0'], '2 (t2_A_v35, G-A1)', 'closure A_j(beta)')
row('A^0 (R-sup)', 'A', [T+'A0_Rsup'], '2 (t2_A_v35, G-A1)', 'closure A_j(beta)')
row('F_g = F/gamma\'', 'A', [T+'Fg'], '2 (t2_A_v35, G-A1)', 'closure A_j(beta)')
row('A rule (R-sup)', 'A', [T+'EA_of_N2'], '2 (t2_A_v35, G-A1)', 'closure A_j(beta)')
row('E_ref4', 'A', [T+'EA_ref4'], '2 (t2_A_v35)', 'domination of Ehat by the 4 N2_ref run (K1, M6)')
row('N2_ref', 'A', ['N2_ref'], '2 (t2_A_v35)', 'Lemma N, M6', 'est')
row('kappa^o_{jk}', 'A', [T+'kappa_o'], '2 (t2_A_v35)', 'Lemma D, closure')
row("gamma'", 'A', [T+'gp'], '2 (t2_A_v35)', 'Lemma D')
row('D^src', 'A', [T+'Dsrc'], '2 (t2_A_v35)', 'Lemma D (D2), C0')
row('D^s (window slope)', 'A', [T+'B0.D_window_slope'], '2 (t2_A_v35, B0)', 'Lemma D (D2), C1')
row('sup_core |D^k d|', 'A', [T+'B0.sup_core', T+'B0'], '2 (t2_A_v35, B0)', "B2', profile difference")
row('rho_J', 'A', [T+'rhoJ'], '2 (t2_A_v35)', 'Lemma D')
row('r_k', 'A', [T+'r_export'], '2 (t2_A_v35)', 'Lemma D')
row('a_c', 'A', [T+'amin_c'], '2 (t2_A_v35)', 'Lemma D')
row('F_k^2', 'A', [T+'F2'], '2 (t2_A_v35)', 'Lemma D')
row('coup_k', 'A', [T+'coup'], '2 (t2_A_v35)', 'Lemma D')
row('C1^0', 'A', [T+'C1_0'], '2 (t2_A_v35)', 'Lemma D')
row('C0 rule', 'A', [T+'C0_rule'], '2 (t2_A_v35)', 'Lemma D')
row('C1 rule', 'A', [T+'C1_rule'], '2 (t2_A_v35)', 'Lemma D')
row('E^A(tau_1)', 'A', [T+'IC', T+'IC.E_tau1'], '2 (t2_A_v35)', 'K6')
row('flux margins at N2_ref (comparison only)', 'A', [T+'core.flux_margin.2'], '2 (t2_A_v35)', 'beta_A paragraph (est outer data)', 'cert at est N2_ref')
row('Nbar rule', 'A', ['lemmaN'], '2', 'Lemma N')
row('segments', 'A', ['segments'], '2', 'K4')
row('A key prefix', 'A', ['res.tight'], '2', 'source list')
# ---- A4 (same assembly at 4 x N2_ref)
row('flux margins beta_A=2 (dominating)', 'A4', [T+'flux.beta.2.margin_FLUX0', T+'flux.beta.2.margin_FLUX2'], '2 (t2_A_v35)', 'beta_A paragraph, K1, table')
row('core margins (dominating)', 'A4', [T+'margins', T+'margins.S0v_vs_0.0153'], '2 (t2_A_v35)', "K1, B2', table")
row('S^v_0, S^v_1 (dominating)', 'A4', [T+'S0', T+'S1'], '2 (t2_A_v35)', "B2', (A) core |z1| (M5), K1")
row('D^c (window constant)', 'A4', [T+'info.0.D_window_const', T+'info.1.D_window_const'], '2 (t2_A_v35)', 'Lemma D (D2)')
# ---- U (Lemma U certificate, v34)
for sym, keys, used in [
  ('abar,bbar,betabar', ['inputs.abar', 'inputs.betabar'], 'setting'),
  ('c_dU', ['launch.c_dU'], 'Lemma dU'), ('w_min', ['launch.wmin'], 'weight, launch, Lemma W'),
  ('C_A^D', ['launch.CA_D'], 'launch (b)'), ('q^D', ['launch.qD'], 'launch (b)'), ('Delta^in', ['launch.Delta_in'], 'launch (b)'),
  ('c_2,c_3', ['launch.c2', 'launch.c3', 'launch.c23_cells', 'launch'], 'launch (a)'),
  ('Phi_j', ['launch.Phi', 'margins.present_rows', 'margins.present_rows.0.Ubar'], 'Step 9, covering'),
  ('Gamma_m(Theta phi)', ['Gamma'], 'Lemma T remark'), ('kappa~, kappa2~', ['kappa', 'kappa2'], 'Lemma T remark'),
  ('Ghat_{jk}', ['Ghat'], 'Step 1'), ('g', ['g'], 'weight'), ('D(56), D(28)', ['D56', 'D28'], 'table'),
  ('theta(28)', ['theta28'], '(U2)'), ('rate margin (R2)', ['rate_margin'], '(U2)'), ('(R<2) margins', ['R_lt2_margin'], '(U2)'),
  ('K_{jk}', ['K'], 'Step 3'), ('footw', ['footw'], 'Step 5'), ('c33', ['c33'], 'Step 4'), ('c223', ['c223'], 'Step 4'),
  ('C_eta', ['C_eta'], 'Step 5'), ('L_max', ['L_max'], 'Step 4'), ('R', ['R_ratio'], 'Step 4'),
  ('sweeps', ['sweeps', 'sweeps.per_yc'], 'Def. sweeps'), ('D_k (late)', ['data_terms_late'], 'Prop. Xterm'),
  ('C0, C1', ['supersolution.C0', 'supersolution.C1'], 'Lemma D'), ('Ehat', ['supersolution.EA'], 'K1'),
  ('Bhat', ['supersolution.Bhat'], 'K2'), ('Nbar^2', ['supersolution.N2', 'supersolution.N2_over_4N2ref', 'supersolution.N2_le_4N2ref'], 'Lemma N, M6'),
  ('K1,K2,K3', ['supersolution.K1_ok', 'supersolution.K2_ok', 'supersolution.K3_ok', 'supersolution.K1_note', 'supersolution.slack_K2', 'supersolution.beta_A'], 'closure'),
  ('M1z checks', ['supersolution.M1z_hyp_check'], '(U1), sec. M1z'), ('box', ['box'], 'sec. zones'),
  ('window margins', ['margins'], 'K3'), ('loop gain', ['loop_gain'], 'Def. loop gain'),
  ('U-fine Bhat', ['ufine', 'ufine.Bhat'], 'Prop. U-fine'), ('S33+foot share', ['S33_foot_over_Bhat2'], 'Step 5')]:
    row(sym, 'U', keys, '4 (t4_cert, M14)', used)
# ---- U1 (u1_ghat.json)
row('Gamma_1(Theta) enclosure', 'U1', ['weight.Theta.Gamma'], '4 (u1_ghat, M4)', 'Lemma T (iii)')
row('phi effect', 'U1', ['phi_effect_rel'], '4 (u1_ghat)', 'Lemma T remark')
row('M1z core check', 'U1', ['M1z_core_check'], '4 (u1_ghat)', 'sec. M1z')
row('outer-foot Y adds', 'U1', ['t5_Yadd'], '4 (u1_ghat)', 'Step 1')
row('level-2 budget slack', 'U1', ['cases.primary.budget_slack'], '4 (u1_ghat)', 'Step 1 reserve (M5)')
row('primary Y_next', 'U1', ['cases.primary.Y_next'], '4 (u1_ghat)', 'Step 8 (M1)')
# ---- Z (zone3 interface)
for sym, keys, used in [
  ('b0,b1,bz (Z1)', ['box_III_380'], 'Step 8'), ('Z2 box', ['box_III_700'], 'Lemma N'), ('n_III', ['n_III', 'n_III_pieces'], 'Lemma N, hyp. (H5)'),
  ('tau_III', ['tau_III'], 'Lemma N'), ('N_k^2 outer ramp', ['DeltaR_detail.N_sq'], 'Step 8'), ('pieces / members', ['DeltaR_detail.pieces'], 'Step 8'),
  ('Delta^R parts', ['DeltaR_detail.parts'], 'Step 8'), ('Delta^{R,b} coefficient', ['DeltaR_detail.DeltaR_b_coefficient'], 'Lemma T (iv)'),
  ('Delta^R', ['DeltaR'], 'Step 8'), ('Delta^{R,b}', ['DeltaR_b'], 'Lemma T (iv)'), ('covering', ['acceptance.covering_asserts', 'acceptance'], 'Step 8, U-fine'),
  ('U-fine K2, K3', ['ufine', 'ufine.K2_at_Bhat_ref', 'ufine.K3_strip_at_Bhat_ref'], 'Prop. U-fine'), ('H6 margins', ['cone_family.H6', 'cone_family.H6.at_Bhat_ref'], 'Prop. U-fine'),
  ('hand-off list', ['handoff_thresholds'], 'hyp. (H3), (H2f)'), ('step-4 reruns', ['step4_reruns'], 'hyp. (H5)'), ('inputs_used', ['inputs_used'], 'Step 8')]:
    row(sym, 'Z', keys, '5 (t5_zone3, M14)', used)
# ---- R (ramp_U.json)
row('K^o_{jk}', 'R', ['K'], '5 (ramp_U, M1)', 'Step 8')
row('Y_next used / contract Y', 'R', ['Y_next_primary_used', 'Y_next_global'], '5 (ramp_U, M1)', 'Step 8')
row('C_xi', 'R', ['CXI'], '5 (ramp_U)', 'Step 8')
row('footw_out, W33', 'R', ['foot_out'], '5 (ramp_U)', 'Steps 1, 8')
row('u1 pin', 'R', ['u1_ghat_sha'], '5 (ramp_U)', 'gaps (O3)')
# ---- S4 (step-4 v33 manifest)
row('(H2) for f - U_1', 'S4', ['results.M9_H2_tightened_fU1'], '6 (step4 v33, B0)', 'Remark C5, K6, hyp.')
row('(H2f)', 'S4', ['results.M9_H2f'], '6 (step4 v33)', 'K6, hyp.')
row('M9 arb rows (final; supersede float rows)', 'S4', ['results.M9_H2f', 'results.M9_H2_tightened_fU1', 'results.M9_H4_far', 'results.M9_H1_centre_map'], '6 (step4 v33, M9)', 'hyp. (H2), (H2f), (H5); gaps (O6)', 'cert')
row('fixed-L cone data', 'S4', ['results.cone_fixedL', 'results.H4_fixedL', 'results.H1_fixedL', 'results.H7_fixedL'], '6 (step4 v33)', 'K6, hyp. (H5)')
# ---- hypotheses / paper (no file), as in the v33 list
C += [
 dict(symbol='H2 (tightened, variable f - U_1)', value='(1e-8,1e-6,1e-3)', producer_task=6, used_in='launch (a), K6', label='hyp (certified in S4 h2u_v33)'),
 dict(symbol='H3', value='(2e-7,1e-5,1e-5)', producer_task=6, used_in='launch (a)', label='hyp'),
 dict(symbol='H2f', value='(1e-6,5e-6,1e-5)', producer_task=6, used_in='K6 (datum of (A))', label='hyp (certified in S4 h2f)'),
 dict(symbol='alpha_1', value='1.5e-4', producer_task=6, used_in='Lemma dU, d', label='hyp'),
 dict(symbol='H4 range', value='graph over the S1 ray on y in [48, 700]', producer_task=6, used_in='G at tau_1', label='hyp'),
 dict(symbol='fixed-L cylinder data (H5, H1, H4, H7)', value='zone3 step4_reruns', producer_task=6, used_in='B4 at tau_1', label='hyp (certified in S4)'),
 dict(symbol="c'_j", value='(0.0520,0.1446,0.4324)', source="step3_proof_v32.tex Prop. B' l.245", used_in='Step 9', label='paper'),
 dict(symbol='(W2) verifier tail lemma', value='tail_cert PASS, closure tail=True (fin12f4)', source='step3nl_v33/V33_VERIFIER.md l.6-7', used_in='Lemma HJ', label='cert (verifier)'),
]
# ---- quoted numbers: [quoted value, file, key path, 'up' (quoted >= value) / 'lo' (quoted <= value)]
N = []
def num(q, f, k, d): N.append([q, f, k, d])
for i, q in enumerate(['1.2018e-6', '1.0441e-4', '0.056694']): num(q, 'U', 'supersolution.EA.%d' % i, 'up')
for i, q in enumerate(['4.8735e-6', '5.5933e-4', '0.36675']): num(q, 'U', 'supersolution.Bhat.%d' % i, 'up')
for i, q in enumerate(['1.1340e-6', '1.6623e-4']): num(q, 'U', 'supersolution.C0.%d' % i, 'up')
for i, q in enumerate(['1.5241e-9', '5.5255e-8']): num(q, 'U', 'supersolution.C1.%d' % i, 'up')
for i, q in enumerate(['0.071472', '0.33012', '5.0406']): num(q, 'U', 'supersolution.N2.%d' % i, 'up')
for i, q in enumerate(['0.18', '0.42', '0.26']): num(q, 'U', 'supersolution.N2_over_4N2ref.%d' % i, 'up')
num('0.28585', 'U', 'supersolution.M1z_hyp_check.zones_I_II.S1_2Bhat', 'up'); num('0.50786', 'U', 'supersolution.M1z_hyp_check.zone_III.b1_plus_pU1', 'up')
num('2.353', 'U', 'margins.late_M1', 'lo'); num('4.892', 'U', 'margins.late_M0', 'lo'); num('12.99', 'U', 'margins.late_M0z', 'lo')
num('2.770', 'U', 'margins.present_M1_min.value', 'lo'); num('2.225', 'U', 'margins.present_M0_min.value', 'lo'); num('14.37', 'U', 'margins.present_M0z_min.value', 'lo')
for i, q in enumerate(['1.23e-7', '0.037169', '0.039035']): num(q, 'U', 'loop_gain.%d.total' % i, 'up')
num('0.037018', 'U', 'loop_gain.2.neighbour', 'up'); num('0.0020172', 'U', 'loop_gain.2.through_A', 'up')
for i, q in enumerate(['4.8723e-6', '2.6266e-4', '0.025842']): num(q, 'U', 'ufine.Bhat.%d' % i, 'up')
num('1.7187e-8', 'U', 'launch.c_dU', 'up'); num('0.73492', 'U', 'launch.wmin', 'lo'); num('2.2429', 'U', 'launch.CA_D', 'up'); num('0.018509', 'U', 'launch.qD', 'up')
for i, q in enumerate(['3.131e-12', '3.287e-12', '3.451e-12']): num(q, 'U', 'launch.Delta_in.%d' % i, 'up')
for i, q in enumerate(['6.434e-9', '1.4217e-8', '3.8404e-8']): num(q, 'U', 'launch.Phi.%d' % i, 'up')
for i, q in enumerate(['4.0000358', '16.000322', '64.00218']): num(q, 'U', 'Gamma.%d' % i, 'up')
num('8.0000358', 'U', 'kappa', 'up'); num('64.00061', 'U', 'kappa2', 'up')
num('51322.6', 'U', 'Ghat.0.0', 'up'); num('33962.9', 'U', 'Ghat.1.1', 'up'); num('348341.1', 'U', 'Ghat.2.2', 'up')
num('4.938', 'U', 'Ghat.1.0', 'up'); num('27.984', 'U', 'Ghat.2.0', 'up'); num('458.53', 'U', 'Ghat.2.1', 'up')
for k, q in {'0.0': '8034.21', '0.1': '724.89', '1.0': '0.08152', '1.1': '6336.94', '1.2': '545.59', '1.3': '0.58700',
             '2.0': '0.87186', '2.1': '528.62', '2.2': '23425.98', '2.3': '555.62'}.items(): num(q, 'U', 'K.' + k, 'up')
num('2.6361', 'U', 'footw.1', 'up'); num('2663.3', 'U', 'footw.2', 'up')
num('108.02', 'U', 'c33', 'up'); num('7.3284', 'U', 'c223', 'up'); num('13.376', 'U', 'C_eta', 'up'); num('3.5352', 'U', 'L_max', 'up'); num('11.593', 'U', 'R_ratio', 'up')
num('5.0892', 'U', 'sweeps.Sigma_in.cover_sup', 'up'); num('0.30846', 'U', 'sweeps.Sigma_E.cover_sup', 'up'); num('1.5243', 'U', 'sweeps.S_X.cover_sup', 'up')
num('3.0081', 'U', 'sweeps.Sigma_X.cover_sup', 'up'); num('4.3350e-3', 'U', 'sweeps.Sigma_rs.2.cover_sup', 'up')
num('1.7330e-6', 'U', 'data_terms_late.D.0', 'up'); num('2.5353e-4', 'U', 'data_terms_late.D.1', 'up')
num('2.7e-3', 'U', 'data_terms_late.slope_share_C1SigmaX_over_Dk.0', 'up'); num('6.6e-4', 'U', 'data_terms_late.slope_share_C1SigmaX_over_Dk.1', 'up')
num('28.231', 'U', 'rate_margin', 'lo'); num('4.432', 'U', 'R_lt2_margin.0', 'lo'); num('6.697', 'U', 'R_lt2_margin.1', 'lo')
num('0.37011', 'U', 'theta28', 'up'); num('1.0368', 'U', 'D56', 'up'); num('1.5876', 'U', 'D28', 'up'); num('1.69e-3', 'U', 'S33_foot_over_Bhat2', 'up')
# A
num('0.042351', 'A', 'HJW.kappa_J', 'up')
for i, q in enumerate(['0.24292', '0.24106', '0.11961']): num(q, 'A', T + 'gp.%d' % i, 'lo')
num('0.97623', 'A', T + 'amin_c', 'lo'); num('5.4995e-7', 'A', T + 'rhoJ', 'up'); num('1.2416e-5', 'A', T + 'kappa_o.2.2', 'up')
for i, q in enumerate(['1.8092e-7', '9.8276e-5']): num(q, 'A', T + 'Dsrc.%d' % i, 'up')
for i, q in enumerate(['8.7399e-10', '5.2599e-8']): num(q, 'A', T + 'B0.D_window_slope.%d' % i, 'up')
for i, q in enumerate(['5.6506e-13', '5.9320e-13']): num(q, 'A', T + 'F2.%d' % i, 'up')
num('1.1568e-12', 'A', T + 'coup.1', 'up')
for i, q in enumerate(['7.2060e-10', '4.3143e-8']): num(q, 'A', T + 'C1_0.%d' % i, 'up')
num('4.23e-5', 'A', T + 'B0.sup_core.0', 'up')
num('3.396', 'A', T + 'core.flux_margin.2.FLUX0_over_K0', 'lo'); num('6.012', 'A', T + 'core.flux_margin.2.FLUX2_over_K2', 'lo')
# A4
num('3.365', 'A4', T + 'flux.beta.2.margin_FLUX0', 'lo'); num('5.955', 'A4', T + 'flux.beta.2.margin_FLUX2', 'lo')
num('7.1878e-3', 'A4', T + 'S0', 'up'); num('0.10430', 'A4', T + 'S1', 'up')
for k, q in [('S0v_vs_0.0153', '2.116'), ('S1v_vs_0.6', '5.748'), ('M1core', '4.409'), ('M0core', '5.450'), ('Ghat22', '1.917'), ('ramp M1', '1.184')]:
    num(q, 'A4', T + 'margins.' + k, 'lo')
for i, q in enumerate(['1.7732e-7', '9.8058e-5']): num(q, 'A4', T + 'info.%d.D_window_const' % i, 'up')
# Z
for i, q in enumerate(['4.4986e-8', '6.7972e-6', '0.011480']): num(q, 'Z', 'DeltaR.%d' % i, 'up')
for i, q in enumerate(['2.66e-12', '3.26e-9', '4.35e-8']): num(q, 'Z', 'DeltaR_b.%d' % i, 'up')
for i, q in enumerate(['8.28e-5', '0.10149', '1.3567', '25350']): num(q, 'Z', 'DeltaR_detail.N_sq.%d' % i, 'up')
num('0.011298', 'Z', 'DeltaR_detail.parts.2.pairs.3', 'up')
num('1.6e-4', 'Z', 'DeltaR_detail.parts.2.L4_foot', 'up'); num('9.3e-6', 'Z', 'DeltaR_detail.parts.2.cubic', 'up'); num('3.7e-7', 'Z', 'DeltaR_detail.parts.2.S33_foot', 'up')
for i, q in enumerate(['0.071445', '0.32707', '3.0445']): num(q, 'Z', 'n_III.%d' % i, 'up')
num('4.0708', 'Z', 'tau_III', 'up'); num('0.10032', 'Z', 'box_III_380.b0', 'up'); num('0.41974', 'Z', 'box_III_380.b1', 'up'); num('0.014291', 'Z', 'box_III_380.bz', 'up')
num('1.495', 'Z', 'box_III_380.margins.b0', 'lo'); num('1.238', 'Z', 'box_III_380.margins.b1', 'lo'); num('11.19', 'Z', 'box_III_380.margins.bz', 'lo')
num('0.22757', 'Z', 'box_III_700.b0', 'up'); num('0.54234', 'Z', 'box_III_700.b1', 'up')
num('0.046956', 'Z', 'ufine.K2_at_Bhat_ref', 'up'); num('0.054140', 'Z', 'ufine.K3_strip_at_Bhat_ref', 'up')
num('1.860', 'Z', 'cone_family.H6.at_Bhat_ref.min_margin', 'lo'); num('2.46', 'Z', 'cone_family.H6.at_Bhat_ref.cone.slope.margin', 'lo')
num('2.54', 'Z', 'cone_family.H6.at_Bhat_ref.cone.K2.margin', 'lo'); num('2.51', 'Z', 'cone_family.H6.at_Bhat_ref.cone.K3.margin', 'lo')
num('2.30', 'Z', 'acceptance.tstar_margin', 'lo')
# R, U1
num('0.0044', 'U1', 't5_Yadd.level1.2', 'up'); num('0.1081', 'U1', 't5_Yadd.level2.2', 'up'); num('0.374', 'U1', 't5_Yadd.level2.3', 'up')
num('4.00000000000000000000054', 'U1', 'weight.Theta.Gamma.0', 'up')  # Gamma_1(Theta) <= 4 + 5.4e-22
for i, q in enumerate(['13436.03', '10526.19', '39482.30']): num(q, 'R', 'K.%d.%d' % (i, i), 'up')
for i, q in enumerate(['457.092', '356.940', '1320.220']): num(q, 'U1', 'cases.primary.Y_next.%d' % i, 'up')
# S4
for j, (b, m) in enumerate([('3.658e-10', '27.3'), ('1.777e-9', '562'), ('1.665e-6', '600')]):
    num(b, 'S4', 'results.M9_H2_tightened_fU1.rows.%d.bound' % j, 'up'); num(m, 'S4', 'results.M9_H2_tightened_fU1.rows.%d.margin' % j, 'lo')
for j, (b, m) in enumerate([('3.1363e-7', '3.188'), ('1.0034e-6', '4.983'), ('2.4351e-6', '4.106')]):
    num(b, 'S4', 'results.M9_H2f.rows.%d.with_error' % j, 'up'); num(m, 'S4', 'results.M9_H2f.rows.%d.margin' % j, 'lo')
num('2.641', 'S4', 'results.cone_fixedL.checks.slope.margin', 'lo'); num('3.037', 'S4', 'results.cone_fixedL.checks.K2.margin', 'lo')
num('27.38', 'S4', 'results.cone_fixedL.checks.N3.margin', 'lo'); num('1119', 'S4', 'results.H4_fixedL.margin', 'lo')
num('2.302', 'S4', 'results.H7_fixedL.margin_tstar_over_ratio', 'lo')
# M9 arb rows (cited in sec:hyp and gaps)
for j, m in enumerate(['27.34', '562.8', '600.9']): num(m, 'S4', 'results.M9_H2_tightened_fU1.rows.%d.margin' % j, 'lo')
for j, m in enumerate(['3.1885', '4.9834', '4.1067']): num(m, 'S4', 'results.M9_H2f.rows.%d.margin' % j, 'lo')
num('15.58', 'S4', 'results.M9_H4_far.min_margin', 'lo')
PIN = {'S4': '8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec'}   # final v33 manifest after N1/N2 and F1/F2 (21 stages, incl. f12_sliver)
out = dict(files={k: dict(path=p, sha256=PIN.get(k, sha(p))) for k, p in FILES.items()}, constants=C, numbers=N,
           note=('keys are full paths from the file root (A, A4 keys start with res.tight. where they live there); array indices as path '
                 'components; ' + "'numbers' lists every quoted bound checked by check_keys_v34.py with its outward direction "
                 "(up: quoted >= certified upper end; lo: quoted <= certified lower end). Labels as in lemmaU_v34.tex."))
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'constants_needed_v34.json')
json.dump(out, open(p + '.tmp', 'w'), indent=1); os.replace(p + '.tmp', p); print('wrote', p, len(C), 'rows,', len(N), 'numbers')
