#!/usr/bin/env python3
"""fix4_r3: curated old/new table of the quantities named in the task (forcing, E-hat, c_dU, N-bar, zone III, every named margin).  Read-only.
usage: python3 -B key_table.py <old_snapshot> <new_snapshot>   -> markdown on stdout.  Values: arb midpoints (radii < 1e-20 relative); 'lo'/'up'
as stored; ratio = new/old.  Labels as in the source files (arbio balls of assembly_rep and the num() exports of t4/t5 are certified arb)."""
import sys, os, json, math
OLD, NEW = sys.argv[1:3]
def v(x):
    if isinstance(x, dict):
        if 'arb' in x and isinstance(x['arb'], list): m, e, r = x['arb']; return math.ldexp(float(int(m)), int(e))
        if 'double' in x: return x['double']
        if 'value' in x: return x['value']
    return float(x)
def get(d, path):
    for k in path: d = d[k]
    return d
J = {}
def load(snap, f):
    k = (snap, f)
    if k not in J: J[k] = json.load(open(os.path.join(snap, f)))
    return J[k]
ROWS = []
def row(name, f, path, label, kind=''):
    a, b = v(get(load(OLD, f), path)), v(get(load(NEW, f), path)); ROWS.append((name, a, b, b/a if a else float('nan'), label, f, '/'.join(map(str, path)), kind))
A = 't2_A/asm_design_N2ref.json'; A4 = 't2_A/asm_design_N2ref4.json'; AI = 't2_A/A_interface.json'; U = 't4_cert/U_cert_r3_pass2.json'; Z = 't5_zone3/zone3_interface.json'
for j in range(3): row('F_j^2 (A, [14,700]) j=%d' % j, A, ['res', 'tight', 'info', j, 'F2'], 'cert')
for j in range(3): row('forcing share of E_j, j=%d' % j, A, ['res', 'tight', 'parts', j, 'forcing'], 'cert')
for j in range(3): row('E_j (A, tight, N2ref) j=%d' % j, A, ['res', 'tight', 'E', j], 'cert')
for j in range(3): row('Ehat = EA_ref j=%d' % j, AI, ['res', 'tight', 'EA_ref', j], 'cert')
for j in range(3): row('EA_ref4 (4 N2ref) j=%d' % j, AI, ['res', 'tight', 'EA_ref4', j], 'cert')
for j in range(3): row('EA0 j=%d' % j, AI, ['res', 'tight', 'EA0', j], 'cert')
for j in range(2): row('C0_ref j=%d' % j, AI, ['res', 'tight', 'C0_ref', j], 'cert')
for j in range(2): row('C1_0 j=%d' % j, AI, ['res', 'tight', 'C1_0', j], 'cert')
for j in range(2): row('C1_ref j=%d' % j, AI, ['res', 'tight', 'C1_ref', j], 'cert')
row('S0v (core sup|w\'|)', AI, ['res', 'tight', 'core', 'S0v'], 'cert'); row('S1v', AI, ['res', 'tight', 'core', 'S1v'], 'cert')
row('K0 at E (flux)', AI, ['res', 'tight', 'core', 'K0_at_E'], 'cert'); row('K2 at E (flux)', AI, ['res', 'tight', 'core', 'K2_at_E'], 'cert')
for s in ('tight', 'v31'):
    for k in load(OLD, A)['res'][s]['margins']: row('margin %s (%s, N2ref)' % (k, s), A, ['res', s, 'margins', k], 'cert', 'margin')
    for k in load(OLD, A4)['res'][s]['margins']: row('margin %s (%s, 4 N2ref)' % (k, s), A4, ['res', s, 'margins', k], 'cert', 'margin')
    for bq in ('5/4', '2'):
        for k in ('margin_FLUX0', 'margin_FLUX2'):
            row('flux %s beta=%s (%s, N2ref)' % (k, bq, s), A, ['res', s, 'flux', 'beta', bq, k], 'cert', 'margin')
            row('flux %s beta=%s (%s, 4 N2ref)' % (k, bq, s), A4, ['res', s, 'flux', 'beta', bq, k], 'cert', 'margin')
for f, tag in (('t2_A/asm_lo_N2ref.json', 'lo'), ('t2_A/asm_hi_N2ref.json', 'hi'), ('t2_A/asm_probe_N2ref.json', 'probe'), ('t2_A/asm_design_N20.json', 'N2=0'), ('t2_A/asm_design_pm3ibp.json', 'pm3ibp')):
    for k in ('S0v_vs_0.0153', 'M1core', 'ramp M1', 'Ghat22'): row('margin %s (tight, %s)' % (k, tag), f, ['res', 'tight', 'margins', k], 'cert', 'margin')
    for k in ('margin_FLUX0', 'margin_FLUX2'): row('flux %s beta=2 (tight, %s)' % (k, tag), f, ['res', 'tight', 'flux', 'beta', '2', k], 'cert', 'margin')
row('c_dU', 't4_cert/u34_primary.json', ['launch', 'c_dU'], 'cert')
for j in range(3): row('Delta_in j=%d' % j, 't4_cert/u34_primary.json', ['launch', 'Delta_in', j], 'cert')
for k in range(3): row('Nbar^2 (Lemma N, beta=2) k=%d' % k, U, ['supersolution', 'N2', k], 'cert')
for k in range(3): row('Nbar^2 / (4 N2ref) k=%d' % k, U, ['supersolution', 'N2_over_4N2ref', k], 'cert', 'ratio<=1')
for k in range(3): row('Bhat (supersolution) k=%d' % k, U, ['supersolution', 'Bhat', k], 'cert')
for k in range(3): row('U-fine Bhat k=%d' % k, U, ['ufine', 'Bhat', k], 'cert')
for k in range(3): row('EA (Lemma U fixed point) k=%d' % k, U, ['supersolution', 'EA', k], 'cert')
for k in range(2): row('C0 (U) k=%d' % k, U, ['supersolution', 'C0', k], 'cert'); row('C1 (U) k=%d' % k, U, ['supersolution', 'C1', k], 'cert')
for k in ('late_M1', 'late_M0', 'late_M0z'): row('margin %s (Lemma U)' % k, U, ['margins', k], 'cert', 'margin')
for k in ('present_M1_min', 'present_M0_min', 'present_M0z_min'): row('margin %s (Lemma U)' % k, U, ['margins', k, 'value'], 'cert', 'margin')
for j in range(2): row('R<2 margin j=%d (Lemma U)' % j, U, ['R_lt2_margin', j], 'cert', 'margin')
row('rate_margin (Lemma U)', U, ['rate_margin'], 'cert', 'margin')
for j in range(3): row('loop gain total j=%d' % j, U, ['loop_gain', j, 'total'], 'cert')
row('M1z check: zone III b1 + pU1 (vs M1z 0.534)', U, ['supersolution', 'M1z_hyp_check', 'zone_III', 'b1_plus_pU1'], 'cert')
row('M1z check: S1(2 Bhat) zones I/II', U, ['supersolution', 'M1z_hyp_check', 'zones_I_II', 'S1_2Bhat'], 'cert')
for k in range(3): row('zone III n_III k=%d' % k, Z, ['n_III', k], 'cert')
for k in range(3): row('zone III DeltaR k=%d' % k, Z, ['DeltaR', k], 'cert')
for k in range(3): row('zone III DeltaR_b k=%d' % k, Z, ['DeltaR_b', k], 'cert')
for k in ('b0', 'b1', 'bz'): row('zone III box_380 %s' % k, Z, ['box_III_380', k], 'cert'); row('zone III box_380 margin %s' % k, Z, ['box_III_380', 'margins', k], 'cert', 'margin')
for k in ('b0', 'b1'): row('zone III box_700 %s' % k, Z, ['box_III_700', k], 'cert')
for fam in ('cone', 'far'):
    for q in ('K2', 'K3', 'slope'):
        row('H6 %s %s margin (at Bhat_ref)' % (fam, q), Z, ['cone_family', 'H6', 'at_Bhat_ref', fam, q, 'margin'], 'cert', 'margin')
        row('H6 %s %s margin (at fine est)' % (fam, q), Z, ['cone_family', 'H6', 'at_fine_est', fam, q, 'margin'], 'est', 'margin')
row('zone III tstar margin', Z, ['acceptance', 'tstar_margin'], 'cert', 'margin'); row('zone III box margin min', Z, ['acceptance', 'box_margin_min'], 'cert', 'margin')
row('zone III DeltaR2 (acceptance)', Z, ['acceptance', 'DeltaR2'], 'cert')
print('| quantity | old (v35) | new (R3) | new/old | label |\n|---|---|---|---|---|')
for n, a, b, r, l, f, p, kind in ROWS: print('| %s | %.6g | %.6g | %.6f | [%s] |' % (n, a, b, r, l))
print('\nsources (new snapshot %s): %s' % (NEW, ', '.join(sorted({f for *_, f, p, k in ROWS}))))
m = [x for x in ROWS if x[7] == 'margin']; worst = min(m, key=lambda x: x[3]); low = min(m, key=lambda x: x[2])
print('\nmargins listed: %d; largest relative decrease: %s x%.6f; smallest new margin: %s = %.6g' % (len(m), worst[0], worst[3], low[0], low[2]))
