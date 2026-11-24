#!/usr/bin/env python3
"""Reviewer's own comparison of run outputs (read-only on the proof tree).  Output compare_runs.json.
 (1) regression: runs/reg_fin12f4 (v4 code, NEW=0) vs step3nl_v33/verifier (fin12f4), JSON-equal apart from timing / new keys
 (2) double run: runs/P vs runs/P_A
 (3) reviewer rerun from a copy: v4_review/rerun/P_R vs runs/P
 (4) box-hypothesis ratios, J constants and the F0*tau_t remark, old/new table
Timing keys removed: secs, runtime_s; 'args' removed for (1) and (3) (paths/hashes differ by construction; hashes of the copied files are compared separately)."""
import json, os, re, hashlib, sys
V4 = '/home/claude/s1spec/step3nl_v4'; V33 = '/home/claude/s1spec/step3nl_v33/verifier'
HERE = os.path.dirname(os.path.abspath(__file__))
def strip(o, drop):
    if isinstance(o, dict): return {k: strip(v, drop) for k, v in o.items() if k not in drop}
    if isinstance(o, list): return [strip(v, drop) for v in o]
    return o
def J(p): return json.load(open(p))
T = ('secs', 'runtime_s')
out = {}
# (1)
reg = {}
for al in ('0.00015', '0.000125'):
    a = J(f'{V33}/rate_table_fin12f4_a{al}.json'); b = J(f'{V4}/runs/reg_fin12f4/rate_table_fin12f4_a{al}.json')
    reg[al] = strip(a, T + ('args', 'box_hypotheses')) == strip(b, T + ('args', 'box_hypotheses'))
    reg[al + '_n_oracle_rows'] = len(b['oracle_rows'])
    reg[al + '_env'] = b.get('args', {}).get('env')
a = J(f'{V33}/modsups_fin12f4/final_alpha_points_5e-8.json'); b = J(f'{V4}/runs/reg_fin12f4/modsups_fin12f4/final_alpha_points_5e-8.json')
reg['final_alpha_points'] = strip(a, T + ('args',)) == strip(b, T + ('args',))
reg['tail_cert_log'] = open(f'{V33}/tail_cert_log.txt').read() == open(f'{V4}/runs/reg_fin12f4/tail_cert_log.txt').read()
out['regression_NEW0_vs_fin12f4'] = reg
# (2), (3)
def cmpdir(A, B, drop):
    r = {}
    for f in ('rate_table_fin13v4_a0.00015.json', 'rate_table_fin13v4_a0.000125.json', 'modsups_fin13v4/final_alpha_points_5e-8.json'):
        pa, pb = os.path.join(A, f), os.path.join(B, f)
        if not (os.path.exists(pa) and os.path.exists(pb)): r[f] = 'missing'; continue
        r[f] = strip(J(pa), drop) == strip(J(pb), drop)
    for f in ('tail_cert_log.txt',):
        pa, pb = os.path.join(A, f), os.path.join(B, f)
        r[f] = os.path.exists(pa) and os.path.exists(pb) and open(pa).read() == open(pb).read()
    mask = lambda t: re.sub(r'runtime [0-9.]+ s', 'runtime T s', re.sub(r'cell [0-9.]+s', 'cell Ts', t))
    for f in ('verify_log_fin13v4_a0.00015.txt', 'verify_log_fin13v4_a0.000125.txt'):
        pa, pb = os.path.join(A, f), os.path.join(B, f)
        r[f] = os.path.exists(pa) and os.path.exists(pb) and mask(open(pa).read()) == mask(open(pb).read())
    return r
out['double_run_P_vs_P_A'] = cmpdir(f'{V4}/runs/P', f'{V4}/runs/P_A', T + ('args',))
PR = os.path.join(HERE, 'rerun', 'P_R')
if os.path.exists(os.path.join(PR, 'modsups_fin13v4', 'final_alpha_points_5e-8.json')):
    out['reviewer_rerun_P_R_vs_P'] = cmpdir(PR, f'{V4}/runs/P', T + ('args',))
    sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
    cp = os.path.join(HERE, 'copy', 'step3nl_v4')
    out['copy_hashes_equal_to_tree'] = {f: sha(os.path.join(cp, f)) == sha(os.path.join(V4, f)) for f in
        ['verifier/' + n for n in ('bb.py', 'closure.py', 'constants.py', 'ival.py', 'mod_sups.py', 'rhs_iv.py', 'tail_cert.py', 'unc.py')] + ['v4cert/out/constants_v4.json', 'cubic_E0fine_0.001_0.01_0.01.json', 'agent_int/planB/planB_ramp.json']}
# (4)
for al in ('0.00015', '0.000125'):
    d = J(f'{V4}/runs/P/rate_table_fin13v4_a{al}.json'); o = J(f'{V33}/rate_table_fin12f4_a{al}.json')
    t = d['tail']; u = d['unc']
    F0 = t['floor_terms']['F0']; tau_t = d['tau_tail_start']
    px = u['parts_Xc']; Xc = u['X_const']
    out['numbers_' + al] = dict(
        box_max_ratio=d['box_hypotheses'], steps=len(d['steps']), peak=d['sqrtL_max'], peak_old=o['sqrtL_max'], rho=t['rho'], rho_old=o['tail']['rho'],
        rate_s_U=t['rate_s_U'], a_star=t['a_star'], sqrt_LF=t['sqrt_LF'], alpha_F=t['alpha_F'], delta_tail=u['delta_tail'], delta_tail_old=o['unc']['delta_tail'],
        J={n: [u[n]['const'], u[n]['slope_per_tau'], u[n]['const'] / o['unc'][n]['const'], u[n]['slope_per_tau'] / o['unc'][n]['slope_per_tau']] for n in ('V0sq', 'V1sq', 'V2sq', 'V4sq')},
        first_cell_U=d['oracle_rows'][0]['U'], first_cell_U_old=o['oracle_rows'][0]['U'],
        F0=F0, tau_t=tau_t, F0_tau_t=F0 * tau_t, X_const=Xc, F0_tau_t_over_Xc=F0 * tau_t / Xc, tail_Lstar_half=px['tail_Lstar_half'],
        Lstar_half_over_F0_tau_t=px['tail_Lstar_half'] / (F0 * tau_t),
        strip_over_sqrtLF=None, env=d['args']['env'], inputs=d['args']['inputs'])
json.dump(out, open(os.path.join(HERE, 'compare_runs.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items()}, indent=1)[:9000])
