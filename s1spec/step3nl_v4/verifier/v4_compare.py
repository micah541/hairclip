#!/usr/bin/env python3
"""fin12f4 (step3nl_v33) vs fin13v4 (step3nl_v4) comparison table for V4_REPORT.md.
usage: python3 v4_compare.py RUNDIR_V4 [OUT.json]   (RUNDIR_V4 contains rate_table_fin13v4_*.json and modsups_fin13v4/final_alpha_points_5e-8.json)"""
import sys, os, json
V33 = '/home/claude/s1spec/step3nl_v33/verifier'
R4 = sys.argv[1]
def load(d, tag, al): return json.load(open(os.path.join(d, f'rate_table_{tag}_a{al}.json')))
out = {}
fa_old = json.load(open(os.path.join(V33, 'modsups_fin12f4', 'final_alpha_points_5e-8.json')))
fa_new = json.load(open(os.path.join(R4, 'modsups_fin13v4', 'final_alpha_points_5e-8.json')))
for key, al in (('1.5e-4', '0.00015'), ('1.25e-4', '0.000125')):
    o = load(V33, 'fin12f4', al); n = load(R4, 'fin13v4', al)
    rows = {}
    def put(nm, a, b): rows[nm] = dict(fin12f4=a, fin13v4=b, ratio=(b / a if a else None))
    put('steps', len(o['steps']), len(n['steps']))
    put('peak sqrt(Lambda)', o['sqrtL_max'], n['sqrtL_max'])
    for k in ('rho', 'rate_s_U', 'a_star', 'L_star', 'sqrt_LF', 'alpha_F', 'kappa'): put('tail ' + k, o['tail'][k], n['tail'][k])
    put('tail F0', o['tail']['floor_terms']['F0'], n['tail']['floor_terms']['F0'])
    put('alpha-condition margin a*/alpha_F', o['tail']['a_star'] / o['tail']['alpha_F'], n['tail']['a_star'] / n['tail']['alpha_F'])
    put('tail rho / rho_min(0.002)', o['tail']['rho'] / 0.002, n['tail']['rho'] / 0.002)
    put('V1 delta_t', o['unc']['delta_tail'], n['unc']['delta_tail'])
    for nm in ('V0sq', 'V1sq', 'V2sq', 'V4sq'):
        put(f'J {nm} const', o['unc'][nm]['const'], n['unc'][nm]['const']); put(f'J {nm} slope', o['unc'][nm]['slope_per_tau'], n['unc'][nm]['slope_per_tau'])
    for k in ('A_max', 'W1_max', 'sI1', 'sI2', 'M1_ramp', 'far_w', 'far_w1'): put('far ' + k, o['far_max'][k], n['far_max'][k])
    put('S_ramp', o['S_ramp_max_along_trajectory'], n['S_ramp_max_along_trajectory'])
    put('sup|(a,b)|', fa_old[key]['sup_ab'], fa_new[key]['sup_ab']); put('sup|beta|', fa_old[key]['sup_beta'], fa_new[key]['sup_beta'])
    put('first cell rate U', o['oracle_rows'][0]['U'], n['oracle_rows'][0]['U'])
    put('max step rate U', max(r['rate_U'] for r in o['steps']), max(r['rate_U'] for r in n['steps']))
    rows['box_hypotheses_max_ratio'] = n.get('box_hypotheses')
    out[key] = rows
if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1, sort_keys=True)
for key, rows in out.items():
    print(f"## alpha1 = {key}")
    print("| quantity | fin12f4 | fin13v4 | ratio |\n|---|---|---|---|")
    for nm, r in rows.items():
        if nm == 'box_hypotheses_max_ratio': continue
        f = lambda v: f"{v:.7g}" if isinstance(v, float) else str(v)
        print(f"| {nm} | {f(r['fin12f4'])} | {f(r['fin13v4'])} | {('%.6f' % r['ratio']) if r['ratio'] is not None else '-'} |")
    print("box hypotheses (max value/bound along the run):", json.dumps(rows['box_hypotheses_max_ratio']))
