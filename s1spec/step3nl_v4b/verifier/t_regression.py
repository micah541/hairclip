#!/usr/bin/env python3
"""T(i) (fin13v4): regression.  The v4 code with NEW=0 (new terms off, fin12f4 typed constants) must reproduce the fin12f4 production outputs
of step3nl_v33/verifier bit for bit, apart from timing fields (secs, runtime_s) and the keys added in v4 (args, box_hypotheses).
Also checks the two cell values quoted in EI_SPEC §8 (tail anchor U, first cell U).
usage: python3 t_regression.py REGDIR  -> REGDIR/t_regression.json"""
import sys, os, json, hashlib
V33 = '/home/claude/s1spec/step3nl_v33/verifier'; R = sys.argv[1]
DROP = {'secs', 'runtime_s', 'args', 'box_hypotheses'}
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k not in DROP}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
res = {}
for al in ('0.00015', '0.000125'):
    a = strip(json.load(open(os.path.join(V33, f'rate_table_fin12f4_a{al}.json')))); b = strip(json.load(open(os.path.join(R, f'rate_table_fin12f4_a{al}.json'))))
    res[f'rate_table_a{al}_equal'] = (a == b)
    if al == '0.00015':
        res['first_cell_U'] = b['oracle_rows'][0]['U']; res['first_cell_U_spec'] = 0.03034405538749633
        ta = [r for r in b['oracle_rows'] if not r['floors']]
        res['tail_anchor_U'] = ta[-1]['U']; res['tail_anchor_U_spec'] = -0.0026930418625762047
fa = strip(json.load(open(os.path.join(V33, 'modsups_fin12f4', 'final_alpha_points_5e-8.json')))); fb = strip(json.load(open(os.path.join(R, 'modsups_fin12f4', 'final_alpha_points_5e-8.json'))))
res['final_alpha_points_equal'] = (fa == fb)
res['tail_cert_log_equal'] = open(os.path.join(V33, 'tail_cert_log.txt')).read() == open(os.path.join(R, 'tail_cert_log.txt')).read()
res['spec_cells_equal'] = (res['first_cell_U'] == res['first_cell_U_spec'] and res['tail_anchor_U'] == res['tail_anchor_U_spec'])
res['PASS'] = all(v for k, v in res.items() if k.endswith('_equal'))
json.dump(res, open(os.path.join(R, 't_regression.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps(res, indent=1))
