"""Print the T6_REPORT tables from residuals.json and crosscheck.json."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'residuals.json')))
S = R['summary']
print('| norm | definition | target | max certified (slab) | margin |\n|---|---|---|---|---|')
rows = [('R2', 'T5'), ('Rs', 'T5'), ('R1', 'T5'), ('RIs', 'T5')] + [(f'RF{j}_opt4n', 'Lemma F, opt4n') for j in range(5)] + \
       [(f'RF{j}_T5', 'Lemma F, T5 ramps') for j in range(5)] + [(f'RFlit{j}', 'T5 literal') for j in range(5)]
for nm, d in rows:
    s = S[nm]
    print(f"| {nm} | {d} | {s['target']:.0e} | {s['max']:.2e} ({s['at_slab']}) | {s['margin']:.0f} |")
print()
for nm in ('Rs_chi', 'RIs3') + tuple(f'RFa{j}_opt4n' for j in range(5)) + tuple(f'RFb{j}_opt4n' for j in range(5)):
    print(nm, f"{S[nm]['max']:.2e} ({S[nm]['at_slab']})")
for k, v in S.items():
    if k.startswith('scalar:') or k.startswith('t7:'): print(k, f"{v['max']:.2e} ({v['at_slab']})")
fails = {k: v['fails'] for k, v in S.items() if v.get('fails')}
print('fails:', fails, ' missing:', R['missing'])
print('cpu total %.0f s, mean %.1f s' % (R['cpu_total_s'], R['cpu_per_slab_mean_s']))
# by slab ranges: early vs late
for nm in ('R2', 'RF3_opt4n', 'RF4_opt4n'):
    v = [(r['norms'][nm], r['slab']) for r in R['per_slab']]
    print(nm, 'slabs 5-9:', ['%.1e' % a for a, s in v if s < 10], ' max slabs >= 30: %.2e' % max(a for a, s in v if s >= 30))
p = os.path.join(HERE, 'crosscheck.json')
if os.path.exists(p):
    C = json.load(open(p))
    keys = ['R2', 'Rs', 'R1', 'RIs', 'RF0', 'RF1', 'RF2', 'RF3', 'RF4']
    print('\n| slab | k | ' + ' | '.join(keys) + ' |\n|' + '---|' * (len(keys) + 2))
    for r in C['rows']:
        print(f"| {r['slab']} | {r['k']:.4g} | " + ' | '.join(f"{r[k]['ratio']:.3f}" for k in keys) + ' |')
    print('all cert >= est:', C['all_cert_ge_est'], 'max ratio', C['max_ratio'])
