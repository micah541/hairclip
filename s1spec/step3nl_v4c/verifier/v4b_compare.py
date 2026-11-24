"""Old (fin13v4, step3nl_v4/runs/P) vs new (fin13v4b, step3nl_v4b/runs/P) numbers for V4B_REPORT.md (read only) -> runs/v4b_compare.{json,md}"""
import json, os, sys
sys.dont_write_bytecode = True
B = '/home/claude/s1spec/step3nl_v4c'; V4 = '/home/claude/s1spec/step3nl_v4'
O = (os.path.join(V4, 'runs', 'P'), 'fin13v4'); N = (os.path.join(B, 'runs', 'P'), 'fin13v4b')
def load(d, tag):
    fa = json.load(open(os.path.join(d, f'modsups_{tag}', 'final_alpha_points_5e-8.json')))
    rt = {k: json.load(open(os.path.join(d, f'rate_table_{tag}_a{a}.json'))) for k, a in (('1.5e-4', '0.00015'), ('1.25e-4', '0.000125'))}
    return fa, rt
fo, ro = load(*O); fn, rn = load(*N)
rows = []; out = {}
for key in ('1.5e-4', '1.25e-4'):
    a, b = fo[key], fn[key]; ta, tb = ro[key]['tail'], rn[key]['tail']
    q = [('steps to the tail', a['steps'], b['steps']), ('peak sqrt(Lambda)', a['peak'], b['peak']), ('tail rho', a['rho'], b['rho']),
         ('a*', a['a_star'], b['a_star']), ('sqrt(L_F)', a['sqrtLF'], b['sqrtLF']), ('alpha_F', a['alphaF'], b['alphaF']),
         ('first cell rate U', ro[key]['steps'][0]['rate_U'], rn[key]['steps'][0]['rate_U']), ('tail anchor U', ta['rate_s_U'], tb['rate_s_U']),
         ('V1 delta_t', a['delta_tail'], b['delta_tail'])]
    for n in ('V0sq', 'V1sq', 'V2sq', 'V4sq'):
        q.append((f'J{n[1]} const', a['J'][n][0], b['J'][n][0])); q.append((f'J{n[1]} slope', a['J'][n][1], b['J'][n][1]))
    q += [('sup|(a,b)|', a['sup_ab'], b['sup_ab']), ('sup|beta|', a['sup_beta'], b['sup_beta'])]
    out[key] = [dict(quantity=n, fin13v4=x, fin13v4b=y, ratio=(y / x if isinstance(x, (int, float)) and x not in (0, None) and isinstance(y, (int, float)) else None)) for n, x, y in q]
    out[key + '_tail_alpha_condition_whole_tail'] = tb.get('alpha_condition_whole_tail')
    out[key + '_box_max_ratio'] = rn[key]['box_hypotheses']['max_ratio']
json.dump(out, open(os.path.join(B, 'runs', 'v4b_compare.json'), 'w'), indent=1)
md = []
for key in ('1.5e-4', '1.25e-4'):
    md.append(f'\n### alpha1 = {key}\n\n| quantity | fin13v4 | fin13v4b | ratio |\n|---|---|---|---|')
    for r in out[key]:
        rr = f"{r['ratio']:.6f}" if r['ratio'] is not None else ''
        md.append(f"| {r['quantity']} | {r['fin13v4']:.7g} | {r['fin13v4b']:.7g} | {rr} |" if isinstance(r['fin13v4'], float) else f"| {r['quantity']} | {r['fin13v4']} | {r['fin13v4b']} | {rr} |")
open(os.path.join(B, 'runs', 'v4b_compare.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
