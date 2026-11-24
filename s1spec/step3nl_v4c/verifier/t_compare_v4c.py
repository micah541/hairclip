"""fin13v4c record checks (read only; writes runs/t_compare_v4c.json):
 (1) regression: runs/reg_fin13v4b (V4C=0, TAG fin13v4b) against the fin13v4b production step3nl_v4b/runs/P: rate tables (both alpha points),
     final_alpha_points, tail_cert log, verify logs; JSON-equal apart from timing keys (runtime_s, secs) and args; logs equal after removing timing.
 (2) double run: runs/P_A against runs/P (fin13v4c), same comparison; and v4bcert/out_A against v4bcert/out (JSON apart from runtime_s/args; pickle bytes);
     production JSON also compared byte for byte apart from runtime_s/secs (args included).
 (3) cert chain reproduction: every fin13v4b cert output in step3nl_v4b/v4bcert/out against v4c v4bcert/out (JSON apart from args/runtime_s; pickle bytes).
 (4) tail_cert and box hypotheses of the production runs (max ratios <= 1, the new CF_* asserts present).
usage: python3 -B t_compare_v4c.py"""
import json, os, re, sys
sys.dont_write_bytecode = True
C = '/home/claude/s1spec/step3nl_v4c'; B = '/home/claude/s1spec/step3nl_v4b'
IGN = {'runtime_s', 'secs', 'args'}
def strip(o, ign=IGN):
    if isinstance(o, dict): return {k: strip(v, ign) for k, v in o.items() if k not in ign}
    if isinstance(o, list): return [strip(v, ign) for v in o]
    return o
def diff(a, b, path=''):
    if type(a) != type(b): return [f'{path}: type {type(a).__name__} != {type(b).__name__}']
    if isinstance(a, dict):
        out = [f'{path}/{k}: missing in second' for k in a if k not in b] + [f'{path}/{k}: missing in first' for k in b if k not in a]
        for k in a:
            if k in b: out += diff(a[k], b[k], f'{path}/{k}')
        return out
    if isinstance(a, list):
        if len(a) != len(b): return [f'{path}: len {len(a)} != {len(b)}']
        out = []
        for i, (x, y) in enumerate(zip(a, b)): out += diff(x, y, f'{path}[{i}]')
        return out
    return [] if a == b else [f'{path}: {a!r} != {b!r}']
TIM = re.compile(r'(cell [0-9.]+s|runtime [0-9.]+ s|\[[0-9]+s\])')
def log_eq(p, q):
    a = [TIM.sub('', l) for l in open(p)]; b = [TIM.sub('', l) for l in open(q)]
    return a == b
def cmp_dir(d1, d2, tag1, tag2, with_args=False):
    res = {}
    for al in ('0.00015', '0.000125'):
        f1 = os.path.join(d1, f'rate_table_{tag1}_a{al}.json'); f2 = os.path.join(d2, f'rate_table_{tag2}_a{al}.json')
        j1, j2 = json.load(open(f1)), json.load(open(f2))
        dd = diff(strip(j1), strip(j2)); res[f'rate_table_a{al}'] = dict(equal=not dd, diffs=dd[:20])
        if with_args:
            dd2 = diff(strip(j1, {'runtime_s', 'secs'}), strip(j2, {'runtime_s', 'secs'})); res[f'rate_table_a{al}_with_args'] = dict(equal=not dd2, diffs=dd2[:20])
        l1 = os.path.join(d1, f'verify_log_{tag1}_a{al}.txt'); l2 = os.path.join(d2, f'verify_log_{tag2}_a{al}.txt')
        res[f'verify_log_a{al}'] = dict(equal=log_eq(l1, l2))
    f1 = os.path.join(d1, f'modsups_{tag1}', 'final_alpha_points_5e-8.json'); f2 = os.path.join(d2, f'modsups_{tag2}', 'final_alpha_points_5e-8.json')
    dd = diff(strip(json.load(open(f1))), strip(json.load(open(f2)))); res['final_alpha_points'] = dict(equal=not dd, diffs=dd[:20])
    res['tail_cert_log'] = dict(equal=open(os.path.join(d1, 'tail_cert_log.txt')).read() == open(os.path.join(d2, 'tail_cert_log.txt')).read())
    return res
def cmp_cert(d1, d2, only_common=False):
    cert = {}
    for f in sorted(os.listdir(d2)):
        p1 = os.path.join(d1, f); p2 = os.path.join(d2, f)
        if not os.path.exists(p1):
            if not only_common: cert[f] = dict(equal=False, note='missing in first')
            continue
        if f.endswith('.json'):
            dd = diff(strip(json.load(open(p1))), strip(json.load(open(p2)))); cert[f] = dict(equal=not dd, byte_identical=open(p1, 'rb').read() == open(p2, 'rb').read(), diffs=dd[:10])
        elif f.endswith('.pkl'): cert[f] = dict(byte_identical=open(p1, 'rb').read() == open(p2, 'rb').read())
        elif f.endswith('.log'): cert[f] = dict(equal_wo_timing=log_eq(p1, p2))
    return cert
out = {}
out['regression_vs_fin13v4b'] = cmp_dir(os.path.join(C, 'runs', 'reg_fin13v4b'), os.path.join(B, 'runs', 'P'), 'fin13v4b', 'fin13v4b')
out['double_run_fin13v4c'] = cmp_dir(os.path.join(C, 'runs', 'P_A'), os.path.join(C, 'runs', 'P'), 'fin13v4c', 'fin13v4c', with_args=False)
out['double_run_cert_chain'] = cmp_cert(os.path.join(C, 'v4bcert', 'out_A'), os.path.join(C, 'v4bcert', 'out'))
out['cert_chain_vs_fin13v4b'] = cmp_cert(os.path.join(C, 'v4bcert', 'out'), os.path.join(B, 'v4bcert', 'out'), only_common=True)
prod = {}
for al in ('0.00015', '0.000125'):
    rt = json.load(open(os.path.join(C, 'runs', 'P', f'rate_table_fin13v4c_a{al}.json')))
    mr = rt['box_hypotheses']['max_ratio']
    prod[al] = dict(tail_ok=rt.get('tail') is not None, box_max_ratio=mr, box_all_le_1=all(v <= 1.0 for v in mr.values()),
                    new_asserts_present=all(k in mr for k in ('CF_w_16_31', 'CF_w1_16_31', 'CF_alpha')), env_V4C=rt['args']['env'].get('V4C'))
tl = open(os.path.join(C, 'runs', 'P', 'tail_cert_log.txt')).read()
prod['tail_cert_log'] = tl
out['production_checks'] = prod
def allok(o):
    if isinstance(o, dict):
        if 'equal' in o or 'byte_identical' in o or 'equal_wo_timing' in o:
            return o.get('equal', True) and o.get('equal_wo_timing', True) and (o.get('byte_identical', True) if 'equal' not in o else True)
        if 'box_all_le_1' in o: return o['box_all_le_1'] and o['tail_ok'] and o['new_asserts_present'] and o['env_V4C'] == '1'
        return all(allok(v) for k, v in o.items() if k != 'tail_cert_log')
    return o is True
out['verdict'] = {k: allok(v) for k, v in out.items()}
json.dump(out, open(os.path.join(C, 'runs', 't_compare_v4c.json'), 'w'), indent=1)
print(json.dumps(out['verdict'], indent=1))
for k, v in out.items():
    if k != 'verdict' and not out['verdict'][k]: print(k, json.dumps(v, indent=1)[:3000])
