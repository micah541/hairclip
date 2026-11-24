"""fin13v4b record checks (read only; writes runs/t_compare_v4b.json):
 (1) regression: runs/reg_fin13v4 (V4B=0, TAG fin13v4) against the closed fin13v4 production step3nl_v4/runs/P: rate tables (both alpha points),
     final_alpha_points, tail_cert log, verify logs; JSON-equal apart from timing keys (runtime_s, secs) and args; logs equal after removing timing.
 (2) double run: runs/P_A against runs/P (fin13v4b), same comparison; and v4bcert/out_A against v4bcert/out (JSON apart from runtime_s; pickle bytes).
usage: python3 -B t_compare_v4b.py"""
import json, os, re, sys, hashlib
sys.dont_write_bytecode = True
B = '/home/claude/s1spec/step3nl_v4b'; V4 = '/home/claude/s1spec/step3nl_v4'
IGN = {'runtime_s', 'secs', 'args'}
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k not in IGN}
    if isinstance(o, list): return [strip(v) for v in o]
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
def cmp_dir(d1, d2, tag1, tag2):
    res = {}
    for al in ('0.00015', '0.000125'):
        f1 = os.path.join(d1, f'rate_table_{tag1}_a{al}.json'); f2 = os.path.join(d2, f'rate_table_{tag2}_a{al}.json')
        dd = diff(strip(json.load(open(f1))), strip(json.load(open(f2))))
        res[f'rate_table_a{al}'] = dict(equal=not dd, diffs=dd[:20])
        l1 = os.path.join(d1, f'verify_log_{tag1}_a{al}.txt'); l2 = os.path.join(d2, f'verify_log_{tag2}_a{al}.txt')
        res[f'verify_log_a{al}'] = dict(equal=log_eq(l1, l2))
    f1 = os.path.join(d1, f'modsups_{tag1}', 'final_alpha_points_5e-8.json'); f2 = os.path.join(d2, f'modsups_{tag2}', 'final_alpha_points_5e-8.json')
    dd = diff(strip(json.load(open(f1))), strip(json.load(open(f2)))); res['final_alpha_points'] = dict(equal=not dd, diffs=dd[:20])
    res['tail_cert_log'] = dict(equal=open(os.path.join(d1, 'tail_cert_log.txt')).read() == open(os.path.join(d2, 'tail_cert_log.txt')).read())
    return res
out = {}
out['regression_vs_fin13v4'] = cmp_dir(os.path.join(B, 'runs', 'reg_fin13v4'), os.path.join(V4, 'runs', 'P'), 'fin13v4', 'fin13v4')
out['double_run_fin13v4b'] = cmp_dir(os.path.join(B, 'runs', 'P_A'), os.path.join(B, 'runs', 'P'), 'fin13v4b', 'fin13v4b')
cert = {}
for f in sorted(os.listdir(os.path.join(B, 'v4bcert', 'out'))):
    p1 = os.path.join(B, 'v4bcert', 'out_A', f); p2 = os.path.join(B, 'v4bcert', 'out', f)
    if not os.path.exists(p1): cert[f] = 'missing in out_A'; continue
    if f.endswith('.json'):
        dd = diff(strip(json.load(open(p1))), strip(json.load(open(p2)))); cert[f] = dict(equal=not dd, byte_identical=open(p1, 'rb').read() == open(p2, 'rb').read(), diffs=dd[:10])
    elif f.endswith('.pkl'): cert[f] = dict(byte_identical=open(p1, 'rb').read() == open(p2, 'rb').read())
    elif f.endswith('.log'): cert[f] = dict(equal_wo_timing=log_eq(p1, p2))
out['double_run_cert_chain'] = cert
def allok(o):
    if isinstance(o, dict):
        if 'equal' in o or 'byte_identical' in o or 'equal_wo_timing' in o:
            return o.get('equal', True) and o.get('equal_wo_timing', True) and (o.get('byte_identical', True) if 'equal' not in o else True)
        return all(allok(v) for v in o.values())
    return o is True
out['verdict'] = {k: allok(v) for k, v in out.items()}
json.dump(out, open(os.path.join(B, 'runs', 't_compare_v4b.json'), 'w'), indent=1)
print(json.dumps(out['verdict'], indent=1))
for k, v in out.items():
    if k != 'verdict' and not out['verdict'][k]: print(k, json.dumps(v, indent=1)[:3000])
