#!/usr/bin/env python3
"""fin13v4: compare pass A and pass B of the final chain (certification outputs v4cert/out_A vs v4cert/out, production runs/P_A vs runs/P).
Files must be byte-identical apart from timing fields: JSON keys 'secs', 'runtime_s' are removed recursively; in text outputs the
'cell <t>s' and 'runtime <t> s' fields are masked.  *.time and *.log files of /usr/bin/time and progress logs are not compared (timing only).
usage: python3 t_double_run.py -> /home/claude/s1spec/step3nl_v4c/runs/t_double_run.json"""
import os, json, re, hashlib
V4 = '/home/claude/s1spec/step3nl_v4c'
PAIRS = [(os.path.join(V4, 'v4cert/out_A'), os.path.join(V4, 'v4cert/out')), (os.path.join(V4, 'runs/P_A'), os.path.join(V4, 'runs/P'))]
SKIP = re.compile(r'\.(time|log)$|chain_status\.txt$')
OUTP = (os.path.join(V4, 'v4cert/out/'), os.path.join(V4, 'runs/P/'))
def strip(o):
    # sha256 pins of upstream chain outputs differ between passes only because those files carry their own runtime_s; each such file is
    # itself compared below (equal apart from timing), so these pins are masked here.
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k not in ('secs', 'runtime_s') and not (isinstance(k, str) and k.startswith(OUTP))}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
def mask(t): return re.sub(r'runtime [0-9.]+ s', 'runtime T s', re.sub(r'cell [0-9.]+s', 'cell Ts', t))
res = {}; ok = True
for A, B in PAIRS:
    fa = sorted(os.path.relpath(os.path.join(r, f), A) for r, _, fs in os.walk(A) for f in fs)
    fb = sorted(os.path.relpath(os.path.join(r, f), B) for r, _, fs in os.walk(B) for f in fs)
    for f in sorted(set(fa) | set(fb)):
        if SKIP.search(f): continue
        pa, pb = os.path.join(A, f), os.path.join(B, f)
        if not (os.path.exists(pa) and os.path.exists(pb)): res[f'{B}:{f}'] = 'missing'; ok = False; continue
        ba, bb = open(pa, 'rb').read(), open(pb, 'rb').read()
        if ba == bb: res[f'{B}:{f}'] = 'identical'; continue
        if f.endswith('.json'):
            eq = strip(json.loads(ba)) == strip(json.loads(bb)); res[f'{B}:{f}'] = 'equal apart from timing' if eq else 'DIFFERENT'
        else:
            eq = mask(ba.decode()) == mask(bb.decode()); res[f'{B}:{f}'] = 'equal apart from timing' if eq else 'DIFFERENT'
        ok &= eq
out = dict(PASS=ok, files=res)
json.dump(out, open(os.path.join(V4, 'runs', 't_double_run.json'), 'w'), indent=1, sort_keys=True)
print('PASS' if ok else 'FAIL', sum(v == 'identical' for v in res.values()), 'identical,', sum(v == 'equal apart from timing' for v in res.values()), 'equal apart from timing,',
      [k for k, v in res.items() if v not in ('identical', 'equal apart from timing')])
