#!/usr/bin/env python3
"""X8: compare two snapshots snap/<a>, snap/<b>: every JSON equal after removing the 'timing' key; every other file byte-identical,
except *.log (timing lines) and *.time (/usr/bin/time output), which are reported only.  usage: cmp_runs.py a b"""
import sys, os, json, hashlib
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix8_flux2_impl/snap'
a, b = sys.argv[1:3]
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k != 'timing'}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
res = {}; bad = []
for root in ('out', 'texts'):
    for f in sorted(os.listdir(os.path.join(D, a, root))):
        pa, pb = os.path.join(D, a, root, f), os.path.join(D, b, root, f)
        if not os.path.exists(pb): res[root + '/' + f] = 'missing in b'; bad.append(f); continue
        if f.endswith('.time') or f.endswith('.timing.json'): res[root + '/' + f] = 'timing file (ignored)'; continue
        ba, bb = open(pa, 'rb').read(), open(pb, 'rb').read()
        if ba == bb: res[root + '/' + f] = 'byte-identical'; continue
        if f.endswith('.json'):
            if strip(json.loads(ba)) == strip(json.loads(bb)): res[root + '/' + f] = 'identical apart from timing'; continue
        if f.endswith('.log'):
            la = [l for l in ba.decode().splitlines() if not any(t in l for t in ('s)', 's,', 'MB)', 'cell', 'Elapsed', 'sha'))]
            lb = [l for l in bb.decode().splitlines() if not any(t in l for t in ('s)', 's,', 'MB)', 'cell', 'Elapsed', 'sha'))]
            res[root + '/' + f] = 'log: equal apart from timing lines' if la == lb else 'log: DIFFERS'; 
            if la != lb: bad.append(f)
            continue
        res[root + '/' + f] = 'DIFFERS'; bad.append(f)
out = dict(a=a, b=b, files=res, ok=not bad, differing=bad)
json.dump(out, open(os.path.join(D, 'cmp_%s_%s.json' % (a, b)), 'w'), indent=1, sort_keys=True)
for k, v in res.items(): print('%-40s %s' % (k, v))
print('OK' if not bad else 'DIFFER: %s' % bad)
sys.exit(0 if not bad else 1)
