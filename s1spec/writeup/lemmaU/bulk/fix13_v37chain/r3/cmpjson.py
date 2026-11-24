#!/usr/bin/env python3
"""fix4_r3: compare JSON outputs key by key, ignoring the record fields (args blocks, timing).  Read-only.
usage: python3 -B cmpjson.py [--ignore k1,k2] [--show N] A1 B1 [A2 B2 ...]     prints, per pair, IDENTICAL or the differing paths."""
import sys, json
IGN = {'args', 'cpu', 'seconds', 'wall', 'wall_s', 'elapsed', 'time'}
SHOW = 12; argv = sys.argv[1:]
while argv and argv[0].startswith('--'):
    if argv[0] == '--ignore': IGN |= set(argv[1].split(',')); argv = argv[2:]
    elif argv[0] == '--show': SHOW = int(argv[1]); argv = argv[2:]
def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k in IGN: continue
            if k not in a or k not in b: out.append((path + '/' + str(k), 'missing in ' + ('A' if k not in a else 'B'), '')); continue
            walk(a[k], b[k], path + '/' + str(k), out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append((path, 'len %d' % len(a), 'len %d' % len(b))); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, '%s[%d]' % (path, i), out)
    elif a != b: out.append((path, a, b))
bad = 0
for i in range(0, len(argv), 2):
    A, B = argv[i], argv[i + 1]; out = []
    walk(json.load(open(A)), json.load(open(B)), '', out)
    print(('IDENTICAL  ' if not out else 'DIFF %4d  ' % len(out)) + A + '  vs  ' + B)
    for d in out[:SHOW]: print('    %s\n        A: %s\n        B: %s' % (d[0], str(d[1])[:160], str(d[2])[:160]))
    bad += bool(out)
print('pairs with differences:', bad, 'of', len(argv)//2)
