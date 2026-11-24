"""[kit file] Compare the outputs of the v3.3 chain of this run (step4/p2/out_v33/...) with the proof tree's (reference/v33_out/..., the 173 non-log outputs of
FINAL_manifest_v33.json).  Tiers:
  EXACT      same bytes
  ARGS       equal after dropping every `args` block and the timing fields (cpu, wall, seconds, ...) and after reading absolute file names modulo the install root
  NUMERIC    same structure and same text, only numbers differ (largest relative difference reported): expected for outputs that depend on the regenerated v3.2 chain
             (the kit's v3.2 files differ from the proof tree's in the last digits, README 'Differences')
  DIFFERENT  structure or text differs
  MISSING    not produced by this run
usage: python verify_v33.py [--json OUT]   -> prints a table; compare() is used by check.py"""
import os, sys, json, re, hashlib, math
KIT = os.path.dirname(os.path.abspath(__file__)); S4 = os.path.join(KIT, 'step4'); REFD = os.path.join(KIT, 'reference', 'v33_out')
TIMING = lambda k: ('cpu' in k or 'wall' in k or k in ('seconds', 'time', 'started', 'elapsed'))
_NPR = re.compile(r'^(?:[A-Za-z]:)?/(?:[^/\s]+/)*?(?=(?:step4|writeup|step3nl_v33|step3nl|step3|lin|s1_cap|s1_tight|reference)/)')
_NUM = re.compile(r'^\s*\[?\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)')


def normp(s):
    if not isinstance(s, str): return s
    s = s.replace('\\', '/')
    if not re.match(r'(?:[A-Za-z]:)?/', s): return s
    return _NPR.sub('@/', s).replace('@/reference/v33_out/', '@/step4/')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()


def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items() if k != 'args' and not TIMING(k)}
    if isinstance(o, list): return [clean(x) for x in o]
    return normp(o)


def walk(a, b, acc, path=''):
    """acc = [n_numeric_diffs, max_rel, structural_problems]"""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b): acc[2].append(path + ' keys'); return
        for k in a: walk(a[k], b[k], acc, path + '/' + str(k))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): acc[2].append(path + ' length'); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, acc, path + '/' + str(i))
    elif isinstance(a, bool) or isinstance(b, bool):
        if a != b: acc[2].append(path + ' bool')
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a != b:
            acc[0] += 1; d = abs(a - b) / max(abs(a), abs(b), 1e-300); acc[1] = max(acc[1], d if math.isfinite(d) else math.inf)
    elif isinstance(a, str) and isinstance(b, str):
        if a != b:
            ma, mb = _NUM.match(a), _NUM.match(b)
            if ma and mb:
                x, y = float(ma.group(1)), float(mb.group(1)); acc[0] += 1; acc[1] = max(acc[1], abs(x - y) / max(abs(x), abs(y), 1e-300))
            else: acc[2].append(path + ' text')
    elif a != b: acc[2].append(path + ' type')


def compare():
    out = []
    if not os.path.isdir(REFD): return out
    for dp, dn, fn in sorted(os.walk(REFD)):
        for f in sorted(fn):
            rp = os.path.join(dp, f); rel = os.path.relpath(rp, REFD).replace(os.sep, '/')
            if rel.endswith('FINAL_manifest_v33.json'): continue                    # the manifest itself is compared separately (stage list, results)
            p = os.path.join(S4, *rel.split('/'))
            if not os.path.exists(p): out.append((rel, 'MISSING', '')); continue
            if sha(p) == sha(rp): out.append((rel, 'EXACT', '')); continue
            if not rel.endswith('.json'): out.append((rel, 'DIFFERENT', 'not JSON')); continue
            try:
                with open(p, encoding='utf-8') as f1, open(rp, encoding='utf-8') as f2: a, b = clean(json.load(f1)), clean(json.load(f2))
            except Exception as e: out.append((rel, 'DIFFERENT', repr(e)[:60])); continue
            if a == b: out.append((rel, 'ARGS', '')); continue
            acc = [0, 0.0, []]; walk(a, b, acc)
            out.append((rel, 'NUMERIC' if not acc[2] else 'DIFFERENT', f'{acc[0]} numbers differ, max rel {acc[1]:.2e}' if not acc[2] else '; '.join(acc[2][:3])))
    return out


if __name__ == '__main__':
    from collections import Counter
    r = compare(); print(dict(Counter(t for _, t, _ in r)))
    for rel, t, n in r:
        if t != 'EXACT': print(f'  {t:9s} {rel}  {n}')
    if '--json' in sys.argv: json.dump(r, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)
