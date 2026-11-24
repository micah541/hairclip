"""[kit file] Compare every output the kit produced with the proof tree's (reference/outputs.json: SHA-256 per file, relative to step4/).
Tiers: EXACT (same bytes) | TIMING (equal after dropping cpu/wall/seconds fields) | PROVENANCE (equal after also dropping args/paths/time stamps: only the
bookkeeping differs) | NPZ (arrays equal) | DIFFERENT (numbers differ) | MISSING (not produced by this run).
usage: python verify_outputs.py [--json OUT]   -> prints a table; returns the list for check.py"""
import os, sys, json, hashlib
KIT = os.path.dirname(os.path.abspath(__file__)); S4 = os.path.join(KIT, 'step4')
T1 = lambda k: ('cpu' in k or 'wall' in k or k in ('seconds', 'time', 'started', 'elapsed'))
T2 = lambda k: T1(k) or k in ('args', 'work_dir', 'base_dir', 'path', 'paths', 'poldir', 'argv', 'env', 'runner', 'cmd', 'sources', 'provenance', 'H1_JSON', 'sha256', 'sha256_of_files', 'code_sha256', 'base_sha256', 'consts', 'amp', 'out', 'sources_sha', 'inputs', 'code', 'T6RES', 'YOUT')


def strip(o, drop):
    if isinstance(o, dict): return {k: strip(v, drop) for k, v in o.items() if not drop(k)}
    if isinstance(o, list): return [strip(x, drop) for x in o]
    if isinstance(o, str): return o.replace('/home/claude/s1spec/step4/', '').replace('\\', '/')
    return o


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()


def compare():
    ref = json.load(open(os.path.join(KIT, 'reference', 'outputs.json'))); out = []
    refroot = os.environ.get('KIT_REFROOT')              # optional: the proof tree's step4 directory (enables the TIMING / PROVENANCE / NPZ tiers); without it only EXACT or DIFFERENT
    for rel, h in sorted(ref.items()):
        p = os.path.join(KIT, 's1_tight', *rel.split('/')[3:]) if rel.startswith('../../s1_tight/') else os.path.join(S4, *rel.split('/'))      # [kit fix 6] the proof tree's s1_tight is KIT/s1_tight here
        if not os.path.exists(p): out.append((rel, 'MISSING', '')); continue
        if sha(p) == h: out.append((rel, 'EXACT', '')); continue
        tier = 'DIFFERENT'; note = ''
        rp = os.path.join(refroot, *rel.split('/')) if refroot else ''
        if rp and os.path.exists(rp):
            try:
                if rel.endswith('.json'):
                    a, b = json.load(open(p)), json.load(open(rp))
                    tier = 'TIMING' if strip(a, T1) == strip(b, T1) else 'PROVENANCE' if strip(a, T2) == strip(b, T2) else 'DIFFERENT'
                elif rel.endswith('.npz'):
                    import numpy as np
                    a, b = np.load(p), np.load(rp)
                    tier = 'NPZ' if sorted(a.files) == sorted(b.files) and all((a[k] == b[k]).all() for k in a.files) else 'DIFFERENT'
            except Exception as e: note = repr(e)[:80]
        out.append((rel, tier, note))
    return out


if __name__ == '__main__':
    r = compare(); from collections import Counter
    c = Counter(t for _, t, _ in r); print(dict(c))
    for rel, t, n in r:
        if t not in ('EXACT',): print(f'  {t:10s} {rel} {n}')
    if '--json' in sys.argv: json.dump(r, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)
