"""Read-only check: every SHA-256 pinned by p2/out/FINAL_manifest.json (stage inputs/outputs at their last write, upstream entries,
script_sha256) against the files now.  Writes nothing; prints a JSON summary; exit code 1 on any mismatch."""
import os, sys, json, hashlib
S4 = '/home/claude/s1spec/step4'; FINAL = f'{S4}/p2/out/FINAL_manifest.json'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
ab = lambda p: p if p.startswith('/') else os.path.join(S4, p)
def verify(path=FINAL):
    m = json.load(open(path)); last = {}
    for P in m['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]
    bad += [u['path'] for u in m['upstream'] if sha(ab(u['path'])) != u['sha256']]
    bad += [q for q, h in m.get('script_sha256', {}).items() if sha(ab(q)) != h]
    return dict(manifest=path, manifest_sha256=sha(path), stages=len(m['stages']), paths=len(last), upstream=len(m['upstream']),
                scripts=len(m.get('script_sha256', {})), total=len(last) + len(m['upstream']), recorded=m['checks'].get('hashes_verified_at_end'),
                mismatches=bad, ok=not bad), m, last
if __name__ == '__main__':
    r, _, _ = verify(); print(json.dumps(r, indent=1)); sys.exit(0 if r['ok'] else 1)
