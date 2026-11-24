"""Read-only: verify the closed v3.2 manifest (verify_v32_manifest.verify) and accept a mismatch only for paths listed in
EXCEPTION_v32_*.json whose current sha256 equals the recorded current_sha256.  Writes nothing; exit 1 on any other mismatch.
usage: python3 -I -B verify_v32_with_exceptions.py"""
import os, sys, json, glob, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from verify_v32_manifest import verify
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
r, m, last = verify()
exc = {}
for p in sorted(glob.glob(os.path.join(HERE, 'EXCEPTION_v32_*.json'))):
    e = json.load(open(p)); exc[e['path']] = (e, os.path.basename(p))
accepted = [q for q in r['mismatches'] if q in exc and sha(exc[q][0]['abs_path']) == exc[q][0]['current_sha256'] and last.get(q) == exc[q][0]['pinned_sha256']]
other = [q for q in r['mismatches'] if q not in accepted]
out = dict(manifest=r['manifest'], manifest_sha256=r['manifest_sha256'], total=r['total'], mismatches=r['mismatches'],
           accepted_by_exception={q: exc[q][1] for q in accepted}, other_mismatches=other, ok=not other)
print(json.dumps(out, indent=1)); sys.exit(0 if out['ok'] else 1)
