"""[v2: STORE (default polished_v2), PLOG (default polish_P11v2.log)] fresh-process reload check: python3 verify_store.py SLAB... ; compares recomputed T5 norms (slab_norms(SP,6)) with the in-run NORMS in polish_P11.log"""
import sys, os, json
from polstore import load_slab, enc
from polish_chain import slab_norms, KEYS
inrun = {}
for l in open(os.environ.get('PLOG', 'polish_P11v2.log')):
    if l.startswith('NORMS '): d = json.loads(l[6:]); inrun[d['slab']] = d
for i in map(int, sys.argv[1:]):
    SP, d = load_slab(os.path.join(os.environ.get('STORE', 'polished_v2'), f'slab_{i:03d}.pkl'))
    assert all(enc(v) == t for row, rt in zip(SP.Yn, d['Yn']) for v, t in zip(row, rt))
    em = max(abs(float(v.mid())) for E, Es in SP.residuals() for v in list(E) + list(Es))
    mx, rows, th = slab_norms(SP, 6)
    diff = {k: (mx[k], inrun[i][k]) for k in KEYS}
    print(f'slab {i} k={d["k"]:.2f} n={d["n"]} N={d["N"]} max|E|mid={em:.2e} (in-run {d["rec"]["maxE"]:.2e}); norms bit-identical: {all(a == b for a, b in diff.values())}', flush=True)
    for k in KEYS: print(f'   {k}: reload {mx[k]!r} in-run {inrun[i][k]!r}')
