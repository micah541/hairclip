"""Portable launcher for the T7 per-slab batches (Windows-safe: subprocesses via sys.executable, no fork, os.path only).
  python run_batch.py --script c2r3.py --base ../t1/polished_v2 --jobs 8 --out out_c2r3      (slabs: the whole store, 0-180)
  python run_batch.py --script nl7.py  --base ../t1/polished_v2 --slabs 0-180 --jobs 8 --out out_nl
Each job runs `<script> --slabs LIST --out DIR` once on its share of the slabs (dealt round-robin, so the costs balance).  Restartable: every script skips slabs whose
JSON already exists.  Afterwards: python merge.py"""
import os, sys, json, argparse, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__))


def store_range(base):
    """'first-last' slab of a polished store (its manifest.json); relative paths are taken from t7/"""
    p = base or os.path.join('..', 't1', 'polished')
    p = p if os.path.isabs(p) else os.path.normpath(os.path.join(HERE, p))
    s = [r['slab'] for r in json.load(open(os.path.join(p, 'manifest.json')))]
    return f'{min(s)}-{max(s)}'


def worker(script, slabs, out, extra):
    subprocess.run([sys.executable, os.path.join(HERE, script), '--slabs', ','.join(str(s) for s in slabs), '--out', out] + extra, check=False)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--script', required=True); ap.add_argument('--slabs', default=None, help='a-b (default: the whole store)')
    ap.add_argument('--jobs', type=int, default=8); ap.add_argument('--out', required=True); ap.add_argument('--base', default=None)
    ap.add_argument('extra', nargs='*')
    A = ap.parse_args(); A.slabs = A.slabs or store_range(A.base); a, b = [int(v) for v in A.slabs.split('-')]
    if A.base: A.extra = ['--base', A.base] + A.extra
    out = A.out if os.path.isabs(A.out) else os.path.join(HERE, A.out); os.makedirs(out, exist_ok=True)
    groups = [list(range(a + j, b + 1, A.jobs)) for j in range(A.jobs)]
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    t0 = time.time()
    procs = [ctx.Process(target=worker, args=(A.script, g, out, A.extra)) for g in groups if g]
    for p in procs: p.start()
    for p in procs: p.join()
    print(f'done {A.script} slabs {a}-{b} in {time.time() - t0:.0f}s wall')


if __name__ == '__main__':
    main()
