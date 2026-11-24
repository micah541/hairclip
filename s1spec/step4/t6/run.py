"""T6 production entry point (portable: Linux/Windows, no fork).

  python run.py --slabs 5-187 --out OUTDIR [--jobs N] [--poldir DIR]

Writes OUTDIR/slab_XXX.json (certified residual norms, scalar defects, T7 corrections, diagnostics) and
OUTDIR/jets/jets_XXX.npz (per-element hulls).  Existing slab files are skipped (restartable).
--jobs N starts N independent worker processes (subprocess, one slab sub-list each, interleaved)."""
import os, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, time, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)


def parse_slabs(spec):
    out = []
    for part in spec.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        elif part: out.append(int(part))
    return out


def work(slabs, outdir, poldir=None):
    import numpy as np
    from slab import run_slab
    os.makedirs(os.path.join(outdir, 'jets'), exist_ok=True)
    for i in slabs:
        fj = os.path.join(outdir, f'slab_{i:03d}.json')
        if os.path.exists(fj): continue
        t0 = time.time(); c0 = time.process_time()
        out, jets = run_slab(i, poldir=poldir)
        out['wall'] = time.time() - t0; out['cpu'] = time.process_time() - c0
        np.savez_compressed(os.path.join(outdir, 'jets', f'jets_{i:03d}.npz'), **jets)
        tmp = fj + '.tmp'
        with open(tmp, 'w') as f: json.dump(out, f, indent=1)
        os.replace(tmp, fj)
        print(f'slab {i}: N={out["N"]} cpu={out["cpu"]:.1f}s R2={out["norms"].get("R2", 0):.2e} '
              f'RF4={out["norms"].get("RF4_T5", 0):.2e}', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slabs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--poldir', default=None)
    a = ap.parse_args()
    slabs = parse_slabs(a.slabs)
    if a.jobs <= 1:
        work(slabs, a.out, a.poldir); return
    procs = []
    for w in range(a.jobs):
        sub = slabs[w::a.jobs]
        if not sub: continue
        cmd = [sys.executable, os.path.abspath(__file__), '--slabs', ','.join(map(str, sub)), '--out', a.out]
        if a.poldir: cmd += ['--poldir', a.poldir]
        log = open(os.path.join(a.out, f'worker{w}.log'), 'a') if os.path.isdir(a.out) else None
        if log is None:
            os.makedirs(a.out, exist_ok=True); log = open(os.path.join(a.out, f'worker{w}.log'), 'a')
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
    for p in procs: p.wait()


if __name__ == '__main__':
    main()
