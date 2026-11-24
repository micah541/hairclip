"""B1: T6 residual of the (SE)-corrected base, slab by slab (slabc.run_slab_c).
  python3 -I runc.py --slabs 2-180 [--out out] [--jobs 2] [--poldir ../../t1/polished_v2]
Writes OUT/slab_XXX.json (T6's run.py format, plus se_direct and an args block) and OUT/jets/jets_XXX.npz.
Restartable: a finished slab is skipped only if its args (base file, code digests, opts, env) match the current ones."""
import os, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, time, subprocess, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
S4 = os.path.normpath(os.path.join(HERE, '..', '..'))
T6DIR = os.path.join(S4, 't6')
CODE = [os.path.join(HERE, 'slabc.py')] + [os.path.join(T6DIR, f) for f in ('base.py', 'engine.py', 'tp.py', 'weights.py', 'slab.py')]
ENVK = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
CMP = ('poldir', 'base_sha256', 'code_sha256', 'env')


def sha(f):
    return hashlib.sha256(open(f, 'rb').read()).hexdigest()


def parse_slabs(spec):
    out = []
    for part in spec.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        elif part: out.append(int(part))
    return out


def provenance(i, poldir):
    pk = os.path.join(poldir, f'slab_{i:03d}.pkl')
    return dict(argv=sys.argv[1:], slab=i, poldir=poldir, base_sha256=sha(pk),
                code_sha256={os.path.relpath(f, S4): sha(f) for f in CODE}, env={k: os.environ[k] for k in ENVK if k in os.environ})


def work(slabs, outdir, poldir):
    import numpy as np
    if HERE not in sys.path: sys.path.insert(0, HERE)
    from slabc import run_slab_c
    os.makedirs(os.path.join(outdir, 'jets'), exist_ok=True)
    for i in slabs:
        fj = os.path.join(outdir, f'slab_{i:03d}.json')
        pv = provenance(i, poldir)
        if os.path.exists(fj):
            old = json.load(open(fj)).get('args') or {}
            if all(old.get(k) == pv[k] for k in CMP): continue
        t0 = time.time(); c0 = time.process_time()
        out, jets = run_slab_c(i, poldir=poldir)
        out['wall'] = time.time() - t0; out['cpu'] = time.process_time() - c0
        out['args'] = pv
        np.savez_compressed(os.path.join(outdir, 'jets', f'jets_{i:03d}.npz'), **jets)
        tmp = fj + '.tmp'
        with open(tmp, 'w') as f: json.dump(out, f, indent=1)
        os.replace(tmp, fj)
        se = out['se_direct'][0]
        print(f"slab {i}: cpu={out['cpu']:.1f}s R1={out['norms']['R1']:.2e} RF0={out['norms']['RF0_opt4n']:.2e} "
              f"RF3={out['norms']['RF3_opt4n']:.2e} a2={se['sharp']['a2']:.1e} a2'={se['sharp']['a2p']:.1e} "
              f"ImW^(ell)={max(s['exact_end']['imW'] for s in out['se_direct']):.1e}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slabs', required=True)
    ap.add_argument('--out', default=os.path.join(HERE, 'out'))
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--poldir', default=os.path.join(S4, 't1', 'polished_v2'))
    a = ap.parse_args()
    a.out = os.path.abspath(a.out); a.poldir = os.path.normpath(os.path.abspath(a.poldir))
    slabs = parse_slabs(a.slabs)
    if a.jobs <= 1:
        work(slabs, a.out, a.poldir); return
    os.makedirs(a.out, exist_ok=True)
    procs = []
    for w in range(a.jobs):
        sub = slabs[w::a.jobs]
        if not sub: continue
        cmd = [sys.executable, '-I', os.path.abspath(__file__), '--slabs', ','.join(map(str, sub)), '--out', a.out, '--poldir', a.poldir]
        log = open(os.path.join(a.out, f'worker{w}.log'), 'a')
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
    for p in procs: p.wait()


if __name__ == '__main__':
    main()
