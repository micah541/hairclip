"""Cross-check of the certified norms against T1's quadrature estimates (resid.py), sampled densely in time.
  python crosscheck.py --cert OUTDIR [--slabs 3,5,30,61,90,122,150,183,187] [--nth 41] [--json crosscheck.json]
Slabs without a certified file in OUTDIR are certified first (run_slab).  Certified (T5-literal names) must be >= estimate."""
import os, sys, json, argparse
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
T1DIR = os.path.normpath(os.path.join(HERE, '..', 't1'))
for p in (T1DIR, HERE):                 # HERE first: T1 also has a run.py
    if p in sys.path: sys.path.remove(p)
    sys.path.insert(0, p)
MAP = {'R2': 'R2', 'R1': 'R1', 'Rs': 'Rs', 'RIs': 'RIs', 'RF0': 'RFlit0', 'RF1': 'RFlit1', 'RF2': 'RFlit2', 'RF3': 'RFlit3', 'RF4': 'RFlit4'}


def estimate(i, nth):
    from polstore import load_slab
    import resid
    SP, d = load_slab(os.path.join(T1DIR, 'polished', f'slab_{i:03d}.pkl'))
    th = sorted(set(list(np.linspace(0, 1, nth)) + [0.0] + list(SP.RD.c)))
    rows = [resid.norms(resid.samples(SP, t))[0] for t in th]
    return {k: max(r[k] for r in rows) for k in MAP}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cert', required=True)
    ap.add_argument('--slabs', default='3,5,30,61,90,122,150,183,187')
    ap.add_argument('--nth', type=int, default=41)
    ap.add_argument('--json', default=os.path.join(HERE, 'crosscheck.json'))
    a = ap.parse_args()
    man = json.load(open(os.path.join(T1DIR, 'polished', 'manifest.json')))
    res = []
    for i in map(int, a.slabs.split(',')):
        fj = os.path.join(a.cert, f'slab_{i:03d}.json')
        if not os.path.exists(fj):
            from run import work
            work([i], a.cert)
        c = json.load(open(fj))['norms']
        est = estimate(i, a.nth)
        t1 = man[i]['norms'] or {}
        row = dict(slab=i, k=man[i]['k0'])
        for k, kc in MAP.items():
            row[k] = dict(cert=c[kc], est_dense=est[k], est_T1=t1.get(k), ratio=c[kc] / est[k] if est[k] > 0 else None,
                          ok=c[kc] >= est[k] and (t1.get(k) is None or c[kc] >= t1[k]))
        res.append(row)
        print(f"slab {i:3d} k={row['k']:.4g}: " + ' '.join(f"{k}:{row[k]['ratio']:.3f}{'' if row[k]['ok'] else '!'}" for k in MAP), flush=True)
    allok = all(r[k]['ok'] for r in res for k in MAP)
    rmax = max(r[k]['ratio'] for r in res for k in MAP)
    json.dump(dict(nth=a.nth, all_cert_ge_est=allok, max_ratio=rmax, rows=res), open(a.json, 'w'), indent=1)
    print('all certified >= estimates:', allok, ' max ratio cert/est:', rmax)


if __name__ == '__main__':
    main()
