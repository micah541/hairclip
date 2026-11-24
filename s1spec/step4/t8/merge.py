"""Assemble T8 certified per-slab constants:  python merge.py [--in constants/slabs] [--out constants/constants.json] [--first a --last b]
Checks that all slabs first..last are present (default: T6's certified range in t6/residuals.json; v2 base: 0..180, v1: 5..187)
and records the T7 inputs used (sources and placeholder flags)."""
import os, sys, json, glob, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', default=os.path.join(HERE, 'constants', 'slabs'))
    ap.add_argument('--out', default=os.path.join(HERE, 'constants', 'constants.json'))
    ap.add_argument('--first', type=int, default=None)
    ap.add_argument('--last', type=int, default=None)
    a = ap.parse_args()
    if a.first is None or a.last is None:
        import cconst
        r6 = json.load(open(os.path.join(cconst.T6DIR, 'residuals.json')))['slabs_range']
        a.first = r6[0] if a.first is None else a.first; a.last = r6[1] if a.last is None else a.last
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(a.inp, 'slab_*.json')))]
    recs = sorted((r for r in recs if a.first <= r['slab'] <= a.last), key=lambda r: r['slab'])
    have = [r['slab'] for r in recs]
    missing = [i for i in range(a.first, a.last + 1) if i not in have]
    import cconst
    t7 = cconst.t7_inputs()
    t7s = {k: v for k, v in t7.items() if k not in ('nVmin', '_cert_status')}
    summ = {}
    for j in range(cconst.JMAX + 1):
        L = [r['lev'][str(j)] for r in recs]
        summ[j] = dict(lam_min=min(l['lam'] for l in L), k_max=max(l['k'][str(j)] for l in L), Dcoef_min=min(l['Dcoef'] for l in L),
                       gam_max=max(l['gam'] for l in L), gr_max=max(l['gr'] for l in L),
                       c_max={i: max(l['c'].get(str(i), 0.0) for l in L) for i in range(j)},
                       cpu=sum(r.get('cpu', 0.0) for r in recs))
    out = dict(design=dict(ramp=[cconst.YA, cconst.YB], dP=cconst.DP, dS=cconst.DS, a=0.01, nu=1e-3, levels=list(range(cconst.JMAX + 1))),
               slabs_range=[a.first, a.last], missing=missing, t7=t7s, summary=summ, per_slab=recs)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    tmp = a.out + '.tmp'
    with open(tmp, 'w') as f: json.dump(out, f)
    os.replace(tmp, a.out)
    print(f'{len(recs)} slabs merged, missing {missing} -> {a.out}')
    for j, s in summ.items():
        print(f"  j={j}: lam_min {s['lam_min']:.4f}  k_jj max {s['k_max']:.4e}  theta_min {s['Dcoef_min']:.3f}  gamma max {s['gam_max']:.3e}")


if __name__ == '__main__':
    main()
