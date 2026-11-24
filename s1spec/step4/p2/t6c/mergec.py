"""B1: assemble the corrected-base T6 run (runc.py) into residuals.json and jets/ in T6's format, rounded outward:
norms, norms_int, t7 and sup_ scalars up by 2 ulps (T6 merge.py's _round_up); range_ pairs widened by 2 ulps (as p2/se_corr.py l.53).
  python3 -I mergec.py [--in out] [--first 2] [--last 180] [--json residuals.json] [--jets jets]
Every slab must come from the same base directory and the same code (args.code_sha256 of runc.py)."""
import os, sys, json, glob, shutil, argparse, math, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
S4 = os.path.normpath(os.path.join(HERE, '..', '..'))
_spec = importlib.util.spec_from_file_location('t6merge', os.path.join(S4, 't6', 'merge.py'))
M6 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(M6)

NOTE = ('T6 run on the (SE)-corrected base itself (review B1): kap^ = kap~ + (a1/ell) B(u) + (a2/ell^2) B\'(u), u = x/ell in [0.6, 0.8], '
        'a1 = -(theta~(ell) + 3pi/2), a2 the root of Im W^(ell) = 0 (Banach per tc, enclosure eps_a2), a1\', a2\' exact/implicit; '
        'p2/t6c/slabc.py (T6 engine, the two window elements cut at u = 0.6, 0.8; slab 11 cut in time at the knot crossing). '
        'Replaces the increments of p2/se_corr.py.')


def w2(v, d):
    return math.nextafter(math.nextafter(v, d), d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', default=os.path.join(HERE, 'out'))
    ap.add_argument('--first', type=int, default=2)
    ap.add_argument('--last', type=int, default=180)
    ap.add_argument('--json', default=os.path.join(HERE, 'residuals.json'))
    ap.add_argument('--jets', default=os.path.join(HERE, 'jets'))
    a = ap.parse_args()
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(a.inp, 'slab_*.json')))]
    recs = sorted((r for r in recs if a.first <= r['slab'] <= a.last), key=lambda r: r['slab'])
    have = [r['slab'] for r in recs]
    missing = [i for i in range(a.first, a.last + 1) if i not in have]
    pv = [r.get('args') for r in recs]
    if not all(pv): raise SystemExit('slab files without args: rerun runc.py')
    com = {k: pv[0][k] for k in ('poldir', 'code_sha256')}
    bad = [r['slab'] for r, p in zip(recs, pv) if any(p[k] != com[k] for k in com)]
    if bad: raise SystemExit(f'slabs {bad[:5]} come from another base dir or code version: rerun runc.py')
    for r in recs:
        M6._round_up(r)
        r['scalars'] = {k: ([w2(v[0], -math.inf), w2(v[1], math.inf)] if k.startswith('range_') and isinstance(v, list) and len(v) == 2 else v)
                        for k, v in r['scalars'].items()}
    summ = {}
    for nm, tg in M6.TARGETS.items():
        vals = [(r['norms'].get(nm, 0.0), r['slab']) for r in recs]
        if not vals: continue
        mx, at = max(vals)
        summ[nm] = dict(max=mx, at_slab=at, target=tg, margin=(tg / mx if mx > 0 else None), fails=[s for v, s in vals if v > tg])
    for nm in sorted(set(k for r in recs for k in r['norms']) - set(M6.TARGETS)):
        mx, at = max((r['norms'].get(nm, 0.0), r['slab']) for r in recs); summ[nm] = dict(max=mx, at_slab=at)
    for nm in sorted(set(k for r in recs for k in r['scalars'] if k.startswith('sup_'))):
        mx, at = max((r['scalars'][nm], r['slab']) for r in recs); summ['scalar:' + nm] = dict(max=mx, at_slab=at)
    for nm in sorted(set(k for r in recs for k in r['t7'])):
        mx, at = max((r['t7'][nm], r['slab']) for r in recs); summ['t7:' + nm] = dict(max=mx, at_slab=at)
    sed = {}
    for k in ('a2', 'a2p', 'h', 'hp', 'a1', 'a1p'):
        mx, at = max((max(s['sharp'][k] for s in r['se_direct']), r['slab']) for r in recs); sed[k] = dict(max=mx, at_slab=at)
    sed['exact_end_imW'] = max(max(s['exact_end']['imW'] for s in r['se_direct']) for r in recs)
    sed['exact_end_turn'] = max(max(s['exact_end']['turn'] for s in r['se_direct']) for r in recs)
    sed['eps_a2_max'] = max(max(s['eps_a2'] for s in r['se_direct']) for r in recs)
    sed['lip_max'] = max(max(s['lip'] for s in r['se_direct']) for r in recs)
    cpu = [r.get('cpu', r['diag']['cpu_total']) for r in recs]
    defs = dict(M6.DEFS, base_corrected=NOTE)
    args = dict(vars(a), argv=sys.argv[1:], inp=os.path.abspath(a.inp), runc=dict(com, env=pv[0]['env']),
                rounding='norms/norms_int/t7/sup_ +2 ulps (t6/merge.py _round_up); range_ pairs widened 2 ulps')
    out = dict(args=args, definitions=defs, targets=M6.TARGETS, slabs_range=[a.first, a.last], missing=missing,
               cpu_total_s=sum(cpu), cpu_per_slab_mean_s=(sum(cpu) / len(cpu) if cpu else None), se_direct=sed,
               summary=summ, per_slab=recs, _ranges_widened=True)
    tmp = a.json + '.tmp'
    with open(tmp, 'w') as f: json.dump(out, f)
    os.replace(tmp, a.json)
    os.makedirs(a.jets, exist_ok=True)
    idx = []
    for r in recs:
        src = os.path.join(a.inp, 'jets', f'jets_{r["slab"]:03d}.npz'); dst = os.path.join(a.jets, os.path.basename(src))
        if os.path.exists(src): shutil.copy2(src, dst)
        if os.path.exists(dst): idx.append(r['slab'])
    with open(os.path.join(a.jets, 'index.json'), 'w') as f:
        json.dump(dict(slabs=idx, base=NOTE,
                       format='as t6/jets (T6 element rows; the two cut elements are the union over their pieces); hulls of the '
                              'corrected base, rounded outward one ulp from the arb bounds (slab.py)'), f, indent=1)
    print(f'{len(recs)} slabs merged, missing {missing}; json -> {a.json}; jets -> {a.jets}')
    for nm in ('R1', 'RIs', 'RF0_opt4n', 'RF1_opt4n', 'RF2_opt4n', 'RF3_opt4n'):
        print(f"  {nm:10s} max {summ[nm]['max']:.3e} (slab {summ[nm]['at_slab']})")
    print('  se_direct', json.dumps(sed))


if __name__ == '__main__':
    main()
