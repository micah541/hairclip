"""F2 analysis: fit (xc, lam, alpha, beta2) on the even part for base and perturbed snapshots, Delta w = w_pert - w_base (common chart), amplification.
usage: python3 f2an.py [tags...]   (default: all p*_{1,2} tags with snapshots)"""
import sys, os, glob, json, time, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
import chart
from chart import hn, hnW
SN = '/home/claude/s1spec/step4/feas/snaps'; CA = '/home/claude/s1spec/step4/feas/logs/f2_cache'
LEV = [10, 30, 100, 400, 1600]

def get(tag, lev):
    f = f'{CA}/{tag}_{lev}.npz'
    if os.path.exists(f):
        z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files} if int(z['ok']) else None
    path = f'{SN}/{tag}_{lev}.npz'
    if not os.path.exists(path): return 'missing'
    try:
        d = chart.analyse(path); d['ok'] = 1
        np.savez(f, **{k: np.asarray(v) for k, v in d.items()})
    except Exception as e:
        np.savez(f, ok=0, err=str(e)); print(f'  analysis failed {tag} {lev}: {e!r}', flush=True); return None
    z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files}

def nrm(d, which, n):
    return hn(d, n) if which == 'rho' else hnW(d, n)

def compare(tag, lev, delta, info):
    b = get('base', lev); p = get(tag, lev)
    if p is None or isinstance(p, str) or b is None or isinstance(b, str): return None
    bx = get('baseX', lev)
    n = min(int(b['n']), int(p['n']), int(bx['n']))
    dl = float(np.log(p['lam'].real / b['lam'].real)); dlx = float(np.log(bx['lam'].real / b['lam'].real))
    sw_ = (bx['w'][:n] - b['w'][:n]) / dlx; sa_ = float(bx['alpha'] - b['alpha']) / dlx
    dw = p['w'][:n] - b['w'][:n] - sw_ * dl            # remove the k_max (T-shift) mismatch: base transported to the perturbed run's lambda
    dwraw = p['w'][:n] - b['w'][:n]
    return dict(lev=lev, n=n, xe=float(xs_[n - 1]), dw_rho=hn(dw, n), dw_W=hnW(dw, n), dalpha=float(p['alpha'] - b['alpha']) - sa_ * dl, dalpha_raw=float(p['alpha'] - b['alpha']), dlnlam=dl, dw_rho_raw=hn(dwraw, n), dlam_rel=float(p['lam'] / b['lam'] - 1),
                dxc_lam=float((p['xc'] - b['xc']) / b['lam'].real), alpha_b=float(b['alpha']), kmax_p=float(p['kmax']), kmax_b=float(b['kmax']),
                dodd=hn(p['odd'][:n] - b['odd'][:n], n), dw=dw)
xs_ = chart.xs

if __name__ == '__main__':
    tags = sys.argv[1:] or sorted({os.path.basename(f).rsplit('_', 1)[0] for f in glob.glob(f'{SN}/p*_*_[12]_*.npz')})
    out = {}
    for lev in LEV: get('base', lev); get('baseX', lev)
    for tag in tags:
        kj = int(tag[1:].split('_')[0]); d = 1e-6 * int(tag[-1]); info = json.load(open(f'/home/claude/s1spec/step4/feas/logs/pert_{tag}.json'))
        rows = []
        for lev in LEV:
            if lev < kj: continue
            r = compare(tag, lev, d, info)
            if r is not None: rows.append(r)
        out[tag] = dict(kj=kj, delta=d, info=info, rows=[{k: v for k, v in r.items() if k != 'dw'} for r in rows])
        np.savez(f'{CA}/dw_{tag}.npz', **{f'dw{r["lev"]}': r['dw'] for r in rows})
        for r in rows:
            print(f"{tag} obs k={r['lev']}: |dw|_rho/d={r['dw_rho']/d:.3e} |dw|_W/d={r['dw_W']/d:.3e} dalpha/d={r['dalpha']/d:+.3e} xe={r['xe']:.1f}", flush=True)
    json.dump(out, open('/home/claude/s1spec/step4/feas/logs/f2_results_partial.json', 'w'), indent=1)
