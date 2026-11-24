"""F2b analysis.  usage: python3 f2an2.py [tj ...]   (default all with snapshots).  Cache: logs/f2_cache."""
import sys, os, glob, json, time, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
import chart
from chart import hn, hnW
D = '/home/claude/s1spec/step4/feas2'; SN = f'{D}/snaps'; CA = f'{D}/logs/f2_cache'
LEV = [100, 400, 1600]; LEVX = {100: 100.3, 400: 401.2, 1600: 1604.8}
TJ = ['4.5', '3.0', '1.5', '0.3', '0']

def get(tag, lev):
    f = f'{CA}/{tag}_{lev:g}.npz'
    if os.path.exists(f):
        z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files} if int(z['ok']) else None
    path = f'{SN}/{tag}_{lev:g}.npz'
    if not os.path.exists(path): return 'missing'
    try:
        t0 = time.time(); d = chart.analyse(path); d['ok'] = 1; np.savez(f, **{k: np.asarray(v) for k, v in d.items()}); print(f'  analysed {tag} {lev:g} [{time.time()-t0:.0f}s]', flush=True)
    except Exception as e:
        np.savez(f, ok=0, err=str(e)); print(f'  analysis FAILED {tag} {lev:g}: {e!r}', flush=True); return None
    z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files}

def compare(ctl, tag, lev):
    b, bx, p = get(ctl, lev), get(ctl, LEVX[lev]), get(tag, lev)
    if any(x is None or isinstance(x, str) for x in (b, bx, p)): return None
    n = min(int(b['n']), int(p['n']), int(bx['n']))
    dl = float(np.log(p['lam'].real / b['lam'].real)); dlx = float(np.log(bx['lam'].real / b['lam'].real))
    sw = (bx['w'][:n] - b['w'][:n]) / dlx; sa = float(bx['alpha'] - b['alpha']) / dlx
    dw = p['w'][:n] - b['w'][:n] - sw * dl
    return dict(lev=lev, n=n, xe=float(chart.xs[n - 1]), dw_rho=hn(dw, n), dw_W=hnW(dw, n), dalpha=float(p['alpha'] - b['alpha']) - sa * dl, dlnlam=dl,
                dw_rho_raw=hn(p['w'][:n] - b['w'][:n], n), alpha_b=float(b['alpha']), kmax_b=float(b['kmax']), kmax_p=float(p['kmax']), dw=dw)

if __name__ == '__main__':
    tjs = sys.argv[1:] or TJ; out = {}
    for tj in tjs:
        ctl = f't{tj}_ctl'
        for lev in LEV: get(ctl, lev); get(ctl, LEVX[lev])
        for sh in ('peak', 'fixed', 'far'):
            for r in (1, 2):
                tag = f't{tj}_{sh}_{r}'
                if not os.path.exists(f'{D}/logs/pert_{tag}.json'): continue
                info = json.load(open(f'{D}/logs/pert_{tag}.json')); rows = []
                for lev in LEV:
                    if not os.path.exists(f'{SN}/{tag}_{lev:g}.npz'): continue
                    c = compare(ctl, tag, lev)
                    if c is not None:
                        rows.append({k: v for k, v in c.items() if k != 'dw'}); np.savez(f'{CA}/dw_{tag}_{lev}.npz', dw=c['dw'])
                out[tag] = dict(tj=tj, shape=sh, delta=info['delta'], info=info, rows=rows)
                for c in rows: print(f"{tag} k={c['lev']}: dw_rho/delta={c['dw_rho']/info['delta']:.3e} dw_W/delta={c['dw_W']/info['delta']:.3e} dalpha/delta={c['dalpha']/info['delta']:+.3e} xe={c['xe']:.1f}", flush=True)
    json.dump(out, open(f'{D}/logs/f2b_results.json', 'w'), indent=1)
