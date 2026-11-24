"""time-step floor: compare conv_{tag}_{k} snapshot fits with the feas base (cdt=1e-2) at k = 100, 400, 1600 modulo (T-shift, x_c, alpha): same transport correction as f2an2.
usage: python3 floorcheck.py tag"""
import sys, numpy as np, json; sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
import chart
from chart import hn, hnW
tag = sys.argv[1]; CAF = '/home/claude/s1spec/step4/feas/logs/f2_cache'; out = {}
for lev in (100, 400, 1600):
    try:
        z = np.load(f'/home/claude/s1spec/step4/feas2/snaps/conv_{tag}_{lev}.npz')
    except Exception as e: print('missing', lev); continue
    b = dict(np.load(f'{CAF}/base_{lev}.npz', allow_pickle=True)); bx = dict(np.load(f'{CAF}/baseX_{lev}.npz', allow_pickle=True))
    p = chart.analyse(f'/home/claude/s1spec/step4/feas2/snaps/conv_{tag}_{lev}.npz')
    n = min(int(b['n']), int(p['n']), int(bx['n'])); dl = float(np.log(p['lam'].real / b['lam'].real)); dlx = float(np.log(bx['lam'].real / b['lam'].real))
    sw = (bx['w'][:n] - b['w'][:n]) / dlx; sa = float(bx['alpha'] - b['alpha']) / dlx; dw = p['w'][:n] - b['w'][:n] - sw * dl
    zb = np.load(f'/home/claude/s1spec/step4/feas/snaps/base_{lev}.npz')
    out[lev] = dict(t_conv=float(z['t']), t_base=float(zb['t']), dT=float(z['t'] - zb['t']), dw_rho=hn(dw, n), dw_W=hnW(dw, n), dalpha=float(p['alpha'] - b['alpha']) - sa * dl, alpha_b=float(b['alpha']), dlnlam=dl)
    print(lev, out[lev], flush=True)
json.dump(out, open(f'/home/claude/s1spec/step4/feas2/logs/floor_{tag}.json', 'w'), indent=1)
