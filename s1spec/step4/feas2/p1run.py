import sys, json, time; sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from p1 import *
import chart
cases = {'k1600': ('/home/claude/s1spec/step4/feas/snaps/base_1600.npz', 'fit'), 'k100': ('/home/claude/s1spec/step4/feas/snaps/base_100.npz', 'fit'),
         'k10': ('/home/claude/s1spec/step4/feas/snaps/base_10.npz', 'tip'), 't4.59': ('snaps/state_t4.59.npz', 'tip'), 't4.5': ('snaps/state_t4.5.npz', 'tip'),
         't2.0': ('snaps/state_t2.0.npz', 'tip'), 't0.5': ('snaps/state_t0.5.npz', 'tip')}
sel = sys.argv[1].split(','); Rcs = [None, 20.0, 15.0, 10.0]
h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
ncon = int(sys.argv[3]) if len(sys.argv) > 3 else 3
SIG = float(sys.argv[4]) if len(sys.argv) > 4 else 35.0
Rcs = [None, 15.0] if SIG <= 40 else [20.0, 15.0, 10.0]
st = Setup(Sigma=SIG, h=h, nmax=82); const = st.ref_constants(); res = {}
for name in sel:
    path, fr = cases[name]; t0 = time.time(); kw = {}
    if fr == 'fit':
        d = chart.analyse(path); kw = dict(xc=d['xc'], lam=d['lam'])
    g = numeric_geometry(path, Smax=40, dsig=0.01 if name in ('k100', 'k1600') else 0.005, **kw); sp = as_splines(g)
    print(f"#### {name} frame={fr} lam={g['lam']:.4e} (geometry {time.time()-t0:.0f}s)", flush=True)
    for Rc in Rcs:
        info, rows = certificate(st, sp, True, Rc, const=const, ncon=ncon); show(f'{name}', info, rows); res[f'{name}_Rc{Rc}_nc{ncon}_S{SIG:g}'] = dict(info=info, rows=rows)
    json.dump(res, open(f'logs/p1_{sys.argv[1].replace(",", "_")}_h{h}_nc{ncon}_S{SIG:g}.json', 'w'), indent=1, default=float)
