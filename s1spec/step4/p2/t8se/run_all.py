"""Run the Lemma F constants over all 188 slabs of the T1 base (theta = 0, 1/2, 1) for several designs -> consts_<tag>.json.
usage: python3 run_all.py  (1 process, ~10 min)"""
import json, time, numpy as np, fb, lfnum as LF
B = fb.base()
tau1 = B.tau_of_k(2.3e4)[0]
LSr = (14.0, 28.0, 0.0)
def Yh(tau, Y1=400.0, ch=0.20): return Y1 * np.exp(-ch * (tau1 - tau)), ch * Y1 * np.exp(-ch * (tau1 - tau))
CFG = {
  'T5':   lambda s: dict(),
  'T5loc': lambda s: dict(om=(8.0, *Yh(s['tau'])), ynear=4 * Yh(s['tau'])[0]),
  'opt':  lambda s: dict(ramps=(LSr, LSr, LSr, (14.0, 22.0, 0.0), (14.0, 22.0, 0.0)), dP=(0.25, 0.25, 0.25, 1.5, 1.5)),
  'nest': lambda s: dict(ramps=((14.0, 20.0, 0.0),) * 3 + ((20.0, 28.0, 0.0),) * 2, nested={3: 1, 4: 2}),
  'opt4n': lambda s: dict(ramps=((14.0, 22.0, 0.0),) * 4 + ((22.0, 30.0, 0.0),), nested={4: 2}, dP=(0.25, 0.25, 0.25, 1.5, 1.5)),
}
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
res = {k: [] for k in CFG}; t0 = time.time()
for i in range(len(B.slabs)):
    for th in (0.0, 0.5, 1.0):
        s = fb.sample(i, th)
        for name, f in CFG.items():
            o = LF.constants(s, **f(s)); o.update(slab=i, theta=th, nV=s['nV'], closure=s['closure'])
            res[name].append(clean(o))
    if i % 10 == 0 or i == len(B.slabs) - 1:
        print(f'slab {i} k={s["k"]:.4g} {time.time()-t0:.0f}s', flush=True)
        for name in CFG: json.dump(dict(tau1=tau1, rows=res[name]), open(f'consts_{name}.json', 'w'))
print('DONE', time.time() - t0)
