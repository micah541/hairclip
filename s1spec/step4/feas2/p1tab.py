import json, glob, numpy as np
fs = sorted(glob.glob('logs/p1_*_S100.json')); R = {}
for f in fs: R.update(json.load(open(f)))
order = [('k1600', 3), ('k100', 3), ('k10', 3), ('t4.59', 3), ('t4.59', 4), ('t4.5', 4), ('t2.0', 4), ('t0.5', 4)]
print('| case | ncon | Rc | d2 / d1 / d0 (sup of tapered coeff diff) | N_close | nu_{N+1} | tail t | pen d2/d1/d0 | lam_min(A_bb) | margin@N_close | N for margin>=0.05 | margin@N=80 | limiting term at N_close |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for name, nc in order:
    for Rc in (20.0, 10.0):
        k = f'{name}_Rc{Rc}_nc{nc}_S100'
        if k not in R: continue
        v = R[k]; rows = v['rows']; info = v['info']; ok = [r for r in rows if r['margin_net'] > 0]
        if not ok: print(f'| {name} | {nc} | {Rc:g} | none | | | | | | | | | |'); continue
        r = ok[0]; r5 = [x for x in rows if x['margin_net'] >= 0.05]; r80 = rows[-1]
        lim = 'tail (nu_{N+1} vs sup-norm penalty)' if (r['t'] < 0.6) else 'cross-term/Galerkin block'
        if r['N'] == 3: lim = 'tail t = nu_4 - omega - small penalty'
        print(f"| {name} | {nc} | {Rc:g} | {info['d2']:.3g} / {info['d1']:.3g} / {info['d0']:.3g} | {r['N']} | {r['nuN1']:.2f} | {r['t']:.3f} | {r['pen_d2']:.2f} / {r['pen_d1']:.2f} / {r['pen_d0']:.2f} | {r['lamA']:.3f} | {r['margin_net']:.3f} | {r5[0]['N'] if r5 else '>80'} | {r80['margin_net']:.3f} | {lim} |")
