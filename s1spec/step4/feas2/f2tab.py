import json, numpy as np
R = json.load(open('logs/f2b_results.json'))
TJ = ['0', '0.3', '1.5', '3.0', '4.5']; SH = ['peak', 'fixed', 'far']
print('| t_j | shape (w, amp, ell_j) | A_rel rho k=100 | 400 | 1600 | A_rel W k=100 | 400 | 1600 | A_phys rho k=1600 [1/len] | dalpha/delta k=100/400/1600 | linearity (d2/d1)/2 at 100/400/1600 |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for tj in TJ:
    for sh in SH:
        a = R.get(f't{tj}_{sh}_1')
        if not a: continue
        i = a['info']; norm = i['amp'] / i['ell']
        rows = {r['lev']: r for r in a['rows']}; b = R.get(f't{tj}_{sh}_2'); rb = {r['lev']: r for r in b['rows']} if b else {}
        d = a['delta']
        f = lambda key, lev, scale=1.0: (f"{rows[lev][key] / d / scale:.2e}" if lev in rows else '-')
        lin = [(f"{rb[l]['dw_rho'] / rows[l]['dw_rho'] / 2:.2f}" if (l in rb and l in rows) else '-') for l in (100, 400, 1600)]
        print(f"| {tj} | {sh} (w={i['w']:.3g}, amp={i['amp']:.3g}, ell={i['ell']:.3g}, sc={i['sc']:.2f}) | {f('dw_rho',100,norm)} | {f('dw_rho',400,norm)} | {f('dw_rho',1600,norm)} | {f('dw_W',100,norm)} | {f('dw_W',400,norm)} | {f('dw_W',1600,norm)} | {f('dw_rho',1600,i['amp'])} | {f('dalpha',100)} / {f('dalpha',400)} / {f('dalpha',1600)} | {' / '.join(lin)} |")
