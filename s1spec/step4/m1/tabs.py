import json, numpy as np, sys
tags = sys.argv[1].split(',')
print("| tag | k_max | N | lambda | lambda*k_max | x_c (lab, real) | alpha | alpha*k_max | alpha/5e-4 | A=\\|w~\\|_rho | A/\\|alpha\\|^3 | sqrt(mu E_I) | sqrt(Lambda) | sup \\|w'\\| on [16.25,57] | S1-chart: first y with \\|f\\|>0.5 | odd-part norm (rho, x<=28) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in tags:
    r = json.load(open(f'res_{t}.json')); i = r['inner']; k = r['kmax']; a = r['alpha']
    print(f"| {t} | {k:.1f} | {r["N"]//16} | {r['lam']:.4e} | {r['lam']*k:.5f} | {r['xc'][0]:.11f} | {a:+.4e} | {a*k:+.4f} | {abs(a)/5e-4:.2f} | {i['A']:.3e} | {i['A']/abs(a)**3:.0f} | {i['sqrtmuEI']:.2e} | {i['sqrtLam']:.2e} | {i['supwp']:.2e} | {r['chart_ext'][0][0]:.0f} | {r['odd_norm']:.1e} |")
