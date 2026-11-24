import json,glob,sys
ks=[1600,2000,3000,4000,5000,6000,7000,8000,10000]
R={k:json.load(open(f'res_main_{k}.json')) for k in ks}
print("k alpha a*k sqrtLam A muEI^.5 B | chart_ext(first|f|>.5, q-range end) | lam")
for k in ks:
    r=R[k]; i=r['inner']
    print(f"{r['kmax']:7.0f} {r['alpha']:+.3e} {r['alpha_k']:.3f} {i['sqrtLam']:.2e} {i['A']:.2e} {i['sqrtmuEI']:.2e} {i['B']:.1e} | {r['chart_ext'][0][0]:.0f} {r['chart_ext'][0][1]:.0f} | lam {r['lam']:.3e} odd {r['odd_norm']:.1e} fit {r['fit_r2']:.2e}")
print("\nshell [14,60] E0 E1 E2 E3 E4")
for k in ks: print(k, ["%.2e"%e for e in R[k]['shell']['14-60']])
for w in ('50-75','75-112','112-169','169-253'):
    print("\nwindow",w)
    for k in ks: print(k, ["%.2e"%e for e in R[k]['shell'][w]])
print("\nslide max over [y,1.5y], y in [50,253] (E0..E4 value@y)")
for k in ks: print(k, {j:("%.2e"%v[0], round(v[1]) if v[1] else None) for j,v in R[k]['slide'].items()})
print("\ncone sup 200:", )
for k in ks: print(k, R[k]['cone_sup']['200'], R[k]['cone_sup']['200-226'] )
print("cone 200-400"); 
for k in ks: print(k, R[k]['cone_sup']['200-400'])
print("far"); 
for k in ks: print(k, R[k]['far'], R[k]['pt'])
