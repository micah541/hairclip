import json, numpy as np, sys
tags = sys.argv[1].split(',')
R = [json.load(open(f'res_{t}.json')) for t in tags]
k = np.array([r['kmax'] for r in R])
def fit(vals, kk):
    p = np.polyfit(np.log(kk), np.log(vals), 1); return -p[0], p
print("exponents p in E ~ k^-p  (fit over k in", k.round(), ")")
items = {}
for w in ('50-75', '75-112', '112-169', '169-253'):
    for j in (0, 1, 2, 3, 4): items[f'E{j}[{w}]'] = np.array([r['shell'][w][j] for r in R])
for j in range(5): items[f'E{j}slide'] = np.array([r['slide'][str(j)][0] for r in R])
thr = {0: 2e-7, 1: 3e-6, 2: 3e-6, 3: 3e-6, 4: 3e-6}
for name, v in items.items():
    j = int(name[1]); sel = k >= 4500
    p, c = fit(v[sel], k[sel]); kc = np.exp((np.log(thr[j]) - c[1]) / c[0])
    print(f"{name:16s} p={p:5.2f}  last={v[-1]:.2e}  k_cross(thr {thr[j]:g})={kc:8.0f}")
