import json, numpy as np, sys
tags = sys.argv[1].split(',')
R = [json.load(open(f'res_{t}.json')) for t in tags]
k = np.array([r['kmax'] for r in R])
c = 40.0
thr = {0: 2e-7, 1: 3e-6, 2: 3e-6, 3: 3e-6, 4: 3e-6}
WINS = {'50-75': (50, 75), '75-112': (75, 112), '112-169': (112, 169), '169-253': (169, 253)}
print("| window | level | exponent p (E ~ k^-p, fit over k>=5000) | E at k=%.0f | threshold | fitted k_cross | closed-form k_cross (c3=40) |" % k[-1])
print("|---|---|---|---|---|---|---|")
for w in ('112-169', '169-253'):
    for j in (0, 1, 2, 3, 4):
        v = np.array([r['shell'][w][j] for r in R]); sel = k >= 4500
        pf = np.polyfit(np.log(k[sel]), np.log(v[sel]), 1); kc = np.exp((np.log(thr[j]) - pf[1]) / pf[0])
        ya, yb = WINS[w]; kcf = (28 * 9**j * c**2 * (yb**5 - ya**5) / (5 * thr[j]))**(1 / 6)
        print(f"| [{ya}, {yb}] | E{j} | {-pf[0]:.2f} | {v[-1]:.2e} | {thr[j]:g} | {kc:.0f} | {kcf:.0f} |")
for j in range(5):
    v = np.array([r['slide'][str(j)][0] for r in R]); sel = k >= 4500
    pf = np.polyfit(np.log(k[sel]), np.log(v[sel]), 1); kc = np.exp((np.log(thr[j]) - pf[1]) / pf[0])
    kcf = (28 * 9**j * c**2 * (380**5 - 253**5) / (5 * thr[j]))**(1 / 6)
    print(f"| sliding max (last window [253, 380]) | E{j} | {-pf[0]:.2f} | {v[-1]:.2e} | {thr[j]:g} | {kc:.0f} | {kcf:.0f} |")
