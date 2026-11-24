import json, numpy as np
TT = 4.598863926611; K0 = 1.5398; b = 0.02 / 0.07
R = json.load(open('logs/p2_all.json'))
def sfix(t, k): return -np.log(TT - t) if t < 4.5905 else 4 * np.log(k / K0)
rows = []
for x in R:
    s_a, s_b = sfix(x['t_a'], x['k_a']), sfix(x['t_b'], x['k_b']); ds = s_b - s_a
    macro = x['t_b'] <= 4.5905 and x['k_a'] < 5.1
    th = x['theta_nc4'] if macro else x['theta_nc3']; co = x['coer_a_nc4'] if macro else x['coer_a_nc3']
    # local clock tau = int dt / lam^4, lam = K0/kappa_tip ; use endpoints average of kappa^4
    dt = x['t_b'] - x['t_a']; dtau = dt * 0.5 * (x['k_a']**4 + x['k_b']**4) / K0**4 if macro else ds * 1.0
    rows.append(dict(t_a=x['t_a'], t_b=x['t_b'], k_a=x['k_a'], k_b=x['k_b'], s_a=s_a, s_b=s_b, ds=ds, theta=th, ncon=4 if macro else 3, coer=co, dtau=dtau,
                     rate=th / ds, Ds_allowed=b * ds / th if th > 0 else np.inf, nrow=th / b))
def phase(r):
    if r['t_b'] <= 0.0101: return 'layer t<0.01'
    if r['t_b'] <= 4.5505: return 'macro 0.01<t<4.55'
    if r['k_a'] < 7: return 'approach k 3.3-7'
    return 'rapid k 7-1810'
out = {}
for r in rows: out.setdefault(phase(r), []).append(r)
print('| t_a -> t_b | k_a -> k_b | s_a -> s_b | ds | ncon | coercivity at a | theta | theta/ds | Delta s allowed (m=0.02) | dtau (local clock) |'); print('|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    print(f"| {r['t_a']:.4g} -> {r['t_b']:.4g} | {r['k_a']:.3g} -> {r['k_b']:.3g} | {r['s_a']:.3f} -> {r['s_b']:.3f} | {r['ds']:.3g} | {r['ncon']} | {r['coer']:.3f} | {r['theta']:.3g} | {r['rate']:.3g} | {r['Ds_allowed']:.3g} | {r['dtau']:.3g} |")
print()
tot = 0
for ph, rs in out.items():
    S = sum(r['theta'] for r in rs); n = S / b
    # greedy count with additivity
    cnt = 1; acc = 0.0
    for r in rs:
        if r['theta'] > b: cnt += int(np.ceil(r['theta'] / b)) - 1; acc = 0
        elif acc + r['theta'] > b: cnt += 1; acc = r['theta']
        else: acc += r['theta']
    print(f"{ph}: sum theta = {S:.3g}, sum ds = {sum(r['ds'] for r in rs):.3g}, certificates >= sum theta/b = {n:.1f}, greedy ~ {cnt}")
# extrapolate rapid beyond 1810 to 3000: theta ~ ratio per step
rap = out['rapid k 7-1810']; th = np.array([r['theta'] for r in rap]); k0 = rap[-1]['k_b']
print('last rows theta:', th[-4:], 'ratio', th[-1] / th[-2]); ex = th[-1] * (np.log(3000 / k0) / (4 * 0.5 * np.log(2) * 1)) * 0.7**0.5
print('extrapolated theta 1810->3000 ~', ex)
# theta' vs k  (rapid): fit theta' = C k^-p
k = np.array([np.sqrt(r['k_a'] * r['k_b']) for r in rap]); rate = np.array([r['rate'] for r in rap]); m = k > 30
p = np.polyfit(np.log(k[m]), np.log(rate[m]), 1); print('rapid rate fit  theta/ds = %.3g * k^%.3f (k>30)' % (np.exp(p[1]), p[0]))
json.dump(rows, open('logs/p2_rows.json', 'w'), indent=1, default=float)
