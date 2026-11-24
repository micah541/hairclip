"""Windows and shells (H2, H3) uniformly on slab 180: 32 t_c balls of radius 1/32, ball m = the half-ball between t_k = -1 + (2k+1)/16
(k = m // 2) and t_k -+ 1/16.  Its parameter hull is the path p(t) in p(t_k) + (t - t_k) pdot (p(t_k): the Krawczyk enclosure of h1.py,
pdot: the verified enclosure of h1dt.py on that ball); Theta from the h1 runs.  Writes out/h1slab/ift_m.json (T9IFT for win.py) and the
win.py environment lines."""
import json, math
O = '/home/claude/s1spec/step4/p2/out'
for m in range(32):
    k, sg = m // 2, (-1 if m % 2 == 0 else 1); tk = -1 + (2 * k + 1) / 16; c = -1 + (2 * m + 1) / 32
    h = json.load(open(f'{O}/t9_h1_p{tk}_coh_s180.json')); d = json.load(open(f'{O}/h1slab/dt_{tk + sg * 0.03125}_0.03125.json'))
    assert abs(d['tc'][0] + d['tc'][1] - 2 * c) < 1e-9
    pd = d['pdot']; hull = {}
    for key, i in (('xc', 0), ('alpha', 2)):
        lo, hi = h[key]; ext = [sg * pd[i][0] / 16, sg * pd[i][1] / 16]
        hull[key] = [math.nextafter(lo + min(0.0, *ext), -math.inf), math.nextafter(hi + max(0.0, *ext), math.inf)]
    ext = [sg * pd[1][0] / 16, sg * pd[1][1] / 16]
    hull['lam'] = [math.nextafter(h['lam'][0] * math.exp(min(0.0, *ext)) * (1 - 1e-15), -math.inf), math.nextafter(h['lam'][1] * math.exp(max(0.0, *ext)) * (1 + 1e-15), math.inf)]
    json.dump(dict(alpha=hull['alpha'], xc=hull['xc'], lam=hull['lam'], Theta=h['Theta'], tc_ball=[c - 1 / 32, c + 1 / 32],
                   args=dict(ball=m, balls=32, radius=1 / 32, tk=tk, h1=f't9_h1_p{tk}_coh_s180.json', dt=f'h1slab/dt_{tk + sg * 0.03125}_0.03125.json')),
              open(f'{O}/h1slab/ift_{m}.json', 'w'), indent=1)
    print(f'T9TC={c},0.03125 T9IFT={O}/h1slab/ift_{m}.json WINOUT={O}/h1slab/win_{m}.json')
