"""T8 float evaluation: amplification table, Ybar requirement, dominant terms, fallbacks.  usage: python3 report.py > report.out
Reads consts_<cfg>.json (run_all.py).  Writes results.json."""
import json, numpy as np, warnings
warnings.filterwarnings('ignore')
import amp, jcas
from fb import X0

def run(cfg, far_from=None, label=None, **kw):
    per, tau1 = amp.load(cfg)
    far = None
    if far_from is not None:
        perg, _ = amp.load(far_from)
        Rg = amp.integrate(perg, tau1, **{k: v for k, v in kw.items() if k != 'gam_key'})
        Omf = float(np.exp(-8.0 * (np.log(4.0) - 0.25)))
        far = (Rg, Omf)
    R = amp.integrate(per, tau1, far=far, **kw)
    h = amp.handoff(R, per, tau1)
    h['label'] = label or cfg; h['R'] = R; h['per'] = per; h['tau1'] = tau1
    return h

def line(h):
    A = h['A']; L = h['lim']
    return (f"| {h['label']} | " + ' | '.join(f"{A[j]:.2e}" for j in range(5)) +
            f" | {L['item3_E3']:.2e} | {L['item4_K4_viaF4']:.2e} | {L['item4_K4_goodtime']:.2e} | {min(L['item2_E0'], L['item2_E1'], L['item2_E2']):.2e} |")

if __name__ == '__main__':
    out = {}
    print('| design | A0 | A1 | A2 | A3 | A4 | Y: item 3 (E3) | Y: K4 via F4 | Y: K4 good time | Y: item 2 |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    rows = [('T5', dict(label='T5 design: Lemma S ramp, global')),
            ('T5loc', dict(far_from='T5', label='T5 design, localized (transport-following Omega)')),
            ('opt', dict(label='levels 3-4 ramp [14,22], delta_P=1.5')),
            ('nest', dict(label='levels 3-4 nested on D_1, D_2')),
            ('opt4n', dict(label='lev 0-3 [14,22], lev 4 nested on D_2')),
            ]
    H = {}
    for cfg, kw in rows:
        h = run(cfg, **kw); H[cfg] = h; print(line(h), flush=True)
    # fallbacks on the best design
    best = 'opt4n'
    for Cell in (100.0, 30.0, 10.0):
        h = run(best, D4pt=Cell, label=f'fallback 1 (time-diff.), D4 <= {Cell:g} Ybar'); print(line(h)); H[f'fb1_{Cell:g}'] = h
    Ws = lambda j, m, kat: 1.0 / (1.0 + (kat + X0)**2)
    h = run(best, Wscale=Ws, label='E_I weight (1+x^2)^-1 (J unchanged, indicative)'); print(line(h)); H['W1'] = h
    # parameter sensitivity on T5 design
    print()
    print('| T5 design, sensitivity | A3 | A4 | Y item 3 | Y K4 good time |')
    print('|---|---|---|---|---|')
    for lab, kw in [('C_x = 2', dict(Cx=2.0)), ('C_x = 10', dict(Cx=10.0)), ('kappa_D = 1 (T5 (Y) as written; not valid)', dict(P=dict(kD=1.0))),
                    ('kappa_D = 0.25', dict(P=dict(kD=0.25))), ('eps = 3e-5', dict(P=dict(eps=3e-5))), ('c3 = 0.2', dict(P=dict(cI=(0.425, 0.424, 0.526, 0.2)))),
                    ('theta0 = 0.05', dict(P=dict(th0=0.05))), ('mu_e = 0.08', dict(P=dict(muE=0.08)))]:
        h = run('T5', label=lab, **kw)
        print(f"| {lab} | {h['A'][3]:.2e} | {h['A'][4]:.2e} | {h['lim']['item3_E3']:.2e} | {h['lim']['item4_K4_goodtime']:.2e} |")
    # dominant terms at levels 3, 4 (T5 design and best)
    print()
    for cfg in ('T5', 'opt4n'):
        h = H[cfg]
        for j in (2, 3, 4):
            sh, J = amp.breakdown(h['R'], h['per'], h['tau1'], j)
            print(f"{cfg} level {j}: " + '; '.join(f"{k} {v:.2e}" for k, v in sh.items()))
        J = h['R']['diag'][3]['J'][-1]
        print(f"   J-cascade at tau1 (lam_min={h['R']['diag'][3]['lam_min'][-1]:.3f}): J0..J4 = " + ', '.join(f"{J[m]:.3e}" for m in range(5)))
        c = jcas.cascade(h['R']['diag'][3]['lam_min'][-1])
        print(f"   J4 shares: D4 part {c['JV4_shares']['D4']:.3e} of {c['JV4_shares']['total']:.3e}; JB={c['JB']:.3g}, JD4={c['JD4']:.3g}, JW1={c['JW1']:.3g}, JW3={c['JW3']:.3g}")
    # time profile of the constants (T5 design): sup over slabs of k_jj, sup of c_ji, sigma
    per = H['T5']['per']
    print()
    for j in range(5):
        kk = [p['lev'][j]['k'][j] for p in per]; lam = [p['lev'][j]['lam'] for p in per]
        cmax = max((max(p['lev'][j]['c'].values()) if p['lev'][j]['c'] else 0.0) for p in per)
        cmax_late = max((max(p['lev'][j]['c'].values()) if p['lev'][j]['c'] else 0.0) for p in per if p['k'] > 850)
        g = [p['lev'][j]['gam'] for p in per]; gc = [p['lev'][j]['gamc'] for p in per]
        print(f"level {j}: k_jj in [{min(kk):.3e}, {max(kk):.3e}] (at y~{per[-1]['lev'][j]['kat']:.1f}), lam in [{min(lam):.4f}, {max(lam):.4f}], "
              f"max c_ji {cmax:.2e} (k>850: {cmax_late:.2e}), gamma_j max {max(g):.2e} (cone y<600: {max(gc):.2e})")
    json.dump({k: dict(A=v['A'], lim=v['lim'], label=v['label'], D2=v['D2_goodtime']) for k, v in H.items()}, open('results.json', 'w'), indent=1)
