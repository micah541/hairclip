"""[kit file] The certified inequalities of STEP4_RESULT (final tail), as (name, value, op, bound) read from the outputs under a step4 tree.
extract(S4) is run on the kit's tree by check.py and, at build time, on the proof tree (reference/ours.json: the 'ours' column).
A value that cannot be read is None (check.py reports it as a FAIL: output missing).  Bounds are the thresholds of the proof, not tolerances."""
import os, re, json, math
up = lambda x: math.nextafter(x, math.inf)


def extract(S4):
    P2 = os.path.join(S4, 'p2'); O = os.path.join(P2, 'out'); T8C = os.path.join(P2, 't8se', 'constants')
    def J(*p):
        try:
            with open(os.path.join(*p), encoding='utf-8') as f: return json.load(f)
        except (OSError, ValueError): return None
    def T(*p):
        try:
            with open(os.path.join(*p), encoding='utf-8', errors='replace') as f: return f.read()
        except OSError: return None
    R = []; I = []
    def c(name, fn, op, bound):
        try: v = fn()
        except Exception: v = None
        R.append((name, v, op, bound))
    def info(name, fn):
        try: v = fn()
        except Exception: v = None
        I.append((name, v))
    sw = J(O, 'switch_s2.json'); yc = J(O, 'ycomp_final.json'); amp = J(T8C, 'amp_A1_w180.json'); kch = J(O, 't9_cone_kch.json'); far = J(O, 't9_far.json')
    cl = J(O, 't9_h1_clam_s180.json'); hs = J(O, 'h1slab', 'summary.json'); ws = J(O, 'h1slab', 'win_summary.json'); h4 = J(O, 'h4loc_A1_w180.json')
    # --- switch, (Y)
    c('switch (A): Y(tau_a) <= typed Y0 = 6.21e-13 (Ybar = 1e-12)', lambda: sw['modulated']['Y_raw'], '<=', 6.21e-13)
    for j, f0 in enumerate((4e-23, 3e-22, 5e-19, 3.4e-15)): c(f'switch: F_{j}(tau_a) <= typed F0[{j}] used by camp', lambda j=j: sw['modulated']['F_raw'][j], '<=', f0)
    c('(Y): comparison level on the direct residuals <= 3.22e-14', lambda: yc['level_max'], '<=', 3.22e-14)
    c('(Y): Ybar = 1e-12 majorant Yc >= typed Y0 = 6.21e-13', lambda: yc['Yc'], '>=', 6.21e-13)
    # --- (F), T8
    lim = lambda k: (lambda: amp['handoff']['lim'][k])
    for k, nm in (('item3_E3', 'item 3 E3 (window, binding)'), ('item3_E2', 'item 3 E2'), ('item3_E1', 'item 3 E1'), ('item3_E0', 'item 3 E0'), ('item4_K4_goodtime', 'item 4 (K4 at the good time)'),
                  ('item2_E0', 'item 2 E0 (shells, own bound)'), ('item2_E1', 'item 2 E1'), ('item2_E2', 'item 2 E2')):
        c(f'(F) {nm}: Ybar limit >= 1e-12', lim(k), '>=', 1e-12)
    c('(F) Ymax (min over items) >= the 1.7e-12 box level of N3', lambda: amp['handoff']['Ymax'], '>=', 1.7e-12)
    c('(F) D2(tau*) at Ybar = 1e-12 <= 9.40e-12 (A1 window = slab 180)', lambda: amp['handoff']['at_Ybar']['D2_goodtime'], '<=', 9.40e-12 * (1 + 1e-3))
    bt = T(O, 'n3boxes_final_B.log') or ''
    m = re.search(r"at the n3 Ybar = [0-9.e+-]+: \|e'\| (\S+) <= (\S+), x\|e''\| (\S+) <= (\S+), x\^2\|e'''\| (\S+) <= (\S+), sqrtF0 ratio (\S+) <= (\S+)", bt)
    for i, nm in ((0, "far core |e'| at Y = 1.7e-12 <= b0 = 1.5e-3"), (2, "far core x|e''| <= b1 = 3e-4"), (4, "far core x^2|e'''| <= b2 = 1e-2"), (6, 'sqrt(F0) ratio <= fac = 1.7')):
        c('boxes: ' + nm, lambda i=i: float(m.group(i + 1).rstrip(',:')), '<=', float(m.group(i + 2).rstrip(',:')) if m else None)
    m2 = re.search(r"T7 Lemma 4b gives (\S+) at Ybar", bt)
    c('boxes: tip |e|, |e\'| (T7 Lemma 4b) at Y = 1.7e-12 <= 1e-5', lambda: float(m2.group(1)), '<=', 1e-5)
    # --- H1
    c('(H1) sqrt(Lambda) over slab 180 incl. C_Lambda Ybar <= 5e-8', lambda: cl['sqrtLam_total_Wp'], '<=', 5e-8)
    c('(H1) |alpha| over slab 180 + C_alpha Ybar <= 1.5e-4', lambda: up(hs['alpha_abs_max'] + cl['C_alpha'] * 1e-12), '<=', 1.5e-4)
    c('(H1) Ybar limit of H1 >= 1e-12', lambda: cl['Ybar_limit_H1'], '>=', 1e-12)
    def thw():
        pts = [J(O, f't9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json') for k in range(16)]
        return (max(p['Theta'][1] for p in pts) - min(p['Theta'][0] for p in pts)) / 2
    c('(H1) Theta half-width over the 16 points <= 1e-10', thw, '<=', 1e-10)
    # --- H2, H3
    for j, thr in enumerate((2e-7, 1e-5, 1e-5, 1e-5)): c(f'(H3) window energy E{j} (base, slab 180) <= threshold', lambda j=j: ws['windows_max'][j], '<=', thr)
    for j, thr in enumerate((2e-7, 1e-4, 0.2)): c(f'(H2) shell energy E{j} on y in [14, 60] (base) <= threshold', lambda j=j: ws['shells'][j], '<=', thr)
    c('(H3) K4 on y in [50, 350] (r_c = 200) <= 0.06', lambda: kch['cone'][0]['K4'], '<=', 0.06)
    c('(H4) graph over the ray: |theta - Theta| <= 0.079', lambda: kch['phi_max'], '<=', 0.079)
    # --- H5
    c('KCH (chart factor, computed) <= 1.188', lambda: kch['KCH'], '<=', 1.188)
    for k in ('K2', 'K3', 'K4'): c(f'(H5) cone {k} <= 0.12 (r_c in [200, 380], base + error at Ybar)', lambda k=k: kch['cone_max'][k], '<=', 0.12)
    c('(H5) cone N3 <= 0.14', lambda: kch['cone_max']['N3'], '<=', 0.14)
    for nm, key in (("|h''|", 'h2'), ("|h'''|", 'h3'), ("|h''''|", 'h4'), ('slope', 'slope'), ("||h'''||_L2", 'L2')):
        c(f'(H5) far cylinders x >= 150: {nm}', lambda key=key: far[key], '<=', far['thr'][key] if far else None)
    wk = lambda d: d['windows']['constants/amp_A1_w180.json']['worst']['ratio']
    c('(H5) Lemma O H4: worst cylinder of the final window, int (h0^(5))^2 / M^2 <= 1', lambda: wk(h4), '<=', 1.0)
    # --- v3.2 candidate
    v = J(O, 't9_v32_checks.json'); fv = J(O, 't9_far_v32.json'); h6 = J(O, 't9_h6_v32.json'); h4v = J(T8C, 'h4loc_v32.json'); lo = J(O, 't9_lemmaO_v32.json')
    for k in ('cone_slope', 'cone_h2', 'cone_N3'): c(f'v3.2: {k} data <= threshold', lambda k=k: v['checks'][k]['value'], '<=', v['checks'][k]['threshold'] if v else None)
    for k in ('h2', 'h3', 'h4', 'slope', 'L2'): c(f'v3.2: far {k} (physical) <= threshold', lambda k=k: fv[k], '<=', fv['thr'][k] if fv else None)
    for k in ('slope_cone', 'slope_far', 'cone_K2', 'cone_K3', 'far_K2', 'far_K3'): c(f'v3.2: closure margin {k} >= 1', lambda k=k: h6['closure_margins'][k], '>=', 1.0)
    c('v3.2: Lemma O t* margin, cone >= 1', lambda: h6['cone']['tstar_margin'], '>=', 1.0)
    c('v3.2: Lemma O t* margin, far (0.15, 0.3) >= 1', lambda: h6['far']['tstar_margin'], '>=', 1.0)
    c('v3.2: Lemma O H4 far at M = 0.15 (centre 376) <= 1', lambda: h6['far']['H4_ratio_at_M'], '<=', 1.0)
    c('v3.2: Lemma O H1 cone sup|h0| <= 0.19', lambda: lo['sup_h0_max']['cone'], '<=', 0.19)
    lf = J(O, 't9_lemmaO_far.json')
    c('v3.2: Lemma O H1 far family (arc by arc): sup|h0| <= H', lambda: float(lf['sup_h0']), '<=', float(lf['H']) if lf else None)
    c('v3.2: Lemma O H1 far family: all internal checks hold (1 = yes)', lambda: 1.0 if lf['ok'] in (True, 'True') and all(lf['checks'].values()) else 0.0, '>=', 1.0)
    # --- informational
    t1c = J(O, 't1clock.json')
    info('t1 (good-time enclosure, theorem 4: [4.5978574775, 4.5978586144])', lambda: t1c['t1'])
    info('Lemma O H6 (v3.1, step 3 as specified) fails: min closure margin of t9_h6.json', lambda: (J(O, 't9_h6.json') or {}).get('min_closure_margin'))
    info('Lemma O H1 (v3.1, informational): number of failing far cylinders', lambda: len((J(O, 't9_lemmaO.json') or {}).get('fails') or []))
    info('far H1 / axis crossing / other stages added by later cleanups', lambda: sorted(f for f in os.listdir(O) if f.startswith(('t9_axis', 't9_far_h1', 't9_far_H1'))) if os.path.isdir(O) else None)
    return R, I
