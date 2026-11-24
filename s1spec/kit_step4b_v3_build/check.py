#!/usr/bin/env python3
"""Kit B v3 check (v3.2 chain: STEP4_RESULT margin table;  v3.3 chain: margins of STEP4_RESULT_v33_addendum sections 1, 5, 7): every certified inequality of STEP4_RESULT's margin table, evaluated on the RECOMPUTED outputs (strict, no tolerance), with our numbers next to them.

usage: python check.py [--quick]      (also writes kit_out/check_report.txt)
PASS  the recomputed value satisfies the inequality (exact comparison of the stored double with the bound, no tolerance)
FAIL  it does not, or the output is missing
INFO  reported only (our value next to yours, bit-identity with the reference)
Our numbers (OURS) are the values of STEP4_RESULT.md / the reference run; yours may differ in the last digits on another BLAS/libm (they are all outward-rounded
certified bounds, so only the inequality matters).  --quick checks the slab 0-3 test run only (PARTIAL)."""
import os, sys, json, math, re
try: sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception: pass

KIT = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(KIT, 'kit_out'); S4 = os.path.join(KIT, 'step4'); REF = os.path.join(KIT, 'reference')
P2 = os.path.join(S4, 'p2'); POUT = os.path.join(P2, 'out'); T8C = os.path.join(P2, 't8se', 'constants'); OV33 = os.path.join(P2, 'out_v33')
QUICK = False
NB = 32
rows = []
def add(status, name, new='', ours='', margin='', note=''): rows.append((status, name, new, ours, margin, note))
def jl(*p):
    try:
        with open(os.path.join(*p), encoding='utf-8') as f: return json.load(f)
    except (OSError, ValueError): return None
def tl(*p):
    try:
        with open(os.path.join(*p), encoding='utf-8', errors='replace') as f: return f.read()
    except OSError: return None
def g(x): return f'{x:.4g}' if isinstance(x, float) else str(x)

DOWNGRADE = [False]          # partial hand-off (KIT_NBALL < 32): FAIL rows depending on the ball hull are shown as INFO


def claim(name, val, op, bound, ours=None, unit=''):
    """strict: val <= bound (op '<=') / val >= bound ('>='); val None -> FAIL (missing output)"""
    if val is None: add('FAIL', name, 'missing', g(ours) if ours is not None else '', '', 'output missing'); return False
    ok = (val <= bound) if op == '<=' else (val >= bound) if op == '>=' else (val < bound) if op == '<' else (val > bound)
    m = (bound / val if op in ('<=', '<') and val > 0 else val / bound if bound else float('nan'))
    if DOWNGRADE[0] and not ok: add('INFO', name + ' [partial: %d of 32 balls]' % NB, g(val) + unit, g(ours) + unit if ours is not None else '', f'x{m:.3g}', f'{op} {g(bound)}'); return ok
    add('PASS' if ok else 'FAIL', name, g(val) + unit, g(ours) + unit if ours is not None else '', f'x{m:.3g}', f'{op} {g(bound)}')
    return ok


# ----------------------------------------------------------------------------------------------------------------- T1 / T6 / (SE)
TARGETS = {'R2': 2e-14, 'R1': 5e-13, 'Rs': 5e-11, 'RIs': 1e-10, 'RF0': 1e-7, 'RF1': 1e-7, 'RF2': 1e-7, 'RF3': 1e-7, 'RF4': 1e-5}

def check_base():
    man = jl(S4, 't1', 'polished_v2', 'manifest.json'); ref = jl(REF, 'polished_v2_manifest.json')
    if man is None: add('FAIL', 'T1 polished store', 'missing'); return
    n_exp = 4 if QUICK else 181
    add('PASS' if len(man) == n_exp else 'FAIL', 'T1: slabs in the polished store', str(len(man)), str(n_exp))
    mx = max(m['maxE'] for m in man); add('PASS' if mx < 1e-26 else 'WARN', 'T1: max polish residual |E| (Newton stop 1e-26)', g(mx), '9.7e-27', '', '< 1e-26 (WARN only: the stop test also accepts dz < 1e-37)')
    norms = {}
    for m in man:
        if m.get('norms'):
            for k, v in m['norms'].items():
                if k in TARGETS: norms[k] = max(norms.get(k, 0), v)
    bad = {k: (v, TARGETS[k]) for k, v in norms.items() if not v < TARGETS[k]}
    add('PASS' if not bad and norms else 'FAIL', 'T1: sampled residual norms R2,R1,Rs,RIs,RF0-RF4 below the T5 targets', ', '.join(f'{k} {v:.1e}' for k, v in norms.items()),
        'R2 2.5e-15, R1 6.9e-15, Rs 7.4e-14, RF3 9.6e-9, RF4 2.6e-7', '', 'strict < targets' + (f'; violated: {bad}' if bad else ''))
    if ref:
        k = min(len(ref), len(man)); same = sum(1 for a, b in zip(man[:k], ref[:k]) if a['maxE'] == b['maxE'] and a.get('norms') == b.get('norms'))
        add('INFO', 'T1: polish residuals / norms bit-identical to the reference (slabs compared)', f'{same}/{k}')
    T6 = jl(S4, 't6', 'residuals.json')
    if T6:
        s = T6['summary']; fails = {k: v['fails'] for k, v in s.items() if isinstance(v, dict) and v.get('fails')}
        mxs = {k: v['max'] for k, v in s.items() if k in ('R2', 'R1', 'Rs', 'RIs', 'RF0_opt4n', 'RF3_opt4n', 'RF4_opt4n')}
        add('PASS' if not fails and not T6.get('missing') else 'FAIL', 'T6: all slabs below the T5 targets (R2, Rs, R1, RIs, RF0-RF4)', ', '.join(f'{k} {v:.2e}' for k, v in mxs.items()),
            'R2 3.1e-15, R1 7.4e-15, RF3 8.6e-9, RF4 2.1e-7', '', f'slabs {T6["slabs_range"]}; failures {fails or "none"}')
    se = jl(POUT, 'se_corr.json')
    if se:
        claim('(SE): |a2| of the bump correction', max(r['a2'] for r in se), '<=', 1.6e-14, 1.57e-14)
        add('INFO', '(SE): max |a1|, |h|, |h\'|', f"{max(r['a1'] for r in se):.2e}, {max(r['h'] for r in se):.2e}, {max(r['hp'] for r in se):.2e}", '6.5e-18, 1.4e-14, 9.0e-13')
    else: add('FAIL', '(SE) correction', 'p2/out/se_corr.json missing')


# ----------------------------------------------------------------------------------------------------------------- T7
def check_t7():
    C = jl(S4, 't7', 'constants.json')
    if C is None: add('FAIL', 'T7 constants', 'step4/t7/constants.json missing'); return
    m = C['meta']; cert = C['certificates']; cl = m['closure']
    add('PASS' if len(cert) == 26 and cl.get('all_certified') else 'FAIL', 'T7: 26 certificates, all certified', f'{len(cert)}, all_certified={cl.get("all_certified")}', '26')
    claim('T7: closure constant rho_c (min over certificates) >= 0.0471', cl['rho_c_min'], '>=', 0.0471, 0.04713)
    claim('T7: min c_I >= omega_eff = 0.045 (closure requirement)', min(c['c_I'] for c in cert), '>=', 0.045, 0.05284)
    mu1 = min(v['E_I']['mu1'] for v in C['slabs'].values() if 'E_I' in v); lam = m.get('lambda_max_T8')
    claim('T7: lambda < 2 mu1 (T8 certified lambda_max < 2 min mu1)', lam, '<', 2 * mu1, 0.2520)
    if lam: add('INFO', 'T7: margin 2 mu1/lambda_max - 1', f'{2 * mu1 / lam - 1:.3f}', '0.190')
    nu = [v['E_I'].get('lam_lt_2mu1_margin') for v in C['slabs'].values() if 'E_I' in v]
    if nu: claim('T7: per-slab lambda<2mu1 margins all positive', min(nu), '>', 0.0, 0.190)


# ----------------------------------------------------------------------------------------------------------------- the claims of STEP4_RESULT
def check_claims():
    sys.path.insert(0, KIT)
    import step4_claims as SC
    R, I = SC.extract(S4)
    ours = jl(REF, 'ours.json') or {}
    for name, val, op, bound in R: claim(name, val, op, bound, ours.get(name))
    for name, val in I: add('INFO', name, g(val) if not isinstance(val, list) else '[' + ', '.join(g(x) for x in val) + ']', g(ours.get('info: ' + name)) if ('info: ' + name) in ours else '')


def check_compare():
    """every reference output (SHA-256 list of the proof tree) against this run's"""
    import verify_outputs as VO
    r = VO.compare(); from collections import Counter
    cnt = Counter(t for _, t, _ in r)
    miss = [x for x, t, _ in r if t == 'MISSING']; dif = [x for x, t, _ in r if t == 'DIFFERENT']
    add('INFO', 'outputs vs proof tree: byte-identical / timing-fields-only / provenance-only / different / missing', ' / '.join(str(cnt.get(k, 0)) for k in ('EXACT', 'TIMING', 'PROVENANCE', 'DIFFERENT', 'MISSING')), f'of {len(r)} reference files')
    if dif: add('INFO', 'outputs that differ from the proof tree (numbers, not only bookkeeping)', f'{len(dif)}: ' + ', '.join(dif[:12]) + (' ...' if len(dif) > 12 else ''), '', '', 'a difference is not a failure by itself: every claim above is recomputed on this run\'s own outputs')
    if miss: add('WARN', 'reference outputs not produced by this run', f'{len(miss)}: ' + ', '.join(miss[:8]) + (' ...' if len(miss) > 8 else ''))
    st = jl(OUT, 'tail_rounds.json')                                        # [kit fix] compared with the `args` block dropped (the two files differ in args only; numbers identical)
    if st: add('PASS' if st.get('n3_b1_unchanged') else 'WARN', 'tail round 1 (T7 without nonlinear closure) and round 2 (republished T7) give the same n3_b1 (args dropped)', str(st.get('n3_b1_unchanged')),
               '', '', 'sha256 of the whole files equal: %s' % st.get('sha256_equal'))
    return r


# ----------------------------------------------------------------------------------------------------------------- v3.3 chain (STEP4_RESULT_v33_addendum sections 1, 5, 7)
def check_v33():
    RM = jl(REF, 'v33_manifest.json'); M = jl(OV33, 'FINAL_manifest_v33.json')
    if RM is None: add('FAIL', 'v3.3: reference/v33_manifest.json', 'missing'); return
    if M is None:
        add('FAIL', 'v3.3 chain: p2/out_v33/FINAL_manifest_v33.json', 'missing', '', '', 'the v33_00 ... v33_20 stages of run_all.py did not finish'); return
    res = M.get('results') or {}; ref = RM['results']
    ok_stages = [s['stage'] for s in M['stages']] == [s['stage'] for s in RM['stages']]
    add('PASS' if ok_stages and len(M['stages']) == 21 else 'FAIL', 'v3.3: the manifest has the 21 stages of the proof tree, same names and order', str(len(M['stages'])), '21')
    add('PASS' if M['checks'].get('all_items_ok') is True else 'FAIL', 'v3.3: every result item ok (driver all_items_ok)', str(M['checks'].get('all_items_ok')), 'True')
    add('PASS' if not M['checks'].get('inputs_with_two_versions') else 'FAIL', 'v3.3: no input read in two versions during the chain', str(M['checks'].get('inputs_with_two_versions')))
    dpin = M['checks'].get('v32_pins_differing_from_proof_tree') or []
    add('INFO', 'v3.3: v3.2 files read by the v3.3 stages that differ (bytes) from the proof tree\'s', f'{len(dpin)} of {len(M.get("reused_v32") or [])}',
        '', '', 'expected: the v3.2 chain is regenerated, numerically equal not byte-equal (README)')
    def g_(d, *ks):
        try:
            for k in ks: d = d[k]
            return d
        except (KeyError, IndexError, TypeError): return None
    def margin_floor(name, val, floor):
        if val is None: return
        ok = val >= floor * 0.999
        add('PASS' if ok else 'WARN', name + f' [margin vs addendum x{floor}]', f'x{val:.4g}', f'x{floor}', '', 'WARN only: the inequality itself is the claim above')
    # ---- (H2f), tightened (H2), (H3) on the arb bases (section 5: arb rows supersede the float rows)
    for j, (thr, fl) in enumerate(zip((1e-6, 5e-6, 1e-5), (3.188, 4.983, 4.106))):
        claim(f'v3.3 (H2f) j={j}: with-error value <= {thr:g}', g_(res, 'M9_H2f', 'rows', j, 'with_error'), '<=', thr, g_(ref, 'M9_H2f', 'rows', j, 'with_error')); margin_floor(f'v3.3 (H2f) j={j}', g_(res, 'M9_H2f', 'rows', j, 'margin'), fl)
    for j, (thr, fl) in enumerate(zip((1e-8, 1e-6, 1e-3), (124.78, 2124.9, 641.5))):
        claim(f"v3.3 tightened (H2), w' j={j}: with-error value <= {thr:g}", g_(res, 'M9_H2_tightened', 'rows', j, 'with_error'), '<=', thr, g_(ref, 'M9_H2_tightened', 'rows', j, 'with_error')); margin_floor(f"v3.3 tightened (H2) w' j={j}", g_(res, 'M9_H2_tightened', 'rows', j, 'margin'), fl)
    for j, (thr, fl) in enumerate(zip((1e-8, 1e-6, 1e-3), (27.34, 562.8, 600.9))):
        claim(f'v3.3 tightened (H2), f - U_1 (B0) j={j}: bound <= {thr:g}', g_(res, 'M9_H2_tightened_fU1', 'rows', j, 'bound'), '<=', thr, g_(ref, 'M9_H2_tightened_fU1', 'rows', j, 'bound')); margin_floor(f'v3.3 tightened (H2) f-U_1 j={j}', g_(res, 'M9_H2_tightened_fU1', 'rows', j, 'margin'), fl)
    for j, (thr, fl) in enumerate(zip((2e-7, 1e-5, 1e-5, 1e-5), (5.700, 79.26, 26.50, 7.799))):
        claim(f'v3.3 (H3) window j={j}: with-error value <= {thr:g}', g_(res, 'M9_H3_H2v32', 'rows', 'win', j, 'with_error_arb'), '<=', thr, g_(ref, 'M9_H3_H2v32', 'rows', 'win', j, 'with_error_arb')); margin_floor(f'v3.3 (H3) window j={j}', g_(res, 'M9_H3_H2v32', 'rows', 'win', j, 'margin_arb'), fl)
    for j, thr in enumerate((2e-7, 1e-4, 0.2)):
        claim(f'v3.3 (H2) v3.2 shell j={j} on the arb base: with-error value <= {thr:g}', g_(res, 'M9_H3_H2v32', 'rows', 'shell', j, 'with_error_arb'), '<=', thr, g_(ref, 'M9_H3_H2v32', 'rows', 'shell', j, 'with_error_arb'))
    claim('v3.3 (H3) Ybar limit, binding (E_3) >= 1e-12', g_(res, 'M9_H3_H2v32', 'binding_Ybar_limit'), '>=', 1e-12, g_(ref, 'M9_H3_H2v32', 'binding_Ybar_limit'))
    claim('v3.3 (H3) Ybar limit, binding (E_3): margin >= 15.75 (addendum x15.75)', (g_(res, 'M9_H3_H2v32', 'binding_Ybar_limit') or 0) / 1e-12, '>=', 15.75, (g_(ref, 'M9_H3_H2v32', 'binding_Ybar_limit') or 0) / 1e-12)
    claim('v3.3 far-family Lemma O H4 at every centre >= 490: min margin >= 15.58 (addendum x15.58)', g_(res, 'M9_H4_far', 'min_margin'), '>=', 15.58, g_(ref, 'M9_H4_far', 'min_margin'))
    # ---- fixed-L cone family (section 1; stages 5-8, 17)
    cc = lambda k: g_(res, 'cone_fixedL', 'checks', k, 'value'); cm = lambda k: g_(res, 'cone_fixedL', 'checks', k, 'margin'); cr = lambda k: g_(ref, 'cone_fixedL', 'checks', k, 'value')
    for k, nm, thr, fl in (('slope', 'slope datum (H5)', 0.08, 2.64), ('K2', 'sup|h_0\'\'| = K2 (H5)', 0.03, 3.03), ('N3', '||h_0\'\'\'||_L2 = N3 (H5)', 0.004, 27.3), ('K3', 'H3 jet K3', 0.12, 2012), ('K4', 'H3 jet K4', 0.12, 143)):
        claim(f'v3.3 fixed-L cone {nm} <= {thr:g}', cc(k), '<=', thr, cr(k)); margin_floor(f'v3.3 fixed-L cone {nm}', cm(k), fl)
    claim('v3.3 fixed-L cone Lemma O H4: int (h_0^(5))^2 <= 0.0144', g_(res, 'H4_fixedL', 'total'), '<=', 0.0144, g_(ref, 'H4_fixedL', 'total')); margin_floor('v3.3 fixed-L cone H4', g_(res, 'H4_fixedL', 'margin'), 1119)
    h1 = res.get('H1_fixedL') or {}; h1m = res.get('M9_H1_centre_map') or {}
    claim('v3.3 fixed-L cone Lemma O H1 (graph): sup|h_0| <= H - C_0 sqrt(t_1)', (h1.get('H_max') - h1.get('C0_sqrt_t1') - h1.get('sup_h0_max')) if h1 else None, '>=', 0.0, None)
    add('PASS' if h1m.get('ok') and h1m.get('centre_map_ok') else 'FAIL', 'v3.3 Lemma O H1 with the certified centre-to-arclength map (M9)', str(h1m.get('ok')), 'True')
    claim('v3.3 fixed-L cone H7: ratio <= needed', g_(res, 'H7_fixedL', 'ratio_upper'), '<=', g_(res, 'H7_fixedL', 'needed') or 0.0, g_(ref, 'H7_fixedL', 'ratio_upper')); margin_floor('v3.3 fixed-L cone H7 t*', g_(res, 'H7_fixedL', 'margin_tstar_over_ratio'), 2.30)
    # ---- N1 chart terms (stage 19), F1/F2 sliver (stage 20)
    claim('v3.3 N1 chart: |e| (ramp) < 2^-19', g_(res, 'N1_chart', 'ramp', 'e'), '<', 2.0 ** -19, g_(ref, 'N1_chart', 'ramp', 'e'))
    sa = g_(res, 'N1_chart', 'ramp', 'sup_add', 0); claim('v3.3 N1 chart: sup|chi^add_0| on the ramp <= 2e-7', sa, '<=', 2e-7, g_(ref, 'N1_chart', 'ramp', 'sup_add', 0))
    f12 = res.get('F12_sliver') or {}
    add('PASS' if f12.get('ok') and not f12.get('quoted_numbers_failing') else 'FAIL', 'v3.3 f12_sliver: all checks pass and no quoted number fails', f'{f12.get("n_checks")} checks; failing quoted numbers: {f12.get("quoted_numbers_failing")}', '326 checks')
    nm = f12.get('named') or {}
    claim('v3.3 f12 ramp M1 (A) >= 1.1847 (quoted x1.1847)', g_(nm, 'ramp_M1', 'new'), '>=', 1.1847, g_(ref, 'F12_sliver', 'named', 'ramp_M1', 'new'))
    claim('v3.3 f12 K6 >= 1.189', g_(nm, 'K6', 'new'), '>=', 1.189, g_(ref, 'F12_sliver', 'named', 'K6', 'new'))
    claim('v3.3 f12 H6 min margin >= 1.860', g_(nm, 'H6', 'new'), '>=', 1.860, g_(ref, 'F12_sliver', 'named', 'H6', 'new'))
    fx = g_(nm, 'flux_beta2_4N2ref', 'new') or [None, None]
    claim('v3.3 f12 flux FLUX0/K0 (beta 2, 4 N2_ref) >= 3.365', fx[0], '>=', 3.365, (g_(ref, 'F12_sliver', 'named', 'flux_beta2_4N2ref', 'new') or [None])[0])
    claim('v3.3 f12 flux FLUX2/K2 (beta 2, 4 N2_ref) >= 5.955', fx[1], '>=', 5.955, (g_(ref, 'F12_sliver', 'named', 'flux_beta2_4N2ref', 'new') or [None, None])[1])
    claim('v3.3 f12 M1\' = 28 sI2\' <= 0.4328 (what (A) reads)', 28 * f12['sI2_prime'] if f12.get('sI2_prime') else None, '<=', 0.4328, 28 * g_(ref, 'F12_sliver', 'sI2_prime') if g_(ref, 'F12_sliver', 'sI2_prime') else None)
    for k in ('ok',):
        for nm_, v in res.items():
            if isinstance(v, dict) and k in v and v[k] is not True: add('FAIL', f'v3.3 result item {nm_}: ok', str(v[k]), 'True')
    # ---- outputs against the proof tree (args dropped)
    import verify_v33 as V3
    r = V3.compare(); from collections import Counter
    cnt = Counter(t for _, t, _ in r)
    add('INFO', 'v3.3 outputs vs proof tree: byte-identical / equal with args+timing dropped / numbers differ / structure differs / missing',
        ' / '.join(str(cnt.get(k, 0)) for k in ('EXACT', 'ARGS', 'NUMERIC', 'DIFFERENT', 'MISSING')), f'of {len(r)} reference files')
    num = [(x, n) for x, t, n in r if t == 'NUMERIC']
    if num:
        w = max(num, key=lambda q: float(q[1].split('max rel ')[1])); add('INFO', 'v3.3 outputs with differing numbers: largest relative difference', w[1].split('max rel ')[1], '', '', w[0])
    dif = [x for x, t, _ in r if t == 'DIFFERENT']; mis = [x for x, t, _ in r if t == 'MISSING']
    if dif: add('WARN', 'v3.3 outputs whose structure/text differs from the proof tree', f'{len(dif)}: ' + ', '.join(dif[:8]))
    if mis: add('WARN', 'v3.3 reference outputs not produced by this run', f'{len(mis)}: ' + ', '.join(mis[:8]))


def main():
    check_base(); check_t7(); check_claims(); check_compare(); check_v33()
    nf = sum(r[0] == 'FAIL' for r in rows); npass = sum(r[0] == 'PASS' for r in rows); nw = sum(r[0] == 'WARN' for r in rows)
    verdict = 'TOTAL PASS' if nf == 0 else 'TOTAL FAIL'
    L = ['=' * 118, 'KIT 4B check: certified inequalities of STEP4_RESULT, recomputed  (strict comparisons, no tolerance)', '=' * 118]
    for st, nm, new, ours, mg, note in rows:
        L.append(f'{st:5s}  {nm}')
        L.append(f'         yours: {new}' + (f'   | ours: {ours}' if ours != '' else '') + (f'   | margin {mg}' if mg else '') + (f'   | {note}' if note else ''))
    L += ['-' * 118, 'Rule: every PASS/FAIL line is an inequality evaluated exactly on the recomputed double (no tolerance); "ours" are the proof tree\'s values.', '',
          '=' * 40 + ' SUMMARY (paste this) ' + '=' * 40, f'{verdict}   [{npass} PASS, {nf} FAIL, {nw} WARN]']
    e = jl(OUT, 'env.json') or {}
    if e: L.append(f"env: Python {e.get('python', '?').split()[0]} | {e.get('platform', '?')} | flint {e.get('versions', {}).get('flint')} numpy {e.get('versions', {}).get('numpy')} | jobs {e.get('jobs')}")
    for st, nm, new, ours, mg, note in rows:
        if st == 'FAIL': L.append(f'FAIL: {nm}: yours {new} {note}')
        if st == 'INFO' and nm.startswith('outputs vs proof tree'): L.append(f'{nm}: {new}  {ours}')
    keep = ('switch (A)', '(Y): comparison', '(F) item 3 E3', '(F) item 4', '(F) Ymax', 'boxes: far core x^2', 'boxes: sqrt(F0)', '(H1) sqrt(Lambda)', '(H1) |alpha|', '(H3) window energy E3', '(H5) cone K4', '(H5) Lemma O H4', 'T7: closure constant', 'v3.2: closure margin cone_K2', 'KCH', 'v3.3 (H2f)', 'v3.3 tightened', 'v3.3 (H3)', 'v3.3 fixed-L', 'v3.3 far-family', 'v3.3 f12', 'v3.3 N1')
    for st, nm, new, ours, mg, note in rows:
        if nm.startswith(keep): L.append(f'  {nm}: yours {new} (ours {ours}) margin {mg}')
    L.append('=' * 103)
    txt = '\n'.join(L); print(txt)
    try:
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, 'check_report.txt'), 'w', encoding='utf-8') as f: f.write(txt + '\n')
    except OSError: pass
    return 0 if nf == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
