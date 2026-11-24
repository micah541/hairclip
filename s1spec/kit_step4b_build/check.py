#!/usr/bin/env python3
"""Kit B check: every certified inequality of STEP4_RESULT's margin table, evaluated on the RECOMPUTED outputs (strict, no tolerance), with our numbers next to them.

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
P2 = os.path.join(S4, 'p2'); POUT = os.path.join(P2, 'out'); T8C = os.path.join(P2, 't8se', 'constants')
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
    st = jl(OUT, 'tail_rounds.json')
    if st: add('PASS' if st.get('n3_b1_unchanged') else 'WARN', 'tail round 1 (T7 without nonlinear closure) and round 2 (republished T7) give the same n3_b1', str(st.get('n3_b1_unchanged')))
    return r


def main():
    check_base(); check_t7(); check_claims(); check_compare()
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
    keep = ('switch (A)', '(Y): comparison', '(F) item 3 E3', '(F) item 4', '(F) Ymax', 'boxes: far core x^2', 'boxes: sqrt(F0)', '(H1) sqrt(Lambda)', '(H1) |alpha|', '(H3) window energy E3', '(H5) cone K4', '(H5) Lemma O H4', 'T7: closure constant', 'v3.2: closure margin cone_K2', 'KCH')
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
