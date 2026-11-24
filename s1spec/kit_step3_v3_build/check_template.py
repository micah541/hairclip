#!/usr/bin/env python3
"""Compare kit_out/ with the reference values of the original runs (reference/) and print a PASS / FAIL / WARN / SKIP table.   (kit_step3 v3, v3.3 chain)

Status:  PASS  claim holds AND the recomputed value reproduces the reference under the rule below
         FAIL  claim violated, or output differs from the reference beyond the rule, or output missing for a part that ran
         WARN  claim holds and the output differs from the reference only in the 9th digit or later (platform float noise), or (old kit rule) a typed decimal lies
               on the wrong side of the recomputed certified value by a rounding amount
         SKIP  part not run (no output in kit_out/)
Rules:
  V1  (v3.3 parts lemmaSA_base, lemmaSA, lemmaU) every JSON output is compared with the shipped reference of the proof tree BYTE FOR BYTE after (a) dropping the
      provenance `args` block (top-level `args` holding `argv`/`inputs`; it records absolute paths, cwd and script hashes) and (b) reading path-valued strings by file
      name only (a few outputs quote the path of an input file).  Identical -> PASS.  Not identical but only numeric leaves differ and max relative difference <= 1e-9
      -> WARN (print it); anything else -> FAIL.
  V2  every quoted number of constants_needed_v34.json (the numbers lemmaU_v34.tex cites; @@NUMBERS_NOTE@@) is re-evaluated, exactly (rational arithmetic), on the
      RERUN outputs: a bound must be >= the certified upper end, a margin <= the certified lower end.
  V3  the shipped reference files themselves carry the sha256 pinned in V33_INDEX.md 4.6 / 4.7 (EXPECT_SHA below); reported even if nothing was run.
  The old rules of kit v2 apply to the parts that did not change (verifier tables, ibp, lemmaO):
#@@RULES_OLD@@
usage: python check.py [--detail] [--help]     (also writes kit_out/check_report.txt)"""
import os, sys, json, re, math
from fractions import Fraction
from decimal import Decimal
import hashlib

#@@HELPERS@@

EXPECT_SHA = @@EXPECT_SHA@@

def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

# ------------------------------------------------------------------------------------------------------------------ V3: pinned reference hashes
def check_refs():
    items = []
    for dp, dn, fs in os.walk(REF):
        for f in sorted(fs):
            r = os.path.relpath(os.path.join(dp, f), KIT).replace(os.sep, '/')
            if f in EXPECT_SHA and not r.startswith(('reference/ibp', 'reference/lemmaO')): items.append((r, f))
    items += [('step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json', 'rate_table_fin12f4_a0.00015.json'),
              ('step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json', 'final_alpha_points_5e-8.json'),
              ('step4/p2/out_v33/b0_dprof.json', 'b0_dprof.json'), ('step4/p2/out_v33/n1_chart.json', 'n1_chart.json'),
              ('seeds/u1_ghat.json', 'u1_ghat.json'), ('seeds/zone3_interface.v34.json', 'zone3_interface.v34.json')]
    bad = []
    for r, f in items:
        p = os.path.join(KIT, *r.split('/'))
        if not os.path.exists(p): bad.append('missing ' + r)
        elif sha_file(p) != EXPECT_SHA[f]: bad.append('%s: %s... vs pinned %s...' % (r, sha_file(p)[:12], EXPECT_SHA[f][:12]))
    add('PASS' if not bad else 'FAIL', 'V3: shipped reference / seed / pinned-input files carry the sha256 pinned in V33_INDEX 4.6-4.7 (A_interface 78fe2f8e, U_cert_v35 960809fd, zone3 45e8d4dc, manifest 8bc47e83, b0_dprof 7fdf378c, n1_chart 341ed1f4, rate table b20e8715, ...)',
        '%d files checked' % len(items) if not bad else '; '.join(bad))

# ------------------------------------------------------------------------------------------------------------------ V1: JSON comparison
SKIPK = ('seconds', 'cpu', 'runtime_s')
PATHLIKE = re.compile(r'^[^\s"\']*[/\\][^\s"\']*\.(json|pkl|py|txt|out|log)$')
def norm(o, top=True):
    if isinstance(o, dict):
        d = {}
        for k, v in o.items():
            if top and k == 'args' and isinstance(v, dict) and ('argv' in v or 'inputs' in v): continue
            if k in SKIPK: continue
            d[k] = norm(v, False)
        return d
    if isinstance(o, list): return [norm(x, False) for x in o]
    if isinstance(o, str) and PATHLIKE.match(o): return re.split(r'[/\\]', o)[-1]
    return o
def cmp_json(label, newp, refp):
    n = jl(newp)
    if n is None: add('FAIL', label, 'output missing: ' + os.path.relpath(newp, KIT)); return None
    r = jl(refp)
    if r is None: add('FAIL', label, 'reference missing: ' + os.path.relpath(refp, KIT)); return n
    a, b = norm(n), norm(r)
    if json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True):
        add('PASS', label, 'identical to the reference apart from the args block'); return n
    mx, bad, worst = tree_diff(a, b)
    if bad == 0 and mx <= 1e-9: add('WARN', label, 'numeric leaves differ, max relative difference %.2e (platform float noise?)' % mx)
    else: add('FAIL', label, ('%d structural/string differences, first: %s' % (bad, worst)) if bad else 'max relative difference %.2e: %s' % (mx, worst))
    return n

# ------------------------------------------------------------------------------------------------------------------ V2: quoted numbers on the rerun outputs
def interval(v):
    if isinstance(v, dict) and 'arb' in v:
        a = v['arb']
        if isinstance(a, list):
            m = Fraction(int(a[0])) * (Fraction(2) ** a[1]); r = Fraction(a[2]); return m - r, m + r
        s = a.strip('[]')
        if '+/-' in s:
            m, r = s.split('+/-'); m = Fraction(Decimal(m.strip())); r = Fraction(Decimal(r.strip())); return m - r, m + r
        x = Fraction(Decimal(s)); return x, x
    if isinstance(v, (int, float)) and not isinstance(v, bool): x = Fraction(v); return x, x
    if isinstance(v, str): x = Fraction(Decimal(v)); return x, x
    raise TypeError(repr(v)[:80])
def resolve(d, path):
    parts = path.split('.') if path else []
    i = 0; cur = d
    while i < len(parts):
        if isinstance(cur, list): cur = cur[int(parts[i])]; i += 1; continue
        if not isinstance(cur, dict): raise KeyError(path)
        for j in range(len(parts), i, -1):
            k = '.'.join(parts[i:j])
            if k in cur: cur = cur[k]; i = j; break
        else: raise KeyError('%s (at %s)' % (path, '.'.join(parts[i:])))
    return cur
CNP = os.path.join(KIT, 'writeup', 'lemmaU', 'bulk', 'text', 'constants_needed_v34.json')
def check_numbers(files, part, label):
    """files: key -> rerun output path"""
    CN = jl(CNP); D = {}; miss = [k for k, p in files.items() if not os.path.exists(p)]
    if miss: add('SKIP', label, 'outputs missing: ' + ', '.join(miss)); return
    for k, p in files.items(): D[k] = jl(p)
    bad = []; n = 0
    for q, f, k, d in CN['numbers']:
        if f not in D: continue
        n += 1
        try: lo, hi = interval(resolve(D[f], k))
        except Exception as e: bad.append('%s %s:%s unresolved (%s)' % (q, f, k, e)); continue
        Q = Fraction(Decimal(q)); ok = Q >= hi if d == 'up' else Q <= lo
        if not ok: bad.append('%s %s:%s certified [%.10g, %.10g]' % (q, f, k, float(lo), float(hi)))
    add('PASS' if not bad else 'FAIL', label, '%d quoted numbers hold on the rerun outputs' % n if not bad else '%d of %d fail: %s' % (len(bad), n, '; '.join(bad[:6])))

# ------------------------------------------------------------------------------------------------------------------ verifier (fin12f4)
def check_verifier():
    V, R = os.path.join(OUT, 'verifier'), os.path.join(REF, 'verifier')
    if not part_present('verifier'): add('SKIP', 'verifier (fin12f4)', 'not run'); return
    tc, tcr = tl(os.path.join(V, 'tail_cert_log.txt')), tl(os.path.join(R, 'tail_cert_log.txt'))
    if tc is None: add('FAIL', 'tail_cert PASS', 'tail_cert_log.txt missing')
    else:
        same = tc.strip() == tcr.strip(); ok = 'TAIL_CERT PASS' in tc
        add('PASS' if ok and same else 'FAIL', 'tail_cert.py: TAIL_CERT PASS (all structural degree conditions)', 'log text identical to reference' if same else 'log differs from reference')
    ref_sum = jl(os.path.join(R, 'final_alpha_points_5e-8.json'))
    # typed inputs of Lemma S (run_F10p_ibp.py of the base chain, fin12f4 slopes x4): quoted decimals
    typed = {'1.5e-4': dict(sqLam='5.3075e-8', rho='2.693e-3', J=['6.0273e-12', '8.4673e-12', '1.8708e-12'], J4='6.2362e-8',
                           slopes=['3.4548e-18', '4.8534e-18', '1.0724e-18'], sJ4='3.5746e-14', abar='2.19103e-7', beta='3.33811e-9')}
    for key, tag in (('1.5e-4', '0.00015'), ('1.25e-4', '0.000125')):
        t = jl(os.path.join(V, f'rate_table_fin12f4_a{tag}.json')); tr = jl(os.path.join(R, f'rate_table_fin12f4_a{tag}.json')); lg = tl(os.path.join(V, f'verify_log_fin12f4_a{tag}.txt'))
        nm = f'verifier closes at alpha1={key} (tag fin12f4)'
        if t is None: add('FAIL', nm, 'rate table missing'); continue
        tail = t.get('tail'); rt = tr['tail']
        peak, rho, floor = t['sqrtL_max'], (tail or {}).get('rho'), (tail or {}).get('sqrt_LF')
        m = re.search(r'alpha-condition OK \(margin ([0-9.e+]+)\)', lg or '')
        closes = tail is not None and rho is not None and rho > 0 and 'tail=True' in (lg or '') and t['steps'] and t['alpha1'] == tr['alpha1'] and t['Ly1'] == tr['Ly1']
        inter = closes and peak <= tr['sqrtL_max'] * (1 + 1e-6) and rho >= rt['rho'] * (1 - 1e-6) and floor <= rt['sqrt_LF'] * (1 + 1e-6) and len(t['steps']) == len(tr['steps'])
        page = True
        if key == '1.5e-4': page = holds(peak, '<=', '5.31e-8') and holds(rho, '>=', '2.693e-3') and holds(floor, '<=', '9.88e-10')
        okr = reldiff(peak, tr['sqrtL_max']) <= 1e-6 and reldiff(rho, rt['rho']) <= 1e-6 and reldiff(floor, rt['sqrt_LF']) <= 1e-6
        ex = (peak == tr['sqrtL_max'] and rho == rt['rho'] and floor == rt['sqrt_LF'])
        det = f"steps {len(t['steps'])}, peak sqrtL {peak:.4e}, rho {rho:.4e} (alpha-cond margin {m.group(1) if m else '?'}), floor sqrt(L_F) {floor:.3e}, tail start tau {t['tau_tail_start']:.3f}" if closes else 'tail not certified'
        st = 'PASS' if (inter and page and okr) else ('FAIL' if not (inter and page) else 'WARN')
        add(st, nm + (' (peak<=5.31e-8, rho>=2.693e-3, floor<=9.88e-10)' if key == '1.5e-4' else ''), det + ('; bit-identical to reference' if ex else ('; within 1e-6 of reference' if okr else '; differs from reference by >1e-6 (claims hold)')))
        u = t['unc']; ur = tr['unc']; okj, exj, mxj, note = repro({k: u[k] for k in ('V0sq', 'V1sq', 'V2sq', 'V4sq')}, {k: ur[k] for k in ('V0sq', 'V1sq', 'V2sq', 'V4sq')}, 1e-6)
        add('PASS' if okj else 'FAIL', f'verifier J_0,1,2,4 and slopes (fin12f4: slopes x4 against fin11) at alpha1={key} reproduce reference', note)
    sm = jl(os.path.join(V, 'final_alpha_points_5e-8.json'))
    if sm is None: add('FAIL', 'final_alpha_points_5e-8.json', 'missing')
    else:
        ok, ex, mx, note = repro(sm, ref_sum, 1e-6, skip=('sup_from_reference',))
        add('PASS' if ok else 'FAIL', 'final_alpha_points_5e-8.json (incl. certified sup|(a,b)|, sup|beta| from mod_sups.py) reproduces reference', note)
        s = sm['1.5e-4']
        add('PASS' if (repr(float(s['sup_ab'])) == '2.19103e-07' and repr(float(s['sup_beta'])) == '3.33811e-09') else 'FAIL',
            'v3.3 expected values: sup|(a,b)| = 2.19103e-07, sup|beta| = 3.33811e-09 at alpha1 = 1.5e-4 (unchanged against fin11)', 'recomputed %r, %r' % (s['sup_ab'], s['sup_beta']))
    t = jl(os.path.join(V, 'rate_table_fin12f4_a0.00015.json'))
    if t:
        u = t['unc']; q = typed['1.5e-4']; sm = sm or ref_sum; s = sm['1.5e-4']
        items = [('sqrtLam (F2 --sqLam)', t['sqrtL_max'], '<=', q['sqLam']), ('rho (--rho)', t['tail']['rho'], '>=', q['rho']),
                 ('abar (--abar) vs certified sup|(a,b)|', s['sup_ab'], '<=', q['abar']), ('beta (--beta) vs certified sup|beta|', s['sup_beta'], '<=', q['beta']),
                 ('J4 (--J4)', u['V4sq']['const'], '<=', q['J4']), ('slope J4 (--sJ4)', u['V4sq']['slope_per_tau'], '<=', q['sJ4'])]
        for i, n in enumerate(('V0sq', 'V1sq', 'V2sq')):
            items.append((f'J{i} (--J)', u[n]['const'], '<=', q['J'][i])); items.append((f'slope J{i} (--slopes)', u[n]['slope_per_tau'], '<=', q['slopes'][i]))
        bad = []
        for n, v, op, qs in items:
            if n.startswith(('abar', 'beta')): ok_ = Fraction(repr(float(v))) <= dec(qs)
            else: ok_ = holds(v, op, qs)
            if not ok_: bad.append(f'{n}: certified {v:.6g} vs typed {qs} ({reldiff(v, float(qs)):.1e} rel)')
        if not bad: add('PASS', 'typed Lemma-S inputs of the base chain (sqLam, rho, J, slopes x4, abar, beta) are outward-rounded vs the certified fin12f4 values', f'{len(items)} inputs checked')
        else:
            worst = max(float(re.search(r'\(([0-9.e+-]+) rel', b).group(1)) for b in bad)
            add('WARN' if worst < 1e-2 else 'FAIL', f'typed Lemma-S inputs below the certified verifier value in {len(bad)} of {len(items)} cases (max {worst:.1e} relative)', '; '.join(bad))

# ------------------------------------------------------------------------------------------------------------------ ibp
#@@IBP@@
# ------------------------------------------------------------------------------------------------------------------ repaired Lemma S(A)
BASEN = ['lemmaS_F10p_M0.52', 'cap_F10p', 'lemmaS_rigcap_M0.52', 'weights_exact', 'assembly_arb']
T2N = ['hjw_rep', 'ic_design', 'ic_lo', 'ic_hi', 'ic_probe', 'a3_signs', 'weights_t2', 'weights_t2_fb', 'ramp_rep', 'ramp_rep_pm3ibp', 'outer_L2',
       'asm_design_N2ref', 'asm_design_N2ref4', 'asm_design_N20', 'asm_lo_N2ref', 'asm_hi_N2ref', 'asm_probe_N2ref', 'asm_design_pm3ibp', 'A_interface']
def check_lemmaSA():
    D, R = os.path.join(OUT, 'lemmaSA'), os.path.join(REF, 'lemmaSA', 'out')
    DB, RB = os.path.join(OUT, 'lemmaSA_base'), os.path.join(REF, 'lemmaSA', 'base')
    if not part_present('lemmaSA_base'): add('SKIP', 'lemmaSA base chain (F10p, capture_F2, weights_exact, assembly_arb)', 'not run')
    else:
        for n in BASEN: cmp_json('SA base chain: ' + n + '.json', os.path.join(DB, n + '.json'), os.path.join(RB, n + '.json'))
        # (the claim flags of assembly_arb.json quote the v3.1 page-1 numbers; some are false at the tightened (H2) in the reference too. They are not claims of v3.3: (A) is assembled by t2_A_v35.)
    if not part_present('lemmaSA'): add('SKIP', 'lemmaSA t2_A_v35 chain (-> A_interface.json)', 'not run')
    else:
        for n in T2N: cmp_json('SA t2_A_v35: ' + n + '.json', os.path.join(D, n + '.json'), os.path.join(R, n + '.json'))
        check_numbers({'A': os.path.join(D, 'A_interface.json'), 'A4': os.path.join(D, 'asm_design_N2ref4.json')}, 'lemmaSA',
                      'SA V2: quoted numbers on A_interface.json (flux margins x3.5 at beta_A = 2, core margins, E^A, D^c, K6 inputs ...) and the 4 x N2_ref run')

# ------------------------------------------------------------------------------------------------------------------ Lemma U + zone III
UN = ['ramp_U', 'u1_ghat', 'u2_K', 'u34_primary', 'U_cert_v35_pass1', 'U_cert_v35', 'U_cert_for_zone3_v35', 'zone3_interface', 'U_cert_v35_pass3']
def check_lemmaU():
    D, R = os.path.join(OUT, 'lemmaU'), os.path.join(REF, 'lemmaU')
    if not part_present('lemmaU'): add('SKIP', 'lemmaU / zone3 (circular chain)', 'not run'); return
    src = tl(os.path.join(D, 'A_interface_source.txt'))
    if src and 'shipped reference' in src: add('WARN', 'U: A_interface.json used by assemble_U was the shipped reference, not a rerun', src.strip())
    out = {}
    for n in UN: out[n] = cmp_json('U: ' + n + '.json', os.path.join(D, n + '.json'), os.path.join(R, n + '.json'))
    if (D and os.path.exists(os.path.join(D, 'zone3_interface_intermediate.json'))):
        add('SKIP', 'U: zone3 intermediate (between pass 1 and pass 2)', 'no reference was kept in the proof tree ([TODO], see BUILD_REPORT.md); pass 2 above is its indirect check')
    # fixed points
    u1, r1 = out.get('u1_ghat'), jl(os.path.join(KIT, 'seeds', 'u1_ghat.json'))
    if u1 and r1:
        same = [x['arb'] for x in u1['cases']['primary']['Y_next']] == [x['arb'] for x in r1['cases']['primary']['Y_next']]
        add('PASS' if same else 'FAIL', 'U: seed consistency, primary Y_next of the recomputed u1_ghat.json equals the pinned seed read by ramp_U.py', 'equal' if same else 'differ')
    p2, p3 = out.get('U_cert_v35'), out.get('U_cert_v35_pass3')
    if p2 and p3:
        a, b = norm(p2), norm(p3)
        add('PASS' if json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True) else 'FAIL', 'U: fixed point of the circular chain, pass 3 == pass 2 apart from args (M14 / V33_INDEX 4.6)', '')
    check_numbers({'U': os.path.join(D, 'U_cert_v35.json'), 'U1': os.path.join(D, 'u1_ghat.json'), 'Z': os.path.join(D, 'zone3_interface.json'), 'R': os.path.join(D, 'ramp_U.json')}, 'lemmaU',
                  'U V2: quoted numbers on U_cert_v35.json (Lemma U), u1_ghat.json, zone3_interface.json (zone III / H6), ramp_U.json')

# ------------------------------------------------------------------------------------------------------------------ text checks
def check_text():
    L = os.path.join(OUT, 'logs')
    a, b = tl(os.path.join(L, 'text_check_keys.log')), tl(os.path.join(L, 'text_check_balance.log'))
    if a is None and b is None: add('SKIP', 'text checks (check_keys_v34.py, check_balance.py)', 'not run'); return
    if a is None: add('FAIL', 'check_keys_v34.py', 'log missing')
    else:
        m = re.search(r'RESULT: (PASS|FAIL[^\n]*)', a)
        det = ' | '.join(l for l in a.splitlines() if l.startswith(('cited keys', 'quoted numbers', 'row keys')))
        add('PASS' if (m and m.group(1) == 'PASS') else 'FAIL', 'check_keys_v34.py: lemmaU_v34.tex keys / sha256 pins / quoted numbers / side conditions against the pinned outputs', det or 'no RESULT line')
    if b is None: add('FAIL', 'check_balance.py', 'log missing')
    else: add('PASS' if b.rstrip().endswith('PASS') else 'FAIL', 'check_balance.py lemmaU_v34.tex (braces, environments, labels/refs)', b.strip().splitlines()[-2] if len(b.strip().splitlines()) > 1 else '')

#@@LEMMAO@@
def main():
    if '--help' in sys.argv or '-h' in sys.argv: print(__doc__); return
    for fn in (check_refs, check_verifier, check_ibp, check_lemmaSA, check_lemmaU, check_lemmaO, check_text):
        try: fn()
        except Exception as e:
            import traceback; add('FAIL', fn.__name__, 'check.py error: ' + repr(e) + ' ' + traceback.format_exc().splitlines()[-3].strip())
    ran = [r for r in rows if r[0] != 'SKIP' and not r[1].startswith('V3:')]
    if not ran: rows.append(('FAIL', 'nothing to check', 'kit_out/ has no outputs: run  python run_all.py  first'))
    info = jl(os.path.join(OUT, 'env.json')) or {}
    lines = ['REPRODUCTION CHECK kit_step3 v3  (python %s, python-flint %s, %s)' % (str(info.get('python', '?')).split()[0], info.get('python_flint', '?'), info.get('platform', '?')), '']
    for st, nm, dt in rows:
        lines.append(f'{st:5s} {nm}' + (f'\n        {dt}' if dt else ''))
    cnt = {k: sum(1 for r in rows if r[0] == k) for k in ('PASS', 'FAIL', 'WARN', 'SKIP')}
    lines += ['', f'TOTAL: {cnt["PASS"]} PASS, {cnt["FAIL"]} FAIL, {cnt["WARN"]} WARN, {cnt["SKIP"]} SKIP', '', 'Rules V1-V3 and the old rules R1-R3: see the docstring of check.py (python check.py --help).']
    txt = '\n'.join(lines); print(txt)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'check_report.txt'), 'w', encoding='utf-8') as f: f.write(txt + '\n')
    sys.exit(1 if cnt['FAIL'] else 0)
if __name__ == '__main__': main()
