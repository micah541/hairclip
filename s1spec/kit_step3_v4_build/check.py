#!/usr/bin/env python3
"""Compare this kit's outputs with the reference outputs of the proof tree and print PASS / FAIL / WARN / SKIP lines and a SUMMARY.   (kit_step3 v4)

Status:  PASS  claim holds and the output reproduces the reference under the rule below
         FAIL  claim violated, output differs beyond the rule, or output missing for a part that ran
         WARN  claim holds, numbers differ only at relative <= 1e-9 (platform float noise); or a guard note (reads outside the kit)
         SKIP  part not run
Rule C (comparison of a JSON output with its reference): equal after normalising PROVENANCE only, namely
  * every `args` block (argv, env, cwd, python/flint versions, input paths with sha256) and every timing field (seconds, cpu, wall, runtime_s, cpu_s) is dropped;
  * path-valued strings and path-valued dict KEYS are read by file name (the proof tree's absolute paths vs this kit's);
  * sha256-valued provenance fields (a key ending in `_sha`, `_sha256`, `sha256`, or a 64-hex value in a path-keyed `inputs`/`pinned`/`outputs` dict) are
    replaced by a placeholder: they hash input FILES whose own provenance (args) legitimately differs in a kit (the kit v3 ramp_U `u1_ghat_sha` FAIL).
  Equal -> PASS; only numbers differ with max relative difference <= 1e-9 -> WARN; anything else -> FAIL.
usage: python check.py      (also writes kit_out/check_report.txt)"""
import os, sys, json, re, math, glob
KIT = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(KIT, 'kit_out'); REF = os.path.join(KIT, 'reference')
NC = os.path.join(KIT, 'writeup', 'lemmaU', 'bulk', 'fix13_v37chain')
rows = []
def add(st, name, detail=''): rows.append((st, name, detail))
def jl(p):
    try:
        with open(p, encoding='utf-8') as f: return json.load(f)
    except Exception: return None
def tl(p):
    try:
        with open(p, encoding='utf-8') as f: return f.read()
    except Exception: return None
TIMING = ('seconds', 'cpu', 'wall', 'runtime_s', 'cpu_s', 'wall_s', 'elapsed', 'time_s', 'secs')
PATHLIKE = re.compile(r'^(?:[A-Za-z]:)?[^\s"\']*[/\\][^\s"\']*\.(json|pkl|py|txt|out|log|tex|md|npz|sh)$')
HEX64 = re.compile(r'^[0-9a-f]{64}$')
ABSP = re.compile(r'^(?:[A-Za-z]:)?[/\\][^\s"\']*$')
NPR = re.compile(r'^(?:[A-Za-z]:)?/(?:[^/\s]+/)*?(?=(?:step4|writeup|step3nl_v4c|step3nl_v4b|step3nl_v33|step3nl|step3|lin|s1_cap|s1_tight|reference)/)')
def base(s): return re.split(r'[/\\]', s)[-1]
def norm(o, key=None):
    if isinstance(o, dict):
        d = {}
        for k, v in o.items():
            if k == 'args' or k in TIMING: continue
            kk = base(k) if isinstance(k, str) and PATHLIKE.match(k) else k
            if isinstance(v, str) and HEX64.match(v) and (str(k).endswith(('_sha', 'sha256')) or kk != k or key in ('inputs', 'pinned', 'outputs')): d[kk] = '<sha256>'; continue
            d[kk] = norm(v, k)
        return d
    if isinstance(o, list): return [norm(x, key) for x in o]
    if isinstance(o, str) and PATHLIKE.match(o): return base(o)
    if isinstance(o, str) and ABSP.match(o): return NPR.sub('@/', o.replace('\\', '/'))      # absolute directory paths: modulo the install root
    return o
def tree_diff(a, b, path='', acc=None):
    acc = acc if acc is not None else [0.0, 0, '']
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b): acc[1] += 1; acc[2] = acc[2] or f'{path}: keys {sorted(set(a) ^ set(b))[:4]}'; return acc
        for k in a: tree_diff(a[k], b[k], f'{path}/{k}', acc)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): acc[1] += 1; acc[2] = acc[2] or f'{path}: length'; return acc
        for i, (x, y) in enumerate(zip(a, b)): tree_diff(x, y, f'{path}/{i}', acc)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if a != b:
            r = abs(a - b) / max(abs(a), abs(b), 1e-300)
            if r > acc[0]: acc[0] = r; acc[2] = acc[2] if acc[1] else f'{path}: {a!r} vs {b!r}'
    elif a != b: acc[1] += 1; acc[2] = acc[2] or f'{path}: {str(a)[:50]!r} vs {str(b)[:50]!r}'
    return acc
def cmp_json(label, newp, refp, required=True):
    n = jl(newp)
    if n is None:
        if required: add('FAIL', label, 'output missing: ' + os.path.relpath(newp, KIT))
        return None
    r = jl(refp)
    if r is None: add('FAIL', label, 'reference missing: ' + os.path.relpath(refp, KIT)); return n
    a, b = norm(n), norm(r)
    if json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True): add('PASS', label, 'equal to the reference (rule C)'); return n
    mx, bad, worst = tree_diff(a, b)
    if bad == 0 and mx <= 1e-9: add('WARN', label, 'only numbers differ, max relative difference %.2e (%s)' % (mx, worst))
    else: add('FAIL', label, ('%d structural/text differences, first %s' % (bad, worst)) if bad else 'max relative difference %.2e at %s' % (mx, worst))
    return n
def cmp_bytes(label, newp, refp):
    if not os.path.exists(newp): add('FAIL', label, 'missing ' + os.path.relpath(newp, KIT)); return
    same = open(newp, 'rb').read() == open(refp, 'rb').read()
    add('PASS' if same else 'FAIL', label, 'byte-identical to the reference' if same else 'bytes differ')
def ran(*logs): return any(os.path.exists(os.path.join(OUT, 'logs', l + '.log')) for l in logs)

def check_verifier():
    V, R = os.path.join(OUT, 'verifier'), os.path.join(REF, 'verifier')
    if not ran('verifier_tail_cert'): add('SKIP', 'verifier fin13v4c', 'not run'); return
    tc = tl(os.path.join(V, 'tail_cert_log.txt')); tcr = tl(os.path.join(R, 'tail_cert_log.txt'))
    add('PASS' if tc and 'TAIL_CERT PASS' in tc and tc.strip() == (tcr or '').strip() else 'FAIL', 'verifier: tail_cert PASS, log identical to the reference', '')
    for a, tg in (('1.5e-4', '0.00015'), ('1.25e-4', '0.000125')):
        t = cmp_json(f'verifier: rate_table_fin13v4c_a{tg}.json (closure alpha1 = {a})', os.path.join(V, f'rate_table_fin13v4c_a{tg}.json'), os.path.join(R, f'rate_table_fin13v4c_a{tg}.json'))
        if t:
            lg = tl(os.path.join(V, f'verify_log_fin13v4c_a{tg}.txt')) or ''
            closes = t.get('tail') is not None and t['tail'].get('rho', 0) > 0 and 'tail=True' in lg
            bh = (t.get('box_hypotheses') or {}).get('max_ratio') or {}
            add('PASS' if closes and bh and max(bh.values()) <= 1 else 'FAIL', f'verifier closes at alpha1 = {a}: tail certified, rho > 0, every box-hypothesis ratio <= 1',
                f"peak sqrtL {t['sqrtL_max']:.6e}, rho {t['tail']['rho']:.6e}, max box ratio {max(bh.values()) if bh else None}")
    s = cmp_json('verifier: mod_sups final_alpha_points_5e-8.json (sup|(a,b)|, sup|beta|, J, rho)', os.path.join(V, 'modsups_fin13v4c', 'final_alpha_points_5e-8.json'), os.path.join(R, 'final_alpha_points_5e-8.json'))
    if s: add('PASS' if s['1.5e-4']['sup_ab'] == 2.14039e-07 and s['1.5e-4']['sup_beta'] == 2.35488e-09 else 'FAIL', 'verifier: sup|(a,b)| = 2.14039e-7, sup|beta| = 2.35488e-9 at alpha1 = 1.5e-4 (v38 values)',
              f"{s['1.5e-4']['sup_ab']!r}, {s['1.5e-4']['sup_beta']!r}")

def check_cert():
    if not ran('cert_strip_cubic_v4c'): add('SKIP', 'cert (strip_cubic_v4c, e0_A2_budget_v4c)', 'not run'); return
    for f in ('strip_cubic_v4c', 'e0_A2_budget_v4c'): cmp_json(f'cert: {f}.json', os.path.join(OUT, 'cert', f + '.json'), os.path.join(REF, 'cert', f + '.json'))

def check_n1():
    if not ran('n1_n1_chart'): add('SKIP', 'n1 chart terms', 'not run'); return
    n = cmp_json('n1: n1_chart_fin13v4c.json', os.path.join(OUT, 'n1', 'n1_chart_fin13v4c.json'), os.path.join(REF, 'chain', 'n1', 'n1_chart_fin13v4c.json'))
    if n: add('PASS' if n['ramp']['e'] < 2 ** -19 and n['ramp']['sup_add'][0] <= 2e-7 else 'FAIL', 'n1: ramp e < 2^-19 and sup chi_add_0 <= 2e-7', f"e {n['ramp']['e']:.3e}, sup_add0 {n['ramp']['sup_add'][0]:.3e}")

T2N = ['hjw_rep', 'ic_design', 'ic_lo', 'ic_hi', 'ic_probe', 'a3_signs', 'weights_t2', 'weights_t2_fb', 'ramp_rep', 'ramp_rep_pm3ibp', 'outer_L2',
       'asm_design_N2ref', 'asm_design_N2ref4', 'asm_design_N20', 'asm_lo_N2ref', 'asm_hi_N2ref', 'asm_probe_N2ref', 'asm_design_pm3ibp', 'A_interface']
def check_lemmaSA():
    if not ran('lemmaSA_base_F10p'): add('SKIP', 'lemmaSA (base chain, R3, t2_A)', 'not run'); return
    for n in ('lemmaS_F10p_M0.52',): cmp_json('SA base: ' + n + '.json', os.path.join(OUT, 'lemmaSA_base', n + '.json'), os.path.join(REF, 'chain', 'base_ibp', n + '.json'))
    for n in ('cap_F10p', 'lemmaS_rigcap_M0.52', 'weights_exact', 'assembly_arb'):
        cmp_json('SA base: ' + n + '.json', os.path.join(OUT, 'lemmaSA_base', n + '.json'), os.path.join(REF, 'chain', 'base_ibp', 'rigor', n + '.json'))
    for n in ('r3_forcing', 'r3_symcheck'): cmp_json('SA R3: ' + n + '.json', os.path.join(OUT, 'lemmaSA_r3', n + '.json'), os.path.join(REF, 'chain', 'r3', n + '.json'))
    out = {n: cmp_json('SA t2_A: ' + n + '.json', os.path.join(OUT, 'lemmaSA', n + '.json'), os.path.join(REF, 'chain', 't2_A', n + '.json')) for n in T2N}
    A = out.get('A_interface')
    if A: add('PASS' if all(A['res'][s]['rule_reproduces_EA_ref'] is True for s in ('tight', 'v31')) else 'FAIL', 'SA: A_interface rule_reproduces_EA_ref (tight, v31)', '')
    for nm in ('asm_design_N2ref', 'asm_design_N2ref4'):
        a = out.get(nm)
        rr = jl(os.path.join(REF, 'chain', 't2_A', nm + '.json'))
        if a and rr:
            fa, fr = sorted(c[0] for c in a['res']['tight']['claims_failing']), sorted(c[0] for c in rr['res']['tight']['claims_failing'])
            add('PASS' if fa == fr else 'FAIL', f'SA: {nm} tight: failing claims are exactly those of the reference (v3.1 page-1 quotes and the v3.1 S0_B quote; not v3.3+ claims)', '%d: %s' % (len(fa), fa[:3]))

UN = ['u1_ghat', 'u2_K', 'u34_primary', 'U_cert_v36_pass1', 'U_cert_for_zone3_v36_pass1', 'U_cert_v36_pass2', 'U_cert_for_zone3_v36', 'U_cert_v36_pass3']
def check_lemmaU():
    if not ran('lemmaU_u1_ghat'): add('SKIP', 'lemmaU / zone III', 'not run'); return
    o = {n: cmp_json('U: ' + n + '.json', os.path.join(OUT, 'lemmaU', n + '.json'), os.path.join(REF, 'chain', 't4_cert', n + '.json')) for n in UN}
    for n in ('zone3_interface_pass1', 'zone3_interface'): o[n] = cmp_json('zone III: ' + n + '.json', os.path.join(OUT, 'lemmaU', n + '.json'), os.path.join(REF, 'chain', 't5_zone3', n + '.json'))
    cmp_bytes('zone III: ramp_U.json (an input of chain v37, the v35 byte copy; not regenerated)', os.path.join(NC, 't5_zone3', 'out', 'ramp_U.json'), os.path.join(REF, 'chain', 't5_zone3', 'ramp_U.json'))
    p2, p3 = o.get('U_cert_v36_pass2'), o.get('U_cert_v36_pass3')
    if p2 and p3: add('PASS' if json.dumps(norm(p2), sort_keys=True) == json.dumps(norm(p3), sort_keys=True) else 'FAIL', 'U: fixed point, pass 3 == pass 2 (rule C)', '')
    if p3:
        ss = p3['supersolution']; ok = all(ss[k] is True for k in ('K1_ok', 'K2_ok', 'K3_ok', 'N2_le_4N2ref')) and p3['ufine']['ok'] is True
        add('PASS' if ok else 'FAIL', 'U: supersolution K1, K2, K3, N2 <= 4 N2_ref, U-fine requirement', '')
    z = o.get('zone3_interface')
    if z: add('PASS' if z['cone_family']['H6']['at_Bhat_ref']['min_margin'] > 1 else 'FAIL', 'zone III: H6 min margin > 1', str(z['cone_family']['H6']['at_Bhat_ref']['min_margin']))

FX = ['chart4_shell', 'flux2_cert', 'flux2_coef', 'flux2_data', 'flux2_forms', 'flux2_sanity', 'gen_texts', 'inputs']
def check_flux():
    if not ran('flux_flux2_forms'): add('SKIP', 'flux lemma v2 (v37)', 'not run'); return
    o = {n: cmp_json('flux: ' + n + '.json', os.path.join(OUT, 'flux', n + '.json'), os.path.join(REF, 'flux', 'out', n + '.json')) for n in FX}
    cmp_bytes('flux: flux2_forms.pkl', os.path.join(OUT, 'flux', 'flux2_forms.pkl'), os.path.join(REF, 'flux', 'out', 'flux2_forms.pkl'))
    for f in sorted(os.listdir(os.path.join(REF, 'flux', 'texts'))): cmp_bytes('flux texts: ' + f, os.path.join(OUT, 'flux', f), os.path.join(REF, 'flux', 'texts', f))
    c = o.get('flux2_cert')
    if c:
        mc = c['res'].get('margin_combined'); mv = mc.get('lo', mc.get('up')) if isinstance(mc, dict) else mc
        add('PASS' if c['verdict'] == 'PASS' and not c['failing'] else 'FAIL', 'flux: flux2_cert verdict PASS, no failing condition (C8 = U-1 by hand)', f'combined margin {mv}')

def check_proofs():
    if not ran('proofs_U1_sym'): add('SKIP', 'proofs (G1_check, U1_sym, GS_boot)', 'not run'); return
    for n in ('G1_check', 'U1_sym', 'GS_boot_v4c'):
        x = cmp_json('proofs: ' + n + '.json', os.path.join(OUT, 'proofs', n + '.json'), os.path.join(REF, 'proofs', n + '.json'))
        if x: add('PASS' if (x.get('all_pass', x.get('ok')) is True) else 'FAIL', f'proofs: {n} all checks pass', '')

def check_old(part, files):
    d = os.path.join(OUT, part)
    if part == 'lemmaO' and os.path.exists(os.path.join(d, 'lemmaO_v3_consts.json')) and os.path.exists(os.path.join(d, 'far_margin2.json')):   # the reference file is the merged one (kit v2/v3 rule)
        m = dict(jl(os.path.join(d, 'lemmaO_v3_consts.json')), far_margin2=jl(os.path.join(d, 'far_margin2.json')))
        json.dump(m, open(os.path.join(d, 'lemmaO_v3_consts_merged.json'), 'w'), indent=1); files = ['lemmaO_v3_consts_merged.json']
        cmp_json('lemmaO: lemmaO_v3_consts.json + far_margin2.json (merged, as the reference)', os.path.join(d, files[0]), os.path.join(REF, 'lemmaO', 'lemmaO_v3_consts.json')); return
    if not os.path.isdir(d) or not os.listdir(d): add('SKIP', part, 'not run'); return
    for f in files:
        r = os.path.join(REF, part, f)
        if os.path.exists(r): cmp_json(f'{part}: {f}', os.path.join(d, f), r)
        else: add('PASS' if os.path.exists(os.path.join(d, f)) else 'FAIL', f'{part}: {f} produced', '(no reference shipped)')

def check_guard():
    g = glob.glob(os.path.join(OUT, 'guard', 'reads_*.txt'))
    reads = sorted({l.split('\t')[0] for p in g for l in open(p, encoding='utf-8') if l.strip()})
    add('PASS' if not reads else 'WARN', 'guard: files read outside the kit and outside the Python installation / OS directories', f'{len(reads)}' + (': ' + ', '.join(reads[:8]) if reads else ''))
    errs = [l for l in glob.glob(os.path.join(OUT, 'logs', '*.log')) if 'kit guard:' in (tl(l) or '')]
    add('PASS' if not errs else 'FAIL', 'guard: no step tried to write / chdir / start a child outside the kit', ', '.join(os.path.basename(e) for e in errs))

def main():
    if '-h' in sys.argv or '--help' in sys.argv: print(__doc__); return
    for fn in (check_verifier, check_cert, check_n1, check_lemmaSA, check_lemmaU, check_flux, check_proofs,
               lambda: check_old('ibp', ['verify_symbolic.json', 'final_5e-4.json']), lambda: check_old('lemmaO', ['lemmaO_v3_consts.json', 'far_margin2.json']), check_guard):
        try: fn()
        except Exception as e:
            import traceback; add('FAIL', getattr(fn, '__name__', 'check'), 'check.py error: %r %s' % (e, traceback.format_exc().splitlines()[-3].strip()))
    info = jl(os.path.join(OUT, 'env.json')) or {}
    lines = ['REPRODUCTION CHECK kit_step3 v4  (python %s, python-flint %s, %s)' % (str(info.get('python', '?')).split()[0], info.get('python_flint', '?'), info.get('platform', '?')), '']
    for st, nm, dt in rows: lines.append(f'{st:5s} {nm}' + (f'\n        {dt}' if dt else ''))
    cnt = {k: sum(1 for r in rows if r[0] == k) for k in ('PASS', 'FAIL', 'WARN', 'SKIP')}
    lines += ['', 'SUMMARY: %d PASS, %d FAIL, %d WARN, %d SKIP' % (cnt['PASS'], cnt['FAIL'], cnt['WARN'], cnt['SKIP'])]
    txt = '\n'.join(lines); print(txt); os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'check_report.txt'), 'w', encoding='utf-8') as f: f.write(txt + '\n')
    sys.exit(1 if cnt['FAIL'] else 0)
if __name__ == '__main__': main()
