#!/usr/bin/env python3
"""Reproduction kit for the step-3 certificates, v3.3 chain (kit_step3 v3).
usage:  python run_all.py [part ...] [--list]
parts (default order): smoke verifier ibp lemmaSA lemmaU lemmaO text        ('zone3' is an alias of 'lemmaU': the two are one circular chain)
  smoke    versions, sha256 of every shipped input (inputs/MANIFEST.sha256), syntax, imports of the verifier / ibp / lemmaO modules   (seconds)
  verifier fin12f4 verifier: tail_cert, closure 1.5e-4 and 1.25e-4 (FLUX0=1.64e-20, FLUX2=8.0e-20), mod_sups
  ibp      modulation lemma (unchanged from kit v2)
  lemmaSA  repaired Lemma S(A): base chain (F10p -> capture_F2 -> weights_exact -> assembly_arb), then the t2_A_v35 chain (-> A_interface.json)
  lemmaU   Lemma U and zone III (circular chain, see README): ramp_U, u1, u2, u34, assemble pass 1, zone3, assemble pass 2, zone3, assemble pass 3
  lemmaO   Lemma O (unchanged from kit v2)
  text     text checks of lemmaU_v34.tex: check_keys_v34.py and check_balance.py (they read the shipped reference outputs)
Each step is a separate Python process (sys.executable -B), run one at a time, with an explicit environment.  Logs and outputs: kit_out/."""
import os, sys, subprocess, time, json, platform, shutil, hashlib, ast

KIT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(KIT, 'kit_out')
VER = os.path.join(KIT, 'step3nl', 'verifier')
IBP = os.path.join(KIT, 'step3nl', 'agent_mod', 'ibp')
LO = os.path.join(KIT, 'writeup', 'lemmaO')
BULK = os.path.join(KIT, 'writeup', 'lemmaU', 'bulk')
T2 = os.path.join(BULK, 't2_A_v35'); BASE = os.path.join(T2, 'base_ibp')
T4 = os.path.join(BULK, 't4_cert'); T5 = os.path.join(BULK, 't5_zone3'); TXT = os.path.join(BULK, 'text')
S4O = os.path.join(KIT, 'step4', 'p2', 'out_v33')
REF = os.path.join(KIT, 'reference')

# environment of the final verifier runs (VERIFIER.md "update 5"; MOD_ALPHA 25, eta 0.01, theta 0.05, S0V 0.0153, S1V 0.6, strip terms in by default) + the v3.3 flux floors (run_v33.sh)
VER_ENV = dict(THETA='0.05', ADD='abcd', FIX='1', S0V='0.0153', S1V='0.6', DELTA='0.2', MOD='ibp', MODA='25', ETA='0.01', FLUX0='1.64e-20', FLUX2='8.0e-20')

def base_env(extra=None):
    e = dict(os.environ)
    e.update(PYTHONHASHSEED='0', PYTHONUTF8='1', PYTHONIOENCODING='utf-8', PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    e.pop('PYTHONPATH', None)
    if extra: e.update(extra)
    return e

TIMES = {}
FAILED = []

def step(part, name, script, args, cwd, env=None, outputs=(), dest=None):
    """run `python -B script args` in cwd; stdout+stderr -> kit_out/logs/<part>_<name>.log.  `outputs`: files created (relative to cwd unless absolute) copied to
    kit_out/<dest or part>/ ; an item may be a pair (src, new_name)."""
    dest = dest or part
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True); os.makedirs(os.path.join(OUT, dest), exist_ok=True)
    log = os.path.join(OUT, 'logs', f'{part}_{name}.log')
    cmd = [sys.executable, '-B', script] + list(args)
    print(f'[{part}] {name}: {" ".join(os.path.basename(c) if i == 2 else c for i, c in enumerate(cmd[1:], 1))}  (cwd {os.path.relpath(cwd, KIT)})', flush=True)
    t0 = time.time()
    with open(log, 'w', encoding='utf-8') as lf:
        p = subprocess.run(cmd, cwd=cwd, env=base_env(env), stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    TIMES[f'{part}/{name}'] = dict(wall_s=round(dt, 1), returncode=p.returncode)
    ok = p.returncode == 0
    for o in outputs:
        o, newname = (o if isinstance(o, tuple) else (o, None))
        src = o if os.path.isabs(o) else os.path.join(cwd, o)
        dst = os.path.join(OUT, dest, newname or os.path.basename(src))
        if os.path.exists(src):
            if os.path.abspath(src) != os.path.abspath(dst): shutil.copy(src, dst)
        elif ok: ok = False; print(f'   MISSING output {o}')
    print(f'   -> {"ok" if ok else "FAILED (rc=%d)" % p.returncode}  {dt:.0f}s wall   log: {os.path.relpath(log, KIT)}', flush=True)
    if not ok: FAILED.append(f'{part}/{name}')
    old = {}
    try: old = json.load(open(os.path.join(OUT, 'timing.json')))
    except Exception: pass
    old.update(TIMES)
    with open(os.path.join(OUT, 'timing.json'), 'w') as f: json.dump(old, f, indent=1)
    return ok

def rm(*paths):
    for p in paths:
        if os.path.isdir(p):
            for f in os.listdir(p):
                if f.endswith('.json') or f.endswith('.tmp'): os.remove(os.path.join(p, f))
        elif os.path.exists(p): os.remove(p)

# ---------------------------------------------------------------------------------------------------------------------------------- smoke
def part_smoke():
    ok = True
    # 1. sha256 of every shipped input
    bad = []
    for l in open(os.path.join(KIT, 'inputs', 'MANIFEST.sha256'), encoding='utf-8'):
        if l.startswith('#') or not l.strip(): continue
        h, r = l.rstrip('\n').split('  ', 1); p = os.path.join(KIT, *r.split('/'))
        if not os.path.exists(p): bad.append(('missing', r))
        elif hashlib.sha256(open(p, 'rb').read()).hexdigest() != h: bad.append(('sha256 differs', r))
    print('   input manifest:', 'ok' if not bad else bad)
    ok &= not bad
    # 2. imports of the cheap parts (as kit v2) and syntax of every .py
    code = {
     'verifier': (VER, 'import constants, bb, rhs_iv, unc, closure, tail_cert, ival; print("verifier imports ok")', VER_ENV),
     'ibp': (IBP, 'import ibp_gen, ibp_consts, ibp_boxes, simplex; print("ibp imports ok")', {}),
     'lemmaO': (LO, 'import lemmaO_v3_interval, lemmaO_consts_v2; print("lemmaO imports ok")', {}),
    }
    for part, (cwd, src, env) in code.items():
        ok &= step('smoke', part, '-c', ['import sys; sys.path.insert(0, "."); ' + src], cwd, env)
    badsyn = []
    for dp, _, fs in os.walk(KIT):
        if os.sep + 'kit_out' in dp or '__pycache__' in dp: continue
        for f in fs:
            if f.endswith('.py'):
                try: ast.parse(open(os.path.join(dp, f), encoding='utf-8').read())
                except SyntaxError as e: badsyn.append((f, e))
    print('   syntax check:', 'ok' if not badsyn else badsyn)
    info = dict(python=sys.version, platform=platform.platform())
    try:
        import flint, numpy, scipy, sympy, mpmath
        info.update(python_flint=flint.__version__, numpy=numpy.__version__, scipy=scipy.__version__, sympy=sympy.__version__, mpmath=mpmath.__version__)
        if flint.__version__ != '0.9.0': print('   WARNING: python-flint %s, the reference runs used 0.9.0 (results may differ in the last digits)' % flint.__version__)
        if sympy.__version__ != '1.14.0': print('   WARNING: sympy %s, the shipped pickles were written with 1.14.0' % sympy.__version__)
    except Exception as e: info['import_error'] = repr(e); ok = False
    os.makedirs(OUT, exist_ok=True); json.dump(info, open(os.path.join(OUT, 'env.json'), 'w'), indent=1); print('   ', info)
    if not ok or badsyn: FAILED.append('smoke')

# ---------------------------------------------------------------------------------------------------------------------------------- verifier
def part_verifier():
    d = os.path.join(OUT, 'verifier'); os.makedirs(d, exist_ok=True)
    step('verifier', 'tail_cert', os.path.join(VER, 'tail_cert.py'), [], d, VER_ENV, ['tail_cert_log.txt'])
    for a in ('1.5e-4', '1.25e-4'):
        tag = '0.00015' if a == '1.5e-4' else '0.000125'
        step('verifier', f'closure_{a}', os.path.join(VER, 'closure.py'), [a], d, dict(VER_ENV, TAG='fin12f4', LY1SQ='5e-8'),
             [f'rate_table_fin12f4_a{tag}.json', f'verify_log_fin12f4_a{tag}.txt'])
    step('verifier', 'mod_sups_summary', os.path.join(VER, 'mod_sups.py'), ['fin12f4', d], d, dict(VER_ENV, TABLE_DIR=d), ['final_alpha_points_5e-8.json'])

# ---------------------------------------------------------------------------------------------------------------------------------- ibp, lemmaO (unchanged)
def part_ibp():
    d = os.path.join(OUT, 'ibp'); os.makedirs(d, exist_ok=True)
    step('ibp', 'symbolic', os.path.join(IBP, 'ibp_sym.py'), [], d, None, ['verify_symbolic.json'])
    step('ibp', 'run_5e-4', os.path.join(IBP, 'ibp_run.py'), ['5e-4', '0.02', '8'], d, None, ['run_0.0005_0.02_n8.json'])
    step('ibp', 'final_5e-4', os.path.join(IBP, 'ibp_final.py'), ['run_0.0005_0.02_n8.json', 'final_5e-4.json'], d, None, ['final_5e-4.json'])

def part_lemmaO():
    d = os.path.join(OUT, 'lemmaO'); os.makedirs(d, exist_ok=True)
    step('lemmaO', 'crosscheck_v2', os.path.join(LO, 'lemmaO_v3_interval.py'), [], d)
    step('lemmaO', 'consts_v3', os.path.join(LO, 'lemmaO_v3_interval.py'), ['run'], d, None, ['lemmaO_v3_consts.json'])
    step('lemmaO', 'far_margin2', os.path.join(LO, 'lemmaO_far_margin2.py'), [], d, None, ['far_margin2.json'])

# ---------------------------------------------------------------------------------------------------------------------------------- repaired Lemma S(A)
BASE_OUTS = [('lemmaS_F10p_M0.52.json', BASE), ('cap_F10p.json', os.path.join(BASE, 'rigor')), ('lemmaS_rigcap_M0.52.json', os.path.join(BASE, 'rigor')),
             ('weights_exact.json', os.path.join(BASE, 'rigor')), ('assembly_arb.json', os.path.join(BASE, 'rigor'))]

def part_lemmaSA():
    """(1) base chain = base_ibp/run_chain_v33.sh of the proof tree: run_F10p_ibp -> capture_F2 -> weights_exact -> assembly_arb (each chdirs to base_ibp itself).
       (2) run_t2.sh of the proof tree (t2_A_v35): hjw -> ic x4 -> signs -> weights x2 -> ramp x2 -> outer -> assembly x7 -> export."""
    for n, d in BASE_OUTS: rm(os.path.join(d, n))                       # never use stale outputs
    rm(os.path.join(T2, 'out'))
    os.makedirs(os.path.join(T2, 'out'), exist_ok=True)
    ok = step('lemmaSA', 'base_F10p', os.path.join(BASE, 'run_F10p_ibp.py'), [], BASE, None, ['lemmaS_F10p_M0.52.json'], dest='lemmaSA_base')
    ok = ok and step('lemmaSA', 'base_capture_F2', os.path.join(BASE, 'rigor', 'capture_F2.py'), [], BASE, None, [os.path.join(BASE, 'rigor', 'cap_F10p.json'), os.path.join(BASE, 'rigor', 'lemmaS_rigcap_M0.52.json')], dest='lemmaSA_base')
    ok = ok and step('lemmaSA', 'base_weights_exact', os.path.join(BASE, 'rigor', 'weights_exact.py'), [], BASE, None, [os.path.join(BASE, 'rigor', 'weights_exact.json')], dest='lemmaSA_base')
    ok = ok and step('lemmaSA', 'base_assembly_arb', os.path.join(BASE, 'rigor', 'assembly_arb.py'), [], BASE, None, [os.path.join(BASE, 'rigor', 'assembly_arb.json')], dest='lemmaSA_base')
    if not ok: print('   lemmaSA: stopping (a base-chain step failed)'); return
    R = lambda n, s, a, outs: step('lemmaSA', n, os.path.join(T2, 'rigor', s), a, T2, None, [os.path.join('out', o) for o in outs], dest='lemmaSA')
    steps = [('hjw', 'hjw_rep.py', [], ['hjw_rep.json'])]
    for nm, h in (('design', '1e-6,5e-6,1e-5'), ('lo', '8e-7,5e-6,1e-5'), ('hi', '1.2e-6,5e-6,1e-5'), ('probe', '3.5e-7,1.1e-6,2.7e-6')):
        steps.append(('ic_' + nm, 'ic_rep.py', ['--h2f', h, '--name', nm], ['ic_%s.json' % nm]))
    steps += [('a3_signs', 'a3_signs.py', [], ['a3_signs.json']),
              ('weights_t2', 'weights_t2.py', [], ['weights_t2.json']),
              ('weights_t2_fb', 'weights_t2.py', ['--seg-box', '20/100', '60/100', '--tag', 'fb'], ['weights_t2_fb.json']),
              ('ramp_rep', 'ramp_rep.py', ['--N', '1200'], ['ramp_rep.json']),
              ('ramp_rep_pm3ibp', 'ramp_rep.py', ['--N', '1200', '--pm3-ibp', '--out', 'out/ramp_rep_pm3ibp.json'], ['ramp_rep_pm3ibp.json']),
              ('outer_L2', 'outer_L2.py', [], ['outer_L2.json'])]
    dprof = os.path.join(S4O, 'b0_dprof.json'); chart = os.path.join(S4O, 'n1_chart.json')
    A = ['--hjw', 'out/hjw_rep.json', '--dprof', dprof, '--chart', chart]
    steps += [('asm_design_N2ref', 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2', '0.1', '0.2', '5.0', '--out', 'out/asm_design_N2ref.json'], ['asm_design_N2ref.json']),
              ('asm_design_N2ref4', 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2', '0.4', '0.8', '20.0', '--out', 'out/asm_design_N2ref4.json'], ['asm_design_N2ref4.json']),
              ('asm_design_N20', 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2', '0', '0', '0', '--out', 'out/asm_design_N20.json'], ['asm_design_N20.json'])]
    for s in ('lo', 'hi', 'probe'):
        steps.append(('asm_%s_N2ref' % s, 'assembly_rep.py', A + ['--ic', 'out/ic_%s.json' % s, '--h3key', s, '--N2', '0.1', '0.2', '5.0', '--out', 'out/asm_%s_N2ref.json' % s], ['asm_%s_N2ref.json' % s]))
    steps += [('asm_design_pm3ibp', 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2', '0.1', '0.2', '5.0', '--ramp', 'out/ramp_rep_pm3ibp.json', '--out', 'out/asm_design_pm3ibp.json'], ['asm_design_pm3ibp.json']),
              ('export_t2', 'export_t2.py', [], ['A_interface.json'])]
    for n, s, a, o in steps:
        if not R(n, s, a, o): print('   lemmaSA: stopping (step %s failed)' % n); return

# ---------------------------------------------------------------------------------------------------------------------------------- Lemma U + zone III
def part_lemmaU():
    """The chain of V33_INDEX 4.6 / 4.4.  Circular dependences: u1_ghat reads t5_zone3/out/ramp_U.json, ramp_U reads the primary Y_next of t4_cert/out/u1_ghat.json
    (seed = the pinned u1_ghat.json); assemble_U pass 1 reads a previous zone3 interface (seed = the pinned v34 file); zone3 reads the Bhat of an assembled U_cert."""
    o4, o5 = os.path.join(T4, 'out'), os.path.join(T5, 'out')
    os.makedirs(o4, exist_ok=True); os.makedirs(o5, exist_ok=True); rm(o4, o5)
    shutil.copy(os.path.join(KIT, 'seeds', 'u1_ghat.json'), os.path.join(o4, 'u1_ghat.json'))
    shutil.copy(os.path.join(KIT, 'seeds', 'zone3_interface.v34.json'), os.path.join(o5, 'zone3_interface.json'))
    aint = os.path.join(T2, 'out', 'A_interface.json'); src = 'rerun of part lemmaSA'
    if not os.path.exists(aint):
        aint = os.path.join(REF, 'lemmaSA', 'out', 'A_interface.json'); src = 'shipped reference (part lemmaSA not run)'
    print(f'   A_interface for assemble_U: {os.path.relpath(aint, KIT)} ({src})', flush=True)
    P = 'lemmaU'
    run = lambda n, s, a, cwd, outs: step(P, n, os.path.join(cwd, s), a, cwd, None, outs)
    seq = [('ramp_U', 'ramp_U.py', ['--fresh'], T5, [os.path.join(o5, 'ramp_U.json')]),
           ('u1_ghat', 'u1_ghat.py', [], T4, [os.path.join(o4, 'u1_ghat.json')]),
           ('u2_K', 'u2_K.py', ['--fresh'], T4, [os.path.join(o4, 'u2_K.json')]),
           ('u34_primary', 'u34_sweeps.py', ['--case', 'primary'], T4, [os.path.join(o4, 'u34_primary.json')])]
    asm = lambda out: ['--A', aint, '--zone3', os.path.join('..', 't5_zone3', 'out', 'zone3_interface.json'), '--out', out]
    uc = os.path.join(o4, 'U_cert_v35.json'); fz = os.path.join('out', 'U_cert_for_zone3_v35.json')
    z3 = lambda: ['--ucert', os.path.join('..', 't4_cert', 'out', 'U_cert_for_zone3_v35.json'), '--ucert_full', uc]
    seq += [('assemble_pass1', 'assemble_U.py', asm('U_cert_v35.json'), T4, [(uc, 'U_cert_v35_pass1.json')]),
            ('for_zone3_pass1', 'make_ucert_for_zone3.py', [os.path.join('out', 'U_cert_v35.json'), fz], T4, []),
            ('zone3_intermediate', 'zone3.py', z3(), T5, [(os.path.join(o5, 'zone3_interface.json'), 'zone3_interface_intermediate.json')]),
            ('assemble_pass2', 'assemble_U.py', asm('U_cert_v35.json'), T4, [uc]),
            ('for_zone3_pass2', 'make_ucert_for_zone3.py', [os.path.join('out', 'U_cert_v35.json'), fz], T4, [os.path.join(o4, 'U_cert_for_zone3_v35.json')]),
            ('zone3_final', 'zone3.py', z3(), T5, [os.path.join(o5, 'zone3_interface.json')]),
            ('assemble_pass3', 'assemble_U.py', asm('U_cert_v35_pass3.json'), T4, [os.path.join(o4, 'U_cert_v35_pass3.json')])]
    for n, s, a, cwd, outs in seq:
        if not run(n, s, a, cwd, outs): print('   lemmaU: stopping (step %s failed)' % n); break
    rm(os.path.join(o4, 'u2_K.partial.json'), os.path.join(o5, 'ramp_U.partial.json'))
    with open(os.path.join(OUT, P, 'A_interface_source.txt'), 'w') as f: f.write(os.path.relpath(aint, KIT).replace(os.sep, '/') + '  # ' + src + '\n')

# ---------------------------------------------------------------------------------------------------------------------------------- text checks
def part_text():
    os.makedirs(os.path.join(OUT, 'text'), exist_ok=True)
    step('text', 'check_keys', os.path.join(TXT, 'check_keys_v34.py'), [], TXT)
    step('text', 'check_balance', os.path.join(TXT, 'check_balance.py'), ['lemmaU_v34.tex'], TXT)

PARTS = dict(smoke=part_smoke, verifier=part_verifier, ibp=part_ibp, lemmaSA=part_lemmaSA, lemmaU=part_lemmaU, lemmaO=part_lemmaO, text=part_text)
ALIAS = dict(zone3='lemmaU')

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if '--list' in sys.argv or '-h' in sys.argv or '--help' in sys.argv: print(__doc__); sys.exit(0)
    args = [ALIAS.get(a, a) for a in args]
    unknown = [a for a in args if a not in PARTS]
    if unknown: print('unknown part(s):', unknown, '- choose from', list(PARTS), '(alias: zone3 = lemmaU)'); sys.exit(2)
    todo = []
    for a in (args or list(PARTS)):
        if a not in todo: todo.append(a)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    for p in todo: PARTS[p]()
    print(f'\nfinished {todo} in {time.time() - t0:.0f}s wall; failed steps: {FAILED or "none"}')
    print('next:  python check.py   (writes kit_out/check_report.txt; paste its text back)')
    sys.exit(1 if FAILED else 0)
