#!/usr/bin/env python3
"""Reproduction kit for the step-3 certificates on verifier fin13v4c, (A) chain v37, flux lemma v2 (v37) and the interface / proof checks (kit_step3 v4).
usage:  python run_all.py [part ...] [--list]
parts (default order): smoke verifier cert n1 lemmaSA lemmaU flux proofs ibp lemmaO
  smoke     versions, sha256 of every shipped input (inputs/MANIFEST.sha256), syntax of every .py                                  (seconds)
  verifier  fin13v4c: tail_cert, closure alpha1 = 1.5e-4 and 1.25e-4, mod_sups (environment of step3nl_v4c/verifier/run_v4c.sh)      (~12 min)
  cert      the two new fin13v4c certificates, rerun from the shipped cert outputs: strip_cubic_v4c (D-3, H8), e0_A2_budget_v4c (A-1)   (~10 s)
  n1        N1 chart terms on the fin13v4c sups (fix13_v37chain/n1/n1_chart.py)                                                   (~10 s)
  lemmaSA   (A): base chain (run_F10p_ibp, capture_F2, weights_exact, assembly_arb), R3 constants, t2_A chain (-> A_interface.json)    (~4 min)
  lemmaU    Lemma U + zone III as run_chain_v37.sh: u1, u34, u2, assemble pass 1 / zone3 / pass 2 / zone3 / pass 3; snapshot -> runs/prod2  (~6 min)
  flux      flux lemma v2 reassembly (flux/code: forms, coef, chart4_shell, inputs, data, assemble, sanity, texts) on this run's chain   (~4 min)
  proofs    G1_check (level-2 export, on this run's chain + flux), U1_sym (sympy identities), GS_boot (GS strict improvement, fin13v4c)  (~15 s)
  ibp, lemmaO  unchanged parts of kit v2 / v3
The chain parts read the SHIPPED pinned verifier run (step3nl_v4c/runs/P), as the proof tree did; part `verifier` recomputes it and check.py compares.
Each step is a separate process (sys.executable -B), one at a time, explicit environment.  EVERY process goes through kit_guard.py: a write, a working
directory or a child process outside this kit raises an error (see kit_guard.py).  Logs and copies of the outputs: kit_out/."""
import os, sys, subprocess, time, json, platform, shutil, hashlib, ast
KIT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, KIT); import kit_guard                                   # noqa: E402  installs the write guard in this process and routes children through it
OUT = os.path.join(KIT, 'kit_out')
VC = os.path.join(KIT, 'step3nl_v4c'); VER = os.path.join(VC, 'verifier'); CERT = os.path.join(VC, 'v4bcert')
BULK = os.path.join(KIT, 'writeup', 'lemmaU', 'bulk'); NC = os.path.join(BULK, 'fix13_v37chain')
T2 = os.path.join(NC, 't2_A'); BASE = os.path.join(T2, 'base_ibp'); T4 = os.path.join(NC, 't4_cert'); T5 = os.path.join(NC, 't5_zone3'); R3 = os.path.join(NC, 'r3'); N1 = os.path.join(NC, 'n1')
FLUX = os.path.join(NC, 'flux'); PROD = os.path.join(NC, 'runs', 'prod2')
IBP = os.path.join(KIT, 'step3nl', 'agent_mod', 'ibp'); LO = os.path.join(KIT, 'writeup', 'lemmaO'); V38 = os.path.join(KIT, 'step4', 'p2', 'v38')
REF = os.path.join(KIT, 'reference')
VER_ENV = dict(NEW='1', V4B='1', V4C='1', THETA='0.05', ADD='abcd', FIX='1', S0V='0.0153', S1V='0.6', DELTA='0.2', MOD='ibp', MODA='25', ETA='0.01', FLUX0='1.64e-20', FLUX2='8.0e-20')
CHAIN_ENV = dict(V36_VER='fin13v4c', V36_R3='1', V36_R2='1')

def base_env(extra=None):
    e = {k: v for k, v in os.environ.items() if not k.startswith(('PYTHON', 'V36_', 'F2_', 'V4B', 'TAG', 'LY1SQ', 'TABLE_DIR'))}
    e.update(PYTHONHASHSEED='0', PYTHONUTF8='1', PYTHONIOENCODING='utf-8', PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if extra: e.update(extra)
    return e
TIMES = {}; FAILED = []
def step(part, name, script, args, cwd, env=None, outputs=(), dest=None):
    """python -B <script> <args> in cwd (through kit_guard.py); log kit_out/logs/<part>_<name>.log; outputs copied to kit_out/<dest or part>/"""
    dest = dest or part
    assert os.path.realpath(cwd).startswith(os.path.realpath(KIT) + os.sep), ('working directory outside the kit', cwd)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True); os.makedirs(os.path.join(OUT, dest), exist_ok=True)
    log = os.path.join(OUT, 'logs', f'{part}_{name}.log')
    cmd = [sys.executable, '-B', script] + list(args)
    print(f'[{part}] {name}: {os.path.basename(script)} {" ".join(os.path.relpath(a, KIT) if os.path.isabs(a) else a for a in args)}  (cwd {os.path.relpath(cwd, KIT)})', flush=True)
    t0 = time.time()
    with open(log, 'w', encoding='utf-8') as lf:
        p = subprocess.run(cmd, cwd=cwd, env=base_env(env), stdout=lf, stderr=subprocess.STDOUT)
    dt = time.time() - t0; ok = p.returncode == 0
    for o in outputs:
        o, nn = (o if isinstance(o, tuple) else (o, None))
        src = o if os.path.isabs(o) else os.path.join(cwd, o); dst = os.path.join(OUT, dest, nn or os.path.basename(src))
        if os.path.exists(src):
            if os.path.abspath(src) != os.path.abspath(dst): shutil.copy(src, dst)
        elif ok: ok = False; print(f'   MISSING output {o}')
    TIMES[f'{part}/{name}'] = dict(wall_s=round(dt, 1), returncode=p.returncode)
    print(f'   -> {"ok" if ok else "FAILED (rc=%d)" % p.returncode}  {dt:.0f}s   log: {os.path.relpath(log, KIT)}', flush=True)
    if not ok: FAILED.append(f'{part}/{name}')
    old = {}
    try: old = json.load(open(os.path.join(OUT, 'timing.json')))
    except Exception: pass
    old.update(TIMES); json.dump(old, open(os.path.join(OUT, 'timing.json'), 'w'), indent=1)
    return ok
def clear(d, keep=()):
    if os.path.isdir(d):
        for f in os.listdir(d):
            p = os.path.join(d, f)
            if os.path.isfile(p) and f not in keep: os.remove(p)

def part_smoke():
    bad = []
    for l in open(os.path.join(KIT, 'inputs', 'MANIFEST.sha256'), encoding='utf-8'):
        if l.startswith('#') or not l.strip(): continue
        h, r = l.rstrip('\n').split('  ', 1); p = os.path.join(KIT, *r.split('/'))
        if not os.path.exists(p): bad.append(('missing', r))
        elif hashlib.sha256(open(p, 'rb').read()).hexdigest() != h: bad.append(('sha256 differs', r))
    print('   input manifest:', 'ok' if not bad else bad[:10])
    badsyn = []
    for dp, _, fs in os.walk(KIT):
        if os.sep + 'kit_out' in dp: continue
        for f in fs:
            if f.endswith('.py'):
                try: ast.parse(open(os.path.join(dp, f), encoding='utf-8').read())
                except SyntaxError as e: badsyn.append((f, str(e)))
    print('   syntax:', 'ok' if not badsyn else badsyn)
    info = dict(python=sys.version, platform=platform.platform())
    try:
        import flint, numpy, sympy, mpmath
        info.update(python_flint=flint.__version__, numpy=numpy.__version__, sympy=sympy.__version__, mpmath=mpmath.__version__)
    except Exception as e: info['import_error'] = repr(e)
    os.makedirs(OUT, exist_ok=True); json.dump(dict(info, input_manifest_bad=bad, syntax_bad=badsyn), open(os.path.join(OUT, 'env.json'), 'w'), indent=1); print('  ', info)
    if bad or badsyn or 'import_error' in info: FAILED.append('smoke')

def part_verifier():
    d = os.path.join(OUT, 'verifier'); os.makedirs(d, exist_ok=True); clear(d)
    step('verifier', 'tail_cert', os.path.join(VER, 'tail_cert.py'), [], d, VER_ENV, ['tail_cert_log.txt'])
    for a, tg in (('1.5e-4', '0.00015'), ('1.25e-4', '0.000125')):
        step('verifier', f'closure_{a}', os.path.join(VER, 'closure.py'), [a], d, dict(VER_ENV, TAG='fin13v4c', LY1SQ='5e-8'),
             [f'rate_table_fin13v4c_a{tg}.json', f'verify_log_fin13v4c_a{tg}.txt'])
    os.makedirs(os.path.join(d, 'modsups_fin13v4c'), exist_ok=True)
    step('verifier', 'mod_sups', os.path.join(VER, 'mod_sups.py'), ['fin13v4c', os.path.join(d, 'modsups_fin13v4c')], d, VER_ENV,
         [os.path.join(d, 'modsups_fin13v4c', 'final_alpha_points_5e-8.json')])

def part_cert():
    w = os.path.join(OUT, 'cert_work'); os.makedirs(w, exist_ok=True); clear(w)
    for f in os.listdir(os.path.join(CERT, 'out')): shutil.copy(os.path.join(CERT, 'out', f), os.path.join(w, f))
    for f in ('strip_cubic_v4c.json', 'e0_A2_budget_v4c.json'):
        if os.path.exists(os.path.join(w, f)): os.remove(os.path.join(w, f))
    for s in ('strip_cubic_v4c', 'e0_A2_budget_v4c'):
        step('cert', s, os.path.join(CERT, s + '.py'), [], CERT, dict(V4BOUT=w), [os.path.join(w, s + '.json')])

def part_n1():
    os.makedirs(os.path.join(N1, 'out'), exist_ok=True); clear(os.path.join(N1, 'out'))
    step('n1', 'n1_chart', os.path.join(N1, 'n1_chart.py'), [], N1, CHAIN_ENV, [os.path.join(N1, 'out', 'n1_chart_fin13v4c.json')])

BASE_OUTS = [os.path.join(BASE, 'lemmaS_F10p_M0.52.json')] + [os.path.join(BASE, 'rigor', f) for f in ('cap_F10p.json', 'lemmaS_rigcap_M0.52.json', 'weights_exact.json', 'assembly_arb.json')]
def part_lemmaSA():
    if not os.path.exists(os.path.join(N1, 'out', 'n1_chart_fin13v4c.json')): print('   lemmaSA needs part n1 first'); FAILED.append('lemmaSA/needs_n1'); return
    for p in BASE_OUTS:
        if os.path.exists(p): os.remove(p)
    os.makedirs(os.path.join(T2, 'out'), exist_ok=True); clear(os.path.join(T2, 'out'))
    for f in ('r3_forcing.json', 'r3_symcheck.json'):
        if os.path.exists(os.path.join(R3, f)): os.remove(os.path.join(R3, f))
    E = CHAIN_ENV
    ok = step('lemmaSA', 'base_F10p', os.path.join(BASE, 'run_F10p_ibp.py'), [], BASE, E, [BASE_OUTS[0]], 'lemmaSA_base')
    for s, o in (('capture_F2', BASE_OUTS[1:3]), ('weights_exact', [BASE_OUTS[3]]), ('assembly_arb', [BASE_OUTS[4]])):
        ok = ok and step('lemmaSA', 'base_' + s, os.path.join(BASE, 'rigor', s + '.py'), [], BASE, E, o, 'lemmaSA_base')
    ok = ok and step('lemmaSA', 'r3_forcing', os.path.join(R3, 'r3_forcing.py'), [], R3, E, [os.path.join(R3, 'r3_forcing.json')], 'lemmaSA_r3')
    ok = ok and step('lemmaSA', 'r3_symcheck', os.path.join(R3, 'r3_symcheck.py'), [], R3, E, [os.path.join(R3, 'r3_symcheck.json')], 'lemmaSA_r3')
    if not ok: print('   lemmaSA: stopping (base / r3 failed)'); return
    R = lambda n, s, a, outs: step('lemmaSA', n, os.path.join(T2, 'rigor', s), a, T2, E, [os.path.join('out', o) for o in outs])
    st = [('hjw', 'hjw_rep.py', [], ['hjw_rep.json'])]
    for nm, h in (('design', '1e-6,5e-6,1e-5'), ('lo', '8e-7,5e-6,1e-5'), ('hi', '1.2e-6,5e-6,1e-5'), ('probe', '3.5e-7,1.1e-6,2.7e-6')):
        st.append(('ic_' + nm, 'ic_rep.py', ['--h2f', h, '--name', nm], ['ic_%s.json' % nm]))
    st += [('a3_signs', 'a3_signs.py', [], ['a3_signs.json']), ('weights_t2', 'weights_t2.py', [], ['weights_t2.json']),
           ('weights_t2_fb', 'weights_t2.py', ['--seg-box', '20/100', '60/100', '--tag', 'fb'], ['weights_t2_fb.json']),
           ('ramp_rep', 'ramp_rep.py', ['--N', '1200'], ['ramp_rep.json']),
           ('ramp_rep_pm3ibp', 'ramp_rep.py', ['--N', '1200', '--pm3-ibp', '--out', 'out/ramp_rep_pm3ibp.json'], ['ramp_rep_pm3ibp.json']),
           ('outer_L2', 'outer_L2.py', [], ['outer_L2.json'])]
    A = ['--hjw', 'out/hjw_rep.json', '--dprof', os.path.join(KIT, 'step4', 'p2', 'out_v33', 'b0_dprof.json'), '--chart', os.path.join(N1, 'out', 'n1_chart_fin13v4c.json')]
    for nm, n2 in (('design_N2ref', ['0.1', '0.2', '5.0']), ('design_N2ref4', ['0.4', '0.8', '20.0']), ('design_N20', ['0', '0', '0'])):
        st.append(('asm_' + nm, 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2'] + n2 + ['--out', f'out/asm_{nm}.json'], [f'asm_{nm}.json']))
    for s in ('lo', 'hi', 'probe'):
        st.append((f'asm_{s}_N2ref', 'assembly_rep.py', A + ['--ic', f'out/ic_{s}.json', '--h3key', s, '--N2', '0.1', '0.2', '5.0', '--out', f'out/asm_{s}_N2ref.json'], [f'asm_{s}_N2ref.json']))
    st += [('asm_design_pm3ibp', 'assembly_rep.py', A + ['--ic', 'out/ic_design.json', '--h3key', 'design', '--N2', '0.1', '0.2', '5.0', '--ramp', 'out/ramp_rep_pm3ibp.json', '--out', 'out/asm_design_pm3ibp.json'], ['asm_design_pm3ibp.json']),
           ('export_t2', 'export_t2.py', [], ['A_interface.json'])]
    for n, s, a, o in st:
        if not R(n, s, a, o): print('   lemmaSA: stopping (step %s failed)' % n); return

def part_lemmaU():
    o4, o5 = os.path.join(T4, 'out'), os.path.join(T5, 'out'); os.makedirs(o4, exist_ok=True)
    clear(o4); clear(o5, keep=('ramp_U.json',))                             # run_chain_v37.sh: ramp_U.json is kept (an input; the v35 byte copy, review B3)
    AI = os.path.join(T2, 'out', 'A_interface.json')
    if not os.path.exists(AI): print('   lemmaU needs part lemmaSA first'); FAILED.append('lemmaU/needs_lemmaSA'); return
    Z3V35 = os.path.join(BULK, 't5_zone3', 'out', 'zone3_interface.json'); E = CHAIN_ENV; P = 'lemmaU'
    run = lambda n, cwd, s, a, outs: step(P, n, os.path.join(cwd, s), a, cwd, E, outs)
    seq = [('u1_ghat', T4, 'u1_ghat.py', [], [os.path.join(o4, 'u1_ghat.json')]),
           ('u34_primary', T4, 'u34_sweeps.py', ['--case', 'primary'], [os.path.join(o4, 'u34_primary.json')]),
           ('u2_K', T4, 'u2_K.py', ['--fresh'], [os.path.join(o4, 'u2_K.json')])]
    for n, cwd, s, a, outs in seq:
        if not run(n, cwd, s, a, outs): print('   lemmaU: stopping'); return
    for f in os.listdir(o4):
        if f.endswith('.partial.json'): os.remove(os.path.join(o4, f))
    p1, p2, p3 = (os.path.join(o4, f'U_cert_v36_pass{i}.json') for i in (1, 2, 3))
    seq = [('assemble_pass1', T4, 'assemble_U.py', ['--A', AI, '--zone3', Z3V35, '--out', 'U_cert_v36_pass1.json'], [p1]),
           ('for_zone3_pass1', T4, 'make_ucert_for_zone3.py', ['out/U_cert_v36_pass1.json', 'out/U_cert_for_zone3_v36_pass1.json'], [os.path.join(o4, 'U_cert_for_zone3_v36_pass1.json')]),
           ('zone3_pass1', T5, 'zone3.py', ['--ucert', '../t4_cert/out/U_cert_for_zone3_v36_pass1.json', '--ucert_full', p1], [(os.path.join(o5, 'zone3_interface.json'), 'zone3_interface_pass1.json')]),
           ('assemble_pass2', T4, 'assemble_U.py', ['--A', AI, '--zone3', os.path.join(o5, 'zone3_interface_pass1.json'), '--out', 'U_cert_v36_pass2.json'], [p2]),
           ('for_zone3', T4, 'make_ucert_for_zone3.py', ['out/U_cert_v36_pass2.json', 'out/U_cert_for_zone3_v36.json'], [os.path.join(o4, 'U_cert_for_zone3_v36.json')]),
           ('zone3', T5, 'zone3.py', ['--ucert', '../t4_cert/out/U_cert_for_zone3_v36.json', '--ucert_full', p2], [os.path.join(o5, 'zone3_interface.json')]),
           ('assemble_pass3', T4, 'assemble_U.py', ['--A', AI, '--zone3', os.path.join(o5, 'zone3_interface.json'), '--out', 'U_cert_v36_pass3.json'], [p3])]
    for n, cwd, s, a, outs in seq:
        if not run(n, cwd, s, a, outs): print('   lemmaU: stopping'); return
        if n == 'zone3_pass1': shutil.copy(os.path.join(o5, 'zone3_interface.json'), os.path.join(o5, 'zone3_interface_pass1.json'))
    # snapshot = run_chain_v37.sh's copy to runs/<name> (read by flux, F2_V36RUN, and by G1_check)
    if os.path.isdir(PROD): shutil.rmtree(PROD)
    for d, src in (('t2_A', os.path.join(T2, 'out')), ('t4_cert', o4), ('t5_zone3', o5)):
        os.makedirs(os.path.join(PROD, d)); [shutil.copy(os.path.join(src, f), os.path.join(PROD, d, f)) for f in os.listdir(src) if os.path.isfile(os.path.join(src, f))]
    os.makedirs(os.path.join(PROD, 'base_ibp', 'rigor')); shutil.copy(BASE_OUTS[0], os.path.join(PROD, 'base_ibp'))
    for p in BASE_OUTS[1:]: shutil.copy(p, os.path.join(PROD, 'base_ibp', 'rigor'))
    os.makedirs(os.path.join(PROD, 'r3')); [shutil.copy(os.path.join(R3, f), os.path.join(PROD, 'r3')) for f in ('r3_forcing.json', 'r3_symcheck.json')]
    os.makedirs(os.path.join(PROD, 'n1')); shutil.copy(os.path.join(N1, 'out', 'n1_chart_fin13v4c.json'), os.path.join(PROD, 'n1'))
    if step(P, 'v36cfg_record', os.path.join(NC, 'v36cfg.py'), [], NC, E): shutil.copy(os.path.join(OUT, 'logs', f'{P}_v36cfg_record.log'), os.path.join(PROD, 'v36cfg.json'))   # run_chain_v37.sh: python3 v36cfg.py > runs/<name>/v36cfg.json
    print('   snapshot of this run\'s chain ->', os.path.relpath(PROD, KIT))

def part_flux():
    if not os.path.isdir(PROD): print('   flux needs part lemmaU first'); FAILED.append('flux/needs_lemmaU'); return
    for d in ('out', 'texts'): os.makedirs(os.path.join(FLUX, d), exist_ok=True); clear(os.path.join(FLUX, d))
    E = dict(F2_VTAG='fin13v4c', F2_V36RUN=PROD)
    for s in ('flux2_forms', 'flux2_coef', 'chart4_shell', 'inputs_manifest', 'flux2_data', 'flux2_assemble', 'flux2_sanity', 'gen_texts'):
        outs = [os.path.join(FLUX, 'out', f) for f in {'inputs_manifest': ['inputs.json'], 'flux2_assemble': ['flux2_cert.json'], 'flux2_forms': ['flux2_forms.json', 'flux2_forms.pkl']}.get(s, [s + '.json'])]
        if not step('flux', s, os.path.join(FLUX, 'code', s + '.py'), [], os.path.join(FLUX, 'code'), E, outs): print('   flux: stopping'); return
    for f in os.listdir(os.path.join(FLUX, 'texts')): shutil.copy(os.path.join(FLUX, 'texts', f), os.path.join(OUT, 'flux', f))

def part_proofs():
    od = os.path.join(KIT, 'step4', 'p2', 'out_v38', 'iface'); os.makedirs(od, exist_ok=True); clear(od)
    for s, n in (('G1_check_v38', 'G1_check'), ('U1_sym_v38', 'U1_sym'), ('GS_boot_v4c_v38', 'GS_boot_v4c')):
        if s == 'G1_check_v38' and not (os.path.isdir(PROD) and os.path.exists(os.path.join(FLUX, 'out', 'flux2_data.json'))):
            print('   G1_check needs parts lemmaU and flux first'); FAILED.append('proofs/G1_needs_chain'); continue
        step('proofs', n, os.path.join(V38, s + '.py'), [], V38, None, [os.path.join(od, n + '.json')])

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

PARTS = dict(smoke=part_smoke, verifier=part_verifier, cert=part_cert, n1=part_n1, lemmaSA=part_lemmaSA, lemmaU=part_lemmaU, flux=part_flux, proofs=part_proofs,
             ibp=part_ibp, lemmaO=part_lemmaO)
if __name__ == '__main__':
    if '--list' in sys.argv or '-h' in sys.argv or '--help' in sys.argv: print(__doc__); sys.exit(0)
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    unknown = [a for a in args if a not in PARTS]
    if unknown: print('unknown part(s):', unknown, '- choose from', list(PARTS)); sys.exit(2)
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    for p in (args or list(PARTS)): PARTS[p]()
    print(f'\nfinished in {time.time() - t0:.0f}s; failed steps: {FAILED or "none"}')
    print('next:  python check.py   (writes kit_out/check_report.txt; paste its text back)')
    sys.exit(1 if FAILED else 0)
