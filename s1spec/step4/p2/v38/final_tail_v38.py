"""v38 driver: step-4 re-pin on the step-3 objects of Nov 20 2026 (approved by Micah on Nov 20).  Writes p2/out_v38/FINAL_manifest_v38.json.
New file; final_tail_v33.py (pinned by FINAL_manifest_v33.json) and every v33 script are unchanged.  The closed v3.2 manifest p2/out/FINAL_manifest.json
(e01cb4da...) and the v3.3 manifest p2/out_v33/FINAL_manifest_v33.json (8bc47e83..., 21 stages) are verified at the start and at the end and are
never written.
Objects pinned (see p2/v38/V38_LOG.md):
  verifier fin13v4c, production run step3nl_v4c/runs/P only (fix15_review REVIEW_v4c_interface.md 1.4: runs/P_A and runs/reg_fin13v4b each carry one
    unverifiable pin and are NOT pinned), with its cert inputs v4bcert/out (strip_cubic_v4c.json = D-3 / H8, e0_A2_budget_v4c.json = A-1);
  (A) chain v37 writeup/lemmaU/bulk/fix13_v37chain (t2_A, t4_cert, t5_zone3, n1; in-place outputs == runs/prod2, checked) and flux v37 (flux/out, texts);
  interface certificates C1_ke, GS_tau1 (fix14_interface) and GS_boot on v4c (fix15_review/c); proof checks G1_check, U1_sym (fix10_proofs).
verify_v32 consults p2/v33/EXCEPTION_v32_*.json exactly as p2/v33/verify_v32_with_exceptions.py does: a mismatch of the closed manifest is accepted only
for a path listed in an exception record whose current sha256 equals the record's current_sha256 and whose pinned sha256 equals the record's
pinned_sha256 (the Nov 15 B8 incident: s1_tight/consistency/check_outA5.json); every other mismatch stops the driver.
Stages (indices as in v33 for 0-20):
  0-7, 10-17  CARRIED from FINAL_manifest_v33.json by hash (record copied; every input/output sha256 re-checked against the disk; not rerun);
  8  ratio_v38      H7 time ratio on the fin13v4c mod_sups (abar);              9  notes_v38   covering notes on the v37 zone3_interface;
  18 h4far_v38      far-family H4 (zone3_interface v37);                         19 n1_chart_v38 chart terms on the fin13v4c rate table;
  20 f12_sliver_v38 F1/F2 standalone check on fin13v4c (+ strip cubic, all rate-table pins) and chain v37 (writes only out_v38/f12/);
  21 flux_v38       flux v37 reproduction (F2_OUT/F2_TX scratch under out_v38/flux) and pin;
  22 C1_ke_v38, 23 GS_tau1_v38, 24 GS_boot_v4c_v38, 25 G1_check_v38 (chain v37 + flux v37), 26 U1_sym_v38: reruns into out_v38/iface;
  27 iface_compare_v38: reruns vs reference outputs;
  28 verify_v4c_runs_with_exceptions: runs/P pins all match; runs/P_A, runs/reg_fin13v4b (NOT pinned) accepted only through p2/v38/EXCEPTION_v4c_PA_reg_e2v4b.json.
Records: each stage record pins cwd, cmd, env, env_child (complete child environment, built from scratch), env_read, inputs (declared + every absolute
path in the args block of the stage's JSON outputs, with the sha256 the stage recorded asserted equal to the disk), outputs (sha256), and the stage
log separately (log: logs of some stages print wall times, so a log hash is a timing field like 'seconds').
usage: nice -n 19 python3 -I -B final_tail_v38.py [--from STAGE] [--until STAGE] [--final PATH]
  --from: rerun from that stage (index or name prefix); earlier stage records come from out_v38/FINAL_manifest_v38.partial.json if it exists, else from
  the final manifest, after checking that their definitions are unchanged and every recorded sha256 matches the file now.  --until: stop after that stage
  (partial file kept).  Progress is written atomically to the partial file after each stage; the final manifest is written atomically at the end."""
import os, sys, json, hashlib, subprocess, time, glob
S4 = '/home/claude/s1spec/step4'; P2 = f'{S4}/p2'; O = f'{P2}/out'; OV33 = f'{P2}/out_v33'; OV = f'{P2}/out_v38'; V33 = f'{P2}/v33'; V = f'{P2}/v38'
T8 = f'{P2}/t8se'; T6C = f'{P2}/t6c'
V32MAN = f'{O}/FINAL_manifest.json'; V32MAN_SHA = 'e01cb4da71ba825318b5395f053595d6e980e84c1ac286b9d72ff1e4376a8c94'
V33MAN = f'{OV33}/FINAL_manifest_v33.json'; V33MAN_SHA = '8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec'
FINAL = sys.argv[sys.argv.index('--final') + 1] if '--final' in sys.argv else f'{OV}/FINAL_manifest_v38.json'
PARTIAL = FINAL[:-5] + '.partial.json'
for d in ('logs', 'f12', 'flux', 'iface'): os.makedirs(f'{OV}/{d}', exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
ab = lambda p: p if p.startswith('/') else os.path.join(S4, p)
rel = lambda p: os.path.relpath(p, S4) if p.startswith(S4 + '/') else p
def write_atomic(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1); fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, path)

def exceptions():
    exc = {}
    for p in sorted(glob.glob(f'{V33}/EXCEPTION_v32_*.json')):
        e = json.load(open(p)); exc[e['path']] = dict(e, record=rel(p), record_sha256=sha(p))
    return exc
def verify_v32():
    """every sha256 pinned by FINAL_manifest.json (stage paths at their last write, upstream, scripts) against the files now; a mismatch is accepted
    only through an EXCEPTION_v32_*.json record (same rule as verify_v32_with_exceptions.py)"""
    m = json.load(open(V32MAN)); last = {}
    for P in m['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]
    bad += [u['path'] for u in m['upstream'] if sha(ab(u['path'])) != u['sha256']]
    bad += [q for q, h in m.get('script_sha256', {}).items() if sha(ab(q)) != h]
    exc = exceptions()
    acc = [q for q in bad if q in exc and sha(exc[q]['abs_path']) == exc[q]['current_sha256'] and last.get(q) == exc[q]['pinned_sha256']]
    other = [q for q in bad if q not in acc]
    return dict(manifest=rel(V32MAN), manifest_sha256=sha(V32MAN), manifest_sha256_expected=V32MAN_SHA, stages=len(m['stages']), paths=len(last),
                upstream=len(m['upstream']), scripts=len(m.get('script_sha256', {})), mismatches=bad,
                accepted_by_exception={q: dict(record=exc[q]['record'], record_sha256=exc[q]['record_sha256'], pinned_sha256=exc[q]['pinned_sha256'],
                                               current_sha256=exc[q]['current_sha256']) for q in acc},
                other_mismatches=other, ok=not other and sha(V32MAN) == V32MAN_SHA), last, {u['path']: u['sha256'] for u in m['upstream']}
# Known, documented difference to FINAL_manifest_v33.json (found by this driver on Nov 20; not an exception record, the path is an input of v33 stage 20
# f12_sliver only, which v38 replaces): fix2/F34_REPORT.md line 6 -- after the v33 manifest was written, make_constants_v34.py regenerated
# constants_needed_v34.json so that its own S4 pin names the final v33 manifest (8bc47e83...) instead of a0d348ec... (a circular pin: the file cannot
# be in the manifest it pins).  The pinned bytes are kept as constants_needed_v34_pre_f34.json; the two files differ in that one sha256 line only.
V33_KNOWN = {'/home/claude/s1spec/writeup/lemmaU/bulk/fix1/lemmaU_rev/constants_needed_v34.json': dict(
    pinned_sha256='0f3c745546a4dbf9ec2dc12e4e3b5571da001ce121039e719e1e2be23c09d37c', current_sha256='e70e70f0d0f31be7e0c011f1cd80dbd859d54657b25aa07ebd87b8837f3f139f',
    pinned_bytes_kept_as='/home/claude/s1spec/writeup/lemmaU/bulk/fix1/lemmaU_rev/constants_needed_v34_pre_f34.json', source='writeup/lemmaU/bulk/fix2/F34_REPORT.md line 6',
    only_line_changed='"sha256": "a0d348ec..." -> "sha256": "8bc47e83..." (its pin of the step-4 v33 manifest)', used_by='v33 stage 20 f12_sliver only (replaced in v38)')}
def _known_ok(q, h):
    k = V33_KNOWN.get(q)
    if not k or h != k['pinned_sha256'] or sha(q) != k['current_sha256'] or sha(k['pinned_bytes_kept_as']) != k['pinned_sha256']: return False
    a, b = open(k['pinned_bytes_kept_as']).read().split('\n'), open(q).read().split('\n')
    d = [(x, y) for x, y in zip(a, b) if x != y]
    return len(a) == len(b) and len(d) == 1 and d[0][0].strip().startswith('"sha256": "a0d348ec') and d[0][1].strip() == '"sha256": "' + V33MAN_SHA + '"'
def verify_v33():
    """every sha256 pinned by FINAL_manifest_v33.json (stage paths at their last write, reused v3.2 inputs, scripts) against the files now;
    the one known, documented difference V33_KNOWN is accepted only if it is exactly as recorded there"""
    m = json.load(open(V33MAN)); last = {}
    for P in m['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]
    bad += [q for q, h in m['reused_v32'] if sha(ab(q)) != h]
    bad += [q for q, h in m['script_sha256'].items() if sha(ab(q)) != h]
    known = [q for q in bad if _known_ok(ab(q), last.get(q))]
    other = [q for q in bad if q not in known]
    return dict(manifest=rel(V33MAN), manifest_sha256=sha(V33MAN), manifest_sha256_expected=V33MAN_SHA, stages=len(m['stages']), paths=len(last),
                reused_v32=len(m['reused_v32']), scripts=len(m['script_sha256']), mismatches=bad, known_documented={q: V33_KNOWN[ab(q)] for q in known},
                other_mismatches=other, ok=not other and sha(V33MAN) == V33MAN_SHA), m

BU = '/home/claude/s1spec/writeup/lemmaU/bulk'; NC = f'{BU}/fix13_v37chain'; VC = '/home/claude/s1spec/step3nl_v4c'
RT = f'{VC}/runs/P/rate_table_fin13v4c_a0.00015.json'; RT125 = f'{VC}/runs/P/rate_table_fin13v4c_a0.000125.json'
MS = f'{VC}/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json'; SCUB = f'{VC}/v4bcert/out/strip_cubic_v4c.json'; A2B = f'{VC}/v4bcert/out/e0_A2_budget_v4c.json'
LOV32 = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json'; Z3 = f'{NC}/t5_zone3/out/zone3_interface.json'
AMP = f'{T8}/constants/amp_A1_w180.json'; T9K = f'{O}/t9_cone_kch.json'; N1C = f'{NC}/n1/out/n1_chart_fin13v4c.json'
JETS = lambda *sl: [f'{T6C}/jets/jets_{s:03d}.npz' for s in sl]
T2, T4, T5 = f'{NC}/t2_A', f'{NC}/t4_cert', f'{NC}/t5_zone3'
_ASM = ['design_N2ref', 'design_N2ref4', 'design_N20', 'lo_N2ref', 'hi_N2ref', 'probe_N2ref', 'design_pm3ibp']
OF = f'{OV}/f12'; OI = f'{OV}/iface'; OFX = f'{OV}/flux'
F12_IN = ([f'{V}/f12_sliver_v38.py', f'{NC}/n1/n1_chart.py', N1C, f'{OV}/n1_chart.json', RT, RT125, MS, SCUB, A2B, f'{VC}/agent_sob/sob_constants_P2.json',
           f'{VC}/verifier/constants.py', f'{P2}/t9h/s1chart.py', f'{NC}/v36cfg.py', f'{NC}/r3/r3_forcing.py',
           f'{T2}/rigor/assembly_rep.py', f'{T2}/rigor/arlib.py', f'{T2}/rigor/export_t2.py', f'{T2}/base_ibp/rigor/assembly_arb.py',
           f'{T4}/assemble_U.py', f'{T4}/t4lib.py', f'{T4}/make_ucert_for_zone3.py', f'{T5}/zone3.py', f'{T5}/z3lib.py']
          + [f'{T2}/out/asm_{n}.json' for n in _ASM] + [f'{T2}/out/A_interface.json', f'{T4}/out/U_cert_v36_pass2.json', f'{T4}/out/U_cert_v36_pass3.json',
           f'{T4}/out/U_cert_for_zone3_v36.json', Z3, f'{BU}/fix1/lemmaU_rev/constants_needed_v34.json', f'{BU}/fix12_texts/tools/gen_numbers_A_v37.py'])
F12_OUT = ([f'{OF}/f12_sliver.json', f'{OF}/repro/n1_chart.json', f'{OF}/new/n1_chart.json'] + [f'{OF}/{d}/asm_{n}.json' for n in _ASM for d in ('repro', 'new')]
           + [f'{OF}/repro/A_interface.json', f'{OF}/new/A_interface.json', f'{OF}/repro/U_cert.json', f'{OF}/new/U_cert.json', f'{OF}/new/U_cert_pass2.json',
              f'{OF}/repro/U_cert_for_zone3.json', f'{OF}/new/U_cert_for_zone3.json', f'{OF}/repro/zone3_interface.json', f'{OF}/new/zone3_interface.json'])
FLUXCODE = f'{NC}/flux/code'
FLUX_STAGES = ['flux2_forms', 'flux2_coef', 'chart4_shell', 'inputs_manifest', 'flux2_data', 'flux2_assemble', 'flux2_sanity', 'gen_texts']
FLUX_PIN = ([f'{NC}/flux/out/{f}.json' for f in ('chart4_shell', 'flux2_cert', 'flux2_coef', 'flux2_data', 'flux2_forms', 'flux2_sanity', 'gen_texts', 'inputs')]
            + [f'{NC}/flux/out/flux2_forms.pkl'] + [f'{NC}/flux/texts/{f}' for f in sorted(os.listdir(f'{NC}/flux/texts'))])
FLUX_OUT = ([f'{OFX}/flux_v38.json'] + [f'{OFX}/out/{f}.json' for f in ('chart4_shell', 'flux2_cert', 'flux2_coef', 'flux2_data', 'flux2_forms', 'flux2_sanity', 'gen_texts', 'inputs')]
            + [f'{OFX}/out/flux2_forms.pkl'] + [f'{OFX}/texts/{os.path.basename(f)}' for f in FLUX_PIN if '/texts/' in f])
IREF = dict(C1_ke=[f'{BU}/fix15_review/c/C1_ke.json', f'{BU}/fix14_interface/C1_ke.json'], GS_tau1=[f'{BU}/fix15_review/c/GS_tau1.json', f'{BU}/fix14_interface/GS_tau1.json'],
            GS_boot_v4c=[f'{BU}/fix15_review/c/GS_boot_v4c.json'], G1_check=[f'{BU}/fix10_proofs/G1_check.json'], U1_sym=[f'{BU}/fix10_proofs/U1_sym.json'])
ISRC = dict(C1_ke=f'{BU}/fix15_review/c/C1_ke.py', GS_tau1=f'{BU}/fix15_review/c/GS_tau1.py', GS_boot_v4c=f'{BU}/fix15_review/c/GS_boot_v4c.py',
            G1_check=f'{BU}/fix10_proofs/G1_check.py', U1_sym=f'{BU}/fix10_proofs/U1_sym.py')

def st(name, cwd, cmd, env, inp, outp, log): return dict(stage=name, cwd=cwd, cmd=cmd, env=env, inputs=inp, outputs=outp, log=f'{OV}/logs/{log}')
def carry(i): return dict(carry=i)
BASE_ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}
READ_VARS = sorted(set('ALIN BASE_T1 BLOOP BSUB CAB CERT_TAG COHERENT CUBFINE EPSI EXE FAR_M FORM_MULT H1JOBS H1NSUB H1ONLY H1REF H1TAG H4OUT HSTEP '
                       'LO_NS MODA MODE2 N3JSON NLFAC OMEGA_Q PDOT POS_C0_MILLI QSHIFT REORG S1PROF S1TT SLAB0 STIPPOINT SUBW T6DIR T6RES T8CONSTS '
                       'T9IFT T9JSON T9OUT T9PIECE T9RATIO T9TC T9WIN TOPLOSS TREF V33_VALIDATE WINOUT YBAR YOUT ZSHIFT '
                       'V36_VER V36_R3 V36_R2 F2_OUT F2_TX F2_VTAG F2_V36RUN F2_NCELL F2_DRYRUN '
                       'PYTHONHASHSEED PYTHONDONTWRITEBYTECODE PYTHONPATH PYTHONHOME OMP_NUM_THREADS LC_ALL LANG PATH'.split()))
def child_env(S): return dict(BASE_ENV, **S['env'])
def env_read(S): E = child_env(S); return {k: E.get(k) for k in READ_VARS}
CHAIN_ENV = dict(V36_VER='fin13v4c', V36_R3='1', V36_R2='1')
STAGES = [carry(i) for i in range(8)] + [
    st('fixed-L cone: H7 time ratio on fin13v4c (ratio_v38.py)', V, ['ratio_v38.py'], {}, [f'{V}/ratio_v38.py', MS, LOV32], [f'{OV}/ratio_v38.json'], 'ratio_v38.log'),
    st('covering notes on the v37 zone3_interface: far family subset, t9_h6_v32 superseded (notes_v38.py)', V, ['notes_v38.py'], {},
       [f'{V}/notes_v38.py', Z3, f'{O}/t9_lemmaO_far.json', f'{O}/t9_far_v32.json', f'{T8}/constants/h4loc_v32.json'], [f'{OV}/notes_v38.json'], 'notes_v38.log'),
    ] + [carry(i) for i in range(10, 18)] + [
    st('M9: far-family H4, every centre at radius >= 490 to the far end, zone3 v37 (h4far_v38.py)', T8, ['../v38/h4far_v38.py', 'constants/amp_A1_w180.json'],
       dict(T6DIR='../t6c', T9JSON='../out/t9_cone_kch.json'),
       [f'{V}/h4far_v38.py', f'{T8}/h4loc.py', f'{P2}/t9chart.py', AMP, T9K, LOV32, f'{T8}/constants/constants_b1.json', f'{T8}/constants/h4loc_v32.json', Z3] + JETS(180),
       [f'{OV}/h4far_v38.json'], 'h4far_v38.log'),
    st('N1: chart terms, S1 chart vs ray graph, on (A)\'s ramp and core, verifier fin13v4c (n1_chart_v38.py)', V, ['n1_chart_v38.py'], {},
       [f'{V}/n1_chart_v38.py', RT, '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
        f'{P2}/t9h/s1chart.py', '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/step3/hs_uforms.pkl', '/home/claude/s1spec/lin/lintail2.py', '/home/claude/s1_cap/tail.py'],
       [f'{OV}/n1_chart.json'], 'n1_chart_v38.log'),
    st('f12_sliver: F1 sliver G_i, F2 sliver sup|w\'\'|, on fin13v4c (strip cubic, rate-table pins) and chain v37 (f12_sliver_v38.py)', V, ['f12_sliver_v38.py'], CHAIN_ENV,
       F12_IN, F12_OUT, 'f12_sliver_v38.log'),
    st('flux v37: reproduction (F2_OUT/F2_TX scratch) and pin of fix13_v37chain/flux/out, texts (flux_v38.py)', V, ['flux_v38.py'], {},
       [f'{V}/flux_v38.py', f'{FLUXCODE}/f2common.py'] + [f'{FLUXCODE}/{s}.py' for s in FLUX_STAGES] + FLUX_PIN, FLUX_OUT, 'flux_v38.log'),
    st('C-1: sup|k~ e| < 1/2 on step-4 P2 boxes (C1_ke_v38.py)', V, ['C1_ke_v38.py'], {}, [f'{V}/C1_ke_v38.py', ISRC['C1_ke']] + IREF['C1_ke'], [f'{OI}/C1_ke.json'], 'C1_ke_v38.log'),
    st('(GS) at tau_1: tip arc a normal graph over S1 on |x| <= 60.25 (GS_tau1_v38.py)', V, ['GS_tau1_v38.py'], {}, [f'{V}/GS_tau1_v38.py', ISRC['GS_tau1']] + IREF['GS_tau1'],
       [f'{OI}/GS_tau1.json'], 'GS_tau1_v38.log'),
    st('(GS) strict improvement under the bootstrap, verifier fin13v4c (GS_boot_v4c_v38.py)', V, ['GS_boot_v4c_v38.py'], {}, [f'{V}/GS_boot_v4c_v38.py', ISRC['GS_boot_v4c'], RT]
       + IREF['GS_boot_v4c'], [f'{OI}/GS_boot_v4c.json'], 'GS_boot_v4c_v38.log'),
    st('G1: level-2 export check on chain v37 + flux v37 (G1_check_v38.py)', V, ['G1_check_v38.py'], {}, [f'{V}/G1_check_v38.py', ISRC['G1_check']] + IREF['G1_check'],
       [f'{OI}/G1_check.json'], 'G1_check_v38.log'),
    st('U-1: exact sympy identities (U1_sym_v38.py)', V, ['U1_sym_v38.py'], {}, [f'{V}/U1_sym_v38.py', ISRC['U1_sym']] + IREF['U1_sym'], [f'{OI}/U1_sym.json'], 'U1_sym_v38.log'),
    st('interface/proof certificates: reruns vs references (iface_compare_v38.py)', V, ['iface_compare_v38.py'], {},
       [f'{V}/iface_compare_v38.py'] + [f'{OI}/{n}.json' for n in IREF] + [r for v in IREF.values() for r in v], [f'{OI}/iface_compare_v38.json'], 'iface_compare_v38.log'),
    st('verifier fin13v4c runs: runs/P pins; runs/P_A, runs/reg_fin13v4b with EXCEPTION_v4c_PA_reg_e2v4b.json, not pinned (verify_v4c_runs_with_exceptions.py)', V,
       ['verify_v4c_runs_with_exceptions.py'], {}, [f'{V}/verify_v4c_runs_with_exceptions.py', f'{V}/EXCEPTION_v4c_PA_reg_e2v4b.json', RT, RT125, MS], [f'{OV}/v4c_runs_pins.json'], 'v4c_runs_pins.log'),
]
READONLY_STEP4 = [f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py', f'{P2}/t9h/s1chart.py', f'{P2}/t9chart.py', f'{P2}/curvjets.py',
                  f'{S4}/t6/reader.py', f'{P2}/t9.py', f'{T8}/h4loc.py', f'{P2}/t9_lemmaO.py']
def stage_name(S, M33): return S['stage'] if 'stage' in S else M33['stages'][S['carry']]['stage']
def args_inputs(p):
    """absolute input paths with the sha256 a stage recorded in the args block of its JSON output p"""
    try: x = json.load(open(p))
    except Exception: return {}
    a = x.get('args') if isinstance(x, dict) else None
    if not isinstance(a, dict): return {}
    inp = a.get('inputs', {}); out = {}
    if isinstance(inp, list): inp = {d['path']: d['sha256'] for d in inp if isinstance(d, dict) and 'path' in d}
    for k, v in (inp.items() if isinstance(inp, dict) else []):
        if isinstance(v, dict) and 'path' in v: k, v = v['path'], v.get('sha256')
        if isinstance(k, str) and k.startswith('/') and isinstance(v, str) and len(v) == 64: out[k] = v
    return out

if __name__ == '__main__':
    v32_start, V32LAST, V32UP = verify_v32(); assert v32_start['ok'], ('v3.2 manifest does not verify (exceptions consulted)', v32_start['other_mismatches'])
    v33_start, M33 = verify_v33(); assert v33_start['ok'], ('v3.3 manifest does not verify', v33_start['other_mismatches'][:20])
    print(f"v3.2 manifest verified: {v32_start['stages']} stages, {v32_start['paths']} paths, {v32_start['upstream']} upstream, {v32_start['scripts']} scripts; "
          f"accepted by exception: {sorted(v32_start['accepted_by_exception'])}", flush=True)
    print(f"v3.3 manifest verified: {v33_start['stages']} stages, {v33_start['paths']} paths, {v33_start['reused_v32']} reused v3.2, {v33_start['scripts']} scripts; "
          f"known documented difference: {sorted(v33_start['known_documented'])}", flush=True)
    man = dict(note=__doc__, v32_manifest=dict(path=rel(V32MAN), sha256=v32_start['manifest_sha256']), v33_manifest=dict(path=rel(V33MAN), sha256=v33_start['manifest_sha256']),
               stages=[], checks=dict(v32_manifest_at_start=v32_start, v33_manifest_at_start=v33_start),
               environment=dict(base_env=BASE_ENV, read_vars=READ_VARS, interpreter=sys.executable, flags=['-I', '-B'],
                                note='children get child_env = base_env + stage env only; -I implies -E: PYTHON* variables are ignored by the child interpreter '
                                     '(PYTHONHASHSEED=0 is set and recorded; byte identity is checked by two clean passes)',
                                parent_values_not_passed={k: os.environ[k] for k in READ_VARS if k in os.environ and k not in ('PATH', 'LANG', 'LC_ALL')}))
    i0, i1 = 0, len(STAGES) - 1
    def stage_index(key, start=0):
        return int(key) if key.isdigit() else next(i for i, S in enumerate(STAGES) if i >= start and stage_name(S, M33).startswith(key))
    def same_def(S, P):
        if 'carry' in S: return P.get('carried_from', {}).get('index') == S['carry'] and P['stage'] == M33['stages'][S['carry']]['stage']
        return (P['stage'] == S['stage'] and P['cmd'] == ['python3', '-I', '-B'] + S['cmd'] and P['env'] == S['env'] and P['cwd'] == rel(S['cwd'])
                and P.get('env_child') == child_env(S) and P.get('env_read') == env_read(S)
                and [f['path'] for f in P['inputs'] if not f.get('from_args')] == [rel(q) for q in S['inputs']] and [f['path'] for f in P['outputs']] == [rel(q) for q in S['outputs']])
    if '--from' in sys.argv:
        i0 = stage_index(sys.argv[sys.argv.index('--from') + 1]); cont = os.path.exists(PARTIAL)
        prev = json.load(open(PARTIAL if cont else FINAL)); last = {}
        assert i0 <= len(prev['stages'])
        for S, P in zip(STAGES[:i0], prev['stages'][:i0]):
            assert same_def(S, P), ('stage definition changed', stage_name(S, M33))
            for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
            for f in P['outputs']: last[f['path']] = f['sha256']
        bad = [q for q, h in last.items() if sha(ab(q)) != h]; assert not bad, bad
        man['stages'] = prev['stages'][:i0]
        man['checks'].update({k: v for k, v in prev['checks'].items() if k in ('segments',) and cont})
        man['checks'].setdefault('segments', []).append([i0, None])
        print(f'from stage {i0}: stages 0-{i0 - 1} reused ({len(last)} paths verified)', flush=True)
    if '--until' in sys.argv: i1 = stage_index(sys.argv[sys.argv.index('--until') + 1], i0)
    for k, S in enumerate(STAGES[i0:i1 + 1], start=i0):
        if 'carry' in S:                                                # carried by hash from FINAL_manifest_v33.json
            P = json.loads(json.dumps(M33['stages'][S['carry']]))
            bad = [f['path'] for f in P['inputs'] + P['outputs'] if sha(ab(f['path'])) != f['sha256']]
            assert not bad, ('carried stage: hash mismatch', P['stage'], bad)
            P['carried_from'] = dict(manifest=rel(V33MAN), manifest_sha256=V33MAN_SHA, index=S['carry'], rerun=False,
                                     verified=f"{len(P['inputs'])} inputs + {len(P['outputs'])} outputs re-hashed equal")
            man['stages'].append(P); write_atomic(PARTIAL, man)
            print(f"[{k}] carried from v33 stage {S['carry']}: {P['stage']} ({len(P['inputs'])} in, {len(P['outputs'])} out, hashes equal)", flush=True)
            continue
        ins = []
        for p in S['inputs']:
            h = sha(p); assert h != 'MISSING', p
            rp = rel(p); pinned = V32LAST.get(rp) or V32UP.get(rp)
            if pinned is not None: assert pinned == h, ('v3.2 pinned input changed', rp)
            ins.append(dict(path=rp, sha256=h, v32_pinned=pinned is not None))
        for p in S['outputs'] + [S['log']]:
            assert p.startswith(OV + '/'), ('output outside out_v38', p)
            if os.path.exists(p): os.remove(p)                                 # no stale output can be pinned
        env = child_env(S); t0 = time.time()
        with open(S['log'], 'w') as lf:
            r = subprocess.run([sys.executable, '-I', '-B'] + S['cmd'], cwd=S['cwd'], env=env, stdout=lf, stderr=subprocess.STDOUT)
        print(f"[{k}] {S['stage']}: rc {r.returncode} ({time.time() - t0:.0f}s)", flush=True)
        assert r.returncode == 0, (S['stage'], open(S['log']).read()[-3000:])
        outs = [dict(path=rel(p), sha256=sha(p)) for p in S['outputs']]
        assert all(o['sha256'] != 'MISSING' for o in outs), [o['path'] for o in outs if o['sha256'] == 'MISSING']
        declared = {i['path'] for i in ins}; outset = {ab(o['path']) for o in outs}; extra = {}
        for p in S['outputs']:
            for q, h in args_inputs(p).items():
                if q in outset or q.startswith(OV + '/f12/') or q.startswith(OV + '/flux/') or rel(q) in declared or q in extra: continue
                hq = sha(q); assert hq == h, ('stage recorded an input sha256 that differs from the disk now', S['stage'], q, h, hq)
                extra[q] = hq
        for q in sorted(extra):
            rp = rel(q); pinned = V32LAST.get(rp) or V32UP.get(rp)
            if pinned is not None: assert pinned == extra[q], ('v3.2 pinned input changed', rp)
            ins.append(dict(path=rp, sha256=extra[q], v32_pinned=pinned is not None, from_args=True))
        man['stages'].append(dict(stage=S['stage'], cwd=rel(S['cwd']), cmd=['python3', '-I', '-B'] + S['cmd'], env=S['env'], env_child=env, env_read=env_read(S),
                                  inputs=ins, outputs=outs, log=dict(path=rel(S['log']), sha256=sha(S['log']), note='timing field (some logs print wall times)'),
                                  seconds=round(time.time() - t0, 1)))
        write_atomic(PARTIAL, man)
    if man['checks'].get('segments'): man['checks']['segments'][-1][1] = i1
    if i1 < len(STAGES) - 1:
        man['checks']['in_progress'] = dict(next=i1 + 1, next_stage=stage_name(STAGES[i1 + 1], M33)); write_atomic(PARTIAL, man)
        print(f'stopped after stage {i1}; continue with --from {i1 + 1}'); sys.exit(0)
    assert [P['stage'] for P in man['stages']] == [stage_name(S, M33) for S in STAGES], 'incomplete stage list'
    seen = {}
    for P in man['stages']:
        for i in P['inputs']: seen.setdefault(i['path'], set()).add(i['sha256'])
    man['checks']['inputs_with_two_versions'] = sorted(p for p, v in seen.items() if len(v) > 1); assert not man['checks']['inputs_with_two_versions']
    man['reused_v32'] = sorted({i['path']: i['sha256'] for P in man['stages'] for i in P['inputs'] if i['v32_pinned']}.items())
    man['script_sha256'] = {rel(f): sha(f) for f in sorted(glob.glob(f'{V}/*.py')) + sorted(glob.glob(f'{V33}/*.py')) + READONLY_STEP4}
    J = lambda p: json.load(open(p))
    res = {k: v for k, v in M33['results'].items() if k not in ('H7_fixedL', 'M9_H4_far', 'N1_chart', 'F12_sliver')}
    res = {k: dict(v, carried_from_v33=True) for k, v in res.items()}
    ra = J(f'{OV}/ratio_v38.json'); res['H7_fixedL'] = dict({k: ra[k] for k in ('abar', 'ratio_upper', 'needed', 'ratio_le_needed', 't_star_lower', 'margin_tstar_over_ratio')}, ok=ra['ratio_le_needed'])
    nt = J(f'{OV}/notes_v38.json'); res['notes'] = dict(far_family_subset=nt['far_family']['subset'], H6_task5_min_margin=nt['h6']['task5_min_margin_at_Bhat_ref'], ok=nt['far_family']['subset'])
    h4f = J(f'{OV}/h4far_v38.json'); res['M9_H4_far'] = dict(min_margin=h4f['min_margin'], worst=h4f['worst'], n_rows=h4f['n_rows'], XSTART=h4f['XSTART'], ok=h4f['all_ok'])
    n1 = J(f'{OV}/n1_chart.json'); res['N1_chart'] = dict(ramp={k: n1['ramp'][k] for k in ('sup_add', 'X_W', 'beta22', 'beta32', 'beta33', 'Dbeta22', 'eps', 'e', 'rho_meas')},
                                                          core={k: n1['core'][k] for k in ('sup_add', 'E', 'beta22', 'beta3', 'e', 'rho_meas')}, dB0_slivers=n1['dB0_slivers'],
                                                          ok=n1['ramp']['e'] < 2 ** -19 and n1['ramp']['sup_add'][0] <= 2e-7)
    f12 = J(f'{OF}/f12_sliver.json')
    res['F12_sliver'] = dict(G_old=f12['F1']['G_old'], G_new=f12['F1']['G_new'], sI2=f12['F2']['rows']['1.5e-4']['sI2'], sI2_prime=f12['F2']['sI2_prime_used'],
                             M1_prime=f12['F2']['M1_prime'], M1_read_by_A=f12['M1_read_by_A'],
                             named={k: f12['named'][k] for k in ('ramp_M1', 'K6', 'H6', 'S0v', 'flux_beta2_4N2ref', 'EA_ref', 'EA_ref4', 'Ehat_U')},
                             n_checks=len(f12['checks']), n_pass=sum(c['ok'] for c in f12['checks']), macro_changes=f12['macro_changes'],
                             rounding_flips=f12['rounding_flips'], ok=f12['all_ok'])
    fx = J(f'{OFX}/flux_v38.json'); res['flux_v37'] = dict(combined_margin_lower=fx['cert']['combined_margin_lower'], verdict=fx['cert']['verdict'],
                                                           conditions={k: v.get('margin') for k, v in fx['cert']['conditions'].items()}, B_v37_1=fx['B_v37_1'],
                                                           n_checks=len(fx['checks']), ok=fx['all_ok'])
    vr = J(f'{OV}/v4c_runs_pins.json'); res['v4c_runs_pins'] = dict(n_checks=len(vr['checks']), exception=vr['exception'], ok=vr['all_ok'],
                                                                  accepted={k: sum(len(r['accepted_by_exception']) for r in v.values()) for k, v in vr['runs'].items()})
    ic = J(f'{OI}/iface_compare_v38.json'); res['interface_proofs'] = dict(n_checks=len(ic['checks']), ok=ic['all_ok'])
    for n in ('C1_ke', 'GS_tau1', 'GS_boot_v4c', 'G1_check', 'U1_sym'):
        x = J(f'{OI}/{n}.json'); res['interface_proofs'][n] = {k: x[k] for k in ('ok', 'all_pass', 'sup_ke', 'sup_ku', 'sup_uprime', 'XS') if k in x}
    man['results'] = res; man['checks']['all_items_ok'] = all(v.get('ok', False) for v in res.values()); assert man['checks']['all_items_ok'], {k: v.get('ok') for k, v in res.items()}
    last = {}
    for P in man['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]; assert not bad, bad
    v32_end, _, _ = verify_v32(); man['checks']['v32_manifest_at_end'] = v32_end; assert v32_end['ok'], v32_end['other_mismatches']
    v33_end, _ = verify_v33(); man['checks']['v33_manifest_at_end'] = v33_end; assert v33_end['ok'], v33_end['other_mismatches'][:20]
    man['checks']['hashes_verified_at_end'] = dict(paths=len(last), scripts=len(man['script_sha256']), reused_v32=len(man['reused_v32']),
                                                   carried_stages=sum(1 for P in man['stages'] if 'carried_from' in P), rerun_or_new_stages=sum(1 for P in man['stages'] if 'carried_from' not in P))
    write_atomic(FINAL, man)
    if os.path.exists(PARTIAL): os.remove(PARTIAL)
    print('results:', json.dumps({k: v for k, v in res.items() if not v.get('carried_from_v33')})[:4000])
    print('checks:', json.dumps({k: v for k, v in man['checks'].items() if not k.startswith('v3')}))
    print(f"v3.2 manifest at end: ok = {v32_end['ok']} (accepted by exception: {sorted(v32_end['accepted_by_exception'])}); v3.3 manifest at end: ok = {v33_end['ok']}")
    print(f'manifest: {FINAL} ({len(man["stages"])} stages)')
