"""v3.3 task 6 driver: the NEW step-4 hand-off stages that step 3 v3.3 (Lemma U) needs.  Writes p2/out_v33/FINAL_manifest_v33.json.
Step 4 v3.2 is closed and pinned by p2/out/FINAL_manifest.json (34 stages, 523 paths + 97 upstream + 26 scripts).  This driver never writes
under p2/out/ or into any existing step-4 file; every v3.2 output it reads is checked against the sha256 that FINAL_manifest.json records
for it (last write), and the whole v3.2 manifest is verified at the start and again at the end.  New scripts: p2/v33/; outputs: p2/out_v33/.
Stages (one process at a time; run the driver under nice -n 19; every child is python3 -I -B, so no bytecode is written anywhere):
  (H2f), v3.3 hand-off threshold on [60, 700]: h2f_ift (widened parameter balls), h2f_balls (32 t_c balls, win.py machinery to x = 730),
         t9_h2f (t9.py's error-transfer constants on [60, 700]), h2f_camp (true curve at Ybar = 1e-12; margins);
  tightened (H2) = (1e-8, 1e-6, 1e-3): h2t_v33 (re-assembly of pinned outputs);
  fixed-L cone family (L = 50, X = 3, centres [300, 380], cylinders [r_c - 150, r_c + 150]): t9_cone_v33 (jets, (H5), H3), h4loc_v33 (H4),
         t9_lemmaO_v33 (H1 graph, C0 sqrt t = 0.00227), ratio_v33 (H7);
  notes_v33: far family (centres >= 490) is a subset of the certified v3.2 family; t9_h6_v32 superseded by task 5's H6.
  B0 fix (review REVIEW_v33.md B0): b0_dprof (certified bounds on d = (f - U) - w', the difference of the v3.3 and v3.2 deviation variables,
         on [14, 60], at |alpha| <= alpha(t_1) and at alpha_1), h2u_v33 (the tightened (H2) for f - U_1 on [14, 60], from h2t_v33 + b0_dprof).
  M9 fix (review REVIEW_v33.md M9; stages 12-18, new files only, the earlier stages unchanged): m9_balls -> win_arb (arb cell half-widths,
         exact cell breakpoints, arb sums: (H3) windows, (H2) shells, (H2f) on 32 t_c balls x 2 parameter-ball sets), h3_arb ((H3), v3.2 (H2)
         on the arb bases), h2t_v33_arb, h2u_v33_arb, h2f_camp_arb (arb-base copies of h2t_v33, h2u_v33, h2f_camp), t9_lemmaO_v33_arb
         (Lemma O H1 with the certified centre-to-arclength map), h4far_v33 (Lemma O H4 for every far centre, radius >= 490 to the far end).
  N1 fix (re-review REREVIEW_v33.md N1; stage 19, new file): n1_chart (certified chart terms between the verifier's S1-chart variable and
         (A)'s ray-graph variable on (A)'s ramp [14, 28] and core [28, 60]; lemma writeup/lemmaU/bulk/fix2/n12/chart_lemma.tex).
  N2 fix (re-review N2): every stage runs with an EXPLICIT environment built from scratch (BASE_ENV + the stage's own env), never the
         parent's; each stage record pins env_child (the complete child environment) and env_read (every variable the stage code reads,
         found by grepping os.environ/getenv in the stage scripts and the modules they import: READ_VARS; null = unset).  Note: the child
         interpreter runs with -I, which implies -E, so it IGNORES all PYTHON* variables (PYTHONHASHSEED included); they are set and recorded
         for completeness and for grandchildren, but str hashing in the child is randomized; byte-identity is checked by two clean passes.
  F1/F2 fix (review t8_review/CHECK_N1_N5.md F1, F2; stage 20, new file, standalone pinned check): f12_sliver (certified sliver sups G_i on the
         extended ramp sliver y in [14 - 2^-19, 16.25 - X0] and sup|w''| <= sI2' on x in [16.249, 16.25]; scratch recomputation, from pinned
         inputs only, of n1_chart, the 7 (A) assembly runs, A_interface, Lemma U (2 passes) and zone III with the corrected constants, every
         harness validated by reproducing the pinned file; all claims, margins and quoted numbers re-checked).  Writes only under out_v33/f12/;
         no existing pinned output changes.  Report p2/v33/F12_REPORT.md.
usage: nice -n 19 python3 -I -B final_tail_v33.py [--from STAGE] [--until STAGE]
  --from: rerun from that stage (index or name prefix); the earlier stages are taken from out_v33/FINAL_manifest_v33.partial.json if it
  exists (continuing that pass), else from FINAL_manifest_v33.json, after checking that their definitions are unchanged and that every
  recorded sha256 matches the file now.  --until: stop after that stage (partial file kept).  Progress after each stage is written
  atomically to the partial file; the final manifest is written atomically only after all stages finished and every hash was checked."""
import os, sys, json, hashlib, subprocess, time, glob
S4 = '/home/claude/s1spec/step4'; P2 = f'{S4}/p2'; O = f'{P2}/out'; OV = f'{P2}/out_v33'; V = f'{P2}/v33'; T8 = f'{P2}/t8se'; T6C = f'{P2}/t6c'
V32MAN = f'{O}/FINAL_manifest.json'
FINAL, PARTIAL = f'{OV}/FINAL_manifest_v33.json', f'{OV}/FINAL_manifest_v33.partial.json'
os.makedirs(f'{OV}/logs', exist_ok=True); os.makedirs(f'{OV}/h2f', exist_ok=True); os.makedirs(f'{OV}/m9', exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
ab = lambda p: p if p.startswith('/') else os.path.join(S4, p)
rel = lambda p: os.path.relpath(p, S4) if p.startswith(S4 + '/') else p
def write_atomic(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1); fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, path)

def verify_v32():
    """every sha256 pinned by FINAL_manifest.json (stage paths at their last write, upstream, scripts) against the files now"""
    m = json.load(open(V32MAN)); last = {}
    for P in m['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]
    bad += [u['path'] for u in m['upstream'] if sha(ab(u['path'])) != u['sha256']]
    bad += [q for q, h in m.get('script_sha256', {}).items() if sha(ab(q)) != h]
    return dict(manifest=rel(V32MAN), manifest_sha256=sha(V32MAN), stages=len(m['stages']), paths=len(last), upstream=len(m['upstream']),
                scripts=len(m.get('script_sha256', {})), mismatches=bad, ok=not bad), last, {u['path']: u['sha256'] for u in m['upstream']}

AMP = f'{T8}/constants/amp_A1_w180.json'; N3 = f'{T8}/constants/n3_b1.json'; T7C = f'{S4}/t7/constants.json'; CLAM = f'{O}/t9_h1_clam_s180.json'
T9K = f'{O}/t9_cone_kch.json'; WSUM = f'{O}/h1slab/win_summary.json'; YC = f'{O}/ycomp_final.json'
IFT = [f'{O}/h1slab/ift_{m}.json' for m in range(32)]; IFTW = [f'{OV}/h2f/ift_{m}.json' for m in range(32)]
WINS = [f'{OV}/h2f/win_{m}.json' for m in range(32)]
OM9 = f'{OV}/m9'; M9WIN = [f'{OM9}/win_arb_{S_}_{m}.json' for S_ in ('v32', 'wide') for m in range(32)]
JETS = lambda *sl: [f'{T6C}/jets/jets_{s:03d}.npz' for s in sl]
PKL = f'{S4}/t1/polished_v2/slab_180.pkl'
MS = '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'
LOV32 = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json'
Z3 = '/home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json'
BU = '/home/claude/s1spec/writeup/lemmaU/bulk'; T2V = f'{BU}/fix2/n12/t2_A_v35'; OF = f'{OV}/f12'
_ASM = ['design_N2ref', 'design_N2ref4', 'design_N20', 'lo_N2ref', 'hi_N2ref', 'probe_N2ref', 'design_pm3ibp']
F12_IN = ([f'{V}/f12_sliver.py', f'{V}/n1_chart.py', f'{OV}/n1_chart.json', '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json', MS,
           '/home/claude/s1spec/step3nl_v33/agent_sob/sob_constants_P2.json', '/home/claude/s1spec/step3nl_v33/verifier/constants.py', f'{P2}/t9h/s1chart.py',
           '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
           '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/step3/hs_uforms.pkl', '/home/claude/s1spec/lin/lintail2.py', '/home/claude/s1_cap/tail.py',
           f'{T2V}/rigor/assembly_rep.py', f'{T2V}/rigor/arlib.py', f'{T2V}/rigor/export_t2.py', f'{T2V}/base_ibp/rigor/assembly_arb.py', f'{T2V}/base_ibp/rigor/arbio.py',
           f'{T2V}/base_ibp/rigor/weights_exact.py', f'{T2V}/base_ibp/rigor/weights_exact.json', f'{T2V}/base_ibp/rigor/cap_F10p.json', f'{T2V}/base_ibp/rigor/ramp_tables.json',
           f'{T2V}/base_ibp/rigor/tables_exact.pkl', f'{T2V}/base_ibp/lemmaS_F10p_M0.52.json', f'{T2V}/base_ibp/lemmaS_cflux.json', f'{T2V}/base_ibp/rigor/assembly_arb.json',
           f'{T2V}/base_v31H2/rigor/assembly_arb.json', '/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.json']
          + [f'{T2V}/out/{f}.json' for f in ('weights_t2', 'weights_t2_fb', 'ramp_rep', 'ramp_rep_pm3ibp', 'outer_L2', 'hjw_rep', 'ic_design', 'ic_lo', 'ic_hi', 'ic_probe', 'a3_signs')]
          + [f'{T2V}/out/asm_{n}.json' for n in _ASM] + [f'{T2V}/out/A_interface.json', f'{OV}/b0_dprof.json', f'{OV}/h2u_v33_arb.json',
           f'{BU}/t4_cert/assemble_U.py', f'{BU}/t4_cert/t4lib.py', f'{BU}/t4_cert/out/u1_ghat.json', f'{BU}/t4_cert/out/u2_K.json', f'{BU}/t4_cert/out/u34_primary.json',
           f'{BU}/t4_cert/out/U_cert_v35.json', f'{BU}/t4_cert/out/U_cert_for_zone3_v35.json', f'{BU}/t5_zone3/zone3.py', f'{BU}/t5_zone3/z3lib.py',
           f'{BU}/t5_zone3/out/ramp_U.json', Z3, LOV32, '/home/claude/s1spec/writeup/lemmaO_v32/lemmaR_consts.json',
           f'{BU}/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json', f'{BU}/fix1/lemmaU_rev/constants_needed_v34.json', f'{BU}/fix2/n35_gen_numbers.py'])
F12_OUT = ([f'{OF}/f12_sliver.json', f'{OF}/repro/n1_chart.json', f'{OF}/new/n1_chart.json'] + [f'{OF}/{d}/asm_{n}.json' for n in _ASM for d in ('repro', 'new')]
           + [f'{OF}/repro/A_interface.json', f'{OF}/new/A_interface.json', f'{OF}/repro/U_cert.json', f'{OF}/new/U_cert.json', f'{OF}/new/U_cert_pass2.json',
              f'{OF}/new/U_cert_for_zone3.json', f'{OF}/repro/zone3_interface.json', f'{OF}/new/zone3_interface.json'])
def st(name, cwd, cmd, env, inp, outp, log): return dict(stage=name, cwd=cwd, cmd=cmd, env=env, inputs=inp, outputs=outp, log=f'{OV}/logs/{log}')
# ---- N2: explicit child environment (nothing inherited from the parent shell)
BASE_ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}
READ_VARS = sorted(set('ALIN BASE_T1 BLOOP BSUB CAB CERT_TAG COHERENT CUBFINE EPSI EXE FAR_M FORM_MULT H1JOBS H1NSUB H1ONLY H1REF H1TAG H4OUT HSTEP '
                       'LO_NS MODA MODE2 N3JSON NLFAC OMEGA_Q PDOT POS_C0_MILLI QSHIFT REORG S1PROF S1TT SLAB0 STIPPOINT SUBW T6DIR T6RES T8CONSTS '
                       'T9IFT T9JSON T9OUT T9PIECE T9RATIO T9TC T9WIN TOPLOSS TREF V33_VALIDATE WINOUT YBAR YOUT ZSHIFT '
                       'PYTHONHASHSEED PYTHONDONTWRITEBYTECODE PYTHONPATH PYTHONHOME OMP_NUM_THREADS LC_ALL LANG PATH'.split()))
def child_env(S): return dict(BASE_ENV, **S['env'])
def env_read(S): E = child_env(S); return {k: E.get(k) for k in READ_VARS}
STAGES = [
    st('h2f: widened parameter balls (h2f_ift.py)', V, ['h2f_ift.py'], {}, [f'{V}/h2f_ift.py', CLAM, N3] + IFT, IFTW + [f'{OV}/h2f/ift_widen.json'], 'h2f_ift.log'),
    st('h2f: base sums on [60, 700], 32 t_c balls (h2f_balls.py -> h2f_win.py)', V, ['h2f_balls.py'], {},
       [f'{V}/h2f_balls.py', f'{V}/h2f_win.py', f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py', f'{P2}/t9h/s1chart.py', PKL] + IFTW,
       WINS + [f'{OV}/h2f/h2f_base.json'], 'h2f_balls.log'),
    st('h2f: error-transfer constants on [60, 700] (t9_h2f.py)', P2, ['v33/t9_h2f.py', 't8se/constants/amp_A1_w180.json'], dict(T6DIR=T6C, N3JSON=N3),
       [f'{V}/t9_h2f.py', AMP, N3, T7C, CLAM, f'{P2}/t9chart.py'] + IFT + JETS(177, 178, 179, 180), [f'{OV}/t9_h2f.json'], 't9_h2f.log'),
    st('h2f: (H2f) for the true curve at Ybar = 1e-12 (h2f_camp.py)', V, ['h2f_camp.py'], {}, [f'{V}/h2f_camp.py', f'{OV}/h2f/h2f_base.json', f'{OV}/t9_h2f.json', AMP, YC],
       [f'{OV}/h2f_result.json'], 'h2f_camp.log'),
    st('tightened (H2) shells (h2t_v33.py)', V, ['h2t_v33.py'], {}, [f'{V}/h2t_v33.py', WSUM, T9K, AMP], [f'{OV}/h2t_v33.json'], 'h2t_v33.log'),
    st('fixed-L cone: jets, (H5), H3 (t9_cone_v33.py)', P2, ['v33/t9_cone_v33.py', 't8se/constants/amp_A1_w180.json'], dict(T6DIR=T6C, N3JSON=N3),
       [f'{V}/t9_cone_v33.py', AMP, N3, f'{P2}/curvjets.py'] + JETS(177, 178, 179, 180), [f'{OV}/t9_cone_v33.json'], 't9_cone_v33.log'),
    st('fixed-L cone: H4 (h4loc_v33.py)', T8, ['../v33/h4loc_v33.py', 'constants/amp_A1_w180.json'], dict(T6DIR='../t6c', H4OUT=f'{OV}/h4loc_v33.json', T9JSON='../out/t9_cone_kch.json'),
       [f'{V}/h4loc_v33.py', AMP, f'{T6C}/residuals.json', T9K, f'{P2}/t9chart.py'] + JETS(180), [f'{OV}/h4loc_v33.json'], 'h4loc_v33.log'),
    st('fixed-L cone: Lemma O H1 graph (t9_lemmaO_v33.py, LO_NS=48)', P2, ['v33/t9_lemmaO_v33.py', 't8se/constants/amp_A1_w180.json', '180'], dict(T6DIR=T6C, LO_NS='48'),
       [f'{V}/t9_lemmaO_v33.py', AMP, f'{P2}/t9h/basegeo.py', PKL] + JETS(180), [f'{OV}/t9_lemmaO_v33.json'], 't9_lemmaO_v33.log'),
    st('fixed-L cone: H7 time ratio (ratio_v33.py)', V, ['ratio_v33.py'], {}, [f'{V}/ratio_v33.py', MS, LOV32], [f'{OV}/ratio_v33.json'], 'ratio_v33.log'),
    st('covering notes: far family subset, t9_h6_v32 superseded (notes_v33.py)', V, ['notes_v33.py'], {},
       [f'{V}/notes_v33.py', Z3, f'{O}/t9_lemmaO_far.json', f'{O}/t9_far_v32.json', f'{T8}/constants/h4loc_v32.json'], [f'{OV}/notes_v33.json'], 'notes_v33.log'),
    st('B0: profile difference d = (f - U) - w\' on [14, 60] (b0_dprof.py)', V, ['b0_dprof.py'], {},
       [f'{V}/b0_dprof.py', '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
        f'{P2}/t9h/s1chart.py'] + IFTW, [f'{OV}/b0_dprof.json'], 'b0_dprof.log'),
    st('B0: tightened (H2) for f - U_1 on [14, 60] (h2u_v33.py)', V, ['h2u_v33.py'], {}, [f'{V}/h2u_v33.py', f'{OV}/h2t_v33.json', f'{OV}/b0_dprof.json'],
       [f'{OV}/h2u_v33.json'], 'h2u_v33.log'),
    # ---- M9 (review REVIEW_v33.md M9) ----
    st('M9: arb window/shell/(H2f) base sums, 32 t_c balls x 2 ball sets (m9_balls.py -> win_arb.py)', V, ['m9_balls.py'], {},
       [f'{V}/m9_balls.py', f'{V}/win_arb.py', f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py', f'{P2}/t9h/s1chart.py',
        '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', PKL] + IFT + IFTW + [WSUM, f'{OV}/h2f/h2f_base.json'] + WINS + [f'{O}/h1slab/win_{m}.json' for m in range(32)],
       M9WIN + [f'{OM9}/win_arb_summary.json'], 'm9_balls.log'),
    st('M9: (H3) windows and v3.2 (H2) shells on the arb bases (h3_arb.py)', V, ['h3_arb.py'], {}, [f'{V}/h3_arb.py', f'{OM9}/win_arb_summary.json', WSUM, T9K, AMP],
       [f'{OV}/h3_arb.json'], 'h3_arb.log'),
    st('M9: tightened (H2) shells on the arb bases (h2t_v33_arb.py)', V, ['h2t_v33_arb.py'], {}, [f'{V}/h2t_v33_arb.py', f'{OM9}/win_arb_summary.json', T9K, AMP],
       [f'{OV}/h2t_v33_arb.json'], 'h2t_v33_arb.log'),
    st('M9: tightened (H2) for f - U_1 on the arb bases (h2u_v33_arb.py)', V, ['h2u_v33_arb.py'], {}, [f'{V}/h2u_v33_arb.py', f'{OV}/h2t_v33_arb.json', f'{OV}/b0_dprof.json'],
       [f'{OV}/h2u_v33_arb.json'], 'h2u_v33_arb.log'),
    st('M9: (H2f) on the arb bases (h2f_camp_arb.py)', V, ['h2f_camp_arb.py'], {}, [f'{V}/h2f_camp_arb.py', f'{OM9}/win_arb_summary.json', f'{OV}/t9_h2f.json', AMP, YC],
       [f'{OV}/h2f_result_arb.json'], 'h2f_camp_arb.log'),
    st('M9: Lemma O H1 with the certified centre map (t9_lemmaO_v33_arb.py, LO_NS=48)', P2, ['v33/t9_lemmaO_v33_arb.py', 't8se/constants/amp_A1_w180.json', '180'],
       dict(T6DIR=T6C, LO_NS='48'), [f'{V}/t9_lemmaO_v33_arb.py', AMP, f'{P2}/t9h/basegeo.py', PKL, T9K] + JETS(180), [f'{OV}/t9_lemmaO_v33_arb.json'], 't9_lemmaO_v33_arb.log'),
    st('M9: far-family H4, every centre at radius >= 490 to the far end (h4far_v33.py)', T8, ['../v33/h4far_v33.py', 'constants/amp_A1_w180.json'],
       dict(T6DIR='../t6c', T9JSON='../out/t9_cone_kch.json'),
       [f'{V}/h4far_v33.py', f'{T8}/h4loc.py', f'{P2}/t9chart.py', AMP, T9K, LOV32, f'{T8}/constants/constants_b1.json', f'{T8}/constants/h4loc_v32.json', Z3] + JETS(180),
       [f'{OV}/h4far_v33.json'], 'h4far_v33.log'),
    # ---- N1 (re-review REREVIEW_v33.md N1) ----
    st('N1: chart terms, S1 chart vs ray graph, on (A)\'s ramp and core (n1_chart.py)', V, ['n1_chart.py'], {},
       [f'{V}/n1_chart.py', '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json', '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py',
        '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py', f'{P2}/t9h/s1chart.py', '/home/claude/s1spec/step3/hstail.py',
        '/home/claude/s1spec/step3/hs_uforms.pkl', '/home/claude/s1spec/lin/lintail2.py', '/home/claude/s1_cap/tail.py'],
       [f'{OV}/n1_chart.json'], 'n1_chart.log'),
    # ---- F1, F2 (review t8_review/CHECK_N1_N5.md; standalone pinned check, writes only under out_v33/f12/) ----
    st('f12_sliver: F1 sliver G_i, F2 sliver sup|w\'\'|, pinned claims re-checked with the corrected constants (f12_sliver.py)', V, ['f12_sliver.py'], {},
       F12_IN, F12_OUT, 'f12_sliver.log'),
]
READONLY_STEP4 = [f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py', f'{P2}/t9h/s1chart.py', f'{P2}/t9chart.py', f'{P2}/curvjets.py',
                  f'{S4}/t6/reader.py', f'{P2}/t9.py', f'{T8}/h4loc.py', f'{P2}/t9_lemmaO.py']

if __name__ == '__main__':
    v32_start, V32LAST, V32UP = verify_v32(); assert v32_start['ok'], ('v3.2 manifest does not verify', v32_start['mismatches'])
    print(f"v3.2 manifest verified: {v32_start['stages']} stages, {v32_start['paths']} paths, {v32_start['upstream']} upstream, {v32_start['scripts']} scripts", flush=True)
    man = dict(note=__doc__, v32_manifest=dict(path=rel(V32MAN), sha256=v32_start['manifest_sha256']), stages=[], checks=dict(v32_manifest_at_start=v32_start),
               environment=dict(base_env=BASE_ENV, read_vars=READ_VARS, interpreter=sys.executable, flags=['-I', '-B'],
                                note='children get child_env = base_env + stage env only (N2); -I implies -E: PYTHON* variables are ignored by the child interpreter',
                                parent_values_not_passed={k: os.environ[k] for k in READ_VARS if k in os.environ and k not in ('PATH', 'LANG', 'LC_ALL')}))
    i0, i1 = 0, len(STAGES) - 1
    def stage_index(key, start=0):
        return int(key) if key.isdigit() else next(i for i, S in enumerate(STAGES) if i >= start and S['stage'].startswith(key))
    if '--from' in sys.argv:
        i0 = stage_index(sys.argv[sys.argv.index('--from') + 1]); cont = os.path.exists(PARTIAL)
        prev = json.load(open(PARTIAL if cont else FINAL)); last = {}
        assert i0 <= len(prev['stages'])
        for S, P in zip(STAGES[:i0], prev['stages'][:i0]):
            assert P['stage'] == S['stage'] and P['cmd'] == ['python3', '-I', '-B'] + S['cmd'] and P['env'] == S['env'] and P['cwd'] == rel(S['cwd']), S['stage']
            assert P.get('env_child') == child_env(S) and P.get('env_read') == env_read(S), ('reused stage ran with another environment', S['stage'])
            assert [f['path'] for f in P['inputs']] == [rel(q) for q in S['inputs']] and [f['path'] for f in P['outputs']] == [rel(q) for q in S['outputs'] + [S['log']]], S['stage']
            for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
            for f in P['outputs']: last[f['path']] = f['sha256']
        bad = [q for q, h in last.items() if sha(ab(q)) != h]; assert not bad, bad
        man['stages'] = prev['stages'][:i0]
        man['checks'].update({k: v for k, v in prev['checks'].items() if k in ('segments',) and cont})
        man['checks'].setdefault('segments', []).append([i0, None])
        print(f'from stage {i0}: stages 0-{i0 - 1} reused ({len(last)} paths verified)', flush=True)
    if '--until' in sys.argv: i1 = stage_index(sys.argv[sys.argv.index('--until') + 1], i0)
    for S in STAGES[i0:i1 + 1]:
        ins = []
        for p in S['inputs']:
            h = sha(p); assert h != 'MISSING', p
            rp = rel(p); pinned = V32LAST.get(rp) or V32UP.get(rp)
            if pinned is not None: assert pinned == h, ('v3.2 pinned input changed', rp)
            ins.append(dict(path=rp, sha256=h, v32_pinned=pinned is not None))
        env = child_env(S); t0 = time.time()                       # N2: explicit, recorded; never the parent's environment
        with open(S['log'], 'w') as lf:
            r = subprocess.run([sys.executable, '-I', '-B'] + S['cmd'], cwd=S['cwd'], env=env, stdout=lf, stderr=subprocess.STDOUT)
        print(f"{S['stage']}: rc {r.returncode} ({time.time() - t0:.0f}s)", flush=True)
        assert r.returncode == 0, (S['stage'], open(S['log']).read()[-2000:])
        for p in S['outputs']: assert p.startswith(OV + '/'), ('output outside out_v33', p)
        outs = [dict(path=rel(p), sha256=sha(p)) for p in S['outputs'] + [S['log']]]
        man['stages'].append(dict(stage=S['stage'], cwd=rel(S['cwd']), cmd=['python3', '-I', '-B'] + S['cmd'], env=S['env'], env_child=env, env_read=env_read(S),
                                  inputs=ins, outputs=outs, seconds=round(time.time() - t0, 1)))
        write_atomic(PARTIAL, man)
    if man['checks'].get('segments'): man['checks']['segments'][-1][1] = i1
    if i1 < len(STAGES) - 1:
        man['checks']['in_progress'] = dict(next=i1 + 1, next_stage=STAGES[i1 + 1]['stage']); write_atomic(PARTIAL, man)
        print(f'stopped after stage {i1}; continue with --from {i1 + 1}'); sys.exit(0)
    assert [P['stage'] for P in man['stages']] == [S['stage'] for S in STAGES], 'incomplete stage list'
    seen = {}
    for P in man['stages']:
        for i in P['inputs']: seen.setdefault(i['path'], set()).add(i['sha256'])
    man['checks']['inputs_with_two_versions'] = sorted(p for p, v in seen.items() if len(v) > 1); assert not man['checks']['inputs_with_two_versions']
    man['reused_v32'] = sorted({i['path']: i['sha256'] for P in man['stages'] for i in P['inputs'] if i['v32_pinned']}.items())
    man['script_sha256'] = {rel(f): sha(f) for f in sorted(glob.glob(f'{V}/*.py')) + READONLY_STEP4}
    res = {}
    def J(p): return json.load(open(p))
    h2f = J(f'{OV}/h2f_result.json'); res['H2f'] = dict(thresholds=h2f['thresholds'], rows=[{k: r[k] for k in ('j', 'base', 'with_error', 'margin', 'Ybar_limit')} for r in h2f['rows']], ok=h2f['all_ok'])
    h2t = J(f'{OV}/h2t_v33.json'); res['H2_tightened'] = dict(rows=[{k: r[k] for k in ('j', 'base', 'with_error', 'margin')} for r in h2t['rows']], ok=h2t['all_ok'])
    h2u = J(f'{OV}/h2u_v33.json'); res['H2_tightened_fU1'] = dict(rows=[{k: r[k] for k in ('j', 'E_wprime_with_error', 'E_d', 'bound', 'margin')} for r in h2u['rows']],
                                                                 ok=h2u['all_ok'] and h2u['all_margin_ge_2'])
    cn = J(f'{OV}/t9_cone_v33.json'); res['cone_fixedL'] = dict(checks=cn['checks'], ok=cn['all_ok'])
    h4 = J(f'{OV}/h4loc_v33.json'); res['H4_fixedL'] = h4['check']
    lo = J(f'{OV}/t9_lemmaO_v33.json'); res['H1_fixedL'] = dict(ok=lo['H1_all'], sup_h0_max=lo['sup_h0_max'], H_max=lo['H_max'], C0_sqrt_t1=lo['C0_sqrt_t1'])
    ra = J(f'{OV}/ratio_v33.json'); res['H7_fixedL'] = dict({k: ra[k] for k in ('ratio_upper', 'needed', 'ratio_le_needed', 't_star_lower', 'margin_tstar_over_ratio')}, ok=ra['ratio_le_needed'])
    wa = J(f'{OM9}/win_arb_summary.json'); res['M9_base_sums'] = dict(combined=wa['combined'], compare_float={k: v for k, v in wa['compare_float'].items() if k != 'slab'}, ok=True)
    h3 = J(f'{OV}/h3_arb.json'); res['M9_H3_H2v32'] = dict(rows={k: [{q: r[q] for q in ('j', 'with_error_arb', 'margin_arb', 'margin_float', 'Ybar_limit_arb')} for r in v] for k, v in h3['rows'].items()},
                                                        binding_Ybar_limit=h3['binding_Ybar_limit_arb'], ok=h3['all_ok'])
    h2ta = J(f'{OV}/h2t_v33_arb.json'); res['M9_H2_tightened'] = dict(rows=[{k: r[k] for k in ('j', 'base', 'with_error', 'margin')} for r in h2ta['rows']], ok=h2ta['all_ok'])
    h2ua = J(f'{OV}/h2u_v33_arb.json'); res['M9_H2_tightened_fU1'] = dict(rows=[{k: r[k] for k in ('j', 'E_wprime_with_error', 'E_d', 'bound', 'margin')} for r in h2ua['rows']],
                                                                     ok=h2ua['all_ok'] and h2ua['all_margin_ge_2'])
    h2fa = J(f'{OV}/h2f_result_arb.json'); res['M9_H2f'] = dict(rows=[{k: r[k] for k in ('j', 'base', 'with_error', 'margin', 'Ybar_limit')} for r in h2fa['rows']], ok=h2fa['all_ok'])
    loa = J(f'{OV}/t9_lemmaO_v33_arb.json'); res['M9_H1_centre_map'] = dict(sup_h0_max=loa['sup_h0_max'], H_max=loa['H_max'], centre_map_ok=loa['centre_map_ok'], ok=loa['H1_all'] and loa['centre_map_ok'])
    h4f = J(f'{OV}/h4far_v33.json'); res['M9_H4_far'] = dict(min_margin=h4f['min_margin'], worst=h4f['worst'], n_rows=h4f['n_rows'], XSTART=h4f['XSTART'], ok=h4f['all_ok'])
    n1 = J(f'{OV}/n1_chart.json'); res['N1_chart'] = dict(ramp={k: n1['ramp'][k] for k in ('sup_add', 'X_W', 'beta22', 'beta32', 'beta33', 'Dbeta22', 'eps', 'e', 'rho_meas')},
                                                          core={k: n1['core'][k] for k in ('sup_add', 'E', 'beta22', 'beta3', 'e', 'rho_meas')}, dB0_slivers=n1['dB0_slivers'],
                                                          ok=n1['ramp']['e'] < 2 ** -19 and n1['ramp']['sup_add'][0] <= 2e-7)
    f12 = J(f'{OF}/f12_sliver.json')
    res['F12_sliver'] = dict(G_old=f12['F1']['G_old'], G_new=f12['F1']['G_new'], sI2=f12['F2']['rows']['1.5e-4']['sI2'], sI2_prime=f12['F2']['sI2_prime_used'],
                             named={k: f12['named'][k] for k in ('ramp_M1', 'K6', 'H6', 'S0v', 'flux_beta2_4N2ref', 'EA_ref', 'EA_ref4', 'Ehat_U')},
                             n_checks=len(f12['checks']), quoted_numbers_failing=f12['quoted_numbers_failing'], macro_changes=f12['macro_changes'],
                             rounding_flips=f12['rounding_flips'], ok=f12['all_ok'])
    man['results'] = res; man['checks']['all_items_ok'] = all(v.get('ok', False) for v in res.values())
    last = {}
    for P in man['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(ab(q)) != h]; assert not bad, bad
    v32_end, _, _ = verify_v32(); man['checks']['v32_manifest_at_end'] = v32_end; assert v32_end['ok'], v32_end['mismatches']
    man['checks']['hashes_verified_at_end'] = dict(paths=len(last), scripts=len(man['script_sha256']), reused_v32=len(man['reused_v32']))
    write_atomic(FINAL, man)
    if os.path.exists(PARTIAL): os.remove(PARTIAL)
    print('results:', json.dumps(res)[:3000]); print('checks:', json.dumps({k: v for k, v in man['checks'].items() if k not in ('v32_manifest_at_start', 'v32_manifest_at_end')}))
    print(f"v3.2 manifest at end: ok = {v32_end['ok']} ({v32_end['paths']} paths, {v32_end['upstream']} upstream, {v32_end['scripts']} scripts)")
