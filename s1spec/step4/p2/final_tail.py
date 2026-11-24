"""Final hand-off tail on the B1 data, one version of every input (coordinator, Nov 12).  Writes p2/out/FINAL_manifest.json.
Data sets: T6 = p2/t6c (residuals.json, jets/; the (SE)-corrected base, T6's engine, B1); T7 = t7/constants.json (current);
T8 = constants_b1.json (cconst + merge on p2/t6c and T7 slabs >= 2, rebuilt here); base windows T9WIN = out/h1slab/win_summary.json.
Superseded and read by no stage: p2/t6/*, step4/t6/residuals.json, se_corr.json, constants_fin.json, n3.json, amp_fin*.
Every camp run reads the final t9_cone_kch.json: the T8 iteration and the window camp run twice (passes A, B); pass B's t9 must
reproduce pass A's byte for byte (t9 reads only the amp's per-slab aF, bF and at_Ybar values, which do not depend on T9JSON).
Cleanup (Nov 12): Lemma O H1 for the far family (t9_lemmaO_far.py, certified pass/fail, arc by arc) and the informational axis crossing
(t9_axis_crossing.py) are stages; h4loc derives its chart factor (t9chart Taylor mode; reads t9_cone_kch.json); t1clock includes
lam/lam~ and x_NL (reads t7/constants.json and switch_s2.json).  Follow-up: h1.py's rotation-mode term is certified (E_I part by cells, M over
the Krawczyk box); the 16 h1 exact points (h1_points.py) and the H1 slab summary (h1_summary.py) are stages before clam, no longer upstream.
The B8 check (s1_tight/consistency/check_outA5.py) is a stage too (the last one: no stage reads it), no longer upstream.
usage: python3 -I final_tail.py [--from STAGE] [--until STAGE]   (STAGE: an index, or the start of a stage name; --from: rerun from that
       stage, the earlier ones taken from the current FINAL_manifest.json after checking that their definitions are unchanged and that the
       last recorded SHA-256 of every path they read or wrote matches the file now; --until: stop after that stage, leaving a partial
       progress file marked checks.in_progress.  FINAL_manifest.json is written only when every stage has finished and every recorded
       SHA-256 has been checked against the files, atomically (temp file + rename); until then the progress after each stage is in
       out/FINAL_manifest.partial.json (also atomic).  If that file exists, --from continues its pass (stopped by --until, or killed): the
       earlier stages are taken from it, and the manifest keeps the pass's rerun_from and lists its segments; delete it to start a new pass.
       The H1 stage is itself resumable point by point (t9h/h1_points.py).)"""
import os, sys, json, hashlib, subprocess, time, glob
S4 = '/home/claude/s1spec/step4'; P2 = f'{S4}/p2'; T8 = f'{P2}/t8se'; O = f'{P2}/out'; T6C = f'{P2}/t6c'
F0 = '4e-23,3e-22,5e-19,3.4e-15'; YB = '1e-12'
S1ENV = dict(COHERENT='1', S1PROF='/home/claude/s1_tight/ext_tip/prof_tip.txt', S1TT='10.249077651120778,10.249077651121361')
CAMPENV = dict(T9JSON='../out/t9_cone_kch.json', T9WIN='../out/h1slab/win_summary.json', T6DIR='../t6c', OMP_NUM_THREADS='1')
JETS = lambda *sl: [f'{T6C}/jets/jets_{s:03d}.npz' for s in sl]
T7C = f'{S4}/t7/constants.json'; CB1 = f'{T8}/constants/constants_b1.json'; AMP = f'{T8}/constants/amp_A1_w180.json'
LOV32 = f'{S4}/../writeup/lemmaO_v32/lemmaO_v32_consts.json'                  # step-3 v3.2 Lemma O constants (read only)
H1PTS = [f'{O}/t9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json' for k in range(16)]   # the 16 h1 exact points of slab 180
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
FINAL, PARTIAL = f'{O}/FINAL_manifest.json', f'{O}/FINAL_manifest.partial.json'
B8D = '/home/claude/s1_tight'                                  # step-1 tight S1 enclosures and their B8 consistency check (read only, except its JSON)
def write_atomic(path, obj):                                   # temp file + rename: a reader never sees a half-written manifest
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1); fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, path)
rel = lambda p: os.path.relpath(p, S4)


def st(name, cwd, cmd, env, inp, outp, log):
    return dict(stage=name, cwd=cwd, cmd=cmd, env=env, inputs=inp, outputs=outp, log=log)


def t8_iteration():
    e = dict(CAMPENV, SUBW='1.46484375e-4')                   # B1's T8 iteration (W' sub-window for its own report; per-slab data only used here)
    return [
        st('T8 camp (no N3)', T8, ['camp.py', '--consts', 'constants/constants_b1.json', '--F0', F0, '--out', 'constants/amp_b1.json'], e,
           [CB1, T7C, f'{T6C}/residuals.json', f'{O}/t9_cone_kch.json', f'{O}/h1slab/win_summary.json'], [f'{T8}/constants/amp_b1.json'], f'{T8}/constants/amp_b1.log'),
        st('T8 n3far it0', T8, ['n3far.py', '--slabs', '2-180', '--amp', 'constants/amp_b1.json', '--consts', 'constants/constants_b1.json', '--out', 'constants/n3_b1_it0.json'],
           dict(T9JSON=CAMPENV['T9JSON'], T6DIR='../t6c', OMP_NUM_THREADS='1'), [f'{T8}/constants/amp_b1.json', CB1, T7C, f'{T6C}/residuals.json'],
           [f'{T8}/constants/n3_b1_it0.json'], f'{T8}/constants/n3_b1_it0.log'),
        st('T8 camp --n3 it0', T8, ['camp.py', '--consts', 'constants/constants_b1.json', '--F0', F0, '--n3', 'constants/n3_b1_it0.json', '--out', 'constants/amp_b1_n3_it0.json'], e,
           [CB1, f'{T8}/constants/n3_b1_it0.json', T7C, f'{T6C}/residuals.json', f'{O}/t9_cone_kch.json', f'{O}/h1slab/win_summary.json'], [f'{T8}/constants/amp_b1_n3_it0.json'],
           f'{T8}/constants/amp_b1_n3_it0.log'),
        st('T8 n3far (final boxes)', T8, ['n3far.py', '--slabs', '2-180', '--amp', 'constants/amp_b1_n3_it0.json', '--consts', 'constants/constants_b1.json', '--out', 'constants/n3_b1.json'],
           dict(T9JSON=CAMPENV['T9JSON'], T6DIR='../t6c', OMP_NUM_THREADS='1'), [f'{T8}/constants/amp_b1_n3_it0.json', CB1, T7C, f'{T6C}/residuals.json'],
           [f'{T8}/constants/n3_b1.json'], f'{T8}/constants/n3_b1.log')]


def window_and_t9(tag):
    return [
        st(f'camp --n3 n3_b1 --window 180-180 (T9WIN) [{tag}]', T8, ['camp.py', '--consts', 'constants/constants_b1.json', '--F0', F0, '--n3', 'constants/n3_b1.json',
           '--window', '180-180', '--out', 'constants/amp_A1_w180.json'], dict(CAMPENV),
           [CB1, f'{T8}/constants/n3_b1.json', T7C, f'{T6C}/residuals.json', f'{O}/t9_cone_kch.json', f'{O}/h1slab/win_summary.json'], [AMP], f'{O}/final_camp_w180_{tag}.log'),
        st(f'n3boxes [{tag}]', T8, ['n3boxes.py', 'constants/n3_b1.json', 'constants/amp_A1_w180.json', 'constants/constants_b1.json', YB], {},
           [f'{T8}/constants/n3_b1.json', AMP, CB1, f'{T8}/constants/amp_b1_n3_it0.json', T7C], [], f'{O}/n3boxes_final_{tag}.log'),
        st(f't9.py (cone jets, KCH, c0, c0tip, shells, frame) [{tag}]', P2, ['t9.py', 't8se/constants/amp_A1_w180.json'],
           dict(T6DIR=T6C, N3JSON=f'{T8}/constants/n3_b1.json'),
           [AMP, f'{T8}/constants/n3_b1.json', T7C, f'{O}/t9_h1_clam_s180.json', f'{P2}/t9chart.py', f'{P2}/curvjets.py'] + JETS(177, 178, 179, 180)
           + sorted(glob.glob(f'{O}/h1slab/ift_*.json')), [f'{O}/t9_cone_kch.json'], f'{O}/t9_final_{tag}.log')]


STAGES = [
    st('T8 cconst (per slab, T7 slabs >= 2 only)', T8, ['cconst.py', '--slabs', '2-180', '--out', 'constants/slabs_b1', '--jobs', '2'],
       dict(TOPLOSS='0.0196', T6DIR='../t6c', T7SLAB0='2', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
       [f'{T6C}/residuals.json', T7C, f'{S4}/t7/t8req.json'] + JETS(*range(2, 181)), [f'{T8}/constants/slabs_b1/slab_{i:03d}.json' for i in range(2, 181)],
       f'{O}/cconst_final.log'),
    st('T8 merge', T8, ['merge.py', '--in', 'constants/slabs_b1', '--first', '2', '--last', '180', '--out', 'constants/constants_b1.json'],
       dict(T6DIR='../t6c', T7SLAB0='2', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
       [f'{T8}/constants/slabs_b1/slab_{i:03d}.json' for i in range(2, 181)] + [T7C, f'{T6C}/residuals.json'], [CB1], f'{O}/merge_final.log'),
    st('switch (P1 -> P2 at slab 2), T6 = t6c', P2, ['switch.py', '../t3/p1/ckpt/slab_00096.json', '2', '0.873067570143915'], dict(T6DIR=T6C),
       [f'{S4}/t3/p1/ckpt/slab_00096.json', f'{T6C}/residuals.json', T7C, f'{S4}/t7/run_v2/out_farg/farg_002.json'] + JETS(2),
       [f'{O}/switch_s2.json'], f'{O}/switch_final.log'),
    st('ycomp (Y0 = 6.21e-13, SLAB0 = 2)', P2, ['ycomp.py', '6.21e-13'], dict(T6RES=f'{T6C}/residuals.json', SLAB0='2', YOUT=f'{O}/ycomp_final.json'),
       [f'{T6C}/residuals.json', T7C], [f'{O}/ycomp_final.json'], f'{O}/ycomp_final.log'),
    st('H1 exact points (16 x h1.py p<t_k>: Krawczyk, Lambda, certified rotation-mode term)', f'{P2}/t9h', ['h1_points.py'], dict(S1ENV, H1JOBS='2'),
       [f'{S4}/t1/polished_v2/slab_180.pkl', S1ENV['S1PROF'], '/home/claude/s1_tight/ext813/ext813_prof.txt', '/home/claude/s1spec/step3/outA5.txt',
        '/home/claude/s1spec/step3nl/agent_mod/lemma/lyap_extra.json', f'{P2}/t9h/h1.py', f'{P2}/t9h/s1chart2.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py',
        '/home/claude/s1spec/step3nl/agent_mod/ibp/ibp_boxes.py', '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3/hstail.py'],
       H1PTS, f'{O}/h1_points_final.log'),
    st('H1 slab summary (h1_summary.py)', f'{P2}/t9h', ['h1_summary.py'], {}, H1PTS + sorted(glob.glob(f'{O}/h1slab/dt_*_0.03125.json')),
       [f'{O}/h1slab/summary.json'], f'{O}/h1_summary_final.log'),
    st('clam (error transfer; mu_e from the B1 constants)', f'{P2}/t9h', ['clam.py'], dict(S1ENV, T8CONSTS=CB1),
       [CB1, f'{O}/h1slab/summary.json'] + [f'{O}/t9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json' for k in range(16)], [f'{O}/t9_h1_clam_s180.json'], f'{O}/clam_final.log'),
] + [dict(x, pass_='A') for x in t8_iteration() + window_and_t9('A')] + t8_iteration() + window_and_t9('B') + [
    st('t9_far.py', P2, ['t9_far.py', 't8se/constants/amp_A1_w180.json', 't8se/constants/constants_b1.json'], dict(T6DIR=T6C), [AMP, CB1] + JETS(180), [f'{O}/t9_far.json'], f'{O}/t9_far_final.log'),
    st('h4loc on the final amp', T8, ['h4loc.py', 'constants/amp_A1_w180.json'], dict(T6DIR='../t6c', H4OUT='../out/h4loc_A1_w180.json', T9JSON='../out/t9_cone_kch.json'),
       [AMP, f'{T6C}/residuals.json', f'{O}/t9_cone_kch.json', f'{P2}/t9chart.py'] + JETS(180), [f'{O}/h4loc_A1_w180.json'], f'{O}/h4loc_final.log'),
    st('Lemma O H1/H8 (t9_lemmaO.py, LO_NS=48)', P2, ['t9_lemmaO.py', 't8se/constants/amp_A1_w180.json', '180'], dict(T6DIR=T6C, LO_NS='48'), [AMP] + JETS(180),
       [f'{O}/t9_lemmaO.json'], f'{O}/t9_lemmaO_final.log'),
    st('Lemma O H6 (t9_h6.py)', P2, ['t9_h6.py'], {}, [f'{O}/t9_cone_kch.json', f'{O}/t9_far.json', f'{S4}/../writeup/lemmaO/lemmaO_v3_consts.json'], [f'{O}/t9_h6.json'], f'{O}/t9_h6_final.log'),
    st('v3.2: t9_far.py --v32', P2, ['t9_far.py', 't8se/constants/amp_A1_w180.json', 't8se/constants/constants_b1.json', '--v32'], dict(T6DIR=T6C), [AMP, CB1] + JETS(180),
       [f'{O}/t9_far_v32.json'], f'{O}/t9_far_v32_final.log'),
    st('v3.2: h4loc FAR_M=0.15', T8, ['h4loc.py', 'constants/amp_A1_w180.json'], dict(T6DIR='../t6c', FAR_M='0.15', H4OUT='constants/h4loc_v32.json', T9JSON='../out/t9_cone_kch.json'),
       [AMP, f'{T6C}/residuals.json', f'{O}/t9_cone_kch.json', f'{P2}/t9chart.py'] + JETS(180), [f'{T8}/constants/h4loc_v32.json'], f'{O}/h4loc_v32_final.log'),
    st('v3.2: t9_v32_checks.py', P2, ['t9_v32_checks.py'], {}, [f'{O}/t9_cone_kch.json'], [f'{O}/t9_v32_checks.json'], f'{O}/t9_v32_checks_final.log'),
    st('v3.2: t9_h6_v32.py', P2, ['t9_h6_v32.py'], {}, [f'{O}/t9_cone_kch.json', f'{O}/t9_far_v32.json', f'{T8}/constants/h4loc_v32.json', f'{P2}/h6/fable/h6f_cases.json',
       f'{S4}/../writeup/lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json'], [f'{O}/t9_h6_v32.json'], f'{O}/t9_h6_v32_final.log'),
    st('v3.2: t9_lemmaO.py --v32', P2, ['t9_lemmaO.py', 't8se/constants/amp_A1_w180.json', '180', '--v32'], dict(T6DIR=T6C, LO_NS='48'), [AMP] + JETS(180),
       [f'{O}/t9_lemmaO_v32.json'], f'{O}/t9_lemmaO_v32_final.log'),
    st('Lemma O H1 far family (t9_lemmaO_far.py, v3.2, arc by arc)', P2, ['t9_lemmaO_far.py', 't8se/constants/amp_A1_w180.json', 'out/t9_far_v32.json',
       '../../writeup/lemmaO_v32/lemmaO_v32_consts.json'], dict(T6DIR=T6C, T9JSON=f'{O}/t9_cone_kch.json'),
       [AMP, f'{O}/t9_far_v32.json', LOV32, f'{O}/t9_cone_kch.json'] + JETS(180), [f'{O}/t9_lemmaO_far.json'], f'{O}/t9_lemmaO_far_final.log'),
    st('axis crossing (informational; t9_axis_crossing.py)', P2, ['t9_axis_crossing.py', 't8se/constants/amp_A1_w180.json', '180'], dict(AX_NT='16', AX_NG='32'),
       [AMP, f'{P2}/t9h/basegeo.py', f'{S4}/t1/polished_v2/manifest.json', f'{S4}/t1/polished_v2/slab_180.pkl'], [f'{O}/t9_axis_crossing.json'],
       f'{O}/t9_axis_crossing_final.log'),
    st('t1clock', P2, ['t1clock.py', 't8se/constants/constants_b1.json', '180', YB], dict(T6DIR=T6C, SWITCHJSON=f'{O}/switch_s2.json'),
       [CB1, f'{T6C}/residuals.json', T7C, f'{O}/switch_s2.json'] + JETS(180), [f'{O}/t1clock.json'], f'{O}/t1clock_final.log'),
    st('S1 hull containment (B8; check_outA5.py; read by no stage, so last)', f'{B8D}/consistency', ['-B', 'check_outA5.py'], {},
       [f'{B8D}/consistency/check_outA5.py', '/home/claude/s1spec/step3/outA5.txt', '/home/claude/s1spec/step3/inA5.txt', f'{B8D}/ext_tip/prof_tip.txt',
        f'{B8D}/ext_tip/ext_tip_check.json', f'{B8D}/prof/prof_tight.txt', f'{B8D}/ext813/ext813_prof.txt', f'{B8D}/ext813/ext813_points.txt',
        '/home/claude/s1spec/step3/gen_hs2.py', '/home/claude/s1spec/step3/hstail.py'], [f'{B8D}/consistency/check_outA5.json'], f'{O}/check_outA5_final.log'),
]

if __name__ == '__main__':
    man = dict(note=__doc__, data=dict(T6=rel(T6C), T7=rel(T7C), T8=rel(CB1), T9WIN='p2/out/h1slab/win_summary.json'),
               upstream=[], stages=[], checks={})
    A1 = {f['path']: f for f in json.load(open(f'{O}/A1_manifest.json'))['files']}
    for p, f in A1.items():                                    # unchanged upstream: S1 hulls and the H1 chain (base polynomials only)
        if any(k in f['stage'] for k in ('h1dt', 'win_slab', 'win.py', 'windows/shells')):     # B8, the h1 exact points and the summary: stages
            pp = p if p.startswith('/') else os.path.join(S4, p); assert sha(pp) == f['sha256'], p
            man['upstream'].append(dict(stage=f['stage'], path=p, sha256=f['sha256']))
    t9A = None; i0 = 0; i1 = len(STAGES) - 1
    def stage_index(key, start=0):
        return int(key) if key.isdigit() else next(i for i, S in enumerate(STAGES) if i >= start and S['stage'].startswith(key))
    if '--from' in sys.argv:                                   # rerun from a stage: reuse the earlier stages of the current manifest, verified
        i0 = stage_index(sys.argv[sys.argv.index('--from') + 1])
        cont = os.path.exists(PARTIAL)                         # continue the pass in progress, or start a new one from the complete manifest
        prev = json.load(open(PARTIAL if cont else FINAL)); last = {}
        assert i0 <= len(prev['stages']), ('--from beyond the stages recorded', i0, len(prev['stages']))
        for S, P in zip(STAGES[:i0], prev['stages'][:i0]):
            assert P['stage'].startswith(S['stage']) and P['cmd'] == ['python3', '-I'] + S['cmd'] and P['env'] == S['env'] and P['cwd'] == rel(S['cwd']), S['stage']
            assert [f['path'] for f in P['inputs']] == [rel(q) for q in S['inputs']] and [f['path'] for f in P['outputs']] == [rel(q) for q in S['outputs'] + [S['log']]], S['stage']
            for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
            for f in P['outputs']: last[f['path']] = f['sha256']
        bad = [q for q, h in last.items() if sha(os.path.join(S4, q)) != h]; assert not bad, bad
        man['stages'] = prev['stages'][:i0]
        man['checks'] = {k: v for k, v in prev['checks'].items()
                         if k not in ('inputs_with_two_versions_in_pass_B', 'in_progress') and (cont or k not in ('rerun_from', 'segments'))}
        if not cont and 'rerun_from' in prev['checks']:        # a new pass: keep how the reused stages were produced
            man['checks']['earlier_reruns'] = man['checks'].get('earlier_reruns', []) + [dict(rerun_from=prev['checks']['rerun_from'], segments=prev['checks'].get('segments'))]
        if 'rerun_from' not in man['checks']:
            man['checks']['rerun_from'] = dict(stage=STAGES[i0]['stage'], index=i0, reused_stages=i0, reused_paths_verified=len(last))
        elif 'segments' not in man['checks']:                  # a killed pass of the earlier code: its finished part is the first segment
            man['checks']['segments'] = [[man['checks']['rerun_from']['index'], i0 - 1, 'killed after this stage']]
        if cont or 'segments' in man['checks']: man['checks'].setdefault('segments', []).append([i0, None])
        for P in man['stages']:                                # the t9 fixed point across segments: pass A's t9 output, as recorded
            if P['stage'].startswith('t9.py') and P['stage'].endswith('[A]'): t9A = next(f['sha256'] for f in P['outputs'] if f['path'] == rel(f'{O}/t9_cone_kch.json'))
        print(f"rerun from stage {i0} ({STAGES[i0]['stage']}): stages 0-{i0 - 1} reused, {len(last)} paths verified", flush=True)
    if '--until' in sys.argv: i1 = stage_index(sys.argv[sys.argv.index('--until') + 1], i0)
    for S in STAGES[i0:i1 + 1]:
        ins = [dict(path=rel(p), sha256=sha(p)) for p in S['inputs']]
        assert all(i['sha256'] != 'MISSING' for i in ins), [i['path'] for i in ins if i['sha256'] == 'MISSING']
        env = dict(os.environ, **S['env']); t0 = time.time()
        with open(S['log'], 'w') as lf:
            r = subprocess.run([sys.executable, '-I'] + S['cmd'], cwd=S['cwd'], env=env, stdout=lf, stderr=subprocess.STDOUT)
        print(f"{S['stage']}: rc {r.returncode} ({time.time() - t0:.0f}s)", flush=True)
        assert r.returncode == 0, S['stage']
        outs = [dict(path=rel(p), sha256=sha(p)) for p in S['outputs'] + [S['log']]]
        man['stages'].append(dict(stage=S['stage'] + (' [pass A]' if S.get('pass_') == 'A' and not S['stage'].endswith('[A]') else ''), pass_A=S.get('pass_') == 'A', cwd=rel(S['cwd']), cmd=['python3', '-I'] + S['cmd'], env=S['env'], inputs=ins, outputs=outs,
                                  seconds=round(time.time() - t0, 1)))
        if S['stage'].startswith('t9.py') and S['stage'].endswith('[A]'): t9A = sha(f'{O}/t9_cone_kch.json')
        if S['stage'].startswith('t9.py') and S['stage'].endswith('[B]'):
            man['checks']['t9_fixed_point'] = (t9A == sha(f'{O}/t9_cone_kch.json')); assert man['checks']['t9_fixed_point']
        if S['stage'] == 'T8 merge':
            cb = json.load(open(CB1))['args']['cconst']; man['checks']['constants_b1_residuals_sha'] = cb['residuals_sha256'] == sha(f'{T6C}/residuals.json')
            man['checks']['constants_b1_T7SLAB0'] = cb.get('env', {}).get('T7SLAB0'); assert man['checks']['constants_b1_residuals_sha']
        if S['stage'].startswith('switch'):
            sw = json.load(open(f'{O}/switch_s2.json'))['modulated']
            man['checks']['switch_vs_typed'] = dict(Y=sw['Y_raw'], Y_typed=6.21e-13, F=sw['F_raw'], F_typed=[float(v) for v in F0.split(',')],
                                                    ok=sw['Y_raw'] <= 6.21e-13 and all(a <= float(b) for a, b in zip(sw['F_raw'], F0.split(','))))
            assert man['checks']['switch_vs_typed']['ok'], man['checks']['switch_vs_typed']
        write_atomic(PARTIAL, man)                             # progress only; FINAL_manifest.json is written at the end
    if man['checks'].get('segments'): man['checks']['segments'][-1][1] = i1
    if i1 < len(STAGES) - 1:                                   # --until: a partial manifest, continued by --from i1 + 1
        man['checks']['in_progress'] = dict(next=i1 + 1, next_stage=STAGES[i1 + 1]['stage'])
        write_atomic(PARTIAL, man)
        print(f'stopped after stage {i1}; continue with --from {i1 + 1} (progress in {rel(PARTIAL)})'); sys.exit(0)
    # one version per input: every path read by several stages has one sha256 at the last pass
    seen = {}
    for S in man['stages']:
        if S['pass_A']: continue
        for i in S['inputs']: seen.setdefault(i['path'], set()).add(i['sha256'])
    man['checks']['inputs_with_two_versions_in_pass_B'] = sorted(p for p, v in seen.items() if len(v) > 1)
    man['script_sha256'] = {rel(f): sha(f) for f in [__file__, f'{P2}/t9.py', f'{P2}/t9chart.py', f'{P2}/t9chart_check.py', f'{T8}/camp.py', f'{T8}/cconst.py', f'{T8}/merge.py', f'{P2}/t9_far.py',
                                                     f'{P2}/t9_lemmaO.py', f'{P2}/t1clock.py', f'{P2}/switch.py', f'{P2}/t9h/clam.py', f'{T8}/h4loc.py', f'{T8}/n3far.py',
                                                     f'{T8}/n3boxes.py', f'{P2}/ycomp.py', f'{P2}/t9_h6.py', f'{P2}/t9_h6_v32.py', f'{P2}/t9_v32_checks.py',
                                                     f'{P2}/t9_lemmaO_far.py', f'{P2}/t9_axis_crossing.py', f'{P2}/curvjets.py', f'{P2}/t9h/basegeo.py',
                                                     f'{P2}/t9h/h1.py', f'{P2}/t9h/h1_points.py', f'{P2}/t9h/h1_summary.py']}
    names = [S['stage'] + (' [pass A]' if S.get('pass_') == 'A' and not S['stage'].endswith('[A]') else '') for S in STAGES]
    assert [P['stage'] for P in man['stages']] == names, 'incomplete stage list'
    last = {}                                                  # every recorded SHA-256 against the files now (the last write of each path)
    for P in man['stages']:
        for f in P['inputs']: last.setdefault(f['path'], f['sha256'])
        for f in P['outputs']: last[f['path']] = f['sha256']
    bad = [q for q, h in last.items() if sha(os.path.join(S4, q)) != h]
    bad += [u['path'] for u in man['upstream'] if sha(u['path'] if u['path'].startswith('/') else os.path.join(S4, u['path'])) != u['sha256']]
    assert not bad, bad
    man['checks']['hashes_verified_at_end'] = dict(paths=len(last), upstream=len(man['upstream']), scripts=len(man['script_sha256']))
    write_atomic(FINAL, man)
    if os.path.exists(PARTIAL): os.remove(PARTIAL)
    print('checks:', json.dumps(man['checks']))
