"""Final hand-off tail on the B1 data, one version of every input (coordinator, Nov 12).  Writes p2/out/FINAL_manifest.json.
Data sets: T6 = p2/t6c (residuals.json, jets/; the (SE)-corrected base, T6's engine, B1); T7 = t7/constants.json (current);
T8 = constants_b1.json (cconst + merge on p2/t6c and T7 slabs >= 2, rebuilt here); base windows T9WIN = out/h1slab/win_summary.json.
Superseded and read by no stage: p2/t6/*, step4/t6/residuals.json, se_corr.json, constants_fin.json, n3.json, amp_fin*.
Every camp run reads the final t9_cone_kch.json: the T8 iteration and the window camp run twice (passes A, B); pass B's t9 must
reproduce pass A's byte for byte (t9 reads only the amp's per-slab aF, bF and at_Ybar values, which do not depend on T9JSON).
usage: python3 -I final_tail.py"""
import os, sys, json, hashlib, subprocess, time, glob
S4 = '/home/claude/s1spec/step4'; P2 = f'{S4}/p2'; T8 = f'{P2}/t8se'; O = f'{P2}/out'; T6C = f'{P2}/t6c'
F0 = '4e-23,3e-22,5e-19,3.4e-15'; YB = '1e-12'
S1ENV = dict(COHERENT='1', S1PROF='/home/claude/s1_tight/ext_tip/prof_tip.txt', S1TT='10.249077651120778,10.249077651121361')
CAMPENV = dict(T9JSON='../out/t9_cone_kch.json', T9WIN='../out/h1slab/win_summary.json', T6DIR='../t6c', OMP_NUM_THREADS='1')
JETS = lambda *sl: [f'{T6C}/jets/jets_{s:03d}.npz' for s in sl]
T7C = f'{S4}/t7/constants.json'; CB1 = f'{T8}/constants/constants_b1.json'; AMP = f'{T8}/constants/amp_A1_w180.json'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
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
    st('clam (error transfer; mu_e from the B1 constants)', f'{P2}/t9h', ['clam.py'], dict(S1ENV, T8CONSTS=CB1),
       [CB1, f'{O}/h1slab/summary.json'] + [f'{O}/t9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json' for k in range(16)], [f'{O}/t9_h1_clam_s180.json'], f'{O}/clam_final.log'),
] + [dict(x, pass_='A') for x in t8_iteration() + window_and_t9('A')] + t8_iteration() + window_and_t9('B') + [
    st('t9_far.py', P2, ['t9_far.py', 't8se/constants/amp_A1_w180.json', 't8se/constants/constants_b1.json'], dict(T6DIR=T6C), [AMP, CB1] + JETS(180), [f'{O}/t9_far.json'], f'{O}/t9_far_final.log'),
    st('h4loc on the final amp', T8, ['h4loc.py', 'constants/amp_A1_w180.json'], dict(T6DIR='../t6c', H4OUT='../out/h4loc_A1_w180.json'),
       [AMP, f'{T6C}/residuals.json'] + JETS(180), [f'{O}/h4loc_A1_w180.json'], f'{O}/h4loc_final.log'),
    st('Lemma O H1/H8 (t9_lemmaO.py, LO_NS=48)', P2, ['t9_lemmaO.py', 't8se/constants/amp_A1_w180.json', '180'], dict(T6DIR=T6C, LO_NS='48'), [AMP] + JETS(180),
       [f'{O}/t9_lemmaO.json'], f'{O}/t9_lemmaO_final.log'),
    st('Lemma O H6 (t9_h6.py)', P2, ['t9_h6.py'], {}, [f'{O}/t9_cone_kch.json', f'{O}/t9_far.json', f'{S4}/../writeup/lemmaO/lemmaO_v3_consts.json'], [f'{O}/t9_h6.json'], f'{O}/t9_h6_final.log'),
    st('v3.2: t9_far.py --v32', P2, ['t9_far.py', 't8se/constants/amp_A1_w180.json', 't8se/constants/constants_b1.json', '--v32'], dict(T6DIR=T6C), [AMP, CB1] + JETS(180),
       [f'{O}/t9_far_v32.json'], f'{O}/t9_far_v32_final.log'),
    st('v3.2: h4loc FAR_M=0.15', T8, ['h4loc.py', 'constants/amp_A1_w180.json'], dict(T6DIR='../t6c', FAR_M='0.15', H4OUT='constants/h4loc_v32.json'),
       [AMP, f'{T6C}/residuals.json'] + JETS(180), [f'{T8}/constants/h4loc_v32.json'], f'{O}/h4loc_v32_final.log'),
    st('v3.2: t9_v32_checks.py', P2, ['t9_v32_checks.py'], {}, [f'{O}/t9_cone_kch.json'], [f'{O}/t9_v32_checks.json'], f'{O}/t9_v32_checks_final.log'),
    st('v3.2: t9_h6_v32.py', P2, ['t9_h6_v32.py'], {}, [f'{O}/t9_cone_kch.json', f'{O}/t9_far_v32.json', f'{T8}/constants/h4loc_v32.json', f'{P2}/h6/fable/h6f_cases.json',
       f'{S4}/../writeup/lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json'], [f'{O}/t9_h6_v32.json'], f'{O}/t9_h6_v32_final.log'),
    st('v3.2: t9_lemmaO.py --v32', P2, ['t9_lemmaO.py', 't8se/constants/amp_A1_w180.json', '180', '--v32'], dict(T6DIR=T6C, LO_NS='48'), [AMP] + JETS(180),
       [f'{O}/t9_lemmaO_v32.json'], f'{O}/t9_lemmaO_v32_final.log'),
    st('t1clock', P2, ['t1clock.py', 't8se/constants/constants_b1.json', '180', YB], dict(T6DIR=T6C), [CB1, f'{T6C}/residuals.json'] + JETS(180), [f'{O}/t1clock.json'],
       f'{O}/t1clock_final.log'),
]

if __name__ == '__main__':
    man = dict(note=__doc__, data=dict(T6=rel(T6C), T7=rel(T7C), T8=rel(CB1), T9WIN='p2/out/h1slab/win_summary.json'),
               upstream=[], stages=[], checks={})
    A1 = {f['path']: f for f in json.load(open(f'{O}/A1_manifest.json'))['files']}
    for p, f in A1.items():                                    # unchanged upstream: S1 hulls and the H1 chain (base polynomials only)
        if any(k in f['stage'] for k in ('S1 hull', 'Krawczyk', 'h1dt', 'H1 slab', 'win_slab', 'win.py', 'windows/shells')):
            pp = p if p.startswith('/') else os.path.join(S4, p); assert sha(pp) == f['sha256'], p
            man['upstream'].append(dict(stage=f['stage'], path=p, sha256=f['sha256']))
    t9A = None
    for S in STAGES:
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
        json.dump(man, open(f'{O}/FINAL_manifest.json', 'w'), indent=1)
    # one version per input: every path read by several stages has one sha256 at the last pass
    seen = {}
    for S in man['stages']:
        if S['pass_A']: continue
        for i in S['inputs']: seen.setdefault(i['path'], set()).add(i['sha256'])
    man['checks']['inputs_with_two_versions_in_pass_B'] = sorted(p for p, v in seen.items() if len(v) > 1)
    man['script_sha256'] = {rel(f): sha(f) for f in [__file__, f'{P2}/t9.py', f'{P2}/t9chart.py', f'{P2}/t9chart_check.py', f'{T8}/camp.py', f'{T8}/cconst.py', f'{T8}/merge.py', f'{P2}/t9_far.py',
                                                     f'{P2}/t9_lemmaO.py', f'{P2}/t1clock.py', f'{P2}/switch.py', f'{P2}/t9h/clam.py', f'{T8}/h4loc.py', f'{T8}/n3far.py',
                                                     f'{T8}/n3boxes.py', f'{P2}/ycomp.py', f'{P2}/t9_h6.py', f'{P2}/t9_h6_v32.py', f'{P2}/t9_v32_checks.py']}
    json.dump(man, open(f'{O}/FINAL_manifest.json', 'w'), indent=1)
    print('checks:', json.dumps(man['checks']))
