#!/usr/bin/env python3
"""Stage X0 (FLUX2_DESIGN.md sec. 10): pin every input of flux lemma v2 by sha256.  Refuses to run if fix6_v36/V36_REPORT.md is missing.
The final v36 snapshot directory is given by env F2_V36RUN (a directory with t2_A/out-like files, e.g. fix6_v36/runs/prod1, or fix6_v36
itself, whose outputs live in t2_A/out and n1/out); it must be the one V36_REPORT.md names as final (the path string is asserted to occur
in the report).  Checks: rate-table env (FLUX0, FLUX2, S0V, S1V) = (1.64e-20, 8.0e-20, 0.0153, 0.6); v36 asm/A_interface use the
fin13v4b rate table (sha equal to the pinned one); the level-2 export block exists.  Output out/inputs.json."""
import sys, os, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
T = FC.Timer()
REP = FC.CHAIN_REPORT
DRY = os.environ.get('F2_DRYRUN') == '1' and 'F2_OUT' in os.environ      # dry run (scratch output only): report not required
assert DRY or os.path.exists(REP), 'V36_REPORT.md missing: X0 refuses to run'
run = os.environ.get('F2_V36RUN')
assert run and os.path.isdir(run), 'set F2_V36RUN to the final v36 snapshot directory'
if not DRY:
    rep = open(REP).read()
    rel = os.path.relpath(run, FC.N36)
    assert (run in rep) or (rel in rep) or ('runs/' + os.path.basename(run) in rep), 'F2_V36RUN not named in V36_REPORT.md'


def pick(*cands):
    for c in cands:
        p = os.path.join(run, c)
        if os.path.exists(p): return p
    raise FileNotFoundError(cands)


v36 = dict(asm_design_N2ref4=pick('t2_A/asm_design_N2ref4.json', 't2_A/out/asm_design_N2ref4.json'),
           A_interface=pick('t2_A/A_interface.json', 't2_A/out/A_interface.json'),
           hjw_rep=pick('t2_A/hjw_rep.json', 't2_A/out/hjw_rep.json'),
           n1_chart=pick('n1/n1_chart_%s.json' % FC.VTAG, 'n1/out/n1_chart_%s.json' % FC.VTAG),
           lemmaS_cflux=os.path.join(FC.N36, 't2_A', 'base_ibp', 'lemmaS_cflux.json'),
           v36cfg=pick('v36cfg.json') if os.path.exists(os.path.join(run, 'v36cfg.json')) else os.path.join(FC.N36, 'v36cfg.py'),
           **({} if DRY else dict(report=REP)))
RT = {a: os.path.join(FC.V4B, 'runs', 'P', 'rate_table_%s_a%s.json' % (FC.VTAG, a)) for a in ('0.00015', '0.000125')}
proof = [os.path.join(FC.V4B, 'v4cert', 'out', 'constants_v4.json'), os.path.join(FC.V4B, 'v4cert', 'out', 'c2_t3d.json'), os.path.join(FC.V4B, 'v4cert', 'c2_t3d.py'),
         os.path.join(FC.V4B, 'graph_struct.pkl'), os.path.join(FC.V4B, 'verifier', 'closure.py'), os.path.join(FC.V4B, 'verifier', 'constants.py'),
         os.path.join(FC.V4B, 'verifier', 'rhs_iv.py'), os.path.join(FC.V4B, 'verifier', 'run_v4b.sh' if FC.VTAG == 'fin13v4b' else 'run_v4c.sh'),
         '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/lin/lintail2.py', '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py',
         '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py', '/home/claude/s1spec/step4/p2/t9h/s1chart.py',
         os.path.join(FC.F1, 'flux_forms.pkl'), os.path.join(FC.F1, 'flux_cert.json'), '/home/claude/s1spec/writeup/lemmaS/lemmaS_cflux.json',
         os.path.join(FC.DESIGN, 'FLUX2_DESIGN.md'), os.path.join(FC.DESIGN, 'flux2_lemma.tex')]
CHK = {}


def ok(n, v):
    CHK[n] = bool(v); print('%-78s %s' % (n, bool(v)), flush=True); assert v, n


for a_, p in RT.items():
    e = json.load(open(p))['args']['env']
    ok('rate table %s: env FLUX0/FLUX2/S0V/S1V as expected' % a_, (e['FLUX0'], e['FLUX2'], e['S0V'], e['S1V']) == ('1.64e-20', '8.0e-20', '0.0153', '0.6'))
asm = json.load(open(v36['asm_design_N2ref4']))
ok('v36 asm has the level-2 export block export["2"] with C0, C1', 'C0' in asm['res']['tight']['export']['2'] and 'C1' in asm['res']['tight']['export']['2'])
ins = asm['args']['inputs'] if 'inputs' in asm['args'] else {}
rt_sha = FC.sha(RT['0.00015'])
ok('v36 asm args pin the %s rate table with the same sha' % FC.VTAG, any(os.path.basename(k) == os.path.basename(RT['0.00015']) and v == rt_sha for k, v in ins.items()))
hj = json.load(open(v36['hjw_rep']))
ok('v36 hjw_rep reads the %s rate table (same sha)' % FC.VTAG, hj.get('verifier_sha256') == rt_sha)
pinned = {p: FC.sha(p) for p in list(v36.values()) + list(RT.values()) + proof}
out = dict(stage='X0 inputs_manifest', v36_run=run, v36=v36, rate_tables=RT, pinned=pinned, checks=CHK, all_pass=all(CHK.values()),
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py')], dict(F2_V36RUN=run)))
out['timing'] = dict(runtime_s=T.s())
h = FC.dump(os.path.join(FC.OUT, 'inputs.json'), out)
print('X0 done: %d inputs pinned; sha %s' % (len(pinned), h))
