"""step-4 v38 stage 28: the input pins of the verifier fin13v4c runs, with the documented exception p2/v38/EXCEPTION_v4c_PA_reg_e2v4b.json.
runs/P (pinned by the v38 manifest): every pin of its two rate tables and of mod_sups must match the disk (no exception applies).
runs/P_A, runs/reg_fin13v4b (NOT pinned): a mismatching pin is accepted only if it is the exception record's path with the recorded pinned sha256 and the file
now has the recorded current sha256; then their numbers are compared (rate tables and mod_sups apart from args and timing fields runtime_s, secs): P_A == P, reg_fin13v4b ==
step3nl_v4b/runs/P.  Read-only; writes p2/out_v38/v4c_runs_pins.json.   usage: python3 -I -B verify_v4c_runs_with_exceptions.py"""
import os, sys, json, glob, hashlib
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else 'MISSING'
HERE = '/home/claude/s1spec/step4/p2/v38'; R = '/home/claude/s1spec/step3nl_v4c/runs'; RB = '/home/claude/s1spec/step3nl_v4b/runs/P'
EXC = f'{HERE}/EXCEPTION_v4c_PA_reg_e2v4b.json'; E = json.load(open(EXC))
def recs(d): return sorted(glob.glob(f'{d}/rate_table_*.json') + glob.glob(f'{d}/modsups_*/final_alpha_points_5e-8.json'))
def strip(x):
    if isinstance(x, dict): return {k: strip(v) for k, v in x.items() if k not in ('args', 'runtime_s', 'secs')}
    if isinstance(x, list): return [strip(v) for v in x]
    return x
out = dict(item='pins of the verifier fin13v4c runs (runs/P pinned; runs/P_A, runs/reg_fin13v4b with the documented exception)', exception=dict(path=EXC, sha256=sha(EXC)), runs={}, checks=[])
def check(n, ok, d=None): out['checks'].append(dict(check=n, ok=bool(ok), detail=d)); print(('PASS ' if ok else 'FAIL ') + n, '' if d is None else json.dumps(d)[:300])
for run in ('P', 'P_A', 'reg_fin13v4b'):
    rows = {}
    for f in recs(f'{R}/{run}'):
        pins = json.load(open(f))['args']['inputs']; bad = {p: dict(pinned=h, disk=sha(p)) for p, h in pins.items() if sha(p) != h}
        acc = {p: v for p, v in bad.items() if run != 'P' and p == E['path'] and v['pinned'] == E['pinned_sha256'] and v['disk'] == E['current_sha256'] and 'runs/' + os.path.relpath(f, R) in E['affected_records']}
        rows[os.path.relpath(f, R)] = dict(sha256=sha(f), pins=len(pins), mismatches=bad, accepted_by_exception=sorted(acc), other=sorted(set(bad) - set(acc)))
    out['runs'][run] = rows
    check(f'runs/{run}: every pin matches the disk' + (' except the documented exception' if run != 'P' else ''), rows and all(not r['other'] for r in rows.values()),
          {k: len(v['accepted_by_exception']) for k, v in rows.items()})
check('runs/P: no pin needs the exception', all(not r['mismatches'] for r in out['runs']['P'].values()))
for a, b, lab in (('P_A', f'{R}/P', 'runs/P_A == runs/P'), ('reg_fin13v4b', RB, 'runs/reg_fin13v4b == step3nl_v4b/runs/P')):
    pairs = [(f, os.path.join(b, os.path.relpath(f, f'{R}/{a}'))) for f in recs(f'{R}/{a}')]
    eq = {os.path.relpath(f, R): os.path.exists(g) and strip(json.load(open(f))) == strip(json.load(open(g))) for f, g in pairs}
    check(f'numbers: {lab} (rate tables and mod_sups apart from args and timing: runtime_s, oracle_rows secs)', pairs and all(eq.values()), eq)
out['all_ok'] = all(c['ok'] for c in out['checks'])
# runs/P_A and runs/reg_fin13v4b are NOT pinned: their record hashes are listed under 'runs' only, not as inputs of this stage
out['args'] = dict(argv=sys.argv, env={}, inputs={p: sha(p) for p in [os.path.abspath(__file__), EXC] + recs(f'{R}/P')}, not_pinned_read=sorted(f for run in ('P_A', 'reg_fin13v4b') for f in recs(f'{R}/{run}')) + recs(RB))
tmp = '/home/claude/s1spec/step4/p2/out_v38/v4c_runs_pins.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, tmp[:-4])
print('all ok:', out['all_ok']); assert out['all_ok']
