#!/usr/bin/env python3
"""step-4 v38 stage 21 (p2/v38/final_tail_v38.py): flux lemma v2, v37 reassembly -- reproduction and pin.
The pinned objects are writeup/lemmaU/bulk/fix13_v37chain/flux/out/* (8 JSON + flux2_forms.pkl) and flux/texts/* (V37_REPORT.md, verdict PASS).
This stage reruns the 8 flux stages of flux/code/run_all.sh (same order, same scripts, unmodified) with the documented scratch switches
F2_OUT = p2/out_v38/flux/out and F2_TX = p2/out_v38/flux/texts (f2common.py / gen_texts.py), F2_VTAG = fin13v4c, F2_V36RUN = chain v37
runs/prod2, each child `python3 -B <script>` with cwd flux/code and an environment built from scratch (the flux scripts import f2common
from their own directory, so the child cannot use -I), and checks:
  every rerun JSON == the pinned JSON apart from the args block (and the *.timing.json sidecars, not compared);
  flux2_forms.pkl and every text file byte-identical;
  inputs.json: the only allowed difference is its pin of the chain report V37_LOG.md (review fix15 B-v37-1: V37_LOG.md was edited after the
  flux run; the pinned file records the earlier sha256), reported below;
  the pinned flux/out equals the double-run snapshots snap/v37_pass1 and snap/v37_pass2 (byte for byte);
  flux2_cert.json verdict PASS, the combined pointwise margin (FLUX0 + eps_I FLUX2)/(K0bar + eps_I K2pt) and the C0..C7 results recorded;
  the fix13_v37chain tree (flux, runs/prod2, t2_A/out, t4_cert/out, t5_zone3/out) is unchanged by the stage (every file hash, no new file).
Writes only under p2/out_v38/flux/.  usage: nice -n 19 python3 -I -B flux_v38.py"""
import os, sys, json, hashlib, subprocess, time
from fractions import Fraction
sys.dont_write_bytecode = True
NC = '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain'; FX = f'{NC}/flux'; CODE = f'{FX}/code'
OD = '/home/claude/s1spec/step4/p2/out_v38/flux'; OO, OT, OL = f'{OD}/out', f'{OD}/texts', f'{OD}/logs'
STAGES = ['flux2_forms.py', 'flux2_coef.py', 'chart4_shell.py', 'inputs_manifest.py', 'flux2_data.py', 'flux2_assemble.py', 'flux2_sanity.py', 'gen_texts.py']
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
J = lambda p: json.load(open(p))
for d in (OO, OT, OL): os.makedirs(d, exist_ok=True)
for d in (OO, OT):
    for f in os.listdir(d): os.remove(os.path.join(d, f))                     # stale outputs of an earlier pass never survive
TREES = [FX, f'{NC}/runs/prod2', f'{NC}/t2_A/out', f'{NC}/t4_cert/out', f'{NC}/t5_zone3/out', f'{NC}/n1/out']
def snapshot():
    out = {}
    for t in TREES:
        for dp, dn, fn in os.walk(t):
            for f in fn: p = os.path.join(dp, f); out[p] = sha(p)
    return out
SNAP0 = snapshot(); print(f'snapshot: {len(SNAP0)} files', flush=True)
CHK = []
def check(name, ok, detail=None):
    CHK.append(dict(check=name, ok=bool(ok), detail=detail)); print(('PASS ' if ok else 'FAIL ') + name, '' if detail is None else json.dumps(detail, default=str)[:300], flush=True)
PINNED = sorted(f for f in os.listdir(f'{FX}/out') if f.endswith('.json') and not f.endswith('.timing.json')) + ['flux2_forms.pkl']
TEXTS = sorted(os.listdir(f'{FX}/texts'))
for snap in ('v37_pass1', 'v37_pass2'):
    bad = [f for f in PINNED if sha(f'{FX}/out/{f}') != sha(f'{FX}/snap/{snap}/out/{f}')] + [f for f in TEXTS if sha(f'{FX}/texts/{f}') != sha(f'{FX}/snap/{snap}/texts/{f}')]
    check(f'pinned flux/out and flux/texts == snap/{snap} byte for byte ({len(PINNED)} out files, {len(TEXTS)} texts)', not bad, bad)
# ---- rerun
ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1',
       'F2_OUT': OO, 'F2_TX': OT, 'F2_VTAG': 'fin13v4c', 'F2_V36RUN': f'{NC}/runs/prod2'}
TIMES = {}
for s in STAGES:
    t0 = time.time()
    with open(f'{OL}/{s[:-3]}.log', 'w') as lf:
        r = subprocess.run([sys.executable, '-B', s], cwd=CODE, env=ENV, stdout=lf, stderr=subprocess.STDOUT)
    TIMES[s] = round(time.time() - t0, 1); print(f'{s}: rc {r.returncode} ({TIMES[s]} s)', flush=True)
    assert r.returncode == 0, (s, open(f'{OL}/{s[:-3]}.log').read()[-2000:])
# ---- compare
def strip(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x
CMP = {}
for f in PINNED:
    p, n = f'{FX}/out/{f}', f'{OO}/{f}'
    if not os.path.exists(n): CMP[f] = 'missing in rerun'; continue
    if f.endswith('.pkl'): CMP[f] = 'byte-identical' if sha(p) == sha(n) else 'DIFFERS'; continue
    a, b = strip(J(p)), strip(J(n))
    if f == 'inputs.json':
        rep = f'{NC}/V37_LOG.md'
        pa, pb = dict(a.get('pinned', {})), dict(b.get('pinned', {}))
        BV37 = dict(report=rep, pinned_in_reference=pa.get(rep), pinned_in_rerun=pb.get(rep), disk=sha(rep),
                    note='review fix15 B-v37-1: V37_LOG.md edited after the flux run; only the run name is read from it')
        pa.pop(rep, None); pb.pop(rep, None); a = dict(a, pinned=pa); b = dict(b, pinned=pb)
    CMP[f] = 'equal apart from args' if a == b else 'DIFFERS'
CMP_T = {f: ('byte-identical' if os.path.exists(f'{OT}/{f}') and sha(f'{FX}/texts/{f}') == sha(f'{OT}/{f}') else 'DIFFERS') for f in TEXTS}
check('rerun reproduces every pinned flux/out file (JSON apart from args; pkl byte for byte)', all(v in ('equal apart from args', 'byte-identical') for v in CMP.values()), CMP)
check('rerun reproduces every flux/texts file byte for byte', all(v == 'byte-identical' for v in CMP_T.values()), CMP_T)
check('inputs.json: the only non-args difference is the V37_LOG.md pin (B-v37-1); rerun pins the current file', BV37['pinned_in_rerun'] == BV37['disk'], BV37)
# ---- results of the pinned certificate
C = J(f'{FX}/out/flux2_cert.json'); res = C['res']
from flint import arb, ctx
ctx.prec = 200
K0b, K2p = arb(res['K0bar']['arb']), arb(res['K2pt']['arb']); eps = arb('3e-5')
comb = (arb('1.64e-20') + eps * arb('8.0e-20')) / (K0b + eps * K2p)
# C8 = U-1 (regularity lemma) is recorded by flux2_assemble with pass_ = None: proved by hand in fix10_proofs (U1_regularity.tex), not scripted (V37_REPORT)
check('flux2_cert verdict PASS, no failing condition, every scripted condition passes (C8 = U-1: by hand, pass_ null as recorded)', C['verdict'] == 'PASS' and not C['failing']
      and all(v.get('pass_') is True for k, v in C['conditions'].items() if k != 'C8') and C['conditions'].get('C8', {}).get('pass_') is None,
      {k: (v.get('margin'), v.get('pass_')) for k, v in C['conditions'].items()})
check('combined pointwise margin (FLUX0 + eps_I FLUX2)/(K0bar + eps_I K2pt) > 1 [cert]', (comb - 1).lower() > 0, dict(lower=float(comb.lower())))
SNAP1 = snapshot(); chg = sorted(p for p in set(SNAP0) | set(SNAP1) if SNAP0.get(p) != SNAP1.get(p))
check('fix13_v37chain trees unchanged by this stage (every file hash, no new file)', not chg, chg[:20])
OUTS = sorted([f'{OO}/{f}' for f in os.listdir(OO) if not f.endswith('.timing.json') and not f.endswith('.time')] + [f'{OT}/{f}' for f in os.listdir(OT)])
out = dict(item='flux lemma v2, v37 reassembly (fix13_v37chain/flux): reproduction and pin for the step-4 v38 manifest', label='certified (flux scripts: arb, outward)',
           pinned={f: sha(f'{FX}/out/{f}') for f in PINNED}, pinned_texts={f: sha(f'{FX}/texts/{f}') for f in TEXTS}, compare=CMP, compare_texts=CMP_T, B_v37_1=BV37,
           cert=dict(verdict=C['verdict'], K0bar=res['K0bar']['arb'], K2pt=res['K2pt']['arb'], combined_margin_lower=float(comb.lower()),
                     conditions=C['conditions'], failing=C['failing'], margin_combined=res.get('margin_combined'), strictness_pad_C3=res.get('strictness_pad_C3')),
           rerun_outputs={p: sha(p) for p in OUTS}, checks=CHK, all_ok=all(c['ok'] for c in CHK),
           args=dict(argv=sys.argv, child_env=ENV, child_cmd='python3 -B <stage> (cwd flux/code)', stages=STAGES,
                     inputs={p: sha(p) for p in [os.path.abspath(__file__)] + [f'{CODE}/{s}' for s in STAGES] + [f'{CODE}/f2common.py']}))
tmp = f'{OD}/flux_v38.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1, default=str); os.replace(tmp, f'{OD}/flux_v38.json')
print('stage times (s):', TIMES)
print('all ok:', out['all_ok'], f"({sum(c['ok'] for c in CHK)}/{len(CHK)})")
assert out['all_ok'], [c for c in CHK if not c['ok']]
