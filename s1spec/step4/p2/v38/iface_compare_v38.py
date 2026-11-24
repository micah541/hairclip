#!/usr/bin/env python3
"""step-4 v38 stage 27 (p2/v38/final_tail_v38.py): interface and proof certificates -- comparison of the v38 reruns with the reference outputs.
Stages 22-26 reran (as v38 copies that differ only in the output path, the timing field, and for G1 the input run) :
  C1_ke      (fix14_interface/C1_ke.py; reference fix15_review/c/C1_ke.json, also fix14_interface/C1_ke.json)        sup|k~ e| < 1/2 on P2's boxes
  GS_tau1    (fix14_interface/GS_tau1.py; reference fix15_review/c/GS_tau1.json, fix14_interface/GS_tau1.json)       (GS) at tau_1
  GS_boot    (fix15_review/c/GS_boot_v4c.py on fin13v4c production; reference fix15_review/c/GS_boot_v4c.json)       (GS) strict improvement
  G1_check   (fix10_proofs/G1_check.py, repointed from chain v36 to chain v37 + flux v37; reference fix10_proofs/G1_check.json = v36)
  U1_sym     (fix10_proofs/U1_sym.py; reference fix10_proofs/U1_sym.json)                                               exact sympy identities
Checks: every rerun has ok/all_pass true; C1_ke, GS_tau1, GS_boot, U1_sym equal their references apart from the args block; G1_check passes
every check on v37 (its values are v37 values, so it differs from the v36 reference: the changed values are listed); each reference output's
args pins of its inputs are compared with the disk (reported).  Writes p2/out_v38/iface/iface_compare_v38.json.
usage: python3 -I -B iface_compare_v38.py"""
import os, sys, json, hashlib
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
J = lambda p: json.load(open(p))
B = '/home/claude/s1spec/writeup/lemmaU/bulk'; OI = '/home/claude/s1spec/step4/p2/out_v38/iface'
REF = dict(C1_ke=[f'{B}/fix15_review/c/C1_ke.json', f'{B}/fix14_interface/C1_ke.json'], GS_tau1=[f'{B}/fix15_review/c/GS_tau1.json', f'{B}/fix14_interface/GS_tau1.json'],
           GS_boot_v4c=[f'{B}/fix15_review/c/GS_boot_v4c.json'], G1_check=[f'{B}/fix10_proofs/G1_check.json'], U1_sym=[f'{B}/fix10_proofs/U1_sym.json'])
def strip(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x
def pins(x):
    a = x.get('args', {}); inp = a.get('inputs', {})
    if isinstance(inp, list): inp = {d['path']: d['sha256'] for d in inp}
    return {p: h for p, h in inp.items() if isinstance(p, str) and p.startswith('/')}
CHK = []; ROWS = {}
def check(name, ok, detail=None):
    CHK.append(dict(check=name, ok=bool(ok), detail=detail)); print(('PASS ' if ok else 'FAIL ') + name, '' if detail is None else json.dumps(detail, default=str)[:300], flush=True)
for n, refs in REF.items():
    new = J(f'{OI}/{n}.json'); okk = new.get('ok', new.get('all_pass'))
    check(f'{n}: rerun ok / all_pass is true', okk is True)
    row = dict(rerun=f'{OI}/{n}.json', rerun_sha256=sha(f'{OI}/{n}.json'), refs={})
    for r in refs:
        R = J(r); pr = pins(R); badp = {p: h for p, h in pr.items() if not os.path.exists(p) or sha(p) != h}
        row['refs'][r] = dict(sha256=sha(r), equal_apart_from_args=strip(R) == strip(new), ref_input_pins=len(pr), ref_pins_not_matching_disk=sorted(badp))
    pn = pins(new); badn = {p: h for p, h in pn.items() if not os.path.exists(p) or sha(p) != h}
    check(f'{n}: every input pin of the rerun matches the disk ({len(pn)} pins)', not badn, sorted(badn))
    if n == 'G1_check':
        ref = J(refs[0])
        row['v36_to_v37_changed_values'] = {k: dict(v36=ref['values'].get(k), v37=v) for k, v in new['values'].items() if ref['values'].get(k) != v}
        row['checks_v37'] = len(new['checks']) if isinstance(new.get('checks'), (list, dict)) else None
        check('G1_check on chain v37: all checks pass (values are v37, differ from the v36 reference by design)', new['all_pass'] is True,
              dict(n_changed=len(row['v36_to_v37_changed_values'])))
        check('G1_check: the rerun reads chain v37 runs/prod2 and flux v37 flux2_data.json only', all(('fix13_v37chain' in p) for p in pn if '/writeup/lemmaU/' in p and not p.endswith('_v38.py')), sorted(pn))
    else:
        check(f'{n}: rerun equals the primary reference {os.path.relpath(refs[0], B)} apart from args', row['refs'][refs[0]]['equal_apart_from_args'])
    ROWS[n] = row
out = dict(item='interface (C-1, GS) and proof (G1, U-1) certificates for the step-4 v38 manifest: reruns vs references', label='comparison (the certificates are [cert], arb / exact sympy)',
           rows=ROWS, checks=CHK, all_ok=all(c['ok'] for c in CHK),
           args=dict(argv=sys.argv, env={}, inputs={p: sha(p) for p in [os.path.abspath(__file__)] + [f'{OI}/{n}.json' for n in REF] + [r for v in REF.values() for r in v]}))
tmp = f'{OI}/iface_compare_v38.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, f'{OI}/iface_compare_v38.json')
print('all ok:', out['all_ok'], f"({sum(c['ok'] for c in CHK)}/{len(CHK)})")
assert out['all_ok'], [c for c in CHK if not c['ok']]
