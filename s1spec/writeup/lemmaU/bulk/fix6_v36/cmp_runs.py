#!/usr/bin/env python3
"""fix6_v36: regression comparison of a v36 run snapshot against v35 (pinned outputs) or a fix4_r3 run.  Read-only.
usage: python3 -I -B cmp_runs.py <snapshot dir> v35|r3run2|<other v36 snapshot dir> [--json OUT]
Ignored: 'args' blocks and timing keys (as fix4_r3/cmpjson.py).  Remaining differences are classified:
  PATH  = both values are strings and differ only by directory prefixes / file names of the chain copies (fix6_v36 vs t2_A_v35/fix4_r3, v36 vs v35/r3 names);
  OTHER = anything else (numbers, structure, other strings).
A regression PASSES iff no OTHER difference exists in any pair."""
import sys, json, os, re
B = '/home/claude/s1spec/writeup/lemmaU/bulk'
IGN = {'args', 'cpu', 'seconds', 'wall', 'wall_s', 'elapsed', 'time', 'runtime_s'}
T2 = ['A_interface', 'a3_signs', 'asm_design_N20', 'asm_design_N2ref', 'asm_design_N2ref4', 'asm_design_pm3ibp', 'asm_hi_N2ref', 'asm_lo_N2ref', 'asm_probe_N2ref',
      'hjw_rep', 'ic_design', 'ic_hi', 'ic_lo', 'ic_probe', 'outer_L2', 'ramp_rep', 'ramp_rep_pm3ibp', 'weights_t2', 'weights_t2_fb']
BASE = ['lemmaS_F10p_M0.52.json', 'rigor/cap_F10p.json', 'rigor/lemmaS_rigcap_M0.52.json', 'rigor/weights_exact.json', 'rigor/assembly_arb.json']

def pairs(snap, ref):
    P = []
    if ref == 'v35':
        o2 = B + '/fix2/n12/t2_A_v35/out'; o4 = B + '/t4_cert/out'; o5 = B + '/t5_zone3/out'; ob = B + '/fix2/n12/t2_A_v35/base_ibp'
        for n in T2: P.append((f'{snap}/t2_A/{n}.json', f'{o2}/{n}.json'))
        for n in ('u1_ghat', 'u2_K', 'u34_primary'): P.append((f'{snap}/t4_cert/{n}.json', f'{o4}/{n}.json'))
        for k in (1, 2, 3): P.append((f'{snap}/t4_cert/U_cert_v36_pass{k}.json', f'{o4}/U_cert_v35_pass{k}.json'))
        P.append((f'{snap}/t4_cert/U_cert_for_zone3_v36.json', f'{o4}/U_cert_for_zone3_v35.json'))
        P.append((f'{snap}/t5_zone3/zone3_interface.json', f'{o5}/zone3_interface.json'))
        for n in BASE: P.append((f'{snap}/base_ibp/{n}', f'{ob}/{n}'))
    else:
        r = ref if ref.startswith('/') else B + '/fix4_r3/runs/' + ref
        nm = (lambda s: s.replace('v36', 'r3')) if '/fix4_r3/' in r else (lambda s: s)
        for n in T2: P.append((f'{snap}/t2_A/{n}.json', f'{r}/t2_A/{n}.json'))
        for n in ('u1_ghat', 'u2_K', 'u34_primary', 'U_cert_v36_pass1', 'U_cert_v36_pass2', 'U_cert_v36_pass3', 'U_cert_for_zone3_v36', 'U_cert_for_zone3_v36_pass1'):
            P.append((f'{snap}/t4_cert/{n}.json', f'{r}/t4_cert/{nm(n)}.json'))
        for n in ('zone3_interface', 'zone3_interface_pass1'): P.append((f'{snap}/t5_zone3/{n}.json', f'{r}/t5_zone3/{n}.json'))
        if '/fix4_r3/' in r: P.append((f'{snap}/r3/r3_forcing.json', B + '/fix4_r3/r3/r3_forcing.json'))
        else:
            for n in BASE: P.append((f'{snap}/base_ibp/{n}', f'{r}/base_ibp/{n}'))
            for q in ('r3/r3_forcing.json', 'r3/r3_symcheck.json', 'n1/n1_chart_fin13v4b.json'):
                if os.path.exists(f'{snap}/{q}'): P.append((f'{snap}/{q}', f'{r}/{q}'))
    return P

NORM = [(re.escape(B) + r'/fix6_v36/t2_A', 'T2'), (re.escape(B) + r'/fix2/n12/t2_A_v35', 'T2'), (re.escape(B) + r'/fix4_r3/t2_A', 'T2'),
        (re.escape(B) + r'/fix6_v36/t4_cert', 'T4'), (re.escape(B) + r'/fix4_r3/t4_cert', 'T4'), (re.escape(B) + r'/t4_cert', 'T4'),
        (re.escape(B) + r'/fix6_v36/t5_zone3', 'T5'), (re.escape(B) + r'/fix4_r3/t5_zone3', 'T5'), (re.escape(B) + r'/t5_zone3', 'T5'),
        (re.escape(B) + r'/fix6_v36/r3', 'R3'), (re.escape(B) + r'/fix4_r3/r3', 'R3'),
        (r'\.\./t4_cert', 'T4'), (r'U_cert(_for_zone3)?_(v36|v35|r3)', r'U_cert\1_X'), (r'_v36_pass', '_X_pass'), (r'U_cert_X_pass2\.json', 'U_cert_X.json')]   # v35 named its pass-2 Lemma U output U_cert_v35.json
def norm(s):
    for a, b in NORM: s = re.sub(a, b, s)
    return s

def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k in IGN: continue
            if k not in a or k not in b: out.append(('OTHER', path + '/' + str(k), 'missing in ' + ('A' if k not in a else 'B'), '')); continue
            walk(a[k], b[k], path + '/' + str(k), out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append(('OTHER', path, 'len %d' % len(a), 'len %d' % len(b))); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, '%s[%d]' % (path, i), out)
    elif a != b:
        kind = 'PATH' if isinstance(a, str) and isinstance(b, str) and norm(a) == norm(b) else 'OTHER'
        out.append((kind, path, a, b))

if __name__ == '__main__':
    snap, ref = sys.argv[1], sys.argv[2]; res = []; nother = 0
    for A, Bf in pairs(snap, ref):
        if not os.path.exists(A) or not os.path.exists(Bf): res.append(dict(a=A, b=Bf, status='MISSING')); nother += 1; continue
        out = []; walk(json.load(open(A)), json.load(open(Bf)), '', out)
        np_ = sum(1 for o in out if o[0] == 'PATH'); no = sum(1 for o in out if o[0] == 'OTHER'); nother += bool(no)
        st = 'IDENTICAL' if not out else ('PATH-ONLY (%d)' % np_ if not no else 'DIFF other=%d path=%d' % (no, np_))
        res.append(dict(a=A, b=Bf, status=st, other=[o[1:] for o in out if o[0] == 'OTHER'][:20], path=[o[1] for o in out if o[0] == 'PATH'][:20]))
        print('%-22s %s' % (st, os.path.relpath(A, snap)))
        for o in out:
            if o[0] == 'OTHER': print('     ', o[1], '| A:', str(o[2])[:110], '| B:', str(o[3])[:110])
    print('REGRESSION', 'PASS' if nother == 0 else 'FAIL (%d files with non-path differences)' % nother)
    if '--json' in sys.argv:
        json.dump(dict(snap=snap, ref=ref, pass_=nother == 0, files=res), open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1, default=str)
