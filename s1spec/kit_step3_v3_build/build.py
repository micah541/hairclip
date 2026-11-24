#!/usr/bin/env python3
"""Builder for the reproduction kit "kit_step3 v3" (step-3 certificates of the v3.3 chain).

  python3 build.py [--n35assert]

Creates /home/claude/s1spec/kit_step3_v3/ (tree, wiped first), /home/claude/s1spec/kit_step3_v3.zip and
/home/claude/s1spec/kit_step3_v3_build/BUILD_REPORT.md.  Reads (never writes) the proof tree.  Runs no certified computation.

Spec: writeup/lemmaU/bulk/t7_assembly/KITS_TODO.md sections A, D.1, E.2 (E.2 wins over D.1 wins over A) and V33_INDEX.md 4.6 / 4.4.
Convention (as kit_step3 v2 and kit_step4b): every script that had a hard-coded /home/claude/... path gets a `_kit_root()` header (walks up to the
KIT_ROOT marker file) and `(_S1 + '/rel/path')` literals, each patched line marked `[kit patch]`; the patch list goes to BUILD_REPORT.md.

--n35assert : use the *_n35assert.py variants of assembly_rep.py, ic_rep.py and assemble_U.py (E.2 literally).  DEFAULT is the on-disk
              scripts that produced the pinned outputs (F34: the assert variants change the string IC.step4_H2_range and the script hashes,
              so the pinned A_interface 78fe2f8e... would not be reproduced up to `args`).
"""
import os, re, sys, ast, json, shutil, hashlib, zipfile, tokenize, io, time

S = '/home/claude/s1spec'
B = S + '/writeup/lemmaU/bulk'
T2 = B + '/fix2/n12/t2_A_v35'
P1 = S + '/writeup/lemmaU/phase1'
P2A = S + '/writeup/lemmaU/phase2a'
V = S + '/step3nl_v33/verifier'
S4O = S + '/step4/p2/out_v33'
OLD = S + '/kit_step3'
HERE = os.path.dirname(os.path.abspath(__file__))
K = S + '/kit_step3_v3'
ZIP = S + '/kit_step3_v3.zip'
N35 = '--n35assert' in sys.argv

# ------------------------------------------------------------------------------------------------ expected pins (V33_INDEX 4.6 / 4.7, KITS_TODO E.2)
EXPECT = {   # path in proof tree -> sha256 (full) ; prefixes are from V33_INDEX, full values from the SHA256SUMS of the pinned out/ dirs
    T2 + '/out/A_interface.json': '78fe2f8e38a340ba515fbc33cbab60cdd85231095af6ccd192e5f6707f457e8c',
    T2 + '/out/asm_design_N2ref4.json': '46bae149800904875e5b510d02de7843dc791f890d53b9a3e45338e63a19aa16',
    B + '/t4_cert/out/u1_ghat.json': '5a06d24cae8591f508f9f7499c7d2b4355e06cd5b9d8d2419a9a0108c5638466',
    B + '/t4_cert/out/u2_K.json': '10ae3c383f7098284e5ae16f4355db9a2645f64531f62c6565c5824213f73261',
    B + '/t4_cert/out/u34_primary.json': 'abaf73724b84a0845b549d74a32050e5efafc8f246360ea2d8c388dd0dd70333',
    B + '/t4_cert/out/U_cert_v35.json': '960809fdc18dd67c1770fa5b77b2658081d0c81cb24f1b246f48ca4a583d5f54',
    B + '/t4_cert/out/U_cert_for_zone3_v35.json': '22c13e4403a1b41ab09ab15a9d8852f153832a97742895757fdb1ce3fd0ca19f',
    B + '/t5_zone3/out/zone3_interface.json': '45e8d4dcaff2bd75fcea725f80e98625af70a40709fc2cf5d3df78b5becb3ea3',
    B + '/t5_zone3/out/ramp_U.json': '55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de',
    S4O + '/FINAL_manifest_v33.json': '8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec',
    S4O + '/b0_dprof.json': '7fdf378c5f59f098263aec8f39bfff7f365e0d25cbbcd154c3f1010bd4a2cd8e',
    S4O + '/n1_chart.json': '341ed1f442ac64114bea6943cf37466a1515361bfcf838b41a7f7d860bf1ece0',
    V + '/rate_table_fin12f4_a0.00015.json': 'b20e87153440144b0537afa4cd4b9dff2d788ffd2ce2163626faf973b13d2ad1',
    V + '/rate_table_fin12f4_a0.000125.json': '79b60067f0f935eec8807a7641ba0af0b514755f7ac0ae553488f1a0e36fef65',
    V + '/modsups_fin12f4/final_alpha_points_5e-8.json': 'd5b14797e9' ,   # prefix only: completed below
    V + '/verify_log_fin12f4_a0.00015.txt': 'f028ed4a1275563de73dee145934f5ce49dac856da1582c5f53f8592246b9a60',
    V + '/verify_log_fin12f4_a0.000125.txt': '303f0b38d06315a5c3ff9ecb2b736b09f751bf4094227370260989783ee627f2',
    V + '/tail_cert_log.txt': 'd2237124555e0f64802797ce28a8ad5457427839eb933513d1e7a3aaaaaa945f',
    B + '/t4_cert/out/U_cert_v35_pass3.json': '35682a55eacea310da544bd8095fdd840966a59b4cd8ac820bb7bbd0602f5f7d',
    B + '/t4_cert/out/U_cert_v35_pass1.json': 'a28538def0cb69a707d2d8c17db5698c385027f6c5009647c6a83413cbc7e0b8',
    B + '/t5_zone3/out/zone3_interface.v34.json': 'a381b1f22248a3adc62fc085086e170d75d1388dd086ea6fc1aff9de7ec704c0',
}
SCRIPT_EXPECT = {   # source scripts (pre-patch) : sha256 prefix, from t4_cert/out/SHA256SUMS_scripts_m14, t5_zone3/SHA256SUMS_scripts, V33_INDEX 4.6
    B + '/t4_cert/u2_K.py': '059a2734', B + '/t4_cert/t4lib.py': '85ddd5f3', B + '/t4_cert/u1_ghat.py': 'daca142d', B + '/t4_cert/u34_sweeps.py': '8d3eaf29',
    B + '/t4_cert/make_ucert_for_zone3.py': '1e0a25c3', B + '/t5_zone3/zone3.py': 'b9f23489', B + '/t5_zone3/ramp_U.py': '0df7ce69', B + '/t5_zone3/z3lib.py': 'fa59483e',
    T2 + '/rigor/export_t2.py': '067e8c4d',
}
if N35:
    SCRIPT_EXPECT[B + '/t4_cert/assemble_U_n35assert.py'] = '6ee50b4d'; SCRIPT_EXPECT[T2 + '/rigor/assembly_rep_n35assert.py'] = '4aef3de6'; SCRIPT_EXPECT[T2 + '/rigor/ic_rep_n35assert.py'] = '59b66e1a'
else:
    SCRIPT_EXPECT[B + '/t4_cert/assemble_U.py'] = 'd14baeea'; SCRIPT_EXPECT[T2 + '/rigor/assembly_rep.py'] = '77d4abe6'

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

REPORT = dict(files=[], patches=[], notes=[], todo=[], asserts=[])

def note(s): REPORT['notes'].append(s)
def todo(s): REPORT['todo'].append(s); print('[TODO]', s)

# ------------------------------------------------------------------------------------------------ preconditions
for p, h in list(EXPECT.items()) + [(p, h) for p, h in SCRIPT_EXPECT.items()]:
    if not os.path.exists(p): sys.exit('missing source: ' + p)
    a = sha(p)
    if not a.startswith(h[:8]): sys.exit('HASH MISMATCH %s: %s vs expected %s' % (p, a[:16], h[:16]))
    REPORT['asserts'].append((p.replace(S + '/', ''), a))
EXPECT[V + '/modsups_fin12f4/final_alpha_points_5e-8.json'] = sha(V + '/modsups_fin12f4/final_alpha_points_5e-8.json')
if os.path.exists(K): shutil.rmtree(K)
os.makedirs(K)

# ------------------------------------------------------------------------------------------------ helpers
def cp(src, dst, binary=True):
    dst = K + '/' + dst
    os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copyfile(src, dst); return dst

def cpdir(src, dst, ignore=('__pycache__',)):
    shutil.copytree(src, K + '/' + dst, ignore=shutil.ignore_patterns(*ignore))

# path map (proof tree path relative to /home/claude/s1spec -> kit path); everything not listed keeps its relative position
MAP = [('/writeup/lemmaU/bulk/fix2/n12/t2_A_v35', '/writeup/lemmaU/bulk/t2_A_v35'),
       ('/writeup/lemmaU/bulk/fix1/lemmaU_rev', '/writeup/lemmaU/bulk/text')]
LIT = re.compile(r"""(?<![\w])(['"])(/home/claude(?:/s1spec|/s1_tight|/s1_cap)?)([^'"\n]*)\1""")
HDR = ("import os as _kos\n"
       "def _kit_root(_d=_kos.path.dirname(_kos.path.abspath(__file__))):\n"
       "    while not _kos.path.exists(_kos.path.join(_d, 'KIT_ROOT')) and _d != _kos.path.dirname(_d): _d = _kos.path.dirname(_d)\n"
       "    return _d\n"
       "_S1 = _kit_root().replace('\\\\', '/')   # [kit patch] the kit root (was an absolute path of the proof tree)\n")

def docstring_lines(tree):
    rng = []
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.body and isinstance(n.body[0], ast.Expr) \
                and isinstance(n.body[0].value, ast.Constant) and isinstance(n.body[0].value.value, str):
            rng.append((n.body[0].lineno, n.body[0].end_lineno))
    return rng

def patch_py(src, dst, extra=()):
    """copy src to kit path dst with the [kit patch] changes; extra = [(old, new)] exact replacements (asserted once each)"""
    text = open(src, encoding='utf-8').read()
    rel = dst
    assert 'from __future__' not in text, src
    patched = []
    for old, new in extra:
        assert text.count(old) == 1, (src, old[:60], text.count(old))
        text = text.replace(old, new); patched.append(('manual', old.strip()[:100], new.strip()[:100]))
    lines = text.split('\n'); used = False
    for i, l in enumerate(lines):
        def rp(m):
            nonlocal used
            root, rest = m.group(2), m.group(3)
            if root == '/home/claude/s1spec':
                for a, b in MAP:
                    if rest.startswith(a): rest = b + rest[len(a):]; break
                used = True
                new = '_S1' if rest == '' else "(_S1 + '%s')" % rest
            else:      # /home/claude/s1_tight etc. (not expected in this kit)
                used = True; new = "(_S1 + '%s%s')" % (root[len('/home/claude'):], rest)
            patched.append((i + 1, m.group(0), new)); return new
        l2 = LIT.sub(rp, l)
        if l2 != l and not l2.rstrip().endswith('\\'): l2 += '  # [kit patch]'
        lines[i] = l2
    text = '\n'.join(lines)
    # remaining occurrences must be in comments or docstrings only; neutralise them
    tree = ast.parse(text)
    dl = docstring_lines(tree)
    for i, l in enumerate(text.split('\n')):
        if '/home/claude' in l:
            ln = i + 1
            code = l.split('#')[0]
            if '/home/claude' in code and not any(a <= ln <= b for a, b in dl):
                sys.exit('unpatched code reference in %s:%d: %s' % (src, ln, l.strip()[:120]))
    text = text.replace('/home/claude/s1spec', '<proof-tree>').replace('/home/claude', '<proof-tree>')
    if used or any(('_S1' in l) for l in text.split('\n')):
        # insert header after docstring / leading comment lines
        tree = ast.parse(text); ls = text.split('\n'); at = 0
        if tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(tree.body[0].value, ast.Constant) and isinstance(tree.body[0].value.value, str):
            at = tree.body[0].end_lineno
        else:
            while at < len(ls) and ls[at].startswith('#'): at += 1
        ls[at:at] = HDR.rstrip('\n').split('\n'); text = '\n'.join(ls)
    ast.parse(text)
    out = K + '/' + dst
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w', encoding='utf-8', newline='\n') as f: f.write(text)
    for p in patched: REPORT['patches'].append((dst,) + tuple(p))
    return out

def cp_py(src, dst, extra=()):
    return patch_py(src, dst, extra)

# ================================================================================================ 1. unchanged old parts (ibp, lemmaO, verifier code, support tree)
for d in ('step3', 'step3nl', 'lin', 's1_cap', 'cap', 's1_tight'):
    cpdir(OLD + '/' + d, d)
cpdir(OLD + '/writeup/lemmaO', 'writeup/lemmaO')
cp(OLD + '/KIT_ROOT', 'KIT_ROOT')
os.remove(K + '/step3nl/verifier/rate_table_fin11_a0.00015.json')       # v3.1 table: not used by v3.3
for d in ('ibp', 'lemmaO'): cpdir(OLD + '/reference/' + d, 'reference/' + d)
# the v3.3 verifier has byte-identical code (checked here)
for f in ('bb', 'closure', 'constants', 'ival', 'mod_sups', 'rhs_iv', 'tail_cert', 'unc'):
    a = open(OLD + '/step3nl/verifier/%s.py' % f, 'rb').read(); b = open(V + '/%s.py' % f, 'rb').read()
    assert a == b, 'verifier code of kit_step3 v2 differs from step3nl_v33/verifier: ' + f
note('verifier code (bb, closure, constants, ival, mod_sups, rhs_iv, tail_cert, unc) is byte-identical in kit_step3 v2 and step3nl_v33/verifier (asserted); only the environment (FLUX0=1.64e-20, FLUX2=8.0e-20, TAG=fin12f4) differs.')

# ================================================================================================ 2. verifier fin12f4 pinned inputs + references
cp(V + '/rate_table_fin12f4_a0.00015.json', 'step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json')
cp(V + '/modsups_fin12f4/final_alpha_points_5e-8.json', 'step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json')
for f in ('rate_table_fin12f4_a0.00015.json', 'rate_table_fin12f4_a0.000125.json', 'verify_log_fin12f4_a0.00015.txt', 'verify_log_fin12f4_a0.000125.txt', 'tail_cert_log.txt'):
    cp(V + '/' + f, 'reference/verifier/' + f)
cp(V + '/modsups_fin12f4/final_alpha_points_5e-8.json', 'reference/verifier/final_alpha_points_5e-8.json')

# ================================================================================================ 3. repaired Lemma S(A): t2_A_v35 chain
TK = 'writeup/lemmaU/bulk/t2_A_v35'
rig = ['a3_signs', 'arlib', 'export_t2', 'hjw_rep', 'outer_L2', 'ramp_rep', 'weights_rep', 'weights_t2']
for f in rig: cp_py(T2 + '/rigor/%s.py' % f, TK + '/rigor/%s.py' % f)
cp_py(T2 + '/rigor/' + ('assembly_rep_n35assert.py' if N35 else 'assembly_rep.py'), TK + '/rigor/assembly_rep.py')
cp_py(T2 + '/rigor/' + ('ic_rep_n35assert.py' if N35 else 'ic_rep.py'), TK + '/rigor/ic_rep.py')
for f in ('lemmaS_B3', 'lemmaS_F1', 'lemmaS_F2', 'lemmaS_L4terms', 'lemmaS_interval', 'lemmaS_prop1_bounds', 'lemmaS_prop1_sym', 'run_F1', 'run_F10p_ibp'):
    cp_py(T2 + '/base_ibp/%s.py' % f, TK + '/base_ibp/%s.py' % f)
for f in ('arbio', 'assembly_arb', 'capture_F2', 'weights_exact'):
    cp_py(T2 + '/base_ibp/rigor/%s.py' % f, TK + '/base_ibp/rigor/%s.py' % f)
for f in ('prop1_forms.pkl', 'prop1_forms3.pkl', 'lemmaS_cflux.json', 'lemmaS_Bbox_1e-05_a0.00015_dv0.15.json', 'lemmaS_F4B_F10p_M0.52.json', 'rigor_assembly_arb.log'):
    cp(T2 + '/base_ibp/' + f, TK + '/base_ibp/' + f)
cp(T2 + '/base_ibp/rigor/ramp_tables.json', TK + '/base_ibp/rigor/ramp_tables.json')
cp(T2 + '/base_ibp/rigor/tables_exact.pkl', TK + '/base_ibp/rigor/tables_exact.pkl')
cp(T2 + '/base_v31H2/rigor/assembly_arb.json', TK + '/base_v31H2/rigor/assembly_arb.json')
cp(T2 + '/base_v31H2/rigor_assembly_arb.log', TK + '/base_v31H2/rigor_assembly_arb.log')
# v3.1 reference of assembly_rep (V31) and the shared exact tables
cp(S + '/writeup/lemmaS/rigor/assembly_arb.json', 'writeup/lemmaS/rigor/assembly_arb.json')
cp(S + '/writeup/lemmaS/rigor/assembly_arb.out', 'writeup/lemmaS/rigor/assembly_arb.out')
cp(S + '/writeup/lemmaS/rigor/tables_exact.pkl', 'writeup/lemmaS/rigor/tables_exact.pkl')
assert sha(K + '/writeup/lemmaS/rigor/tables_exact.pkl') == sha(K + '/' + TK + '/base_ibp/rigor/tables_exact.pkl')
# step-4 inputs (pinned)
cp(S4O + '/b0_dprof.json', 'step4/p2/out_v33/b0_dprof.json')
cp(S4O + '/n1_chart.json', 'step4/p2/out_v33/n1_chart.json')
# references of the (A) chain
T2OUT = ['A_interface', 'a3_signs', 'asm_design_N20', 'asm_design_N2ref', 'asm_design_N2ref4', 'asm_design_pm3ibp', 'asm_hi_N2ref', 'asm_lo_N2ref', 'asm_probe_N2ref',
         'hjw_rep', 'ic_design', 'ic_hi', 'ic_lo', 'ic_probe', 'outer_L2', 'ramp_rep', 'ramp_rep_pm3ibp', 'weights_t2', 'weights_t2_fb']
for f in T2OUT: cp(T2 + '/out/%s.json' % f, 'reference/lemmaSA/out/%s.json' % f)
for f, g in (('lemmaS_F10p_M0.52.json', 'lemmaS_F10p_M0.52.json'), ('rigor/cap_F10p.json', 'cap_F10p.json'), ('rigor/lemmaS_rigcap_M0.52.json', 'lemmaS_rigcap_M0.52.json'),
             ('rigor/weights_exact.json', 'weights_exact.json'), ('rigor/assembly_arb.json', 'assembly_arb.json')):
    cp(T2 + '/base_ibp/' + f, 'reference/lemmaSA/base/' + g)
nfiles_t2 = len(T2OUT)
note('run_t2.sh writes 19 JSON files plus A_interface.json = %d files in out/ (listed in KITS_TODO as "20"; the 20th would be SHA256SUMS, which is not shipped as it hashes the args blocks).' % nfiles_t2)

# ================================================================================================ 4. Lemma U + zone III
U = 'writeup/lemmaU'
t4 = ['t4lib', 'u1_ghat', 'u2_K', 'u34_sweeps', 'make_ucert_for_zone3']
cp_py(B + '/t4_cert/t4lib.py', U + '/bulk/t4_cert/t4lib.py', extra=[("T2 = ROOT + '/bulk/t2_A'", "T2 = ROOT + '/bulk/t2_A_v35'")])   # [kit patch] t2_A (v3.1 interface dir) -> t2_A_v35; only used by the unshipped u5_enum and the provisional branch
for f in t4[1:]: cp_py(B + '/t4_cert/%s.py' % f, U + '/bulk/t4_cert/%s.py' % f)
cp_py(B + '/t4_cert/' + ('assemble_U_n35assert.py' if N35 else 'assemble_U.py'), U + '/bulk/t4_cert/assemble_U.py')
for f in ('ramp_U', 'z3lib'): cp_py(B + '/t5_zone3/%s.py' % f, U + '/bulk/t5_zone3/%s.py' % f)
cp_py(B + '/t5_zone3/zone3.py', U + '/bulk/t5_zone3/zone3.py',
      extra=[("default='/home/claude/s1spec/writeup/lemmaU/bulk/t4_cert/out/U_cert_v34.json'", "default='/home/claude/s1spec/writeup/lemmaU/bulk/t4_cert/out/U_cert_v35.json'")])   # stale v34 default of --ucert_full
for f in ('ghat_arb', 'weight_consts_arb', 'arbio_local'): cp_py(P1 + '/ghat/%s.py' % f, U + '/phase1/ghat/%s.py' % f)
cp_py(P1 + '/ramp/ulib.py', U + '/phase1/ramp/ulib.py')
cp(P2A + '/U_ramp/ghat_design_arb.json', U + '/phase2a/U_ramp/ghat_design_arb.json')
cp(S + '/writeup/lemmaO_v32/lemmaO_v32_consts.json', 'writeup/lemmaO_v32/lemmaO_v32_consts.json')
cp(S + '/writeup/lemmaO_v32/lemmaR_consts.json', 'writeup/lemmaO_v32/lemmaR_consts.json')
cp(B + '/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json', U + '/bulk/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json')
# seeds of the circular chain (see README "Order")
cp(B + '/t4_cert/out/u1_ghat.json', 'seeds/u1_ghat.json')
cp(B + '/t5_zone3/out/zone3_interface.v34.json', 'seeds/zone3_interface.v34.json')
# references
for f in ('ramp_U', 'u1_ghat', 'u2_K', 'u34_primary', 'U_cert_v35', 'U_cert_v35_pass1', 'U_cert_v35_pass3', 'U_cert_for_zone3_v35', 'zone3_interface'):
    src = (B + '/t5_zone3/out/%s.json' if f in ('ramp_U', 'zone3_interface') else B + '/t4_cert/out/%s.json') % f
    cp(src, 'reference/lemmaU/%s.json' % f)
cp(S4O + '/FINAL_manifest_v33.json', 'reference/step4/FINAL_manifest_v33.json')

# ================================================================================================ 5. text checks
TXT = U + '/bulk/text'
cp(B + '/fix1/lemmaU_rev/lemmaU_v34.tex', TXT + '/lemmaU_v34.tex')
cp_py(B + '/t3_proof/check_balance.py', TXT + '/check_balance.py',
      extra=[("src = open(sys.argv[1]).read()", "src = open(sys.argv[1], encoding='utf-8').read()")])      # [kit patch] explicit encoding (Windows)
cp_py(B + '/fix1/lemmaU_rev/check_keys_v34.py', TXT + '/check_keys_v34.py', extra=[
    ("CN = json.load(open(os.path.join(HERE, 'constants_needed_v34.json')))",
     "CN = json.load(open(os.path.join(HERE, 'constants_needed_v34.json'), encoding='utf-8'))\n"
     "_KR = os.path.normpath(os.path.join(HERE, '..', '..', '..', '..'))   # [kit patch] kit root = HERE/../../../..\n"
     "for _f in CN['files'].values(): _f['path'] = os.path.join(_KR, _f['path'])   # [kit patch] paths in constants_needed_v34.json are kit-relative (reference/...)"),
    ("TEX = open(os.path.join(HERE, 'lemmaU_v34.tex')).read()", "TEX = open(os.path.join(HERE, 'lemmaU_v34.tex'), encoding='utf-8').read()"),
    ("D[k] = json.load(open(f['path']))", "D[k] = json.load(open(f['path'], encoding='utf-8'))"),
    ("pins = D['U']['args']['inputs']", "pins = {os.path.basename(k): v for k, v in D['U']['args']['inputs'].items()}   # [kit patch] pins compared by file name (the recorded paths are those of the proof tree)"),
    ("if pins.get(CN['files'][f]['path']) != CN['files'][f]['sha256']:", "if pins.get(os.path.basename(CN['files'][f]['path'])) != CN['files'][f]['sha256']:"),
    ("P3 = json.load(open(p3))", "P3 = json.load(open(p3, encoding='utf-8'))"),
    ("pin3 = {os.path.normpath(os.path.join(P3['args']['cwd'], k)): v for k, v in P3['args']['inputs'].items()}", "pin3 = {os.path.basename(k): v for k, v in P3['args']['inputs'].items()}"),
    ("pin3.get(CN['files'][f]['path']) == CN['files'][f]['sha256']", "pin3.get(os.path.basename(CN['files'][f]['path'])) == CN['files'][f]['sha256']"),
])
KREL = {'A': 'reference/lemmaSA/out/A_interface.json', 'A4': 'reference/lemmaSA/out/asm_design_N2ref4.json', 'U': 'reference/lemmaU/U_cert_v35.json',
        'U1': 'reference/lemmaU/u1_ghat.json', 'Z': 'reference/lemmaU/zone3_interface.json', 'R': 'reference/lemmaU/ramp_U.json', 'S4': 'reference/step4/FINAL_manifest_v33.json'}
CNj = json.load(open(B + '/fix1/lemmaU_rev/constants_needed_v34.json', encoding='utf-8'))
for k, f in CNj['files'].items():
    assert os.path.basename(f['path']) == os.path.basename(KREL[k]), (k, f['path'])
    assert sha(K + '/' + KREL[k]) == f['sha256'], ('constants_needed_v34.json pins %s but the shipped reference differs' % k)
    f['path'] = KREL[k]
CNj['note'] = CNj['note'] + ' [kit patch: the "path" entries are kit-relative references to the pinned outputs of the proof tree; sha256 unchanged]'
assert '/home/claude' not in json.dumps(CNj['constants']) and '/home/claude' not in json.dumps(CNj['numbers'])
os.makedirs(K + '/' + TXT, exist_ok=True)
json.dump(CNj, open(K + '/' + TXT + '/constants_needed_v34.json', 'w', encoding='utf-8', newline='\n'), indent=1)
NUMBERS = CNj['numbers']
note('lemmaU_v34.tex sha256 %s, constants_needed_v34.json (source) sha256 %s: V33_INDEX 4.6 quotes 5131776d... / 0f3c7455... (before the F34 edits); the kit ships the current files.' % (sha(B + '/fix1/lemmaU_rev/lemmaU_v34.tex')[:12], sha(B + '/fix1/lemmaU_rev/constants_needed_v34.json')[:12]))

# ================================================================================================ 6. drivers
EXP_JSON = {os.path.relpath(p, S): h for p, h in EXPECT.items()}
shutil.copyfile(HERE + '/run_all.py', K + '/run_all.py')
tmpl = open(HERE + '/check_template.py', encoding='utf-8').read()
old = open(OLD + '/check.py', encoding='utf-8').read()
helpers = old[old.index('KIT = os.path.dirname'):old.index('# ------', old.index('def need('))]
ibp_blk = old[old.index('def check_ibp():'):old.index('# ------', old.index('def check_ibp():'))]
lo_blk = old[old.index('def check_lemmaO():'):old.index('def main():')]
doc_rules = old[old.index('Rules (R1-R3), printed again at the end:'):old.index('usage: python check.py')]
tmpl = tmpl.replace('#@@HELPERS@@', helpers).replace('#@@IBP@@', ibp_blk).replace('#@@LEMMAO@@', lo_blk).replace('#@@RULES_OLD@@', doc_rules.replace('"""', "'''"))
short = {k: v for k, v in EXP_JSON.items()}
tmpl = tmpl.replace('@@EXPECT_SHA@@', json.dumps({os.path.basename(k): v for k, v in short.items()}, indent=1, sort_keys=True))
tmpl = tmpl.replace('@@NUMBERS_NOTE@@', '%d quoted numbers (constants_needed_v34.json)' % len(NUMBERS))
with open(K + '/check.py', 'w', encoding='utf-8', newline='\n') as f: f.write(tmpl)
ast.parse(tmpl)
with open(K + '/requirements.txt', 'w', newline='\n') as f:
    f.write('# versions the kit was tested with (CPython 3.13 / 3.14, Linux); wheels for cp313 win_amd64 exist for all of them\n'
            'python-flint==0.9.0\nsympy==1.14.0\nnumpy==2.5.4\nscipy==1.18.1\nmpmath==1.3.0\n')
for p in ('run_all.py',): ast.parse(open(K + '/' + p, encoding='utf-8').read())
shutil.copyfile(HERE + '/README_template.md', K + '/README.md')

# ---- scrub '/home/claude' from comments/docstrings of the old-kit scripts (before hashing)
def scrub_py(p):
    """neutralise '/home/claude' in comments and docstrings of an (old-kit) script; refuse if it occurs in code"""
    text = open(p, encoding='utf-8').read()
    if '/home/claude' not in text: return False
    tree = ast.parse(text); dl = docstring_lines(tree)
    for i, l in enumerate(text.split('\n')):
        if '/home/claude' in l and '/home/claude' in l.split('#')[0] and not any(a <= i + 1 <= b for a, b in dl):
            sys.exit('code reference to /home/claude in %s:%d: %s' % (p, i + 1, l.strip()[:120]))
    text = text.replace('/home/claude/s1spec', '<proof-tree>').replace('/home/claude', '<proof-tree>')
    ast.parse(text)
    with open(p, 'w', encoding='utf-8', newline='\n') as f: f.write(text)
    return True
nscrub = 0
for dp, dn, fs in os.walk(K):
    for f in fs:
        if f.endswith('.py'): nscrub += scrub_py(os.path.join(dp, f))
note('%d old-kit / patched scripts had "/home/claude" only inside comments or docstrings (e.g. the old "[kit patch] was the absolute path ..." remarks); replaced by "<proof-tree>" so that the grep proof is literal. Only comments/docstrings changed: build.py refuses if the string occurs in executable code.' % nscrub)

# ================================================================================================ 7. input manifest (all shipped files except reference/, drivers and outputs)
man = []
for dp, dn, fs in os.walk(K):
    dn[:] = sorted(d for d in dn if d != '__pycache__')
    for f in sorted(fs):
        p = os.path.join(dp, f); r = os.path.relpath(p, K).replace(os.sep, '/')
        if r.startswith(('reference/', 'inputs/')) or r in ('README.md', 'run_all.py', 'check.py'): continue
        man.append('%s  %s' % (sha(p), r))
os.makedirs(K + '/inputs', exist_ok=True)
with open(K + '/inputs/MANIFEST.sha256', 'w', newline='\n') as f:
    f.write('# SHA-256 of every shipped source and data input (kit-patched scripts included); run_all.py smoke verifies them.\n' + '\n'.join(man) + '\n')
refman = []
for dp, dn, fs in os.walk(K + '/reference'):
    for f in sorted(fs):
        p = os.path.join(dp, f); refman.append('%s  %s' % (sha(p), os.path.relpath(p, K).replace(os.sep, '/')))
refman.sort(key=lambda s: s.split('  ')[1])
with open(K + '/inputs/REFERENCE.sha256', 'w', newline='\n') as f:
    f.write('# SHA-256 of the shipped reference outputs (pinned outputs of the proof tree)\n' + '\n'.join(refman) + '\n')

# ================================================================================================ 8. grep proof, syntax
bad = []
for dp, dn, fs in os.walk(K):
    for f in fs:
        p = os.path.join(dp, f); r = os.path.relpath(p, K)
        if f.endswith(('.py', '.sh', '.md', '.txt', '.bat', '.ps1', '.sha256')) or f == 'KIT_ROOT':
            t = open(p, encoding='utf-8', errors='replace').read()
            if '/home/claude' in t: bad.append(r)
        if f.endswith('.py'):
            ast.parse(open(p, encoding='utf-8').read())
if bad: sys.exit('references to /home/claude remain in: %s' % bad)
data_hits = []
for dp, dn, fs in os.walk(K):
    for f in fs:
        if f.endswith(('.json', '.out', '.log', '.pkl')):
            p = os.path.join(dp, f)
            if b'/home/claude' in open(p, 'rb').read(): data_hits.append(os.path.relpath(p, K))
note('grep proof: no .py/.sh/.md/.txt/.sha256 file of the kit contains "/home/claude" (checked by build.py). %d data files (pinned reference outputs / seeds / pinned logs) contain the string inside provenance fields (`args` blocks, `src`/`path` strings of the original run); they are only parsed as data, never opened by that path (check.py compares path-valued strings by file name).' % len(data_hits))
REPORT['data_hits'] = sorted(data_hits)

# ================================================================================================ 9. zip + report
if os.path.exists(ZIP): os.remove(ZIP)
names = []
for dp, dn, fs in os.walk(K):
    dn[:] = sorted(dn)
    for f in sorted(fs): names.append(os.path.join(dp, f))
names.sort(key=lambda p: os.path.relpath(p, K))
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in names:
        zi = zipfile.ZipInfo('kit_step3_v3/' + os.path.relpath(p, K).replace(os.sep, '/'), date_time=(2026, 11, 14, 12, 0, 0))
        zi.compress_type = zipfile.ZIP_DEFLATED; zi.external_attr = 0o644 << 16
        z.writestr(zi, open(p, 'rb').read())
total = sum(os.path.getsize(p) for p in names)
REPORT['zip'] = (os.path.getsize(ZIP), sha(ZIP)); REPORT['tree_bytes'] = total; REPORT['nfiles'] = len(names)
REPORT['files'] = [(os.path.relpath(p, K).replace(os.sep, '/'), os.path.getsize(p), sha(p)) for p in names]

TODOS = open(HERE + '/todo_list.txt', encoding='utf-8').read().strip().split('\n') if os.path.exists(HERE + '/todo_list.txt') else []
with open(HERE + '/BUILD_REPORT.md', 'w', encoding='utf-8', newline='\n') as f:
    w = f.write
    w('# BUILD_REPORT: kit_step3 v3\n\nBuilt by `kit_step3_v3_build/build.py` (%s) from the proof tree `<proof-tree>` = /home/claude/s1spec. Nothing was run: no certified computation, only hashing, `ast.parse` and zip.\n\n' % time.strftime('%Y-%m-%d'))
    w('Variant: %s\n\n' % ('`--n35assert` (assert variants of assembly_rep, ic_rep, assemble_U)' if N35 else 'DEFAULT (on-disk scripts that produced the pinned outputs; see build.py docstring)'))
    w('Tree `kit_step3_v3/`: %d files, %d bytes. Zip `kit_step3_v3.zip`: %d bytes, sha256 `%s`.\n\n' % (len(names), total, REPORT['zip'][0], REPORT['zip'][1]))
    w('## Pinned source files verified at build time (sha256 of the proof-tree file)\n\n| file | sha256 |\n|---|---|\n')
    for p, h in REPORT['asserts']: w('| `%s` | `%s` |\n' % (p, h))
    w('\n## Notes\n\n')
    for n in REPORT['notes']: w('- %s\n' % n)
    w('\n## [TODO] and unresolved items (stated plainly, nothing guessed)\n\n')
    for t in TODOS + REPORT['todo']: w('- %s\n' % t)
    w('\n## Patches (`[kit patch]`), %d entries\n\nFormat: kit file, line (in the proof-tree source, before header insertion) or `manual`, old text, new text. The `_kit_root()` header is inserted after the module docstring of every file that has a `_S1` reference.\n\n' % len(REPORT['patches']))
    w('| kit file | line | old | new |\n|---|---|---|---|\n')
    for p in REPORT['patches']:
        w('| `%s` | %s | `%s` | `%s` |\n' % (p[0], p[1], str(p[2]).replace('|', '/').replace('`', "'"), str(p[3]).replace('|', '/').replace('`', "'")))
    w('\n## Data files that contain the string /home/claude (provenance fields only)\n\n')
    for d in REPORT['data_hits']: w('- `%s`\n' % d)
    w('\n## Every file of the kit (path, bytes, sha256)\n\n| path | bytes | sha256 |\n|---|---|---|\n')
    for r, n, h in REPORT['files']: w('| `%s` | %d | `%s` |\n' % (r, n, h))
print('built', K, len(names), 'files', total, 'bytes; zip', REPORT['zip'])
