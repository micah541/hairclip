"""[v4: kit_step4b v3 + the 29-stage step-4 v38 chain (v38_block.py), write guard kit_guard.py; writes kit_step4b_v4]
assemble kit_step4b v3 from the proof tree (not shipped): kit B (34-stage v3.2 chain, with the six fixes of the first end-to-end run, kit4b_run_wsl1/RUN_NOTES.md,
plus kit fix 7) + the 21 stages of the v3.3 step-4 chain.  Writes /home/claude/s1spec/kit_step4b_v3 (tree), kit_step4b_v3.zip (deterministic) and BUILD_REPORT.md next to this file.
Never writes outside those three; reads the proof tree read-only.   usage: python3 build.py"""
import os, re, shutil, sys, ast, json, hashlib, zipfile, tokenize, io
S = '/home/claude/s1spec'; K = S + '/kit_step4b_v4'; BD = os.path.dirname(os.path.abspath(__file__)); ZIP = S + '/kit_step4b_v4.zip'   # v4
sys.dont_write_bytecode = True
FM33P = S + '/step4/p2/out_v33/FINAL_manifest_v33.json'
FM33_SHA = '8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec'
if hashlib.sha256(open(FM33P, 'rb').read()).hexdigest() != FM33_SHA: sys.exit('FINAL_manifest_v33.json is not the 21-stage manifest 8bc47e83...: refusing to build')
_FMP = os.environ.get('KIT_FM', S + '/step4/p2/out/FINAL_manifest.json')
_fm = json.load(open(_FMP)); _pa = [i for i, t in enumerate(_fm['stages']) if t.get('pass_A')]
if len(_fm['stages']) < 34 or _pa != list(range(_pa[0], _pa[0] + len(_pa))):
    sys.exit('FINAL_manifest is not a complete final manifest (%d stages; pass_A at %s): refusing to rebuild' % (len(_fm['stages']), _pa))
if os.path.exists(K): shutil.rmtree(K)
os.makedirs(K)
def cp(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copy2(src, dst)
def cpdir_py(sd, dd, extra=(), skip=()):
    for f in sorted(os.listdir(sd)):
        p = os.path.join(sd, f)
        if os.path.isfile(p) and (f.endswith('.py') or f in extra) and f not in skip: cp(p, os.path.join(dd, f))
# --- step-3 support tree from kit_step3 (already patched, tested on Windows)
K3 = S + '/kit_step3'
for d in ('step3', 'step3nl', 'lin', 's1_cap', 'cap'):
    shutil.copytree(K3 + '/' + d, K + '/' + d, ignore=shutil.ignore_patterns('__pycache__'))
cp(K3 + '/KIT_ROOT', K + '/KIT_ROOT')
# --- step4 sources
S4 = S + '/step4'; D4 = K + '/step4'
cpdir_py(S4 + '/t1', D4 + '/t1')
cpdir_py(S4 + '/t6', D4 + '/t6')
cpdir_py(S4 + '/t7', D4 + '/t7', extra=('t8req.json', 'phi_store.json', 'fixed_coarse.json', 'eicon.json', 'cx.json', 'closure.json', 'c2float.json', 'c1float_clamped_0.1.json'))
shutil.copytree(S4 + '/t7/nlsym', D4 + '/t7/nlsym', ignore=shutil.ignore_patterns('__pycache__'))     # nonlinear-closure form generators + ramp forms
cpdir_py(S4 + '/t8', D4 + '/t8', extra=('forms.json', 'lf_forms.pkl'))
cpdir_py(S4 + '/p2', D4 + '/p2')
cpdir_py(S4 + '/p2/t8se', D4 + '/p2/t8se', extra=('forms.json', 'lf_forms.pkl'))
cpdir_py(S4 + '/p2/t9h', D4 + '/p2/t9h')
cpdir_py(S4 + '/p2/t6c', D4 + '/p2/t6c')                                                  # B1: direct T6 residual of the corrected base
for rel_ in ('p2/h6/fable/h6f_cases.json', 'p2/h6/H6_PLAN.md'):
    if os.path.exists(S4 + '/' + rel_): cp(S4 + '/' + rel_, D4 + '/' + rel_)
for rel_ in ('writeup/lemmaO/lemmaO_v3_consts.json', 'writeup/lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json', 'writeup/lemmaO_v32/lemmaO_v32_consts.json'): cp(S + '/' + rel_, K + '/' + rel_)
for k_ in range(16):                                                                         # the dt balls of the proof tree were computed from these earlier point enclosures
    f_ = f't9_h1_p{-1 + (2 * k_ + 1) / 16}_coh_s180.json'; cp(S4 + '/p2/out/h1slab/pre_u4/' + f_, K + '/reference/h1_pre_u4/' + f_)
cp('/home/claude/s1_tight/consistency/check_outA5.py', K + '/s1_tight/consistency/check_outA5.py')      # B8
cp(os.path.join(BD, 'h1slab_summaries.py'), D4 + '/p2/t9h/h1slab_summaries.py')
# P1 reference checkpoints 95, 96 (inputs)
for k in (95, 96):
    for e in ('json', 'npz'): cp(f'{S4}/t3/p1/ckpt/slab_{k:05d}.{e}', f'{D4}/t3/p1/ckpt/slab_{k:05d}.{e}')
# S1 data
SD = '/home/claude/s1_tight'
cp(SD + '/prof/prof_tight.txt', K + '/s1_tight/prof/prof_tight.txt')
cp(SD + '/ext813/ext813_prof.txt', K + '/s1_tight/ext813/ext813_prof.txt')
cp(SD + '/ext_tip/prof_tip.txt', K + '/s1_tight/ext_tip/prof_tip.txt')
cp(SD + '/ext_tip/ext_tip_check.json', K + '/s1_tight/ext_tip/ext_tip_check.json')
cp(SD + '/ext813/ext813_points.txt', K + '/s1_tight/ext813/ext813_points.txt')
cp(S + '/step3/inA5.txt', K + '/step3/inA5.txt')
for f in ('bracket2.json', 'laplace.py', 'tail2.py'): cp('/home/claude/s1spec/kit_step3/s1_tight/' + f, K + '/s1_tight/' + f)
cp(S4 + '/t1/polished_v2/manifest.json', K + '/reference/polished_v2_manifest.json')
# --- path patching of step4 sources
HDR = ("import os as _kos\n"
       "def _kit_root(_d=_kos.path.dirname(_kos.path.abspath(__file__))):\n"
       "    while not _kos.path.exists(_kos.path.join(_d, 'KIT_ROOT')) and _d != _kos.path.dirname(_d): _d = _kos.path.dirname(_d)\n"
       "    return _d\n"
       "_S1 = _kit_root().replace('\\\\', '/')   # [kit patch] was the absolute path /home/claude/s1spec\n")
pat = re.compile(r"""(f?)(['"])/home/claude/(s1spec|s1_tight|s1_cap)(/[^'"]*|)\2""")
pat2 = re.compile(r"""(f?)(['"])/home/claude/(s1spec|s1_tight|s1_cap)(/[^'"]*|)(['"])?""")
def repl(m):
    f, q, root, rest = m.group(1), m.group(2), m.group(3), m.group(4)
    sub = {'s1spec': '', 's1_tight': '/s1_tight', 's1_cap': '/s1_cap'}[root]
    if f: return f"f{q}{{_S1}}{sub}{rest}{q}"
    return f"(_S1 + {q}{sub}{rest}{q})"
npatch = 0
for dp, dn, fs in [t for r_ in (D4, K + '/s1_tight/consistency') for t in os.walk(r_)]:
    for f in fs:
        if not f.endswith('.py'): continue
        p = os.path.join(dp, f); s = open(p, encoding='utf-8').read()
        if '/home/claude/' not in s: continue
        s2, n = pat.subn(repl, s)
        # f-strings with extra quote kinds are handled; any remaining non-docstring occurrences:
        s2 = HDR + s2
        try: ast.parse(s2)
        except SyntaxError as e: print('SYNTAX', p, e); sys.exit(1)
        rem = [l for l in s2.splitlines() if '/home/claude/' in l and not l.lstrip().startswith(('#', '"""')) ]
        open(p, 'w', encoding='utf-8').write(s2); npatch += 1
        if rem: print('REMAINING in', os.path.relpath(p, K), ':', [r.strip()[:110] for r in rem][:3])
print('patched', npatch, 'files')
p = D4 + '/p2/t8se/cconst.py'; t = open(p, encoding='utf-8').read()
t = t.replace("T7DIR = os.path.normpath(os.path.join(HERE, '..', 't7'))", "T7DIR = os.path.normpath(os.path.join(HERE, '..', '..', 't7'))   # [kit patch] p2/t7 was a symlink to t7")
assert "'..', '..', 't7'" in t; open(p, 'w', encoding='utf-8').write(t)
p = D4 + '/p2/t9h/win_slab.py'; t = open(p, encoding='utf-8').read()                          # point enclosures as the proof tree's dt/ift/win files were built (see run_all.py st_h1dt)
t = t.replace("O = (_S1 + '/step4/p2/out')\n", "O = (_S1 + '/step4/p2/out')\nimport os\nPTS = O if os.environ.get('KIT_DT_POINTS') == 'final' else (_S1 + '/reference/h1_pre_u4')   # [kit patch]\n", 1)
t = t.replace("json.load(open(f'{O}/t9_h1_p{tk}_coh_s180.json'))", "json.load(open(f'{PTS}/t9_h1_p{tk}_coh_s180.json'))", 1); assert 'PTS}/t9_h1_p' in t and 'KIT_DT_POINTS' in t
open(p, 'w', encoding='utf-8').write(t)
# ====================================================================================================================================================
# v3.3 step-4 chain (21 stages) + bulk inputs.  (spliced into build.py)
# ====================================================================================================================================================
FM33 = json.load(open(FM33P)); PATCHES = []          # PATCHES: (file, what) -> BUILD_REPORT.md
V33S = S4 + '/p2/v33'; V33D = D4 + '/p2/v33'; BULK = S + '/writeup/lemmaU/bulk'; KB = K + '/writeup/lemmaU/bulk'
def normrel(p):
    """a manifest path of the proof tree in the kit's normal form: relative to step4, or '@/' + relative to the kit root"""
    for a, b in (('/home/claude/s1spec/step4/', ''), ('/home/claude/s1spec/', '@/'), ('/home/claude/s1_cap/', '@/s1_cap/'), ('/home/claude/s1_tight/', '@/s1_tight/')):
        if p.startswith(a): return b + p[len(a):]
    return p
def kitpath(n): return K + '/' + n[2:] if n.startswith('@/') else D4 + '/' + n
def sanitize_nonexec(s, comments=True):
    """comments and docstrings that still name the proof tree's absolute path get neutral wording (not executed; keeps `grep /home/claude` meaningful)"""
    lines = s.splitlines(keepends=True); starts = [0]
    for l in lines: starts.append(starts[-1] + len(l))
    toks = list(tokenize.generate_tokens(io.StringIO(s).readline)); edits = []
    skip = (tokenize.NL, tokenize.COMMENT)
    for i, t in enumerate(toks):
        if '/home/claude' not in t.string: continue
        if t.type == tokenize.COMMENT: ok = comments
        elif t.type == tokenize.STRING:
            j = i - 1
            while j >= 0 and toks[j].type in skip: j -= 1
            k = i + 1
            ok = (j < 0 or toks[j].type in (tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT)) and toks[k].type == tokenize.NEWLINE
        else: ok = False
        if ok: edits.append((starts[t.start[0] - 1] + t.start[1], starts[t.end[0] - 1] + t.end[1], t.string))
    for a, b, txt in sorted(edits, reverse=True): s = s[:a] + txt.replace('/home/claude/s1spec', '<proof-tree>').replace('/home/claude', '<proof-tree-host>') + s[b:]
    return s
def code_hits(s):
    """string tokens / names left that contain the proof tree's absolute path (code, not comment/docstring)"""
    return [t.string[:100] for t in tokenize.generate_tokens(io.StringIO(s).readline) if '/home/claude' in t.string and t.type not in (tokenize.COMMENT,)]
def patch_py(dst, rel, pre=(), hdr=True):
    """[kit patch] one .py file of the kit tree: asserted single-site source edits `pre` (old, new[, count]), then the generic /home/claude -> KIT_ROOT rewrite"""
    s = open(dst, encoding='utf-8').read()
    for pr in pre:
        a, b = pr[0], pr[1]; cnt = pr[2] if len(pr) > 2 else 1
        assert s.count(a) == cnt, ('kit-patch target', rel, s.count(a), a[:90]); s = s.replace(a, b); PATCHES.append((rel, pr[3] if len(pr) > 3 else b.strip().splitlines()[0][:110]))
    s, n = pat.subn(repl, s)
    s = sanitize_nonexec(s)
    left = [h for h in code_hits(s)]
    if left: print('REMAINING absolute path in code:', rel, left[:4]); sys.exit(1)
    if hdr: s = HDR + s
    ast.parse(s)
    open(dst, 'w', encoding='utf-8').write(s)
    return n

# ---------------------------------------------------------------------------------------------------------------------------------- stage scripts
stage_scripts = sorted(f for f in os.listdir(V33S) if f.endswith('.py') and not f.startswith('final_tail_v33') and f != 'verify_v32_manifest.py')
for f in stage_scripts: cp(V33S + '/' + f, V33D + '/' + f)
cp(V33S + '/final_tail_v33.py', V33D + '/final_tail_v33.py')
WPRE = ("sys.executable, '-I', '-B'", "sys.executable, '-I', '-B', '-X', 'utf8'", 1)      # children keep the reference flags (-I -B); -X utf8 so that text files are UTF-8 on Windows too
PRE = {
 'h2f_balls.py': [WPRE + ("h2f_balls: child interpreter flags +-X utf8",)],
 'm9_balls.py': [WPRE + ("m9_balls: child interpreter flags +-X utf8",)],
 'f12_sliver.py': [
  ("PY = [sys.executable, '-I', '-B']; ME = os.path.abspath(__file__)", "PY = [sys.executable, '-I', '-B', '-X', 'utf8']; ME = os.path.abspath(__file__)", 1, "f12: child interpreter flags +-X utf8"),
  ("N1PY = f'{P2}/v33/n1_chart.py'; N1PIN = f'{OV}/n1_chart.json'",
   "N1PY = f'{P2}/v33/n1_chart.py'; PINV = '/home/claude/s1spec/reference/v33_out/p2/out_v33'; N1PIN = f'{PINV}/n1_chart.json'   # [kit patch] the PINNED n1_chart.json of the proof tree (shipped reference), see _KP/NORMP",
   1, "f12: pinned inputs b0_dprof/n1_chart/h2u_v33_arb are the proof tree's shipped copies (reference/v33_out), not this run's"),
  ("def strip_args(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x",
   "_NPR = re.compile(r'^(?:[A-Za-z]:)?/(?:[^/\\s]+/)*?(?=(?:step4|writeup|step3nl_v33|step3nl|step3|lin|s1_cap|s1_tight|reference)/)')\n"
   "def NORMP(s):\n"
   "    \"\"\"[kit patch] absolute file names inside pinned JSON are provenance of the proof tree; compare them modulo the install root ('@/...')\"\"\"\n"
   "    if not isinstance(s, str): return s\n"
   "    s = s.replace('\\\\', '/')\n"
   "    if not re.match(r'(?:[A-Za-z]:)?/', s): return s\n"
   "    return _NPR.sub('@/', s).replace('@/reference/v33_out/', '@/step4/')\n"
   "def _KP(p):\n"
   "    n = NORMP(p)\n"
   "    if n.startswith('@/step4/p2/out_v33/'): return PINV + n[len('@/step4/p2/out_v33'):]\n"
   "    return _S1 + n[1:] if n.startswith('@/') else p\n"
   "def _nt(x):\n"
   "    if isinstance(x, dict): return {k: _nt(v) for k, v in x.items()}\n"
   "    if isinstance(x, list): return [_nt(v) for v in x]\n"
   "    return NORMP(x)\n"
   "def strip_args(x): return _nt({k: v for k, v in x.items() if k != 'args'}) if isinstance(x, dict) else x",
   1, "f12: strip_args / leaves compare file names modulo the install root (NORMP)"),
  ("    else: yield pre, x\ndef rnd(x, n, up):", "    else: yield pre, (NORMP(x) if isinstance(x, str) else x)\ndef rnd(x, n, up):", 1, "f12: leaves() normalises file names"),
  ("    i = av.index('--chart'); assert av[i + 1] == N1PIN; o = av.index('--out'); assert av[o + 1] == f'out/asm_{nm}.json'\n    rep = list(av); rep[o + 1] = f'{REP}/asm_{nm}.json'\n",
   "    i = av.index('--chart'); assert NORMP(av[i + 1]).endswith('/p2/out_v33/n1_chart.json'); o = av.index('--out'); assert av[o + 1] == f'out/asm_{nm}.json'\n"
   "    d_ = av.index('--dprof'); assert NORMP(av[d_ + 1]).endswith('/p2/out_v33/b0_dprof.json'); av = list(av); av[d_ + 1] = f'{PINV}/b0_dprof.json'; av[i + 1] = N1PIN\n"
   "    rep = list(av); rep[o + 1] = f'{REP}/asm_{nm}.json'\n", 1, "f12: the pinned asm argv points --dprof/--chart at the shipped pins"),
  ("assert J(UPIN)['args']['argv'] == ['assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', 'U_cert_v35.json']",
   "assert [NORMP(a) for a in J(UPIN)['args']['argv']] == [NORMP(a) for a in ['assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', 'U_cert_v35.json']]", 1, "f12: argv assertion modulo install root"),
  ("assert J(ZPIN)['args']['argv'] == ['zone3.py', '--ucert', '../t4_cert/out/U_cert_for_zone3_v35.json', '--ucert_full', UPIN]",
   "assert [NORMP(a) for a in J(ZPIN)['args']['argv']] == [NORMP(a) for a in ['zone3.py', '--ucert', '../t4_cert/out/U_cert_for_zone3_v35.json', '--ucert_full', UPIN]]", 1, "f12: argv assertion modulo install root"),
  ("g = dict(__name__='n35_gen_f12', _A=Apath, _A4=A4path, _U=Upath)", "g = dict(__name__='n35_gen_f12', __file__=GENPY, _A=Apath, _A4=A4path, _U=Upath)", 1, "f12: exec of n35_gen_numbers.py gets __file__ (its KIT_ROOT header)"),
  ("    if k not in FILES: FILES[k] = J(f['path'])", "    if k not in FILES: FILES[k] = J(_KP(f['path']))", 1, "f12: files named in constants_needed_v34.json are read from the kit (pins for step-4 v33 outputs)"),
 ],
}
for f in stage_scripts + ['final_tail_v33.py']:
    if f == 'final_tail_v33.py': continue
    n = patch_py(V33D + '/' + f, 'step4/p2/v33/' + f, PRE.get(f, ()))

# ---------------------------------------------------------------------------------------------------------------------------------- the kit driver (final_tail_v33.py)
DPRE = [
 ("V32MAN = f'{O}/FINAL_manifest.json'", "V32MAN = '/home/claude/s1spec/reference/final_manifest.json'; PINV = '/home/claude/s1spec/reference/v33_out/p2/out_v33'; REFM33 = '/home/claude/s1spec/reference/v33_manifest.json'",
  1, "driver: V32MAN is the slim reference manifest; PINV = shipped pins"),
 ("ab = lambda p: p if p.startswith('/') else os.path.join(S4, p)",
  "ab = lambda p: os.path.join(_S1, p[2:]) if p.startswith('@/') else p if os.path.isabs(p) else os.path.join(S4, p)", 1, "driver: ab() understands '@/' (kit-root relative) paths"),
 ("rel = lambda p: os.path.relpath(p, S4) if p.startswith(S4 + '/') else p",
  "def rel(p):\n    p = p.replace('\\\\', '/')\n    if p.startswith(S4 + '/'): return os.path.relpath(p, S4).replace('\\\\', '/')\n"
  "    return '@/' + os.path.relpath(p, _S1).replace('\\\\', '/') if p.startswith(_S1 + '/') else p", 1, "driver: rel() records kit-root paths as '@/...' (same normal form as reference/v33_manifest.json)"),
 ("BASE_ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}",
  "BASE_ENV = {'PATH': os.pathsep.join([os.path.dirname(sys.executable)] + (['/usr/bin', '/bin'] if os.name != 'nt' else [os.path.join(os.environ.get('SystemRoot', 'C:\\\\Windows'), 'System32')])),\n"
  "            'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}\n"
  "if os.name == 'nt': BASE_ENV.update({k: os.environ[k] for k in ('SystemRoot', 'SYSTEMROOT', 'TEMP', 'TMP') if k in os.environ})   # [kit patch] Windows needs these in a clean child environment [TODO: untested on Windows]",
  1, "driver: portable explicit child PATH (os.pathsep; Windows SystemRoot/TEMP)"),
 ("P['cmd'] == ['python3', '-I', '-B'] + S['cmd']", "P['cmd'] == ['python3', '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver: interpreter flags + -X utf8"),
 ("cmd=['python3', '-I', '-B'] + S['cmd']", "cmd=['python3', '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver: interpreter flags + -X utf8"),
 ("subprocess.run([sys.executable, '-I', '-B'] + S['cmd']", "subprocess.run([sys.executable, '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver: interpreter flags + -X utf8"),
 ("flags=['-I', '-B']", "flags=['-I', '-B', '-X', 'utf8']", 1, "driver: recorded flags"),
 ("F12_IN = ([f'{V}/f12_sliver.py', f'{V}/n1_chart.py', f'{OV}/n1_chart.json',", "F12_IN = ([f'{V}/f12_sliver.py', f'{V}/n1_chart.py', f'{PINV}/n1_chart.json',", 1, "driver: f12_sliver reads the shipped PINNED n1_chart.json"),
 ("[f'{T2V}/out/A_interface.json', f'{OV}/b0_dprof.json', f'{OV}/h2u_v33_arb.json',", "[f'{T2V}/out/A_interface.json', f'{PINV}/b0_dprof.json', f'{PINV}/h2u_v33_arb.json', f'{PINV}/FINAL_manifest_v33.json',", 1,
  "driver: f12_sliver reads the shipped PINNED b0_dprof.json, h2u_v33_arb.json and the reference manifest"),
 ("if __name__ == '__main__':\n    v32_start, V32LAST, V32UP = verify_v32();",
  "if __name__ == '__main__':\n    if '--list' in sys.argv:                                                  # [kit patch] list the stages and exit\n"
  "        for i_, S_ in enumerate(STAGES): print(f\"{i_:2d}  {S_['stage']}\")\n        sys.exit(0)\n    v32_start, V32LAST, V32UP = verify_v32();", 1, "driver: --list"),
 ("    if '--from' in sys.argv:\n        i0 = stage_index(", "    if '--from' in sys.argv and stage_index(sys.argv[sys.argv.index('--from') + 1]) > 0:   # [kit patch] --from 0 starts a fresh manifest (no previous one exists in the kit)\n        i0 = stage_index(", 1,
  "driver: --from 0 = fresh start"),
 ("            if pinned is not None: assert pinned == h, ('v3.2 pinned input changed', rp)\n", "            # [kit patch] the v3.2 chain is regenerated in the kit (numerically equal, not byte-equal): pins are recorded, not asserted\n", 1,
  "driver: v3.2 pins recorded, not asserted"),
 ("    man['reused_v32'] = sorted({i['path']: i['sha256'] for P in man['stages'] for i in P['inputs'] if i['v32_pinned']}.items())\n",
  "    man['reused_v32'] = sorted({i['path']: i['sha256'] for P in man['stages'] for i in P['inputs'] if i['v32_pinned']}.items())\n"
  "    man['checks']['v32_pins_differing_from_proof_tree'] = sorted(p_ for p_, h_ in man['reused_v32'] if V32LAST.get(p_) != h_)   # [kit patch] informational: numerically equal, not byte-equal\n", 1,
  "driver: lists the v3.2 inputs whose sha256 differs from the proof tree's"),
]
s = open(V33D + '/final_tail_v33.py', encoding='utf-8').read()
i0_ = s.index('def verify_v32():'); i1_ = s.index('\nAMP = ')
NEWV = '''def verify_v32():
    """[kit patch] the proof tree verified its closed v3.2 manifest (34 stages, 523 paths) byte for byte.  This kit REGENERATES the v3.2 chain (numerically equal, not byte-equal), so
    that check is not applicable; instead the v3.2 files that the v3.3 stages read are listed with the sha256 the proof tree had (reference/v33_manifest.json) and the stages record
    whether this run's file differs (checks.v32_pins_differing_from_proof_tree).  A missing input is still a hard error (assert in the stage loop)."""
    rm = json.load(open(REFM33)); pins = {}
    for P in rm['stages']:
        for f in P['inputs']:
            if f['v32_pinned']: pins.setdefault(f['path'], f['sha256'])
    return dict(manifest=rel(V32MAN), manifest_sha256=sha(V32MAN), stages=rm['v32_manifest']['stages'], paths=len(pins), upstream=0, scripts=0, mismatches=[], ok=True,
                kit_patch='v3.2 manifest verification replaced by the pin list of reference/v33_manifest.json'), pins, {}
'''
s = s[:i0_] + NEWV + s[i1_:]; PATCHES.append(('step4/p2/v33/final_tail_v33.py', 'driver: verify_v32() replaced (pin list instead of byte verification of the closed v3.2 manifest)'))
open(V33D + '/final_tail_v33.py', 'w', encoding='utf-8').write(s)
n = patch_py(V33D + '/final_tail_v33.py', 'step4/p2/v33/final_tail_v33.py', DPRE)

# ---------------------------------------------------------------------------------------------------------------------------------- bulk inputs (f12_sliver, N1, zone3, ...)
EXC_DIR = {'probe', 'prov', 'ref_phase2a', 'v3.0_backup', '__pycache__'}
EXC_F = re.compile(r'(_pre_|\.bak|\.pdf$|\.pyc$|\.tex$|\.md$|_pre\.|n35assert|\.partial|\.tmp$|^run_.*\.log$)')
def copy_tree_sel(src, dst):
    n = 0
    for dp, dn, fn in os.walk(src):
        dn[:] = sorted(d for d in dn if d not in EXC_DIR and not d.startswith('out_'))
        for f in sorted(fn):
            if EXC_F.search(f): continue
            cp(os.path.join(dp, f), os.path.join(dst, os.path.relpath(os.path.join(dp, f), src))); n += 1
    return n
T2S = BULK + '/fix2/n12/t2_A_v35'
nT2 = copy_tree_sel(T2S, KB + '/fix2/n12/t2_A_v35')
for f in ('assemble_U.py', 't4lib.py', 'out/u1_ghat.json', 'out/u2_K.json', 'out/u34_primary.json', 'out/U_cert_v35.json', 'out/U_cert_for_zone3_v35.json'): cp(BULK + '/t4_cert/' + f, KB + '/t4_cert/' + f)
for f in ('zone3.py', 'z3lib.py', 'out/ramp_U.json', 'out/zone3_interface.json'): cp(BULK + '/t5_zone3/' + f, KB + '/t5_zone3/' + f)
cp(BULK + '/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json', KB + '/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json')
cp(BULK + '/fix1/lemmaU_rev/constants_needed_v34.json', KB + '/fix1/lemmaU_rev/constants_needed_v34.json')
cp(BULK + '/fix2/n35_gen_numbers.py', KB + '/fix2/n35_gen_numbers.py')
for f in ('assembly_arb.json', 'assembly_arb.out'): cp(S + '/writeup/lemmaS/rigor/' + f, K + '/writeup/lemmaS/rigor/' + f)
cp(S + '/writeup/lemmaO_v32/lemmaR_consts.json', K + '/writeup/lemmaO_v32/lemmaR_consts.json')
V33V = S + '/step3nl_v33'
for f in ('verifier/constants.py', 'verifier/rate_table_fin12f4_a0.00015.json', 'verifier/modsups_fin12f4/final_alpha_points_5e-8.json', 'agent_sob/sob_constants_P2.json'): cp(V33V + '/' + f, K + '/step3nl_v33/' + f)
# patch the .py of the bulk trees
BPRE = {'writeup/lemmaU/bulk/fix2/n35_gen_numbers.py': [
    ("H=J('/home/claude/s1spec/step4/p2/out_v33/h2u_v33_arb.json'); D=J('/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json')",
     "H=J('/home/claude/s1spec/reference/v33_out/p2/out_v33/h2u_v33_arb.json'); D=J('/home/claude/s1spec/reference/v33_out/p2/out_v33/b0_dprof.json')   # [kit patch] the PINNED step-4 v33 outputs", 1,
     "n35_gen_numbers: reads the shipped pinned h2u_v33_arb.json / b0_dprof.json")]}
nb = 0
for dp, dn, fn in os.walk(K + '/writeup'):
    for f in fn:
        if not f.endswith('.py'): continue
        p = os.path.join(dp, f); rel_ = os.path.relpath(p, K).replace(os.sep, '/')
        if '/home/claude' in open(p, encoding='utf-8').read() or rel_ in BPRE: nb += patch_py(p, rel_, BPRE.get(rel_, ())) is not None
print('bulk: %d files of t2_A_v35, patched %d .py files under writeup/' % (nT2, nb))
# SDS: the kit-step3 scripts the v33 stages import (c4fun, psi3fun, h2fun, hstail, lintail2, s1_cap/tail) are already the patched kit_step3 copies (step3nl, step3, lin, s1_cap)
for f in ('step3nl/agent_int/v3/c4fun.py', 'step3nl/psi3fun.py', 'step3nl/h2fun.py', 'step3/hstail.py', 'step3/hs_uforms.pkl', 'lin/lintail2.py', 's1_cap/tail.py'): assert os.path.exists(K + '/' + f), f

# ---------------------------------------------------------------------------------------------------------------------------------- reference data of the v33 chain
RV = K + '/reference/v33_out'
outs33 = []
for st_ in FM33['stages']:
    for o in st_['outputs']:
        if not o['path'].endswith('.log'): outs33.append(o['path'])
for p_ in outs33: cp(S4 + '/' + p_, RV + '/' + p_)
cp(FM33P, RV + '/p2/out_v33/FINAL_manifest_v33.json')                      # also the f12 pin 'S4' of constants_needed_v34.json
slim = dict(note='slim copy of p2/out_v33/FINAL_manifest_v33.json of the proof tree (sha256 %s): the 21 stages with their inputs/outputs (sha256), results and checks.  Paths: relative to step4, or @/ = kit root.' % FM33_SHA,
            manifest_sha256=FM33_SHA, v32_manifest=dict(FM33['v32_manifest'], stages=34, paths=523), script_sha256={normrel(k): v for k, v in FM33['script_sha256'].items()},
            results=FM33['results'], checks={k: v for k, v in FM33['checks'].items() if k in ('all_items_ok', 'hashes_verified_at_end', 'inputs_with_two_versions')}, stages=[])
for st_ in FM33['stages']:
    slim['stages'].append(dict(stage=st_['stage'], cwd=normrel(st_['cwd']), cmd=st_['cmd'], env=st_['env'], seconds=st_['seconds'],
                               inputs=[dict(path=normrel(i['path']), sha256=i['sha256'], v32_pinned=i['v32_pinned']) for i in st_['inputs']],
                               outputs=[dict(path=normrel(o['path']), sha256=o['sha256']) for o in st_['outputs']]))
json.dump(slim, open(K + '/reference/v33_manifest.json', 'w'), indent=1)
print('v33 reference: %d output files, slim manifest with %d stages' % (len(outs33), len(slim['stages'])))

import glob
exec(open(BD + '/v38_block.py', encoding='utf-8').read())                     # [v4] the v38 stages
# --- docstrings (not comments: the [kit patch] header keeps its provenance comment) of the files inherited from kit B that still name the proof tree's path: neutral wording
nds = 0
for dp, dn, fs_ in os.walk(K):
    for f_ in fs_:
        if not f_.endswith('.py'): continue
        p_ = os.path.join(dp, f_); t_ = open(p_, encoding='utf-8').read()
        if '/home/claude' not in t_: continue
        t2_ = sanitize_nonexec(t_, comments=False)
        if t2_ != t_: ast.parse(t2_); open(p_, 'w', encoding='utf-8').write(t2_); nds += 1
print('docstrings neutralised in', nds, 'inherited files')
man = []
for root in ('step4/t3/p1/ckpt', 's1_tight', 'step3/outA5.txt', 'step3/hs_uforms.pkl', 'step3/lyap_consts2.json', 'step3nl/agent_mod/lemma/lyap_extra.json',
             'step4/t7/phi_store.json', 'step4/t7/fixed_coarse.json', 'step4/p2/h6/fable/h6f_cases.json', 'writeup', 'reference/h1_pre_u4',
             'step3nl_v33', 'reference/v33_out', 'step3nl_v4c', 'reference/v38_out', 'step4/p2/v38/EXCEPTION_v4c_PA_reg_e2v4b.json', 'step3nl/agent_int/v3/c4fun.py', 'step3nl/psi3fun.py', 'step3nl/h2fun.py', 'step3/hstail.py', 'lin/lintail2.py', 's1_cap/tail.py'):
    fp = os.path.join(K, root)
    fl = [os.path.join(dp, f) for dp, _, fs in os.walk(fp) for f in fs] if os.path.isdir(fp) else [fp]
    for f in sorted(fl):
        h = hashlib.sha256(open(f, 'rb').read()).hexdigest(); man.append(f'{h}  ' + os.path.relpath(f, K).replace(os.sep, '/'))
os.makedirs(K + '/inputs', exist_ok=True)
open(K + '/inputs/MANIFEST.sha256', 'w').write('# SHA-256 of the shipped data inputs (cannot be regenerated on Windows / belong to earlier steps)\n' + '\n'.join(man) + '\n')
# --- reference data from the final manifest (FINAL_manifest.json + A1_manifest.json): stage list for run_all.py, output hashes for check.py
import hashlib
O_ = S4 + '/p2/out'; FM = json.load(open(_FMP)); A1 = json.load(open(O_ + '/A1_manifest.json'))
pas = [i for i, st_ in enumerate(FM['stages']) if st_['pass_A']]; assert pas == list(range(pas[0], pas[0] + len(pas))), pas
npass = len(pas); stages = []
for i, st_ in enumerate(FM['stages']):
    g = 'pre' if i < pas[0] else 'A' if i < pas[0] + npass else 'B' if i < pas[0] + 2 * npass else 'post'
    assert st_['cmd'][:2] == ['python3', '-I']
    stages.append(dict(stage=st_['stage'], group=g, cwd=st_['cwd'], cmd=st_['cmd'][2:], env=st_['env'], log=st_['outputs'][-1]['path'], outputs=[o['path'] for o in st_['outputs'][:-1]]))
outs = {}
for st_ in FM['stages']:
    for o in st_['outputs']:
        if not o['path'].endswith('.log'): outs[o['path']] = o['sha256']
    for i_ in st_['inputs']:
        if i_['path'].startswith(('p2/t6c/', 'p2/out/h1slab/', 'p2/out/t9_h1')) or i_['path'] in ('t7/constants.json', 't7/run_v2/out_farg/farg_002.json'): outs.setdefault(i_['path'], i_['sha256'])
for f_ in A1['files']:
    if f_['path'].startswith('p2/out/h1slab/') and f_['path'].endswith('.json'): outs.setdefault(f_['path'], f_['sha256'])
os.makedirs(K + '/reference', exist_ok=True)
json.dump(dict(note='slim copy of p2/out/FINAL_manifest.json (the final tail of the proof tree); run_all.py executes these stages', stages=stages,
               h1env=dict(S1TT='10.249077651120778,10.249077651121361'), script_sha256=FM['script_sha256'], checks=FM['checks']), open(K + '/reference/final_manifest.json', 'w'), indent=1)
json.dump(outs, open(K + '/reference/outputs.json', 'w'), indent=1)
for rel_ in ('p2/out/h4loc_A1_w180.json',):
    pass
print('reference: %d tail stages (%d per pass), %d reference output hashes' % (len(stages), npass, len(outs)))
sys.path.insert(0, BD)
import step4_claims as SC
R_, I_ = SC.extract(S4); ours = {n: v for n, v, _, _ in R_}; ours.update({'info: ' + n: v for n, v in I_})
json.dump(ours, open(K + '/reference/ours.json', 'w'), indent=1)
print('ours.json: %d values; claims failing on the proof tree itself: %s' % (len(R_), [n for n, v, op, b in R_ if v is None or b is None or not ((v <= b) if op == '<=' else (v >= b))]))
for f in ('run_all.py', 'check.py', 'README.md', 'requirements.txt', 'step4_claims.py', 'verify_outputs.py', 'verify_v33.py'):
    cp(os.path.join(BD, f), K + '/' + f)
cp(S + '/kit_common_v4/kit_guard.py', K + '/kit_guard.py')   # [v4] write guard


# ====================================================================================================================================================
# post: audits, deterministic zip, BUILD_REPORT.md
# ====================================================================================================================================================
for dp, dn, fn in os.walk(K):
    for d in list(dn):
        if d == '__pycache__': shutil.rmtree(os.path.join(dp, d)); dn.remove(d)
if os.path.exists(K + '/kit_out'): shutil.rmtree(K + '/kit_out')
ALL = sorted(os.path.relpath(os.path.join(dp, f), K).replace(os.sep, '/') for dp, _, fs in os.walk(K) for f in fs)
# --- audit 1: /home/claude in the kit (code string tokens vs comments/docstrings vs pinned data)
aud = dict(code=[], comment_or_doc=[], data=[], binary=0)
for rp in ALL:
    p = K + '/' + rp
    try: t = open(p, encoding='utf-8').read()
    except UnicodeDecodeError: aud['binary'] += 1; continue
    if '/home/claude' not in t: continue
    if rp.endswith('.py'):
        try:
            hits = [(tk.type, tk.string) for tk in tokenize.generate_tokens(io.StringIO(t).readline) if '/home/claude' in tk.string]
        except Exception: hits = [(0, '?')]
        for ty, st_ in hits: (aud['comment_or_doc'] if ty == tokenize.COMMENT else aud['code']).append((rp, st_[:90].replace('\n', ' ')))
    else: aud['data'].append(rp)
code_ok = {'run_all.py': 'tr(): rewrite table for provenance strings of the proof tree manifest', 'verify_outputs.py': 'strip(): provenance prefix of reference JSON'}
code_bad = [(rp, st_) for rp, st_ in aud['code'] if rp not in code_ok]
print('audit /home/claude: code-string hits outside the 2 allowed rewrite tables: %d; comment hits: %d (HDR convention); data files: %d' % (len(code_bad), len(aud['comment_or_doc']), len(aud['data'])))
for x in code_bad[:20]: print('   CODE HIT', x)
# --- audit 2: local imports of every .py under step4/p2/v33 and writeup/ resolve inside the kit
std = set(sys.stdlib_module_names) | {'flint', 'sympy', 'numpy', 'scipy', 'mpmath'}
miss_imp = []
for rp in ALL:
    if not (rp.endswith('.py') and (rp.startswith('step4/p2/v33/') or rp.startswith('writeup/lemmaU/') or rp.startswith('step3nl_v33/'))): continue
    tree = ast.parse(open(K + '/' + rp, encoding='utf-8').read()); here = os.path.dirname(rp)
    cand_dirs = [here, 'writeup/lemmaU/bulk/fix2/n12/t2_A_v35/base_ibp', 'writeup/lemmaU/bulk/fix2/n12/t2_A_v35/base_ibp/rigor', 'writeup/lemmaU/bulk/fix2/n12/t2_A_v35/rigor', 'step4/p2/t9h', 'step4/p2', 'step4/t6',
                 'step3nl/agent_int/v3', 'step3nl', 'step3', 'lin', 's1_cap', 'step4/p2/t8se', 'step4/p2/v33']
    for nd in ast.walk(tree):
        names = [a.name.split('.')[0] for a in nd.names] if isinstance(nd, ast.Import) else [nd.module.split('.')[0]] if isinstance(nd, ast.ImportFrom) and nd.module and nd.level == 0 else []
        for nm in names:
            if nm in std or nm.startswith('_k'): continue
            if not any(os.path.exists(K + '/' + d + '/' + nm + '.py') for d in cand_dirs): miss_imp.append((rp, nm))
print('audit imports: unresolved local imports: %s' % sorted(set(miss_imp))[:20])
# --- audit 3: every input of the reference v33 stages is shipped, produced by an earlier stage, or a v3.2-chain output of the WSL test run (read-only look at kit4b_run_wsl1)
prod = set(); unresolved = []; chain = []; WSL = S + '/kit4b_run_wsl1'
for st_ in slim['stages']:
    for i_ in st_['inputs']:
        p_ = i_['path']
        if p_ in prod: continue
        if os.path.exists(kitpath(p_)): continue
        if p_.startswith('@/'): unresolved.append((st_['stage'][:30], p_)); continue
        (chain if os.path.exists(WSL + '/step4/' + p_) else unresolved).append((st_['stage'][:30], p_) if not os.path.exists(WSL + '/step4/' + p_) else p_)
    for o_ in st_['outputs']: prod.add(o_['path'])
print('audit inputs: %d kit-chain (v3.2) inputs found in the WSL test run tree, unresolved: %s' % (len(set(chain)), unresolved[:10]))
# --- audit 4: path literals rewritten to the kit (existence; outputs and logs excluded)
lit = re.compile(r"""_S1 \+ ['"](/[^'"]*)['"]|\{_S1\}(/[^'"{}]*)""")
absent = {}
for rp in ALL:
    if not (rp.endswith('.py') and (rp.startswith(('step4/p2/v33/', 'writeup/lemmaU/', 'step3nl_v33/')))): continue
    for m in lit.finditer(open(K + '/' + rp, encoding='utf-8').read()):
        q = m.group(1) or m.group(2)
        if '*' in q or q.endswith('/') and False: continue
        if not os.path.exists(K + q): absent.setdefault(q, set()).add(rp)
print('audit literals: %d kit-path literals do not exist in the tree (outputs of the chain expected):' % len(absent))
for q, v in sorted(absent.items()): print('   ', q, '<-', sorted(v)[:2])
AUDIT = dict(code_hits_bad=code_bad, comment_hits=len(aud['comment_or_doc']), data_files_with_provenance=len(aud['data']), unresolved_imports=sorted(set(miss_imp)), unresolved_inputs=unresolved,
             chain_inputs=len(set(chain)), absent_literals={q: sorted(v) for q, v in absent.items()}, code_allowed=code_ok, code_all=aud['code'])
json.dump(AUDIT, open(BD + '/audit.json', 'w'), indent=1)
# --- deterministic zip
if os.path.exists(ZIP): os.remove(ZIP)
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for rp in ALL:
        zi = zipfile.ZipInfo(rp, date_time=(2026, 11, 20, 0, 0, 0)); zi.compress_type = zipfile.ZIP_DEFLATED; zi.external_attr = 0o644 << 16
        with open(K + '/' + rp, 'rb') as fh: z.writestr(zi, fh.read(), compresslevel=9)
zsha = hashlib.sha256(open(ZIP, 'rb').read()).hexdigest()
print('zip', os.path.getsize(ZIP), zsha)
# --- report data for BUILD_REPORT.md
FILES = {rp: (hashlib.sha256(open(K + '/' + rp, 'rb').read()).hexdigest(), os.path.getsize(K + '/' + rp)) for rp in ALL}
json.dump(dict(files=FILES, patches=PATCHES, zip=dict(path=ZIP, sha256=zsha, size=os.path.getsize(ZIP)), nfiles=len(ALL)), open(BD + '/build_data.json', 'w'), indent=1)
print('files', len(ALL), 'patches', len(PATCHES))

# --- relation to kit B (kit_step4b, read-only)
V1 = S + '/kit_step4b'
v1 = {os.path.relpath(os.path.join(dp, f), V1).replace(os.sep, '/') for dp, _, fs in os.walk(V1) for f in fs if '__pycache__' not in dp}
only_new = sorted(set(ALL) - v1); only_old = sorted(v1 - set(ALL))
changed = sorted(rp for rp in set(ALL) & v1 if hashlib.sha256(open(K + '/' + rp, 'rb').read()).hexdigest() != hashlib.sha256(open(V1 + '/' + rp, 'rb').read()).hexdigest())
print('vs kit B: %d files only in v3, %d only in kit B, %d changed' % (len(only_new), len(only_old), len(changed)))
RELDATA = dict(only_new=only_new, only_old=only_old, changed=changed)
json.dump(RELDATA, open(BD + '/vs_kitB.json', 'w'), indent=1)

exec(open(BD + '/report_v4.py', encoding='utf-8').read())   # [v4] BUILD_REPORT.md
