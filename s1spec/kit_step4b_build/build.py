"""assemble kit_step4b from the proof tree (not shipped)"""
import os, re, shutil, sys, ast, json
S = '/home/claude/s1spec'; K = S + '/kit_step4b'
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
cp(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'h1slab_summaries.py'), D4 + '/p2/t9h/h1slab_summaries.py')
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
import hashlib
man = []
for root in ('step4/t3/p1/ckpt', 's1_tight', 'step3/outA5.txt', 'step3/hs_uforms.pkl', 'step3/lyap_consts2.json', 'step3nl/agent_mod/lemma/lyap_extra.json',
             'step4/t7/phi_store.json', 'step4/t7/fixed_coarse.json', 'step4/p2/h6/fable/h6f_cases.json', 'writeup', 'reference/h1_pre_u4'):
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step4_claims as SC
R_, I_ = SC.extract(S4); ours = {n: v for n, v, _, _ in R_}; ours.update({'info: ' + n: v for n, v in I_})
json.dump(ours, open(K + '/reference/ours.json', 'w'), indent=1)
print('ours.json: %d values; claims failing on the proof tree itself: %s' % (len(R_), [n for n, v, op, b in R_ if v is None or b is None or not ((v <= b) if op == '<=' else (v >= b))]))
for f in ('run_all.py', 'check.py', 'README.md', 'requirements.txt', 'step4_claims.py', 'verify_outputs.py'):
    if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), f)): cp(os.path.join(os.path.dirname(os.path.abspath(__file__)), f), K + '/' + f)

