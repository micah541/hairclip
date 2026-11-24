#!/usr/bin/env python3
"""Builder for the reproduction kit "kit_step3 v4" (step-3 certificates on verifier fin13v4c, (A) chain v37, flux v37, interface / proof checks).

  python3 build.py

Creates /home/claude/s1spec/kit_step3_v4/ (wiped first), /home/claude/s1spec/kit_step3_v4.zip (deterministic) and kit_step3_v4_build/BUILD_REPORT.md.
Reads the proof tree read-only; runs no certified computation (hashing, ast.parse, zip only).  The objects are those pinned by the step-4 v38 manifest
step4/p2/out_v38/FINAL_manifest_v38.json (approved by Micah on Nov 20 2026): verifier fin13v4c production run step3nl_v4c/runs/P (NOT runs/P_A, NOT
runs/reg_fin13v4b), (A) chain v37 writeup/lemmaU/bulk/fix13_v37chain (runs/prod2), flux v37 (flux/out, flux/texts), C-1/GS/G1/U-1 checks.

Layout: the kit mirrors the proof tree (paths relative to /home/claude/s1spec; /home/claude/s1_tight -> s1_tight, /home/claude/s1_cap -> s1_cap).  Every
shipped .py gets the `_kit_root()` header and every absolute path literal of the proof tree becomes `(_S1 + '/rel')` ([kit patch], generic, listed in
BUILD_REPORT.md); comments and docstrings get '<proof-tree>'.  Support trees step3, step3nl, lin, s1_cap, cap, s1_tight, writeup/lemmaO and the ibp/lemmaO
parts are taken from kit_step3 v3 (already patched, run on Windows in kit v2); everything of the v4 objects is taken fresh from the proof tree.
Every process the kit starts goes through kit_guard.py (writes outside the kit raise PermissionError; Nov 15 incident)."""
import os, re, sys, ast, json, shutil, hashlib, zipfile, tokenize, io, time
sys.dont_write_bytecode = True
S = '/home/claude/s1spec'; B = S + '/writeup/lemmaU/bulk'; NC = B + '/fix13_v37chain'; VC = S + '/step3nl_v4c'
OLD = S + '/kit_step3_v3'; HERE = os.path.dirname(os.path.abspath(__file__)); K = S + '/kit_step3_v4'; ZIP = S + '/kit_step3_v4.zip'
M38 = S + '/step4/p2/out_v38/FINAL_manifest_v38.json'
REPORT = dict(patches=[], notes=[], todo=[], asserts=[])
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
def note(s): REPORT['notes'].append(s)
def todo(s): REPORT['todo'].append(s)

# ---------------------------------------------------------------------------------------------------- preconditions: the objects are the pinned ones
man = json.load(open(M38)); PIN = {}
for P in man['stages']:
    for f in P['inputs'] + P['outputs']: PIN.setdefault(f['path'] if f['path'].startswith('/') else S + '/step4/' + f['path'], f['sha256'])
def pinned(p):
    """p must be pinned by FINAL_manifest_v38.json with its current sha256 (or, for files the v38 stages do not read, is only hashed)"""
    h = sha(p)
    if p in PIN: assert PIN[p] == h, ('pinned file changed since the v38 manifest', p)
    REPORT['asserts'].append((p.replace(S + '/', ''), h, p in PIN)); return h
for p in [VC + '/runs/P/rate_table_fin13v4c_a0.00015.json', VC + '/runs/P/rate_table_fin13v4c_a0.000125.json', VC + '/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json',
          VC + '/v4bcert/out/strip_cubic_v4c.json', VC + '/v4bcert/out/e0_A2_budget_v4c.json', NC + '/t2_A/out/A_interface.json', NC + '/t4_cert/out/U_cert_v36_pass3.json',
          NC + '/t5_zone3/out/zone3_interface.json', NC + '/n1/out/n1_chart_fin13v4c.json']:
    pinned(p)
if os.path.exists(K): shutil.rmtree(K)
os.makedirs(K)

# ---------------------------------------------------------------------------------------------------- patching
HDR = ("import os as _kos\n"
       "def _kit_root(_d=_kos.path.dirname(_kos.path.abspath(__file__))):\n"
       "    while not _kos.path.exists(_kos.path.join(_d, 'KIT_ROOT')) and _d != _kos.path.dirname(_d): _d = _kos.path.dirname(_d)\n"
       "    return _d\n"
       "_S1 = _kit_root().replace('\\\\', '/')   # [kit patch] the kit root (was an absolute path of the proof tree)\n")
PAT = re.compile(r"""(f?)(['"])/home/claude/(s1spec|s1_tight|s1_cap)(/[^'"\n]*|)\2""")
def _repl(m):
    f, q, root, rest = m.groups(); sub = {'s1spec': '', 's1_tight': '/s1_tight', 's1_cap': '/s1_cap'}[root]
    return f"f{q}{{_S1}}{sub}{rest}{q}" if f else f"(_S1 + {q}{sub}{rest}{q})"
def _nonexec(s):
    """comments / docstrings naming the proof tree's path get '<proof-tree>' (not executed)"""
    lines = s.splitlines(keepends=True); starts = [0]
    for l in lines: starts.append(starts[-1] + len(l))
    toks = list(tokenize.generate_tokens(io.StringIO(s).readline)); edits = []
    for i, t in enumerate(toks):
        if '/home/claude' not in t.string: continue
        ok = t.type == tokenize.COMMENT
        if t.type == tokenize.STRING:
            j = i - 1
            while j >= 0 and toks[j].type in (tokenize.NL, tokenize.COMMENT): j -= 1
            ok = (j < 0 or toks[j].type in (tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT)) and toks[i + 1].type == tokenize.NEWLINE
        if not ok: sys.exit('absolute path left in CODE: %s' % t.string[:100])
        edits.append((starts[t.start[0] - 1] + t.start[1], starts[t.end[0] - 1] + t.end[1], t.string))
    for a, b, txt in sorted(edits, reverse=True): s = s[:a] + txt.replace('/home/claude/s1spec', '<proof-tree>').replace('/home/claude', '<proof-tree-host>') + s[b:]
    return s
def patch_text(s, rel, pre=()):
    for pr in pre:
        a, b_ = pr[0], pr[1]
        assert s.count(a) == 1, ('kit-patch target', rel, s.count(a), a[:90]); s = s.replace(a, b_); REPORT['patches'].append((rel, 'manual', a.strip()[:90], b_.strip()[:110]))
    n0 = len(PAT.findall(s)); s = PAT.sub(_repl, s)
    if n0: REPORT['patches'].append((rel, 'generic', '%d absolute path literal(s)' % n0, "(_S1 + '/...')"))
    s = _nonexec(s)
    if n0 and 'def _kit_root' not in s:                                       # header only where a literal was rewritten (files without one stay byte-identical)
        tree = ast.parse(s); ls = s.split('\n'); at = 0
        if tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(getattr(tree.body[0], 'value', None), ast.Constant) and isinstance(tree.body[0].value.value, str):
            at = tree.body[0].end_lineno
        else:
            while at < len(ls) and (ls[at].startswith('#') or ls[at].startswith('from __future__')): at += 1
        ls[at:at] = HDR.rstrip('\n').split('\n'); s = '\n'.join(ls)
    ast.parse(s); return s
def cp(src, rel, pre=()):
    dst = os.path.join(K, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
    if src.endswith('.py'):
        with open(dst, 'w', encoding='utf-8', newline='\n') as f: f.write(patch_text(open(src, encoding='utf-8').read(), rel, pre))
    else: shutil.copyfile(src, dst)
    return dst
def rel_of(p):
    for a, b_ in (('/home/claude/s1spec/', ''), ('/home/claude/s1_tight/', 's1_tight/'), ('/home/claude/s1_cap/', 's1_cap/')):
        if p.startswith(a): return b_ + p[len(a):]
    raise ValueError(p)
def cpm(src, pre=()): return cp(src, rel_of(src), pre)          # mirror position
EXC_F = re.compile(r'(_pre_|\.bak|\.pdf$|\.pyc$|\.tex$|\.md$|_pre\.|n35assert|\.partial|\.tmp$|\.sh$|\.time$|\.out$)')
def cptree(src, skipdirs=(), skipfiles=()):
    n = 0
    for dp, dn, fn in os.walk(src):
        dn[:] = sorted(d for d in dn if d not in ('__pycache__',) + tuple(skipdirs) and os.path.relpath(os.path.join(dp, d), src) not in skipdirs)
        for f in sorted(fn):
            p = os.path.join(dp, f)
            if EXC_F.search(f) or os.path.relpath(p, src) in skipfiles: continue
            cpm(p); n += 1
    return n

# ---------------------------------------------------------------------------------------------------- 1. support trees from kit_step3 v3 (patched, tested)
for d in ('step3', 'step3nl', 'lin', 's1_cap', 'cap', 's1_tight', 'writeup/lemmaO', 'writeup/lemmaS', 'writeup/lemmaU/phase1', 'writeup/lemmaU/phase2a', 'writeup/lemmaO_v32',
          'writeup/lemmaU/bulk/t1_statements', 'reference/ibp', 'reference/lemmaO'):
    shutil.copytree(OLD + '/' + d, K + '/' + d, ignore=shutil.ignore_patterns('__pycache__'))
shutil.copyfile(OLD + '/KIT_ROOT', K + '/KIT_ROOT')
shutil.rmtree(K + '/step3nl/verifier')                                     # the v3.3 (fin12f4) verifier copy of kit v2/v3 is not used by v4
for f in ('step3nl/graph_eq.pkl',): cpm(S + '/' + f)                       # missing in kit v3 (lemmaS_symbolic of base_ibp)
for f in ('/home/claude/s1_tight/prof/prof_tight.txt', '/home/claude/s1_tight/ext813/ext813_prof.txt'): cpm(f)
note('support trees step3, step3nl (minus its v3.3 verifier), lin, s1_cap, cap, s1_tight, writeup/{lemmaO, lemmaS, lemmaO_v32, lemmaU/phase1, lemmaU/phase2a, lemmaU/bulk/t1_statements}, reference/{ibp, lemmaO} are the patched copies of kit_step3 v3 (unchanged); added step3nl/graph_eq.pkl, s1_tight/prof/prof_tight.txt, s1_tight/ext813/ext813_prof.txt.')

# ---------------------------------------------------------------------------------------------------- 2. verifier fin13v4c (code, cert inputs) + references
for f in ('bb', 'closure', 'constants', 'ival', 'mod_sups', 'rhs_iv', 'tail_cert', 'unc'): cpm(f'{VC}/verifier/{f}.py')
rt = json.load(open(VC + '/runs/P/rate_table_fin13v4c_a0.00015.json'))
for p, h in rt['args']['inputs'].items():                                  # the 12 pins of the production rate table: code above + the cert inputs
    assert sha(p) == h, ('rate-table pin does not match the disk', p)
    if not p.endswith('.py'): cpm(p)
cpm(VC + '/agent_int/planB/planB_ramp.json')
for f in ('graph_struct.pkl',): cpm(VC + '/' + f)
cpm(VC + '/agent_sob/sob_constants_P2.json'); cpm(VC + '/agent_int/v3/c4fun.py')
for f in os.listdir(VC + '/runs/P'):                                        # the pinned production run (consumed by the chain, as in the proof tree)
    if os.path.isfile(VC + '/runs/P/' + f): cpm(VC + '/runs/P/' + f)
for f in os.listdir(VC + '/runs/P/modsups_fin13v4c'): cpm(VC + '/runs/P/modsups_fin13v4c/' + f)
for f in os.listdir(VC + '/runs/P'):
    if os.path.isfile(VC + '/runs/P/' + f) and not f.endswith('.out'): cp(VC + '/runs/P/' + f, 'reference/verifier/' + f)
cp(VC + '/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json', 'reference/verifier/final_alpha_points_5e-8.json')
# cert stages D-3 (strip cubic) and A-1 (E0 A^2 budget): rerun from the pinned cert outputs
for f in ('strip_cubic_v4c', 'e0_A2_budget_v4c', 'v4bcommon'): cpm(f'{VC}/v4bcert/{f}.py')
for f in sorted(os.listdir(VC + '/v4bcert/out')):
    if f.endswith(('.json', '.pkl')): cpm(f'{VC}/v4bcert/out/{f}')
cpm(VC + '/agent_mod/bloop/bloop.json'); cpm(B + '/fix11_review/d_cubic_strip.json'); cpm(VC + '/v4cert/out/constants_v4.json')
for f in ('strip_cubic_v4c.json', 'e0_A2_budget_v4c.json'): cp(f'{VC}/v4bcert/out/{f}', 'reference/cert/' + f)
for p in json.load(open(VC + '/v4bcert/out/strip_cubic_v4c.json'))['args']['inputs'] if isinstance(json.load(open(VC + '/v4bcert/out/strip_cubic_v4c.json'))['args']['inputs'], dict) else []:
    if p.startswith('/home/claude') and os.path.isfile(p) and not os.path.exists(os.path.join(K, rel_of(p))): cpm(p)

# ---------------------------------------------------------------------------------------------------- 3. (A) chain v37: code + pinned inputs (no in-place outputs)
cpm(NC + '/v36cfg.py'); cpm(NC + '/V37_LOG.md') if False else cp(NC + '/V37_LOG.md', rel_of(NC + '/V37_LOG.md'))   # read by the flux inputs stage (chain report)
cpm(NC + '/n1/n1_chart.py')
for f in sorted(os.listdir(NC + '/r3')):
    if f.endswith('.py') and '_pre_' not in f: cpm(NC + '/r3/' + f)
BASE_OUTS = {'base_ibp/lemmaS_F10p_M0.52.json', 'base_ibp/rigor/cap_F10p.json', 'base_ibp/rigor/lemmaS_rigcap_M0.52.json', 'base_ibp/rigor/weights_exact.json', 'base_ibp/rigor/assembly_arb.json'}
nT2 = cptree(NC + '/t2_A', skipdirs=('out', 'prov', 'ref_phase2a', 'probe'), skipfiles=BASE_OUTS)
for f in sorted(os.listdir(NC + '/t4_cert')):
    if f.endswith('.py') and '_pre_' not in f and 'n35assert' not in f: cpm(NC + '/t4_cert/' + f)
for f in sorted(os.listdir(NC + '/t5_zone3')):
    if f.endswith('.py') and '_pre_' not in f: cpm(NC + '/t5_zone3/' + f)
cpm(NC + '/t5_zone3/out/ramp_U.json')                                       # kept by run_chain_v37.sh (the v35 byte copy, review B3): an input, not regenerated
cpm(B + '/t5_zone3/out/zone3_interface.json')                               # v35 zone3 interface = initial input of assemble pass 1 (run_chain_v37.sh Z3V35)
for p in ('/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json',): cpm(p)
for f in sorted(os.listdir(S + '/step4/p2/t9h')):
    if f.endswith('.py') and '_pre_' not in f: cpm(S + '/step4/p2/t9h/' + f)
note('(A) chain v37: t2_A copied without out/, prov/, ref_phase2a/, probe/ and without the 5 base-chain outputs (regenerated); %d files.' % nT2)
# references of the chain (runs/prod2 = the in-place outputs, byte for byte, checked by the v38 stage 20)
for dp, dn, fn in os.walk(NC + '/runs/prod2'):
    for f in fn:
        p = os.path.join(dp, f); r = os.path.relpath(p, NC + '/runs/prod2')
        if not f.endswith(('.log', 'SHA256SUMS', 'stages.txt', 'env.txt')): cp(p, 'reference/chain/' + r)

# ---------------------------------------------------------------------------------------------------- 4. flux v37
for f in sorted(os.listdir(NC + '/flux/code')):
    if f.endswith('.py') and '_pre_' not in f and f != 'cmp_runs.py': cpm(NC + '/flux/code/' + f)
for f in ('flux_forms.pkl', 'flux_cert.json'): cpm(B + '/fix3/flux/' + f)
cpm(S + '/writeup/lemmaS/lemmaS_cflux.json')
_pin = json.load(open(NC + '/flux/out/inputs.json')).get('pinned') or {}
for p in _pin:                                                              # files pinned by the flux inputs stage
    if isinstance(p, str) and p.startswith('/home/claude') and os.path.isfile(p) and not any(x in p for x in ('/runs/prod2/', '/flux/out/', '/out_v38/', '/step4/p2/out', '/t2_A/out/', '/t4_cert/out/', '/t5_zone3/out/', '/n1/out/')) \
            and not os.path.exists(os.path.join(K, rel_of(p))): cpm(p); note('flux pinned input added (inputs.json pinned): %s' % rel_of(p))
for jf in sorted(os.listdir(NC + '/flux/out')):                          # every recorded input of the pinned flux outputs (args blocks) that the kit does not ship yet
    if not jf.endswith('.json') or jf.endswith('.timing.json'): continue
    inp = (json.load(open(NC + '/flux/out/' + jf)).get('args') or {}).get('inputs') or {}
    for p in (inp if isinstance(inp, dict) else []):
        if not (isinstance(p, str) and p.startswith('/home/claude') and os.path.isfile(p)): continue
        if any(x in p for x in ('/runs/prod2/', '/flux/out/', '/out_v38/', '/step4/', '/t2_A/out/', '/t4_cert/out/', '/t5_zone3/out/', '/n1/out/')): continue
        if not os.path.exists(os.path.join(K, rel_of(p))): cpm(p); note('flux input added from the args block of flux/out/%s: %s' % (jf, rel_of(p)))                                # read by flux2_coef / inputs_manifest (OLDCF regression)
for f in ('flux2_lemma.tex', 'FLUX2_DESIGN.md'): cp(B + '/fix7_flux2/' + f, rel_of(B + '/fix7_flux2/' + f))
for f in sorted(os.listdir(NC + '/flux/out')):
    if f.endswith(('.json', '.pkl')) and not f.endswith('.timing.json'): cp(NC + '/flux/out/' + f, 'reference/flux/out/' + f)
for f in sorted(os.listdir(NC + '/flux/texts')): cp(NC + '/flux/texts/' + f, 'reference/flux/texts/' + f)

# ---------------------------------------------------------------------------------------------------- 5. interface / proof checks (v38 copies; outputs inside the kit)
V38 = S + '/step4/p2/v38'
for f in ('G1_check_v38', 'U1_sym_v38', 'GS_boot_v4c_v38'): cpm(f'{V38}/{f}.py')
for f in ('G1_check', 'U1_sym', 'GS_boot_v4c'): cp(S + f'/step4/p2/out_v38/iface/{f}.json', f'reference/proofs/{f}.json')
cp(S + '/step4/p2/out_v38/FINAL_manifest_v38.json', 'reference/step4/FINAL_manifest_v38.json')

# ---------------------------------------------------------------------------------------------------- 6. drivers, guard, README
for f in ('run_all.py', 'check.py', 'README.md'): shutil.copyfile(os.path.join(HERE, f), os.path.join(K, f))
shutil.copyfile(S + '/kit_common_v4/kit_guard.py', K + '/kit_guard.py')
with open(K + '/requirements.txt', 'w', newline='\n') as f:
    f.write('# versions the reference runs used (CPython 3.14.4, Linux/WSL)\npython-flint==0.9.0\nsympy==1.14.0\nnumpy==2.5.3\nscipy\nmpmath==1.3.0\n')

# ---------------------------------------------------------------------------------------------------- 7. checks: no /home/claude in code; syntax; manifests
bad = []
for dp, dn, fs in os.walk(K):
    for f in fs:
        p = os.path.join(dp, f)
        if f.endswith('.py'):
            t = open(p, encoding='utf-8').read(); ast.parse(t)
            if '/home/claude' in t: bad.append(os.path.relpath(p, K))
if bad: sys.exit('references to /home/claude remain in: %s' % bad)
data_hits = sorted(os.path.relpath(os.path.join(dp, f), K) for dp, dn, fs in os.walk(K) for f in fs if not f.endswith('.py') and b'/home/claude' in open(os.path.join(dp, f), 'rb').read())
note('grep proof: no .py file of the kit contains "/home/claude" (build.py refuses otherwise). %d data files contain it inside provenance fields (args blocks, recorded paths); they are parsed as data only.' % len(data_hits))
mf = []
for dp, dn, fs in os.walk(K):
    dn[:] = sorted(dn)
    for f in sorted(fs):
        p = os.path.join(dp, f); r = os.path.relpath(p, K).replace(os.sep, '/')
        if r.startswith(('reference/', 'inputs/')) or r in ('README.md', 'run_all.py', 'check.py'): continue
        mf.append('%s  %s' % (sha(p), r))
os.makedirs(K + '/inputs', exist_ok=True)
open(K + '/inputs/MANIFEST.sha256', 'w', newline='\n').write('# SHA-256 of every shipped source and data input (kit-patched scripts included)\n' + '\n'.join(mf) + '\n')
rf = sorted('%s  %s' % (sha(os.path.join(dp, f)), os.path.relpath(os.path.join(dp, f), K).replace(os.sep, '/')) for dp, dn, fs in os.walk(K + '/reference') for f in fs)
open(K + '/inputs/REFERENCE.sha256', 'w', newline='\n').write('# SHA-256 of the shipped reference outputs (pinned outputs of the proof tree)\n' + '\n'.join(rf) + '\n')

# ---------------------------------------------------------------------------------------------------- 8. zip + report
names = sorted((os.path.join(dp, f) for dp, dn, fs in os.walk(K) for f in fs), key=lambda p: os.path.relpath(p, K))
if os.path.exists(ZIP): os.remove(ZIP)
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in names:
        zi = zipfile.ZipInfo('kit_step3_v4/' + os.path.relpath(p, K).replace(os.sep, '/'), date_time=(2026, 11, 20, 12, 0, 0))
        zi.compress_type = zipfile.ZIP_DEFLATED; zi.external_attr = 0o644 << 16
        z.writestr(zi, open(p, 'rb').read())
total = sum(os.path.getsize(p) for p in names); zsz, zsha = os.path.getsize(ZIP), sha(ZIP)
with open(HERE + '/BUILD_REPORT.md', 'w', encoding='utf-8', newline='\n') as f:
    w = f.write
    w('# BUILD_REPORT: kit_step3 v4\n\nBuilt by `kit_step3_v4_build/build.py` (%s) from the proof tree `<proof-tree>` = /home/claude/s1spec. Nothing certified was run by the build.\n\n' % time.strftime('%Y-%m-%d %H:%M'))
    w('**Zip `kit_step3_v4.zip`: %d bytes, sha256 `%s`.** Tree `kit_step3_v4/`: %d files, %d bytes.\n\n' % (zsz, zsha, len(names), total))
    w('Step-4 manifest the objects are pinned by: `step4/p2/out_v38/FINAL_manifest_v38.json` sha256 `%s`.\n\n' % sha(M38))
    w('## Source files asserted at build time\n\n| file | sha256 | pinned by v38 manifest |\n|---|---|---|\n')
    for p, h, i in REPORT['asserts']: w('| `%s` | `%s` | %s |\n' % (p, h, 'yes' if i else 'hashed only'))
    w('\n## Notes\n\n'); [w('- %s\n' % n) for n in REPORT['notes']]
    extra = open(HERE + '/todo_list.txt', encoding='utf-8').read().strip().split('\n') if os.path.exists(HERE + '/todo_list.txt') else []
    w('\n## [TODO] / limits (stated plainly)\n\n'); [w('- %s\n' % t) for t in extra + REPORT['todo']]
    w('\n## Patches (%d entries)\n\n| kit file | kind | old | new |\n|---|---|---|---|\n' % len(REPORT['patches']))
    for p in REPORT['patches']: w('| `%s` | %s | `%s` | `%s` |\n' % tuple(str(x).replace('|', '/').replace('`', "'") for x in p))
    w('\n## Data files containing /home/claude (provenance only)\n\n'); [w('- `%s`\n' % d) for d in data_hits]
    w('\n## Every file (path, bytes, sha256)\n\n| path | bytes | sha256 |\n|---|---|---|\n')
    for p in names: w('| `%s` | %d | `%s` |\n' % (os.path.relpath(p, K).replace(os.sep, '/'), os.path.getsize(p), sha(p)))
print('built', K, len(names), 'files', total, 'bytes; zip', zsz, zsha)
