# check_balance_v44.py -- (v44: copy of fix24_texts_v43/tools/check_balance_v43.py made by tools/copy_tools_v44.py; on code lines fix24_texts_v43 -> fix27_texts_v44 and v43 -> v44; comments below are the v43 history, unchanged.)
# check_balance_v43.py -- (v43: copy of fix22_texts_v42/tools/check_balance_v42.py; on code lines fix22_texts_v42 -> fix24_texts_v43 and v42 -> v43; comments below are the v42 history, unchanged; v43 additions marked 'v43'.)
# check_balance_v42.py -- (v42: copy of fix21_texts_v41/tools/check_balance_v41.py with v42 names; NAMEDIRS also include fix21_texts_v41/,
#  writeup/lean_bridge/ (+ snap_final/) and the Lean repository CWFlow (read-only, for the existence check (d) of new-text paths).)
# brace, \begin/\end, $ parity, \ref/\label checks for the text-v40 files (LaTeX is not installed).
# (v42: copy of fix19_texts_v40/tools/check_balance_v40.py with v42 names; NAMEDIRS also include fix19_texts_v40/ and fix20_delta_check/;
#  check (d) now covers the text written new for v42, tools/v42_inserts.json.)
# (v40: copy of fix17_texts_v39/tools/check_balance_v39.py with v40 names; the path roots also include fix17_texts_v39/ and fix18_final_review/.)
# Successor of fix16_texts_v38/tools/check_balance_v38.py (unchanged logic) plus: (a) a cross-document pass, in which every \ref of a
# standalone file must resolve to a \label somewhere in the set (unresolved anywhere = FAIL); (b) duplicate labels within a file and across
# the \input-expanded step3_v42_main.tex; (c) every \input target exists.  Writes tools/balance_v42.json.
import re, os, json
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/'
FILES = ['bootstrap_v44.tex', 'ei_lemma_v4c_v44.tex', 'E0_LEMMA_v4c_v44.tex', 'flux2_lemma_v44.tex', 'lemmaS_A_v44.tex', 'numbers_A_v44.tex',
         'continuation_v44.tex', 'lemmaU_v44.tex', 'main_v3_v44.tex', 'erratum_v44.tex', 'step3_v44_main.tex', 'b12_preamble_v44.tex',
         'chart_lemma_body_v44.tex', 'ei_constants_v44.tex', 'C1_identification_v44.tex', 'GS_lemma_v44.tex', 'U1_regularity_v44.tex',
         'G1_export_k2_v44.tex', 'step4_handoff_v44.tex']
def strip(t):
    out = []
    for line in t.split('\n'):
        line = re.sub(r'\\\\', '  ', line); line = re.sub(r'\\[{}%$]', '  ', line)
        i = line.find('%'); out.append(line if i < 0 else line[:i])
    return '\n'.join(out)
missing_inputs = []
def expand(name, seen=()):
    p = D + (name if name.endswith('.tex') else name + '.tex')
    if not os.path.exists(p): missing_inputs.append(name); return None
    t = open(p).read()
    def sub(m):
        e = expand(m.group(1), seen + (name,))
        return e if e is not None else m.group(0)
    return re.sub(r'^(?!%)\\input\{([^}]*)\}', lambda m: sub(m), t, flags=re.M)
def check(t):
    s = strip(t); r = {}
    d = 0; mn = 0
    for ch in s:
        d += (ch == '{') - (ch == '}'); mn = min(mn, d)
    r['brace_final'] = d; r['brace_min'] = mn
    st = []; bad = []
    for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', s):
        if m.group(1) == 'begin': st.append(m.group(2))
        elif not st or st.pop() != m.group(2): bad.append(m.group(2))
    r['env_unclosed'] = st; r['env_mismatch'] = bad
    r['dollar_even'] = (s.count('$') % 2 == 0)
    labs = re.findall(r'\\label\{([^}]*)\}', s)
    labels = set(labs); refs = set(re.findall(r'\\(?:ref|eqref|pageref)\{([^}]*)\}', s))
    r['refs_missing'] = sorted(refs - labels)
    r['duplicate_labels'] = sorted({l for l in labs if labs.count(l) > 1})
    r['ok'] = d == 0 and mn >= 0 and not st and not bad and r['dollar_even'] and not r['duplicate_labels']
    return r, labels, refs
res = {}; ALL = set(); REFS = {}
CTRL = {f: len(re.findall(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', open(D + f, encoding='utf-8', errors='replace').read())) for f in FILES if os.path.exists(D + f)}
for f in FILES:
    if not os.path.exists(D + f): res[f] = {'missing': True, 'ok': False}; continue
    r, l, rf = check(open(D + f).read()); res[f] = r; ALL |= l; REFS[f] = rf
    for m in re.finditer(r'^(?!%)\\input\{([^}]*)\}', open(D + f).read(), flags=re.M):
        if not os.path.exists(D + m.group(1) + ('' if m.group(1).endswith('.tex') else '.tex')): missing_inputs.append(f + ' -> ' + m.group(1))
t = expand('step3_v44_main.tex'); r, l, rf = check(t); r['ok'] = r['ok'] and not r['refs_missing']
res['step3_v44_main.tex [expanded]'] = r
cross = {f: sorted(rf - ALL) for f, rf in REFS.items() if rf - ALL}
res['control_chars'] = {'per_file': {f: n for f, n in CTRL.items() if n}, 'ok': not any(CTRL.values())}
res['cross_document'] = {'labels_in_set': len(ALL), 'unresolved_anywhere': cross, 'missing_inputs': sorted(set(missing_inputs)),
                         'ok': not cross and not missing_inputs}
json.dump(res, open(D + 'tools/balance_v44.json', 'w'), indent=1)
for k, v in res.items():
    print(('PASS' if v.get('ok') else 'FAIL'), k, {a: b for a, b in v.items() if a != 'ok' and b not in ([], 0, True, {})})
# (d) every \fn{...}/\pth{...} path with a '/' or a known file suffix in the text written NEW for v39 (tools/v42_inserts.json) must exist
#     under one of the roots below (backticked paths of the .md inserts as well).  Pure names without a directory are looked up too.
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
ROOTS = ['/home/claude/', '/home/claude/s1spec/', '/home/claude/s1spec/step4/', '/home/claude/s1spec/step4/p2/', '/home/claude/s1spec/step3nl_v4c/',
         '/home/claude/s1spec/writeup/', '/home/claude/s1spec/writeup/lemmaU/bulk/', D, '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/']
NAMEDIRS = ['/home/claude/s1spec/', '/home/claude/s1spec/step4/p2/out_v38/', '/home/claude/s1spec/step4/p2/out_v38/iface/', '/home/claude/s1spec/step4/p2/v38/',
            '/home/claude/s1spec/step4/p2/v33/', D, '/home/claude/s1spec/writeup/step4_proof/', '/home/claude/s1spec/step3nl_v4c/v4bcert/out/',
            '/home/claude/s1spec/step3nl_v4c/runs/', '/home/claude/s1spec/step4/p2/out_v33/', B + 'fix1/lemmaU_rev/', B + 'fix10_proofs/', B + 'fix16_texts_v38/',
            B + 'fix13_v37chain/flux/out/', B + 'fix13_v37chain/n1/out/', B + 'fix13_v37chain/t5_zone3/out/',
            B + 'fix17_texts_v39/', B + 'fix19_texts_v40/', B + 'fix21_texts_v41/', B + 'fix18_final_review/', B + 'fix20_delta_check/', B + 'fix22_texts_v42/', B + 'fix23_delta_check/', B + 'fix24_texts_v43/', B + 'fix25_check_v43/', B + 'fix26_latex_v43/', '/home/claude/s1spec/writeup/lean_bridge/snap_s3/',  # v43: + fix22, fix23, snap_s3; v44: + fix24, fix25, fix26
           
            '/home/claude/s1spec/writeup/lean_bridge/', '/home/claude/s1spec/writeup/lean_bridge/snap_final/', '/mnt/c/Users/micahw/',
            '/mnt/c/Users/micahw/agents/lean/CWFlow/', '/mnt/c/Users/micahw/agents/lean/CWFlow/CWFlow/', '/mnt/c/Users/micahw/agents/lean/CWFlow/CWFlow/Curve/', B + 'fix13_v37chain/t2_A/out/', '/home/claude/s1spec/writeup/lit_P1/', B + 'fix1/b12_text/', '/home/claude/s1spec/step3nl_v4c/verifier/', '/home/claude/s1spec/step4/p2/v33/']
# names that a new sentence quotes ON PURPOSE although they do not exist (the wrong name that text v38 printed, erratum O4)
ALLOWED_ABSENT = {'runs/reg_fin13v4c': 'erratum O4 quotes the wrong directory name printed by text v38',
                  'leanprover/lean4:v4.33.1': 'v44: Lean toolchain name (content of CWFlow/lean-toolchain), not a path'}
ins = json.load(open(D + 'tools/v44_inserts.json'))
paths = set()
for it in ins:
    for m in re.finditer(r'\\(?:fn|pth)\{([^}]*)\}|`([^`\s]+)`', it['text']):
        p = (m.group(1) or m.group(2)).replace('\\_', '_')
        if re.fullmatch(r'[0-9a-f]{8,64}(\.\.\.)?', p) or ' ' in p or '$' in p or p.startswith('--'): continue
        if '/' in p or re.search(r'\.(tex|json|md|py|zip|pkl)$', p): paths.add(p.rstrip('/').rstrip('.'))
def exists(p):
    if p in ALLOWED_ABSENT: return True
    if p.endswith('/*'): p = p[:-2]
    if p.startswith('/'): return os.path.exists(p)
    return any(os.path.exists(r + p) for r in ROOTS + NAMEDIRS)
bad = sorted(p for p in paths if not exists(p))
res['new_text_paths'] = {'checked': len(paths), 'not_found': bad, 'ok': not bad}
json.dump(res, open(D + 'tools/balance_v44.json', 'w'), indent=1)
print(('PASS' if not bad else 'FAIL'), 'new_text_paths', len(paths), 'checked; not found:', bad)
