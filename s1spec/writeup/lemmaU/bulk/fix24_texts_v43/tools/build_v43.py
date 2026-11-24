# build_v43.py -- builds the step-3 text set v43 in fix24_texts_v43/ from the v42 set in fix22_texts_v42/ (read-only).
# (v43: follows fix22_texts_v42/tools/build_v42.py.)  Work: the findings of the delta check
# writeup/lemmaU/bulk/fix23_delta_check/DELTA_CHECK_v42.md, should-fix S-1..S-4 and wording W-1..W-9.
# The replacements are in tools/edits_v43.txt (each asserted to occur exactly once).  Text only: no number is computed here.
# S-1: the compile records of the S3 files are writeup/lean_bridge/raw_{Periodic2Circle,MaxShift,All}.log if all three exist on disk
#      at build time; otherwise the text points to IMPL_LOG.md "S3 links" and says the per-file logs are pending.
# Erratum: Parts 1-6 verbatim, Part 7 (section R) with four marked corrections (edits), Part 8 (section S) from tools/erratum_part8_v43.tex;
# records addition from tools/records_errata_add_v43.md.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v43.py   (after gen_records_v43.py, gen_numbers_A_v43.py, gen_ei_constants_v43.py)
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b43_common as C
from b43_common import Builder, D, D9
R = json.load(open(D + 'tools/records_v43.json')); P = R['pins']
assert (P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']) == (992, 646, 629, 55, 38), P
# ---- S-1: compile records of the S3 files
LB = '/home/claude/s1spec/writeup/lean_bridge/'
LOGS = ['raw_Periodic2Circle.log', 'raw_MaxShift.log', 'raw_All.log']
HAVE_LOGS = all(os.path.isfile(LB + f) for f in LOGS)
if HAVE_LOGS:
    S3TEX = ('text v43 (delta check fix23 S-1): the compile records of the S3 files \\fn{Periodic2Circle}, \\fn{MaxShift} and \\fn{All} are the raw logs '
             '\\fn{writeup/lean_bridge/raw_Periodic2Circle.log}, \\fn{writeup/lean_bridge/raw_MaxShift.log} and \\fn{writeup/lean_bridge/raw_All.log}; '
             'text v42 had no raw log for them')
    S3ERR = 'the raw logs \\fn{writeup/lean_bridge/raw_Periodic2Circle.log}, \\fn{raw_MaxShift.log}, \\fn{raw_All.log}'
    S3MD = '`writeup/lean_bridge/raw_Periodic2Circle.log`, `raw_MaxShift.log`, `raw_All.log` (cited in U-1 section 1.2).'
else:
    S3TEX = ('text v43 (delta check fix23 S-1): for the S3 files \\fn{Periodic2Circle}, \\fn{MaxShift} and \\fn{All} no raw compile log was kept; '
             'compile logs: see \\fn{writeup/lean_bridge/IMPL_LOG.md} ``S3 links\'\'; per-file logs pending')
    S3ERR = 'no raw log was kept; compile logs: see \\fn{writeup/lean_bridge/IMPL_LOG.md} ``S3 links\'\'; per-file logs pending'
    S3MD = 'no raw log kept; see `writeup/lean_bridge/IMPL_LOG.md` "S3 links"; per-file logs pending (expected next to the other raw logs in writeup/lean_bridge/, named raw_<file>.log).'
print('S-1 compile logs on disk:', HAVE_LOGS)
FILL = {'@S3LOGS@': S3TEX, '@S3LOGS_ERR@': S3ERR, '@S3LOGS_MD@': S3MD}
def sub(s):
    for k, v in FILL.items(): s = s.replace(k, v)
    assert not re.search(r'@[A-Z0-9_]+@', s), s[:200]
    return s
# ---- parse edits_v43.txt
EDITS = {}; ORDER = []
cur = None; mode = None; buf = []; old = None; kind = None
for line in open(D + 'tools/edits_v43.txt').read().split('\n'):
    if line.startswith('## ') and mode is None: continue
    m = re.match(r'@@ (FILE|REP|REPQ|NEW|END)(?: (.*))?$', line)
    if m:
        tag, arg = m.group(1), m.group(2)
        if tag == 'FILE': cur = arg; EDITS[cur] = []; ORDER.append(cur); mode = None
        elif tag in ('REP', 'REPQ'): kind = (tag, int(arg) if arg else None); mode = 'old'; buf = []
        elif tag == 'NEW': old = '\n'.join(buf); mode = 'new'; buf = []
        elif tag == 'END': EDITS[cur].append((kind, old, '\n'.join(buf))); mode = None
        continue
    if mode: buf.append(line)
HDR = ('%% %s -- text v43 (fix24_texts_v43, 23 Nov 2026), built by fix24_texts_v43/tools/build_v43.py from fix22_texts_v42/%s (unchanged):\n'
       '%% the findings of the delta check fix23_delta_check/DELTA_CHECK_v42.md (S-1..S-4, W-1..W-9), where they touch this file.  Not compiled.')
TEX = ['bootstrap_v42.tex', 'flux2_lemma_v42.tex', 'lemmaU_v42.tex', 'U1_regularity_v42.tex', 'continuation_v42.tex', 'GS_lemma_v42.tex',
       'step3_v42_main.tex', 'main_v3_v42.tex', 'step4_handoff_v42.tex', 'lemmaS_A_v42.tex', 'G1_export_k2_v42.tex',
       'C1_identification_v42.tex', 'E0_LEMMA_v4c_v42.tex', 'b12_preamble_v42.tex', 'chart_lemma_body_v42.tex', 'ei_lemma_v4c_v42.tex']
assert set(ORDER) <= set(TEX) | {'erratum_v42.tex'}, set(ORDER) - set(TEX)
out = []; nedits = 0; per = {}
def apply(b, name):
    global nedits
    for (tag, k), o, n in EDITS.get(name, []):
        n = sub(n)
        if k is None: b.rep(o, n, new=(tag == 'REP'))
        else: b.rep_all(o, n, k, new=(tag == 'REP'))
        nedits += 1
    per[name] = b.n
for name in TEX:
    b = Builder(name, header=HDR % (C.v43name(name), name))
    apply(b, name)
    out.append(b.write())
# ---- the step-4 patch note (markdown): retag, plus a v43 section
b = Builder('STEP4_PATCH_NOTE_v42.md')
add = ('\n\n## Text v43 (23 Nov 2026; delta check fix23_delta_check/DELTA_CHECK_v42.md)\n'
       '- No edit is added or changed. The proposed edits still need the authors\' go before a new step-4 version is written.\n')
b.t = b.t.rstrip('\n') + add
C.INSERTS.append({'file': b.name, 'text': add})
out.append(b.write())
# ---- erratum: Parts 1-6 verbatim; Part 7 with four marked corrections; Part 8 new
p8 = sub(open(D + 'tools/erratum_part8_v43.tex').read())
b = Builder('erratum_v42.tex', retag_names=False, header='% erratum_v43.tex -- text v43: erratum_v42.tex (fix22_texts_v42, unchanged) with Parts 1-6 verbatim, Part 7 (section R) with four corrections marked [text v43], plus Part 8 (section S, text v43).  Built by fix24_texts_v43/tools/build_v43.py.  Not compiled.')
apply(b, 'erratum_v42.tex')
b.rep('\\end{document}', p8.rstrip('\n') + '\n\\end{document}')
out.append(b.write())
md = open(D9 + 'ERRATUM_v42.md').read()
open(D + 'ERRATUM_v43.md', 'w').write('# ERRATUM_v43 (text v43): ERRATUM_v42.md verbatim, then a pointer to Part 8\n\nPart 8 (text v43: the delta check fix23_delta_check/DELTA_CHECK_v42.md of text v42) is section S of `erratum_v43.tex`; its source is `tools/erratum_part8_v43.tex`. Part 7 (section R) of `erratum_v43.tex` carries four corrections marked [text v43] (R2, R3(b), R4, R6).\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v43.md')
rec = open(D9 + 'RECORDS_ERRATA_v42.md').read(); rec43 = sub(open(D + 'tools/records_errata_add_v43.md').read())
open(D + 'RECORDS_ERRATA_v43.md', 'w').write('# RECORDS_ERRATA_v43: RECORDS_ERRATA_v42.md verbatim, then the v43 section at the end\n\n' + rec + '\n\n' + rec43)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v43.md', 'text': rec43})
out.append(D + 'RECORDS_ERRATA_v43.md')
C.save_inserts()
# ---- sanity: no 0.1475 quote outside the history files; no v42 set-file name left outside them
for f in out:
    t = re.sub(r'fix22_texts_v42/[\w.\\-]+', ' ', open(f).read()); bn = os.path.basename(f)   # provenance paths into fix22 are intended
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
    left = sorted(set(m.group(0) for m in re.finditer(r'(?<![A-Za-z0-9])(?:%s)(?:_|\\_)v42\.(?:tex|md|json)|step3(?:_|\\_)v42(?:_|\\_)main' % '|'.join(map(re.escape, C.STEMS)), t)))
    assert not left, (bn, left)
# ---- v43 content checks (delta check fix23)
NONHIST = [f for f in out if os.path.basename(f).endswith('.tex') and not os.path.basename(f).startswith('erratum')]
ALLT = {os.path.basename(f): open(f).read() for f in out}
for f in NONHIST:
    t = ALLT[os.path.basename(f)]; bn = os.path.basename(f)
    if '\\LeanSthree' in t.replace('\\providecommand{\\LeanSthree}', ''):
        if bn == 'step3_v43_main.tex': assert '\\providecommand{\\LeanSthree}' in ALLT['b12_preamble_v43.tex'] and '\\input{b12_preamble_v43}' in t, bn
        else: assert t.count('\\providecommand{\\LeanSthree}') == 1, bn
    # S-2: nobody says All.lean was never compiled, except as a quote of what text v42 said
    for line in t.split('\n'):
        if 'never compiled' in line and 'All.lean' in line: assert 'text v42 said' in line, (bn, line[:160])
    # S-4: no limit form of the L^2 criterion
    assert not re.search(r'\\int\\kappa\^2\s*ds\s*\\to\s*\\infty', t.replace("``$\\int\\kappa^2ds\\to\\infty$''", '')), bn
U1 = ALLT['U1_regularity_v43.tex']
assert '\\tD' not in U1                                            # W-9
assert 'T_{\\max}-t_1)$ \\LeanSthree\\ (review W7 S3b)' in U1     # W-3
assert ('raw_MaxShift.log' in U1) == HAVE_LOGS                     # S-1
# W-7: every file that names snap_final/SHA256SUMS also names snap_s3/SHA256SUMS
for bn, t in ALLT.items():
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')) or bn.endswith('.md'): continue
    if 'snap_final/SHA256SUMS' in t: assert 'snap_s3/SHA256SUMS' in t, bn
# W-1: every "review W7 ... sound" sentence that sits next to (PC)/(MS)/(GE) or the S3 names in a status line also cites fix23
for bn in ('step3_v43_main.tex', 'main_v3_v43.tex', 'continuation_v43.tex', 'C1_identification_v43.tex', 'U1_regularity_v43.tex'):
    assert 'fix23_delta_check/DELTA_CHECK_v42.md' in ALLT[bn], bn
assert 'its citation; an independent check of the S3 theorems' not in ALLT['erratum_v43.tex']   # W-1 (R6)
assert 'These theorems have not been reviewed independently' not in ALLT['erratum_v43.tex']
assert 'text v42 states it as \\fn{rst_glue}' not in ALLT['erratum_v43.tex']   # S-3
assert set(k for k, v in per.items() if v) == {'U1_regularity_v42.tex', 'C1_identification_v42.tex', 'continuation_v42.tex', 'step3_v42_main.tex',
        'main_v3_v42.tex', 'lemmaU_v42.tex', 'bootstrap_v42.tex', 'b12_preamble_v42.tex', 'erratum_v42.tex'}, per
print('built %d files, %d edits' % (len(out), nedits))
for k in TEX + ['erratum_v42.tex']: print('  edits in %s: %d' % (k, per[k]))
