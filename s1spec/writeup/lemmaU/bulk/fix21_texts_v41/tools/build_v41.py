# build_v41.py -- builds the step-3 text set v41 in fix21_texts_v41/ from the v40 set in fix19_texts_v40/ (read-only).
# (v41: follows fix19_texts_v40/tools/build_v40.py.)  Work order: the wording items W-1..W-5 and the cosmetic retag item of the delta
# check fix20_delta_check/DELTA_CHECK_v40.md, section 4:
#   W-1 order inside the closure proof (bootstrap: Closed step and Open Step 1: GS, G -> U-1 -> the rest; Q continuous where GS, G hold);
#   W-2 Theorem U-1 (c), (d): the inequalities hold at every tau wherever the hypotheses of their lemmas hold;
#   W-3 the alpha-barrier on a march step (bootstrap Step 2; flux lemma Prop. F:prop:abs);
#   W-4 flux lemma Lemma F:lem:data: |alpha| <= 2.5e-4 (box (B-a)) instead of alpha_1;
#   W-5 lemmaU implication step, main_v3 bootstrap paragraph: B1(a,b) excepted from "strictly better";
#   cosmetic: lemmaS_A title / date retag.
# The replacements are in tools/edits_v41.txt (plain text, each asserted to occur the stated number of times).  Text only: no number is
# computed here.  Erratum Part 6 (section Q) from tools/erratum_part6_v41.tex; records addition from tools/records_errata_add_v41.md.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v41.py   (after gen_records_v41.py, gen_numbers_A_v41.py, gen_ei_constants_v41.py)
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b41_common as C
from b41_common import Builder, D, D9
R = json.load(open(D + 'tools/records_v41.json')); P = R['pins']
assert (P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']) == (992, 646, 629, 55, 38), P
PH = {'@DC@': r'\fn{fix20_delta_check/DELTA_CHECK_v40.md}'}
def sub(s):
    for k, v in PH.items(): s = s.replace(k, v)
    assert not re.search(r'@[A-Z0-9]+@', s), s[:200]
    return s
# ---- parse edits_v41.txt
EDITS = {}; ORDER = []
cur = None; mode = None; buf = []; old = None; kind = None
for line in open(D + 'tools/edits_v41.txt').read().split('\n'):
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
HDR = ('%% %s -- text v41 (fix21_texts_v41, 21 Nov 2026), built by fix21_texts_v41/tools/build_v41.py from fix19_texts_v40/%s (unchanged):\n'
       '%% wording items W-1..W-5 and the retag item of the delta check fix20_delta_check/DELTA_CHECK_v40.md where they touch this file.  Not compiled.')
TEX = ['bootstrap_v40.tex', 'flux2_lemma_v40.tex', 'lemmaU_v40.tex', 'U1_regularity_v40.tex', 'continuation_v40.tex', 'GS_lemma_v40.tex',
       'step3_v40_main.tex', 'main_v3_v40.tex', 'step4_handoff_v40.tex', 'lemmaS_A_v40.tex', 'G1_export_k2_v40.tex',
       'C1_identification_v40.tex', 'E0_LEMMA_v4c_v40.tex', 'b12_preamble_v40.tex', 'chart_lemma_body_v40.tex', 'ei_lemma_v4c_v40.tex']
assert set(ORDER) <= set(TEX), set(ORDER) - set(TEX)
out = []; nedits = 0; per = {}
for name in TEX:
    b = Builder(name, header=HDR % (C.v41name(name), name))
    for (tag, k), o, n in EDITS.get(name, []):
        n = sub(n)
        if k is None: b.rep(o, n, new=(tag == 'REP'))
        else: b.rep_all(o, n, k, new=(tag == 'REP'))
        nedits += 1
    per[name] = b.n
    out.append(b.write())
# ---- the step-4 patch note (markdown): retag, plus one v41 line (no edit added)
b = Builder('STEP4_PATCH_NOTE_v40.md')
add = '\n\n## Text v41 (21 Nov 2026; delta check fix20)\n- No edit is added or changed (the text-v41 items are wording in step-3 files only).\n'
b.t = b.t.rstrip('\n') + add
C.INSERTS.append({'file': b.name, 'text': add})
out.append(b.write())
# ---- erratum: Parts 1-5 verbatim (history: no retag of the body), Part 6 new
p6 = open(D + 'tools/erratum_part6_v41.tex').read()
b = Builder('erratum_v40.tex', retag_names=False, header='% erratum_v41.tex -- text v41: erratum_v40.tex (fix19_texts_v40, unchanged) verbatim as Parts 1-5, plus Part 6 (section Q, text v41).  Built by fix21_texts_v41/tools/build_v41.py.  Not compiled.')
b.rep('\\end{document}', p6.rstrip('\n') + '\n\\end{document}')
out.append(b.write())
md = open(D9 + 'ERRATUM_v40.md').read()
open(D + 'ERRATUM_v41.md', 'w').write('# ERRATUM_v41 (text v41): ERRATUM_v40.md verbatim, then a pointer to Part 6\n\nPart 6 (text v41, delta check fix20) is section Q of `erratum_v41.tex`; its source is `tools/erratum_part6_v41.tex`.\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v41.md')
rec = open(D9 + 'RECORDS_ERRATA_v40.md').read(); rec41 = open(D + 'tools/records_errata_add_v41.md').read()
open(D + 'RECORDS_ERRATA_v41.md', 'w').write('# RECORDS_ERRATA_v41: RECORDS_ERRATA_v40.md verbatim, then the v41 section at the end\n\n' + rec + '\n\n' + rec41)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v41.md', 'text': rec41})
out.append(D + 'RECORDS_ERRATA_v41.md')
C.save_inserts()
# ---- sanity: no 0.1475 quote outside the history files; no v40 set-file name left outside them
for f in out:
    t = re.sub(r'fix19_texts_v40/[\w.\\-]+', ' ', open(f).read()); bn = os.path.basename(f)   # provenance paths into fix19 are intended
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
    left = sorted(set(m.group(0) for m in re.finditer(r'(?<![A-Za-z0-9])(?:%s)(?:_|\\_)v40\.(?:tex|md|json)|step3(?:_|\\_)v40(?:_|\\_)main' % '|'.join(map(re.escape, C.STEMS)), t)))
    assert not left, (bn, left)
# ---- the files touched by edits are exactly those of the delta-check items
assert {k for k, v in per.items() if v} == {'bootstrap_v40.tex', 'flux2_lemma_v40.tex', 'U1_regularity_v40.tex', 'lemmaU_v40.tex', 'main_v3_v40.tex', 'lemmaS_A_v40.tex'}, per
print('built %d files, %d edits' % (len(out), nedits))
for k in TEX: print('  edits in %s: %d' % (k, per[k]))
