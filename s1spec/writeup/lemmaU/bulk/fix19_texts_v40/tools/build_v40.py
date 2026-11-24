# build_v40.py -- builds the step-3 text set v40 in fix19_texts_v40/ from the v39 set in fix17_texts_v39/ (read-only).
# Work order: the findings of the final review fix18_final_review/FINAL_REVIEW_v39.md:
#   F-1 closure proof (B1(a,b) have no gap: re-derived by first touch in the Open step; U-1 regularity without B1(b));
#   F-2 macros (numbers file by gen_numbers_A_v40.py; \input of numbers_A_v40 in lemmaU, main_v3 and the flux lemma standalone preamble);
#   F-3 symmetry of gamma(t_1) in (H6); F-4 one meaning of T_max (a length); F-5 DKS Theorem 3.1 cited by its statement, (D1) a remark,
#   author note on "maximal solution"; F-6 n^III level 0; F-7 (H4) range certified by stage 12; R-1 pin counts; R-2, R-3 records.
# The replacements are in tools/edits_v40.txt (plain text, each asserted to occur the stated number of times).  Text only: no number is
# computed here; record values come from tools/records_v40.json (gen_records_v40.py, read from disk).
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v40.py   (after gen_records_v40.py, gen_numbers_A_v40.py, gen_ei_constants_v40.py)
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b40_common as C
from b40_common import Builder, D, D9
R = json.load(open(D + 'tools/records_v40.json')); P = R['pins']
assert (P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']) == (992, 646, 629, 55, 38), P
PH = {'@PINS@': '%d pin entries on %d distinct files: %d stage paths and %d scripts, %d of them also stage paths' % (
          P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']),
      '@REV@': r'\fn{fix18_final_review/FINAL_REVIEW_v39.md}'}
PH['@FINREV@'] = (r'an independent check of the text-v40 changes (the final review of text v39, ' + PH['@REV@'] + r', is done: ``closes with fixes, '
                  r"none blocking'', modulo the external and author items; its fixes F-1 to F-7, R-1 to R-3 are applied in text v40)")
PH['@H4CERT@'] = (r'stage 12 of \fn{step4/p2/out_v38/FINAL_manifest_v38.json} (M9 arb base sums, \fn{step4/p2/v33/win_arb.py}, item (e) of its header) '
                  r'asserts, for every $t_c$ ball and both ball sets, that the cells used are consecutive along the arc and that $\cos\varphi>0$ on every '
                  r'cell of the ranges $y\in[14,60]$, $[50,380.25]$ and $[60,700]$, i.e.\ $y$ is strictly increasing along the arc there; so the arc is '
                  r'a graph over the $S_1$ ray on $y\in[14,700]\supset[48,700]$ \tC')
def sub(s):
    for k, v in PH.items(): s = s.replace(k, v)
    assert '@' not in re.sub(r'@@', '', s) or not re.search(r'@[A-Z0-9]+@', s), s[:200]
    return s
# ---- parse edits_v40.txt
EDITS = {}; ORDER = []
cur = None; mode = None; buf = []; old = None; kind = None
for line in open(D + 'tools/edits_v40.txt').read().split('\n'):
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
HDR = ('%% %s -- text v40 (fix19_texts_v40, 21 Nov 2026), built by fix19_texts_v40/tools/build_v40.py from fix17_texts_v39/%s (unchanged):\n'
       '%% fixes of the final review fix18_final_review/FINAL_REVIEW_v39.md (F-1..F-7, R-1..R-3) where they touch this file.  Not compiled.')
TEX = ['bootstrap_v39.tex', 'flux2_lemma_v39.tex', 'lemmaU_v39.tex', 'U1_regularity_v39.tex', 'continuation_v39.tex', 'GS_lemma_v39.tex',
       'step3_v39_main.tex', 'main_v3_v39.tex', 'step4_handoff_v39.tex', 'lemmaS_A_v39.tex', 'G1_export_k2_v39.tex',
       'C1_identification_v39.tex', 'E0_LEMMA_v4c_v39.tex', 'b12_preamble_v39.tex', 'chart_lemma_body_v39.tex', 'ei_lemma_v4c_v39.tex']
assert set(ORDER) <= set(TEX), set(ORDER) - set(TEX)
out = []; nedits = 0
for name in TEX:
    b = Builder(name, header=HDR % (C.v40name(name), name))
    for (tag, k), o, n in EDITS.get(name, []):
        n = sub(n)
        if k is None: b.rep(o, n, new=(tag == 'REP'))
        else: b.rep_all(o, n, k, new=(tag == 'REP'))
        nedits += 1
    out.append(b.write())
# ---- the step-4 patch note (markdown): retag, plus one v40 section (no edit added)
b = Builder('STEP4_PATCH_NOTE_v39.md')
add = ('\n\n## Text v40 (21 Nov 2026; review fix18)\n- No edit is added or changed. Two notes for a future step-4 version: (F-7) the (H4) range '
       '$y\\in[14,700]$ is certified by the coverage assertion of stage 12 of `FINAL_manifest_v38.json` (`p2/v33/win_arb.py` item (e)), so the '
       'hand-off list can state the extended range as certified; (F-3) the hand-off list should state the symmetry of $\\gamma(t_1)$ about the '
       '$x$-axis next to (H6), as `step4_handoff_v40.tex` item (3) does.\n')
b.t = b.t.rstrip('\n') + add
C.INSERTS.append({'file': b.name, 'text': add})
out.append(b.write())
# ---- erratum: Parts 1-4 verbatim (history: no retag of the body), Part 5 new
p5 = open(D + 'tools/erratum_part5_v40.tex').read()
b = Builder('erratum_v39.tex', retag_names=False, header='% erratum_v40.tex -- text v40: erratum_v39.tex (fix17_texts_v39, unchanged) verbatim as Parts 1-4, plus Part 5 (section P, text v40).  Built by fix19_texts_v40/tools/build_v40.py.  Not compiled.')
b.rep('\\end{document}', p5.rstrip('\n') + '\n\\end{document}')
out.append(b.write())
md = open(D9 + 'ERRATUM_v39.md').read()
open(D + 'ERRATUM_v40.md', 'w').write('# ERRATUM_v40 (text v40): ERRATUM_v39.md verbatim, then a pointer to Part 5\n\nPart 5 (text v40, review fix18) is section P of `erratum_v40.tex`; its source is `tools/erratum_part5_v40.tex`.\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v40.md')
rec = open(D9 + 'RECORDS_ERRATA_v39.md').read(); rec40 = open(D + 'tools/records_errata_add_v40.md').read()
open(D + 'RECORDS_ERRATA_v40.md', 'w').write('# RECORDS_ERRATA_v40: RECORDS_ERRATA_v39.md verbatim, then the v40 section at the end\n\n' + rec + '\n\n' + rec40)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v40.md', 'text': rec40})
out.append(D + 'RECORDS_ERRATA_v40.md')
C.save_inserts()
# ---- sanity: no 0.1475 quote outside the history files; no v39 set-file name left outside them
for f in out:
    t = re.sub(r'fix17_texts_v39/[\w.\\-]+', ' ', open(f).read()); bn = os.path.basename(f)   # provenance paths into fix17 are intended
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
    left = sorted(set(m.group(0) for m in re.finditer(r'(?<![A-Za-z0-9])(?:%s)(?:_|\\_)v39\.(?:tex|md|json)|step3(?:_|\\_)v39(?:_|\\_)main' % '|'.join(map(re.escape, C.STEMS)), t)))
    assert not left, (bn, left)
print('built %d files, %d edits' % (len(out), nedits))
