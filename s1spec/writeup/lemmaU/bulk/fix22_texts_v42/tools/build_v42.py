# build_v42.py -- builds the step-3 text set v42 in fix22_texts_v42/ from the v41 set in fix21_texts_v41/ (read-only).
# (v42: follows fix21_texts_v41/tools/build_v41.py.)  Work order:
#   (1) Lean replaces DKS: the twelve text changes of writeup/lean_bridge/REVIEW_W7.md section 3 (U1 sections 1.2-2 and gaps; C-1 labels,
#       (I4), (I5), Lemma C1:lem:max, Theorem (v), gaps; continuation, bootstrap, step3_main, ei_lemma, GS, lemmaU, main_v3, step4_handoff,
#       b12_preamble; erratum Part 7; records section);
#   (2) the bridges S3a, S3b of review W7, formalized after it (lean_bridge/snap_s3), cited through ONE macro per file, \LeanSthree (prints [Lean]);
#   (3) the v41 leftovers (step3_main title via b42_common.TITLE_REPS; step3_main open list);
#   (4) the open-item lists (no Lean build / n=1 / DKS items; full lake build and commit added; S3a/S3b proved in Lean).
# The replacements are in tools/edits_v42.txt (plain text, each asserted to occur exactly once).  Text only: no number is computed here.
# Erratum Part 7 (section R) from tools/erratum_part7_v42.tex; records addition from tools/records_errata_add_v42.md.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v42.py   (after gen_records_v42.py, gen_numbers_A_v42.py, gen_ei_constants_v42.py)
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b42_common as C
from b42_common import Builder, D, D9
R = json.load(open(D + 'tools/records_v42.json')); P = R['pins']
assert (P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']) == (992, 646, 629, 55, 38), P
def sub(s):
    assert not re.search(r'@[A-Z0-9]+@', s), s[:200]
    return s
# ---- parse edits_v42.txt
EDITS = {}; ORDER = []
cur = None; mode = None; buf = []; old = None; kind = None
for line in open(D + 'tools/edits_v42.txt').read().split('\n'):
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
HDR = ('%% %s -- text v42 (fix22_texts_v42, 23 Nov 2026), built by fix22_texts_v42/tools/build_v42.py from fix21_texts_v41/%s (unchanged):\n'
       '%% Lean replaces DKS (writeup/lean_bridge/REVIEW_W7.md section 3) and the v41 leftovers, where they touch this file.  Not compiled.')
TEX = ['bootstrap_v41.tex', 'flux2_lemma_v41.tex', 'lemmaU_v41.tex', 'U1_regularity_v41.tex', 'continuation_v41.tex', 'GS_lemma_v41.tex',
       'step3_v41_main.tex', 'main_v3_v41.tex', 'step4_handoff_v41.tex', 'lemmaS_A_v41.tex', 'G1_export_k2_v41.tex',
       'C1_identification_v41.tex', 'E0_LEMMA_v4c_v41.tex', 'b12_preamble_v41.tex', 'chart_lemma_body_v41.tex', 'ei_lemma_v4c_v41.tex']
assert set(ORDER) <= set(TEX), set(ORDER) - set(TEX)
out = []; nedits = 0; per = {}
for name in TEX:
    b = Builder(name, header=HDR % (C.v42name(name), name))
    for (tag, k), o, n in EDITS.get(name, []):
        n = sub(n)
        if k is None: b.rep(o, n, new=(tag == 'REP'))
        else: b.rep_all(o, n, k, new=(tag == 'REP'))
        nedits += 1
    per[name] = b.n
    out.append(b.write())
# ---- the step-4 patch note (markdown): retag; the two sentences that describe (P1') / C-1 inputs are updated (they are proposals, not
#      applied to the closed step4_proof.tex), plus a v42 section
b = Builder('STEP4_PATCH_NOTE_v41.md')
b.rep("which is DKS 2002 Thm 3.1 + U-1 + C-1(v) and keeps its open items (Lean build, DKS re-read).",
      "which since text v42 is the compiled Lean theorem `cdf_exists_maximal` + U-1 + C-1(v) and keeps its open items (full `lake build CWFlow.Curve.All`, commit of `CWFlow/Curve/`, DKS is no longer an input; the bridges S3a, S3b are Lean theorems).")
b.rep("(Lean (L1), (L2) read not built; DKS (D1); Lemma (P2-a);",
      "(compiled Lean theorems `cdf_exists_maximal`, `cdf_symmetric` of `CWFlow/Curve/`, commit pending; data bridge by `exists_isLagImmersion_of_periodic`; Lemma (P2-a);")
add = ('\n\n## Text v42 (23 Nov 2026; Lean replaces DKS, review writeup/lean_bridge/REVIEW_W7.md)\n'
       '- No edit is added. Two proposed sentences were updated: the note on l.789/l.791 ((P1\') is now the compiled Lean theorem `cdf_exists_maximal` + U-1 + C-1(v)) '
       'and the proposed P1/continuation sentence (C-1 rests on compiled Lean theorems, commit pending, instead of Lean read from source and DKS (D1)). '
       'These edits still need the authors\' go before a new step-4 version is written.\n')
b.t = b.t.rstrip('\n') + add
C.INSERTS.append({'file': b.name, 'text': add})
out.append(b.write())
# ---- erratum: Parts 1-6 verbatim (history: no retag of the body), Part 7 new
p7 = open(D + 'tools/erratum_part7_v42.tex').read()
b = Builder('erratum_v41.tex', retag_names=False, header='% erratum_v42.tex -- text v42: erratum_v41.tex (fix21_texts_v41, unchanged) verbatim as Parts 1-6, plus Part 7 (section R, text v42).  Built by fix22_texts_v42/tools/build_v42.py.  Not compiled.')
b.rep('\\end{document}', p7.rstrip('\n') + '\n\\end{document}')
out.append(b.write())
md = open(D9 + 'ERRATUM_v41.md').read()
open(D + 'ERRATUM_v42.md', 'w').write('# ERRATUM_v42 (text v42): ERRATUM_v41.md verbatim, then a pointer to Part 7\n\nPart 7 (text v42: the Lean theorems replace DKS, review lean_bridge/REVIEW_W7.md) is section R of `erratum_v42.tex`; its source is `tools/erratum_part7_v42.tex`.\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v42.md')
rec = open(D9 + 'RECORDS_ERRATA_v41.md').read(); rec42 = open(D + 'tools/records_errata_add_v42.md').read()
open(D + 'RECORDS_ERRATA_v42.md', 'w').write('# RECORDS_ERRATA_v42: RECORDS_ERRATA_v41.md verbatim, then the v42 section at the end\n\n' + rec + '\n\n' + rec42)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v42.md', 'text': rec42})
out.append(D + 'RECORDS_ERRATA_v42.md')
C.save_inserts()
# ---- sanity: no 0.1475 quote outside the history files; no v41 set-file name left outside them
for f in out:
    t = re.sub(r'fix21_texts_v41/[\w.\\-]+', ' ', open(f).read()); bn = os.path.basename(f)   # provenance paths into fix21 are intended
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
    left = sorted(set(m.group(0) for m in re.finditer(r'(?<![A-Za-z0-9])(?:%s)(?:_|\\_)v41\.(?:tex|md|json)|step3(?:_|\\_)v41(?:_|\\_)main' % '|'.join(map(re.escape, C.STEMS)), t)))
    # v42: the history sentences "text v41 ... step3_v41_main.tex said" in erratum Part 7 are excluded above; elsewhere none may remain
    assert not left, (bn, left)
# ---- v42 content checks (Lean replaces DKS)
NONHIST = [f for f in out if os.path.basename(f).endswith('.tex') and not os.path.basename(f).startswith('erratum')]
for f in NONHIST:
    t = open(f).read(); bn = os.path.basename(f)
    # the one switch: every file that uses \LeanSthree defines it (providecommand) exactly once, or gets it from b12_preamble in step3_main
    if '\\LeanSthree' in t.replace('\\providecommand{\\LeanSthree}', ''):
        if bn == 'step3_v42_main.tex': assert '\\providecommand{\\LeanSthree}' in open(D + 'b12_preamble_v42.tex').read() and '\\input{b12_preamble_v42}' in t, bn
        else: assert t.count('\\providecommand{\\LeanSthree}') == 1, bn
    assert '\\providecommand{\\LeanSthree}' not in t or t.count('\\providecommand{\\LeanSthree}') == 1, bn
# d72b744 / "read from source" may survive only in sentences that say what text v41 did
for f in NONHIST:
    for line in open(f).read().split('\n'):
        if line.lstrip().startswith('%'): continue
        for pat in ('d72b744', 'read from source', 'not elaborated', 'DKS re-read', 're-read DKS'):
            if pat in line:
                assert re.search(r'[Tt]ext v4[01]|text v39|v41|no longer|void|Historical|historical|v39 text|v38 status|Text v38|Text v39|Text v40', line), (os.path.basename(f), pat, line[:160])
# ---- the files touched by edits
assert {k for k, v in per.items() if v} == {'U1_regularity_v41.tex', 'C1_identification_v41.tex', 'continuation_v41.tex', 'bootstrap_v41.tex',
        'step3_v41_main.tex', 'b12_preamble_v41.tex', 'step4_handoff_v41.tex', 'ei_lemma_v4c_v41.tex', 'GS_lemma_v41.tex', 'lemmaU_v41.tex',
        'main_v3_v41.tex'}, per
print('built %d files, %d edits' % (len(out), nedits))
for k in TEX: print('  edits in %s: %d' % (k, per[k]))
