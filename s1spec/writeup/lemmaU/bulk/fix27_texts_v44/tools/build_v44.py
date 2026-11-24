# build_v44.py -- builds the step-3 text set v44 in fix27_texts_v44/ from the v43 set in fix24_texts_v43/ (read-only).
# (v44: follows fix24_texts_v43/tools/build_v43.py.)  Work: (A) the check writeup/lemmaU/bulk/fix25_check_v43/CHECK_v43.md of text v43
# (F-1 records, F-2..F-6 wording, the Note: Windows runs of both v4 kits); (B) the LaTeX defects of fix26_latex_v43/LATEX_v43.md
# (\key defined where used; short titles for titles with \fn; fig_curves.pdf beside main_v3); (C) records/open lists.
# The replacements are in tools/edits_v44.txt (each asserted to occur exactly once).  Text only: no number is computed here.
# Erratum: Parts 1-7 verbatim, Part 8 (section S) with one marked correction (F-2), Part 9 (section T) from tools/erratum_part9_v44.tex;
# records addition from tools/records_errata_add_v44.md (and one marked correction in the v43 section, F-2).
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v44.py   (after gen_records_v44.py, gen_numbers_A_v44.py, gen_ei_constants_v44.py)
import sys, os, json, re, shutil, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b44_common as C
from b44_common import Builder, D, D9
R = json.load(open(D + 'tools/records_v44.json')); P = R['pins']
assert (P['entries'], P['distinct_files_normalised'], P['stage_paths_distinct'], P['scripts'], P['scripts_also_stage_paths']) == (992, 646, 629, 55, 38), P
assert R['kit3_win_summary'] == '72 PASS, 6 FAIL, 3 WARN' and R['kit4b_win_summary'] == 'TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]'
# ---- (B3) the figure of main_v3: byte copy of writeup/main/fig_curves.pdf
FIG = '/home/claude/s1spec/writeup/main/fig_curves.pdf'
shutil.copyfile(FIG, D + 'fig_curves.pdf')
assert hashlib.sha256(open(D + 'fig_curves.pdf', 'rb').read()).hexdigest() == R['args']['inputs']['fig_curves']['sha256']
def sub(s):
    assert not re.search(r'@[A-Z0-9_]+@', s), s[:200]
    return s
# ---- parse edits_v44.txt
EDITS = {}; ORDER = []
cur = None; mode = None; buf = []; old = None; kind = None
for line in open(D + 'tools/edits_v44.txt').read().split('\n'):
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
HDR = ('%% %s -- text v44 (fix27_texts_v44, 24 Nov 2026), built by fix27_texts_v44/tools/build_v44.py from fix24_texts_v43/%s (unchanged):\n'
       '%% the check fix25_check_v43/CHECK_v43.md (F-1..F-6, Note) and the LaTeX defects of fix26_latex_v43/LATEX_v43.md, where they touch this file.  Text v44 not yet compiled.')
TEX = ['bootstrap_v43.tex', 'flux2_lemma_v43.tex', 'lemmaU_v43.tex', 'U1_regularity_v43.tex', 'continuation_v43.tex', 'GS_lemma_v43.tex',
       'step3_v43_main.tex', 'main_v3_v43.tex', 'step4_handoff_v43.tex', 'lemmaS_A_v43.tex', 'G1_export_k2_v43.tex',
       'C1_identification_v43.tex', 'E0_LEMMA_v4c_v43.tex', 'b12_preamble_v43.tex', 'chart_lemma_body_v43.tex', 'ei_lemma_v4c_v43.tex']
assert set(ORDER) <= set(TEX) | {'erratum_v43.tex'}, set(ORDER) - set(TEX)
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
    b = Builder(name, header=HDR % (C.v44name(name), name))
    apply(b, name)
    out.append(b.write())
# ---- the step-4 patch note (markdown): retag, plus a v44 section
b = Builder('STEP4_PATCH_NOTE_v43.md')
add = ('\n\n## Text v44 (24 Nov 2026; check fix25_check_v43/CHECK_v43.md, LaTeX compile fix26_latex_v43/LATEX_v43.md)\n'
       '- No edit is added or changed. The proposed edits still need the authors\' go before a new step-4 version is written.\n'
       '- When the proposed kit item is applied, it should also record the Windows runs of 24 Nov 2026 (`kit3v4_run_win1/RUN_NOTES.md`, '
       '`kit4bv4_run_win1/RUN_NOTES.md`: both reproduce after Windows-only fixes in the run copies) and the kit v5 that carries those fixes '
       '(pending; no hash yet), as `step4_handoff_v44.tex` section "Kits" does.\n')
b.t = b.t.rstrip('\n') + add
C.INSERTS.append({'file': b.name, 'text': add})
out.append(b.write())
# ---- erratum: Parts 1-7 verbatim; Part 8 with one marked correction (F-2); Part 9 new
p9 = sub(open(D + 'tools/erratum_part9_v44.tex').read())
b = Builder('erratum_v43.tex', retag_names=False, header='% erratum_v44.tex -- text v44: erratum_v43.tex (fix24_texts_v43, unchanged) with Parts 1-7 verbatim, Part 8 (section S) with one correction marked [Text v44, check fix25 F-2:], \\key defined (LaTeX compile fix26, defect 1), plus Part 9 (section T, text v44).  Built by fix27_texts_v44/tools/build_v44.py.  Text v44 not yet compiled.')
apply(b, 'erratum_v43.tex')
b.rep('\\end{document}', p9.rstrip('\n') + '\n\\end{document}')
out.append(b.write())
md = open(D9 + 'ERRATUM_v43.md').read()
open(D + 'ERRATUM_v44.md', 'w').write('# ERRATUM_v44 (text v44): ERRATUM_v43.md verbatim, then a pointer to Part 9\n\nPart 9 (text v44: the check fix25_check_v43/CHECK_v43.md of text v43 and the LaTeX compile fix26_latex_v43/LATEX_v43.md) is section T of `erratum_v44.tex`; its source is `tools/erratum_part9_v44.tex`. Part 8 (section S) of `erratum_v44.tex` carries one correction marked [Text v44, check fix25 F-2:].\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v44.md')
rec = open(D9 + 'RECORDS_ERRATA_v43.md').read(); rec44 = sub(open(D + 'tools/records_errata_add_v44.md').read())
o = '- No step-3 or step-4 record changed; nothing was rerun or re-pinned; no Lean file was built or edited. The records items above are unchanged.\n- (L-4'
assert rec.count(o) == 1
nw = ('- No step-3 or step-4 record changed; nothing was rerun or re-pinned; no Lean file was edited [text v44, check fix25 F-2: the three S3 files were '
      'recompiled per file, sources unchanged, to make the logs of L-4, which end with return code 0]. The records items above are unchanged.\n- (L-4')
rec = rec.replace(o, nw)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v44.md', 'text': nw})
open(D + 'RECORDS_ERRATA_v44.md', 'w').write('# RECORDS_ERRATA_v44: RECORDS_ERRATA_v43.md verbatim (one marked correction in its v43 section), then the v44 section at the end\n\n' + rec + '\n\n' + rec44)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v44.md', 'text': rec44})
out.append(D + 'RECORDS_ERRATA_v44.md')
C.save_inserts()
# ---- sanity: no 0.1475 quote outside the history files; no v43 set-file name left outside them
for f in out:
    t = re.sub(r'fix24_texts_v43/[\w.\\-]+', ' ', open(f).read()); bn = os.path.basename(f)   # provenance paths into fix24 are intended
    if bn.startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
    left = sorted(set(m.group(0) for m in re.finditer(r'(?<![A-Za-z0-9])(?:%s)(?:_|\\_)v43\.(?:tex|md|json)|step3(?:_|\\_)v43(?:_|\\_)main' % '|'.join(map(re.escape, C.STEMS)), t)))
    assert not left, (bn, left)
# ---- v44 content checks
ALLT = {os.path.basename(f): open(f).read() for f in out}
NONHIST = [bn for bn in ALLT if bn.endswith('.tex') and not bn.startswith('erratum')]
U1 = ALLT['U1_regularity_v44.tex']
assert 'anywhere in the import closure' not in U1                                   # F-3
assert 'For $X:\\mathbb R\\supset U\\to\\mathbb C$ an immersion' not in U1           # F-5
assert 'tools/check_numbers_v44.py} of this text set' in U1
c1 = ALLT['C1_identification_v44.tex']; i = c1.index('\\item \\textbf{Lean} (text v42): the theorems used are compiled')
assert 'fix23_delta_check/DELTA_CHECK_v42.md} Part~A' in c1[i:i + 600]               # F-4
assert 'commit pending; \\fn{cdf_exists_maximal}; the regularity lemma' not in ALLT['main_v3_v44.tex']   # F-6
assert 'no Lean file was built or edited. Parts 1--6' not in ALLT['erratum_v44.tex']                    # F-2
assert 'fix17_texts_v39/tools/gen_numbers_A' not in ALLT['lemmaS_A_v44.tex']
# Note: no non-history text says that the Windows runs are open or were not made
for bn in NONHIST:
    t = ALLT[bn]
    for ph in ("Windows runs of both v4 kits are open", "the authors' independent Windows runs are open", "No Windows run of either v4 kit has been made",
               "\\textbf{Open:} the authors' independent Windows runs", "the authors' Windows runs are open", "the authors' Windows runs of the two v4 kits; LaTeX compilation; the step-4",
               "LaTeX compile} of every v43 file", "(a$'$) the authors' Windows runs"):
        assert ph not in t, (bn, ph)
    assert 'kit v5' not in t or 'pending' in t, bn
# (B1) every standalone document that uses \key defines it (also checked by check_macros_v44.py); (B2) titles: check_titles_v44.py
for bn in ('GS_lemma_v44.tex', 'erratum_v44.tex', 'flux2_lemma_v44.tex'):
    assert ALLT[bn].count('\\providecommand{\\key}[1]{\\texttt{\\detokenize{#1}}}') == 1, bn
assert '\\section[Lemma U (summary)]' in ALLT['step3_v44_main.tex']
# every status-bearing file cites the check fix25 or the compile fix26
for bn in ('step3_v44_main.tex', 'main_v3_v44.tex', 'erratum_v44.tex'):
    assert 'fix25_check_v43/CHECK_v43.md' in ALLT[bn] and 'fix26_latex_v43/LATEX_v43.md' in ALLT[bn], bn
assert set(k for k, v in per.items() if v) == {'GS_lemma_v43.tex', 'flux2_lemma_v43.tex', 'step3_v43_main.tex', 'ei_lemma_v4c_v43.tex', 'G1_export_k2_v43.tex',
        'main_v3_v43.tex', 'U1_regularity_v43.tex', 'C1_identification_v43.tex', 'lemmaS_A_v43.tex', 'step4_handoff_v43.tex', 'erratum_v43.tex'}, per
print('built %d files, %d edits; fig_curves.pdf copied' % (len(out), nedits))
for k in TEX + ['erratum_v43.tex']: print('  edits in %s: %d' % (k, per[k]))
