# LaTeX compile of text set v43 (Nov 24, 2026)
Compiler: MiKTeX on Micah's Windows (MiKTeX-pdfTeX 4.26 (MiKTeX 26.2)); pdflatex + bibtex + 2 pdflatex, nonstopmode, in scratch copy C:\Users\micahw\kitruns\latex_v43 (fix24_texts_v43 untouched). PDFs in pdf/ (SHA256SUMS), logs in logs/.

**Result: all 15 standalone documents compile once 3 source defects are worked around; the assembled step3_v43_main compiles to 55 pp with 0 errors, 0 undefined references, 0 undefined citations.**

## Source defects (to fix in v44; all already present in v42)
1. `\key` used but not defined in GS_lemma_v43.tex (l.70), erratum_v43.tex (l.646), flux2_lemma_v43.tex (l.37, 4 uses) -> fatal. Worked around with `\providecommand{\key}` on the command line. (check_macros_v40.py missed it.)
2. step3_v43_main.tex:73 section title contains `\fn{fix24_texts_v43/lemmaU_v43.tex}`; the underscore breaks in the .toc ("Missing $") -> fatal on pass 2. Fix: optional short title `\section[...]{...}` (used in scratch copy step3_v43_main_scratch.tex).
3. main_v3_v43.tex needs fig_curves.pdf (from writeup/main/) next to it.

## Per document (pages; undefined refs when standalone)
C1_identification 8; E0_LEMMA 3; G1_export_k2 6; GS_lemma 6; U1_regularity 10; bootstrap 6 (7 refs); continuation 8 (2); ei_lemma_v4c 10 (1); erratum 26; flux2_lemma 9 (2); lemmaS_A 19; lemmaU 27 (3); main_v3 21; step4_handoff 4; step3_v43_main 55 (0).
Every standalone undefined reference is a label defined in another file of the set (checked) and resolves in the assembly.
Overfull-box warnings (cosmetic): 384 in total across logs.
Not done: main.bib has no Lean entry (main_v3 cites none yet); visual proofreading of PDFs.
