# LaTeX compile of text set v44 (Nov 24, 2026)

Compiler: MiKTeX on Micah's Windows machine (`MiKTeX-pdfTeX 4.26 (MiKTeX 26.2)`), called from WSL through `cmd.exe`, one document at a time,
nice 19. Recipe `latex_logs/run_v44.bat`: `pdflatex -interaction=nonstopmode -halt-on-error <doc>.tex`, then `bibtex <doc>` if the .aux has
`\bibdata`, then pdflatex twice more. **No command-line workaround**: no `\providecommand` on the command line, no scratch edit of any source.

Scratch directory: `C:\Users\micahw\kitruns\latex_v44` (new, empty before; never a WSL-UNC path). Inputs copied there: every `*.tex` of
`fix27_texts_v44/` and `fig_curves.pdf` (byte-identical to the set's files, checked with `cmp` after the compile), and
`writeup/main/main.bib` (for main_v3's `\bibliography{main}`; same sha256 as the source, 1c1b5601...). fix27_texts_v44 itself was not
compiled in place. PDFs: `pdf/` (`pdf/SHA256SUMS`); logs: `latex_logs/` (+ `main_v3_v44.blg`, `version.txt`). PDFs carry timestamps, so
they are not byte-reproducible.

**Result: all 15 standalone documents of text v44 compile (exit code 0 on the final pass, PDF written, no LaTeX error, no undefined
citation, no "Rerun" warning). The assembled `step3_v44_main` has 56 pages, 0 errors, 0 undefined references, 0 undefined citations.**
The three defects of the v43 compile (`fix26_latex_v43/LATEX_v43.md`) are gone: `\key` is defined where used; the table of contents of
step3_main reads back (short title of the Lemma U section); main_v3 finds `fig_curves.pdf` beside it.

| document | pages | errors | undefined refs (standalone) | overfull boxes |
|---|---|---|---|---|
| step3_v44_main (assembly) | 56 | 0 | 0 | 97 |
| lemmaU_v44 | 27 | 0 | 3 (BS:def, BS:prop:close, BS:sec:order) | 71 |
| main_v3_v44 | 21 | 0 | 0 | 4 |
| U1_regularity_v44 | 11 | 0 | 0 | 19 |
| C1_identification_v44 | 8 | 0 | 0 | 11 |
| GS_lemma_v44 | 6 | 0 | 0 | 10 |
| G1_export_k2_v44 | 6 | 0 | 0 | 7 |
| E0_LEMMA_v4c_v44 | 3 | 0 | 0 | 6 |
| erratum_v44 | 27 | 0 | 0 | 64 |
| step4_handoff_v44 | 4 | 0 | 0 | 14 |
| bootstrap_v44 | 6 | 0 | 7 (C:def:G, C:lem:G, C:prop:global, EI:sec:hyp, F:prop:abs, F:rem:boot, M:thm:W) | 8 |
| continuation_v44 | 8 | 0 | 2 (BS:def, BS:prop:close) | 12 |
| ei_lemma_v4c_v44 | 10 | 0 | 1 (BS:def) | 13 |
| flux2_lemma_v44 | 9 | 0 | 2 (BS:def, BS:prop:close) | 13 |
| lemmaS_A_v44 | 19 | 0 | 0 | 36 |

- Every standalone undefined reference is a label defined in another file of the set (bootstrap_v44, continuation_v44, ei_lemma_v4c_v44,
  flux2_lemma_v44, step3_v44_main; checked by grep) and resolves in the assembly; the counts equal those of the v43 compile.
- "errors" = log lines of `! LaTeX Error`, `! Undefined control sequence`, `! Missing`, `! Package ... Error`, `Emergency stop`, TeX capacity;
  none. (Four log lines starting with `!` in step3_main, U1_regularity and lemmaS_A are text dumps of overfull boxes, not errors; v43 has the same.)
- Overfull-box warnings (cosmetic): 385 in total (v43: 384). Other warnings: one float specifier `h` -> `ht`, five "No \author given".
- Page changes against v43: step3_main 56 (55), U1_regularity 11 (10), erratum 27 (26); the added v44 text.
- Not done: visual proofreading of the PDFs; main.bib still has no Lean entry (main_v3 cites none).
- This compile is recorded here only; the v44 texts themselves still say "text v44 not yet compiled" (they were built before it).
