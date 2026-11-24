# LATEX_v3: MiKTeX compile of step4_proof_v3.tex (24 Nov 2026)

Engine: MiKTeX-pdfTeX 3.141592653-2.6-1.40.29 (MiKTeX 26.2), Micah's Windows machine, run from WSL via cmd.exe
(`nice -n 19 cmd.exe /c run.bat step4_proof_v3`), one pdflatex at a time (a kit run was active at below-normal priority).
Scratch dir (new): `C:\Users\micahw\kitruns\latex_step4_v3\` (`run.bat` = the v2 recipe, copied; `halt_logs\`).
Recipe: pdflatex `-interaction=nonstopmode -halt-on-error`, bibtex only if `\bibdata` in .aux (not the case: inline
thebibliography), pdflatex x2. (A first attempt via `start /belownormal` returned "Access is denied" before pdflatex
ran and produced no files; the plain invocation of the v2 recipe was used.)

| document | mode | result | pages | `!` errors | undefined refs | undefined cites | LaTeX warnings | overfull hbox | underfull |
|---|---|---|---|---|---|---|---|---|---|
| step4_proof_v3.tex (sha256 e8683090...) | `-halt-on-error` | PDF written (`step4_proof_v3.pdf`, sha256 ad2b1d9279db63b0...) | 44 | 0 | 0 | 0 | 0 | 0 | 59 |

v2 under -halt-on-error was fatal (the literal `\n`, removed in v3, edit S14). Log: `latex_logs/v3_halt.log`
(final pass; `Output written on step4_proof_v3.pdf (44 pages, 943376 bytes)`; rerunfilecheck: `.out` unchanged).
The new labels resolve (`ss:v33status` = 10.1, p. 35; `ss:mainv32` = 10.2, p. 36).
The PDF was not proofread visually (no PDF renderer in WSL); log checks only.
