# LATEX_v2: MiKTeX compile of step4_proof_v2.tex and baseline step4_proof.tex (24 Nov 2026)

Engine: MiKTeX-pdfTeX 3.141592653-2.6-1.40.29 (MiKTeX 26.2), Micah's Windows machine, run from WSL via cmd.exe, nice 19.
Scratch dir (new): `C:\Users\micahw\kitruns\latex_step4_v2\` (`base\`, `v2\`, `halt_logs\`, `run.bat`, `run_nohalt.bat`).
Recipe (as fix26_latex_v43/logs/runp.bat): pdflatex, bibtex only if `\bibdata` in .aux (not the case: inline thebibliography), pdflatex x2.
Logs copied to `latex_logs/` here.

| document | mode | result | pages | `!` errors | undefined refs | undefined cites | overfull hbox |
|---|---|---|---|---|---|---|---|
| step4_proof.tex (closed, scratch copy, sha 6b4c183a...) | `-halt-on-error` | **fatal, no PDF** | - | 1 | - | - | - |
| step4_proof.tex | nonstopmode, no halt | PDF written | 39 | 1 (same) | 0 | 0 | 0 |
| step4_proof_v2.tex (sha a972c349...) | `-halt-on-error` | **fatal, no PDF** | - | 1 (same) | - | - | - |
| step4_proof_v2.tex | nonstopmode, no halt | PDF written (`step4_proof_v2.pdf`, sha 81a43c2e...) | 40 | 1 (same) | 0 | 0 | 0 |

The one error is PRE-EXISTING in the closed source: a literal `\n` after "This item is closed." at the end of Remark rem:t1
(closed l.668, v2 l.679) -> "Undefined control sequence". TeX drops it; output is otherwise unaffected. The closed PDF of Nov 12
(TeX Live 2023) was built with the same error (closed `step4_proof.log` l.729; 39 pages). It is NOT fixed in v2 (not in the note);
a one-character fix (delete `\n`) is proposed for a later version.
An overfull line (101.8 pt) caused by the full kit sha256 in `\texttt` was removed by splitting each hash into 16-char `\fn` chunks
with `\allowbreak` (P9). No other warnings of type undefined reference/citation or "Rerun".
The PDF was not proofread visually (no PDF text tool in WSL); log checks only.
