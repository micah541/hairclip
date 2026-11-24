# PATCH_REPORT_v2: step4_proof_v2.tex (24 Nov 2026)

New file `step4_proof_v2.tex` = closed `writeup/step4_proof/step4_proof.tex` (sha256 6b4c183a...) + edits P1-P10 of
`lemmaU/bulk/fix27_texts_v44/STEP4_PATCH_NOTE_v44.md` (sha256 8bb6bfe3...), approved by Micah Warren on 24 Nov 2026.
Built by `tools/apply_patch_v2.py` (asserts base and note sha256; each old string occurs exactly once in the base and in the
working text); record `apply_patch_v2.json` (args, per-edit old/new sha256, base line). Output sha256
a972c349e75265731527ff54942bfb03a15414f9c4030094b826874e13b6006c; two runs byte-identical. Plain diff: `tools/diff_v2.txt`.
Step 3 is NOT closed; this patch changes only step-4 wording/citations; no number of step 4 changed.

## Edits
| id | closed line | status | notes |
|---|---|---|---|
| P1 summary box | 53 | adapted | note wording; `step4_handoff_v44.tex` given with its directory `writeup/lemmaU/bulk/fix27_texts_v44/` (the step-4 doc otherwise has no way to find it). |
| P2 summary table | 89 | as note | |
| P3 Thm 4(c) (H6) + new (GS) | 119 | adapted | "symmetric" -> "symmetric about the $x$-axis" (note v40 F-3 / handoff item (3)); C-1 file given with directory; (GS) typeset as `\item[(GS)]` in the (H\arabic*) list; manifest given with path `step4/p2/out_v38/`; tag `\tC` added (stage 23 is [cert]). Numbers checked against `out_v38/iface/GS_tau1.json`: sup\|ku\| 0.0065169067 <= 0.0065170, sup\|u'\| 0.0541816344 <= 0.054182; the certificate covers \|x\| <= 60.25 (superset of 60). |
| P4 status (b) | 130 | adapted | "commit pending" replaced: commit 3130aef39829f02c126aeb166fbb9583f7670108 of the authors' local repository CWFlow (public location pending); full `lake build CWFlow.Curve.All` 24 Nov 2026, 10316 jobs, no errors, no `sorry`, log `writeup/lean_bridge/full_build_Curve_All_3130aef.log` (sha256 6f448b64...; tail "Build completed successfully (10316 jobs)", 0 matches for "error"/"sorry"). `cdf_symmetric` is in `CWFlow/Curve/N1.lean` (port), so "of CWFlow/Curve/" is correct. |
| P5 status (c) | 131 | adapted | "[cert]" -> `\tC`; GS lemma file with directory. Labels exist: GS-1 = `GS:lem:tau1` in GS_lemma_v44.tex, C-1(iv) in C1_identification_v44.tex. |
| P6 P2 status | 336 | as note | old string incl. trailing newline, as in the note. |
| P7 P2 Step 7 | 370 | adapted | the note's "Proposition C1:prop:P2" is a LaTeX label; written as the Proposition "identification on P2" (label `C1:prop:P2`) so the PDF reads sensibly. |
| P8 (H6) section | 672 | as note | |
| P9 open item O4 (Kit B) | 798 | adapted | note text + its v44 addendum (Windows runs, kit v5) + instructed status: kit_step4b_v4 Windows 159/0/1 (kit4bv4_run_win1), kit_step3_v4 Windows (kit3v4_run_win1: 72 PASS / 6 FAIL CRLF-only / 3 WARN), kit_step3_v5 Windows 78/0/3 (kit3v5_run_win1), kit v5 hashes from `kit_v5_build/KIT_V5_REPORT.md`; kit_step4b_v5 Windows run "in progress, no result yet" marked `\OPEN`. The old string was extended to include the following sentence "Until then the P2 and T9 chain is reproduced only by ... sandbox", which would be false after the edit; it is removed. Hashes split into 16-char `\fn` chunks (overfull line otherwise). |
| P10 appendix Kit B | 874 | adapted | "As P9" -> "See open item O4 (\S\ref{s:open})"; old string extended by "It will cover" -> "It covers" (kit 4B v4 README: base store, T6/t6c, T7, final tail, H1 stages). |

## Checks run
- base sha256 before/after = 6b4c183a... ; all 11 files of `writeup/step4_proof/` sha256 identical before/after (closed_before vs after).
- apply script: 10 exact-once assertions in base and working text; double run byte-identical.
- numbers in P3/P9 checked against GS_tau1.json, C1_ke.json (stage 22: sup\|k~e\| 0.0017428), manifest_v38 stage list (22 = C-1, 23 = GS tau1; sha be909bbf...), KIT_V5_REPORT.md, RUN_NOTES.md of the three Windows runs.
- MiKTeX compile of base and v2: see `LATEX_v2.md` (v2: 40 pp, 0 undefined refs/cites, 1 pre-existing error shared with the closed file).
- Nothing outside `writeup/step4_proof_v2/` and `C:\Users\micahw\kitruns\latex_step4_v2\` changed (find -newer start stamp over /home/claude, excluding these; see final section).

## Not done / left as is (outside the note; for the authors)
1. Status box l.55 still says "kit~B: pending" in its Open list, and the status-box header still says "version 2, Nov~12, 2026": now inconsistent with P9 / with this file. Not in the note, not changed.
2. Hand-off table row l.146 "(H6) & smooth, closed, immersed & Theorem~4(b) & \tA": still [A]; not in the note.
3. l.789 "$T$ is the maximal time by ... (P1) \tA" and the dependency paragraph l.791 ("literature facts (P1)"): deliberately unchanged per the note (step 3's (P1')).
4. Theorem 4(c) still says "v3.1 as amended by v3.2"; the note's v40 items (F-7 (H4) range [14,700] certified; v3.3 items in Thm 4(c)) are notes for a future version, not edits.
5. Pre-existing literal `\n` at closed l.668 (v2 l.679): compile error under -halt-on-error in both files; not fixed.
6. The PDF has not been proofread visually.
7. Lean: no public remote; cited as local commit, "public location pending". kit_step4b_v5 Windows result not available.
8. The closed dir contains an unrelated older `step4_proof_v2.tex` (Nov 12 draft); the name collision is noted in the new file's header.
