#!/usr/bin/env python3
"""apply_sync_v3.py: build step4_proof_v3.tex (status-synced version 3 of the step-4 write-up) from
writeup/step4_proof_v2/step4_proof_v2.tex (patch v2, approved 24 Nov 2026).  Version 3 approved by Micah Warren
on 24 Nov 2026 (status sync, new files only).

Reads the v2 file only (never writes it).  Asserts the base sha256, asserts that every old string occurs exactly
once in the base and once in the working text, and writes step4_proof_v3.tex + apply_sync_v3.json (args block,
per-edit record).  Every edit and its source is listed in SYNC_REPORT_v3.md.
Usage: nice -n 19 python3 -I apply_sync_v3.py
"""
import hashlib, json, os, sys

BASE = "/home/claude/s1spec/writeup/step4_proof_v2/step4_proof_v2.tex"
BASE_SHA = "a972c349e75265731527ff54942bfb03a15414f9c4030094b826874e13b6006c"
OUTD = "/home/claude/s1spec/writeup/step4_proof_v3"
OUT = os.path.join(OUTD, "step4_proof_v3.tex")
REC = os.path.join(OUTD, "apply_sync_v3.json")
F27 = "writeup/lemmaU/bulk/fix27_texts_v44/"
# sources read for the new text (hashed into the record; not edited)
SOURCES = [
    "/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/main_v3_v44.tex",
    "/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/step3_v44_main.tex",
    "/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/step4_handoff_v44.tex",
    "/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/C1_identification_v44.tex",
    "/home/claude/s1spec/step4/STEP4_RESULT.md",
    "/home/claude/s1spec/step4/p2/out/FINAL_manifest.json",
    "/home/claude/s1spec/step4/p2/out/FINAL_manifest_cleanup_pass2.json",
    "/home/claude/s1spec/step4/p2/out_v33/FINAL_manifest_v33.json",
    "/home/claude/s1spec/step4/p2/out_v38/FINAL_manifest_v38.json",
    "/home/claude/s1spec/step4/p2/out_v38/REPIN_REPORT.md",
    "/home/claude/s1spec/step4/p2/out_v38/iface/C1_ke.json",
    "/home/claude/s1spec/step4/p2/v33/EXCEPTION_v32_B8_check_outA5.json",
    "/home/claude/s1spec/writeup/lean_bridge/full_build_Curve_All_3130aef.log",
    "/home/claude/s1spec/kit4bv4_run_wsl1/RUN_NOTES.md",
    "/home/claude/s1spec/kit3v4_run_wsl5/RUN_NOTES.md",
    "/home/claude/s1spec/kit3v4_run_win1/RUN_NOTES.md",
    "/home/claude/s1spec/kit4bv4_run_win1/RUN_NOTES.md",
    "/home/claude/s1spec/kit3v5_run_win1/RUN_NOTES.md",
    "/home/claude/s1spec/kit_step4a_v3.zip",
]


def hx(h):
    """a 64-hex sha256 as four \\fn chunks joined by \\allowbreak"""
    assert len(h) == 64
    return r"\allowbreak".join(r"\fn{" + h[i:i + 16] + "}" for i in range(0, 64, 16))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


assert sha(BASE) == BASE_SHA, "base sha256 mismatch"
base = open(BASE, encoding="utf-8").read()

V33STATUS = (
    r"\subsection{Current status: the v3.3 line of step~3 (Nov~24, 2026)}\label{ss:v33status}" "\n"
    r"\textbf{Step~3 is not closed.} The current step-3 texts are the v3.3 line, text v44 (\fn{" + F27 + r"}): the "
    r"assembled step~3 \fn{step3_v44_main.tex} and the main-theorem document \fn{main_v3_v44.tex}. In this line "
    r"Hypothesis~(F) of version~2 (and the further hypothesis (O), the outer cutoff of Lemma~S(A) on $[350,700]$) is "
    r"\emph{no longer assumed}: the far-field estimate is Lemma~U, a uniformly local Lemma~S(B) whose energies are "
    r"centred at every physical radius and keep the edge and ramp terms (\fn{lemmaU_v44.tex}), and the outer cutoff is "
    r"treated in the repaired Lemma~S(A) (\fn{lemmaS_A_v44.tex}; erratum \fn{erratum_v44.tex}). The main theorem of "
    r"\fn{main_v3_v44.tex} states no hypothesis about the far field, but it is conditional on the open step-3 items: the "
    r"authors' confirmations ((P3$'$), (P2$'$) and the frame conventions, Lemma~O arc by arc, the single-point "
    r"statement); the public location of the Lean commit \fn{3130aef39829f02c126aeb166fbb9583f7670108} of the authors' "
    r"local repository \fn{CWFlow} (built in full on 24~Nov~2026: \fn{lake build CWFlow.Curve.All}, 10316 jobs, no "
    r"errors, no \fn{sorry}; log \fn{writeup/lean_bridge/full_build_Curve_All_3130aef.log}); the \texttt{[read]} inputs "
    r"and \texttt{[hand]} steps listed in the step-3 texts; the records items (\fn{RECORDS_ERRATA_v44.md}); and a final "
    r"review. Its values in (i) and (iv) are those of the verifier \fn{fin13v4c}: $\bar a\le2.14039\cdot10^{-7}$, tail "
    r"rate $\varrho\ge2.6792\cdot10^{-3}$, floor $\sqrt{\Lam_F}\le9.8972\cdot10^{-10}$, eventually "
    r"$\abs\alpha\le6.8824\cdot10^{-7}$; the values of \S\ref{ss:mainv32} ($2.19103\cdot10^{-7}$, $2.69\cdot10^{-3}$, "
    r"$4.94\cdot10^{-10}$, $3.43\cdot10^{-7}$) belong to the v3.2 line and do not apply to v3.3. In the v3.3 line, "
    r"(P1$'$) (local well-posedness with uniqueness, the maximal flow and the $L^2$ blow-up criterion) rests on compiled "
    r"Lean theorems of \fn{CWFlow/Curve/} (\fn{cdf_exists_maximal} and others) together with the regularity lemma U-1 "
    r"(\fn{U1_regularity_v44.tex}), and no longer on the literature fact (P1) \tA\ used in \S\ref{ss:mainv32}." "\n\n"
    r"\emph{What step~4 delivers to the v3.3 line} (\fn{step4_handoff_v44.tex}; the values of Theorem~4 are unchanged):" "\n"
    r"\begin{itemize}" "\n"
    r"\item (H1)--(H5) of Theorem~4(c), with the v3.2 values;" "\n"
    r"\item the tightened (H2) in the variable $f-U_1$ on $y\in[14,60]$: $(3.658\cdot10^{-10},1.777\cdot10^{-9},"
    r"1.665\cdot10^{-6})$ against $(10^{-8},10^{-6},10^{-3})$, margins $\times27.34$, $\times562.8$, $\times600.9$; and "
    r"(H2f) on $[60,700]$, with error $(3.137\cdot10^{-7},1.004\cdot10^{-6},2.436\cdot10^{-6})$ against "
    r"$(10^{-6},5\cdot10^{-6},10^{-5})$, margins $\times3.188$, $\times4.983$, $\times4.106$ \tC;" "\n"
    r"\item the fixed-$L$ cone family (jets, (H5) data, Lemma~O H4 and H1) and its time ratio H7: ratio "
    r"$\le1.6000014\cdot10^{-7}\le1.6001\cdot10^{-7}$, $t^*\ge3.68334\cdot10^{-7}$, margin $\times2.302$ \tC;" "\n"
    r"\item the far-family Lemma~O H4 at every centre of radius $\ge490\lambda_1$, minimum margin $\times15.58$ \tC;" "\n"
    r"\item (H4) on $y\in[14,700]$ (step~3 v3.3 uses $[48,700]$; Theorem~4(c) states $[48,350]$) \tC;" "\n"
    r"\item (H6) by Theorem C-1 of \fn{C1_identification_v44.tex} (\S\ref{ss:H6reg}) \tA, with "
    r"$\sup|\tilde\kappa e|\le0.0017429$ \tC; and (GS) of Theorem~4(c) \tC." "\n"
    r"\end{itemize}" "\n"
    r"These items are pinned by \fn{FINAL_manifest_v33.json} and \fn{FINAL_manifest_v38.json} (\S\ref{app:repro}). "
    r"Step~3 recomputes the Lemma~O H6 closure with the Lemma~U strip input (smallest factor $\times1.86$, far $K_3$); "
    r"this supersedes the v3.2 closure rows of Table~\ref{tab:v32}. Step~4's own results depend neither on (F) nor on "
    r"Lemma~U." "\n\n"
    r"\subsection{The main theorem of version 2 (v3.2 line; history)}\label{ss:mainv32}" "\n"
    r"\emph{This subsection is kept as written in version~2 (Nov~12); it is superseded by \S\ref{ss:v33status}.} "
)

V33_MANIFEST_BULLET = (
    r"\item \emph{H1 follow-up rerun (O6(c)).} \fn{final_tail.py --from 'H1 exact points'} was run twice (stages 0--3 "
    r"reused and verified); pass~2 reproduces pass~1 byte for byte apart from the stages' seconds, the \fn{cpu} field of "
    r"the 16 H1 outputs and timing lines in 6 logs (\fn{checks.h1rot_second_pass}; pass manifests "
    r"\fn{FINAL_manifest_h1rot_pass1.json}, \fn{FINAL_manifest_h1rot_pass2.json})." "\n"
    r"\item \emph{On top of it (v3.3 line).} \fn{p2/out_v33/FINAL_manifest_v33.json} (21 stages, 366 paths, SHA-256 "
    r"\fn{8bc47e83e370fec3...}; driver \fn{final_tail_v33.py}) and \fn{p2/out_v38/FINAL_manifest_v38.json} (29 stages, "
    r"992 pin entries on 646 distinct files, SHA-256 " + hx("be909bbf42543ab017bb9c81ba58227c80fc4f010cc1b50a66a515b848f90418")
    + r"; driver \fn{p2/v38/final_tail_v38.py}; report \fn{p2/out_v38/REPIN_REPORT.md}; two passes byte-identical apart "
    r"from timing). The v38 driver verifies the closed manifest at start and end with exactly one mismatch, "
    r"\fn{s1_tight/consistency/check_outA5.json} (rewritten on Nov~15 by a kit test run that escaped its folder; same "
    r"script and inputs, result unchanged: all rows inside), accepted only through the documented exception "
    r"\fn{p2/v33/EXCEPTION_v32_B8_check_outA5.json}; and manifest\_v33 with one known difference, "
    r"\fn{constants_needed_v34.json} (a circular pin). The closed manifest itself is unchanged."
)

# (id, old, new, source)
EDITS = [
    ("D1 date",
     r"\date{Version 2, November 12, 2026}",
     r"\date{Version 3 (status sync), November 24, 2026; the mathematics is that of version 2, November 12, 2026}",
     "task"),
    ("S1 status box: version, manifests",
     r"\textbf{Status} (version 2, Nov~12, 2026; it supersedes draft v1 of Nov~11, kept as \fn{step4_proof_v1.tex}). "
     r"Every step-4 computation is done in \fn{arb} with outward rounding and no sampling. The final tail is pinned by "
     r"one manifest, \fn{step4/p2/out/FINAL_manifest.json}, which records every stage with its command, environment and "
     r"the SHA-256 of its inputs and outputs; a second pass reproduces it byte for byte.",
     r"\textbf{Status} (version 3, Nov~24, 2026: a status sync of version~2 of Nov~12, 2026, whose file "
     r"\fn{writeup/step4_proof/step4_proof.tex} is closed and unchanged, with the approved patch v2 of Nov~24 "
     r"(\fn{writeup/step4_proof_v2/}); no certified step-4 bound changed; stale records (hash values, stage numbers) "
     r"are corrected; version~2 superseded draft v1 of Nov~11, kept as \fn{step4_proof_v1.tex}). Every step-4 "
     r"computation is done in \fn{arb} with outward rounding and no sampling. The v3.2 final tail is pinned by the "
     r"closed manifest \fn{step4/p2/out/FINAL_manifest.json} (34 stages, 523 paths, 26 scripts, 97 upstream entries; "
     r"SHA-256 \fn{e01cb4da71ba8253...}), which records every stage with its command, environment and the SHA-256 of its "
     r"inputs and outputs; a second pass reproduces it byte for byte. On top of it, the v3.3 stages are pinned by "
     r"\fn{step4/p2/out_v33/FINAL_manifest_v33.json} (21 stages, \fn{8bc47e83...}) and the current chain by "
     r"\fn{step4/p2/out_v38/FINAL_manifest_v38.json} (29 stages, 992 pin entries, \fn{be909bbf...}; driver "
     r"\fn{step4/p2/v38/final_tail_v38.py}), which verifies the closed manifest with one documented exception "
     r"(\fn{EXCEPTION_v32_B8_check_outA5.json}) and manifest\_v33 with one known difference (\S\ref{app:repro}).",
     "manifests on disk; REPIN_REPORT.md; step4_handoff_v44.tex sec. 1"),
    ("S2 status box: current status + history label",
     r"The amendment v3.2 has been approved by the authors and is built in \fn{writeup/step3_proof_v32/}, "
     r"\fn{writeup/lemmaO_v32/} and \fn{writeup/step3_amend_v32/}. It proves Lemma~O H6 in step~3. \emph{Step-3 dependency.}",
     r"\emph{Current status} (\S\ref{ss:v33status}): Hypothesis~(F) below is \emph{no longer assumed}. In the v3.3 line "
     r"of step~3 (texts v44, \fn{" + F27 + r"}: \fn{step3_v44_main.tex}, \fn{main_v3_v44.tex}) the far-field estimate is "
     r"Lemma~U, a uniformly local Lemma~S(B) with the edge terms kept, and the outer cutoff is treated in the repaired "
     r"Lemma~S(A). \textbf{Step~3 is not closed}: the main theorem of \fn{main_v3_v44.tex} is conditional on the open "
     r"step-3 items (\S\ref{ss:v33status}). \emph{History (version~2, Nov~12; v3.2 line).} The amendment v3.2 has been "
     r"approved by the authors and is built in \fn{writeup/step3_proof_v32/}, \fn{writeup/lemmaO_v32/} and "
     r"\fn{writeup/step3_amend_v32/}. It proves Lemma~O H6 in step~3. \emph{Step-3 dependency.}",
     "main_v3_v44.tex Remark rem:hist, s:open; step3_v44_main.tex status box"),
    ("S2b status box: main theorem (past tense)",
     r"The main theorem is therefore stated \emph{conditional on Hypothesis~(F)} (\fn{writeup/main/main.pdf}).",
     r"The main theorem of version~2 was therefore stated \emph{conditional on Hypothesis~(F)} (\fn{writeup/main/main.pdf}).",
     "history"),
    ("S2c status box: fix being evaluated",
     r"A fix, a uniformly local Lemma~S(B), is being evaluated (go/no-go tests). Under v3.1 alone, Lemma~O H6 is not met (\S\ref{ss:H6v31}).",
     r"A fix, a uniformly local Lemma~S(B), was then being evaluated; it became Lemma~U of the v3.3 line (above). Under "
     r"v3.1 alone, Lemma~O H6 is not met (\S\ref{ss:H6v31}).",
     "main_v3_v44.tex ss:gap"),
    ("S3 status box: Open list",
     r"\textbf{Open} (\S\ref{s:open}): Hypothesis~(F) of step~3 (the Lemma~S(B) gap; the main theorem is conditional on it); "
     r"the true curve's $k_{\max}$ at $t_1$ (only the base value is certified); kit~B: pending. \emph{Closed in this version:} "
     r"Lemma~O H1 for the far family, Lemma~R, the script items O6a--c.",
     r"\textbf{Open} (\S\ref{s:open}, updated Nov~24): step~3 is not closed (its open items are listed in "
     r"\S\ref{ss:v33status}; Hypothesis~(F) is no longer one of them); the true curve's $k_{\max}$ at $t_1$ (only the base "
     r"value is certified; step~3 does not use it); the Windows run of \fn{kit_step4b_v5.zip} (in progress, no result yet). "
     r"The WSL test runs and the Windows runs of the v4 kits and of \fn{kit_step3_v5.zip} are done (O4). "
     r"\emph{Closed in version~2:} Lemma~O H1 for the far family, Lemma~R, the script items O6a--c (O6(c) by the follow-up "
     r"of Nov~12).",
     "RUN_NOTES.md of the kit runs; STEP4_RESULT.md O6(c)"),
    ("S4 abstract: v3.3 items",
     r"checks the hypotheses of step~3 as amended by v3.2, including the Lemma~O data per cylinder (\S\ref{s:T9}--\ref{s:H6}).",
     r"checks the hypotheses of step~3 as amended by v3.2, including the Lemma~O data per cylinder "
     r"(\S\ref{s:T9}--\ref{s:H6}); further stages certify the items added by the v3.3 line of step~3 "
     r"(\S\ref{ss:v33status}; step~3 is not closed).",
     "step4_handoff_v44.tex"),
    ("S5 tag legend: [A] includes written hand proofs",
     r"\tA\ analytic argument, standard or from the literature;",
     r"\tA\ analytic argument, standard or from the literature, or a written hand proof in a named project file (not "
     r"machine-checked; it may use compiled Lean theorems, cited by name and commit);",
     "C1_identification_v44.tex (proof of C-1 is [hand])"),
    ("S6 conventions: v3.3 named explicitly",
     r"``Step~3'' means \fn{step3_proof.tex} v3.1 as amended by v3.2 (\S\ref{s:H6}) unless v3.1 is named.",
     r"``Step~3'' means \fn{step3_proof.tex} v3.1 as amended by v3.2 (\S\ref{s:H6}) unless v3.1 is named; the v3.3 line "
     r"of step~3 (texts v44) is named explicitly where it is meant (\S\ref{ss:v33status}).",
     "task"),
    ("S7 chain table: H6 row",
     r"step-4 closure $\times1.69$, same condition & \S\ref{s:H6}\\\bottomrule",
     r"step-4 closure $\times1.69$, same condition (v3.2 line, history). In v3.3 the strip input is Lemma~U and step~3 "
     r"recomputes H6 (smallest factor $\times1.86$; \S\ref{ss:v33status}) & \S\ref{s:H6}\\\bottomrule",
     "main_v3_v44.tex (Theorem H6 v3.3: smallest factor x1.86); manifest_v38 stage 9 (x1.86087)"),
    ("S8 'Under v3.1 alone' paragraph",
     r"Theorem~4(c) delivers the v3.2 data, which do not depend on (F). The v3.1 values are recorded in \S\ref{ss:LOv31}.",
     r"Theorem~4(c) delivers the v3.2 data, which do not depend on (F). The v3.1 values are recorded in \S\ref{ss:LOv31}. "
     r"In the v3.3 line (F) is no longer assumed, and the hand-off list has further items certified by step~4 "
     r"(\S\ref{ss:v33status}).",
     "step4_handoff_v44.tex"),
    ("S9 status item (c): strip statement",
     r"The strip statement is part of Hypothesis~(F), so H6 in step~3, and step~4's certified closure margins for it "
     r"(Table~\ref{tab:v32}), are conditional on (F) (\S\ref{ss:v32status}).",
     r"In the v3.2 line the strip statement is part of Hypothesis~(F), so H6 in step~3, and step~4's certified closure "
     r"margins for it (Table~\ref{tab:v32}), are conditional on (F) (\S\ref{ss:v32status}). In the v3.3 line the strip "
     r"statement is Lemma~U, and step~3 recomputes H6 (\S\ref{ss:v33status}).",
     "main_v3_v44.tex"),
    ("S10 'Consequently the main theorem holds once ...'",
     r"Consequently the main theorem (\S\ref{s:main}) holds once the far-family H1 item is closed and Lemma~R is proved.",
     r"The two items on which version~1 made the main theorem depend, the far-family H1 item and Lemma~R, are closed "
     r"(O1, O2). The main theorem's present status is in \S\ref{ss:v33status}: step~3 is not closed.",
     "O1, O2 of this document; CLAUDE.md status"),
    ("S11 hand-off table: (H6) row",
     r"(H6) & smooth, closed, immersed & Theorem~4(b) & \tA\\\bottomrule",
     r"(H6) & smooth, closed, immersed, symmetric about the $x$-axis & Theorem C-1(iv) of \fn{C1_identification_v44.tex} "
     r"(hand proof on compiled Lean theorems, \S\ref{ss:H6reg}); $\sup|\tilde\kappa e|\le0.0017429$, stage 22 of "
     r"\fn{FINAL_manifest_v38.json} \tC & \tA\\\bottomrule",
     "C1_identification_v44.tex Thm C1:thm (iv), proof tagged [hand]; out_v38/iface/C1_ke.json sup_ke 0.0017428478"),
    ("S12 margins table: v3.2 step-3 side header",
     r"\multicolumn{5}{@{}l}{\emph{Step-3 side of v3.2 (for reference; \fn{t9_h6_v32.json}); the closure rows are conditional on Hypothesis~(F)}}\\",
     r"\multicolumn{5}{@{}l}{\emph{Step-3 side of v3.2 (history; \fn{t9_h6_v32.json}); closure rows conditional on (F); see \S\ref{ss:v33status}}}\\",
     "history"),
    ("S13 rotation-mode norm (follow-up O6(c))",
     r"The rotation-mode norm still uses generous hand bounds for the frame-mode derivatives on $[10.25,60]$ (\fn{h1.py}, "
     r"block ``E\_I part''); they enter only through $\operatorname{rad}C$, a contribution $\le3.1\cdot10^{-11}$ to "
     r"$\sqrt{\Lam(t_k)}$, against the threshold $5\cdot10^{-8}$.",
     r"The rotation-mode norm enters only through $\operatorname{rad}C$, a contribution $\le3.1\cdot10^{-11}$ to "
     r"$\sqrt{\Lam(t_k)}$, against the threshold $5\cdot10^{-8}$. Since the follow-up O6(c) of Nov~12 "
     r"(\fn{STEP4_RESULT.md}; stages ``H1 exact points'' and ``H1 slab summary'' of \fn{FINAL_manifest.json}), the $E_I$ "
     r"part of the mode norms on $[10.25,60]$ is summed by cells, with the $S_1$ rows for $\sigma\le13$, the step-1 tail "
     r"envelopes beyond and C4 jets for $\hat\psi_3+2\alpha\hat h_2$, in place of the typed hand bounds of the first "
     r"version (one of them, $\abs{sc'}\le10^{-3}$, was false near $x=10.25$, sup $1.204\cdot10^{-3}$, though its "
     r"integral still over-covered). The maximal rotation-mode norm is $17.626$ (it was $17.642$), and the quoted bounds "
     r"stand. The allowance $1.01$ on the $C$-term is unchanged.",
     "STEP4_RESULT.md O6(c) follow-up; out/h1slab/h1_p*.log max rot_mode_norm 17.62543, C_rad 1.72069e-12"),
    ("S14 (H6) regularity: literal \\n (fatal under -halt-on-error)",
     "This item is closed.\\n\\end{remark}",
     "This item is closed.\n\\end{remark}",
     "LATEX_v2.md"),
    ("S15 Lemma O hypotheses: H6 'conditional on Lemma R'",
     r"v3.2 proves H6 in step~3, conditional on Lemma~R (\S\ref{s:H6});",
     r"v3.2 proves H6 in step~3 using Lemma~R, which is proved (\S\ref{s:H6}; in the v3.3 line H6 rests on Lemma~U, "
     r"\S\ref{ss:v33status});",
     "O2 of this document"),
    ("S16 Lemma O table: H6 row",
     r"conditional on Hypothesis~(F) (Lemma~R is proved, \S\ref{ss:v32status}); not met under v3.1 (\S\ref{ss:H6v31})}\\",
     r"conditional on Hypothesis~(F) (Lemma~R is proved, \S\ref{ss:v32status}); not met under v3.1 (\S\ref{ss:H6v31}); "
     r"in v3.3 the strip input is Lemma~U and (F) is not assumed (\S\ref{ss:v33status})}\\",
     "main_v3_v44.tex"),
    ("S17 far H1: stage number",
     r"\fn{p2/out/t9_lemmaO_far.json}, stage 28 of the final tail, \fn{ok = true}). Data of v3.2",
     r"\fn{p2/out/t9_lemmaO_far.json}, stage 30 of the final tail, counted from 0, \fn{ok = true}). Data of v3.2",
     "FINAL_manifest.json stage index 30 (28 in the 31-stage cleanup manifest)"),
    ("S18 axis crossing: stage number",
     r"\fn{p2/out/t9_axis_crossing.json}, stage 29 of the final tail)",
     r"\fn{p2/out/t9_axis_crossing.json}, stage 31 of the final tail)",
     "FINAL_manifest.json stage index 31"),
    ("S19 v32status title + note",
     r"\subsection{The step-3 dependency: Lemma R proved, Hypothesis (F)}\label{ss:v32status}",
     r"\subsection{The step-3 dependency in the v3.2 line: Lemma R proved, Hypothesis (F) (history)}\label{ss:v32status}" "\n"
     r"\emph{This subsection records the status of version~2 (Nov~12). In the v3.3 line (F) is no longer assumed; see "
     r"\S\ref{ss:v33status}.}" "\n",
     "main_v3_v44.tex"),
    ("S20 v32status: Consequence",
     r"\emph{Consequence.} The main theorem (\S\ref{s:main}) is stated \textbf{conditional on Hypothesis~(F)} "
     r"(\fn{writeup/main/main.pdf}). (This is not the item (F) of the bootstrap hypotheses of Theorem~P2, \S\ref{s:P2}.) "
     r"A fix, a uniformly local Lemma~S(B), is being evaluated (go/no-go tests); until it is in, (F) stays a hypothesis.",
     r"\emph{Consequence (version~2).} The main theorem of version~2 (\S\ref{ss:mainv32}) is stated \textbf{conditional "
     r"on Hypothesis~(F)} (\fn{writeup/main/main.pdf}). (This is not the item (F) of the bootstrap hypotheses of "
     r"Theorem~P2, \S\ref{s:P2}.) A fix, a uniformly local Lemma~S(B), was then being evaluated. It has since been "
     r"carried out as Lemma~U of the v3.3 line, where (F) is no longer assumed (\S\ref{ss:v33status}).",
     "main_v3_v44.tex Remark rem:hist"),
    ("S21 v32status: other caveats",
     r"\emph{Other step-3 caveats} (\fn{step3_proof_v32.tex} \S6, final-rigor-pass items, not step-4 items):",
     r"\emph{Other step-3 caveats of v3.2} (\fn{step3_proof_v32.tex} \S6, final-rigor-pass items, not step-4 items; for "
     r"the v3.3 line see \S\ref{ss:v33status}):",
     "history"),
    ("S22 main theorem section: new status subsection",
     r"\section{The main theorem}\label{s:main}" "\n" r"We combine the four steps.",
     r"\section{The main theorem}\label{s:main}" "\n" + V33STATUS + r"We combine the four steps.",
     "main_v3_v44.tex (Main Theorem, s:open, s:step4); step3_v44_main.tex (What remains); step4_handoff_v44.tex; "
     "manifest_v38 results; lean build log; task list of author items"),
    ("S23 main theorem name",
     r"\begin{named}{Main Theorem (conditional on Hypothesis~(F))}",
     r"\begin{named}{Main Theorem of version 2 (v3.2 line; conditional on Hypothesis~(F); history)}",
     "history"),
    ("S24 open item O2",
     r"\item \textbf{Step-3 dependency: Lemma~R proved; Hypothesis~(F) open} (\S\ref{ss:v32status}).",
     r"\item \textbf{Step-3 dependency} (updated Nov~24). \emph{Now:} Hypothesis~(F) is no longer assumed; in the v3.3 "
     r"line it is replaced by Lemma~U and the repaired Lemma~S(A), but \textbf{step~3 is not closed}; its open items, and "
     r"what step~4 delivers to it, are in \S\ref{ss:v33status}. \emph{History (version~2, \S\ref{ss:v32status}): Lemma~R "
     r"proved; Hypothesis~(F) open.}",
     "main_v3_v44.tex"),
    ("S25 O2: main theorem conditional (past)",
     r"The main theorem is therefore conditional on Hypothesis~(F) (\fn{writeup/main/main.pdf}).",
     r"The main theorem of version~2 was therefore conditional on Hypothesis~(F) (\fn{writeup/main/main.pdf}).",
     "history"),
    ("S26 O2: fix being evaluated",
     r"A fix, a uniformly local Lemma~S(B), is being evaluated (go/no-go tests). Under v3.1 alone, H6 is not met",
     r"A fix, a uniformly local Lemma~S(B), was then being evaluated; it became Lemma~U. Under v3.1 alone, H6 is not met",
     "main_v3_v44.tex"),
    ("S27 O1: stage number",
     r"is a pass/fail check arc by arc (stage 28 of the final tail, \fn{ok = true})",
     r"is a pass/fail check arc by arc (stage 30 of the final tail, counted from 0, \fn{ok = true})",
     "FINAL_manifest.json"),
    ("S28 O4 header + WSL step-3 kit run",
     r"\item \textbf{Kit B}: \fn{kit_step4b_v4.zip}",
     r"\item \textbf{Kits} (kit~B and the step-3 kits; updated Nov~24): \fn{kit_step4b_v4.zip}",
     "task"),
    ("S29 O4: WSL runs",
     r"test run in WSL \fn{kit4bv4_run_wsl1}: \texttt{TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]}. ",
     r"test run in WSL \fn{s1spec/kit4bv4_run_wsl1}: \texttt{TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]}; step-3 kit "
     r"\fn{kit_step3_v4.zip} (sha256 \fn{0f899eb2...}), test run in WSL \fn{s1spec/kit3v4_run_wsl5}: 80 PASS, 0 FAIL, "
     r"1 WARN (float noise). The WSL runs used the installation that holds the proof tree; they test the kits and are "
     r"not independent reproductions. ",
     "kit4bv4_run_wsl1/RUN_NOTES.md; kit3v4_run_wsl5/RUN_NOTES.md; step4_handoff_v44.tex sec. Kits"),
    ("S30 O4: kit4b v5 run dir",
     r"Windows run in progress, no result yet \OPEN{kit\_step4b v5 on Windows}.",
     r"Windows run in progress (run copy \fn{kitruns/kit4bv5_win1} on the authors' Windows machine, started Nov~24), "
     r"no result yet \OPEN{kit\_step4b v5 on Windows}.",
     "directory listing of the Windows kitruns folder"),
    ("S31 O6 header",
     r"\item \textbf{Script items found while writing: closed, except (c)}",
     r"\item \textbf{Script items found while writing: closed}",
     "STEP4_RESULT.md O6(c) follow-up"),
    ("S32 O6 (d): stage numbers",
     r"and it and the far-H1 stage are stages 28 and 29 of \fn{FINAL_manifest.json}, which now has 31 stages (\S\ref{app:repro}).",
     r"and the far-H1 stage and it are stages 30 and 31 of \fn{FINAL_manifest.json}, counted from 0 (28 and 29 in the "
     r"31-stage manifest of the cleanup; the final manifest has 34 stages, \S\ref{app:repro}).",
     "FINAL_manifest.json; FINAL_manifest_cleanup_pass2.json"),
    ("S33 O6 (c): closed",
     r"(c) \emph{Remaining, not covered by the cleanup:} the rotation-mode norm in \fn{h1.py} uses generous hand bounds, "
     r"a contribution $\le3.1\cdot10^{-11}$ to $\sqrt\Lam$ (\S\ref{ss:alpha}); it changes no conclusion.",
     r"(c) \emph{Closed by the follow-up O6(c) of Nov~12, which version~2 did not record:} the rotation-mode norm in "
     r"\fn{h1.py} is summed by cells instead of hand bounds (\S\ref{ss:alpha}); its contribution to $\sqrt\Lam$ stays "
     r"$\le3.1\cdot10^{-11}$; the 16 H1 points and the H1 summary became stages of the final tail.",
     "STEP4_RESULT.md O6(c) follow-up; FINAL_manifest.json"),
    ("S34 AI statement facts",
     r"independent Windows runs of step 3 and of kit A;",
     r"independent Windows runs of step 3 and of kit A, and the Windows runs of Nov~24 (step-3 kits v4 and v5, kit~B v4; O4);",
     "RUN_NOTES.md of the Windows runs"),
    ("S35 typed allowances",
     r"Typed allowances remain only where stated: \fn{h1.py}'s rotation-mode hand bounds (\S\ref{ss:alpha}) and",
     r"Typed allowances remain only where stated: the allowance $1.01$ on \fn{h1.py}'s $C$-term (\S\ref{ss:alpha}; its "
     r"rotation-mode hand bounds were replaced by cell sums in the follow-up O6(c)) and",
     "STEP4_RESULT.md O6(c) follow-up"),
    ("S36 manifest appendix: scripts/upstream",
     r"plus the SHA-256 of the scripts and of 117 upstream files.",
     r"plus the SHA-256 of the 26 scripts and of 97 upstream entries; at the end it re-checks all 523 pinned paths "
     r"(\fn{checks.hashes_verified_at_end}). (Version~2 described the 31-stage manifest of the cleanup, with 117 upstream "
     r"files; the final closed manifest of Nov~12 has 34 stages.)",
     "FINAL_manifest.json checks.hashes_verified_at_end; FINAL_manifest_cleanup_pass2.json"),
    ("S37 manifest appendix: stage list",
     r"\item \emph{The 31 stages}: \fn{cconst} (179 slabs), \fn{merge}, \fn{switch}, \fn{ycomp}, \fn{clam};",
     r"\item \emph{The 34 stages} (counted from 0): \fn{cconst} (179 slabs), \fn{merge}, \fn{switch}, \fn{ycomp}, the 16 "
     r"H1 exact points (\fn{t9h/h1_points.py}) and the H1 slab summary (\fn{t9h/h1_summary.py}), \fn{clam};",
     "FINAL_manifest.json stages"),
    ("S38 manifest appendix: stage list end",
     r"the far-family Lemma~O H1 (\fn{t9_lemmaO_far.py}, stage 28), the axis crossing (\fn{t9_axis_crossing.py}, stage 29) and \fn{t1clock}.",
     r"the far-family Lemma~O H1 (\fn{t9_lemmaO_far.py}, stage 30), the axis crossing (\fn{t9_axis_crossing.py}, stage "
     r"31), \fn{t1clock}, and the B8 containment check (stage 33, last, since no stage reads it).",
     "FINAL_manifest.json stages"),
    ("S39 manifest appendix: upstream",
     r"\item \emph{Upstream} (hashes in the manifest; not rerun by the tail): the B8 containment check, the 16 Krawczyk "
     r"points, the 32 derivative balls, the 32 window balls and their summaries of the H1 stage (\fn{t9h}).",
     r"\item \emph{Upstream} (hashes in the manifest; not rerun by the tail): the 32 derivative balls and the 32 window "
     r"balls of the H1 stage (\fn{t9h}; 97 entries). (In version~2's 31-stage manifest the B8 check, the 16 Krawczyk "
     r"points and the H1 summary were upstream too.)",
     "FINAL_manifest.json upstream (97); STEP4_RESULT.md"),
    ("S40 manifest appendix: h1rot rerun and v33/v38 manifests",
     r"(\fn{checks.cleanup_second_pass})." "\n" r"\end{itemize}",
     r"(\fn{checks.cleanup_second_pass})." "\n" + V33_MANIFEST_BULLET + "\n" r"\end{itemize}",
     "FINAL_manifest.json checks.h1rot_second_pass; REPIN_REPORT.md; EXCEPTION_v32_B8_check_outA5.json"),
    ("S41 hash table caption",
     r"\caption{SHA-256 (first 16 hex digits) of the certified outputs used here, as recorded in \fn{FINAL_manifest.json} "
     r"and checked against the files on Nov~12.",
     r"\caption{SHA-256 (first 16 hex digits) of the certified outputs used here, as pinned in the closed "
     r"\fn{FINAL_manifest.json} (34 stages) and re-checked against the files on Nov~24, 2026 (version~3). For five entries "
     r"version~2 printed the values of the 31-stage cleanup manifest; they are corrected here. \fn{t2/gamma0_final.json} "
     r"is not pinned by the manifest (it is a kit~A output); \fn{check_outA5.json} differs on disk from its pin by the "
     r"documented exception (\S\ref{app:repro}).",
     "sha256 of the files; FINAL_manifest.json pins; FINAL_manifest_cleanup_pass2.json"),
    ("S42 hash: FINAL_manifest.json", r"\fn{cd37f42214d97561}", r"\fn{e01cb4da71ba8253}", "sha256 on disk"),
    ("S43 hash: h1slab/summary.json", r"\fn{5813699c98e56dcf}", r"\fn{657c751a03411e39}", "manifest pin = disk"),
    ("S44 hash: final_tail.py", r"\fn{534e4ce4e2bd7fc0}", r"\fn{54b2d6428add1903}", "manifest script_sha256 = disk"),
    ("S45 hash: t9_h1_clam_s180.json", r"\fn{0a759b5b786d38bf}", r"\fn{ac2670ee3a4e45cc}", "manifest pin = disk"),
    ("S46 hash: check_outA5.json (pin)", r"\fn{43a4ae8b56ed562f}", r"\fn{73843eff7ab481cb}",
     "manifest pin (disk 549619b793d447df by the B8 exception)"),
    ("S47 hash table: v33/v38 manifests",
     r" & & \fn{p2/out/t9_axis_crossing.json} & \fn{ce63ff1087e9347e}\\\bottomrule",
     r"\fn{p2/out_v33/FINAL_manifest_v33.json} & \fn{8bc47e83e370fec3} & \fn{p2/out/t9_axis_crossing.json} & "
     r"\fn{ce63ff1087e9347e}\\" "\n" r"\fn{p2/out_v38/FINAL_manifest_v38.json} & \fn{be909bbf42543ab0} & & \\\bottomrule",
     "sha256 on disk"),
    ("S48 step-3 kit paragraph",
     r"reports \textbf{91 PASS, 0 FAIL, 0 WARN}, with the verifier bit-identical to the reference.",
     r"reports \textbf{91 PASS, 0 FAIL, 0 WARN}, with the verifier bit-identical to the reference. This kit covers step~3 "
     r"as of v3.1; the kits of the v3.3 line are \fn{kit_step3_v4.zip} and \fn{kit_step3_v5.zip} (O4).",
     "O4"),
    ("S49 kit A v3",
     r"The plain-overlap rule was the wrong criterion for intermediate slabs; the transported comparison should replace it in the kit's checker.",
     r"The plain-overlap rule was the wrong criterion for intermediate slabs. Kit~A v3 (\fn{s1spec/kit_step4a_v3.zip}, "
     r"SHA-256 " + hx("5eecc2e1888bd4abe7632d104a614b4a7b4f4b6db8e23a83ba180b654b53db82") + r") drops it: its checker "
     r"compares $L$ and $X_0$ at a common time (WARN only) and keeps the other rules (every slab validated, $k_{\max}$ and "
     r"final-slab overlap, final $k_{\max}\ge28$); a run finished with the earlier kit is re-checked by copying the new "
     r"\fn{check.py} (kit README).",
     "kit_step4a_v3.zip README.md (v3 changes)"),
    ("S50 commands: closed tail + v38",
     r"python3 -I final_tail.py     # the whole final tail on the B1 data -> out/FINAL_manifest.json",
     r"python3 -I final_tail.py     # the whole final tail on the B1 data -> out/FINAL_manifest.json" "\n"
     r"#   (closed v3.2 tail; it rewrites B8's check_outA5.json: see the exception record first)" "\n"
     r"cd v38; nice -n 19 python3 -I -B final_tail_v38.py; cd ..   # re-pin v38 -> out_v38/FINAL_manifest_v38.json",
     "final_tail_v38.py usage line; EXCEPTION_v32_B8_check_outA5.json"),
    ("S51 commands: upstream comment",
     r"# upstream H1 stage (hashes in the manifest's upstream list)",
     r"# upstream H1 stage: derivative and window balls (hashes in the manifest's upstream list);" "\n"
     r"# the 16 h1.py points, the H1 summary and B8 are also stages of final_tail.py (follow-up O6(c))",
     "FINAL_manifest.json"),
    ("S52 notes on sources",
     r"while this document follows the later step-3 status (Lemma~R proved, main theorem conditional on Hypothesis~(F));",
     r"while version~2 followed the later step-3 status of Nov~12 (Lemma~R proved, main theorem conditional on "
     r"Hypothesis~(F)) and version~3 follows the v3.3 line (\S\ref{ss:v33status}); \fn{STEP4_RESULT.md} predates the "
     r"v3.3 stages and the re-pin v38;",
     "STEP4_RESULT.md (no mention of v33/v38)"),
]

text = base
recs = []
for eid, old, new, src in EDITS:
    nb = base.count(old)
    nt = text.count(old)
    assert nb == 1, f"{eid}: old string occurs {nb} times in base"
    assert nt == 1, f"{eid}: old string occurs {nt} times in working text"
    line = base[: base.index(old)].count("\n") + 1
    text = text.replace(old, new, 1)
    assert text.count(new) >= 1
    recs.append({"id": eid, "base_line_v2": line, "source": src,
                 "old_sha256": hashlib.sha256(old.encode()).hexdigest(),
                 "new_sha256": hashlib.sha256(new.encode()).hexdigest()})

# post-conditions
assert "\\n\\end{remark}" not in text
assert "being evaluated (go/no-go" not in text
assert "kit~B: pending" not in text
assert "holds once the far-family" not in text
assert "cd37f42214d97561" not in text

HEADER = (
    "% step4_proof_v3.tex -- step-4 write-up, version 3 (status sync; new file).  The closed\n"
    "% s1spec/writeup/step4_proof/ and s1spec/writeup/step4_proof_v2/ are NOT edited and stay byte-identical.\n"
    f"% Base: {BASE}\n%   sha256 {BASE_SHA}\n"
    "% Approval: Micah Warren approved this status-synced version 3 on 24 Nov 2026.\n"
    "% Built: 24 Nov 2026 by tools/apply_sync_v3.py (exact-match assertions); record apply_sync_v3.json;\n"
    "% every edit, its old/new text and its source: SYNC_REPORT_v3.md.  Step 3 is NOT closed.\n"
    "% (The v2 header below is kept.)\n"
)
text = HEADER + text
os.makedirs(OUTD, exist_ok=True)
with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)
assert sha(BASE) == BASE_SHA
rec = {"args": {"argv": sys.argv, "env": {"python": sys.version},
                "inputs": {BASE: BASE_SHA, **{p: sha(p) for p in SOURCES}}},
       "edits": recs, "output": {OUT: sha(OUT)}}
with open(REC, "w") as f:
    json.dump(rec, f, indent=1)
print("OK", len(recs), "edits;", OUT, sha(OUT))
