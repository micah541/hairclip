#!/usr/bin/env python3
"""apply_patch_v2.py: build step4_proof_v2.tex from the CLOSED step4_proof.tex and the
10 edits P1..P10 of STEP4_PATCH_NOTE_v44.md (approved by Micah Warren, 24 Nov 2026).

Reads the closed file only (never writes it).  Asserts the base and note sha256, asserts that
every old string occurs exactly once in the base (and once in the working text when applied),
and writes step4_proof_v2.tex + apply_patch_v2.json (args block, per-edit record).
Usage: python3 -I apply_patch_v2.py
"""
import hashlib, json, os, sys, datetime

BASE = "/home/claude/s1spec/writeup/step4_proof/step4_proof.tex"
NOTE = "/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/STEP4_PATCH_NOTE_v44.md"
OUTD = "/home/claude/s1spec/writeup/step4_proof_v2"
OUT = os.path.join(OUTD, "step4_proof_v2.tex")
REC = os.path.join(OUTD, "apply_patch_v2.json")
BASE_SHA = "6b4c183ab14a339b25472ff4e0a1c566d907c244059b35ad62299dfb2f15fb4f"
NOTE_SHA = "8bb6bfe3d343295b4e6da8f4657a955f240ad7dd444ce9dcd16af38067fa7ebc"
F27 = "writeup/lemmaU/bulk/fix27_texts_v44/"


def hx(h):
    """a 64-hex sha256 as four \\fn chunks joined by \\allowbreak (the full hash in one \\texttt overflows)"""
    assert len(h) == 64
    return r"\allowbreak".join(r"\fn{" + h[i:i + 16] + "}" for i in range(0, 64, 16))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


assert sha(BASE) == BASE_SHA, "base sha256 mismatch"
assert sha(NOTE) == NOTE_SHA, "note sha256 mismatch"
base = open(BASE, encoding="utf-8").read()

LEAN = (r"compiled Lean theorems \fn{cdf_exists_maximal}, \fn{cdf_symmetric} of \fn{CWFlow/Curve/}, "
        r"commit \fn{3130aef39829f02c126aeb166fbb9583f7670108} of the authors' local repository \fn{CWFlow} "
        r"(public location pending); the full \fn{lake build CWFlow.Curve.All} at that commit completed on "
        r"24~Nov~2026 (10316 jobs, no errors, no \fn{sorry}; log "
        r"\fn{writeup/lean_bridge/full_build_Curve_All_3130aef.log})")

# (id, old, new, mode) ; mode: 'note' = as the note, 'adapted' = see PATCH_REPORT_v2.md
EDITS = [
    ("P1",
     r"Theorem~4 delivers the hand-off conditions (H1)--(H6) of step~3 v3.1 \emph{as amended by v3.2}",
     r"Theorem~4 delivers the hand-off conditions (H1)--(H6) and (GS) of step~3 v3.1 \emph{as amended by v3.2}, "
     r"with the v3.3 items of \fn{" + F27 + r"step4_handoff_v44.tex} item~(1)",
     "adapted"),
    ("P2",
     r"hand-off (T9) & (H1)--(H4), (H6) for every $\tau^*$ in slab 180",
     r"hand-off (T9) & (H1)--(H4), (H6) by C-1, (GS) for every $\tau^*$ in slab 180",
     "note"),
    ("P3",
     r"\item $\gamma(t_1)$ is a smooth closed immersed curve.",
     r"\item $\gamma(t_1)$ is a smooth closed immersed curve, symmetric about the $x$-axis, the time-$t_1$ curve "
     r"of the maximal smooth CDF flow from $\Gamma_0$ (Theorem C-1 of \fn{" + F27 + r"C1_identification_v44.tex})."
     "\n"
     r"\item[(GS)] on $|x|\le60$ the tip arc of $\Gamma_1$ in the step-3 frame is a normal graph $u$ over $S_1$ "
     r"with $\sup|k u|\le0.0065170$, $\sup|u'|\le0.054182$ (stage 23 of "
     r"\fn{step4/p2/out_v38/FINAL_manifest_v38.json}) \tC.",
     "adapted"),
    ("P4",
     r"Continuation uses the literature fact (P1) of step~3 \tA.",
     r"The step-4 solution is the smooth CDF flow from $\Gamma_0$, and it continues past $t(27.15)$, by "
     r"Theorem C-1 of \fn{C1_identification_v44.tex} (" + LEAN + r"; data bridge by "
     r"\fn{exists_isLagImmersion_of_periodic}; Lemma (P2-a); \fn{C1_ke.json}, stage 22 of manifest\_v38).",
     "adapted"),
    ("P5",
     r"(H6) follows from (b) \tA.",
     r"(H6) is Theorem C-1(iv) of \fn{C1_identification_v44.tex}; (GS) is Lemma GS-1 of "
     r"\fn{" + F27 + r"GS_lemma_v44.tex} \tC.",
     "adapted"),
    ("P6",
     "continuation uses the literature fact (P1) of step~3 \\tA.\n",
     r"continuation: Theorem C-1 and Lemma (P2-a) of \fn{C1_identification_v44.tex} "
     r"(the a-priori form of Theorem P2)." "\n",
     "note"),
    ("P7",
     r"this is the continuation fact (P1) of step 3, from the literature \tA",
     r"by Lemma (P2-a) and the Proposition ``identification on P2'' (label \fn{C1:prop:P2}) of "
     r"\fn{C1_identification_v44.tex}",
     "adapted"),
    ("P8",
     r"it is smooth for $t>0$ by parabolic regularity, and closed and immersed as the CDF evolution of a closed immersed curve \tA",
     r"it is the smooth symmetric CDF flow from $\Gamma_0$ on $[0,t(27.15)]$, including $t=0$ "
     r"(Theorem C-1 (i)--(iv) of \fn{C1_identification_v44.tex}; no parabolic smoothing is used, since "
     r"$\Gamma_0$ is a trigonometric curve)",
     "note"),
    ("P9",
     r"\item \textbf{Kit B} (P2 and the hand-off on the user's machine): \textbf{still pending}. Until then the P2 "
     r"and T9 chain is reproduced only by the two passes of the final tail on the sandbox (\S\ref{app:repro}).",
     (r"\item \textbf{Kit B}: \fn{kit_step4b_v4.zip} (sha256 "
     + hx("ab47b5416e8a6d58d5267e2a93c154bf979809005de0abb4e992dd989df7f0fa") + r") test run in WSL "
     r"\fn{kit4bv4_run_wsl1}: \texttt{TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]}. "
     r"Windows runs on the authors' Windows~11 machine (24~Nov~2026, run copies): "
     r"\fn{s1spec/kit4bv4_run_win1/RUN_NOTES.md}, \texttt{TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]} after four "
     r"Windows-only fixes in the run copy; the step-3 kit \fn{kit_step3_v4.zip}, "
     r"\fn{s1spec/kit3v4_run_win1/RUN_NOTES.md}: reproduces (72 PASS; 6 FAIL, all CRLF-only differences of "
     r"text outputs; 3 WARN, last-bit float differences; two Windows-only fixes in the run copy). "
     r"Kit v5 carries these fixes: \fn{kit_step3_v5.zip} (sha256 "
     + hx("f0ea28743d340af8a440b4280e7ccc7ed71b5306635c577da0d568b619ee413e") + r"), Windows run "
     r"\fn{s1spec/kit3v5_run_win1}: 78 PASS, 0 FAIL, 3 WARN (last-bit float), no run-copy patch; "
     r"\fn{kit_step4b_v5.zip} (sha256 "
     + hx("49c4b7227069090d0b3573e6c6a10d42d77e94a437b3306101c3e303dbdf8ea2") + r"): Windows run in progress, "
     r"no result yet \OPEN{kit\_step4b v5 on Windows}."),
     "adapted"),
    ("P10",
     r"\paragraph{Kit B (P2 and the hand-off).} \textbf{Pending.} It will cover",
     r"\paragraph{Kit B (P2 and the hand-off).} See open item O4 (\S\ref{s:open}). It covers",
     "adapted"),
]

text = base
recs = []
for eid, old, new, mode in EDITS:
    nb = base.count(old)
    nt = text.count(old)
    assert nb == 1, f"{eid}: old string occurs {nb} times in base"
    assert nt == 1, f"{eid}: old string occurs {nt} times in working text"
    line = base[: base.index(old)].count("\n") + 1
    text = text.replace(old, new, 1)
    assert text.count(new) >= 1
    recs.append({"id": eid, "mode": mode, "base_line": line,
                 "old_sha256": hashlib.sha256(old.encode()).hexdigest(),
                 "new_sha256": hashlib.sha256(new.encode()).hexdigest()})

HEADER = (
    "% step4_proof_v2.tex -- step-4 write-up, version with the approved patch v2 (new file; the closed\n"
    "% s1spec/writeup/step4_proof/step4_proof.tex is NOT edited and stays byte-identical).\n"
    f"% Base: {BASE}\n%   sha256 {BASE_SHA}\n"
    f"% Patch note: {NOTE}\n%   sha256 {NOTE_SHA} (edits P1..P10)\n"
    "% Approval: Micah Warren approved the 10 edits on 24 Nov 2026.\n"
    "% Built: 24 Nov 2026 by tools/apply_patch_v2.py (exact-match assertions); record apply_patch_v2.json;\n"
    "% adaptations of the note's wording are listed in PATCH_REPORT_v2.md.  Step 3 is NOT closed.\n"
    "% (Unrelated: step4_proof/step4_proof_v2.tex in the closed dir is an older Nov 12 draft, not this file.)\n"
)
text = HEADER + text
os.makedirs(OUTD, exist_ok=True)
with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)
assert sha(BASE) == BASE_SHA
rec = {"args": {"argv": sys.argv, "env": {"python": sys.version},
                "inputs": {BASE: BASE_SHA, NOTE: NOTE_SHA}},
       "edits": recs, "output": {OUT: sha(OUT)}}
with open(REC, "w") as f:
    json.dump(rec, f, indent=1)
print("OK", len(recs), "edits;", OUT, sha(OUT))
