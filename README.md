# Hair clip: a computer-assisted blow-up proof for curve diffusion flow

Micah Warren (University of Oregon).

**Status: work in progress. The proof is NOT complete; see "Status" below.**
**License:** code under the MIT License (`LICENSE`); the mathematical texts (.tex, .pdf, .md write-ups) are all rights reserved by the authors until publication.

## What this is

We want a computer-assisted proof that an explicit, embedded, symmetric planar "hair clip"
curve Γ0 develops a finite-time singularity under curve diffusion flow (CDF: normal velocity
−κ_ss). The singularity is modelled on the non-embedded one-loop self-similar shrinker **S1**. The
evolving curve loses embeddedness early (t ≈ 0.003), as S1 itself self-crosses.

The proof has four steps:

| Step | Content | Main locations |
|---|---|---|
| 1 | Existence of S1 with tight validated enclosures | `s1_cap/`, `s1_tight/` (hull check `s1_tight/consistency/check_outA5.py`) |
| 2 | Linear stability of S1 modulo symmetries (independent result; step 3 does not use it) | `s1spec/writeup/step2_proof.tex/.pdf`, step-2 scripts at the top of `s1spec/` (`README_step2.md`), `s1spec/lin/` |
| 3 | Nonlinear trapping near S1 in the rescaled frame | verifiers `s1spec/step3nl*/` (final: **`s1spec/step3nl_v4c/`**, which also reads `step3nl/`, `step3nl_v4/`, `step3nl_v4b/`, `step3/`); lemma chain `s1spec/writeup/lemmaU/bulk/` (final (A) chain `fix13_v37chain/`, flux lemma `fix8_flux2_impl/`, proofs `fix10_proofs/`, interface `fix14_interface/`); **final text set `s1spec/writeup/lemmaU/bulk/fix27_texts_v44/`** (PDFs in `pdf/`; assembled step 3 is `step3_v44_main`) |
| 4 | Validated integration from Γ0 to the hand-off time t₁ | `s1spec/step4/` (authoritative summary `STEP4_RESULT.md`; manifests below); write-up `s1spec/writeup/step4_proof_v3/` |

Main-theorem drafts: `s1spec/writeup/main/main.tex/.pdf` (older; conditional on Hypothesis (F)) and
`s1spec/writeup/lemmaU/bulk/fix27_texts_v44/main_v3_v44.tex/.pdf` (current draft).

The reviews are kept next to
their subjects (`*REVIEW*.md`, `CHECK*.md`, `DELTA_CHECK*.md`, `AUDIT*.md`).

## Status (Nov 24, 2026)

- Steps 1, 2 and 4 are complete as certified computations with pinned manifests.
- **Step 3 is NOT closed.** The certified chain (verifier fin13v4c, (A) chain v37, flux lemma v2)
  and the text set v44 have passed several independent reviews with no blocking item, but these remain open:
  - author confirmations: (P3′); Lemma O checked arc by arc; the single-point statement;
    (P2′)/frames;
  - proofreading of the v44 PDFs, and a final text round (Lean commit hash and build log,
    confirmations, a bibliography entry for the Lean code);
  - the final independent review.
- The step-4 write-up patch (`fix27_texts_v44/STEP4_PATCH_NOTE_v44.md`) has not been applied.
- Lean: the CDF well-posedness, continuation and L² blow-up criterion parts are formalised in Lean
  (see `lean/`). They compiled in the original repository; **the vendored copy here has not been built.**

Every number in the documents is labelled *certified* (interval arithmetic, python-flint `arb`,
outward rounding) or *estimated*.

## Layout and paths

The directory structure mirrors the original working tree. **The code uses absolute paths under
`/home/claude/`**: scripts expect this repository's contents to be at `/home/claude/` (for example
`/home/claude/s1spec/step4/...`). To rerun scripts in place, put a copy or symlink of the tree there,
or use the self-contained kits below, which carry their own paths.

This repository is a *curated* copy (about 20,700 files). See `MANIFEST.tsv` (sha256, size, repo
path, source path, reason) for the full list. The rules were:
- included: all files pinned by the three step-4 manifests (`s1spec/step4/p2/out/FINAL_manifest.json`,
  `out_v33/FINAL_manifest_v33.json`, `out_v38/FINAL_manifest_v38.json`; 1055 files across steps 1–4),
  and all code, text, PDFs and small data in the step directories;
- excluded: large non-pinned intermediates (step-4 `.npz/.pkl` snapshots and JSON over 1 MB, other
  binaries over 1 MB, text over 5 MB), compiled executables, the kit trees and kit run directories
  (their `RUN_NOTES.md` and logs are kept), the third-party CAPD library, and a third-party arXiv source.
  Some exploratory scripts may therefore point to files that are not here.

Final pinned manifests (sha256):
- step 4, v3.1/v3.2 (34 stages): `s1spec/step4/p2/out/FINAL_manifest.json` e01cb4da…
- v3.3 (21 stages): `s1spec/step4/p2/out_v33/FINAL_manifest_v33.json` 8bc47e83…
- **v38 (29 stages, current; pins the step-3 v4c/v37 objects)**: `s1spec/step4/p2/out_v38/FINAL_manifest_v38.json` be909bbf…
  (driver `s1spec/step4/p2/v38/final_tail_v38.py`; documented exceptions `EXCEPTION_*.json` there and in `p2/v33/`).

## Reproducing (kits)

The kits are self-contained reproductions, needing Python 3 and python-flint. See the README in each zip.

| Kit | File | sha256 | Runtime | Recorded results |
|---|---|---|---|---|
| Step 3 | `s1spec/kit_step3_v5.zip` (16 MB) | f0ea2874…413e | ~35–40 min | WSL 80 PASS / 0 FAIL / 1 WARN; Windows 78 PASS / 0 FAIL / 3 WARN (float last bit) |
| Step 4 B | `s1spec/kit_step4b_v5.zip` (44 MB) | 49c4b722…8ea2 | ~5–6.5 h (`run_all.py --jobs 4..6`) | v4 kit: WSL and Windows 159 PASS / 0 FAIL / 1 WARN; v5 not yet run end to end on Windows |
| Step 4 A | `s1spec/kit_step4a_v3.zip` | 5eecc2e1…db82 | short | Γ0 bit-identical, 96/96 slabs (Windows) |

Usage: unzip into an empty directory, then `python run_all.py` (step 4 B: `--jobs N`), then `python check.py`.
Run notes: `s1spec/kit*_run*/RUN_NOTES.md`, `s1spec/kit4bv5_test_wsl1/`, `s1spec/kit_v5_build/KIT_V5_REPORT.md`.

## Lean

`lean/` contains the 16 files `CWFlow/Curve/*.lean`, the n = 1 CDF results used by step 3
(maximal flow and L² energy blow-up criterion, extension under bounded energy, uniqueness,
symmetry, shift/glue), together with the part of the authors' CWFlow library they import. See
`lean/README.md`, including the note about which library files this would make public.
Earlier drafts and logs: `s1spec/writeup/lean_n1/`, `s1spec/writeup/lean_bridge/`.

## Related work

G. Wheeler, "Pinchoff by surface diffusion", arXiv 2608.21882, a computer-assisted embedded torus
pinch-off for surfaces. It is complementary to this work, which treats planar curves, an explicit
initial curve and a non-embedded shrinker. The bibliography is in `s1spec/writeup/main/main.bib`.
