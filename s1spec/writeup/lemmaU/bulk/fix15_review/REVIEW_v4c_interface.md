# REVIEW_v4c_interface: independent review of verifier fin13v4c, the (A) chain v37 with the flux reassembly, the interface lemmas C-1 and (GS), and the coherence of the step-3 proof set (fix round 15)

Nov 20, 2026. Fable, independent reviewer. The proof tree was read-only; everything written is in this directory (`fix15_review/`).

- Paths: `V` = `/home/claude/s1spec/step3nl_v4c`; `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`; `N` = `B/fix13_v37chain`; `I` = `B/fix14_interface`.
- Labels: **[cert]** arb with outward rounding on certified inputs taken as given (mine unless a file is named); **[sym]** exact sympy; **[num]** float check; **[hand]** my derivation, not machine-checked; **[read]** verified by reading the cited file:line; **[est]** estimate.
- Parts B and C were done by forks of this reviewer (same model, same context). The review was interrupted once by the usage limit. Part B's fork had finished its file; I reran its two scripts (same results) and spot-checked its diff, protected-tree and wiring claims. Part C was relaunched from the partial files and I spot-checked its main finding (C1-1) against the two cited texts. Parts A and D are mine.
- Runs: one process at a time (shared lock), nice 19, PYTHONHASHSEED=0, `python3 -B`; peak RSS 173 MB.
- **Not reviewed:** the inherited mathematics of v35/v36 beyond the v37 diff; the full text of `fix12_texts/` (only the parts needed for Part D); Lean was not built and DKS 2002 was not re-read; the chains were not rerun.

## 0. Verdicts

| Part | Verdict | One line |
|---|---|---|
| **A** verifier fin13v4c | **sound modulo items** | D-3, D-2a, A-1, A-2 are implemented correctly; H8 is a valid bound (coded coefficients are 2.8 to 17 times above my recomputation [cert]). The records issue was repaired adequately for the production run `runs/P`; one statement in the report about it is wrong (A-1r). |
| **B** (A) chain v37 + flux | **sound modulo items** | Wiring to `V/runs/P` is correct, regressions and double runs hold, no margin decreased; only records and text items. |
| **C** C-1 identification | **sound modulo items** | The chain Lean flow → step-4 gauge → P1 → P2 → t₁ is correct and not circular. One piece of text is owed: the a-priori form (P2-a) of Theorem P2 is a reading, not a written lemma (C1-1). |
| **C** (GS) | **sound modulo items** | τ₁ certificate, strict improvement, open/closed and the endpoint at 60 are correct. Minor records and wording items. |
| **D** coherence | **sound modulo items** | No circularity in the bootstrap order. Every hypothesis is either a bootstrap member that is strictly improved, or a listed input. The texts do not yet state this in one place and still carry fin13v4b / v36 / X_S = 60¼ (D-1 to D-5). |

**Step 3 is not closed.** No new blocking mathematical item was found. What remains is text (one lemma to write out, C1-1; the v38 text set), records, the re-pin and kits, the author items, and the final review.

## 1. Part A: verifier fin13v4c

### 1.1 What I did
- Diffed `V/verifier/{constants,rhs_iv,closure}.py` against their `_pre_v4c` backups; the backups are byte-identical to `step3nl_v4b/verifier/`. Three hunks only: read of `strip_cubic_v4c.json` with key-set and box asserts (`constants.py:145-159`), term H8 (`rhs_iv.py:273-282`), the D-2a asserts and args (`closure.py:96-97, 306-307`). Each is gated by `SC_V4C`.
- Read `v4bcert/strip_cubic_v4c.py`, `e0_A2_budget_v4c.py`, `tc_weighted.py`, the diff `E0_LEMMA_v4b.tex` → `E0_LEMMA_v4c.tex`.
- Recomputed the strip cubic coefficients with my own code: `a/a_strip_cubic.py` → `a/a_strip_cubic.json`, 0 fails.
- Checked every sha256 pin of the six rate tables against the disk, compared the three cert passes file by file, and compared the rate tables (production, double run, regression) with my own stripper.
- `step3nl_v4b` (652 files), `fix6_v36` (1039), `fix8_flux2_impl` (243) verify against the snapshots in `V/audit_v4c/`.

### 1.2 D-3, term H8: correct
- **Derivation [hand].** `tc_weighted.pkl` holds C_l = coefficients of T_c^X/ρ, so the strip form is ∫ρ·R₃^X/ρ. R₃^X = ∫₀¹(1−t) Σ ∂²C_l/∂U_i∂U_j (u₀+tw) w^(i)w^(j)w^(l) dt; the ½, the factor 2 for i ≠ j, the max over cells per (l,i,j) and the sum over (l,i,j) per monomial are right (`strip_cubic_v4c.py:70-81`).
- **Hölder assembly [hand].** On the strip w′ = v. ∫ρ|v|² ≤ RW·V0², ∫ρ|v||v′| ≤ RW·V0V1, ∫ρ|v′|² ≤ RW·V1², ∫ρ|v′| ≤ IRW·V1, ∫ρ ≤ IR, with the sup taken on the lowest-order factor (sI1r = raw E_I sup of |w′|, S0s = sup|w| on [16.25, 34.25]). `rhs_iv.py:278-280` is this term by term. Every term has x-degree 3; `tail_cert` passes with an unchanged degree table [read].
- **Hypotheses.** |α| ≤ 2.5e-4, |w|,|w′|,|w″| ≤ 0.02 on the strip: asserted (E0_alpha, U3_strip_S0..S2; worst ratio 0.773) and `constants.py:155` asserts that the JSON box matches.
- **Recomputation [cert]** (300 cells per unit σ instead of 200, x-cells as arb hulls, χ′ and χ″ in factored form):

  | monomial | coded | mine | coded/mine |
  |---|---|---|---|
  | (0,1,2) | 0.567028 | 0.204483 | 2.77 |
  | (0,2,1) | 2.571937 | 0.496073 | 5.18 |
  | (0,3,0) | 0.765674 | 0.136340 | 5.62 |
  | (1,0,2) | 0.283185 | 0.050979 | 5.56 |
  | (1,1,1) | 3.541830 | 0.423531 | 8.36 |
  | (1,2,0) | 1.067966 | 0.124175 | 8.60 |
  | (2,0,1) | 2.446e-17 | 2.418e-18 | 10.1 |
  | (2,1,0) | 5.587e-18 | 3.207e-19 | 17.4 |
  | (3,0,0) | 4.02e-37 | 4.94e-38 | 8.14 |

  All coded values cover mine. The coded values are loose because the stage evaluates the cutoff polynomial in expanded form on intervals; they are valid upper bounds. Pointwise float test of the monomial bound on 400 random states in the box: worst ratio 0.044 [num].
- No other use of `ADD_STRIP_CUBIC` remains on the V4C=1 path (`rhs_iv.py:282` is the regression branch only) [read].

### 1.3 D-2a, A-1, A-2
- **D-2a:** `closure.py:96-97` asserts far_w ≤ η0, far_w1 ≤ η1, α ≤ α_max with the three values read from `cubic_E0fine_0.001_0.01_0.01.json`. far_w1 = max(E0 bound, E_I bound) is conservative. Worst ratios 0.059 / 0.116 / 0.16 [read, `runs/P`]. Asserted at every march step and at the tail anchor.
- **A-1:** `e0_A2_budget_v4c.py` adds (d_ab + d_lead)(‖(g_A,g_B)‖ + L_C) = 0.179382 per αA² and the third-order loop piece 0.002641. Need 38.7232, margin **×1.65695** [cert, their JSON]; my fix-11 value was ×1.65695 (`fix11_review/a_recompute.json`). The stage asserts the absorption and is in the chain; the verifier constants are unchanged, as intended. Item (iv) of `E0_LEMMA_v4c.tex` carries both pieces.
- **A-2:** the whole-tail α-condition is reported as a ratio, ×1554.46 / ×1512.16 [read].
- A-3, A-4 (text) are done in `E0_LEMMA_v4c.tex:43, 112-115`.

### 1.4 The records issue (out_A overwritten, cert chain rerun)
What happened, from timestamps and pins [read]:
1. Cert pass A ran 09:21–09:24:21 into `v4bcert/out`.
2. `verifier/{constants,rhs_iv,closure}.py` were last written at 09:24:52, after pass A and before any verifier run.
3. Regression and production pass A (`runs/reg_fin13v4b`, `runs/P_A`) consumed pass A's cert outputs.
4. The driver found a stale `out_A` (copied from v4b), so it deleted pass A's outputs instead of keeping them; pass B was written to `out` and production `runs/P` consumed it.
5. The stale copy was renamed `out_A_from_v4b`; a third cert pass was written to `out_A`.

Assessment:
- **Production `runs/P` is clean [cert/read].** All 12 pins in each of its two rate tables match the disk. `out` and the third pass `out_A` are byte-identical except `tc_weighted.json` (runtime), `e0_identities.json` (its pin of that file), two logs and `chain_status.txt`. So the consumed cert constants are reproduced by an independent pass.
- **P = P_A and regression = fin13v4b** hold for every number (my comparison, args and timing ignored).
- **A-1r: one statement is wrong.** V4C_REPORT (double-run row) and V37_REPORT say P_A and P "pin the same sha256 of every cert input". They do not: `runs/P_A` and `runs/reg_fin13v4b` pin `e2v4b_constants.json` as 37386e13…, which matches no file on disk (now 7d6ada65…). Cause: `e2_consts.py:24` pins the text of `verifier/constants.py`, which was edited after pass A (step 2 above); the difference propagates through the pin of `e2_consts.json`. The numeric content was the same (the rate tables are equal). V4C_LOG step 4 ("renamed before any run") is true for verifier runs only.
- The third cert pass has no driver record (`runs/certC.out` is empty; the stages write to a hard-coded `out/`, so directories must have been swapped by hand).
- Verdict on the handling: adequate for the object to be pinned (`runs/P` with `v4bcert/out`). `runs/P_A` and `runs/reg_fin13v4b` each carry one unverifiable pin and should not go into a manifest as they are.

### 1.5 Findings, part A
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| A-1r | `V/V4C_REPORT.md` (double-run row, record note); `N/V37_REPORT.md` (checks table); `V/runs/P_A/*`, `V/runs/reg_fin13v4b/*` args | minor (records) | P_A and the regression pin an `e2v4b_constants.json` (37386e13…) that no longer exists; the "same sha256" statement is false for that file. Numbers equal. | Correct the two reports. Either rerun P_A and the regression on the current `out` (about 12 min each) or record an exception note; pin only `runs/P`. |
| A-2r | `V/runs/certC.out`, `V/run_prod_v4c.sh` | minor (records) | Third cert pass has no script or log of how `out_A` was produced. | One line in V4C_LOG, or make the pass a driver stage with `V4BOUT`. |
| A-3r | `V/SHA256SUMS_v4c` | minor (records) | `V4C_LOG.md` changed after the checksum file; the file also lists itself. | Regenerate. |
| A-4 | `V/v4bcert/strip_cubic_v4c.py:64-65` | minor rigor (negligible) | x-cells are built from float midpoint and radius; adjacent cells can leave a gap of about one ulp. My recomputation with arb hulls is far below the coded values, so no effect. | Use arb hulls at the next touch. |
| A-5 | `V/v4bcert/strip_cubic_v4c.py:52-55` | info | Expanded-polynomial cutoff jets make the C3S constants 2.8 to 17 times larger than needed. | Optional. |
| A-6 | `V/e0_lemma/E0_LEMMA_v4c.tex:115` | text | R-3 ratios quoted are the fin13v4b values (labelled so). | Use ×1554.46 / ×1512.16. |
| A-7 | inherited | minor ([read] input) | Monotonicity of W on [31.25, 34.25] behind RW, IRW (fix-11 A-5) is still not re-derived; it now also carries H8. | Keep in the [read] list or certify. |

## 2. Part B: (A) chain v37 and flux reassembly (fork; `b/PART_B.md`, `b/b_check.py` 29 checks 0 fail, `b/b_flux.py` 11 checks 1 fail = B-v37-1; rerun by me with equal results)

**Sound modulo items.** No numerical or wiring defect.

- **Diff [read].** 27 changed scripts; every `*_pre_v37` backup is byte-identical to its fix6_v36 or fix8 original. Hunks: path repoints, the `fin13v4c` branch of `v36cfg.py:31-37`, `KJ_PAD` for fin13v4c (`v36cfg.py:49`), `F2_VTAG` in `flux/code/f2common.py`, tag-parametrised names, args. No mathematical line changed. I confirmed the `v36cfg.py` diff myself.
- **Stale outputs.** `run_chain_v37.sh:21-30` deletes the copied outputs before each run. In-place outputs equal `runs/prod2` byte for byte; no prod2 JSON mentions step3nl_v4b, fix6_v36 or fin13v4b. Two files equal v36's: `ramp_U.json` (the known v35 byte copy, B3) and `weights_exact.json` (verifier-independent). The out_A-type problem did not occur here.
- **Wiring [cert].** Consumed J0, J1, J2, J4, slopes and peak ≥ raw `V/runs/P` values (excess 2.6e-7 to 2.2e-5 relative); ρ = 2.679e-3 ≤ raw; ā, β̄ equal `mod_sups`. The 1.5e-4 row dominates the 1.25e-4 row in all 12 consumed quantities. All 259 pins in 131 args blocks match the disk and point to `V/runs/P` (pass B).
- **Numbers.** κ_J certified 0.04242965327 reproduced, pad 1e-9 present [cert]. c_dU = 1.46342737776e-8, 2e-8/c_dU = ×1.366654 [cert]. [read, prod2]: ramp M1 1.184464; Ĝ22 1.916398; late/present M1 2.353209 / 2.77032; H6 1.860867; zone III 1.238888; C0_2 / C1_2 = 114.76328 / 0.06981372 (N2_ref) and 116.35923 / 0.26059988 (4·N2_ref). Minimum over 133 margin rows 1.18446; minimum ratio v37/v36 1.000000.
- **Regressions.** `cmp_runs` rerun: reg_v36 against `fix6_v36/runs/prod2` = 25 identical + 12 path-only; prod1 = prod2 for every JSON.
- **Flux.** `flux/out` = both snapshots (9 files, bytes); regression snapshot = fix8 apart from args/timing; FLUX0 = 1.64e-20, FLUX2 = 8.0e-20 unchanged. Recomputed [cert]: FLUX0/K̄0 ×4.80971; combined ×4.80453; C1 ×9.85e6; C6 ×141.32. [read]: C′ 2.65189e-26; C3 pad ×1547.6 / ×2236.4; C5 ×1554.46 / ×1512.16; C0–C7 and C9 pass; C8 unset.

| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| B-v37-1 | `N/flux/out/inputs.json`; `N/SHA256SUMS_v37` | minor (records) | `V37_LOG.md` was edited after the flux run and after the checksum file; its flux pin and checksum line no longer match. Only the run name is taken from it. | Regenerate the sums; pin a non-log file at the next flux touch. |
| B-v37-2 | `N/flux/out/flux2_cert.json` C8; `N/flux/texts/*` | text | C8 still "OPEN"; texts carry fix-8 wording. | Next text round. |
| B-v37-3 | `N/V37_REPORT.md` | text | C0_2 "114.763", "116.359", C1_2 "0.0698137" are rounded down. | Round up. |
| B-v37-4 | `N/run_chain_v37.sh:2-8` | text | Header comments describe fix6_v36 / fin13v4b. | Cosmetic. |

Fix-11 items in v37: B1, B2, B3, B4, B5, B7 and C-7 are **still open** (unchanged code); B6 is superseded by B-v37-3.

## 3. Part C: interface lemmas (fork; `c/PART_C.md`, reruns and `c/c_frame_shift.py`)

### 3.1 C-1 (`I/C1_identification.tex`): sound modulo items
- **Gauge and identities.** The P1 equation in step 4's gauge was re-derived by hand and equals `step4_proof.tex:227-230`. The 26 identities of `C1_sym.py` (rerun: 26/26 [sym]) are the ones the text uses. The parity and Parseval steps hold for infinitely many modes [hand].
- **Uniqueness.** The Grönwall lemma is correct [hand]: L² norm of the difference; it needs only sup‖v‖∞ and sup‖w‖_{L²}; the periodic integrations by parts are justified in H⁵. No gauge term arises because both solutions are in the same gauge. Step 4's Galerkin lemma proves exactly C⁰H⁵∩C¹L² [read].
- **No circularity past t₁.** The ∫k² bound is used only where the identification is already proved (P1), or inside a contradiction on T* (P2).
- **P2 as an a-priori theorem.** Step 4's proof is a-priori for a given smooth solution; its only existence input is the "[A] literature (P1)" sentence (`step4_proof.tex:370`), which Lean + DKS replace legitimately. But see C1-1.
- **`C1_ke.json`:** sup|κ̃e| = 1.74285e-3 [cert], slabs 2..180; reproduced.

### 3.2 (GS) (`I/GS_lemma.tex`): sound modulo items
- Tube lemma, τ₁ certificate (all of slab 180; 6.517e-3 / 0.05418 [cert], reproduced), open/closed argument, and the endpoint at x = 60 are correct. §6's objections to the v37 sketch are right.
- **GS_boot rerun on fin13v4c production** (`c/GS_boot_v4c.py`, input `V/runs/P/rate_table_fin13v4c_a0.00015.json`, the three D-2a asserts added to the required list) [cert]:
  - sup|ku| ≤ **5.3902e-4**, sup|u′| ≤ **0.040935**: unchanged from the quoted 5.40e-4 / 0.0410.
  - |w′| on [8, 10.25] ≤ 1.414421e-3 (was 1.414426e-3); sqrtL_max 5.308499e-8.
  - All box ratios ≤ 1 (largest 0.789); the new asserts ≤ 0.116.
  - Difference from the original `GS_boot.json`: 5th to 6th digits of ratios and sqrtL-derived values, and paths.

### 3.3 Findings, part C
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| **C1-1** | `I/C1_identification.tex:36-40, 111`; `step4_proof.tex:370` | **major (text owed; no numerical risk)** | (P2-a), the dichotomy "τ_e = 27.15 or t(τ) → T′", is a reading of step 4's short Step 7, not a proved lemma. Missing: limits of (λ, c, t) at τ_e; closedness and openness of the global normal-graph representation over the base (same pattern as GS-3/GS-4, with a ≥ ½ from `C1_ke`); the implicit function theorem at the limit. | Write it out, about 15 lines, in the next version of the C-1 text or in step 4's P2 proof. |
| C1-2 | `I/C1_sym.py:72-76` | minor | Y4's inputs (d_t s, L_t) are typed in, so the drift formula is [hand] + [sym]. The hand derivation confirms it. | Relabel. |
| C1-3 | `I/C1_identification.tex:105` | text | Lemma kP2 mixes x and y; conclusion right. | Edit. |
| C1-4 | `I/C1_identification.tex:42, 44` | inherited | Lean read, not built; n = 1 not elaborated; DKS not re-read. | Author items. |
| GS-1 | `I/GS_lemma.tex:57`, `I/GS_tau1.py:48` | minor (records) | The cited "true-curve shift ≤ 3.3e-11 (step4_proof)" is found nowhere. Bound from the stored matrix: ≤ 1.65e-11 [cert on 5-digit ball strings, `c/c_frame_shift.json`]; the 1e-9 widening covers it ×60. | Cite the computed bound. |
| GS-2 | `I/GS_boot.py:61, 66-70` | minor (negligible) | A 1e-6 sliver at x₁₃ uses the tail envelope slightly outside its range; k ≈ 2e-7 there. | Align the cell. |
| GS-3 | `I/GS_lemma.tex:71-83, 122` | minor | Lemma GS-2 also uses S0V = 0.0153 (the B2′ member), not "D1's boxes only". | State it; see D-2. |
| GS-4 | `I/GS_lemma.tex:71, 129` | minor (inherited) | Boxes beyond the tail anchor rest on the monotonicity statement of VERIFIER.md, read, not audited. `flux2_lemma_v37.tex:172` gives the argument for Λ ≤ L* + 2C′. | Audit once with the final text. |
| GS-5 | `I/GS_lemma.tex:71, 83`; `I/GS_boot.py:31-32` | text | Cites fin13v4b and `closure.py` l.143/162 (v4c: l.145/164). | Repoint; ship `GS_boot_v4c`. |
| GS-6 | `c/GS_boot_v4c.json` sqrtL_max | minor (negligible) | The script uses √L_max, while B1 is Λ ≤ L + 2C′ (relative 1e-11 on Λ). Margins are ×24 and more. | Multiply by (1 + 2C′/L) at the next touch. |

## 4. Part D: coherence of U-1, C-1, (GS), G1, flux lemma v2 and the E_I lemma

### 4.1 Dependency table (what each uses, what each proves)
| lemma | uses | proves |
|---|---|---|
| **C-1** | step 4: Prop Γ0, the Galerkin-limit lemma on the P1 slabs, P2 in the form (P2-a), `C1_ke.json`; Lean (L1), (L2) at n = 1; DKS (D1); U-1's gluing, maximal-flow, symmetry and blow-up lemmas (which need only a smooth initial immersion) | γ is a smooth symmetric immersed flow on [0, t(27.15)] ⊃ [0, t₁]; step 4's enclosures enclose γ; (S4.1), (S4.2); T_max = T* − t₁ |
| **(GS-1)** | C-1 (i), (iii); P2 boxes of slab 180; frame hull; S₁ chart | GS at τ₁ on \|x\| ≤ 60¼ (part of (S4.3)) |
| **(P2′)** in U-1 | (S4.1)–(S4.3); the closed set 𝔅 = GS + modulation box + B1; Neumann bound 0.0466 [read] | C¹ modulation while P ∈ 𝔅; τ_max = τ_mod |
| **Theorem U-1** | (P2′); P(τ) ∈ 𝔅 and (G) on [τ₁, τ′] | C¹ energies and identities at every τ, for the verifier, (A), Lemma U and the flux lemma |
| **E_I lemma** (D1) | U-1; boxes (B-a)–(B-f), (E0), (E2) (derived from B1 and asserted); S0V, S1V from B2′; **(I3′)** | Y = Λ − Z ≤ L strictly; H-J-W; ā, β̄ |
| **flux lemma v2** | U-1 (c); (A)'s data at β_A·Ê (B2′); BF = (BΦ) doubled exports, \|(a,b)\| ≤ 1e-6; OR boxes | (I3′); B1 in the form Λ ≤ L + 2C′; (F2)–(F6) strictly |
| **G1** (Lemma D, k = 2) | Lemma S(A) at j = 2; H-J-W on windows; U-1 (d) | exports C0₂, C1₂ without the factor 2: (BΦ) strictly |
| **(GS-2..4)** | B1 with its asserted boxes; S0V (B2′); ψ̂₃, ĥ₂, k_{S1} on cells; smooth flow; P(τ) continuous from (P2′)(b) | GS strictly (5.40e-4, 0.0410), open and closed |

### 4.2 Circularity check [hand]
- **C-1 is upstream of everything** and uses no step-3 statement. U-1's four flow lemmas are used by C-1 with F = Γ0; they do not use (S4.1).
- **The loop verifier → (A) → G1 → flux → verifier** is a simultaneous bootstrap, not a circle: the verifier consumes (I3′) and S0V at the hypothesis level (doubled exports, β_A·Ê, 0.0153); (A) consumes only the verifier's conclusions and "does not use (BΦ)" (`flux2_lemma_v37.tex:194`); each member is returned with a certified factor (×2 exports, ×2.12 S0V, pad and cushion for B1).
- **(P2′) ↔ GS:** (P2′) needs GS in 𝔅; GS-4 needs continuity of P(τ) from (P2′)(b). This is a joint continuation on [τ₁, τ_mod), because the strict improvement GS-2 does not use (P2′).
- **GS-2 uses hypothesis-level members only** (B1 boxes, S0V) plus GS at the same τ (needed for w to exist). The boxes are Sobolev consequences of Λ ≤ barrier, not of the energy inequality.
- **B1 with 2C′:** the box asserts cover Λ ≤ L + 2C′ on the march by the 1e-3 cushion (C2) and in the tail by `flux2_lemma_v37.tex:172`.
- **Regularity of the E_I identities.** All integrations by parts in `ei_lemma_v4.tex` hold at U-1's regularity (v‴ Lipschitz, v⁗ bounded with jumps on a finite set, v = v′ = v″ = v‴ = 0 at 10.25 and 60), with one exception in wording: see D-4.
- **Inputs that remain inputs** (proved by none of the six): Lean build and n = 1 instance; DKS (D1); (P3′) (used only for smoothness away from p*, not in the bootstrap); step 4's P1/P2 statements with their own labels and (P2-a) (C1-1); the graph equation, the Neumann bound and the C⁴ matching of ψ̂₃, ĥ₂ [read]; the chart lemma at k = 4 [hand, sympy structure check]; W monotone on the strip [read]; the tail monotonicity statement of VERIFIER.md [read].

### 4.3 Findings, part D
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| D-1 | `B/fix10_proofs/U1_regularity.tex:93, 99-101, 111`; `B/fix12_texts/continuation_v37.tex:36, 41` | medium (text owed; no mathematical gap) | U-1 and the continuation text pose (GS) on \|x\| ≤ 60¼ with the closed set as domain of the modulation, and Lemma C:lem:GS still has the invalid sketch. As written, U-1's hypothesis is in no bootstrap: GS-2 improves it on \|x\| ≤ 60 only. The fix-14 edits (X_S := 60, 𝒪 through GS°, replacement lemma) make the set fit; I checked that U-1's proofs need nothing on (60, 60¼] [hand]. | Issue U1 and continuation in new versions with the edits of `GS_lemma.tex` §6. |
| D-2 | `B/fix12_texts/lemmaU_v37.tex:388-395, 422-441`; `step3_v37_main.tex:68-90`; `GS_lemma.tex:122` | medium (text) | The bootstrap set exists in two versions. `lemmaU_v37` defines BS_U = {B1, B2′, B3U, B4, G} and proves closure for it; BF and GS are added only by reference in `step3_v37_main`. The order table there has "GS: B1 boxes; U-1", omitting S0V (B2′) and the joint continuation with (P2′); the closed/open step for (BΦ) is only a remark. | One statement of the extended set, one order table with the rows of §4.1, one closure proof that covers BF and GS. |
| D-3 | `B/fix12_texts/ei_lemma_v4.tex:1-4, 17, 33-44, 128, 147`; `flux2_lemma_v37.tex:24, 27`; `continuation_v37.tex:19, 65`; `lemmaU_v37.tex` | medium (text) | The text set describes fin13v4b / fix6_v36 / fix8. Owed for fin13v4c: H8 and the D-2a hypotheses in the E_I/E0 hypothesis list; delete "ADD_STRIP_CUBIC not audited"; A-1 budget ×1.657; v37 numbers; (S4.1) by C-1 instead of step 4's [A] sentence; U-1 gap items 1 and 4 closed. | Text set v38, generated from the v37 chain and `V/runs/P`. |
| D-4 | `B/fix12_texts/ei_lemma_v4.tex:114` (item (i), level 2) | minor (one sentence) | "(L₁v)″ = M₁v″ + e₁v′ + e₀v, so ⟨v″, (L₁v)″⟩_W ≤ …" uses six derivatives of v, while U-1 gives v⁗ bounded with jumps at B ∩ [10.25, 28]. The inequality ⟨D₂v, L₁v⟩_W ≤ −½V₄² − μ₂V₂² + E₁V₁V₂ + E₀V₀V₂ is still true: both sides are continuous on H⁴ with v, v′, v″, v‴ vanishing at the ends, and smooth functions are dense there [hand]. Fix-11 item C-6 (list each integration by parts with its regularity) is otherwise met in substance, not in form. | Add the density sentence and a short list of the integrations by parts with the regularity used. |
| D-5 | `B/fix10_proofs/U1_regularity.tex:111` | minor (text) | 𝔅 is stated with √Λ ≤ √Λ_bar, while B1 is now Λ ≤ L + 2C′. | State B1 in the flux-lemma form. |

## 5. Updated closing list for "step 3 closed"

Numbering follows §5 of `fix11_review/REVIEW_v4b_v36_proofs.md`.

**A. Flux lemma v2**
1. U-1: reviewed. ~~(1a) (GS) lemma with script~~ **done and reviewed (sound modulo GS-1..GS-6)**; ship `GS_boot_v4c`. Remaining: (1b) text: U1 in a new version with X_S = 60, 𝒪 through GS°, C-4, B1 form, (S4.1) by C-1 (D-1, D-5).
2. G1: reviewed. Remaining as before: text edits of `lemmaS_A` and `lemmaU` Lemma D; per-cell arb assert (C-7, still open in v37); PROOFS_REPORT count.
3. Flux text: retag to v37 numbers (B-v37-2, D-3).
4. Pin `N/flux/out/*` (not fix8) in a manifest and kit; fix B-v37-1 first.

**B. Verifier**
5. `ei_lemma_v4.tex`: written, read here for coherence only. Remaining: retag to fin13v4c with H8, D-2a, A-1 (D-3); the density sentence and the list of integrations by parts (D-4); a line-by-line review with the final text set (item 14).
6. v4b delta: closed by fin13v4c (A-1, A-3, A-4 done; A-2 done).
7. Audits: D-2a **done** (asserted). Remaining: D-1a pins; DS0..DS3 and the K_ρ/K_W Sobolev constants unaudited; W monotone on the strip [read] (A-7); tail monotonicity statement (GS-4).
- ~~N1 verifier pass with a certified strip cubic term~~ **done: fin13v4c, reviewed, sound modulo records items.** Remaining: A-1r (correct the reports; rerun or except P_A and the regression), A-2r, A-3r.

**C. (A) chain, Lemma U, zone III**
8. ~~Chain on the new tag~~ **done: v37, reviewed, sound modulo items.** Remaining records: B1, B2, B3, B4 (chain-side asserts), B-v37-1.
9. Texts: B5, B7, B-v37-3, B-v37-4; one bootstrap set and one closure proof including BF and GS (D-2).

**D. Step 4 and kits (Micah's go)**
10. Re-pin stages 8, 9, 18, 19, 20 on fin13v4c / v37 (`V/runs/P`, `V/v4bcert/out` including `strip_cubic_v4c.json` and `e0_A2_budget_v4c.json`, `N/runs/prod2`, `N/flux/out`), new final_tail driver reading the EXCEPTION record, new manifest. **The hold reason (N1) is gone**; it needs Micah's go. Add to the manifest: `I/C1_ke.json`, `I/GS_tau1.json`, `GS_boot_v4c.json`. Step-4 wording: (H6) by C-1, GS named in the hand-off list.
11. Kits rebuilt and run; Windows runs.

**E. Author items**
12. Lean: `lake build` and an n = 1 instantiation file (C-3 / C1-4). DKS pp. 1234–35 re-read, including the uniform bounds on ∂_s^m κ used in C-1. (P3′), Lemma O arc by arc, single-point statement: unchanged.
- ~~N2 step-4 to step-3 interface~~ **done and reviewed: C-1 sound modulo items.** Remaining: **C1-1, write out (P2-a)** (about 15 lines; the one piece of mathematics still owed in text); C1-2, C1-3; GS-1 citation.

**F. Whole document**
13. LaTeX compile: unchanged (now also `E0_LEMMA_v4c.tex`, `C1_identification.tex`, `GS_lemma.tex`).
14. Final independent review of the assembled step 3 on the v38 text set.

## 6. Files in this directory
- `REVIEW_v4c_interface.md`: this file.
- `a/a_strip_cubic.py` → `a/a_strip_cubic.json`: part A, independent strip cubic coefficients [cert] and pointwise test [num].
- `b/PART_B.md`, `b/b_check.py` → `b_check.json`, `b/b_flux.py` → `b_flux.json`, `b/cmp_reg_vs_v36prod2.txt`, rerun logs, `b/first_pass/` (outputs of the first run, equal).
- `c/PART_C.md`, `c/GS_boot_v4c.py` → `GS_boot_v4c.json`, `c/c_frame_shift.py` → `c_frame_shift.json`, and rerun copies of `C1_sym.py`, `C1_num.py`, `C1_ke.py`, `GS_tau1.py` with outputs.
- `SHA256SUMS`.

Proof-tree check after the runs: no file under `step3`, `step3nl`, `step3nl_v33`, `step3nl_v4`, `step3nl_v4b`, `step3nl_v4c`, `step4`, `s1_tight`, `fix6_v36`, `fix8_flux2_impl`, `fix10_proofs`, `fix11_review`, `fix12_texts`, `fix13_v37chain`, `fix14_interface` is newer than the start of this review; no `__pycache__` was written; `fix14_interface/SHA256SUMS` verifies.
