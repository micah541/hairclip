# FLUX2_IMPL_REPORT: flux lemma v2 implemented (fix round 8)

Nov 20, 2026. This implements the design approved by Micah on Nov 20: `B/fix7_flux2/FLUX2_DESIGN.md`, `flux2_lemma.tex` and `probe/`.
- Paths: `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix8_flux2_impl`; `B` = `.../lemmaU/bulk`.
- Log: `D/FLUX2_LOG.md`.
- Labels: [cert] = python-flint arb with outward rounding, constants taken from raw JSON and rounded up (lower bounds down); [sym] = exact sympy; [hand] = derivation by hand, not machine-checked; [est] = estimate; [read] = verified by reading code.

## Verdict
**Flux lemma v2 holds with the present FLUX0 = 1.64e-20 and FLUX2 = 8.0e-20, so no verifier rerun is needed. Every scripted condition C0–C7 passes at both α points of fin13v4b.**
- No stop rule fired:
  - C0_2 = 116.38 against the kill value ≈ 4e4;
  - the combined pointwise margin is ×4.80 (≥ ×2);
  - the χ₄ structure assertion passed;
  - FLUX0/FLUX2 are unchanged.
- **This does not close step 3.** Still open (§4): the regularity lemma U-1, the written k = 2 export proof (G1), the [hand] verifier-side readings, and an independent review.

## 1. Certified results (`out/flux2_cert.json`, verdict PASS)
**Inputs:** the (A) chain v36 `fix6_v36/runs/prod2`, which is the final production run and byte-identical to prod1. Setting 4·N2_ref, β_A = 2, exports doubled to the (BΦ) level.

| quantity | value | label |
|---|---|---|
| K̄0 (level 0, pointwise) = Σκ⁽⁰⁾N̄² + K0^nl + K0^w | **3.4099e-21** (lin 3.40954e-21, nl 3.09e-25, w 9.9e-70) | [cert] |
| FLUX0 / K̄0 | ×4.8096 | [cert] |
| K2^pt (level 2, pointwise part; λ1 = λ2 = 1) | 1.3900e-19; ε_I·K2^pt = 4.17e-24 | [cert] |
| **(FLUX0 + ε_I FLUX2)/(K̄0 + ε_I K2^pt)** | **×4.8044** | [cert] |
| integrated part: ∫_I φ2 ≤ C″ + S″T | C″ = 1.1051e-20, S″ = 4.3963e-23 | [cert] |
| C′ = με_I C″, S′ = με_I S″ | **C′ = 2.6523e-26**, S′ = 1.0551e-28 | [cert] |
| slack s = F0 − μK̄0 − με_I K2^pt (F0 = 1.312192e-21) | **s ≥ 1.0391e-21** (0.79·F0) | [cert] |
| strictness pad (C3): min over steps of pad·L/(U⁺·2C′) | ×1547 (α 1.5e-4), ×2236 (α 1.25e-4) | [cert] |

- K2^pt alone exceeds FLUX2 (×1.74). A separate bound K2 ≤ FLUX2 is therefore not claimed. It is not needed, because the verifier uses only the sum (design §1, §9.3).
- **Constants behind these numbers** [cert]:
  - κ⁽⁰⁾ = (8.46429, 4.17394, 2.26745). This is F1's κ⁽⁰⁾ with j ∈ [1 ± 3.3655e-3], |j′| ≤ 1.872e-3, and C_eff carrying the modulation ball. F1 itself is reproduced to 1e-12.
  - 𝔞₀..₃ = (10.397, 10.418, 8.973, 4.433); x₃ = 2.9312; z_i = sup|ζ⁽ⁱ⁾|, i = 1..4: (0.72917, 0.83482, 1.94446, 10.37038).
  - N̄ = (1.9925e-11, 3.2577e-12, 1.4620e-12) and M̄ = (1.7211e-2, 2.8139e-3, 1.2628e-3), both at β_A = 2 with d and chart included.
  - ∫_I M_3² ≤ 1.5798e-5 + 2.810e-7·T; ∫_I M_4² ≤ 3.6765e-3 + 8.462e-6·T; ∫_I s_2² ≤ 1.995e-8 + 1.263e-3·T.
  - W(57) ≤ 1.3403e-18.

## 2. Conditions C0–C9

| | condition | status | certified margin (lower bound) |
|---|---|---|---|
| C0 | s > 0 | PASS | s ≥ 0.79·F0 |
| C1 | S′ ≤ s | PASS | ×9.8e6 |
| C2 | L(τ) + 2C′ ≤ Ly_hi on every march cell | PASS | ×5.0e7 / ×5.3e7 (α 1.25e-4 / 1.5e-4) |
| C3 | U_k⁺·2C′ < 1e-9·max(1, \|U_k\|)·L_k | PASS | ×2236 / ×1547 |
| C4 | 2C′ ≤ η·L_F, η = 0.01, L_F ≥ 2F0/\|rate_s_U\| | PASS | ×1.7e5 / ×1.8e5 |
| C5 | whole-tail α-condition with rhs·(1 + 2C′/L_F) | PASS | ×1512 / ×1555 |
| C6 | C′ ≤ (κ_J − κ_J^cert)·X_const (v36 pad 1e-9; certified reserve 9.99999e-10) | PASS | ×141 |
| C7 | consumers of "Λ ≤ barrier" | PASS (scripted parts) | see below |
| C8 | regularity U-1 | **OPEN** | — |
| C9 | bootstrap order: all window constants at the (BΦ) level, i.e. exports doubled and E at β_A = 2 | by construction | — |

**C7 detail** (`flux2_cert.json` `res.C7`):
- **Scripted checks:**
  - The mod_sups tail row at L* + 2C′ stays under the 6-digit-rounded sup_ab/sup_beta: reserve ×2.5e8 or more.
  - Tail states lie inside the last march cell, L* + 2C′ ≤ Ly_hi_last: ×5e7 or more. This covers far_max, S_ramp and sqrtL_max, because those are maxima over march steps (closure.py:153, 308–310 [read]).
- **Consumers covered by C2:** the March consumers, which read Ly_hi.
- **Consumers covered by C6:** the J's.
- **Floor consumers:** these need a relative reserve of C′/L_F = 2.9e-8. The only uncushioned ones are PHI_S and PHI_C in assembly_arb.f4b (raw sqrt_LF and √L*). Per V36_REPORT §5, they feed only I3 and v3.1-only rows.
- **Flux lemma v2 does not use I3**: the integrated part is recomputed here from the exports. So no cushion is required, and none is applied.

## 3. Stages, checks, reproducibility
Each stage writes an `args` block (argv, env, every input with sha256). Timing goes to `*.timing.json` sidecars.

| stage | script → output | key checks |
|---|---|---|
| X1 | `flux2_forms.py` → `flux2_forms.{pkl,json}` | 35 [sym] checks all PASS. Checks run on: the §2.1 identity from the v-equation; the level-0 reduction with generic j (Euler check; no j″; only 4WZZ′J·g g‴ moved, with a vanishing boundary coefficient; affine in (j, j′); F1 reproduced exactly); D2v; the N_low structure from graph_struct.pkl and the 12-monomial N_low′; the conversions; the 5/2 trace residue. |
| X2 | `flux2_coef.py` → `flux2_coef.json` | OR sups ≤ the ZB values. F1 regression for Γ and κ⁽⁰⁾. W(57) regression: the old value is reproduced bitwise by a mirror of its code. C_eff > 13. |
| X3 | `chart4_shell.py` → `chart4_shell.json` | **The χ₄ assertion passes: b₂² = (u‴)² is the only quadratic high-jet monomial; there is no b₅ (u⁽⁶⁾); χ₄ = 0 for a ray.** The k ≤ 3 structure is re-asserted. The k ≤ 3 sups on the shell are ≤ the n1_chart core values. A float spot check by the geometric route agrees to 16 digits (diagnostic only). |
| X0 | `inputs_manifest.py` → `inputs.json` | 27 inputs pinned. Refuses to run without V36_REPORT. Rate-table env checked. v36 pins the same fin13v4b table. |
| X4 | `flux2_data.py` → `flux2_data.json` | Stop rule on C0_2. The f4b formula reproduces the pipeline's K0_w at β = 2 (rel ≤ 1e-9). N̄² is ≤ 1.05× the old formula. |
| X5 | `flux2_assemble.py` → `flux2_cert.json` | C0–C7. Exit code ≠ 0 on any failure. |
| X6 | `flux2_sanity.py` → `flux2_sanity.json` | Float only, not part of the proof: the reduced form equals the defining integral (3 tests); the bounds hold; the reflected-Z simulation behaves as expected. |
| X7 | `gen_texts.py` → `texts/` | Numbers are rounded in the safe direction. Brace and \begin/\end balance pass. Not compiled (no LaTeX). |

- **Double run:** two consecutive production passes (`snap/run4`, `snap/run5`, `cmp_run4_run5.json`). All JSON, pkl and tex outputs are byte-identical; logs differ only in timing lines.
- **Earlier runs** (`snap/run1`–`run3`): the outputs were already identical, apart from the sha of upstream timing fields. Those fields were then moved to sidecars.
- **Resources:** one process at a time, nice 19, PYTHONHASHSEED=0, `python3 -B`. Peak RSS is 135 MB (X3). A full pass takes about 4 min.
- **Hashes:** `SHA256SUMS`.

**Deviations from the design** (each is a valid alternative and is recorded in `texts/flux2_lemma_v2.tex` §impl):
1. **Coefficient sups on OR.** These are read from the pinned C2 output `c2_t3d.json`, sub-region [57, 60] (same code path) rather than re-run. Each is asserted ≤ its ZB value.
2. **W(57) ≤ 1.3403e-18** for the verifier's actual weight (S7 blend). The v3.1 value 1.43716e-18 was for a quintic blend. The (A)-chain flux block, which is records only, still uses the old value.
3. **Conversions in Minkowski form,** squared with fixed rational weights. The factor σ⁻¹ in ‖g″²‖ is kept.
4. **ā = 1e-6 is a typed hypothesis.** It is ≥ 2.19103e-7, so it covers every level in use. Sensitivity is about 1e-6 relative.
5. **Chart boxes on the shell** use the standing boxes B-a, B-d and B-e: α ≤ 2.5e-4, |w| ≤ ½, |w′| ≤ S0V, |w″| ≤ S1V/54.75.

## 4. What remains unjustified or open
1. **Regularity lemma U-1 (C8).** Λ absolutely continuous; g⁗(τ,·) ∈ L²(OR) a.e.; the level-2 identity with a C³ cutoff; Agmon for g″ a.e. This has been open since EI_SPEC §3, and the whole lemma depends on it.
2. **The level-2 export of (A), unit G1.** The certified computation exists (v36). I wrote the statement and a proof sketch in `texts/lemmaD_k2_sketch.tex` [hand]. Still open: a written, reviewed verification that every level-2 `Tr` source holds in window form. The code asserts this by construction, but there is no text.
3. **Verifier-side readings for C2–C5 [read/hand, not re-derived from bb.py]:**
   - "a cell bounds tot/Λ for all Λ ∈ [Ly_lo, Ly_hi]";
   - the first-touch arguments with Y = Λ − Z;
   - in C5, β nondecreasing and of degree ≤ 1 in Λ.
   A Fable review should check these against closure.py and bb.py.
4. **[hand] steps not machine-checked:**
   - chart lemma (d) at k = 4 ("same proof"; the structure itself is [sym] and the bounds are [cert]);
   - the weight reading "Ê_k bounds ∫(28/y) z_k² on the core", taken from the f4b rule;
   - the (BΦ) closed/open argument.
5. **Texts.** Not compiled; to be written: the E_I lemma hypothesis (I3) → (I3′), B1 → Λ ≤ L + 2C′, and the tail-floor quotes with factor (1 + 2C′/L_F)^½.
6. **No independent review yet.** The design calls for one Fable pass done from the scripts.

## 5. Protected trees
These are unchanged by sha256 before and after this work (`snap/before_*.sha` vs `snap/after_*.sha`): step3nl_v4b, fix3, fix4_r3, fix5, fix7_flux2 and the step-4 manifests. fix6_v36 is unchanged from its final state; see FLUX2_LOG for the final check. No new `__pycache__` was created. Nothing outside `D` was written.

## Files
- `code/`: X0–X8 scripts, `f2common.py`, `run_stage.sh`, `run_all.sh`, `cmp_runs.py`.
- `out/`: stage JSONs, the pkl, logs and timing sidecars.
- `texts/`: `flux2_lemma_v2.tex`, `absorption_lemma_v2.tex`, `erratum_flux2.tex`, `BPhi_continuation_v36.tex`, `lemmaD_k2_sketch.tex`, `numbers_flux2.json`.
- `snap/`: before/after shas and run snapshots.
- `FLUX2_LOG.md`, `SHA256SUMS`.
