# Task 4: arb certification of Lemma U. All acceptance rows pass except (R<2) ≥ 5 (×4.43, report only). No stop trigger.

Opus, Nov 13 2026. Contract: `CONTRACTS.md` Task 4 as amended by `CONTRACTS_r2.md`. Only `t4_cert/` was written. Other tasks' code was imported read-only (`-B`).

**Inputs**
- (A): `t2_A/out/A_interface.json` (task 2, certified), setting `tight` (design H2f), β_A = 5/4.
- Zone III: `t5_zone3/out/zone3_interface.json` (task 5; it is itself marked provisional pending this file).

**Labels**
- **[cert]**: arb with outward rounding, project script.
- **[paper]**: a hand step, stated in the script docstrings.
- **[hyp]**: a hand-off threshold.

## Headline (primary box (0.15, 0.52), U_max 0.09, M0z 0.16, data56)

| quantity | value | threshold | status |
|---|---|---|---|
| late M1 (vs 0.52) | **×2.358** [cert] | ≥ 1.5 (stop), report < 1.8 | pass |
| late M0 (vs 0.15) | ×5.166 [cert] | ≥ 2.0 | pass |
| present M1 min, y_c ∈ [56, 168.67] | ×2.775 at [56, 60] [cert] | ≥ 2.0 | pass |
| present M0 min | **×2.226** at [159, 168.67] [cert] | ≥ 2.0, report < 2.1 | pass |
| present / late M0z (vs 0.16) | ×14.99 / ×14.96 | — | pass |
| loop gain G₀, G₁, G₂ | 1.2e-7, **0.0372**, **0.0390** [cert] | ≤ 0.5 (stop) | pass |
| — of which through neighbours / through (A) | G₂: 0.0370 / 0.0020 | — | — |
| rate margin (R2) | ×28.23 [cert] | ≥ 10 | pass |
| (R<2) margins k = 0, 1 | **×4.43**, ×6.70 [cert] | ≥ 5 | **fails at k = 0 (report)**, see §3 |
| θ(28) | 0.3701 | ≤ 3/5 | pass |
| φ effect on Γ_m / Ĝ / K | ≤ 3.4e-5 / ≤ 3.2e-5 / ≤ 1.2e-5 relative | ≤ 1e-3 | pass |
| S33 + foot (incl. the (2,2,3) term) at level 2 | 1.7e-3 of B̂₂ | ≤ 1e-2 | pass |
| enumerator cross-check | no class missing; max weight difference 1.2% (A ≤ P2) | 5% | pass |
| slope part C1·Σ_X / 𝓓_k | 2.2e-3 (k = 0), 1.4e-4 (k = 1) | report | — |
| N̄² / (4 × N2_ref) (coordinator's check) | 0.179, 0.413, 0.252 | ≤ 1 | pass |
| N̄² / N2_ref | 0.71, **1.65**, **1.007** | — | above N2_ref at levels 1 and 2; Ê is evaluated at the actual N̄² |

**Supersolution** [cert]
- B̂ = (4.860e-6, 5.567e-4, 0.36560).
- Ê^A = (6.879e-7, 1.0390e-4, 0.05659).
- C0 = (7.07e-7, 1.654e-4); C1 = (7.6e-10, 8.4e-9).
- (K1) A-part, (K2) and (K3) all hold in arb. Core strictness and fluxes are task 2's rows; they are copied into `U_cert.json`.

**U-fine B̂ for task 5:** (4.856e-6, 2.626e-4, **0.02578**) [cert]. Requirement: B̂₂ ≤ 0.106, B̂₀ ≤ 5e-6, B̂₁ ≤ 1e-3. All are met. The thinnest is B̂₀, at ×1.03.

**Fallback box (0.20, 0.60), M0z 0.20** [cert; Δ^R taken from task 5's primary-box run]

| quantity | value |
|---|---|
| late M1 / M0 | ×2.635 / ×6.79 |
| present M1 / M0 | ×3.08 / ×2.97 |
| G₂ | 0.039 |
| rate margin | ×13.80 |
| g | 768479 (K certified on [6.2e5, 7.7e5]) |
| θ(28) | 0.556 |

**Late B₂ = 0.3656, by part:**

| part | value |
|---|---|
| launch | 0.1405 |
| inflow | 0.0720 |
| X-terms | 0.1401 |
| outer Δ^R | 0.01148 |
| rampE | 6.7e-4 |
| foot+S33+(2,2,3) | 6.2e-4 |
| rs | 2.5e-4 |
| couplings | 1.3e-6 |
| b | 4.9e-7 |
| forcing Φ | 3.8e-8 |

## 1. Files (`t4_cert/`)

| file | content |
|---|---|
| `t4lib.py` | helpers (arb export with string + outward double + label, args block, verifier sups) |
| `u1_ghat.py` → `out/u1_ghat.json` | U1: Γ_m, κ̃, κ̃₂ of Θφ; Ĝ_jk at both boxes; g, D(56), D(28), θ(28), (R2), (R<2) |
| `u2_K.py` → `out/u2_K{,_nophi,_fallback}.json` | U2: K_jk with φ (branch and bound, g a dimension), folded classes, footw, c33, c223, C_η, L_max, R |
| `rampphi.py` | the φ ramp geometry (same formulas as u2), used by u5 |
| `u34_sweeps.py` → `out/u34_{primary,fallback}.json` | U3: Σ_in, Σ_E, S_X, Σ_X, Σ_rs on y_c ∈ [56, 253]. U4: C_A^D, q, q^D, w_min, c₂, c₃ per present cell, Δ^in, c_δU, Φ_j, Sobolev |
| `u5_enum.py` → `out/u5_enum.json` | U5: cross-check against (A)'s enumerator (t2_A/rigor/ramp_rep.py) |
| `assemble_U.py` → `out/U_cert.json`, `out/U_cert_fallback.json`, `out/U_cert_for_zone3.json` | U6 |
| `run_all.sh [A_interface path]` | the whole chain, about 5 min; writes `out/SHA256SUMS` |

- Every JSON has an `args` block.
- Two consecutive runs of the whole chain are byte-identical (§5).
- Peak RSS is 80 MB; one process at a time; `nice 19`.

## 2. Certified constants

**U1** [cert]
- Γ(Θφ) = (4.0000357, 16.000321, 64.002179); κ̃ = 8.0000357; κ̃₂ = 64.000607. These come from P1's rerun plus exact Leibniz/triangle bounds with β = K_b e^{−s} ≤ K_b/28.
- Primary Ĝ_jj = (51322.5, 33962.9, 348341.0). This includes φ, task 5's outer-foot Y adds and M1z = 0.534 (§3).
- The contract convention (no φ, M1z = M1) reproduces 338258.77 (asserted).
- g = 361147; D(56) = 1.03672; D(28) = 1.58756.

**U2** [cert]: K table (K_aug, folded cubic classes included):

| j \ k | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| 0 | 8034.2 | 724.88 | | |
| 1 | 0.0815 | 6336.9 | 545.58 | 0.587 |
| 2 | 0.872 | 528.62 | 23426.0 | 555.62 |

- footw = (0, 2.636, 2663.3).
- c33 = 108.02: the exact foot W33·C_η/3 = 35.67, plus (2,2,2) = 72.35.
- c223 = 7.328 (Cauchy–Schwarz route; §3).
- C_η = 13.375; L_max = 3.5351; R = 11.593.
- g-range [3e5, 6.4e5] for the primary, which also covers the contract fallback's 6.29e5. The honest fallback (g = 7.68e5) has its own run on [6.2e5, 7.7e5].

**U3** [cert]: sweeps at y_c = 56. Every sweep is nonincreasing in y_c (paper: W increases on σ ≤ l/2), so the cover sup over [56, 253] is the value at 56. The grid check is in the JSON.

| sweep | value |
|---|---|
| Σ_in | 5.0892 |
| Σ_E | 0.30846 |
| S_X | 1.5242 |
| Σ_X | 3.0080 |
| Σ_rs² | 4.33e-3 |
| Σ_rs⁰, Σ_rs¹ | 0 under (R<2), asserted |

The x-tail is bounded analytically with W ≤ e^{4(σ + a)}.

**U4** [cert]

| constant | value |
|---|---|
| w_min | 0.734921596 (= W(0) = W(l)) |
| C_A^D | 2.24284 |
| q | 0.018508 |
| q^D | 0.018508 |
| **c_δU** | **1.718604e-8** (from the verifier's sup_β = 3.33811e-9 and sup_a = 2.19103e-7; below task 5's placeholder 2e-8) |
| Δ^in | (3.13, 3.29, 3.45)e-12 |
| ∫W e^{−2σ} | 1.28799 |
| Φ_late | (≤ 1.4e-8, …, 3.8e-8) |
| Sobolev | sup u² ≤ ‖u‖²/l + 2‖u‖‖u′‖, ‖u‖² ≤ 𝔈/w_min |

c₂ and c₃ are given per present cell.

**U5** [cert]
- Direct route (A's code with every pair kept, no IBP): every column and cubic class matches P2's.
  - The pair columns differ by at most 1.2%, always A ≤ P2: A sums its entries before the sup, while P2 folds per key after the sup.
  - The cubic classes are identical.
  - Nothing is unmatched.
- (A)'s native route J creates two extra classes: the pair (0,2) at level 0 and the pair (2,4) at level 2. Both are artefacts of its integration by parts of the pairs (a, j+1), which brings in the dissipation index j+2. Lemma U keeps those pairs as data (column k = j+1 of K), so it creates neither. The list is in `out/u5_enum.json` under `native`.

## 3. Things the parent must know

1. **Convention fix: sup|z₁| ≠ M1.**
   - The Ĝ pipeline (ghat_arb, ghat_design_arb, P2's folding, task 5's ramp_U) uses the box value M1 = 0.52 both for |f_s| in the A-coefficients and for the factors z₁ of the identity.
   - But z₁ = y f_y − p U₁, so the box gives only |z₁| ≤ 0.52 + p U₁. I introduced a separate M1z = 0.534, justified region by region:
     - core [28, 60]: 0.52 + p U₁(60) + p c_δU 60^p = 0.5333, asserted;
     - zones I/II: S₁(2B̂) = 0.285, asserted;
     - zone III: b₁ + p U₁(379.5) = 0.508, asserted.
   - Cost: Ĝ₂₂ +2.9%; rate margin ×29.07 → ×28.23.
   - Task 5's outer ramp is unaffected locally, because |z₁| ≤ 0.508 < 0.52 there. Its global Y was built with 0.52; the effect on Δ^R is below 1e-3 relative [est]. Task 3 should state the M1z hypothesis.
2. **(R<2) margin at level 0 is ×4.43 < 5** (acceptance row; not a stop trigger).
   - The inequality holds strictly. It enters only the coupling surplus (couplings 1.3e-6 of B₂).
   - It could be raised by taking g above the (R2) minimum: at g = 6e5 it is about ×5.9 [est], but D(28) also grows. I did not change it.
3. **(2,2,3) by Cauchy–Schwarz, not by Lemma C's integration by parts:** T223 = c223·(∬εΘz₂⁴)^{1/2}(∬εΘz₃²)^{1/2}.
   - The L⁴ factor is bounded by the same sup-free Sobolev form as S33.
   - This avoids a point term (the m = 3 jump of (η²)‴ at y = 56). The class has only m = 1 in fact (recorded), so the integration-by-parts route would also be valid.
   - T223 is about 3e-4 of B₂.
4. **K_S.** With the recertified K the honest multiplier is K_S = 1. The ×1.33 of P2 compared against feas_P. This is why late M1 is ×2.36 rather than the [rev] ×2.28.
5. **(H2) and (H3) for present centres** are typed as the hand-off thresholds [hyp]: (1e-8, 1e-6, 1e-3) and (2e-7, 1e-5, 1e-5). α₁ = 1.5e-4.
6. **zone3 loader.**
   - `zone3.py --ucert` reads `supersolution.Bhat`, but that key is the global crude B̂ in the contract. So the U-fine B̂ is written to `out/U_cert_for_zone3.json` under that key. Task 5 should rerun with `--ucert .../t4_cert/out/U_cert_for_zone3.json`.
   - I did not run zone3.py, because it writes into t5_zone3/.
   - At B̂₂ = 0.0258 the H6 margins should exceed task 5's fine-estimate row (≥ ×1.86 at B̂₂ = 0.033).
7. **Task 5's Δ^R** is multiplied by D(253)φ(253) (conservative, ≤ 1 + 9e-5), and DeltaR_b is added.

## 4. What is still provisional or not mine
- **The v3.3 interface (`t2_A_v33/out/A_interface.json`) is not yet a v3.3 result.** It is byte-identical to t2_A's (sha256 4f86e3cb…: β_A = 5/4, H2f₀ = 5e-7). Rerun with `sh run_all.sh <v33 path>`, or only the assembly with `--A <path> --out U_cert_v33.json` (under 1 s), once it changes.
  - β_A enters only task 2's flux rows. Lemma U uses Ê, C0 and C1 at β = 2 (Lemma N).
  - A larger H2f₀ raises Ê₀ only, which is negligible here: the level-0 margins are not binding.
- `zone3_interface.json` is marked provisional by task 5 (it waits for this U-fine B̂). Δ^R, n_III and b₁ are certified given that file.
- (H2f), (H2), (H3) and α₁ are hand-off thresholds [hyp]; task 6 certifies them.
- c'_j of Prop. B′ (0.0520, 0.1446, 0.4324) are taken from v3.2 [paper input].
- (K1)'s core-strictness and flux rows are task 2's; they are copied, not recomputed.

## 5. Reproducibility
`run_all.sh` was run twice; `out/SHA256SUMS` is identical between the two runs (no timing fields in the JSONs).
