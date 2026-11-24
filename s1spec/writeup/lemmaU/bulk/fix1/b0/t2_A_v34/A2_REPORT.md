# Task 2 (r2): (A) with R5′ in L² form, H-J-W, the (H2f) datum, point terms and localisation

**All acceptance rows of CONTRACTS_r2 pass at the `design` setting. No stop trigger trips.** Two changes need a look:
- The ramp coefficient box is tightened to the core box. Its certified margin is ×1.18 (A5).
- The R5′ addition to E₂ᴬ is 2.1e-3 at 4 × N2_ref: below the 0.01 trigger, but above Fable's 1e-6 to 1e-4 estimate.

Nov 13, 2026. Contract: `CONTRACTS.md` Task 2 as amended by `CONTRACTS_r2.md`. Everything is in `bulk/t2_A/`; `phase2a/` is untouched. `T2_REPORT.md` is the r1 stop report.

Labels: **[cert]** arb with outward rounding, this task. **[est]** reference value, not a bound. **[hyp]** hypothesis, not certified here.

## Rerun
- `./run_t2.sh` runs 9 stages in about 40 s at `nice 19`, with 93 MB peak RSS.
- Output is `out/A_interface.json` plus 18 stage JSONs, every one with an args block. `out/SHA256SUMS` holds the hashes.
- **Two consecutive full runs are byte-identical** (19 JSON).

## Results [cert unless labelled]

All runs use κ_J applied. The outer data are at N2_ref = (0.1, 0.2, 5.0) [est] unless stated. The datum is H2 + H2f, with `design` H2f = (5e-7, 5e-6, 1e-5) [hyp].

| | tight (H2) | v3.1 (H2) |
|---|---|---|
| Ê (at N2_ref) | 6.877e-7, 1.0388e-4, **0.05658** | 8.777e-7, 2.0288e-4, 0.25558 |
| EA0 (no outer source) | 6.839e-7, 1.0384e-4, 0.05607 | 8.739e-7, 2.0284e-4, 0.25507 |
| Ê at 4 × N2_ref | 6.991e-7, 1.0401e-4, 0.05814 | — |
| γ′ (lower) | 0.242926, 0.241064, 0.119615 | same |
| C0 (z₂, z₃) | 7.0652e-7, 1.6544e-4 | 8.623e-7, 2.4663e-4 |
| C1_0 → C1_ref | 3.83e-12 → 7.62e-10; 1.56e-11 → 8.44e-9 | same |
| S0v / S1v | 6.177e-3 / 0.1034 | 7.724e-3 / 0.1773 |
| S0v against 0.0153 | **×2.477** | ×1.981 |
| core M1, M0 | **×4.443**, ×5.737 | ×2.723, ×5.309 |
| core Ĝ₂₂ (vs 153664) | 80141.6, **×1.917** | same |
| K0(Ê), K2(Ê) | 1.370e-21, 4.002e-21 | 1.789e-21, 6.151e-21 |
| FLUX0/(β·K0), FLUX2/(β·K2), β = 5/4 | **×2.394, ×3.998** | ×1.834, ×2.601 |
| same, β = 2 | ×1.496, ×2.499 | ×1.146, ×1.626 |
| FLUX0_needed, FLUX2_needed (4K) | 5.48e-21, 1.60e-20 | 7.15e-21, 2.46e-20 |

**κ°** (o_j = Σ_k κ°_jk N̄_k² per unit τ):

| level | N̄₀² | N̄₁² | N̄₂² |
|---|---|---|---|
| 0 | 5.84e-9 | 1.70e-9 | 0 |
| 1 | 4.98e-11 | 3.44e-8 | 6.79e-10 |
| 2 | 6.23e-9 | 3.91e-8 | **1.242e-5** |

**(H2f) settings** at the tightened (H2) (A4/A4b):

| setting (H2f) | Ê₀ | FLUX0/K0 at 5/4 · Ê, at 2Ê | FLUX2/K2 at 5/4 · Ê, at 2Ê | S0v (vs 0.0153) | FLUX0_needed |
|---|---|---|---|---|---|
| design (5e-7, 5e-6, 1e-5) | 6.88e-7 | ×2.394, ×1.496 | ×3.998, ×2.499 | ×2.477 | 5.48e-21 |
| lo (4e-7, …) | 5.88e-7 | ×2.789, ×1.743 | ×4.557, ×2.848 | ×2.581 | 4.70e-21 |
| hi (6e-7, …) | 7.88e-7 | ×2.097, ×1.311 | ×3.561, ×2.226 | ×2.390 | 6.26e-21 |
| probe (3.5e-7, 1.1e-6, 2.7e-6) | 5.38e-7 | ×3.043, ×1.902 | ×4.922, ×3.076 | ×2.666 | 4.31e-21 |

- Parts of Ê₀ at the design setting: hand-off 5.1e-7, inflow 1.12e-7, ramp 2.97e-9, TOP 5.17e-8, CLS 7.59e-9 (A3), outer 3.81e-9.
- At the hi setting the FLUX0 margin at β = 5/4 is ×2.10, which is below the 2.2 report line.
- E₂ᴬ fell from 0.0726 to 0.0566. The ramp-rate cost dropped from 0.0287 to 0.0117 because of the A5 repair (smaller ramp Ĝ₂₂, 95e3 against 239e3).

### Acceptance (CONTRACTS_r2)

| quantity | threshold | value | status |
|---|---|---|---|
| FLUX0/(β_A K0), FLUX2/(β_A K2), β_A = 5/4, design | ≥ 2.0 (report < 2.2; stop < 1.5) | ×2.394, ×3.998 | pass |
| S0v against 0.0153 | ≥ 2.0 | ×2.477 | pass |
| core M1 | ≥ 2.0 | ×4.443 | pass |
| E₂ᴬ at N2_ref | ≤ 0.1 (report > 0.08) | 0.0566 | pass |
| ô₀/γ′₀ at 4 × N2_ref | ≤ 1e-7 | 1.52e-8 (3.8e-9 at N2_ref) | pass |
| R5′ addition to E₂ᴬ at 4 × N2_ref | ≤ 0.01 (stop trigger) | 2.08e-3 | pass |
| every class on O closed by T1–T7, P | yes (report if not) | yes (asserted) | pass |
| H-J-W | κ_J ≤ 0.13 | 0.042348 | pass (closed) |
| core Ĝ₂₂ margin after A5 | ≥ 1.5 (report < 1.7) | ×1.917 | pass |
| retained r₀, r₁ on [28, 350] | ≥ 1.249 a_c | 1.249 (all zones) | pass |
| segment conditions, both boxes | all true | true for (0.15, 0.52) and (0.20, 0.60) | pass |
| reruns | byte-identical | yes | pass |

- The r1 row "K at 2Ê, margin ≥ 2" is superseded. At β = 2 the margins are ×1.50 and ×2.50, so FLUX0 fails ×2 at β = 2 (as R2 predicted).
- At the v3.1 (H2), S0v vs 0.0153 is ×1.98 and the flux margin at β = 5/4 is ×1.83. Both are for the (O) note only.

## What each item did

- **A1, H-J-W** (`hjw_rep.py`; closed in r1). Verifier `rate_table_fin11_a0.00015.json`.
  - (W1) fails at every cell with P_k ≥ 0. (W2) holds; it relies on the verifier's tail lemma, which was not re-run.
  - (W3) needs F0 added to c_k: **erratum E-V1**, `closure.py:182`, `rhs_iv.py:198`, `unc.py:17–27`.
  - κ_J = 0.042347 (W1) + 4.8e-7 (W3) = **0.042348**. Every J_m enters (A) as J_m(1 + κ_J).

- **A2, R5′** (`outer_L2.py`, new; `outer_rep.py` is no longer used).
  - The class enumerator of `ramp_rep.py`, mirrored to η_o = 1 − S. Lemma P: case (iii) at levels 0 and 1 (no integration by parts on O); case (i) at level 2 (top pair and (2,3) integrated by parts; net point term −g_c z₂(350)² ≤ 0 dropped).
  - **Class list on O and tools**, all asserted:

| level | class | tool |
|---|---|---|
| 0 | (0,1), (1,1) | T2 |
| 1 | (0,1), (1,1), (1,2), (2,2) | T2 |
| 1 | (1,3) | T3 |
| 2 | (0,2), (1,2), (2,2) | T2 |
| 2 | (2,4) | T3 |
| 2 | (2,3) | envelope e = εω t(1−t)⁴, then AM-GM, then T5 on θ_* |
| 2 | (2,2,2) | T6 (θ = K·e), then T2 plus the (2,3) chain |
| 2 | (2,2,3) a′, b | T7 (κ continuous), then (2,2,2) |

  - Further sources on O:
    - the interpolation remainder ½Ỹ_o εω |(η_o²)″ − 2κ(η_o²)′|: T2, 2.07e-7 per unit N̄₂² at level 2;
    - the L⁴/S33 outer foot: its base halves go into the O rate; the data halves cost 9.9e-8 per unit N̄₂² (`weights_t2.py`).
  - Reserve δ° = ¼ a_o; used 0.2496 at level 2.
  - **κ°₂₂ = 1.24e-5 is 82% the (2,3)/(2,2,2) chain**, made of the envelope constant K, then T6 and T5.
  - No pointwise bound on z₂ or z₃ is used, and no datum on z₃ or z₄. No `lam` and no `eps4`.
  - The level-0 outer source is ô₀/γ′₀ = 3.8e-9 at N2_ref and 1.5e-8 at 4 × N2_ref.

- **A3, point terms** (`a3_signs.py`, `ramp_rep.py` route).
  - On the stored forms (`lemmaS_F1.div_form`), the coefficient of (2,2,0) on (j+1, j+1) is +4 and that of (3,2,0) on (j, j+1) is +2, at j = 0, 1, 2. The stored form is the dE_j/dτ contribution: (n,n) is −2.
  - So Lemma P holds as stated: case (i) gives −g_c z_j² ≤ 0, and case (ii) leaves +g_c z_j².
  - Levels 0 and 1: the m = 3 monomial of (j, j+1) is now bounded directly by data (space-time Cauchy–Schwarz with conv_j and conv_{j+1}). Level 2 keeps case (i).
  - g_c(28) = 5.862e-4. dE = (+6.41e-9, +3.00e-6, +9.7e-9), i.e. +0.9% of E₀ᴬ and +2.9% of E₁ᴬ.

- **A4, initial datum** (`ic_rep.py --h2f`). E^A_j(τ₁) ≤ H2_j + H2f_j, with no H3 and no zone-III term. Step 4 certifies (H2) on [14, 60] only, so **(H2f) is [hyp]**: task 6 must certify it.

- **A4b, flux at the bootstrap level.** K0 and K2 are linear in E; this was asserted by recomputing f4b at β·Ê. Values are in the tables above.

- **A5, interface at y = 28 and segments** (`weights_t2.py`, new; `weights_rep.py` is kept for reference).
  - **Convention found.** v3.1 and phase 2a take the interpolation, L⁴ and S33 identities separately on the ramp (box (0.05, 0.9)) and on the core (box (0.11, 0.52)), each with its own constants, and drop the boundary terms at 28. This is not covered by L(i):
    - Y_ramp[3] = 1156.7 against Y_core[3] = 644.2;
    - there is no single μ, so L(iii) fails.
  - **L(ii) with a transition on [28, 56] is impossible.** It would make the core Ĝ₂₂ at y = 28 equal 2.4e5, above 153664.
  - **Repair:**
    1. The ramp coefficient box is set to the core box (0.11, 0.52), with z₀ ≤ 1/50 and z₁ ≤ 0.52 on the ramp. This is valid because the verifier certifies |y f_y| ≤ 0.4396 and |z₁| ≤ 0.4328 there. **The new ramp M1 margin is ×1.185** (it was ×2.05 against 0.9); ramp z₁ is ×1.20.
    2. One C^{1,1} Ỹ_j covers [14, 700]. It rises on [40, 50] and on [175, 250]; every cell is checked to satisfy Ỹ ≥ coefficient and ε·Ĝ ≤ ε(28)·Ĝ_core.
    3. The L⁴/S33 identities are localised: piece P1 on [14, 350] with core constants (L(iii)), and piece P2 on O, joined by a partition of unity on [250, 350]. The S33 weight rises C¹ on [50, 60].
    4. Ramp-foot halves now act on all of P1.
  - Results:
    - core Ĝ_jj = 4347.7, 5492.1, 80141.6 (phase 2a: 4347.7, 5288, 79364.6);
    - γ′ = (0.24293, 0.24106, 0.11962), against (0.24292, 0.24139, 0.12087) before;
    - ramp Ĝ₂₂ = 95e3 with B_c = 1.37 (slack 0.326 a_min for the R1/R2 Young steps).
  - **Segments.** The box (0.15, 0.52) holds on (60, 350] with z₀ ≤ 0.161 (B3U(c) plus δU) and z₁ ≤ 0.61. The fallback (0.20, 0.60), with z₁ ≤ 0.69, also holds.
    - Worst zone ratios ε·Ĝ/(ε(28)·Ĝ_core): 0.27 on [40, 50]; 0.15 on [50, 60]; 0.08 on the plateau; 0.03 on [175, 250]; 0.02 on [250, 350].
    - O passes the (S2) condition with B_c lowered by δ° = ¼.
  - **(O1) can be relaxed from 0.11 to 0.15** on (60, 350]: the segment checks are certified at (0.15, 0.52), and also at (0.20, 0.60).

- **A6, export.** `out/A_interface.json`, with the r1 keys plus the r2 additions:
  - `IC` = {H2, H2f, E_tau1, setting};
  - `core` gains beta_A, K0/K2 at Ê and at β·Ê, FLUX*_needed and `S0v_margin_vs_0p0153`;
  - `H2f_settings`, `v31_only`, `classes_outer`, `outer_other`, `lemmaN`.
  - Each number is exported as an arb string and as an outward double, with a label.
  - The rule Ê = EA0 + (o + 2·coupling)/γ′ reproduces the assembled EA_ref (asserted). The coupling term is below 1e-4 of o/γ′.

- **A7, audit.**
  - (i) `R['S0']` and `R['S1']` bound sup|v_A| and sup|y v_A′| on **[28, 60]**: conv = (60/28)^γ and ℓ = ln(60/28), with η_A = 1 there. This contains the verifier's [34.25, 60].
  - (ii) `S0_B (y_E = 50) ≤ S0V` is 0.01478 at design, ×1.035 (close to 1). It fails at `hi` (0.0153) and at v3.1 (H2) (0.0185). **No statement of `STATEMENTS_r2.tex` consumes it**: a grep for S0_B, (B), exit and y_E finds only D1's S0^v and S1^v on [34.25, 60], which are R['S0'] and R['S1']. It is the (B)-exit of v3.1, replaced by Lemma U. It is kept under `v31_only`.

## New errors (file:line)
- **Verifier:** `step3nl/verifier/closure.py:182` omits F0·D_k on march cells (E-V1, already in R2_DECISIONS).
- **v3.1, interface at 28:**
  - `writeup/lemmaS/rigor/weights_exact.py:246, 254` (also `base_ibp` copy): `core()` and `ramp()` build separate L⁴/S33 constants and separate YC[j+1].
  - `writeup/lemmaS/rigor/ramp_arb.py:29, 57`: the ramp interpolation remainder uses the ramp-only YCr[j+1] while the core uses its own. The boundary terms at 28 (interpolation jump and [Θ z₁ z₂³]) are dropped.
  - Same in phase 2a: `phase2a/A_repair/rigor/ramp_rep.py:187`, `weights_rep.py` (`core`/`ramp`).
- **Phase 2a segment check:** `phase2a/A_repair/rigor/weights_rep.py:149` runs [220, 350] with sup|z₀| = 1/50 and sup|z₁| = 0.52. Under B3U these are 0.161 and 0.61. It is repaired here; it also passes with them.
- **Phase 2a `A_REPAIR.md` R1** ("c, c′ vanish at t = 0, 1") is false as written. Level 2 is still correct by Lemma P(i).
- **Not re-audited:** the (B) windows (`window()`, `weights_exact.py:277`) use the same per-region convention. They are v3.1-only, feeding `win M0/M1`, K2/K3 and the (B) exit, and r2 does not consume them.

## Still [est] or [hyp]
- **[hyp]:**
  - (H2f) = (5e-7, 5e-6, 1e-5); task 6 must certify it on step-4 data.
  - The outer data N̄_k², supplied by Lemma N (B3U at 2B̂, Z3).
  - The box H-O1 (½, 1) on O (Z2).
  - The segment box (0.15, 0.52) and z₀ ≤ 0.161 on (60, 350] (B3U(b, c), Z1).
  - The core box B2′ (bootstrap).
- **[est]:** N2_ref = (0.1, 0.2, 5.0), a reference only; the export is affine in N̄².
- **Taken as given, not re-run:**
  - the verifier's tail lemma (W2);
  - the verifier inputs FLUX0 = 4.1e-21 and FLUX2 = 2.0e-20;
  - the ramp box margins, which come from the verifier's M1_ramp = 0.43277.

## Files
- **Code** (`rigor/`):
  - new: `hjw_rep.py`, `ic_rep.py` (r1 version kept as `ic_rep_pre_r2.py`), `a3_signs.py`, `weights_t2.py`, `outer_L2.py`, `export_t2.py`;
  - changed: `ramp_rep.py` (backup `_pre_t2`), `assembly_rep.py` (backup `_pre_t2`).
- **Driver:** `run_t2.sh`.
- **Outputs:** `out/*.json`, `out/*.log`, `out/SHA256SUMS`.
- **Records:**
  - `probe/` holds the r1 stop runs and their hashes;
  - `ref_phase2a/` holds the phase-2a outputs, a rerun copy, and the old `run_chain.sh`.
