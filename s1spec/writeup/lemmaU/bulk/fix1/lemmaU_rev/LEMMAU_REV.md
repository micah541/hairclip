# Lemma U proof text, revision v34 (fix round 1)

Nov 14, 2026. I wrote only in `fix1/lemmaU_rev/`. `t3_proof/lemmaU.tex` is unchanged (sha256 1080d051…).

**Result.**
- `lemmaU_v34.tex` passes the balance check (`t3_proof/check_balance.py`: PASS, 36 labels, 26 refs).
- `check_keys_v34.py` passes every check except one: the S4 hash. The step-4 manifest was rewritten by a running M9 driver while I worked (see P6).
  - 6 of 7 source hashes verified; S4 is pinned at 8202d2c6…, the manifest the text was written against;
  - all row keys resolve, including on the new manifest;
  - all 131 distinct `\key{}` citations are covered;
  - 177 quoted numbers are rounded outward against the certified values (exact rational comparison);
  - the side conditions hold (see the last section).
- No margin fails.
- One hypothesis-level wording problem was found (B2′, item P1 below) and fixed in this text.

## Files

| file | content |
|---|---|
| `lemmaU_v34.tex` | revised proof text (copy of `t3_proof/lemmaU.tex`, edited) |
| `constants_needed_v34.json` | same format as `t3_proof/constants_needed.json`. Changes: (1) keys are full paths from the file root; (2) 7 sources with sha256; (3) 103 rows; (4) a new `numbers` list of 177 quoted values, each with its outward direction |
| `make_constants_v34.py` | writes the JSON |
| `check_keys_v34.py` | read-only checker (exit status 1 on failure). Run: `nice -n 19 python3 -I check_keys_v34.py` |

Sources and their sha256:

| tag | file | sha256 |
|---|---|---|
| SA | `fix1/b0/t2_A_v34/out/A_interface.json` | 34832f43… |
| SA₄ | `fix1/b0/t2_A_v34/out/asm_design_N2ref4.json` (pinned in `t2_A_v34/out/SHA256SUMS`) | 6ad175c1… |
| SU | `t4_cert/out/U_cert_v34.json` | 271430d7… |
| (cited by name) | `t4_cert/out/u1_ghat.json` | 5a06d24c… |
| SZ | `t5_zone3/out/zone3_interface.json` | a381b1f2… |
| (cited by name) | `t5_zone3/out/ramp_U.json` | 55317ce5… |
| S4 | `step4/p2/out_v33/FINAL_manifest_v33.json` | 8202d2c6… |

## Changes, with line numbers in `lemmaU_v34.tex`

### Header and sources
- l.14–15: the title says revision v34; the date line was updated.
- l.18: the purpose paragraph now points to `lemmaS_A_v33.tex` and `continuation_v33.tex`.
- l.22–29: the source list was replaced with SA, SA₄, SU (+u1), SZ (+ramp_U) and S4, with the M14 hashes. The superseded files are named. There is a pointer to the constants file.

### §1 Setting
- l.36: ā and β̄ are read from the fin12f4 file (M4). zone3 reads ā from the verifier (M3).
- l.47–48: Definition (C5) now contains w′ and the profile difference d (B0).
- l.62, remark "why (C5)". **B12 §4 item 1.** The old claim "(H2) already stated for f − U₁" is replaced:
  - (H2) is certified for f − U₁ by `h2u_v33` (shells of w′ plus d by Minkowski);
  - the bound is (3.658e-10, 1.777e-9, 1.665e-6), with margins ×27.3, ×562, ×600 [cert, S4 `results.H2_tightened_fU1`].
- l.64: a new remark, "where d enters". d enters (A) only: its ramp data, S0V/S1V, the fluxes and the ramp boxes. Lemma U consumes only f − U outputs and uses δU. The d-terms are rates, so they go into C1.
- l.69: a_W = 3/10 is contained exactly in the certified covers (M3).

### §2–4 Lemma T, launch, H-J-W, Lemma D
- l.94: assemble_U uses the Γ₁ enclosure (old G3, resolved).
- l.124: launch (a). H₂ is in the variable f − U₁ and is certified by h2u_v33.
- l.146: (W2) was rerun for fin12f4: tail_cert PASS and closure tail=True (`V33_VERIFIER.md` l.6–7). The statement "not rerun [hyp]" is removed (M5).
- l.157–166, **Lemma D, window form** (B12 §4 item 2):
  - (D2) now reads ∫𝔰 ≤ D^c + D^s·T, with D^c ≤ D^src;
  - C1 contains D^s, not ρ_J·D^src;
  - a new paragraph, "Why the slope is not ρ_J·D^src", covers the η-split and the d-rates.
- l.169, Lemma D inputs (all v34 values):
  - F², coup, D^src (1.8122e-7, 9.8267e-5);
  - D^s (9.1216e-10, 5.1442e-8) (`B0.D_window_slope`) and D^c ≤ (1.7747e-7, 9.8053e-5) (SA₄);
  - C1⁰ (7.5191e-10, 4.2195e-8), with the reason for the increase;
  - C0 (1.1345e-6, 1.6621e-4) and C1 (1.5554e-9, 5.4307e-8).
- l.198: late 𝒟 is (1.7339e-6, 2.5351e-4); the slope shares are 2.7e-3 and 6.5e-4.

### §5 Proposition U proof
- l.263–265, Step 0 split (M5 typo): η²χ²(Θφ)^(m) replaces η²(χ²Θφ)^(m). rem^χ is now written out, and the exactness of the split is stated: the ramps are disjoint, and η = 1 on supp χ′ and χ = 1 on supp η′.
- l.277–286, **Step 1 reserve argument (M5; B12 §4 item 8).** The "used locally / disjoint supports" argument is replaced:
  - level-2 slack is 0 (`u1_ghat.json cases.primary.budget_slack`);
  - L⁴ identity with inner and outer foot terms;
  - X_in² + X_out² ≤ X² (disjoint supports), then Cauchy–Schwarz in ℝ², giving X ≤ M₁ᶻ(3N₃ + κ̃N₂) + (b_in²F_in + b_out²F_out)^{1/2};
  - one Young step, (a² + b²)^{1/2}·N₄ ≤ δ_r·N₄² + c_in²F_in/(4δ_r) + c_out²F_out/(4δ_r). These are exactly the two data terms the scripts charge (u1 footw; ramp_U l.86–88, c₁²/(4·dd)).
  - Two separate Young steps would overdraw the budget by a_min/10. The ramp_U docstring gives the right conclusion with the wrong reason; that is noted.
  - The N₃ parts and the L⁴ classes carry no N₄, so they are charged separately (√(a² + b²) ≤ a + b).
  - The outer z₁ factor uses M₁ = 0.52, which is valid because the zone-III bound is 0.50786.
- l.286: outer-foot Y add is +0.1081 (was 0.108, which is not outward).
- l.292: Ĝ₀₀ ≤ 51322.6 (was 51322.5 < 51322.514, not outward).
- l.302: c₃₃ ≤ 108.02.
- l.319: Step 8 pairs cite `ramp_U.json K`.
- l.327: Step 8 interpolation remainder (M1). K°_jj uses the same primary Y_next (457.092, 356.940, 1320.220) as Ĝ, asserted to dominate the contract Y. K°_jj ≤ (13436.03, 10526.19, 39482.30).
- l.333:
  - ΔR is (4.4986e-8, **6.7972e-6**, 0.011480);
  - the level-2 parts are 0.011298 and 9.3e-6 (were 0.011297 and 9.2e-6, not outward).
- l.344: S₁(2B̂) ≤ 0.28583 (was 0.28528, below the certified value). Also in (U1) l.226 and in the table.
- l.346, **new paragraph** (M5 last bullet; B12 §4 item 5): (A)'s core |z₁ᴬ| ≤ 0.52 follows from E^A ≤ 2Ê by Sobolev: |z₁ᴬ| ≤ √2·S₁ᵛ(Ê) ≤ √2·0.10429 ≤ 0.1475 [cert; SA₄ `S1` at `EA_ref4` ≥ Ê]. Lemma U's own 0.5333 ≤ M₁ᶻ is unchanged.

### §7 Bootstrap
- l.364–366, §7.1 standing facts (B12 §4 item 3):
  - (P1) is now in L² form (Chou / Wheeler; to be confirmed by the authors);
  - v3.2's wording ("curvature and its first derivatives") is not supplied and not needed;
  - the graph property moves into the bootstrap as G (Lemma G of `continuation_v33.tex`).
- l.375: n^III ≤ (0.071445, 0.32707, 3.0445), computed with the certified c_δU (M3).
- l.377: N̄² ≤ (0.071472, 0.33012, 5.0403) ≤ 4·N2_ref. This domination is why (A)'s margins are quoted from SA₄ (M6).
- l.383: **B2′ restated** as sup|v_A| + sup_core|d| ≤ 0.0153 (problem P1 below).
- l.386: new bootstrap item G.
- l.388: an explanation of the B2′ change: the derived S₀ᵛ + sup|d| has margin ×2.116, and the chart-identification caveat is stated.
- l.390: the β_A paragraph.
  - Fluxes use (√K(2E^A) + √K(E_d))².
  - Margins ×3.363 / ×5.953 at 4·N2_ref [SA₄] (were ×3.505 / ×6.196).
  - The N2_ref values ×3.395 / ×6.010 are given as a comparison only.
- l.395–402, order table:
  - D1 consumes bounds on w′;
  - D3′ consumes (H2) and (H2f) in f − U₁, and d;
  - new row G.
- l.406–411: the closure rule 𝒜_j(β) is the **G-A1 (R-sup) rule**, with A⁰ (`A0_Rsup`), F_g (`Fg`) and coupling coefficient 1. The text notes that M2 restored the coupling term in assemble_U.
- l.420–422, closure proposition (B12 §4 item 6):
  - (K6) now includes G: (H4) on [48, 700] plus the verifier zone, and (H5);
  - the conclusion is BS_U including G on [τ₁, τ_max).
- l.424, **new pointer** (B12 §4 item 4): τ_max = ∞ by Proposition "global existence in rescaled time" of `continuation_v33.tex`, which also gives (W1)–(W5).
- l.425, certified values. All are M14 numbers.
  - Ê = (1.2021e-6, 1.0440e-4, 0.056686).
  - (A) rows at 4·N2_ref: S₀ᵛ 7.1881e-3, ×2.116, ×5.748, M1 ×4.409, M0 ×5.449, Ĝ₂₂ ×1.917, ramp M1 ×1.184. A note says that no script asserts the SA₄ rows (M6).
  - B̂ = (4.8735e-6, 5.5927e-4, 0.36670).
  - late M1 ×2.353, M0 ×4.892; present M1 ×2.770, M0 ×2.225.
  - M0z ≥ ×12.99 (late). The old "×13.04" was the present-cell value and above the late value.
  - **(K6) ×1.190.**
  - The S4 certification of (H2), (H2f) (×3.188, ×4.983, ×4.106) and the fixed-L data. The statement "pending task 6" is removed.
- l.427: the implication step derives G strictly.
- l.429: "Nonempty" cites the f − U₁ certification of (H2) and (H2f), and G.
- l.431: "Closed" covers G by Lemma G.
- l.439: loop gain ≤ (1.23e-7, 0.037169, 0.039035). The level-2 parts are 0.037018 and 0.0020172 (were 0.037017 and 0.0020171, not outward).
- l.444: U-fine B̂ is (4.8723e-6, 2.6266e-4, 0.025840); K₂ ≤ 0.046955 and K₃ ≤ 0.054139. The H6 margins are unchanged.

### §9–11 Table, hypotheses, gaps
- l.451–501: the table was rewritten with the v34 values and keys. New rows: (A)'s core |z₁ᴬ|, d, the budget slack and primary Y, K°_jj, (K6), and S4.
- l.503–514, "What remains a hypothesis":
  - (H2) and (H2f) are *certified* in S4, with values and margins;
  - the (H2f)₀ = 5e-7 label in SZ `handoff_thresholds` is stale (M9);
  - the step-4 Taylor-radius caveat (M9) is stated;
  - **(H4) on y ∈ [48, 700]** (B12 §4 item 6; review M7);
  - the fixed-L cone data are certified (slope ×2.641, h₀″ ×3.037, n₃ ×27.38, H4 ×1119, H1, H7 ×2.302);
  - **far-family coverage** (M9, inherited): centres in (490, 602) and above 602;
  - (P1) is in L² form;
  - (W2) is removed from the hypotheses.
- l.516–535, gaps:
  - old G1–G3 and G5–G7 are resolved, with reasons;
  - G8 records the M5 repair;
  - new open items O1–O7 (see below); O7, at l.535, covers the M9 rerun of S4.

## Problems found (stated plainly)

- **P1. The variable in B2′ (hypothesis level; fixed here, no margin moves).**
  - The old B2′ read "sup|v_A| ≤ 0.0153", with v_A = f − U.
  - After B0, the verifier (D1) consumes B2′ as a bound on w′, the verifier's S0V. From sup|v_A| ≤ 0.0153 one gets only sup|w′| ≤ 0.0153 + 4.25e-5, which does not meet S0V.
  - (A)'s *derived* value is already S₀ᵛ + sup_core|d| ≤ 0.0153/2.116. So the fix is to state B2′ as sup|v_A| + sup_core|d| ≤ 0.0153, which implies both uses. lemmaU_v34 does this.
  - **`lemmaS_A_v33.tex` (Lemma z1, boxes table) and `step3_v33_main.tex` l.66 still state B2′ in the old form and should be aligned.**
- **P2. Non-outward roundings in the old text.**
  - Ĝ₀₀ 51322.5 (true 51322.514).
  - The ΔR level-2 parts 0.011297 and 9.2e-6.
  - The outer Y add 0.108 (0.10802).
  - The loop-gain parts 0.037017 and 0.0020171.
  - All are corrected. None is a hypothesis-level quantity.
- **P3. M6 is only checked by my script.** No project script asserts "(A) margins at 4·N2_ref together with N̄² ≤ 4·N2_ref and E_ref4 ≥ Ê". assemble_U only records `N2_le_4N2ref`. `check_keys_v34.py` checks all three.
- **P4. Stale provenance pin.** `ramp_U.json` pins u1_ghat 4b86cc91… (pre-M4). The current file is 5a06d24c…, and only its args block differs.
- **P5. Stale record in SZ.** `zone3_interface.json handoff_thresholds` still lists (H2f)₀ = 5e-7 with a "probe/est" label. The certified 3.14e-7 passes both 5e-7 and 1e-6.
- **P6. The step-4 manifest changed during this task.**
  - At 08:19 today a running `final_tail_v33.py --from 9` (review M9) rewrote `FINAL_manifest_v33.json`: 19 stages, no longer 8202d2c6…. A `.partial.json` and the pass is still active.
  - The rows I cite (`results.H2f`, `H2_tightened_fU1`, `cone_fixedL`, `H4/H1/H7_fixedL`) are unchanged in the new file, and all quoted numbers pass against it.
  - The new arb rows `M9_H2f` and `M9_H2_tightened_fU1` agree with the quoted digits.
  - `M9_H4_far` reports far-family H4 at every centre ≥ 490 (min ×15.58). If it stands, it closes the far-family coverage caveat.
  - Recorded as gap O7 (l.535). **When M9 finishes:** set PIN in `make_constants_v34.py` to the new hash, rerun both scripts, update l.27, cite the M9 rows, and drop the caveat in §hyp.
- **Not re-checked, inherited** (listed as O5 and in §hyp):
  - the chart identification of w′ at levels ≥ 1;
  - far-family coverage;
  - (P1) for immersed curves;
  - the zone3 far-K₃ four-piece list.

## Stale numbers in `fix1/b12_text/numbers_A_v33.tex` (not edited)

Against the M14 sources:

| macro | in file | now | source |
|---|---|---|---|
| header hashes | A 04947b7f…, U 4c3f23bc… | A 34832f43…, U 271430d7… | M14 |
| `\EhV` | (1.1953e-6, 1.0433e-4, 0.056683) | (1.2021e-6, 1.0440e-4, 0.056686); the old value is no longer an upper bound | U supersolution.EA |
| `\EAzeroV` | (1.1913e-6, 1.0427e-4, 0.056160), comment says "EA0" | EA0 = (1.1980e-6, 1.0434e-4, 0.056163). The old numbers match the new `A0_Rsup`, so relabel or update | A EA0 / A0_Rsup |
| `\coupV` | (0, 1.1506e-12) | (0, 1.1571e-12) | A coup |
| `\CzeroV` | (1.1290e-6, 1.6616e-4) | (1.1345e-6, 1.6621e-4) | U supersolution.C0 |
| `\NbarV` | (0.071473, 0.33012, 5.0400) | (0.071472, 0.33012, 5.0403); level 2 is below the certified 5.040298 | U supersolution.N2 |
| `\SzerovV` | 7.177e-3 | 7.1881e-3 | A4 S0 |
| `\SonevV` | 0.10428 | 0.10429 (0.1042886). √2·S₁ᵛ = 0.14749 < 0.148 still holds | A4 S1 |
| `\MarFluxZeroV`, `\MarFluxTwoV` | 3.382, 5.983 | 3.363, 5.953 | A4 flux.beta.2 |
| `\MarSzeroV` | 2.119 | 2.116 | A4 margins |
| `\MarSoneV` | 5.749 | 5.748 | A4 margins |
| `\MarMzeroV` | 5.453 | 5.449 | A4 margins |
| `\MarRampMzeroV` | 7.021 | 7.020 (sI1 0.001161, M3) | A core.margins |
| `\MarFluxRefV` | 3.414, 6.041 | 3.395, 6.010 | A core.flux_margin.2 |
| `\BhV` | (4.8734e-6, 5.5898e-4, 0.36664) | (4.8735e-6, 5.5927e-4, 0.36670) | U supersolution.Bhat |
| `\UmarV` | late M1 ×2.354, late M0 ×4.895, present M1 ×2.771 | ×2.353, ×4.892, ×2.770 (present M0 ×2.225 unchanged) | U margins |
| `\BhFineV` | (4.8719e-6, 2.6265e-4, 0.025836) | (4.8723e-6, 2.6266e-4, 0.025840) | U ufine.Bhat |
| `\KsixV` | 1.183 | 1.190 | U EA₀ / 1.01e-6 |
| `\GAoneV` | (6.8e-9, 6.6e-8, 3.2e-6) [est] | G-A1 is now applied; measured (6.7e-9, 6.5e-8, 3.0e-6). The sentences in lemmaS_A_v33 §(a) and G-A1, and in step3_v33_main's status box, saying the values are "script values until rerun", are obsolete | M14 |

Unchanged and still valid:
- `\gpV`, `\GcoreV`, `\GcoreOffV`, `\kappaJV`, `\DsrcV`, `\rhoJV`, `\acV`, `\FtwoV`, `\ConeZeroV`, `\ConeV`, `\kapOtwotwoV`;
- `\MarMoneV` (4.409), `\MarGV`, `\MarRampMoneV`, `\MarRampZoneV`;
- the d tables, the datum macros and the verifier values.

Not caused by M14, but not outward:
- `\PhiV` = 3.319e-5. The certified Φ is 3.31955e-5 to 3.31967e-5 (A4 `info[j].Phi`), so it should read 3.3197e-5.
- `\dsupRampV` (9.07e-6, …) is the b0_dprof value. A's `B0.sup_ramp` (9.2716e-6, …) includes the 2e-7 chart allowance. This is consistent but worth one comment.

## Side conditions checked by `check_keys_v34.py`

| check | result |
|---|---|
| E_ref4 ≥ Ê componentwise; SA₄'s E = SA's EA_ref4 | holds |
| N̄² ≤ 4·N2_ref | holds; `N2_le_4N2ref` = true |
| √2·S₁ᵛ ≤ 0.1475 | 0.147486 |
| D^c ≤ D^src | holds |
| SA₄ and SA D^s agree | yes |
| (K6) margin | 1.19013 |
| B2′ margin (S₀ᵛ + sup_core\|d\| against 0.0153) | 2.11604 |
| U_cert pins of A, Z, u1 | match |
| S4 rows (H2f, H2_tightened_fU1, cone_fixedL, H4/H1/H7_fixedL) | all `ok` |
