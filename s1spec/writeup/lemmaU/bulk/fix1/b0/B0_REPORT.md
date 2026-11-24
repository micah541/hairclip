# Fix B0: the deviation variable of (A) versus the variable of the shells and the verifier

Nov 14, 2026. Finding B0 of `t8_review/REVIEW_v33.md`. Labels: [cert] = python-flint arb, outward rounding; [hyp] = hypothesis; [est] = estimate.

**Result.** Route (b): (A) stays in f − U throughout, and a certified bound on the profile difference is inserted at every interface where the two variables meet. The (H2) threshold is unchanged: the tightened (H2) = (1e-8, 1e-6, 1e-3) holds for f − U₁ on [14, 60] with margin ×27.3 [cert]. The chain was rerun t2 → t4 → zone3. K1–K3 hold, and no stop trigger was hit.

## 1. Which variable each component uses

Notation:
- w′ = f − αψ̂₃′(x) − α²ĥ₂′(x), with x = y − σ_tip. This is the v3.2 shell variable and the verifier variable.
- U(α, y) = 1.3056 α y^p, with p = 1.0246.
- d := (f − U) − w′ = α(ψ̂₃′ − 1.3056 y^p) + α²ĥ₂′.

| Component | Variable | Where it touches the other variable |
|---|---|---|
| Verifier (`step3nl_v33`) | w′ in the S1 chart; x-derivatives; W = (1+x²)⁻² | J_m and √Λ are for w′. Its far-field inputs FLUX0/FLUX2 (I3) and S0V/S1V are also for w′. |
| (A), energy identity | z^A_k = ∂_s^k(f − U(α(τ), y)) on [14, 700] | The forcing r₀ = 1.4(β + 1.02αā) is the forcing of f − U(α(τ), ·). The core box f = U + z^A and the outputs E^A, C0/C1 are also in this variable. |
| (A), ramp data [14, 28] | verifier J_m (w′) | **mismatch:** z^A = (w′-part) + D^k d at every τ |
| (A), datum (K6) | (H2) on [14, 60] (step 4: w′) + (H2f) on [60, 700] (f − U₁) | **mismatch** on [14, 60] at τ₁ |
| (A) → verifier | E^A used as w′ energies in K0, K2, I3; S0v, S1v used as sup\|w′\|, sup\|y w″\| | **mismatch** at every τ |
| Lemma U (`t4_cert`) | z = ∂_s^k(f − U₁) | Uses E^A, C0, C1 through δU only, which is consistent. Launch (a) uses (H2) on [28, 60] for f − U₁: **mismatch**, closed by stage `h2u_v33`. Launch (b), (K6) and Lemma dU are unchanged. |
| (H2f), (H3) windows (step 4 win mode) | f − U₁ | consistent |
| Step-4 shells (H2) | w′ (`win.py` shell mode, with a level-0 chart ball of ±2e-7) | see the datum row |

**Can (A) run in w′ on its core?** Not cheaply. Its forcing and its core box are derived for f − U(α(τ), y). (A) in w′ would need a new forcing (the profile residual of ψ̂₃ and ĥ₂) and a new core box, and the interface at y = 28–56 with Lemma U would still need d.

So (A) keeps f − U throughout. d is added in exactly four places:
1. the datum at τ₁ (step 4, new stage);
2. the ramp J conversion, at every τ;
3. the verifier fluxes;
4. S0V/S1V (and the ramp boxes).

Lemma U needs no change, because its (H2) is now certified for f − U₁ with the same threshold.

## 2. The certified d tables (new step-4 v33 stage `b0_dprof`)

Script: `step4/p2/v33/b0_dprof.py`. It uses the certified C⁴ profiles from `step3nl/agent_int/v3/c4fun.py` (the module that win.py and ift.py use) and σ_tip from `t9h/s1chart.py`.

Method:
- 2944 cells of width 1/64 on y ∈ [14, 60].
- The level-k terms of ψ̂₃′ − 1.3056 y^p are bounded in first-order form about the cell centre, to avoid the cancellation between the two terms.
- x-derivatives come from y^m ∂^m = D(D−1)…(D−m+1).

Results [cert]:
- Coefficients:
  - √E[g] = (3.97, 4.70, 5.57, 8.16)e-3, for g = ψ̂₃′ − 1.3056 y^p;
  - √E[h] = (706, 1446, 2964, 6074), for h = ĥ₂′.
  - These are energies ∫(28/y)(D^k ·)² ds over [14, 60]. The α² ĥ₂′ part dominates.
- At |α| ≤ α(t₁) = 1.17387e-4 (the maximum over the 32 widened slab-180 balls):
  - E_k[d] on [14, 60] = (1.035e-10, 4.19e-10, 1.72e-9, 7.13e-9).
- At |α| ≤ α₁ = 1.5e-4 [hyp, v3.3 bootstrap]:
  - E_k on [28, 60] = (2.47e-10, 1.01e-9, 4.16e-9, 1.73e-8);
  - sup on the ramp = (9.07e-6, 1.82e-5, 3.69e-5, 7.50e-5);
  - sup on the core = (4.22e-5, 8.56e-5, 1.74e-4, 3.56e-4);
  - ramp W-norms ∫W(∂^m d)² = (1.23e-15, 1.21e-17, 3.17e-19, 1.69e-20).

## 3. (H2) for f − U₁ (new step-4 v33 stage `h2u_v33`)

At t₁, U₁ = U(α(τ₁), ·). Minkowski in L²((28/y) ds), level by level, gives:

√G_k[f − U₁] ≤ √E_k[w′] + √E_k[d]

Here E_k[w′] is h2t_v33's with-error value (the v3.2 shells plus their error transfer at Ȳ = 1e-12, true curve, all of slab 180).

| k | E[w′] (with error) | E[d] | bound | threshold | margin |
|---|---|---|---|---|---|
| 0 | 8.014e-11 | 1.035e-10 | 3.658e-10 | 1e-8 | **×27.3** |
| 1 | 4.706e-10 | 4.186e-10 | 1.777e-9 | 1e-6 | ×562.8 |
| 2 | 1.5588e-6 | 1.720e-9 | 1.665e-6 | 1e-3 | ×600.9 |

All values are [cert]; N35: rounded outward from the arb rows (`h2u_v33_arb.json`, `h2t_v33_arb.json`), the earlier digits 3.657e-10, 1.719e-9, 1.664e-6, ×563, ×601 were not outward. **The threshold is unchanged.** The reviewer's guess of about 1e-7 was too high: on [14, 60], α²ĥ₂′ is at most 2.5e-5.

## 4. What changed in (A) (`fix1/b0/t2_A_v34`, a copy; t2_A_v33 is untouched)

- **`rigor/assembly_rep.py`** (backup `*_pre_b0.py`), new option `--dprof`:
  - Every quadratic ramp data term t splits as t ≤ (1+η)t_w + (1+1/η)t_d. η = √(t_d/t_w) per term, taken as a point value; this is the optimum, which equals exact Minkowski for the pointwise E^A.
  - Cubic terms (Φ: ∫∫εω|z₂|³, as used by the phi and ibp223 classes and by S33) split as (1+η)²t_w + (1+1/η)²t_d.
  - Cauchy–Schwarz products use the AM–GM window form.
  - t_d is a rate per unit τ, since |α| ≤ α₁ holds at every τ. In the pointwise E^A it contributes t_d/γ′_j. In the export windows it contributes t_d·T, which goes into C1 (C0 keeps the damped D).
  - Fluxes: K_i ≤ (√K_i(E) + √K_i(E_d))². This is valid because K_i = Σ k_l N_l with k_l ≥ 0 (asserted) and N_l linear in E. The same holds at β·E.
  - I3: η = 1/8 (I3 is not consumed downstream).
  - S0V/S1V and H-B2 checks use S0v + sup_core|d| and S1v + sup_core|Dd|.
  - Ramp boxes: z₁ uses 0.4328 + sup_ramp|Dd|. The typed 2e-5 is now asserted to be at least the certified 1.82e-5. M0 adds sup_ramp|d|.
  - The step-4 level-0 chart allowance of 2e-7 is included.
- **`run_t2.sh`** passes `--dprof`.
- **`ic_rep.py`** only has its text changed.
- **`export_t2.py`** writes a `B0` block.
- Two consecutive `run_t2.sh` runs gave byte-identical outputs (`out/SHA256SUMS`).

## 5. New certified numbers (v34 against v33; all [cert])

**(A), at the tightened (H2) and the (H2f) design setting**

| Quantity | v33 | v34 | Trigger |
|---|---|---|---|
| E^A_ref | (1.18766e-6, 1.03883e-4, 0.0565847) | (1.19506e-6, 1.04309e-4, 0.0566785) | |
| S0v vs 0.0153 (verifier: S0v + sup\|d\| = 7.199e-3) | ×2.1436 | **×2.1254** | ≥2.0, OK |
| S1v vs 0.6 | (not checked with d) | ×5.789 | |
| Core M1 / M0 | ×4.443 / ×5.464 | ×4.437 / ×5.459 | |
| Ramp M1 / z1 / M0 | ×1.185 / ×1.201 / ×7.030 | ×1.185 / ×1.201 / ×7.021 | |
| K0(E), K2(E) | 2.339e-21, 6.456e-21 | 2.422e-21, 6.674e-21 | |
| FLUX0/K0, FLUX2/K2 at β_A = 2 | ×3.505, ×6.196 | **×3.414, ×6.041** | ≥2.0, OK |
| FLUX0/(2·K0(E)) | ×3.505 | ×3.386 | ≥2.0, OK |
| C0 (k = 0, 1) | (1.1166e-6, 1.6544e-4) | (1.1287e-6, 1.6614e-4) | |
| C1_0 (k = 0, 1) | (3.89e-12, 4.91e-11) | (7.52e-10, 4.22e-8) | |

The C1_0 increase is the d-rate at α₁.

**Lemma U (`U_cert_v34.json`)**

| Quantity | v33 | v34 |
|---|---|---|
| Ê | (1.18788e-6, 1.03902e-4, 0.0565883) | (1.19529e-6, 1.04327e-4, 0.0566827) |
| B̂ | (4.8731e-6, 5.5667e-4, 0.36560) | (4.8733e-6, 5.5898e-4, 0.36664) |
| C1 (supersolution) | (8.07e-10, 1.216e-8) | (1.555e-9, 5.431e-8) |
| Slope share C1Σ_X/D_k | (0.0014, 0.00014) | (0.0027, 0.00064) |
| K1, K2, K3 | true | true |
| Late M1 / M0 / M0z | ×2.358 / ×4.900 / ×13.05 | **×2.354** / ×4.895 / ×13.01 |
| Present M1 / M0 | ×2.775 / ×2.225 | ×2.771 / ×2.225 |
| Loop gain | (1.22e-7, 0.037169, 0.039034) | (1.22e-7, 0.037169, 0.039034) |
| U-fine B̂ | (4.8716e-6, 2.6258e-4, 0.025780) | (4.8719e-6, 2.6265e-4, 0.025836) |
| (K6): (H2 + H2f)₀ against Ê₀ | ×1.176 | ×1.183 |

The late-centre M1 trigger is ≥1.5, so it is met.

(K6) now rests on (H2) for f − U₁ from `h2u_v33`. The flux rows are at E^A_ref (N2ref). Ê/E_ref − 1 = 2e-4, so FLUX0/K0(2Ê) ≈ ×3.413 [est from linearity].

**zone3 (rerun with `--ucert U_cert_for_zone3_v34.json`)**
- H6 minimum margin ×1.8609 (unchanged, far K3).
- Cone K3 ×2.511 (was ×2.512); cone slope ×2.463; t* ×2.302.
- U-fine K2 0.046953 and K3-strip 0.054134 (against 0.18).
- ΔR, ΔR_b, n_III and box_III are identical to v33, so U_cert_v34 (assembled against the v33 zone3 interface, now saved as `zone3_interface.v33.json`) is at its fixed point. A second zone3 run was byte-identical.

**Stop triggers: none hit.** FLUX0/K0 at β_A = 2 is ×3.41 (≥2.0). Late M1 is ×2.35 (≥1.5). S0v is ×2.125 against 0.0153 (≥2.0). No design change was needed beyond the route.

## 6. Records and hashes

**Closed step 4.** `p2/out/FINAL_manifest.json` (e01cb4da…) verifies: 523 paths, 97 upstream, 26 scripts, no mismatches. This was checked before, at the start and end of every driver pass, and after.

**v33 step-4 manifest.** `p2/out_v33/FINAL_manifest_v33.json` is now **12 stages**, 161 paths, 26 scripts, sha256 `8202d2c6c4962cae23a8051a3830d2908b601417b3c9c95041b95485884a1d86`.
- Two consecutive full passes are identical apart from `seconds`.
- It supersedes 52af64e6… (10 stages). No backup of that file was kept, but its stage records 0–8 carry over unchanged.
- Compared with 52af64e6:
  - `h2f/h2f_base.json` changed only in its `reused_balls` field (the balls were reused). The numbers are identical, which changes the args hash of `h2f_result.json`.
  - `notes_v33` was re-pinned to the new `zone3_interface.json`; the old pin 9611888c was already stale before this fix.

| File | sha256 |
|---|---|
| `step4/p2/v33/b0_dprof.py` | 26a220698d5f3d77831639c611bb6601f84ccc14672075f5dbfcaf5ef16379fa |
| `step4/p2/v33/h2u_v33.py` | 9fc567361d4009f661636992cfe905451f70f7cf17f034037985dfdc6179fd6b |
| `step4/p2/v33/final_tail_v33.py` (backup `final_tail_v33_pre_b0.py`, afe914d0…) | b58e9b75c430feb1dc1d6fc5ddb9f4595c44fe116cbb17f1745d9d1bcc40d378 |
| `step4/p2/out_v33/b0_dprof.json` | 7fdf378c5f59f098263aec8f39bfff7f365e0d25cbbcd154c3f1010bd4a2cd8e |
| `step4/p2/out_v33/h2u_v33.json` | 009c9d1197f0cd6d7bf79ae1cf6cead0f444d74c5bcb096546f8c6b374d36c2a |
| `fix1/b0/t2_A_v34/rigor/assembly_rep.py` | 196f8767ef8f0945d6fa40fc1f7facbe32ed5433366c56aafb99bcf0e0789bd5 |
| `fix1/b0/t2_A_v34/out/A_interface.json` | 04947b7fbc79fe28a67de10cb5a3f8bf249bc68b4cbd29b8c1a574a9cf16313e |
| `t4_cert/out/U_cert_v34.json` | 4c3f23bcd8761bcbb2cc72a1aa0787714d4201f36aa3613817e043be5a5e7668 |
| `t4_cert/out/U_cert_for_zone3_v34.json` | ce0798b05eb6b274409a0342660c1c4242085f809c6d5934b5166cf0f478f219 |
| `t5_zone3/out/zone3_interface.json` (new) | 8152c65076c649365edf0997cb0b92b12d3e4ab49c6c0ce6821fe7b09e51f195 |
| `t5_zone3/out/zone3_interface.v33.json` (saved copy) | 0d6338d0f8580cd0cc3c2b0f829b412f8ca0be37f077173023fe447fe382c754 |

Other record changes:
- `t5_zone3/out/SHA256SUMS` was regenerated; the old file is kept as `SHA256SUMS.v33`, and `zone3.v33.log` is kept.
- t4_cert's `SHA256SUMS` is unchanged and still verifies. The new files are `U_cert_v34.json`, `U_cert_for_zone3_v34.json` and `assemble_v34.log`.

## 7. Open points (not B0, stated plainly)

- **Chart identification.** The step-4 allowance of 2e-7 is included at level 0 only. Two further identifications are inherited from v3.1/v3.2, and B0 does not address them:
  - the verifier's S1-chart w′ (a normal graph over S1, in x) with the ray-graph quantity f − αψ̂₃′ − α²ĥ₂′;
  - the chart terms at levels ≥ 1 between them.
  The reviewer should look at this.
- **d at every τ uses α₁ uniformly.** This is valid but crude: α decays to the floor of about 7e-7. It is the reason C1_0 grew.
- **Typed constants (review M3).** The typed 0.00116 is unchanged. The typed 2e-5 is replaced in the margins by the certified sup|Dd| = 1.82e-5.
- **M2 is not fixed here.** assemble_U still drops the coupling term.

## 8. Documents whose quoted numbers must be updated

- **`t3_proof/lemmaU.tex`**:
  - l.60: the remark that (H2) "already stated for f − U₁" is now true via `h2u_v33` (×27.3);
  - add d to (A)'s ramp data, the fluxes and S0V/S1V;
  - (K1)–(K6) values: Ê, B̂, S0v ×2.125, core M1 ×4.437, M0 ×5.459, late ×2.354 / ×4.895, present ×2.771 / ×2.225, M0z, K6 ×1.18, loop gain, U-fine B̂, U-fine K2/K3, cone K3 ×2.511;
  - C0/C1 in Lemma D;
  - the U-fine proposition.
- **`t1_statements/STATEMENTS_r2.tex`**: Remark IC, the datum of (A), and the (H2) variable.
- **`t7_assembly/main_v3.tex`** (and `main_v3_CHANGES.md`): (A) and U numbers, and the (H2) variable sentence.
- **`t7_assembly/V33_INDEX.md`**, **`STEP4_RESULT_v33_addendum.md`**, **`step4/p2/v33/T6_REPORT.md`**: the v33 manifest is now 12 stages, sha 8202d2c6…, 161 paths, 26 scripts, with the new stages `b0_dprof` and `h2u_v33`.
- **`t7_assembly/erratum/ERRATUM.md`**, **`BULK.md`**, **`KITS_TODO.md`**: kit B/v33 must ship `b0_dprof.py` and `h2u_v33.py`, the c4fun/psi3fun/h2fun inputs, and t2_A_v34.
- **The proof text in `fix1/b12_text/`** (another agent): it must cite the v34 interface (`fix1/b0/t2_A_v34/out/A_interface.json`) and `U_cert_v34.json`, and it must describe the d-terms.
- **`t4_cert/run_all.sh`** default, and any citation of `t2_A_v33/out/A_interface.json` or `U_cert_v33.json` as current.
