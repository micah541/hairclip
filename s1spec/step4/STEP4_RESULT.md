# STEP 4 RESULT (base v2): certified chain Γ0 → switch → P2 → hand-off on slab 180

Status (Nov 12). **Theorem 4 is not yet proved.** A1–A3 and B8 are closed. B1 and B5 are closed (Nov 12). The final tail runs on B1's direct residuals and reads one version of every input (`p2/out/FINAL_manifest.json`), including T7's constants republished on the direct residuals; a second pass reproduces it. Lemma O H6 is closed in step 3 v3.2 by the covering argument plus Lemma R (proved, certified; `writeup/lemmaO_v32/lemmaR.pdf`), except for the strip input. That input belongs to the step-3 gap **Hypothesis (F)**: Lemma S(B) drops the window-edge terms (`writeup/step3_gap_v32/gap_v32.pdf`, `writeup/step3_gap_v32/REFEREE_check.md`). The main theorem is conditional on (F) (`writeup/main/main.pdf`). Step 4's own results do not depend on (F). Other open items: §3. Everything below is arb with outward rounding and no sampling, unless it is marked otherwise. Every quoted bound is rounded outward.

The hand-off good-time window is now **all of slab 180** (t_c ∈ [−1, 1], Δτ = 0.15, k ∈ [2.27e4, 2.38e4]). The old W′ (1.46e-4 long) failed Lemma O H4 by ×5–×280 (A1). All hand-off items are certified uniformly for τ* anywhere in the slab. The thinnest certified margins are:
- the far box |e| (√F₀ ratio against fac 1.7), ×1.69, and the level-3 far box x²|e‴|, ×1.99 (at Ȳ = 1e-12, against the N3′ boxes certified at Y = 1.7e-12; `n3boxes`: 1.69999995 and 1.9940, so ×1.69 outward, not ×1.70);
- |α|, ×1.27.

## 1. Chain
| link | result | margin |
|---|---|---|
| Γ0 (T2) | 3762 sine modes; b₂ by Krawczyk (half-width 8.1e-24); embedded | — |
| P1 (T3), slabs 0–96 | t₉₆ = 4.597850539485830831 ± 7.8e-19; k ∈ [29.352874545325, 29.352874545328] (`t3/p1/status.json` kmax_enclosure [29.352874545325992, 29.352874545327012], outward) | — |
| base v2 + (SE) | slabs 0–180, (SE) exact up to \|a₂\| ≤ 1.6e-14; T6 hulls and ranges widened outward on read; residual summary regenerated. **B1:** T6's engine on the (SE)-corrected base itself (`p2/t6c`: residuals.json, jets), the one T6 data set of the final tail | — |
| switch (A), base slab 2 | t_c,a = 982986695892350/2⁵⁰, so τ_a lies in slab 2 and P2 starts there (slabs 0–1 are unused). **Y(τ_a) ≤ 6.21e-13**; F_j(τ_a) ≤ 3.75e-23, 2.45e-22, 4.72e-19, 3.19e-15 (rerun on `p2/t6c`: unchanged; the typed Y0 = 6.21e-13 and F0 hold) | ×1.6 |
| T7 | Nonlinear closure (Lemma 4c, A2) certified on all 25 intervals, min ϱ_c ≥ 0.04711; DN = D4CUB = 0. Tip box (Lemma 4b, A3): sup\|e\| ≤ 1.6031e4·Λ_e^½ and sup\|e′\| ≤ 1.6352e4·Λ_e^½ on [0, 34.25]. Republished Nov 12 with se7 on `p2/t6c`: min ϱ_c 0.04711, (Y) 3.21e-14, Y_cert 6.21e-13, x_NL ≤ 1.95e-15 against its slack. Slabs 0–1 keep the se_corr widths (flagged; outside P2) and no final stage reads them | — |
| (Y) | level ≤ 3.22e-14 on the direct residuals (was 1.33e-13); ycomp asserts Ȳ = 1e-12 ≥ Y(τ_a) | — |
| (F), T8, window = slab 180 | N3′ and 𝒬 per slab (`n3far`, b₀, b₁, b₂ = 1.5e-3, 3e-4, 1e-2, fac 1.7, Y = 1.7e-12) on `constants_b1`, `n3_b1` (`p2/t6c`); `camp --n3 --window 180-180`, with the slab-wide base windows. B_j = 7.61e-7, 1.76e-9, 5.11e-8, 7.82e-6. D₂(τ*) ≤ 9.40e-12. Ȳ limits (computed KCH, c0, c0tip; B5): **item 3 E3 1.575e-11**; item 4 2.89e-11; item 2 (own bound, Λ_e route) E2 3.63e-10 | ×15.7 |
| boxes at Ȳ = 1e-12 (`n3boxes`, `boxes.py`) | \|e′\|, x\|e″\|, x²\|e‴\| against the N3′ boxes, and \|e\| (√F₀ ratio against fac 1.7). All four hold at Y = 1.7e-12 (2.38e-4, 1.08e-4, 8.34e-3, 1.00001). Tip \|e\|, \|e′\| ≤ 1.61e-8, 1.64e-8 against 1e-5 | ×6.42, ×4.69, ×1.99, **×1.69** (6.4287, 4.6917, 1.9940, 1.69999995) |

## 2. Hand-off items, uniform for τ* ∈ slab 180 (`p2/t9h`, `p2/t9.py`, `p2/t9_far.py`)
- **S1 chart.** prof_tip (conditioned, nested in prof_tight) with the conditioned τ_tip, plus ext813 and tail envelopes; the tip offset enters as one scalar. Θ = 1.8617364190148 ± 5.2e-11.
- **H1 (α, √Λ).**
  - 16 exact points t_k = −1 + (2k+1)/16. At each, a Krawczyk test (review B3) proves that the three constraints have exactly one zero, and Λ is evaluated there with the exact u⁗ cell hull (B5): √Λ(t_k) ≤ 2.01e-9.
  - 32 balls of radius 1/32 (`h1dt.py`) bound sup‖∂_t w‖_Λ ≤ 3.33e-7. On each ball, dw/dt = u_t + Σ ṗ_i m_i, with exact first variations and ṗ from the IFT by defect correction. The parameter box is verified a priori.
  - So √Λ(t) ≤ √Λ(t_k) + |t − t_k| sup‖w_t‖ ≤ 2.29e-8 on the whole slab. This includes C_ΛȲ ≤ 2.42e-11 (C_Λ ≤ 24.2, with T7's tip constant).
- **H2, H3.** Windows and shells use 32 t_c balls of radius 1/32, with the parameter paths p(t_k) + (t − t_k)ṗ.

| item | base (arb, slab 180) | with error (Ȳ = 1e-12) | margin |
|---|---|---|---|
| (H1) \|α\| | ≤ 1.1739e-4 | + C_αȲ ≤ 7.3e-13 | ×1.27 |
| (H1) √Λ | ≤ 2.282e-8 | ≤ 2.29e-8; holds for Ȳ ≤ 1.12e-9 (clam's earlier 1.58e-9 took the wrong end of the C_Λ ball) | ×2.18 |
| (H2) shells [14, 60] | 8.01e-11, 1.34e-10, 2.94e-10 | Ȳ ≤ 3.63e-10 (E2; own bound by Λ_e) | ×363 |
| (H3) windows E0–E3 | 3.48e-8, 1.27e-7, 3.76e-7, 9.71e-7 (against 2e-7, 1e-5, 1e-5, 1e-5) | Ȳ ≤ 1.575e-11 (E3) | ×15.7 |
| (H3) K4, y ∈ [50, 350] | 3.1633e-4 | 5.60e-3 against 0.06 | ×10.7 |
| (H4) graph over the ray | \|θ − Θ\| ≤ 0.078 | + 5.5e-4 | — |
| (H5) cone, r_c ∈ [200, 380] | exact curvature jets (`curvjets.py`) | K2, K3, K4 ≤ 1.89e-2, 4.78e-4, 1.23e-2 against 0.12; N3 ≤ 1.17e-3 against 0.14; slope ≤ 0.0573 | ×6.3, ×251, ×9.7, ×119 |
| (H5) far cylinders, x ≥ 150 | slab-180 hulls, exact jets | \|h″\|, \|h‴\|, \|h⁗\| ≤ 1.80e-4, 2.60e-8, 6.97e-9; slope ≤ 0.0416; ‖h‴‖_L² ≤ 5.52e-7 | ×48, ×2698, ×89, ×5.0, ×1509 |
| (H5) Lemma O H4, ∫(h₀⁽⁵⁾)² | base ≤ 6.8e-6 | worst cone r_c 370–380 ≤ 0.2393 M²; far ≤ 0.00397 M² (`h4loc`; chart factor computed, ≤ 1.0103, cleanup O6b) | ×4.18 |

## 3. Gaps
- **Open:**
  - **Lemma O H6.** Lemma O H6 is closed in step 3 v3.2 by the covering argument plus Lemma R (proved, certified; `writeup/lemmaO_v32/lemmaR.pdf`), except for the strip input. That input belongs to the step-3 gap **Hypothesis (F)**: Lemma S(B) drops the window-edge terms (`writeup/step3_gap_v32/gap_v32.pdf`, `writeup/step3_gap_v32/REFEREE_check.md`). The main theorem is conditional on (F) (`writeup/main/main.pdf`). Step 4's own results do not depend on (F).
  - *For the record, v3.1:* Lemma O H6 fails as specified there. The amendment plan was `p2/h6/H6_PLAN.md`: cone K = 0.25, far (M, K) = (0.1, 0.2), and the |x| ≤ 1 covering. Its thinnest margins are ×1.13 (cone K₃) and ×1.14 (cone K₂). The details below are from `t9_h6.py`; step 3 is not changed here.
    - Cone family by overlapping centres: at a cylinder's inner end L_c/L_p = 4, so the recovered bound is 4·(0.0189 + 0.0395) = 0.234 against K₂ = 0.18 (×0.77). It closes only for L_c/L_p ≤ 3.08.
    - Cone cylinders beyond the last cone centre (r > 380) take their bound from the far family: (95/113)·0.374 = 0.315 against 0.18 (×0.57).
    - Both far checks pass: far family by overlap 0.374 against 1.5 (×4.0); far inner strip from the cone family 0.176 against 1.5 (×8.5).
    - Minimal fixes in step 3, each with its margin:
      - (a) Re-certify the cone Lemma O at K₂ = K₃ ≥ 0.25 (×1.07); t* must stay ≥ 1.6e-7 (it is 3.70e-7 at 0.18).
      - (b) Re-certify the far family with C₂√t ≤ 0.193 (now 0.353), e.g. at smaller M and K, since the far data are ≤ 0.021.
      - (c) Extending the cone family instead fails Lemma O H4 beyond r ≈ 482.
  - **H7 and H9** are internal to step 3: the time ratio holds (cone 1.6e-7 vs 3.70e-7, far 6.13e-9 vs 1.23e-8), and regularity holds by the parabolic flow and Lemma O's conclusion.
  - **Theorem 4(a):** the k_max of the true curve needs a pointwise tip bound on e″.
- **Closed:**
  - **Lemma O H1, arc by arc (cleanup O1, Nov 12).** H8 is not required.
    - Cone family: sup|h₀| ≤ 0.1872, H ≤ 0.1895 in all 18 cylinders (`t9_lemmaO.py --v32`).
    - Far family (v3.2: X = 2, L = 113, (M, K) = (0.15, 0.3)): `p2/t9_lemmaO_far.py` → `p2/out/t9_lemmaO_far.json`, a stage of the final tail, **ok = true**. For every centre on the arm (and its mirror) from the first far centre to the far end, on base arclength ±G = ±227.125 around the centre the true tangent turns by ≤ K₀·J_max·G ≤ 0.04095 rad (K₀ ≤ h₂ = 1.797e-4 from `t9_far_v32.json`; J = |γ′| ∈ [0.99667, 1.00333] from k₀, |e|, |e′|). So ξ is increasing and |ξ(±G)| ≥ 226.18 > XL = 226: the arc through the centre is a graph over [−X, X] in its own tangent frame. sup|h₀| ≤ (X²/2)·L·K₀·sec³ ≤ 0.04072; with C₀√t₁ ≤ 9.212e-4 (`lemmaO_v32_consts.json`, far_v32), H = 0.0417 and sup|h₀| ≤ H − C₀√t₁ (slack 6.7e-5).
    - The arc stays where t9_far's bounds hold (x ≥ 148) for every centre at base arclength ≥ 375.125. That covers "centres ≥ 376" in hand-off units read as arclength (≥ 375.998 in base units) or as radius (points with x < 375.125 have r ≤ 373.013). The frame scale λ_f ≤ 1 makes the base-unit computation conservative.
    - The old script's separation test (H8's content; it failed at the self-crossing) is dropped (backup `p2/t9_lemmaO_far_pre_cleanup.py`).
    - *Informational:* the hand-off curve self-intersects. `p2/t9_axis_crossing.py` → `p2/out/t9_axis_crossing.json` is now a stage of the final tail (the earlier JSON had no script on disk). On each of 16 t_c-balls covering slab 180, Im of the true arm changes sign on element 93, s ∈ [29/32, 30/32], x ∈ [3.269e4, 3.398e4], Re ∈ [2.235e4, 2.316e4], with |e| ≤ 0.525 (the earlier JSON had 1.24, from older data) and certified |Im| ≥ |e| + 7.9 at both points. This is the expected loss of embeddedness (STEP4_DESIGN §2). Step 3 applies Lemma O arc by arc and needs only H1–H7 and H9 (step3_proof l.97, l.105, l.111), so H8 is not a requirement.
  - A1, A2, A3; B2, B3, B4, B6, B7, B9.
  - B8 (a stage of the final tail since Nov 12): `/home/claude/s1_tight/consistency/check_outA5.py` shows all 1025 + 1025 + 501 rows of prof_tip, prof_tight and ext813 lie inside the outA5/ibp_boxes hulls, with positive margins. Theorem 4 is stated for an S1 in the certified set.
  - The writer's rounding and stale-summary items.
  - **B1 (T8, Nov 12):** the residual of the (SE)-corrected base is computed directly by T6's engine (`p2/t6c`; T8_REPORT, LEMMA_F §6.1); se_corr is not used.
  - **Cleanup of the proof's script items O6a–c (Nov 12).** Backups `*_pre_cleanup.py`; manifest before: `p2/out/FINAL_manifest_pre_cleanup.json`.
    - O6a, `t1clock.py`: dt = λ⁴(1 − ϑ)dτ with the true λ. Per slab, |ln(λ/λ̃)| ≤ ρ_s = ρ_a + Σ Δτ(x̄_s + θ_s c̄_s), since d ln(λ/λ̃)/dτ = a − ã = x₁ + ϑc. Here ρ_a = μ₁ + s·c̄₂ ≤ 4.14e-13 from the switch's IFT (`switch_s2.json`), x̄_s = C_xȲ + x_r,s + x_NL with x_r,s computed per slab as in n3far (max 2.57e-16; the typed 2.0e-16 was below it) and x_NL ≤ 1.95e-15, and θ_s = x̄_s/n_V,s. The clock factor is now per slab, plus a ±s·λ⁴ start term ≤ 1.1e-17. Result: t₁ ∈ [4.597857477538968, 4.597858614388128] (was [4.597857477538703, 4.597858614388391]; the per-slab θ_s tightens more than λ/λ̃ widens), ⊂ [4.5978574775, 4.5978586144] as quoted. |ln(λ/λ̃)| ≤ 3.38e-8 on slab 180 (the hand bound was 2.32e-7), and θ_max = 3.3641e-8 ≤ 3.37e-8. The base k_max hull is unchanged.
    - O6b, `h4loc.py`: the typed chart factor 1.001² is replaced by a computed one. The true curve's h⁽⁵⁾ over the cylinder's own centre tangent is expanded in arb Taylor mode (as `t9chart.py`/KCH) at the same base point. The e⁽⁵⁾ coefficient is A5 ≤ 1.0103, and the lower-order terms e…e⁗ are charged explicitly (≤ 5.4e-12 pointwise in h⁽⁵⁾). The weight is dξ = M_t dx, and the angle includes the true-vs-base centre tangent. The hand-off frame (λ_f ≤ 1, |x_c|) and the arc's excess over XL fit the supports' slack 1 (asserted, ≤ 0.43). Chart factor A5·√M_t ≤ 1.0103 (cone r_c 370–380) and 1.0067 (far). H4: worst cone 0.2297 → 0.2393 M² (×4.35 → ×4.18); far at M = 0.15 0.1712 → 0.1761 M² (×5.84 → ×5.68); far at M = 1 0.00386 → 0.00397 M².
    - O6c: the far-H1 stage and the axis-crossing stage are in `final_tail.py` and `FINAL_manifest.json` (now stages 30, 31). h4loc's stages now also pin `t9_cone_kch.json` and `t9chart.py`, and t1clock's pin `t7/constants.json` and `switch_s2.json`. `final_tail.py --from 'h4loc on the final amp'` was run twice; stages 0–19 were reused (442 paths verified). Pass 2 reproduces all 11 non-log outputs of stages 20–30 byte for byte. The only differences are the stages' seconds and basegeo's build-time lines in three logs (`checks.cleanup_second_pass`).
    - O6(c), `t9h/h1.py` rotation-mode term (follow-up, Nov 12). The E_I part of the mode norms on [10.25, 60] is now summed by cells, replacing the typed |g′| ≤ 1.01, 1e-3, 1e-3, 1.6x^1.03 and the factor 50. It uses the S1 rows (with dT) for σ ≤ 13, the step-1 tail envelopes beyond, and C4 jets for ψ̂₃ + 2αĥ₂. The typed |sc′| ≤ 1e-3 was false near x = 10.25 (sup 1.204e-3), though its integral still over-covered. dp now uses the Jacobian over the Krawczyk box, which contains the C-ball path (asserted). The 1.01 allowance on the C-term is unchanged.
      - Moved: rotation-mode norm max 17.642 → 17.626 (al-mode norm 7.7414 → 7.7216), √Λ(t_k) by −2.8e-14, base sup 2.2812816e-8 → 2.2812788e-8, C_Λ 24.175 → 24.161, with error 2.283699e-8 → 2.283695e-8. The quoted ≤ 2.282e-8, ≤ 2.29e-8 and C_Λ ≤ 24.2 stand.
      - |α| did not move: alpha_abs_max, P_hull and C_α are unchanged, and the Krawczyk/IFT and Λ fields of all 16 points are identical.
      - `t9h/clam.py` took the upper end of (5e-8 − base)/C_Λ for the H1 Ȳ limit, i.e. divided by the lower end of the wide C_Λ ball. Fixed (backup `clam_pre_cleanup.py`): 1.586e-9 → 1.125e-9.
      - The 16 h1 points (`t9h/h1_points.py`, two in parallel) and the H1 summary (`t9h/h1_summary.py`; no script for it was on disk, and the new one reproduces the old `summary.json` byte for byte) are now stages 4–5 of the final tail. `final_tail.py --from 'H1 exact points'` was run twice. Pass 2 was killed by an environment restart after stage 12 and continued with `--from 13`. It reproduces pass 1, apart from stage seconds, the h1 outputs' `cpu` field and timing lines in logs. Every output after clam is byte-identical to before the follow-up (t9_cone_kch.json included).
    - Restart-safe tail (Nov 12). `final_tail.py` writes `p2/out/FINAL_manifest.json` only after every stage has finished and every recorded SHA-256 has been checked against the files, atomically (temp file + rename). Progress after each stage goes to `p2/out/FINAL_manifest.partial.json`, and `--from N` continues a stopped or killed pass from it; `--until N` stops a pass. The H1 stage resumes point by point (`t9h/h1_points.py`; marker `out/h1slab/h1_points.pass`; `H1ONLY` for chunks).
    - B8 as a stage (Nov 12). The packaging agent's smoke test reran `check_outA5.py` in place at 17:09, which changed only the `seconds` field of `check_outA5.json`. The check is now the last stage of the final tail and is pinned normally (no stage reads it). Its results are unchanged: all 1025 + 1025 + 501 rows inside.
    - FINAL_manifest (Nov 12, final): 34 stages, the 31 of the first cleanup plus the H1 points, the H1 summary and B8. Every pinned hash matches: 523 paths, 26 scripts and 97 upstream entries (h1dt balls and windows). `checks.earlier_reruns` records how stages 4–32 were produced.
  - **Final tail, one version of every input (Nov 12).** `p2/final_tail.py` writes `p2/out/FINAL_manifest.json`: every stage in order, with cmd, cwd, env and the SHA-256 of its inputs and outputs. Data sets: T6 = `p2/t6c` (residuals and jets), T7 = `t7/constants.json`, T8 = `constants_b1`/`n3_b1`, base windows = `h1slab/win_summary.json`.
    - Settled versions: switch, ycomp, camp, n3far, h4loc, t9, t9_far, Lemma O and t1clock read `p2/t6c` only (`switch.py`, `t9.py`, `t9_far.py`, `t9_lemmaO.py` and `t1clock.py` now take T6DIR). clam takes μ_e from constants_b1 (T8CONSTS). No final stage reads `p2/t6/*`, `step4/t6/residuals.json`, se_corr.json, constants_fin.json or n3.json.
    - T8's flag (a): every camp run records T9WIN and reads `h1slab/win_summary.json` and the final `t9_cone_kch.json`. Both passes are in the manifest, and t9 reproduces bit for bit between them.
    - T8's flag (b): `t6c/residuals.json` is the version that constants_b1 records (SHA-256 checked).
    - Rebuilt on t6c:
      - switch: Y(τ_a) ≤ 6.2032e-13 and F_j unchanged, so the typed Y0 and F0 hold;
      - ycomp: level 3.22e-14;
      - clam: identical;
      - T8's camp/n3far iteration with B1's args: identical boxes;
      - then the window camp, t9, t9_far, h4loc, Lemma O, H6, the v3.2 scripts and t1clock.
    - Kept from A1: the h1dt derivative balls and the windows (`t9h`); since the follow-up the 16 h1 exact points and the H1 summary are stages. They read only the base polynomials, which equal the (SE)-corrected base for x < 0.6ℓ̃. Their hashes are in the manifest's upstream list.
    - Changes against the A1 data: Ȳ ≤ 1.575e-11 (E3), cone K4 ×9.7, H3 K4 ×10.7, H4 ×4.35 (×4.18 after cleanup O6b), far h₄ ×89. KCH, c0, c0tip and the shells are unchanged to 6 digits; t₁ is unchanged.
    - T7 republished (Nov 12; se7 on `p2/t6c`). Rerun: ycomp, cconst/merge, then every stage that reads `t7/constants.json` (switch, the T8 iteration, camp, n3boxes, t9) and the rest of `final_tail.py`. constants_b1 changes only in c_I, μ₂ and C_x's source text. The final numbers move in the 4th–5th digit (E3 limit 1.5755e-11 → 1.5757e-11, cone K4 1.2297e-2 → 1.2296e-2), so the rounded values quoted here stand.
    - T7 slabs 0–1 (se_corr widths): no final stage reads them.
      - cconst and merge take T7's worst cases over slabs ≥ 2 only (new env T7SLAB0 = 2, recorded in their provenance). Before this, cconst's min/max ran over all slabs (conservative).
      - ycomp builds its farg table for m ≥ SLAB0 only (edited; rerun byte-identical). camp, n3far and the switch index slabs ≥ 2.
      - n3boxes and t9 read only T7's meta; h4loc, t9_far, Lemma O, the axis crossing and clam do not read `t7/constants.json`. t1clock (cleanup) reads its slabs ≥ 2 only (C3 m̄ for x_r), the certificates' c_φ and meta.modulation_nl (x_NL).
    - Second pass (cconst's slab files moved aside, so every stage recomputes): all 199 non-log outputs reproduce byte for byte, except the cpu/wall timing fields. Those occur in the 179 cconst slab files, constants_b1.json and switch_s2.json, and the three are identical once the fields are removed (checked). The pass-1 manifest is `p2/out/FINAL_manifest_pass1.json`.
  - **B5, remainder (Nov 12).** `p2/t9.py` computes KCH, c0 and c0tip with no typed allowance or factor (`p2/t9chart.py`: arb Taylor mode, no sympy at run time, 8 s). The `u4_cell` hull is exact (closed earlier).
    - Exact split: D^j f_t(y) − D^j f_b(y) = [Q^t_j − Q^b_j](X_t(y)) + [Q^b_j(X_t(y)) − Q^b_j(X_b(y))], Q_j = D^j of the slope at one base point x; true curve F_t = (a sinφ + e′cosφ)/(a cosφ − e′sinφ), a = 1 − κ̃e, Y_t = Y_b − e sinφ, D_t = (Y_t/Y_t′)d/dx.
    - First bracket: exact linear part plus a remainder charged to η_i = y₀^{i−1}e^{(i)} with the n3 boxes (\|e′\|, x\|e″\|, x²\|e‴\| ≤ 1.5e-3, 3e-4, 1e-2; \|e\| ≤ c_eȲ + √F₀(y³/84)^½ at Y = 1.7e-12, c_e = 1.6031e4 from T7 Lemma 4b). Second bracket ≤ sup\|dQ^b_j/dx\|·\|e\|·tan φ*. e-terms per window bin with sqrt(28/(3y_w)).
    - Check (`p2/t9chart_check.py`): the linear part equals D^j(sec²φ e′ + sec²φ tanφ κ̃ e) (sympy, exact agreement); the remainder bound holds on 200 random jets (worst 0.99 of it).
    - KCH 1.1512 → **1.1878** (linear 1.1836 × computed nonlinear factor 1.0036, against the typed 1.0001). Most of the rise is the e-term's own w₁, w₂ terms at levels 2–3 (15κ̃y²tanφ sec⁴φ and 10κ̃²y³ sec⁵φ in the level-3 w₁ coefficient, against 12 and 7 without it), which the c0 allowance did not cover.
    - c0 = 2.69e-3, 2.20e-2, 4.71e-2, 9.87e-2 (was 0.002, 2.27e-2, 4.43e-2, 6.30e-2; the level-0 0.002 was too small); c0tip = 0.111, 0.114, 0.118, 0.125 per unit Ȳ (new; `t9_cone_kch.json` ≤ 0.11040, 0.11321, 0.11707, 0.12423).
    - Effect (camp, window 180): on the final B1 data (`amp_A1_w180.json`) Ȳ ≤ 1.575e-11 (E3, binding), ×15.7; item 4 2.89e-11. (On the A1 data at the time: 1.517e-11, ×15.1, and 2.77e-11; before B5, 1.567e-11 and 2.87e-11.) b₀ for the cone jets now comes from `n3.json` (1.5e-3; it was `n3_fin.json`, 7e-4, outside the A1 chain): K3 4.90e-4 → 4.92e-4, N3 1.200e-3 → 1.204e-3. α, √Λ, H4 and t₁ are unchanged.
    - y_ray/y_T8 is now certified by a tip anchor: y_ray(x) = p(0)cosΘ + x − ∫₀ˣ(1 − cos φ), with the integral bounded on the arrangement of the cells' x-hulls, and q(0) = 0, θ(0) = −π/2 by symmetry (asserted against the tip hulls). The former enclosure q − x ∈ [q_lo − a, q_hi − b] assumed q − x static over the slab. The slab-wide hulls do not imply that, and q − x drifts by about 0.013 per slab at x = 380. KCH is 0.03% lower with the anchor (`T9RATIO=qx` keeps the old enclosure).
  - **Shells [14, 60], own bound (Nov 12; B5 follow-up).** Item 2 no longer reuses the windows' KCH. T8's F_m cannot carry the shells: their ramp weights vanish at y_T8 = 14, and y_h − y_T8 is not certified to the 1e-4 that would be needed there. `t9.py` (block `shell`) bounds the shells by Λ_e alone, per unit Ȳ, and camp's item 2 reads it.
    - T5/T7 norms: v = ζ_∞e′ with ζ_∞ = S((x − 10.25)/6) on [10.25, 16.25]; V_m = ‖v^{(m)}‖_W, W = (1 + x²)^{-2}; Λ_e ≥ μ_e(V₀² + ε_IV₂²); V₁ ≤ k_{V1}V₀ + V₂/2. So V₀, V₁, V₂ ≤ 5.28, 486, 966 per unit Ȳ, including the ζ corrections on the sliver x ∈ [15.94, 16.25] (ζ ≥ 0.99977).
    - Same exact split and Taylor mode as the windows, weight 28/y_h² dy_h with the cut-off (≤ 1) dropped, as for the base in `win.py`. |e| from T7 (Lemma 4b on x ≤ 34.249, T7's far bound beyond). Remainder boxes by Agmon on ψ_c(1 + x²)^{-1}v^{(m)}.
    - Result: √E_j^sh[err] ≤ 1.832e3, 1.014e7, 1.231e9 · Ȳ (j = 0, 1, 2). Item 2 limits are 2.39e-7, 9.85e-10, 3.63e-10 (was 5.57e-8, 1.34e-7, 3.08e-7 with the windows' factor), so ×363 at Ȳ = 1e-12. Not binding.
  - **Hand-off frame in the chart (Nov 12; B5 follow-up).** y_h = (y_ray − x_c cosΘ)/λ enters D_h = ((Y − c)/M)d/dx and the weight 28λM/(Y − c)² dx, for the windows and the shells.
    - (x_c, λ) is the union of the 32 slab-180 IFT balls, widened by the true-curve shift 2·max(C_tr, C_sc)·Ȳ. Here C_tr = 9.66 and C_sc = 5.18 come from clam's transfer M^{-T} (stored M, f_norms, tail t).
    - c ∈ [−3.013e-5, −2.878e-5] and λ ∈ [0.99999490, 0.99999511] (`t9_cone_kch.json` frame [0.99999490736, 0.99999510645], rounded outward). KCH 1.187780 → 1.187777; window limits unchanged at 4 digits (E3 1.518e-11 on the A1 data; final 1.575e-11).

## 4. Files and commands (final tail: `p2/final_tail.py` → `p2/out/FINAL_manifest.json`; the A1 commands below are kept for the upstream stages)
```
cd /home/claude/s1spec/step4/p2
python3 -I final_tail.py   # final tail on the B1 data, one version of every input: switch ... far H1, axis crossing, t1clock; writes out/FINAL_manifest.json
python3 -I final_tail.py --from 'h4loc on the final amp'   # cleanup rerun: stages 0-19 reused after checking their SHA-256
python3 -I se_corr.py; python3 -I switch.py ../t3/p1/ckpt/slab_00096.json 2 0.873067570143915
T6RES=$PWD/t6/residuals.json SLAB0=2 YOUT=$PWD/out/ycomp_fin.json python3 -I ycomp.py 6.21e-13
# T8: cconst/merge -> constants_fin; camp -> amp_fin; n3far (it0) -> camp --n3 -> amp_fin_n3_it0; n3far --amp amp_fin_n3_it0 -> n3.json
cd t9h; export COHERENT=1 S1PROF=/home/claude/s1_tight/ext_tip/prof_tip.txt S1TT=10.249077651120778,10.249077651121361
./run_h1slab.sh 0 7 & ./run_h1slab.sh 8 15      # 16 x h1.py p<t_k> (Krawczyk + Lambda), 32 x h1dt.py <tc> 1/32 (TREF = t_k)
python3 -I win_slab.py > ../out/h1slab/win_cmds.txt; ./run_winslab.sh 0 15 & ./run_winslab.sh 16 31; python3 -I clam.py
cd ../t8se; T9JSON=../out/t9_cone_kch.json T9WIN=../out/h1slab/win_summary.json python3 -I camp.py --consts constants/constants_fin.json \
   --F0 4e-23,3e-22,5e-19,3.4e-15 --n3 constants/n3.json --window 180-180 --out constants/amp_A1_w180.json
python3 -I n3boxes.py constants/n3.json constants/amp_A1_w180.json constants/constants_fin.json 1e-12
H4OUT=../out/h4loc_A1_w180.json python3 -I h4loc.py constants/amp_A1_w180.json
cd ..; A=t8se/constants/amp_A1_w180.json; C=t8se/constants/constants_fin.json
python3 -I boxes.py $A $C 1e-12 t8se/constants/n3.json; python3 -I t9.py $A; python3 -I t9_far.py $A $C; python3 -I t1clock.py $C 180 1e-12
# t9 <-> camp: t9's KCH, c0, c0tip use only the amp's per-slab aF0, bF0 (not its limits); run t9, camp, t9 (the second t9 output is identical)
LO_NS=48 python3 -I t9_lemmaO.py $A 180; python3 -I t9_h6.py; python3 -I t9chart_check.py   # the last: offline sympy check of t9chart.py
```

## 5. Theorem 4 (draft)
**Theorem 4.** Let Γ0 be the closed curve X(α) = X0 + (L0/2π)∫₀^α e^{iθ}, with θ = −π/2 + α + Σ_{n≤3762} b_n sin nα, and with the exact data of `t2/gamma0_final.json`:
- b_n exact doubles for n ≠ 2;
- b₂ the unique closing value in its certified interval;
- L0 = 0x1.9e9c2d62f0874p+6.

Let S1 be any shrinker in the certified step-1 set (B8). Let γ be the curve diffusion flow from Γ0. Then there is t₁ ∈ [4.5978574775, 4.5978586144] (time since Γ0; the good time τ* in slab 180), with λ₁ and c₁, such that:
- (a) γ is smooth and axially symmetric on [0, t₁], and the base curve's maximal curvature at t₁ lies in [2.279e4, 2.376e4] **[certified]**;
- (b) Γ₁ = (γ(t₁) − c₁e₁)/λ₁ satisfies the step-3 hand-off hypotheses:
  - (H1) **[certified: |α| ≤ 1.1739e-4, √Λ ≤ 2.29e-8]**;
  - (H2), (H3), (H4) **[certified]**;
  - (H5) **[certified: Lemma O H1–H5 for both families (H1 arc by arc, cleanup O1; H8 is not required). Lemma O H6: closed in v3.2 by the covering argument plus Lemma R, except for the strip input, which belongs to the step-3 gap Hypothesis (F) (§3)]**;
  - (H6) **[step 3; to be stated]**.

Hence, by step 3, γ develops a finite-time singularity modelled on S1. **This is conditional on Hypothesis (F) of step 3: Lemma S(B) drops the window-edge terms (`writeup/step3_gap_v32/gap_v32.pdf`, `REFEREE_check.md` there). The main theorem is stated conditionally on (F) (`writeup/main/main.pdf`). Theorem 4 itself (step 4's results) does not depend on (F).**

## 6. v3.2 hand-off (step-3 amendment v3.2)
Candidate: far (M, K) from (1, 1.5) to (0.15, 0.3); cone unchanged (M = 0.12, K = 0.18); six H6-data thresholds added to (H5). Design: `p2/h6/fable/DESIGN.md` §4, §7. All v3.1 outputs and rows above are unchanged; the v3.2 outputs sit beside them (manifest `p2/out/v32_manifest.json`, SHA-256; **superseded by `p2/out/FINAL_manifest.json`**, whose "v3.2: ..." stages and the far-H1 stage pin these outputs). Every margin below is a certified lower bound (python-flint `arb`, outward rounding). **No check fails.** Margin = threshold / bound.

| item | bound | against | margin | source |
|---|---|---|---|---|
| **cone slope (data), thinnest** | 0.05724 | 0.08 | **×1.40** | `t9_v32_checks.json` |
| cone sup\|h₀″\| (K2 data) | 0.01883 | 0.03 | ×1.59 | `t9_v32_checks.json` |
| cone ‖h₀‴‖_L2 (N3 data) | 1.170e-3 | 0.004 | ×3.41 | `t9_v32_checks.json` |
| far slope (data) | 0.04153 | 0.08 | ×1.92 | `t9_far_v32.json` |
| far h₂ (physical, λ₁ units) | 1.797e-4 | 3.53e-4 | ×1.96 | `t9_far_v32.json` |
| far h₃ (physical) | 2.595e-8 | 1.17e-5 | ×450.9 | `t9_far_v32.json` |
| far h₄ (physical) | 6.964e-9 | 1.03e-7 | ×14.79 | `t9_far_v32.json` |
| far ‖h₃‖_L2 (physical) | 5.516e-7 | 3.32e-6 | ×6.01 | `t9_far_v32.json` |
| far h₂, N3 in cylinder units (vs 0.04, 0.004) | 0.02031, 6.626e-4 | | ×1.96, ×6.03 | `t9_h6_v32.json` |
| far jets vs M = 0.15: slope, h₃, h₄ | 0.0416, 3.313e-4, 0.01005 | 0.15 | ×3.61, ×452.8, ×14.92 | `t9_h6_v32.json` |
| closure: cone slope (strip is worst; covered 0.08753) | 0.13532 | 0.25 | ×1.85 | `t9_h6_v32.json` |
| closure: far slope | 0.08753 | 0.25 | ×2.86 | `t9_h6_v32.json` |
| closure: cone K2 (r_c = 380, r_p = 150; strip 0.04969) | 0.10664 | 0.18 | ×1.69 (thinnest closure item) | `t9_h6_v32.json` |
| closure: cone K3 | 0.09348 | 0.18 | ×1.93 | `t9_h6_v32.json` |
| closure: far K2 | 0.13125 | 0.3 | ×2.29 | `t9_h6_v32.json` |
| closure: far K3 | 0.12549 | 0.3 | ×2.39 | `t9_h6_v32.json` |
| Lemma O t\* cone (unchanged, grid 14) | 3.683e-7 | 1.6001e-7 | ×2.30 | `t9_h6_v32.json` |
| Lemma O t\* far at (0.15, 0.3) | 9.511e-8 | 6.134e-9 | ×15.5 | `t9_h6_v32.json` |
| Lemma O H4 far at M = 0.15 (centre 376; chart computed, cleanup O6b) | 0.1761 M² | 1 | ×5.68 | `t8se/constants/h4loc_v32.json` |
| Lemma O H4 cone (r_c 370-380; chart computed) | 0.2393 M² | 1 | ×4.18 | `h4loc_v32.json` |
| Lemma O H1 far: sup\|h₀\| (cylinder units; H = 0.0417, C₀√t₁ ≤ 9.212e-4) | 0.04072 | H − C₀√t₁ = 0.04077 | slack 6.7e-5 (graph: turning ≤ 0.04095 rad, \|ξ(±G)\| ≥ 226.18 > 226) | `t9_lemmaO_far.json` |
| time-equivalent margin of the closure (cone K2) | | | ×3.86 in time | `t9_h6_v32.json` |

- **Thinnest margin: cone slope data, ×1.40** (0.05724 against the new threshold 0.08). Among the closure items the thinnest is cone K2, ×1.69. The cone slope and h″ thresholds are design choices that trade against the closure margins (DESIGN §6).
- **Lemma O H1 (cone), `t9_lemmaO_v32.json`** with C₀√t = 0.0023: sup|h₀| ≤ 0.1872, H ≤ 0.1895 (v3.1, C₀ = 0.0021: H ≤ 0.1893). The 18 cone cylinders pass. The H8 distance is ≥ 0.136 L against 0.1 (×1.36, unchanged from v3.1). H8 is not a requirement (§3). The far rows of `t9_lemmaO.py` are identical to v3.1 and informational; the far family's H1 is certified by `t9_lemmaO_far.py` (§3, closed; row above).
- Far thresholds in λ₁ units are the DESIGN values (3.53e-4, 1.17e-5, 1.03e-7, 3.32e-6, slope 0.08), each asserted ≤ the exact value from M2F/113, M/113², M/113³, N3F/113^1.5.
- Final data (Nov 12): the v3.2 data rows come from the final tail on the B1 data (`final_tail.py`, stages "v3.2: ..."); the closure rows and the Lemma O t\* rows depend only on design values and are unchanged. `t9_h6.py`/`t9_h6.json` (v3.1) are left in place and still show the v3.1 failure.
- Reruns: `python3 -I final_tail.py` (from `p2`; T6DIR = p2/t6c, T8 = constants_b1, n3_b1); `--from STAGE` reruns from a stage and reuses the earlier ones after checking their hashes.
