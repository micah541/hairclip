# STEP 4 RESULT (base v2): certified chain Γ0 → switch → P2 → hand-off on slab 180

Status (Nov 12). **Theorem 4 is not yet proved.** A1–A3 and B8 are closed. Open (§3): B1, part of B5, Lemma O H6 (step-3 amendment), and Lemma O H6, which needs a change in step 3. Everything below is arb with outward rounding and no sampling, unless it is marked otherwise. Every quoted bound is rounded outward.

The hand-off good-time window is now **all of slab 180** (t_c ∈ [−1, 1], Δτ = 0.15, k ∈ [2.27e4, 2.38e4]). The old W′ (1.46e-4 long) failed Lemma O H4 by ×5–×280 (A1). All hand-off items are certified uniformly for τ* anywhere in the slab. The thinnest certified margins are:
- the level-3 far box x²|e‴|, ×1.68 (Ȳ ≤ 1.7e-12);
- |α|, ×1.27.

## 1. Chain
| link | result | margin |
|---|---|---|
| Γ0 (T2) | 3762 sine modes; b₂ by Krawczyk (half-width 8.1e-24); embedded | — |
| P1 (T3), slabs 0–96 | t₉₆ = 4.597850539485830831 ± 7.8e-19; k ∈ [29.352874545326, 29.352874545327] | — |
| base v2 + (SE) | slabs 0–180, (SE) exact up to \|a₂\| ≤ 1.6e-14; T6 hulls and ranges widened outward on read; residual summary regenerated | — |
| switch (A), base slab 2 | t_c,a = 982986695892350/2⁵⁰, so τ_a lies in slab 2 and P2 starts there (slabs 0–1 are unused). **Y(τ_a) ≤ 6.21e-13**; F_j(τ_a) ≤ 3.75e-23, 2.45e-22, 4.72e-19, 3.19e-15 | ×1.6 |
| T7 | Nonlinear closure (Lemma 4c, A2) certified on all 25 intervals, min ϱ_c ≥ 0.04711; DN = D4CUB = 0. Tip box (Lemma 4b, A3): sup\|e\| ≤ 1.6031e4·Λ_e^½ and sup\|e′\| ≤ 1.6352e4·Λ_e^½ on [0, 34.25] | — |
| (Y) | level ≤ 1.33e-13; ycomp asserts Ȳ = 1e-12 ≥ Y(τ_a) | — |
| (F), T8, window = slab 180 | N3′ and 𝒬 per slab (`n3far`, b₀, b₁, b₂ = 1.5e-3, 3e-4, 1e-2, fac 1.7, Y = 1.7e-12); `camp --n3 --window 180-180`, with the slab-wide base windows. D₂(τ*) ≤ 1.27e-11. Ȳ limits: **item 3 E3 1.56e-11**; item 4 2.86e-11; item 2 E0 5.7e-8 | ×15.6 |
| boxes at Ȳ = 1e-12 (`n3boxes`, `boxes.py`) | \|e′\|, x\|e″\|, x²\|e‴\| ≤ 5.43e-4, 7.57e-5, 5.92e-3, against the N3′ boxes; \|e\|: √F₀ ratio ≤ 1.00001 against fac 1.7. All four hold at Y = 1.7e-12. Tip \|e\|, \|e′\| ≤ 1.61e-8, 1.64e-8 against 1e-5 | ×2.76, ×3.96, **×1.68**, ×1.69 |

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
| (H1) √Λ | ≤ 2.282e-8 | ≤ 2.29e-8; holds for Ȳ ≤ 1.58e-9 | ×2.18 |
| (H2) shells [14, 60] | 8.01e-11, 1.34e-10, 2.94e-10 | Ȳ ≤ 5.7e-8 | ×5.7e4 |
| (H3) windows E0–E3 | 3.48e-8, 1.27e-7, 3.76e-7, 9.71e-7 (against 2e-7, 1e-5, 1e-5, 1e-5) | Ȳ ≤ 1.56e-11 (E3) | ×15.6 |
| (H3) K4, y ∈ [50, 350] | 3.1633e-4 | 6.53e-3 against 0.06 | ×9.1 |
| (H4) graph over the ray | \|θ − Θ\| ≤ 0.078 | + 5.5e-4 | — |
| (H5) cone, r_c ∈ [200, 380] | exact curvature jets (`curvjets.py`) | K2, K3, K4 ≤ 1.89e-2, 4.90e-4, 1.44e-2 against 0.12; N3 ≤ 1.21e-3 against 0.14; slope ≤ 0.0573 | ×6.3, ×244, ×8.3, ×115 |
| (H5) far cylinders, x ≥ 150 | slab-180 hulls, exact jets | \|h″\|, \|h‴\|, \|h⁗\| ≤ 1.84e-4, 2.81e-8, 8.39e-9; slope ≤ 0.0428; ‖h‴‖_L² ≤ 5.98e-7 | ×47, ×2490, ×73, ×4.9, ×1390 |
| (H5) Lemma O H4, ∫(h₀⁽⁵⁾)² | base ≤ 6.8e-6 | worst cone r_c 370–380 ≤ 0.305 M²; far ≤ 0.0052 M² (`h4loc`) | ×3.28 |

## 3. Gaps
- **Open:**
  - **B1** (T8): the `se_corr.py` majorants are not derived.
  - **B5, remainder.**
    - The c0 reparametrisation allowance (+0.01 per level, 0.002 at level 0) and the KCH factor 1.0001 are still typed. The linear reparametrisation term is sec²φ tanφ κ̃·e. A computed bound needs its D^j and the nonlinear graph-variation remainder, which is not yet done.
    - The `u4_cell` hull is now exact (closed).
  - **Lemma O H1, arc by arc** (H8 is not required):
    - Cone family: sup|h₀| ≤ 0.19 (`t9_lemmaO.py`).
    - Far family: the tangent turns ≤ 0.092 rad over an arc, and sup|h₀| ≤ (X²/2)·L·h₂ = 0.042 (`t9_lemmaO_far.json`, t9_far bounds).
    - *Informational:* the hand-off curve self-intersects (`p2/out/t9_axis_crossing.json`, axis crossing near x ≈ 3.27e4). This is the expected loss of embeddedness (STEP4_DESIGN §2). Step 3 applies Lemma O arc by arc and needs only H1–H7 and H9 (step3_proof l.97, l.105, l.111), so H8 is not a requirement.
  - **Lemma O H6 fails as specified.** The amendment plan is in `p2/h6/H6_PLAN.md`: cone K = 0.25, far (M, K) = (0.1, 0.2), and the |x| ≤ 1 covering. Its thinnest margins are ×1.13 (cone K₃) and ×1.14 (cone K₂). The details below are from `t9_h6.py`; step 3 is not changed here.
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
  - A1, A2, A3; B2, B3, B4, B6, B7, B9.
  - B8: `/home/claude/s1_tight/consistency/check_outA5.py` shows all 1025 + 1025 + 501 rows of prof_tip, prof_tight and ext813 lie inside the outA5/ibp_boxes hulls, with positive margins. Theorem 4 is stated for an S1 in the certified set.
  - The writer's rounding and stale-summary items.

## 4. Files and commands (nice 10; manifest `p2/out/A1_manifest.json`, SHA-256, in run order)
```
cd /home/claude/s1spec/step4/p2
python3 -I se_corr.py; python3 -I switch.py ../t3/p1/ckpt/slab_00096.json 2 0.873067570143915
T6RES=$PWD/t6/residuals.json SLAB0=2 YOUT=$PWD/out/ycomp_fin.json python3 -I ycomp.py 6.21e-13
# T8: cconst/merge -> constants_fin; camp -> amp_fin; n3far (it0) -> camp --n3 -> amp_fin_n3_it0; n3far --amp amp_fin_n3_it0 -> n3.json
cd t9h; export COHERENT=1 S1PROF=/home/claude/s1_tight/ext_tip/prof_tip.txt S1TT=10.249077651120778,10.249077651121361
./run_h1slab.sh 0 7 & ./run_h1slab.sh 8 15      # 16 x h1.py p<t_k> (Krawczyk + Lambda), 32 x h1dt.py <tc> 1/32 (TREF = t_k)
python3 -I win_slab.py > ../out/h1slab/win_cmds.txt; ./run_winslab.sh 0 15 & ./run_winslab.sh 16 31; python3 -I clam.py
cd ../t8se; T9JSON=../out/t9_cone_kch.json T9WIN=../out/h1slab/win_summary.json python3 -I camp.py --consts constants/constants_fin.json \
   --F0 4e-23,3e-22,5e-19,3.4e-15 --n3 constants/n3.json --window 180-180 --out constants/amp_A1_w180.json
python3 -I n3boxes.py constants/n3.json constants/amp_A1_w180.json constants/constants_fin.json 1e-12; python3 -I h4loc.py constants/amp_A1_w180.json
cd ..; A=t8se/constants/amp_A1_w180.json; C=t8se/constants/constants_fin.json
python3 -I boxes.py $A $C 1e-12 t8se/constants/n3.json; python3 -I t9.py $A; python3 -I t9_far.py $A $C; python3 -I t1clock.py $C 180 1e-12
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
  - (H5) **[certified: Lemma O H2–H5 for both families, and H1, H8 for the cone family. Open: H6. H6 fails as specified (§3); it needs a change in step 3]**;
  - (H6) **[step 3; to be stated]**.

Hence, by step 3, γ develops a finite-time singularity modelled on S1. **This is conditional on B1, the B5 remainder (c0 and KCH), Lemma O H6 (step-3 amendment), and Lemma O H6, which needs a change in step 3 (§3).**

## 6. v3.2 hand-off (candidate; pending the user's approval of the step-3 amendment)
Candidate: far (M, K) from (1, 1.5) to (0.15, 0.3); cone unchanged (M = 0.12, K = 0.18); six H6-data thresholds added to (H5). Design: `p2/h6/fable/DESIGN.md` §4, §7. All v3.1 outputs and rows above are unchanged; the v3.2 outputs sit beside them (manifest `p2/out/v32_manifest.json`, SHA-256). Every margin below is a certified lower bound (python-flint `arb`, outward rounding). **No check fails.** Margin = threshold / bound.

| item | bound | against | margin | source |
|---|---|---|---|---|
| **cone slope (data), thinnest** | 0.05724 | 0.08 | **×1.40** | `t9_v32_checks.json` |
| cone sup\|h₀″\| (K2 data) | 0.01883 | 0.03 | ×1.59 | `t9_v32_checks.json` |
| cone ‖h₀‴‖_L2 (N3 data) | 1.200e-3 | 0.004 | ×3.33 | `t9_v32_checks.json` |
| far slope (data) | 0.04279 | 0.08 | ×1.87 | `t9_far_v32.json` |
| far h₂ (physical, λ₁ units) | 1.839e-4 | 3.53e-4 | ×1.92 | `t9_far_v32.json` |
| far h₃ (physical) | 2.808e-8 | 1.17e-5 | ×416.6 | `t9_far_v32.json` |
| far h₄ (physical) | 8.380e-9 | 1.03e-7 | ×12.29 | `t9_far_v32.json` |
| far ‖h₃‖_L2 (physical) | 5.970e-7 | 3.32e-6 | ×5.56 | `t9_far_v32.json` |
| far h₂, N3 in cylinder units (vs 0.04, 0.004) | 0.02078, 7.172e-4 | | ×1.93, ×5.58 | `t9_h6_v32.json` |
| far jets vs M = 0.15: slope, h₃, h₄ | 0.0428, 3.586e-4, 0.01209 | 0.15 | ×3.51, ×418.3, ×12.40 | `t9_h6_v32.json` |
| closure: cone slope (strip is worst; covered 0.08753) | 0.13532 | 0.25 | ×1.85 | `t9_h6_v32.json` |
| closure: far slope | 0.08753 | 0.25 | ×2.86 | `t9_h6_v32.json` |
| closure: cone K2 (r_c = 380, r_p = 150; strip 0.04969) | 0.10664 | 0.18 | ×1.69 (thinnest closure item) | `t9_h6_v32.json` |
| closure: cone K3 | 0.09348 | 0.18 | ×1.93 | `t9_h6_v32.json` |
| closure: far K2 | 0.13125 | 0.3 | ×2.29 | `t9_h6_v32.json` |
| closure: far K3 | 0.12549 | 0.3 | ×2.39 | `t9_h6_v32.json` |
| Lemma O t\* cone (unchanged, grid 14) | 3.683e-7 | 1.6001e-7 | ×2.30 | `t9_h6_v32.json` |
| Lemma O t\* far at (0.15, 0.3) | 9.511e-8 | 6.134e-9 | ×15.5 | `t9_h6_v32.json` |
| Lemma O H4 far at M = 0.15 (centre 376) | 0.2286 M² | 1 | ×4.37 | `t8se/constants/h4loc_v32.json` |
| Lemma O H4 cone (unchanged, r_c 370-380) | 0.3044 M² | 1 | ×3.29 | `h4loc_v32.json` |
| time-equivalent margin of the closure (cone K2) | | | ×3.86 in time | `t9_h6_v32.json` |

- **Thinnest margin: cone slope data, ×1.40** (0.05724 against the new threshold 0.08). Among the closure items the thinnest is cone K2, ×1.69. The cone slope and h″ thresholds are design choices that trade against the closure margins (DESIGN §6).
- **Lemma O H1 (cone), `t9_lemmaO_v32.json`** with C₀√t = 0.0023: sup|h₀| ≤ 0.1872, H ≤ 0.1895 (v3.1, C₀ = 0.0021: H ≤ 0.1893). The 18 cone cylinders pass. The H8 distance is ≥ 0.136 L against 0.1 (×1.36, unchanged from v3.1). H8 is not a requirement (§3). The far rows are identical to v3.1 and informational (`t9_lemmaO_far.py` unchanged).
- Far thresholds in λ₁ units are the DESIGN values (3.53e-4, 1.17e-5, 1.03e-7, 3.32e-6, slope 0.08), each asserted ≤ the exact value from M2F/113, M/113², M/113³, N3F/113^1.5.
- Not changed: the chain, T7, T8, H1-H4, cone jets, `t9.py`, `camp.py`, `t1clock.py`, all step-3 files. `t9_h6.py`/`t9_h6.json` (v3.1) are left in place and still show the v3.1 failure.
- Reruns (nice 10, from `p2`): `python3 -I t9_far.py --v32 t8se/constants/amp_A1_w180.json t8se/constants/constants_fin.json`; `cd t8se; FAR_M=0.15 H4OUT=constants/h4loc_v32.json python3 -I h4loc.py constants/amp_A1_w180.json`; `python3 -I t9_v32_checks.py`; `python3 -I t9_h6_v32.py`; `LO_NS=48 python3 -I t9_lemmaO.py --v32 t8se/constants/amp_A1_w180.json 180`.
