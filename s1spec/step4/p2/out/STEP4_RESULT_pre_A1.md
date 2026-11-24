# STEP 4 RESULT (base v2): certified chain Γ0 → switch → P2 → hand-off sub-window W′

Status (Nov 11, 15:37 PT). Everything below is computed in arb with outward rounding and no sampling, unless it is marked otherwise.

The chain is certified from Γ0 to the hand-off sub-window W′, which is slab 180 with t_c ∈ [1 − 1/512, 1] (Δτ = 1.46e-4, k ≈ 2.37e4). This includes α by the step-3 IFT and the base window and shell energies.

**Status withdrawn pending the independent review (`STEP4_REVIEW.md`): Theorem 4 is not yet proved; A1 (Lemma O H4 at t₁), A2 (nonlinear terms in (Y)) and A3 (tip pointwise box) are open.** H1 gives √Λ ≤ 2.30e-8 over all of W′, against 5e-8 (§3, G1), and T8 has closed G4. The thinnest certified margins are:
- the far box at level 3, ×1.85 in Ȳ (the box needs Ȳ ≤ 1.85e-12);
- the cone K4 (H5), ×1.61.

The v1 report is kept as `p2/out/STEP4_RESULT_v1_stopped.md`.

## 1. Chain
| link | result | margin |
|---|---|---|
| Γ0 (T2) | 3762 sine modes; b₂ by Krawczyk (half-width 8.1e-24); embedded | — |
| P1 (T3), slabs 0–96 | t₉₆ = 4.597850539485830831 ± 7.8e-19 (time since Γ0); k ∈ [29.352874545326, 29.352874545327] | — |
| base v2 (T1/T6) + (SE) §2(ii) | slabs 0–180. The bump on u ∈ [0.6, 0.8] makes (SE) exact (\|a₂\| ≤ 1.6e-14); certified increments are added to all far T6 norms. The jet hulls are widened outward with `nextafter` (rounding fix, rerun) | exact |
| switch (A) | P1 slab 96 ↔ base slab 2, t_c,a = 982986695892350/2⁵⁰. IFT q = 4.2e-9. **Y(τ_a) ≤ 6.20e-13**; F_j(τ_a) ≤ 3.7e-23, 2.4e-22, 4.7e-19, 3.2e-15 | ×1.6 |
| T7 (republished, gaps 9 and 10 closed) | 26 certificates, ϱ_c ≥ 0.0471. The SE-corrected far constants and farg are now read from `t7/constants.json` (far n moved ≤ 8.7e-7 relative). Lemma 4a removes the e″ tip box: its cubic term m₂₃\|e″\|³ (m₂₃ ≤ 3.0003) costs ≤ 4.0e-9 in rate, inside DN = 1e-6, and ≤ 3.0e-9 of εD₄². That slice is reserved: the residual share is now ε(1/16 − 3e-9) | — |
| (Y) | level ≤ 1.33e-13 → **Ȳ = 1e-12** | — |
| (F), T8 on W′ | Lemma N3′ and 𝒬 are certified per slab at levels 0–3 (T8, `n3far.py`) and applied by `camp.py --n3`. The top loss 0.0196 is kept in every θ_j, so the top term is counted twice (conservative). KCH = 1.152 (arb). With the C5 good time inside W′, D₂(τ*) ≤ 1.23e-8. Item limits Ȳ ≤: item 4 **5.048e-12**; item 3 E3 1.79e-11; item 2 E0 5.9e-8 | ×5.0 |
| boxes at Ȳ (with N3′) | tip \|e\|, \|e′\| ≤ 1.25·2e3Ȳ/√μ_e = 1.32e-8 (T8's certified box) against 1e-5. Far core: \|e′\| ≤ 5.43e-4, x\|e″\| ≤ 7.54e-5, **x²\|e‴\| ≤ 5.69e-3** against 1e-2; it is 9.84e-3 at Ȳ = 1.85e-12, so the box closes for every Ȳ ≤ 1.85e-12. Lemma 4a gives S_W ≤ 1.14e-6, which proves the e″ part of the box on [16.25, 34.25]. \|κ̃e\| ≤ 2.39e-3 (consistent with the top loss) | ×18, ×132, **×1.85 in Ȳ** |

## 2. Hand-off items at τ* ∈ W′ (T9, `p2/t9h/`, `p2/t9.py`)
- **S1 chart.** It uses only step-1 objects: prof_tip (σ ≤ 8) and ext813 (σ ∈ [8, 13]). prof_tip is the S1 agent's drop-in for prof_tight, conditioned on M = 0 (`/home/claude/s1_tight/ext_tip`). Its boxes are nested in prof_tight's, with θ width ≤ 3.2e-12 on x ∈ [0, 2] (it was 1.76e-9), and the conditioned τ_tip has width 5.8e-13. Correlation in τ_tip is kept. Beyond σ = 13 it uses tail envelopes; Θ = THINF comes from the identity ∫_σ^∞k = p/σ + 4m/σ² + … .
  - **Θ = 1.8617364190148 ± 5.2e-11**; it was ± 1.7e-9 with prof_tight and ± 2.3e-7 with envelopes only.
- **The step-3 IFT for α** uses the constraints ⟨χw, f_j⟩_ρ = 0 on the 3425 step-3 IBP hull boxes, with f₁, f₂ and f₃ (adjoint w*).
  - M = [⟨χψ̂₃, f₃⟩ = 0.030 ± 0.004; f₃ is decoupled from the frame modes to 1.5e-4]. Two linear passes are made on each of the 32 t_c-balls that cover W′, and at the midpoint.
  - **α = −1.1296366e-4 ± 3.0e-9 over W′** (−1.1296360e-4 ± 1.8e-10 at the midpoint). So **|α| ≤ 1.1297e-4** against 1.5e-4 (×1.33). The error adds at most C_αȲ, with C_α ≤ 0.72.
  - With envelopes only (before ext813): ±0.47%.

| item | base (arb, W′) | with error (Ȳ = 1e-12) | margin |
|---|---|---|---|
| (H1) α | −1.1296366e-4 ± 3.0e-9 over W′ | + C_αȲ ≤ 7.2e-13 | ×1.33 |
| (H1) √Λ | **≤ 2.291e-8 over W′** (32 t_c-balls). The worst ball has A ≤ 2.25e-8, √ε B ≤ 4.4e-9, √μ V₀ ≤ 2.9e-10, and rotation mode 1.7e-12 × 17.64. The midpoint alone gives 6.72e-9 | + C_ΛȲ ≤ 2.4e-11 (C_Λ ≤ 24.1); holds for Ȳ ≤ 1.58e-9 | **×2.18** |
| (H2) shells [14, 60] | 2.8e-12, 3.7e-12, 8.9e-11 | Ȳ ≤ 5.9e-8 | ×5.9e4 |
| (H3) windows E0–E3 | 1.44e-8, 5.19e-8, 1.70e-7, 4.70e-7 | Ȳ ≤ 1.79e-11 (E3) | ×17.9 |
| (H3) K4, y ∈ [50, 350] | 3.2e-4 (hulls) | 3.31e-2 against 0.06 | ×1.81 |
| (H4) graph over the ray | \|θ − Θ\| ≤ 0.079 | + 5.4e-4 | — |
| (H5) cone, r_c ∈ [200, 380] | — | K2, K3, K4 ≤ 1.89e-2, 4.86e-4, **7.45e-2** against 0.12; N3 ≤ 1.19e-3 against 0.14 | ×6.3, ×246, **×1.61**, ×117 |
| (H5) far cylinders, x ≥ 150 (`p2/t9_far.py`) | \|κ̃\|, \|κ̃′\|, \|κ̃″\| ≤ 1.76e-4, 1.3e-8, 5.3e-11 (slab 180 hulls) | \|h″\|, \|h‴\|, \|h⁗\| ≤ 1.79e-4, 2.62e-8, 4.22e-8; slope ≤ 0.0416; ‖h‴‖_L² ≤ 5.56e-7 | ×49, ×2679, ×14.7, ×5.0, ×1499 (scripted, arb) |

## 3. Gaps
Open: A1–A3 of `STEP4_REVIEW.md` (and its B items); see the review. Closed in this round:
- **G1 (√Λ, H1) is closed.**
  - The tip offset C = −π/2 − θ(tip) + dT·k_tip goes through the projection as one scalar (`COHERENT=1`). Λ is evaluated at C's midpoint, and ‖Πq‖_Λ ≤ 17.64 bounds the rest (q is the rotation mode). With prof_tip and the conditioned τ_tip, rad(C) = 1.7e-12.
  - W′ is covered by 32 t_c-balls of radius 2⁻¹⁵, each with its own two-pass IFT and Λ cells. This bounds √Λ at every τ ∈ W′, including the C5 good time, so no separate timing term is needed. The ≤ 1e-11 timing estimate itself is not certified; that would need a ∂_τw pass.
  - The drift of (xc, λ, α) across each ball now limits the bound (A ≈ 2.25e-8; the midpoint alone gives 6.72e-9).
  - The validated error enters through the same projection (T5_DESIGN l.161, `t9h/clam.py`): √Λ[w] ≤ √Λ[w̃] + C_ΛȲ, with C_Λ ≤ 24.1, and |α − α̃| ≤ C_αȲ, with C_α ≤ 0.72. The bound uses:
    - Λ^{S3}[e] ≤ (0.08/μ_e)Λ_e, with μ_e = 0.036;
    - Cauchy–Schwarz on ⟨χe, f_j⟩_ρ, with ‖f_j‖_ρ ≤ 1.84, 1.19 and 0.0196;
    - T8's tip box on the cutoff tail.
- **G4 (Lemma N3′ and 𝒬) is closed by T8.** `t8se/n3far.py` certifies N3′ at levels 0–3 and 𝒬 per slab, and `camp.py --n3` applies them. Both were rerun on this pipeline (`n3_fin.json`, `amp_n3_fin.json`) and reproduce T8's numbers.
- **Closed earlier:**
  - G2: arb base windows, shells and K4 at W′;
  - G5: item 5 scripted (`p2/t9_far.py`);
  - KCH (arb);
  - gap 9: Lemma 4a;
  - gap 10: SE far constants;
  - the rounding bug: outward `nextafter` hulls throughout p2.

## 4. Files and commands (nice 10)
```
cd /home/claude/s1spec/step4/p2
python3 -I se_corr.py; python3 -I switch.py ../t3/p1/ckpt/slab_00096.json 2 0.873067570143915
T6RES=$PWD/t6/residuals.json SLAB0=2 YOUT=$PWD/out/ycomp_fin.json python3 -I ycomp.py
cd t8se; TOPLOSS=0.0196 python3 -I cconst.py --slabs 2-180 --out constants/slabs_fin; python3 -I merge.py --in constants/slabs_fin --first 2 --last 180 --out constants/constants_fin.json
SUBW=1.46484375e-4 T9JSON=../out/t9_cone_kch.json python3 -I camp.py --consts constants/constants_fin.json --F0 4e-23,3e-22,5e-19,3.4e-15 --out constants/amp_fin.json
python3 -I n3far.py --slabs 2-180 --out constants/n3_fin.json --amp constants/amp_fin.json --consts constants/constants_fin.json --Y 1.85e-12
SUBW=1.46484375e-4 T9JSON=../out/t9_cone_kch.json python3 -I camp.py --consts constants/constants_fin.json --F0 4e-23,3e-22,5e-19,3.4e-15 --n3 constants/n3_fin.json --out constants/amp_n3_fin.json
cd ..; P="t8se/constants/amp_n3_fin.json t8se/constants/constants_fin.json"
python3 -I boxes.py $P 1e-12; python3 -I boxes.py $P 1.85e-12; python3 -I t9.py t8se/constants/amp_n3_fin.json; python3 -I t9_far.py $P
cd t9h; export COHERENT=1 S1PROF=/home/claude/s1_tight/ext_tip/prof_tip.txt S1TT=10.249077651120778,10.249077651121361
H1TAG=_tip python3 -I h1.py m                                            # midpoint -> out/t9_h1_m_coh_tip.json
for K in $(seq 0 31); do H1NSUB=32 H1TAG=_tip32 python3 -I h1.py $K; done   # W' cover -> out/t9_h1_K_coh_tip32.json
python3 -I clam.py                                                       # C_Lambda, C_alpha -> out/t9_h1_clam.json
T9IFT=../out/t9_ift_tight.json WINOUT=../out/t9_win_tight.json python3 -I win.py      # windows/shells at W'
```
- KCH and c0 depend only on the base, so camp can read `t9_cone_kch.json` before t9.py reruns. t9.py now writes the c0 allowance itself (+0.01 per level and 0.002 at level 0, as used since v2); `c0_chain` is the bare chain-rule term.

## 5. Theorem 4 (draft)
**Theorem 4.** Let Γ0 be the closed curve X(α) = X0 + (L0/2π)∫₀^α e^{iθ}, with θ = −π/2 + α + Σ_{n≤3762} b_n sin nα, and with the exact data of `t2/gamma0_final.json`:
- b_n exact doubles for n ≠ 2;
- b₂ the unique closing value in its certified interval;
- L0 = 0x1.9e9c2d62f0874p+6.

Let γ be the curve diffusion flow from Γ0. Then there are t₁ ∈ [4.5978580262, 4.5978580263] (time since Γ0), λ₁ and c₁ such that:
- (a) γ is smooth and axially symmetric on [0, t₁], and k_max(γ(t₁)) ∈ [2.36e4, 2.38e4] **[certified]**;
- (b) Γ₁ = (γ(t₁) − c₁e₁)/λ₁ satisfies the step-3 hand-off hypotheses:
  - (H1) |α| ≤ 1.5e-4 and √Λ ≤ 5e-8 **[certified: 1.1297e-4 and 2.30e-8 over W′]**;
  - (H2) **[certified]**;
  - (H3) **[certified]**;
  - (H4) **[certified]**;
  - (H5) **[certified, scripted]**;
  - the nonlinear terms at levels 0–3 use Lemma N3′ and 𝒬 **[certified, T8 `n3far.py`]**.

Hence, by step 3, γ develops a finite-time singularity modelled on S1. **Not yet proved: conditional on A1–A3 of `STEP4_REVIEW.md`.**
