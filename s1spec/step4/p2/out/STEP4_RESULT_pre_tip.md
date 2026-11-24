# STEP 4 RESULT (base v2): certified chain Γ0 → switch → P2 → hand-off sub-window W′

Status (Nov 11, 14:51 PT). Everything below is computed in arb with outward rounding and no sampling, unless it is marked otherwise.

The chain is certified from Γ0 to the hand-off sub-window W′, which is slab 180 with t_c ∈ [1 − 1/512, 1] (Δτ = 1.46e-4, k ≈ 2.37e4). This includes α by the step-3 IFT and the base window and shell energies.

**One hand-off item is still open: √Λ ≤ 5e-8 (H1).** The certified bound is 8.18e-8 after fix (a) (§3, G1). The thinnest certified margins are:
- the far box at level 3, ×1.8;
- the cone K4 (H5), ×1.67.

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
| (F), T8 on W′ | top loss 0.0196 in every θ_j; N3 sensitivity ×100 changes the result by −2%; KCH = 1.152 (arb). With the C5 good time inside W′, D₂(τ*) ≤ 1.19e-8. Item limits Ȳ ≤: item 4 **5.11e-12**; item 3 E3 1.80e-11; item 2 E0 5.9e-8 | ×5.1 |
| boxes at Ȳ | tip 1.1e-8 against 1e-5. Far core: \|e′\| ≤ 5.4e-4, x\|e″\| ≤ 7.5e-5, **x²\|e‴\| ≤ 5.5e-3** against 1e-2 (needs Ȳ ≤ 1.9e-12). Lemma 4a gives S_W ≤ 1.14e-6, which proves the e″ part of the box on [16.25, 34.25]. \|κ̃e\| ≤ 2.39e-3 (consistent with the top loss) | ×18, ×133, **×1.8** |

## 2. Hand-off items at τ* ∈ W′ (T9, `p2/t9h/`, `p2/t9.py`)
- **S1 chart.** It uses only step-1 objects: prof_tight (σ ≤ 8) and ext813 (σ ∈ [8, 13]). Correlation in τ_tip is kept. Beyond σ = 13 it uses tail envelopes; Θ = THINF comes from the identity ∫_σ^∞k = p/σ + 4m/σ² + … .
  - **Θ = 1.8617364190 ± 1.7e-9**; it was ± 2.3e-7 with envelopes only.
- **The step-3 IFT for α** uses the constraints ⟨χw, f_j⟩_ρ = 0 on the 3425 step-3 IBP hull boxes, with f₁, f₂ and f₃ (adjoint w*).
  - M = [⟨χψ̂₃, f₃⟩ = 0.030 ± 0.004; f₃ is decoupled from the frame modes to 1.5e-4]. Two linear passes are made at the W′ midpoint.
  - **α = −1.129637e-4 ± 3.0e-9**, ± 5.2e-9 over W′. So **|α| ≤ 1.1297e-4** against 1.5e-4 (×1.33).
  - With envelopes only (before ext813): ±0.47%.

| item | base (arb, W′) | with error (Ȳ = 1e-12) | margin |
|---|---|---|---|
| (H1) α | −1.12964e-4 ± 5.2e-9 | + C_αȲ | ×1.33 |
| (H1) √Λ | **≤ 8.18e-8** = 5.31e-8 at the tip-offset midpoint (A ≤ 5.29e-8, √ε B ≤ 5.2e-9, √μ V₀ ≤ 2.8e-10) + rad(C) 1.61e-9 × ‖Πq‖_Λ ≤ 17.64 (rotation mode) | — | **open: needs 5e-8** |
| (H2) shells [14, 60] | 2.8e-12, 3.7e-12, 8.9e-11 | Ȳ ≤ 5.9e-8 | ×5.9e4 |
| (H3) windows E0–E3 | 1.44e-8, 5.19e-8, 1.70e-7, 4.70e-7 | Ȳ ≤ 1.80e-11 (E3) | ×18 |
| (H3) K4, y ∈ [50, 350] | 3.2e-4 (hulls) | 3.20e-2 against 0.06 | ×1.88 |
| (H4) graph over the ray | \|θ − Θ\| ≤ 0.079 | + 5.4e-4 | — |
| (H5) cone, r_c ∈ [200, 380] | — | K2, K3, K4 ≤ 1.9e-2, 4.8e-4, **7.2e-2** against 0.12; N3 ≤ 1.2e-3 against 0.14 | ×6.4, ×250, **×1.67**, ×118 |
| (H5) far cylinders, x ≥ 150 (`p2/t9_far.py`) | \|κ̃\|, \|κ̃′\|, \|κ̃″\| ≤ 1.76e-4, 1.3e-8, 5.3e-11 (slab 180 hulls) | \|h″\|, \|h‴\|, \|h⁗\| ≤ 1.79e-4, 2.60e-8, 4.08e-8; slope ≤ 0.0416; ‖h‴‖_L² ≤ 5.5e-7 | ×49, ×2690, ×15, ×5.0, ×1500 (scripted, arb) |

## 3. Gaps
- **G1. √Λ ≤ 5e-8 (H1) is open; the bound is 8.18e-8 (it was 1.09e-7).** Fix (a) is implemented (`COHERENT=1`, `p2/t9h/h1.py m`, `out/t9_h1_m_coh.json`).
  - The tip-angle offset C = −π/2 − θ(tip) + dT·k_tip is one scalar, radius 1.61e-9. The IFT and Λ are evaluated at its midpoint, which gives α = −1.12964e-4 ± 8.9e-10, xc ± 2.32e-8 and λ ± 8.9e-9.
  - The C-dependence is bounded through the projection: d√Λ/dC ≤ ‖Π q‖_Λ ≤ 17.64 (q is the rotation mode; dp = M⁻ᵀ⟨χq, f⟩).
  - Both parts come from prof_tight's pointwise tip-end widths (full widths, θ, p and q: 4.7–5.3e-10 at x = 2.25, rising to θ 1.76e-9 and q 2.69e-9 at the tip). The bound has to drop by a factor of 1.64; if it scales with these widths, θ widths ≲ 2e-10 would give ≈ 1e-8 (an estimate, not run). The S1 agent's tighter tip-end profile (`/home/claude/s1_tight/ext_tip`) has not landed, and only a float diagnostic is there; rerunning on it means swapping the table in `s1chart2.py` and rerunning `h1.py m` (about 1 min).
  - The W′ τ-Lipschitz term, 7.3e-5·sup‖w_τ‖_Λ, is still to be added. It is expected to be ≤ 1e-11 with a coherent w_τ, but this is not certified.
- **G4. Lemma N3: lower-order terms and 𝒬 (written lemma, with T8).** Its numerical part is in the run.
- Closed this round:
  - **G2**: arb base windows, shells and K4 at W′.
  - **G5**: item 5 is now a scripted arb composition (`p2/t9_far.py`).
  - **KCH** (arb).
  - **gap 9**: Lemma 4a.
  - **gap 10**: SE far constants.
  - **the rounding bug**: se_corr hulls via `nextafter`; outward U/ub in all p2 scripts; se_corr, cconst/merge/camp and the switch were rerun, with the switch unchanged.

## 4. Files and commands (nice 10)
```
cd /home/claude/s1spec/step4/p2
python3 -I se_corr.py; python3 -I switch.py ../t3/p1/ckpt/slab_00096.json 2 0.873067570143915
T6RES=$PWD/t6/residuals.json SLAB0=2 YOUT=$PWD/out/ycomp_fin.json python3 -I ycomp.py
cd t8se; TOPLOSS=0.0196 python3 -I cconst.py --slabs 2-180 --out constants/slabs_fin; python3 -I merge.py --in constants/slabs_fin --first 2 --last 180 --out constants/constants_fin.json
SUBW=1.46484375e-4 T9JSON=../out/t9_cone_kch.json python3 -I camp.py --consts constants/constants_fin.json --F0 4e-23,3e-22,5e-19,3.4e-15 --out constants/amp_fin.json
cd ..; python3 -I boxes.py t8se/constants/amp_fin.json t8se/constants/constants_fin.json 1e-12; python3 -I t9.py t8se/constants/amp_fin.json
python3 -I t9_far.py t8se/constants/amp_fin.json t8se/constants/constants_fin.json   # item 5 -> out/t9_far.json
cd t9h; COHERENT=1 python3 -I h1.py m   # H1: IFT (alpha) + Lambda, tip offset as one scalar -> out/t9_h1_m_coh.json
T9IFT=../out/t9_ift_tight.json WINOUT=../out/t9_win_tight.json python3 -I win.py      # windows/shells at W'
```

## 5. Theorem 4 (draft)
**Theorem 4.** Let Γ0 be the closed curve X(α) = X0 + (L0/2π)∫₀^α e^{iθ}, with θ = −π/2 + α + Σ_{n≤3762} b_n sin nα, and with the exact data of `t2/gamma0_final.json`:
- b_n exact doubles for n ≠ 2;
- b₂ the unique closing value in its certified interval;
- L0 = 0x1.9e9c2d62f0874p+6.

Let γ be the curve diffusion flow from Γ0. Then there are t₁ ∈ [4.5978580262, 4.5978580263] (time since Γ0), λ₁ and c₁ such that:
- (a) γ is smooth and axially symmetric on [0, t₁], and k_max(γ(t₁)) ∈ [2.36e4, 2.38e4] **[certified]**;
- (b) Γ₁ = (γ(t₁) − c₁e₁)/λ₁ satisfies the step-3 hand-off hypotheses:
  - (H1) |α| ≤ 1.5e-4 **[certified]**, and √Λ ≤ 5e-8 **[open: G1, bound 8.18e-8]**;
  - (H2) **[certified]**;
  - (H3) **[certified]**;
  - (H4) **[certified]**;
  - (H5) **[certified, scripted]**;
  - the nonlinear level-3 terms use Lemma N3 **[G4, T8]**.

Hence, by step 3, γ develops a finite-time singularity modelled on S1. **This is conditional on G1 (√Λ) and G4 (T8's N3 lower-order terms and 𝒬).**
