# Lemma U, Phase 2a, agent B: ramp inequality, transport identity, Ĝ design (Nov 13, 2026)

Labels:
- **[cert]**: arb, `ghat_design_arb.py`.
- **[est]**: float model (`U_model2.py`, `top_pair_est.py`) on certified inputs.
- **[paper]**: proved by hand in `U_ramp.tex`. That file is **not compiled** (no LaTeX here); braces, `$`, `\[ \]` and environments are balance-checked by script.
- **[rev]**: the Phase-1 reviewer's numbers.

Nothing outside this directory was edited. The scripts write no bytecode.

## Files and reruns
All reruns use `PYTHONHASHSEED=0 nice -n 19 python3 -B <script>`, one process at a time, each under 110 MB RSS.

| file | content | time |
|---|---|---|
| `U_ramp.tex` | the derivations: Lemma T (transport), Prop. U (time-integrated ramp inequality), Lemmas D, Y, C, the Ĝ design, the model table, open items | |
| `ghat_design_arb.py` → `.json` / `.out` | Ĝ_jj per box in arb, with the design quantities (g, D(56), D(28), rate margins, c_rs, Γ validity) | 2 s |
| `U_model2.py` → `.json` / `.out` | float model. Self-check against `t1_checkpoint.json` (asserted); `--inputs A.json` reads agent A's E^A and D_k | 3 min |
| `top_pair_est.py` → `.json` / `.out` | [est] split of the z_{j+1} column; IBP + Young route for the top pair | 10 s |

Every JSON has an `args` block (argv, env, input sha256).

## 1. Inequality (Prop. U, `U_ramp.tex` §3)
For every centre and every τ ≥ τ_ρ, at levels j = 0, 1, 2:

𝔈_j(τ) ≤ ( √(𝔈_j(τ_ρ) + 𝓘_j + 𝓡^E_j + 𝓡^X_j + 𝓡^rs_j + 𝓕_j + 𝓒_j + 𝓑_j + Δ^R_j) + Φ_j )²

The terms:
- **inflow** 𝓘_j = c⁺ Σ_in E^k_j;
- **ramp energy terms** 𝓡^E_j = Σ_{k≤j} K_jk ε₂₈ Σ_E E^k_k;
- **X-terms** 𝓡^X_j = Σ_{k>j} K_jk 𝓧_k, where:
  - 𝓧₁ = ε₂₈ Σ_E E^k₁;
  - 𝓧₂ = min(energy route, √(S_X D₀) Minkowski route);
  - 𝓧₃ = (√(S_X D₁) + δU part)²;
- **ramp share of the rate** 𝓡^rs (only in the data56 design);
- **foot terms and S33** 𝓕_j;
- **couplings** 𝓒_j, from Ĝ_jk with k < j (ramp as data, bulk via the level-k damping surplus);
- **b-coupling** 𝓑₂ = b̄ 𝔅₁ / (c⁻K_b − 5b̄);
- **outer ramp** Δ^R (Lemma R);
- **frozen-profile forcing** Φ_j ≤ 4e-8 [est].

Data:
- E^k_j = (√E^A_j + √δU_j)²: (A)'s energies, sup in time;
- **D_k ≥ ∬_{τ_ρ}^∞ ∫_{28}^{56} εω (z^A_{k+2})²**: (A)'s exported dissipation **time integrals**, symbolic. No rate model remains.

The sweep constants are joint (y, τ) sweeps with the pointwise damping factor D(y):
- Σ_in, Σ_E and Σ_rs are time integrals of sup_ramp[profile · Ψ/ω];
- S_X = sup_{τ, ramp} Ψ/ω.

They replace the Phase-1 products G_w · I_W · sup(y/28) · sup(·).

Hypotheses:
- (R2): (1 − 4ā) g ≥ Ĝ₂₂ (1 + g/y_R⁴) φ, with y_R = 56;
- (R<2): (1 − 4ā) g > Ĝ_kk D(28) φ for k = 0, 1;
- P2's K, certified for g ∈ [3e5, 6e5];
- B3U.

**Top pair.** The (j+1, j+1) pair is **in** Lemma U's K_{j,j+1}: P2 `coef_arrays`, k = j+1, w = 1. So (A)'s gap does not recur here. What Lemma U shares with (A) is the data issue: z_{j+1} on the ramp is (A)'s dissipation. Only level 2 is affected (z₃ comes from D₁): z₁ and z₂ have energy routes.

Optional Lemma Y [est coefficients]:
- mixed pairs (a, 3) go to energies by exact IBP;
- the top pair (3,3) goes by one IBP plus Young against Lemma U's own ramp dissipation (δ a_min taken from B_c, so Ĝ₂₂ grows by at most 1.5/(1.5 − δ));
- this reduces the level-2 X-term from 0.141 to 0.005;
- it creates new cubic classes, from derivatives of the jet parts of A_i^{(r)}, which are **not bounded**.

**Cubic class (2,2,3).** It is put in sup-free S33 form by z₂²z₃ = (z₂³)′/3:
- c₂₂₃ ≤ 60.5 [est], against the placeholder 32.6;
- the jet part of c′ is a quartic term, bounded by the ramp L⁴ bound;
- (2,2,2) stays in c₃₃ (72.3 [cert]).

## 2. Honest X-term factor (Lemma D)
- **Export bound** [paper]: D_k ≤ (E_k(τ₁) + ∫ f_k⁺)/(κ_k a_min).
  - (A)'s certified E^A is H2 plus undiscounted impulse integrals, apart from a coupling term of 2.4e-7 in sup form.
  - So D_k ≤ c_D E^A_k/(κ_k a_min), with c_D ≈ 1 once agent A confirms. The model default c_D = 2 (the reviewer's "2E") is conservative.
- **Factor against the t1 rate model: f_X = 4.14 c_D/κ** [est], computed as 1.4713 c_D/(κ · I_W^in(56) · γ′).

| κ | c_D = 2 | c_D = 1 |
|---|---|---|
| 1.25 | ×6.6 | ×3.3 |
| (A) as it stands (0.3 / 0.2) | ×27.6 / 41.4 | ×13.8 / 20.7 |

- With pointwise D, the factor multiplying S_X is D(56)·φ = 1.036, where the model had G_w = 1.82.

## 3. Transport identity (Lemma T) [paper]
**Exact identity.**

∫∂_τΨ z_j² + 2∫Ψ z_j(−c z_{j+1}) = c∫(η²)′χ²WDφ z_j² − c∫Ψ(4θ + K_b e^{−s}) z_j²

Consequences:
- W and χ_out are exactly co-moving (ds_ρ/dτ = c = ¼ − a), so no W′ or χ′ term survives.
- a enters only through c: no ∫|a| and no sup|W′/W|.
- The D part is exactly the damping −(1 − 4a) g ∫εη²χ²Wφ z².

**Sign conditions:** ā < ¼, g ≥ 0, K_b ≥ 0, η nondecreasing. The inflow (≥ 0) is the only positive term.

**b-term.** The rate is ≤ b̄(C_j ∫e^{−s}Ψz_j² + couplings), with C = (5, 5, 8) and C₂₁ = 1. The reviewer's residual 4·b̄·sup|W′/W|/28 ≈ 1.25e-7 is the right size but the wrong argument:
- the relative rate b̄C/28 is constant in time on the ramp, and Lemma U has no decay;
- so Gronwall diverges on the infinite horizon.

The fix is a new factor in the weight: **φ = exp(K_b e^{−s}), K_b = 1e-3.** It works like D:
- it gives the decay c K_b e^{−s}, ×71 above the b-rate;
- it costs e^{K_b/28} − 1 = 3.6e-5 on the data, and Γ_m changes by ≤ 1e-4 relative;
- the level-2 b-coupling is 8.8e-4·𝔅₁.

The outer-ramp b-part is time-integrable: ≤ 2.5e-8 · sup ∫_outer z_j², added to Δ^R.

## 4. Ĝ design
**Recommendation: data56 with a global box (0.15, 0.52).**
- The ramp share of the rate is data, as in (A); the rate criterion starts at y = 56.
- g = Ĝ₂₂/(1 − 4ā − Ĝ₂₂/56⁴).

| box | Ĝ₂₂ [cert] | global28: G_w, rate margin | data56: D(56), D(28), rate margin |
|---|---|---|---|
| (0.11, 0.52), U .09, M0z .16 | 268125.9 (= reviewer's 268126) | 1.774, ×2.29 | 1.0280, 1.448, ×36.7 |
| **(0.15, 0.52)** | **338258.8** (×1.2616, reviewer: ×1.26) | 2.224, ×1.82 | **1.0356**, 1.570, **×29.07** |
| (0.20, 0.60) | 590852.8 | 25.8, ×1.04 | 1.064, 2.02, ×16.6 (Γ valid) |
| (0.30, 0.80) | 1.84e6 | infeasible | ×5.3, but θ(28) > 3/5 |
| (0.50, 1.00), Lemma-O-like | 6.04e6 | infeasible | ×1.63, D(28) = 26: unusable |

- **The reviewer's "Gw ≈ 1.03"** is D(56) = 1.0356 [cert].
- **What data56 costs:**
  - the data term 𝓡^rs₂ = 2.4e-4 (0.06% of 𝔅₂), with c_rs = Ĝg/28⁸ = 0.314 [cert];
  - no change to the dissipation budget, dfrac, Young steps or Ĝ pipeline;
  - P2's g range covers g = 3.50e5.
- **Main gain: robustness.** The rate margin goes from ×1.8 to ×29, so the box and U conventions stop mattering.
- **Most of the G_w gain does not need data56.** It comes from accounting for D pointwise, which is valid in either design: ladder row 6, global28, gives ×2.39.
- **Localised Ĝ** (paper derivation, §4(b)):
  - partition ζ_I² + ζ_III² = 1, physical, so co-moving with no transport term;
  - extra terms in the L⁴ step, the {2,3,3} IBP and the interpolation, through C_ζ (≈ 80 for a transition of length ln 1.5; Ĝ₂₂ about +35% [est]);
  - zone III by Gronwall, factor ≤ 1.005.
  - Needed only if D5a cannot give a box of (0.20, 0.60) in zone III. A global box up to (0.20, 0.60) costs nothing material under data56.

## 5. U_model2 margins [est], default inputs
Default inputs:
- E^A = certified IBP E^A × (1.41, 1.13, 1.25) = (1.708e-7, 9.237e-5, 0.06006);
- D_k = 2E^A_k/(1.25 a_min), i.e. D₀ = 3.08e-7, D₁ = 1.67e-4;
- certified K, footw and c₃₃; weight (4,4,0.3).

**Self-check: passed, relative deviation 0.** With every new ingredient off, the model reproduces the t1_checkpoint late row (M1 ×2.0189, M0 ×3.834 at KS = 1; M1 ×1.9316 at KS = 2) and the present rows.

| configuration | late M1 (KS 1 / 2) | late M0 | present M1 | present M0 (y_c = 126.6) |
|---|---|---|---|---|
| checkpoint (row 0) | 2.019 / 1.932 | 3.83 | 2.41 | 1.96 |
| **headline: data56, box (0.15, 0.52)** | **2.40 / 2.09** | 5.83 (vs 0.15; 4.3 vs 0.11) | 2.82 / 2.35 | 2.85 (vs 0.15; 2.09 vs 0.11) |
| + Lemma Y, δ = 0.25 | 2.65 / 2.46 | 5.82 | 3.35 / 3.05 | 2.85 |
| D_k with (A) as it stands (κ = 0.3/0.2) | 1.50 | 5.43 | 1.59 | 2.80 |
| D_k with c_D = 1 | 2.62 | 5.89 | 3.19 | 2.85 |
| E^A₂ repair ×1.65 / ×2.5 | 2.31 / 2.16 | | 2.75 / 2.63 | |

- **Late 𝔅₂ in the headline is 0.376:**
  - launch 0.142, inflow 0.076, X-terms 0.141 (all from K₂₃𝓧₃), outer 0.016;
  - ramp 6e-4, rs 2.4e-4, S33/foot 1.8e-4, couplings and b ≤ 2e-6.
- **Ladder** (`.out`): the E^A repair alone gives ×1.86. Adding the honest X-terms with a rate-free accounting but without pointwise D gives ×1.68 to ×1.85. Pointwise D brings it to ×2.39.
- **Kill check 2** (late M1 ≥ ×1.3 at honest KS) is met in every row. The weakest is (A) left with κ = (0.3, 0.2), at ×1.50.

## 6. Open items
1. **Agent A's export:** D₀ and D₁ (with κ, and the split confirming c_D ≈ 1), and the repaired E^A. Rerun with `python3 -B U_model2.py --inputs A.json --tag A1`.
2. **Certify in arb:** the joint sweeps Σ_in, Σ_E, S_X, Σ_rs and the launch constants with pointwise D. These are explicit sups and integrals.
3. **New cubic classes** from the IBPs of Lemma Y and Lemma C (derivatives of the jet parts of A_i^{(r)}). Also c₂₂₃ and the Young coefficients in arb.
4. **Zone III box from D5a** in the S₁ chart. Data56 needs it ≤ (0.15, 0.52), or ≤ (0.20, 0.60) as fallback; otherwise the localised Ĝ constants.
5. **φ in P1/P2:** Γ_m and K change by ≤ 1e-4 relative; rerun them with the φ factor.
6. **P2's K** must be rerun if a box with g outside [3e5, 6e5] is chosen.
7. **Φ_j in arb** (4e-8 [est]).
8. **Still typed in the model:** H2, H3, DT6, ā, b̄, α₁ = 1.5e-4, and the reviewer's repair factors.
