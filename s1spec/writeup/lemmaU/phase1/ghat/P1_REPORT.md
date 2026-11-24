# P1 report: certified weight constants and Ĝ_jj for the Lemma-U weight (Nov 13, 2026)

Labels: **[cert]** = python-flint arb with outward rounding, no sampling. **[est]** = float. **[cond]** = certified arithmetic for a step whose derivation is missing.

Files in this directory:
- `weight_consts_arb.py` → `weight_consts.{json,out}`: Γ_m, κ̃, κ̃₂.
- `ghat_arb.py` → `ghat_cert.{json}`, `ghat_arb.out`: Ĝ, the fixed point, the sanity check and the sensitivities.
- `arbio_local.py`: a verbatim copy of `lemmaS/rigor/arbio.py`.

Rerun with `PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 nice -n 19 python3 <script>`. Each run is one process, takes under 3 s and stays under 100 MB. A rerun gives `ghat_cert.json` byte-identical apart from `wall_s`, and the JSON has an `args` block with input sha256s. Nothing outside `ghat/` was touched.

## 1. Weight constants [cert]

Geometry: s ≥ ln 28 and s_ρ ≥ ln y_L, with y_L = 56 or 150, and all a ∈ [0, 0.3].

Method:
- W^(k)/W is written as an exact integer polynomial P_k(p, q), with p = A/S and q = B/S. p is strictly decreasing in σ and q strictly increasing.
- The damping factor enters through D^(j)/D = (−4)^j θ, with θ ∈ [0, θ_max(σ)] and θ_max(σ) = g/(max(28, y_L e^σ)^4 + g). Every point of this (σ, θ) set is realised by the geometry, and no other point is.
- Γ_m, L and L² + L′ are all **affine in θ**; the θ² terms of L² + L′ cancel exactly. So the supremum over θ is taken at θ = 0 and at θ = θ_max.
- σ ∈ [−12, 12]: adaptive arb cells in (σ, a), about 3000 leaves.
- σ ≤ −12 and σ ≥ 12: by monotonicity, the (p, q) box on each tail is tiny, and one arb evaluation covers the whole tail. This includes σ → ±∞ and s → ∞.
- Each bound is monotone in θ_max, so the bounds computed at θ_a = 3/5 hold for every θ_a ≤ 3/5.

| weight | Γ₁, Γ₂, Γ₃ | κ̃ | κ̃₂ |
|---|---|---|---|
| (4,4,a), all a ∈ [0, 0.3], y_L = 56 and 150 | ≤ 4.0000000000000009, 16.000000000000004, 64.000000000000014 | ≤ 8.0000000000000018 | ≤ 64.000000000000014 |
| (4,2,0.5), y_L = 56 and 150 | same | ≤ 6.0000000000000009 | ≤ 36.00000000093 |

- **The float values hold:** Γ_m = 4^m, κ̃ = ν + 4 and κ̃₂ = κ̃², up to a relative 1e-12 of cell overestimate. They are suprema, approached only as σ → ±∞.
- They do not depend on a, on θ_a (up to 0.6) or on y_L. So the fixed point in g needs no iteration.
- Γ_m is used per order m; the coarser Γ_eff^m is not needed, though here it would give the same numbers.

## 2. Ĝ_jj with the general-weight L⁴ classes [cert]

**What was ported.** This is a line-by-line copy of `weights_exact.py`'s `core` path (ab_exact, At, classify, l4_extra_YC, budget, ghat), with three changes:
- γ^m → Γ_m;
- κ → κ̃ in the `L4` class and in the `S33` IBP class (κ² there is κ̃·κ̃);
- κ² → κ̃₂ in the interpolation step.

**Sign cancellations.** None are assumed. The exact tables have every |coef| ≥ 0 and weight-derivative orders m ≤ 3, both asserted. The tables keep one rational coefficient per monomial (m, i, r), and that coefficient is weight-independent (`om` is a generic sympy Function). So the bound |c·ω^(m)·A| ≤ |c|·Γ_m·Θ·|A| is applied term by term.

**The {2,3,3} class.** It uses `ibp` (no sup of z₂). The level-2 dissipation budget is used exactly (slack 0, an exact rational check), as stated in lemmaS.tex.

| box (M0, M1) | weight | Ĝ₀₀ | Ĝ₁₁ | Ĝ₂₂ | Ĝ₂₂ without L⁴ | L⁴ factor |
|---|---|---|---|---|---|---|
| crude (0.11, 0.52) | (4,4,a) | 49317.054 | 28556.252 | **239750.34** | 78448.6 | ×3.056 |
| fine (0.21, 0.21) | (4,4,a) | 50537.03 | 24412.23 | 109185.6 | 54663.3 | |
| (A) core (0.0167137, 0.1070809) | (4,4,a) | 45352.91 | 19642.31 | 43180.95 | 40056.0 | |
| crude | (4,2,0.5) | 43012.33 | 23839.73 | 196498.0 | 70242.3 | |
| fine | (4,2,0.5) | 44304.27 | 20268.10 | 89645.86 | 48174.3 | |
| (A) core | (4,2,0.5) | 39294.05 | 16015.30 | 36898.52 | 34463.2 | |

All entries are certified upper bounds; the ball radii are about 1e-11. The off-diagonal Ĝ_jk are in the JSON.

## 3. Certified fixed point (θ_a := Ĝ₂₂/28⁴ rounded up to 1e-9; check: 28⁴θ_a ≥ Ĝ₂₂)

| case | θ_a | g | Gw = 1/(1−θ_a) | rate margin 28⁴/Ĝ₂₂ |
|---|---|---|---|---|
| **crude, (4,4,a)** | **0.390056123** | 3.931e5 | **≤ 1.639496** | **≥ 2.56373** |
| crude, (4,2,0.5) | 0.319687791 | 2.888e5 | ≤ 1.469914 | ≥ 3.12805 |
| fine box, own fixed point, (4,4) / (4,2) | 0.17764 / 0.14585 | | 1.2161 / 1.1708 | 5.629 / 6.856 |

**Comparison with t8 [est].** t8 gave crude Ĝ₂₂ = 2.398e5, rate margin 2.564 and Gw 1.639. The certified numbers agree to four or more digits, so t8's float pipeline was accurate.

## 4. Which box governs the rate criterion at y = 28

**What the certified pipeline gives.** Ĝ is a single constant multiplying ∫εΘ z_j². The P-terms are pointwise, but the L⁴ identity, the {2,3,3} integration by parts and the interpolation N_{j+1}² ≤ N_jN_{j+2} + ½κ̃₂N_j² are global integrations over supp Ξ. So M0 and M1 must be sups over the whole of supp Ξ.

For **every** centre ρ, crude or fine, the Lemma-U weight does not vanish anywhere from y = 28 out to the outer ramp. On that set only the crude bootstrap box (0.11, 0.52) is available. **So the crude box governs: Ĝ₂₂ = 2.3975e5, Gw ≤ 1.6395, rate margin ≥ 2.5637.**

**The radius-dependent box: [cond], not certified as a Lemma-U step.**
- Combining (A)'s core box on [28, 60] with the crude box beyond 60 requires g ≥ max(Ĝ_small/(1 − Ĝ_small/28⁴), Ĝ_crude/(1 − Ĝ_crude/60⁴)).
- The **binding term is the crude box at y = 60**, not y = 28. The small-box term needs only g ≥ 4.6e4, against 2.44e5.
- Result: θ_a ≤ 0.28439, **Gw ≤ 1.3974** for (4,4); Gw ≤ 1.3246 for (4,2). This confirms GONOGO's "≈ 1.40".
- **What it still needs:** a localised Ĝ, i.e. a position-dependent box inside the L⁴, IBP and interpolation steps. That means a smooth partition or a weight Ξ·φ(s), with extra κ̃ terms from φ′. That derivation does not exist yet.
- It also inherits (A)'s conditionality: (A)'s core box (0.0167, 0.107) is certified arithmetic for a statement that assumes the crude box on (A)'s support (GONOGO error 6). It is usable inside the continuity argument, in the order (A) → U-crude, but not beyond.

**The fine box has a related problem.** The fine family's weight also reaches y = 28, so its box must hold on [28, 150) as well. T1's crude-family conclusion is y|f′| ≤ 0.30–0.36 [est], which is above 0.21. So the fine-box Ĝ (rate margin 5.6) is legitimate only if (0.21, 0.21) is established from y = 28 outward. Otherwise use the crude Gw = 1.6395, as GONOGO already does downstream.

## 5. Sanity check against the power weight, and the t8 discrepancy

With Γ = 1, κ = 5, κ₂ = 25 on the crude box and **s33 = `assumed`**, the pipeline reproduces `lemmaS/rigor/weights_exact.json`'s core Ĝ for all k ≤ j: Ĝ₂₂ = 65729.7217, maximum relative deviation 5e-59 [cert].

With **s33 = `ibp`** at the same power weight, Ĝ₂₂ = **77720.88** [cert]. That equals t8's 7.772e4.

So the discrepancy is not a t8 error. **The certified (A) value 65729.72 was produced with `--s33 assumed`, i.e. sup|z₂| ≤ M₂* = 3 on the core.** This is the F2.py default, and it is recorded in `lemmaS_F10p_M0.52.json` args `"s33": "assumed"` and in `rigor/weights_exact.py: core(..., s33='assumed')`.

**Possible inconsistency in closed v3.1 Lemma S(A), for the coordinator.**
- lemmaS.tex's hypothesis table says *"H-L4: core cubic/quartic terms incl. {2,3,3} — closed by the L⁴ identity and one integration by parts"*, i.e. the ibp route. The final certified run did not use it.
- With ibp, the core Ĝ₂₂ becomes 7.772e4: the margin against 1.54e5 goes from 2.34 to about 1.98, and γ′₂ falls slightly, about 0.143 → 0.124 [est]. That changes E (A) and E^A of the C1 rerun a little.
- I did not touch any lemmaS file. Lemma U itself uses ibp throughout, with no sup|z₂| assumption.

## 6. Inherited hypotheses and sensitivities [cert]

- **M0z = 1/50** (sup|z₀| = sup|v| as a coefficient), inherited from (A)'s H-B2. Lemma U must supply it on supp Ξ.
  - Raising it to 0.11 or 0.13 changes Ĝ₂₂ only to 239773 / 239779 (rate margin 2.5635). It is negligible.
- **Umax = 1/20** (|U_j| ≤ p^j Umax in the coefficient jets), inherited from (A). It is **not** valid out to 380λ₁, where U ≈ 0.086 (GONOGO error 5).
  - With Umax = 0.09: Ĝ₂₂ = 250918, rate margin ≥ **2.4496**, Gw about 1.69.
  - Lemma U should use the Umax of its actual support; it costs about 5%.
- **Unchanged from (A):** p = 1.0246, dfrac = 0.1, the ±Young splitting, and the box structure V₂…V₇ → U-part only.
- **Not covered here:**
  - ramp constants and C_η^Ξ (P2);
  - the ∂_τΘ / transport term of the moving weight. GONOGO assumes it gives a rate ≤ 0 through the damping factor, and the energy identity for Lemma U must show it.

## 7. Effect on GONOGO

- **Risk 1 is retired at the estimated level.** The certified crude Ĝ₂₂ = 2.3975e5 equals the t8 estimate (2.398e5), not 18% above it.
- The M1-flip threshold (Ĝ₂₂ ≥ about 4.4e5) is ×1.84 away, and the rate margin is ×2.564 [cert].
- The radius-dependent box would bring Gw to 1.397, but it needs a new localisation derivation, and its binding term sits at y = 60, not 28.
