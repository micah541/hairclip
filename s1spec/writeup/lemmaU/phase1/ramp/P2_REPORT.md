# P2 — Lemma-U inner-ramp constants on y ∈ [28, 56], certified (Nov 13, 2026)

Scripts: `ulib.py` (library), `ramp_K_cert.py` (driver). Output: `ramp_K_cert.json` (contains an `args` block with argv, env and input sha256), log `ramp_K_cert.out`.
Rerun: `cd lemmaU/phase1/ramp && PYTHONHASHSEED=0 nice -n 19 python3 -B ramp_K_cert.py --maxeval 150000 --tol 0.003 --fresh`. It is a single process with peak RSS 80 MB and takes 1:48 min. It resumes from `ramp_K_cert.partial.json`.
`run2_g450k/` is the same run with g ∈ [3e5, 4.5e5] (superseded).

Labels:
- **[cert]**: arb, outward. Each value is the maximum of arb upper endpoints over a branch-and-bound cover of the domain. Every quoted sup has a certified lower bound within 0.3% (`bb` block in the JSON).
- **[est]**: float arithmetic on certified numbers and t1's inputs.
- **[paper]**: hand argument.

## Setting
- **Weight:** Θ = W(σ)(1 + g e^{−4s}), with W = 1/(e^{−μ(σ+a)} + 1 + e^{ν(σ−ln1.5−a)}) and (μ, ν) = (4, 4).
- **Domain:** a ∈ [0, 0.3] and g ∈ [3.0e5, 6.0e5] are both branch-and-bound dimensions, so the bounds hold for every value in these ranges and no monotonicity is used. The centres satisfy s_ρ ≥ ln 56 (σ = (t−1)ln2 − w, w ∈ [0, ∞), with an exact σ → −∞ tail cell).
- **Ramp:** η = S(t) with t = (s − ln 28)/ln 2.
- **Why g reaches 6e5:** with my global constants the rate fixed point gives Ĝ₂₂ ≈ 2.77e5, so g ≈ 5.0e5. That is above the 4.5e5 I was given (see "Problems").
- **Tables:** `lemmaS/rigor/tables_exact.pkl` (exact Rational |coef|, read only).
- **A-coefficient bounds:** the symbolic Aterms port of `weights_exact`, with U and box as parameters.
  - Ramp pair coefficients use Ab = the Z-free (jet-free) part, including the pᵏU parts.
  - The ramp box is (M0, M1) = (0.11, 0.52), or (A)'s core box (0.0168, 0.108) for variant `smallA`.
  - U_ramp = 1/80 ≥ 1.3056·α₁·56^p = 0.012109 [cert arith; paper: |α| ≤ α₁, λ(τ) ≤ λ(τ_ρ)].
- **Global Y_{j+1}:** uses the crude box, M0z = 0.11 (sup |z₀| ≤ S0 ≤ 0.11) and U_glob = 0.11.

## Definitions (per unit εΘ, the same normalisation as `feas_P.json`)
- rem_m = Σ_{l=1..m} C(m,l) (η²)^{(l)} Θ^{(m−l)}/Θ. This is the cutoff part of (η²Θ)^{(m)}, enclosed per cell, with |·| taken per m.
- ρ̂ = ½ |(η²)″ + 2L(η²)′|, with L = (ln εΘ)′ = −4 + W′/W − 4θ.
- K_jk = sup [ Σ_pairs w_k Σ_{m≥1} |coef| Ab_(i,r) |rem_m| + [k = j] Y_{j+1} ρ̂ ], where w_k = 1 for the pair (k,k) and ½ for a mixed pair.

## Results [cert]

**Global weight constants** (support s ≥ ln28, centres ≥ ln56; lower bounds are the limits as σ → ±∞):

| constant | upper bound | lower bound |
|---|---|---|
| Γ₁ | 4.00008 | 4 |
| Γ₂ | 16.046 | 16 |
| Γ₃ | 64.183 | 64 |
| κ̃ | 8.00008 | 8 |
| κ̃₂ | 64.181 | 64 |

**On the ramp:**
- sup |W′/W| = 4.
- sup |Θ^{(m)}/Θ| = 3.882, 15.53, 62.13.
- sup |L| = 3.521.
- ρ̂ ≤ 16.78.
- **C_η^Ξ = sup |Lηη′ + η′² + ηη″| ≤ 13.349.** (A)'s power-weight value is 15.84.

**K table.** Primary entry `["4.0"]["K"]`: Lemma-U weight, crude box, Y_{j+1} from the full general-weight Ĝ (L⁴ extras, dfrac 2/25, ramp-foot Young adds), with the cubic cutoff classes folded in by box sups.

| j \ k | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| 0 | 7930.1 | 724.96 | | |
| 1 | 0.0815 | 5952.8 | 545.27 | **0.5716** (new; z₃ data X₃) |
| 2 | 0.872 | 528.68 | 20960.4 | 555.19 |

**Other variants**, under `variants.*.K`:

| variant | K₀₀ | K₁₁ | K₂₂ | K_{j,j+1} (j = 0, 1, 2) |
|---|---|---|---|---|
| U, pair part only (no ρ̂·Y) | 549.7 | 449.0 | 1894.6 | 725.0, 532.0, 537.2 |
| U, Y from L⁴ extras only (dfrac 1/10, t8 convention) | 7930 | 5940 | 20705 | same |
| U, Y without L⁴ | 7930 | 5170 | 11545 | same |
| U, small (A) core box on the ramp, full Y | 7922 | 5896 | 20423 | 707.9, 523.1, 495.2 |
| power weight, honest \|coef\|·4^m, Y without L⁴ | 4980 | 3206 | 6673 | 789.6, 586.8, 587.6 |
| feas_P reproduced (signed γ = −4) | 677 | 907 | 2869 | 789.6, 586.8, 587.6 |

- `K_t1_shape_no_cubic` is the same as the primary entry but without the cubic folding.
- The feas_P reproduction validates the pipeline: it matches `feas_P.json`, whose 300-cell naive grid sits 2–40% above my certified sups.

**Global Ĝ cross-check for agent A** (not my deliverable; my conventions are M0z = U = 0.11 and certified Γ):
- Y₃ = 641.4 without L⁴, 1152.9 with L⁴, 1178.5 with the foot adds.
- Ĝ₂₂ = 9.19e4 without L⁴, **2.658e5** with L⁴ (dfrac 1/10), **2.769e5** with dfrac 2/25 and the foot adds.
- t8's estimate was 2.398e5. At 2.769e5: rate margin 28⁴/Ĝ₂₂ = ×2.22, θ_a = 0.45, Gw = 1.82.
- The dissipation budget holds at every level; the slack at level 2 is 0.

**Foot / S33 weights.** Port of `weights_exact.ramp` with γ^m → Γ_m and κ → κ̃, crude box, dfrac 2/25, η′_max = 1.875/ln2:

| level | footw | W33 | W33·C_η^Ξ/3 | Young adds to the global Y |
|---|---|---|---|---|
| 0 | 0 | | | |
| 1 | 1.9407 | | | Y₂ += 1.94 |
| 2 | 2068.74 | 6.0492 | 26.92 | Y₂ += 47.91, Y₃ += 18.35 |

- These must enter Ĝ; they are included in the Ĝ₂₂ = 2.769e5 above.
- The level-2 budget is used exactly: 0.5 a_min. **The global Ĝ must then use dfrac = 2/25 for its z₄ Young steps, not 1/10.**
- t1/t8 used (A)'s 1.5707 and 2862.5 times the extras ratio 2.75 (i.e. 4.32 and 7872). They were ×2.2 and ×3.8 conservative.

## Comparison and effective KS [est: float arithmetic on certified K and t1's CLEAN-GO inputs]
Effective KS = t1's ramp bracket Σ K_jk ε(28) E_k + K X computed with the honest K, divided by the same bracket with the feas_P K that t1 uses.

| variant | level 0 | level 1 | level 2 |
|---|---|---|---|
| primary (with cubic folding and K₁₃) | 0.874 | **1.331** | 0.965 |
| without cubic folding | 0.874 | 0.788 | 0.937 |

Per-entry ratios, honest against feas_P (primary):
- **K_{j,j+1}: ×0.86, ×0.75, ×0.79.** The damping lowers them and feas_P's grid overestimated.
- **K_jj: ×9.9, ×4.7, ×6.8.** The cause is error 2 plus the L⁴ and foot parts of Y.
- **K_j0: K₁₀ = 0.082, K₂₀ = 0.87, against about 1e-12 in feas_P.** feas_P used U = 1e-12; with U_ramp they are negligible in the bracket.

So GONOGO's claim holds: the K_{j,j+1} dominate, and the K_jj change has no material effect on B₂ (the honest bracket is 3.5% *below* the one used). The level-1 rise to ×1.33 comes from the new K₁₃·X₃ entry; level 1 is not binding in T1. **The worst effective KS is 1.33, against the flip at KS ≈ 14.8: margin ×11.** T1's KS = 1/2/4 rows remain conservative.

## Problems / items for the coordinator
1. **New cubic cutoff-commutator classes**: rem_m × (jet part of A_i^{(r)}). Neither feas_P nor t1 has them. Raw weights per εΘ are under `cubic_cutoff_classes_crude`:
   - Level 1: (1,1,2) = 46.5. The others are about 2.2 or less.
   - Level 2: (1,2,2) = 458.8, (1,2,3) = 64.7, (2,2,2) = 72.3, (2,2,3) = 7.33.
   - The classes with a low index are folded into K via |z₀| ≤ 0.11 and |z₁| ≤ 0.52. This adds the extra entry **K[1][3] = 0.5716**, which **t1/t5 as written do not read**: extend the bracket loop to k = j+2, with data X₃.
   - Still open (genuinely cubic, level 2): **(2,2,2), weight 72.3**. It has the same form as the exact S33 foot term ∫∫εΘ|z₂|³, so add it to the sup-free S33 data: the coefficient goes from 26.9 to 99.2 (×3.7). t1 had 1.17e4·rr ≙ 41, so T1's S33 term becomes about 2.6e-4 [est], still negligible. **(2,2,3), weight 7.33** (εΘ z₂²|z₃|) still needs a bound, for example Cauchy–Schwarz with the cutoff-free L⁴ quantity plus z₃ from (A)'s dissipation. It is not done.
2. **Possible omission in (A) v3.1, flagged and not established.** (A)'s certified ramp commutators use the jet-free Ab (box M2.. = 0): `lemmaS_F2.py:138`, `rigor/ramp_arb.py:76`. Its ramp L⁴ classes carry only the core weight factor. I find no place where the same cutoff × jet terms on (A)'s ramp [14, 28] are bounded. Please check before relying on this. I did not edit anything in lemmaS*.
3. **Range of g.** My Ĝ₂₂ fixed point (≈2.77e5, so g ≈ 5.0e5) exceeds the given [3.0e5, 4.5e5]. I therefore certified on [3.0e5, 6.0e5]:
   - the pair parts peak at g = 3e5;
   - ρ̂ rises from 16.56 to 16.78 compared with the 4.5e5 run.

   Agent A's g ≈ 3.93e5 is inside. **If agent A's certified Ĝ₂₂ comes out near 2.7e5 rather than 2.4e5, Gw ≈ 1.82 (not 1.64) should go into T1.** The difference is my conventions M0z = 0.11 and U = 0.11, against t8's 0.02 and 0.05, plus the dfrac 2/25 and foot adds.
4. **Conventions to reconcile with agent A:** M0z, U_glob and dfrac 2/25 (forced by the ramp S33 foot Young). The same Y_{j+1} must be used in Ĝ and in K_jj. The JSON stores Y without L⁴, with L⁴, and with the foot adds; K_jj is linear in Y with slope ρ̂ ≤ 16.78, so it can be rescaled conservatively.
5. **Scope.** Only the pair/ρ̂ structure that t1 consumes was certified. The space–time conversion (IW, Gw) and the inflow term are t1's. The S33 normalisation (εΘ)_min^{−1/2}·F_max^{1/2} is left to the assembly: on the ramp, ln(εΘ) varies by at most ln2·sup|L| ≤ 2.44, so the max/min ratio is ≤ 11.5.
