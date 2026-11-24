# FLUX2_DESIGN: flux lemma v2 for the redefined outer-ramp functional K^OR (levels 0 and 2)

Fix round 7, unit "flux2 design". Nov 20, 2026. Fable (design lead). Approved by Micah on Nov 20.
**This is a design. Nothing here is implemented, certified or closed.** The proof tree was read-only; everything I wrote is in this directory.

Labels:
- **[cert]**: certified number quoted from an existing output (file named), or computed here in arb with outward rounding (`probe/`);
- **[read]**: verified by reading the cited file:line;
- **[sym]**: exact sympy identity (`probe/id_probe.py`, or the file named);
- **[hand]**: my derivation by hand, not machine-checked;
- **[est]**: estimate, not a bound.

Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `V` = `/home/claude/s1spec/step3nl_v4b/verifier` (tag fin13v4b), `SPEC` = `B/fix5/ei_spec/EI_SPEC.md`, `F1` = `B/fix3/flux/`, `ROUTE` = `B/fix3/ROUTE_DECISION.md`, `N36` = `B/fix6_v36` (the running (A) chain v36; I did not write there).
Companion file: `flux2_lemma.tex` (LaTeX-ready statements and proof sketches; not compiled, no LaTeX here).

## 0. Summary

1. **Verdict: feasible, with large room; no change of FLUX0/FLUX2 expected** [est].
   - Pointwise part: K0 + eps_I·K2^pt ≈ 3.93e-21 against FLUX0 + eps_I·FLUX2 = 1.64024e-20, margin **≈ ×4.17** [est from cert inputs, v35 data; `probe/feas_probe.out`]. It is the level-0 margin of F1 (×4.20 [cert]) less about 0.7% for the new terms.
   - Time-integrated part: 2C′ ≈ 2.3e-25 against the binding limit 8.2e-23, **≈ ×360** [est]; ×39 if the level-2 export of (A) is ten times worse than I guess. S′ ≈ 6e-29 against the slack s ≈ 1.0e-21.
2. **One change of route against ROUTE §1, forced by the redefinition of K^OR.** The "trace" form (endpoint term through g‴(57) by Agmon) does **not** close the new K2^OR: the term ζ·N_low′ contains d110·(g″² + g′g‴), and after any integration by parts that removes g⁗ one is left with ∫W ζ² d110 g″ g‴². With g‴ controlled only in L²(dτ) this term is not integrable in time [hand, §6.2]. So the design uses the level-2 export of (A) **as a space-time L² bound of g⁗ on [57, 60]** (N4), not as a trace.
   - Gain: at level 2 **no integration by parts at all**, hence **no endpoint terms**; the septic cutoff enters only through sup|ζ^(i)|, i ≤ 4.
   - Cost: the chart identity at level 4 on the shell σ ∈ [54.75, 57.76] and D⁴d there (new stage X3, about 4–5 agent-hours). All chart coefficients there carry factors ≤ 1e-40 [est from the cert level ≤ 3 values, which are ≤ 6.6e-14 on the whole core and far smaller on the shell].
   - What (A) must supply is unchanged: the level-2 window export with r_2 = 1/1000, which `N36` is producing now (§12).
3. **Level 0 stays pointwise in τ** and keeps F1's exact reduction. New: the quasilinear factor j = J⁻⁴ inside the commutator, the modulation transport inside the drift ball, the lower-order nonlinear term, three formal items.
4. **The verifier needs no rerun.** The absorption lemma (§7) shows that fin13v4b's certificates, unchanged, give Λ ≤ L(τ) + 2C′. Conditions C0–C9 are to be checked by a script on the stored rate tables.
5. **Plan**: 9 stages, **about 20–22 agent-hours (Opus) plus one Fable review**; compute under 15 min; one process, nice 19, < 400 MB.
6. **Not yet justified** (§11): the level-4 chart structure [hand]; that the v36 level-2 export has the window-affine form (the v36 agent says so in `export_t2.py`; its report is not written yet); the verifier-side reading "rate cell = sup of tot/Λ over Λ ∈ [Ly_lo, Ly_hi]" on which C2/C3 rest [read, not re-derived]; the regularity lemma U-1 (open since SPEC §3).

## 1. What fin13v4b consumes [read]

- `V/run_v4b.sh:8`: FLUX0 = 1.64e-20, FLUX2 = 8.0e-20. `V/constants.py:65` → ABS_FLUX_I0, ABS_FLUX_I2.
- `V/rhs_iv.py:229`: `L = L + MU·FLUX0 + MUEPSI·FLUX2`, only when floors are on. One constant in the degree-0 group of every rate cell.
- `V/closure.py:175`: tail floor F0 = mu·(FLUX0 + eps_I·FLUX2) = 1.312192e-21 [cert, rate table `floor_terms`]; L_F = 2F0/rho.
- `V/closure.py` `unconditional()`: slope of the J's contains F0; the march part does not (the v35/v36 chain adds F0·tau_march in `hjw_rep`, kappa_J).
- mu = 0.08, eps_I = 3e-5, so mu·eps_I = 2.4e-6 and eps_I·FLUX2 = 2.4e-24 = 1.5e-4 of FLUX0.
- **Consequence (as ROUTE §2.1): the verifier uses only the sum.** The hypothesis it needs is, for a.e. τ,

  mu·K0^OR(τ) + mu·eps_I·K2^OR(τ) ≤ F0 − s + φ(τ),   (I3′)

  with a slack s > 0 and a function φ ≥ 0 that is controlled on windows (§7). A separate bound K2 ≤ FLUX2 is not needed and is not claimed (the crude level-2 pointwise bound of §5 is 1.7e-19 > 8.0e-20 [est]; see §9.3).
- No other consumer of FLUX exists in `V/` (grep) [read].

## 2. The objects

### 2.1 Definitions (SPEC §4.10, rewritten) [hand, each step is a [sym] identity of `ei_identities.py` or one line]
OR = [57, 60]. ζ = S7((60 − x)/3), S7(t) = t⁴(35 − 84t + 70t² − 20t³). g = w′, v = ζg. W = e^{−ψ}, ψ′ = φ′ = x/(3(1 + x²)^{1/3}) on x ≥ 38.25.
D2v = (W v″)″/W = v⁗ − 2ψ′v‴ + ψ₂v″, ψ₂ := ψ′² − ψ″.
K_m g := [∂^m, ζ]g: K2g = 2ζ′g′ + ζ″g, K3g = 3ζ′g″ + 3ζ″g′ + ζ‴g, K4g = 4ζ′g‴ + 6ζ″g″ + 4ζ‴g′ + ζ⁗g.

Write **j := 1 − c = J⁻⁴**, J² = (1 − kU0)² + U1². Then −∂⁴ + (c∂³)′ = −(j∂³)′, and

  −[L1, ζ]g − cK4g − c′K3g = **j·K4g + j′·K3g** + k²·K2g + C·ζ′g,  C = 7km + q/4

(F1's H̃ is the case j ≡ 1). With the modulation commutator T9, put **C_eff := C − a·Q_u − b·C_u**. Then

  X0 := j K4g + j′ K3g + k² K2g + C_eff ζ′g + ζ c0′ w + ζ(a E_A′ + b E_B′) w,
  X2 := X0 + ζ N_low′,

  **K0^OR = ∫_OR W v X0 dx − ∫_OR (W ζ v)′ N_low dx,  K2^OR = ∫_OR W D2v · X2 dx.**

N_low = d010 g″ + d100 g′ + d200 g′² + d110 g′g″ + d300 g′³ + d000 (in g-jets: w″ = g′, w‴ = g″), |d000| ≤ A0|w| + A1|g|.
N_low′ as in SPEC §4.5 (twelve monomials), with Dd = d^a + d^b·g′.

Check against the task list:
- "outer-ramp quasilinear commutator remainder" (ROUTE §3 item 1) = j K4g + j′ K3g − K4g. It is in X0, X2.
- "level-2 modulation block on the outer ramp" = the T9 term −ζ′(aQ_u + bC_u)g (in C_eff) and ζ(aE_A′ + bE_B′)w. T8, T10, T11 and the v-parts of T12 are global V-norm terms of the verifier (SPEC §4.9, §4.10) [read `V/rhs_iv.py:158-160, 179-182` through V4_REVIEW §1.2–1.3].
- The top-order quasilinear term (c v‴)′ is global and in the verifier (SPEC §4.3); it is **not** in K^OR.

### 2.2 Sign and size of j on OR
On σ ≥ 54.75: |k| ≤ 3.7e-44 [cert, F1 `flux_cert.json` profile], so J² = 1 + U1² up to 1e-40 and j ∈ (0, 1 + 1e-40]. The code will not use the sign: it takes the ball **j ∈ [1 − CBAR_OR, 1 + CBAR_OR]**, CBAR_OR ≤ 3.3655e-3 at the alpha-level 2.5e-4 [cert, `step3nl_v4b/v4cert/out/constants_v4.json` key `TQ[0.00025][ZB][c]`; zone ZB = [34.25, 60] ⊃ OR], and **|j′| = |c′| ≤ C1_OR := c_a + c_b·s1** with c_a ≤ 6.165e-5, c_b ≤ 0.16513 [cert, same file], s1 := sup_OR|g′| ≤ S1V/54.75 = 0.010959 (box B-e, S1V = 0.6). So C1_OR ≤ 1.872e-3 [cert arithmetic].
The production march uses alpha_up ≤ 1.6e-4 ≤ 2.5e-4 at every step [cert, rate table `steps[].alpha_up`, `ramp_levels`], so the 2.5e-4 tables apply.

### 2.3 Boxes used on OR (all are standing hypotheses of fin13v4b, asserted by `V/closure.py` `box_assert` [read :84-93])
| Box | Use here |
|---|---|
| B-a: alpha ≤ 2.5e-4 | coefficient tables |
| B-c: \|g\| ≤ 0.02 on [16.25, 60] | coefficient tables (U1 box) |
| B-d: \|w\| ≤ 0.5 =: w̄ on [34.25, 60] (asserted as far_w + 25.75·S0V ≤ 0.5) | terms with w on OR; U0 box |
| B-e: S0V = 0.0153, S1V = 0.6 | s0 := sup_OR\|g\| ≤ S0V; s1 := sup_OR\|g′\| ≤ S1V/y_min |
| sup\|(a, b)\| ≤ ā | C_eff ball; ā = 2.14041e-7 [cert, `runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json`] at the conclusion level. **Which level the bootstrap supplies (ā or a typed hypothesis ≥ ā) must be taken from the v36 continuation text**; the flux is insensitive (relative 1e-6) |

## 3. Data on [57, 60]

### 3.1 Notation
σ = x + σ_tip, σ_tip ∈ [−2.249127, −2.249019] [read F1 tex]. OR in σ is the shell S := [57 + σ_tip, 60 + σ_tip] ⊂ [54.7508, 57.7510]. Note 57.7510 > 57.75: the old code uses yw = 57.75 as the upper end (`assembly_arb.py:175-177`). That is covered there by the cruder factor 60/28 for 1/ω [hand], but **the v2 code must use the enclosure of S**, not 57.75.

  M_l(τ)² := ∫_OR (∂_x^l g)² dx,  N_l² := ∫_OR W (∂_x^l g)² dx ≤ W(57)·M_l²,  W(57) ≤ 1.4371630e-18 [cert, `lemmaS_cflux.json` `W57_upper`; to be recomputed in X2].
  s2(τ) := sup_OR|g″|,  ω_W := (∫_OR W dx)^{1/2} ≤ (W(57)/inf_OR ψ′)^{1/2}, inf ψ′ ≥ 1.2826 [cert, probe].

### 3.2 Three classes of data
| Quantity | Control | Source |
|---|---|---|
| M_0, M_1, M_2 | **pointwise in τ**, at the bootstrap level β_A = 2 | (A): Ê_j, the conversions of `assembly_arb.f4b` (N0, N1, N2), d by `b0_dprof`, chart by `n1_chart` (k ≤ 2), per level by Minkowski |
| s0, s1 | pointwise (boxes B-e) | verifier hypothesis; implied by (A) (I2) |
| ∫_I M_3² dτ ≤ m3c + m3s·\|I\| | **window-affine** | (A) level-1 export (C0_1, C1_1), chart k ≤ 3, D³d |
| ∫_I M_4² dτ ≤ m4c + m4s·\|I\| | **window-affine** | (A) **level-2 export (C0_2, C1_2), new in v36**; chart k = 4 and D⁴d on S, **new (stage X3)** |
| ∫_I s2² dτ | window-affine, derived | Agmon on [57, 60]: s2² ≤ M_2²/3 + 2·M_2·M_3 |

"Window-affine" always means: for every window I = [τ_a, τ_a + T] ⊂ [τ_1, τ′].

### 3.3 Conversions [hand]
With D = σ∂_σ: σ³∂³ = D³ − 3D² + 2D and σ⁴∂⁴ = D⁴ − 6D³ + 11D² − 6D, so with dσ = σ·ds and σ ≥ σ_min := inf S:

  M_3² ≤ 3·σ_min⁻⁵ ∫_S [(D³g)² + 9(D²g)² + 4(Dg)²] ds,
  M_4² ≤ 4·σ_min⁻⁷ ∫_S [(D⁴g)² + 36(D³g)² + 121(D²g)² + 36(Dg)²] ds.

The first is the rule behind `assembly_arb.py:178` [read]. The second is new.
D^k g = z^A_k∘Y − D^k d − χ_k (chart lemma (d), `B/fix1/b12_text/chart_lemma_body.tex`), so in L²(S × I; ds dτ), by Minkowski:

  ‖D^k g‖ ≤ (1 + ρ_c)^{1/2}·‖z^A_k‖_{L²(Y(S) × I)} + ‖D^k d‖ + ‖χ_k‖.

- ‖z^A_k‖² over the shell: for k = 3, 4, ≤ (y_max⁵/28)·(C0_{k−2} + C1_{k−2}·T) (the export has the weight εω = 28·y⁻⁵ ≥ 28·y_max⁻⁵ on Y(S), y_max = sup Y(S)); for k ≤ 2, ≤ T·(y_max/28)·β_A·Ê_k (pointwise energies; slope only).
- ‖D^k d‖²: a rate per unit τ (|α| ≤ α_1 at every τ): T·∫_S (D^k d)² ds. k ≤ 3 from `b0_dprof` (core table), k = 4 new.
- ‖χ_k‖: k ≤ 3 by `n1_chart` core values (χ^add as a rate; relative parts β_k2, β_k3 ≤ 2.8e-17 [cert], moved to the left side); k = 4 by §3.4.

### 3.4 Chart identity at level 4 on the shell (new; stage X3) [hand]
Chart lemma (d) holds for every k with the same proof: z^A_k∘Y = D^k(w′ + d) + χ_k, χ_k = (T^k − D^k)b + T^kΦ + p^k(U(σ) − U(Y)), T = (1 + ε_c)D. For k = 4, by weight counting (D^i b has weight i, total weight ≤ 4, and D^iε_c, D^iΦ contain D^i b linearly at top order):

  **χ_4 = A_4 + B_42·b₂ + B_43·b₃ + B_44·b₄ + B_422·b₂²,**  b_i := ∂_σ^i b = u^{(i+1)},

with coefficients depending on (σ, e, δ, k, …, k⁗, u, b, b₁) only. The only product of two jets of order ≥ 2 is b₂² (weights 2 + 2). **The script must assert this structure in sympy** (third derivatives in (b₂, b₃, b₄) vanish except ∂²/∂b₂²; no b₅; χ_4 = 0 for a ray), as `n1_chart.py:79-81` does for k ≤ 3.
Bound: with b_i = ∂^i g + P^{(i)} (P = αψ̂3′ + α²ĥ2′; P^{(4)} needs ψ̂3^{(5)}, ĥ2^{(5)}, available for x ≥ 28 from `c4fun.ext_jets`, orders ≤ 8 [read ROUTE §2.2]),

  ‖χ_4‖_{L²(S×I)} ≤ Ā_4·(|S|_s·T)^{1/2} + B̄_42(‖g″‖ + ‖P″‖) + B̄_43(‖g‴‖ + ‖P‴‖) + B̄_44(‖g⁗‖ + ‖P⁗‖) + B̄_422·(2‖g″²‖ + 2‖P″²‖),
  ‖g″²‖²_{L²(S×I)} = ∫_I∫ g″⁴ ≤ ∫_I s2²·M_2² dτ ≤ sup M_2²·∫_I s2².

The g⁗ term is circular with coefficient B̄_44·(conversion) ≲ 1e-40 [est]: solve (1 − β̄_44)·‖D⁴g‖ ≤ ….
All of Ā, B̄ carry a factor sin δ, k^{(i)}, η or e, each ≤ about 1e-40 on S [est]. Whole-shell boxes (no cells) are enough.

### 3.5 The time-integrated bounds, in the form the lemma uses
From §3.3–3.4, with all cross terms by AM–GM so that every bound is affine in T:

  ∫_I M_3² ≤ m3c + m3s·T,  ∫_I M_4² ≤ m4c + m4s·T,  ∫_I s2² ≤ σ2c + σ2s·T,

  σ2c := M̄_2·m3c/λ_s,  σ2s := M̄_2²/3 + M̄_2·(λ_s + m3s/λ_s)  (from 2 M_2 M_3 ≤ M̄_2(λ_s + M_3²/λ_s), any λ_s > 0).

**Bootstrap level.** The exports enter at the level of a new bootstrap hypothesis (BΦ) (§7.4): X_k(τ_a, T) ≤ 2(C0_k + C1_k·T), k = 1, 2. The pointwise energies enter at β_A·Ê = 2Ê as now.

## 4. Level 0 (pointwise in τ)

**Statement (L0).** For a.e. τ: K0^OR ≤ K̄0 := Σ_{l ≤ 2} κ⁽⁰⁾_l·N_l² + K0^nl + K0^w.

1. **Linear and quasilinear commutator, with modulation** [the F1 reduction with one more generic function].
   Integrand W ζ g·(j K4g + j′ K3g + k² K2g + C_eff ζ′g). Reduce with F1's algorithm (`F1/flux_reduce.py` `reduce_level(0, 2)`): a monomial c·g^{(a)}g^{(b)} with a > 2 is integrated by parts only if c vanishes identically at 57+ and 60− (exact one-sided values of ζ^{(i)}; W, C, k², j, j′ generic and nonzero). Only the monomial 4Wζζ′j·g g‴ qualifies. Result: pairs (0,0), (1,0), (2,0), (2,1), coefficients affine in (j, j′).
   - **No j″ may appear** (j″ contains g″, which has no pointwise bound). With target 2 the single integration by parts differentiates the coefficient once. The script asserts it.
   - Hand check of the same reduction in another form [hand]: with κ := 4Wζζ′j, ∫κ g g‴ = −∫κ′ g g″ + ½∫κ′ g′² (κ and κ′ vanish at both ends), giving pairs (2,0): −κ′ + 6Wζζ″j + 3Wζζ′j′; (1,1): ½κ′; (1,0): 4Wζζ‴j + 3Wζζ″j′ + 2k²Wζζ′; (0,0): Wζζ⁗j + Wζζ‴j′ + k²Wζζ″ + W·C_eff·ζζ′.
   - Sups: Γ_ab = sup|c_ab|/W (a ≠ b), Γ_aa = sup max(c_aa/W, 0), with **j, j′, C_eff as arb balls** (§2.2; C_eff ∈ C ± ā·(sup|Q_u| + sup|C_u|)). The one-sided diagonal uses the drift sign (C_eff ≥ 13.6 > 0, ζζ′ ≤ 0) once, here, and the modulation transport is inside the same ball. This removes the double use that ROUTE §4 pointed out.
   - κ⁽⁰⁾ by weighted AM–GM with fixed rationals as in F1.
   - **Regression:** j ≡ 1, j′ ≡ 0, ā = 0 must reproduce `F1/flux_forms.pkl` level 0 and κ⁽⁰⁾ = (8.4244, 4.1592, 2.2596) [cert F1].
2. **Lower-order nonlinear term** [hand]. (Wζv)′ = W·(ζ²g′ + (2ζζ′ − ψ′ζ²)g), so by Cauchy–Schwarz on OR

   K0^nl := (N_1 + KAP2·N_0)·[(d010 + d110·s1)·N_2 + (d100 + d200·s1 + d300·s1²)·N_1 + A1·N_0 + A0·w̄·ω_W],  KAP2 := sup_OR|2ζζ′ − ψ′ζ²|.
3. **Formal terms with w** [hand]. K0^w := N_0·w̄·ω_W·(sup_OR|c0′| + ā·sup_OR|(E_A′, E_B′)|₂). Size ≤ 1e-50 [est]; carried with certified constants, not dropped.

Sizes at the v35 data, β_A = 2 [est from cert]: N_0 = 2.10e-11, N_1 = 3.46e-12, N_2 = 2.11e-12; item 1 = 3.904e-21·(1 + 0.006) (quasilinear ≤ 2.2e-23, modulation 1e-6 relative); item 2 ≈ 4e-25. K̄0 ≈ 3.93e-21.

## 5. Level 2 (no integration by parts)

**Statement (L2).** For a.e. τ: |K2^OR| ≤ 𝔇(τ)·𝔛(τ), where 𝔇 ≥ ‖D2v‖_{W,OR} and 𝔛 ≥ ‖X2‖_{W,OR}:

  D2v = ζ·g⁗ + (4ζ′ − 2ψ′ζ)·g‴ + (6ζ″ − 6ψ′ζ′ + ψ₂ζ)·g″ + (4ζ‴ − 6ψ′ζ″ + 2ψ₂ζ′)·g′ + (ζ⁗ − 2ψ′ζ‴ + ψ₂ζ″)·g   [hand; sympy in X1]
  𝔇 := 𝔇_pt + 𝔇_int,  𝔇_pt := 𝔞_2 N_2 + 𝔞_1 N_1 + 𝔞_0 N_0,  𝔇_int := N_4 + 𝔞_3 N_3,
  𝔞_i := sup_OR|coefficient of g^{(i)} in D2v|: 𝔞_3 ≤ 4.55, 𝔞_2 ≤ 9.09, 𝔞_1 ≤ 10.48, 𝔞_0 ≤ 11.17 [cert, probe; plain ball evaluation, not sharp].

  𝔛 := 𝔛_pt + 𝔛_int, term by term (z_i := sup_OR|ζ^{(i)}|: z_1 = 35/48, z_4 = 280/27 exact; z_2 ≈ 0.835, z_3 ≈ 1.944 [est], to certify):

| Term of X2 | Bound of its W-norm on OR | Class |
|---|---|---|
| j·4ζ′g‴ | j̄·4z_1·N_3, j̄ := 1 + CBAR_OR | **int** |
| j·(6ζ″g″ + 4ζ‴g′) | j̄·(6z_2·N_2 + 4z_3·N_1) | pt |
| j·ζ⁗g + C_eff·ζ′g | sup_OR(j̄\|ζ⁗\| + \|C_eff ζ′\|)·N_0 (≤ 12.4·N_0 [cert, probe, at j̄ = 1]) | pt |
| j′·K3g | C1_OR·(3z_1·N_2 + 3z_2·N_1 + z_3·N_0) | pt |
| k²·K2g | sup k²·(2z_1·N_1 + z_2·N_0) (≤ 1e-86) | pt |
| ζN_low′: d010·g‴ | d010·N_3 | **int** |
| (Dd010 + d100)·g″ | (d010^a + d010^b·s1 + d100)·N_2 | pt |
| Dd100·g′ | (d100^a + d100^b·s1)·N_1 | pt |
| 2 d200·g′g″ | 2·d200·s1·N_2 | pt |
| Dd200·g′² | (d200^a + d200^b·s1)·s1·N_1 | pt |
| d110·g″² | d110·s2(τ)·N_2 | **int** (s2 is only L²/L⁴ in τ) |
| d110·g′g‴ | d110·s1·N_3 | **int** |
| Dd110·g′g″ | (d110^a + d110^b·s1)·s1·N_2 | pt |
| 3 d300·g′²g″ | 3·d300·s1²·N_2 | pt |
| Dd300·g′³ | (d300^a + d300^b·s1)·s1²·N_1 | pt |
| D(d000) = E1c·g′ + E0c | A1·N_1 + B1·N_0 + B0·w̄·ω_W | pt |
| ζc0′w + ζ(aE_A′ + bE_B′)w | (sup\|c0′\| + ā·sup\|(E_A′, E_B′)\|₂)·w̄·ω_W | pt |

  So 𝔛_int = x_3·N_3 + d110·N̄_2·s2, with x_3 := 4z_1·j̄ + d010 + d110·s1 and N̄_2 := sup_τ N_2.

**Split, pointwise in τ** (any fixed λ_1, λ_2 > 0; default 1):

  𝔇·𝔛 ≤ **K2^pt** + **φ2**,
  K2^pt := 𝔇_pt·𝔛_pt + (λ_1/2)·𝔇_pt² + (λ_2/2)·𝔛_pt²   (a constant: sup over τ at β_A·Ê),
  φ2 := (1/(2λ_1) + 1/2)·𝔛_int² + (1/(2λ_2) + 1/2)·𝔇_int²,
  𝔛_int² ≤ 2·x_3²·N_3² + 2·(d110·N̄_2)²·s2²,  𝔇_int² ≤ 2·N_4² + 2·𝔞_3²·N_3².

**Window bound of φ2** from §3.5 and N_l² ≤ W(57)·M_l²: ∫_I φ2 ≤ C″ + S″·T with

  C″ = a_X·[2x_3²·W57·m3c + 2(d110·N̄_2)²·σ2c] + a_D·[2·W57·m4c + 2𝔞_3²·W57·m3c],  a_X := 1/(2λ_1) + 1/2,  a_D := 1/(2λ_2) + 1/2,

and S″ the same with (m3s, m4s, σ2s).

Why this is enough: every monomial of D2v·X2 has at most two factors that are only L² in time (g‴, g⁗, s2), and s2 always comes with a pointwise N_2. No monomial needs sup|g‴|, a pointwise sup|g″|, or g^{(5)}.

Coefficient sups d010 … B1 on OR: stage X2 recomputes them on [57, 60] with the code path of `v4cert/c2_t3d.py` and asserts each ≤ the ZB value. ZB values at alpha-level 2.5e-4 [cert, `constants_v4.json`]: d010 1.536e-4 (a 1.51e-6, b 3.79e-3); d100 1.37e-6 (a 6.6e-9, b 4.2e-6); d200 3.41e-3 (a 3.8e-6, b 2.5e-3); d110 0.4115 (a 3.79e-3, b 10.153); d300 3.046 (a 8.4e-4, b 2.246); A0 1.5e-22, A1 2.45e-9, B0 3.0e-22, B1 1.14e-10.

## 6. The cutoff and the endpoint terms

### 6.1 The cutoff
ζ = S7((60 − x)/3), septic, C³; ζ^{(i)}(57) = ζ^{(i)}(60) = 0 for i = 1, 2, 3; ζ⁗(57+) = −280/27, ζ⁗(60−) = +280/27 [cert, `F1/flux_reduce.json` `zeta_end_values`]. The v3.1 tables used a quintic (U2_TRACE R3.3a). Every table of v2 is generated from the S7 polynomial on the closed interval; no object of `writeup/lemmaS/lemmaS_cflux*` is used except W57_upper as a regression value.

### 6.2 Endpoint terms: complete list
| Where | Term | How it is treated |
|---|---|---|
| Energy identity, level 2: ½ d/dτ V_2² = ∫ W D2v·v_τ | boundary terms at 10.25 and 60 of the two integrations by parts; interior points 16.25, 57 where ζ⁗ jumps | zero at the ends (v_τ, (v_τ)′ vanish there); no interior term because v ∈ C³ and (W v″)′ is Lipschitz (SPEC §3). Needs the regularity lemma U-1 (open, §11) |
| Level 0, global: ∫ W v ζ N_low′ = −∫ (Wζv)′ N_low | boundary at 10.25, 60 | zero (ζ = 0); then split by zones without further integration by parts (SPEC zone rule) |
| Level 0, OR, the monomial 4Wζζ′j·g g‴ | boundary at 57+ and 60− | **zero**: the coefficient vanishes identically at both ends (ζ′(57) = 0; ζ(60) = 0), machine-checked with exact one-sided values, as in F1 |
| Level 2, OR | — | **none: no integration by parts is done.** ζ⁗ enters only through sup\|ζ⁗\| = 280/27 on the open interval |
| F1's O4 = ∫Wζζ⁗(g g⁗ − 4g′g‴) + … and the point term W(57)·(280/27)·g(57)·g‴(57) | — | **does not arise.** The (4,0) monomial is bounded directly: ≤ (280/27)·N_0·N_4, inside 𝔇·𝔛 |

**Why not the trace route of ROUTE §1 item 3** [hand]. It bounds g‴(57) in L²(dτ) by Agmon from z_3, z_4 and uses no space-time bound of g⁗. For the old functional (commutator of L1 only) that suffices. For the new K2^OR the term ζN_low′ pairs with ζg⁗ and contains d110·(g″² + g′g‴) = d110·(g′g″)′. Removing g⁗ by parts gives, besides point terms at 57 that Agmon does control,

  −(5/2)∫ W ζ² d110·g″·g‴² dx + (terms in which the coefficient is differentiated),

and no rearrangement removes it: (g′g″)′·g⁗ + (5/2)·g″g‴² is a total derivative and no other multiple is [sym, `probe/id_probe.py`]; d110 = 10J⁻⁶U1 has no sign. With g‴ ∈ L²(dτ; L²_x) only, sup_x|g″|·‖g‴‖² ≲ M_2^{1/2}·M_3^{5/2} is not integrable in τ. The alternatives are a pointwise sup|g″| on OR (needs a pointwise level-3 energy of (A)), the verifier's V_3 (loses W(57)^{1/2}, changes the verifier's inequality to a time-dependent rate), or a space-time bound of g⁗. The last one costs only the conversion of §3.4, because the export is the same.
The same d110 structure is why the E_I spec needed S1z, S2z on the inner ramp (SPEC §4.6).

## 7. The time-integrated absorption lemma

### 7.1 Setting
Verifier: for a.e. τ, ½Λ′ ≤ tot_s(α, x) + μ·K0^OR + μ·eps_I·K2^OR, where tot_s is the right-hand side of `rhs_iv.rhs_pn` without the two flux constants (floors off), evaluated at the true state x(τ). The certificates of fin13v4b bound **tot := tot_s + F0**, F0 = μ(FLUX0 + eps_I·FLUX2).
Flux lemma (§4, §5): μK0^OR + μ eps_I K2^OR ≤ μK̄0 + μ eps_I K2^pt + φ, φ := μ·eps_I·φ2 ≥ 0, and ∫_I φ ≤ C′ + S′·|I| on every window, C′ := μ eps_I C″, S′ := μ eps_I S″.
Slack: **s := F0 − μK̄0 − μ eps_I K2^pt**. So

  ½Λ′ ≤ tot(α, x) − s + φ(τ)  a.e.   (★)

### 7.2 Lemma (reflected integral)
Assume (★) on [τ_1, τ′], φ ≥ 0 locally integrable, ∫_σ^τ φ ≤ C′ + S′(τ − σ) for all τ_1 ≤ σ ≤ τ ≤ τ′, and **S′ ≤ s**. Put P(τ) := 2∫_{τ_1}^τ (φ − s) and

  **Z(τ) := P(τ) − min_{σ ∈ [τ_1, τ]} P(σ) = 2·sup_{σ ≤ τ} ∫_σ^τ (φ − s).**

Then: (i) 0 ≤ Z ≤ 2C′; (ii) Z is absolutely continuous and Z′ ≥ 2(φ − s) a.e. (the running minimum of an absolutely continuous function is absolutely continuous and nonincreasing); (iii) **Y := Λ − Z** satisfies ½Y′ ≤ tot(α(τ), x(τ)) a.e., and Y ≤ Λ ≤ Y + 2C′; Z(τ_1) = 0.
Proof: one line each [hand]; written out in `flux2_lemma.tex`.

So Y obeys the verifier's inequality with the **constant** fluxes, but evaluated at a state whose energy is Λ ∈ [Y, Y + 2C′], not Y. The rest is to check that every certificate tolerates that shift.

### 7.3 Why the verifier's structure permits it (each item: what is used, where it is in the code, the condition)
| # | Certificate | Argument with Y [hand] | Condition |
|---|---|---|---|
| C0 | — | s > 0 after all pointwise terms are charged | s = F0 − μK̄0 − μ eps_I K2^pt > 0 |
| C1 | — | Z ≤ 2C′ needs the slope to be absorbed | S′ ≤ s |
| C2 | March cell k: `orc.rate(aq, Ly_lo, Ly_hi)` bounds tot/Λ by U_k for all shapes, α ≤ aq, Λ ∈ [Ly_lo, Ly_hi] (`closure.py:119-121`, `bb.bound_cell`; TAIL_LEMMA §"leaf certificate") | Claim Y ≤ L(τ). Near a first touch Y = L: Λ ≥ Y > Ly_lo and Λ ≤ L + 2C′; the state is in the cell if **L(τ) + 2C′ ≤ Ly_hi_k** on the cell. Ly_hi = max(L_k, g)·(1 + 1e-3) (`closure.py:112`) | 2C′ ≤ 1e-3·Ly1 = 2.5e-18 |
| C3 | Strictness pad: barrier L′ = 2(U_k + 1e-9·max(1, \|U_k\|))·L (`closure.py:120`) | G := Y − L. If U_k ≥ 0: Y′ ≤ 2U_kΛ ≤ 2U_k(Y + 2C′), so G′ ≤ 2U_kG + (4U_kC′ − 2·pad_k·L) and the bracket must be < 0. If U_k < 0: Y′ ≤ 2U_kΛ ≤ 2U_kY (Λ ≥ Y), no condition. Gronwall on the left neighbourhood handles "a.e." | **U_k⁺·2C′ < 1e-9·max(1, \|U_k\|)·L_k for every k.** Binding at k = 0 of the 1.5e-4 table: U_0 = 0.0304299, L_0 ≥ 2.5e-15 [cert, rate table], so **2C′ < 8.2155e-23**; 1.25e-4 table: U_0 = 0.0210584, 2C′ < 1.187e-22 [cert arithmetic] |
| C4 | Tail anchor: tot_s/Λ ≤ −ρ on Λ ∈ [L*, (1 + η)L*], α ≤ a*, η = ETA = 0.01 (`closure.py:163`); joint scaling (TAIL_LEMMA Lemma joint) | Claim Y ≤ ℓ(τ) = max(L_F, L*e^{−ρ(τ − τ_t)}). At a touch Y = ℓ: Λ ∈ [ℓ, ℓ + 2C′] ⊂ [ℓ, (1 + η)ℓ]; scale by (ℓ/L*)^{1/2}; then ½Y′ ≤ tot_s + F0 ≤ −ρΛ + F0 ≤ −ρℓ + F0 (the coefficient of Λ is negative and Λ ≥ ℓ). This is exactly the inequality of the floor paragraph of TAIL_LEMMA.tex; its proof goes through unchanged, with the same L_F = 2F0/ρ | 2C′ ≤ η·L_F = 9.80e-21 [cert: sqrt_LF 9.900007e-10] |
| C5 | Tail α-condition, whole tail (`closure.py:185-190`, margin 177.08) | β is nondecreasing in Λ, and Λ ≤ (1 + 2C′/L_F)·ℓ. The α-barrier κ√ℓ stays as it is | lhs ≥ rhs_whole_tail·(1 + 2C′/L_F); 177.19 against 0.114 [cert] |
| C6 | Unconditional integrals (`closure.unconditional`, TAIL_LEMMA §V1): θ∫NEGV ≤ ½Λ(start) + Σ(r1⁺·Ly_hi + Φ)D + … | Apply to Y: 2θ∫_a^b NEGV ≤ Y(a) − Y(b) + 2∫(r1⁺Λ + Φ + F0), with Λ ≤ Ly_hi (C2), Y(a) ≤ Λ(a), −Y(b) ≤ 2C′. **X_const gains C′; the slope is unchanged.** The same for windows from any τ_a (`hjw_rep`) | J_m gains J_m·C′/X_const, X_const = 3.7483e-15 [cert]. Needed: **κ_J(used by (A)) − κ_J(certified) ≥ C′/X_const** (≈ 3e-11 [est]; up to 3e-10 pessimistic) |
| C7 | Every other consumer of "Λ ≤ barrier" | must use a cushioned value: `farsup`, `sramp_sup`, `box_assert`, `sqrtL_max`, `mod_sups` use Ly_hi [read ROUTE §5 for fin12f4; **to re-read on fin13v4b**]; quoted tail floors (√L_F, α_F) get the factor (1 + 2C′/L_F)^{1/2} ≤ 1 + 1.2e-7 [est]; the v36 chain's consumers of the rate table: `v36cfg.verifier_inputs` rounds J, slopes, peak up to 5 digits and ρ down [read `N36/v36cfg.py`]; `assembly_arb.f4b` reads the raw sqrt_LF into PHI_S [read], used only for v3.1-era quantities (I3, (B)-exit) | list every consumer; each has reserve ≥ 2C′/L (march, ≈ 1e-10) or ≥ C′/L_F (floor, ≈ 1.2e-7), or is unused |
| C8 | Regularity | Λ absolutely continuous; φ ∈ L¹_loc; (★) a.e. | lemma U-1 (open) |
| C9 | Bootstrap order | §7.4 | C1–C6 with the (BΦ)-level constants |

Why this works at all, in one paragraph: (1) both flux inputs enter the verifier additively, as one constant in the degree-0 group, so only μK0 + μ eps_I K2 matters; (2) every certificate is a statement about all states in a cell [Ly_lo, Ly_hi] whose upper end has a cushion of 1e-3 relative, 1e7 times larger than 2C′/Λ; (3) the barrier is strict by a pad of 1e-9; (4) in the tail the inequality is monotone in Λ in the right direction; (5) the J's are computed from the cushioned Ly_hi and move by one additive constant.

Difference from TAIL_LEMMA.tex Lemma "L¹_loc forcing": that lemma uses a discounted integral, replaces the floor and covers the tail only. Here the pointwise slack s pays the slope, the floor L_F is unchanged, and the march is covered. TAIL_LEMMA.tex is not edited; the new lemma is a separate file.

### 7.4 Bootstrap order (C9) [hand]
The window bounds of M_3, M_4 are conclusions of (A) (Lemma D, k = 1, 2). (A) uses the verifier's J's (H-J-W), which now depend on C′. To break the loop, add to the bootstrap set:

  **(BΦ)**: for k = 1, 2 and every window [τ_a, τ_a + T] ⊂ [τ_1, τ′]: X_k(τ_a, T) := ∫_{τ_a}^{τ_a+T}∫_{ln 28}^{ln 350} εω (z^A_{k+2})² ≤ 2·(C0_k + C1_k·T).

- Under BS(τ′) with (BΦ): the flux lemma holds with (C′, S′) computed from the doubled exports and from 2Ê; the absorption lemma gives the verifier's conclusions with Λ ≤ L + 2C′; (A) then gives X_k ≤ C0_k + C1_k·T, strictly better.
- (BΦ) is closed and open in τ′ in the sense of the closure proposition (X_k is continuous in the window end; strict improvement).
- (A)'s inputs do not move: the J's it reads are the verifier's rounded values, valid by C6/C7.
- Lemma D at k = 1 already exists in v35. **v36 does not need (BΦ) as a hypothesis of (A)**; it is a hypothesis of the verifier side only.

### 7.5 What changes in the verifier's statement
- Hypothesis (I3) of the E_I lemma becomes (I3′) of §1 with (s, φ, C′, S′).
- Conclusion B1 becomes Λ(τ) ≤ L(τ) + 2C′ (march) and Λ ≤ ℓ(τ) + 2C′ (tail); α-barriers unchanged.
- No number of the rate tables changes. No rerun.

## 8. Constants to certify

All in python-flint arb, 200 bits, outward rounding to doubles by `nextafter`, PYTHONHASHSEED=0, `python3 -B`. "OR cells" = 3000 cells of [57, 60] with ζ^{(i)} in mean-value form on the clipped cell of u = (60 − x)/3 (F1's method; plain ball evaluation overshoots by 2–9% [cert, probe]). Every sup is rounded **up**; the two lower bounds (inf ψ′, s) are rounded **down**.

| # | Constant | Definition | Method | Inputs | Stage |
|---|---|---|---|---|---|
| 1 | z_1 … z_4 | sup_OR\|ζ^{(i)}\| | OR cells; exact values 35/48 and 280/27 asserted as lower checks | — | X2 |
| 2 | 𝔞_0 … 𝔞_3, KAP2, sup ψ′, sup\|ψ₂\|, inf ψ′ | §5, §4 | OR cells; ψ^{(j)} exact in x and (1 + x²)^{1/3} | — | X2 |
| 3 | W(57), ω_W² | W = e^{−ψ} at 57; ∫_OR W ≤ W(57)/inf ψ′ | W(57) from the verifier's definition of ψ (integral of ψ′ from the switch, as `agent_sob`); regression against `W57_upper` | `V/constants.py`, `agent_sob` | X2 |
| 4 | profile balls on S: k, m, p, q − σ, C, sup\|c0′\|, Q_u, C_u, E_A′, E_B′ | `hstail.profile_bounds(inf S)` and the profile ODE | arb | `step3/hstail.py` | X2 |
| 5 | CBAR_OR, c_a, c_b; d010, d100, d200, d110, d300 with ^a, ^b; A0, A1, B0, B1 on OR | SPEC §4.3–4.5 | the code path of `v4cert/c2_t3d.py` restricted to x ∈ [57, 60], boxes B-a, B-c, B-d; assert ≤ the ZB values of `constants_v4.json` | `step3nl_v4b/v4cert/` (imported read-only), c4fun | X2 |
| 6 | Γ⁽⁰⁾_ab, κ⁽⁰⁾_l | §4 item 1 | reduced forms of X1 compiled to arb; OR cells; j, j′, C_eff balls; one-sided diagonal; fixed rational AM–GM weights | X1 forms, #4, #5, ā | X2 |
| 7 | s0, s1 | S0V; S1V/inf S | arb from the env of the rate table | rate table `args.env` | X2 |
| 8 | Ā_4, B̄_42, B̄_43, B̄_44, B̄_422, β̄_44; sup and ∫_S of D⁴d, P″ … P⁗ | §3.4 | sympy structure (asserted) → arb on S (cells 1/4 or one box); c4fun ext_jets | `n1_chart.py` (model), c4fun, hstail | X3 |
| 9 | N̄_0, N̄_1, N̄_2 (and M̄_l) at β_A = 2 | §3.2 | conversions of `f4b` with the enclosure of S; d, chart per level by Minkowski; **regression:** Σκ_old·N_l² reproduces the pipeline's K0, K2 at β = 2 of the v36 assembly | v36 `A_interface.json`, `asm_design_N2ref4.json`, `b0_dprof.json`, v36 n1_chart | X4 |
| 10 | m3c, m3s, m4c, m4s, σ2c, σ2s | §3.3–3.5, at the (BΦ) level | arb from the exports k = 1, 2 and #8, #9 | v36 export blocks | X4 |
| 11 | K̄0, K2^pt, C″, S″, C′, S′, s | §4, §5, §7 | arb | #1–#10; μ, eps_I, FLUX0, FLUX2 from the rate table args | X5 |
| 12 | C0–C7 margins, both α points | §7.3 | arb on the stored rate tables (every step k), tail block, `unc`, v36 `hjw_rep.json` | `runs/P/rate_table_fin13v4b_a*.json`, v36 hjw | X5 |

Free parameters (any positive rational is valid; fixed in the script, not searched in the certified run): AM–GM weights of level 0; λ_1, λ_2, λ_s.

## 9. Feasibility

Inputs: [cert] v35 chain (`B/fix2/n12/t2_A_v35/out/asm_design_N2ref4.json`: Ê = (1.2130e-6, 1.0451e-4, 5.8246e-2); level-1 export C0_1 = 1.6631e-4, C1_1 = 7.684e-8); [cert] F1 (K0 ≤ 3.9038e-21 at β_A = 2); [cert] fin13v4b rate tables; [cert] `constants_v4.json`. The v36 values will differ slightly (verifier J's +0.08%, R3 forcing, r_2): **all results below are [est]** and must be recomputed on v36. Script and output: `probe/feas_probe.py`, `probe/feas_probe.out`.

### 9.1 Pointwise condition
| Quantity | Value | Label |
|---|---|---|
| K0 linear commutator (F1) | 3.9038e-21 | [cert, v35] |
| + quasilinear (j, j′), modulation, nonlinear, w-terms | +0.6% → K̄0 ≈ 3.93e-21 | [est] |
| 𝔇_pt, 𝔛_pt | 2.90e-10, 2.99e-10 | [est from cert] |
| K2^pt (λ_1 = λ_2 = 1) | 1.73e-19; eps_I·K2^pt = 5.2e-24 | [est] |
| **(FLUX0 + eps_I FLUX2)/(K̄0 + eps_I K2^pt)** | **≈ ×4.17** | [est] |
| slack s | ≈ 9.97e-22 | [est] |

### 9.2 Integrated condition
The level-2 export is not available yet. Guess: C0_2 = (Ê_2 + D_2)/(r_2·a_c) ≈ 2·0.058/(1e-3·0.976) ≈ 119, C1_2/C0_2 as at level 1 [est]. All at the (BΦ) level (exports doubled):
| Quantity | Guess | ×10 worse C0_2, C1_2 | Limit | Label |
|---|---|---|---|---|
| ∫N_3² | 6.8e-23 + 5.1e-26·T | same | | [est from cert] |
| ∫N_4² | 2.2e-20 + 1.0e-23·T | 2.2e-19 + 1.0e-22·T | | [est] |
| C′ | 1.14e-25 | 1.05e-24 | | [est] |
| **2C′ against C3** | 2.3e-25, **×361** | 2.1e-24, ×39 | 8.2155e-23 [cert] | [est] |
| S′ against s | 5.5e-29 | 4.9e-28 | 9.97e-22 | [est] |
| C6: C′/X_const | 3.0e-11 | 2.8e-10 | reserve of κ_J in v36 (unknown) | [est] |

**Kill value for (A):** C3 fails when C0_2 ≳ 4e4 [est], about ×350 above the guess. If C0_2 comes out between 4e3 and 4e4, raise λ_2 (moves cost from C′ to K2^pt, which has room of about 2000×) before anything else.

### 9.3 The separate comparison K2 ≤ FLUX2
Not needed (§1). For the record: the crude K2^pt = 1.7e-19 exceeds FLUX2 = 8.0e-20. With F1's reduced forms for the linear part one gets P2 ≈ 5.93e-20 (F1's 5.886e-20 [cert] plus ≈ 4e-22 [est] for the new pointwise terms), ×1.35. If a referee prefers separate statements, either add that optional stage or state the combined condition; **I recommend the combined condition** and no optional stage.

## 10. Coding plan (Opus; mechanical parts may go to Sonnet)

Directory: `B/fix7_flux2/code/` (scripts), `B/fix7_flux2/out/` (outputs), log `FLUX2_LOG.md` (resumable, one entry per stage), report `FLUX2_REPORT.md`.
Rules: proof tree read-only (imports with `sys.path`, `python3 -B`, PYTHONDONTWRITEBYTECODE=1; check no `__pycache__` appears); one process, nice 19, PYTHONHASHSEED=0, < 800 MB; every stage writes an `args` block (argv, env, every input path with sha256); atomic writes; two consecutive runs byte-identical apart from timing; `SHA256SUMS`; before/after sha snapshot of `step3nl_v4b`, `fix3`, `fix5`, `fix6_v36`.

| Stage | Script → output | Content | Checks | Agent-h | Compute, RAM |
|---|---|---|---|---|---|
| X0 | `inputs_manifest.py` → `inputs.json` | Pin every input by sha256: rate tables (2), `constants_v4.json`, c2_t3d sources, hstail, c4fun, F1 files, **v36** `A_interface.json`, `asm_design_N2ref4.json`, `hjw_rep.json`, n1_chart (fin13v4b), `b0_dprof.json`. Refuse to run if `N36/V36_REPORT.md` is missing | env of the rate tables = (FLUX0, FLUX2, S0V, S1V) expected | 0.5 | s |
| X1 | `flux2_forms.py` → `flux2_forms.pkl`, `.json` | sympy, exact. (a) the identity of §2.1 (commutators, j-form, C_eff) from the v-equation; (b) level-0 reduction with generic j (F1's algorithm), Euler-operator check, assert no j″; (c) D2v coefficients; (d) N_low′ in g-jets from `graph_struct.pkl`, equal to the SPEC §4.5 list | regression: j ≡ 1 reproduces `F1/flux_forms.pkl` level 0 exactly | 2.5 | < 1 min, < 300 MB |
| X2 | `flux2_coef.py` → `flux2_coef.json` | constants #1–#7 | regression: j ≡ 1, ā = 0 gives F1's Γ⁽⁰⁾, κ⁽⁰⁾; OR sups ≤ ZB values; W(57) against `W57_upper`; z_1, z_4 against exact values | 3 | 2–4 min, < 300 MB |
| X3 | `chart4_shell.py` → `chart4_shell.json` | constants #8: χ_4 structure in sympy (asserted), arb bounds on S, D⁴d, P-jets | k ≤ 3 on S reproduces values ≤ the `n1_chart` core values; χ_4 = 0 for a ray; float spot check of the level-4 identity on a synthetic curve (not part of the proof) | 4.5 | 2–5 min (sympy), < 500 MB |
| X4 | `flux2_data.py` → `flux2_data.json` | constants #9, #10 | regression against the v36 pipeline K0, K2 at β = 2 (old κ); m3 against the v36 "honest I3" within the stated rule differences | 2.5 | s |
| X5 | `flux2_assemble.py` → `flux2_cert.json` | constants #11, #12; verdict PASS/FAIL per condition, both α points | all margins as certified lower bounds; fails loudly | 2.5 | s |
| X6 | `flux2_sanity.py` → `flux2_sanity.json` | float, not part of the proof: 30-digit quadrature of K0^OR, K2^OR for 3 test triples (g, w, coefficients) against the bounds and against the level-0 reduced form; a simulated φ with the reflected Z (Z ≤ 2C′, Z′ ≥ 2(φ − s)) | — | 1.5 | s |
| X7 | texts: `flux2_lemma.tex` (numbers filled from `flux2_cert.json` by a generator, as `n35_gen_numbers.py`), `absorption_lemma.tex`, erratum entries (cutoff shape, jet box, drift double use, trace route), SPEC §4.10 pointer, (BΦ) paragraph for the continuation text (new file, v36 suffix) | — | balance check | 2.5 | — |
| X8 | double run, `SHA256SUMS`, tree snapshots, `FLUX2_REPORT.md` | — | byte-identity apart from timing; nothing outside `fix7_flux2` changed | 1 | 2 × total |

**Total: about 20–22 agent-hours**, compute under 15 min per full pass, peak RAM under 500 MB [est].
Then: **Fable review** (one pass, from scripts and texts, not from the report): the identity of §2.1 against the verifier code; the level-4 chart structure; C2–C7 against `closure.py`/`bb.py`; the bootstrap order.
Order: X1–X3 do not need v36 and can start now. X0, X4, X5 need `N36/V36_REPORT.md` and the final v36 outputs.
Stop rules: stop and report if (i) X1 finds a j″ or a non-vanishing boundary coefficient at level 0; (ii) the sympy structure of χ_4 is not as in §3.4; (iii) any of C0–C6 fails; (iv) the combined margin is below ×2.

## 11. What is not yet justified

1. **§2.1 identity** [hand]. Assembled from SPEC's [sym] pieces; the sign bookkeeping (j = 1 − c, C_eff) is mine. Checked here in sympy [sym, `probe/id_probe.py`]: the j-commutator identity ζ·(−(j g‴)′) = −(j v‴)′ + j K4g + j′ K3g, the D2v coefficients of §5, the level-0 hand reduction of §4, and the two conversions of §3.3. Not checked here: the C_eff sign (it follows SPEC T9) and the N_low′ list. X1 redoes all of it from the v-equation.
2. **Level-4 chart structure** (§3.4) [hand, weight count]. If sympy finds another product of high jets, the bound needs one more Hölder step; I do not expect it.
3. **Level-2 export of (A)** in window form: asserted in `N36/t2_A/rigor/export_t2.py:26-36` [read], not yet reported or reviewed. ROUTE's unit G1 (written statement of Lemma D for k = 2) is not in the tree.
4. **The verifier-side reading in C2–C5**: "the cell certificate bounds tot/Λ on all of [Ly_lo, Ly_hi]" and the first-touch proofs. I read `closure.py` and TAIL_LEMMA.tex; I did not re-derive `bb.py`. The Fable review must.
5. **C6/C7 reserves in v36** are unknown until its report exists. If κ_J or a floor consumer has no reserve, the v36 chain needs a rerun with a padded value (minutes of compute, new hashes).
6. **Regularity (U-1)**: Λ absolutely continuous, g⁗(τ, ·) ∈ L²(OR) for a.e. τ, the level-2 identity with a C³ cutoff. Open since SPEC §3; the flux lemma adds nothing new to it, but depends on it.
7. **Bootstrap level of ā** (§2.3) and the wording of (BΦ) inside the closure proposition: text work, to be fixed in X7 against the v36 continuation text.
8. **Old code, wrong-side constant (not used by v2):** `assembly_arb.py:175-178` uses 57.75 as the upper end of y on the shell, while sup Y(S) ≈ 57.7510. For N0 this is covered by the 60/28 factor [hand]; for `fac` (the I3 pipeline) I did not check. v2 replaces these formulas.
9. Every number of §9 is an estimate on v35 data.

## 12. What v2 needs from the (A) chain v36

`N36/V36_REPORT.md` did not exist when this was finalised (checked last at 07:31 on Nov 20: the v35 regression had passed, the R3 regression was running, no production run yet). So the design is against ROUTE's specification and the patched v36 scripts, not against v36 numbers. From `N36/V36_LOG.md` and the patched scripts [read], v36 is producing the items marked "yes".

| # | Quantity | Form | In v36? |
|---|---|---|---|
| A1 | Ê_j, j = 0, 1, 2, at the supersolution, and β_A = 2 | `A_interface.json` rule `EA_of_N2`; `asm_design_N2ref4.json` `E` | yes (as v35) |
| A2 | **Level-2 export**: for every window, ∫∫_{[ln 28, ln 350]} εω (z^A_4)² ≤ C0_2 + C1_2·T, with r_2 = 1/1000, a_c | `export_t2` block `r2_export`; `assembly_rep` `export['2']` | yes (log item 6) |
| A3 | Level-1 export C0_1, C1_1 | `export['1']` | yes |
| A4 | Flux block at β = 2 with d and chart: K0_w, K2_w, the d-part, RC, BC (for the regression of N_l) | `asm_design_N2ref4.json` `flux` | yes (as v35) |
| A5 | n1_chart rerun on the fin13v4b sups (core block: χ^add_k, β, ρ_meas, e, k ≤ 3) | `N36/n1/out/n1_chart_fin13v4b.json` | yes (log item 2) |
| A6 | `hjw_rep.json`: κ_J certified and the value (A) uses; X_const | as v35 | yes; **reserve ≥ 3e-10 to be checked (C6)** |
| A7 | The list of rate-table quantities the chain reads raw (unrounded) | report | **to ask the v36 agent now**, together with A6: if κ_J is consumed with less than 3e-10 of reserve, pad it before the production run rather than after |
| A8 | Margins S^v_0, S^v_1 against S0V, S1V (the boxes B-e this lemma uses) | `margins` | yes |

Not needed from (A): a pointwise level-3 quantity; a trace; fifth derivatives; a change of FLUX0/FLUX2 (v36 asserts the fin13v4b values 1.64e-20, 8.0e-20 [read `v36cfg.py`]).

## 13. Downstream

**If the lemma passes with the present FLUX0, FLUX2 (expected):**
- No verifier rerun; fin13v4b tag, rate tables, J's unchanged. v36 chain unchanged, unless C6/C7 lack reserve (§11 item 5).
- New pinned objects: `fix7_flux2/out/*` (to add to the step-3 kit: `kit_step3` v4 gains X1–X5; about 10 min more run time [est]).
- Texts: the E_I lemma hypothesis (I3) → (I3′); B1 → Λ ≤ L + 2C′; tail floor quotes with the factor (1 + 2C′/L_F)^{1/2}; (BΦ) in the bootstrap set and the closure proposition; Lemma S(A) (e) replaced by this lemma, Lemma D stated for k ≤ 2; erratum entries; numbers file; index.
- Step 4: nothing. The hand-off thresholds do not see the flux.

**If FLUX0 or FLUX2 had to change** (not expected; needed only if the combined margin fell below 1):
1. Verifier: new tag (fin13v4c), `run_v4b.sh` with the new env: tail_cert, both closures, mod_sups (about 15 min compute). F0, L_F, α_F, the J slopes and sJ4 move in proportion to F0; the peak and ρ essentially do not (V33_VERIFIER.md: the ×4 change of fin12f4 moved only slopes, sJ4, F0).
2. (A) chain: new version (v37) with the new verifier inputs and the new thresholds typed in `assembly_rep.py`/`v36cfg.py`: t2 → t4 → zone3; then X0, X4, X5 again. The loop FLUX → J slopes → Ê → K → FLUX is contracting with a tiny gain (slopes enter Ê through window slopes only) [est]; one pass should do.
3. Step 4 (v33-successor manifest): every stage that reads verifier sups or the (A)/U/zone3 interfaces: `n1_chart`, `f12_sliver`, `t9_cone_v33`, `h4loc_v33`, `t9_lemmaO_v33`, `ratio_v33`, `notes_v33`, and the manifest itself. **The exact list is the re-pin list of `N36/V36_REPORT.md`; use that, not this line.** The closed 34-stage manifest is not touched.
4. Kits: kit_step3 (verifier outputs), kit_step4b (v33-successor stages); texts: (W1) floor values, numbers files, main_v3.
All of that needs Micah's go first.

## Files in this directory
- `FLUX2_DESIGN.md`: this file.
- `flux2_lemma.tex`: LaTeX-ready statements (flux lemma v2, data lemma, absorption lemma) with proof sketches. Not compiled.
- `probe/feas_probe.py`, `probe/feas_probe.out`: feasibility probe (arb for the cutoff sups; the rest [est]). Not part of any proof.
- `probe/id_probe.py`, `probe/id_probe.out`: exact sympy checks of five hand identities (all True).
