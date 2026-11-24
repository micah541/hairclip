# SCOPE: the L² blow-up criterion for curve diffusion flow in Lean (replacing DKS Thm 3.1)

Date: 2026-11-23. Scoping only: nothing was built, no file outside this directory was written.
Labels: **[read]** = checked by reading the cited source file today; **[hand]** = hand derivation in this note, not machine-checked; **[GUESS]** = Lean/Mathlib name or estimate not verified by compilation.

## 0. Verdict

**GO.** The criterion is a short, self-contained piece of 1D analysis on top of what the sibling repo already proves. Estimated 1,900–2,900 new Lean lines in about 7 new files, 30–50 agent-hours, no existing file touched.

- The maximal flow, uniqueness, gluing and "bounded |A| ⇒ extends" are already Lean theorems [read]. Only "∫κ² ds bounded ⇒ sup|κ| bounded on [t₀, T)" is missing.
- That missing step needs one energy, E₁ = ∫κ_s² ds, three elementary periodic inequalities, and no Grönwall, no length monotonicity, no Fenchel theorem and none of the repo's n-dimensional interpolation machinery (§2).
- Main risks: the pointwise evolution equations (κ_t, commutator) in Lean, and the chart bookkeeping between `Circle`'s stereographic charts and the parametrisation θ ↦ F(e^{iθ}). Neither is a mathematical risk.
- Side benefit: the same work gives a Lean replacement for DKS (D1) as used in C-1 ("bounded ∫κ² ⇒ the flow extends past T"), and `exists_maximal_flow` + `rst_glue` replace the [hand] Lemmas "maximal smooth solution" and "gluing" of `U1_regularity_v41.tex`.

## 1. Target statements

### 1.1 What exists [read]

Sibling repo `/mnt/c/Users/micahw/agents/lean/CWFlow` (HEAD 8718f95, with about 20 uncommitted modified files in `AK/`, `Bridge/`, `Extension/`: the user is editing it).

- `CWFlow.exists_maximal_flow` (`Flow/MaximalFlow.lean:45`). Hypotheses: `ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E`, `hJi`, `hJJ`, `hdim : finrank ℝ V = 2 * n`, `hF₀ : IsLagImmersion E J F₀`. Conclusion: `Tmax : ℝ≥0∞`, `0 < Tmax`, `Φ` with `Φ 0 = F₀`; `IsFlow E J T Φ` for every `0 < T`, `ofReal T < Tmax`; every flow from `F₀` on `[0,T]` has `ofReal T < Tmax` and equals `Φ` there; and
  `Tmax ≠ ⊤ → ∀ K, ¬ ∀ t ∈ Ico 0 Tmax.toReal, ∀ z, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ K ^ 2`.
- `CWFlow.flow_extends_of_bounded_A` (`Flow/ExtensionTheorem.lean:35`): flow on every `[0,T']`, `T' < T`, with `aNormSq ≤ K²` on `Ico 0 T` ⇒ an `IsFlow` on `[0, T+ε]` agreeing on `Ico 0 T`.
- `rst_glue` (`Flow/Restart.lean:122`): glue `Φ` on `[0, a+s]` and `Θ` on `[0,T₁]` when `Θ t = Φ (a+t)` on an overlap `[0,s]`, `s > 0`. This avoids the seam problem that left `glue_seam` as a `sorry` in N1.lean.
- `aNormSq J X x = gNormSq (ginvZ X x) (cX J X x)` (`Extension/ImmCoord.lean:83`), with `gNormSq h T = Σ_I Σ_J (Π_i h![I i, J i]) T I T J` (`Extension/GNorm.lean:193`), `ginvZ` = entries of `(gX X x)⁻¹`, `cX` = `⟪∂_a∂_b X, J ∂_c X⟫`. Reparametrisation invariance: `aNormSq_comp` (`Extension/ImmReparam.lean:428`); locality: `aNormSq_congr` (`Extension/ImmSmooth.lean:132`).
- n = 1 instance: `/mnt/c/hsl/CWFlow/CWFlow/N1.lean` (= `lean_n1/CWFlowN1_v5.lean`, 814 lines), compiled only in the **older** copy of the repo (commit d72b744, which has no `Extension/`, no `MaximalFlow`). It has `E1`, `Jc`, `Jc_isom`, `Jc_JJ`, `hdim1`, `lag_dim_one`, `isLagImmersion_iff`, `curv`, `dS`, `flowVel_eq_neg_kss` (for `X : ℝ → ℂ`), `IsFlow.shift_gen`, `IsFlow.conj_reflect`, `cdf_symmetric`, and one `sorry` (`glue_seam`).

### 1.2 Instantiation at n = 1: no obstruction found [read, not compiled]

| Hypothesis of `exists_maximal_flow` | n = 1 instance |
|---|---|
| `E : Type u`, `V : Type u` (same universe) | `E1 = EuclideanSpace ℝ (Fin 1)`, `ℂ`, both in `Type` |
| `[NeZero (finrank ℝ E)]` | `instNeZeroFinrankE1` (N1.lean:57) |
| `[ChartedSpace E M] [IsManifold 𝓘(ℝ,E) ∞ M] [CompactSpace M] [T2Space M] [SigmaCompactSpace M]` | `Circle`; Mathlib `Geometry/Manifold/Instances/Sphere.lean:571–574` gives the charted space and `IsManifold (𝓡 1) ω`. N1.lean already instantiates `exists_flow_of_immersion`, which has the same variable block |
| `ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E` | `n = 1`, `ι = ContinuousLinearEquiv.refl ℝ E1` |
| `hJi`, `hJJ` | `Jc_isom`, `Jc_JJ` |
| `hdim : finrank ℝ ℂ = 2 * 1` | `Complex.finrank_real_complex` |
| `IsLagImmersion E1 Jc F₀` | `isLagImmersion_iff` (Lagrangian is automatic: `lag_dim_one`) |

Caveat: N1.lean must first be ported to the sibling repo as a **new** file. Its imports (`Flow/Main`, `UniqueMain`, `NormalFlow`, `GaugeEvolution`, `Scalar/Reparam`, `Scalar/Smooth`, `Flow/DiffeoChart`, `Flow/TanAlgebra`) all exist there, but the sibling repo is 5 commits ahead plus uncommitted edits, so small breakages are possible [GUESS: minor]. The `glue` section (the `sorry`) is dropped in the port; the repo forbids `sorry`, and `rst_glue` supersedes it.

### 1.3 Definitions (new)

```lean
namespace CWFlow.Curve
/-- The 2π-periodic lift of a map on the circle. -/
def lift (F : Circle → ℂ) (θ : ℝ) : ℂ := F (Circle.exp θ)
/-- κ of a parametrised plane curve (as `CWFlowN1.curv`): ⟪c'', i c'⟫ / |c'|³. -/
def curv (c : ℝ → ℂ) (θ : ℝ) : ℝ := ⟪deriv (deriv c) θ, Complex.I * deriv c θ⟫_ℝ / ‖deriv c θ‖ ^ 3
/-- ∫ κ² ds over one period, ds = |c'| dθ. -/
def bendingEnergy (c : ℝ → ℂ) : ℝ := ∫ θ in (0:ℝ)..(2 * π), curv c θ ^ 2 * ‖deriv c θ‖
```

`bendingEnergy (lift (Φ t))` is the paper's ∫_{γ(t₁+t)} κ² ds: an integral over the immersed parametrised curve, one traversal. It is an explicit interval integral, with no measure on `Circle` and no partition of unity.

### 1.4 Main targets

```lean
/-- (T1) Blow-up criterion; replaces DKS Thm 3.1 as used in Lemma "blow-up criterion, the form needed". -/
theorem cdf_exists_maximal {F₀ : Circle → ℂ} (hF₀ : IsLagImmersion E1 Jc F₀) :
    ∃ Tmax : ℝ≥0∞, 0 < Tmax ∧ ∃ Φ : ℝ → Circle → ℂ, Φ 0 = F₀ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E1 Jc T Φ) ∧
      (∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = F₀ →
        ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0:ℝ) T, Ψ t = Φ t) ∧
      (Tmax ≠ ⊤ → ∀ K : ℝ, ∃ t ∈ Ico (0:ℝ) Tmax.toReal, K < bendingEnergy (lift (Φ t)))

/-- (T2) Extension form; replaces DKS (D1) as used in C-1. -/
theorem cdf_extends_of_bounded_energy {T : ℝ} (hT : 0 < T) {Φ : ℝ → Circle → ℂ}
    (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E1 Jc T' Φ) {K t₀ : ℝ} (ht₀ : t₀ ∈ Ioo 0 T)
    (hK : ∀ t ∈ Ico t₀ T, bendingEnergy (lift (Φ t)) ≤ K) :
    ∃ ε > 0, ∃ Ψ : ℝ → Circle → ℂ, IsFlow E1 Jc (T + ε) Ψ ∧ ∀ t ∈ Ico 0 T, Ψ t = Φ t
```

(T1) follows from (T2) and clauses 1–3 of `exists_maximal_flow`: the extended `Ψ` has `Ψ 0 = F₀`, so `ofReal (T+ε) < Tmax = ofReal T`, a contradiction. The fourth clause of `exists_maximal_flow` is not needed. (T2) only assumes the energy bound near `T`, which is what the continuation argument has anyway.

Supporting statements that will be visible in the paper:

```lean
/-- (T3) |A|² in an ι-chart is κ² of the lift. -/
theorem aNormSq_eq_curv_sq {F : Circle → ℂ} (hF : IsLagImmersion E1 Jc F) (z : Circle)
    {y : E1} (hy : y ∈ Extension.chartDomE ι1 z) :
    ∃ θ : ℝ, Circle.exp θ = Extension.chartInv ι1 z y ∧
      Extension.aNormSq Jc (F ∘ Extension.chartInv ι1 z) y = curv (lift F) θ ^ 2
/-- (T4) the lifted flow is CDF: ∂ₜγ = −κ_ss N at interior times. -/
theorem IsFlow.lift_isCDF {T : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ) :
    IsCDF (Ioo 0 T) (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2))
/-- (T5) the analytic core, pure Mathlib. -/
theorem IsCDF.curv_sq_le {a b : ℝ} {γ : ℝ × ℝ → ℂ} (h : IsCDF (Ioo a b) γ) {t₀ K : ℝ}
    (ht₀ : t₀ ∈ Ioo a b) (hK : ∀ t ∈ Ico t₀ b, bendingEnergy (fun θ => γ (t, θ)) ≤ K) :
    ∃ C : ℝ, ∀ t ∈ Ico t₀ b, ∀ θ, curv (fun θ' => γ (t, θ')) θ ^ 2 ≤ C
```

with

```lean
structure IsCDF (I : Set ℝ) (γ : ℝ × ℝ → ℂ) : Prop where
  isOpen : IsOpen I
  smooth : ContDiffOn ℝ ∞ γ (I ×ˢ univ)
  periodic : ∀ t θ, γ (t, θ + 2 * π) = γ (t, θ)
  reg : ∀ t ∈ I, ∀ θ, deriv (fun θ' => γ (t, θ')) θ ≠ 0
  vel : ∀ t ∈ I, ∀ θ, HasDerivAt (fun s => γ (s, θ)) (-(κ_ss) • N) t   -- κ_ss, N as in `flowVel_eq_neg_kss`
```

(T3) at n = 1 [hand]: sums over `Fin 3 → Fin 1` have one term, so `aNormSq = (g⁻¹)³ C²` with `g = ‖X'‖²`, `C = ⟪X'', i X'⟫`, which is κ². The 1×1 inverse is `Matrix.inv_subsingleton` / `Matrix.det_fin_one` [GUESS on the exact simp set]. The transfer from the stereographic chart to θ uses `aNormSq_comp` with τ = `chartFwd ι1 z ∘ Circle.exp` (as a map `E1 → E1` after the identification ℝ ≃ E1); surjectivity of `Circle.exp` means no local inverse of `exp` is needed.

## 2. Proof route

### 2.1 Choice

DKS run an induction over all E_m = ∫|∂_s^m κ|² with Gagliardo–Nirenberg and get the rate (T−t)^{-1/4}. We need neither the rate nor higher m, because `flow_extends_of_bounded_A` already turns a sup bound on |A| into the extension; in the repo that is where the higher-order estimates live.

**Minimal route: bound E₁ = ∫κ_s² ds only.** Then sup κ² follows from 1D Sobolev plus a crude length lower bound. Checked by hand below; every inequality is elementary.

Not needed: length monotonicity (`volM_antitoneOn` exists but is unused), Fenchel's 2π, Grönwall, `Interp/Hamilton*`, `Interp/BallInterp`, `Extension/GlobalIneq`, time regularity at t = 0.

### 2.2 The mathematics [hand]

Notation at fixed t: v = |γ_θ|, ds = v dθ, ∂_s = v⁻¹∂_θ, T = γ_θ/v, N = iT, T_s = κN, N_s = −κT. Flow: γ_t = V N with V = −κ_ss. Norms over one period: A₀² = ∫κ² ds, a_j² = ∫(∂_s^j κ)² ds, L = ∫ds.

**(E) Pointwise evolution**, for any smooth normal velocity V:
1. v_t = −κ V v. From γ_θt = γ_tθ = v(V_s N − κ V T).
2. [∂_t, ∂_s] f = κ V f_s.
3. T_t = V_s N, N_t = −V_s T.
4. κ_t = V_ss + κ² V. From the N-component of T_θt = T_tθ.

With V = −κ_ss: v_t = κκ_ss v, (f_s)_t = (f_t)_s − κκ_ss f_s, κ_t = −κ_ssss − κ²κ_ss.

**(I) Energy identity.** Using (E), periodic integration by parts and κκ_s²κ_ss = κ(κ_s³/3)_s:

  d/dt ∫κ_s² ds = −2∫κ_sss² ds + 2∫κ²κ_ss² ds + (1/3)∫κ_s⁴ ds.

**(P) Periodic inequalities**, f, g smooth and 2π-periodic in θ:
- (P1) ∫ f_s g ds = −∫ f g_s ds.
- (P2) Cauchy–Schwarz in L²(ds).
- (P3) If f has a zero then sup f² ≤ 2 ‖f‖‖f_s‖. Applies to κ_s and κ_ss, which vanish at a maximum of κ, κ_s (Rolle).
- (P4) There is θ₀ with κ(θ₀)² ≤ A₀²/L, hence sup κ² ≤ A₀²/L + 2A₀a₁.
- (P5) L·A₀² ≥ 2. Proof: ∫T ds = γ(2π) − γ(0) = 0, so for any θ₀ there is θ₁ ∈ [θ₀, θ₀+2π] with ⟨T(θ₁), T(θ₀)⟩ ≤ 0, so √2 ≤ |T(θ₁) − T(θ₀)| ≤ ∫|κ| ds ≤ √L·A₀.

So **sup κ² ≤ A₀⁴/2 + 2A₀a₁**.

**(C) Closing the energy inequality.** By (P1), (P2): a₁² ≤ A₀a₂ and a₂² ≤ a₁a₃. Hence, with no rational exponents, a₂³ ≤ A₀a₃² and a₁³ ≤ A₀²a₃. By (P3): 2∫κ²κ_ss² ≤ 4A₀²a₂a₃ and (1/3)∫κ_s⁴ ≤ (2/3)a₁³a₂ ≤ (2/3)A₀²a₂a₃. With p = a₂a₃, x = a₃, c = (14/3)A₀²: p³ ≤ A₀x⁵ and

  d/dt E₁ ≤ −2x² + c·p ≤ C(A₀) := 1 + c³A₀R⁵, R = c³A₀/8

(two cases x ≥ R, x < R; `nlinarith`-level). So with A₀² ≤ K on [t₀, T): E₁(t) ≤ E₁(t₀) + C(K)(t − t₀), and

  sup κ²(t) ≤ K²/2 + 2√K·√(E₁(t₀) + C(K)(T − t₀)) for t ∈ [t₀, T).

Scaling check: under γ ↦ λγ, t ↦ λ⁴t, both sides of (I) scale as λ⁻⁷ and C(A₀) ∝ A₀¹⁴ scales as λ⁻⁷. Consistent.

**(A) Assembly of (T2).** Everything is at interior times, so no one-sided calculus at t = 0.
1. Shift: `Φ(t₀ + ·)` is an `IsFlow` on every `[0, T'−t₀]` (`IsFlow.shift_gen`, N1.lean:660, generic, proved; or `rst_shiftR`).
2. (T4), (T5), (T3): `aNormSq ≤ C` for the shifted flow on `Ico 0 (T−t₀)`.
3. `flow_extends_of_bounded_A` gives `Θ` on `[0, T−t₀+ε]` agreeing with `Φ(t₀+·)` on `[0, T−t₀)`.
4. `rst_glue` with `a = t₀`, any overlap `s ∈ (0, T−t₀)`, gives `Ψ` on `[0, T+ε]`, equal to `Φ` on `[0, t₀+s]`; on `(t₀+s, T)` equality comes from step 3 (or from `flow_unique_of_immersion`).

Recommended before any Lean: a sympy check of (E1)–(E4) and of the pointwise integrand manipulations in (I) (project convention; about 1 h, Sonnet).

## 3. Inventory

### 3.1 Reusable from the sibling repo [read: names and statements]

| Need | Lemma | File |
|---|---|---|
| Maximal flow, uniqueness | `exists_maximal_flow`, `flow_unique_of_immersion` | `Flow/MaximalFlow.lean`, `Flow/UniqueMain.lean` |
| Bounded \|A\| ⇒ extends | `flow_extends_of_bounded_A` | `Flow/ExtensionTheorem.lean` |
| Gluing with overlap, restriction | `rst_glue`, `rst_shiftR`, `IsFlow.mono` | `Flow/Restart.lean`, `Flow/GaugeEvolution.lean:37` |
| Chart-level joint smoothness and velocity of an `IsFlow` | `fc_contDiffOn_chart_of_isFlow`, `fc_vel_of_isFlow`, `IsFlow.isFlowM` | `Extension/FlowClosed.lean:25,41`, `Extension/FlowBridge.lean:176` |
| `velX` ↔ `flowVel`, linear change of domain | `velX_chartInv_eq_flowVel`, `flowVel_eq_velX`, `fb_flowVel_comp`, `flowVel_comp_linearEquiv` | `Extension/FlowBridge.lean`, `Extension/ImmPot.lean:232`, `Scalar/LinearChange.lean:110` |
| Reparametrisation invariance | `aNormSq_comp`, `velX_comp`, `flowVel_comp'`, `flowVel_congr_of_eventuallyEq`, `aNormSq_congr` | `Extension/ImmReparam.lean:428`, `Extension/PotReparam.lean:188`, `Scalar/Smooth.lean:55`, `Flow/TanAlgebra.lean:77`, `Extension/ImmSmooth.lean:132` |
| `IsLagImm ι` ↔ `IsLagImmersion` | `isLagImm_iff_isLagImmersion` | `Extension/FlowBridge.lean:123` |
| Pattern for d/dt under an integral (tube lemma + dominated convergence) | `hasDerivAt_integral_mul_slice` (n-dim, compactly supported weight; a pattern to copy, not directly applicable) | `Extension/SliceSmooth.lean:104` |
| n = 1 velocity and setting | `flowVel_eq_neg_kss`, `Jc*`, `isLagImmersion_iff`, `IsFlow.shift_gen` | old repo `N1.lean` (to port) |

From Mathlib [names read in the pinned Mathlib unless marked]: `contMDiff_circleExp`, `contMDiff_coe_sphere` (Sphere.lean:595, 424); `intervalIntegral.integral_mul_deriv_eq_deriv_mul` (IBP); `intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le` (ParametricIntervalIntegral.lean:97); `Function.Periodic.intervalIntegral_add_eq`; `exists_deriv_eq_zero` (Rolle); `Convex.image_sub_le_mul_sub_of_deriv_le` (MeanValue.lean:354); `ContDiffAt.isSymmSndFDerivAt` (mixed partials); `MeasureTheory.exists_le_setAverage`; `Matrix.inv_subsingleton`, `Matrix.det_fin_one`; `Circle.periodic_exp` [GUESS].

### 3.2 Not reusable, although it looks related

- `Interp/Hamilton*.lean`, `BallInterp`, `ThreeOrders`, `SupBound`: interpolation on ℝⁿ with cutoffs on balls and jets `jet j ξ`; statement shapes do not match periodic arclength integrals. The 1D versions (P1)–(P5) are shorter to prove than to adapt.
- `Extension/IntM.lean` (`intM`, `IntAtlas`, natural densities), `VolM`, `energyM`: integration on M through a partition of unity. Using it for ∫κ² ds would give a statement the paper cannot read off easily; the interval integral over the lift is simpler and explicit.
- `Extension/EvolPointwise`, `LocalEvolution`, `GlobalIneq`, `EnergyBound`: evolution of |∇^k C|² in graph charts *under an assumed bound on |A|*. The present task is to produce that bound.
- `Interp/EnergyGronwall`: not needed (the inequality is E₁' ≤ C).

### 3.3 New

1. Port of N1.lean (sorry-free part).
2. Lift: `Circle.exp` is a smooth immersion in charts; `IsFlow ⇒ IsCDF` (T4); `aNormSq = curv²` (T3).
3. Arclength calculus for a time-dependent periodic immersed curve: joint smoothness of v, T, κ, ∂_s^j κ; (E1)–(E4); d/dt ∫ f ds = ∫ (f_t + κκ_ss f) ds.
4. Periodic inequalities (P1)–(P5) and the real-number lemma (C).
5. Energy identity (I), the bound (T5), assembly (T1), (T2).

## 4. Work packages

New directory `CWFlow/Curve/` in the sibling repo. New files only; `CWFlow.lean` and `Extension/All.lean` are **not** edited (the user adds imports when he wants). W2–W4 import only Mathlib, so they compile fast and cannot be broken by the user's concurrent edits.

| WP | File | Content | Imports | Lines [GUESS] | Agent-h [GUESS] | Model |
|---|---|---|---|---|---|---|
| W0 | (sympy, in `lean_bridge/`) | check (E1)–(E4), integrand algebra of (I) | – | 150 py | 1 | Sonnet |
| W1 | `Curve/N1.lean` | port of N1.lean minus `glue`; `cdf_exists_maximal_A` = `exists_maximal_flow` at n = 1 (instantiation test) | CWFlow | 800 (copied) + 60 | 2–4 | Sonnet |
| W2 | `Curve/Periodic.lean` | `curv`, `bendingEnergy`, weighted periodic IBP (P1), C–S (P2), (P3), (P4), (P5), lemma (C) | Mathlib | 450–650 | 6–9 | Sonnet, Opus for (P5) |
| W3 | `Curve/Arclength.lean` | `IsCDF`; joint smoothness; (E1)–(E4); d/dt of ∫ f ds | Mathlib | 550–850 | 9–15 | Opus |
| W4 | `Curve/Energy.lean` | identity (I), E₁' ≤ C(K), (T5) | W2, W3 | 300–450 | 5–8 | Opus |
| W5 | `Curve/Lift.lean` | `Circle.exp` immersion in charts, (T4), (T3) | W1, `Extension/FlowClosed`, `ImmReparam` | 350–550 | 6–10 | Opus |
| W6 | `Curve/Criterion.lean` | (T2), (T1), `#print axioms` guards | W4, W5, `Flow/MaximalFlow`, `Flow/Restart` | 150–250 | 3–5 | Opus |
| W7 | review | independent read of statements against `U1_regularity`/`continuation`/C-1; text patch list | – | – | 2–3 | Fable |

Total: about 1,900–2,900 new lines (plus the 800-line port), 30–50 agent-hours plus review. W2 and W3 are independent of each other and of W1; W5 needs W1.

### Order and kill points

- **K0 (W0, 1 h):** sympy disagrees with (E) or (I) ⇒ stop and fix the math first.
- **K1 (W1, ≤ 4 h):** the port compiles in the sibling repo and `exists_maximal_flow` instantiates at `Circle`. If instance or universe problems cannot be solved in 4 h, stop and report; everything else depends on this. Expected to pass, since the same variable block already instantiated in the old copy.
- **K2 (W3, at 10 h):** (E4) κ_t not proved ⇒ reassess. Fallback: take the evolution of κ and v as fields of a larger hypothesis structure and prove them later (the theorem would then be conditional, which the repo's rule 2 allows only as explicit `Prop` hypotheses; not acceptable as the final state).
- **K3 (W5, at 8 h):** chart friction with `Circle` (stereographic charts, `𝓡 1` vs `𝓘(ℝ, E1)`, product model `𝓘(ℝ,ℝ).prod 𝓘(ℝ,E1)`). Fallback: state (T3)/(T4) for an arbitrary smooth covering map `p : ℝ → M` with nonvanishing chart derivative, and prove the `Circle.exp` instance separately.
- W2, W4, W6 carry low risk.

### Risks

1. **(E) in Lean (medium).** Mixed partials for ℂ-valued maps on `ℝ × ℝ`, slice derivatives vs `fderiv` on the product, real inner product on ℂ. Tedious, not deep. Mitigation: prove (E) for general normal speed V with a small private API (`∂θ`, `∂t` as `fderiv γ q (0,1)`, `(1,0)`), and never expand κ in θ-derivatives.
2. **Circle charts (medium).** N1.lean deliberately stayed on `ℝ → ℂ` and deferred exactly this bridge.
3. **Port drift (low–medium).** The sibling repo changed after d72b744 (gauge refactor, user's uncommitted edits). If the user edits `Flow/` or `Extension/` files during the work, W1/W5/W6 oleans go stale; W2–W4 are immune.
4. **Build resources (practical).** One Lean at a time via `scripts/compile1.sh` (lock shared with `../SLAG master`, hard memory limit in `../tools/leanjob.py`), `LEAN_NUM_THREADS=4`, no `lake update` (`.lake/packages` is a junction to `../MM/.lake/packages`). Repo files with the full import closure need about 3.5 GB and about 7 min each (N1 report: 7.5 min); Mathlib-only files about 35 s–2 min. This exceeds the 800 MB/process preference recorded for the proof tree, so Lean compiles need Micah's explicit OK and must queue behind his own Lean jobs.
5. **Statement fidelity (low).** `IsFlow` asks joint smoothness on `[0,T] × S¹` including t = 0, and velocity with no tangential part; (T1) is about that class. This is the class U1 already uses (Lemma "n = 1: the Lean velocity is CDF").
6. **Estimates.** Line and hour figures are guesses from the repo's own rate (about 29k lines of `Extension/` + `Interp/` in a few days) and from the N1 sorry-closing effort. Uncertainty about ×1.5 either way.

### What the paper gains, and what it does not

- Lemma "blow-up criterion, the form needed" becomes [Lean] via (T1); Lemmas "maximal smooth solution" and "gluing" become [Lean] via `exists_maximal_flow`, `rst_glue`. The author item A-2 (meaning of "maximal solution" in DKS) disappears.
- C-1's use of DKS (D1) can be replaced by (T2); `C1_identification_v41.tex` §(I5) and the two "DKS continuation steps" would need rewording, and it must be checked that C-1 needs nothing beyond "extends as an `IsFlow`" (the C-4 completion, uniform ‖∂_s^m κ‖_∞, is implied by an `IsFlow` on the closed interval). Not checked here beyond reading l.49 and l.172.
- Still [hand] after this: the identification of the step-4 solution with the Lean flow (C-1), (P2′) modulation, the symmetry lemma's use (Lean `cdf_symmetric` exists in N1).
- "[Lean]" would then mean compiled, with `#print axioms` clean, once W1–W6 build in the sibling repo.

## 5. Recommendation

**Go**, in this order: W0 → W1 (K1) → W2 ∥ W3 → W4 → W5 → W6 → W7. Before starting, Micah should confirm: (a) new directory `CWFlow/Curve/` in the sibling repo is acceptable (or name another place); (b) Lean compiles from agent sessions may use the shared lock and about 3.5 GB; (c) whether N1.lean should be ported by copy (proposed) or whether he prefers to move it himself.

Optional upstream simplification (the user's call, not required): `isFlow_closed_of_bound` uses its bound only on `Ioo 0 T` (`Extension/FlowClosed.lean:91`); the route above does not depend on changing it.
