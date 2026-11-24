/-
CWFlowN1.lean  --  the n = 1 instance of CWFlow (curves in ℂ ≅ ℝ²): curve diffusion flow.

DRAFT, NEVER ELABORATED. Written on 2026-11-21 against HSLflow-Whitney commit
d72b7447edb0c830331a76152559902c23a4092e (CWFlow/, toolchain leanprover/lean4:v4.33.1,
Mathlib v4.33.1).  No Lean toolchain was available when this was written; every
Mathlib lemma name below that is not quoted from the CWFlow source is marked [GUESS].
See INSTRUCTIONS.md (same directory) for where to put it and how to build it.

What is proved here (modulo compiling) and what is `sorry`:
  (0)  setting: E1 = EuclideanSpace ℝ (Fin 1), M = Circle (unit circle in ℂ), V = ℂ, J = (I * ·);
       hJi, hJJ, hdim, and "every immersion of a 1-manifold is Lagrangian"      -- proved
  (i)  existence  `cdf_exists`                     -- proved from `exists_flow_of_immersion`
  (ii) uniqueness `cdf_unique`                     -- proved from `flow_unique_of_immersion`
  (ii') geometric uniqueness up to reparametrization `cdf_unique_up_to_reparam`
                                                   -- proved from `exists_flow_of_normal_family`
                                                      + `flow_unique_of_immersion`
  (iii) n = 1 velocity identity  flowVel = -k_ss N  `flowVel_eq_neg_kss`   -- SORRY (S1)
  (iv) symmetry: conj ∘ Φ_t ∘ ι = Φ_t            `cdf_symmetric`  -- proved from (ii) and
       the closure lemma `IsFlow.conj_reflect`, which is SORRY (S3); it rests on
       `flowVel_conj` (SORRY (S2)) and the repo's `flowVel_comp`.
  (v)  optional extras used as hand arguments in U1_regularity.tex: time shift (SORRY (S4)),
       gluing (SORRY (S5)).  `IsFlow.mono` is already in the repo.

All `sorry`s are tagged `-- SORRY (Sk)`; `grep -n "SORRY" CWFlowN1.lean` lists them.
-/
import CWFlow.Flow.Main               -- exists_flow_of_immersion
import CWFlow.Flow.UniqueMain         -- flow_unique_of_immersion
import CWFlow.Flow.NormalFlow         -- exists_flow_of_normal_family (imports GraphUnique: IsFlow)
import CWFlow.Flow.GaugeEvolution     -- IsFlow.mono
import CWFlow.Scalar.Reparam          -- flowVel_comp (used only in the sketch of S3)
import Mathlib.Geometry.Manifold.Instances.Sphere   -- [GUESS: module path] charts on `Circle`

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold ComplexConjugate

noncomputable section

namespace CWFlowN1

open CWFlow

/-! ## 0. The setting -/

/-- The model space of a 1-manifold.  Mathlib charts `Circle` on `EuclideanSpace ℝ (Fin 1)`
(`𝓡 1` is notation for `𝓘(ℝ, EuclideanSpace ℝ (Fin 1))`). -/
abbrev E1 := EuclideanSpace ℝ (Fin 1)

-- `[NeZero (Module.finrank ℝ E)]` is the only dimension hypothesis of the two main theorems:
--   Main.lean l.22-23:  variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
--                         [NeZero (Module.finrank ℝ E)]
instance instNeZeroFinrankE1 : NeZero (Module.finrank ℝ E1) :=
  ⟨by simp [E1]⟩   -- [GUESS] uses simp lemma `finrank_euclideanSpace_fin`; else
                   -- `⟨by rw [finrank_euclideanSpace_fin]; exact one_ne_zero⟩`

/-! Instance probes.  The manifold hypotheses of the main theorems are (Main.lean l.24-25):
    `[TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
     [CompactSpace M] [T2Space M] [SigmaCompactSpace M]`.
  If one of these probes fails, see INSTRUCTIONS.md, "likely failures", item 2. -/
example : ChartedSpace E1 Circle := inferInstance
example : IsManifold 𝓘(ℝ, E1) ∞ Circle := inferInstance
  -- [GUESS] Mathlib may only provide `IsManifold (𝓡 1) ω Circle`; then add
  -- `instance : IsManifold 𝓘(ℝ, E1) ∞ Circle := IsManifold.of_le (n := ω) le_top`
example : CompactSpace Circle := inferInstance
example : T2Space Circle := inferInstance
example : SigmaCompactSpace Circle := inferInstance
example : CompleteSpace E1 := inferInstance   -- needed by `exists_flow_of_normal_family`

/-- `J` = multiplication by `i` on `ℂ`, as an `ℝ`-linear map. -/
def Jc : ℂ →L[ℝ] ℂ := ContinuousLinearMap.mul ℝ ℂ Complex.I

@[simp] theorem Jc_apply (v : ℂ) : Jc v = Complex.I * v := by
  simp [Jc]   -- [GUESS] `ContinuousLinearMap.mul_apply'` is simp; `rfl` may also work

/-- `hJi` of the main theorems: `J` is an isometry.  Proof by real polarization. -/
theorem Jc_isom (u v : ℂ) : ⟪Jc u, Jc v⟫ = ⟪u, v⟫ := by
  simp only [Jc_apply]
  -- [GUESS: lemma name] ⟪x, y⟫_ℝ = (‖x + y‖ * ‖x + y‖ - ‖x‖ * ‖x‖ - ‖y‖ * ‖y‖) / 2
  rw [real_inner_eq_norm_add_mul_self_sub_norm_mul_self_sub_norm_mul_self_div_two,
    real_inner_eq_norm_add_mul_self_sub_norm_mul_self_sub_norm_mul_self_div_two, ← mul_add]
  simp   -- uses ‖I‖ = 1 (`Complex.norm_I`)
  -- fallback: `simp [Complex.inner, map_mul, Complex.conj_I]; ring_nf; simp [Complex.I_sq]`

/-- `hJJ` of the main theorems: `J² = -1`. -/
theorem Jc_JJ (v : ℂ) : Jc (Jc v) = -v := by
  simp only [Jc_apply]
  rw [← mul_assoc, Complex.I_mul_I]
  ring

/-- `hdim` of the main theorems: `dim V = 2 dim E`, here `2 = 2 · 1`. -/
theorem hdim1 : Module.finrank ℝ ℂ = 2 * Module.finrank ℝ E1 := by
  simp [E1, Complex.finrank_real_complex]

/-- `⟪J v, v⟫ = 0` (from `hJi`, `hJJ` only). -/
theorem inner_Jc_self (v : ℂ) : ⟪Jc v, v⟫ = 0 := by
  have h1 : ⟪Jc v, v⟫ = ⟪Jc (Jc v), Jc v⟫ := (Jc_isom _ _).symm
  rw [Jc_JJ, inner_neg_left, real_inner_comm] at h1
  -- h1 : ⟪v, Jc v⟫ = -⟪v, Jc v⟫   (rw of `real_inner_comm` hits the first inner product)
  have h2 : ⟪v, Jc v⟫ = 0 := by linarith
  rw [real_inner_comm]
  exact h2

/-- In dimension one every linear map is Lagrangian: `⟪J (L u), L w⟫ = 0`. -/
theorem lag_dim_one (L : E1 →L[ℝ] ℂ) (u w : E1) : ⟪Jc (L u), L w⟫ = 0 := by
  by_cases hu : u = 0
  · simp [hu]
  · -- [GUESS: lemma name] finrank_eq_one_iff_of_nonzero' (v) (nz : v ≠ 0) :
    --   finrank K V = 1 ↔ ∀ w : V, ∃ c : K, c • v = w
    obtain ⟨c, hc⟩ := (finrank_eq_one_iff_of_nonzero' (K := ℝ) u hu).mp
      (by simp [finrank_euclideanSpace]) w
    rw [← hc, map_smul, real_inner_smul_right, inner_Jc_self, mul_zero]

/-- A smooth immersed closed curve `F : S¹ → ℂ`, chartwise (the first two clauses of
`IsLagImmersion`; quoted from Weinstein/Global.lean l.53-58):
```
structure IsLagImmersion (J : V →L[ℝ] V) (F : M → V) : Prop where
  smooth : ∀ z : M, ContDiffOn ℝ ∞ (chartRep E F z) (chart E z).target
  inj : ∀ z : M, ∀ y ∈ (chart E z).target, Function.Injective (fderiv ℝ (chartRep E F z) y)
  lag : ∀ z : M, ∀ y ∈ (chart E z).target, ∀ u w,
    ⟪J (fderiv ℝ (chartRep E F z) y u), fderiv ℝ (chartRep E F z) y w⟫ = 0
```
-/
structure IsCurveImmersion (F : Circle → ℂ) : Prop where
  smooth : ∀ z : Circle, ContDiffOn ℝ ∞ (chartRep E1 F z) (chart E1 z).target
  inj : ∀ z : Circle, ∀ y ∈ (chart E1 z).target, Function.Injective (fderiv ℝ (chartRep E1 F z) y)

theorem IsCurveImmersion.isLag {F : Circle → ℂ} (h : IsCurveImmersion F) :
    IsLagImmersion E1 Jc F :=
  ⟨h.smooth, h.inj, fun _ _ _ u w => lag_dim_one _ u w⟩

/-- At `n = 1` "Lagrangian immersion" is just "immersion". -/
theorem isLagImmersion_iff {F : Circle → ℂ} : IsLagImmersion E1 Jc F ↔ IsCurveImmersion F :=
  ⟨fun h => ⟨h.smooth, h.inj⟩, fun h => h.isLag⟩

/-! ## (i) Existence

Source, Flow/Main.lean l.29-36:
```
theorem exists_flow_of_immersion {J : V →L[ℝ] V} {F : M → V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) :
    ∃ T > 0, ∃ Φ : ℝ → M → V, Φ 0 = F ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ) ∧
      (∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)) ∧
      ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
        HasDerivWithinAt (fun s => Φ s x) (flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t
```
The three clauses are exactly the fields of `IsFlow` (Flow/GraphUnique.lean l.42-46):
```
structure IsFlow (J : V →L[ℝ] V) (T : ℝ) (Φ : ℝ → M → V) : Prop where
  smooth : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ)
  imm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)
  vel : ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
    HasDerivWithinAt (fun s => Φ s x) (flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t
```
-/

/-- **(i) Short-time existence of CDF** from a smooth immersed closed curve. -/
theorem cdf_exists {F : Circle → ℂ} (hF : IsCurveImmersion F) :
    ∃ T > 0, ∃ Φ : ℝ → Circle → ℂ, Φ 0 = F ∧ IsFlow E1 Jc T Φ := by
  obtain ⟨T, hT, Φ, h0, hsm, himm, hvel⟩ :=
    exists_flow_of_immersion (E := E1) Jc_isom Jc_JJ hF.isLag hdim1
  exact ⟨T, hT, Φ, h0, ⟨hsm, himm, hvel⟩⟩

/-! ## (ii) Uniqueness

Source, Flow/UniqueMain.lean l.28-32:
```
theorem flow_unique_of_immersion {J : V →L[ℝ] V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {T : ℝ} {Φ₁ Φ₂ : ℝ → M → V}
    (hΦ₁ : IsFlow E J T Φ₁) (hΦ₂ : IsFlow E J T Φ₂) (h0 : Φ₁ 0 = Φ₂ 0) :
    ∀ t ∈ Icc (0 : ℝ) T, Φ₁ t = Φ₂ t
```
Note: this is equality of the maps, not "up to reparametrization": an `IsFlow` has velocity
exactly `flowVel` (purely normal at n = 1, by (iii)), so there is no gauge freedom left. -/

/-- **(ii) Uniqueness of CDF** among `IsFlow` solutions (purely normal velocity). -/
theorem cdf_unique {T : ℝ} {Φ₁ Φ₂ : ℝ → Circle → ℂ}
    (h₁ : IsFlow E1 Jc T Φ₁) (h₂ : IsFlow E1 Jc T Φ₂) (h0 : Φ₁ 0 = Φ₂ 0) :
    ∀ t ∈ Icc (0 : ℝ) T, Φ₁ t = Φ₂ t :=
  flow_unique_of_immersion Jc_isom Jc_JJ hdim1 h₁ h₂ h0

/-! ## (ii') Geometric uniqueness: up to reparametrization

Source, Flow/NormalFlow.lean l.282-291:
```
theorem exists_flow_of_normal_family {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ} (hT : 0 < T)
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ α x : M, x ∈ (chart E α).source → ∃ c : E,
      dtW (flowRep E G α) T (t, chart E α x) =
        flowVel J (chartRep E (G t) α) (chart E α x) +
          fderiv ℝ (chartRep E (G t) α) (chart E α x) c) :
    ∃ χ : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M), (∀ x, χ 0 x = x) ∧
      IsFlow E J T (fun t x => G t (χ t x))
```
(`flowRep E Φ β p = Φ p.1 ((chart E β).symm p.2)`, GraphUnique.lean l.50;
 `dtW w T p = derivWithin (fun s => w (s, p.2)) (Icc 0 T) p.1`, Scalar/GraphTime.lean l.24.) -/

/-- **(ii') Any smooth family of immersed curves moving by `−k_ss N` plus a tangential term,
from the same initial curve as an `IsFlow` `Φ`, is `Φ` up to a family of diffeomorphisms
`χ t` of `S¹` with `χ 0 = id`:** `G t ∘ χ t = Φ t`. -/
theorem cdf_unique_up_to_reparam {T : ℝ} (hT : 0 < T) {G Φ : ℝ → Circle → ℂ}
    (hsm : ∀ α : Circle, ContDiffOn ℝ ∞ (flowRep E1 G α) (Icc 0 T ×ˢ (chart E1 α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsCurveImmersion (G t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ α x : Circle, x ∈ (chart E1 α).source → ∃ c : E1,
      dtW (flowRep E1 G α) T (t, chart E1 α x) =
        flowVel Jc (chartRep E1 (G t) α) (chart E1 α x) +
          fderiv ℝ (chartRep E1 (G t) α) (chart E1 α x) c)
    (hΦ : IsFlow E1 Jc T Φ) (h0 : G 0 = Φ 0) :
    ∃ χ : ℝ → (Circle ≃ₘ⟮𝓘(ℝ, E1), 𝓘(ℝ, E1)⟯ Circle), (∀ x, χ 0 x = x) ∧
      ∀ t ∈ Icc (0 : ℝ) T, (fun x => G t (χ t x)) = Φ t := by
  obtain ⟨χ, hχ0, hflow⟩ :=
    exists_flow_of_normal_family hT hsm (fun t ht => (himm t ht).isLag) hvel
  refine ⟨χ, hχ0, fun t ht => ?_⟩
  have h00 : (fun t x => G t (χ t x)) 0 = Φ 0 := by
    funext x
    show G 0 (χ 0 x) = Φ 0 x
    rw [hχ0]
    exact congrFun h0 x
  exact cdf_unique hflow hΦ h00 t ht

/-! ## (iii) At n = 1 the velocity `J∇div(JH)` is `−k_ss N`

Source, Scalar/Pointwise.lean l.81-90:
```
def meanCurv (X : E → V) (y : E) : V := HOf (fderiv ℝ X y) (fderiv ℝ (fderiv ℝ X) y)
def flowPot (J : V →L[ℝ] V) (X : E → V) (y : E) : ℝ :=
  trOf (fderiv ℝ X y) fun a c =>
    ⟪fderiv ℝ (fun z => J (meanCurv X z)) y a, fderiv ℝ X y c⟫
def flowVel (J : V →L[ℝ] V) (X : E → V) (y : E) : V :=
  J (fderiv ℝ X y (ginvOf (fderiv ℝ X y) (fderiv ℝ (flowPot J X) y)))
```
Hand computation (also U1_sym.py S1, exact sympy): with `T = X'/|X'|`, `N = J T`,
`k = ⟪X'', J X'⟫/|X'|³`, `∂_s = |X'|⁻¹ ∂_θ`:  `H = k N`, `JH = −k T`, `flowPot = −k_s`,
`flowVel = J (T · (−k_ss)) = −k_ss N`.  No tangential component.

The statement is made for `E = ℝ` (a curve `X : ℝ → ℂ`), where `deriv` is available.  It transfers
to the charts `E1 = EuclideanSpace ℝ (Fin 1)` by the repo's invariance lemmas
(Scalar/LinearChange.lean `flowVel_comp_linearEquiv`, Scalar/Reparam.lean `flowVel_comp`). -/

/-- Signed curvature `k = ⟪X'', i X'⟫ / |X'|³` of a plane curve. -/
def curv (X : ℝ → ℂ) (θ : ℝ) : ℝ := ⟪deriv (deriv X) θ, Complex.I * deriv X θ⟫ / ‖deriv X θ‖ ^ 3

/-- Arclength derivative `∂_s f = f' / |X'|`. -/
def dS (X : ℝ → ℂ) (f : ℝ → ℝ) (θ : ℝ) : ℝ := deriv f θ / ‖deriv X θ‖

/-- **(iii) n = 1: `flowVel = −k_ss N`, `N = i T`.**  NOT PROVED (needs a jet computation of
`HOf`, `ginvOf`, `trOf` at `E = ℝ`, i.e. `trOf L B = B 1 1 / ‖L 1‖²`, plus three `fderiv`
chain rules); the identity itself is checked exactly in U1_sym.py (S1). -/
theorem flowVel_eq_neg_kss {X : ℝ → ℂ} (hX : ContDiff ℝ ∞ X) (hreg : ∀ θ, deriv X θ ≠ 0) (θ : ℝ) :
    flowVel Jc X θ =
      -((dS X (dS X (curv X)) θ : ℝ) : ℂ) * (Complex.I * (deriv X θ / (‖deriv X θ‖ : ℂ))) := by
  sorry -- SORRY (S1): n = 1 velocity identity (statement checked by sympy, not in Lean)

/-! ## (iv) Symmetry under reflection, via uniqueness

`σ = conj` on `ℂ` (an `ℝ`-linear isometry with `σ J = −J σ`), `ι p = p⁻¹ = conj p` on `Circle`.
Neither the closure of `IsFlow` under `σ ∘ · ∘ ι` nor the pointwise identity `flowVel_conj`
is in the repo (U1_regularity.tex already says "it is not a Lean lemma"). -/

/-- `flowVel (σ ∘ X) = σ (flowVel X)` for `σ = conj` (J appears twice, so the sign cancels:
`H' = σH`, `JH' = −σJH`, `flowPot' = −flowPot`, `flowVel' = J σ L g⁻¹ D(−pot) = σ flowVel`).
NOT PROVED; U1_sym.py S2 checks it at n = 1. -/
theorem flowVel_conj {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
    {X : E → ℂ} {U : Set E} (hU : IsOpen U) (hX : ContDiffOn ℝ ∞ X U)
    (hinj : ∀ y ∈ U, Function.Injective (fderiv ℝ X y)) {y : E} (hy : y ∈ U) :
    flowVel Jc (fun z => conj (X z)) y = conj (flowVel Jc X y) := by
  sorry -- SORRY (S2): conjugation equivariance of flowVel (sympy S2; not in Lean)

/-- Closure of `IsFlow` under `Φ_t ↦ conj ∘ Φ_t ∘ ι`, `ι p = p⁻¹`.  NOT PROVED.
Sketch: (smooth) `conj` is a continuous ℝ-linear map and `(t, p) ↦ (t, p⁻¹)` is smooth
(`Circle` is a Lie group, [GUESS] `contMDiff_inv`); (imm) chart reps are
`conj ∘ chartRep (Φ t) (ι z) ∘ τ`, `τ = chart (ι z) ∘ ι ∘ (chart z).symm` with invertible
derivative, so smooth/injective; `lag` is automatic (`lag_dim_one`); (vel) `HasDerivWithinAt`
composes with the CLM `conj`, and the velocity matches by `flowVel_conj` (S2) and the repo's
`flowVel_comp` (Scalar/Reparam.lean l.118-120):
```
theorem flowVel_comp {z : E} (hz : z ∈ N)
    (hpot : DifferentiableAt ℝ (flowPot J X) (τ z)) :
    flowVel J (X ∘ τ) z = flowVel J X (τ z)
```
(`flowVel_comp` needs only an invertible derivative of `τ`, so orientation reversal is fine.) -/
theorem IsFlow.conj_reflect {T : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ) :
    IsFlow E1 Jc T (fun t p => conj (Φ t p⁻¹)) := by
  sorry -- SORRY (S3): closure of IsFlow under conj ∘ · ∘ ι (plumbing + S2)

/-- **(iv) A flow from a reflection-symmetric curve stays symmetric:**
`conj (Φ t p⁻¹) = Φ t p` for all `t ∈ [0, T]`.  Proved from (ii) and (S3). -/
theorem cdf_symmetric {T : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ)
    (hsym : ∀ p : Circle, conj (Φ 0 p⁻¹) = Φ 0 p) :
    ∀ t ∈ Icc (0 : ℝ) T, ∀ p : Circle, conj (Φ t p⁻¹) = Φ t p := by
  intro t ht p
  have h0 : (fun t p => conj (Φ t p⁻¹)) 0 = Φ 0 := funext hsym
  exact congrFun (cdf_unique (IsFlow.conj_reflect h) h h0 t ht) p

/-! ## (v) Optional: the hand arguments of U1_regularity.tex §2 as Lean targets

Restriction is in the repo (Flow/GaugeEvolution.lean l.34-35):
```
theorem IsFlow.mono {J : V →L[ℝ] V} {T δ : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) (hδ : δ ≤ T) :
    IsFlow E J δ Φ
```
Time shift and gluing are not. -/

example {T δ : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ) (hδ : δ ≤ T) : IsFlow E1 Jc δ Φ :=
  h.mono hδ

/-- Time shift (U1 "time shift", C1 Thm (v)).  NOT PROVED; expected easy
(`ContMDiffOn.comp` with `(t, x) ↦ (t + a, x)`; [GUESS] `HasDerivWithinAt.comp_add_const`). -/
theorem IsFlow.shift {T a : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ) (ha : 0 ≤ a)
    (haT : a ≤ T) : IsFlow E1 Jc (T - a) (fun t => Φ (t + a)) := by
  sorry -- SORRY (S4): time shift

/-- Gluing (U1 Lemma "gluing").  NOT PROVED and NOT cheap: joint C^∞ at the seam `t = a` needs
all one-sided `∂ₜ^m` to agree, which U1 gets from the equation by induction on `m`. -/
theorem IsFlow.glue {a b : ℝ} {Φ₁ Φ₂ : ℝ → Circle → ℂ} (h₁ : IsFlow E1 Jc a Φ₁)
    (h₂ : IsFlow E1 Jc b Φ₂) (ha : 0 ≤ a) (hb : 0 ≤ b) (h : Φ₂ 0 = Φ₁ a) :
    IsFlow E1 Jc (a + b) (fun t => if t ≤ a then Φ₁ t else Φ₂ (t - a)) := by
  sorry -- SORRY (S5): gluing

/-! ## Axiom audit (expected after a successful build)
`cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam` should depend only on
`[propext, Classical.choice, Quot.sound]`; `cdf_symmetric` additionally on `sorryAx` (through S3). -/
#print axioms cdf_exists
#print axioms cdf_unique
#print axioms cdf_unique_up_to_reparam
#print axioms cdf_symmetric

end CWFlowN1
