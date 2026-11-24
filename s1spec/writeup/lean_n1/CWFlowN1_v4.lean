/-
CWFlowN1.lean  --  the n = 1 instance of CWFlow (curves in ℂ ≅ ℝ²): curve diffusion flow.

STATUS 2026-11-23 (v4): builds with Lean v4.33.1 / Mathlib in C:\hsl\CWFlow; S1-S4 closed, only S5
(gluing) is `sorry`.  See N1_SORRY_REPORT.md.  Original draft note follows.
DRAFT (v1), written on 2026-11-21 against HSLflow-Whitney commit
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
  (iii) n = 1 velocity identity  flowVel = -k_ss N  `flowVel_eq_neg_kss`   -- (S1) PROVED 2026-11-23
  (iv) symmetry: conj ∘ Φ_t ∘ ι = Φ_t            `cdf_symmetric`  -- proved from (ii) and
       the closure lemma `IsFlow.conj_reflect` ((S3) PROVED 2026-11-23); it rests on
       `flowVel_conj` (S2, PROVED 2026-11-23, no hypotheses needed) and the repo's `flowVel_comp`.
  (v)  optional extras used as hand arguments in U1_regularity.tex: time shift (S4, PROVED 2026-11-23),
       gluing (SORRY (S5)).  `IsFlow.mono` is already in the repo.

All `sorry`s are tagged `-- SORRY (Sk)`; `grep -n "SORRY" CWFlowN1.lean` lists them.
-/
import CWFlow.Flow.Main               -- exists_flow_of_immersion
import CWFlow.Flow.UniqueMain         -- flow_unique_of_immersion
import CWFlow.Flow.NormalFlow         -- exists_flow_of_normal_family (imports GraphUnique: IsFlow)
import CWFlow.Flow.GaugeEvolution     -- IsFlow.mono
import CWFlow.Scalar.Reparam          -- flowVel_comp (used in S3)
import CWFlow.Scalar.Smooth           -- contDiffOn_flowPot (used in S3)
import CWFlow.Flow.DiffeoChart        -- diffeoRep, isLagImmersion_comp_diffeo (used in S3)
import CWFlow.Flow.TanAlgebra         -- flowVel_congr_of_eventuallyEq (used in S3)
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

/-! ### Dimension-one jet computations for (iii) -/

section Dim1
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem clm_real_apply {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (L : ℝ →L[ℝ] F) (a : ℝ) :
    L a = a • L 1 := by
  rw [← map_smul, smul_eq_mul, mul_one]

theorem injective_of_apply_one {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0) : Function.Injective L := by
  intro a b hab
  rw [clm_real_apply L a, clm_real_apply L b] at hab
  exact smul_left_injective ℝ hL hab

/-- `g⁻¹ ξ = ξ(1) / |L 1|²` in dimension one. -/
theorem ginvOf_dim1 {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0) (ξ : ℝ →L[ℝ] ℝ) :
    ginvOf L ξ = ξ 1 / ‖L 1‖ ^ 2 := by
  have hinj := injective_of_apply_one hL
  apply gOf_injective hinj
  rw [gOf_ginvOf hinj]
  ext
  have hn : ‖L 1‖ ≠ 0 := norm_ne_zero_iff.mpr hL
  rw [gOf_apply, clm_real_apply L (ξ 1 / ‖L 1‖ ^ 2), real_inner_smul_left,
    real_inner_self_eq_norm_sq]
  field_simp

/-- The trace in dimension one: `tr_g B = B(1,1) / |L 1|²` for `B` bilinear. -/
theorem trOf_dim1 {W : Type*} [AddCommGroup W] [Module ℝ W] {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0)
    (B : ℝ → ℝ → W) (hB : ∀ a c, B a c = (a * c) • B 1 1) :
    trOf L B = (‖L 1‖ ^ 2)⁻¹ • B 1 1 := by
  have hs : ∑ i, bas ℝ i * (dualBas i 1 / ‖L 1‖ ^ 2) = (‖L 1‖ ^ 2)⁻¹ := by
    have h1 := sum_dualBas_smul (E := ℝ) 1
    simp only [smul_eq_mul] at h1
    calc ∑ i, bas ℝ i * (dualBas i 1 / ‖L 1‖ ^ 2)
        = (∑ i, dualBas i 1 * bas ℝ i) / ‖L 1‖ ^ 2 := by
          rw [Finset.sum_div]; refine Finset.sum_congr rfl fun i _ => ?_; ring
      _ = (‖L 1‖ ^ 2)⁻¹ := by rw [h1, one_div]
  simp only [trOf, ginvOf_dim1 hL]
  simp_rw [hB _ (_ / _)]
  rw [← Finset.sum_smul, hs]

theorem fderiv_real_apply {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (f : ℝ → F) (x a : ℝ) :
    fderiv ℝ f x a = a • deriv f x := by
  rw [clm_real_apply]; rfl


/-- The plane identity behind `JH = -k T`. -/
theorem I_mul_perp (a v : ℂ) (hv : v ≠ 0) :
    Complex.I * ((‖v‖ ^ 2)⁻¹ • (a - (⟪a, v⟫ / ‖v‖ ^ 2) • v)) =
      -(((⟪a, Complex.I * v⟫ / ‖v‖ ^ 3) / ‖v‖) • v) := by
  have hn : ‖v‖ ≠ 0 := norm_ne_zero_iff.mpr hv
  have e3 : (⟪a, Complex.I * v⟫ / ‖v‖ ^ 3) / ‖v‖ = ⟪a, Complex.I * v⟫ / (‖v‖ ^ 2) ^ 2 := by
    field_simp
  rw [e3, Complex.sq_norm, Complex.normSq_apply]
  have hs : v.re * v.re + v.im * v.im ≠ 0 := by
    rw [← Complex.normSq_apply]; exact (Complex.normSq_pos.mpr hv).ne'
  have hs2 : v.re ^ 2 + v.im ^ 2 ≠ 0 := by rw [sq, sq]; exact hs
  have hs3 : v.im ^ 2 + v.re ^ 2 ≠ 0 := by rw [add_comm]; exact hs2
  apply Complex.ext <;>
  · simp only [Complex.real_smul, Complex.inner, Complex.mul_re, Complex.mul_im, Complex.I_re,
      Complex.I_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.neg_re, Complex.neg_im, Complex.conj_re, Complex.conj_im, Complex.inv_re,
      Complex.inv_im, Complex.div_re, Complex.div_im]
    field_simp
    ring

section Curve
variable {X : ℝ → ℂ} (hX : ContDiff ℝ ∞ X) (hreg : ∀ θ, deriv X θ ≠ 0)
include hX

theorem contDiff_deriv' : ContDiff ℝ ∞ (deriv X) := (contDiff_infty_iff_deriv.mp hX).2

theorem contDiff_deriv2 : ContDiff ℝ ∞ (deriv (deriv X)) :=
  (contDiff_infty_iff_deriv.mp (contDiff_deriv' hX)).2

theorem fderiv_fderiv_real (w b c : ℝ) :
    fderiv ℝ (fderiv ℝ X) w b c = (b * c) • deriv (deriv X) w := by
  have hd : DifferentiableAt ℝ (fderiv ℝ X) w :=
    ((hX.fderiv_right (m := ∞) le_rfl).differentiable (by simp)) w
  have h1 : deriv (fderiv ℝ X) w 1 = deriv (deriv X) w := by
    have := deriv_clm_apply (u := fun _ : ℝ => (1 : ℝ)) hd (differentiableAt_const _)
    simp only [deriv_const, map_zero, add_zero] at this
    rw [← this]; rfl
  rw [fderiv_real_apply, ContinuousLinearMap.smul_apply, clm_real_apply (deriv (fderiv ℝ X) w) c,
    h1, smul_smul]

include hreg

theorem meanCurv_real (w : ℝ) :
    meanCurv X w = (‖deriv X w‖ ^ 2)⁻¹ •
      (deriv (deriv X) w - (⟪deriv (deriv X) w, deriv X w⟫ / ‖deriv X w‖ ^ 2) • deriv X w) := by
  have hL : fderiv ℝ X w 1 = deriv X w := rfl
  have hL0 : fderiv ℝ X w 1 ≠ 0 := hreg w
  rw [meanCurv, HOf, trOf_dim1 hL0]
  · rw [fderiv_fderiv_real hX, one_mul, one_smul, hL, perpOf, ginvOf_dim1 hL0,
      clm_real_apply (fderiv ℝ X w), hL]
    rfl
  · intro b c
    rw [fderiv_fderiv_real hX, fderiv_fderiv_real hX, one_mul, one_smul, perpOf_smul]

theorem JmeanCurv_real :
    (fun w => Jc (meanCurv X w)) = fun w => -((curv X w / ‖deriv X w‖) • deriv X w) := by
  funext w
  rw [Jc_apply, meanCurv_real hX hreg, I_mul_perp _ _ (hreg w)]
  rfl

theorem differentiable_curv : Differentiable ℝ (curv X) := by
  intro w
  have hv := (contDiff_deriv' hX).differentiable (by simp)
  have ha := (contDiff_deriv2 hX).differentiable (by simp)
  exact ((ha w).inner ℝ ((hv w).const_mul Complex.I)).div (((hv w).norm ℝ (hreg w)).pow 3)
    (pow_ne_zero 3 (norm_ne_zero_iff.mpr (hreg w)))

theorem flowPot_real (w : ℝ) : flowPot Jc X w = -(deriv (curv X) w / ‖deriv X w‖) := by
  have hL : fderiv ℝ X w 1 = deriv X w := rfl
  have hL0 : fderiv ℝ X w 1 ≠ 0 := hreg w
  set D := fderiv ℝ (fun z => Jc (meanCurv X z)) w with hDdef
  rw [flowPot, ← hDdef, trOf_dim1 hL0]
  swap
  · intro b c
    rw [clm_real_apply D b, clm_real_apply (fderiv ℝ X w) c, real_inner_smul_left,
      real_inner_smul_right, smul_eq_mul]
    ring
  -- derivative of `JH = -(k/n) • v`
  have hv := (contDiff_deriv' hX).differentiable (by simp)
  have hvd : HasDerivAt (deriv X) (deriv (deriv X) w) w := (hv w).hasDerivAt
  have hk : HasDerivAt (curv X) (deriv (curv X) w) w := (differentiable_curv hX hreg w).hasDerivAt
  have hnd : DifferentiableAt ℝ (fun y => ‖deriv X y‖) w := (hv w).norm ℝ (hreg w)
  have hn : HasDerivAt (fun y => ‖deriv X y‖) (deriv (fun y => ‖deriv X y‖) w) w := hnd.hasDerivAt
  have hn0 : ‖deriv X w‖ ≠ 0 := norm_ne_zero_iff.mpr (hreg w)
  have hnn : deriv (fun y => ‖deriv X y‖) w * ‖deriv X w‖ = ⟪deriv X w, deriv (deriv X) w⟫ := by
    have h1 := hvd.norm_sq
    have h2 := hn.pow 2
    have := h1.unique h2
    simp only [Nat.cast_ofNat, Nat.add_one_sub_one, pow_one] at this
    linarith [this]
  have hρ := hk.div hn hn0
  have hJ := (hρ.smul hvd).neg
  have hD1 : D 1 = -((curv X w / ‖deriv X w‖) • deriv (deriv X) w +
      ((deriv (curv X) w * ‖deriv X w‖ - curv X w * deriv (fun y => ‖deriv X y‖) w) /
        ‖deriv X w‖ ^ 2) • deriv X w) := by
    rw [hDdef, JmeanCurv_real hX hreg]
    exact hJ.deriv
  rw [hD1, hL, smul_eq_mul, inner_neg_left, inner_add_left, real_inner_smul_left,
    real_inner_smul_left, real_inner_self_eq_norm_sq, real_inner_comm (deriv X w), ← hnn]
  field_simp
  ring

theorem flowVel_real (θ : ℝ) :
    flowVel Jc X θ =
      -((dS X (dS X (curv X)) θ : ℝ) : ℂ) * (Complex.I * (deriv X θ / (‖deriv X θ‖ : ℂ))) := by
  have hL0 : fderiv ℝ X θ 1 ≠ 0 := hreg θ
  have hpot : flowPot Jc X = -(dS X (curv X)) := funext fun w => by
    rw [flowPot_real hX hreg]; rfl
  have hn0 : (‖deriv X θ‖ : ℂ) ≠ 0 := by exact_mod_cast norm_ne_zero_iff.mpr (hreg θ)
  rw [flowVel, hpot, ginvOf_dim1 hL0, clm_real_apply (fderiv ℝ X θ)]
  change Jc ((deriv (-(dS X (curv X))) θ / ‖deriv X θ‖ ^ 2) • deriv X θ) = _
  rw [deriv.neg, Jc_apply, Complex.real_smul]
  simp only [dS]
  push_cast
  field_simp

end Curve
end Dim1

/-- **(iii) n = 1: `flowVel = −k_ss N`, `N = i T`.**  NOT PROVED (needs a jet computation of
`HOf`, `ginvOf`, `trOf` at `E = ℝ`, i.e. `trOf L B = B 1 1 / ‖L 1‖²`, plus three `fderiv`
chain rules); the identity itself is checked exactly in U1_sym.py (S1). -/
theorem flowVel_eq_neg_kss {X : ℝ → ℂ} (hX : ContDiff ℝ ∞ X) (hreg : ∀ θ, deriv X θ ≠ 0) (θ : ℝ) :
    flowVel Jc X θ =
      -((dS X (dS X (curv X)) θ : ℝ) : ℂ) * (Complex.I * (deriv X θ / (‖deriv X θ‖ : ℂ))) :=
  -- (S1) closed 2026-11-23 (statement unchanged): `flowVel_real` above
  flowVel_real hX hreg θ

/-! ## (iv) Symmetry under reflection, via uniqueness

`σ = conj` on `ℂ` (an `ℝ`-linear isometry with `σ J = −J σ`), `ι p = p⁻¹ = conj p` on `Circle`.
Neither the closure of `IsFlow` under `σ ∘ · ∘ ι` nor the pointwise identity `flowVel_conj`
is in the repo (U1_regularity.tex already says "it is not a Lean lemma"). -/

section Conj
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
set_option linter.unusedSectionVars false

/-- `conj` as a continuous ℝ-linear map. -/
abbrev cL : ℂ →L[ℝ] ℂ := (Complex.conjLIE : ℂ →L[ℝ] ℂ)

@[simp] theorem cL_apply (w : ℂ) : cL w = conj w := rfl

theorem inner_conj_conj (u v : ℂ) : ⟪conj u, conj v⟫ = ⟪u, v⟫ := by
  have := Complex.conjLIE.inner_map_map u v
  simpa using this

theorem fderiv_conj_comp {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (X : F → ℂ) (y : F) :
    fderiv ℝ (fun z => conj (X z)) y = cL.comp (fderiv ℝ X y) := by
  have h := Complex.conjLIE.comp_fderiv (f := X) (x := y)
  exact h

theorem gOf_conj (L : E →L[ℝ] ℂ) : gOf (cL.comp L) = gOf L := by
  ext a b
  simp

theorem ginvOf_conj (L : E →L[ℝ] ℂ) : ginvOf (cL.comp L) = ginvOf L := by
  simp only [ginvOf, gOf_conj]

theorem trOf_conj {W : Type*} [AddCommMonoid W] (L : E →L[ℝ] ℂ) (B : E → E → W) :
    trOf (cL.comp L) B = trOf L B := by
  simp only [trOf, ginvOf_conj]

theorem perpOf_conj (L : E →L[ℝ] ℂ) (v : ℂ) : perpOf (cL.comp L) (conj v) = conj (perpOf L v) := by
  have h : ((innerSL ℝ : ℂ →L[ℝ] ℂ →L[ℝ] ℝ) (conj v)).comp (cL.comp L) =
      ((innerSL ℝ : ℂ →L[ℝ] ℂ →L[ℝ] ℝ) v).comp L := by
    ext a
    simp
  simp only [perpOf, h, ginvOf_conj, ContinuousLinearMap.comp_apply, cL_apply, map_sub]

theorem meanCurv_conj (X : E → ℂ) (y : E) :
    meanCurv (fun z => conj (X z)) y = conj (meanCurv X y) := by
  have h1 : fderiv ℝ (fun z => conj (X z)) = fun z => cL.comp (fderiv ℝ X z) :=
    funext (fderiv_conj_comp X)
  let A : (E →L[ℝ] ℂ) ≃L[ℝ] (E →L[ℝ] ℂ) :=
    ContinuousLinearEquiv.arrowCongr (ContinuousLinearEquiv.refl ℝ E)
      (Complex.conjLIE.toContinuousLinearEquiv)
  have hA : ∀ L, A L = cL.comp L := fun L => by ext a; simp [A]
  have h2 : fderiv ℝ (fderiv ℝ (fun z => conj (X z))) y =
      (A : (E →L[ℝ] ℂ) →L[ℝ] (E →L[ℝ] ℂ)).comp (fderiv ℝ (fderiv ℝ X) y) := by
    have e : fderiv ℝ (fun z => conj (X z)) = A ∘ fderiv ℝ X := by
      rw [h1]; funext z; simp [hA]
    rw [e]
    exact A.comp_fderiv
  have hQ : ∀ a c, fderiv ℝ (fderiv ℝ (fun z => conj (X z))) y a c =
      conj (fderiv ℝ (fderiv ℝ X) y a c) := fun a c => by
    rw [h2]; simp [hA]
  simp only [meanCurv, HOf, hQ]
  rw [fderiv_conj_comp X y]
  simp only [trOf_conj, perpOf_conj]
  simp only [trOf, map_sum]

theorem flowPot_conj (X : E → ℂ) :
    flowPot Jc (fun z => conj (X z)) = fun y => -flowPot Jc X y := by
  funext y
  have hJ : (fun z => Jc (meanCurv (fun z => conj (X z)) z)) =
      fun z => -(cL (Jc (meanCurv X z))) := by
    funext z
    rw [meanCurv_conj]
    simp only [Jc_apply, cL_apply, map_mul, Complex.conj_I]
    ring
  have hD : fderiv ℝ (fun z => Jc (meanCurv (fun z => conj (X z)) z)) y =
      -(cL.comp (fderiv ℝ (fun z => Jc (meanCurv X z)) y)) := by
    rw [hJ, fderiv_fun_neg]
    congr 1
    exact fderiv_conj_comp (fun z => Jc (meanCurv X z)) y
  simp only [flowPot, hD, fderiv_conj_comp X y, trOf_conj]
  simp only [trOf, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp

theorem flowVel_conj_any (X : E → ℂ) (y : E) :
    flowVel Jc (fun z => conj (X z)) y = conj (flowVel Jc X y) := by
  simp only [flowVel, flowPot_conj, fderiv_fun_neg, fderiv_conj_comp X y, ginvOf_conj, map_neg,
    ContinuousLinearMap.comp_apply, cL_apply, Jc_apply, map_mul, Complex.conj_I]
  ring
end Conj

/-- `flowVel (σ ∘ X) = σ (flowVel X)` for `σ = conj` (J appears twice, so the sign cancels:
`H' = σH`, `JH' = −σJH`, `flowPot' = −flowPot`, `flowVel' = J σ L g⁻¹ D(−pot) = σ flowVel`).
NOT PROVED; U1_sym.py S2 checks it at n = 1. -/
theorem flowVel_conj {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
    {X : E → ℂ} {U : Set E} (hU : IsOpen U) (hX : ContDiffOn ℝ ∞ X U)
    (hinj : ∀ y ∈ U, Function.Injective (fderiv ℝ X y)) {y : E} (hy : y ∈ U) :
    flowVel Jc (fun z => conj (X z)) y = conj (flowVel Jc X y) :=
  -- (S2) closed: the identity holds with no hypotheses at all (`flowVel_conj_any`)
  flowVel_conj_any X y

section FlowOps
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
set_option linter.unusedSectionVars false

/-- `conj ∘ F` is an (automatically Lagrangian) immersion when `F` is. -/
theorem isLagImmersion_conj {F : M → ℂ} (hF : IsLagImmersion E Jc F) :
    IsLagImmersion E Jc (fun p => conj (F p)) where
  smooth z := cL.contDiff.comp_contDiffOn (hF.smooth z)
  inj z y hy := by
    have e : fderiv ℝ (chartRep E (fun p => conj (F p)) z) y = cL.comp (fderiv ℝ (chartRep E F z) y) :=
      fderiv_conj_comp (chartRep E F z) y
    rw [e]
    exact Complex.conjLIE.injective.comp (hF.inj z y hy)
  lag z y hy u w := by
    have e : fderiv ℝ (chartRep E (fun p => conj (F p)) z) y = cL.comp (fderiv ℝ (chartRep E F z) y) :=
      fderiv_conj_comp (chartRep E F z) y
    rw [e]
    have h := hF.lag z y hy u w
    simp only [ContinuousLinearMap.comp_apply, cL_apply, Jc_apply] at h ⊢
    have : Complex.I * conj (fderiv ℝ (chartRep E F z) y u) =
        -conj (Complex.I * fderiv ℝ (chartRep E F z) y u) := by
      simp only [map_mul, Complex.conj_I]; ring
    rw [this, inner_neg_left, inner_conj_conj, h, neg_zero]

/-- A flow followed by `conj` is a flow (for `J = i·`, `conj J = −J conj`). -/
theorem IsFlow.conj {T : ℝ} {Φ : ℝ → M → ℂ} (h : IsFlow E Jc T Φ) :
    IsFlow E Jc T (fun t p => conj (Φ t p)) where
  smooth := cL.contDiff.contMDiff.comp_contMDiffOn h.smooth
  imm t ht := isLagImmersion_conj (h.imm t ht)
  vel t ht x z hx := by
    have h1 := (cL.hasFDerivAt (x := Φ t x)).comp_hasDerivWithinAt t (h.vel t ht x z hx)
    have e : chartRep E (fun p => conj (Φ t p)) z = fun y => conj (chartRep E (Φ t) z y) := rfl
    rw [e, flowVel_conj_any]
    exact h1

variable [IsManifold 𝓘(ℝ, E) ∞ M]

/-- A flow precomposed with a fixed diffeomorphism of `M` is a flow. -/
theorem IsFlow.comp_diffeo {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ)
    (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) : IsFlow E J T (fun t => Φ t ∘ φ) where
  smooth := by
    have hg : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((q.1, φ q.2) : ℝ × M)) :=
      contMDiff_fst.prodMk (φ.contMDiff.comp contMDiff_snd)
    exact h.smooth.comp hg.contMDiffOn (fun q hq => ⟨hq.1, trivial⟩)
  imm t ht := isLagImmersion_comp_diffeo (h.imm t ht) φ
  vel t ht x z hx := by
    set α := φ x
    have hα : φ x ∈ (chart E α).source := mem_extChartAt_source _
    have h1 := h.vel t ht (φ x) α hα
    have hy : chart E z x ∈ diffeoDom E φ z α := mem_diffeoDom hx hα
    have hN := isOpen_diffeoDom (E := E) φ z α
    -- locality: the two chart readings agree near `chart z x`
    have hev : chartRep E (Φ t ∘ φ) z =ᶠ[𝓝 (chart E z x)]
        chartRep E (Φ t) α ∘ diffeoRep E φ z α := by
      filter_upwards [hN.mem_nhds hy] with y hy'
      exact chartRep_comp_diffeo (Φ t) hy'
    have himm := h.imm t ht
    have hmapsT : ∀ y ∈ diffeoDom E φ z α, diffeoRep E φ z α y ∈ (chart E α).target := by
      intro y hy'
      simp only [diffeoRep, Function.comp_apply]
      exact (chart E α).map_source hy'.2
    have hX : ∀ y ∈ diffeoDom E φ z α, ContDiffAt ℝ ∞ (chartRep E (Φ t) α) (diffeoRep E φ z α y) :=
      fun y hy' => (himm.smooth α).contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds (hmapsT y hy'))
    have hinj : ∀ y ∈ diffeoDom E φ z α,
        Function.Injective (fderiv ℝ (chartRep E (Φ t) α) (diffeoRep E φ z α y)) :=
      fun y hy' => himm.inj α _ (hmapsT y hy')
    have hpot : DifferentiableAt ℝ (flowPot J (chartRep E (Φ t) α)) (diffeoRep E φ z α (chart E z x)) := by
      have hc := contDiffOn_flowPot (isOpen_chart_target (E := E) α) (himm.smooth α)
        (fun y hy' => himm.inj α y hy') J
      exact (hc.contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds (hmapsT _ hy))).differentiableAt
        (by simp)
    have hvel : flowVel J (chartRep E (Φ t ∘ φ) z) (chart E z x) =
        flowVel J (chartRep E (Φ t) α) (chart E α (φ x)) := by
      rw [flowVel_congr_of_eventuallyEq hev,
        flowVel_comp hN (contDiffOn_diffeoRep φ z α) (fun y hy' => exists_equiv_fderiv_diffeoRep φ hy')
          hX hinj J hy hpot, diffeoRep_chart hx]
    rw [hvel]
    exact h1

end FlowOps

/-- The reflection `p ↦ p⁻¹ = conj p` of `Circle`, as a diffeomorphism. -/
def invDiffeo : Circle ≃ₘ⟮𝓘(ℝ, E1), 𝓘(ℝ, E1)⟯ Circle where
  toEquiv := Equiv.inv Circle
  contMDiff_toFun := (contMDiff_inv (I := 𝓘(ℝ, E1)) (n := ω) (G := Circle)).of_le le_top
  contMDiff_invFun := (contMDiff_inv (I := 𝓘(ℝ, E1)) (n := ω) (G := Circle)).of_le le_top

@[simp] theorem invDiffeo_apply (p : Circle) : invDiffeo p = p⁻¹ := rfl

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
    IsFlow E1 Jc T (fun t p => conj (Φ t p⁻¹)) :=
  -- (S3) closed 2026-11-23 (statement unchanged)
  IsFlow.conj (IsFlow.comp_diffeo h invDiffeo)

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

/-- Time shift for a general flow (any `E`, `M`, `V`, `J`); only `0 ≤ a` is needed. -/
theorem IsFlow.shift_gen {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V} {T a : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) (ha : 0 ≤ a) :
    IsFlow E J (T - a) (fun t => Φ (t + a)) where
  smooth := by
    have h1 : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => q.1 + a) :=
      ((contDiff_id.add contDiff_const : ContDiff ℝ ∞ (fun x : ℝ => x + a)).contMDiff).comp
        contMDiff_fst
    have hg : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((q.1 + a, q.2) : ℝ × M)) := h1.prodMk contMDiff_snd
    exact h.smooth.comp hg.contMDiffOn
      (fun q hq => ⟨⟨by linarith [hq.1.1], by linarith [hq.1.2]⟩, trivial⟩)
  imm := fun t ht => h.imm (t + a) ⟨by linarith [ht.1], by linarith [ht.2]⟩
  vel := fun t ht x z hx => by
    have h1 := h.vel (t + a) ⟨by linarith [ht.1], by linarith [ht.2]⟩ x z hx
    have hs : HasDerivWithinAt (fun s : ℝ => s + a) 1 (Icc 0 (T - a)) t :=
      ((hasDerivAt_id t).add_const a).hasDerivWithinAt
    have hmap : MapsTo (fun s : ℝ => s + a) (Icc 0 (T - a)) (Icc 0 T) :=
      fun s hs' => ⟨by linarith [hs'.1], by linarith [hs'.2]⟩
    have h2 := h1.scomp t hs hmap
    rw [one_smul] at h2
    exact h2

/-- Time shift (U1 "time shift", C1 Thm (v)).  NOT PROVED; expected easy
(`ContMDiffOn.comp` with `(t, x) ↦ (t + a, x)`; [GUESS] `HasDerivWithinAt.comp_add_const`). -/
theorem IsFlow.shift {T a : ℝ} {Φ : ℝ → Circle → ℂ} (h : IsFlow E1 Jc T Φ) (ha : 0 ≤ a)
    (haT : a ≤ T) : IsFlow E1 Jc (T - a) (fun t => Φ (t + a)) :=
  -- (S4) closed; the hypothesis `haT` is not needed (`IsFlow.shift_gen`)
  IsFlow.shift_gen h ha

/-- Gluing (U1 Lemma "gluing").  NOT PROVED and NOT cheap: joint C^∞ at the seam `t = a` needs
all one-sided `∂ₜ^m` to agree, which U1 gets from the equation by induction on `m`. -/
theorem IsFlow.glue {a b : ℝ} {Φ₁ Φ₂ : ℝ → Circle → ℂ} (h₁ : IsFlow E1 Jc a Φ₁)
    (h₂ : IsFlow E1 Jc b Φ₂) (ha : 0 ≤ a) (hb : 0 ≤ b) (h : Φ₂ 0 = Φ₁ a) :
    IsFlow E1 Jc (a + b) (fun t => if t ≤ a then Φ₁ t else Φ₂ (t - a)) := by
  sorry -- SORRY (S5): gluing

/-! ## Axiom audit (expected after a successful build)
`cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam` should depend only on
`[propext, Classical.choice, Quot.sound]`.  Since v4 (2026-11-23) the same holds for `cdf_symmetric`
(S2, S3 closed), `flowVel_eq_neg_kss` (S1), `IsFlow.shift` (S4); only `IsFlow.glue` (S5) uses `sorryAx`. -/
#print axioms cdf_exists
#print axioms cdf_unique
#print axioms cdf_unique_up_to_reparam
#print axioms cdf_symmetric
#print axioms flowVel_eq_neg_kss
#print axioms flowVel_conj
#print axioms IsFlow.conj_reflect
#print axioms IsFlow.shift
#print axioms IsFlow.glue

end CWFlowN1
