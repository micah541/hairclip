/-
# The `n = 1` instance of `CWFlow`: curve diffusion flow (port of N1.lean)

Port of `CWFlowN1_v5.lean` (compiled in the older copy of the repository, commit d72b744) to this repository,
as a new file. Changes against v5: namespace `CWFlow.Curve`; `curv`, `dS` moved to `CWFlow.Curve.Defs`;
the gluing section (`IsFlow.glue_of_seam`, the `sorry` `glue_seam`, `IsFlow.glue`) dropped
(`CWFlow.rst_glue` in `Flow/Restart.lean` supersedes it); `cdf_exists_maximal_A` added (instantiation of
`exists_maximal_flow` at `n = 1`).

Setting: `E1 = EuclideanSpace ℝ (Fin 1)`, `M = Circle`, `V = ℂ`, `J = (I * ·)`.
Contents: (i) `cdf_exists`, (ii) `cdf_unique`, (ii') `cdf_unique_up_to_reparam`,
(iii) `flowVel_eq_neg_kss` (`J∇div(JH) = −κ_ss N` for curves), (iv) `cdf_symmetric`, (v) `IsFlow.shift_gen`,
and `cdf_exists_maximal_A`.
-/
import CWFlow.Flow.Main               -- exists_flow_of_immersion
import CWFlow.Flow.UniqueMain         -- flow_unique_of_immersion
import CWFlow.Flow.NormalFlow         -- exists_flow_of_normal_family (imports GraphUnique: IsFlow)
import CWFlow.Flow.GaugeEvolution     -- IsFlow.mono
import CWFlow.Scalar.Reparam          -- flowVel_comp (used in S3)
import CWFlow.Scalar.Smooth           -- contDiffOn_flowPot (used in S3)
import CWFlow.Flow.DiffeoChart        -- diffeoRep, isLagImmersion_comp_diffeo (used in S3)
import CWFlow.Flow.TanAlgebra         -- flowVel_congr_of_eventuallyEq (used in S3)
import Mathlib.Geometry.Manifold.Instances.Sphere   -- charts on `Circle`
import CWFlow.Curve.N1Basic
import CWFlow.Curve.CircleChart
import CWFlow.Flow.MaximalFlow       -- exists_maximal_flow

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold ComplexConjugate

noncomputable section

namespace CWFlow.Curve

open CWFlow

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
theorem IsFlow.comp_diffeo_n1 {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
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
  IsFlow.conj (IsFlow.comp_diffeo_n1 h invDiffeo)

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

/-! ## The maximal flow at `n = 1` (instantiation of `exists_maximal_flow`) -/

-- `ι1 : ℝ¹ ≃ E1` (the identity) is defined in `CWFlow.Curve.CircleChart`.

/-- **The maximal curve diffusion flow**, with blow-up of `|A|²` (`exists_maximal_flow` at `n = 1`). -/
theorem cdf_exists_maximal_A {F₀ : Circle → ℂ} (hF₀ : IsLagImmersion E1 Jc F₀) :
    ∃ Tmax : ENNReal, 0 < Tmax ∧ ∃ Φ : ℝ → Circle → ℂ, Φ 0 = F₀ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E1 Jc T Φ) ∧
      (∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = F₀ →
        ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t) ∧
      (Tmax ≠ ⊤ → ∀ K : ℝ, ¬ ∀ t ∈ Ico (0 : ℝ) Tmax.toReal, ∀ z : Circle,
        ∀ y ∈ Extension.chartDomE ι1 z, Extension.aNormSq Jc (Φ t ∘ Extension.chartInv ι1 z) y ≤ K ^ 2) :=
  exists_maximal_flow ι1 Jc_isom Jc_JJ Complex.finrank_real_complex hF₀

/-! ## Axiom audit -/

/-- info: 'CWFlow.Curve.cdf_exists' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_exists
/-- info: 'CWFlow.Curve.cdf_unique' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_unique
/-- info: 'CWFlow.Curve.cdf_unique_up_to_reparam' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_unique_up_to_reparam
/-- info: 'CWFlow.Curve.cdf_symmetric' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_symmetric
/-- info: 'CWFlow.Curve.flowVel_eq_neg_kss' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowVel_eq_neg_kss
/-- info: 'CWFlow.Curve.IsFlow.shift_gen' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsFlow.shift_gen
/-- info: 'CWFlow.Curve.cdf_exists_maximal_A' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_exists_maximal_A

end CWFlow.Curve
