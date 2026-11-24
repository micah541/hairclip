import CWFlow.Flow.LocalInverse
import DifferentialGeometry.Analysis.Calculus.SmoothExtension.BorelHalfLine.Parametric
import Mathlib.LinearAlgebra.Dual.Lemmas

/-!
# The graph gauge, locally in a chart

Let `Ψ : E × E* → V` be smooth with injective derivative on an open set `Q` (a Weinstein map in cotangent
coordinates) and let `X : ℝ × E → V` be smooth on `[0, T] × U` with `X(0, z) = Ψ(z, 0)` (a family of immersions
starting at the zero section, in a chart).

`exists_local_gauge`: near every `y₀ ∈ U` and for short time there are smooth `ζ(t, y) ∈ E`, `ξ(t, y) ∈ E*` with
`X(t, ζ(t, y)) = Ψ(y, ξ(t, y))`, `ζ(0, y) = y`, `ξ(0, y) = 0`: after the reparametrization `ζ_t` the family is a
section in cotangent coordinates. The solution is unique among pairs `(z, ξ')` with `(t, z) ∈ O`, `(y, ξ') ∈ A`
for open sets `O`, `A` that contain the solution; this is what patches the local gauges together.

Proof: invert `Ψ` near `(y₀, 0)`, `Θ = Ψ⁻¹`; then `χ(t, z) = Θ₁(X(t, z))` has `χ(0, z) = z`, so
`(t, z) ↦ (t, χ(t, z))` is invertible near `(0, y₀)` (after extending `X` to negative times), and
`ζ_t = χ_t⁻¹`, `ξ = Θ₂(X(t, ζ))`.
-/

open Set Filter

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] [FiniteDimensional ℝ V]

theorem finrank_dualCLM : Module.finrank ℝ (E →L[ℝ] ℝ) = Module.finrank ℝ E := by
  rw [← (LinearMap.toContinuousLinearMap (𝕜 := ℝ) (E := E) (F' := ℝ)).finrank_eq]
  exact Subspace.dual_finrank_eq

theorem exists_local_gauge (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {Ψ : E × (E →L[ℝ] ℝ) → V} {Q : Set (E × (E →L[ℝ] ℝ))} (hQ : IsOpen Q)
    (hΨ : ContDiffOn ℝ ∞ Ψ Q) (hΨinj : ∀ q ∈ Q, Function.Injective (fderiv ℝ Ψ q))
    {T : ℝ} (hT : 0 < T) {U : Set E} (hU : IsOpen U) {X : ℝ × E → V}
    (hX : ContDiffOn ℝ ∞ X (Icc 0 T ×ˢ U)) (hX0 : ∀ z ∈ U, X (0, z) = Ψ (z, 0))
    (hU0 : ∀ z ∈ U, (z, (0 : E →L[ℝ] ℝ)) ∈ Q) {y₀ : E} (hy₀ : y₀ ∈ U) :
    ∃ ε, 0 < ε ∧ ε ≤ T ∧ ∃ B : Set E, IsOpen B ∧ y₀ ∈ B ∧ B ⊆ U ∧
      ∃ (ζ : ℝ × E → E) (ξ : ℝ × E → (E →L[ℝ] ℝ)) (O : Set (ℝ × E)) (A : Set (E × (E →L[ℝ] ℝ))),
        IsOpen O ∧ IsOpen A ∧ A ⊆ Q ∧
        ContDiffOn ℝ ∞ ζ (Icc 0 ε ×ˢ B) ∧ ContDiffOn ℝ ∞ ξ (Icc 0 ε ×ˢ B) ∧
        (∀ y ∈ B, ζ (0, y) = y ∧ ξ (0, y) = 0) ∧
        (∀ t ∈ Icc (0 : ℝ) ε, ∀ y ∈ B, ζ (t, y) ∈ U ∧ (t, ζ (t, y)) ∈ O ∧ (y, ξ (t, y)) ∈ A ∧
          X (t, ζ (t, y)) = Ψ (y, ξ (t, y))) ∧
        (∀ t ∈ Icc (0 : ℝ) ε, ∀ y ∈ B, ∀ (z : E) (ξ' : E →L[ℝ] ℝ), (t, z) ∈ O → (y, ξ') ∈ A →
          X (t, z) = Ψ (y, ξ') → z = ζ (t, y) ∧ ξ' = ξ (t, y)) := by
  classical
  have hVc : CompleteSpace V := FiniteDimensional.complete ℝ V
  have hdimΨ : Module.finrank ℝ (E × (E →L[ℝ] ℝ)) = Module.finrank ℝ V := by
    rw [Module.finrank_prod, finrank_dualCLM, hdim]
    ring
  -- the local inverse of `Ψ`
  obtain ⟨e, he0, heQ, hecoe, hesm⟩ :=
    exists_openPartialHomeomorph_of_injective_fderiv hdimΨ hQ hΨ hΨinj (hU0 y₀ hy₀)
  -- extension of `X` to all times
  have hunc : Function.uncurry (fun (t : ℝ) (z : E) => X (t, z)) = X := funext fun p => rfl
  obtain ⟨gext, V0, hV0, hgsm, hgeq⟩ := DifferentialGeometry.Analysis.borel_interval_extend_param
    (fun (t : ℝ) (z : E) => X (t, z)) T hT U y₀ (by rw [hU.interior_eq]; exact hy₀)
    (by rw [hunc]; exact hX)
  have hW0open : IsOpen (interior V0 ∩ U) := isOpen_interior.inter hU
  have hy₀W : y₀ ∈ interior V0 ∩ U := ⟨mem_interior_iff_mem_nhds.mpr hV0, hy₀⟩
  have hXt : ContDiffOn ℝ ∞ (Function.uncurry gext) ((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) :=
    hgsm.mono (prod_mono_right (inter_subset_left.trans interior_subset))
  have hXteq : ∀ t ∈ Icc (0 : ℝ) T, ∀ z ∈ interior V0 ∩ U, Function.uncurry gext (t, z) = X (t, z) :=
    fun t ht z hz => hgeq t ht z (interior_subset hz.1)
  have hDopen : IsOpen (((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target) :=
    hXt.continuousOn.isOpen_inter_preimage (isOpen_univ.prod hW0open) e.open_target
  have hκ : ContDiffOn ℝ ∞ (fun p : ℝ × E => (e.symm (Function.uncurry gext p)).1)
      (((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target) :=
    contDiff_fst.comp_contDiffOn (hesm.comp (hXt.mono inter_subset_left) fun p hp => hp.2)
  have hK : ContDiffOn ℝ ∞ (fun p : ℝ × E => (p.1, (e.symm (Function.uncurry gext p)).1))
      (((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target) :=
    contDiffOn_fst.prodMk hκ
  -- at time zero
  have h0 : ∀ z ∈ interior V0 ∩ U, (z, (0 : E →L[ℝ] ℝ)) ∈ e.source →
      (0, z) ∈ ((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target ∧
        e.symm (Function.uncurry gext (0, z)) = (z, 0) := by
    intro z hz hzs
    have h1 : Function.uncurry gext (0, z) = e (z, 0) := by
      rw [hXteq 0 (left_mem_Icc.mpr hT.le) z hz, hX0 z hz.2, hecoe]
    refine ⟨⟨⟨mem_univ _, hz⟩, ?_⟩, ?_⟩
    · rw [mem_preimage, h1]
      exact e.map_source hzs
    · rw [h1]
      exact e.left_inv hzs
  have hsrc0 : IsOpen {z : E | (z, (0 : E →L[ℝ] ℝ)) ∈ e.source} :=
    e.open_source.preimage (continuous_id.prodMk continuous_const)
  have hp0 := (h0 y₀ hy₀W he0).1
  -- the derivative of `K` at `(0, y₀)` is injective
  have hκd : DifferentiableAt ℝ (fun p : ℝ × E => (e.symm (Function.uncurry gext p)).1) (0, y₀) :=
    (hκ.contDiffAt (hDopen.mem_nhds hp0)).differentiableAt (by simp)
  have hκb : ∀ b : E,
      fderiv ℝ (fun p : ℝ × E => (e.symm (Function.uncurry gext p)).1) (0, y₀) (0, b) = b := by
    intro b
    have hι : HasFDerivAt (fun z : E => ((0 : ℝ), z))
        ((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E)) y₀ :=
      (hasFDerivAt_const (0 : ℝ) y₀).prodMk (hasFDerivAt_id y₀)
    have hcomp := HasFDerivAt.comp (g := fun p : ℝ × E => (e.symm (Function.uncurry gext p)).1)
      (f := fun z : E => ((0 : ℝ), z)) y₀ hκd.hasFDerivAt hι
    have hid : ((fun p : ℝ × E => (e.symm (Function.uncurry gext p)).1) ∘ fun z : E => ((0 : ℝ), z))
        =ᶠ[𝓝 y₀] id := by
      filter_upwards [hW0open.mem_nhds hy₀W, hsrc0.mem_nhds he0] with z hz hzs
      simp only [Function.comp_apply, id]
      rw [(h0 z hz hzs).2]
    have hc := (hcomp.congr_of_eventuallyEq hid.symm).unique (hasFDerivAt_id y₀)
    have := congrArg (fun L : E →L[ℝ] E => L b) hc
    simpa using this
  have hKinj : Function.Injective
      (fderiv ℝ (fun p : ℝ × E => (p.1, (e.symm (Function.uncurry gext p)).1)) (0, y₀)) := by
    rw [DifferentiableAt.fderiv_prodMk differentiableAt_fst hκd, injective_iff_map_eq_zero]
    rintro ⟨τ, b⟩ hv
    have h1 : τ = 0 := by
      have := congrArg Prod.fst hv
      simpa [fderiv_fst] using this
    have h2 := congrArg Prod.snd hv
    simp only [ContinuousLinearMap.prod_apply, Prod.snd_zero] at h2
    rw [h1, hκb b] at h2
    rw [h1, h2]
    rfl
  -- the open set where the derivative of `K` is injective
  have hQk : IsOpen ((((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target) ∩
      (fderiv ℝ (fun p : ℝ × E => (p.1, (e.symm (Function.uncurry gext p)).1))) ⁻¹'
        {L : ℝ × E →L[ℝ] ℝ × E | Function.Injective L}) :=
    (hK.continuousOn_fderiv_of_isOpen hDopen (by simp)).isOpen_inter_preimage hDopen
      ContinuousLinearMap.isOpen_injective
  obtain ⟨k, hk0, hkQ, hkcoe, hksm⟩ := exists_openPartialHomeomorph_of_injective_fderiv rfl hQk
    (hK.mono inter_subset_left) (fun p hp => hp.2) (q₀ := ((0 : ℝ), y₀)) ⟨hp0, hKinj⟩
  -- facts about points of the target of `k`
  have hfacts : ∀ p ∈ k.target,
      (p.1, (k.symm p).2) ∈ k.source ∧ (e.symm (Function.uncurry gext (p.1, (k.symm p).2))).1 = p.2 := by
    intro p hp
    have hr : k.symm p ∈ k.source := k.map_target hp
    have hkr : k (k.symm p) = p := k.right_inv hp
    rw [hkcoe] at hkr
    have h1 : (k.symm p).1 = p.1 := by
      have h := congrArg Prod.fst hkr
      exact h
    have hrr : k.symm p = (p.1, (k.symm p).2) := Prod.ext h1 rfl
    refine ⟨hrr ▸ hr, ?_⟩
    have h2 := congrArg Prod.snd hkr
    rw [← hrr]
    exact h2
  have hk0' : k ((0 : ℝ), y₀) = (0, y₀) := by
    rw [hkcoe]
    exact Prod.ext rfl (congrArg Prod.fst (h0 y₀ hy₀W he0).2)
  have hk0t : ((0 : ℝ), y₀) ∈ k.target := hk0' ▸ k.map_source hk0
  obtain ⟨u, v, hu, hv, h0u, hy₀v, huv⟩ := isOpen_prod_iff.mp k.open_target 0 y₀ hk0t
  obtain ⟨ε', hε', hεu⟩ := Metric.isOpen_iff.mp hu 0 h0u
  have hsrck : IsOpen {y : E | ((0 : ℝ), y) ∈ k.source} :=
    k.open_source.preimage (continuous_const.prodMk continuous_id)
  have hε : 0 < min (ε' / 2) T := lt_min (by positivity) hT
  have hsub : Icc (0 : ℝ) (min (ε' / 2) T) ×ˢ
      (v ∩ (interior V0 ∩ U) ∩ {z : E | (z, (0 : E →L[ℝ] ℝ)) ∈ e.source} ∩
        {y : E | ((0 : ℝ), y) ∈ k.source}) ⊆ k.target := by
    rintro ⟨t, y⟩ ⟨ht, hy⟩
    refine huv ⟨hεu ?_, hy.1.1.1⟩
    rw [Real.ball_eq_Ioo]
    have h1 : t ≤ ε' / 2 := ht.2.trans (min_le_left _ _)
    constructor <;> linarith [ht.1]
  have hD : ∀ p ∈ k.target, (p.1, (k.symm p).2) ∈
      ((univ : Set ℝ) ×ˢ (interior V0 ∩ U)) ∩ Function.uncurry gext ⁻¹' e.target :=
    fun p hp => (hkQ (hfacts p hp).1).1
  have hζ : ContDiffOn ℝ ∞ (fun p : ℝ × E => (k.symm p).2) k.target :=
    contDiff_snd.comp_contDiffOn hksm
  have hpair : ContDiffOn ℝ ∞ (fun p : ℝ × E => (p.1, (k.symm p).2)) k.target :=
    contDiffOn_fst.prodMk hζ
  have hξ : ContDiffOn ℝ ∞
      (fun p : ℝ × E => (e.symm (Function.uncurry gext (p.1, (k.symm p).2))).2) k.target :=
    contDiff_snd.comp_contDiffOn
      (hesm.comp ((hXt.mono inter_subset_left).comp hpair fun p hp => (hD p hp))
        fun p hp => (hD p hp).2)
  refine ⟨min (ε' / 2) T, hε, min_le_right _ _,
    v ∩ (interior V0 ∩ U) ∩ {z : E | (z, (0 : E →L[ℝ] ℝ)) ∈ e.source} ∩
      {y : E | ((0 : ℝ), y) ∈ k.source},
    ((hv.inter hW0open).inter hsrc0).inter hsrck, ⟨⟨⟨hy₀v, hy₀W⟩, he0⟩, hk0⟩,
    fun y hy => hy.1.1.2.2,
    fun p => (k.symm p).2, fun p => (e.symm (Function.uncurry gext (p.1, (k.symm p).2))).2,
    k.source, e.source, k.open_source, e.open_source, heQ, hζ.mono hsub, hξ.mono hsub, ?_, ?_, ?_⟩
  · -- time zero
    intro y hy
    have hz := h0 y hy.1.1.2 hy.1.2
    have hky : k ((0 : ℝ), y) = (0, y) := by
      rw [hkcoe]
      exact Prod.ext rfl (congrArg Prod.fst hz.2)
    have hsy : k.symm ((0 : ℝ), y) = (0, y) := by
      have := k.left_inv hy.2
      rwa [hky] at this
    refine ⟨by simp only [hsy], ?_⟩
    simp only [hsy]
    exact congrArg Prod.snd hz.2
  · -- the equation
    intro t ht y hy
    have hp : (t, y) ∈ k.target := hsub ⟨ht, hy⟩
    have hf := hfacts (t, y) hp
    have hd := hD (t, y) hp
    have htT : t ∈ Icc (0 : ℝ) T := ⟨ht.1, ht.2.trans (min_le_right _ _)⟩
    have hsymm : e.symm (Function.uncurry gext (t, (k.symm (t, y)).2)) =
        (y, (e.symm (Function.uncurry gext (t, (k.symm (t, y)).2))).2) := Prod.ext hf.2 rfl
    have hmem : (y, (e.symm (Function.uncurry gext (t, (k.symm (t, y)).2))).2) ∈ e.source :=
      hsymm ▸ e.map_target hd.2
    refine ⟨hd.1.2.2, hf.1, hmem, ?_⟩
    have hr := e.right_inv hd.2
    rw [hsymm, hecoe] at hr
    exact (hXteq t htT _ hd.1.2).symm.trans hr.symm
  · -- uniqueness
    intro t ht y hy z ξ' hzO hA hEq
    have htT : t ∈ Icc (0 : ℝ) T := ⟨ht.1, ht.2.trans (min_le_right _ _)⟩
    have hzD := (hkQ hzO).1
    have h1 : e.symm (Function.uncurry gext (t, z)) = (y, ξ') := by
      rw [hXteq t htT z hzD.1.2, hEq, ← hecoe]
      exact e.left_inv hA
    have hkz : k (t, z) = (t, y) := by
      rw [hkcoe]
      exact Prod.ext rfl (congrArg Prod.fst h1)
    have hsz : k.symm (t, y) = (t, z) := by
      have := k.left_inv hzO
      rwa [hkz] at this
    refine ⟨by simp only [hsz], ?_⟩
    simp only [hsz]
    exact (congrArg Prod.snd h1).symm

end CWFlow
