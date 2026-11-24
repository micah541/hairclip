import CWFlow.Weinstein.Global

/-!
# A diffeomorphism of `M` read in charts

For a diffeomorphism `φ` of `M`, the chart representative `diffeoRep E φ z α = chart α ∘ φ ∘ (chart z)⁻¹`
is smooth with invertible derivative on its (open) domain, and composing a Lagrangian immersion with `φ`
gives a Lagrangian immersion.
-/

open Set Function Filter

open scoped Manifold Topology ContDiff InnerProductSpace RealInnerProductSpace

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) in
/-- A self-map `φ` of `M` read from the chart at `z` to the chart at `α`. -/
def diffeoRep (φ : M → M) (z α : M) : E → E := chart E α ∘ φ ∘ (chart E z).symm

variable (E) in
/-- The points of the chart at `z` that `φ` maps into the chart at `α`. -/
def diffeoDom (φ : M → M) (z α : M) : Set E :=
  (chart E z).target ∩ (chart E z).symm ⁻¹' (φ ⁻¹' (chart E α).source)

theorem isOpen_diffeoDom (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) (z α : M) :
    IsOpen (diffeoDom E φ z α) := by
  have h1 : ContinuousOn (chart E z).symm (chart E z).target := continuousOn_extChartAt_symm z
  have h2 : IsOpen (φ ⁻¹' (chart E α).source) :=
    (isOpen_extChartAt_source α).preimage φ.continuous
  exact h1.isOpen_inter_preimage (isOpen_chart_target z) h2

theorem contDiffOn_diffeoRep (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) (z α : M) :
    ContDiffOn ℝ ∞ (diffeoRep E φ z α) (diffeoDom E φ z α) := by
  have h1 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E z).symm (chart E z).target :=
    contMDiffOn_extChartAt_symm z
  have h2 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E α) (chartAt E α).source :=
    contMDiffOn_extChartAt
  have hφ : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (φ ∘ (chart E z).symm) (diffeoDom E φ z α) :=
    φ.contMDiff.comp_contMDiffOn (h1.mono inter_subset_left)
  have h3 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E α ∘ (φ ∘ (chart E z).symm))
      (diffeoDom E φ z α) := by
    refine h2.comp hφ fun y hy => ?_
    have := hy.2
    simpa [extChartAt_source] using this
  exact contMDiffOn_iff_contDiffOn.mp h3

theorem mem_diffeoDom {φ : M → M} {z α x : M} (hz : x ∈ (chart E z).source)
    (hα : φ x ∈ (chart E α).source) : chart E z x ∈ diffeoDom E φ z α := by
  refine ⟨(chart E z).map_source hz, ?_⟩
  simp only [mem_preimage]
  rw [(chart E z).left_inv hz]
  exact hα

theorem diffeoRep_chart {φ : M → M} {z α x : M} (hz : x ∈ (chart E z).source) :
    diffeoRep E φ z α (chart E z x) = chart E α (φ x) := by
  simp only [diffeoRep, Function.comp_apply]
  rw [(chart E z).left_inv hz]

theorem diffeoRep_symm_diffeoRep (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) {z α : M} {y : E}
    (hy : y ∈ diffeoDom E φ z α) :
    diffeoRep E φ.symm α z (diffeoRep E φ z α y) = y := by
  have hx : (chart E z).symm y ∈ (chart E z).source := (chart E z).map_target hy.1
  have hφx : φ ((chart E z).symm y) ∈ (chart E α).source := hy.2
  simp only [diffeoRep, Function.comp_apply]
  rw [(chart E α).left_inv hφx]
  have : φ.symm (φ ((chart E z).symm y)) = (chart E z).symm y := φ.symm_apply_apply _
  rw [this]
  exact (chart E z).right_inv hy.1

theorem diffeoRep_mem_diffeoDom (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) {z α : M} {y : E}
    (hy : y ∈ diffeoDom E φ z α) :
    diffeoRep E φ z α y ∈ diffeoDom E φ.symm α z := by
  have hx : (chart E z).symm y ∈ (chart E z).source := (chart E z).map_target hy.1
  have hφx : φ ((chart E z).symm y) ∈ (chart E α).source := hy.2
  refine ⟨(chart E α).map_source hφx, ?_⟩
  simp only [diffeoRep, Function.comp_apply, mem_preimage]
  rw [(chart E α).left_inv hφx]
  have : φ.symm (φ ((chart E z).symm y)) = (chart E z).symm y := φ.symm_apply_apply _
  rw [this]
  exact hx

theorem exists_equiv_fderiv_diffeoRep (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) {z α : M} {y : E}
    (hy : y ∈ diffeoDom E φ z α) :
    ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (diffeoRep E φ z α) y := by
  have hτ := ((contDiffOn_diffeoRep (E := E) φ z α).contDiffAt
    ((isOpen_diffeoDom (E := E) φ z α).mem_nhds hy)).differentiableAt (by simp)
  have hy' := diffeoRep_mem_diffeoDom φ hy
  have hτ' := ((contDiffOn_diffeoRep (E := E) φ.symm α z).contDiffAt
    ((isOpen_diffeoDom (E := E) φ.symm α z).mem_nhds hy')).differentiableAt (by simp)
  have hid : diffeoRep E φ.symm α z ∘ diffeoRep E φ z α =ᶠ[𝓝 y] id := by
    filter_upwards [(isOpen_diffeoDom (E := E) φ z α).mem_nhds hy] with y' hy'' using
      diffeoRep_symm_diffeoRep φ hy''
  have hcomp := fderiv_comp y hτ' hτ
  rw [hid.fderiv_eq, fderiv_id] at hcomp
  have hinj : Function.Injective (fderiv ℝ (diffeoRep E φ z α) y) := by
    intro a b hab
    have := congrArg (fderiv ℝ (diffeoRep E φ.symm α z) (diffeoRep E φ z α y)) hab
    rw [← ContinuousLinearMap.comp_apply, ← ContinuousLinearMap.comp_apply, ← hcomp] at this
    simpa using this
  exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (diffeoRep E φ z α) y : E →ₗ[ℝ] E) hinj
    rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩

theorem chartRep_comp_diffeo (G : M → V) {φ : M → M} {z α : M} {y : E}
    (hy : y ∈ diffeoDom E φ z α) :
    chartRep E (G ∘ φ) z y = chartRep E G α (diffeoRep E φ z α y) := by
  have hφx : φ ((chart E z).symm y) ∈ (chart E α).source := hy.2
  simp only [chartRep, diffeoRep, Function.comp_apply]
  rw [(chart E α).left_inv hφx]

theorem isLagImmersion_comp_diffeo {J : V →L[ℝ] V} {G : M → V} (hG : IsLagImmersion E J G)
    (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) : IsLagImmersion E J (G ∘ φ) := by
  -- local facts at a point
  have key : ∀ (z : M) (y : E), y ∈ (chart E z).target →
      y ∈ diffeoDom E φ z (φ ((chart E z).symm y)) ∧
      diffeoRep E φ z (φ ((chart E z).symm y)) y ∈
        (chart E (φ ((chart E z).symm y))).target := by
    intro z y hy
    have hx : (chart E z).symm y ∈ (chart E z).source := (chart E z).map_target hy
    have hα : φ ((chart E z).symm y) ∈ (chart E (φ ((chart E z).symm y))).source :=
      mem_extChartAt_source _
    refine ⟨⟨hy, hα⟩, ?_⟩
    have := diffeoRep_chart (E := E) (φ := φ) (z := z) (α := φ ((chart E z).symm y)) hx
    rw [(chart E z).right_inv hy] at this
    rw [this]
    exact (chart E _).map_source hα
  have hdiff : ∀ (z : M) (y : E), y ∈ (chart E z).target →
      ContDiffAt ℝ ∞ (chartRep E (G ∘ φ) z) y ∧
      fderiv ℝ (chartRep E (G ∘ φ) z) y =
        (fderiv ℝ (chartRep E G (φ ((chart E z).symm y)))
          (diffeoRep E φ z (φ ((chart E z).symm y)) y)).comp
        (fderiv ℝ (diffeoRep E φ z (φ ((chart E z).symm y))) y) := by
    intro z y hy
    obtain ⟨hd, ht⟩ := key z y hy
    set α := φ ((chart E z).symm y) with hα
    have h1 : ContDiffAt ℝ ∞ (diffeoRep E φ z α) y :=
      (contDiffOn_diffeoRep (E := E) φ z α).contDiffAt
        ((isOpen_diffeoDom (E := E) φ z α).mem_nhds hd)
    have h2 : ContDiffAt ℝ ∞ (chartRep E G α) (diffeoRep E φ z α y) :=
      (hG.smooth α).contDiffAt ((isOpen_chart_target α).mem_nhds ht)
    have hev : chartRep E (G ∘ φ) z =ᶠ[𝓝 y] chartRep E G α ∘ diffeoRep E φ z α := by
      filter_upwards [(isOpen_diffeoDom (E := E) φ z α).mem_nhds hd] with y' hy' using
        chartRep_comp_diffeo G hy'
    refine ⟨(h2.comp y h1).congr_of_eventuallyEq hev, ?_⟩
    rw [hev.fderiv_eq]
    exact fderiv_comp y (h2.differentiableAt (by simp)) (h1.differentiableAt (by simp))
  refine ⟨?_, ?_, ?_⟩
  · intro z y hy
    exact (hdiff z y hy).1.contDiffWithinAt
  · intro z y hy
    obtain ⟨hd, ht⟩ := key z y hy
    rw [(hdiff z y hy).2]
    obtain ⟨e, he⟩ := exists_equiv_fderiv_diffeoRep (E := E) φ hd
    rw [← he]
    exact (hG.inj _ _ ht).comp e.injective
  · intro z y hy u w
    obtain ⟨hd, ht⟩ := key z y hy
    rw [(hdiff z y hy).2]
    simpa using hG.lag _ _ ht _ _

end CWFlow
