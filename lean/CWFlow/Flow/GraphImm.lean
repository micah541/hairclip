import CWFlow.Bridge.GraphChart

/-!
# The graph of an exact form over a Lagrangian immersion, as a map on `M`

Let `F : M → V` be a Lagrangian immersion and `f : M → ℝ` smooth. The chart graphs `chartGraph f α` of
`Bridge/GraphChart.lean` agree on overlaps, so they define a map `graphImm f : M → V`, at every point `x` that is
*covered*: `x` lies in some chart `α` with the graph of `D f_α` at `φ_α x` in the Weinstein domain of `α`.

* `graphImm_eq`: `graphImm f x = chartGraph f α (φ_α x)` for every chart `α` covering `x`;
* `isLagImmersion_graphImm`: if every point is covered, `graphImm f` is a Lagrangian immersion;
* `graphImm_zero`: `graphImm 0 = F`.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) (J : V →L[ℝ] V) (F : M → V)

/-- The chart at `α` covers `x`: `x` is in the chart and the graph of `D f_α` at `φ_α x` lies in the Weinstein
domain of the chart. -/
def Covers (f : M → ℝ) (α x : M) : Prop :=
  x ∈ (chart E α).source ∧ chart E α x ∈ graphDom E J F f α

open Classical in
/-- The graph of `df` over the immersion `F`, as a map on `M` (equal to `F` at points that no chart covers). -/
def graphImm (f : M → ℝ) (x : M) : V :=
  if h : ∃ α, Covers E J F f α x then chartGraph E J F f h.choose (chart E h.choose x) else F x

variable {E J F}

theorem transition_chart {β α x : M} (hβ : x ∈ (chart E β).source) :
    transition E β α (chart E β x) = chart E α x := by
  simp only [transition, Function.comp_apply]
  rw [(chart E β).left_inv hβ]

theorem chart_mem_overlap {β α x : M} (hβ : x ∈ (chart E β).source) (hα : x ∈ (chart E α).source) :
    chart E β x ∈ overlap E β α := by
  refine ⟨(chart E β).map_source hβ, ?_⟩
  rw [mem_preimage, (chart E β).left_inv hβ]
  exact hα

variable [T2Space M]
variable {f : M → ℝ} (hf : ∀ α : M, ContDiffOn ℝ ∞ (chartRep E f α) (chart E α).target)
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hf hJi hJJ hF

/-- **Chart independence**: `graphImm f` is the chart graph in every chart covering the point. -/
theorem graphImm_eq {α x : M} (h : Covers E J F f α x) :
    graphImm E J F f x = chartGraph E J F f α (chart E α x) := by
  have hex : ∃ α, Covers E J F f α x := ⟨α, h⟩
  rw [graphImm, dif_pos hex]
  have hβ := hex.choose_spec
  have ht := transition_chart (α := α) hβ.1
  have hc := chartGraph_transition hf hJi hJJ hF (chart_mem_overlap hβ.1 h.1) hβ.2
    (by rw [ht]; exact h.2)
  rw [ht] at hc
  exact hc.symm

/-- Near a point of the chart at `z` covered by the chart at `α`, `graphImm f` in the chart at `z` is the chart
graph at `α` composed with the transition. -/
theorem chartRep_graphImm_eventuallyEq {z α : M} {y : E} (hy : y ∈ overlap E z α)
    (hα : transition E z α y ∈ graphDom E J F f α) :
    chartRep E (graphImm E J F f) z =ᶠ[𝓝 y] chartGraph E J F f α ∘ transition E z α := by
  have hN : IsOpen (overlap E z α ∩ transition E z α ⁻¹' graphDom E J F f α) :=
    (contDiffOn_transition (E := E) z α).continuousOn.isOpen_inter_preimage
      (isOpen_overlap (E := E) z α) (isOpen_graphDom' hf hJi hJJ hF α)
  filter_upwards [hN.mem_nhds ⟨hy, hα⟩] with w hw
  exact graphImm_eq hf hJi hJJ hF (α := α) (x := (chart E z).symm w) ⟨hw.1.2, hw.2⟩

theorem contDiffAt_transition {z α : M} {y : E} (hy : y ∈ overlap E z α) :
    ContDiffAt ℝ ∞ (transition E z α) y :=
  (contDiffOn_transition (E := E) z α).contDiffAt ((isOpen_overlap (E := E) z α).mem_nhds hy)

theorem contDiffAt_chartRep_graphImm {z α : M} {y : E} (hy : y ∈ overlap E z α)
    (hα : transition E z α y ∈ graphDom E J F f α) :
    ContDiffAt ℝ ∞ (chartRep E (graphImm E J F f) z) y :=
  ((contDiffAt_chartGraph hf hJi hJJ hF hα).comp y
    (contDiffAt_transition hf hJi hJJ hF hy)).congr_of_eventuallyEq
    (chartRep_graphImm_eventuallyEq hf hJi hJJ hF hy hα)

theorem fderiv_chartRep_graphImm {z α : M} {y : E} (hy : y ∈ overlap E z α)
    (hα : transition E z α y ∈ graphDom E J F f α) :
    fderiv ℝ (chartRep E (graphImm E J F f) z) y =
      (fderiv ℝ (chartGraph E J F f α) (transition E z α y)).comp (fderiv ℝ (transition E z α) y) := by
  rw [(chartRep_graphImm_eventuallyEq hf hJi hJJ hF hy hα).fderiv_eq]
  exact fderiv_comp y ((contDiffAt_chartGraph hf hJi hJJ hF hα).differentiableAt (by simp))
    ((contDiffAt_transition hf hJi hJJ hF hy).differentiableAt (by simp))

/-- The chart graphs are Lagrangian. -/
theorem lagrangian_fderiv_chartGraph {α : M} {y : E} (hy : y ∈ graphDom E J F f α) (a b : E) :
    ⟪J (fderiv ℝ (chartGraph E J F f α) y a), fderiv ℝ (chartGraph E J F f α) y b⟫ = 0 :=
  graphFu_lagrangian
    (((contDiffOn_chartPsi hJi hJJ hF α).contDiffAt
      ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)).differentiableAt (by simp))
    (((hf α).contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds hy.1)).of_le two_le_infty)
    (chartPsi_symplectic hJi hJJ hF α _ hy.2) a b

/-- **The graph of `df` is a Lagrangian immersion**, if every point is covered by a chart. -/
theorem isLagImmersion_graphImm (hcov : ∀ x : M, ∃ α, Covers E J F f α x) :
    IsLagImmersion E J (graphImm E J F f) := by
  have key : ∀ z : M, ∀ y ∈ (chart E z).target, ∃ α, y ∈ overlap E z α ∧
      transition E z α y ∈ graphDom E J F f α := by
    intro z y hy
    obtain ⟨α, hα⟩ := hcov ((chart E z).symm y)
    exact ⟨α, ⟨hy, hα.1⟩, hα.2⟩
  refine ⟨fun z y hy => ?_, fun z y hy => ?_, fun z y hy a b => ?_⟩
  · obtain ⟨α, hyα, hα⟩ := key z y hy
    exact (contDiffAt_chartRep_graphImm hf hJi hJJ hF hyα hα).contDiffWithinAt
  · obtain ⟨α, hyα, hα⟩ := key z y hy
    obtain ⟨e, he⟩ := exists_equiv_fderiv_transition (E := E) hyα
    rw [fderiv_chartRep_graphImm hf hJi hJJ hF hyα hα, ← he]
    exact (injective_fderiv_chartGraph hf hJi hJJ hF hα).comp e.injective
  · obtain ⟨α, hyα, hα⟩ := key z y hy
    rw [fderiv_chartRep_graphImm hf hJi hJJ hF hyα hα]
    exact lagrangian_fderiv_chartGraph hf hJi hJJ hF hα _ _

omit hf in
/-- The graph of `d0` is the immersion itself. -/
theorem graphImm_zero : graphImm E J F (0 : M → ℝ) = F := by
  funext x
  unfold graphImm
  split_ifs with h
  · have hβ := h.choose_spec
    have h0 : fderiv ℝ (chartRep E (0 : M → ℝ) h.choose) (chart E h.choose x) = 0 :=
      fderiv_const_apply (0 : ℝ)
    have hmem := hβ.2.2
    simp only [graphMap] at hmem
    rw [h0] at hmem
    simp only [chartGraph, graphFu, Function.comp_apply, graphMap]
    rw [h0, ((weinstein_global hJi hJJ hF).1 h.choose).1.zero_section _ hmem]
    simp only [chartRep, Function.comp_apply]
    rw [(chart E h.choose).left_inv hβ.1]
  · rfl

end CWFlow
