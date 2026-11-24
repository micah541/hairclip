import CWFlow.Weinstein.Global
import CWFlow.Scalar.JetFunction
import CWFlow.Scalar.Reparam

/-!
# Graphs of exact forms over the immersion, chart by chart

Let `F : M → V` be a Lagrangian immersion and `f : M → ℝ` smooth. In the chart at `α`, the graph of
`df` over `F` is `chartGraph f α = chartPsi α ∘ (y ↦ (y, D f_α(y)))`, `f_α = f ∘ φ_α⁻¹`, defined where
`(y, D f_α(y))` lies in the Weinstein domain of the chart (`graphDom f α`, open).

* `graphMap_transition`: the cotangent lift of a chart transition maps the graph of `D f_α` to the graph
  of `D f_β`;
* `chartGraph_transition`: hence (Weinstein compatibility, `chartPsi_compat`) the chart graphs agree;
* `flowPot_chartGraph_transition`: hence (`flowPot_comp`) `div(JH)` of the chart graphs agree: the
  operator `f ↦ −div(JH)[F_f]` is chart independent;
* `neg_flowPot_chartGraph`: in each chart it is `Pjet` of the coordinate 4-jet of `f_α`.
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

/-- The graph of `df` over the immersion, in the chart at `α`. -/
def chartGraph (f : M → ℝ) (α : M) : E → V := graphFu (chartPsi E J F α) (chartRep E f α)

/-- The points of the chart at `α` where the graph of `df` lies in the Weinstein domain. -/
def graphDom (f : M → ℝ) (α : M) : Set E :=
  {y | y ∈ (chart E α).target ∧ graphMap (chartRep E f α) y ∈ chartDom E J F α}

variable {E J F}

/-- **Locality of `div(JH)`**: it only depends on the germ of the immersion. -/
theorem flowPot_congr_of_eventuallyEq {X₁ X₂ : E → V} {y : E} (h : X₁ =ᶠ[𝓝 y] X₂) :
    flowPot J X₁ y = flowPot J X₂ y := by
  have h1 : fderiv ℝ X₁ =ᶠ[𝓝 y] fderiv ℝ X₂ := h.fderiv
  have h2 : fderiv ℝ (fderiv ℝ X₁) =ᶠ[𝓝 y] fderiv ℝ (fderiv ℝ X₂) := h1.fderiv
  have h3 : (fun z => J (meanCurv X₁ z)) =ᶠ[𝓝 y] fun z => J (meanCurv X₂ z) := by
    filter_upwards [h1, h2] with z hz1 hz2
    simp only [meanCurv, hz1, hz2]
  simp only [flowPot, h1.self_of_nhds, h3.fderiv_eq]

variable {f : M → ℝ} (hf : ∀ α : M, ContDiffOn ℝ ∞ (chartRep E f α) (chart E α).target)
include hf

/-- The cotangent lift of a chart transition maps the graph of `D f_α` to the graph of `D f_β`. -/
theorem graphMap_transition {α β : M} {y : E} (hy : y ∈ overlap E α β) :
    cotLift (transition E α β) (graphMap (chartRep E f α) y) =
      graphMap (chartRep E f β) (transition E α β y) := by
  obtain ⟨e, he⟩ := exists_equiv_fderiv_transition (E := E) hy
  have hτd : DifferentiableAt ℝ (transition E α β) y :=
    ((contDiffOn_transition (E := E) α β).contDiffAt
      ((isOpen_overlap (E := E) α β).mem_nhds hy)).differentiableAt (by simp)
  have hfβ : DifferentiableAt ℝ (chartRep E f β) (transition E α β y) :=
    ((hf β).contDiffAt ((isOpen_chart_target (E := E) β).mem_nhds
      (mapsTo_transition (E := E) α β hy))).differentiableAt (by simp)
  have heq : chartRep E f α =ᶠ[𝓝 y] chartRep E f β ∘ transition E α β := by
    filter_upwards [(isOpen_overlap (E := E) α β).mem_nhds hy] with w hw
    exact chartRep_transition f hw
  have hD : fderiv ℝ (chartRep E f α) y =
      (fderiv ℝ (chartRep E f β) (transition E α β y)).comp (fderiv ℝ (transition E α β) y) := by
    rw [heq.fderiv_eq, fderiv_comp y hfβ hτd]
  simp only [graphMap, cotLift]
  refine Prod.ext rfl ?_
  simp only
  rw [hD, ← he, ContinuousLinearMap.inverse_equiv]
  ext v
  simp

theorem isOpen_graphDom {α : M} (hdom : IsOpen (chartDom E J F α)) : IsOpen (graphDom E J F f α) := by
  have hS := isOpen_chart_target (E := E) α
  have hc : ContinuousOn (graphMap (chartRep E f α)) (chart E α).target :=
    continuousOn_id.prodMk ((hf α).continuousOn_fderiv_of_isOpen hS (by simp))
  exact hc.isOpen_inter_preimage hS hdom

/-- The graph point of `D f_α` is a smooth function of `y` on the chart target. -/
theorem contDiffAt_graphMap {α : M} {y : E} (hy : y ∈ (chart E α).target) :
    ContDiffAt ℝ ∞ (graphMap (chartRep E f α)) y :=
  (contDiffOn_id.prodMk ((hf α).fderiv_of_isOpen (isOpen_chart_target (E := E) α) (m := ∞)
    le_rfl)).contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds hy)

variable [T2Space M]

variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

omit hf in
/-- The chart Weinstein maps pull `ω₀` back to the canonical form, in inner-product form. -/
theorem chartPsi_symplectic (α : M) :
    ∀ q ∈ chartDom E J F α, ∀ a b,
      ⟪J (fderiv ℝ (chartPsi E J F α) q a), fderiv ℝ (chartPsi E J F α) q b⟫ = omegaCanVal a b := by
  intro q hq a b
  rw [← pullback2_omega0 hJi hJJ]
  exact ((weinstein_global hJi hJJ hF).1 α).1.symplectic q hq a b

omit hf in
theorem isOpen_chartDom (α : M) : IsOpen (chartDom E J F α) :=
  ((weinstein_global hJi hJJ hF).1 α).1.isOpen

omit hf in
theorem contDiffOn_chartPsi (α : M) : ContDiffOn ℝ ∞ (chartPsi E J F α) (chartDom E J F α) :=
  ((weinstein_global hJi hJJ hF).1 α).1.contDiffOn

theorem isOpen_graphDom' (α : M) : IsOpen (graphDom E J F f α) :=
  isOpen_graphDom hf (isOpen_chartDom hJi hJJ hF α)

theorem contDiffAt_chartGraph {α : M} {y : E} (hy : y ∈ graphDom E J F f α) :
    ContDiffAt ℝ ∞ (chartGraph E J F f α) y :=
  ((contDiffOn_chartPsi hJi hJJ hF α).contDiffAt
    ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)).comp y (contDiffAt_graphMap hf hy.1)

theorem injective_fderiv_chartGraph {α : M} {y : E} (hy : y ∈ graphDom E J F f α) :
    Function.Injective (fderiv ℝ (chartGraph E J F f α) y) :=
  injective_fderiv_graphFu
    (((contDiffOn_chartPsi hJi hJJ hF α).contDiffAt
      ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)).differentiableAt (by simp))
    (((hf α).contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds hy.1)).of_le two_le_infty)
    (chartPsi_symplectic hJi hJJ hF α _ hy.2)

/-- **Compatibility of the chart graphs.** -/
theorem chartGraph_transition {α β : M} {y : E} (hy : y ∈ overlap E α β)
    (hyα : y ∈ graphDom E J F f α) (hyβ : transition E α β y ∈ graphDom E J F f β) :
    chartGraph E J F f β (transition E α β y) = chartGraph E J F f α y := by
  have hq' : cotLift (transition E α β) (graphMap (chartRep E f α) y) ∈ chartDom E J F β := by
    rw [graphMap_transition hf hy]
    exact hyβ.2
  simp only [chartGraph, graphFu, Function.comp_apply]
  rw [← graphMap_transition hf hy]
  exact chartPsi_compat hF hyα.2 hy hq'

/-- **Chart independence of `div(JH)` on graphs.** -/
theorem flowPot_chartGraph_transition {α β : M} {y : E} (hy : y ∈ overlap E α β)
    (hyα : y ∈ graphDom E J F f α) (hyβ : transition E α β y ∈ graphDom E J F f β) :
    flowPot J (chartGraph E J F f α) y = flowPot J (chartGraph E J F f β) (transition E α β y) := by
  set τ := transition E α β
  have hτc : ContinuousOn τ (overlap E α β) := (contDiffOn_transition (E := E) α β).continuousOn
  set N := (overlap E α β ∩ graphDom E J F f α) ∩ τ ⁻¹' graphDom E J F f β with hNdef
  have hN : IsOpen N :=
    (hτc.mono inter_subset_left).isOpen_inter_preimage
      ((isOpen_overlap (E := E) α β).inter (isOpen_graphDom' hf hJi hJJ hF α))
      (isOpen_graphDom' hf hJi hJJ hF β)
  have hyN : y ∈ N := ⟨⟨hy, hyα⟩, hyβ⟩
  have heq : chartGraph E J F f α =ᶠ[𝓝 y] chartGraph E J F f β ∘ τ := by
    filter_upwards [hN.mem_nhds hyN] with w hw
    exact (chartGraph_transition hf hJi hJJ hF hw.1.1 hw.1.2 hw.2).symm
  rw [flowPot_congr_of_eventuallyEq heq]
  exact flowPot_comp hN ((contDiffOn_transition (E := E) α β).mono fun w hw => hw.1.1)
    (fun _ hw => exists_equiv_fderiv_transition (E := E) hw.1.1)
    (fun w hw => contDiffAt_chartGraph hf hJi hJJ hF hw.2)
    (fun w hw => injective_fderiv_chartGraph hf hJi hJJ hF hw.2) J hyN

/-- **The operator in a chart is `Pjet` of the coordinate 4-jet.** -/
theorem neg_flowPot_chartGraph {α : M} {y : E} (hy : y ∈ graphDom E J F f α) :
    -flowPot J (chartGraph E J F f α) y = Pjet J (chartPsi E J F α) (y, jetOf (chartRep E f α) y) :=
  neg_flowPot_graphFu_eq_Pjet (isOpen_chartDom hJi hJJ hF α) (contDiffOn_chartPsi hJi hJJ hF α)
    (chartPsi_symplectic hJi hJJ hF α) (isOpen_graphDom' hf hJi hJJ hF α)
    ((hf α).mono fun w hw => hw.1) (fun w hw => hw.2) hy

end CWFlow
