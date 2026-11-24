import CWFlow.Bridge.Existence
import CWFlow.Flow.GraphImm
import CWFlow.Flow.Tangential
import CWFlow.Scalar.GraphTime

/-!
# The family of immersions of a scalar solution, and its velocity

Let `u` solve the scalar equation (`(cwData g hJi hJJ hF).IsSolution g u₀ T u`). The family of immersions is
`t ↦ graphImm (fnOf g (u t))`, the graphs of `d u(t)` over `F`.

* `exists_good`, `covers_of_solution`: every point is covered by a chart at every time;
* `isLagImmersion_solution`: each `graphImm (fnOf g (u t))` is a Lagrangian immersion;
* `contDiffOn_wOf`: the solution is jointly smooth in a chart;
* `hasDerivWithinAt_graphImm`: the velocity of the family at `x`, in a good chart `α`, is the vertical vector
  `graphVel = DΨ_α(0, −D div(JH))`;
* `graphVel_eq`: `graphVel = J∇div(JH) + DF_t (tanCoord)`: the velocity is the flow velocity plus a tangential
  vector (Chen–Warren, Claim 3.3, pointwise).
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow.Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

/-- A time-dependent function in the chart at `α`. -/
def wOf (u : ℝ → SmoothCcTensor g 0 0) (α : M) (q : ℝ × E) : ℝ := chartRep E (fnOf g (u q.1)) α q.2

/-- **A jointly smooth family is jointly smooth in each chart.** -/
theorem contDiffOn_wOf {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0} (hu : JointlySmoothOn g 0 0 (Icc 0 T) u)
    (α : M) : ContDiffOn ℝ ∞ (wOf g u α) (Icc 0 T ×ˢ (chart E α).target) := by
  have h1 : ContDiffOn ℝ ∞ (fun q : ℝ × E => scal E (chartRd g α (u q.1) q.2))
      (Icc 0 T ×ˢ (chart E α).target) := by
    intro q hq
    have h := Vendor.chartRepr_euclid_jointContDiffWithinAt (I := 𝓘(ℝ, E)) g 0 0 u (Icc 0 T) α hu
      hq.1 hq.2
    exact (scal E).contDiff.contDiffAt.comp_contDiffWithinAt q h
  exact h1.congr fun q hq => (scalRd_eq_chartRep g α (u q.1) hq.2).symm

variable {J : V →L[ℝ] V} {F : M → V}

variable (J F) in
/-- The velocity of the graph family at the chart point `y`, in the chart at `α`: the vertical vector
`DΨ_α(0, −D div(JH))`. -/
def graphVel (f : M → ℝ) (α : M) (y : E) : V :=
  fderiv ℝ (chartPsi E J F α) (graphMap (chartRep E f α) y)
    (0, -fderiv ℝ (flowPot J (chartGraph E J F f α)) y)

variable (J F) in
/-- The tangential part of the velocity of the graph family, in the coordinates of the chart at `α`. -/
def tanCoord (f : M → ℝ) (α : M) (y : E) : E :=
  tanPart (fderiv ℝ (chartGraph E J F f α) y)
    (graphVel J F f α y - flowVel J (chartGraph E J F f α) y)

variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

omit [CompactSpace M] [SigmaCompactSpace M] [NeZero (Module.finrank ℝ E)] [CompleteSpace E] in
/-- **The velocity of a graph is the flow velocity plus a tangential vector.** -/
theorem graphVel_eq [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {f : M → ℝ} (hf : ∀ α : M, ContDiffOn ℝ ∞ (chartRep E f α) (chart E α).target) {α : M} {y : E}
    (hy : y ∈ graphDom E J F f α) :
    graphVel J F f α y =
      flowVel J (chartGraph E J F f α) y + fderiv ℝ (chartGraph E J F f α) y (tanCoord J F f α y) := by
  have h := graphFu_vertical_sub_flowVel hJi hJJ hdim
    (((hf α).contDiffAt ((isOpen_chart_target (E := E) α).mem_nhds hy.1)).of_le two_le_infty)
    (((contDiffOn_chartPsi hJi hJJ hF α).contDiffAt
      ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)).differentiableAt (by simp))
    (chartPsi_symplectic hJi hJJ hF α _ hy.2)
    (-fderiv ℝ (flowPot J (chartGraph E J F f α)) y) rfl
  exact eq_add_of_sub_eq' h

variable {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
variable (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
include hu

/-- At every time, every point has a good patch. -/
theorem exists_good {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) :
    ∃ i : PIdx E M, x ∈ patchU E i.1 ∧ good g J F i.1 x (covJet g 0 0 (2 * 2) (u t) x) := by
  obtain ⟨i, hi⟩ := (patchρ E M).exists_pos_of_mem (mem_univ x)
  have hx : x ∈ tsupport (patchρ E M i) := subset_tsupport _ (Function.mem_support.2 hi.ne')
  exact ⟨i, tsupport_patchρ_subset i hx, hu.jets_mem t ht x i hx⟩

/-- At every time, every point is covered by a chart. -/
theorem covers_of_solution {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) :
    ∃ α, Covers E J F (fnOf g (u t)) α x := by
  obtain ⟨i, hi, hgi⟩ := exists_good g hJi hJJ hF hu ht x
  exact ⟨i.1, hi.1, mem_graphDom_of_good g J F (u t) hi hgi⟩

/-- **Each time slice of the family is a Lagrangian immersion.** -/
theorem isLagImmersion_solution {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) :
    IsLagImmersion E J (graphImm E J F (fnOf g (u t))) :=
  isLagImmersion_graphImm (contDiffOn_chartRep_fnOf g (u t)) hJi hJJ hF
    (covers_of_solution g hJi hJJ hF hu ht)

/-- Goodness of the jet of a function in the chart at `α` is membership of the chart graph in the Weinstein
domain. -/
theorem good_of_mem_graphDom (v : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ patchU E α)
    (h : chart E α x ∈ graphDom E J F (fnOf g v) α) :
    good g J F α x (covJet g 0 0 (2 * 2) v x) := by
  unfold good
  rw [cJet_covJet g v hx]
  exact h.2

/-- **The time derivative of the solution in a chart is `−div(JH)` of the chart graph**, near a good point. -/
theorem dtW_wOf_eventuallyEq {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ patchU E α)
    (hgα : good g J F α x (covJet g 0 0 (2 * 2) (u t) x)) :
    (fun z => dtW (wOf g u α) T (t, z)) =ᶠ[𝓝 (chart E α x)]
      fun z => -flowPot J (chartGraph E J F (fnOf g (u t)) α) z := by
  have hy := mem_graphDom_of_good g J F (u t) hα hgα
  have hopen : IsOpen (graphDom E J F (fnOf g (u t)) α ∩ (chart E α).symm ⁻¹' patchU E α) :=
    ((continuousOn_extChartAt_symm (I := 𝓘(ℝ, E)) α).mono fun z hz => hz.1).isOpen_inter_preimage
      (isOpen_graphDom' (contDiffOn_chartRep_fnOf g (u t)) hJi hJJ hF α) (isOpen_patchU α)
  have hmem : chart E α x ∈ graphDom E J F (fnOf g (u t)) α ∩ (chart E α).symm ⁻¹' patchU E α := by
    refine ⟨hy, ?_⟩
    rw [mem_preimage, (chart E α).left_inv hα.1]
    exact hα
  filter_upwards [hopen.mem_nhds hmem] with z hz
  have hzr : chart E α ((chart E α).symm z) = z := (chart E α).right_inv hz.1.1
  have hg' : good g J F α ((chart E α).symm z) (covJet g 0 0 (2 * 2) (u t) ((chart E α).symm z)) :=
    good_of_mem_graphDom g hJi hJJ hF hu (u t) hz.2 (by rw [hzr]; exact hz.1)
  have hd := hasDerivWithinAt_fnOf g hJi hJJ hF hu ht hz.2 hg'
  rw [hzr] at hd
  exact hd.derivWithin (uniqueDiffOn_Icc hu.pos t ht)

/-- **The velocity of the family of immersions** at a good point: the vertical vector `graphVel`. -/
theorem hasDerivWithinAt_graphImm {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ patchU E α)
    (hgα : good g J F α x (covJet g 0 0 (2 * 2) (u t) x)) :
    HasDerivWithinAt (fun s => graphImm E J F (fnOf g (u s)) x)
      (graphVel J F (fnOf g (u t)) α (chart E α x)) (Icc 0 T) t := by
  have hy := mem_graphDom_of_good g J F (u t) hα hgα
  have hw := contDiffOn_wOf g hu.smooth α
  have hΨ : DifferentiableAt ℝ (chartPsi E J F α)
      (graphMap (fun z => wOf g u α (t, z)) (chart E α x)) :=
    ((contDiffOn_chartPsi hJi hJJ hF α).contDiffAt
      ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)).differentiableAt (by simp)
  have hA := hasDerivWithinAt_graphFu_time hu.pos (isOpen_chart_target (E := E) α) hw ht hy.1 hΨ
  have hη : fderiv ℝ (fun z => dtW (wOf g u α) T (t, z)) (chart E α x) =
      -fderiv ℝ (flowPot J (chartGraph E J F (fnOf g (u t)) α)) (chart E α x) := by
    rw [(dtW_wOf_eventuallyEq g hJi hJJ hF hu ht hα hgα).fderiv_eq]
    exact fderiv_neg
  rw [hη] at hA
  have hc : ContinuousWithinAt (fun s => graphMap (fun z => wOf g u α (s, z)) (chart E α x))
      (Icc 0 T) t :=
    continuousWithinAt_const.prodMk
      (hasDerivWithinAt_slice_fderiv hu.pos (isOpen_chart_target (E := E) α) hw ht
        hy.1).continuousWithinAt
  have hev : (fun s => graphImm E J F (fnOf g (u s)) x) =ᶠ[𝓝[Icc 0 T] t]
      fun s => graphFu (chartPsi E J F α) (fun z => wOf g u α (s, z)) (chart E α x) := by
    filter_upwards [hc.preimage_mem_nhdsWithin ((isOpen_chartDom hJi hJJ hF α).mem_nhds hy.2)]
      with s hs
    exact graphImm_eq (contDiffOn_chartRep_fnOf g (u s)) hJi hJJ hF ⟨hα.1, hy.1, hs⟩
  exact hA.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ht)

end CWFlow.Bridge
