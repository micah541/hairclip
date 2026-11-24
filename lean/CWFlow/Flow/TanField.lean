import CWFlow.Flow.JointFamily
import CWFlow.Flow.TanAlgebra

/-!
# The tangential part of the velocity, as a vector field on `M`

For a solution `u` of the scalar equation, with family of immersions `G_t = graphImm (fnOf g (u t))`:

* `eRep α (t, y)`: minus the tangential coordinates, `tanPart (D G_t) (−∂ₜ G_t)`, in the chart at `α`, jointly smooth on
  `Icc 0 T ×ˢ (chart E α).target` (`contDiffOn_eRep`);
* `eRep_transition`: they transform as a vector under chart transitions;
* `tanField t x`: the tangent vector at `x` (read in the chart at `x` itself); its reading in any chart `α` is
  `eRep α` (`trivAt_tanField`);
* `contMDiffOn_tanField`: it is a smooth section of `TM` over `Icc 0 T × M`.
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

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) (J : V →L[ℝ] V) (F : M → V)

/-- Minus the tangential coordinates of the velocity of the family, in the chart at `α`. -/
def eRep (u : ℝ → SmoothCcTensor g 0 0) (T : ℝ) (α : M) (p : ℝ × E) : E :=
  tanPart (fderiv ℝ (fun z => famRep g J F u α (p.1, z)) p.2) (-dtW (famRep g J F u α) T p)

/-- Minus the tangential part of the velocity of the family, as a tangent vector at `x`: the field whose flow
reparametrizes the family into a solution of the flow. -/
def tanField (u : ℝ → SmoothCcTensor g 0 0) (T : ℝ) (t : ℝ) (x : M) : TangentSpace 𝓘(ℝ, E) x :=
  (show E from eRep g J F u T x (t, chart E x x))

variable {J F}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
variable {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
variable (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
include hu

theorem contDiffOn_eRep (α : M) :
    ContDiffOn ℝ ∞ (eRep g J F u T α) (Icc 0 T ×ˢ (chart E α).target) := by
  have hG := contDiffOn_famRep g hJi hJJ hF hu α
  have hsl := contDiffOn_slice_fderiv hu.pos (isOpen_chart_target (E := E) α) hG
  have hdt := contDiffOn_dtW hu.pos (isOpen_chart_target (E := E) α) hG
  intro p hp
  have hinj : Function.Injective (fderiv ℝ (fun z => famRep g J F u α (p.1, z)) p.2) :=
    (isLagImmersion_solution g hJi hJJ hF hu hp.1).inj α p.2 hp.2
  exact (contDiffAt_tanPart hinj _).comp_contDiffWithinAt p ((hsl p hp).prodMk (hdt p hp).neg)

/-- **The tangential coordinates transform as a vector.** -/
theorem eRep_transition {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α β x : M} (hα : x ∈ (chart E α).source)
    (hβ : x ∈ (chart E β).source) :
    fderiv ℝ (transition E α β) (chart E α x) (eRep g J F u T α (t, chart E α x)) =
      eRep g J F u T β (t, chart E β x) := by
  have hG := isLagImmersion_solution g hJi hJJ hF hu ht
  have hy : chart E α x ∈ overlap E α β := chart_mem_overlap hα hβ
  have hτ : transition E α β (chart E α x) = chart E β x := transition_chart hα
  obtain ⟨e, he⟩ := exists_equiv_fderiv_transition (E := E) hy
  have heq : chartRep E (graphImm E J F (fnOf g (u t))) α =ᶠ[𝓝 (chart E α x)]
      chartRep E (graphImm E J F (fnOf g (u t))) β ∘ transition E α β := by
    filter_upwards [(isOpen_overlap (E := E) α β).mem_nhds hy] with w hw
    exact chartRep_transition _ hw
  have hdβ : DifferentiableAt ℝ (chartRep E (graphImm E J F (fnOf g (u t))) β)
      (transition E α β (chart E α x)) :=
    ((hG.smooth β).contDiffAt ((isOpen_chart_target (E := E) β).mem_nhds
      (mapsTo_transition (E := E) α β hy))).differentiableAt (by simp)
  have hdτ : DifferentiableAt ℝ (transition E α β) (chart E α x) :=
    ((contDiffOn_transition (E := E) α β).contDiffAt
      ((isOpen_overlap (E := E) α β).mem_nhds hy)).differentiableAt (by simp)
  have hL : fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x) =
      (fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) β) (chart E β x)).comp
        (e : E →L[ℝ] E) := by
    rw [heq.fderiv_eq, fderiv_comp _ hdβ hdτ, hτ, he]
  have hv : dtW (famRep g J F u α) T (t, chart E α x) = dtW (famRep g J F u β) T (t, chart E β x) := by
    show derivWithin (fun s => famRep g J F u α (s, chart E α x)) (Icc 0 T) t =
      derivWithin (fun s => famRep g J F u β (s, chart E β x)) (Icc 0 T) t
    congr 1
    funext s
    simp only [famRep, chartRep, Function.comp_apply]
    rw [(chart E α).left_inv hα, (chart E β).left_inv hβ]
  have key : eRep g J F u T α (t, chart E α x) = e.symm (eRep g J F u T β (t, chart E β x)) := by
    show tanPart (fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x))
        (-dtW (famRep g J F u α) T (t, chart E α x)) =
      e.symm (tanPart (fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) β) (chart E β x))
        (-dtW (famRep g J F u β) T (t, chart E β x)))
    rw [hL, hv]
    exact tanPart_comp_equiv (hG.inj β _ ((chart E β).map_source hβ)) e _
  rw [← he, key]
  exact e.apply_symm_apply _

/-- **The reading of the tangential field in the chart at `α`.** -/
theorem trivAt_tanField {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ (chart E α).source) :
    (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).continuousLinearMapAt ℝ x (tanField g J F u T t x) =
      eRep g J F u T α (t, chart E α x) := by
  have hs : x ∈ (chartAt E α).source := by rwa [extChartAt_source] at hα
  rw [TangentBundle.continuousLinearMapAt_trivializationAt_eq_core hs,
    tangentBundleCore_coordChange_achart, ModelWithCorners.range_eq_univ, fderivWithin_univ]
  exact eRep_transition g hJi hJJ hF hu ht (mem_extChartAt_source (I := 𝓘(ℝ, E)) x) hα

/-- **The tangential field is a smooth section of `TM` over `Icc 0 T × M`.** -/
theorem contMDiffOn_tanField :
    ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, E).prod 𝓘(ℝ, E)) ∞
      (fun q : ℝ × M => (TotalSpace.mk' E q.2 (tanField g J F u T q.1 q.2) : TangentBundle 𝓘(ℝ, E) M))
      (Icc 0 T ×ˢ univ) := by
  rintro ⟨t, x⟩ ⟨ht, -⟩
  rw [contMDiffWithinAt_totalSpace]
  refine ⟨contMDiffWithinAt_snd, ?_⟩
  have hsrc : ∀ q : ℝ × M, q ∈ Icc 0 T ×ˢ (chart E x).source → q.2 ∈ (chartAt E x).source := by
    intro q hq
    have h := hq.2
    rwa [extChartAt_source] at h
  have h1 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
      (fun q : ℝ × M => (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source) (t, x) :=
    contMDiffWithinAt_fst.prodMk_space
      ((contMDiffOn_extChartAt (I := 𝓘(ℝ, E)) (n := ∞) (x := x)).comp contMDiffOn_snd hsrc (t, x)
        ⟨ht, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩)
  have h2 : ContMDiffWithinAt 𝓘(ℝ, ℝ × E) 𝓘(ℝ, E) ∞ (eRep g J F u T x)
      (Icc 0 T ×ˢ (chart E x).target) (t, chart E x x) :=
    (contDiffOn_eRep g hJi hJJ hF hu x (t, chart E x x)
      ⟨ht, mem_extChartAt_target (I := 𝓘(ℝ, E)) x⟩).contMDiffWithinAt
  have h3 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞
      (fun q : ℝ × M => eRep g J F u T x (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source) (t, x) := by
    have hm : MapsTo (fun q : ℝ × M => (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source)
        (Icc 0 T ×ˢ (chart E x).target) := fun q hq => ⟨hq.1, (chart E x).map_source hq.2⟩
    exact ContMDiffWithinAt.comp (g := eRep g J F u T x)
      (f := fun q : ℝ × M => (q.1, chart E x q.2)) (t, x) h2 h1 hm
  have hev : (fun q : ℝ × M => (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x
        (TotalSpace.mk' E q.2 (tanField g J F u T q.1 q.2) : TangentBundle 𝓘(ℝ, E) M)).2)
      =ᶠ[𝓝[Icc 0 T ×ˢ (chart E x).source] (t, x)]
      fun q : ℝ × M => eRep g J F u T x (q.1, chart E x q.2) := by
    filter_upwards [self_mem_nhdsWithin] with q hq
    have hb : q.2 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x).baseSet := by
      rw [TangentBundle.trivializationAt_baseSet]
      exact hsrc q hq
    rw [← trivAt_tanField g hJi hJJ hF hu hq.1 hq.2]
    show _ = (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x).linearMapAt ℝ q.2 _
    rw [Trivialization.coe_linearMapAt_of_mem _ hb]
  refine (h3.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin
    ⟨ht, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩)).mono_of_mem_nhdsWithin ?_
  exact mem_nhdsWithin.2 ⟨univ ×ˢ (chart E x).source,
    isOpen_univ.prod (isOpen_extChartAt_source (I := 𝓘(ℝ, E)) x),
    ⟨mem_univ _, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩, fun p hp => ⟨hp.2.1, hp.1.2⟩⟩

end CWFlow.Bridge
