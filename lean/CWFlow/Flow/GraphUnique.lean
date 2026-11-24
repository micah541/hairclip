import CWFlow.Flow.FlowExistence

/-!
# Uniqueness of the flow among graphs over `F`

`IsFlow E J T Φ`: `Φ` is a jointly smooth family of Lagrangian immersions on `[0, T]` moving by `J∇div(JH)`.

Let `u` solve the scalar equation and `G_t = graphImm (fnOf g (u t))`. Suppose a flow `Φ` is a graph up to a
reparametrization: `Φ_t ∘ ψ_t = G_t` for a jointly smooth family of maps `ψ_t : M → M`, `ψ_0 = id`. Then:

* `fderivWithin_psiRep_eRep`: `ψ` solves the transport equation `∂ₜψ + Dψ(X) = 0`, `X = tanField`;
* `psi_comp_chi`: hence `ψ_t ∘ χ_t = id` for every flow `χ` of `X` (so `ψ_t` is the inverse of a diffeomorphism);
* `flow_eq_graph_reparam`: `Φ_t = G_t ∘ χ_t`, the solution constructed in `exists_flow`;
* `flow_unique_of_graph`: **two flows from `F` that are graphs of scalar solutions up to reparametrization
  coincide** (the scalar solutions agree by the Mantegazza–Martinazzi uniqueness theorem).
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

namespace CWFlow

section Defs

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) in
/-- **A solution of the flow `∂ₜF = J∇div(JH)` on `[0, T]`**: a family of Lagrangian immersions, jointly smooth
in `(t, x)`, whose velocity at every point, read in any chart, is `J∇div(JH)` of the current immersion
(one-sided at the endpoints). -/
structure IsFlow (J : V →L[ℝ] V) (T : ℝ) (Φ : ℝ → M → V) : Prop where
  smooth : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ)
  imm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)
  vel : ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
    HasDerivWithinAt (fun s => Φ s x) (flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t

variable (E) in
/-- A family `Φ : ℝ → M → V` in the chart at `β`. -/
def flowRep (Φ : ℝ → M → V) (β : M) (p : ℝ × E) : V := Φ p.1 ((chart E β).symm p.2)

variable (E) in
/-- A family of maps `ψ : ℝ → M → M` read from the chart at `α` to the chart at `β`. -/
def psiRep (ψ : ℝ → M → M) (α β : M) (q : ℝ × E) : E := chart E β (ψ q.1 ((chart E α).symm q.2))

end Defs

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-! ### Chart readings of jointly smooth families -/

theorem contMDiffOn_chartPre (α : M) (T : ℝ) :
    ContMDiffOn 𝓘(ℝ, ℝ × E) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
      (fun p : ℝ × E => ((p.1, (chart E α).symm p.2) : ℝ × M)) (Icc 0 T ×ˢ (chart E α).target) := by
  have h1 : ContMDiffOn 𝓘(ℝ, ℝ × E) 𝓘(ℝ, ℝ) ∞ (fun p : ℝ × E => p.1)
      (Icc 0 T ×ˢ (chart E α).target) := contDiff_fst.contMDiff.contMDiffOn
  have h2 : ContMDiffOn 𝓘(ℝ, ℝ × E) 𝓘(ℝ, E) ∞ (fun p : ℝ × E => (chart E α).symm p.2)
      (Icc 0 T ×ˢ (chart E α).target) :=
    ContMDiffOn.comp (g := (chart E α).symm) (f := fun p : ℝ × E => p.2)
      (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E)) α) contDiff_snd.contMDiff.contMDiffOn
      fun p hp => hp.2
  exact h1.prodMk h2

theorem contDiffOn_flowRep {Φ : ℝ → M → V} {T : ℝ}
    (hΦ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2)
      (Icc 0 T ×ˢ univ)) (β : M) :
    ContDiffOn ℝ ∞ (flowRep E Φ β) (Icc 0 T ×ˢ (chart E β).target) := by
  have h3 := ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2) hΦ (contMDiffOn_chartPre (E := E) β T)
    (fun p hp => ⟨hp.1, mem_univ _⟩)
  exact contMDiffOn_iff_contDiffOn.mp h3

section Psi

variable {ψ : ℝ → M → M} {T : ℝ}
variable (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
  (Icc 0 T ×ˢ univ))
include hψ

theorem contMDiffOn_psiPre (α : M) :
    ContMDiffOn 𝓘(ℝ, ℝ × E) 𝓘(ℝ, E) ∞ (fun p : ℝ × E => ψ p.1 ((chart E α).symm p.2))
      (Icc 0 T ×ˢ (chart E α).target) :=
  ContMDiffOn.comp (g := fun q : ℝ × M => ψ q.1 q.2) hψ (contMDiffOn_chartPre (E := E) α T)
    (fun p hp => ⟨hp.1, mem_univ _⟩)

theorem contDiffWithinAt_psiRep {α β : M} {t : ℝ} {z : E} (ht : t ∈ Icc (0 : ℝ) T)
    (hz : z ∈ (chart E α).target) (hβ : ψ t ((chart E α).symm z) ∈ (chart E β).source) :
    ContDiffWithinAt ℝ ∞ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, z) := by
  have hβ' : ψ t ((chart E α).symm z) ∈ (chartAt E β).source := by rwa [extChartAt_source] at hβ
  have h4 : ContMDiffAt 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E β) (ψ t ((chart E α).symm z)) :=
    contMDiffAt_extChartAt' hβ'
  have h5 := ContMDiffAt.comp_contMDiffWithinAt (g := chart E β)
    (f := fun p : ℝ × E => ψ p.1 ((chart E α).symm p.2)) (t, z) h4
    (contMDiffOn_psiPre hψ α (t, z) ⟨ht, hz⟩)
  exact contMDiffWithinAt_iff_contDiffWithinAt.mp h5

theorem contMDiff_slice {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (ψ t) := by
  intro x
  have hc : ContMDiffAt 𝓘(ℝ, E) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞ (fun y : M => ((t, y) : ℝ × M)) x :=
    contMDiffAt_const.prodMk contMDiffAt_id
  have h := ContMDiffWithinAt.comp (g := fun q : ℝ × M => ψ q.1 q.2)
    (f := fun y : M => ((t, y) : ℝ × M)) (s := univ) x (hψ (t, x) ⟨ht, mem_univ _⟩)
    hc.contMDiffWithinAt (fun y _ => ⟨ht, mem_univ _⟩)
  exact contMDiffWithinAt_univ.mp h

end Psi

theorem isOpen_diffeoDom' {φ : M → M} (hφ : Continuous φ) (z α : M) : IsOpen (diffeoDom E φ z α) := by
  have h1 : ContinuousOn (chart E z).symm (chart E z).target := continuousOn_extChartAt_symm z
  have h2 : IsOpen (φ ⁻¹' (chart E α).source) := (isOpen_extChartAt_source α).preimage hφ
  exact h1.isOpen_inter_preimage (isOpen_chart_target z) h2

theorem contDiffOn_diffeoRep' {φ : M → M} (hφ : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ φ) (z α : M) :
    ContDiffOn ℝ ∞ (diffeoRep E φ z α) (diffeoDom E φ z α) := by
  have h1 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E z).symm (chart E z).target :=
    contMDiffOn_extChartAt_symm z
  have h2 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E α) (chartAt E α).source :=
    contMDiffOn_extChartAt
  have hφ' : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (φ ∘ (chart E z).symm) (diffeoDom E φ z α) :=
    hφ.comp_contMDiffOn (h1.mono inter_subset_left)
  have h3 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (chart E α ∘ (φ ∘ (chart E z).symm))
      (diffeoDom E φ z α) := by
    refine h2.comp hφ' fun y hy => ?_
    have := hy.2
    simpa [extChartAt_source] using this
  exact contMDiffOn_iff_contDiffOn.mp h3

/-! ### A curve whose chart readings have zero derivative is constant -/

omit [NeZero (Module.finrank ℝ E)] [CompleteSpace E] [IsManifold 𝓘(ℝ, E) ∞ M] [CompactSpace M]
  [SigmaCompactSpace M] in
theorem eq_of_chart_hasDerivWithinAt_zero {c : ℝ → M} {T : ℝ} (hc : ContinuousOn c (Icc 0 T))
    (hd : ∀ t ∈ Icc (0 : ℝ) T, ∀ β : M, c t ∈ (chart E β).source →
      HasDerivWithinAt (fun s => chart E β (c s)) 0 (Icc 0 T) t) :
    ∀ t ∈ Icc (0 : ℝ) T, c t = c 0 := by
  have hsub : Icc (0 : ℝ) T ⊆ {s | c s = c 0} := by
    refine IsClosed.Icc_subset_of_forall_mem_nhdsWithin ?_ rfl ?_
    · have h := hc.preimage_isClosed_of_isClosed isClosed_Icc (isClosed_singleton (x := c 0))
      rw [inter_comm] at h
      exact h
    · rintro x ⟨hx, hxI⟩
      have hxI' : x ∈ Icc (0 : ℝ) T := ⟨hxI.1, hxI.2.le⟩
      have hsrc : c x ∈ (chart E (c x)).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) _
      have hmem : c ⁻¹' (chart E (c x)).source ∈ 𝓝[Icc 0 T] x :=
        (hc x hxI').preimage_mem_nhdsWithin
          ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) (c x)).mem_nhds hsrc)
      obtain ⟨ε, hε, hεsub⟩ : ∃ ε > 0, Ioo (x - ε) (x + ε) ∩ Icc 0 T ⊆ c ⁻¹' (chart E (c x)).source := by
        rw [mem_nhdsWithin_iff_exists_mem_nhds_inter] at hmem
        obtain ⟨U, hU, hUs⟩ := hmem
        obtain ⟨ε, hε, hεU⟩ := Metric.mem_nhds_iff.mp hU
        refine ⟨ε, hε, fun s hs => hUs ⟨hεU ?_, hs.2⟩⟩
        rw [Real.ball_eq_Ioo]
        exact hs.1
      have hbr : min T (x + ε / 2) ≤ x + ε / 2 := min_le_right _ _
      have hxb : x < min T (x + ε / 2) := lt_min hxI.2 (by linarith)
      have hbT : min T (x + ε / 2) ≤ T := min_le_left _ _
      have hIcc : Icc x (min T (x + ε / 2)) ⊆ Icc 0 T := Icc_subset_Icc hxI.1 hbT
      have hsrcs : ∀ s ∈ Icc x (min T (x + ε / 2)), c s ∈ (chart E (c x)).source := fun s hs =>
        hεsub ⟨⟨by linarith [hs.1], by linarith [hs.2]⟩, hIcc hs⟩
      have hcont : ContinuousOn (fun s => chart E (c x) (c s)) (Icc x (min T (x + ε / 2))) :=
        ContinuousOn.comp (g := chart E (c x)) (f := c)
          (continuousOn_extChartAt (I := 𝓘(ℝ, E)) (c x)) (hc.mono hIcc) hsrcs
      have hconst := constant_of_has_deriv_right_zero hcont (fun s hs => by
        have hsI : s ∈ Icc (0 : ℝ) T := hIcc ⟨hs.1, hs.2.le⟩
        have h := hd s hsI (c x) (hsrcs s ⟨hs.1, hs.2.le⟩)
        have hsT : s < T := lt_of_lt_of_le hs.2 hbT
        exact h.mono_of_mem_nhdsWithin
          (mem_of_superset (Icc_mem_nhdsGE hsT) (Icc_subset_Icc hsI.1 le_rfl)))
      rw [mem_nhdsGT_iff_exists_Ioc_subset]
      refine ⟨min T (x + ε / 2), hxb, fun s hs => ?_⟩
      have hs' : s ∈ Icc x (min T (x + ε / 2)) := ⟨hs.1.le, hs.2⟩
      have h1 : chart E (c x) (c s) = chart E (c x) (c x) := hconst s hs'
      have h2 : c s = c x := (chart E (c x)).injOn (hsrcs s hs') hsrc h1
      show c s = c 0
      rw [h2, hx]
  exact fun t ht => hsub ht

/-! ### The transport equation for the reparametrization -/

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
variable (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)

omit [NeZero (Module.finrank ℝ E)] [CompleteSpace E] [CompactSpace M] [SigmaCompactSpace M]
  [FiniteDimensional ℝ V] in
/-- The field `tanField` at times in `[0, T]` only depends on the family there. -/
theorem tanField_congr {T : ℝ} {u₁ u₂ : ℝ → SmoothCcTensor g 0 0}
    (h : ∀ t ∈ Icc (0 : ℝ) T, u₁ t = u₂ t) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) :
    tanField g J F u₁ T t x = tanField g J F u₂ T t x := by
  have hfam : ∀ α : M, ∀ s ∈ Icc (0 : ℝ) T, ∀ z : E,
      famRep g J F u₁ α (s, z) = famRep g J F u₂ α (s, z) := by
    intro α s hs z
    simp only [famRep]
    rw [h s hs]
  have hdt : ∀ (α : M) (z : E),
      dtW (famRep g J F u₁ α) T (t, z) = dtW (famRep g J F u₂ α) T (t, z) := fun α z =>
    derivWithin_congr (fun s hs => hfam α s hs z) (hfam α t ht z)
  have hfun : (fun z => famRep g J F u₁ x (t, z)) = fun z => famRep g J F u₂ x (t, z) :=
    funext fun z => hfam x t ht z
  show tanPart (fderiv ℝ (fun z => famRep g J F u₁ x (t, z)) (chart E x x))
      (-dtW (famRep g J F u₁ x) T (t, chart E x x)) =
    tanPart (fderiv ℝ (fun z => famRep g J F u₂ x (t, z)) (chart E x x))
      (-dtW (famRep g J F u₂ x) T (t, chart E x x))
  rw [hfun, hdt]

variable {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
variable (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
variable {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) {ψ : ℝ → M → M}
variable (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
  (Icc 0 T ×ˢ univ))
variable (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ t (ψ t x) = graphImm E J F (fnOf g (u t)) x)
include hu hΦ hψ hG

/-- The graph family in the chart at `α` is the flow in the chart at `β` composed with `ψ_t`. -/
theorem chartRep_graph_eventuallyEq {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α β : M} {w : E}
    (hw : w ∈ diffeoDom E (ψ t) α β) :
    chartRep E (graphImm E J F (fnOf g (u t))) α =ᶠ[𝓝 w]
      chartRep E (Φ t) β ∘ diffeoRep E (ψ t) α β := by
  filter_upwards [(isOpen_diffeoDom' (contMDiff_slice hψ ht).continuous α β).mem_nhds hw] with w' hw'
  have h := chartRep_comp_diffeo (Φ t) hw'
  rw [Function.comp_apply, ← h]
  simp only [chartRep, Function.comp_apply]
  exact (hG t ht _).symm

/-- `ψ_t` is a local diffeomorphism: its composition with an immersion is an immersion. -/
theorem exists_equiv_fderiv_psi {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α β : M} {w : E}
    (hw : w ∈ diffeoDom E (ψ t) α β) :
    ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (diffeoRep E (ψ t) α β) w := by
  have hψt := contMDiff_slice hψ ht
  have hρ : DifferentiableAt ℝ (diffeoRep E (ψ t) α β) w :=
    ((contDiffOn_diffeoRep' hψt α β).contDiffAt
      ((isOpen_diffeoDom' hψt.continuous α β).mem_nhds hw)).differentiableAt (by simp)
  have hGt := isLagImmersion_solution g hJi hJJ hF hu ht
  have hmem : diffeoRep E (ψ t) α β w ∈ (chart E β).target := (chart E β).map_source hw.2
  have hX : DifferentiableAt ℝ (chartRep E (Φ t) β) (diffeoRep E (ψ t) α β w) :=
    (((hΦ.imm t ht).smooth β).contDiffAt
      ((isOpen_chart_target (E := E) β).mem_nhds hmem)).differentiableAt (by simp)
  have heq := chartRep_graph_eventuallyEq g hJi hJJ hF hu hΦ hψ hG ht hw
  have hcomp := fderiv_comp w hX hρ
  have hinj : Function.Injective (fderiv ℝ (diffeoRep E (ψ t) α β) w) := by
    intro a b hab
    apply hGt.inj α w hw.1
    rw [heq.fderiv_eq, hcomp, ContinuousLinearMap.comp_apply, ContinuousLinearMap.comp_apply, hab]
  exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (diffeoRep E (ψ t) α β) w : E →ₗ[ℝ] E) hinj
    rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩

include hdim

/-- **The transport equation** `∂ₜψ + Dψ(X) = 0`, `X = tanField`, in charts: the derivative of `ψ` in the
direction `(1, X)` vanishes. -/
theorem fderivWithin_psiRep_eRep {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α y β : M} (hα : y ∈ patchU E α)
    (hgα : good g J F α y (covJet g 0 0 (2 * 2) (u t) y)) (hβ : ψ t y ∈ (chart E β).source) :
    fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)
      (1, eRep g J F u T α (t, chart E α y)) = 0 := by
  have hzt : chart E α y ∈ (chart E α).target := (chart E α).map_source hα.1
  have hsymm : (chart E α).symm (chart E α y) = y := (chart E α).left_inv hα.1
  have hβz : ψ t ((chart E α).symm (chart E α y)) ∈ (chart E β).source := by rw [hsymm]; exact hβ
  have hceq : psiRep E ψ α β (t, chart E α y) = chart E β (ψ t y) := by
    simp only [psiRep]
    rw [hsymm]
  have hct : chart E β (ψ t y) ∈ (chart E β).target := (chart E β).map_source hβ
  have hP := contDiffWithinAt_psiRep hψ ht hzt hβz
  have hUD : UniqueDiffWithinAt ℝ (Icc (0 : ℝ) T ×ˢ (chart E α).target) (t, chart E α y) :=
    ((uniqueDiffOn_Icc hu.pos).prod (isOpen_chart_target (E := E) α).uniqueDiffOn) _ ⟨ht, hzt⟩
  have hPd : HasFDerivWithinAt (psiRep E ψ α β)
      (fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y))
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y) :=
    (hP.differentiableWithinAt (by simp)).hasFDerivWithinAt
  have hN : {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E β).source} ∈
      𝓝[Icc 0 T ×ˢ (chart E α).target] (t, chart E α y) :=
    (contMDiffOn_psiPre hψ α (t, chart E α y) ⟨ht, hzt⟩).continuousWithinAt.preimage_mem_nhdsWithin
      ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) β).mem_nhds hβz)
  have hΦrep := contDiffOn_flowRep hΦ.smooth β
  have hΦd := hasFDerivWithinAt_prod hu.pos (isOpen_chart_target (E := E) β) hΦrep ht hct
  rw [← hceq] at hΦd
  have hκ : HasFDerivWithinAt (fun q : ℝ × E => (q.1, psiRep E ψ α β q))
      ((ContinuousLinearMap.fst ℝ ℝ E).prod
        (fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)))
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y) := hasFDerivWithinAt_fst.prodMk hPd
  have hmaps : MapsTo (fun q : ℝ × E => (q.1, psiRep E ψ α β q))
      ((Icc 0 T ×ˢ (chart E α).target) ∩
        {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E β).source})
      (Icc 0 T ×ˢ (chart E β).target) := fun q hq => ⟨hq.1.1, (chart E β).map_source hq.2⟩
  have h1 := HasFDerivWithinAt.comp (g := flowRep E Φ β)
    (f := fun q : ℝ × E => (q.1, psiRep E ψ α β q)) (t, chart E α y) hΦd
    (hκ.mono inter_subset_left) hmaps
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin hN)
  have hfam := contDiffOn_famRep g hJi hJJ hF hu α
  have hfd := hasFDerivWithinAt_prod hu.pos (isOpen_chart_target (E := E) α) hfam ht hzt
  have hev : famRep g J F u α =ᶠ[𝓝[Icc 0 T ×ˢ (chart E α).target] (t, chart E α y)]
      flowRep E Φ β ∘ fun q : ℝ × E => (q.1, psiRep E ψ α β q) := by
    filter_upwards [inter_mem self_mem_nhdsWithin hN] with q hq
    simp only [famRep, chartRep, flowRep, psiRep, Function.comp_apply]
    rw [(chart E β).left_inv hq.2]
    exact (hG q.1 hq.1.1 _).symm
  have h3 := h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨ht, hzt⟩)
  have hEq := hUD.eq hfd h3
  have hEq1 := congrArg (fun L : ℝ × E →L[ℝ] V => L (1, eRep g J F u T α (t, chart E α y))) hEq
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.prod_apply,
    ContinuousLinearMap.coe_fst'] at hEq1
  rw [fderivWithin_prod_apply hu.pos (isOpen_chart_target (E := E) α) hfam ht hzt, hceq,
    fderivWithin_prod_apply hu.pos (isOpen_chart_target (E := E) β) hΦrep ht hct, one_smul,
    one_smul] at hEq1
  -- the left side is the flow velocity of the graph
  have hL : dtW (famRep g J F u α) T (t, chart E α y) +
      fderiv ℝ (fun w => famRep g J F u α (t, w)) (chart E α y)
        (eRep g J F u T α (t, chart E α y)) =
      flowVel J (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α y) :=
    dtW_add_fderiv_eRep g hJi hJJ hF hdim hu ht hα hgα
  -- the time derivative of the flow
  have hdt : dtW (flowRep E Φ β) T (t, chart E β (ψ t y)) =
      flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)) := by
    have hfun : (fun s => flowRep E Φ β (s, chart E β (ψ t y))) = fun s => Φ s (ψ t y) := by
      funext s
      simp only [flowRep]
      rw [(chart E β).left_inv hβ]
    show derivWithin (fun s => flowRep E Φ β (s, chart E β (ψ t y))) (Icc 0 T) t = _
    rw [hfun]
    exact (hΦ.vel t ht (ψ t y) β hβ).derivWithin (uniqueDiffOn_Icc hu.pos t ht)
  -- invariance of the velocity
  have hψt := contMDiff_slice hψ ht
  have hmem : chart E α y ∈ diffeoDom E (ψ t) α β := mem_diffeoDom hα.1 hβ
  have hfv : flowVel J (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α y) =
      flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)) := by
    have h := flowVel_comp' (isOpen_chart_target (E := E) β) ((hΦ.imm t ht).smooth β)
      ((hΦ.imm t ht).inj β) J (isOpen_diffeoDom' hψt.continuous α β) (contDiffOn_diffeoRep' hψt α β)
      (fun w hw => exists_equiv_fderiv_psi g hJi hJJ hF hu hΦ hψ hG ht hw)
      (fun w hw => (chart E β).map_source hw.2) hmem
    rw [diffeoRep_chart hα.1] at h
    rw [flowVel_congr_of_eventuallyEq
      (chartRep_graph_eventuallyEq g hJi hJJ hF hu hΦ hψ hG ht hmem), h]
  rw [hL, hfv, hdt] at hEq1
  have h0 : fderiv ℝ (chartRep E (Φ t) β) (chart E β (ψ t y))
      (fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)
        (1, eRep g J F u T α (t, chart E α y))) = 0 :=
    add_left_cancel (a := flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)))
      (by rw [add_zero]; exact hEq1.symm)
  exact (hΦ.imm t ht).inj β _ hct (h0.trans (map_zero _).symm)

/-- Along a curve of the field `tanField`, the chart readings of `ψ` have zero derivative. -/
theorem hasDerivWithinAt_chart_psi_chi {χ : ℝ → M → M} {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {x : M}
    (hχ : HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u T t (χ t x))))
    {β : M} (hβ : ψ t (χ t x) ∈ (chart E β).source) :
    HasDerivWithinAt (fun s => chart E β (ψ s (χ s x))) 0 (Icc 0 T) t := by
  obtain ⟨i, hi, hgi⟩ := exists_good g hJi hJJ hF hu ht (χ t x)
  have hs : χ t x ∈ (chartAt E i.1).source := by
    have h := hi.1
    rwa [extChartAt_source] at h
  have hch := (mdifferentiableAt_extChartAt (I := 𝓘(ℝ, E)) hs).hasMFDerivAt
  have hcomp := hch.comp t hχ
  have hfd : HasFDerivAt (fun s => chart E i.1 (χ s x)) (_ : ℝ →L[ℝ] E) t :=
    hasMFDerivAt_iff_hasFDerivAt.mp hcomp
  have hread := trivAt_tanField g hJi hJJ hF hu ht hi.1
  rw [TangentBundle.continuousLinearMapAt_trivializationAt hs] at hread
  have hd : HasDerivAt (fun s => chart E i.1 (χ s x))
      (eRep g J F u T i.1 (t, chart E i.1 (χ t x))) t := by
    refine hfd.hasDerivAt.congr_deriv ?_
    rw [← hread]
    exact congrArg _ (one_smul ℝ _)
  have hc : HasDerivAt (fun s => (s, chart E i.1 (χ s x)))
      ((1 : ℝ), eRep g J F u T i.1 (t, chart E i.1 (χ t x))) t := (hasDerivAt_id t).prodMk hd
  have hN : {s | χ s x ∈ (chart E i.1).source} ∈ 𝓝 t :=
    hχ.continuousAt.preimage_mem_nhds ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) i.1).mem_nhds hi.1)
  have hzt : chart E i.1 (χ t x) ∈ (chart E i.1).target := (chart E i.1).map_source hi.1
  have hβz : ψ t ((chart E i.1).symm (chart E i.1 (χ t x))) ∈ (chart E β).source := by
    rw [(chart E i.1).left_inv hi.1]
    exact hβ
  have hP := contDiffWithinAt_psiRep hψ (α := i.1) (β := β) ht hzt hβz
  have hPd := (hP.differentiableWithinAt (by simp)).hasFDerivWithinAt
  have hm : MapsTo (fun s => (s, chart E i.1 (χ s x)))
      (Icc 0 T ∩ {s | χ s x ∈ (chart E i.1).source}) (Icc 0 T ×ˢ (chart E i.1).target) :=
    fun s hs => ⟨hs.1, (chart E i.1).map_source hs.2⟩
  have h1 := HasFDerivWithinAt.comp_hasDerivWithinAt (f := fun s => (s, chart E i.1 (χ s x))) t hPd
    (hc.hasDerivWithinAt (s := Icc 0 T ∩ {s | χ s x ∈ (chart E i.1).source})) hm
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin (nhdsWithin_le_nhds hN))
  have hev : (fun s => chart E β (ψ s (χ s x))) =ᶠ[𝓝[Icc 0 T] t]
      psiRep E ψ i.1 β ∘ fun s => (s, chart E i.1 (χ s x)) := by
    filter_upwards [nhdsWithin_le_nhds hN] with s hs
    simp only [psiRep, Function.comp_apply]
    rw [(chart E i.1).left_inv hs]
  refine (h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ht)).congr_deriv ?_
  exact fderivWithin_psiRep_eRep g hJi hJJ hF hdim hu hΦ hψ hG ht hi hgi hβ

/-- **`ψ_t` is the inverse of the reparametrization `χ_t`**: `ψ_t ∘ χ_t = id` for every flow `χ` of the field
`tanField` starting at the identity. -/
theorem psi_comp_chi {χ : ℝ → M → M} (hψ0 : ∀ x, ψ 0 x = x) (hχ0 : ∀ x, χ 0 x = x)
    (hχ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u T t (χ t x))))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) : ψ t (χ t x) = x := by
  have hc : ContinuousOn (fun s => ψ s (χ s x)) (Icc 0 T) := by
    intro s hs
    have h1 : ContinuousWithinAt (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 T ×ˢ univ) (s, χ s x) :=
      (hψ (s, χ s x) ⟨hs, mem_univ _⟩).continuousWithinAt
    have h2 : ContinuousWithinAt (fun r => ((r, χ r x) : ℝ × M)) (Icc 0 T) s :=
      (continuousAt_id.prodMk (hχ s hs x).continuousAt).continuousWithinAt
    exact ContinuousWithinAt.comp (g := fun q : ℝ × M => ψ q.1 q.2)
      (f := fun r => ((r, χ r x) : ℝ × M)) h1 h2 (fun r hr => ⟨hr, mem_univ _⟩)
  have h := eq_of_chart_hasDerivWithinAt_zero (E := E) hc (fun s hs β hβ =>
    hasDerivWithinAt_chart_psi_chi g hJi hJJ hF hdim hu hΦ hψ hG hs (hχ s hs x) hβ) t ht
  have h' : ψ t (χ t x) = ψ 0 (χ 0 x) := h
  rw [h', hχ0, hψ0]

/-- **A flow that is a graph up to reparametrization is the constructed solution** `G_t ∘ χ_t`. -/
theorem flow_eq_graph_reparam {χ : ℝ → M → M} (hψ0 : ∀ x, ψ 0 x = x) (hχ0 : ∀ x, χ 0 x = x)
    (hχ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u T t (χ t x))))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) :
    Φ t x = graphImm E J F (fnOf g (u t)) (χ t x) := by
  rw [← hG t ht (χ t x), psi_comp_chi g hJi hJJ hF hdim hu hΦ hψ hG hψ0 hχ0 hχ ht x]

omit hu hΦ hψ hG

/-- **Uniqueness of the flow among graphs over `F`.** Two solutions of the flow that, after reparametrizations
`ψ₁`, `ψ₂` starting at the identity, are the graphs of solutions `u₁`, `u₂` of the scalar equation with the same
initial datum, coincide. -/
theorem flow_unique_of_graph {u₁ u₂ : ℝ → SmoothCcTensor g 0 0}
    (hu₁ : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u₁)
    (hu₂ : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u₂)
    {Φ₁ Φ₂ : ℝ → M → V} (hΦ₁ : IsFlow E J T Φ₁) (hΦ₂ : IsFlow E J T Φ₂) {ψ₁ ψ₂ : ℝ → M → M}
    (hψ₁ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ₁ q.1 q.2)
      (Icc 0 T ×ˢ univ))
    (hψ₂ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ₂ q.1 q.2)
      (Icc 0 T ×ˢ univ))
    (h₁0 : ∀ x, ψ₁ 0 x = x) (h₂0 : ∀ x, ψ₂ 0 x = x)
    (hG₁ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ₁ t (ψ₁ t x) = graphImm E J F (fnOf g (u₁ t)) x)
    (hG₂ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ₂ t (ψ₂ t x) = graphImm E J F (fnOf g (u₂ t)) x) :
    ∀ t ∈ Icc (0 : ℝ) T, Φ₁ t = Φ₂ t := by
  have huu := scalar_solution_unique g hJi hJJ hF u₀ T u₁ u₂ hu₁ hu₂
  obtain ⟨Y, hYsm, hYeq⟩ := DifferentialGeometry.Analysis.ODE.seeley_time_extend (I := 𝓘(ℝ, E))
    (tanField g J F u₁ T) T hu₁.pos (contMDiffOn_tanField g hJi hJJ hF hu₁)
  obtain ⟨χ, h0, -, -, hode⟩ := exists_reparam_flow Y hYsm hu₁.pos.le
  have hχ₁ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u₁ T t (χ t x))) := fun t ht x => by
    have h := hode t (Icc_subset_Ioo_reparam ht) x
    rwa [hYeq t ht] at h
  have hχ₂ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u₂ T t (χ t x))) := fun t ht x => by
    have h := hχ₁ t ht x
    rwa [tanField_congr g huu ht] at h
  intro t ht
  funext x
  rw [flow_eq_graph_reparam g hJi hJJ hF hdim hu₁ hΦ₁ hψ₁ hG₁ h₁0 h0 hχ₁ ht x,
    flow_eq_graph_reparam g hJi hJJ hF hdim hu₂ hΦ₂ hψ₂ hG₂ h₂0 h0 hχ₂ ht x, huu t ht]

/-- info: 'CWFlow.Bridge.flow_unique_of_graph' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flow_unique_of_graph

end Bridge

end CWFlow
