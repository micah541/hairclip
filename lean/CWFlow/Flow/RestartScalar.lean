import CWFlow.Flow.FlowExistence
import CWFlow.Bridge.ScalarTensor
import CWFlow.Flow.GraphUnique
import ShortTime.Vendor.ParametricJets
import ShortTime.Nonlinear.ContinuousDependenceCov

/-!
# Restarting a flow: the scalar side

* `exists_flow_of_solution`: `exists_flow` for an arbitrary datum `u₀` of the scalar equation. A solution `u` of
  the scalar equation on `[0, T]` with datum `u₀` gives a flow of Lagrangian immersions on `[0, T]` starting from
  the graph `graphImm (fnOf g u₀)`, of the form `t ↦ G_t ∘ χ_t` with `G_t = graphImm (fnOf g (u t))`.
* `tendsto_scalarFam_zero`: a family `w_s` of functions, jointly smooth on `[0, δ] × M` and vanishing at `s = 0`,
  tends to `0` as `s → 0⁺` in every covariant `C^j` norm (the form of the hypothesis of the continuous dependence
  part of Mantegazza–Martinazzi, `ShortTime.mantegazza_martinazzi`). This is the joint continuity of the fibre
  norm along a jointly smooth family (`ShortTime.Vendor.joint_jet_small`), applied to `s = t²` so that the base
  point `0` is interior.
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
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
variable (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)

include hdim

/-- Short-time existence of the flow from the graph of any datum of the scalar equation. -/
theorem exists_flow_of_solution {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
    (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u) :
    ∃ χ : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M), (∀ x, χ 0 x = x) ∧
      IsFlow E J T (fun t x => graphImm E J F (fnOf g (u t)) (χ t x)) := by
  have hT : 0 < T := hu.pos
  obtain ⟨Y, hYsm, hYeq⟩ := DifferentialGeometry.Analysis.ODE.seeley_time_extend (I := 𝓘(ℝ, E))
    (tanField g J F u T) T hT (contMDiffOn_tanField g hJi hJJ hF hu)
  obtain ⟨χ, h0, hχsm, -, hode⟩ := exists_reparam_flow Y hYsm hT.le
  refine ⟨χ, h0, ?_, ?_, ?_⟩
  · exact contMDiffOn_reparam g hJi hJJ hF hu (χ := fun t x => χ t x) hχsm
  · intro t ht
    exact isLagImmersion_comp_diffeo (isLagImmersion_solution g hJi hJJ hF hu ht) (χ t)
  · intro t ht x z hz
    obtain ⟨i, hi, hgi⟩ := exists_good g hJi hJJ hF hu ht (χ t x)
    have hχ := hode t (Icc_subset_Ioo_reparam ht) x
    rw [hYeq t ht] at hχ
    have h := hasDerivWithinAt_reparam g hJi hJJ hF hdim hu ht hχ hi hgi
    have hG := isLagImmersion_solution g hJi hJJ hF hu ht
    have hmem : chart E z x ∈ diffeoDom E (χ t) z i.1 := mem_diffeoDom hz hi.1
    have heq : chartRep E (fun x => graphImm E J F (fnOf g (u t)) (χ t x)) z =ᶠ[𝓝 (chart E z x)]
        chartRep E (graphImm E J F (fnOf g (u t))) i.1 ∘ diffeoRep E (χ t) z i.1 := by
      filter_upwards [(isOpen_diffeoDom (χ t) z i.1).mem_nhds hmem] with w hw
      exact chartRep_comp_diffeo (graphImm E J F (fnOf g (u t))) hw
    have hfv := flowVel_comp' (isOpen_chart_target (E := E) i.1) (hG.smooth i.1) (hG.inj i.1) J
      (isOpen_diffeoDom (χ t) z i.1) (contDiffOn_diffeoRep (χ t) z i.1)
      (fun w hw => exists_equiv_fderiv_diffeoRep (χ t) hw)
      (fun w hw => (chart E i.1).map_source hw.2) hmem
    rw [diffeoRep_chart hz] at hfv
    show HasDerivWithinAt _ (flowVel J (chartRep E (fun x => graphImm E J F (fnOf g (u t)) (χ t x)) z)
      (chart E z x)) _ _
    rw [flowVel_congr_of_eventuallyEq heq, hfv]
    exact h

/-- info: 'CWFlow.Bridge.exists_flow_of_solution' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_flow_of_solution

omit hdim

/-- `scalarFam` only depends on the value of the family at the time. -/
theorem rs_scalarFam_congr {f f' : ℝ → M → ℝ} {t t' : ℝ} (h : f t = f' t') :
    scalarFam g f t = scalarFam g f' t' := by
  unfold scalarFam
  rw [h]

/-- A jointly smooth family of functions vanishing at time `0` tends to `0` in every covariant `C^j` norm. -/
theorem tendsto_scalarFam_zero {w : ℝ → M → ℝ} {δ : ℝ} (hδ : 0 < δ)
    (hw : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => w q.1 q.2) (Icc 0 δ ×ˢ univ))
    (hw0 : ∀ x, w 0 x = 0) :
    ∀ j : ℕ, ∀ ε > 0, ∀ᶠ s in 𝓝[>] (0 : ℝ), ∀ x : M,
      fiberNorm g 0 (0 + j) x ((iteratedCovGrad g 0 0 j (scalarFam g w s)).toSection x) ≤ ε := by
  intro j ε hε
  obtain ⟨a, ha, haδ⟩ : ∃ a : ℝ, 0 < a ∧ a ^ 2 = δ :=
    ⟨Real.sqrt δ, Real.sqrt_pos.mpr hδ, Real.sq_sqrt hδ.le⟩
  obtain ⟨w', hw'def⟩ : ∃ w' : ℝ → M → ℝ, w' = fun t x => w (t ^ 2) x := ⟨_, rfl⟩
  have hw'a : ∀ t x, w' t x = w (t ^ 2) x := fun t x => by rw [hw'def]
  have hw' : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => w' q.1 q.2)
      (Icc (-a) a ×ˢ univ) := by
    have hf : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((q.1 ^ 2, q.2) : ℝ × M)) :=
      ((contDiff_id.pow 2).contMDiff.comp contMDiff_fst).prodMk contMDiff_snd
    refine (ContMDiffOn.comp (g := fun q : ℝ × M => w q.1 q.2) hw hf.contMDiffOn ?_).congr ?_
    · intro q hq
      refine ⟨⟨sq_nonneg _, ?_⟩, mem_univ _⟩
      rw [← haδ]
      nlinarith [mul_nonneg (sub_nonneg.2 hq.1.2) (show 0 ≤ q.1 + a by linarith [hq.1.1])]
    · intro q _
      exact hw'a q.1 q.2
  have hjoint := jointlySmoothOn_scalarFam g (f := w') (S := Icc (-a) a) hw'
  have hzero : scalarFam g w' 0 = 0 := by
    have hs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (w' 0) :=
      contMDiff_slice_fn hw' (t := 0) ⟨by linarith, ha.le⟩
    rw [scalarFam_eq g hs]
    apply SmoothCcTensor.ext
    apply ContMDiffSection.ext
    intro x
    change (scalarTensor g (w' 0) hs).toSection x = (0 : SmoothCcTensor g 0 0).toSection x
    have h0 : w' 0 x = 0 := by rw [hw'a]; simpa using hw0 x
    rw [scalarTensor_toSection, h0, zero_smul, SmoothCcTensor.toSection_zero]
    rfl
  have hS : Icc (-a) a ∈ 𝓝 (0 : ℝ) := Icc_mem_nhds (by linarith) ha
  have hε2 : 0 < ε ^ 2 := by positivity
  have h := ShortTime.Vendor.joint_jet_small (I := 𝓘(ℝ, E)) g 0 0 j (scalarFam g w') hS hzero hjoint hε2
  have hsq : Tendsto Real.sqrt (𝓝[>] (0 : ℝ)) (𝓝 0) :=
    (Real.continuous_sqrt.tendsto' 0 0 Real.sqrt_zero).mono_left nhdsWithin_le_nhds
  filter_upwards [hsq.eventually h, self_mem_nhdsWithin] with s hs hpos x
  have hpos' : (0 : ℝ) < s := hpos
  have hc : scalarFam g w' (Real.sqrt s) = scalarFam g w s :=
    rs_scalarFam_congr g (by funext y; rw [hw'a, Real.sq_sqrt hpos'.le])
  have h1 := hs j le_rfl x
  rw [hc] at h1
  exact (fiberNorm_le_iff hε.le).2 h1.le

/-- info: 'CWFlow.Bridge.tendsto_scalarFam_zero' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms tendsto_scalarFam_zero

end CWFlow.Bridge
