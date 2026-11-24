import CWFlow.Extension.FlowM
import CWFlow.Extension.ImmPot
import CWFlow.Flow.GraphUnique
import CWFlow.Scalar.LinearChange

/-!
# Bridge between the two formalisms for flows of Lagrangian immersions

The old formalism (`CWFlow.IsLagImmersion`, `CWFlow.IsFlow`: charts `extChartAt 𝓘(ℝ, E) z` with values in `E`,
closed time interval, basis-free velocity `flowVel`) and the new one (`IsLagImm`, `IsFlowM`: charts read on
`ℝⁿ` through a fixed `ι : ℝⁿ ≃L[ℝ] E`, open time set, coordinate velocity `velX`).

* `isLagImm_iff_isLagImmersion`: the two immersion notions agree.
* `velX_chartInv_eq_flowVel`: `velX = flowVel` (through `flowVel_eq_velX` and a local-in-space version of
  the linear invariance `flowVel_comp_linearEquiv`).
* `IsFlow.isFlowM`: an `IsFlow` on `[0, T]` is an `IsFlowM` on `(0, T)`.
* `isFlow_of_chart`: the converse, from chartwise data on `[0, T]` (one-sided velocity).
-/

open Set Filter
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

/-! ### Local-in-space invariance of the velocity under a linear equivalence -/

section LinearLocal

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {E E' : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup E'] [NormedSpace ℝ E'] [FiniteDimensional ℝ E']

theorem fb_meanCurv_comp (A : E' ≃L[ℝ] E) {X : E → V} {S : Set E} (hS : IsOpen S)
    (hX : ContDiffOn ℝ ∞ X S) (hinj : ∀ y ∈ S, Function.Injective (fderiv ℝ X y))
    {y : E'} (hy : A y ∈ S) : meanCurv (X ∘ A) y = meanCurv X (A y) := by
  have hD : DifferentiableAt ℝ (fderiv ℝ X) (A y) :=
    ((hX.contDiffAt (hS.mem_nhds hy)).fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)
  rw [meanCurv, meanCurv, fderiv_comp_linearEquiv A X y]
  exact HOf_comp_linearEquiv (hinj _ hy) A _ _ fun a c => fderiv_fderiv_comp_linearEquiv A hD a c

theorem fb_flowPot_comp (A : E' ≃L[ℝ] E) {X : E → V} {S : Set E} (hS : IsOpen S)
    (hX : ContDiffOn ℝ ∞ X S) (hinj : ∀ y ∈ S, Function.Injective (fderiv ℝ X y))
    (J : V →L[ℝ] V) {y : E'} (hy : A y ∈ S) : flowPot J (X ∘ A) y = flowPot J X (A y) := by
  have hL := fderiv_comp_linearEquiv A X y
  have hev : (fun w => J (meanCurv (X ∘ A) w)) =ᶠ[𝓝 y] (fun v => J (meanCurv X v)) ∘ A := by
    filter_upwards [A.continuous.continuousAt.preimage_mem_nhds (hS.mem_nhds hy)] with w hw
    simp only [Function.comp_apply, fb_meanCurv_comp A hS hX hinj hw]
  have hM : fderiv ℝ (fun w => J (meanCurv (X ∘ A) w)) y =
      (fderiv ℝ (fun v => J (meanCurv X v)) (A y)).comp (A : E' →L[ℝ] E) := by
    rw [hev.fderiv_eq, fderiv_comp y (differentiableAt_J_meanCurv J
      (hX.contDiffAt (hS.mem_nhds hy)) (hinj _ hy)) A.differentiableAt, A.fderiv]
  set M := fderiv ℝ (fun v => J (meanCurv X v)) (A y)
  set L' := fderiv ℝ X (A y)
  let B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ := LinearMap.mk₂ ℝ (fun a c => ⟪M a, L' c⟫)
    (fun a a' c => by simp [inner_add_left]) (fun r a c => by simp [inner_smul_left])
    (fun a c c' => by simp [inner_add_right]) (fun r a c => by simp [inner_smul_right])
  have h := trOf_comp_linearEquiv (hinj (A y) hy) A B
  simp only [flowPot, hL, hM]
  exact h

theorem fb_flowVel_comp (A : E' ≃L[ℝ] E) {X : E → V} {S : Set E} (hS : IsOpen S)
    (hX : ContDiffOn ℝ ∞ X S) (hinj : ∀ y ∈ S, Function.Injective (fderiv ℝ X y))
    (J : V →L[ℝ] V) {y : E'} (hy : A y ∈ S) : flowVel J (X ∘ A) y = flowVel J X (A y) := by
  have hpot : DifferentiableAt ℝ (flowPot J X) (A y) :=
    ((contDiffOn_flowPot hS hX hinj J).contDiffAt (hS.mem_nhds hy)).differentiableAt (by simp)
  have hL := fderiv_comp_linearEquiv A X y
  have hD : fderiv ℝ (flowPot J (X ∘ A)) y =
      (fderiv ℝ (flowPot J X) (A y)).comp (A : E' →L[ℝ] E) := by
    have hev : flowPot J (X ∘ A) =ᶠ[𝓝 y] flowPot J X ∘ A := by
      filter_upwards [A.continuous.continuousAt.preimage_mem_nhds (hS.mem_nhds hy)] with w hw
      exact fb_flowPot_comp A hS hX hinj J hw
    rw [hev.fderiv_eq, fderiv_comp y hpot A.differentiableAt, A.fderiv]
  simp only [flowVel, hL, hD, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [apply_ginvOf_comp_linearEquiv (hinj _ hy) A]
  congr 3
  ext a
  simp

end LinearLocal

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

/-! ### Composition with `ι` -/

theorem fb_isOpen_chartDomE (z : M) : IsOpen (chartDomE ι z) :=
  (isOpen_chart_target (E := E) z).preimage ι.continuous

theorem fb_contDiffOn_iff {X : E → V} {S : Set E} :
    ContDiffOn ℝ ∞ (X ∘ ι) (ι ⁻¹' S) ↔ ContDiffOn ℝ ∞ X S := by
  constructor
  · intro h
    have := h.comp (ι.symm.contDiff.contDiffOn (s := S)) (fun e he => by simpa using he)
    refine this.congr fun e _ => ?_
    simp
  · intro h
    exact h.comp ι.contDiff.contDiffOn (fun y hy => hy)

theorem fb_inj_iff {X : E → V} {y : EuclideanSpace ℝ (Fin n)} :
    Function.Injective (fderiv ℝ (X ∘ ι) y) ↔ Function.Injective (fderiv ℝ X (ι y)) := by
  rw [fderiv_comp_linearEquiv]
  constructor
  · intro h a b hab
    have h' : (fderiv ℝ X (ι y) ∘L (ι : EuclideanSpace ℝ (Fin n) →L[ℝ] E)) (ι.symm a) =
        (fderiv ℝ X (ι y) ∘L (ι : EuclideanSpace ℝ (Fin n) →L[ℝ] E)) (ι.symm b) := by
      simpa using hab
    simpa using congrArg ι (h h')
  · intro h a b hab
    exact ι.injective (h hab)

theorem fb_lag_iff {J : V →L[ℝ] V} {X : E → V} {y : EuclideanSpace ℝ (Fin n)} :
    (∀ u w, ⟪J (fderiv ℝ (X ∘ ι) y u), fderiv ℝ (X ∘ ι) y w⟫ = 0) ↔
      ∀ u w, ⟪J (fderiv ℝ X (ι y) u), fderiv ℝ X (ι y) w⟫ = 0 := by
  rw [fderiv_comp_linearEquiv]
  constructor
  · intro h u w
    simpa using h (ι.symm u) (ι.symm w)
  · intro h u w
    exact h _ _

theorem isLagImm_iff_isLagImmersion {J : V →L[ℝ] V} {F : M → V} :
    IsLagImm ι J F ↔ IsLagImmersion E J F := by
  have hc : ∀ z : M, F ∘ chartInv ι z = chartRep E F z ∘ ι := fun z => rfl
  constructor
  · intro h
    refine ⟨fun z => ?_, fun z y hy => ?_, fun z y hy u w => ?_⟩
    · have := h.smooth z
      rw [hc] at this
      exact (fb_contDiffOn_iff ι).1 this
    · have := h.inj z (ι.symm y) (by simpa [chartDomE] using hy)
      rw [hc, fb_inj_iff] at this
      simpa using this
    · have := h.lag z (ι.symm y) (by simpa [chartDomE] using hy)
      rw [hc, fb_lag_iff] at this
      simpa using this u w
  · intro h
    refine ⟨fun z => ?_, fun z y hy => ?_, fun z y hy u w => ?_⟩
    · rw [hc]
      exact (fb_contDiffOn_iff ι).2 (h.smooth z)
    · rw [hc, fb_inj_iff]
      exact h.inj z _ hy
    · rw [hc]
      exact (fb_lag_iff ι).2 (fun u w => h.lag z _ hy u w) u w

theorem velX_chartInv_eq_flowVel {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} (hF : IsLagImmersion E J F) (z : M) {y : EuclideanSpace ℝ (Fin n)}
    (hy : y ∈ chartDomE ι z) :
    velX J (F ∘ chartInv ι z) y = flowVel J (chartRep E F z) (ι y) := by
  have hI := (isLagImm_iff_isLagImmersion ι).2 hF
  have h1 := flowVel_eq_velX (J := J) (X := F ∘ chartInv ι z) hJi hJJ (fb_isOpen_chartDomE ι z)
    (hI.smooth z) (hI.inj z) (hI.lag z) hy
  rw [← h1]
  exact fb_flowVel_comp ι (isOpen_chart_target (E := E) z) (hF.smooth z) (hF.inj z) J hy

/-! ### Flows -/

theorem fb_chartInv_mem (z : M) {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z) :
    chartInv ι z y ∈ (chart E z).source ∧ chart E z (chartInv ι z y) = ι y :=
  ⟨(chart E z).map_target hy, (chart E z).right_inv hy⟩

theorem fb_contMDiffOn_chartPre (α : M) (T : ℝ) :
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

theorem IsFlow.isFlowM {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {T : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) : IsFlowM ι J (Ioo 0 T) Φ := by
  refine ⟨isOpen_Ioo, fun z => ?_, fun t ht => ?_, fun t ht z y hy => ?_⟩
  · have h3 := ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2) hΦ.smooth
      (fb_contMDiffOn_chartPre (E := E) z T) (fun p hp => ⟨hp.1, mem_univ _⟩)
    have h4 : ContDiffOn ℝ ∞ (fun p : ℝ × E => Φ p.1 ((chart E z).symm p.2))
        (Icc 0 T ×ˢ (chart E z).target) := contMDiffOn_iff_contDiffOn.mp h3
    have hmap : ContDiff ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => ((q.1, ι q.2) : ℝ × E)) :=
      contDiff_fst.prodMk (ι.contDiff.comp contDiff_snd)
    exact h4.comp hmap.contDiffOn fun q hq => ⟨Ioo_subset_Icc_self hq.1, hq.2⟩
  · exact (isLagImm_iff_isLagImmersion ι).2 (hΦ.imm t (Ioo_subset_Icc_self ht))
  · obtain ⟨hs, hc⟩ := fb_chartInv_mem ι z hy
    have h := hΦ.vel t (Ioo_subset_Icc_self ht) (chartInv ι z y) z hs
    rw [hc, ← velX_chartInv_eq_flowVel ι hJi hJJ (hΦ.imm t (Ioo_subset_Icc_self ht)) z hy] at h
    exact h.hasDerivAt (Icc_mem_nhds ht.1 ht.2)

theorem isFlow_of_chart {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V}
    (hs : ∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (chartInv ι z q.2))
      (Icc 0 T ×ˢ chartDomE ι z))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImm ι J (Φ t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ z : M, ∀ y ∈ chartDomE ι z,
      HasDerivWithinAt (fun s => Φ s (chartInv ι z y)) (velX J (Φ t ∘ chartInv ι z) y) (Icc 0 T) t) :
    IsFlow E J T Φ := by
  refine ⟨?_, fun t ht => (isLagImm_iff_isLagImmersion ι).1 (himm t ht), ?_⟩
  · rintro ⟨t, x⟩ hq
    have ht : t ∈ Icc (0 : ℝ) T := hq.1
    -- the family read in the chart at `x`, on `ℝ × E`
    have hmap : ContDiff ℝ ∞ (fun q : ℝ × E => ((q.1, ι.symm q.2) : ℝ × EuclideanSpace ℝ (Fin n))) :=
      contDiff_fst.prodMk (ι.symm.contDiff.comp contDiff_snd)
    have c2' : ContDiffOn ℝ ∞ (fun q : ℝ × E => Φ q.1 ((chart E x).symm q.2))
        (Icc 0 T ×ˢ (chart E x).target) := by
      refine ((hs x).comp hmap.contDiffOn (fun q hq => ⟨hq.1, ?_⟩)).congr fun q _ => ?_
      · simpa [chartDomE] using hq.2
      · simp [chartInv]
    have c2 : ContMDiffOn 𝓘(ℝ, ℝ × E) 𝓘(ℝ, V) ∞ (fun q : ℝ × E => Φ q.1 ((chart E x).symm q.2))
        (Icc 0 T ×ˢ (chart E x).target) := contMDiffOn_iff_contDiffOn.mpr c2'
    have c1 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
        (fun q : ℝ × M => ((q.1, chart E x q.2) : ℝ × E)) (Icc 0 T ×ˢ (chart E x).source) := by
      have ca : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => q.1)
          (Icc 0 T ×ˢ (chart E x).source) := contMDiffOn_fst
      have cb : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => chart E x q.2)
          (Icc 0 T ×ˢ (chart E x).source) :=
        ContMDiffOn.comp (g := chart E x) (f := fun q : ℝ × M => q.2)
          (t := (chart E x).source)
          (by rw [extChartAt_source]; exact contMDiffOn_extChartAt) contMDiffOn_snd
          (fun q hq => hq.2)
      exact ca.prodMk_space cb
    have hxs : x ∈ (chart E x).source := mem_extChartAt_source x
    have c3 := ContMDiffOn.comp (g := fun q : ℝ × E => Φ q.1 ((chart E x).symm q.2))
      (f := fun q : ℝ × M => ((q.1, chart E x q.2) : ℝ × E)) c2 c1
      (fun q hq => ⟨hq.1, (chart E x).map_source hq.2⟩)
    have c5 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞
        (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ (chart E x).source) (t, x) := by
      refine (c3 (t, x) ⟨ht, hxs⟩).congr (fun q hq' => ?_) ?_
      · show Φ q.1 q.2 = Φ q.1 ((chart E x).symm (chart E x q.2))
        rw [(chart E x).left_inv hq'.2]
      · show Φ t x = Φ t ((chart E x).symm (chart E x x))
        rw [(chart E x).left_inv hxs]
    exact c5.mono_of_mem_nhdsWithin (mem_nhdsWithin.2 ⟨univ ×ˢ (chart E x).source,
      isOpen_univ.prod (isOpen_extChartAt_source x), ⟨mem_univ _, hxs⟩,
      fun p hp => ⟨hp.2.1, hp.1.2⟩⟩)
  · intro t ht x z hx
    have hy : ι.symm (chart E z x) ∈ chartDomE ι z := by
      simpa [chartDomE] using (chart E z).map_source hx
    have hι : ι (ι.symm (chart E z x)) = chart E z x := by simp
    have hxx : chartInv ι z (ι.symm (chart E z x)) = x := by
      simp only [chartInv, hι]
      exact (chart E z).left_inv hx
    have h := hvel t ht z _ hy
    rw [hxx, velX_chartInv_eq_flowVel ι hJi hJJ ((isLagImm_iff_isLagImmersion ι).1 (himm t ht)) z hy,
      hι] at h
    exact h

/-- info: 'CWFlow.Extension.isLagImm_iff_isLagImmersion' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isLagImm_iff_isLagImmersion

/-- info: 'CWFlow.Extension.velX_chartInv_eq_flowVel' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms velX_chartInv_eq_flowVel

/-- info: 'CWFlow.Extension.IsFlow.isFlowM' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsFlow.isFlowM

/-- info: 'CWFlow.Extension.isFlow_of_chart' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isFlow_of_chart

end CWFlow.Extension
