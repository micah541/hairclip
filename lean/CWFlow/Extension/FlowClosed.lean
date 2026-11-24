import CWFlow.Extension.LimitScaled
import CWFlow.Extension.FlowBridge

/-!
# A flow with bounded second fundamental form extends to the closed interval

If `Φ` is a flow (`IsFlow E J T' Φ` on `[0, T']` for every `T' < T`) whose second fundamental form is
bounded on `[0, T)`, then it extends to a flow `Φbar` on `[0, T]` (`isFlow_closed_of_bound`), in the
original `CWFlow.IsFlow` format.

Proof: on `(0, T)` the flow is an `IsFlowM` (`IsFlow.isFlowM`, glued over `T'`);
`exists_limit_flow_of_bound` extends it to `[t₁, T]`; the two pieces (`Φ` on `[0, T']`, the extension on
`[t₁, T]`, with `t₁ < T' < T`) are glued in the chartwise format and `isFlow_of_chart` converts to `IsFlow`.
-/

open Set Filter
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem fc_mid {T t : ℝ} (ht : t ∈ Ioo 0 T) : (t + T) / 2 ∈ Ioo 0 T ∧ t < (t + T) / 2 := by
  obtain ⟨h1, h2⟩ := ht
  refine ⟨⟨by linarith, by linarith⟩, by linarith⟩

theorem fc_contDiffOn_chart_of_isFlow {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V}
    (hΦ : IsFlow E J T Φ) (z : M) :
    ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (chartInv ι z q.2))
      (Icc 0 T ×ˢ chartDomE ι z) := by
  have h3 := ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2) hΦ.smooth
    (fb_contMDiffOn_chartPre (E := E) z T) (fun p hp => ⟨hp.1, mem_univ _⟩)
  have h4 : ContDiffOn ℝ ∞ (fun p : ℝ × E => Φ p.1 ((chart E z).symm p.2))
      (Icc 0 T ×ˢ (chart E z).target) := contMDiffOn_iff_contDiffOn.mp h3
  have hmap : ContDiff ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => ((q.1, ι q.2) : ℝ × E)) :=
    contDiff_fst.prodMk (ι.contDiff.comp contDiff_snd)
  exact h4.comp hmap.contDiffOn fun q hq => ⟨hq.1, hq.2⟩

theorem fc_vel_of_isFlow {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {T : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (z : M)
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z) :
    HasDerivWithinAt (fun s => Φ s (chartInv ι z y)) (velX J (Φ t ∘ chartInv ι z) y) (Icc 0 T) t := by
  obtain ⟨hs, hc⟩ := fb_chartInv_mem ι z hy
  have h := hΦ.vel t ht (chartInv ι z y) z hs
  rw [hc, ← velX_chartInv_eq_flowVel ι hJi hJJ (hΦ.imm t ht) z hy] at h
  exact h

theorem fc_prod_mem {X : Type*} [TopologicalSpace X] {A B : Set ℝ} {D : Set X} {t : ℝ} {y : X}
    (h : B ∈ 𝓝[A] t) : B ×ˢ D ∈ 𝓝[A ×ˢ D] (t, y) := by
  rw [nhdsWithin_prod_eq]
  exact prod_mem_prod h self_mem_nhdsWithin

theorem fc_Icc_mem_left {T T' t : ℝ} (ht : t < T') : Icc 0 T' ∈ 𝓝[Icc 0 T] t :=
  mem_nhdsWithin.2 ⟨Iio T', isOpen_Iio, ht, fun _ hs => ⟨hs.2.1, hs.1.le⟩⟩

theorem fc_Icc_mem_right {T t₁ t : ℝ} (ht : t₁ < t) : Icc t₁ T ∈ 𝓝[Icc 0 T] t :=
  mem_nhdsWithin.2 ⟨Ioi t₁, isOpen_Ioi, ht, fun _ hs => ⟨hs.1.le, hs.2.2⟩⟩

theorem isFlow_closed_of_bound {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E J T' Φ) {K : ℝ}
    (hK : ∀ t ∈ Ico 0 T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ K ^ 2) :
    ∃ Φbar : ℝ → M → V, IsFlow E J T Φbar ∧ ∀ t ∈ Ico 0 T, Φbar t = Φ t := by
  classical
  have hM : IsFlowM ι J (Ioo 0 T) Φ := by
    refine ⟨isOpen_Ioo, fun z => ?_, fun t ht => ?_, fun t ht z y hy => ?_⟩
    · rintro ⟨t, y⟩ ⟨ht, hy⟩
      obtain ⟨hm, hlt⟩ := fc_mid ht
      have h := (IsFlow.isFlowM ι hJi hJJ (hΦ _ hm)).smooth z
      have hmem : Ioo 0 ((t + T) / 2) ×ˢ chartDomE ι z ∈ 𝓝 (t, y) :=
        (isOpen_Ioo.prod (fb_isOpen_chartDomE ι z)).mem_nhds ⟨⟨ht.1, hlt⟩, hy⟩
      exact ((h.contDiffAt (by simpa using hmem)).contDiffWithinAt)
    · obtain ⟨hm, hlt⟩ := fc_mid ht
      exact (IsFlow.isFlowM ι hJi hJJ (hΦ _ hm)).imm t ⟨ht.1, hlt⟩
    · obtain ⟨hm, hlt⟩ := fc_mid ht
      exact (IsFlow.isFlowM ι hJi hJJ (hΦ _ hm)).vel t ⟨ht.1, hlt⟩ z y hy
  obtain ⟨t₁, ht₁, Φ₀, hag, hsm, hlag, hvel⟩ :=
    exists_limit_flow_of_bound ι hJi hJJ hdim hT hM (K := K)
      (fun t ht z y hy => hK t ⟨ht.1.le, ht.2⟩ z y hy)
  let Φbar : ℝ → M → V := fun t p => if t < T then Φ t p else Φ₀ t p
  have hEq : ∀ t ∈ Icc t₁ T, Φbar t = Φ₀ t := by
    intro t ht
    by_cases htT : t < T
    · funext p
      simp only [Φbar, if_pos htT]
      rw [hag t ⟨lt_of_lt_of_le ht₁.1 ht.1, htT⟩]
    · funext p
      simp only [Φbar, if_neg htT]
  have hLow : ∀ t, t < T → Φbar t = Φ t := fun t ht => by
    funext p
    simp only [Φbar, if_pos ht]
  set T' : ℝ := (t₁ + T) / 2 with hT'
  have hT'm : T' ∈ Ioo 0 T := (fc_mid ht₁).1
  have ht₁T' : t₁ < T' := (fc_mid ht₁).2
  have hF := hΦ T' hT'm
  have hT'T : T' < T := hT'm.2
  have hs : ∀ z : M, ContDiffOn ℝ ∞
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
      (Icc 0 T ×ˢ chartDomE ι z) := by
    intro z
    rintro ⟨t, y⟩ ⟨ht, hy⟩
    have glue : ∀ (g : ℝ × EuclideanSpace ℝ (Fin n) → V)
        (S : Set (ℝ × EuclideanSpace ℝ (Fin n))),
        ContDiffWithinAt ℝ ∞ g S (t, y) → S ∈ 𝓝[Icc 0 T ×ˢ chartDomE ι z] (t, y) →
        (∀ q ∈ S, Φbar q.1 (chartInv ι z q.2) = g q) →
        ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
          (Icc 0 T ×ˢ chartDomE ι z) (t, y) := by
      intro g S hg hS hfg
      refine (hg.mono_of_mem_nhdsWithin hS).congr_of_eventuallyEq ?_
        (hfg _ (mem_of_mem_nhdsWithin (show (t, y) ∈ Icc 0 T ×ˢ chartDomE ι z from ⟨ht, hy⟩) hS))
      filter_upwards [hS] with q hq using hfg q hq
    by_cases htT' : t < T'
    · refine glue _ (Icc 0 T' ×ˢ chartDomE ι z)
        (fc_contDiffOn_chart_of_isFlow ι hF z (t, y)
          (show (t, y) ∈ Icc 0 T' ×ˢ chartDomE ι z from ⟨⟨ht.1, htT'.le⟩, hy⟩))
        (fc_prod_mem (fc_Icc_mem_left htT')) ?_
      rintro ⟨s, w⟩ ⟨hs1, _⟩
      simp only [hLow s (lt_of_le_of_lt hs1.2 hT'T)]
    · have ht1 : t₁ < t := lt_of_lt_of_le ht₁T' (not_lt.1 htT')
      refine glue _ (Icc t₁ T ×ˢ chartDomE ι z)
        (hsm z (t, y) (show (t, y) ∈ Icc t₁ T ×ˢ chartDomE ι z from ⟨⟨ht1.le, ht.2⟩, hy⟩))
        (fc_prod_mem (fc_Icc_mem_right ht1)) ?_
      rintro ⟨s, w⟩ ⟨hs1, _⟩
      simp only [hEq s hs1]
  have himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImm ι J (Φbar t) := by
    intro t ht
    by_cases htT' : t < T'
    · rw [hLow t (lt_trans htT' hT'T)]
      exact (isLagImm_iff_isLagImmersion ι).2 (hF.imm t ⟨ht.1, htT'.le⟩)
    · have ht1 : t₁ ≤ t := (not_lt.1 htT').trans' ht₁T'.le
      rw [hEq t ⟨ht1, ht.2⟩]
      exact hlag t ⟨ht1, ht.2⟩
  have hvel' : ∀ t ∈ Icc (0 : ℝ) T, ∀ z : M, ∀ y ∈ chartDomE ι z,
      HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar t ∘ chartInv ι z) y)
        (Icc 0 T) t := by
    intro t ht z y hy
    by_cases htT' : t < T'
    · have h1 := (fc_vel_of_isFlow ι hJi hJJ hF ⟨ht.1, htT'.le⟩ z hy).mono_of_mem_nhdsWithin
        (fc_Icc_mem_left (T := T) htT')
      rw [hLow t (lt_trans htT' hT'T)]
      refine h1.congr_of_eventuallyEq ?_ (by rw [hLow t (lt_trans htT' hT'T)])
      filter_upwards [fc_Icc_mem_left (T := T) htT'] with s hs
      rw [hLow s (lt_of_le_of_lt hs.2 hT'T)]
    · have ht1 : t₁ < t := lt_of_lt_of_le ht₁T' (not_lt.1 htT')
      have h1 := (hvel t ⟨ht1.le, ht.2⟩ z y hy).mono_of_mem_nhdsWithin
        (fc_Icc_mem_right (T := T) ht1)
      rw [hEq t ⟨ht1.le, ht.2⟩]
      refine h1.congr_of_eventuallyEq ?_ (by rw [hEq t ⟨ht1.le, ht.2⟩])
      filter_upwards [fc_Icc_mem_right (T := T) ht1] with s hs
      rw [hEq s hs]
  exact ⟨Φbar, isFlow_of_chart ι hJi hJJ hT hs himm hvel', fun t ht => hLow t ht.2⟩

/-- info: 'CWFlow.Extension.isFlow_closed_of_bound' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isFlow_closed_of_bound

end CWFlow.Extension
