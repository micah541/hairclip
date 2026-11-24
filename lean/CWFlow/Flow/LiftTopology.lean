import Mathlib.Analysis.Calculus.MeanValue

/-!
# Two topological lemmas for lifting a family through a local diffeomorphism

* `eqOn_of_lift_of_locallyInjective`: two continuous lifts on `[a, b]` of the same path through a map `p` that is
  locally injective along the first lift, with the same starting point, coincide.
* `exists_pos_forall_mem_of_isCompact`: a tube lemma in time: if `f(0, x)` lies in an open set `O` for all `x` in
  a compact set `K`, then `f(t, x) ∈ O` for all `x ∈ K` and all small `t ≥ 0`.
-/

open Set Filter

open scoped Topology

namespace CWFlow

/-- **Uniqueness of lifts.** Let `γ₁, γ₂ : [a, b] → Q` be continuous with `p ∘ γ₁ = p ∘ γ₂` and `γ₁ a = γ₂ a`.
If `p` is injective near every point of `γ₁`, then `γ₁ = γ₂` on `[a, b]`. -/
theorem eqOn_of_lift_of_locallyInjective {Q B : Type*} [TopologicalSpace Q] [T2Space Q] {p : Q → B}
    {γ₁ γ₂ : ℝ → Q} {a b : ℝ} (h₁ : ContinuousOn γ₁ (Icc a b)) (h₂ : ContinuousOn γ₂ (Icc a b))
    (hinj : ∀ t ∈ Icc a b, ∃ N ∈ 𝓝 (γ₁ t), InjOn p N)
    (hp : ∀ t ∈ Icc a b, p (γ₁ t) = p (γ₂ t)) (ha : γ₁ a = γ₂ a) : EqOn γ₁ γ₂ (Icc a b) := by
  have hsub : Icc a b ⊆ {t | γ₁ t = γ₂ t} := by
    refine IsClosed.Icc_subset_of_forall_mem_nhdsWithin ?_ ha ?_
    · have hc : ContinuousOn (fun t => (γ₁ t, γ₂ t)) (Icc a b) := h₁.prodMk h₂
      have h := hc.preimage_isClosed_of_isClosed isClosed_Icc (isClosed_diagonal (X := Q))
      rw [inter_comm] at h
      exact h
    · rintro x ⟨hx, hxI⟩
      have hxI' : x ∈ Icc a b := ⟨hxI.1, hxI.2.le⟩
      obtain ⟨N, hN, hNinj⟩ := hinj x hxI'
      have hN₂ : N ∈ 𝓝 (γ₂ x) := by
        have hx' : γ₁ x = γ₂ x := hx
        rwa [← hx']
      have hI : Icc a b ∈ 𝓝[>] x :=
        mem_of_superset (Ioc_mem_nhdsGT hxI.2) fun t ht => ⟨hxI.1.trans ht.1.le, ht.2⟩
      have hle : 𝓝[>] x ≤ 𝓝[Icc a b] x := nhdsWithin_le_of_mem hI
      have e₁ : ∀ᶠ t in 𝓝[>] x, γ₁ t ∈ N := hle ((h₁ x hxI').preimage_mem_nhdsWithin hN)
      have e₂ : ∀ᶠ t in 𝓝[>] x, γ₂ t ∈ N := hle ((h₂ x hxI').preimage_mem_nhdsWithin hN₂)
      filter_upwards [e₁, e₂, hI] with t ht₁ ht₂ htI
      exact hNinj ht₁ ht₂ (hp t htI)
  exact fun t ht => hsub ht

/-- **A tube lemma in time.** -/
theorem exists_pos_forall_mem_of_isCompact {X Q : Type*} [TopologicalSpace X] [TopologicalSpace Q]
    {K : Set X} (hK : IsCompact K) {T : ℝ} (hT : 0 < T) {f : ℝ × X → Q}
    (hf : ContinuousOn f (Icc 0 T ×ˢ K)) {O : Set Q} (hO : IsOpen O) (h0 : ∀ x ∈ K, f (0, x) ∈ O) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∀ t ∈ Icc (0 : ℝ) δ, ∀ x ∈ K, f (t, x) ∈ O := by
  have hev : ∀ᶠ t in 𝓝 (0 : ℝ), ∀ x ∈ K, x ∈ K → t ∈ Icc (0 : ℝ) T → f (t, x) ∈ O := by
    refine hK.eventually_forall_of_forall_eventually fun x hx => ?_
    have hmem : f ⁻¹' O ∈ 𝓝[Icc 0 T ×ˢ K] ((0 : ℝ), x) :=
      (hf (0, x) ⟨left_mem_Icc.mpr hT.le, hx⟩).preimage_mem_nhdsWithin (hO.mem_nhds (h0 x hx))
    rw [mem_nhdsWithin_iff_exists_mem_nhds_inter] at hmem
    obtain ⟨U, hU, hUs⟩ := hmem
    filter_upwards [hU] with z hz
    intro hzK hzI
    exact hUs ⟨hz, hzI, hzK⟩
  obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.mp hev
  refine ⟨min (ε / 2) T, lt_min (by positivity) hT, min_le_right _ _, fun t ht x hx => ?_⟩
  have htε : t ∈ Metric.ball (0 : ℝ) ε := by
    rw [Real.ball_eq_Ioo]
    have h1 : t ≤ ε / 2 := ht.2.trans (min_le_left _ _)
    constructor <;> linarith [ht.1]
  exact hball htε x hx hx ⟨ht.1, ht.2.trans (min_le_right _ _)⟩

end CWFlow

namespace CWFlow

/-- **Uniqueness of continuous solution curves** under a local uniqueness hypothesis: if near every point of
`γ₁` the two curves agree whenever both are near that point, and `γ₁ a = γ₂ a`, then `γ₁ = γ₂` on `[a, b]`. -/
theorem eqOn_of_locally_unique {Q : Type*} [TopologicalSpace Q] [T2Space Q] {γ₁ γ₂ : ℝ → Q} {a b : ℝ}
    (h₁ : ContinuousOn γ₁ (Icc a b)) (h₂ : ContinuousOn γ₂ (Icc a b))
    (huniq : ∀ t ∈ Icc a b, ∃ N ∈ 𝓝 (γ₁ t), ∀ s ∈ Icc a b, γ₁ s ∈ N → γ₂ s ∈ N → γ₁ s = γ₂ s)
    (ha : γ₁ a = γ₂ a) : EqOn γ₁ γ₂ (Icc a b) := by
  have hsub : Icc a b ⊆ {t | γ₁ t = γ₂ t} := by
    refine IsClosed.Icc_subset_of_forall_mem_nhdsWithin ?_ ha ?_
    · have hc : ContinuousOn (fun t => (γ₁ t, γ₂ t)) (Icc a b) := h₁.prodMk h₂
      have h := hc.preimage_isClosed_of_isClosed isClosed_Icc (isClosed_diagonal (X := Q))
      rw [inter_comm] at h
      exact h
    · rintro x ⟨hx, hxI⟩
      have hxI' : x ∈ Icc a b := ⟨hxI.1, hxI.2.le⟩
      obtain ⟨N, hN, hNu⟩ := huniq x hxI'
      have hN₂ : N ∈ 𝓝 (γ₂ x) := by
        have hx' : γ₁ x = γ₂ x := hx
        rwa [← hx']
      have hI : Icc a b ∈ 𝓝[>] x :=
        mem_of_superset (Ioc_mem_nhdsGT hxI.2) fun t ht => ⟨hxI.1.trans ht.1.le, ht.2⟩
      have hle : 𝓝[>] x ≤ 𝓝[Icc a b] x := nhdsWithin_le_of_mem hI
      have e₁ : ∀ᶠ t in 𝓝[>] x, γ₁ t ∈ N := hle ((h₁ x hxI').preimage_mem_nhdsWithin hN)
      have e₂ : ∀ᶠ t in 𝓝[>] x, γ₂ t ∈ N := hle ((h₂ x hxI').preimage_mem_nhdsWithin hN₂)
      filter_upwards [e₁, e₂, hI] with t ht₁ ht₂ htI
      exact hNu t htI ht₁ ht₂
  exact fun t ht => hsub ht

end CWFlow
