import CWFlow.Interp.RadialLift

/-!
# Lifting a ray under an a priori bound

Let `f` be a local homeomorphism on an open set `U` of a Hausdorff space with `f p = 0`.  If every partial lift
of the segment `s ↦ s • v` starting at `p` stays in a fixed compact subset `K ⊆ U`, then the whole segment lifts
(`exists_ray_lift`).  The proof is an open/closed argument on the supremum `b*` of the times up to which a lift
exists: lifts are unique (`eqOn_of_lift`), hence glue to a lift on `[0, b*)`; this stays in `K`, so it has a
cluster point in `U`, and a local lift around the cluster point extends it past `b*`.
-/

open Set Filter Metric

open scoped Topology

namespace CWFlow.Interp

variable {X E : Type*} [TopologicalSpace X] [T2Space X] [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- Two lifts of the ray on `[a, b]` that agree at `a` agree on `[a, b]`. -/
theorem rl_agree {f : X → E} {U : Set X} (hf : IsLocalHomeomorphOn f U) {v : E} {γ₁ γ₂ : ℝ → X} {a b : ℝ}
    (h₁ : ContinuousOn γ₁ (Icc a b)) (h₂ : ContinuousOn γ₂ (Icc a b)) (h0 : γ₁ a = γ₂ a)
    (l₁ : ∀ s ∈ Icc a b, γ₁ s ∈ U ∧ f (γ₁ s) = s • v) (l₂ : ∀ s ∈ Icc a b, f (γ₂ s) = s • v) :
    EqOn γ₁ γ₂ (Icc a b) :=
  eqOn_of_lift (p := f) h₁ h₂ (fun t ht => injOn_nhds_of_isLocalHomeomorphOn hf (l₁ t ht).1)
    (fun t ht => by rw [(l₁ t ht).2, l₂ t ht]) h0

/-- Local lifting of the ray around a point `q` over `s₀ • v`. -/
theorem rl_local {f : X → E} {U : Set X} (hU : IsOpen U) (hf : IsLocalHomeomorphOn f U) {q : X} (hq : q ∈ U)
    {v : E} {s₀ : ℝ} (hfq : f q = s₀ • v) :
    ∃ ε > 0, ∃ N ∈ 𝓝 q, ∃ δ : ℝ → X, ContinuousOn δ (Icc (s₀ - ε) (s₀ + ε)) ∧ δ s₀ = q ∧
      (∀ s ∈ Icc (s₀ - ε) (s₀ + ε), δ s ∈ U ∧ f (δ s) = s • v) ∧
      ∀ s ∈ Icc (s₀ - ε) (s₀ + ε), ∀ x ∈ N, f x = s • v → x = δ s := by
  obtain ⟨e, hqe, rfl⟩ := hf q hq
  have hs₀ : s₀ • v ∈ e.target := by
    rw [← hfq]; exact e.map_source hqe
  have hq' : e.symm (s₀ • v) = q := by
    rw [← hfq]; exact e.left_inv hqe
  have hcont : Continuous fun s : ℝ => s • v := continuous_id.smul continuous_const
  have hA : (fun s : ℝ => s • v) ⁻¹' e.target ∈ 𝓝 s₀ := hcont.continuousAt.preimage_mem_nhds
    (e.open_target.mem_nhds hs₀)
  have hcs : ContinuousAt (fun s : ℝ => e.symm (s • v)) s₀ :=
    (e.continuousAt_symm hs₀).comp (f := fun s : ℝ => s • v) hcont.continuousAt
  have hB : (fun s : ℝ => e.symm (s • v)) ⁻¹' U ∈ 𝓝 s₀ :=
    hcs.preimage_mem_nhds (hU.mem_nhds (by rw [hq']; exact hq))
  obtain ⟨ε, hε, hεW⟩ := Metric.mem_nhds_iff.mp (inter_mem hA hB)
  have hsub : ∀ s ∈ Icc (s₀ - ε / 2) (s₀ + ε / 2), s ∈ ball s₀ ε := by
    intro s hs
    rw [Metric.mem_ball, Real.dist_eq, abs_sub_lt_iff]
    constructor <;> linarith [hs.1, hs.2]
  refine ⟨ε / 2, half_pos hε, e.source, e.open_source.mem_nhds hqe, fun s => e.symm (s • v), ?_, ?_, ?_, ?_⟩
  · exact e.continuousOn_symm.comp hcont.continuousOn fun s hs => (hεW (hsub s hs)).1
  · simpa using hq'
  · intro s hs
    have h := hεW (hsub s hs)
    exact ⟨h.2, e.right_inv h.1⟩
  · intro s hs x hx hxs
    show x = e.symm (s • v)
    rw [← hxs]
    exact (e.left_inv hx).symm

/-- The extension step: a lift on `[0, b)` staying in a compact subset of `U` extends past `b`. -/
theorem rl_extend {f : X → E} {U : Set X} (hU : IsOpen U) (hf : IsLocalHomeomorphOn f U) {v : E} {K : Set X}
    (hK : IsCompact K) (hKU : K ⊆ U) {b : ℝ} (hb : 0 < b) {γ : ℝ → X} (hc : ContinuousOn γ (Ico 0 b))
    (hl : ∀ s ∈ Ico 0 b, γ s ∈ U ∧ f (γ s) = s • v) (hKγ : ∀ s ∈ Ico 0 b, γ s ∈ K) :
    ∃ ε > 0, ∃ γ' : ℝ → X, ContinuousOn γ' (Icc 0 (b + ε)) ∧ (∀ s ∈ Ico 0 b, γ' s = γ s) ∧
      ∀ s ∈ Icc 0 (b + ε), γ' s ∈ U ∧ f (γ' s) = s • v := by
  have hev : Ico 0 b ∈ 𝓝[<] b := Ico_mem_nhdsLT hb
  have hle : map γ (𝓝[<] b) ≤ 𝓟 K := by
    rw [le_principal_iff]
    exact Filter.mem_map.2 (mem_of_superset hev fun s hs => hKγ s hs)
  obtain ⟨q, hqK, hq⟩ := hK.exists_clusterPt hle
  have hqU : q ∈ U := hKU hqK
  have hfreq : ∀ N ∈ 𝓝 q, ∃ᶠ s in 𝓝[<] b, γ s ∈ N := fun N hN => mapClusterPt_iff_frequently.1 hq N hN
  have hfq : f q = b • v := by
    by_contra hne
    obtain ⟨A, B, hA, hB, hqA, hbB, hAB⟩ := t2_separation hne
    have hfc : ContinuousAt f q := hf.continuousAt hqU
    have h1 : ∃ᶠ s in 𝓝[<] b, γ s ∈ f ⁻¹' A :=
      hfreq (f ⁻¹' A) (hfc.preimage_mem_nhds (hA.mem_nhds hqA))
    have ht : Tendsto (fun s : ℝ => s • v) (𝓝 b) (𝓝 (b • v)) :=
      (continuous_id.smul continuous_const).tendsto b
    have h2 : ∀ᶠ s in 𝓝[<] b, s • v ∈ B :=
      (ht.mono_left nhdsWithin_le_nhds).eventually (hB.mem_nhds hbB)
    obtain ⟨s, hs1, hs2, hs3⟩ := (h1.and_eventually (h2.and hev)).exists
    have hs1' : f (γ s) ∈ A := hs1
    rw [(hl s hs3).2] at hs1'
    exact Set.disjoint_left.1 hAB hs1' hs2
  obtain ⟨ε, hε, N, hN, δ, hδc, hδb, hδl, hδN⟩ := rl_local hU hf hqU hfq
  have hIoo : Ioo (b - ε) b ∈ 𝓝[<] b := Ioo_mem_nhdsLT (by linarith)
  obtain ⟨s1, hs1N, hs1I, hs1K⟩ := ((hfreq N hN).and_eventually (inter_mem hIoo hev)).exists
  have hs1b : s1 ∈ Icc (b - ε) (b + ε) := ⟨hs1I.1.le, by linarith [hs1I.2]⟩
  have h0 : γ s1 = δ s1 := hδN s1 hs1b (γ s1) hs1N (hl s1 hs1K).2
  have hγδ : ∀ s ∈ Ico s1 b, γ s = δ s := by
    intro s hs
    have hsb : Icc s1 s ⊆ Ico 0 b := fun t ht => ⟨hs1K.1.trans ht.1, lt_of_le_of_lt ht.2 hs.2⟩
    have hsd : Icc s1 s ⊆ Icc (b - ε) (b + ε) := fun t ht =>
      ⟨hs1I.1.le.trans ht.1, by linarith [ht.2, hs.2]⟩
    exact rl_agree hf (hc.mono hsb) (hδc.mono hsd) h0 (fun t ht => hl t (hsb ht))
      (fun t ht => (hδl t (hsd ht)).2) ⟨hs.1, le_rfl⟩
  refine ⟨ε, hε, fun s => if s ≤ s1 then γ s else δ s, ?_, ?_, ?_⟩
  · have hsub : Icc (0 : ℝ) (b + ε) ⊆ Icc 0 s1 ∪ Icc s1 (b + ε) := by
      intro s hs
      by_cases h : s ≤ s1
      · exact Or.inl ⟨hs.1, h⟩
      · exact Or.inr ⟨(not_le.1 h).le, hs.2⟩
    refine ContinuousOn.mono (ContinuousOn.union_of_isClosed ?_ ?_ isClosed_Icc isClosed_Icc) hsub
    · refine (hc.mono fun t ht => ⟨ht.1, lt_of_le_of_lt ht.2 hs1K.2⟩).congr fun s hs => ?_
      simp only [if_pos hs.2]
    · refine (hδc.mono fun t ht => ⟨by linarith [ht.1, hs1I.1], ht.2⟩).congr fun s hs => ?_
      by_cases h : s ≤ s1
      · have : s = s1 := le_antisymm h hs.1
        subst this
        simp only [if_pos h]
        exact h0
      · simp only [if_neg h]
  · intro s hs
    by_cases h : s ≤ s1
    · simp only [if_pos h]
    · simp only [if_neg h]
      exact (hγδ s ⟨(not_le.1 h).le, hs.2⟩).symm
  · intro s hs
    by_cases h : s ≤ s1
    · simp only [if_pos h]
      exact hl s ⟨hs.1, lt_of_le_of_lt h hs1K.2⟩
    · simp only [if_neg h]
      exact hδl s ⟨by linarith [not_le.1 h, hs1I.1], hs.2⟩

/-- **Ray lifting under an a priori bound.** Let `f` be a local homeomorphism on the open set `U` of a Hausdorff
space, `f p = 0`. If every partial lift of the segment `s ↦ s • v` starting at `p` stays in a fixed compact subset
`K` of `U`, then the whole segment `s ∈ [0, 1]` lifts. -/
theorem exists_ray_lift {X E : Type*} [TopologicalSpace X] [T2Space X] [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : X → E} {U : Set X} {p : X} (hU : IsOpen U) (hf : IsLocalHomeomorphOn f U) (hp : p ∈ U)
    (hfp : f p = 0) {v : E} {K : Set X} (hK : IsCompact K) (hKU : K ⊆ U)
    (hapriori : ∀ b ∈ Ioc (0 : ℝ) 1, ∀ γ : ℝ → X, ContinuousOn γ (Ico 0 b) → γ 0 = p →
      (∀ s ∈ Ico 0 b, γ s ∈ U ∧ f (γ s) = s • v) → ∀ s ∈ Ico 0 b, γ s ∈ K) :
    ∃ γ : ℝ → X, ContinuousOn γ (Icc 0 1) ∧ γ 0 = p ∧
      ∀ s ∈ Icc (0 : ℝ) 1, γ s ∈ U ∧ f (γ s) = s • v := by
  classical
  have : Nonempty X := ⟨p⟩
  set S : Set ℝ := {b | b ∈ Icc (0 : ℝ) 1 ∧ ∃ γ : ℝ → X, ContinuousOn γ (Icc 0 b) ∧ γ 0 = p ∧
    ∀ s ∈ Icc (0 : ℝ) b, γ s ∈ U ∧ f (γ s) = s • v} with hSdef
  have h0S : (0 : ℝ) ∈ S := by
    refine ⟨⟨le_rfl, zero_le_one⟩, fun _ => p, continuousOn_const, rfl, fun s hs => ?_⟩
    have : s = 0 := le_antisymm hs.2 hs.1
    rw [this]
    exact ⟨hp, by simp [hfp]⟩
  have hSne : S.Nonempty := ⟨0, h0S⟩
  have hbdd : BddAbove S := ⟨1, fun b hb => hb.1.2⟩
  have hSmem : ∀ b ∈ S, ∃ γ : ℝ → X, ContinuousOn γ (Icc 0 b) ∧ γ 0 = p ∧
      ∀ s ∈ Icc (0 : ℝ) b, γ s ∈ U ∧ f (γ s) = s • v := fun b hb => hb.2
  choose! g hgc hg0 hgl using hSmem
  have hagree : ∀ b₁ ∈ S, ∀ b₂ ∈ S, ∀ s ∈ Icc (0 : ℝ) (min b₁ b₂), g b₁ s = g b₂ s := by
    intro b₁ h₁ b₂ h₂ s hs
    have m1 := min_le_left b₁ b₂
    have m2 := min_le_right b₁ b₂
    exact rl_agree hf ((hgc b₁ h₁).mono (Icc_subset_Icc_right m1))
      ((hgc b₂ h₂).mono (Icc_subset_Icc_right m2)) (by rw [hg0 b₁ h₁, hg0 b₂ h₂])
      (fun t ht => hgl b₁ h₁ t ⟨ht.1, ht.2.trans m1⟩) (fun t ht => (hgl b₂ h₂ t ⟨ht.1, ht.2.trans m2⟩).2) hs
  -- the glued lift
  let γ : ℝ → X := fun s => if h : ∃ b ∈ S, s ≤ b then g h.choose s else p
  have hγ : ∀ b ∈ S, ∀ s ∈ Icc (0 : ℝ) b, γ s = g b s := by
    intro b hb s hs
    have h : ∃ b ∈ S, s ≤ b := ⟨b, hb, hs.2⟩
    show (if h : ∃ b ∈ S, s ≤ b then g h.choose s else p) = g b s
    rw [dif_pos h]
    obtain ⟨hc1, hc2⟩ := h.choose_spec
    exact hagree _ hc1 b hb s ⟨hs.1, le_min hc2 hs.2⟩
  have hb0 : 0 ≤ sSup S := le_csSup hbdd h0S
  have hb1 : sSup S ≤ 1 := csSup_le hSne fun b hb => hb.1.2
  have hγ0 : γ 0 = p := by rw [hγ 0 h0S 0 ⟨le_rfl, le_rfl⟩, hg0 0 h0S]
  rcases eq_or_lt_of_le hb0 with hz | hpos
  · -- `sSup S = 0`: local lift at `p`
    exfalso
    obtain ⟨ε, hε, N, hN, δ, hδc, hδp, hδl, -⟩ := rl_local hU hf hp (v := v) (s₀ := 0) (by simp [hfp])
    have hm : min ε 1 ∈ S := by
      have hsub : Icc (0 : ℝ) (min ε 1) ⊆ Icc (0 - ε) (0 + ε) := fun s hs =>
        ⟨by linarith [hs.1], by linarith [hs.2, min_le_left ε 1]⟩
      exact ⟨⟨(le_min hε.le zero_le_one), min_le_right _ _⟩, δ, hδc.mono hsub, hδp, fun s hs => hδl s (hsub hs)⟩
    have := le_csSup hbdd hm
    have h2 : 0 < min ε 1 := lt_min hε one_pos
    linarith
  · set b := sSup S with hbdef
    have hlt : ∀ s ∈ Ico (0 : ℝ) b, ∃ b' ∈ S, s < b' := fun s hs => exists_lt_of_lt_csSup hSne hs.2
    have hlift : ∀ s ∈ Ico (0 : ℝ) b, γ s ∈ U ∧ f (γ s) = s • v := by
      intro s hs
      obtain ⟨b', hb', hsb'⟩ := hlt s hs
      rw [hγ b' hb' s ⟨hs.1, hsb'.le⟩]
      exact hgl b' hb' s ⟨hs.1, hsb'.le⟩
    have hcont : ContinuousOn γ (Ico 0 b) := by
      intro t ht
      obtain ⟨b', hb', htb'⟩ := hlt t ht
      have h1 : ContinuousWithinAt γ (Icc 0 b') t := by
        refine ((hgc b' hb') t ⟨ht.1, htb'.le⟩).congr (fun y hy => hγ b' hb' y hy) (hγ b' hb' t ⟨ht.1, htb'.le⟩)
      refine h1.mono_of_mem_nhdsWithin ?_
      refine mem_nhdsWithin.2 ⟨Iio b', isOpen_Iio, htb', fun x hx => ⟨hx.2.1, hx.1.le⟩⟩
    have hKb := hapriori b ⟨hpos, hb1⟩ γ hcont hγ0 hlift
    obtain ⟨ε, hε, γ', hγ'c, hγ'eq, hγ'l⟩ := rl_extend hU hf hK hKU hpos hcont hlift hKb
    have hγ'0 : γ' 0 = p := by rw [hγ'eq 0 ⟨le_rfl, hpos⟩, hγ0]
    have hm : min (b + ε) 1 ∈ S := by
      have hsub : Icc (0 : ℝ) (min (b + ε) 1) ⊆ Icc 0 (b + ε) := Icc_subset_Icc_right (min_le_left _ _)
      exact ⟨⟨le_min (by linarith) zero_le_one, min_le_right _ _⟩, γ', hγ'c.mono hsub, hγ'0,
        fun s hs => hγ'l s (hsub hs)⟩
    have hle := le_csSup hbdd hm
    have hb1' : b = 1 := by
      by_contra hne
      have hb_lt : b < 1 := lt_of_le_of_ne hb1 hne
      have : b < min (b + ε) 1 := lt_min (by linarith) hb_lt
      linarith
    have hmin : min (b + ε) 1 = 1 := min_eq_right (by linarith)
    rw [hmin] at hm
    exact hm.2

end CWFlow.Interp

/-- info: 'CWFlow.Interp.exists_ray_lift' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms CWFlow.Interp.exists_ray_lift
