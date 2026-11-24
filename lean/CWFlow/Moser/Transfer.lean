import CWFlow.Moser.Natural

/-!
# Transfer of Moser curves between charts

Let `T` relate two sets of Moser data as in `Moser/Natural.lean`. Let `γ₁` be a Moser curve (data 1) from `q ∈ N`
and `γ₂` a Moser curve (data 2) from `T q`. Neither curve is assumed to stay in `N`. The only assumption is a
closure condition: if the curves correspond on `[0, t)`, then `γ₁ t ∈ N`. Then they correspond on all of `[0, 1]`:
`γ₁ t ∈ N` and `γ₂ t = T (γ₁ t)` (`moserCurve_transfer`). On a manifold the closure condition comes from the
Hausdorff property, and the argument is real induction on `[0, 1]`.
-/

open Filter Set Function

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] [FiniteDimensional ℝ W]

variable {lam₁ β₁ lam₂ β₂ : W → W →L[ℝ] ℝ} {O₁ O₂ N : Set W} {T : W → W}
variable (hO₂ : IsOpen O₂) (hlam₂ : ContDiffOn ℝ ∞ lam₂ O₂) (hβ₂ : ContDiffOn ℝ ∞ β₂ O₂)
  (hN : IsOpen N) (hTN : MapsTo T N O₂) (hT : ContDiffOn ℝ ∞ T N)
  (hTinv : ∀ q ∈ N, ∃ e : W ≃L[ℝ] W, (e : W →L[ℝ] W) = fderiv ℝ T q)
  (hlamT : ∀ q ∈ N, lam₁ q = pb1form T lam₂ q) (hβT : ∀ q ∈ N, β₁ q = pb1form T β₂ q)
include hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT

/-- **Transfer of Moser curves.** -/
theorem moserCurve_transfer {q : W} (hqN : q ∈ N) {γ₁ γ₂ : ℝ → W}
    (h₁ : IsMoserCurve lam₁ β₁ O₁ q γ₁) (h₂ : IsMoserCurve lam₂ β₂ O₂ (T q) γ₂)
    (hcl : ∀ t ∈ Ioc (0 : ℝ) 1, (∀ r ∈ Ico (0 : ℝ) t, γ₁ r ∈ N ∧ γ₂ r = T (γ₁ r)) → γ₁ t ∈ N) :
    ∀ t ∈ Icc (0 : ℝ) 1, γ₁ t ∈ N ∧ γ₂ t = T (γ₁ t) := by
  have hTc : ∀ p ∈ N, ContinuousAt T p := fun p hp =>
    (hT.contDiffAt (hN.mem_nhds hp)).continuousAt
  set s : Set ℝ := {t | ∀ r ∈ Icc (0 : ℝ) t, γ₁ r ∈ N ∧ γ₂ r = T (γ₁ r)} with hs
  have h0 : (0 : ℝ) ∈ s := by
    intro r hr
    have hr0 : r = 0 := le_antisymm hr.2 hr.1
    subst hr0
    exact ⟨h₁.start ▸ hqN, by rw [h₂.start, h₁.start]⟩
  -- `s ∩ [0, 1]` is closed
  have hclosed : IsClosed (s ∩ Icc 0 1) := by
    refine isClosed_of_closure_subset fun t ht => ?_
    have htI : t ∈ Icc (0 : ℝ) 1 := (isClosed_Icc.closure_subset_iff.2 inter_subset_right) ht
    have hbefore : ∀ r ∈ Ico (0 : ℝ) t, γ₁ r ∈ N ∧ γ₂ r = T (γ₁ r) := by
      intro r hr
      obtain ⟨u, hu, hus⟩ := mem_closure_iff.1 ht (Ioi r) isOpen_Ioi hr.2
      exact hus.1 r ⟨hr.1, le_of_lt hu⟩
    refine ⟨?_, htI⟩
    rcases htI.1.eq_or_lt with ht0 | ht0
    · rw [← ht0]
      exact h0
    have hN1 : γ₁ t ∈ N := hcl t ⟨ht0, htI.2⟩ hbefore
    have heq : γ₂ t = T (γ₁ t) := by
      have : (𝓝[Ico 0 t] t).NeBot := mem_closure_iff_nhdsWithin_neBot.1
        (by rw [closure_Ico ht0.ne]; exact right_mem_Icc.2 ht0.le)
      have hsub : Ico (0 : ℝ) t ⊆ Icc 0 1 := fun r hr => ⟨hr.1, hr.2.le.trans htI.2⟩
      have hA : Tendsto γ₂ (𝓝[Ico 0 t] t) (𝓝 (γ₂ t)) := ((h₂.continuousOn t htI).mono hsub).tendsto
      have hB : Tendsto (T ∘ γ₁) (𝓝[Ico 0 t] t) (𝓝 (T (γ₁ t))) :=
        (((hTc _ hN1).comp_continuousWithinAt (h₁.continuousOn t htI)).mono hsub).tendsto
      exact tendsto_nhds_unique_of_eventuallyEq hA hB
        (eventually_nhdsWithin_of_forall fun r hr => (hbefore r hr).2)
    intro r hr
    rcases hr.2.lt_or_eq with hrt | hrt
    · exact hbefore r ⟨hr.1, hrt⟩
    · rw [hrt]
      exact ⟨hN1, heq⟩
  -- from a good time, one can go a bit further
  have hgt : ∀ x ∈ s ∩ Ico (0 : ℝ) 1, (s ∩ Ioc x 1).Nonempty := by
    rintro x ⟨hxs, hx0, hx1⟩
    have hxI : x ∈ Icc (0 : ℝ) 1 := ⟨hx0, hx1.le⟩
    obtain ⟨hxN, hxeq⟩ := hxs x ⟨hx0, le_rfl⟩
    have hmem : γ₁ ⁻¹' N ∈ 𝓝[Icc 0 1] x :=
      (h₁.continuousOn x hxI).preimage_mem_nhdsWithin (hN.mem_nhds hxN)
    obtain ⟨ε, hε, hball⟩ := Metric.mem_nhdsWithin_iff.1 hmem
    set y := min (x + ε / 2) 1 with hy
    have hyle : y ≤ x + ε / 2 := min_le_left _ _
    have hy1 : y ≤ 1 := min_le_right _ _
    have hxy : x < y := lt_min (by linarith) hx1
    have hrI : ∀ r ∈ Icc x y, r ∈ Icc (0 : ℝ) 1 := fun r hr => ⟨hx0.trans hr.1, hr.2.trans hy1⟩
    have hseg : ∀ r ∈ Icc x y, γ₁ r ∈ N := fun r hr => hball ⟨by
      rw [Metric.mem_ball, Real.dist_eq, abs_lt]
      constructor <;> linarith [hr.1, hr.2], hrI r hr⟩
    have hcurve₁ : ∀ r ∈ Icc x y,
        HasDerivWithinAt (T ∘ γ₁) (moserField lam₂ β₂ r (T (γ₁ r))) (Icc 0 1) r := fun r hr => by
      have hTd := (hT.contDiffAt (hN.mem_nhds (hseg r hr))).differentiableAt (by simp)
      rw [moserField_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT (hseg r hr)
        (h₁.good r (hrI r hr))]
      exact hTd.hasFDerivAt.comp_hasDerivWithinAt r (h₁.deriv r (hrI r hr))
    have hcont₁ : ContinuousOn (T ∘ γ₁) (Icc x y) := fun r hr =>
      (hcurve₁ r hr).continuousWithinAt.mono fun r' hr' => hrI r' hr'
    have hC₁ : IsCompact ((fun r => (r, (T ∘ γ₁) r)) '' Icc x y) :=
      isCompact_Icc.image_of_continuousOn (continuousOn_id.prodMk hcont₁)
    have hC₁G : (fun r => (r, (T ∘ γ₁) r)) '' Icc x y ⊆ moserGood lam₂ β₂ O₂ :=
      image_subset_iff.2 fun r hr =>
        good_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT (hseg r hr) (h₁.good r (hrI r hr))
    obtain ⟨f, hf, hc, L, hCL, -, hfX, -⟩ := exists_cutoff hO₂ hlam₂ hβ₂
      (h₂.isCompact_graph.union hC₁) (union_subset h₂.graph_subset hC₁G)
    obtain ⟨K, hK⟩ := lipschitzWith_slice hf hc
    have hIco : ∀ r ∈ Ico x y, r ∈ Ico (0 : ℝ) 1 := fun r hr =>
      ⟨hx0.trans hr.1, lt_of_lt_of_le hr.2 hy1⟩
    have huniq := ODE_solution_unique_of_mem_Icc_right (v := f) (s := fun _ => univ) (K := K)
      (a := x) (b := y) (f := γ₂) (g := T ∘ γ₁) (fun t _ => (hK t).lipschitzOnWith)
      (h₂.continuousOn.mono fun r hr => hrI r hr)
      (fun t ht => by
        rw [hfX t _ (interior_subset (hCL (Or.inl ⟨t, Ico_subset_Icc_self (hIco t ht), rfl⟩)))]
        exact (h₂.deriv t (Ico_subset_Icc_self (hIco t ht))).mono_of_mem_nhdsWithin
          (Icc_mem_nhdsWithin_Ici_of_mem (hIco t ht)))
      (fun _ _ => mem_univ _) hcont₁
      (fun t ht => by
        rw [hfX t _ (interior_subset (hCL (Or.inr ⟨t, Ico_subset_Icc_self ht, rfl⟩)))]
        exact (hcurve₁ t (Ico_subset_Icc_self ht)).mono_of_mem_nhdsWithin
          (Icc_mem_nhdsWithin_Ici_of_mem (hIco t ht)))
      (fun _ _ => mem_univ _) hxeq
    refine ⟨y, fun r hr => ?_, hxy, hy1⟩
    rcases le_or_gt r x with hrx | hrx
    · exact hxs r ⟨hr.1, hrx⟩
    · exact ⟨hseg r ⟨hrx.le, hr.2⟩, huniq ⟨hrx.le, hr.2⟩⟩
  have h1 := hclosed.mem_of_ge_of_forall_exists_gt h0 zero_le_one hgt
  exact fun t ht => h1 t ht

end CWFlow
