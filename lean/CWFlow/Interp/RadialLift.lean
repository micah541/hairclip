import Mathlib.Topology.IsLocalHomeomorph
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Analysis.Convex.PathConnected
import Mathlib.Analysis.Normed.Module.Convex

/-!
# The radial lifting lemma

Let `f : X → E` be a local homeomorphism on an open set `U` of a Hausdorff space, mapping `U` into the ball
`B = B(0, r)`, and let `p ∈ U` with `f p = 0`.  Suppose every ray `s ↦ s • v` (`v ∈ B`, `s ∈ [0,1]`) lifts to a
path in `U` starting at `p`.  Then the endpoint map `σ v = (lift of the ray of v)(1)` is a continuous section of `f`
over `B`, and its image is exactly the connected component of `p` in `U` (`radial_lift`).

The proof is: uniqueness of lifts (`eqOn_of_lift`), a real induction on the time showing that the lifts depend
continuously on `v` (`tendsto_lift`), and an open/closed argument for the image.
-/

open Set Filter Metric

open scoped Topology

namespace CWFlow.Interp

/-- **Uniqueness of lifts.** Let `γ₁, γ₂ : [a, b] → Q` be continuous with `p ∘ γ₁ = p ∘ γ₂` and `γ₁ a = γ₂ a`.
If `p` is injective near every point of `γ₁`, then `γ₁ = γ₂` on `[a, b]`. -/
theorem eqOn_of_lift {Q B : Type*} [TopologicalSpace Q] [T2Space Q] {p : Q → B}
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

variable {X E : Type*} [TopologicalSpace X] [T2Space X] [NormedAddCommGroup E] [NormedSpace ℝ E]

omit [T2Space X] [NormedSpace ℝ E] in
/-- A local homeomorphism on `U` is injective near each point of `U`. -/
theorem injOn_nhds_of_isLocalHomeomorphOn {f : X → E} {U : Set X} (hf : IsLocalHomeomorphOn f U) {x : X}
    (hx : x ∈ U) : ∃ N ∈ 𝓝 x, InjOn f N := by
  obtain ⟨e, hxe, rfl⟩ := hf x hx
  exact ⟨e.source, e.open_source.mem_nhds hxe, e.injOn⟩

/-- **Continuous dependence of the lifts.** Let `γ v` be a lift in `U` of the ray `s ↦ s • v` for each `v` in the
ball, starting at `p`.  Then for every `s ∈ [0,1]`, `v ↦ γ v s` is continuous at each `v₀` of the ball. -/
theorem tendsto_lift {f : X → E} {U : Set X} {p : X} {r : ℝ} (hf : IsLocalHomeomorphOn f U)
    {γ : E → ℝ → X} (hc : ∀ v ∈ ball (0 : E) r, ContinuousOn (γ v) (Icc 0 1))
    (h0 : ∀ v ∈ ball (0 : E) r, γ v 0 = p)
    (hs : ∀ v ∈ ball (0 : E) r, ∀ s ∈ Icc (0 : ℝ) 1, γ v s ∈ U ∧ f (γ v s) = s • v)
    {v₀ : E} (hv₀ : v₀ ∈ ball (0 : E) r) :
    ∀ s ∈ Icc (0 : ℝ) 1, Tendsto (fun v => γ v s) (𝓝 v₀) (𝓝 (γ v₀ s)) := by
  have hball : ball (0 : E) r ∈ 𝓝 v₀ := isOpen_ball.mem_nhds hv₀
  set P : ℝ → Prop := fun s => Tendsto (fun v => γ v s) (𝓝 v₀) (𝓝 (γ v₀ s)) with hP
  -- the local step
  have L : ∀ s₀ ∈ Icc (0 : ℝ) 1, ∃ δ > 0, ∀ s₁ ∈ Icc (0 : ℝ) 1, s₀ - δ < s₁ → P s₁ →
      ∀ s ∈ Icc (0 : ℝ) 1, s₁ ≤ s → s < s₀ + δ → P s := by
    intro s₀ hs₀
    obtain ⟨e, hxe, rfl⟩ := hf (γ v₀ s₀) (hs v₀ hv₀ s₀ hs₀).1
    have hA : ∀ᶠ s in 𝓝[Icc (0 : ℝ) 1] s₀, γ v₀ s ∈ e.source :=
      (hc v₀ hv₀ s₀ hs₀).preimage_mem_nhdsWithin (e.open_source.mem_nhds hxe)
    obtain ⟨ε, hε, hεA⟩ := Metric.mem_nhdsWithin_iff.mp hA
    have hs₀v : s₀ • v₀ ∈ e.target := by
      rw [← (hs v₀ hv₀ s₀ hs₀).2]
      exact e.map_source hxe
    have hB : {q : ℝ × E | q.1 • q.2 ∈ e.target} ∈ 𝓝 (s₀, v₀) :=
      (continuous_fst.smul continuous_snd).continuousAt.preimage_mem_nhds
        (e.open_target.mem_nhds hs₀v)
    obtain ⟨u, hu, w, hw, huw⟩ := mem_nhds_prod_iff.mp hB
    obtain ⟨δ₁, hδ₁, hδu⟩ := Metric.mem_nhds_iff.mp hu
    have hδ : 0 < min δ₁ ε := lt_min hδ₁ hε
    refine ⟨min δ₁ ε, hδ, ?_⟩
    intro s₁ hs₁ hs₁l hP1 s hsI hs₁s hsu
    -- membership in the interval of radius `δ`
    have hin : ∀ t, s₁ ≤ t → t ≤ s → t ∈ ball s₀ δ₁ ∩ ball s₀ ε := by
      intro t h1 h2
      have a1 : s₀ - min δ₁ ε < t := lt_of_lt_of_le hs₁l h1
      have a2 : t < s₀ + min δ₁ ε := lt_of_le_of_lt h2 hsu
      have m1 := min_le_left δ₁ ε
      have m2 := min_le_right δ₁ ε
      constructor
      · rw [Metric.mem_ball, Real.dist_eq, abs_sub_lt_iff]; constructor <;> linarith
      · rw [Metric.mem_ball, Real.dist_eq, abs_sub_lt_iff]; constructor <;> linarith
    have hs₁Icc : s₁ ∈ Icc s₁ s := ⟨le_rfl, hs₁s⟩
    have hsIcc : s ∈ Icc s₁ s := ⟨hs₁s, le_rfl⟩
    have hsrc : ∀ t ∈ Icc s₁ s, γ v₀ t ∈ e.source := by
      intro t ht
      have := hin t ht.1 ht.2
      have htI : t ∈ Icc (0 : ℝ) 1 := ⟨hs₁.1.trans ht.1, ht.2.trans hsI.2⟩
      exact hεA ⟨this.2, htI⟩
    have hG1 : ∀ᶠ v in 𝓝 v₀, γ v s₁ ∈ e.source :=
      hP1.eventually (e.open_source.mem_nhds (hsrc s₁ hs₁Icc))
    have key : ∀ᶠ v in 𝓝 v₀, γ v s = e.symm (s • v) := by
      filter_upwards [hG1, hw, hball] with v hv1 hvw hvb
      have htgt : ∀ t ∈ Icc s₁ s, t • v ∈ e.target := fun t ht =>
        huw (mk_mem_prod (hδu (hin t ht.1 ht.2).1) hvw)
      have hstart : e.symm (s₁ • v) = γ v s₁ := by
        symm
        refine (e.eq_symm_apply hv1 (htgt s₁ hs₁Icc)).2 ?_
        exact (hs v hvb s₁ hs₁).2
      have hcont₁ : ContinuousOn (fun t : ℝ => e.symm (t • v)) (Icc s₁ s) :=
        e.continuousOn_symm.comp (continuous_id.smul continuous_const).continuousOn htgt
      have hcont₂ : ContinuousOn (γ v) (Icc s₁ s) :=
        (hc v hvb).mono fun t ht => ⟨hs₁.1.trans ht.1, ht.2.trans hsI.2⟩
      have := eqOn_of_lift (p := e) hcont₁ hcont₂
        (fun t ht => ⟨e.source, e.open_source.mem_nhds (e.map_target (htgt t ht)), e.injOn⟩)
        (fun t ht => by
          rw [e.right_inv (htgt t ht)]
          exact ((hs v hvb t ⟨hs₁.1.trans ht.1, ht.2.trans hsI.2⟩).2).symm)
        hstart hsIcc
      exact this.symm
    have hcs : ContinuousAt (fun v : E => e.symm (s • v)) v₀ := by
      have hsv₀ : s • v₀ ∈ e.target := huw (mk_mem_prod (hδu (hin s hs₁s le_rfl).1) (mem_of_mem_nhds hw))
      exact (e.continuousAt_symm hsv₀).comp (continuous_const.smul continuous_id).continuousAt
    have hval : e.symm (s • v₀) = γ v₀ s := by
      have hsrc' := hsrc s hsIcc
      have := (hs v₀ hv₀ s hsI).2
      rw [← this]
      exact e.left_inv hsrc'
    refine Tendsto.congr' (key.mono fun v hv => hv.symm) ?_
    rw [← hval]
    exact hcs
  -- the real induction
  have hP0 : P 0 := by
    have : (fun v => γ v 0) =ᶠ[𝓝 v₀] fun _ => γ v₀ 0 := by
      filter_upwards [hball] with v hv
      rw [h0 v hv, h0 v₀ hv₀]
    exact tendsto_const_nhds.congr' this.symm
  let T : Set ℝ := {s | s ∈ Icc (0 : ℝ) 1 ∧ ∀ s' ∈ Icc (0 : ℝ) s, P s'}
  have h0T : (0 : ℝ) ∈ T := ⟨⟨le_rfl, zero_le_one⟩, fun s' hs' => by
    rw [le_antisymm hs'.2 hs'.1]; exact hP0⟩
  have hbdd : BddAbove T := ⟨1, fun s hs => hs.1.2⟩
  set c := sSup T with hcdef
  have hc0 : 0 ≤ c := le_csSup hbdd h0T
  have hc1 : c ≤ 1 := csSup_le ⟨0, h0T⟩ fun s hs => hs.1.2
  obtain ⟨δ, hδ, hLδ⟩ := L c ⟨hc0, hc1⟩
  obtain ⟨s₁, hs₁T, hs₁c⟩ := exists_lt_of_lt_csSup ⟨0, h0T⟩ (show c - δ < c by linarith)
  have hs₁le : s₁ ≤ c := le_csSup hbdd hs₁T
  have hall : ∀ s ∈ Icc (0 : ℝ) 1, s < c + δ → P s := by
    intro s hs hsc
    by_cases hle : s ≤ s₁
    · exact hs₁T.2 s ⟨hs.1, hle⟩
    · exact hLδ s₁ hs₁T.1 hs₁c (hs₁T.2 s₁ ⟨hs₁T.1.1, le_rfl⟩) s hs (not_le.mp hle).le hsc
  have hc_eq : c = 1 := by
    by_contra hne
    have hlt : c < 1 := lt_of_le_of_ne hc1 hne
    have hmT : min (c + δ / 2) 1 ∈ T := by
      refine ⟨⟨le_min (by linarith) zero_le_one, min_le_right _ _⟩, fun s' hs' => ?_⟩
      refine hall s' ⟨hs'.1, hs'.2.trans (min_le_right _ _)⟩ ?_
      have := min_le_left (c + δ / 2) 1
      linarith [hs'.2]
    have := le_csSup hbdd hmT
    have h2 : c < min (c + δ / 2) 1 := lt_min (by linarith) hlt
    linarith
  intro s hs
  exact hall s hs (by linarith [hs.2, hc_eq, hδ])

/-- **The radial lifting lemma.** Let `f` be a local homeomorphism on the open set `U` (in a Hausdorff space) with
`f p = 0` and `f U ⊆ B(0, r)`.  If every ray `s ↦ s • v` (`v ∈ B(0, r)`, `s ∈ [0, 1]`) lifts to a path in `U` starting
at `p`, then there is a continuous section `σ` of `f` over the ball with `σ 0 = p` whose image is the connected
component of `p` in `U`. -/
theorem radial_lift {X E : Type*} [TopologicalSpace X] [T2Space X]
    [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : X → E} {U : Set X} {p : X} {r : ℝ}
    (hU : IsOpen U) (hf : IsLocalHomeomorphOn f U) (hp : p ∈ U) (hfp : f p = 0)
    (hmaps : Set.MapsTo f U (Metric.ball (0 : E) r))
    (hlift : ∀ v ∈ Metric.ball (0 : E) r, ∃ γ : ℝ → X, ContinuousOn γ (Set.Icc 0 1) ∧ γ 0 = p ∧
      ∀ s ∈ Set.Icc (0 : ℝ) 1, γ s ∈ U ∧ f (γ s) = s • v) :
    ∃ σ : E → X, ContinuousOn σ (Metric.ball (0 : E) r) ∧ σ 0 = p ∧
      (∀ v ∈ Metric.ball (0 : E) r, σ v ∈ U ∧ f (σ v) = v) ∧
      σ '' Metric.ball (0 : E) r = connectedComponentIn U p := by
  have h0b : (0 : E) ∈ ball (0 : E) r := by
    have := hmaps hp
    rwa [hfp] at this
  have : Nonempty X := ⟨p⟩
  choose! γ hγc hγ0 hγs using hlift
  have hcont : ∀ v₀ ∈ ball (0 : E) r, ContinuousAt (fun v => γ v 1) v₀ := fun v₀ hv₀ =>
    tendsto_lift hf hγc hγ0 hγs hv₀ 1 ⟨zero_le_one, le_rfl⟩
  have hsec : ∀ v ∈ ball (0 : E) r, γ v 1 ∈ U ∧ f (γ v 1) = v := fun v hv => by
    have := hγs v hv 1 ⟨zero_le_one, le_rfl⟩
    rwa [one_smul] at this
  have hσ0 : γ 0 1 = p := by
    have := eqOn_of_lift (p := f) (γ₁ := fun _ => p) (γ₂ := γ 0) (a := 0) (b := 1)
      continuousOn_const (hγc 0 h0b)
      (fun t _ => injOn_nhds_of_isLocalHomeomorphOn hf hp)
      (fun t ht => by rw [hfp, (hγs 0 h0b t ht).2, smul_zero]) (hγ0 0 h0b).symm
      ⟨zero_le_one, le_rfl⟩
    exact this.symm
  refine ⟨fun v => γ v 1, fun v hv => (hcont v hv).continuousWithinAt, hσ0, hsec, ?_⟩
  set σ : E → X := fun v => γ v 1 with hσ
  set A : Set X := σ '' ball (0 : E) r with hA
  have hAU : A ⊆ U := by
    rintro _ ⟨v, hv, rfl⟩
    exact (hsec v hv).1
  have hp_A : p ∈ A := ⟨0, h0b, hσ0⟩
  apply Subset.antisymm
  · exact ((convex_ball (0 : E) r).isPreconnected.image σ
      (fun v hv => (hcont v hv).continuousWithinAt)).subset_connectedComponentIn hp_A hAU
  · -- A is open
    have hAopen : IsOpen A := by
      rw [isOpen_iff_mem_nhds]
      rintro _ ⟨v₀, hv₀, rfl⟩
      obtain ⟨e, hxe, hfe⟩ := hf (σ v₀) (hsec v₀ hv₀).1
      have hpre : σ ⁻¹' e.source ∩ ball (0 : E) r ∈ 𝓝 v₀ :=
        inter_mem ((hcont v₀ hv₀).preimage_mem_nhds (e.open_source.mem_nhds hxe))
          (isOpen_ball.mem_nhds hv₀)
      obtain ⟨N, hNsub, hNopen, hv₀N⟩ := _root_.mem_nhds_iff.mp hpre
      have hW : IsOpen (e.source ∩ e ⁻¹' N) :=
        e.continuousOn.isOpen_inter_preimage e.open_source hNopen
      refine mem_of_superset (hW.mem_nhds ⟨hxe, ?_⟩) ?_
      · show e (σ v₀) ∈ N
        rw [← hfe, (hsec v₀ hv₀).2]; exact hv₀N
      · rintro y ⟨hy, hyN⟩
        have hyb := hNsub hyN
        refine ⟨e y, hyb.2, ?_⟩
        have h1 := hsec (e y) hyb.2
        have h2 : σ (e y) ∈ e.source := hyb.1
        have : e (σ (e y)) = e y := by
          have h3 : f (σ (e y)) = e y := h1.2
          rw [hfe] at h3; exact h3
        exact e.injOn h2 hy this
    -- A is closed relative to U
    have hAclosed : closure A ∩ connectedComponentIn U p ⊆ A := by
      rintro x ⟨hxcl, hxC⟩
      have hxU : x ∈ U := connectedComponentIn_subset U p hxC
      obtain ⟨e, hxe, hfe⟩ := hf x hxU
      have hfx : f x ∈ ball (0 : E) r := hmaps hxU
      have hfc : ContinuousAt f x := by rw [hfe]; exact e.continuousAt hxe
      have hlim : Tendsto (fun y => σ (f y)) (𝓝[A] x) (𝓝 (σ (f x))) :=
        (((hcont (f x) hfx).comp hfc).mono_left nhdsWithin_le_nhds)
      have hid : Tendsto (fun y : X => y) (𝓝[A] x) (𝓝 x) := tendsto_id.mono_left nhdsWithin_le_nhds
      have : (𝓝[A] x).NeBot := mem_closure_iff_nhdsWithin_neBot.mp hxcl
      have heq : (fun y => σ (f y)) =ᶠ[𝓝[A] x] fun y => y := by
        filter_upwards [self_mem_nhdsWithin] with y hy
        obtain ⟨v, hv, rfl⟩ := hy
        show σ (f (σ v)) = σ v
        rw [(hsec v hv).2]
      have hxeq : σ (f x) = x := tendsto_nhds_unique (hlim.congr' heq) hid
      exact ⟨f x, hfx, hxeq⟩
    exact isPreconnected_connectedComponentIn.subset_of_closure_inter_subset hAopen
      ⟨p, mem_connectedComponentIn hp, hp_A⟩ hAclosed

end CWFlow.Interp

/-- info: 'CWFlow.Interp.radial_lift' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms CWFlow.Interp.radial_lift
