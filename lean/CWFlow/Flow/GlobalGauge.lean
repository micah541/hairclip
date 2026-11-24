import CWFlow.Flow.GraphUnique
import CWFlow.Flow.FlowS
import CWFlow.Flow.Gauge
import CWFlow.Flow.LiftTopology

/-!
# The global gauge of a flow

`exists_global_gauge`: let `Φ` be a flow (`IsFlow`) with `Φ 0 = F`. For a short time `δ` there is a jointly
smooth family of maps `ψ_t : M → M`, `ψ_0 = id`, such that `Φ_t ∘ ψ_t` is a section in the Weinstein coordinates
of `F`: every point has a chart `α` and a smooth `ξ(t, y)` on a neighbourhood `B` of the point with
`Φ_t(ψ_t(φ_α⁻¹ y)) = Ψ_α(y, ξ(t, y))`, `ξ(0, ·) = 0`.

Proof: the local gauges of `exists_local_gauge` (one for the chart at each point) agree on overlaps, because
both solve the same equation `Φ_s(p) = Ψ(x, η)` for `(p, η)` (`chartPsi_compat`), start at `(x, 0)`, and the
solution is locally unique (`eqOn_of_locally_unique`). A finite subcover gives the uniform time.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
variable {J : V →L[ℝ] V} {F : M → V}

/-- The local reparametrization read off a chart gauge `ζ`. -/
private abbrev gl (ζ : M → ℝ × E → E) (α : M) (t : ℝ) (x : M) : M :=
  (chart E α).symm (ζ α (t, chart E α x))

/-- The domain in `M` of the gauge of the chart at `α`. -/
private abbrev gD (B : M → Set E) (α : M) : Set M := (chart E α).source ∩ chart E α ⁻¹' B α

/-- Covectors of the chart at `α'` at `x`, read in the chart at `α`. -/
private abbrev gL (α' α x : M) (η : E →L[ℝ] ℝ) : E →L[ℝ] ℝ :=
  η.comp (ContinuousLinearMap.inverse (fderiv ℝ (transition E α' α) (chart E α' x)))

/-- **The global gauge.** A flow from `F` is, for a short time and after a reparametrization `ψ_t` of `M`, a
section in the Weinstein coordinates of `F`. Here `Φ` has velocity `κ • flowVel`; only its
smoothness is used. -/
theorem exists_global_gaugeS {κ : ℝ} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hF : IsLagImmersion E J F) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlowS E J κ T Φ) (hΦ0 : Φ 0 = F) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∃ ψ : ℝ → M → M,
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ x, ψ 0 x = x) ∧
      ∀ x : M, ∃ (α : M) (B : Set E) (ξ : ℝ × E → (E →L[ℝ] ℝ)), IsOpen B ∧ B ⊆ (chart E α).target ∧
        x ∈ (chart E α).source ∧ chart E α x ∈ B ∧ ContDiffOn ℝ ∞ ξ (Icc 0 δ ×ˢ B) ∧
        (∀ y ∈ B, ξ (0, y) = 0) ∧
        ∀ t ∈ Icc (0 : ℝ) δ, ∀ y ∈ B, ψ t ((chart E α).symm y) ∈ (chart E α).source ∧
          (y, ξ (t, y)) ∈ chartDom E J F α ∧
          Φ t (ψ t ((chart E α).symm y)) = chartPsi E J F α (y, ξ (t, y)) := by
  classical
  have hW := weinstein_global (E := E) (M := M) hJi hJJ hF
  have hX0 : ∀ α : M, ∀ z ∈ (chart E α).target,
      flowRep E Φ α (0, z) = chartPsi E J F α (z, 0) := by
    intro α z hz
    rw [(hW.1 α).1.zero_section z ((hW.1 α).2 z hz)]
    show Φ 0 ((chart E α).symm z) = F ((chart E α).symm z)
    rw [hΦ0]
  have hloc := fun α : M => exists_local_gauge (E := E) hdim (isOpen_chartDom hJi hJJ hF α)
    (contDiffOn_chartPsi hJi hJJ hF α)
    (fun q hq => injective_fderiv_of_symplectic (chartPsi_symplectic hJi hJJ hF α q hq)) hT
    (isOpen_chart_target (E := E) α) (X := flowRep E Φ α) (contDiffOn_flowRep hΦ.smooth α) (hX0 α)
    (fun z hz => (hW.1 α).2 z hz) (y₀ := chart E α α) (mem_extChartAt_target (I := 𝓘(ℝ, E)) α)
  choose ε hε hεT B hBo hB0 hBU ζ ξ O A hO hA hAQ hζ hξ h0 hsol huniq using hloc
  have hDo : ∀ α, IsOpen (gD B α) := fun α =>
    (continuousOn_extChartAt (I := 𝓘(ℝ, E)) α).isOpen_inter_preimage
      (isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α) (hBo α)
  have hDα : ∀ α, α ∈ gD B α := fun α => ⟨mem_extChartAt_source (I := 𝓘(ℝ, E)) α, hB0 α⟩
  have hloc0 : ∀ β x, x ∈ gD B β → gl ζ β 0 x = x := by
    intro β x hx
    show (chart E β).symm (ζ β (0, chart E β x)) = x
    rw [(h0 β _ hx.2).1, (chart E β).left_inv hx.1]
  -- two local gauges agree on the overlap of their domains
  have hagree : ∀ α α' x, x ∈ gD B α → x ∈ gD B α' → ∀ m : ℝ, m ≤ ε α → m ≤ ε α' →
      ∀ t ∈ Icc (0 : ℝ) m, gl ζ α t x = gl ζ α' t x := by
    intro α α' x hx hx' m hmα hmα' t ht
    have hcl : ∀ β, x ∈ gD B β → m ≤ ε β →
        ContinuousOn (fun s : ℝ => gl ζ β s x) (Icc 0 m) ∧
        ContinuousOn (fun s : ℝ => ξ β (s, chart E β x)) (Icc 0 m) := by
      intro β hβ hmβ
      have hin : MapsTo (fun s : ℝ => (s, chart E β x)) (Icc 0 m) (Icc 0 (ε β) ×ˢ B β) :=
        fun s hs => ⟨⟨hs.1, hs.2.trans hmβ⟩, hβ.2⟩
      have hpair : ContinuousOn (fun s : ℝ => (s, chart E β x)) (Icc 0 m) :=
        (continuous_id.prodMk continuous_const).continuousOn
      have hz := (hζ β).continuousOn.comp hpair hin
      refine ⟨?_, (hξ β).continuousOn.comp hpair hin⟩
      exact (continuousOn_extChartAt_symm (I := 𝓘(ℝ, E)) β).comp hz
        (fun s hs => (hsol β s ⟨hs.1, hs.2.trans hmβ⟩ _ hβ.2).1)
    have hL : Continuous (gL (E := E) α' α x) := continuous_id.clm_comp continuous_const
    have hc₁ : ContinuousOn (fun s : ℝ => ((s, gl ζ α s x, ξ α (s, chart E α x)) :
        ℝ × M × (E →L[ℝ] ℝ))) (Icc 0 m) :=
      continuousOn_id.prodMk ((hcl α hx hmα).1.prodMk (hcl α hx hmα).2)
    have hc₂ : ContinuousOn (fun s : ℝ => ((s, gl ζ α' s x, gL α' α x (ξ α' (s, chart E α' x))) :
        ℝ × M × (E →L[ℝ] ℝ))) (Icc 0 m) :=
      continuousOn_id.prodMk ((hcl α' hx' hmα').1.prodMk (hL.comp_continuousOn (hcl α' hx' hmα').2))
    have key := eqOn_of_locally_unique hc₁ hc₂ ?_ ?_ ht
    · have h := congrArg (fun p : ℝ × M × (E →L[ℝ] ℝ) => p.2.1) key
      exact h
    · intro t₀ ht₀
      have h1 := hsol α t₀ ⟨ht₀.1, ht₀.2.trans hmα⟩ _ hx.2
      have hp₀ : gl ζ α t₀ x ∈ (chart E α).source := (chart E α).map_target h1.1
      have hcp₀ : chart E α (gl ζ α t₀ x) = ζ α (t₀, chart E α x) := (chart E α).right_inv h1.1
      refine ⟨{p : ℝ × M × (E →L[ℝ] ℝ) | p.2.1 ∈ (chart E α).source ∧
        (p.1, chart E α p.2.1) ∈ O α ∧ (chart E α x, p.2.2) ∈ A α}, ?_, ?_⟩
      · have e1 : ∀ᶠ p : ℝ × M × (E →L[ℝ] ℝ) in 𝓝 (t₀, gl ζ α t₀ x, ξ α (t₀, chart E α x)),
            p.2.1 ∈ (chart E α).source :=
          (continuous_fst.comp continuous_snd).continuousAt.eventually_mem
            ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α).mem_nhds hp₀)
        have hch : ContinuousAt (chart E α) (gl ζ α t₀ x) :=
          (continuousOn_extChartAt (I := 𝓘(ℝ, E)) α).continuousAt
            ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α).mem_nhds hp₀)
        have hc : ContinuousAt (fun p : ℝ × M × (E →L[ℝ] ℝ) => (p.1, chart E α p.2.1))
            (t₀, gl ζ α t₀ x, ξ α (t₀, chart E α x)) :=
          continuousAt_fst.prodMk (ContinuousAt.comp (g := chart E α)
            (f := fun p : ℝ × M × (E →L[ℝ] ℝ) => p.2.1) hch
            (continuous_fst.comp continuous_snd).continuousAt)
        have e2 : ∀ᶠ p : ℝ × M × (E →L[ℝ] ℝ) in 𝓝 (t₀, gl ζ α t₀ x, ξ α (t₀, chart E α x)),
            (p.1, chart E α p.2.1) ∈ O α := by
          refine hc.eventually_mem ((hO α).mem_nhds ?_)
          show (t₀, chart E α (gl ζ α t₀ x)) ∈ O α
          rw [hcp₀]
          exact h1.2.1
        have e3 : ∀ᶠ p : ℝ × M × (E →L[ℝ] ℝ) in 𝓝 (t₀, gl ζ α t₀ x, ξ α (t₀, chart E α x)),
            (chart E α x, p.2.2) ∈ A α :=
          (continuous_const.prodMk (continuous_snd.comp continuous_snd)).continuousAt.eventually_mem
            ((hA α).mem_nhds h1.2.2.1)
        filter_upwards [e1, e2, e3] with p a b c using ⟨a, b, c⟩
      · rintro s hs - ⟨hs2src, hs2O, hs2A⟩
        have hs2src' : gl ζ α' s x ∈ (chart E α).source := hs2src
        have hs2O' : (s, chart E α (gl ζ α' s x)) ∈ O α := hs2O
        have hs2A' : (chart E α x, gL α' α x (ξ α' (s, chart E α' x))) ∈ A α := hs2A
        have hsα : s ∈ Icc (0 : ℝ) (ε α) := ⟨hs.1, hs.2.trans hmα⟩
        have g2 := hsol α' s ⟨hs.1, hs.2.trans hmα'⟩ _ hx'.2
        have hcot : cotLift (transition E α' α) (chart E α' x, ξ α' (s, chart E α' x)) =
            (chart E α x, gL α' α x (ξ α' (s, chart E α' x))) :=
          Prod.ext (transition_chart hx'.1) rfl
        have hchain : flowRep E Φ α (s, chart E α (gl ζ α' s x)) =
            chartPsi E J F α (chart E α x, gL α' α x (ξ α' (s, chart E α' x))) := by
          have a1 : flowRep E Φ α (s, chart E α (gl ζ α' s x)) = Φ s (gl ζ α' s x) := by
            show Φ s ((chart E α).symm (chart E α (gl ζ α' s x))) = _
            rw [(chart E α).left_inv hs2src']
          have a2 : Φ s (gl ζ α' s x) =
              chartPsi E J F α' (chart E α' x, ξ α' (s, chart E α' x)) := g2.2.2.2
          have a3 := chartPsi_compat hF (z := α') (w := α)
            (q := (chart E α' x, ξ α' (s, chart E α' x))) (hAQ α' g2.2.2.1)
            (chart_mem_overlap hx'.1 hx.1) (by rw [hcot]; exact hAQ α hs2A')
          rw [a1, a2, ← a3, hcot]
        obtain ⟨u1, u2⟩ := huniq α s hsα _ hx.2 _ _ hs2O' hs2A' hchain
        have hpp : gl ζ α s x = gl ζ α' s x := by
          show (chart E α).symm (ζ α (s, chart E α x)) = _
          rw [← u1, (chart E α).left_inv hs2src']
        show ((s, gl ζ α s x, ξ α (s, chart E α x)) : ℝ × M × (E →L[ℝ] ℝ)) =
          (s, gl ζ α' s x, gL α' α x (ξ α' (s, chart E α' x)))
        rw [hpp, u2]
    · show ((0, gl ζ α 0 x, ξ α (0, chart E α x)) : ℝ × M × (E →L[ℝ] ℝ)) =
        (0, gl ζ α' 0 x, gL α' α x (ξ α' (0, chart E α' x)))
      rw [hloc0 α x hx, hloc0 α' x hx', (h0 α _ hx.2).2, (h0 α' _ hx'.2).2]
      show ((0, x, 0) : ℝ × M × (E →L[ℝ] ℝ)) = (0, x, (0 : E →L[ℝ] ℝ).comp _)
      rw [ContinuousLinearMap.zero_comp]
  -- a finite subcover and a uniform time
  have hmin : ∀ S : Finset M, ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∀ k ∈ S, δ ≤ ε k := by
    intro S
    induction S using Finset.induction_on with
    | empty => exact ⟨T, hT, le_rfl, by simp⟩
    | insert a s ha ih =>
      obtain ⟨δ, h1, h2, h3⟩ := ih
      refine ⟨min δ (ε a), lt_min h1 (hε a), (min_le_left _ _).trans h2, ?_⟩
      intro k hk
      rcases Finset.mem_insert.1 hk with hka | hk
      · rw [hka]
        exact min_le_right _ _
      · exact (min_le_left _ _).trans (h3 k hk)
  obtain ⟨S, hS⟩ := isCompact_univ.elim_finite_subcover (gD B) hDo
    (fun x _ => mem_iUnion.2 ⟨x, hDα x⟩)
  obtain ⟨δ, hδ, hδT, hδS⟩ := hmin S
  have hcov : ∀ x : M, ∃ k, k ∈ S ∧ x ∈ gD B k := fun x => by
    have h := hS (mem_univ x)
    simp only [mem_iUnion] at h
    obtain ⟨k, hk, hxk⟩ := h
    exact ⟨k, hk, hxk⟩
  choose k hkS hkD using hcov
  have hψeq : ∀ α ∈ S, ∀ x ∈ gD B α, ∀ t ∈ Icc (0 : ℝ) δ, gl ζ (k x) t x = gl ζ α t x :=
    fun α hα x hx t ht => hagree (k x) α x (hkD x) hx δ (hδS _ (hkS x)) (hδS α hα) t ht
  have hsm : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞
      (fun q : ℝ × M => gl ζ (k q.2) q.1 q.2) (Icc 0 δ ×ˢ univ) := by
    rintro ⟨t, x⟩ hq
    have hδα : δ ≤ ε (k x) := hδS _ (hkS x)
    have c1 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
        (fun q : ℝ × M => ((q.1, chart E (k x) q.2) : ℝ × E)) (Icc 0 δ ×ˢ gD B (k x)) := by
      have ca : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => q.1)
          (Icc 0 δ ×ˢ gD B (k x)) := contMDiffOn_fst
      have cb : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => chart E (k x) q.2)
          (Icc 0 δ ×ˢ gD B (k x)) :=
        ContMDiffOn.comp (g := chart E (k x)) (f := fun q : ℝ × M => q.2)
          (t := (chart E (k x)).source)
          (by rw [extChartAt_source]; exact contMDiffOn_extChartAt) contMDiffOn_snd
          (fun q hq => hq.2.1)
      exact ca.prodMk_space cb
    have c2 : ContMDiffOn 𝓘(ℝ, ℝ × E) 𝓘(ℝ, E) ∞ (ζ (k x)) (Icc 0 δ ×ˢ B (k x)) :=
      contMDiffOn_iff_contDiffOn.mpr
        ((hζ (k x)).mono (prod_mono (Icc_subset_Icc_right hδα) subset_rfl))
    have c3 := ContMDiffOn.comp (g := ζ (k x))
      (f := fun q : ℝ × M => ((q.1, chart E (k x) q.2) : ℝ × E)) c2 c1 (fun q hq => ⟨hq.1, hq.2.2⟩)
    have c4 := ContMDiffOn.comp (g := (chart E (k x)).symm)
      (f := fun q : ℝ × M => ζ (k x) (q.1, chart E (k x) q.2))
      (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E)) (k x)) c3
      (fun q hq => (hsol (k x) q.1 ⟨hq.1.1, hq.1.2.trans hδα⟩ _ hq.2.2).1)
    have c5 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞
        (fun q : ℝ × M => gl ζ (k q.2) q.1 q.2) (Icc 0 δ ×ˢ gD B (k x)) (t, x) :=
      (c4 (t, x) ⟨hq.1, hkD x⟩).congr
        (fun q hq' => hψeq (k x) (hkS x) q.2 hq'.2 q.1 hq'.1)
        (hψeq (k x) (hkS x) x (hkD x) t hq.1)
    refine c5.mono_of_mem_nhdsWithin (mem_nhdsWithin.2 ⟨univ ×ˢ gD B (k x),
      isOpen_univ.prod (hDo (k x)), ⟨mem_univ _, hkD x⟩, fun p hp => ⟨hp.2.1, hp.1.2⟩⟩)
  refine ⟨δ, hδ, hδT, fun t x => gl ζ (k x) t x, hsm, fun x => hloc0 (k x) x (hkD x), fun x => ?_⟩
  refine ⟨k x, B (k x), ξ (k x), hBo _, hBU _, (hkD x).1, (hkD x).2,
    (hξ _).mono (prod_mono (Icc_subset_Icc_right (hδS _ (hkS x))) subset_rfl),
    fun y hy => (h0 _ y hy).2, fun t ht y hy => ?_⟩
  have g := hsol (k x) t ⟨ht.1, ht.2.trans (hδS _ (hkS x))⟩ y hy
  have hyt : y ∈ (chart E (k x)).target := hBU _ hy
  have hyD : (chart E (k x)).symm y ∈ gD B (k x) := by
    refine ⟨(chart E (k x)).map_target hyt, ?_⟩
    show chart E (k x) ((chart E (k x)).symm y) ∈ B (k x)
    rw [(chart E (k x)).right_inv hyt]
    exact hy
  have e : gl ζ (k ((chart E (k x)).symm y)) t ((chart E (k x)).symm y) =
      (chart E (k x)).symm (ζ (k x) (t, y)) := by
    rw [hψeq (k x) (hkS x) _ hyD t ht]
    show (chart E (k x)).symm (ζ (k x) (t, chart E (k x) ((chart E (k x)).symm y))) = _
    rw [(chart E (k x)).right_inv hyt]
  refine ⟨?_, hAQ _ g.2.2.1, ?_⟩
  · show gl ζ (k ((chart E (k x)).symm y)) t ((chart E (k x)).symm y) ∈ (chart E (k x)).source
    rw [e]
    exact (chart E (k x)).map_target g.1
  · show Φ t (gl ζ (k ((chart E (k x)).symm y)) t ((chart E (k x)).symm y)) = _
    rw [e]
    exact g.2.2.2

/-- **The global gauge.** A flow from `F` is, for a short time and after a reparametrization `ψ_t` of `M`, a
section in the Weinstein coordinates of `F`. -/
theorem exists_global_gauge (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hF : IsLagImmersion E J F) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) (hΦ0 : Φ 0 = F) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∃ ψ : ℝ → M → M,
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ x, ψ 0 x = x) ∧
      ∀ x : M, ∃ (α : M) (B : Set E) (ξ : ℝ × E → (E →L[ℝ] ℝ)), IsOpen B ∧ B ⊆ (chart E α).target ∧
        x ∈ (chart E α).source ∧ chart E α x ∈ B ∧ ContDiffOn ℝ ∞ ξ (Icc 0 δ ×ˢ B) ∧
        (∀ y ∈ B, ξ (0, y) = 0) ∧
        ∀ t ∈ Icc (0 : ℝ) δ, ∀ y ∈ B, ψ t ((chart E α).symm y) ∈ (chart E α).source ∧
          (y, ξ (t, y)) ∈ chartDom E J F α ∧
          Φ t (ψ t ((chart E α).symm y)) = chartPsi E J F α (y, ξ (t, y)) := by
  exact exists_global_gaugeS (κ := 1) hJi hJJ hF hdim hT hΦ.toS hΦ0

end Bridge

end CWFlow
