import CWFlow.Extension.TimeGraphFlow

/-!
# A uniform version of the time graph construction

A smooth family `Y_t` of Lagrangian immersions of a ball, whose horizontal maps `χ_t = horC ∘ Y_t` are
`C¹`-close to the identity for all `t` in an open interval `(c, b)`, is on the whole interval the graph of
`∇φ(t, ·)` over the ball of radius `5`, after the reparametrization `ξ_t = χ_t⁻¹`; `φ` is jointly smooth
and `∂ₜφ = potX` (`exists_time_graph_unif`). The qualitative `exists_param_inverse` is replaced by the
given uniform bounds (`tu_inverse`), and the rest follows `exists_time_graph` and
`exists_time_graph_flow`.
-/

open Set Metric
open CWFlow.Interp
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-- Uniform inverse with a time parameter. -/
theorem tu_inverse {c b : ℝ} (χ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (hχ : ContDiffOn ℝ ∞ χ (Ioo c b ×ˢ ball (0 : 𝔼) 6))
    (hbd : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : 𝔼) (11 / 2),
      ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      ‖χ (t, y) - y‖ ≤ 1 / 8) :
    ∃ ξ : ℝ × 𝔼 → 𝔼, ContDiffOn ℝ ∞ ξ (Ioo c b ×ˢ ball (0 : 𝔼) 5) ∧
      ∀ t ∈ Ioo c b,
        InjOn (fun y => χ (t, y)) (ball (0 : 𝔼) (11 / 2)) ∧
        (∀ x ∈ ball (0 : 𝔼) 5, ξ (t, x) ∈ ball (0 : 𝔼) (11 / 2) ∧ χ (t, ξ (t, x)) = x) ∧
        (∀ y ∈ ball (0 : 𝔼) (11 / 2), χ (t, y) ∈ ball (0 : 𝔼) 5 → ξ (t, χ (t, y)) = y) := by
  classical
  have hU : IsOpen (Ioo c b ×ˢ ball (0 : 𝔼) 6) := isOpen_Ioo.prod isOpen_ball
  have hχc : ∀ p ∈ Ioo c b ×ˢ ball (0 : 𝔼) 6, ContDiffAt ℝ ∞ χ p :=
    fun p hp => hχ.contDiffAt (hU.mem_nhds hp)
  have hdiff : ∀ p ∈ Ioo c b ×ˢ ball (0 : 𝔼) 6, DifferentiableAt ℝ χ p :=
    fun p hp => (hχc p hp).differentiableAt (by simp)
  have h56 : ball (0 : 𝔼) (11 / 2) ⊆ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num)
  have per : ∀ t ∈ Ioo c b, InjOn (fun y => χ (t, y)) (ball (0 : 𝔼) (11 / 2)) ∧
      (∀ x ∈ ball (0 : 𝔼) 5, ∃ y ∈ ball (0 : 𝔼) (11 / 2), χ (t, y) = x) := by
    intro t ht
    have hdy : ∀ y ∈ ball (0 : 𝔼) (11 / 2), DifferentiableAt ℝ (fun y' => χ (t, y')) y :=
      fun y hy => (pi_slice_hasFDerivAt (hdiff (t, y) ⟨ht, h56 hy⟩)).differentiableAt
    have hDle : ∀ y ∈ ball (0 : 𝔼) (11 / 2),
        ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ 𝔼‖ ≤ 1 / 4 :=
      fun y hy => ((hbd t ht y hy).1).trans (by norm_num)
    have happ := pi_approx (g := fun y' => χ (t, y')) (convex_ball (0 : 𝔼) (11 / 2)) hdy hDle
    refine ⟨?_, ?_⟩
    · intro x hx y hy hxy
      have h := happ x hx y hy
      have h' : ‖x - y‖ ≤ 1 / 4 * ‖x - y‖ := by
        simp only at hxy
        rw [hxy, sub_self, zero_sub, norm_neg] at h
        exact h
      have : ‖x - y‖ = 0 := by linarith [norm_nonneg (x - y)]
      exact sub_eq_zero.1 (norm_eq_zero.1 this)
    · intro x hx
      have hxR' : x ∈ ball (0 : 𝔼) (11 / 2) := ball_subset_ball (by norm_num) hx
      have hsub : closedBall x (1 / 4) ⊆ ball (0 : 𝔼) (11 / 2) := by
        intro z hz
        rw [mem_closedBall, dist_eq_norm] at hz
        rw [mem_ball_zero_iff]
        rw [mem_ball_zero_iff] at hx
        have : ‖z‖ ≤ ‖x‖ + ‖z - x‖ := by
          calc ‖z‖ = ‖x + (z - x)‖ := by rw [add_sub_cancel]
            _ ≤ ‖x‖ + ‖z - x‖ := norm_add_le _ _
        linarith
      have hs := pi_surj (f := fun y' => χ (t, y')) isOpen_ball happ (b := x)
        (ρ := 1 / 4) (by norm_num) hsub
      have hxm : x ∈ closedBall (χ (t, x)) (3 / 4 * (1 / 4)) := by
        rw [mem_closedBall, dist_comm, dist_eq_norm]
        have := (hbd t ht x hxR').2
        linarith
      exact hs hxm |> fun ⟨y, hy, hyx⟩ => ⟨y, hy, hyx⟩
  set ξ : ℝ × 𝔼 → 𝔼 := fun p => Function.invFunOn (fun y => χ (p.1, y)) (ball (0 : 𝔼) (11 / 2)) p.2
    with hξ
  have hpos : ∀ t ∈ Ioo c b, ∀ x ∈ ball (0 : 𝔼) 5,
      ξ (t, x) ∈ ball (0 : 𝔼) (11 / 2) ∧ χ (t, ξ (t, x)) = x := fun t ht x hx =>
    Function.invFunOn_pos ((per t ht).2 x hx)
  have hleft : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : 𝔼) (11 / 2), ξ (t, χ (t, y)) = y :=
    fun t ht y hy => (per t ht).1.leftInvOn_invFunOn hy
  refine ⟨ξ, ?_, fun t ht => ⟨(per t ht).1, hpos t ht, fun y hy _ => hleft t ht y hy⟩⟩
  intro p hp
  obtain ⟨t, x⟩ := p
  obtain ⟨ht, hx⟩ := hp
  obtain ⟨y, hy, hχy⟩ := (per t ht).2 x hx
  have hdy : ∀ y ∈ ball (0 : 𝔼) (11 / 2), DifferentiableAt ℝ (fun y' => χ (t, y')) y :=
    fun y hy => (pi_slice_hasFDerivAt (hdiff (t, y) ⟨ht, h56 hy⟩)).differentiableAt
  have hD : ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ 𝔼‖ ≤ 1 / 4 :=
    ((hbd t ht y hy).1).trans (by norm_num)
  obtain ⟨g, hg, hgy, hginv⟩ := pi_local_inverse (hχc (t, y) ⟨ht, h56 hy⟩) (pi_inj_of_le hD)
  rw [hχy] at hg hgy hginv
  have hnb : ∀ᶠ z in 𝓝 (t, x), z ∈ Ioo c b ×ˢ ball (0 : 𝔼) 5 :=
    (isOpen_Ioo.prod isOpen_ball).mem_nhds ⟨ht, hx⟩
  have hgb : ∀ᶠ z in 𝓝 (t, x), g z ∈ ball (0 : 𝔼) (11 / 2) :=
    hg.continuousAt.eventually_mem (isOpen_ball.mem_nhds (by rw [hgy]; exact hy))
  have hat : ContDiffAt ℝ ∞ ξ (t, x) := by
    refine hg.congr_of_eventuallyEq ?_
    filter_upwards [hnb, hginv, hgb] with z hz hz2 hz3
    have h1 := hleft z.1 hz.1 (g z) hz3
    rw [hz2] at h1
    exact h1
  exact hat.contDiffWithinAt

theorem exists_time_graph_unif {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) (y₀ : V)
    {c b : ℝ} (hcb : c < b) {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Ioo c b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6))
    (hinj : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    (hχ : ∀ t ∈ Ioo c b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
      ‖fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      ‖horC U y₀ (Y (t, y)) - y‖ ≤ 1 / 8) :
    ∃ (ξ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
      (φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
      ContDiffOn ℝ ∞ ξ (Ioo c b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5) ∧
      ContDiffOn ℝ ∞ φ (Ioo c b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5) ∧
      ∀ t ∈ Ioo c b,
        InjOn (fun y => horC U y₀ (Y (t, y))) (ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2)) ∧
        (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
          ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) ∧ horC U y₀ (Y (t, ξ (t, x))) = x ∧
          Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ (t, x')) x ∧
          HasDerivAt (fun s => φ (s, x)) (potX J (fun y => Y (t, y)) (ξ (t, x))) t) ∧
        (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
          horC U y₀ (Y (t, y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 → ξ (t, horC U y₀ (Y (t, y))) = y) := by
  classical
  have hhor : ContDiff ℝ ∞ (horC U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
      (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hver : ContDiff ℝ ∞ (verC J U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => verCL J U (v - y₀)) :=
      (verCL J U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hχ' : ContDiffOn ℝ ∞ (fun p : ℝ × 𝔼 => horC U y₀ (Y p))
      (Ioo c b ×ˢ ball (0 : 𝔼) 6) := hhor.comp_contDiffOn hY
  obtain ⟨ξ, hξ, hper⟩ := tu_inverse (fun p : ℝ × 𝔼 => horC U y₀ (Y p)) hχ' hχ
  have h56 : ball (0 : 𝔼) (11 / 2) ⊆ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num)
  set S : Set (ℝ × 𝔼) := Ioo c b ×ˢ ball (0 : 𝔼) 5 with hS
  have hSo : IsOpen S := isOpen_Ioo.prod isOpen_ball
  have hmaps : MapsTo (fun p : ℝ × 𝔼 => (p.1, ξ p)) S (Ioo c b ×ˢ ball (0 : 𝔼) 6) := by
    intro p hp
    exact ⟨hp.1, h56 ((hper p.1 hp.1).2.1 p.2 hp.2).1⟩
  set v : ℝ × 𝔼 → 𝔼 := fun p => verC J U y₀ (Y (p.1, ξ p)) with hv
  have hvC : ContDiffOn ℝ ∞ v S :=
    hver.comp_contDiffOn (hY.comp (contDiffOn_fst.prodMk hξ) hmaps)
  have hdec : ∀ t ∈ Ioo c b, ∀ x ∈ ball (0 : 𝔼) 5,
      Y (t, ξ (t, x)) = y₀ + U x + J (U (v (t, x))) := by
    intro t ht x hx
    have h := decomp_horC_verC (y₀ := y₀) hJi hJJ hU hdim (Y (t, ξ (t, x)))
    have h2 : horC U y₀ (Y (t, ξ (t, x))) = x := ((hper t ht).2.1 x hx).2
    rw [h2] at h
    exact h
  set φ₀ : ℝ × 𝔼 → ℝ := fun q => ∫ σ in (0 : ℝ)..1, ⟪v (q.1, σ • q.2), q.2⟫ with hφ
  have hφC : ContDiffOn ℝ ∞ φ₀ S := tg_pot_smooth isOpen_Ioo (convex_Ioo _ _) hvC
  have hgrad_all : ∀ t ∈ Ioo c b, ∀ x ∈ ball (0 : 𝔼) 5,
      Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ₀ (t, x')) x := by
    intro t ht x hx
    have hpos := (hper t ht).2.1
    have hvt : ContDiffOn ℝ ∞ (fun x => v (t, x)) (ball (0 : 𝔼) 5) := tg_slice hvC ht
    have hξt : ContDiffOn ℝ ∞ (fun x => ξ (t, x)) (ball (0 : 𝔼) 5) := tg_slice hξ ht
    have hYt : ContDiffOn ℝ ∞ (fun y => Y (t, y)) (ball (0 : 𝔼) 6) := tg_slice hY ht
    have hsym : ∀ x ∈ ball (0 : 𝔼) 5, ∀ a b : 𝔼,
        inner ℝ (fderiv ℝ (fun x => v (t, x)) x a) b =
          inner ℝ a (fderiv ℝ (fun x => v (t, x)) x b) := by
      intro x hx
      have hxi := (hpos x hx).1
      refine tg_symm hJi hJJ hU (y₀ := y₀) (X := fun y => Y (t, y)) (s := fun x => ξ (t, x))
        ?_ ?_ ?_ ?_ ?_
      · exact (hYt.contDiffAt (isOpen_ball.mem_nhds (h56 hxi))).differentiableAt (by simp)
      · exact (hξt.contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
      · exact (hvt.contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
      · filter_upwards [isOpen_ball.mem_nhds hx] with z hz using hdec t ht z hz
      · exact hlag t ht _ (h56 hxi)
    have hφt : ∀ x ∈ ball (0 : 𝔼) 5,
        HasFDerivAt (fun y => φ₀ (t, y)) (innerSL ℝ (v (t, x))) x := by
      intro x hx
      have := bp_hasFDerivAt_pot (x₀ := (0 : 𝔼)) hvt hsym hx
      have e : (fun y => φ₀ (t, y)) =
          fun y => ∫ σ in (0 : ℝ)..1, inner ℝ (v (t, 0 + σ • (y - 0))) (y - 0) := by
        funext y
        simp [hφ]
      rw [e]
      exact this
    have hgrad : ∀ x ∈ ball (0 : 𝔼) 5, gradE (fun x' => φ₀ (t, x')) x = v (t, x) := by
      intro x hx
      ext c
      show fderiv ℝ (fun x' => φ₀ (t, x')) x (eb c) = v (t, x) c
      rw [(hφt x hx).fderiv]
      simp [eb, EuclideanSpace.inner_single_right]
    rw [hdec t ht x hx]
    simp [graphX, hgrad x hx]
  have hξR : ∀ t ∈ Ioo c b, ∀ x ∈ ball (0 : 𝔼) 5, ξ (t, x) ∈ ball (0 : 𝔼) 6 :=
    fun t ht x hx => h56 ((hper t ht).2.1 x hx).1
  have h0 : ∀ t ∈ Ioo c b, φ₀ (t, 0) = 0 := by
    intro t _
    simp [hφ]
  have hdt := tf_dt_phi hJi hJJ hU hY hinj hlag hvel (isOpen_Ioo (a := c) (b := b))
    (Subset.refl _) hξ hξR hφC hgrad_all h0 (by norm_num : (0 : ℝ) < 5)
  -- the time-dependent constant
  have hO : IsOpen (Ioo c b ×ˢ ball (0 : 𝔼) 6) := isOpen_Ioo.prod isOpen_ball
  have hP := contDiffOn_potX_family (J := J) hO hY (fun q hq => hinj q.1 hq.1 q.2 hq.2)
  have hp0 : ContDiffOn ℝ ∞ (fun s => potX J (fun y => Y (s, y)) (ξ (s, 0))) (Ioo c b) := by
    have h1a : ContDiffOn ℝ ∞ (fun s : ℝ => ((s, (0 : 𝔼)) : ℝ × 𝔼)) (Ioo c b) :=
      contDiffOn_id.prodMk (contDiffOn_const (c := (0 : 𝔼)))
    have h1b : ContDiffOn ℝ ∞ (fun s : ℝ => ξ (s, 0)) (Ioo c b) :=
      hξ.comp h1a (fun s hs => ⟨hs, mem_ball_self (by norm_num)⟩)
    have h1 : ContDiffOn ℝ ∞ (fun s : ℝ => ((s, ξ (s, 0)) : ℝ × 𝔼)) (Ioo c b) :=
      contDiffOn_id.prodMk h1b
    have hmaps' : MapsTo (fun s : ℝ => ((s, ξ (s, 0)) : ℝ × 𝔼)) (Ioo c b)
        (Ioo c b ×ˢ ball (0 : 𝔼) 6) :=
      fun s hs => ⟨hs, hξR s hs 0 (mem_ball_self (by norm_num))⟩
    have hc := hP.comp h1 hmaps'
    refine hc.congr fun s hs => ?_
    simp only [Function.comp_apply]
  have ht₀ : (c + b) / 2 ∈ Ioo c b := ⟨by linarith, by linarith⟩
  obtain ⟨hcd, hcC⟩ := tf_primitive ht₀ hp0
  refine ⟨ξ, fun q => φ₀ q + ∫ s in (c + b) / 2..q.1, potX J (fun y => Y (s, y)) (ξ (s, 0)),
    hξ, ?_, ?_⟩
  · exact hφC.add (hcC.comp contDiffOn_fst (fun q hq => hq.1))
  intro t ht
  obtain ⟨hinj', hpos, hleft⟩ := hper t ht
  refine ⟨hinj', fun x hx => ⟨(hpos x hx).1, (hpos x hx).2, ?_, ?_⟩, hleft⟩
  · rw [hgrad_all t ht x hx]
    simp only [graphX, gradE, fderiv_add_const]
  · have h1 := hdt t ht x hx
    have h2 := h1.add (hcd t ht)
    exact h2.congr_deriv (by ring)

/-- info: 'CWFlow.Extension.exists_time_graph_unif' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_time_graph_unif

end CWFlow.Extension
