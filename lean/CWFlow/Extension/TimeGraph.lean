import CWFlow.Extension.GraphLocal
import CWFlow.Interp.ParamInverse
import ShortTime.Analysis.ParamIntegral

/-!
# A smooth family of Lagrangian immersions is a family of graphs

A smooth family `Y_t` of Lagrangian immersions of a ball, which at time `t₀` is parametrized as a graph
over the Lagrangian plane (`horC ∘ Y_{t₀} = id`), is for `t` near `t₀` a graph, after the
reparametrization `ξ_t = (horC ∘ Y_t)⁻¹`, of the gradient of a function `φ(t, ·)` which is jointly smooth
in `(t, x)` (`exists_time_graph`). The potential is `φ(t, x) = ∫₀¹ ⟪v(t, σ x), x⟫ dσ` and its joint
smoothness comes from `ShortTime.contDiffOn_integral_Icc` (MM; Mathlib-only imports).
-/

open Set Metric
open CWFlow.Interp
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-- A slice `x ↦ g (t, x)` of a smooth function on a product is smooth. -/
theorem tg_slice {G : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G] {g : ℝ × 𝔼 → G}
    {A : Set ℝ} {B : Set 𝔼} (hg : ContDiffOn ℝ ∞ g (A ×ˢ B)) {t : ℝ} (ht : t ∈ A) :
    ContDiffOn ℝ ∞ (fun x => g (t, x)) B :=
  hg.comp (contDiffOn_const.prodMk contDiffOn_id) (fun x hx => ⟨ht, hx⟩)

/-- The derivative of the vertical coordinate of a reparametrized Lagrangian immersion is symmetric. -/
theorem tg_symm (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) {y₀ : V} {X : 𝔼 → V} {s v : 𝔼 → 𝔼} {x : 𝔼}
    (hXd : DifferentiableAt ℝ X (s x)) (hsd : DifferentiableAt ℝ s x)
    (hvd : DifferentiableAt ℝ v x)
    (hdec : ∀ᶠ z in 𝓝 x, X (s z) = y₀ + U z + J (U (v z)))
    (hl : ∀ u w, ⟪J (fderiv ℝ X (s x) u), fderiv ℝ X (s x) w⟫ = 0) :
    ∀ a b : 𝔼, ⟪fderiv ℝ v x a, b⟫ = ⟪a, fderiv ℝ v x b⟫ := by
  intro a b
  have h1 : HasFDerivAt (fun z => X (s z)) ((fderiv ℝ X (s x)).comp (fderiv ℝ s x)) x :=
    hXd.hasFDerivAt.comp x hsd.hasFDerivAt
  have h2 : HasFDerivAt (fun z => y₀ + U z + J (U (v z)))
      (U.toContinuousLinearMap + (J.comp U.toContinuousLinearMap).comp (fderiv ℝ v x)) x :=
    ((U.toContinuousLinearMap.hasFDerivAt (x := x)).const_add y₀).add
      ((J.comp U.toContinuousLinearMap).hasFDerivAt.comp x hvd.hasFDerivAt)
  have h3 : HasFDerivAt (fun z => X (s z))
      (U.toContinuousLinearMap + (J.comp U.toContinuousLinearMap).comp (fderiv ℝ v x)) x :=
    h2.congr_of_eventuallyEq hdec
  have h4 := h1.unique h3
  have hl' := hl (fderiv ℝ s x a) (fderiv ℝ s x b)
  have e : ∀ w, fderiv ℝ X (s x) (fderiv ℝ s x w) = U w + J (U (fderiv ℝ v x w)) := fun w => by
    have := congrArg (fun L => L w) h4
    simpa using this
  rw [e, e] at hl'
  simp only [map_add, hJJ, inner_add_left, inner_add_right, inner_neg_left, gb_JUU hU,
    gb_JUJU hJi, gb_UJU hU, gb_UU, neg_zero, add_zero, zero_add] at hl'
  linarith

/-- The potential `∫₀¹ ⟪v (t, σ x), x⟫ dσ` of a smooth family of fields is jointly smooth. -/
theorem tg_pot_smooth {A : Set ℝ} {r : ℝ} (hA : IsOpen A) (hAc : Convex ℝ A) {v : ℝ × 𝔼 → 𝔼}
    (hvC : ContDiffOn ℝ ∞ v (A ×ˢ ball (0 : 𝔼) r)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ∫ σ in (0 : ℝ)..1, ⟪v (q.1, σ • q.2), q.2⟫)
      (A ×ˢ ball (0 : 𝔼) r) := by
  refine ShortTime.contDiffOn_integral_Icc
    (f := fun (q : ℝ × 𝔼) (σ : ℝ) => ⟪v (q.1, σ • q.2), q.2⟫) (hAc.prod (convex_ball _ _))
    (hA.prod isOpen_ball).uniqueDiffOn ?_
  have hg : ContDiff ℝ ∞ (fun z : (ℝ × 𝔼) × ℝ => ((z.1.1, z.2 • z.1.2) : ℝ × 𝔼)) := by
    fun_prop
  have hm : MapsTo (fun z : (ℝ × 𝔼) × ℝ => ((z.1.1, z.2 • z.1.2) : ℝ × 𝔼))
      ((A ×ˢ ball (0 : 𝔼) r) ×ˢ Icc (0 : ℝ) 1) (A ×ˢ ball (0 : 𝔼) r) := by
    intro z hz
    have h1 : z.1.2 ∈ ball (0 : 𝔼) r := hz.1.2
    have h2 := bp_mem_ball (x₀ := (0 : 𝔼)) (R := r) (t := z.2) h1 hz.2
    refine ⟨hz.1.1, ?_⟩
    simpa using h2
  have h1 := hvC.comp hg.contDiffOn hm
  have h2 : ContDiff ℝ ∞ (fun z : (ℝ × 𝔼) × ℝ => z.1.2) := by fun_prop
  exact h1.inner ℝ h2.contDiffOn

theorem exists_time_graph {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) (y₀ : V)
    {a b t₀ : ℝ} (ht₀ : t₀ ∈ Ioo a b) {R : ℝ} {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) R))
    (hlag : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hgraph : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, horC U y₀ (Y (t₀, y)) = y)
    {r R' : ℝ} (hr : 0 < r) (hrR' : r < R') (hR'R : R' < R) :
    ∃ δ > 0, Ioo (t₀ - δ) (t₀ + δ) ⊆ Ioo a b ∧
      ∃ (ξ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
        (φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
        ContDiffOn ℝ ∞ ξ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) ∧
        ContDiffOn ℝ ∞ φ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) ∧
        ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ),
          InjOn (fun y => horC U y₀ (Y (t, y))) (ball (0 : EuclideanSpace ℝ (Fin n)) R') ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R',
            Function.Injective (fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y) ∧
            ‖fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4) ∧
          (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r,
            ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' ∧ horC U y₀ (Y (t, ξ (t, x))) = x ∧
            Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ (t, x')) x) ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R',
            horC U y₀ (Y (t, y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r →
              ξ (t, horC U y₀ (Y (t, y))) = y) ∧
          φ (t, 0) = 0 := by
  classical
  have hhor : ContDiff ℝ ∞ (horC U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
      (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hver : ContDiff ℝ ∞ (verC J U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => verCL J U (v - y₀)) :=
      (verCL J U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hχ : ContDiffOn ℝ ∞ (fun p : ℝ × 𝔼 => horC U y₀ (Y p))
      (Ioo a b ×ˢ ball (0 : 𝔼) R) := hhor.comp_contDiffOn hY
  obtain ⟨δ, hδ, hsub, ξ, hξ, hper⟩ := exists_param_inverse (χ := fun p : ℝ × 𝔼 => horC U y₀ (Y p))
    ht₀ hχ hgraph hr hrR' hR'R
  set I : Set ℝ := Ioo (t₀ - δ) (t₀ + δ) with hI
  set S : Set (ℝ × 𝔼) := I ×ˢ ball (0 : 𝔼) r with hS
  have hSo : IsOpen S := isOpen_Ioo.prod isOpen_ball
  have hmaps : MapsTo (fun p : ℝ × 𝔼 => (p.1, ξ p)) S (Ioo a b ×ˢ ball (0 : 𝔼) R) := by
    intro p hp
    exact ⟨hsub hp.1, ball_subset_ball hR'R.le ((hper p.1 hp.1).2.2.1 p.2 hp.2).1⟩
  set v : ℝ × 𝔼 → 𝔼 := fun p => verC J U y₀ (Y (p.1, ξ p)) with hv
  have hvC : ContDiffOn ℝ ∞ v S :=
    hver.comp_contDiffOn (hY.comp (contDiffOn_fst.prodMk hξ) hmaps)
  have hdec : ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r,
      Y (t, ξ (t, x)) = y₀ + U x + J (U (v (t, x))) := by
    intro t ht x hx
    have h := decomp_horC_verC (y₀ := y₀) hJi hJJ hU hdim (Y (t, ξ (t, x)))
    have h2 : horC U y₀ (Y (t, ξ (t, x))) = x := ((hper t ht).2.2.1 x hx).2
    rw [h2] at h
    exact h
  -- the potential
  set φ : ℝ × 𝔼 → ℝ := fun q => ∫ σ in (0 : ℝ)..1, ⟪v (q.1, σ • q.2), q.2⟫ with hφ
  have hφC : ContDiffOn ℝ ∞ φ S :=
    tg_pot_smooth isOpen_Ioo (convex_Ioo _ _) hvC
  refine ⟨δ, hδ, hsub, ξ, φ, hξ, hφC, ?_⟩
  intro t ht
  obtain ⟨hinj, hder, hpos, hleft⟩ := hper t ht
  have htab : t ∈ Ioo a b := hsub ht
  -- the slice of `v`
  have hvt : ContDiffOn ℝ ∞ (fun x => v (t, x)) (ball (0 : 𝔼) r) := tg_slice hvC ht
  have hξt : ContDiffOn ℝ ∞ (fun x => ξ (t, x)) (ball (0 : 𝔼) r) := tg_slice hξ ht
  have hYt : ContDiffOn ℝ ∞ (fun y => Y (t, y)) (ball (0 : 𝔼) R) := tg_slice hY htab
  have hsym : ∀ x ∈ ball (0 : 𝔼) r, ∀ a b : 𝔼,
      inner ℝ (fderiv ℝ (fun x => v (t, x)) x a) b = inner ℝ a (fderiv ℝ (fun x => v (t, x)) x b) := by
    intro x hx
    have hxi := (hpos x hx).1
    refine tg_symm hJi hJJ hU (y₀ := y₀) (X := fun y => Y (t, y)) (s := fun x => ξ (t, x))
      ?_ ?_ ?_ ?_ ?_
    · exact (hYt.contDiffAt (isOpen_ball.mem_nhds (ball_subset_ball hR'R.le hxi))).differentiableAt
        (by simp)
    · exact (hξt.contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
    · exact (hvt.contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
    · filter_upwards [isOpen_ball.mem_nhds hx] with z hz using hdec t ht z hz
    · exact hlag t htab _ (ball_subset_ball hR'R.le hxi)
  have hφt : ∀ x ∈ ball (0 : 𝔼) r, HasFDerivAt (fun y => φ (t, y)) (innerSL ℝ (v (t, x))) x := by
    intro x hx
    have := bp_hasFDerivAt_pot (x₀ := (0 : 𝔼)) hvt hsym hx
    have e : (fun y => φ (t, y)) =
        fun y => ∫ σ in (0 : ℝ)..1, inner ℝ (v (t, 0 + σ • (y - 0))) (y - 0) := by
      funext y
      simp [hφ]
    rw [e]
    exact this
  have hgrad : ∀ x ∈ ball (0 : 𝔼) r, gradE (fun x' => φ (t, x')) x = v (t, x) := by
    intro x hx
    ext c
    show fderiv ℝ (fun x' => φ (t, x')) x (eb c) = v (t, x) c
    rw [(hφt x hx).fderiv]
    simp [eb, EuclideanSpace.inner_single_right]
  refine ⟨hinj, fun y hy => ⟨?_, ?_⟩, fun x hx => ⟨(hpos x hx).1, (hpos x hx).2, ?_⟩, hleft, ?_⟩
  · exact (hder y hy).1
  · exact (hder y hy).2
  · rw [hdec t ht x hx]
    simp [graphX, hgrad x hx]
  · simp [hφ]

/-- info: 'CWFlow.Extension.exists_time_graph' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_time_graph

end CWFlow.Extension
