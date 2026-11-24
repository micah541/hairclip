import CWFlow.Extension.ChartGraphEvol
import CWFlow.Extension.LocalEvolution
import CWFlow.Extension.EnergyCompare
import Mathlib.MeasureTheory.Function.Jacobian

/-!
# The derivative of a cut-off energy in a fixed graph chart

For a flow `Φ` read in the graph chart at time `t₀`, the cut-off energy `∫ ρ² |∇^k C|²_g dV` is, for `t` near
`t₀`, a change of variables `χ_t` of the corresponding integral over the graph function `Ψ_t`
(`cd_cov`), and splits as the integral `∫ ρ² eSym k (jp Ψ_t)` handled by the local evolution inequality plus a
reparametrization term, whose derivative at `t₀` is a gauge term (`exists_const_chart_energy_deriv`).
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

section Lemmas

variable {n : ℕ}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

theorem cd_det_ne (L : 𝔼 →L[ℝ] 𝔼) (h : Function.Injective L) : L.det ≠ 0 := by
  unfold ContinuousLinearMap.det
  intro h0
  have := LinearMap.det_eq_zero_iff_ker_ne_bot.1 h0
  exact this (LinearMap.ker_eq_bot.2 h)

/-- The change of variables by `χ`. -/
theorem cd_cov (χ ξ' : 𝔼 → 𝔼) (ρ D E : 𝔼 → ℝ)
    (hχd : ∀ y ∈ ball (0 : 𝔼) (11 / 2), HasFDerivAt χ (fderiv ℝ χ y) y)
    (hχi : InjOn χ (ball (0 : 𝔼) (11 / 2)))
    (hpos : ∀ x ∈ ball (0 : 𝔼) 5, ξ' x ∈ ball (0 : 𝔼) (11 / 2) ∧ χ (ξ' x) = x)
    (hleft : ∀ y ∈ ball (0 : 𝔼) (11 / 2), χ y ∈ ball (0 : 𝔼) 5 → ξ' (χ y) = y)
    (hρ : ∀ y, y ∉ ball (0 : 𝔼) 2 → ρ y = 0)
    (hmap : ∀ y ∈ ball (0 : 𝔼) 2, χ y ∈ ball (0 : 𝔼) 5)
    (hD : ∀ y ∈ ball (0 : 𝔼) 2, D y = E (χ y) * |(fderiv ℝ χ y).det|) :
    ∫ y in ball (0 : 𝔼) 6, ρ y ^ 2 * D y = ∫ x in ball (0 : 𝔼) 5, ρ (ξ' x) ^ 2 * E x := by
  classical
  have h1 : ∫ y in ball (0 : 𝔼) 6, ρ y ^ 2 * D y = ∫ y in ball (0 : 𝔼) (11 / 2), ρ y ^ 2 * D y := by
    refine setIntegral_eq_of_subset_of_forall_sdiff_eq_zero measurableSet_ball
      (ball_subset_ball (by norm_num)) fun y hy => ?_
    have : y ∉ ball (0 : 𝔼) 2 := fun h2 => hy.2 (ball_subset_ball (by norm_num) h2)
    rw [hρ y this]; simp
  set g : 𝔼 → ℝ := (ball (0 : 𝔼) 5).indicator (fun x => ρ (ξ' x) ^ 2 * E x) with hg
  have h2 := integral_image_eq_integral_abs_det_fderiv_smul (volume : Measure 𝔼)
    (measurableSet_ball (x := (0 : 𝔼)) (ε := 11 / 2)) (fun y hy => (hχd y hy).hasFDerivWithinAt) hχi g
  have h3 : ∫ y in ball (0 : 𝔼) (11 / 2), ρ y ^ 2 * D y =
      ∫ y in ball (0 : 𝔼) (11 / 2), |(fderiv ℝ χ y).det| • g (χ y) := by
    refine setIntegral_congr_fun measurableSet_ball fun y hy => ?_
    simp only [smul_eq_mul, hg]
    by_cases h2y : y ∈ ball (0 : 𝔼) 2
    · rw [indicator_of_mem (hmap y h2y), hleft y hy (hmap y h2y), hD y h2y]; ring
    · rw [hρ y h2y]
      by_cases h5 : χ y ∈ ball (0 : 𝔼) 5
      · rw [indicator_of_mem h5, hleft y hy h5, hρ y h2y]; simp
      · rw [indicator_of_notMem h5]; simp
  have h4 : ∫ x in χ '' ball (0 : 𝔼) (11 / 2), g x = ∫ x, g x := by
    refine setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => ?_
    have : x ∉ ball (0 : 𝔼) 5 := fun h5 => hx ⟨ξ' x, (hpos x h5).1, (hpos x h5).2⟩
    simp [hg, this]
  have h5 : ∫ x, g x = ∫ x in ball (0 : 𝔼) 5, ρ (ξ' x) ^ 2 * E x :=
    integral_indicator measurableSet_ball
  rw [h1, h3, ← h2, h4, h5]

theorem cd_contDiff_horC {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (U : 𝔼 →ₗᵢ[ℝ] V) (c : V) : ContDiff ℝ ∞ (horC U c) := by
  have : horC U c = fun v : V => horCL U (v - c) := funext fun v => rfl
  rw [this]
  exact (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)

/-- The density of `Y = graph ∘ χ` near a point of the ball of radius `2`. -/
theorem cd_pointwise {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : 𝔼 →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (c : V) {Ψ' : 𝔼 → ℝ} (hΨ : ContDiff ℝ ∞ Ψ')
    {Y : 𝔼 → V} {ξ' : 𝔼 → 𝔼} (hY : ContDiffOn ℝ ∞ Y (ball (0 : 𝔼) (11 / 2)))
    (hder : ∀ y ∈ ball (0 : 𝔼) (11 / 2),
      Function.Injective (fderiv ℝ (fun y' => horC U c (Y y')) y))
    (hgr : ∀ x ∈ ball (0 : 𝔼) 4, Y (ξ' x) = graphX J U c Ψ' x)
    (hleft : ∀ y ∈ ball (0 : 𝔼) (11 / 2), horC U c (Y y) ∈ ball (0 : 𝔼) 5 →
      ξ' (horC U c (Y y)) = y)
    (hmap : ∀ y ∈ ball (0 : 𝔼) 2, horC U c (Y y) ∈ ball (0 : 𝔼) 4)
    {N : ℕ} (k : ℕ) (hN : k + 3 ≤ N) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 2) :
    eDens J Y k y = eSym k (jp N Ψ' (horC U c (Y y))) *
      |(fderiv ℝ (fun y' => horC U c (Y y')) y).det| := by
  set χ : 𝔼 → 𝔼 := fun y' => horC U c (Y y') with hχ
  have h211 : ball (0 : 𝔼) 2 ⊆ ball (0 : 𝔼) (11 / 2) := ball_subset_ball (by norm_num)
  have hev : Y =ᶠ[𝓝 y] (graphX J U c Ψ') ∘ χ := by
    filter_upwards [isOpen_ball.mem_nhds hy] with y' hy'
    have h4 := hmap y' hy'
    have h5 : χ y' ∈ ball (0 : 𝔼) 5 := ball_subset_ball (by norm_num) h4
    have := hgr (χ y') h4
    rw [hleft y' (h211 hy') h5] at this
    simpa using this
  have hχs : ContDiffOn ℝ ∞ χ (ball (0 : 𝔼) 2) :=
    ((cd_contDiff_horC U c).comp_contDiffOn (hY.mono h211))
  have hdet : ∀ z ∈ ball (0 : 𝔼) 2, (jacM χ z).det ≠ 0 := fun z hz => by
    rw [← im_det_eq_jacM]
    exact cd_det_ne _ (hder z (h211 hz))
  rw [eDens_congr hev k]
  rw [eDens_comp (Ω := Set.univ) (X := graphX J U c Ψ') isOpen_univ
    (contDiff_graphX hΨ).contDiffOn (fun x _ => injective_fderiv_graphX hJi hJJ hU hΨ x)
    (fun x _ u w => lagrangian_graphX hJi hJJ hU hΨ x u w) isOpen_ball hχs
    (fun _ _ => mem_univ _) hdet k hy, eDens_graphX hJi hJJ hU hΨ hN, ← im_det_eq_jacM]

/-- The time derivative of the reparametrization at the initial time. -/
theorem cd_dtxi {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {U : 𝔼 →ₗᵢ[ℝ] V} (c : V) {T : Set ℝ} (hT : IsOpen T) {t₀ : ℝ} (ht₀ : t₀ ∈ T)
    {Y : ℝ × 𝔼 → V} (hY : ContDiffOn ℝ ∞ Y (T ×ˢ ball (0 : 𝔼) 6))
    {ξ : ℝ × 𝔼 → 𝔼} (hξ : ContDiffOn ℝ ∞ ξ (T ×ˢ ball (0 : 𝔼) 5))
    (hξ0 : ∀ x ∈ ball (0 : 𝔼) 5, ξ (t₀, x) = x)
    (hpos : ∀ t ∈ T, ∀ x ∈ ball (0 : 𝔼) 5,
      ξ (t, x) ∈ ball (0 : 𝔼) (11 / 2) ∧ horC U c (Y (t, ξ (t, x))) = x)
    (hvel : ∀ y ∈ ball (0 : 𝔼) 6,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t₀, y')) y) t₀)
    (hgraph : ∀ y ∈ ball (0 : 𝔼) 6, horC U c (Y (t₀, y)) = y)
    {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) :
    fderiv ℝ ξ (t₀, x) (1, 0) = -horCL U (velX J (fun y' => Y (t₀, y')) x) := by
  have hx6 : x ∈ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num) hx
  have hOξ : IsOpen (T ×ˢ ball (0 : 𝔼) 5) := hT.prod isOpen_ball
  have hOY : IsOpen (T ×ˢ ball (0 : 𝔼) 6) := hT.prod isOpen_ball
  have hm5 : (t₀, x) ∈ T ×ˢ ball (0 : 𝔼) 5 := ⟨ht₀, hx⟩
  have hξc : ContDiffAt ℝ ∞ ξ (t₀, x) := hξ.contDiffAt (hOξ.mem_nhds hm5)
  have hξd : HasFDerivAt ξ (fderiv ℝ ξ (t₀, x)) (t₀, x) :=
    (hξc.differentiableAt (by simp)).hasFDerivAt
  set v : 𝔼 := fderiv ℝ ξ (t₀, x) (1, 0) with hv
  have hm6 : (t₀, x) ∈ T ×ˢ ball (0 : 𝔼) 6 := ⟨ht₀, hx6⟩
  have hYc : ContDiffAt ℝ ∞ Y (t₀, x) := hY.contDiffAt (hOY.mem_nhds hm6)
  have hYd : HasFDerivAt Y (fderiv ℝ Y (t₀, x)) (t₀, x) :=
    (hYc.differentiableAt (by simp)).hasFDerivAt
  -- the curve `s ↦ (s, ξ (s, x))`
  have hc1 : HasDerivAt (fun s : ℝ => (s, x)) ((1 : ℝ), (0 : 𝔼)) t₀ :=
    (hasDerivAt_id t₀).prodMk (hasDerivAt_const t₀ x)
  have hξs : HasDerivAt (fun s : ℝ => ξ (s, x)) v t₀ := hξd.comp_hasDerivAt t₀ hc1
  have hc2 : HasDerivAt (fun s : ℝ => (s, ξ (s, x))) ((1 : ℝ), v) t₀ :=
    (hasDerivAt_id t₀).prodMk hξs
  have hYd' : HasFDerivAt Y (fderiv ℝ Y (t₀, x)) (t₀, ξ (t₀, x)) := by
    rw [hξ0 x hx]; exact hYd
  have hγ : HasDerivAt (fun s : ℝ => Y (s, ξ (s, x))) (fderiv ℝ Y (t₀, x) (1, v)) t₀ :=
    hYd'.comp_hasDerivAt t₀ hc2
  have hhor : HasDerivAt (fun s : ℝ => horC U c (Y (s, ξ (s, x))))
      (horCL U (fderiv ℝ Y (t₀, x) (1, v))) t₀ :=
    (hasFDerivAt_horC (U := U) (y₀ := c) _).comp_hasDerivAt t₀ hγ
  have hconst : HasDerivAt (fun s : ℝ => horC U c (Y (s, ξ (s, x)))) 0 t₀ := by
    refine (hasDerivAt_const t₀ x).congr_of_eventuallyEq ?_
    filter_upwards [hT.mem_nhds ht₀] with s hs using (hpos s hs x hx).2
  have h0 : horCL U (fderiv ℝ Y (t₀, x) (1, v)) = 0 := hhor.unique hconst
  -- the two partial derivatives
  have hpt : fderiv ℝ Y (t₀, x) (1, 0) = velX J (fun y' => Y (t₀, y')) x := by
    have h1 : HasDerivAt (fun s : ℝ => Y (s, x)) (fderiv ℝ Y (t₀, x) (1, 0)) t₀ :=
      hYd.comp_hasDerivAt t₀ hc1
    exact h1.unique (hvel x hx6)
  have hspace : HasFDerivAt (fun y : 𝔼 => Y (t₀, y))
      ((fderiv ℝ Y (t₀, x)).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼)) x :=
    hYd.comp x (hasFDerivAt_prodMk_right t₀ x)
  have hsid : HasFDerivAt (fun y : 𝔼 => horC U c (Y (t₀, y)))
      ((horCL U).comp ((fderiv ℝ Y (t₀, x)).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼))) x :=
    (hasFDerivAt_horC (U := U) (y₀ := c) _).comp x hspace
  have hid : HasFDerivAt (fun y : 𝔼 => horC U c (Y (t₀, y))) (ContinuousLinearMap.id ℝ 𝔼) x := by
    refine (hasFDerivAt_id x).congr_of_eventuallyEq ?_
    filter_upwards [isOpen_ball.mem_nhds hx6] with y hy using hgraph y hy
  have hlin := hsid.unique hid
  have hvv : horCL U (fderiv ℝ Y (t₀, x) (0, v)) = v := by
    have := congrArg (fun L : 𝔼 →L[ℝ] 𝔼 => L v) hlin
    simpa using this
  have hsplit : fderiv ℝ Y (t₀, x) (1, v) =
      fderiv ℝ Y (t₀, x) (1, 0) + fderiv ℝ Y (t₀, x) (0, v) := by
    rw [← map_add]; simp
  rw [hsplit, map_add, hpt, hvv] at h0
  exact eq_neg_of_add_eq_zero_right h0

/-- A chart map close to the identity at the origin maps the ball of radius `2` into the ball of
radius `11/4`. -/
theorem cd_close {χ : 𝔼 → 𝔼} (hd : ∀ y ∈ ball (0 : 𝔼) (11 / 2), DifferentiableAt ℝ χ y)
    (hb : ∀ y ∈ ball (0 : 𝔼) (11 / 2), ‖fderiv ℝ χ y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4)
    (h0 : ‖χ 0‖ < 1 / 4) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 2) : ‖χ y‖ < 11 / 4 := by
  have h211 : ball (0 : 𝔼) 2 ⊆ ball (0 : 𝔼) (11 / 2) := ball_subset_ball (by norm_num)
  have := (convex_ball (0 : 𝔼) (11 / 2)).norm_image_sub_le_of_norm_fderiv_le' (f := χ)
    (φ := ContinuousLinearMap.id ℝ 𝔼) hd hb (mem_ball_self (by norm_num)) (h211 hy)
  simp only [sub_zero, ContinuousLinearMap.id_apply] at this
  have hy2 : ‖y‖ < 2 := by simpa using hy
  have h1 : ‖χ y‖ ≤ ‖χ y - χ 0 - y‖ + ‖χ 0‖ + ‖y‖ := by
    calc ‖χ y‖ = ‖(χ y - χ 0 - y) + χ 0 + y‖ := by congr 1; abel
      _ ≤ ‖(χ y - χ 0 - y) + χ 0‖ + ‖y‖ := norm_add_le _ _
      _ ≤ ‖χ y - χ 0 - y‖ + ‖χ 0‖ + ‖y‖ := by gcongr; exact norm_add_le _ _
  nlinarith [norm_nonneg y]

theorem cd_hasFDerivAt_of_contDiffOn {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    {f : 𝔼 → V} {s : Set 𝔼} (hs : IsOpen s) (hf : ContDiffOn ℝ ∞ f s) {y : 𝔼} (hy : y ∈ s) :
    HasFDerivAt f (fderiv ℝ f y) y :=
  ((hf.contDiffAt (hs.mem_nhds hy)).differentiableAt (by simp)).hasFDerivAt

/-- The smoothness of the reparametrization term. -/
theorem cd_G_smooth {ρ : 𝔼 → ℝ} (hρ : ContDiff ℝ ∞ ρ) {T : Set ℝ} {ξ : ℝ × 𝔼 → 𝔼}
    (hξ : ContDiffOn ℝ ∞ ξ (T ×ˢ ball (0 : 𝔼) 5)) {eF : ℝ × 𝔼 → ℝ}
    (heF : ContDiffOn ℝ ∞ eF (T ×ˢ ball (0 : 𝔼) 3)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eF q) (T ×ˢ ball (0 : 𝔼) 3) := by
  have h1 : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ρ (ξ q)) (T ×ˢ ball (0 : 𝔼) 3) :=
    hρ.comp_contDiffOn (hξ.mono (prod_mono subset_rfl (ball_subset_ball (by norm_num))))
  have h2 : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ρ q.2) (T ×ˢ ball (0 : 𝔼) 3) :=
    hρ.comp_contDiffOn contDiffOn_snd
  exact ((h1.pow 2).sub (h2.pow 2)).mul heF

/-- The time derivative of the reparametrization term at the initial time. -/
theorem cd_dG {ρ : 𝔼 → ℝ} (hρ : ContDiff ℝ ∞ ρ) {ξ : ℝ × 𝔼 → 𝔼} {eF : ℝ × 𝔼 → ℝ} {t₀ : ℝ}
    {x : 𝔼} (hξ0 : ξ (t₀, x) = x) (hξd : DifferentiableAt ℝ ξ (t₀, x))
    (hed : DifferentiableAt ℝ eF (t₀, x))
    (hGd : DifferentiableAt ℝ (fun q : ℝ × 𝔼 => (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eF q) (t₀, x)) :
    fderiv ℝ (fun q : ℝ × 𝔼 => (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eF q) (t₀, x) (1, 0) =
      fderiv ℝ (fun y => ρ y ^ 2) x (fderiv ℝ ξ (t₀, x) (1, 0)) * eF (t₀, x) := by
  have hc1 : HasDerivAt (fun s : ℝ => (s, x)) ((1 : ℝ), (0 : 𝔼)) t₀ :=
    (hasDerivAt_id t₀).prodMk (hasDerivAt_const t₀ x)
  have h1 := hGd.hasFDerivAt.comp_hasDerivAt t₀ hc1
  have hξs : HasDerivAt (fun s : ℝ => ξ (s, x)) (fderiv ℝ ξ (t₀, x) (1, 0)) t₀ :=
    hξd.hasFDerivAt.comp_hasDerivAt t₀ hc1
  have he : HasDerivAt (fun s : ℝ => eF (s, x)) (fderiv ℝ eF (t₀, x) (1, 0)) t₀ :=
    hed.hasFDerivAt.comp_hasDerivAt t₀ hc1
  have hρ2d : Differentiable ℝ (fun y : 𝔼 => ρ y ^ 2) := (hρ.pow 2).differentiable (by simp)
  have h2 : HasFDerivAt (fun y : 𝔼 => ρ y ^ 2) (fderiv ℝ (fun y : 𝔼 => ρ y ^ 2) x) (ξ (t₀, x)) := by
    rw [hξ0]; exact (hρ2d x).hasFDerivAt
  have hρξ := h2.comp_hasDerivAt t₀ hξs
  have hb := hρξ.sub_const (ρ x ^ 2)
  have hm := hb.mul he
  have := h1.unique hm
  rw [this]
  simp [hξ0]

theorem cd_fderiv_sq_zero {ρ : 𝔼 → ℝ} {x : 𝔼} (hx : x ∉ tsupport ρ) :
    fderiv ℝ (fun y => ρ y ^ 2) x = 0 := by
  have hev : (fun y => ρ y ^ 2) =ᶠ[𝓝 x] fun _ => (0 : ℝ) := by
    filter_upwards [(isClosed_tsupport ρ).isOpen_compl.mem_nhds hx] with y hy
    rw [image_eq_zero_of_notMem_tsupport hy]; simp
  rw [hev.fderiv_eq]; simp

theorem cd_match (L : 𝔼 →L[ℝ] ℝ) (v : 𝔼) (w : Fin n → ℝ) (hv : ∀ c, v c = -w c) (e : ℝ) :
    L v * e = -(e * ∑ c, L (EuclideanSpace.single c 1) * w c) := by
  have hrepr : v = ∑ c, v c • EuclideanSpace.single c (1 : ℝ) := by
    have := (EuclideanSpace.basisFun (Fin n) ℝ).sum_repr v
    simpa [EuclideanSpace.basisFun_apply] using this.symm
  conv_lhs => rw [hrepr]
  simp only [map_sum, map_smul, smul_eq_mul, Finset.sum_mul, Finset.mul_sum, hv]
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun c _ => ?_
  ring

/-- Splitting of the integral. -/
theorem cd_split (ρ ρ' e : 𝔼 → ℝ) (ξ' : 𝔼 → 𝔼)
    (hρ : ∀ y, y ∉ ball (0 : 𝔼) 2 → ρ y = 0)
    (hρ'1 : ∀ x ∈ closedBall (0 : 𝔼) (11 / 4), ρ' x = 1)
    (hρ'0 : ∀ x, x ∉ ball (0 : 𝔼) 3 → ρ' x = 0)
    (hvan : ∀ x ∈ ball (0 : 𝔼) 5, x ∉ closedBall (0 : 𝔼) (11 / 4) → ρ (ξ' x) = 0)
    (hi1 : IntegrableOn (fun x => ρ x ^ 2 * e x) (ball (0 : 𝔼) 5))
    (hi2 : IntegrableOn (fun x => ρ' x * ((ρ (ξ' x) ^ 2 - ρ x ^ 2) * e x)) (ball (0 : 𝔼) 5)) :
    ∫ x in ball (0 : 𝔼) 5, ρ (ξ' x) ^ 2 * e x =
      (∫ x, ρ x ^ 2 * e x) + ∫ x, ρ' x * ((ρ (ξ' x) ^ 2 - ρ x ^ 2) * e x) := by
  have h5 : ∀ x ∈ ball (0 : 𝔼) 5, ρ (ξ' x) ^ 2 * e x =
      ρ x ^ 2 * e x + ρ' x * ((ρ (ξ' x) ^ 2 - ρ x ^ 2) * e x) := by
    intro x hx
    by_cases hc : x ∈ closedBall (0 : 𝔼) (11 / 4)
    · rw [hρ'1 x hc]; ring
    · have hx2 : x ∉ ball (0 : 𝔼) 2 := fun h => hc (by
        have h1 : ‖x‖ < 2 := by simpa using h
        rw [mem_closedBall, dist_zero_right]; linarith)
      rw [hρ x hx2, hvan x hx hc]; simp
  rw [setIntegral_congr_fun measurableSet_ball h5, integral_add hi1 hi2]
  congr 1
  · refine setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => ?_
    have : x ∉ ball (0 : 𝔼) 2 := fun h => hx (ball_subset_ball (by norm_num) h)
    rw [hρ x this]; simp
  · refine setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => ?_
    have : x ∉ ball (0 : 𝔼) 3 := fun h => hx (ball_subset_ball (by norm_num) h)
    rw [hρ'0 x this]; simp

theorem cd_jet_eq {f g : 𝔼 → ℝ} (h : ∀ x ∈ ball (0 : 𝔼) 4, f x = g x) (m : ℕ) {x : 𝔼}
    (hx : x ∈ ball (0 : 𝔼) 4) : jet m f x = jet m g x := by
  refine jet_congr_of_eventuallyEq ?_ m
  filter_upwards [isOpen_ball.mem_nhds hx] with y hy using h y hy

theorem cd_energy3_congr {f g : 𝔼 → ℝ} (h : ∀ x ∈ ball (0 : 𝔼) 4, f x = g x) (m : ℕ) :
    energy3 m 3 f = energy3 m 3 g := by
  unfold energy3
  refine Finset.sum_congr rfl fun i _ => ?_
  refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
  rw [cd_jet_eq h _ (ball_subset_ball (by norm_num) hx)]

theorem cd_energy3_nonneg (m : ℕ) (f : 𝔼 → ℝ) : 0 ≤ energy3 m 3 f :=
  Finset.sum_nonneg fun i _ => setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _

theorem cd_integrals_congr {ρ f g : 𝔼 → ℝ} (hρ : ∀ y, y ∉ ball (0 : 𝔼) 2 → ρ y = 0)
    (h : ∀ x ∈ ball (0 : 𝔼) 4, f x = g x) (m : ℕ) :
    (∫ x, ρ x ^ 2 * ‖jet m f x‖ ^ 2) = ∫ x, ρ x ^ 2 * ‖jet m g x‖ ^ 2 := by
  refine integral_congr_ae (Filter.Eventually.of_forall fun x => ?_)
  by_cases hx : x ∈ ball (0 : 𝔼) 2
  · simp only
    rw [cd_jet_eq h m (ball_subset_ball (by norm_num) hx)]
  · simp only [hρ x hx]; simp

theorem cd_ball_congr {f g : 𝔼 → ℝ} (h : ∀ x ∈ ball (0 : 𝔼) 4, f x = g x) (m : ℕ) :
    (∫ x in ball (0 : 𝔼) 2, ‖jet m f x‖ ^ 2) = ∫ x in ball (0 : 𝔼) 2, ‖jet m g x‖ ^ 2 := by
  refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
  rw [cd_jet_eq h m (ball_subset_ball (by norm_num) hx)]

end Lemmas

/-- The derivative at `t₀` of the cut-off energy of a flow, computed in a fixed graph chart at time `t₀`,
with an upper bound: the local evolution inequality plus the gauge term of the reparametrization. -/
theorem exists_const_chart_energy_deriv {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n) (k : ℕ) (hk : 1 ≤ k) (Λ : ℝ) {ε : ℝ} (hε : 0 < ε)
    {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ {I : Set ℝ} {Φ : ℝ → M → V}, IsFlowM ι J I Φ → ∀ {t₀ : ℝ}, t₀ ∈ I →
      ∀ {p : M} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {σ : EuclideanSpace ℝ (Fin n) → M}
        {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ}, IsGraphChart ι J (Φ t₀) c₀ p U σ φ₀ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ₀ x‖ ≤ K₃) →
      ∀ {ρ : EuclideanSpace ℝ (Fin n) → ℝ}, IsCutoff ρ Λ →
      ∃ D : ℝ,
        HasDerivAt (fun s => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
          ρ x ^ 2 * eDens J (Φ s ∘ σ) k x) D t₀ ∧
        D ≤ -(2 * ((1 + (2 : ℝ) ^ 2)⁻¹) ^ (k + 5)) * (∫ x, ρ x ^ 2 * ‖jet (k + 5) φ₀ x‖ ^ 2)
          + ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ₀ x‖ ^ 2)
          + ε * energy3 (k + 2) 3 φ₀ + C * (1 + energy3 (k + 1) 3 φ₀) := by
  classical
  have hN4 : 4 ≤ k + 7 := by omega
  have hε2 : 0 < ε / 2 := by positivity
  obtain ⟨C₁, hC₁0, hC₁⟩ := exists_const_local_evolution (n := n) (N := k + 7) (k := k) hk le_rfl
    (K := 2) (by norm_num) Λ hε2
  obtain ⟨C₂, hC₂0, hC₂⟩ := exists_const_gauge_term (n := n) (N := k + 7) (k := k) hk le_rfl
    (K := 2) Λ hε2
  refine ⟨C₁ + C₂, add_nonneg hC₁0 hC₂0, ?_⟩
  intro I Φ hΦ t₀ ht₀ p U σ φ₀ hch hK ρ hρ
  obtain ⟨δ, hδ, hδI, ξ, Ψ, hGE, hξs, hΨ0, hξ0, hper⟩ :=
    exists_chart_graphEvol ι hJi hJJ hdim (k + 7) hN4 hc₀ hc₀1 hK₃ hΦ ht₀ hch hK
  have ht₀T : t₀ ∈ Ioo (t₀ - δ) (t₀ + δ) := ⟨by linarith, by linarith⟩
  have hYI : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2))
      (I ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    cge_contDiffOn_family ι hΦ.smooth hch.param hΦ.isOpen
  have hYT : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2))
      (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    hYI.mono (prod_mono hδI subset_rfl)
  -- smallness of the chart map at the origin
  obtain ⟨δ', hδ', hδ'δ, hsmall⟩ : ∃ δ' : ℝ, 0 < δ' ∧ δ' ≤ δ ∧ ∀ t ∈ Ioo (t₀ - δ') (t₀ + δ'),
      ‖horC U (Φ t₀ p) (Φ t (σ 0))‖ < 1 / 4 := by
    have hd := cge_hasDerivAt_family ι hΦ hch.param ht₀ (y := 0) (mem_ball_self (by norm_num))
    have hcont : ContinuousAt (fun s => horC U (Φ t₀ p) (Φ s (σ 0))) t₀ :=
      (cd_contDiff_horC U (Φ t₀ p)).continuous.continuousAt.comp hd.continuousAt
    have hval : horC U (Φ t₀ p) (Φ t₀ (σ 0)) = 0 :=
      (hch.graph 0 (mem_ball_self (by norm_num))).2.1
    have hev : ∀ᶠ s in 𝓝 t₀, ‖horC U (Φ t₀ p) (Φ s (σ 0))‖ < 1 / 4 :=
      hcont.norm.eventually_lt continuousAt_const (by rw [hval]; norm_num)
    obtain ⟨e, he, hball⟩ := Metric.eventually_nhds_iff.1 hev
    refine ⟨min δ e, lt_min hδ he, min_le_left _ _, fun t ht => hball ?_⟩
    rw [Real.dist_eq, abs_lt]
    constructor <;> linarith [ht.1, ht.2, min_le_right δ e]
  have hT'T : Ioo (t₀ - δ') (t₀ + δ') ⊆ Ioo (t₀ - δ) (t₀ + δ) :=
    Ioo_subset_Ioo (by linarith) (by linarith)
  have hT'mem : Ioo (t₀ - δ') (t₀ + δ') ∈ 𝓝 t₀ := Ioo_mem_nhds (by linarith) (by linarith)
  -- the local evolution inequality
  obtain ⟨D₁, hD₁, hD₁le⟩ := hC₁ (Ioo (t₀ - δ) (t₀ + δ)) Ψ hGE ρ hρ t₀ ht₀T
  -- a cutoff equal to one on the ball of radius `11/4`
  obtain ⟨ρ', hρ'c, hρ'cs, hρ'1, hρ'ts⟩ : ∃ ρ' : EuclideanSpace ℝ (Fin n) → ℝ, Continuous ρ' ∧
      HasCompactSupport ρ' ∧ (∀ x ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) (11 / 4), ρ' x = 1) ∧
      tsupport ρ' ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 3 := by
    let b : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨11 / 4, 29 / 10, by norm_num, by norm_num⟩
    refine ⟨b, b.continuous, b.hasCompactSupport, fun x hx => b.one_of_mem_closedBall hx, ?_⟩
    rw [b.tsupport_eq]
    exact closedBall_subset_ball (by show (29 : ℝ) / 10 < 3; norm_num)
  have hρ'0 : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 3 → ρ' x = 0 := fun x hx =>
    image_eq_zero_of_notMem_tsupport (fun h => hx (hρ'ts h))
  have hρ0 : ∀ y, y ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → ρ y = 0 := fun y hy =>
    image_eq_zero_of_notMem_tsupport (fun h => hy (hρ.supp h))
  have hOT3 : IsOpen (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    isOpen_Ioo.prod isOpen_ball
  have hOT5 : IsOpen (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5) :=
    isOpen_Ioo.prod isOpen_ball
  -- smoothness of the energy density along the graph family
  have hA := isOpen_admSet (n := n) (N := k + 7) (2 : ℝ)
  have hlow := admSet_low (n := n) (N := k + 7) (2 : ℝ)
  have heF : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2))
      (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    (contDiffOn_symbol_slice (isSymbol_eSym hA hlow k (by omega)) hGE.isOpen hGE.smooth).mono
      (fun q hq => ⟨⟨hq.1, mem_univ _⟩, hGE.adm q.1 hq.1 q.2 hq.2⟩)
  have hG := cd_G_smooth hρ.smooth hξs heF
  have hB : HasDerivAt (fun s : ℝ => ∫ x, ρ' x * ((ρ (ξ (s, x)) ^ 2 - ρ x ^ 2) *
        eSym k (jp (k + 7) (fun y => Ψ (s, y)) x)))
      (∫ x, ρ' x * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2)) (t₀, x) (1, 0)) t₀ :=
    hasDerivAt_integral_mul_slice hOT3 hG hρ'c hρ'cs (fun x hx => ⟨ht₀T, hρ'ts hx⟩)
  have hΨ₀ : ContDiff ℝ ∞ (fun y => Ψ (t₀, y)) := contDiff_slice hGE.isOpen hGE.smooth ht₀T
  -- the value of the derivative of the reparametrization term
  have hpt : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2)) (t₀, x) (1, 0) =
      -(eSym k (jp (k + 7) (fun y => Ψ (t₀, y)) x) * ∑ c : Fin n,
        fderiv ℝ (fun y => ρ y ^ 2) x (EuclideanSpace.single c 1) *
          vSym (jp (k + 7) (fun y => Ψ (t₀, y)) x) c) := by
    intro x hx
    have hx5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) hx
    have hm5 : (t₀, x) ∈ Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
      ⟨ht₀T, hx5⟩
    have hm3 : (t₀, x) ∈ Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3 :=
      ⟨ht₀T, hx⟩
    have hξd : DifferentiableAt ℝ ξ (t₀, x) :=
      (hξs.contDiffAt (hOT5.mem_nhds hm5)).differentiableAt (by simp)
    have hed : DifferentiableAt ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2)) (t₀, x) :=
      (heF.contDiffAt (hOT3.mem_nhds hm3)).differentiableAt (by simp)
    have hGd : DifferentiableAt ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2)) (t₀, x) :=
      (hG.contDiffAt (hOT3.mem_nhds hm3)).differentiableAt (by simp)
    rw [cd_dG hρ.smooth (hξ0 x hx5) hξd hed hGd]
    have hdt := cd_dtxi (J := J) (U := U) (Φ t₀ p) isOpen_Ioo ht₀T
      (Y := fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2)) hYT hξs hξ0
      (fun t ht x hx => ((hper t ht).2.2.1 x hx))
      (fun y hy => cge_hasDerivAt_family ι hΦ hch.param ht₀ hy)
      (fun y hy => (hch.graph y hy).2.1) hx5
    have hvel : velX J (fun y' => Φ t₀ (σ y')) x =
        velX J (graphX J U (Φ t₀ p) (fun y => Ψ (t₀, y))) x := by
      refine cge_velX_congr ?_
      filter_upwards [isOpen_ball.mem_nhds (ball_subset_ball (by norm_num : (3 : ℝ) ≤ 4) hx)]
        with y hy
      have hy6 : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := ball_subset_ball (by norm_num) hy
      rw [(hch.graph y hy6).1]
      refine cge_graphX_congr ?_
      filter_upwards [isOpen_ball.mem_nhds hy] with z hz using (hΨ0 z hz).symm
    have hhor := horCL_velX_graphX hJi hJJ hch.frame (Φ t₀ p) hΨ₀ (by omega : 5 ≤ k + 7) x
      ⟨2, hGE.adm t₀ ht₀T x hx⟩
    refine cd_match _ _ (fun c => vSym (jp (k + 7) (fun y => Ψ (t₀, y)) x) c) (fun c => ?_) _
    rw [hdt, hvel]
    simp [hhor c]
  have hBd : (∫ x, ρ' x * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        (ρ (ξ q) ^ 2 - ρ q.2 ^ 2) * eSym k (jp (k + 7) (fun y => Ψ (q.1, y)) q.2)) (t₀, x) (1, 0)) =
      -(∫ x, eSym k (jp (k + 7) (fun y => Ψ (t₀, y)) x) * ∑ c : Fin n,
        fderiv ℝ (fun y => ρ y ^ 2) x (EuclideanSpace.single c 1) *
          vSym (jp (k + 7) (fun y => Ψ (t₀, y)) x) c) := by
    rw [← integral_neg]
    refine integral_congr_ae (Filter.Eventually.of_forall fun x => ?_)
    by_cases hx3 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3
    · simp only
      rw [hpt x hx3]
      by_cases hxs : x ∈ tsupport ρ
      · rw [hρ'1 x (ball_subset_closedBall.trans (closedBall_subset_closedBall (by norm_num))
          (hρ.supp hxs))]
        ring
      · rw [cd_fderiv_sq_zero hxs]; simp
    · simp only
      rw [hρ'0 x hx3]
      have : x ∉ tsupport ρ := fun h => hx3 (ball_subset_ball (by norm_num) (hρ.supp h))
      rw [cd_fderiv_sq_zero this]; simp
  -- the gauge term estimate
  have hg := hC₂ ρ hρ (fun y => Ψ (t₀, y)) hΨ₀ (fun x hx => hGE.adm t₀ ht₀T x hx)
  -- the identity of the integrals for `t` near `t₀`
  have hρ2c : Continuous fun x : EuclideanSpace ℝ (Fin n) => ρ x ^ 2 := hρ.smooth.continuous.pow 2
  have hρ2cs : HasCompactSupport fun x : EuclideanSpace ℝ (Fin n) => ρ x ^ 2 :=
    (ibb_hasCompactSupport hρ).comp_left (g := fun a : ℝ => a ^ 2) (by simp)
  have hts2 : tsupport (fun x : EuclideanSpace ℝ (Fin n) => ρ x ^ 2) ⊆ tsupport ρ :=
    closure_mono (fun x hx => by simpa using hx)
  have hev : ∀ s ∈ Ioo (t₀ - δ') (t₀ + δ'),
      (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, ρ x ^ 2 * eDens J (Φ s ∘ σ) k x) =
        (∫ x, ρ x ^ 2 * eSym k (jp (k + 7) (fun y => Ψ (s, y)) x)) +
        ∫ x, ρ' x * ((ρ (ξ (s, x)) ^ 2 - ρ x ^ 2) * eSym k (jp (k + 7) (fun y => Ψ (s, y)) x)) := by
    intro s hs
    have hsT := hT'T hs
    obtain ⟨hinj, hder, hpos, hgr, hleft⟩ := hper s hsT
    have hΨs : ContDiff ℝ ∞ (fun y => Ψ (s, y)) := contDiff_slice hGE.isOpen hGE.smooth hsT
    have hYs : ContDiffOn ℝ ∞ (fun y => Φ s (σ y)) (ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
      hYT.comp (contDiff_const.prodMk contDiff_id).contDiffOn (fun y hy => ⟨hsT, hy⟩)
    have h116 : ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 6 :=
      ball_subset_ball (by norm_num)
    have hYs11 := hYs.mono h116
    have hχsm : ContDiffOn ℝ ∞ (fun y => horC U (Φ t₀ p) (Φ s (σ y)))
        (ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2)) :=
      (cd_contDiff_horC U (Φ t₀ p)).comp_contDiffOn hYs11
    have hχd : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
        HasFDerivAt (fun y' => horC U (Φ t₀ p) (Φ s (σ y')))
          (fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ s (σ y'))) y) y := fun y hy =>
      cd_hasFDerivAt_of_contDiffOn isOpen_ball hχsm hy
    have hclose : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
        ‖horC U (Φ t₀ p) (Φ s (σ y))‖ < 11 / 4 := fun y hy =>
      cd_close (χ := fun y' => horC U (Φ t₀ p) (Φ s (σ y'))) (fun y hy => (hχd y hy).differentiableAt)
        (fun y hy => (hder y hy).2) (hsmall s hs) hy
    have hmap4 : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
        horC U (Φ t₀ p) (Φ s (σ y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 4 := fun y hy => by
      rw [mem_ball_zero_iff]; linarith [hclose y hy]
    have hmap5 : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
        horC U (Φ t₀ p) (Φ s (σ y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := fun y hy => by
      rw [mem_ball_zero_iff]; linarith [hclose y hy]
    have hD : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
        eDens J (Φ s ∘ σ) k y = eSym k (jp (k + 7) (fun y => Ψ (s, y)) (horC U (Φ t₀ p) (Φ s (σ y)))) *
          |(fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ s (σ y'))) y).det| := fun y hy =>
      cd_pointwise hJi hJJ hch.frame (Φ t₀ p) hΨs (Y := fun y => Φ s (σ y))
        (ξ' := fun x => ξ (s, x)) hYs11 (fun y hy => (hder y hy).1) hgr hleft hmap4 k
        (N := k + 7) (by omega) hy
    have hcov := cd_cov (fun y => horC U (Φ t₀ p) (Φ s (σ y))) (fun x => ξ (s, x)) ρ
      (fun y => eDens J (Φ s ∘ σ) k y)
      (fun x => eSym k (jp (k + 7) (fun y => Ψ (s, y)) x))
      hχd hinj (fun x hx => hpos x hx) (fun y hy h5 => hleft y hy h5) hρ0 hmap5 hD
    have hvan : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
        x ∉ closedBall (0 : EuclideanSpace ℝ (Fin n)) (11 / 4) → ρ (ξ (s, x)) = 0 := by
      intro x hx hxc
      by_contra hne
      have hmem : ξ (s, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2 := hρ.supp (subset_tsupport ρ hne)
      have h1 := hclose _ hmem
      rw [(hpos x hx).2] at h1
      exact hxc (by rw [mem_closedBall, dist_zero_right]; exact h1.le)
    have hi1 : IntegrableOn (fun x => ρ x ^ 2 * eSym k (jp (k + 7) (fun y => Ψ (s, y)) x))
        (ball (0 : EuclideanSpace ℝ (Fin n)) 5) := by
      have hc : Continuous fun x => ρ x ^ 2 * eSym k (jp (k + 7) (fun y => Ψ (s, y)) x) :=
        slc_continuous_mul (ρ := fun x => ρ x ^ 2) hOT3 heF.continuousOn hρ2c
          (fun x hx => ⟨hsT, ball_subset_ball (by norm_num) (hρ.supp (hts2 hx))⟩)
      exact (hc.integrable_of_hasCompactSupport hρ2cs.mul_right).integrableOn
    have hi2 : IntegrableOn (fun x => ρ' x * ((ρ (ξ (s, x)) ^ 2 - ρ x ^ 2) *
        eSym k (jp (k + 7) (fun y => Ψ (s, y)) x))) (ball (0 : EuclideanSpace ℝ (Fin n)) 5) := by
      have hc : Continuous fun x => ρ' x * ((ρ (ξ (s, x)) ^ 2 - ρ x ^ 2) *
          eSym k (jp (k + 7) (fun y => Ψ (s, y)) x)) :=
        slc_continuous_mul hOT3 hG.continuousOn hρ'c (fun x hx => ⟨hsT, hρ'ts hx⟩)
      exact (hc.integrable_of_hasCompactSupport hρ'cs.mul_right).integrableOn
    exact hcov.trans (cd_split ρ ρ' (fun x => eSym k (jp (k + 7) (fun y => Ψ (s, y)) x))
      (fun x => ξ (s, x)) hρ0 hρ'1 hρ'0 hvan hi1 hi2)
  have hHD := (hD₁.add hB).congr_of_eventuallyEq
    (show (fun s => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, ρ x ^ 2 * eDens J (Φ s ∘ σ) k x)
      =ᶠ[𝓝 t₀] _ by filter_upwards [hT'mem] with s hs using hev s hs)
  refine ⟨_, hHD, ?_⟩
  -- the bound
  have e1 := cd_integrals_congr (f := fun y => Ψ (t₀, y)) (g := φ₀) hρ0 hΨ0 (k + 5)
  have e2 := cd_ball_congr (f := fun y => Ψ (t₀, y)) (g := φ₀) hΨ0 (k + 5)
  have e3 := cd_energy3_congr (f := fun y => Ψ (t₀, y)) (g := φ₀) hΨ0 (k + 2)
  have e4 := cd_energy3_congr (f := fun y => Ψ (t₀, y)) (g := φ₀) hΨ0 (k + 1)
  have hI2 : 0 ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ₀ x‖ ^ 2 :=
    setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _
  have hE1 := cd_energy3_nonneg (n := n) (k + 1) φ₀
  rw [e1, e2, e3, e4] at hD₁le
  rw [e3] at hg
  have hBdle := (abs_le.1 hg).1
  rw [hBd]
  nlinarith [mul_nonneg hε.le hI2, mul_nonneg hC₂0 hE1]

/-- info: 'CWFlow.Extension.exists_const_chart_energy_deriv' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_chart_energy_deriv

end CWFlow.Extension
