import CWFlow.Extension.GraphChart
import CWFlow.Extension.EnergyCompare
import CWFlow.Extension.ChartEnergy
import CWFlow.Extension.Slope

/-!
# Bounds in a graph chart

For a graph chart `φ` of a Lagrangian immersion `F` of a compact manifold:

* `graphChart_jet3_le`: `|D³φ| ≤ (1 + c₀²)^{3/2} c₁` on the ball of radius `6`, where `c₁` bounds the
  second fundamental form of `F` in all charts;
* `graphChart_adm`: the jets of `φ` of order `≤ N` are admissible (all slots of order `≤ 3` of norm `< 2`) on
  the ball of radius `5`;
* `exists_const_energy3_le_energyM`: the `L²` norms on `ball 0 3` of `D^3 φ, …, D^{m+3} φ` are bounded by
  `C (1 + Σ_{l ≤ m} E_l(F))` (the chart part of Lemma 4.13 of Chen–Warren, combined with the comparison of
  chart integrals with the global energy).
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

local notation "𝔼" => EuclideanSpace ℝ (Fin _)

/-! ### Low jets -/

theorem cb_norm_jet_one {n : ℕ} (f : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    ‖jet 1 f x‖ = ‖gradE f x‖ := by
  have h : ‖jet 1 f x‖ ^ 2 = ‖gradE f x‖ ^ 2 := by
    rw [norm_jet_sq, EuclideanSpace.norm_sq_eq]
    refine Fintype.sum_equiv (Equiv.funUnique (Fin 1) (Fin n)) _ _ fun I => ?_
    rw [iteratedFDeriv_one_apply]
    rfl
  exact (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp h

theorem cb_norm_fderiv_le {n : ℕ} (f : EuclideanSpace ℝ (Fin n) → ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : ‖fderiv ℝ f x‖ ≤ ‖gradE f x‖ := by
  rw [← cb_norm_jet_one, ← norm_iteratedFDeriv_one (𝕜 := ℝ) (x := x) f]
  exact norm_iteratedFDeriv_le_norm_jet 1 f x

theorem cb_jet_bounds {n : ℕ} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ₀ (ball (0 : EuclideanSpace ℝ (Fin n)) 6))
    (h0 : φ₀ 0 = 0)
    (hgr : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 2 φ₀ x‖ ≤ c₀ ∧ ‖gradE φ₀ x‖ ≤ c₀ * ‖x‖)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ₀ x‖ ≤ K₃) :
    ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, ∀ m : ℕ, m ≤ 3 → ‖jet m φ₀ x‖ < 2 := by
  intro x hx m hm
  have hx6 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := ball_subset_ball (by norm_num) hx
  have hx3 : ‖x‖ < 5 := by simpa using hx
  have hx3' : x ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_closedBall hx
  interval_cases m
  · rw [norm_jet_zero, Real.norm_eq_abs]
    have hb : ∀ y ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) 5, ‖fderiv ℝ φ₀ y‖ ≤ 5 / 100 := by
      intro y hy
      have hy6 : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := closedBall_subset_ball (by norm_num) hy
      have hy3 : ‖y‖ ≤ 5 := by simpa using hy
      calc ‖fderiv ℝ φ₀ y‖ ≤ ‖gradE φ₀ y‖ := cb_norm_fderiv_le _ _
        _ ≤ c₀ * ‖y‖ := (hgr y hy6).2
        _ ≤ 1 / 100 * 5 := by nlinarith [norm_nonneg y]
        _ = 5 / 100 := by norm_num
    have hd : ∀ y ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) 5, DifferentiableAt ℝ φ₀ y := fun y hy =>
      ((hφ.contDiffAt (isOpen_ball.mem_nhds
        (closedBall_subset_ball (by norm_num) hy))).differentiableAt (by simp))
    have := (convex_closedBall (0 : EuclideanSpace ℝ (Fin n)) 5).norm_image_sub_le_of_norm_fderiv_le hd hb
      (mem_closedBall_self (by norm_num)) hx3'
    rw [h0, sub_zero, sub_zero, Real.norm_eq_abs] at this
    nlinarith [norm_nonneg x]
  · rw [cb_norm_jet_one]
    have := (hgr x hx6).2
    nlinarith [norm_nonneg x]
  · have := (hgr x hx6).1
    linarith
  · have := hK x hx6
    linarith

theorem graphChart_adm {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
    {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) (N : ℕ) :
    ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, jp N φ x ∈ admSet n N 2 := by
  intro x hx m hm
  have := cb_jet_bounds hc₀ hc₀1 hK₃ hch.smooth hch.zero
    (fun x hx => ⟨(hch.graph x hx).2.2.1, (hch.graph x hx).2.2.2⟩) hK x hx m hm
  exact this

/-! ### The third derivative -/

theorem cb_exists_global_eventually {n : ℕ} {B : Set (EuclideanSpace ℝ (Fin n))} (hB : IsOpen B)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ B) {x : EuclideanSpace ℝ (Fin n)}
    (hx : x ∈ B) : ∃ ψ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ ψ ∧ φ =ᶠ[𝓝 x] ψ := by
  obtain ⟨R, hR, hRB⟩ : ∃ R : ℝ, 0 < R ∧
      closedBall (0 : EuclideanSpace ℝ (Fin n)) R ⊆ (fun v => x + v) ⁻¹' B := by
    have hO : IsOpen ((fun v : EuclideanSpace ℝ (Fin n) => x + v) ⁻¹' B) :=
      hB.preimage (continuous_const.add continuous_id)
    obtain ⟨ε, hε, hεB⟩ := Metric.isOpen_iff.1 hO 0 (by simpa using hx)
    exact ⟨ε / 2, by positivity, (closedBall_subset_ball (by linarith)).trans hεB⟩
  have hO : IsOpen ((fun v : EuclideanSpace ℝ (Fin n) => x + v) ⁻¹' B) :=
    hB.preimage (continuous_const.add continuous_id)
  have hφ' : ContDiffOn ℝ ∞ (fun v => φ (x + v)) ((fun v => x + v) ⁻¹' B) :=
    hφ.comp (contDiff_const.add contDiff_id).contDiffOn fun v hv => hv
  obtain ⟨R', hR', g, hg, -, hgφ⟩ := exists_global_of_contDiffOn hO hR hRB hφ'
  refine ⟨fun y => g (y - x), hg.comp (contDiff_id.sub contDiff_const), ?_⟩
  have hball : ball x R' ∈ 𝓝 x := ball_mem_nhds x (by linarith)
  filter_upwards [hball] with y hy
  have : y - x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' := by
    rw [mem_ball, dist_eq_norm] at hy ⊢
    simpa using hy
  rw [hgφ _ this]; simp

theorem cb_aNormSq_graph {n : ℕ} {M : Type*} {V : Type*} [NormedAddCommGroup V]
    [InnerProductSpace ℝ V] {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    {F : M → V} {σ : EuclideanSpace ℝ (Fin n) → M} {B : Set (EuclideanSpace ℝ (Fin n))}
    (hB : IsOpen B) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ B) {y₀ : V}
    (hgraph : ∀ x ∈ B, F (σ x) = graphX J U y₀ φ x) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ B) :
    aNormSq J (F ∘ σ) x = gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x) := by
  obtain ⟨ψ, hψc, hev⟩ := cb_exists_global_eventually hB hφ hx
  have hgrad : graphX J U y₀ φ =ᶠ[𝓝 x] graphX J U y₀ ψ := by
    have hev' : ∀ᶠ y in 𝓝 x, φ =ᶠ[𝓝 y] ψ := hev.eventually_nhds
    filter_upwards [hev'] with y hy
    unfold graphX gradE
    rw [hy.fderiv_eq]
  have hFσ : (F ∘ σ) =ᶠ[𝓝 x] graphX J U y₀ ψ := by
    filter_upwards [hB.mem_nhds hx, hgrad] with y hy hy2
    rw [Function.comp_apply, hgraph y hy, hy2]
  rw [aNormSq_congr hFσ, aNormSq_graphX hJi hJJ hU hψc, ← jet_congr_of_eventuallyEq hev 2,
    ← jet_congr_of_eventuallyEq hev 3]

theorem cb_aNormSq_param_le {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {F : M → V} (hF : IsLagImm ι J F)
    {c₁ : ℝ}
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ B) : aNormSq J (F ∘ σ) x ≤ c₁ ^ 2 := by
  set z := σ x with hz
  have hW := im_isOpen_W ι hσ z
  have hxW : x ∈ B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source := im_mem_W (σ := σ) hx
  have hev : (F ∘ σ) =ᶠ[𝓝 x] (F ∘ chartInv ι z) ∘ (chartFwd ι z ∘ σ) := by
    filter_upwards [hW.mem_nhds hxW] with w hw
    simp [Function.comp_def, im_chartInv_chartFwd ι hw.2]
  rw [aNormSq_congr hev,
    aNormSq_comp (im_isOpen_chartDomE ι z) (hF.smooth z)
      (fun y hy => hF.inj z y hy) (fun y hy => hF.lag z y hy) hW (hσ.trans_smooth z)
      (fun w hw => im_chartFwd_mem ι hw.2) (hσ.trans_det z) hxW]
  exact hA z _ (im_chartFwd_mem ι hxW.2)

theorem graphChart_jet3_le {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} (hF : IsLagImm ι J F) {c₀ c₁ : ℝ} (hc₀ : 0 ≤ c₀) (hc₁ : 0 ≤ c₁)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ) :
    ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ := by
  intro x hx
  have h1 := cb_aNormSq_param_le ι hF hA hch.param hx
  have h2 := cb_aNormSq_graph hJi hJJ hch.frame isOpen_ball hch.smooth
    (fun x hx => (hch.graph x hx).1) hx
  rw [h2] at h1
  have h3 := norm_jet3_sq_le φ x
  have hj2 := (hch.graph x hx).2.2.1
  have hj0 := norm_nonneg (jet 2 φ x)
  have h4 : ‖jet 3 φ x‖ ^ 2 ≤ ((1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁) ^ 2 := by
    have hrp : ((1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) ^ 2 = (1 + c₀ ^ 2) ^ 3 := by
      rw [← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]
      norm_num
    rw [mul_pow, hrp]
    calc ‖jet 3 φ x‖ ^ 2 ≤ (1 + ‖jet 2 φ x‖ ^ 2) ^ 3 *
          gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x) := h3
      _ ≤ (1 + ‖jet 2 φ x‖ ^ 2) ^ 3 * c₁ ^ 2 :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ ≤ (1 + c₀ ^ 2) ^ 3 * c₁ ^ 2 := by gcongr
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).1 h4

/-! ### The energy bound -/

theorem cb_jp_congr {n : ℕ} (N : ℕ) {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ}
    {x : EuclideanSpace ℝ (Fin n)} (h : φ =ᶠ[𝓝 x] ψ) : jp N φ x = jp N ψ x := by
  unfold jp
  congr 1
  funext m
  exact jet_congr_of_eventuallyEq h _

theorem exists_const_energy3_le_energyM {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M)
    (m : ℕ) {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ {F : M → V}, IsLagImm ι J F → ∀ {p U σ φ},
      IsGraphChart ι J F c₀ p U σ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) →
      energy3 m 3 φ ≤ C * (1 + ∑ l ∈ Finset.range (m + 1), energyM ι A J l F) := by
  obtain ⟨C₀, hC₀, hC⟩ := exists_const_energy3_le_eSym (n := n) m (N := m + 3) le_rfl (ρ := 3)
    (R := 4) (by norm_num) (by norm_num) (K := 2) (by norm_num)
  refine ⟨C₀, hC₀, fun {F} hF {p U σ φ} hch hK => ?_⟩
  obtain ⟨R', hR', ψ, hψ, -, hψφ⟩ := exists_global_of_contDiffOn (U := ball (0 : EuclideanSpace ℝ (Fin n)) 6)
    isOpen_ball (R := 5) (by norm_num) (closedBall_subset_ball (by norm_num)) hch.smooth
  have hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', φ =ᶠ[𝓝 x] ψ := fun x hx => by
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy using (hψφ y hy).symm
  have h5 : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' :=
    fun x hx => ball_subset_ball hR'.le hx
  have hadm : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (4 + 1), jp (m + 3) ψ x ∈ admSet n (m + 3) 2 := by
    intro x hx
    have hx5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
      ball_subset_ball (by norm_num : (4 : ℝ) + 1 ≤ 5) hx
    rw [← cb_jp_congr _ (hev x (h5 x hx5))]
    exact graphChart_adm ι hJi hJJ hc₀ hc₀1 hK₃ hch hK (m + 3) x hx5
  have h1 := hC ψ hψ hadm
  have hE : energy3 m 3 φ = energy3 m 3 ψ := by
    unfold energy3
    refine Finset.sum_congr rfl fun i _ => ?_
    refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
    have hx5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) hx
    rw [jet_congr_of_eventuallyEq (hev x (h5 x hx5))]
  have h2 : ∑ l ∈ Finset.range (m + 1),
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 4, eSym l (jp (m + 3) ψ x) ≤
      ∑ l ∈ Finset.range (m + 1), energyM ι A J l F := by
    refine Finset.sum_le_sum fun l hl => ?_
    have hl' : l ≤ m := Nat.lt_succ_iff.1 (Finset.mem_range.1 hl)
    refine le_of_eq_of_le ?_ (integral_eSym_le_energyM ι A hJi hJJ hch.frame hF hch.param hch.smooth
      (fun x hx => (hch.graph x hx).1) (R := 4) (by norm_num) l (m + 3) (by omega))
    refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
    have hx5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) hx
    rw [cb_jp_congr _ (hev x (h5 x hx5))]
  rw [hE]
  exact h1.trans (mul_le_mul_of_nonneg_left (by linarith) hC₀)

/-- info: 'CWFlow.Extension.graphChart_jet3_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms graphChart_jet3_le

/-- info: 'CWFlow.Extension.graphChart_adm' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms graphChart_adm

/-- info: 'CWFlow.Extension.exists_const_energy3_le_energyM' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_energy3_le_energyM

end CWFlow.Extension
