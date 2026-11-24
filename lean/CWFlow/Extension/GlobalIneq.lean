import CWFlow.Extension.GlobalSetup
import CWFlow.Extension.ChartDeriv
import CWFlow.Extension.ChartAbsorb
import CWFlow.Interp.EnergyGronwall

/-!
# The global differential inequality for the energies

Along a flow of Lagrangian immersions with `|A| ≤ c₁` and a volume bound, at every time `t₀`,
`E_k' ≤ C (1 + Σ_{l ≤ k} E_l)` (Chen–Warren, Prop. 4.2).  The proof assembles, over the charts of
`exists_global_setup`, the chart inequalities (G1) (`exists_const_chart_energy_deriv`), (G2), (G3)
(`exists_const_eSym_top_le`, `exists_const_eSym_next_le`) with Lemma 4.13
(`exists_const_energy3_le_energyM`), and concludes with `absorb_le`.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

/-! ### Real-number lemmas -/

theorem gi_sum_le {N N₀ : ℕ} (hN : N ≤ N₀) (f : Fin N → ℝ) {b : ℝ} (hb : 0 ≤ b)
    (h : ∀ j, f j ≤ b) : ∑ j, f j ≤ N₀ * b := by
  calc ∑ j, f j ≤ ∑ _j : Fin N, b := Finset.sum_le_sum fun j _ => h j
    _ = N * b := by simp
    _ ≤ N₀ * b := mul_le_mul_of_nonneg_right (by exact_mod_cast hN) hb

theorem gi_sum_lin {N N₀ : ℕ} (hN : N ≤ N₀) (D 𝔇 : Fin N → ℝ) (l : ℝ) {b : ℝ} (hb : 0 ≤ b)
    (h : ∀ j, D j ≤ l * 𝔇 j + b) : ∑ j, D j ≤ l * ∑ j, 𝔇 j + N₀ * b := by
  have := gi_sum_le hN (fun j => D j - l * 𝔇 j) hb (fun j => by linarith [h j])
  simp only [Finset.sum_sub_distrib, ← Finset.mul_sum] at this
  linarith

theorem gi_alg1 {lam 𝔇 Ds ε ε' cD Cc N₀ C2 S Y X : ℝ} (hε : 0 < ε)
    (hε' : ε' * (2 * N₀ * C2 + 1) = ε) (hε'0 : 0 ≤ ε') (hS : 0 ≤ S) (hX : 0 ≤ X) (hY : 0 ≤ Y)
    (h : Ds ≤ -lam * 𝔇 + N₀ * (2 * ε' * C2 * (1 + S + Y + X) + cD * (Cc * (1 + S + Y)))) :
    Ds ≤ -lam * 𝔇 + ε * (X + Y) + (ε + cD * N₀ * Cc) * (1 + S + Y) := by
  have hk : 2 * ε' * N₀ * C2 ≤ ε := by nlinarith
  have hW : 0 ≤ 1 + S + Y + X := by linarith
  nlinarith [mul_le_mul_of_nonneg_right hk hW, mul_nonneg hε.le hY]

theorem gi_alg2 {A₂ N₀ Cc 𝔇 Xs S Y : ℝ} (hA₂ : 0 ≤ A₂) (hN₀ : 0 ≤ N₀) (hCc : 0 ≤ Cc)
    (h𝔇 : 0 ≤ 𝔇) (hS : 0 ≤ S) (hY : 0 ≤ Y)
    (h : Xs ≤ A₂ * 𝔇 + N₀ * (A₂ * (Cc * (1 + S + Y)))) :
    Xs ≤ (A₂ * N₀ * Cc + A₂ + 1) * 𝔇 + (A₂ * N₀ * Cc + A₂ + 1) * (1 + Y + S) := by
  have h0 : 0 ≤ A₂ * N₀ * Cc := by positivity
  nlinarith [mul_nonneg (by linarith : 0 ≤ A₂ * N₀ * Cc + 1) h𝔇,
    mul_nonneg (by linarith : 0 ≤ A₂ + 1) (by linarith : 0 ≤ 1 + Y + S)]

theorem gi_alg3 {δ ε c S X Y : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ ε / 2) (hδ2 : δ ≤ 1 / 2) (hc : 0 ≤ c)
    (hS : 0 ≤ S) (hX : 0 ≤ X) (hY : 0 ≤ Y) (h : Y ≤ δ * (1 + S + Y + X) + c * (1 + S)) :
    Y ≤ ε * X + 2 * (1 + c) * (1 + S) := by
  have hS1 : 0 ≤ 1 + S := by linarith
  nlinarith [mul_le_mul_of_nonneg_right hδ hX, mul_le_mul_of_nonneg_right hδ2 hY,
    mul_le_mul_of_nonneg_right hδ2 hS1, mul_nonneg hc hS1]

theorem gi_per1 {Dj dj B2 e1 e2 lam ε' cD C1 C2 S Y X : ℝ}
    (h1 : Dj ≤ -lam * dj + ε' * B2 + ε' * e2 + cD * (1 + e1)) (hB : B2 ≤ e2)
    (he1 : e1 ≤ C1 * (1 + (S + Y))) (he2 : e2 ≤ C2 * (1 + (S + Y + X)))
    (hε' : 0 ≤ ε') (hcD : 0 ≤ cD) (hC1 : 0 ≤ C1) (hS : 0 ≤ S) (hY : 0 ≤ Y) :
    Dj ≤ -lam * dj + (2 * ε' * C2 * (1 + S + Y + X) + cD * ((1 + C1) * (1 + S + Y))) := by
  have h5 := mul_le_mul_of_nonneg_left (hB.trans he2) hε'
  have h7 := mul_le_mul_of_nonneg_left he2 hε'
  have h3 : C1 * (1 + (S + Y)) ≤ (1 + C1) * (1 + S + Y) := by nlinarith
  have h8 := mul_le_mul_of_nonneg_left (by linarith : 1 + e1 ≤ (1 + C1) * (1 + S + Y)) hcD
  nlinarith

theorem gi_per2 {Xj dj e1 A₂ C1 S Y : ℝ} (h1 : Xj ≤ A₂ * dj + A₂ * (1 + e1))
    (he1 : e1 ≤ C1 * (1 + (S + Y))) (hA₂ : 0 ≤ A₂) (hC1 : 0 ≤ C1) (hS : 0 ≤ S) (hY : 0 ≤ Y) :
    Xj ≤ A₂ * dj + A₂ * ((1 + C1) * (1 + S + Y)) := by
  have h3 : 1 + e1 ≤ (1 + C1) * (1 + S + Y) := by nlinarith
  have := mul_le_mul_of_nonneg_left h3 hA₂
  linarith

theorem gi_per3 {Yj B3 e0 e2 e c C0 C2 S Y X : ℝ}
    (h1 : Yj ≤ e * B3 + c * (1 + e0)) (hB : B3 ≤ e2) (he2 : e2 ≤ C2 * (1 + (S + Y + X)))
    (he0 : e0 ≤ C0 * (1 + S)) (he : 0 ≤ e) (hc : 0 ≤ c) (hC0 : 0 ≤ C0) (hS : 0 ≤ S) :
    Yj ≤ e * (C2 * (1 + S + Y + X)) + c * ((1 + C0) * (1 + S)) := by
  have h3 := mul_le_mul_of_nonneg_left (hB.trans he2) he
  have h6 : 1 + e0 ≤ (1 + C0) * (1 + S) := by nlinarith
  have h7 := mul_le_mul_of_nonneg_left h6 hc
  nlinarith

theorem gi_alg3' {ε e N₀ C2 c Cc S X Y : ℝ} (hε : 0 < ε) (he : 0 ≤ e)
    (heq : e * (N₀ * C2 + 1) = min (ε / 2) (1 / 2)) (hN₀ : 0 ≤ N₀) (hC2 : 0 ≤ C2) (hc : 0 ≤ c)
    (hCc : 0 ≤ Cc) (hS : 0 ≤ S) (hX : 0 ≤ X) (hY : 0 ≤ Y)
    (hs : Y ≤ N₀ * (e * (C2 * (1 + S + Y + X)) + c * (Cc * (1 + S)))) :
    Y ≤ ε * X + 2 * (1 + c * N₀ * Cc) * (1 + S) := by
  have hδ : e * N₀ * C2 ≤ min (ε / 2) (1 / 2) := by nlinarith
  exact gi_alg3 (δ := e * N₀ * C2) (c := c * N₀ * Cc) (by positivity)
    (hδ.trans (min_le_left _ _)) (hδ.trans (min_le_right _ _)) (by positivity) hS hX hY
    (by nlinarith [hs])

/-! ### Chart lemmas -/

section Chart

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- A global smooth extension of a graph chart function, admissible on `ball 0 5`. -/
theorem gi_ext (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) (N : ℕ) :
    ∃ ψ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ ψ ∧
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ) ∧
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, jp N ψ x ∈ admSet n N 2 := by
  obtain ⟨R', hR', ψ, hψ, -, hψφ⟩ := exists_global_of_contDiffOn
    (U := ball (0 : EuclideanSpace ℝ (Fin n)) 6)
    isOpen_ball (R := 5) (by norm_num) (closedBall_subset_ball (by norm_num)) hch.smooth
  have hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ := fun x hx => by
    filter_upwards [isOpen_ball.mem_nhds (ball_subset_ball hR'.le hx)] with y hy
    exact (hψφ y hy).symm
  refine ⟨ψ, hψ, hev, fun x hx => ?_⟩
  rw [← cb_jp_congr _ (hev x hx)]
  exact graphChart_adm ι hJi hJJ hc₀ hc₀1 hK₃ hch hK N x hx

theorem gi_cut_congr {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ)
    {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ) (m : ℕ) :
    ∫ x, ρ x ^ 2 * ‖jet m ψ x‖ ^ 2 = ∫ x, ρ x ^ 2 * ‖jet m φ x‖ ^ 2 := by
  refine integral_congr_ae (Filter.Eventually.of_forall fun x => ?_)
  by_cases h2 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2
  · have h5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) h2
    show ρ x ^ 2 * ‖jet m ψ x‖ ^ 2 = ρ x ^ 2 * ‖jet m φ x‖ ^ 2
    rw [jet_congr_of_eventuallyEq (hev x h5)]
  · show ρ x ^ 2 * ‖jet m ψ x‖ ^ 2 = ρ x ^ 2 * ‖jet m φ x‖ ^ 2
    rw [ca_rho_zero hρ h2]; ring

theorem gi_energy3_congr {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ} {R : ℝ} (hR : R ≤ 5)
    (hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ) (m : ℕ) :
    energy3 m R ψ = energy3 m R φ := by
  unfold energy3
  refine Finset.sum_congr rfl fun i _ => ?_
  refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
  have h5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball hR hx
  show ‖jet (3 + i) ψ x‖ ^ 2 = ‖jet (3 + i) φ x‖ ^ 2
  rw [jet_congr_of_eventuallyEq (hev x h5)]

theorem gi_setIntegral_congr {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ} {R : ℝ} (hR : R ≤ 5)
    (hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ) (m : ℕ) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet m ψ x‖ ^ 2 =
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet m φ x‖ ^ 2 := by
  refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
  have h5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball hR hx
  show ‖jet m ψ x‖ ^ 2 = ‖jet m φ x‖ ^ 2
  rw [jet_congr_of_eventuallyEq (hev x h5)]

/-- The chart integral of the energy density is the integral of the symbol. -/
theorem gi_chart_int (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ : ℝ} {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ) {Λ : ℝ}
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ) {ψ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, φ =ᶠ[𝓝 x] ψ) {N l : ℕ} (hN : l + 3 ≤ N) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, ρ x ^ 2 * eDens J (F ∘ σ) l x =
      ∫ x, ρ x ^ 2 * eSym l (jp N ψ x) := by
  have hz : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → ρ x = 0 := fun x hx =>
    ca_rho_zero hρ hx
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero (s := ball (0 : EuclideanSpace ℝ (Fin n)) 6)
    (f := fun x => ρ x ^ 2 * eSym l (jp N ψ x))
    (fun x hx => by
      show ρ x ^ 2 * eSym l (jp N ψ x) = 0
      rw [hz x (fun h => hx (ball_subset_ball (by norm_num) h))]; ring)]
  refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
  by_cases h2 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2
  · have h5 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) h2
    show ρ x ^ 2 * eDens J (F ∘ σ) l x = ρ x ^ 2 * eSym l (jp N ψ x)
    rw [eDens_eq_eSym_of_graph hJi hJJ hch.frame isOpen_ball hch.smooth
      (fun x hx => (hch.graph x hx).1) hN hx, cb_jp_congr _ (hev x h5)]
  · show ρ x ^ 2 * eDens J (F ∘ σ) l x = ρ x ^ 2 * eSym l (jp N ψ x)
    rw [hz x h2]; ring

/-- A top-order integral over a ball of radius `R ≤ 3` is bounded by `energy3`. -/
theorem gi_chart_ball (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) (k : ℕ) {R : ℝ}
    (hR : R ≤ 3) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet (k + 5) φ x‖ ^ 2 ≤ energy3 (k + 2) 3 φ := by
  obtain ⟨ψ, hψ, hev, -⟩ := gi_ext ι hJi hJJ hc₀ hc₀1 hK₃ hch hK 0
  rw [← gi_setIntegral_congr (by linarith) hev, ← gi_energy3_congr (by norm_num) hev]
  have hc : Continuous fun x => ‖jet (k + 5) ψ x‖ ^ 2 :=
    (contDiff_jet hψ (k + 5)).continuous.norm.pow 2
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet (k + 5) ψ x‖ ^ 2
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) ψ x‖ ^ 2 :=
        setIntegral_mono_set (integrableOn_ball_of_continuous hc 3)
          (Filter.Eventually.of_forall fun x => sq_nonneg _)
          (ball_subset_ball hR).eventuallyLE
    _ ≤ energy3 (k + 2) 3 ψ := by
        unfold energy3
        have := Finset.single_le_sum (f := fun i => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3,
          ‖jet (3 + i) ψ x‖ ^ 2) (s := Finset.range (k + 2 + 1))
          (fun i _ => setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _)
          (Finset.self_mem_range_succ (k + 2))
        have e : k + 5 = 3 + (k + 2) := by omega
        rw [e]
        exact this

/-- (G2) at the level of a chart of the flow's time slice. -/
theorem gi_chart_top (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) (k : ℕ) {Λ A₂ : ℝ}
    (hA₂ : ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ ψ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ ψ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp (k + 7) ψ x ∈ admSet n (k + 7) 2) →
      ∫ x, ρ x ^ 2 * eSym (k + 2) (jp (k + 7) ψ x) ≤
        A₂ * (∫ x, ρ x ^ 2 * ‖jet (k + 5) ψ x‖ ^ 2) + A₂ * (1 + energy3 (k + 1) 3 ψ))
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, ρ x ^ 2 * eDens J (F ∘ σ) (k + 2) x ≤
      A₂ * (∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2) + A₂ * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨ψ, hψ, hev, hadm⟩ := gi_ext ι hJi hJJ hc₀ hc₀1 hK₃ hch hK (k + 7)
  rw [gi_chart_int ι hJi hJJ hch hρ hev (N := k + 7) (by omega)]
  have h := hA₂ ρ hρ ψ hψ (fun x hx => hadm x (ball_subset_ball (by norm_num) hx))
  rwa [gi_cut_congr hρ hev, gi_energy3_congr (by norm_num) hev] at h

/-- (G3) at the level of a chart of the flow's time slice. -/
theorem gi_chart_next (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {F : M → V} {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {p U σ φ} (hch : IsGraphChart ι J F c₀ p U σ φ)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) (k : ℕ) {Λ C₃ ε : ℝ}
    (hA₃ : ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ ψ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ ψ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp (k + 7) ψ x ∈ admSet n (k + 7) 2) →
      ∫ x, ρ x ^ 2 * eSym (k + 1) (jp (k + 7) ψ x) ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) ψ x‖ ^ 2) +
          C₃ * (1 + energy3 k 3 ψ))
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, ρ x ^ 2 * eDens J (F ∘ σ) (k + 1) x ≤
      ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) φ x‖ ^ 2) +
        C₃ * (1 + energy3 k 3 φ) := by
  obtain ⟨ψ, hψ, hev, hadm⟩ := gi_ext ι hJi hJJ hc₀ hc₀1 hK₃ hch hK (k + 7)
  rw [gi_chart_int ι hJi hJJ hch hρ hev (N := k + 7) (by omega)]
  have h := hA₃ ρ hρ ψ hψ (fun x hx => hadm x (ball_subset_ball (by norm_num) hx))
  rwa [gi_setIntegral_congr (by norm_num) hev, gi_energy3_congr (by norm_num) hev] at h

end Chart

/-- The global differential inequality `E_k' ≤ C (1 + Σ_{l ≤ k} E_l)` (Chen–Warren Prop. 4.2). -/
theorem exists_const_energy_ineq {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n) (k : ℕ) (hk : 1 ≤ k)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hK₃ : (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ ≤ 1) (V₀ : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ {I : Set ℝ} {Φ : ℝ → M → V}, IsFlowM ι J I Φ →
      (∀ t ∈ I, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2) →
      ∀ {t₀ : ℝ}, t₀ ∈ I → volM ι A (Φ t₀) ≤ V₀ →
      ∃ D : ℝ, HasDerivAt (fun s => energyM ι A J k (Φ s)) D t₀ ∧
        D ≤ C * (1 + ∑ l ∈ Finset.range (k + 1), energyM ι A J l (Φ t₀)) := by
  classical
  obtain ⟨N₀, Λ, hsetup⟩ := exists_global_setup ι A hJi hJJ hdim hc₀ hc₀1 hc₁ hsmall V₀
  choose C4 hC4nn hC4 using fun m : ℕ =>
    exists_const_energy3_le_energyM ι hJi hJJ A m hc₀.le hc₀1 hK₃
  choose! cD hcD0 hcD using fun (ε : ℝ) (hε : 0 < ε) =>
    exists_const_chart_energy_deriv (M := M) ι hJi hJJ hdim k hk Λ hε hc₀.le hc₀1 hK₃
  obtain ⟨A₂, hA₂0, hA₂⟩ := exists_const_eSym_top_le (n := n) (N := k + 7) (k := k) le_rfl
    (K := 2) (by norm_num) Λ
  choose! c3 hc30 hc3 using fun (ε : ℝ) (hε : 0 < ε) =>
    exists_const_eSym_next_le (n := n) (N := k + 7) (k := k) le_rfl (K := 2) (by norm_num) Λ hε
  obtain ⟨lam, hlamdef⟩ : ∃ lam : ℝ, lam = 2 * ((1 + (2 : ℝ) ^ 2)⁻¹) ^ (k + 5) := ⟨_, rfl⟩
  have hlam : 0 < lam := by rw [hlamdef]; positivity
  have hN₀ : (0 : ℝ) ≤ N₀ := Nat.cast_nonneg _
  have hCc1 : 0 ≤ 1 + C4 (k + 1) := by linarith [hC4nn (k + 1)]
  have hCc0 : 0 ≤ 1 + C4 k := by linarith [hC4nn k]
  obtain ⟨Aa, hAa⟩ : ∃ Aa : ℝ, Aa = A₂ * N₀ * (1 + C4 (k + 1)) + A₂ + 1 := ⟨_, rfl⟩
  have hApos : 0 < Aa := by rw [hAa]; positivity
  -- the functions `c₁`, `c₃` of the absorption lemma
  obtain ⟨c1f, hc1f⟩ : ∃ c1f : ℝ → ℝ, ∀ ε, c1f ε =
      ε + cD (ε / (2 * N₀ * C4 (k + 2) + 1)) * N₀ * (1 + C4 (k + 1)) := ⟨_, fun _ => rfl⟩
  obtain ⟨ε2, hε2⟩ : ∃ ε2 : ℝ → ℝ, ∀ ε, ε2 ε = min (ε / 2) (1 / 2) / (N₀ * C4 (k + 2) + 1) :=
    ⟨_, fun _ => rfl⟩
  obtain ⟨c3f, hc3f⟩ : ∃ c3f : ℝ → ℝ, ∀ ε,
      c3f ε = 2 * (1 + c3 (ε2 ε) * N₀ * (1 + C4 k)) := ⟨_, fun _ => rfl⟩
  obtain ⟨Cabs, hCabs0, hCabs⟩ := absorb_le (lam := lam) (A := Aa) hlam hApos c1f c3f
  refine ⟨Cabs, hCabs0, ?_⟩
  intro I Φ hΦ hbd t₀ ht₀ hV
  obtain ⟨N, p, U, σ, φ, ρ, hN, hch, hjet, hcut, hrep⟩ := hsetup hΦ hbd ht₀ hV
  have himm := hΦ.imm t₀ ht₀
  obtain ⟨ρh, hρh⟩ : ∃ ρh : Fin N → EuclideanSpace ℝ (Fin n) → ℝ, ∀ j x,
      ρh j x = if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0 :=
    ⟨fun j x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0,
      fun _ _ => rfl⟩
  have hcut' : ∀ j, IsCutoff (ρh j) Λ := fun j => by
    have : ρh j = fun x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0 :=
      funext (hρh j)
    rw [this]; exact hcut j
  have hrep' : ∀ s ∈ I, ∀ l : ℕ, energyM ι A J l (Φ s) =
      ∑ j, ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ s ∘ σ j) l x := by
    intro s hs l
    rw [hrep s hs l]
    refine Finset.sum_congr rfl fun j _ => ?_
    refine setIntegral_congr_fun measurableSet_ball fun x _ => ?_
    simp only [hρh]
  have hE0 : ∀ (l : ℕ) (s : ℝ), 0 ≤ energyM ι A J l (Φ s) := fun l s => by
    unfold energyM
    exact intM_densOf_nonneg ι A (fun _ => zero_le_one) (fun X y => ec_eDens_nonneg J X l y)
  -- derivatives
  have hder : ∀ j : Fin N, ∃ Dj : ℝ, HasDerivAt (fun s => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ s ∘ σ j) k x) Dj t₀ ∧
      ∀ ε : ℝ, 0 < ε → Dj ≤ -lam * (∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) +
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) (φ j) x‖ ^ 2) +
        ε * energy3 (k + 2) 3 (φ j) + cD ε * (1 + energy3 (k + 1) 3 (φ j)) := by
    intro j
    obtain ⟨D₁, hD₁, -⟩ := hcD 1 one_pos hΦ ht₀ (hch j) (hjet j) (hcut' j)
    refine ⟨D₁, hD₁, fun ε hε => ?_⟩
    obtain ⟨D', hD', hb⟩ := hcD ε hε hΦ ht₀ (hch j) (hjet j) (hcut' j)
    rw [hD₁.unique hD', hlamdef]
    exact hb
  choose D hDd hDb using hder
  have hderiv : HasDerivAt (fun s => energyM ι A J k (Φ s)) (∑ j, D j) t₀ := by
    have h := HasDerivAt.fun_sum (u := Finset.univ) (fun j _ => hDd j)
    refine h.congr_of_eventuallyEq ?_
    filter_upwards [hΦ.isOpen.mem_nhds ht₀] with s hs
    exact hrep' s hs k
  refine ⟨∑ j, D j, hderiv, ?_⟩
  -- the quantities at time `t₀`
  obtain ⟨S, hS⟩ : ∃ S : ℝ, S = ∑ l ∈ Finset.range (k + 1), energyM ι A J l (Φ t₀) := ⟨_, rfl⟩
  obtain ⟨Y, hY⟩ : ∃ Y : ℝ, Y = energyM ι A J (k + 1) (Φ t₀) := ⟨_, rfl⟩
  obtain ⟨X, hX⟩ : ∃ X : ℝ, X = energyM ι A J (k + 2) (Φ t₀) := ⟨_, rfl⟩
  rw [← hS]
  have hS0 : 0 ≤ S := by rw [hS]; exact Finset.sum_nonneg fun l _ => hE0 l t₀
  have hY0 : 0 ≤ Y := by rw [hY]; exact hE0 _ _
  have hX0 : 0 ≤ X := by rw [hX]; exact hE0 _ _
  have hsum1 : ∑ l ∈ Finset.range (k + 1 + 1), energyM ι A J l (Φ t₀) = S + Y := by
    rw [Finset.sum_range_succ, hS, hY]
  have hsum2 : ∑ l ∈ Finset.range (k + 2 + 1), energyM ι A J l (Φ t₀) = S + Y + X := by
    rw [Finset.sum_range_succ, hsum1, hX]
  -- Lemma 4.13
  have he0 : ∀ j, energy3 k 3 (φ j) ≤ C4 k * (1 + S) := fun j => by
    have := hC4 k himm (hch j) (hjet j)
    rwa [← hS] at this
  have he1 : ∀ j, energy3 (k + 1) 3 (φ j) ≤ C4 (k + 1) * (1 + (S + Y)) := fun j => by
    have := hC4 (k + 1) himm (hch j) (hjet j)
    rwa [hsum1] at this
  have he2 : ∀ j, energy3 (k + 2) 3 (φ j) ≤ C4 (k + 2) * (1 + (S + Y + X)) := fun j => by
    have := hC4 (k + 2) himm (hch j) (hjet j)
    rwa [hsum2] at this
  have h𝔇j : ∀ j, 0 ≤ ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2 := fun j =>
    integral_nonneg fun x => mul_nonneg (sq_nonneg _) (sq_nonneg _)
  have h𝔇 : 0 ≤ ∑ j, ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2 :=
    Finset.sum_nonneg fun j _ => h𝔇j j
  -- (G1)
  have hG1 : ∀ ε > 0, 0 ≤ c1f ε ∧ ∑ j, D j ≤ -lam * (∑ j, ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) +
      ε * (X + Y) + c1f ε * (1 + S + Y) := by
    intro ε hε
    have hd : 0 < 2 * (N₀ : ℝ) * C4 (k + 2) + 1 := by
      have := hC4nn (k + 2)
      positivity
    obtain ⟨ε', hε'def⟩ : ∃ e : ℝ, e = ε / (2 * (N₀ : ℝ) * C4 (k + 2) + 1) := ⟨_, rfl⟩
    have hε' : 0 < ε' := by rw [hε'def]; exact div_pos hε hd
    have hε'eq : ε' * (2 * N₀ * C4 (k + 2) + 1) = ε := by
      rw [hε'def]; exact div_mul_cancel₀ _ hd.ne'
    rw [hc1f, ← hε'def]
    refine ⟨?_, ?_⟩
    · have := hcD0 _ hε'
      positivity
    · 
      have hper : ∀ j, D j ≤ -lam * (∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) +
          (2 * ε' * C4 (k + 2) * (1 + S + Y + X) +
            cD ε' * ((1 + C4 (k + 1)) * (1 + S + Y))) := by
        intro j
        exact gi_per1 (hDb j ε' hε') (gi_chart_ball ι hJi hJJ hc₀.le hc₀1 hK₃ (hch j) (hjet j) k
          (R := 2) (by norm_num)) (he1 j) (he2 j) hε'.le (hcD0 _ hε') (hC4nn _) hS0 hY0
      have hb0 : 0 ≤ 2 * ε' * C4 (k + 2) * (1 + S + Y + X) +
          cD ε' * ((1 + C4 (k + 1)) * (1 + S + Y)) := by
        have := hcD0 _ hε'
        have := hC4nn (k + 2)
        positivity
      have hs := gi_sum_lin hN D (fun j => ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) (-lam) hb0
        hper
      exact gi_alg1 hε hε'eq hε'.le hS0 hX0 hY0 hs
  -- (G2)
  have hG2 : X ≤ Aa * (∑ j, ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) + Aa * (1 + Y + S) := by
    have hper : ∀ j, ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ t₀ ∘ σ j) (k + 2) x ≤
        A₂ * (∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) +
          A₂ * ((1 + C4 (k + 1)) * (1 + S + Y)) := by
      intro j
      exact gi_per2 (gi_chart_top ι hJi hJJ hc₀.le hc₀1 hK₃ (hch j) (hjet j) k hA₂ (hcut' j))
        (he1 j) hA₂0 (hC4nn _) hS0 hY0
    have hb0 : 0 ≤ A₂ * ((1 + C4 (k + 1)) * (1 + S + Y)) :=
      mul_nonneg hA₂0 (mul_nonneg hCc1 (by linarith))
    have hs := gi_sum_lin hN (fun j => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ t₀ ∘ σ j) (k + 2) x)
      (fun j => ∫ x, ρh j x ^ 2 * ‖jet (k + 5) (φ j) x‖ ^ 2) A₂ hb0 hper
    rw [← hrep' t₀ ht₀ (k + 2), ← hX] at hs
    rw [hAa]
    exact gi_alg2 hA₂0 hN₀ hCc1 h𝔇 hS0 hY0 hs
  -- (G3)
  have hG3 : ∀ ε > 0, 0 ≤ c3f ε ∧ Y ≤ ε * X + c3f ε * (1 + S) := by
    intro ε hε
    have hd : 0 < (N₀ : ℝ) * C4 (k + 2) + 1 := by
      have := mul_nonneg hN₀ (hC4nn (k + 2))
      linarith
    have hm : 0 < min (ε / 2) (1 / 2) := lt_min (by positivity) (by norm_num)
    obtain ⟨e, hedef⟩ : ∃ e : ℝ, e = ε2 ε := ⟨_, rfl⟩
    have he : 0 < e := by rw [hedef, hε2]; exact div_pos hm hd
    have heq : e * (N₀ * C4 (k + 2) + 1) = min (ε / 2) (1 / 2) := by
      rw [hedef, hε2]; exact div_mul_cancel₀ _ hd.ne'
    have hc30' := hc30 e he
    refine ⟨by
      rw [hc3f, ← hedef]
      exact mul_nonneg zero_le_two (add_nonneg zero_le_one
        (mul_nonneg (mul_nonneg hc30' hN₀) hCc0)), ?_⟩
    have hper : ∀ j, ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ t₀ ∘ σ j) (k + 1) x ≤
        e * (C4 (k + 2) * (1 + S + Y + X)) + c3 e * ((1 + C4 k) * (1 + S)) := by
      intro j
      exact gi_per3 (gi_chart_next ι hJi hJJ hc₀.le hc₀1 hK₃ (hch j) (hjet j) k (hc3 e he)
        (hcut' j)) (gi_chart_ball ι hJi hJJ hc₀.le hc₀1 hK₃ (hch j) (hjet j) k (R := 3) le_rfl)
        (he2 j) (he0 j) he.le hc30' (hC4nn _) hS0
    have hb0 : 0 ≤ e * (C4 (k + 2) * (1 + S + Y + X)) + c3 e * ((1 + C4 k) * (1 + S)) :=
      add_nonneg (mul_nonneg he.le (mul_nonneg (hC4nn _) (by linarith)))
        (mul_nonneg hc30' (mul_nonneg hCc0 (by linarith)))
    have hs := gi_sum_le hN (fun j => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ρh j x ^ 2 * eDens J (Φ t₀ ∘ σ j) (k + 1) x) hb0 hper
    rw [← hrep' t₀ ht₀ (k + 1), ← hY] at hs
    rw [hc3f, ← hedef]
    exact gi_alg3' hε he.le heq hN₀ (hC4nn _) hc30' hCc0 hS0 hX0 hY0 hs
  have hfin := hCabs (∑ j, D j) X Y S _ h𝔇 hS0 hX0 hY0 hG1 hG2 hG3
  exact hfin

/-- info: 'CWFlow.Extension.exists_const_energy_ineq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_energy_ineq

end CWFlow.Extension
