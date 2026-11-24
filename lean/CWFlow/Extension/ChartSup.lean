import CWFlow.Extension.ChartBounds
import CWFlow.Interp.SupBall

/-!
# Pointwise bounds for graph functions from bounded energies

If the global energies `E_l(F)`, `l ≤ m + n`, of a Lagrangian immersion are bounded by `B l`, then in every
graph chart (with the standard slope and third-derivative bounds) `‖D^m φ‖ ≤ C` on the ball of radius `2`,
with `C` depending only on `m, n, B, c₀, K₃`.  The proof combines the local sup bound `exists_const_sup_jet`
with the chart `L²` bound `exists_const_energy3_le_energyM` (Chen–Warren, Lemma 4.13).
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem cs_integral_sq_low {n : ℕ} {ψ : EuclideanSpace ℝ (Fin n) → ℝ} {j : ℕ}
    (h : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet j ψ y‖ < 2) :
    ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet j ψ y‖ ^ 2 ≤
      4 * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) 3)).toReal := by
  have := norm_setIntegral_le_of_norm_le_const_ae (μ := volume)
    (s := ball (0 : EuclideanSpace ℝ (Fin n)) 3) (f := fun y => ‖jet j ψ y‖ ^ 2) (C := 4)
    measure_ball_lt_top
    (ae_restrict_of_forall_mem measurableSet_ball fun y hy => by
      have h1 := h y hy
      have h2 := norm_nonneg (jet j ψ y)
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      nlinarith)
  refine (le_abs_self _).trans ?_
  rw [← Real.norm_eq_abs]
  simpa [measureReal_def] using this

theorem exists_const_chart_sup {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1) (m : ℕ) (B : ℕ → ℝ) :
    ∃ C : ℝ, ∀ {F : M → V}, IsLagImm ι J F →
      (∀ l ≤ m + n, energyM ι A J l F ≤ B l) →
      ∀ {p U σ φ}, IsGraphChart ι J F c₀ p U σ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ x‖ ≤ K₃) →
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet m φ x‖ ≤ C := by
  obtain ⟨C_S, hCS0, hCS⟩ := exists_const_sup_jet (n := n) m (ρ := 2) (R := 3) (by norm_num) (by norm_num)
  choose C4 hC4nn hC4 using fun q : ℕ => exists_const_energy3_le_energyM ι hJi hJJ A q hc₀ hc₀1 hK₃
  set S : ℝ := ∑ l ∈ Finset.range (m + n + 1), max (B l) 0 with hS
  set vol3 : ℝ := (volume (ball (0 : EuclideanSpace ℝ (Fin n)) 3)).toReal with hvol
  set D : ℝ := 4 * vol3 + ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) with hD
  refine ⟨C_S * (1 + (n + 1) * D), fun {F} hF hE {p U σ φ} hch hK x hx => ?_⟩
  have hS0 : 0 ≤ S := Finset.sum_nonneg fun l _ => le_max_right _ _
  have hvol0 : 0 ≤ vol3 := ENNReal.toReal_nonneg
  have hDsum : 0 ≤ ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) :=
    Finset.sum_nonneg fun q _ => mul_nonneg (hC4nn q) (by linarith)
  have hD0 : 0 ≤ D := by rw [hD]; positivity
  obtain ⟨R', hR', ψ, hψ, -, hψφ⟩ := exists_global_of_contDiffOn (U := ball (0 : EuclideanSpace ℝ (Fin n)) 6)
    isOpen_ball (R := 5) (by norm_num) (closedBall_subset_ball (by norm_num)) hch.smooth
  have hev : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', φ =ᶠ[𝓝 y] ψ := fun y hy => by
    filter_upwards [isOpen_ball.mem_nhds hy] with z hz using (hψφ z hz).symm
  have h3 : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' :=
    fun y hy => ball_subset_ball (by linarith) hy
  have hint : ∀ j, ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet j ψ y‖ ^ 2 =
      ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet j φ y‖ ^ 2 := by
    intro j
    refine setIntegral_congr_fun measurableSet_ball fun y hy => ?_
    rw [jet_congr_of_eventuallyEq (hev y (h3 y hy)) j]
  have hbd : ∀ j ≤ m + n, ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet j ψ y‖ ^ 2 ≤ D := by
    intro j hj
    rw [hint j]
    by_cases hj3 : j < 3
    · have := cs_integral_sq_low (n := n) (ψ := φ) (j := j) fun y hy => by
        have hy5 : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := ball_subset_ball (by norm_num) hy
        exact cb_jet_bounds hc₀ hc₀1 hK₃ hch.smooth hch.zero
          (fun x hx => ⟨(hch.graph x hx).2.2.1, (hch.graph x hx).2.2.2⟩) hK y hy5 j (by omega)
      rw [hD]
      linarith
    · obtain ⟨q, rfl⟩ : ∃ q, j = 3 + q := ⟨j - 3, by omega⟩
      have hq : q ≤ m + n := by omega
      have h1 : ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (3 + q) φ y‖ ^ 2 ≤ energy3 q 3 φ :=
        Finset.single_le_sum (f := fun i => ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3,
          ‖jet (3 + i) φ y‖ ^ 2) (fun i _ => setIntegral_nonneg measurableSet_ball fun y _ => by positivity)
          (Finset.mem_range.2 (Nat.lt_succ_self q))
      have h2 := hC4 q hF hch hK
      have h3' : ∑ l ∈ Finset.range (q + 1), energyM ι A J l F ≤ S := by
        calc ∑ l ∈ Finset.range (q + 1), energyM ι A J l F
            ≤ ∑ l ∈ Finset.range (q + 1), max (B l) 0 :=
              Finset.sum_le_sum fun l hl => (hE l (by
                have := Finset.mem_range.1 hl; omega)).trans (le_max_left _ _)
          _ ≤ S := Finset.sum_le_sum_of_subset_of_nonneg
              (Finset.range_subset_range.2 (by omega)) fun l _ _ => le_max_right _ _
      have h4 : energy3 q 3 φ ≤ C4 q * (1 + S) :=
        h2.trans (mul_le_mul_of_nonneg_left (by linarith) (hC4nn q))
      have h5 : C4 q * (1 + S) ≤ ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) :=
        Finset.single_le_sum (f := fun q => C4 q * (1 + S))
          (fun q _ => mul_nonneg (hC4nn q) (by linarith)) (Finset.mem_range.2 (by omega))
      rw [hD]
      linarith
  have hx' := hCS ψ hψ x hx
  have hxR : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' := ball_subset_ball (by linarith) hx
  have hsum : ∑ i ∈ Finset.range (n + 1),
      ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (m + i) ψ y‖ ^ 2 ≤ (n + 1) * D := by
    have := Finset.sum_le_card_nsmul (Finset.range (n + 1))
      (fun i => ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (m + i) ψ y‖ ^ 2) D
      fun i hi => hbd (m + i) (by have := Finset.mem_range.1 hi; omega)
    simpa using this
  calc ‖jet m φ x‖ = ‖jet m ψ x‖ := by rw [jet_congr_of_eventuallyEq (hev x hxR) m]
    _ ≤ _ := hx'
    _ ≤ C_S * (1 + (n + 1) * D) := mul_le_mul_of_nonneg_left (by linarith) hCS0

/-- info: 'CWFlow.Extension.exists_const_chart_sup' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_chart_sup

end CWFlow.Extension
