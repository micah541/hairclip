import CWFlow.Extension.PieceFlow
import CWFlow.Extension.GraphPiece
import CWFlow.Extension.GraphTower
import CWFlow.Extension.LimitChartAux
import CWFlow.Interp.OdeTower
import CWFlow.Interp.ClosureExt

/-!
# The limit at `T` in a graph chart

Let `Φ` be a flow of Lagrangian immersions of `M` on `(t₀, T)` with a graph chart `(p, U, σ, φ₀)` of
`Φ t₀` that persists up to `T` (`hpers`, the conclusion of `exists_persist`), and with uniformly
bounded energies and second fundamental form.  Then `Y (t, y) = Φ t (σ y)` has all derivatives
bounded on `(t₁, T) × B₂` for some `t₁ < T`, hence extends smoothly to `[t₁, T] × B̄₂`
(`exists_limit_chart`).

Proof (blueprint, "limit at `T`"):

1. `exists_piece_flow` describes `Φ t ∘ σ` on `(t₀, T)` by one graph function `φ` over `B₅`, with
   `∂ₜφ = potSym (jp 4 φ_t)` and `χ = horC ∘ Φ ∘ σ` obeying `∂ₜχ = vSym (jp 5 φ_t) ∘ χ`;
2. `exists_const_piece_sup` bounds all jets of order `≥ 1` of `φ_t` on `B₃` uniformly in `t`;
3. (`lc_shift_bounds`) the potential does not see constants, so `∂ₜφ` is bounded; subtracting
   `φ(t₂, 0)` for `t₂` close to `T` makes `φ̃` small, and `graph_iteratedFDeriv_bounded` bounds all
   derivatives of `φ̃` and of `vSym ∘ jp 5 φ̃` on `(t₂, T) × B₃`;
4. (`lc_ode_comp`) `ode_iteratedFDeriv_bounded` bounds the derivatives of `χ`, and the
   composition `Y = graphX φ̃ ∘ (t, χ)` has bounded derivatives;
5. `exists_contDiffOn_closure_of_bounded` extends to the closure.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem exists_limit_chart {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁) (hc₁1 : c₁ ≤ 1 / 2) (B : ℕ → ℝ)
    {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {t₀ T : ℝ} (hT : t₀ < T) (hI : Ioo t₀ T ⊆ I)
    (hA : ∀ t ∈ Ioo t₀ T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2)
    (hE : ∀ t ∈ Ioo t₀ T, ∀ l, energyM ι A J l (Φ t) ≤ B l)
    {p : M} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {σ : EuclideanSpace ℝ (Fin n) → M}
    {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ} (hch : IsGraphChart ι J (Φ t₀) c₀ p U σ φ₀)
    (hpers : ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
      ‖fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      ‖horC U (Φ t₀ p) (Φ t (σ y)) - y‖ ≤ 1 / 8 ∧
      ‖verC J U (Φ t₀ p) (Φ t (σ y))‖ ≤ 1 / 2 ∧
      ‖slopeZ J U (fderiv ℝ (fun y' => Φ t (σ y')) y)‖ ≤ 2 * c₀) :
    ∃ t₁ ∈ Ioo t₀ T, ∃ Ybar : ℝ × EuclideanSpace ℝ (Fin n) → V,
      ContDiffOn ℝ ∞ Ybar (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) ∧
      EqOn Ybar (fun q => Φ q.1 (σ q.2)) (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) := by
  classical
  obtain ⟨ξ, φ, hφ, hmain⟩ := exists_piece_flow ι hJi hJJ hdim hc₀ hc₀1 hc₁ hΦ hT hI hA hch hpers
  have hB3 : ∀ x₀ ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      ball x₀ 2 ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 5 := by
    intro x₀ hx₀ x hx
    rw [mem_ball_zero_iff] at hx₀ ⊢
    rw [mem_ball, dist_eq_norm] at hx
    calc ‖x‖ = ‖(x - x₀) + x₀‖ := by simp
      _ ≤ ‖x - x₀‖ + ‖x₀‖ := norm_add_le _ _
      _ < 5 := by linarith
  have hsl : ∀ t ∈ Ioo t₀ T, ContDiffOn ℝ ∞ (fun y => φ (t, y))
      (ball (0 : EuclideanSpace ℝ (Fin n)) 5) := fun t ht => tg_slice hφ ht
  have hpar : ∀ t ∈ Ioo t₀ T,
      IsParam ι (ball (0 : EuclideanSpace ℝ (Fin n)) 5) (fun x => σ (ξ (t, x))) :=
    fun t ht => (hmain t ht).1
  -- Step 1: uniform bounds of the jets of order `≥ 1` on `B₃`
  have hCm : ∀ m : ℕ, ∃ C : ℝ, 1 ≤ m → ∀ t ∈ Ioo t₀ T,
      ∀ x₀ ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet m (fun y => φ (t, y)) x₀‖ ≤ C := by
    intro m
    by_cases hm : 1 ≤ m
    · obtain ⟨C, hC⟩ := exists_const_piece_sup ι A hJi hJJ hm B
      refine ⟨C, fun _ t ht x₀ hx₀ => ?_⟩
      refine hC (F := Φ t) (y₀ := Φ t₀ p) (hΦ.imm t (hI ht)) (fun l _ => hE t ht l) hch.frame
        (lc_isParam_mono ι (hpar t ht) isOpen_ball (hB3 x₀ hx₀)) ((hsl t ht).mono (hB3 x₀ hx₀))
        (fun x hx => ?_)
      obtain ⟨-, -, hgr, -, hgd, hj2, hj3⟩ := (hmain t ht).2.1 x (hB3 x₀ hx₀ hx)
      exact ⟨hgr, (by linarith), by linarith, by linarith⟩
    · exact ⟨0, fun h => absurd h hm⟩
  -- Step 2: shift and bound
  obtain ⟨t₂, ht₂, c₂, hlow, hGb, hVb⟩ := lc_shift_bounds hJi hJJ hch.frame hT hφ
    (fun t ht x hx => (hmain t ht).2.1 x hx |>.2.2.2.1)
    (fun t ht x hx => (hmain t ht).2.1 x hx |>.2.2.2.2.1)
    (fun t ht x hx => by
      have := (hmain t ht).2.1 x hx |>.2.2.2.2.2.1
      linarith)
    (fun t ht x hx => by
      have := (hmain t ht).2.1 x hx |>.2.2.2.2.2.2
      linarith)
    hCm
  have hsubO : ∀ t ∈ Ioo t₂ T, t ∈ Ioo t₀ T := fun t ht => ⟨by linarith [ht.1, ht₂.1], ht.2⟩
  set t₁ : ℝ := (t₂ + T) / 2 with ht₁def
  have h21 : t₂ < t₁ := by rw [ht₁def]; linarith [ht₂.2]
  have h1T : t₁ < T := by rw [ht₁def]; linarith [ht₂.2]
  have ht₁ : t₁ ∈ Ioo t₀ T := ⟨by linarith [ht₂.1], h1T⟩
  -- Step 3: the horizontal coordinate
  have hY : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2))
      (Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    (cge_contDiffOn_family ι hΦ.smooth hch.param hΦ.isOpen).mono (prod_mono hI subset_rfl)
  have hχ6 : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      horC U (Φ t₀ p) (Φ q.1 (σ q.2))) (Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    (horCL U).contDiff.comp_contDiffOn (hY.sub contDiffOn_const)
  have hχ5 : ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      horC U (Φ t₀ p) (Φ t (σ y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3 := by
    intro t ht y hy
    have h1 := (hpers t ht y (ball_subset_ball (by norm_num) hy)).2.1
    rw [mem_ball_zero_iff] at hy ⊢
    calc ‖horC U (Φ t₀ p) (Φ t (σ y))‖
        = ‖(horC U (Φ t₀ p) (Φ t (σ y)) - y) + y‖ := by simp
      _ ≤ ‖horC U (Φ t₀ p) (Φ t (σ y)) - y‖ + ‖y‖ := norm_add_le _ _
      _ < 3 := by linarith
  have hS2 : Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2 ⊆
      Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6 :=
    prod_mono (Ioo_subset_Ioo_left ht₂.1.le) (ball_subset_ball (by norm_num))
  have hS3 : Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3 ⊆
      Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
    prod_mono (Ioo_subset_Ioo_left ht₂.1.le) (ball_subset_ball (by norm_num))
  have hF : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => φ q - c₂)
      (Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    (hφ.sub contDiffOn_const).mono hS3
  have hsb := lc_ode_comp (J := J) (U := U) (y₀ := Φ t₀ p) h21 h1T
    (F := fun q : ℝ × EuclideanSpace ℝ (Fin n) => φ q - c₂) hF hlow hGb hVb
    (χ := fun q : ℝ × EuclideanSpace ℝ (Fin n) => horC U (Φ t₀ p) (Φ q.1 (σ q.2)))
    (hχ6.mono hS2)
    (fun q hq => ⟨hq.1, hχ5 q.1 (hsubO q.1 hq.1) q.2 hq.2⟩)
    (fun q hq c => by
      have hqO := hsubO q.1 hq.1
      have hy6 : q.2 ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) :=
        ball_subset_ball (by norm_num) hq.2
      have hx5 : horC U (Φ t₀ p) (Φ q.1 (σ q.2)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
        ball_subset_ball (by norm_num) (hχ5 q.1 hqO q.2 hq.2)
      have h := ((hmain q.1 hqO).2.2 q.2 hy6 hx5).2 c
      exact h.congr_deriv (lc_vSym_sub_const hJi hJJ hch.frame isOpen_ball (hsl q.1 hqO) hx5 c₂
        (le_refl 5) c).symm)
    (fun j => by
      have hf : ContDiffOn ℝ ∞ (fun y' => horC U (Φ t₀ p) (Φ t₁ (σ y')))
          (ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2)) :=
        tg_slice (g := fun q : ℝ × EuclideanSpace ℝ (Fin n) => horC U (Φ t₀ p) (Φ q.1 (σ q.2)))
          (hχ6.mono (prod_mono subset_rfl (ball_subset_ball (by norm_num)))) ht₁
      obtain ⟨C, hC⟩ := lc_exists_bound_jet (W := EuclideanSpace ℝ (Fin n)) isOpen_ball
        (R := 2) (by norm_num) (closedBall_subset_ball (by norm_num)) hf j
      exact ⟨C, fun y hy => hC y (ball_subset_closedBall hy)⟩)
  -- Step 4: the extension to the closure
  obtain ⟨Ybar, hYb, hYeq⟩ := exists_contDiffOn_closure_of_bounded
    (isOpen_Ioo.prod isOpen_ball) ((convex_Ioo t₁ T).prod (convex_ball _ _)) hsb.1 hsb.2
  have hcl : closure (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) =
      Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2 := by
    rw [closure_prod_eq, closure_Ioo (ne_of_lt h1T), closure_ball _ (by norm_num)]
  refine ⟨t₁, ht₁, Ybar, hcl ▸ hYb, fun q hq => ?_⟩
  have hqO : q.1 ∈ Ioo t₀ T := ⟨by linarith [ht₁.1, hq.1.1], hq.1.2⟩
  have hy6 : q.2 ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) :=
    ball_subset_ball (by norm_num) hq.2
  have hx5 : horC U (Φ t₀ p) (Φ q.1 (σ q.2)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
    ball_subset_ball (by norm_num) (hχ5 q.1 hqO q.2 hq.2)
  obtain ⟨hξ, -⟩ := (hmain q.1 hqO).2.2 q.2 hy6 hx5
  obtain ⟨-, -, hgraph, -⟩ := (hmain q.1 hqO).2.1 _ hx5
  calc Ybar q
      = graphX J U (Φ t₀ p) (fun y => φ (q.1, y) - c₂) (horC U (Φ t₀ p) (Φ q.1 (σ q.2))) := hYeq hq
    _ = graphX J U (Φ t₀ p) (fun y => φ (q.1, y)) (horC U (Φ t₀ p) (Φ q.1 (σ q.2))) := by
        rw [lc_graphX_sub_const]
    _ = Φ q.1 (σ (ξ (q.1, horC U (Φ t₀ p) (Φ q.1 (σ q.2))))) := hgraph.symm
    _ = Φ q.1 (σ q.2) := by rw [hξ]

/-- info: 'CWFlow.Extension.exists_limit_chart' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_limit_chart

end CWFlow.Extension
