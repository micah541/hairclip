import CWFlow.Extension.Partition
import CWFlow.Extension.ChartBounds
import CWFlow.Extension.VolM

/-!
# Global setup at a time `t₀`

At a time `t₀` we take an adapted cover by graph charts (number of charts bounded through the
volume) and a partition of unity `Σ ρ_j² = 1` subordinate to it, whose pieces are uniform cutoffs in
graph coordinates; the energies at every time `s` are then sums of chart integrals over these
fixed charts.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem gs_tsupport_sq_subset {M : Type*} [TopologicalSpace M] (f : M → ℝ) :
    tsupport (fun q => f q ^ 2) ⊆ tsupport f :=
  closure_mono fun q hq => by
    simp only [Function.mem_support] at hq ⊢
    exact fun h => hq (by simp [h])

open Classical in
theorem exists_global_setup {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀) (V₀ : ℝ) :
    ∃ (N₀ : ℕ) (Λ : ℝ), ∀ {I : Set ℝ} {Φ : ℝ → M → V}, IsFlowM ι J I Φ →
      (∀ t ∈ I, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2) →
      ∀ {t₀ : ℝ}, t₀ ∈ I → volM ι A (Φ t₀) ≤ V₀ →
      ∃ (N : ℕ) (p : Fin N → M) (U : Fin N → (EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V))
        (σ : Fin N → EuclideanSpace ℝ (Fin n) → M) (φ : Fin N → EuclideanSpace ℝ (Fin n) → ℝ)
        (ρ : Fin N → M → ℝ),
        N ≤ N₀ ∧
        (∀ j, IsGraphChart ι J (Φ t₀) c₀ (p j) (U j) (σ j) (φ j)) ∧
        (∀ j, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
          ‖jet 3 (φ j) x‖ ≤ (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁) ∧
        (∀ j, IsCutoff (fun x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0) Λ) ∧
        ∀ s ∈ I, ∀ l : ℕ, energyM ι A J l (Φ s) =
          ∑ j, ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
            (if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0) ^ 2 *
              eDens J (Φ s ∘ σ j) l x := by
  set v : ℝ := (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal with hv
  have hvpos : 0 < v :=
    ENNReal.toReal_pos (Metric.measure_ball_pos _ _ (by norm_num)).ne'
      (measure_ball_lt_top).ne
  obtain ⟨Λ, hΛ⟩ := exists_partition_cutoff (n := n) (E := E) (M := M) ι hJi hJJ ⌊V₀ / v⌋₊
    (K₃ := (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁) hc₀.le (by linarith)
  refine ⟨⌊V₀ / v⌋₊, Λ, ?_⟩
  intro I Φ hΦ hA t₀ ht₀ hV
  obtain ⟨N, p, U, σ, φ, hvol, hch, hcov⟩ :=
    exists_graphChart_cover ι A hJi hJJ hdim (hΦ.imm t₀ ht₀) hc₀ (by linarith) hc₁ hsmall
      (hA t₀ ht₀)
  have hN : N ≤ ⌊V₀ / v⌋₊ := by
    apply Nat.le_floor
    rw [le_div_iff₀ hvpos]
    exact hvol.trans hV
  have hjet := fun j => graphChart_jet3_le ι hJi hJJ (hΦ.imm t₀ ht₀) hc₀.le hc₁ (hA t₀ ht₀) (hch j)
  obtain ⟨ρ, hρs, hρ1, hρsupp, hρc⟩ := hΛ (hΦ.imm t₀ ht₀) N hN p U σ φ hch (fun j => hjet j) hcov
  refine ⟨N, p, U, σ, φ, ρ, hN, hch, hjet, hρc, ?_⟩
  intro s hs l
  rw [energyM_eq_sum ι A (hΦ.imm s hs) Finset.univ (fun j q => ρ j q ^ 2)
    (fun j _ => (hρs j).pow 2) (fun q => hρ1 q) (fun _ => ball 0 6) σ
    (fun j _ => (hch j).param)
    (fun j _ => (gs_tsupport_sq_subset (ρ j)).trans
      ((hρsupp j).trans (image_mono (ball_subset_ball (by norm_num)))))
    l]
  refine Finset.sum_congr rfl fun j _ => setIntegral_congr_fun measurableSet_ball ?_
  intro x hx
  simp [hx]

/-- info: 'CWFlow.Extension.exists_global_setup' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_global_setup

end CWFlow.Extension
