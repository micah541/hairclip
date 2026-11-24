import CWFlow.Extension.LimitChart
import CWFlow.Extension.LimitT
import CWFlow.Extension.LimitMAux
import CWFlow.Extension.Persist
import CWFlow.Extension.GraphChart

/-!
# The flow extends to the final time

Let `Φ` be a flow of Lagrangian immersions of the compact manifold `M` on `(a, T)` with uniformly
bounded second fundamental form and energies.  Then `Φ` extends to a flow on `[t₁, T]` for some
`t₁ < T` (`exists_limit_flow`):

1. `exists_persist` gives `τ₀`; pick `t₀ ∈ (T − τ₀, T)` and cover `M` by finitely many graph charts of
   `Φ t₀` (`exists_graphChart_cover`), each of which persists on `(t₀, T)`;
2. in each chart `exists_limit_chart` extends `Y (t, y) = Φ t (σ y)` smoothly to `[t₁, T] × B̄₂`;
3. `Φ̄ T m := lim_{t → T⁻} Φ t m`; in the chart `σ`, `Φ̄ T ∘ σ = Y (T, ·)`;
4. `limit_chart_at_T` gives the bound `‖Dχ − 1‖ ≤ 1/8`, the Lagrangian condition and the velocity at `T`;
5. `lm_local_T` (in `LimitMAux`) transfers all this to an arbitrary chart of `M` through the smooth
   transition map `κ_j ∘ chartInv ι z`.
-/

open Metric Set Filter
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem lm_exists_bound {N : ℕ} (f : Fin N → ℝ) {c T : ℝ} (hc : c < T) (hf : ∀ j, f j < T) :
    ∃ m : ℝ, c ≤ m ∧ m < T ∧ ∀ j, f j ≤ m := by
  classical
  let S : Finset ℝ := insert c (Finset.univ.image f)
  have hne : S.Nonempty := ⟨c, Finset.mem_insert_self _ _⟩
  refine ⟨S.max' hne, Finset.le_max' S c (Finset.mem_insert_self _ _), ?_,
    fun j => Finset.le_max' S _ (by simp [S])⟩
  have h := Finset.max'_mem S hne
  simp only [S, Finset.mem_insert, Finset.mem_image, Finset.mem_univ, true_and] at h
  rcases h with h | ⟨j, hj⟩
  · rw [h]; exact hc
  · rw [← hj]; exact hf j

theorem lm_exists_ext {α V : Type*} [NormedAddCommGroup V] (Φ : ℝ → α → V) (T : ℝ) :
    ∃ Φbar : ℝ → α → V, (∀ t, t ≠ T → Φbar t = Φ t) ∧
      ∀ m v, Filter.Tendsto (fun s => Φ s m) (𝓝[<] T) (𝓝 v) → Φbar T m = v := by
  classical
  refine ⟨fun t m => if t = T then limUnder (𝓝[<] T) (fun s => Φ s m) else Φ t m,
    fun t ht => funext fun m => if_neg ht, fun m v hv => ?_⟩
  simp only [if_true]
  exact hv.limUnder_eq

/-- Assembly of the extended flow from the facts at the final time. -/
theorem lm_glue {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {a T t₁ : ℝ} {Φ Φbar : ℝ → M → V}
    (hΦ : IsFlowM ι J (Ioo a T) Φ) (ht₁ : t₁ ∈ Ioo a T) (hP1 : ∀ t, t ≠ T → Φbar t = Φ t)
    (hAtT : ∀ z : M, ∀ y ∈ chartDomE ι z,
      ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z) (T, y) ∧
      ContDiffAt ℝ ∞ (Φbar T ∘ chartInv ι z) y ∧
      Function.Injective (fderiv ℝ (Φbar T ∘ chartInv ι z) y) ∧
      (∀ u w, ⟪J (fderiv ℝ (Φbar T ∘ chartInv ι z) y u),
        fderiv ℝ (Φbar T ∘ chartInv ι z) y w⟫ = 0) ∧
      HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar T ∘ chartInv ι z) y)
        (Icc t₁ T) T) :
    ∃ t₁ ∈ Ioo a T, ∃ Φbar : ℝ → M → V,
      (∀ t ∈ Ioo a T, Φbar t = Φ t) ∧
      (∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z)) ∧
      (∀ t ∈ Icc t₁ T, IsLagImm ι J (Φbar t)) ∧
      (∀ t ∈ Icc t₁ T, ∀ z : M, ∀ y ∈ chartDomE ι z,
        HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar t ∘ chartInv ι z) y)
          (Icc t₁ T) t) := by
  classical
  have hIoo' : ∀ t, t ∈ Icc t₁ T → t ≠ T → t ∈ Ioo a T := fun t ht h =>
    ⟨by linarith [ht.1, ht₁.1], lt_of_le_of_ne ht.2 h⟩
  refine ⟨t₁, ht₁, Φbar,
    fun t ht => hP1 t ht.2.ne, ?_, ?_, ?_⟩
  · intro z
    rintro ⟨t, y⟩ ⟨ht, hy⟩
    rcases eq_or_ne t T with h | h
    · have := (hAtT z y hy).1
      rw [h]
      exact this
    · have ht' := hIoo' t ht h
      have hS : ContDiffWithinAt ℝ ∞
          (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (chartInv ι z q.2))
          (Icc t₁ T ×ˢ chartDomE ι z) (t, y) :=
        ((hΦ.smooth z).contDiffAt ((isOpen_Ioo.prod (im_isOpen_chartDomE ι z)).mem_nhds
          ⟨ht', hy⟩)).contDiffWithinAt
      refine hS.congr_of_eventuallyEq ?_ (by simp only [hP1 t h])
      have hev : ∀ᶠ q : ℝ × EuclideanSpace ℝ (Fin n) in nhds (t, y), q.1 ≠ T :=
        (continuous_fst.tendsto (t, y)).eventually (isOpen_ne.mem_nhds h)
      filter_upwards [nhdsWithin_le_nhds hev] with q hq
      simp only [hP1 q.1 hq]
  · intro t ht
    rcases eq_or_ne t T with h | h
    · rw [h]
      exact ⟨fun z y hy => (hAtT z y hy).2.1.contDiffWithinAt, fun z y hy => (hAtT z y hy).2.2.1,
        fun z y hy => (hAtT z y hy).2.2.2.1⟩
    · rw [hP1 t h]
      exact hΦ.imm t (hIoo' t ht h)
  · intro t ht z y hy
    rcases eq_or_ne t T with h | h
    · rw [h]
      exact (hAtT z y hy).2.2.2.2
    · have ht' := hIoo' t ht h
      have hd := hΦ.vel t ht' z y hy
      have hd' : HasDerivAt (fun s => Φbar s (chartInv ι z y))
          (velX J (Φ t ∘ chartInv ι z) y) t := by
        refine hd.congr_of_eventuallyEq ?_
        filter_upwards [isOpen_ne.mem_nhds h] with s hs
        rw [hP1 s hs]
      rw [hP1 t h]
      exact hd'.hasDerivWithinAt


theorem exists_limit_flow {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hK₃ : (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ ≤ 1) (B : ℕ → ℝ)
    {a T : ℝ} (haT : a < T) {Φ : ℝ → M → V} (hΦ : IsFlowM ι J (Ioo a T) Φ)
    (hA : ∀ t ∈ Ioo a T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2)
    (hE : ∀ t ∈ Ioo a T, ∀ l, energyM ι A J l (Φ t) ≤ B l) :
    ∃ t₁ ∈ Ioo a T, ∃ Φbar : ℝ → M → V,
      (∀ t ∈ Ioo a T, Φbar t = Φ t) ∧
      (∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z)) ∧
      (∀ t ∈ Icc t₁ T, IsLagImm ι J (Φbar t)) ∧
      (∀ t ∈ Icc t₁ T, ∀ z : M, ∀ y ∈ chartDomE ι z,
        HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar t ∘ chartInv ι z) y)
          (Icc t₁ T) t) := by
  classical
  have hc₁1 : c₁ ≤ 1 / 2 := by
    have h1 : (1 : ℝ) ≤ (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) :=
      Real.one_le_rpow (by nlinarith [sq_nonneg c₀]) (by norm_num)
    nlinarith [mul_le_mul_of_nonneg_left h1 hc₁]
  obtain ⟨τ₀, hτ₀, hpersist⟩ := exists_persist ι A hJi hJJ hdim hc₀ hc₀1 hc₁ hsmall hK₃ B
  obtain ⟨t₀, ht₀, hTt₀⟩ : ∃ t₀ ∈ Ioo a T, T - t₀ < τ₀ := by
    refine ⟨max (T - τ₀ / 2) ((a + T) / 2), ⟨?_, ?_⟩, ?_⟩
    · exact lt_of_lt_of_le (by linarith) (le_max_right _ _)
    · exact max_lt (by linarith) (by linarith)
    · have := le_max_left (T - τ₀ / 2) ((a + T) / 2)
      linarith
  have ht₀T : t₀ < T := ht₀.2
  have hIoo : Ioo t₀ T ⊆ Ioo a T := Ioo_subset_Ioo_left ht₀.1.le
  obtain ⟨N, p, U, σ, φ, -, hch, hcov⟩ :=
    exists_graphChart_cover ι A hJi hJJ hdim (hΦ.imm t₀ ht₀) hc₀ (by linarith) hc₁ hsmall (hA t₀ ht₀)
  have hpers : ∀ j, ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
      ‖fderiv ℝ (fun y' => horC (U j) (Φ t₀ (p j)) (Φ t (σ j y'))) y -
        ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      ‖horC (U j) (Φ t₀ (p j)) (Φ t (σ j y)) - y‖ ≤ 1 / 8 ∧
      ‖verC J (U j) (Φ t₀ (p j)) (Φ t (σ j y))‖ ≤ 1 / 2 ∧
      ‖slopeZ J (U j) (fderiv ℝ (fun y' => Φ t (σ j y')) y)‖ ≤ 2 * c₀ := by
    intro j t ht y hy
    have h := hpersist hΦ hA hE ht₀ (hch j) t ⟨ht.1.le, by linarith [ht.2, hTt₀]⟩ ht.2 y hy
    exact ⟨h.1, h.2.1, h.2.2.1, h.2.2.2.2⟩
  have hlim : ∀ j, ∃ t₁ ∈ Ioo t₀ T, ∃ Ybar : ℝ × EuclideanSpace ℝ (Fin n) → V,
      ContDiffOn ℝ ∞ Ybar (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) ∧
      EqOn Ybar (fun q => Φ q.1 (σ j q.2)) (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    fun j => exists_limit_chart ι A hJi hJJ hdim hc₀.le hc₀1 hc₁ hc₁1 B hΦ ht₀T hIoo
      (fun t ht z => hA t (hIoo ht) z) (fun t ht => hE t (hIoo ht)) (hch j) (hpers j)
  choose t₁j ht₁j Ybar hYb hYeq using hlim
  obtain ⟨m, hm0, hmT, hmj⟩ := lm_exists_bound t₁j ht₀T (fun j => (ht₁j j).2)
  set t₁ : ℝ := (m + T) / 2 with ht₁def
  have ht₁T : t₁ < T := by rw [ht₁def]; linarith
  have ht₁0 : t₀ < t₁ := by rw [ht₁def]; linarith
  have ht₁j' : ∀ j, t₁j j < t₁ := fun j => by have := hmj j; rw [ht₁def]; linarith
  have hTt : T ∈ Icc t₁ T := right_mem_Icc.mpr ht₁T.le
  have hYb' : ∀ j, ContDiffOn ℝ ∞ (Ybar j) (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    fun j => (hYb j).mono (prod_mono (Icc_subset_Icc_left (ht₁j' j).le) subset_rfl)
  have hYeq' : ∀ j, EqOn (Ybar j) (fun q => Φ q.1 (σ j q.2))
      (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    fun j => (hYeq j).mono (prod_mono (Ioo_subset_Ioo_left (ht₁j' j).le) subset_rfl)
  -- the extension
  obtain ⟨Φbar, hP1, hlimΦ⟩ := lm_exists_ext (V := V) Φ T
  have hP2 : ∀ j, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, Φbar T (σ j y) = Ybar j (T, y) := by
    intro j y hy
    have hmem : ((T, y) : ℝ × EuclideanSpace ℝ (Fin n)) ∈
        Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2 :=
      ⟨hTt, ball_subset_closedBall hy⟩
    have hmaps : MapsTo (fun s : ℝ => ((s, y) : ℝ × EuclideanSpace ℝ (Fin n))) (Icc t₁ T)
        (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) :=
      fun s hs => ⟨hs, ball_subset_closedBall hy⟩
    have hcont : ContinuousWithinAt (fun s : ℝ => Ybar j (s, y)) (Icc t₁ T) T :=
      ContinuousWithinAt.comp (g := Ybar j) (f := fun s : ℝ => ((s, y) : ℝ × EuclideanSpace ℝ (Fin n)))
        ((hYb' j).continuousOn (T, y) hmem)
        (continuous_id.prodMk continuous_const).continuousWithinAt hmaps
    have htend : Filter.Tendsto (fun s : ℝ => Ybar j (s, y)) (𝓝[<] T) (𝓝 (Ybar j (T, y))) := by
      rw [← nhdsWithin_Ioo_eq_nhdsLT ht₁T]
      exact (hcont.mono Ioo_subset_Icc_self).tendsto
    refine hlimΦ _ _ ?_
    refine Filter.Tendsto.congr' ?_ htend
    filter_upwards [Ioo_mem_nhdsLT ht₁T] with s hs
    exact hYeq' j (show ((s, y) : ℝ × EuclideanSpace ℝ (Fin n)) ∈ _ from ⟨hs, hy⟩)
  have hYF : ∀ j, ∀ s ∈ Icc t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      Φbar s (σ j y) = Ybar j (s, y) := by
    intro j s hs y hy
    rcases eq_or_ne s T with h | h
    · rw [h]; exact hP2 j y hy
    · rw [hP1 s h]
      exact (hYeq j (show ((s, y) : ℝ × EuclideanSpace ℝ (Fin n)) ∈ _ from
        ⟨⟨lt_of_lt_of_le (ht₁j' j) hs.1, lt_of_le_of_ne hs.2 h⟩, hy⟩)).symm
  -- the facts at the final time, in the graph charts
  have hTf : ∀ j, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖horCL (U j) ∘L fderiv ℝ (fun y' => Ybar j (T, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      (∀ u w, ⟪J (fderiv ℝ (fun y' => Ybar j (T, y')) y u),
        fderiv ℝ (fun y' => Ybar j (T, y')) y w⟫ = 0) ∧
      HasDerivWithinAt (fun s => Ybar j (s, y)) (velX J (fun y' => Ybar j (T, y')) y) (Icc t₁ T) T :=
    fun j => lm_chart_T ι hΦ ht₀.1 ht₁0 ht₁T (hch j).param
      (fun t ht y hy => (hpers j t ht y hy).1) (hYb' j) (hYeq' j)
  have hAtT : ∀ z : M, ∀ y ∈ chartDomE ι z,
      ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z) (T, y) ∧
      ContDiffAt ℝ ∞ (Φbar T ∘ chartInv ι z) y ∧
      Function.Injective (fderiv ℝ (Φbar T ∘ chartInv ι z) y) ∧
      (∀ u w, ⟪J (fderiv ℝ (Φbar T ∘ chartInv ι z) y u),
        fderiv ℝ (Φbar T ∘ chartInv ι z) y w⟫ = 0) ∧
      HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar T ∘ chartInv ι z) y)
        (Icc t₁ T) T := by
    intro z y hy
    obtain ⟨j, x₀, hx₀, hσx⟩ := hcov (chartInv ι z y)
    exact lm_local_T ι (F := Φbar) (hch j).param ht₁T (hYb' j) (hYF j) (horCL (U j))
      (fun y hy => (hTf j y hy).1) (fun y hy => (hTf j y hy).2.1) (fun y hy => (hTf j y hy).2.2)
      hy (ball_subset_ball (by norm_num) hx₀) hσx
  exact lm_glue ι hΦ ⟨by linarith [ht₁0, ht₀.1], ht₁T⟩ hP1 hAtT

/-- info: 'CWFlow.Extension.exists_limit_flow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_limit_flow

end CWFlow.Extension
