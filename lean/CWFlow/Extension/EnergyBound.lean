import CWFlow.Extension.GlobalIneq

/-!
# Boundedness of the energies up to the final time

Along a flow `Φ` of Lagrangian immersions of a compact manifold on a finite interval with
`|A| ≤ c₁`, all energies `E_k = ∫ |∇^k C|² dV` stay bounded on `[t₁, T)`.
`E_0 ≤ c₁² Vol` (`energyM_zero_le`), the volume is nonincreasing (`volM_antitoneOn`), and the
global differential inequality `E_k' ≤ C_k (1 + Σ_{l ≤ k} E_l)` (`exists_const_energy_ineq`) gives
the bounds by the Gronwall induction `bounded_of_deriv_le_sum`.
-/

open Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

/-- `E_0 = ∫ |A|² dV ≤ c₁² Vol`.  (Deviation: the immersion hypothesis `hF` is added, needed for
integrability of the chart integrands.) -/
theorem energyM_zero_le {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    {F : M → V} (hF : IsLagImm ι J F) {c₁ : ℝ}
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2) :
    energyM ι A J 0 F ≤ c₁ ^ 2 * volM ι A F := by
  unfold energyM volM intM
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun z hz => ?_
  have hsupp : tsupport (A.θ z) ⊆ chartInv ι z '' chartDomE ι z := fun p hp =>
    ⟨chartFwd ι z p, im_chartFwd_mem ι (A.supp z hz hp), im_chartInv_chartFwd ι (A.supp z hz hp)⟩
  have hmeas := (im_isOpen_chartDomE ι z).measurableSet
  have I1 := im_integrableOn_param ι hF (A.smooth z hz) (isNatDens_densOf_eDens (n := n) J 0)
    (isParam_chart ι z) hsupp
  have I2 := im_integrableOn_param ι hF (A.smooth z hz) (isNatDens_densOf_volX (n := n) J)
    (isParam_chart ι z) hsupp
  have I1' : IntegrableOn (fun y => A.θ z (chartInv ι z y) *
      densOf (fun X => eDens J X 0) (F ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
      (chartDomE ι z) := by
    refine I1.congr_fun (fun y _ => ?_) hmeas
    simp [densOf]
  have I2' : IntegrableOn (fun y => c₁ ^ 2 * (A.θ z (chartInv ι z y) *
      densOf volX (F ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y))
      (chartDomE ι z) := by
    refine IntegrableOn.congr_fun (I2.const_mul (c₁ ^ 2)) (fun y _ => ?_) hmeas
    simp [densOf]
  rw [← integral_const_mul]
  refine setIntegral_mono_on I1' I2' hmeas fun y hy => ?_
  have hθ := A.nonneg z (chartInv ι z y)
  have hv : 0 ≤ volX (F ∘ chartInv ι z) y := Real.sqrt_nonneg _
  have he : eDens J (F ∘ chartInv ι z) 0 y = aNormSq J (F ∘ chartInv ι z) y *
      volX (F ∘ chartInv ι z) y := rfl
  simp only [densOf, Function.comp_apply, one_mul, he]
  have := mul_le_mul_of_nonneg_right (hA z y hy) hv
  nlinarith [mul_le_mul_of_nonneg_left this hθ]

theorem eb_volM_nonneg {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) (F : M → V) :
    0 ≤ volM ι A F :=
  intM_densOf_nonneg ι A (fun _ => zero_le_one) fun X y => Real.sqrt_nonneg _

theorem eb_energyM_nonneg {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) (J : V →L[ℝ] V) (k : ℕ)
    (F : M → V) : 0 ≤ energyM ι A J k F :=
  intM_densOf_nonneg ι A (fun _ => zero_le_one) fun X y => ec_eDens_nonneg J X k y

theorem energyM_bounded {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hK₃ : (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ ≤ 1)
    {a T : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J (Ioo a T) Φ)
    (hA : ∀ t ∈ Ioo a T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {t₁ : ℝ} (ht₁ : t₁ ∈ Ioo a T) :
    ∀ k, ∃ B : ℝ, ∀ t ∈ Ico t₁ T, energyM ι A J k (Φ t) ≤ B := by
  classical
  set V₀ : ℝ := volM ι A (Φ t₁) with hV₀
  have hmem : ∀ t ∈ Ico t₁ T, t ∈ Ioo a T := fun t ht => ⟨by linarith [ht₁.1, ht.1], ht.2⟩
  have hvol : ∀ t ∈ Ico t₁ T, volM ι A (Φ t) ≤ V₀ := fun t ht =>
    volM_antitoneOn ι A hJi hJJ hΦ ht₁ (hmem t ht) ht.1
  have hex : ∀ k : ℕ, ∃ C : ℝ, 0 ≤ C ∧ (1 ≤ k → ∀ t ∈ Ioo a T, volM ι A (Φ t) ≤ V₀ →
      ∃ D : ℝ, HasDerivAt (fun s => energyM ι A J k (Φ s)) D t ∧
        D ≤ C * (1 + ∑ l ∈ Finset.range (k + 1), energyM ι A J l (Φ t))) := by
    intro k
    by_cases hk : 1 ≤ k
    · obtain ⟨C, hC0, hC⟩ := exists_const_energy_ineq ι A hJi hJJ hdim k hk hc₀ hc₀1 hc₁ hsmall
        hK₃ V₀
      exact ⟨C, hC0, fun _ t ht hv => hC hΦ hA ht hv⟩
    · exact ⟨0, le_rfl, fun h => absurd h hk⟩
  choose C hC0 hC using hex
  have hE0 : ∀ t ∈ Ico t₁ T, energyM ι A J 0 (Φ t) ≤ c₁ ^ 2 * V₀ := fun t ht =>
    (energyM_zero_le ι A (hΦ.imm t (hmem t ht)) (hA t (hmem t ht))).trans
      (mul_le_mul_of_nonneg_left (hvol t ht) (sq_nonneg _))
  have hV₀nn : 0 ≤ V₀ := eb_volM_nonneg ι A _
  let E' : ℕ → ℝ → ℝ := fun k t => if k = 0 then c₁ ^ 2 * V₀ else energyM ι A J k (Φ t)
  have hE'pos : ∀ k, k ≠ 0 → E' k = fun t => energyM ι A J k (Φ t) := fun k hk => by
    funext t; simp [E', hk]
  have hmain := bounded_of_deriv_le_sum (a := t₁) (T := T) ht₁.2 E' (fun k t => deriv (E' k) t) C hC0
    ⟨c₁ ^ 2 * V₀, fun t _ => by simp [E']⟩
    (fun k t ht => by
      by_cases hk : k = 0
      · simp only [E', hk, if_true]; exact mul_nonneg (sq_nonneg _) hV₀nn
      · simp only [E', hk, if_false]; exact eb_energyM_nonneg ι A J k _)
    (fun k => by
      by_cases hk : k = 0
      · subst hk; simp only [E', if_true]; exact continuousOn_const
      · intro t ht
        obtain ⟨D, hD, -⟩ := hC k (Nat.one_le_iff_ne_zero.2 hk) t (hmem t ht) (hvol t ht)
        rw [hE'pos k hk]
        exact hD.continuousAt.continuousWithinAt)
    (fun k hk t ht => by
      have hk0 : k ≠ 0 := by omega
      have ht' : t ∈ Ioo a T := ⟨by linarith [ht₁.1, ht.1], ht.2⟩
      obtain ⟨D, hD, -⟩ := hC k hk t ht' (hvol t ⟨ht.1.le, ht.2⟩)
      rw [hE'pos k hk0]
      rw [← hD.deriv] at hD
      exact hD)
    (fun k hk t ht => by
      have hk0 : k ≠ 0 := by omega
      have ht' : t ∈ Ioo a T := ⟨by linarith [ht₁.1, ht.1], ht.2⟩
      obtain ⟨D, hD, hle⟩ := hC k hk t ht' (hvol t ⟨ht.1.le, ht.2⟩)
      have hdd : deriv (E' k) t = D := by rw [hE'pos k hk0]; exact hD.deriv
      rw [hdd]
      refine hle.trans (mul_le_mul_of_nonneg_left ?_ (hC0 k))
      have hsum : ∑ l ∈ Finset.range (k + 1), energyM ι A J l (Φ t) ≤
          ∑ l ∈ Finset.range (k + 1), E' l t := by
        refine Finset.sum_le_sum fun l _ => ?_
        by_cases hl : l = 0
        · subst hl; simpa [E'] using hE0 t ⟨ht.1.le, ht.2⟩
        · simp [E', hl]
      linarith)
  intro k
  by_cases hk : k = 0
  · subst hk
    exact ⟨c₁ ^ 2 * V₀, hE0⟩
  · obtain ⟨B, hB⟩ := hmain k
    exact ⟨B, fun t ht => by simpa [E', hk] using hB t ht⟩

/-- info: 'CWFlow.Extension.energyM_zero_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms energyM_zero_le

/-- info: 'CWFlow.Extension.energyM_bounded' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms energyM_bounded

end CWFlow.Extension
