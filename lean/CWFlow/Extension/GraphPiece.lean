import CWFlow.Extension.ChartSup

/-!
# Pointwise bounds for a graph piece not centred at the origin

`exists_const_chart_sup` bounds all jets of order `≥ 1` of a graph chart centred at the origin of its frame.
Here the same is proved for an arbitrary Lagrangian graph piece over `ball x₀ 2`, with the slope bound
`‖∇φ‖ ≤ 1/2` and `‖D²φ‖, ‖D³φ‖ < 2` (the value of `φ` is arbitrary):

* `integral_eSym_le_energyM_ball`: the integral of the energy symbol over `ball x₀ R` of a graph
  parametrization defined on `ball x₀ R₂` (`R < R₂`) is bounded by the global energy;
* `exists_const_piece_sup`: `‖D^m φ (x₀)‖ ≤ C` with `C` depending only on `m, n, B`.

The proof translates `x₀` to the origin and subtracts the constant `φ(x₀)` (neither is seen by the energy
density or by jets of order `≥ 1`), after which the Euclidean lemmas centred at `0` apply.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

/-! ### Translation invariance of the coordinate energy density -/

section Shift

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem gp_fderiv_shift (X : EuclideanSpace ℝ (Fin n) → V) (a : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (fun x => X (a + x)) = fun x => fderiv ℝ X (a + x) := by
  funext x; exact fderiv_comp_add_left a

theorem gp_dX_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) (b : Fin n) :
    dX (fun x => X (a + x)) x b = dX X (a + x) b := by
  unfold dX; rw [gp_fderiv_shift]

theorem gp_ddX_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) (b c : Fin n) :
    ddX (fun x => X (a + x)) x b c = ddX X (a + x) b c := by
  unfold ddX
  rw [gp_fderiv_shift, fderiv_comp_add_left (f := fderiv ℝ X) a]

theorem gp_gX_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) :
    gX (fun x => X (a + x)) x = gX X (a + x) := by
  ext b c
  unfold gX
  rw [gp_dX_shift, gp_dX_shift]

theorem gp_ginvZ_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) :
    ginvZ (fun x => X (a + x)) x = ginvZ X (a + x) := by
  unfold ginvZ
  rw [gp_gX_shift]

theorem gp_volX_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) :
    volX (fun x => X (a + x)) x = volX X (a + x) := by
  unfold volX
  rw [gp_gX_shift]

theorem gp_gamX_shift (X : EuclideanSpace ℝ (Fin n) → V) (a x : EuclideanSpace ℝ (Fin n)) :
    gamX (fun x => X (a + x)) x = gamX X (a + x) := by
  unfold gamX
  simp only [gp_gX_shift, gp_ddX_shift, gp_dX_shift]

theorem gp_cX_shift (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (a x : EuclideanSpace ℝ (Fin n)) :
    cX J (fun x => X (a + x)) x = cX J X (a + x) := by
  unfold cX
  simp only [gp_ddX_shift, gp_dX_shift]

theorem gp_SX_shift (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (a : EuclideanSpace ℝ (Fin n)) (k : ℕ) :
    SX J (fun x => X (a + x)) k = fun x => SX J X k (a + x) := by
  induction k with
  | zero => funext x; exact gp_cX_shift J X a x
  | succ k ih =>
    funext x
    show SX J _ (k + 1) x = SX J X (k + 1) (a + x)
    rw [SX_succ, SX_succ, ih, gp_gamX_shift, fderiv_comp_add_left (f := SX J X k) a]

theorem gp_eDens_shift (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (a : EuclideanSpace ℝ (Fin n)) (k : ℕ) (x : EuclideanSpace ℝ (Fin n)) :
    eDens J (fun x => X (a + x)) k x = eDens J X k (a + x) := by
  unfold eDens
  rw [gp_SX_shift, gp_ginvZ_shift, gp_volX_shift]

end Shift

/-! ### Jets under translation and subtraction of a constant -/

theorem gp_jet_shift {n : ℕ} (f : EuclideanSpace ℝ (Fin n) → ℝ) (a : EuclideanSpace ℝ (Fin n)) (j : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) : jet j (fun v => f (a + v)) x = jet j f (a + x) := by
  unfold jet
  rw [iteratedFDeriv_comp_add_left]

theorem gp_jet_sub_const {n : ℕ} {g : EuclideanSpace ℝ (Fin n) → ℝ} (hg : ContDiff ℝ ∞ g) (c : ℝ)
    (j : ℕ) (hj : 1 ≤ j) (x : EuclideanSpace ℝ (Fin n)) : jet j (fun v => g v - c) x = jet j g x := by
  obtain ⟨i, rfl⟩ : ∃ i, j = i + 1 := ⟨j - 1, by omega⟩
  have hgi : ContDiffAt ℝ ((i + 1 : ℕ) : WithTop ℕ∞) g x :=
    (hg.of_le (by exact_mod_cast le_top)).contDiffAt
  have hc : ContDiffAt ℝ ((i + 1 : ℕ) : WithTop ℕ∞) (fun _ : EuclideanSpace ℝ (Fin n) => c) x :=
    contDiffAt_const
  have h := iteratedFDeriv_sub_apply hgi hc
  unfold jet
  have h1 : (fun v => g v - c) = g - fun _ => c := rfl
  rw [h1, h, iteratedFDeriv_succ_const]
  simp

theorem gp_integral_translate {n : ℕ} (f : EuclideanSpace ℝ (Fin n) → ℝ)
    (x₀ : EuclideanSpace ℝ (Fin n)) (r : ℝ) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) r, f (x₀ + x) = ∫ y in ball x₀ r, f y := by
  have h := (measurePreserving_add_left (volume : Measure (EuclideanSpace ℝ (Fin n))) x₀
    ).setIntegral_preimage_emb (Homeomorph.addLeft x₀).measurableEmbedding f (ball x₀ r)
  have hp : (fun x : EuclideanSpace ℝ (Fin n) => x₀ + x) ⁻¹' ball x₀ r = ball 0 r := by
    ext x
    simp [mem_ball, dist_eq_norm]
  rw [hp] at h
  exact h

theorem gp_integral_sq_low {n : ℕ} {ψ : EuclideanSpace ℝ (Fin n) → ℝ} {j : ℕ} {r : ℝ}
    (h : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r, ‖jet j ψ y‖ < 2) :
    ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) r, ‖jet j ψ y‖ ^ 2 ≤
      4 * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) r)).toReal := by
  have := norm_setIntegral_le_of_norm_le_const_ae (μ := volume)
    (s := ball (0 : EuclideanSpace ℝ (Fin n)) r) (f := fun y => ‖jet j ψ y‖ ^ 2) (C := 4)
    measure_ball_lt_top
    (ae_restrict_of_forall_mem measurableSet_ball fun y hy => by
      have h1 := h y hy
      have h2 := norm_nonneg (jet j ψ y)
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      nlinarith)
  refine (le_abs_self _).trans ?_
  rw [← Real.norm_eq_abs]
  simpa [measureReal_def] using this

/-! ### Chart integrals of the energy symbol, centre `x₀` -/

/-- A chart integral of the energy symbol over `ball x₀ R` is bounded by the global energy, for a graph
parametrization on the larger ball `ball x₀ R₂`. -/
theorem integral_eSym_le_energyM_ball {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}
    {y₀ : V} {F : M → V} [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M)
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (hF : IsLagImm ι J F) {x₀ : EuclideanSpace ℝ (Fin n)} {R R₂ : ℝ}
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι (ball x₀ R₂) σ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball x₀ R₂))
    (hgraph : ∀ x ∈ ball x₀ R₂, F (σ x) = graphX J U y₀ φ x)
    (hR : R < R₂) (l N : ℕ) (hN : l + 3 ≤ N) :
    ∫ x in ball x₀ R, eSym l (jp N φ x) ≤ energyM ι A J l F := by
  have hnd := isNatDens_densOf_eDens (n := n) J l
  have hd0 : ∀ (X : EuclideanSpace ℝ (Fin n) → V) y, 0 ≤ eDens J X l y :=
    fun X y => ec_eDens_nonneg J X l y
  by_cases hR0 : R ≤ 0
  · rw [Metric.ball_eq_empty.2 hR0, Measure.restrict_empty, integral_zero_measure]
    have h0 : 0 ≤ intM ι A F (fun _ => (1 : ℝ)) (densOf fun X => eDens J X l) :=
      intM_densOf_nonneg ι A (fun _ => zero_le_one) hd0
    exact h0
  have hR0 : 0 < R := not_le.1 hR0
  -- the cutoff
  set R₁ : ℝ := (R + R₂) / 2 with hR₁
  have h1 : R < R₁ := by rw [hR₁]; linarith
  have h2 : R₁ < R₂ := by rw [hR₁]; linarith
  let bump : ContDiffBump x₀ := ⟨R, R₁, hR0, h1⟩
  set η : EuclideanSpace ℝ (Fin n) → ℝ := fun x => bump x with hηdef
  have hη : ContDiff ℝ ∞ η := bump.contDiff
  have hη0 : ∀ x, x ∉ closedBall x₀ R₁ → η x = 0 := by
    intro x hx
    have : x ∉ Function.support η := by
      rw [hηdef, bump.support_eq]
      exact fun hm => hx (ball_subset_closedBall hm)
    simpa using this
  have hη1 : ∀ x ∈ closedBall x₀ R, η x = 1 :=
    fun x hx => bump.one_of_mem_closedBall hx
  have hηn : ∀ x, 0 ≤ η x := fun x => bump.nonneg
  have hηl : ∀ x, η x ≤ 1 := fun x => bump.le_one
  have hKB : closedBall x₀ R₁ ⊆ ball x₀ R₂ := closedBall_subset_ball h2
  obtain ⟨κ, hκc, hκ⟩ := hσ.inv
  set h : M → ℝ := ec_ext (ball x₀ R₂) σ η κ with hh
  have hhs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h :=
    ec_contMDiff_ext ι hσ (isCompact_closedBall _ _) hKB hη hη0 hκc hκ
  have hσx : ∀ x ∈ ball x₀ R₂, h (σ x) = η x := by
    intro x hx
    rw [hh, ec_ext_of_mem ⟨x, hx, rfl⟩, hκ x hx]
  have hh1 : ∀ p, h p ≤ 1 := by
    intro p
    by_cases hp : p ∈ σ '' ball x₀ R₂
    · rw [hh, ec_ext_of_mem hp]; exact hηl _
    · rw [hh, ec_ext_of_notMem hp]; exact zero_le_one
  have hsupp : tsupport h ⊆ σ '' ball x₀ R₂ := by
    have hKc : IsClosed (σ '' closedBall x₀ R₁) :=
      ((isCompact_closedBall _ _).image_of_continuousOn (hσ.cont.mono hKB)).isClosed
    have hsub : Function.support h ⊆ σ '' closedBall x₀ R₁ := by
      intro p hp
      by_cases hpB : p ∈ σ '' ball x₀ R₂
      · obtain ⟨x, hx, rfl⟩ := hpB
        refine ⟨x, ?_, rfl⟩
        by_contra hxK
        apply hp
        rw [hσx x hx]; exact hη0 x hxK
      · exact absurd (by rw [hh, ec_ext_of_notMem hpB]) (Function.mem_support.1 hp)
    exact (closure_minimal hsub hKc).trans (Set.image_mono hKB)
  -- the integrals
  have hFp := hF.comp_param ι hσ
  have hcont : ContinuousOn (eDens J (F ∘ σ) l) (ball x₀ R₂) :=
    (contDiffOn_eDens isOpen_ball hFp.1 hFp.2.1 l).continuousOn
  have hfc : ContinuousOn (fun x => η x * eDens J (F ∘ σ) l x) (ball x₀ R₂) :=
    hη.continuous.continuousOn.mul hcont
  have hfi : IntegrableOn (fun x => η x * eDens J (F ∘ σ) l x) (ball x₀ R₂) :=
    im_integrableOn_of_vanish isOpen_ball hfc (isCompact_closedBall _ _) hKB
      (fun y _ hy => by rw [hη0 y hy, zero_mul])
  have hint : intM ι A F h (densOf fun X => eDens J X l) =
      ∫ x in ball x₀ R₂, η x * eDens J (F ∘ σ) l x := by
    rw [intM_eq_param ι A hF hhs hnd hσ hsupp]
    refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
    simp only [densOf, Function.comp_apply]
    rw [hσx x hx]
  have hmono := intM_densOf_mono ι A hF hnd hd0 hhs (h₂ := fun _ => 1) contMDiff_const hh1
  have hballR : ball x₀ R ⊆ ball x₀ R₂ := ball_subset_ball hR.le
  calc ∫ x in ball x₀ R, eSym l (jp N φ x)
      = ∫ x in ball x₀ R, η x * eDens J (F ∘ σ) l x := by
        refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
        rw [hη1 x (ball_subset_closedBall hx), one_mul,
          eDens_eq_eSym_of_graph hJi hJJ hU isOpen_ball hφ hgraph hN (hballR hx)]
    _ ≤ ∫ x in ball x₀ R₂, η x * eDens J (F ∘ σ) l x :=
        setIntegral_mono_set hfi (Filter.Eventually.of_forall fun x => mul_nonneg (hηn x) (hd0 _ _))
          hballR.eventuallyLE
    _ = intM ι A F h (densOf fun X => eDens J X l) := hint.symm
    _ ≤ intM ι A F (fun _ => 1) (densOf fun X => eDens J X l) := hmono
    _ = energyM ι A J l F := rfl

/-! ### The pointwise bound -/

theorem exists_const_piece_sup {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {m : ℕ} (hm : 1 ≤ m) (B : ℕ → ℝ) :
    ∃ C : ℝ, ∀ {F : M → V}, IsLagImm ι J F →
      (∀ l ≤ m + n, energyM ι A J l F ≤ B l) →
      ∀ {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
        {φ : EuclideanSpace ℝ (Fin n) → ℝ} {x₀ : EuclideanSpace ℝ (Fin n)},
      IsLagFrame J U → IsParam ι (ball x₀ 2) σ → ContDiffOn ℝ ∞ φ (ball x₀ 2) →
      (∀ x ∈ ball x₀ 2, F (σ x) = graphX J U y₀ φ x ∧
        ‖gradE φ x‖ ≤ 1 / 2 ∧ ‖jet 2 φ x‖ < 2 ∧ ‖jet 3 φ x‖ < 2) →
      ‖jet m φ x₀‖ ≤ C := by
  obtain ⟨C_S, hCS0, hCS⟩ := exists_const_sup_jet (n := n) m (ρ := 1 / 8) (R := 1 / 4)
    (by norm_num) (by norm_num)
  choose C4 hC4nn hC4 using fun q : ℕ => exists_const_energy3_le_eSym (n := n) q (N := q + 3) le_rfl
    (ρ := 1 / 4) (R := 1 / 2) (by norm_num) (by norm_num) (K := 2) (by norm_num)
  set S : ℝ := ∑ l ∈ Finset.range (m + n + 1), max (B l) 0 with hS
  set vol : ℝ := (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4))).toReal with hvol
  set D : ℝ := 4 * vol + ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) with hD
  refine ⟨C_S * (1 + (n + 1) * D), fun {F} hF hE {U y₀ σ φ x₀} hU hσ hφ hgr => ?_⟩
  have hS0 : 0 ≤ S := Finset.sum_nonneg fun l _ => le_max_right _ _
  have hvol0 : 0 ≤ vol := ENNReal.toReal_nonneg
  have hDsum : 0 ≤ ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) :=
    Finset.sum_nonneg fun q _ => mul_nonneg (hC4nn q) (by linarith)
  have hD0 : 0 ≤ D := by rw [hD]; positivity
  have hmem : ∀ (r : ℝ) (v : EuclideanSpace ℝ (Fin n)), x₀ + v ∈ ball x₀ r ↔ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r := by
    intro r v
    simp [mem_ball, dist_eq_norm]
  have hgraph : ∀ x ∈ ball x₀ 2, F (σ x) = graphX J U y₀ φ x := fun x hx => (hgr x hx).1
  -- the global extension, translated to the origin
  have hφ' : ContDiffOn ℝ ∞ (fun v => φ (x₀ + v)) (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    hφ.comp (contDiff_const.add contDiff_id).contDiffOn fun v hv => (hmem 2 v).2 hv
  obtain ⟨R', hR', g, hg, -, hgφ⟩ := exists_global_of_contDiffOn
    (U := ball (0 : EuclideanSpace ℝ (Fin n)) 2) isOpen_ball (R := 7 / 4) (by norm_num)
    (closedBall_subset_ball (by norm_num)) hφ'
  obtain ⟨ψ, hψdef⟩ : ∃ ψ : EuclideanSpace ℝ (Fin n) → ℝ, ψ = fun v => g v - φ x₀ := ⟨_, rfl⟩
  have hψ : ContDiff ℝ ∞ ψ := by rw [hψdef]; exact hg.sub contDiff_const
  have hR'ball : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4),
      v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' := fun v hv => ball_subset_ball hR'.le hv
  have h2ball : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4), x₀ + v ∈ ball x₀ 2 :=
    fun v hv => (hmem 2 v).2 (ball_subset_ball (by norm_num) hv)
  have hev : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', g =ᶠ[𝓝 v] fun v => φ (x₀ + v) := by
    intro v hv
    filter_upwards [isOpen_ball.mem_nhds hv] with z hz using hgφ z hz
  have hfd : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', fderiv ℝ ψ v = fderiv ℝ φ (x₀ + v) := by
    intro v hv
    have h1 : fderiv ℝ ψ v = fderiv ℝ g v := by rw [hψdef]; exact fderiv_sub_const _
    rw [h1, (hev v hv).fderiv_eq]
    exact fderiv_comp_add_left x₀
  have hjet : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', ∀ j : ℕ, 1 ≤ j →
      jet j ψ v = jet j φ (x₀ + v) := by
    intro v hv j hj
    rw [hψdef, gp_jet_sub_const hg _ j hj v, jet_congr_of_eventuallyEq (hev v hv) j,
      gp_jet_shift φ x₀ j v]
  have hψ0 : ψ 0 = 0 := by
    rw [hψdef]
    have := hgφ 0 (mem_ball_self (by linarith))
    simp only [add_zero] at this
    simp [this]
  -- admissibility
  have hgrad : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4), ‖jet 1 ψ y‖ ≤ 1 / 2 := by
    intro y hy
    rw [hjet y (hR'ball y hy) 1 le_rfl, cb_norm_jet_one]
    exact (hgr (x₀ + y) (h2ball y hy)).2.1
  have hlow : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (3 / 2), ∀ j : ℕ, j ≤ 3 → ‖jet j ψ v‖ < 2 := by
    intro v hv j hj
    have hv7 : v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4) := ball_subset_ball (by norm_num) hv
    have hvn : ‖v‖ < 3 / 2 := by simpa using hv
    interval_cases j
    · rw [norm_jet_zero, Real.norm_eq_abs]
      have hb : ∀ y ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2), ‖fderiv ℝ ψ y‖ ≤ 1 / 2 := by
        intro y hy
        have hy' : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4) :=
          closedBall_subset_ball (by norm_num) hy
        calc ‖fderiv ℝ ψ y‖ ≤ ‖gradE ψ y‖ := cb_norm_fderiv_le _ _
          _ = ‖jet 1 ψ y‖ := (cb_norm_jet_one _ _).symm
          _ ≤ 1 / 2 := hgrad y hy'
      have hd : ∀ y ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2), DifferentiableAt ℝ ψ y :=
        fun y _ => hψ.contDiffAt.differentiableAt (by simp)
      have := (convex_closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2)).norm_image_sub_le_of_norm_fderiv_le
        hd hb (mem_closedBall_self (by norm_num)) (ball_subset_closedBall hv)
      rw [hψ0, sub_zero, sub_zero, Real.norm_eq_abs] at this
      nlinarith [norm_nonneg v]
    · have := hgrad v hv7
      linarith
    · rw [hjet v (hR'ball v hv7) 2 (by norm_num)]
      exact (hgr (x₀ + v) (h2ball v hv7)).2.2.1
    · rw [hjet v (hR'ball v hv7) 3 (by norm_num)]
      exact (hgr (x₀ + v) (h2ball v hv7)).2.2.2
  -- the energy symbol does not see the translation or the constant
  have hpt : ∀ (q l : ℕ), l ≤ q → ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4),
      eSym l (jp (q + 3) ψ x) = eSym l (jp (q + 3) φ (x₀ + x)) := by
    intro q l hlq x hx
    have hgraph' : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4),
        F (σ (x₀ + v)) = graphX J U (y₀ + U x₀) ψ v := by
      intro v hv
      rw [hgraph (x₀ + v) (h2ball v hv)]
      have hgg : gradE ψ v = gradE φ (x₀ + v) := by
        unfold gradE; rw [hfd v (hR'ball v hv)]
      unfold graphX
      rw [hgg, map_add]
      abel
    have e1 := eDens_eq_eSym_of_graph (J := J) (U := U) (y₀ := y₀ + U x₀) (F := F)
      (σ := fun v => σ (x₀ + v)) (B := ball (0 : EuclideanSpace ℝ (Fin n)) (7 / 4)) hJi hJJ hU
      isOpen_ball hψ.contDiffOn hgraph' (N := q + 3) (l := l) (by omega) hx
    have e2 := eDens_eq_eSym_of_graph (J := J) (U := U) (y₀ := y₀) (F := F) (σ := σ)
      (B := ball x₀ 2) hJi hJJ hU isOpen_ball hφ hgraph (N := q + 3) (l := l) (by omega)
      (h2ball x hx)
    have e3 := gp_eDens_shift J (F ∘ σ) x₀ l x
    calc eSym l (jp (q + 3) ψ x) = eDens J (F ∘ fun v => σ (x₀ + v)) l x := e1.symm
      _ = eDens J (F ∘ σ) l (x₀ + x) := e3
      _ = _ := e2
  have hlocal : ∀ q l : ℕ, l ≤ q →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2), eSym l (jp (q + 3) ψ x) ≤ energyM ι A J l F := by
    intro q l hlq
    calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2), eSym l (jp (q + 3) ψ x)
        = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2), eSym l (jp (q + 3) φ (x₀ + x)) :=
          setIntegral_congr_fun measurableSet_ball fun x hx =>
            hpt q l hlq x (ball_subset_ball (by norm_num) hx)
      _ = ∫ y in ball x₀ (1 / 2), eSym l (jp (q + 3) φ y) :=
          gp_integral_translate (fun y => eSym l (jp (q + 3) φ y)) x₀ (1 / 2)
      _ ≤ energyM ι A J l F :=
          integral_eSym_le_energyM_ball ι A hJi hJJ hU hF hσ hφ hgraph (by norm_num : (1 / 2 : ℝ) < 2)
            l (q + 3) (by omega)
  have hbd : ∀ j ≤ m + n, ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4), ‖jet j ψ y‖ ^ 2 ≤ D := by
    intro j hj
    by_cases hj3 : j < 3
    · have := gp_integral_sq_low (n := n) (ψ := ψ) (j := j) (r := 1 / 4) fun y hy =>
        hlow y (ball_subset_ball (by norm_num) hy) j (by omega)
      rw [hD]
      linarith
    · obtain ⟨q, rfl⟩ : ∃ q, j = 3 + q := ⟨j - 3, by omega⟩
      have hq : q ≤ m + n := by omega
      have h1 : ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4), ‖jet (3 + q) ψ y‖ ^ 2 ≤
          energy3 q (1 / 4) ψ :=
        Finset.single_le_sum (f := fun i => ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4),
          ‖jet (3 + i) ψ y‖ ^ 2) (fun i _ => setIntegral_nonneg measurableSet_ball fun y _ => by positivity)
          (Finset.mem_range.2 (Nat.lt_succ_self q))
      have hadm : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2 + 1),
          jp (q + 3) ψ x ∈ admSet n (q + 3) 2 := by
        intro x hx k hk
        exact hlow x (ball_subset_ball (by norm_num) hx) k hk
      have h2 := hC4 q ψ hψ hadm
      have h3' : ∑ l ∈ Finset.range (q + 1),
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2), eSym l (jp (q + 3) ψ x) ≤ S := by
        calc ∑ l ∈ Finset.range (q + 1),
              ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 2), eSym l (jp (q + 3) ψ x)
            ≤ ∑ l ∈ Finset.range (q + 1), max (B l) 0 :=
              Finset.sum_le_sum fun l hl => by
                have hlq : l ≤ q := by have := Finset.mem_range.1 hl; omega
                exact (hlocal q l hlq).trans ((hE l (by omega)).trans (le_max_left _ _))
          _ ≤ S := Finset.sum_le_sum_of_subset_of_nonneg
              (Finset.range_subset_range.2 (by omega)) fun l _ _ => le_max_right _ _
      have h4 : energy3 q (1 / 4) ψ ≤ C4 q * (1 + S) :=
        h2.trans (mul_le_mul_of_nonneg_left (by linarith) (hC4nn q))
      have h5 : C4 q * (1 + S) ≤ ∑ q ∈ Finset.range (m + n + 1), C4 q * (1 + S) :=
        Finset.single_le_sum (f := fun q => C4 q * (1 + S))
          (fun q _ => mul_nonneg (hC4nn q) (by linarith)) (Finset.mem_range.2 (by omega))
      rw [hD]
      linarith
  have h0mem : (0 : EuclideanSpace ℝ (Fin n)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) :=
    mem_ball_self (by norm_num)
  have hx' := hCS ψ hψ 0 h0mem
  have hsum : ∑ i ∈ Finset.range (n + 1),
      ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4), ‖jet (m + i) ψ y‖ ^ 2 ≤ (n + 1) * D := by
    have := Finset.sum_le_card_nsmul (Finset.range (n + 1))
      (fun i => ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 4), ‖jet (m + i) ψ y‖ ^ 2) D
      fun i hi => hbd (m + i) (by have := Finset.mem_range.1 hi; omega)
    simpa using this
  calc ‖jet m φ x₀‖ = ‖jet m ψ 0‖ := by
        rw [hjet 0 (mem_ball_self (by linarith)) m hm, add_zero]
    _ ≤ _ := hx'
    _ ≤ C_S * (1 + (n + 1) * D) := mul_le_mul_of_nonneg_left (by linarith) hCS0

/-- info: 'CWFlow.Extension.integral_eSym_le_energyM_ball' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_eSym_le_energyM_ball

/-- info: 'CWFlow.Extension.exists_const_piece_sup' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_piece_sup

end CWFlow.Extension
