import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.MeasureTheory.Measure.Haar.InnerProductSpace
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.MeasureTheory.Integral.Pi
import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Tactic

/-!
# Sup bound by the integral of the top derivative

For `u ∈ C_c^∞(ℝⁿ; V)`, `‖u x‖ ≤ ∫ ‖Dⁿu(y)‖ dy`: the value at `x` is bounded by the integral
of the mixed partial `∂₁⋯∂ₙ u`, which is a value of `Dⁿu` on standard basis vectors.
The proof is by induction on `n` on `Fin n → ℝ` (one-dimensional fundamental theorem of calculus
in the first variable, then Fubini), transferred to `EuclideanSpace ℝ (Fin n)` at the end.
-/

open MeasureTheory Set
open scoped ContDiff

namespace CWFlow.Interp

section OneDim

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] [CompleteSpace V]

theorem cons_eq_add (t : ℝ) (x' : Fin n → ℝ) :
    (Fin.cons t x' : Fin (n + 1) → ℝ) = Fin.cons 0 x' + t • (Pi.single 0 1 : Fin (n + 1) → ℝ) := by
  ext i
  refine Fin.cases ?_ (fun j => ?_) i <;> simp

theorem contDiff_cons_slice (x' : Fin n → ℝ) :
    ContDiff ℝ ∞ (fun t : ℝ => (Fin.cons t x' : Fin (n + 1) → ℝ)) := by
  rw [show (fun t : ℝ => (Fin.cons t x' : Fin (n + 1) → ℝ)) =
    fun t => Fin.cons 0 x' + t • (Pi.single 0 1 : Fin (n + 1) → ℝ) from funext fun t => cons_eq_add t x']
  exact contDiff_const.add (contDiff_id.smul contDiff_const)

theorem hasDerivAt_cons_slice (x' : Fin n → ℝ) (t : ℝ) :
    HasDerivAt (fun t : ℝ => (Fin.cons t x' : Fin (n + 1) → ℝ)) (Pi.single 0 1) t := by
  rw [show (fun t : ℝ => (Fin.cons t x' : Fin (n + 1) → ℝ)) =
    fun t => Fin.cons 0 x' + t • (Pi.single 0 1 : Fin (n + 1) → ℝ) from funext fun t => cons_eq_add t x']
  simpa using ((hasDerivAt_id t).smul_const (Pi.single 0 1 : Fin (n + 1) → ℝ)).const_add
    (Fin.cons 0 x' : Fin (n + 1) → ℝ)

theorem norm_cons_le_integral {u : (Fin (n + 1) → ℝ) → V} (hu : ContDiff ℝ ∞ u)
    (hc : HasCompactSupport u) (x' : Fin n → ℝ) (t : ℝ) :
    ‖u (Fin.cons t x')‖ ≤
      ∫ τ, ‖fderiv ℝ u (Fin.cons τ x') (Pi.single 0 1)‖ := by
  set s : ℝ → V := fun t => u (Fin.cons t x') with hs
  have hs1 : ContDiff ℝ 1 s :=
    (hu.of_le (by exact_mod_cast le_top)).comp (contDiff_cons_slice x' |>.of_le (by exact_mod_cast le_top))
  have hsc : HasCompactSupport s := by
    refine HasCompactSupport.intro (hc.image (continuous_apply 0)) fun t ht => ?_
    by_contra h
    apply ht
    refine ⟨Fin.cons t x', subset_tsupport _ h, by simp⟩
  have hds : ∀ τ, deriv s τ = fderiv ℝ u (Fin.cons τ x') (Pi.single 0 1) := by
    intro τ
    exact ((hu.differentiable (by simp)).differentiableAt.hasFDerivAt.comp_hasDerivAt τ
      (hasDerivAt_cons_slice x' τ)).deriv
  have hint : Integrable (fun τ => ‖deriv s τ‖) := by
    have h1 : Continuous (deriv s) := hs1.continuous_deriv le_rfl
    exact (h1.norm).integrable_of_hasCompactSupport hsc.deriv.norm
  calc ‖s t‖ = ‖∫ τ in Iic t, deriv s τ‖ := by rw [hsc.integral_Iic_deriv_eq hs1]
    _ ≤ ∫ τ in Iic t, ‖deriv s τ‖ := norm_integral_le_integral_norm _
    _ ≤ ∫ τ, ‖deriv s τ‖ :=
        setIntegral_le_integral hint (Filter.Eventually.of_forall fun _ => norm_nonneg _)
    _ = _ := by simp_rw [hds]

end OneDim

section Step

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- The embedding `y ↦ (0, y)` of `Fin n → ℝ` into `Fin (n+1) → ℝ` as a continuous linear map. -/
noncomputable def consL (n : ℕ) : (Fin n → ℝ) →L[ℝ] (Fin (n + 1) → ℝ) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun y => Fin.cons 0 y
      map_add' := fun a b => by
        ext i; refine Fin.cases ?_ (fun j => ?_) i <;> simp
      map_smul' := fun c a => by
        ext i; refine Fin.cases ?_ (fun j => ?_) i <;> simp }

theorem consL_apply (y : Fin n → ℝ) : consL n y = (Fin.cons 0 y : Fin (n + 1) → ℝ) := rfl

theorem cons_eq_add_consL (τ : ℝ) (y : Fin n → ℝ) :
    (Fin.cons τ y : Fin (n + 1) → ℝ) = Fin.cons τ 0 + consL n y := by
  ext i; refine Fin.cases ?_ (fun j => ?_) i <;> simp [consL_apply]

theorem iteratedFDeriv_deriv0 {u : (Fin (n + 1) → ℝ) → V} (hu : ContDiff ℝ ∞ u)
    (z : Fin (n + 1) → ℝ) (m : Fin n → Fin (n + 1) → ℝ) :
    iteratedFDeriv ℝ n (fun z => fderiv ℝ u z (Pi.single 0 1)) z m =
      iteratedFDeriv ℝ (n + 1) u z (Fin.snoc (α := fun _ => Fin (n + 1) → ℝ) m (Pi.single 0 1)) := by
  have h1 : (fun z => fderiv ℝ u z (Pi.single 0 1)) =
      (ContinuousLinearMap.apply ℝ V (Pi.single 0 1 : Fin (n + 1) → ℝ)) ∘ (fderiv ℝ u) := rfl
  have hd : ContDiff ℝ ∞ (fderiv ℝ u) := hu.fderiv_right (by simp)
  rw [h1, ContinuousLinearMap.iteratedFDeriv_comp_left _
    ((hd.of_le (by exact_mod_cast le_top)).contDiffAt) le_rfl,
    iteratedFDeriv_succ_apply_right]
  simp

theorem iteratedFDeriv_slice {u : (Fin (n + 1) → ℝ) → V} (hu : ContDiff ℝ ∞ u)
    (τ : ℝ) (y : Fin n → ℝ) (m : Fin n → Fin n → ℝ) :
    iteratedFDeriv ℝ n (fun y' => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) y m =
      iteratedFDeriv ℝ (n + 1) u (Fin.cons τ y)
        (Fin.snoc (α := fun _ => Fin (n + 1) → ℝ) (fun k => (Fin.cons 0 (m k) : Fin (n + 1) → ℝ))
          (Pi.single 0 1)) := by
  set g : (Fin (n + 1) → ℝ) → V := fun z => fderiv ℝ u z (Pi.single 0 1) with hg
  have hd : ContDiff ℝ ∞ (fderiv ℝ u) := hu.fderiv_right (by simp)
  have hgs : ContDiff ℝ ∞ g := hd.clm_apply contDiff_const
  have h1 : (fun y' => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) =
      (fun z => g ((Fin.cons τ 0 : Fin (n + 1) → ℝ) + z)) ∘ (consL n) := by
    funext y'; simp only [Function.comp, hg, ← cons_eq_add_consL]
  have hG : ContDiff ℝ ∞ (fun z => g ((Fin.cons τ 0 : Fin (n + 1) → ℝ) + z)) :=
    hgs.comp (contDiff_const.add contDiff_id)
  rw [h1, ContinuousLinearMap.iteratedFDeriv_comp_right _ (hG.of_le (by exact_mod_cast le_top))
    _ le_rfl]
  simp only [ContinuousMultilinearMap.compContinuousLinearMap_apply]
  rw [iteratedFDeriv_comp_add_left, ← cons_eq_add_consL, ← iteratedFDeriv_deriv0 hu]
  rfl

end Step

section Induction

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] [CompleteSpace V]

/-- The induction statement: there is a list `m` of standard basis vectors of `Fin n → ℝ` such
that `‖u x‖ ≤ ∫ ‖Dⁿu(y) m‖` for all smooth compactly supported `u`. -/
def SupStmt (V : Type*) [NormedAddCommGroup V] [NormedSpace ℝ V] (n : ℕ) : Prop :=
  ∃ m : Fin n → Fin n → ℝ, (∀ k, ∃ j, m k = Pi.single j 1) ∧
    ∀ u : (Fin n → ℝ) → V, ContDiff ℝ ∞ u → HasCompactSupport u → ∀ x,
      ‖u x‖ ≤ ∫ y, ‖iteratedFDeriv ℝ n u y m‖

theorem supStmt_zero : SupStmt V 0 := by
  refine ⟨fun k => k.elim0, fun k => k.elim0, fun u _ _ x => ?_⟩
  have hv : (volume : Measure (Fin 0 → ℝ)) univ = 1 := by
    rw [volume_pi, Measure.pi_of_empty]; simp
  have hx : x = default := Subsingleton.elim _ _
  simp only [iteratedFDeriv_zero_apply]
  rw [integral_unique]
  rw [hx, Measure.real, hv]
  simp
  exact le_of_eq (congrArg (fun y => ‖u y‖) (Subsingleton.elim _ _))

theorem cons_zero_single (j : Fin n) :
    (Fin.cons 0 (Pi.single j 1 : Fin n → ℝ) : Fin (n + 1) → ℝ) = Pi.single j.succ 1 := by
  ext i; refine Fin.cases ?_ (fun l => ?_) i <;> simp [Pi.single_apply]

theorem contDiff_slice_deriv {u : (Fin (n + 1) → ℝ) → V} (hu : ContDiff ℝ ∞ u) (τ : ℝ) :
    ContDiff ℝ ∞ (fun y' : Fin n → ℝ => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) := by
  have hd : ContDiff ℝ ∞ (fderiv ℝ u) := hu.fderiv_right (by simp)
  have hgs : ContDiff ℝ ∞ (fun z => fderiv ℝ u z (Pi.single 0 1)) :=
    hd.clm_apply contDiff_const
  have : (fun y' : Fin n → ℝ => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) =
      (fun z => fderiv ℝ u z (Pi.single 0 1)) ∘ (fun y' => Fin.cons τ 0 + consL n y') := by
    funext y'; simp only [Function.comp, ← cons_eq_add_consL]
  rw [this]
  exact hgs.comp (contDiff_const.add (consL n).contDiff)

theorem hasCompactSupport_slice_deriv {u : (Fin (n + 1) → ℝ) → V} (hc : HasCompactSupport u)
    (τ : ℝ) :
    HasCompactSupport (fun y' : Fin n → ℝ => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) := by
  refine HasCompactSupport.intro (K := (fun z : Fin (n + 1) → ℝ => Fin.tail z) '' tsupport u)
    (hc.image (continuous_pi fun i => continuous_apply _)) fun y hy => ?_
  by_contra h
  apply hy
  refine ⟨Fin.cons τ y, ?_, by simp⟩
  have hne : fderiv ℝ u (Fin.cons τ y) ≠ 0 := fun h0 => h (by simp [h0])
  exact support_fderiv_subset ℝ hne

theorem supStmt_succ (h : SupStmt V n) : SupStmt V (n + 1) := by
  obtain ⟨m, hm, H⟩ := h
  refine ⟨Fin.snoc (α := fun _ => Fin (n + 1) → ℝ) (fun k => (Fin.cons 0 (m k) : Fin (n + 1) → ℝ))
    (Pi.single 0 1), ?_, ?_⟩
  · intro k
    refine Fin.lastCases ⟨0, by simp⟩ (fun k => ?_) k
    obtain ⟨j, hj⟩ := hm k
    exact ⟨j.succ, by simp [hj, cons_zero_single]⟩
  intro u hu hc x
  set M : Fin (n + 1) → Fin (n + 1) → ℝ :=
    Fin.snoc (α := fun _ => Fin (n + 1) → ℝ) (fun k => (Fin.cons 0 (m k) : Fin (n + 1) → ℝ))
      (Pi.single 0 1) with hM
  set F : (Fin (n + 1) → ℝ) → ℝ := fun z => ‖iteratedFDeriv ℝ (n + 1) u z M‖ with hF
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) 0 with he
  have hFc : Continuous F :=
    ((ContinuousMultilinearMap.apply ℝ (fun _ => Fin (n + 1) → ℝ) V M).continuous.comp
      (hu.continuous_iteratedFDeriv (by exact_mod_cast le_top))).norm
  have hFs : HasCompactSupport F :=
    hc.mono' fun z hz => support_iteratedFDeriv_subset (𝕜 := ℝ) (n + 1) (by
      intro h0; apply hz; simp [hF, h0])
  have hint : Integrable F := hFc.integrable_of_hasCompactSupport hFs
  have hmp : MeasurePreserving e.symm volume volume :=
    (volume_preserving_piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) 0).symm
  have hG : Integrable (fun p : ℝ × (Fin n → ℝ) => F (e.symm p))
      ((volume : Measure ℝ).prod (volume : Measure (Fin n → ℝ))) :=
    hmp.integrable_comp_of_integrable hint
  have hsymm : ∀ (t : ℝ) (y : Fin n → ℝ), e.symm (t, y) = Fin.cons t y := by
    intro t y; simp [he]; rfl
  have htot : ∫ z, F z = ∫ t : ℝ, ∫ y : Fin n → ℝ, F (Fin.cons t y) := by
    have h2 : ∫ z, F z = ∫ p, F (e.symm p) ∂((volume : Measure ℝ).prod (volume : Measure (Fin n → ℝ))) :=
      (hmp.integral_comp' F).symm
    rw [h2, integral_prod _ hG]
    simp only [hsymm]
  have hd : ContDiff ℝ ∞ (fderiv ℝ u) := hu.fderiv_right (by simp)
  have h1 := norm_cons_le_integral hu hc (Fin.tail x) (x 0)
  rw [Fin.cons_self_tail] at h1
  refine h1.trans ?_
  rw [htot]
  apply integral_mono
  · have hcont : Continuous fun τ : ℝ => ‖fderiv ℝ u (Fin.cons τ (Fin.tail x)) (Pi.single 0 1)‖ :=
      ((hd.continuous.comp (contDiff_cons_slice (Fin.tail x)).continuous).clm_apply
        continuous_const).norm
    refine hcont.integrable_of_hasCompactSupport ?_
    refine HasCompactSupport.intro (hc.image (continuous_apply 0)) fun t ht => ?_
    by_contra h
    apply ht
    have hne : fderiv ℝ u (Fin.cons t (Fin.tail x)) ≠ 0 := fun h0 => h (by simp [h0])
    exact ⟨Fin.cons t (Fin.tail x), support_fderiv_subset ℝ hne, by simp⟩
  · have := hG.integral_prod_left
    simp only [hsymm] at this
    exact this
  · intro τ
    have := H (fun y' => fderiv ℝ u (Fin.cons τ y') (Pi.single 0 1)) (contDiff_slice_deriv hu τ)
      (hasCompactSupport_slice_deriv hc τ) (Fin.tail x)
    simp only [iteratedFDeriv_slice hu τ] at this
    exact this

end Induction

theorem supStmt (V : Type*) [NormedAddCommGroup V] [NormedSpace ℝ V] [CompleteSpace V] (n : ℕ) :
    SupStmt V n :=
  Nat.rec supStmt_zero (fun _ h => supStmt_succ h) n

theorem norm_euclid_single (n : ℕ) (j : Fin n) :
    ‖(EuclideanSpace.equiv (Fin n) ℝ).symm (Pi.single j 1)‖ = 1 := by
  have : (EuclideanSpace.equiv (Fin n) ℝ).symm (Pi.single j 1) = EuclideanSpace.single j 1 := by
    ext i
    simp [EuclideanSpace.equiv, Pi.single_apply, EuclideanSpace.single_apply]
  rw [this]; simp

/-- For `u ∈ C_c^∞(ℝⁿ; V)`, `‖u x‖ ≤ ∫ ‖Dⁿu(y)‖ dy`. -/
theorem norm_le_integral_norm_iteratedFDeriv {n : ℕ} {V : Type*} [NormedAddCommGroup V]
    [NormedSpace ℝ V] [FiniteDimensional ℝ V]
    {u : EuclideanSpace ℝ (Fin n) → V} (hu : ContDiff ℝ ∞ u) (hc : HasCompactSupport u)
    (x : EuclideanSpace ℝ (Fin n)) :
    ‖u x‖ ≤ ∫ y, ‖iteratedFDeriv ℝ n u y‖ := by
  have : CompleteSpace V := FiniteDimensional.complete ℝ V
  obtain ⟨m, hm, H⟩ := supStmt V n
  set ψ := EuclideanSpace.equiv (Fin n) ℝ with hψ
  set G : EuclideanSpace ℝ (Fin n) → ℝ := fun z => ‖iteratedFDeriv ℝ n u z‖ with hG
  have hũ : ContDiff ℝ ∞ (u ∘ ψ.symm) := hu.comp ψ.symm.contDiff
  have hũc : HasCompactSupport (u ∘ ψ.symm) := hc.comp_homeomorph ψ.symm.toHomeomorph
  have h := H (u ∘ ψ.symm) hũ hũc (ψ x)
  simp only [Function.comp, ContinuousLinearEquiv.symm_apply_apply] at h
  have hGi : Integrable G := by
    refine (hu.continuous_iteratedFDeriv (by exact_mod_cast le_top)).norm.integrable_of_hasCompactSupport
      ?_
    exact (hc.iteratedFDeriv n).norm
  have hmp : MeasurePreserving (MeasurableEquiv.toLp 2 (Fin n → ℝ)) volume volume :=
    (EuclideanSpace.volume_preserving_symm_measurableEquiv_toLp (Fin n)).symm
  have hGi2 : Integrable (fun y : Fin n → ℝ => G (ψ.symm y)) :=
    hmp.integrable_comp_of_integrable hGi
  have hint : ∫ y : Fin n → ℝ, G (ψ.symm y) = ∫ z, G z := hmp.integral_comp' G
  refine h.trans ?_
  rw [← hint]
  refine integral_mono_of_nonneg (Filter.Eventually.of_forall fun _ => norm_nonneg _) hGi2
    (Filter.Eventually.of_forall fun y => ?_)
  have h1 : iteratedFDeriv ℝ n (u ∘ ψ.symm) y =
      (iteratedFDeriv ℝ n u (ψ.symm y)).compContinuousLinearMap fun _ => (ψ.symm : _ →L[ℝ] _) :=
    ContinuousLinearMap.iteratedFDeriv_comp_right (ψ.symm : _ →L[ℝ] _)
      (hu.of_le (by exact_mod_cast le_top)) y le_rfl
  show ‖(iteratedFDeriv ℝ n (u ∘ ψ.symm) y) m‖ ≤ ‖iteratedFDeriv ℝ n u (ψ.symm y)‖
  rw [h1]
  simp only [ContinuousMultilinearMap.compContinuousLinearMap_apply]
  refine (ContinuousMultilinearMap.le_opNorm _ _).trans ?_
  have : ∏ k, ‖((ψ.symm : (Fin n → ℝ) →L[ℝ] EuclideanSpace ℝ (Fin n)) (m k))‖ = 1 := by
    refine Finset.prod_eq_one fun k _ => ?_
    obtain ⟨j, hj⟩ := hm k
    rw [hj]; exact norm_euclid_single n j
  rw [this, mul_one]

/-- info: 'CWFlow.Interp.norm_le_integral_norm_iteratedFDeriv' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms norm_le_integral_norm_iteratedFDeriv

end CWFlow.Interp
