import CWFlow.Scalar.Smooth

/-!
# Invariance under a linear change of the domain space

For a continuous linear equivalence `A : E' ≃L E` between two (possibly different) domain spaces, the flow
quantities of `X ∘ A` at `y` are those of `X` at `A y`. Since `A` is linear, `D²(X ∘ A) = D²X(A·, A·)`
with no first-order correction. This is the analogue of `Invariance.lean`/`Reparam.lean` with `S = 0`.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E E' : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup E'] [NormedSpace ℝ E'] [FiniteDimensional ℝ E']

theorem fderiv_comp_linearEquiv (A : E' ≃L[ℝ] E) (X : E → V) (y : E') :
    fderiv ℝ (X ∘ A) y = (fderiv ℝ X (A y)).comp (A : E' →L[ℝ] E) :=
  A.comp_right_fderiv

/-- `D²(X ∘ A)(a, c) = D²X(A a, A c)`. -/
theorem fderiv_fderiv_comp_linearEquiv (A : E' ≃L[ℝ] E) {X : E → V} {y : E'}
    (hX : DifferentiableAt ℝ (fderiv ℝ X) (A y)) (a c : E') :
    fderiv ℝ (fderiv ℝ (X ∘ A)) y a c = fderiv ℝ (fderiv ℝ X) (A y) (A a) (A c) := by
  let Φ : (E →L[ℝ] V) →L[ℝ] (E' →L[ℝ] V) :=
    (ContinuousLinearMap.compL ℝ E' E V).flip (A : E' →L[ℝ] E)
  have e : fderiv ℝ (X ∘ A) = Φ ∘ (fderiv ℝ X ∘ A) :=
    funext fun w => fderiv_comp_linearEquiv A X w
  have h := Φ.hasFDerivAt.comp y (hX.hasFDerivAt.comp y A.hasFDerivAt)
  rw [e, h.fderiv]
  simp [Φ]

variable {L : E →L[ℝ] V}

theorem apply_ginvOf_comp_linearEquiv (hL : Function.Injective L) (A : E' ≃L[ℝ] E)
    (ξ : E' →L[ℝ] ℝ) :
    A (ginvOf (L.comp (A : E' →L[ℝ] E)) ξ) = ginvOf L (ξ.comp (A.symm : E →L[ℝ] E')) := by
  apply gOf_injective hL
  rw [gOf_ginvOf hL]
  ext a
  have h := congrArg (fun f => f (A.symm a))
    (gOf_ginvOf (L := L.comp (A : E' →L[ℝ] E)) (hL.comp A.injective) ξ)
  simpa using h

theorem trOf_comp_linearEquiv (hL : Function.Injective L) (A : E' ≃L[ℝ] E) {W : Type*}
    [AddCommGroup W] [Module ℝ W] (B : E →ₗ[ℝ] E →ₗ[ℝ] W) :
    trOf (L.comp (A : E' →L[ℝ] E)) (fun a c => B (A a) (A c)) = trOf L (fun a c => B a c) := by
  simp only [trOf, apply_ginvOf_comp_linearEquiv hL A]
  refine sum_frame_eq (fun i => A (bas E' i)) (fun i => (dualBas i).comp (A.symm : E →L[ℝ] E'))
    (fun v => ?_) B (ginvOf L : (E →L[ℝ] ℝ) →ₗ[ℝ] E)
  have := congrArg A (sum_dualBas_smul (A.symm v))
  simpa [map_sum, map_smul] using this

theorem perpOf_comp_linearEquiv (hL : Function.Injective L) (A : E' ≃L[ℝ] E) (v : V) :
    perpOf (L.comp (A : E' →L[ℝ] E)) v = perpOf L v := by
  simp only [perpOf, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [apply_ginvOf_comp_linearEquiv hL A]
  congr 3
  ext a
  simp

/-- `H` for `L ∘ A` and `Q(A·, A·)` is `H` for `L` and `Q`. -/
theorem HOf_comp_linearEquiv (hL : Function.Injective L) (A : E' ≃L[ℝ] E) (Q : Hess V E)
    (Q' : Hess V E') (hQ : ∀ a c, Q' a c = Q (A a) (A c)) :
    HOf (L.comp (A : E' →L[ℝ] E)) Q' = HOf L Q := by
  have hB : ∀ a c, perpOf (L.comp (A : E' →L[ℝ] E)) (Q' a c) =
      ((Q : E →L[ℝ] E →L[ℝ] V).toLinearMap₁₂.compr₂ (perpLin L)) (A a) (A c) := fun a c => by
    rw [perpOf_comp_linearEquiv hL A, hQ]
    rfl
  have hB' : ∀ a c, perpOf L (Q a c) =
      ((Q : E →L[ℝ] E →L[ℝ] V).toLinearMap₁₂.compr₂ (perpLin L)) a c := fun a c => rfl
  simp only [HOf, hB, hB']
  exact trOf_comp_linearEquiv hL A _

variable (A : E' ≃L[ℝ] E) {X : E → V} (hX : ContDiff ℝ ∞ X)
  (hinj : ∀ z, Function.Injective (fderiv ℝ X z))
include hX hinj

theorem meanCurv_comp_linearEquiv (y : E') : meanCurv (X ∘ A) y = meanCurv X (A y) := by
  have hD : DifferentiableAt ℝ (fderiv ℝ X) (A y) :=
    (hX.fderiv_right (m := ∞) le_rfl).differentiable (by simp) _
  rw [meanCurv, meanCurv, fderiv_comp_linearEquiv A X y]
  exact HOf_comp_linearEquiv (hinj _) A _ _ fun a c => fderiv_fderiv_comp_linearEquiv A hD a c

theorem flowPot_comp_linearEquiv (J : V →L[ℝ] V) (y : E') :
    flowPot J (X ∘ A) y = flowPot J X (A y) := by
  have hL := fderiv_comp_linearEquiv A X y
  have hev : (fun w => J (meanCurv (X ∘ A) w)) = (fun v => J (meanCurv X v)) ∘ A := by
    funext w
    simp only [Function.comp_apply, meanCurv_comp_linearEquiv A hX hinj w]
  have hM : fderiv ℝ (fun w => J (meanCurv (X ∘ A) w)) y =
      (fderiv ℝ (fun v => J (meanCurv X v)) (A y)).comp (A : E' →L[ℝ] E) := by
    rw [hev, fderiv_comp y (differentiableAt_J_meanCurv J hX.contDiffAt (hinj _))
      A.differentiableAt, A.fderiv]
  set M := fderiv ℝ (fun v => J (meanCurv X v)) (A y)
  set L' := fderiv ℝ X (A y)
  let B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ := LinearMap.mk₂ ℝ (fun a c => ⟪M a, L' c⟫)
    (fun a a' c => by simp [inner_add_left]) (fun r a c => by simp [inner_smul_left])
    (fun a c c' => by simp [inner_add_right]) (fun r a c => by simp [inner_smul_right])
  have h := trOf_comp_linearEquiv (hinj (A y)) A B
  simp only [flowPot, hL, hM]
  exact h

theorem flowVel_comp_linearEquiv (J : V →L[ℝ] V) (y : E') :
    flowVel J (X ∘ A) y = flowVel J X (A y) := by
  have hpot : DifferentiableAt ℝ (flowPot J X) (A y) :=
    ((contDiffOn_flowPot isOpen_univ hX.contDiffOn (fun z _ => hinj z) J).contDiffAt
      Filter.univ_mem).differentiableAt (by simp)
  have hL := fderiv_comp_linearEquiv A X y
  have hD : fderiv ℝ (flowPot J (X ∘ A)) y =
      (fderiv ℝ (flowPot J X) (A y)).comp (A : E' →L[ℝ] E) := by
    have hev : flowPot J (X ∘ A) = flowPot J X ∘ A :=
      funext fun w => flowPot_comp_linearEquiv A hX hinj J w
    rw [hev, fderiv_comp y hpot A.differentiableAt, A.fderiv]
  simp only [flowVel, hL, hD, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [apply_ginvOf_comp_linearEquiv (hinj _) A]
  congr 3
  ext a
  simp

theorem shrinker_comp_linearEquiv (J : V →L[ℝ] V) {c : ℝ} {y : E'}
    (h : flowVel J X (A y) = -c • perpOf (fderiv ℝ X (A y)) (X (A y))) :
    flowVel J (X ∘ A) y = -c • perpOf (fderiv ℝ (X ∘ A) y) ((X ∘ A) y) := by
  rw [flowVel_comp_linearEquiv A hX hinj J y, h, fderiv_comp_linearEquiv A X y,
    perpOf_comp_linearEquiv (hinj _) A]
  rfl

end CWFlow
