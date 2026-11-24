import CWFlow.Scalar.Pointwise

/-!
# Invariance of the jet geometry under reparametrization (algebraic part)

For an invertible `T : E ≃L E` and injective `L : E →L V`:

* the trace does not depend on the frame: `Σᵢ B(cᵢ, w δᵢ) = Σₗ B(bₗ, w εₗ)` for any frame `(cᵢ, δᵢ)`
  with `Σᵢ δᵢ(v) cᵢ = v` (`sum_frame_eq`);
* `T (g_{L∘T})⁻¹ ξ = g_L⁻¹ (ξ ∘ T⁻¹)` (`apply_ginvOf_comp`);
* `tr_{L∘T} B(T·, T·) = tr_L B` (`trOf_comp_equiv`);
* `perp_{L∘T} = perp_L` (`perpOf_comp_equiv`), and `perp_L (L a) = 0` (`perpOf_apply_self`);
* hence `H` is unchanged when `Q` is replaced by `Q(T·, T·) + L S(·, ·)` (`HOf_comp_equiv`); this is
  how `D²(X ∘ τ) = D²X(Dτ·, Dτ·) + DX(D²τ)` enters.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

theorem sum_dualBas_smul (v : E) : ∑ i, dualBas i v • bas E i = v :=
  (bas E).sum_repr v

theorem functional_expand (δ : E →L[ℝ] ℝ) : δ = ∑ l, δ (bas E l) • dualBas l := by
  ext v
  conv_lhs => rw [← sum_dualBas_smul v]
  simp [map_sum, mul_comm]

/-- The trace is independent of the frame. -/
theorem sum_frame_eq {ι : Type*} [Fintype ι] (c : ι → E) (δ : ι → E →L[ℝ] ℝ)
    (hcd : ∀ v, ∑ i, δ i v • c i = v) {W : Type*} [AddCommGroup W] [Module ℝ W]
    (B : E →ₗ[ℝ] E →ₗ[ℝ] W) (w : (E →L[ℝ] ℝ) →ₗ[ℝ] E) :
    ∑ i, B (c i) (w (δ i)) = ∑ l, B (bas E l) (w (dualBas l)) := by
  have h1 : ∀ l, B (bas E l) (w (dualBas l)) = ∑ i, δ i (bas E l) • B (c i) (w (dualBas l)) := by
    intro l
    conv_lhs => rw [← hcd (bas E l)]
    simp [map_sum, map_smul]
  have h2 : ∀ i, B (c i) (w (δ i)) = ∑ l, δ i (bas E l) • B (c i) (w (dualBas l)) := by
    intro i
    conv_lhs => rw [functional_expand (δ i)]
    simp [map_sum, map_smul]
  simp only [h1, h2]
  exact Finset.sum_comm

variable {L : E →L[ℝ] V}

theorem gOf_comp (T : E →L[ℝ] E) (a b : E) : gOf (L.comp T) a b = gOf L (T a) (T b) := rfl

theorem injective_comp_equiv (hL : Function.Injective L) (T : E ≃L[ℝ] E) :
    Function.Injective (L.comp (T : E →L[ℝ] E)) :=
  hL.comp T.injective

/-- `T (g_{L∘T})⁻¹ ξ = g_L⁻¹ (ξ ∘ T⁻¹)`. -/
theorem apply_ginvOf_comp (hL : Function.Injective L) (T : E ≃L[ℝ] E) (ξ : E →L[ℝ] ℝ) :
    T (ginvOf (L.comp (T : E →L[ℝ] E)) ξ) = ginvOf L (ξ.comp (T.symm : E →L[ℝ] E)) := by
  apply gOf_injective hL
  rw [gOf_ginvOf hL]
  ext a
  have h := congrArg (fun f => f (T.symm a))
    (gOf_ginvOf (injective_comp_equiv hL T) ξ)
  simp only [gOf_comp] at h
  simpa using h

/-- **Frame invariance of the trace:** `tr_{L∘T} B(T·, T·) = tr_L B`. -/
theorem trOf_comp_equiv (hL : Function.Injective L) (T : E ≃L[ℝ] E) {W : Type*} [AddCommGroup W]
    [Module ℝ W] (B : E →ₗ[ℝ] E →ₗ[ℝ] W) :
    trOf (L.comp (T : E →L[ℝ] E)) (fun a c => B (T a) (T c)) = trOf L (fun a c => B a c) := by
  simp only [trOf, apply_ginvOf_comp hL T]
  refine sum_frame_eq (fun i => T (bas E i)) (fun i => (dualBas i).comp (T.symm : E →L[ℝ] E))
    (fun v => ?_) B (ginvOf L : (E →L[ℝ] ℝ) →ₗ[ℝ] E)
  have := congrArg T (sum_dualBas_smul (T.symm v))
  simpa [map_sum, map_smul] using this

/-- The normal projection depends only on the image of `L`. -/
theorem perpOf_comp_equiv (hL : Function.Injective L) (T : E ≃L[ℝ] E) (v : V) :
    perpOf (L.comp (T : E →L[ℝ] E)) v = perpOf L v := by
  simp only [perpOf, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [apply_ginvOf_comp hL T]
  congr 3
  ext a
  simp

/-- `perp_L` kills the tangent space. -/
theorem perpOf_apply_self (hL : Function.Injective L) (a : E) : perpOf L (L a) = 0 := by
  have h : ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) (L a)).comp L = gOf L a := rfl
  have hg : ginvOf L (gOf L a) = a := by
    apply gOf_injective hL
    rw [gOf_ginvOf hL]
  rw [perpOf, h, hg, sub_self]

theorem perpOf_add (v w : V) : perpOf L (v + w) = perpOf L v + perpOf L w := by
  simp only [perpOf, map_add, ContinuousLinearMap.add_comp]
  abel

theorem perpOf_smul (r : ℝ) (v : V) : perpOf L (r • v) = r • perpOf L v := by
  simp only [perpOf, map_smul, ContinuousLinearMap.smul_comp, smul_sub]

/-- `perp_L` as a linear map. -/
def perpLin (L : E →L[ℝ] V) : V →ₗ[ℝ] V where
  toFun := perpOf L
  map_add' := perpOf_add
  map_smul' := perpOf_smul

/-- **Invariance of the mean curvature:** replacing `Q` by `Q(T·, T·) + L S(·, ·)` and `L` by `L ∘ T`
does not change `H`. -/
theorem HOf_comp_equiv (hL : Function.Injective L) (T : E ≃L[ℝ] E) (Q Q' : Hess V E)
    (S : E → E → E) (hQ : ∀ a c, Q' a c = Q (T a) (T c) + L (S a c)) :
    HOf (L.comp (T : E →L[ℝ] E)) Q' = HOf L Q := by
  have hB : ∀ a c, perpOf (L.comp (T : E →L[ℝ] E)) (Q' a c) =
      ((Q : E →L[ℝ] E →L[ℝ] V).toLinearMap₁₂.compr₂ (perpLin L)) (T a) (T c) := fun a c => by
    rw [perpOf_comp_equiv hL T, hQ, perpOf_add, perpOf_apply_self hL, add_zero]
    rfl
  have hB' : ∀ a c, perpOf L (Q a c) =
      ((Q : E →L[ℝ] E →L[ℝ] V).toLinearMap₁₂.compr₂ (perpLin L)) a c := fun a c => rfl
  simp only [HOf, hB, hB']
  exact trOf_comp_equiv hL T _

end CWFlow
