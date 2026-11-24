import CWFlow.Weinstein.Lift

/-!
# The Weinstein interface

Everything downstream of the Weinstein neighbourhood theorem (the scalar equation, Claim 3.3, uniqueness) uses
only the structure `IsLocalWeinstein Ω F U Ψ`, for a 2-form field `Ω` on the target. The target enters
nowhere else, so the almost Kähler case is obtained by supplying a different instance.

* `Ω` is a field of 2-forms on a normed space `V` (a chart of the ambient manifold). For `ℂⁿ` it is the constant
  form `omega0 J`, `ω₀(a, b) = ⟨Ja, b⟩` (`omega0_apply`).
* `IsLocalWeinstein Ω F U Ψ`: on an open `U ⊆ E × E*` (cotangent coordinates over a chart of `L`), `Ψ` is smooth,
  `Ψ(x, 0) = F x`, and `Ψ*Ω = ω_can`, `ω_can(a, b) = a₂ b₁ − b₂ a₁`.

The intended split of the existence proof:
(1) a target-specific first guess `Φ₀` with `Φ₀(x, 0) = F x` and a primitive `β` of `Φ₀*Ω − ω_can` vanishing on
    the zero section (for `ℂⁿ`: `Phi0`, `beta` with `A = liftA F`, `Weinstein/Coord.lean`; for almost Kähler: the
    exponential map and a fibrewise homotopy operator);
(2) the target-independent Moser step (`Moser/Theorem.lean`, `Weinstein/Existence.lean`), producing
    `Ψ = Φ₀ ∘ ψ₁`; the `ℂⁿ` instance is `exists_isLocalWeinstein_Cn` (`Weinstein/Cn.lean`).
-/

open ContinuousAlternatingMap

open scoped InnerProductSpace RealInnerProductSpace ContDiff

noncomputable section

namespace CWFlow

/-! ### Pullback of 2-forms, and the canonical form -/

section Pullback

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- `(Φ*Ω)(q) = Ω(Φ q) ∘ (DΦ(q) × DΦ(q))`. -/
def pullback2 (Ω : V → V [⋀^Fin 2]→L[ℝ] ℝ) (Φ : W → V) (q : W) : W [⋀^Fin 2]→L[ℝ] ℝ :=
  (Ω (Φ q)).compContinuousLinearMap (fderiv ℝ Φ q)

@[simp] theorem pullback2_apply (Ω : V → V [⋀^Fin 2]→L[ℝ] ℝ) (Φ : W → V) (q a b : W) :
    pullback2 Ω Φ q ![a, b] = Ω (Φ q) ![fderiv ℝ Φ q a, fderiv ℝ Φ q b] := by
  simp only [pullback2, compContinuousLinearMap_apply]
  congr 1
  ext i
  fin_cases i <;> rfl

end Pullback

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The canonical value `ω_can(a, b) = a₂ b₁ − b₂ a₁` on cotangent coordinates `E × E*`. -/
def omegaCanVal (a b : E × (E →L[ℝ] ℝ)) : ℝ := a.2 b.1 - b.2 a.1

/-! ### The interface -/

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- A local Weinstein map over a chart: `Ψ` is smooth on the open set `U` of cotangent coordinates, restricts to
`F` on the zero section, and pulls `Ω` back to the canonical form. -/
structure IsLocalWeinstein (Ω : V → V [⋀^Fin 2]→L[ℝ] ℝ) (F : E → V) (U : Set (E × (E →L[ℝ] ℝ)))
    (Ψ : E × (E →L[ℝ] ℝ) → V) : Prop where
  isOpen : IsOpen U
  contDiffOn : ContDiffOn ℝ ∞ Ψ U
  zero_section : ∀ x, (x, (0 : E →L[ℝ] ℝ)) ∈ U → Ψ (x, 0) = F x
  symplectic : ∀ q ∈ U, ∀ a b, pullback2 Ω Ψ q ![a, b] = omegaCanVal a b

/-! ### The `ℂⁿ` symplectic form -/

section Cn

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- `ω₀ = dλ₀`, `λ₀(z) = ½⟨Jz, ·⟩`, as a (constant) field of 2-forms. -/
def omega0 (J : V →L[ℝ] V) (z : V) : V [⋀^Fin 2]→L[ℝ] ℝ := extDeriv (toAlt1 (lam0L J)) z

theorem omega0_apply {J : V →L[ℝ] V} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v, J (J v) = -v) (z a b : V) : omega0 J z ![a, b] = ⟪J a, b⟫ := by
  have hd : DifferentiableAt ℝ (lam0L J) z := (lam0L J).differentiableAt
  rw [omega0, extDeriv_toAlt1 hd, dOne, (lam0L J).fderiv]
  simp only [lam0L, ContinuousLinearMap.coe_smul', ContinuousLinearMap.coe_comp', Pi.smul_apply,
    Function.comp_apply, innerSL_apply_apply, smul_eq_mul]
  rw [inner_J_swap hJi hJJ b a]
  ring

/-- For `ℂⁿ`, `Φ*ω₀(a, b) = ⟨J DΦ a, DΦ b⟩`. -/
theorem pullback2_omega0 {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {J : V →L[ℝ] V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (Φ : W → V) (q a b : W) :
    pullback2 (omega0 J) Φ q ![a, b] = ⟪J (fderiv ℝ Φ q a), fderiv ℝ Φ q b⟫ := by
  rw [pullback2_apply, omega0_apply hJi hJJ]

end Cn

end CWFlow
