import CWFlow.Weinstein.Interface
import CWFlow.Moser.Theorem

/-!
# The Weinstein step, independent of the target

On cotangent coordinates `E × E*`, `λ_can(x, ξ)(a) = ξ a₁` (`lamCan`; defined as a function rather than as a
continuous linear map of `(x, ξ)`, since `fderiv` of a coerced `→L`-valued `→L` hits the instance mismatch noted in
`Moser/Forms.lean`), and `dλ_can = ω_can` is nondegenerate
(`extDeriv_lamCan`, `curryLeft_lamCan_injective`).

`weinsteinMap Φ₀ β O = Φ₀ ∘ moserEnd λ_can β O` is the canonical candidate, on `weinsteinDom β O = moserDom λ_can β O`.
`isLocalWeinstein_weinsteinMap`: let `Ω` be any field of 2-forms on `V`.
Suppose:
* `Φ₀` and `β` are smooth on an open `O ⊆ E × E*`;
* `Φ₀(x, 0) = F x` and `β(x, 0) = 0` on the zero section;
* `Φ₀*Ω = ω_can + dβ` on `O`, and `Φ₀*Ω = ω_can` on the zero section.

Then `weinsteinMap` is a local Weinstein map on `weinsteinDom`, which contains the zero section over `O`. No
compactness is needed. The target enters only through the pair `(Φ₀, β)`: for `ℂⁿ`, see `Weinstein/Cn.lean`.
-/

open ContinuousAlternatingMap Filter Set

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The Liouville form `λ_can(x, ξ)(a) = ξ a₁`. -/
def lamCan (q : E × (E →L[ℝ] ℝ)) : E × (E →L[ℝ] ℝ) →L[ℝ] ℝ :=
  q.2.comp (ContinuousLinearMap.fst ℝ E (E →L[ℝ] ℝ))

@[simp] theorem lamCan_apply (q a : E × (E →L[ℝ] ℝ)) : lamCan q a = q.2 a.1 := rfl

theorem contDiff_lamCan : ContDiff ℝ ∞ (lamCan (E := E)) :=
  contDiff_snd.clm_comp contDiff_const

theorem fderiv_lamCan_apply (q a b : E × (E →L[ℝ] ℝ)) : fderiv ℝ lamCan q a b = a.2 b.1 := by
  have h := (hasFDerivAt_snd (p := q)).clm_comp
    (hasFDerivAt_const (ContinuousLinearMap.fst ℝ E (E →L[ℝ] ℝ)) q)
  have e : (lamCan (E := E)) = fun q => q.2.comp (ContinuousLinearMap.fst ℝ E (E →L[ℝ] ℝ)) := rfl
  rw [e, h.fderiv]
  simp

theorem dOne_lamCan (q a b : E × (E →L[ℝ] ℝ)) : dOne lamCan q a b = omegaCanVal a b := by
  rw [dOne, fderiv_lamCan_apply, fderiv_lamCan_apply]
  rfl

/-- `dλ_can = ω_can`. -/
theorem extDeriv_lamCan (q a b : E × (E →L[ℝ] ℝ)) :
    extDeriv (toAlt1 (lamCan (E := E))) q ![a, b] = omegaCanVal a b := by
  rw [extDeriv_toAlt1 (contDiff_lamCan.differentiable (by simp) q), dOne_lamCan]

/-- `ω_can` is nondegenerate. -/
theorem curryLeft_lamCan_injective (q : E × (E →L[ℝ] ℝ)) :
    Function.Injective (extDeriv (toAlt1 (lamCan (E := E))) q).curryLeft := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h : ∀ b, omegaCanVal v b = 0 := fun b => by
    have h := congrArg (fun f => f ![b]) hv
    simp only [curryLeft_apply_apply] at h
    have e : (Matrix.vecCons v ![b] : Fin 2 → E × (E →L[ℝ] ℝ)) = ![v, b] := rfl
    rw [← extDeriv_lamCan q]
    simpa [e] using h
  have h2 : v.2 = 0 := by
    ext u
    simpa [omegaCanVal] using h (u, 0)
  have h1 : v.1 = 0 := by
    refine SeparatingDual.eq_zero_of_forall_dual_eq_zero (R := ℝ) fun f => ?_
    simpa [omegaCanVal, h2] using h (0, f)
  exact Prod.ext h1 h2

variable [FiniteDimensional ℝ E]
variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- The canonical Weinstein map built from `(Φ₀, β)`. -/
def weinsteinMap (Φ₀ : E × (E →L[ℝ] ℝ) → V) (β : E × (E →L[ℝ] ℝ) → (E × (E →L[ℝ] ℝ)) →L[ℝ] ℝ)
    (O : Set (E × (E →L[ℝ] ℝ))) : E × (E →L[ℝ] ℝ) → V :=
  Φ₀ ∘ moserEnd lamCan β O

/-- **The Weinstein step for an arbitrary target form `Ω`**, given a first guess `Φ₀` and a primitive `β` of
`Φ₀*Ω − ω_can` vanishing on the zero section. -/
theorem isLocalWeinstein_weinsteinMap {Ω : V → V [⋀^Fin 2]→L[ℝ] ℝ} {F : E → V}
    {Φ₀ : E × (E →L[ℝ] ℝ) → V} {β : E × (E →L[ℝ] ℝ) → (E × (E →L[ℝ] ℝ)) →L[ℝ] ℝ}
    {O : Set (E × (E →L[ℝ] ℝ))} (hO : IsOpen O) (hΦ₀ : ContDiffOn ℝ ∞ Φ₀ O)
    (hβ : ContDiffOn ℝ ∞ β O) (hΦF : ∀ x, (x, 0) ∈ O → Φ₀ (x, 0) = F x)
    (hβ0 : ∀ x, (x, 0) ∈ O → β (x, 0) = 0)
    (hprim : ∀ q ∈ O, ∀ a b, pullback2 Ω Φ₀ q ![a, b] = omegaCanVal a b + dOne β q a b)
    (hsymp : ∀ x, (x, 0) ∈ O → ∀ a b, pullback2 Ω Φ₀ (x, 0) ![a, b] = omegaCanVal a b) :
    IsLocalWeinstein Ω F (moserDom lamCan β O) (weinsteinMap Φ₀ β O) ∧
      ∀ x, (x, 0) ∈ O → (x, (0 : E →L[ℝ] ℝ)) ∈ moserDom lamCan β O := by
  have hlam : ContDiffOn ℝ ∞ (lamCan (E := E)) O := contDiff_lamCan.contDiffOn
  have hβd : ∀ y ∈ O, DifferentiableAt ℝ β y := fun y hy =>
    (hβ.contDiffAt (hO.mem_nhds hy)).differentiableAt (by simp)
  have hdβ0 : ∀ x, (x, 0) ∈ O → extDeriv (toAlt1 β) (x, 0) = 0 := by
    intro x hx
    ext v
    have hv : v = ![v 0, v 1] := by funext i; fin_cases i <;> rfl
    have h := hprim _ hx (v 0) (v 1)
    rw [hsymp x hx] at h
    have h0 : dOne β (x, 0) (v 0) (v 1) = 0 := by linarith
    rw [hv, extDeriv_toAlt1 (hβd _ hx), h0]
    simp
  have hzero : ∀ x, (x, 0) ∈ O → (x, (0 : E →L[ℝ] ℝ)) ∈ moserDom lamCan β O := fun x hx =>
    mem_moserDom_of_zero hO hlam hβ (hβ0 x hx) fun t _ =>
      ⟨hx, iotaT_invertible t (hdβ0 x hx) (curryLeft_lamCan_injective _)⟩
  have hΩ1 : ∀ y ∈ O, ∀ c d, omegaT lamCan β 1 y ![c, d] = omegaCanVal c d + dOne β y c d := by
    intro y hy c d
    simp only [omegaT, ContinuousAlternatingMap.add_apply, ContinuousAlternatingMap.smul_apply,
      one_smul]
    rw [extDeriv_lamCan, extDeriv_toAlt1 (hβd y hy)]
  have hD := isOpen_moserDom hO hlam hβ
  have hE := contDiffOn_moserEnd hO hlam hβ
  have hmaps : MapsTo (moserEnd lamCan β O) (moserDom lamCan β O) O := fun q hq =>
    moserEnd_mem hO hlam hβ hq
  refine ⟨⟨hD, hΦ₀.comp hE hmaps, fun x hx => ?_, fun q hq a b => ?_⟩, hzero⟩
  · have hxO := moserDom_subset hx
    rw [weinsteinMap, Function.comp_apply, moserEnd_of_zero hO hlam hβ (hβ0 x hxO) hx, hΦF x hxO]
  · have hΦd : DifferentiableAt ℝ Φ₀ (moserEnd lamCan β O q) :=
      (hΦ₀.contDiffAt (hO.mem_nhds (hmaps hq))).differentiableAt (by simp)
    have hψd : DifferentiableAt ℝ (moserEnd lamCan β O) q :=
      ((hE.contDiffAt (hD.mem_nhds hq)).differentiableAt (by simp))
    rw [weinsteinMap, pullback2_apply, fderiv_comp q hΦd hψd]
    simp only [ContinuousLinearMap.comp_apply, Function.comp_apply]
    rw [← pullback2_apply, hprim _ (hmaps hq), ← hΩ1 _ (hmaps hq),
      moserEnd_symplectic hO hlam hβ hq, extDeriv_lamCan]

end CWFlow
