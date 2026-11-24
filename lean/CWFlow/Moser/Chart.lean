import Mathlib

/-!
# Moser's argument in a chart

For a field of bilinear forms `Ω : E → E →L E →L ℝ` (a 2-form in a chart) and a vector field `X`:

* `liePair Ω X p a b = DΩ(p)[X p](a, b) + Ω p (DX a) b + Ω p a (DX b)` is the Lie derivative `L_X Ω`;
* **Cartan for closed 2-forms**: if `Ω` is antisymmetric and closed (`DΩ[a](b,c)` has vanishing cyclic sum),
  then `L_X Ω (a, b) = d(ι_X Ω)(a, b) = D(ι_X Ω)[a](b) − D(ι_X Ω)[b](a)`;
* **the pullback along a flow**: if `ψ_s` is the flow of `X_s` and `Dψ_s` solves the variational equation,
  then `d/ds (ψ_s* Ω_s)(v, w) = (Ω̇_s + L_{X_s} Ω_s)(Dψ_s v, Dψ_s w)`.

So if `ι_{X_s} Ω_s = −β` with `dβ = Ω̇_s`, the pullback `ψ_s* Ω_s` is constant in `s`.
-/

open scoped Topology

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The Lie derivative `L_X Ω` of a field of bilinear forms, at `p`. -/
def liePair (Ω : E → E →L[ℝ] E →L[ℝ] ℝ) (X : E → E) (p a b : E) : ℝ :=
  fderiv ℝ Ω p (X p) a b + Ω p (fderiv ℝ X p a) b + Ω p a (fderiv ℝ X p b)

/-- The exterior derivative of a 1-form field `β` in a chart: `dβ(a, b) = Dβ[a](b) − Dβ[b](a)`. -/
def dOne (β : E → E →L[ℝ] ℝ) (p a b : E) : ℝ := fderiv ℝ β p a b - fderiv ℝ β p b a

/-- The derivative of `q ↦ Ω q (X q)` (the 1-form `ι_X Ω`). -/
theorem fderiv_iota_apply {Ω : E → E →L[ℝ] E →L[ℝ] ℝ} {X : E → E} {p : E}
    (hΩ : DifferentiableAt ℝ Ω p) (hX : DifferentiableAt ℝ X p) (a b : E) :
    fderiv ℝ (fun q => Ω q (X q)) p a b = Ω p (fderiv ℝ X p a) b + fderiv ℝ Ω p a (X p) b := by
  rw [(hΩ.hasFDerivAt.clm_apply hX.hasFDerivAt).fderiv]
  rfl

/-- **Cartan's formula for closed 2-forms**: `L_X Ω = d(ι_X Ω)`. -/
theorem liePair_eq_dOne {Ω : E → E →L[ℝ] E →L[ℝ] ℝ} {X : E → E} {p : E}
    (hΩ : DifferentiableAt ℝ Ω p) (hX : DifferentiableAt ℝ X p)
    (hanti : ∀ a b, Ω p a b = -Ω p b a)
    (hDanti : ∀ c a b, fderiv ℝ Ω p c a b = -fderiv ℝ Ω p c b a)
    (hclosed : ∀ a b c, fderiv ℝ Ω p a b c + fderiv ℝ Ω p b c a + fderiv ℝ Ω p c a b = 0)
    (a b : E) :
    liePair Ω X p a b = dOne (fun q => Ω q (X q)) p a b := by
  rw [liePair, dOne, fderiv_iota_apply hΩ hX, fderiv_iota_apply hΩ hX]
  have h1 := hclosed a b (X p)
  have h2 := hDanti b (X p) a
  have h3 := hDanti a b (X p)
  have h4 := hanti (fderiv ℝ X p b) a
  linarith

/-- The derivative of an antisymmetric field of forms is antisymmetric. -/
theorem fderiv_anti {Ω : E → E →L[ℝ] E →L[ℝ] ℝ} {p : E} (hΩ : DifferentiableAt ℝ Ω p)
    (hanti : ∀ q a b, Ω q a b = -Ω q b a) (c a b : E) :
    fderiv ℝ Ω p c a b = -fderiv ℝ Ω p c b a := by
  have hab : (fun q => Ω q a b) = fun q => -Ω q b a := funext fun q => hanti q a b
  have e1 : fderiv ℝ (fun q => Ω q a b) p c = fderiv ℝ Ω p c a b := by
    rw [((hΩ.hasFDerivAt.clm_apply (hasFDerivAt_const a p)).clm_apply
      (hasFDerivAt_const b p)).fderiv]
    simp
  have e2 : fderiv ℝ (fun q => Ω q b a) p c = fderiv ℝ Ω p c b a := by
    rw [((hΩ.hasFDerivAt.clm_apply (hasFDerivAt_const b p)).clm_apply
      (hasFDerivAt_const a p)).fderiv]
    simp
  rw [← e1, ← e2, hab, show (fun q => -Ω q b a) = -(fun q => Ω q b a) from rfl, fderiv_neg]
  rfl

/-- **The pullback along a flow.** Let `ψ_s z` solve `∂ₛψ = X_s(ψ)`, let `Dψ_s` solve the variational
equation `∂ₛDψ = DX_s(ψ) ∘ Dψ`, and let `s ↦ Ω_s(ψ_s z)` have derivative `Ω̇ + DΩ[X]`. Then
`d/ds Ω_s(ψ_s z)(Dψ_s v, Dψ_s w) = (Ω̇ + L_X Ω)(Dψ v, Dψ w)`. -/
theorem hasDerivAt_pullback {Ω Ωdot : ℝ → E → E →L[ℝ] E →L[ℝ] ℝ} {X : ℝ → E → E}
    {ψ : ℝ → E} {Dψ : ℝ → E →L[ℝ] E} {s : ℝ}
    (hcurve : HasDerivAt (fun r => Ω r (ψ r)) (Ωdot s (ψ s) + fderiv ℝ (Ω s) (ψ s) (X s (ψ s))) s)
    (hDψ : HasDerivAt Dψ ((fderiv ℝ (X s) (ψ s)).comp (Dψ s)) s) (v w : E) :
    HasDerivAt (fun r => Ω r (ψ r) (Dψ r v) (Dψ r w))
      (Ωdot s (ψ s) (Dψ s v) (Dψ s w) + liePair (Ω s) (X s) (ψ s) (Dψ s v) (Dψ s w)) s := by
  have hv : HasDerivAt (fun r => Dψ r v) (fderiv ℝ (X s) (ψ s) (Dψ s v)) s := by
    have := hDψ.clm_apply (hasDerivAt_const s v)
    simpa using this
  have hw : HasDerivAt (fun r => Dψ r w) (fderiv ℝ (X s) (ψ s) (Dψ s w)) s := by
    have := hDψ.clm_apply (hasDerivAt_const s w)
    simpa using this
  have h := (hcurve.clm_apply hv).clm_apply hw
  refine h.congr_deriv ?_
  rw [liePair]
  simp only [ContinuousLinearMap.coe_add', Pi.add_apply]
  ring

/-- **Moser.** If along the flow the pulled-back derivative vanishes (`Ω̇ + L_X Ω = 0` at `ψ_s z` for all
`s ∈ [0, 1]`), then `ψ_s* Ω_s` is constant on `[0, 1]`. -/
theorem pullback_const {Ω Ωdot : ℝ → E → E →L[ℝ] E →L[ℝ] ℝ} {X : ℝ → E → E}
    {ψ : ℝ → E} {Dψ : ℝ → E →L[ℝ] E}
    (hcurve : ∀ s, HasDerivAt (fun r => Ω r (ψ r))
      (Ωdot s (ψ s) + fderiv ℝ (Ω s) (ψ s) (X s (ψ s))) s)
    (hDψ : ∀ s, HasDerivAt Dψ ((fderiv ℝ (X s) (ψ s)).comp (Dψ s)) s)
    (hmoser : ∀ s a b, Ωdot s (ψ s) a b + liePair (Ω s) (X s) (ψ s) a b = 0) (v w : E) (s : ℝ) :
    Ω s (ψ s) (Dψ s v) (Dψ s w) = Ω 0 (ψ 0) (Dψ 0 v) (Dψ 0 w) := by
  have h : ∀ r, HasDerivAt (fun r => Ω r (ψ r) (Dψ r v) (Dψ r w)) 0 r := fun r => by
    have := hasDerivAt_pullback (hcurve r) (hDψ r) v w
    rwa [hmoser] at this
  exact is_const_of_deriv_eq_zero (fun r => (h r).differentiableAt) (fun r => (h r).deriv) s 0

end CWFlow
