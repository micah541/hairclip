import CWFlow.Moser.Forms
import CWFlow.Weinstein.Coord

/-!
# From linear 1-forms to Mathlib's alternating forms

`toAlt1 β q = ofSubsingleton (β q)` turns a field of linear functionals `β : W → W →L ℝ` into a 1-form
`W → W [⋀^Fin 1]→L ℝ`, with `toAlt1 β q ![b] = β q b` and `d(toAlt1 β)(a, b) = dOne β a b`. With this,
`dOne_beta` and `pb2_Phi0_zero` (stated with `dOne`) feed `moser_condition` (stated with `extDeriv`).
-/

open ContinuousAlternatingMap Filter

open scoped Topology

noncomputable section

namespace CWFlow

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- A field of linear functionals as a 1-form. -/
def toAlt1 (β : W → W →L[ℝ] ℝ) (q : W) : W [⋀^Fin 1]→L[ℝ] ℝ :=
  ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1) (β q)

@[simp] theorem toAlt1_apply (β : W → W →L[ℝ] ℝ) (q : W) (v : Fin 1 → W) :
    toAlt1 β q v = β q (v 0) := rfl

theorem toAlt1_eq (β : W → W →L[ℝ] ℝ) :
    toAlt1 β = fun q => (ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toContinuousLinearEquiv
      (β q) := rfl

theorem fderiv_toAlt1_apply {β : W → W →L[ℝ] ℝ} {p : W} (hβ : DifferentiableAt ℝ β p) (a b : W) :
    fderiv ℝ (toAlt1 β) p a ![b] = fderiv ℝ β p a b := by
  set L := (ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toContinuousLinearEquiv
  have h : HasFDerivAt (fun q => L (β q))
      ((L : (W →L[ℝ] ℝ) →L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ)).comp (fderiv ℝ β p)) p :=
    L.hasFDerivAt.comp p hβ.hasFDerivAt
  rw [toAlt1_eq, h.fderiv]
  rfl

theorem differentiableAt_toAlt1 {β : W → W →L[ℝ] ℝ} {p : W} (hβ : DifferentiableAt ℝ β p) :
    DifferentiableAt ℝ (toAlt1 β) p := by
  rw [toAlt1_eq]
  exact ((ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toContinuousLinearEquiv
    |>.hasFDerivAt).comp p hβ.hasFDerivAt |>.differentiableAt

theorem contDiffAt_toAlt1 {β : W → W →L[ℝ] ℝ} {p : W} {k : WithTop ℕ∞} (hβ : ContDiffAt ℝ k β p) :
    ContDiffAt ℝ k (toAlt1 β) p := by
  rw [toAlt1_eq]
  exact ((ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toContinuousLinearEquiv
    |>.contDiff).contDiffAt.comp p hβ

/-- `d(toAlt1 β)(a, b) = dOne β a b`. -/
theorem extDeriv_toAlt1 {β : W → W →L[ℝ] ℝ} {p : W} (hβ : DifferentiableAt ℝ β p) (a b : W) :
    extDeriv (toAlt1 β) p ![a, b] = dOne β p a b := by
  rw [extDeriv_one_apply (differentiableAt_toAlt1 hβ), fderiv_toAlt1_apply hβ,
    fderiv_toAlt1_apply hβ, dOne]

end CWFlow
