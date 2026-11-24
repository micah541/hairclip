import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.Complex.Circle
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.InnerProductSpace.PiL2

/-!
# Curves in `ℂ`: basic definitions (Mathlib only)

For a parametrised plane curve `c : ℝ → ℂ`:

* `curv c θ = ⟪c'', i c'⟫ / |c'|³`, the signed curvature (the same definition as in the `n = 1` instance
  `CWFlow.Curve.flowVel_eq_neg_kss`);
* `dS c f θ = f' / |c'|`, the arclength derivative of a scalar function;
* `lift F θ = F (e^{iθ})`, the `2π`-periodic lift of a map on the circle;
* `bendingEnergy c = ∫₀^{2π} κ² |c'| dθ`, i.e. `∫ κ² ds` over one period.
-/

open scoped RealInnerProductSpace

noncomputable section

namespace CWFlow.Curve

/-- The model space of a 1-manifold: Mathlib charts `Circle` on `EuclideanSpace ℝ (Fin 1)`
(`𝓡 1` is notation for `𝓘(ℝ, EuclideanSpace ℝ (Fin 1))`). -/
abbrev E1 := EuclideanSpace ℝ (Fin 1)

/-- Signed curvature `κ = ⟪c'', i c'⟫ / |c'|³` of a plane curve. -/
def curv (c : ℝ → ℂ) (θ : ℝ) : ℝ := ⟪deriv (deriv c) θ, Complex.I * deriv c θ⟫ / ‖deriv c θ‖ ^ 3

/-- Arclength derivative `∂_s f = f' / |c'|`. -/
def dS (c : ℝ → ℂ) (f : ℝ → ℝ) (θ : ℝ) : ℝ := deriv f θ / ‖deriv c θ‖

/-- The `2π`-periodic lift of a map on the circle. -/
def lift (F : Circle → ℂ) (θ : ℝ) : ℂ := F (Circle.exp θ)

/-- `∫ κ² ds` over one period, `ds = |c'| dθ`. -/
def bendingEnergy (c : ℝ → ℂ) : ℝ := ∫ θ in (0 : ℝ)..(2 * Real.pi), curv c θ ^ 2 * ‖deriv c θ‖

end CWFlow.Curve
