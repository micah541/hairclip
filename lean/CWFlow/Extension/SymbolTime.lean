import CWFlow.Extension.SymbolDeriv
import CWFlow.Extension.JetTime

/-!
# Time derivative of a symbol along jets

For `Φ` smooth on an open set `Ω ⊆ ℝ × ℝⁿ` and a symbol `𝔭`, the function `s ↦ 𝔭 (jp N Φ(s,·) x)` is
differentiable with derivative `D𝔭 (0, jets of ∂ₜΦ)` (`hasDerivAt_symbol_slice`), and the derivative of a
symbol in the pure jet directions is bounded slot by slot (`IsSymbol.fderiv_jets_le`).
-/

open scoped ContDiff
open CWFlow.Interp

namespace CWFlow.Extension

variable {n N : ℕ}

/-- The time derivative of a symbol along the jets of the slices of `Φ`. -/
theorem hasDerivAt_symbol_slice {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {Ω : Set (ℝ × EuclideanSpace ℝ (Fin n))}
    (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω) {t : ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (htx : (t, x) ∈ Ω) (hp : jp N (fun y => Φ (t, y)) x ∈ A) :
    HasDerivAt (fun s : ℝ => 𝔭 (jp N (fun y => Φ (s, y)) x))
      (fderiv ℝ 𝔭 (jp N (fun y => Φ (t, y)) x)
        ((0 : EuclideanSpace ℝ (Fin n)), fun m : Fin (N + 1) =>
          CWFlow.Interp.jet (m : ℕ) (fun y => dt Φ (t, y)) x)) t := by
  have hcurve : HasDerivAt (fun s : ℝ => jp N (fun y => Φ (s, y)) x)
      ((0 : EuclideanSpace ℝ (Fin n)), fun m : Fin (N + 1) =>
        CWFlow.Interp.jet (m : ℕ) (fun y => dt Φ (t, y)) x) t := by
    refine (hasDerivAt_const t x).prodMk ?_
    exact hasDerivAt_pi.2 fun m => hasDerivAt_jet_slice hΩ hΦ (m : ℕ) htx
  have hd : DifferentiableAt ℝ 𝔭 (jp N (fun y => Φ (t, y)) x) :=
    (h.smooth.contDiffAt (h.isOpen.mem_nhds hp)).differentiableAt (by simp)
  exact hd.hasFDerivAt.comp_hasDerivAt t hcurve

/-- The derivative of a symbol in the jet directions, slot by slot. -/
theorem IsSymbol.fderiv_jets_le {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p ∈ A, ∀ ζ' : J n N,
      ‖fderiv ℝ 𝔭 p ((0 : EuclideanSpace ℝ (Fin n)), ζ')‖ ≤
        C * ∑ m : Fin (N + 1), OmZ ((w : ℤ) - wt (some m : Slot N)) ℓ (aOf p) * ‖ζ' m‖ := by
  obtain ⟨C, hC⟩ := h.fderiv_le
  refine ⟨max C 0, le_max_right _ _, fun p hp ζ' => ?_⟩
  refine (hC p hp ((0 : EuclideanSpace ℝ (Fin n)), ζ')).trans ?_
  rw [Fintype.sum_option]
  have h0 : comp (none : Slot N) ((0 : EuclideanSpace ℝ (Fin n)), ζ') = 0 := by simp [comp]
  have hs : 0 ≤ ∑ m : Fin (N + 1), OmZ ((w : ℤ) - wt (some m : Slot N)) ℓ (aOf p) * ‖ζ' m‖ :=
    Finset.sum_nonneg fun m _ => mul_nonneg (OmZ_nonneg (aOf_nonneg p)) (norm_nonneg _)
  simp only [h0, mul_zero, zero_add]
  exact mul_le_mul_of_nonneg_right (le_max_left _ _) hs

/-- info: 'CWFlow.Extension.hasDerivAt_symbol_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_symbol_slice

/-- info: 'CWFlow.Extension.IsSymbol.fderiv_jets_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.fderiv_jets_le

end CWFlow.Extension
