import CWFlow.Extension.JetSymm
import CWFlow.Interp.HamiltonStep
import Mathlib.Analysis.Calculus.LineDeriv.IntegrationByParts

/-!
# Integration by parts for coordinate jets

For `φ` smooth and `w` compactly supported of class `C¹` (resp. `C²`), the coordinate jets satisfy
`∫ (∂ₐ F) w = -∫ F ∂ₐ w` for `F = jet m φ · I`; combined with the symmetry of jets
(`jet_comp_perm`) this moves any index of a jet coordinate onto the weight.
-/

open scoped ContDiff
open CWFlow.Interp MeasureTheory

namespace CWFlow.Extension

variable {n : ℕ}

private theorem jibp_ibp {F w : EuclideanSpace ℝ (Fin n) → ℝ} (hF : ContDiff ℝ 1 F)
    (hw : ContDiff ℝ 1 w) (hc : HasCompactSupport w) (a : Fin n) :
    ∫ x, fderiv ℝ F x (EuclideanSpace.single a 1) * w x =
      -∫ x, F x * fderiv ℝ w x (EuclideanSpace.single a 1) := by
  have h1 : Integrable (fun x => fderiv ℝ F x (EuclideanSpace.single a 1) * w x) :=
    (((hF.continuous_fderiv (by norm_num)).clm_apply continuous_const).mul
      hw.continuous).integrable_of_hasCompactSupport hc.mul_left
  have h2 : Integrable (fun x => F x * fderiv ℝ w x (EuclideanSpace.single a 1)) :=
    (hF.continuous.mul ((hw.continuous_fderiv (by norm_num)).clm_apply
      continuous_const)).integrable_of_hasCompactSupport
      ((hc.fderiv_apply ℝ (EuclideanSpace.single a (1 : ℝ))).mul_left)
  have h3 : Integrable (fun x => F x * w x) :=
    (hF.continuous.mul hw.continuous).integrable_of_hasCompactSupport hc.mul_left
  have := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume)
    (f := F) (g := w) (v := EuclideanSpace.single a (1 : ℝ)) h1 h2 h3
    (fun x _ => hF.differentiable (by norm_num) x) (fun x _ => hw.differentiable (by norm_num) x)
  linarith

private theorem jibp_coord_contDiff {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (m : ℕ) (I : Fin m → Fin n) : ContDiff ℝ ∞ (fun y => jet m φ y I) :=
  (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin m → Fin n => ℝ) I).contDiff.comp (contDiff_jet hφ m)

/-- One integration by parts in the first index. -/
theorem integral_jet_cons_mul {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (a : Fin n) (I : Fin m → Fin n) {w : EuclideanSpace ℝ (Fin n) → ℝ} (hw : ContDiff ℝ 1 w)
    (hc : HasCompactSupport w) :
    ∫ x, jet (m + 1) φ x (Fin.cons a I) * w x =
      -∫ x, jet m φ x I * fderiv ℝ w x (EuclideanSpace.single a 1) := by
  have := jibp_ibp ((jibp_coord_contDiff hφ m I).of_le (by simp)) hw hc a
  simp_rw [jet_cons_eq_fderiv hφ]
  exact this

/-- Two integrations by parts in the first two indices. -/
theorem integral_jet_cons_cons_mul {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (a c : Fin n) (I : Fin m → Fin n) {w : EuclideanSpace ℝ (Fin n) → ℝ} (hw : ContDiff ℝ 2 w)
    (hc : HasCompactSupport w) :
    ∫ x, jet (m + 2) φ x (Fin.cons a (Fin.cons c I)) * w x =
      ∫ x, jet m φ x I *
        fderiv ℝ (fun y => fderiv ℝ w y (EuclideanSpace.single c 1)) x
          (EuclideanSpace.single a 1) := by
  have hw1 : ContDiff ℝ 1 w := hw.of_le (by norm_num)
  have hg : ContDiff ℝ 1 (fun y => fderiv ℝ w y (EuclideanSpace.single c (1 : ℝ))) :=
    ((hw.fderiv_right (m := 1) (by norm_num)).clm_apply contDiff_const)
  have hgc : HasCompactSupport (fun y => fderiv ℝ w y (EuclideanSpace.single c (1 : ℝ))) :=
    hc.fderiv_apply ℝ _
  have e : ∀ x, jet (m + 2) φ x (Fin.cons a (Fin.cons c I)) =
      jet (m + 2) φ x (Fin.cons c (Fin.cons a I)) := fun x => jet_swap_first hφ m x a c I
  simp_rw [e]
  rw [integral_jet_cons_mul hφ (m + 1) c (Fin.cons a I) hw1 hc,
    integral_jet_cons_mul hφ m a I hg hgc, neg_neg]

/-- Moving an index to the front does not change a jet coordinate. -/
theorem jet_insertNth_eq_cons {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (i : Fin (m + 1)) (a : Fin n) (I : Fin m → Fin n) :
    jet (m + 1) φ x (Fin.insertNth i a I) = jet (m + 1) φ x (Fin.cons a I) := by
  rw [← Fin.cons_comp_cycleRange a I i]
  exact jet_comp_perm hφ (m + 1) x i.cycleRange _

/-- Appending instead of prepending: the last index can be moved to the front. -/
theorem jet_snoc_eq_cons {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (a : Fin n) (I : Fin m → Fin n) :
    jet (m + 1) φ x (Fin.snoc I a) = jet (m + 1) φ x (Fin.cons a I) := by
  rw [← jet_insertNth_eq_cons hφ m x (Fin.last m) a I, Fin.insertNth_last']

/-- One integration by parts in an arbitrary index position. -/
theorem integral_jet_insertNth_mul {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (m : ℕ) (i : Fin (m + 1)) (a : Fin n) (I : Fin m → Fin n)
    {w : EuclideanSpace ℝ (Fin n) → ℝ} (hw : ContDiff ℝ 1 w) (hc : HasCompactSupport w) :
    ∫ x, jet (m + 1) φ x (Fin.insertNth i a I) * w x =
      -∫ x, jet m φ x I * fderiv ℝ w x (EuclideanSpace.single a 1) := by
  simp_rw [jet_insertNth_eq_cons hφ]
  exact integral_jet_cons_mul hφ m a I hw hc

/-- info: 'CWFlow.Extension.integral_jet_cons_mul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_jet_cons_mul

/-- info: 'CWFlow.Extension.integral_jet_cons_cons_mul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_jet_cons_cons_mul

/-- info: 'CWFlow.Extension.integral_jet_insertNth_mul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_jet_insertNth_mul

end CWFlow.Extension
