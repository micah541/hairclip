import Mathlib.Analysis.Calculus.ContDiff.CPolynomial
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Analysis.Calculus.ContDiff.FaaDiBruno
import Mathlib.Analysis.Normed.Module.Multilinear.Curry

/-!
# Joint smoothness of multilinear evaluation

`contDiff_cml_apply`: `(f, m) ↦ f m` is smooth, jointly in the continuous multilinear map `f` and the
arguments `m` (by induction on the number of arguments, currying the first one). Consequently the Taylor
composition `taylorComp` used in Faà di Bruno's formula is smooth in the coefficients that it uses
(`contDiff_compAlongOrderedFinpartition`).
-/

open scoped ContDiff

noncomputable section

namespace CWFlow

universe u v

theorem contDiff_cml_apply : ∀ (n : ℕ) {Ei : Fin n → Type u} [∀ i, NormedAddCommGroup (Ei i)]
    [∀ i, NormedSpace ℝ (Ei i)] {G : Type v} [NormedAddCommGroup G] [NormedSpace ℝ G],
    ContDiff ℝ ∞ (fun q : ContinuousMultilinearMap ℝ Ei G × (∀ i, Ei i) => q.1 q.2) := by
  intro n
  induction n with
  | zero =>
    intro Ei _ _ G _ _
    have h : (fun q : ContinuousMultilinearMap ℝ Ei G × (∀ i, Ei i) => q.1 q.2) =
        fun q => ContinuousMultilinearMap.apply ℝ Ei G 0 q.1 := by
      funext q
      simp only [ContinuousMultilinearMap.apply_apply]
      congr 1
      funext i
      exact i.elim0
    rw [h]
    exact (ContinuousMultilinearMap.apply ℝ Ei G 0).contDiff.comp contDiff_fst
  | succ n ih =>
    intro Ei _ _ G _ _
    have hL : ContDiff ℝ ∞ (fun q : ContinuousMultilinearMap ℝ Ei G × (∀ i, Ei i) =>
        continuousMultilinearCurryLeftEquiv ℝ Ei G q.1 (q.2 0)) :=
      ((continuousMultilinearCurryLeftEquiv ℝ Ei G).contDiff.comp contDiff_fst).clm_apply
        ((ContinuousLinearMap.proj (R := ℝ) (φ := Ei) 0).contDiff.comp contDiff_snd)
    have hT : ContDiff ℝ ∞ (fun q : ContinuousMultilinearMap ℝ Ei G × (∀ i, Ei i) =>
        (fun i : Fin n => q.2 i.succ)) :=
      contDiff_pi.2 fun i => (ContinuousLinearMap.proj (R := ℝ) (φ := Ei) i.succ).contDiff.comp
        contDiff_snd
    have h := (ih (Ei := fun i : Fin n => Ei i.succ) (G := G)).comp (hL.prodMk hT)
    convert h using 1
    funext q
    simp only [Function.comp_apply, continuousMultilinearCurryLeftEquiv_apply]
    congr 1
    exact (Fin.cons_self_tail q.2).symm

end CWFlow
