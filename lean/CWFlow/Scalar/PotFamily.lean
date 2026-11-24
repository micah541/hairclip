import CWFlow.Scalar.Smooth
import CWFlow.Scalar.GraphTime

/-!
# The potential `div(JH)` of a jointly smooth family

`contDiffOn_flowPot_family`: for `w(s, z)` jointly smooth on `[0, T] × U` with immersed slices,
`(s, y) ↦ flowPot J (w(s, ·)) y` is jointly smooth on `[0, T] × U`.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

theorem contDiffOn_flowPot_family (J : V →L[ℝ] V) {T : ℝ} (hT : 0 < T) {U : Set E} (hU : IsOpen U)
    {w : ℝ × E → V} (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U))
    (hinj : ∀ p ∈ Icc (0 : ℝ) T ×ˢ U, Function.Injective (fderiv ℝ (fun z => w (p.1, z)) p.2)) :
    ContDiffOn ℝ ∞ (fun p : ℝ × E => flowPot J (fun z => w (p.1, z)) p.2) (Icc 0 T ×ˢ U) := by
  have h1 : ContDiffOn ℝ ∞ (fun p : ℝ × E => fderiv ℝ (fun z => w (p.1, z)) p.2) (Icc 0 T ×ˢ U) :=
    contDiffOn_slice_fderiv hT hU hw
  have h2 : ContDiffOn ℝ ∞ (fun p : ℝ × E => fderiv ℝ (fderiv ℝ (fun z => w (p.1, z))) p.2)
      (Icc 0 T ×ˢ U) := contDiffOn_slice_fderiv hT hU h1
  have hm : ContDiffOn ℝ ∞ (fun p : ℝ × E => J (meanCurv (fun z => w (p.1, z)) p.2))
      (Icc 0 T ×ˢ U) := by
    intro p hp
    have hH := contDiffAt_HOfP (V := V) (E := E)
      (p := (fderiv ℝ (fun z => w (p.1, z)) p.2, fderiv ℝ (fderiv ℝ (fun z => w (p.1, z))) p.2))
      (hinj p hp)
    have hc := ContDiffAt.comp_contDiffWithinAt (g := HOfP)
      (f := fun q : ℝ × E => (fderiv ℝ (fun z => w (q.1, z)) q.2,
        fderiv ℝ (fderiv ℝ (fun z => w (q.1, z))) q.2)) p hH ((h1 p hp).prodMk (h2 p hp))
    exact J.contDiff.contDiffAt.comp_contDiffWithinAt p hc
  have h3 : ContDiffOn ℝ ∞ (fun p : ℝ × E =>
      fderiv ℝ (fun z => J (meanCurv (fun z' => w (p.1, z')) z)) p.2) (Icc 0 T ×ˢ U) :=
    contDiffOn_slice_fderiv hT hU hm
  intro p hp
  have hG : ContDiffWithinAt ℝ ∞ (fun q : ℝ × E => ginvOf (fderiv ℝ (fun z => w (q.1, z)) q.2))
      (Icc 0 T ×ˢ U) p :=
    ContDiffAt.comp_contDiffWithinAt (g := (ginvOf : (E →L[ℝ] V) → (E →L[ℝ] ℝ) →L[ℝ] E))
      (f := fun q : ℝ × E => fderiv ℝ (fun z => w (q.1, z)) q.2) p
      (contDiffAt_ginvOf (hinj p hp)) (h1 p hp)
  have hterm : ∀ i, ContDiffWithinAt ℝ ∞ (fun q : ℝ × E =>
      ⟪fderiv ℝ (fun z => J (meanCurv (fun z' => w (q.1, z')) z)) q.2 (bas E i),
        fderiv ℝ (fun z => w (q.1, z)) q.2
          (ginvOf (fderiv ℝ (fun z => w (q.1, z)) q.2) (dualBas i))⟫) (Icc 0 T ×ˢ U) p := fun i =>
    ((h3 p hp).clm_apply contDiffWithinAt_const).inner ℝ
      ((h1 p hp).clm_apply (hG.clm_apply contDiffWithinAt_const))
  have h := ContDiffWithinAt.sum (fun i (_ : i ∈ Finset.univ) => hterm i)
  refine h.congr (fun q _ => ?_) ?_ <;> simp [flowPot, trOf]

end CWFlow
