import CWFlow.Scalar.Smooth
import CWFlow.Scalar.Multilinear

/-!
# The geometry in terms of flat jets

With second and third derivatives as multilinear maps `T₂ = D²X(y) ∈ E[×2]→L V`, `T₃ = D³X(y)`:

* `HOfM L T₂ = tr_g (T₂)^⊥`, and `meanCurv X w = HOfM (DX w) (D²X w)` for every `w` (`meanCurv_eq_HOfM`);
* `PotM J L T₂ T₃ = tr_g ⟨J DH(T₂(a, ·), T₃(a, ·, ·)), L·⟩`, and `flowPot J X y = PotM J (DX y) (D²X y) (D³X y)`
  (`flowPot_eq_potM`);
* `PotM` is smooth where `L` is injective (`contDiffAt_PotMP`).

Flat (multilinear) jets avoid the nested `E →L E →L V` types and their instance problems.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- The mean curvature from a flat second derivative. -/
def HOfM (L : E →L[ℝ] V) (T₂ : E [×2]→L[ℝ] V) : V := trOf L fun a c => perpOf L (T₂ ![a, c])

def HOfMP (p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V)) : V := HOfM p.1 p.2

/-- `T₂(a, ·)` as a linear map. -/
def slot1 (T₂ : E [×2]→L[ℝ] V) (a : E) : E →L[ℝ] V :=
  continuousMultilinearCurryFin1 ℝ E V (T₂.curryLeft a)

@[simp] theorem slot1_apply (T₂ : E [×2]→L[ℝ] V) (a b : E) : slot1 T₂ a b = T₂ ![a, b] := by
  simp only [slot1, continuousMultilinearCurryFin1_apply, ContinuousMultilinearMap.curryLeft_apply]
  congr 1
  ext i
  fin_cases i <;> rfl

/-- The potential from flat jets. -/
def PotM (J : V →L[ℝ] V) (L : E →L[ℝ] V) (T₂ : E [×2]→L[ℝ] V) (T₃ : E [×3]→L[ℝ] V) : ℝ :=
  trOf L fun a c => ⟪J (fderiv ℝ HOfMP (L, T₂) (slot1 T₂ a, T₃.curryLeft a)), L c⟫

theorem meanCurv_eq_HOfM (X : E → V) (w : E) :
    meanCurv X w = HOfM (fderiv ℝ X w) (iteratedFDeriv ℝ 2 X w) := by
  simp only [meanCurv, HOf, HOfM, iteratedFDeriv_two_apply]
  rfl

theorem contDiffAt_HOfMP {p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V)} (hL : Function.Injective p.1) :
    ContDiffAt ℝ ∞ HOfMP p := by
  have hG : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) => ginvOf p.1) p :=
    (contDiffAt_ginvOf hL).comp p contDiffAt_fst
  have hterm : ∀ i, ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) =>
      p.2 ![bas E i, ginvOf p.1 (dualBas i)]) p := fun i => by
    have hargs : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) =>
        (![bas E i, ginvOf p.1 (dualBas i)] : Fin 2 → E)) p := by
      refine contDiffAt_pi.2 fun j => ?_
      fin_cases j
      · exact contDiffAt_const
      · exact hG.clm_apply contDiffAt_const
    exact (contDiff_cml_apply 2 (Ei := fun _ => E) (G := V)).contDiffAt.comp p
      (contDiffAt_snd.prodMk hargs)
  have hperp : ∀ i, ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) =>
      perpOf p.1 (p.2 ![bas E i, ginvOf p.1 (dualBas i)])) p := fun i => by
    have hin : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) =>
        ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) (p.2 ![bas E i, ginvOf p.1 (dualBas i)])).comp p.1) p :=
      ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).contDiff.contDiffAt.comp p (hterm i)).clm_comp
        contDiffAt_fst
    exact (hterm i).sub (contDiffAt_fst.clm_apply (hG.clm_apply hin))
  have h := ContDiffAt.sum (fun i (_ : i ∈ Finset.univ) => hperp i)
  refine h.congr_of_eventuallyEq (Eventually.of_forall fun p' => ?_)
  simp [HOfMP, HOfM, trOf]

/-- **`div(JH)` from the flat 3-jet.** -/
theorem flowPot_eq_potM (J : V →L[ℝ] V) {X : E → V} {S : Set E} (hS : IsOpen S)
    (hX : ContDiffOn ℝ ∞ X S) (hinj : ∀ y ∈ S, Function.Injective (fderiv ℝ X y)) {y : E}
    (hy : y ∈ S) :
    flowPot J X y = PotM J (fderiv ℝ X y) (iteratedFDeriv ℝ 2 X y) (iteratedFDeriv ℝ 3 X y) := by
  have hXy := hX.contDiffAt (hS.mem_nhds hy)
  have h1 := (hXy.fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)
  have h2 := (hXy.iteratedFDeriv_right (m := ∞) (i := 2) (by norm_cast)).differentiableAt (by simp)
  have hjet := h1.hasFDerivAt.prodMk h2.hasFDerivAt
  have hH := (contDiffAt_HOfMP (V := V) (E := E) (p := (fderiv ℝ X y, iteratedFDeriv ℝ 2 X y))
    (hinj y hy)).differentiableAt (by simp)
  have hin := HasFDerivAt.comp (g := HOfMP)
    (f := fun w => (fderiv ℝ X w, iteratedFDeriv ℝ 2 X w)) y hH.hasFDerivAt hjet
  have hmc := HasFDerivAt.comp (g := ⇑J) y J.hasFDerivAt hin
  have e : (fun z => J (meanCurv X z)) =
      ⇑J ∘ (HOfMP ∘ fun w => (fderiv ℝ X w, iteratedFDeriv ℝ 2 X w)) := by
    funext z
    simp only [Function.comp_apply, meanCurv_eq_HOfM]
    rfl
  have hA : ∀ a, fderiv ℝ (fderiv ℝ X) y a = slot1 (iteratedFDeriv ℝ 2 X y) a := fun a => by
    ext b
    rw [slot1_apply, iteratedFDeriv_two_apply]
    rfl
  have hB : ∀ a, fderiv ℝ (iteratedFDeriv ℝ 2 X) y a = (iteratedFDeriv ℝ 3 X y).curryLeft a :=
    fun a => by
      ext m
      rw [ContinuousMultilinearMap.curryLeft_apply, iteratedFDeriv_succ_apply_left]
      simp
  simp only [flowPot, PotM, e, hmc.fderiv, ContinuousLinearMap.comp_apply,
    ContinuousLinearMap.prod_apply, hA, hB]

/-- `PotM` as a function of the flat 3-jet. -/
def PotMP (J : V →L[ℝ] V) (p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V)) : ℝ :=
  PotM J p.1 p.2.1 p.2.2

theorem contDiffAt_PotMP (J : V →L[ℝ] V)
    {p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V)} (hL : Function.Injective p.1) :
    ContDiffAt ℝ ∞ (PotMP J) p := by
  have hDH := (contDiffAt_HOfMP (p := (p.1, p.2.1)) hL).fderiv_right (m := ∞) le_rfl
  have hDH' := ContDiffAt.comp (g := fderiv ℝ (HOfMP (V := V) (E := E)))
    (f := fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V) => (p.1, p.2.1)) p hDH
    (contDiffAt_fst.prodMk (contDiffAt_fst.comp p contDiffAt_snd))
  have hG : ContDiffAt ℝ ∞
      (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V) => ginvOf p.1) p :=
    (contDiffAt_ginvOf hL).comp p contDiffAt_fst
  have hterm : ∀ i, ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V) =>
      ⟪J (fderiv ℝ HOfMP (p.1, p.2.1) (slot1 p.2.1 (bas E i), p.2.2.curryLeft (bas E i))),
        p.1 (ginvOf p.1 (dualBas i))⟫) p := fun i => by
    have h1 : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V) =>
        slot1 p.2.1 (bas E i)) p :=
      (continuousMultilinearCurryFin1 ℝ E V).contDiff.contDiffAt.comp p
        (((continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin 2 => E) V).contDiff.contDiffAt.comp p
          (contDiffAt_fst.comp p contDiffAt_snd)).clm_apply contDiffAt_const)
    have h2 : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × (E [×2]→L[ℝ] V) × (E [×3]→L[ℝ] V) =>
        p.2.2.curryLeft (bas E i)) p :=
      ((continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin 3 => E) V).contDiff.contDiffAt.comp p
        (contDiffAt_snd.comp p contDiffAt_snd)).clm_apply contDiffAt_const
    exact (J.contDiff.contDiffAt.comp p (hDH'.clm_apply (h1.prodMk h2))).inner ℝ
      (contDiffAt_fst.clm_apply (hG.clm_apply contDiffAt_const))
  have h := ContDiffAt.sum (fun i (_ : i ∈ Finset.univ) => hterm i)
  refine h.congr_of_eventuallyEq (Eventually.of_forall fun p' => ?_)
  simp [PotMP, PotM, trOf]

end CWFlow
