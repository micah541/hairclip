import CWFlow.Scalar.Reparam

/-!
# Smoothness of `H`, `div(JH)` and the velocity

For `X` smooth on an open `S` with injective derivative there:
* `meanCurv X` is smooth on `S` (`contDiffOn_meanCurv`), as `HOfP` composed with the 2-jet;
* `flowPot J X = div(JH)` is smooth on `S` (`contDiffOn_flowPot`);
* the velocity is invariant under reparametrization with no extra hypothesis (`flowVel_comp'`).
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {X : E → V} {S : Set E} (hS : IsOpen S) (hX : ContDiffOn ℝ ∞ X S)
  (hinj : ∀ y ∈ S, Function.Injective (fderiv ℝ X y))
include hS hX hinj

theorem contDiffOn_meanCurv : ContDiffOn ℝ ∞ (meanCurv X) S := by
  intro y hy
  have hXy := hX.contDiffAt (hS.mem_nhds hy)
  have h1 := hXy.fderiv_right (m := ∞) le_rfl
  have h2 := h1.fderiv_right (m := ∞) le_rfl
  have hH := contDiffAt_HOfP (V := V) (E := E) (p := (fderiv ℝ X y, fderiv ℝ (fderiv ℝ X) y))
    (hinj y hy)
  have hc := ContDiffAt.comp (g := HOfP) (f := fun v => (fderiv ℝ X v, fderiv ℝ (fderiv ℝ X) v))
    y hH (h1.prodMk h2)
  exact hc.contDiffWithinAt

theorem contDiffOn_flowPot (J : V →L[ℝ] V) : ContDiffOn ℝ ∞ (flowPot J X) S := by
  have hJH : ContDiffOn ℝ ∞ (fun z => J (meanCurv X z)) S :=
    J.contDiff.comp_contDiffOn (contDiffOn_meanCurv hS hX hinj)
  have hD := hJH.fderiv_of_isOpen hS (m := ∞) le_rfl
  have hDX := hX.fderiv_of_isOpen hS (m := ∞) le_rfl
  intro y hy
  have hG : ContDiffAt ℝ ∞ (fun z => ginvOf (fderiv ℝ X z)) y :=
    (contDiffAt_ginvOf (hinj y hy)).comp y (hDX.contDiffAt (hS.mem_nhds hy))
  have hterm : ∀ i, ContDiffAt ℝ ∞ (fun z =>
      ⟪fderiv ℝ (fun w => J (meanCurv X w)) z (bas E i),
        fderiv ℝ X z (ginvOf (fderiv ℝ X z) (dualBas i))⟫) y := fun i =>
    ((hD.contDiffAt (hS.mem_nhds hy)).clm_apply contDiffAt_const).inner ℝ
      ((hDX.contDiffAt (hS.mem_nhds hy)).clm_apply (hG.clm_apply contDiffAt_const))
  have h := ContDiffAt.sum (fun i (_ : i ∈ Finset.univ) => hterm i)
  refine (h.congr_of_eventuallyEq (Eventually.of_forall fun z => ?_)).contDiffWithinAt
  simp [flowPot, trOf]

/-- **Invariance of the velocity `J∇div(JH)` under reparametrization.** -/
theorem flowVel_comp' (J : V →L[ℝ] V) {τ : E → E} {N : Set E} (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hτinv : ∀ z ∈ N, ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ τ z)
    (hmaps : MapsTo τ N S) {z : E} (hz : z ∈ N) :
    flowVel J (X ∘ τ) z = flowVel J X (τ z) :=
  flowVel_comp hN hτ hτinv (fun w hw => hX.contDiffAt (hS.mem_nhds (hmaps hw)))
    (fun w hw => hinj _ (hmaps hw)) J hz
    (((contDiffOn_flowPot hS hX hinj J).contDiffAt (hS.mem_nhds (hmaps hz))).differentiableAt
      (by simp))

end CWFlow
