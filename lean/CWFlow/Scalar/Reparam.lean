import CWFlow.Scalar.Invariance

/-!
# The flow quantities are invariant under reparametrization

Let `τ` be smooth with invertible derivative on an open `N`, and let `X` be smooth with injective derivative
at the points `τ z`, `z ∈ N`. Then, on `N`:

* `meanCurv (X ∘ τ) = meanCurv X ∘ τ` (`meanCurv_comp`), from `D²(X∘τ) = D²X(Dτ·, Dτ·) + DX(D²τ)`
  (`fderiv_fderiv_comp_apply`) and `HOf_comp_equiv`;
* `flowPot (X ∘ τ) = flowPot X ∘ τ` (`flowPot_comp`);
* `flowVel (X ∘ τ) = flowVel X ∘ τ` (`flowVel_comp`).

So `div(JH)` and the velocity `J∇div(JH)` are well defined on a manifold, computed in any chart.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- `z ↦ J H(z)` is differentiable where `X` is smooth with injective derivative. -/
theorem differentiableAt_J_meanCurv (J : V →L[ℝ] V) {X : E → V} {y : E} (hX : ContDiffAt ℝ ∞ X y)
    (hinj : Function.Injective (fderiv ℝ X y)) :
    DifferentiableAt ℝ (fun v => J (meanCurv X v)) y := by
  have h1 := hX.fderiv_right (m := ∞) le_rfl
  have h2 := h1.fderiv_right (m := ∞) le_rfl
  have hjet := h1.prodMk h2
  have hH := contDiffAt_HOfP (V := V) (E := E) (p := (fderiv ℝ X y, fderiv ℝ (fderiv ℝ X) y)) hinj
  have hc := ContDiffAt.comp (g := HOfP) (f := fun v => (fderiv ℝ X v, fderiv ℝ (fderiv ℝ X) v))
    y hH hjet
  have e : (fun v => J (meanCurv X v)) =
      ⇑J ∘ (HOfP ∘ fun v => (fderiv ℝ X v, fderiv ℝ (fderiv ℝ X) v)) := rfl
  rw [e]
  exact (J.contDiff.contDiffAt.comp y hc).differentiableAt (by simp)

variable {X : E → V} {τ : E → E} {N : Set E} (hN : IsOpen N) (hτ : ContDiffOn ℝ ∞ τ N)
  (hτinv : ∀ z ∈ N, ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ τ z)
  (hX : ∀ z ∈ N, ContDiffAt ℝ ∞ X (τ z))
  (hinj : ∀ z ∈ N, Function.Injective (fderiv ℝ X (τ z)))
include hN hτ hX

omit [FiniteDimensional ℝ E] in
theorem fderiv_comp_eventually {z : E} (hz : z ∈ N) :
    fderiv ℝ (X ∘ τ) =ᶠ[𝓝 z] fun w => (fderiv ℝ X (τ w)).comp (fderiv ℝ τ w) := by
  filter_upwards [hN.mem_nhds hz] with w hw
  exact fderiv_comp w ((hX w hw).differentiableAt (by simp))
    ((hτ.contDiffAt (hN.mem_nhds hw)).differentiableAt (by simp))

omit [FiniteDimensional ℝ E] in
/-- `D²(X∘τ)(a, c) = D²X(Dτ a, Dτ c) + DX(D²τ(a, c))`. -/
theorem fderiv_fderiv_comp_apply {z : E} (hz : z ∈ N) (a c : E) :
    fderiv ℝ (fderiv ℝ (X ∘ τ)) z a c =
      fderiv ℝ (fderiv ℝ X) (τ z) (fderiv ℝ τ z a) (fderiv ℝ τ z c) +
        fderiv ℝ X (τ z) (fderiv ℝ (fderiv ℝ τ) z a c) := by
  have hτz := hτ.contDiffAt (hN.mem_nhds hz)
  have hτd := hτz.differentiableAt (by simp)
  have hDX := ((hX z hz).fderiv_right (m := 1) one_add_one_le_infty).differentiableAt (by simp)
  have hDτ := (hτz.fderiv_right (m := 1) one_add_one_le_infty).differentiableAt (by simp)
  have hA := hDX.hasFDerivAt.comp z hτd.hasFDerivAt
  have h : HasFDerivAt (fun w => (fderiv ℝ X (τ w)).comp (fderiv ℝ τ w)) _ z :=
    hA.clm_comp hDτ.hasFDerivAt
  rw [(fderiv_comp_eventually hN hτ hX hz).fderiv_eq, h.fderiv]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_comp', Function.comp_apply,
    ContinuousLinearMap.compL_apply, ContinuousLinearMap.flip_apply]
  abel

include hτinv hinj

theorem meanCurv_comp {z : E} (hz : z ∈ N) : meanCurv (X ∘ τ) z = meanCurv X (τ z) := by
  obtain ⟨e, he⟩ := hτinv z hz
  have hL : fderiv ℝ (X ∘ τ) z = (fderiv ℝ X (τ z)).comp (e : E →L[ℝ] E) := by
    rw [(fderiv_comp_eventually hN hτ hX hz).self_of_nhds, he]
  rw [meanCurv, hL]
  refine HOf_comp_equiv (hinj z hz) e _ _ (fun a c => fderiv ℝ (fderiv ℝ τ) z a c) fun a c => ?_
  rw [fderiv_fderiv_comp_apply hN hτ hX hz, ← he]
  rfl

theorem meanCurv_comp_eventually {z : E} (hz : z ∈ N) :
    meanCurv (X ∘ τ) =ᶠ[𝓝 z] fun w => meanCurv X (τ w) := by
  filter_upwards [hN.mem_nhds hz] with w hw
  exact meanCurv_comp hN hτ hτinv hX hinj hw

variable (J : V →L[ℝ] V)

theorem flowPot_comp {z : E} (hz : z ∈ N) : flowPot J (X ∘ τ) z = flowPot J X (τ z) := by
  obtain ⟨e, he⟩ := hτinv z hz
  have hτd := (hτ.contDiffAt (hN.mem_nhds hz)).differentiableAt (by simp)
  have hL : fderiv ℝ (X ∘ τ) z = (fderiv ℝ X (τ z)).comp (e : E →L[ℝ] E) := by
    rw [(fderiv_comp_eventually hN hτ hX hz).self_of_nhds, he]
  have hev : (fun w => J (meanCurv (X ∘ τ) w)) =ᶠ[𝓝 z] (fun v => J (meanCurv X v)) ∘ τ := by
    filter_upwards [meanCurv_comp_eventually hN hτ hτinv hX hinj hz] with w hw
    simp only [Function.comp_apply, hw]
  have hM : fderiv ℝ (fun w => J (meanCurv (X ∘ τ) w)) z =
      (fderiv ℝ (fun v => J (meanCurv X v)) (τ z)).comp (e : E →L[ℝ] E) := by
    rw [hev.fderiv_eq, fderiv_comp z (differentiableAt_J_meanCurv J (hX z hz) (hinj z hz)) hτd, he]
  set M := fderiv ℝ (fun v => J (meanCurv X v)) (τ z)
  set L := fderiv ℝ X (τ z)
  let B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ := LinearMap.mk₂ ℝ (fun a c => ⟪M a, L c⟫)
    (fun a a' c => by simp [inner_add_left]) (fun r a c => by simp [inner_smul_left])
    (fun a c c' => by simp [inner_add_right]) (fun r a c => by simp [inner_smul_right])
  have h := trOf_comp_equiv (hinj z hz) e B
  simp only [flowPot, hL, hM]
  exact h

theorem flowPot_comp_eventually {z : E} (hz : z ∈ N) :
    flowPot J (X ∘ τ) =ᶠ[𝓝 z] fun w => flowPot J X (τ w) := by
  filter_upwards [hN.mem_nhds hz] with w hw
  exact flowPot_comp hN hτ hτinv hX hinj J hw

/-- **Invariance of the velocity** `J∇div(JH)`. -/
theorem flowVel_comp {z : E} (hz : z ∈ N)
    (hpot : DifferentiableAt ℝ (flowPot J X) (τ z)) :
    flowVel J (X ∘ τ) z = flowVel J X (τ z) := by
  obtain ⟨e, he⟩ := hτinv z hz
  have hτd := (hτ.contDiffAt (hN.mem_nhds hz)).differentiableAt (by simp)
  have hL : fderiv ℝ (X ∘ τ) z = (fderiv ℝ X (τ z)).comp (e : E →L[ℝ] E) := by
    rw [(fderiv_comp_eventually hN hτ hX hz).self_of_nhds, he]
  have hD : fderiv ℝ (flowPot J (X ∘ τ)) z =
      (fderiv ℝ (flowPot J X) (τ z)).comp (e : E →L[ℝ] E) := by
    have hev : flowPot J (X ∘ τ) =ᶠ[𝓝 z] flowPot J X ∘ τ :=
      flowPot_comp_eventually hN hτ hτinv hX hinj J hz
    rw [hev.fderiv_eq, fderiv_comp z hpot hτd, he]
  simp only [flowVel, hL, hD, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [apply_ginvOf_comp (hinj z hz) e]
  congr 3
  ext a
  simp

end CWFlow
