import CWFlow.Weinstein.Lift

/-!
# Immersion geometry as explicit functions of jets

For an immersion `X : E → V` (`E` finite-dimensional, no inner product on `E`), all the quantities in
the flow `∂ₜX = J∇div(JH)` are written as explicit functions of the jets `L = DX`, `Q = D²X`, `R = D³X`
at a point, smooth where `L` is injective:

* `gOf L : E →L E*`, `g(a, b) = ⟨L a, L b⟩`, and `ginvOf L = g⁻¹`;
* `trOf L B = Σᵢ B(bᵢ, g⁻¹ εᵢ) = gⁱʲ B(bᵢ, bⱼ)`, for a basis `bᵢ` of `E` with dual basis `εᵢ`;
* `perpOf L v = v − L g⁻¹⟨v, L·⟩`, the normal part;
* `HOf L Q = tr_g (Q)^⊥`, the mean curvature;
* `PotOf L Q R = tr_g ⟨J DH[(Q·, R·)], L·⟩`, the potential `div(JH)`.

For a map `X`, `meanCurv X = HOf (DX) (D²X)`, `flowPot X = div_X(JH)` (Whitney's `tdiv` formula
`gⁱʲ⟨∂ᵢ(JH), ∂ⱼX⟩`), and `flowVel X = J∇flowPot`. Then `flowPot X y = PotOf (DX y) (D²X y) (D³X y)`
(`flowPot_eq_potOf`), and the jet functions are smooth (`contDiffAt_ginvOf`, `contDiffAt_HOf`).
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

variable (V E) in
/-- The type of second derivatives `D²X(y)`, taken from `fderiv` so that it carries the normed instances
(a written `E →L[ℝ] E →L[ℝ] V` picks the non-defeq-by-unification TVS topology on the inner space). -/
abbrev Hess := type_of% (fderiv ℝ (fun _ : E => (0 : E →L[ℝ] V)) 0)

variable (V E) in
/-- The type of third derivatives `D³X(y)`. -/
abbrev Hess3 := type_of% (fderiv ℝ (fun _ : E => (0 : Hess V E)) 0)

/-- The basis of `E` used for traces. -/
abbrev bas (E : Type*) [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E] :=
  Module.finBasis ℝ E

/-- The dual basis, as continuous functionals. -/
def dualBas (i : Fin (Module.finrank ℝ E)) : E →L[ℝ] ℝ :=
  LinearMap.toContinuousLinearMap ((bas E).coord i)

/-- The metric `g(a, b) = ⟨L a, L b⟩` of a linear map `L : E → V`. -/
def gOf (L : E →L[ℝ] V) : E →L[ℝ] (E →L[ℝ] ℝ) :=
  ((ContinuousLinearMap.compL ℝ E V ℝ).flip L).comp ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).comp L)

@[simp] theorem gOf_apply (L : E →L[ℝ] V) (a b : E) : gOf L a b = ⟪L a, L b⟫ := rfl

theorem metricMap_eq_gOf (X : E → V) (y : E) : metricMap X y = gOf (fderiv ℝ X y) := rfl

/-- `g⁻¹ : E* → E`. -/
def ginvOf (L : E →L[ℝ] V) : (E →L[ℝ] ℝ) →L[ℝ] E := ContinuousLinearMap.inverse (gOf L)

/-- The trace `gⁱʲ B(bᵢ, bⱼ) = Σᵢ B(bᵢ, g⁻¹ εᵢ)`, for any `B : E → E → W`. -/
def trOf {W : Type*} [AddCommMonoid W] (L : E →L[ℝ] V) (B : E → E → W) : W :=
  ∑ i, B (bas E i) (ginvOf L (dualBas i))

/-- The normal part `v − L g⁻¹ ⟨v, L ·⟩`. -/
def perpOf (L : E →L[ℝ] V) (v : V) : V :=
  v - L (ginvOf L (((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) v).comp L))

/-- The mean curvature `tr_g Q^⊥` for `L = DX`, `Q = D²X`. -/
def HOf (L : E →L[ℝ] V) (Q : Hess V E) : V := trOf L fun a c => perpOf L (Q a c)

/-- `HOf` as a function of the pair `(L, Q)`. -/
def HOfP (p : (E →L[ℝ] V) × Hess V E) : V := HOf p.1 p.2

/-- The potential `div(JH) = tr_g ⟨J DH(Q·, R·), L·⟩` from the 3-jet `(L, Q, R)`. -/
def PotOf (J : V →L[ℝ] V) (L : E →L[ℝ] V) (Q : Hess V E) (R : Hess3 V E) : ℝ :=
  trOf L fun a c => ⟪J (fderiv ℝ HOfP (L, Q) (Q a, R a)), L c⟫

/-! ### The flow quantities of a map -/

/-- The mean curvature vector of `X`. -/
def meanCurv (X : E → V) (y : E) : V := HOf (fderiv ℝ X y) (fderiv ℝ (fderiv ℝ X) y)

/-- The potential `div_X(JH) = gⁱʲ ⟨∂ᵢ(JH), ∂ⱼX⟩`. -/
def flowPot (J : V →L[ℝ] V) (X : E → V) (y : E) : ℝ :=
  trOf (fderiv ℝ X y) fun a c =>
    ⟪fderiv ℝ (fun z => J (meanCurv X z)) y a, fderiv ℝ X y c⟫

/-- The velocity `J∇_g div(JH)` of the Hamiltonian gradient flow of volume. -/
def flowVel (J : V →L[ℝ] V) (X : E → V) (y : E) : V :=
  J (fderiv ℝ X y (ginvOf (fderiv ℝ X y) (fderiv ℝ (flowPot J X) y)))

theorem flowVel_eq_liftA (J : V →L[ℝ] V) (X : E → V) (y : E) :
    flowVel J X y = J (liftA X y (fderiv ℝ (flowPot J X) y)) := rfl

/-! ### Smoothness of the jet functions -/

theorem gOf_injective {L : E →L[ℝ] V} (hL : Function.Injective L) : Function.Injective (gOf L) := by
  rw [injective_iff_map_eq_zero]
  intro u hu
  have h := congrArg (fun f => f u) hu
  simp only [gOf_apply, ContinuousLinearMap.zero_apply, inner_self_eq_zero] at h
  exact hL (by rw [h, map_zero])

/-- Where `L` is injective, `gOf L` is an isomorphism. -/
def gEquiv (L : E →L[ℝ] V) (hL : Function.Injective L) : E ≃L[ℝ] (E →L[ℝ] ℝ) :=
  (LinearMap.linearEquivOfInjective (gOf L : E →ₗ[ℝ] (E →L[ℝ] ℝ)) (gOf_injective hL)
    (by
      rw [← (LinearMap.toContinuousLinearMap : Module.Dual ℝ E ≃ₗ[ℝ] (E →L[ℝ] ℝ)).finrank_eq,
        Subspace.dual_finrank_eq])).toContinuousLinearEquiv

theorem coe_gEquiv {L : E →L[ℝ] V} (hL : Function.Injective L) :
    (gEquiv L hL : E →L[ℝ] (E →L[ℝ] ℝ)) = gOf L :=
  ContinuousLinearMap.ext fun _ => rfl

theorem gOf_ginvOf {L : E →L[ℝ] V} (hL : Function.Injective L) (ξ : E →L[ℝ] ℝ) :
    gOf L (ginvOf L ξ) = ξ := by
  rw [ginvOf, ← coe_gEquiv hL, ContinuousLinearMap.inverse_equiv]
  exact (gEquiv L hL).apply_symm_apply ξ

/-- `⟨L g⁻¹ξ, L a⟩ = ξ a`. -/
theorem inner_ginvOf {L : E →L[ℝ] V} (hL : Function.Injective L) (ξ : E →L[ℝ] ℝ) (a : E) :
    ⟪L (ginvOf L ξ), L a⟫ = ξ a := by
  rw [← gOf_apply, gOf_ginvOf hL]

theorem contDiff_gOf : ContDiff ℝ ∞ (gOf : (E →L[ℝ] V) → E →L[ℝ] (E →L[ℝ] ℝ)) := by
  have h1 := (ContinuousLinearMap.contDiff (𝕜 := ℝ) (E := E →L[ℝ] V)
    (F := (V →L[ℝ] ℝ) →L[ℝ] E →L[ℝ] ℝ) (n := ∞) ((ContinuousLinearMap.compL ℝ E V ℝ).flip))
  have h2 := (ContinuousLinearMap.contDiff (𝕜 := ℝ) (E := E →L[ℝ] V)
    (F := E →L[ℝ] V →L[ℝ] ℝ) (n := ∞) ((ContinuousLinearMap.compL ℝ E V (V →L[ℝ] ℝ))
    (innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ)))
  exact h1.clm_comp h2

theorem contDiffAt_ginvOf {L : E →L[ℝ] V} (hL : Function.Injective L) :
    ContDiffAt ℝ ∞ (ginvOf : (E →L[ℝ] V) → (E →L[ℝ] ℝ) →L[ℝ] E) L := by
  have h0 := contDiffAt_map_inverse (𝕜 := ℝ) (n := ∞) (gEquiv L hL)
  rw [coe_gEquiv hL] at h0
  exact h0.comp L contDiff_gOf.contDiffAt

theorem isOpen_injective : IsOpen {L : E →L[ℝ] V | Function.Injective L} := by
  rw [isOpen_iff_mem_nhds]
  intro L hL
  have h := (contDiff_gOf (E := E) (V := V)).continuous.continuousAt (x := L)
  have hmem : gOf L ∈ Set.range ((↑) : (E ≃L[ℝ] (E →L[ℝ] ℝ)) → _) :=
    ⟨gEquiv L hL, coe_gEquiv hL⟩
  filter_upwards [h.preimage_mem_nhds (ContinuousLinearEquiv.isOpen.mem_nhds hmem)] with L' hL'
  obtain ⟨e, he⟩ := hL'
  intro a b hab
  apply e.injective
  change (e : E →L[ℝ] (E →L[ℝ] ℝ)) a = (e : E →L[ℝ] (E →L[ℝ] ℝ)) b
  rw [he]
  ext c
  simp [hab]

theorem contDiffAt_HOfP {p : (E →L[ℝ] V) × Hess V E}
    (hL : Function.Injective p.1) : ContDiffAt ℝ ∞ HOfP p := by
  have hG : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × Hess V E => ginvOf p.1) p :=
    (contDiffAt_ginvOf hL).comp p contDiffAt_fst
  have hterm : ∀ i, ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × Hess V E =>
      p.2 (bas E i) (ginvOf p.1 (dualBas i))) p := fun i =>
    (contDiffAt_snd.clm_apply contDiffAt_const).clm_apply (hG.clm_apply contDiffAt_const)
  have hperp : ∀ i, ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × Hess V E =>
      perpOf p.1 (p.2 (bas E i) (ginvOf p.1 (dualBas i)))) p := fun i => by
    have hin : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × Hess V E =>
        ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) (p.2 (bas E i) (ginvOf p.1 (dualBas i)))).comp p.1) p :=
      ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).contDiff.contDiffAt.comp p (hterm i)).clm_comp
        contDiffAt_fst
    exact (hterm i).sub (contDiffAt_fst.clm_apply (hG.clm_apply hin))
  have h := ContDiffAt.sum (fun i (_ : i ∈ Finset.univ) => hperp i)
  refine h.congr_of_eventuallyEq (Eventually.of_forall fun p' => ?_)
  simp [HOfP, HOf, trOf]

/-! ### The potential as a function of the 3-jet -/

/-- `flowPot X y = PotOf J (DX y) (D²X y) (D³X y)`. -/
theorem flowPot_eq_potOf (J : V →L[ℝ] V) {X : E → V} {y : E} (hX : ContDiffAt ℝ 3 X y)
    (hinj : Function.Injective (fderiv ℝ X y)) :
    flowPot J X y = PotOf J (fderiv ℝ X y) (fderiv ℝ (fderiv ℝ X) y)
      (fderiv ℝ (fderiv ℝ (fderiv ℝ X)) y) := by
  have h1 :=
    (hX.fderiv_right (m := 2) (by norm_num)).differentiableAt (by norm_num)
  have h2 :=
    ((hX.fderiv_right (m := 2) (by norm_num)).fderiv_right (m := 1) (by norm_num)).differentiableAt
      (by norm_num)
  have hjet := h1.hasFDerivAt.prodMk h2.hasFDerivAt
  have hH := (contDiffAt_HOfP (V := V) (E := E) (p := (fderiv ℝ X y, fderiv ℝ (fderiv ℝ X) y))
    hinj).differentiableAt (by simp)
  have hin := HasFDerivAt.comp (g := HOfP)
    (f := fun z => (fderiv ℝ X z, fderiv ℝ (fderiv ℝ X) z)) y hH.hasFDerivAt hjet
  have hmc := HasFDerivAt.comp (g := ⇑J) y J.hasFDerivAt hin
  have e : (fun z => J (meanCurv X z)) =
      ⇑J ∘ (HOfP ∘ fun z => (fderiv ℝ X z, fderiv ℝ (fderiv ℝ X) z)) := rfl
  simp only [flowPot, PotOf, e, hmc.fderiv]
  rfl

end CWFlow
