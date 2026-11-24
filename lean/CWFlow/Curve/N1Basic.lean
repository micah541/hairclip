/-
# The `n = 1` instance of `CWFlow`, light part (port of N1.lean, sections 0, (iii), (iv)-algebra)

Setting `E1`, `J = (I * ·)` (`Jc`), `hJi`, `hJJ`, `hdim`, "every linear map from `E1` is Lagrangian";
the `n = 1` velocity identity `flowVel Jc X = −κ_ss N` (`flowVel_eq_neg_kss`); `flowVel_conj`.
Split from `CWFlow.Curve.N1` so that it imports only `CWFlow.Scalar.*` (smaller import closure).
-/
import CWFlow.Scalar.Invariance
import Mathlib.Geometry.Manifold.Instances.Sphere   -- charts on `Circle`
import CWFlow.Curve.CircleChart   -- E1, ι1, and the instance `IsManifold 𝓘(ℝ, E1) ∞ Circle`

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold ComplexConjugate

noncomputable section

namespace CWFlow.Curve

open CWFlow

/-! ## 0. The setting -/

/- The model space of a 1-manifold.  Mathlib charts `Circle` on `EuclideanSpace ℝ (Fin 1)`
(`𝓡 1` is notation for `𝓘(ℝ, EuclideanSpace ℝ (Fin 1))`). -/
-- `E1 := EuclideanSpace ℝ (Fin 1)` is defined in `CWFlow.Curve.Defs`.

-- `[NeZero (Module.finrank ℝ E)]` is the only dimension hypothesis of the two main theorems:
--   Main.lean l.22-23:  variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
--                         [NeZero (Module.finrank ℝ E)]
instance instNeZeroFinrankE1 : NeZero (Module.finrank ℝ E1) :=
  ⟨by simp [E1]⟩   -- [GUESS] uses simp lemma `finrank_euclideanSpace_fin`; else
                   -- `⟨by rw [finrank_euclideanSpace_fin]; exact one_ne_zero⟩`

/-! Instance probes.  The manifold hypotheses of the main theorems are (Main.lean l.24-25):
    `[TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
     [CompactSpace M] [T2Space M] [SigmaCompactSpace M]`.
  If one of these probes fails, see INSTRUCTIONS.md, "likely failures", item 2. -/
example : ChartedSpace E1 Circle := inferInstance
example : IsManifold 𝓘(ℝ, E1) ∞ Circle := inferInstance
  -- [GUESS] Mathlib may only provide `IsManifold (𝓡 1) ω Circle`; then add
  -- `instance : IsManifold 𝓘(ℝ, E1) ∞ Circle := IsManifold.of_le (n := ω) le_top`
example : CompactSpace Circle := inferInstance
example : T2Space Circle := inferInstance
example : SigmaCompactSpace Circle := inferInstance
example : CompleteSpace E1 := inferInstance   -- needed by `exists_flow_of_normal_family`

/-- `J` = multiplication by `i` on `ℂ`, as an `ℝ`-linear map. -/
def Jc : ℂ →L[ℝ] ℂ := ContinuousLinearMap.mul ℝ ℂ Complex.I

@[simp] theorem Jc_apply (v : ℂ) : Jc v = Complex.I * v := by
  simp [Jc]   -- [GUESS] `ContinuousLinearMap.mul_apply'` is simp; `rfl` may also work

/-- `hJi` of the main theorems: `J` is an isometry.  Proof by real polarization. -/
theorem Jc_isom (u v : ℂ) : ⟪Jc u, Jc v⟫ = ⟪u, v⟫ := by
  simp only [Jc_apply]
  -- [GUESS: lemma name] ⟪x, y⟫_ℝ = (‖x + y‖ * ‖x + y‖ - ‖x‖ * ‖x‖ - ‖y‖ * ‖y‖) / 2
  rw [real_inner_eq_norm_add_mul_self_sub_norm_mul_self_sub_norm_mul_self_div_two,
    real_inner_eq_norm_add_mul_self_sub_norm_mul_self_sub_norm_mul_self_div_two, ← mul_add]
  simp   -- uses ‖I‖ = 1 (`Complex.norm_I`)
  -- fallback: `simp [Complex.inner, map_mul, Complex.conj_I]; ring_nf; simp [Complex.I_sq]`

/-- `hJJ` of the main theorems: `J² = -1`. -/
theorem Jc_JJ (v : ℂ) : Jc (Jc v) = -v := by
  simp only [Jc_apply]
  rw [← mul_assoc, Complex.I_mul_I]
  ring

/-- `hdim` of the main theorems: `dim V = 2 dim E`, here `2 = 2 · 1`. -/
theorem hdim1 : Module.finrank ℝ ℂ = 2 * Module.finrank ℝ E1 := by
  simp [E1, Complex.finrank_real_complex]

/-- `⟪J v, v⟫ = 0` (from `hJi`, `hJJ` only). -/
theorem inner_Jc_self (v : ℂ) : ⟪Jc v, v⟫ = 0 := by
  have h1 : ⟪Jc v, v⟫ = ⟪Jc (Jc v), Jc v⟫ := (Jc_isom _ _).symm
  rw [Jc_JJ, inner_neg_left, real_inner_comm] at h1
  -- h1 : ⟪v, Jc v⟫ = -⟪v, Jc v⟫   (rw of `real_inner_comm` hits the first inner product)
  have h2 : ⟪v, Jc v⟫ = 0 := by linarith
  rw [real_inner_comm]
  exact h2

/-- In dimension one every linear map is Lagrangian: `⟪J (L u), L w⟫ = 0`. -/
theorem lag_dim_one (L : E1 →L[ℝ] ℂ) (u w : E1) : ⟪Jc (L u), L w⟫ = 0 := by
  by_cases hu : u = 0
  · simp [hu]
  · -- [GUESS: lemma name] finrank_eq_one_iff_of_nonzero' (v) (nz : v ≠ 0) :
    --   finrank K V = 1 ↔ ∀ w : V, ∃ c : K, c • v = w
    obtain ⟨c, hc⟩ := (finrank_eq_one_iff_of_nonzero' (K := ℝ) u hu).mp
      (by simp [finrank_euclideanSpace]) w
    rw [← hc, map_smul, real_inner_smul_right, inner_Jc_self, mul_zero]

/-! ## (iii) At n = 1 the velocity `J∇div(JH)` is `−k_ss N`

Source, Scalar/Pointwise.lean l.81-90:
```
def meanCurv (X : E → V) (y : E) : V := HOf (fderiv ℝ X y) (fderiv ℝ (fderiv ℝ X) y)
def flowPot (J : V →L[ℝ] V) (X : E → V) (y : E) : ℝ :=
  trOf (fderiv ℝ X y) fun a c =>
    ⟪fderiv ℝ (fun z => J (meanCurv X z)) y a, fderiv ℝ X y c⟫
def flowVel (J : V →L[ℝ] V) (X : E → V) (y : E) : V :=
  J (fderiv ℝ X y (ginvOf (fderiv ℝ X y) (fderiv ℝ (flowPot J X) y)))
```
Hand computation (also U1_sym.py S1, exact sympy): with `T = X'/|X'|`, `N = J T`,
`k = ⟪X'', J X'⟫/|X'|³`, `∂_s = |X'|⁻¹ ∂_θ`:  `H = k N`, `JH = −k T`, `flowPot = −k_s`,
`flowVel = J (T · (−k_ss)) = −k_ss N`.  No tangential component.

The statement is made for `E = ℝ` (a curve `X : ℝ → ℂ`), where `deriv` is available.  It transfers
to the charts `E1 = EuclideanSpace ℝ (Fin 1)` by the repo's invariance lemmas
(Scalar/LinearChange.lean `flowVel_comp_linearEquiv`, Scalar/Reparam.lean `flowVel_comp`). -/

-- `curv` and `dS` are defined in `CWFlow.Curve.Defs` (same definitions as in the original N1.lean).

/-! ### Dimension-one jet computations for (iii) -/

section Dim1
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem clm_real_apply {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (L : ℝ →L[ℝ] F) (a : ℝ) :
    L a = a • L 1 := by
  rw [← map_smul, smul_eq_mul, mul_one]

theorem injective_of_apply_one {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0) : Function.Injective L := by
  intro a b hab
  rw [clm_real_apply L a, clm_real_apply L b] at hab
  exact smul_left_injective ℝ hL hab

/-- `g⁻¹ ξ = ξ(1) / |L 1|²` in dimension one. -/
theorem ginvOf_dim1 {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0) (ξ : ℝ →L[ℝ] ℝ) :
    ginvOf L ξ = ξ 1 / ‖L 1‖ ^ 2 := by
  have hinj := injective_of_apply_one hL
  apply gOf_injective hinj
  rw [gOf_ginvOf hinj]
  ext
  have hn : ‖L 1‖ ≠ 0 := norm_ne_zero_iff.mpr hL
  rw [gOf_apply, clm_real_apply L (ξ 1 / ‖L 1‖ ^ 2), real_inner_smul_left,
    real_inner_self_eq_norm_sq]
  field_simp

/-- The trace in dimension one: `tr_g B = B(1,1) / |L 1|²` for `B` bilinear. -/
theorem trOf_dim1 {W : Type*} [AddCommGroup W] [Module ℝ W] {L : ℝ →L[ℝ] V} (hL : L 1 ≠ 0)
    (B : ℝ → ℝ → W) (hB : ∀ a c, B a c = (a * c) • B 1 1) :
    trOf L B = (‖L 1‖ ^ 2)⁻¹ • B 1 1 := by
  have hs : ∑ i, bas ℝ i * (dualBas i 1 / ‖L 1‖ ^ 2) = (‖L 1‖ ^ 2)⁻¹ := by
    have h1 := sum_dualBas_smul (E := ℝ) 1
    simp only [smul_eq_mul] at h1
    calc ∑ i, bas ℝ i * (dualBas i 1 / ‖L 1‖ ^ 2)
        = (∑ i, dualBas i 1 * bas ℝ i) / ‖L 1‖ ^ 2 := by
          rw [Finset.sum_div]; refine Finset.sum_congr rfl fun i _ => ?_; ring
      _ = (‖L 1‖ ^ 2)⁻¹ := by rw [h1, one_div]
  simp only [trOf, ginvOf_dim1 hL]
  simp_rw [hB _ (_ / _)]
  rw [← Finset.sum_smul, hs]

theorem fderiv_real_apply {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (f : ℝ → F) (x a : ℝ) :
    fderiv ℝ f x a = a • deriv f x := by
  rw [clm_real_apply]; rfl


/-- The plane identity behind `JH = -k T`. -/
theorem I_mul_perp (a v : ℂ) (hv : v ≠ 0) :
    Complex.I * ((‖v‖ ^ 2)⁻¹ • (a - (⟪a, v⟫ / ‖v‖ ^ 2) • v)) =
      -(((⟪a, Complex.I * v⟫ / ‖v‖ ^ 3) / ‖v‖) • v) := by
  have hn : ‖v‖ ≠ 0 := norm_ne_zero_iff.mpr hv
  have e3 : (⟪a, Complex.I * v⟫ / ‖v‖ ^ 3) / ‖v‖ = ⟪a, Complex.I * v⟫ / (‖v‖ ^ 2) ^ 2 := by
    field_simp
  rw [e3, Complex.sq_norm, Complex.normSq_apply]
  have hs : v.re * v.re + v.im * v.im ≠ 0 := by
    rw [← Complex.normSq_apply]; exact (Complex.normSq_pos.mpr hv).ne'
  have hs2 : v.re ^ 2 + v.im ^ 2 ≠ 0 := by rw [sq, sq]; exact hs
  have hs3 : v.im ^ 2 + v.re ^ 2 ≠ 0 := by rw [add_comm]; exact hs2
  apply Complex.ext <;>
  · simp only [Complex.real_smul, Complex.inner, Complex.mul_re, Complex.mul_im, Complex.I_re,
      Complex.I_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.neg_re, Complex.neg_im, Complex.conj_re, Complex.conj_im, Complex.inv_re,
      Complex.inv_im, Complex.div_re, Complex.div_im]
    field_simp
    ring

section Curve
variable {X : ℝ → ℂ} (hX : ContDiff ℝ ∞ X) (hreg : ∀ θ, deriv X θ ≠ 0)
include hX

theorem contDiff_deriv' : ContDiff ℝ ∞ (deriv X) := (contDiff_infty_iff_deriv.mp hX).2

theorem contDiff_deriv2 : ContDiff ℝ ∞ (deriv (deriv X)) :=
  (contDiff_infty_iff_deriv.mp (contDiff_deriv' hX)).2

theorem fderiv_fderiv_real (w b c : ℝ) :
    fderiv ℝ (fderiv ℝ X) w b c = (b * c) • deriv (deriv X) w := by
  have hd : DifferentiableAt ℝ (fderiv ℝ X) w :=
    ((hX.fderiv_right (m := ∞) le_rfl).differentiable (by simp)) w
  have h1 : deriv (fderiv ℝ X) w 1 = deriv (deriv X) w := by
    have := deriv_clm_apply (u := fun _ : ℝ => (1 : ℝ)) hd (differentiableAt_const _)
    simp only [deriv_const, map_zero, add_zero] at this
    rw [← this]; rfl
  rw [fderiv_real_apply, ContinuousLinearMap.smul_apply, clm_real_apply (deriv (fderiv ℝ X) w) c,
    h1, smul_smul]

include hreg

theorem meanCurv_real (w : ℝ) :
    meanCurv X w = (‖deriv X w‖ ^ 2)⁻¹ •
      (deriv (deriv X) w - (⟪deriv (deriv X) w, deriv X w⟫ / ‖deriv X w‖ ^ 2) • deriv X w) := by
  have hL : fderiv ℝ X w 1 = deriv X w := rfl
  have hL0 : fderiv ℝ X w 1 ≠ 0 := hreg w
  rw [meanCurv, HOf, trOf_dim1 hL0]
  · rw [fderiv_fderiv_real hX, one_mul, one_smul, hL, perpOf, ginvOf_dim1 hL0,
      clm_real_apply (fderiv ℝ X w), hL]
    rfl
  · intro b c
    rw [fderiv_fderiv_real hX, fderiv_fderiv_real hX, one_mul, one_smul, perpOf_smul]

theorem JmeanCurv_real :
    (fun w => Jc (meanCurv X w)) = fun w => -((curv X w / ‖deriv X w‖) • deriv X w) := by
  funext w
  rw [Jc_apply, meanCurv_real hX hreg, I_mul_perp _ _ (hreg w)]
  rfl

theorem differentiable_curv : Differentiable ℝ (curv X) := by
  intro w
  have hv := (contDiff_deriv' hX).differentiable (by simp)
  have ha := (contDiff_deriv2 hX).differentiable (by simp)
  exact ((ha w).inner ℝ ((hv w).const_mul Complex.I)).div (((hv w).norm ℝ (hreg w)).pow 3)
    (pow_ne_zero 3 (norm_ne_zero_iff.mpr (hreg w)))

theorem flowPot_real (w : ℝ) : flowPot Jc X w = -(deriv (curv X) w / ‖deriv X w‖) := by
  have hL : fderiv ℝ X w 1 = deriv X w := rfl
  have hL0 : fderiv ℝ X w 1 ≠ 0 := hreg w
  set D := fderiv ℝ (fun z => Jc (meanCurv X z)) w with hDdef
  rw [flowPot, ← hDdef, trOf_dim1 hL0]
  swap
  · intro b c
    rw [clm_real_apply D b, clm_real_apply (fderiv ℝ X w) c, real_inner_smul_left,
      real_inner_smul_right, smul_eq_mul]
    ring
  -- derivative of `JH = -(k/n) • v`
  have hv := (contDiff_deriv' hX).differentiable (by simp)
  have hvd : HasDerivAt (deriv X) (deriv (deriv X) w) w := (hv w).hasDerivAt
  have hk : HasDerivAt (curv X) (deriv (curv X) w) w := (differentiable_curv hX hreg w).hasDerivAt
  have hnd : DifferentiableAt ℝ (fun y => ‖deriv X y‖) w := (hv w).norm ℝ (hreg w)
  have hn : HasDerivAt (fun y => ‖deriv X y‖) (deriv (fun y => ‖deriv X y‖) w) w := hnd.hasDerivAt
  have hn0 : ‖deriv X w‖ ≠ 0 := norm_ne_zero_iff.mpr (hreg w)
  have hnn : deriv (fun y => ‖deriv X y‖) w * ‖deriv X w‖ = ⟪deriv X w, deriv (deriv X) w⟫ := by
    have h1 := hvd.norm_sq
    have h2 := hn.pow 2
    have := h1.unique h2
    simp only [Nat.cast_ofNat, Nat.add_one_sub_one, pow_one] at this
    linarith [this]
  have hρ := hk.div hn hn0
  have hJ := (hρ.smul hvd).neg
  have hD1 : D 1 = -((curv X w / ‖deriv X w‖) • deriv (deriv X) w +
      ((deriv (curv X) w * ‖deriv X w‖ - curv X w * deriv (fun y => ‖deriv X y‖) w) /
        ‖deriv X w‖ ^ 2) • deriv X w) := by
    rw [hDdef, JmeanCurv_real hX hreg]
    exact hJ.deriv
  rw [hD1, hL, smul_eq_mul, inner_neg_left, inner_add_left, real_inner_smul_left,
    real_inner_smul_left, real_inner_self_eq_norm_sq, real_inner_comm (deriv X w), ← hnn]
  field_simp
  ring

theorem flowVel_real (θ : ℝ) :
    flowVel Jc X θ =
      -((dS X (dS X (curv X)) θ : ℝ) : ℂ) * (Complex.I * (deriv X θ / (‖deriv X θ‖ : ℂ))) := by
  have hL0 : fderiv ℝ X θ 1 ≠ 0 := hreg θ
  have hpot : flowPot Jc X = -(dS X (curv X)) := funext fun w => by
    rw [flowPot_real hX hreg]; rfl
  have hn0 : (‖deriv X θ‖ : ℂ) ≠ 0 := by exact_mod_cast norm_ne_zero_iff.mpr (hreg θ)
  rw [flowVel, hpot, ginvOf_dim1 hL0, clm_real_apply (fderiv ℝ X θ)]
  change Jc ((deriv (-(dS X (curv X))) θ / ‖deriv X θ‖ ^ 2) • deriv X θ) = _
  rw [deriv.neg, Jc_apply, Complex.real_smul]
  simp only [dS]
  push_cast
  field_simp

end Curve
end Dim1

/-- **(iii) n = 1: `flowVel = −k_ss N`, `N = i T`** (proved in `flowVel_real`). -/
theorem flowVel_eq_neg_kss {X : ℝ → ℂ} (hX : ContDiff ℝ ∞ X) (hreg : ∀ θ, deriv X θ ≠ 0) (θ : ℝ) :
    flowVel Jc X θ =
      -((dS X (dS X (curv X)) θ : ℝ) : ℂ) * (Complex.I * (deriv X θ / (‖deriv X θ‖ : ℂ))) :=
  -- (S1) closed 2026-11-23 (statement unchanged): `flowVel_real` above
  flowVel_real hX hreg θ

/-! ## (iv) Symmetry under reflection, via uniqueness

`σ = conj` on `ℂ` (an `ℝ`-linear isometry with `σ J = −J σ`), `ι p = p⁻¹ = conj p` on `Circle`.
Neither the closure of `IsFlow` under `σ ∘ · ∘ ι` nor the pointwise identity `flowVel_conj`
is in the repo (U1_regularity.tex already says "it is not a Lean lemma"). -/

section Conj
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
set_option linter.unusedSectionVars false

/-- `conj` as a continuous ℝ-linear map. -/
abbrev cL : ℂ →L[ℝ] ℂ := (Complex.conjLIE : ℂ →L[ℝ] ℂ)

@[simp] theorem cL_apply (w : ℂ) : cL w = conj w := rfl

theorem inner_conj_conj (u v : ℂ) : ⟪conj u, conj v⟫ = ⟪u, v⟫ := by
  have := Complex.conjLIE.inner_map_map u v
  simpa using this

theorem fderiv_conj_comp {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (X : F → ℂ) (y : F) :
    fderiv ℝ (fun z => conj (X z)) y = cL.comp (fderiv ℝ X y) := by
  have h := Complex.conjLIE.comp_fderiv (f := X) (x := y)
  exact h

theorem gOf_conj (L : E →L[ℝ] ℂ) : gOf (cL.comp L) = gOf L := by
  ext a b
  simp

theorem ginvOf_conj (L : E →L[ℝ] ℂ) : ginvOf (cL.comp L) = ginvOf L := by
  simp only [ginvOf, gOf_conj]

theorem trOf_conj {W : Type*} [AddCommMonoid W] (L : E →L[ℝ] ℂ) (B : E → E → W) :
    trOf (cL.comp L) B = trOf L B := by
  simp only [trOf, ginvOf_conj]

theorem perpOf_conj (L : E →L[ℝ] ℂ) (v : ℂ) : perpOf (cL.comp L) (conj v) = conj (perpOf L v) := by
  have h : ((innerSL ℝ : ℂ →L[ℝ] ℂ →L[ℝ] ℝ) (conj v)).comp (cL.comp L) =
      ((innerSL ℝ : ℂ →L[ℝ] ℂ →L[ℝ] ℝ) v).comp L := by
    ext a
    simp
  simp only [perpOf, h, ginvOf_conj, ContinuousLinearMap.comp_apply, cL_apply, map_sub]

theorem meanCurv_conj (X : E → ℂ) (y : E) :
    meanCurv (fun z => conj (X z)) y = conj (meanCurv X y) := by
  have h1 : fderiv ℝ (fun z => conj (X z)) = fun z => cL.comp (fderiv ℝ X z) :=
    funext (fderiv_conj_comp X)
  let A : (E →L[ℝ] ℂ) ≃L[ℝ] (E →L[ℝ] ℂ) :=
    ContinuousLinearEquiv.arrowCongr (ContinuousLinearEquiv.refl ℝ E)
      (Complex.conjLIE.toContinuousLinearEquiv)
  have hA : ∀ L, A L = cL.comp L := fun L => by ext a; simp [A]
  have h2 : fderiv ℝ (fderiv ℝ (fun z => conj (X z))) y =
      (A : (E →L[ℝ] ℂ) →L[ℝ] (E →L[ℝ] ℂ)).comp (fderiv ℝ (fderiv ℝ X) y) := by
    have e : fderiv ℝ (fun z => conj (X z)) = A ∘ fderiv ℝ X := by
      rw [h1]; funext z; simp [hA]
    rw [e]
    exact A.comp_fderiv
  have hQ : ∀ a c, fderiv ℝ (fderiv ℝ (fun z => conj (X z))) y a c =
      conj (fderiv ℝ (fderiv ℝ X) y a c) := fun a c => by
    rw [h2]; simp [hA]
  simp only [meanCurv, HOf, hQ]
  rw [fderiv_conj_comp X y]
  simp only [trOf_conj, perpOf_conj]
  simp only [trOf, map_sum]

theorem flowPot_conj (X : E → ℂ) :
    flowPot Jc (fun z => conj (X z)) = fun y => -flowPot Jc X y := by
  funext y
  have hJ : (fun z => Jc (meanCurv (fun z => conj (X z)) z)) =
      fun z => -(cL (Jc (meanCurv X z))) := by
    funext z
    rw [meanCurv_conj]
    simp only [Jc_apply, cL_apply, map_mul, Complex.conj_I]
    ring
  have hD : fderiv ℝ (fun z => Jc (meanCurv (fun z => conj (X z)) z)) y =
      -(cL.comp (fderiv ℝ (fun z => Jc (meanCurv X z)) y)) := by
    rw [hJ, fderiv_fun_neg]
    congr 1
    exact fderiv_conj_comp (fun z => Jc (meanCurv X z)) y
  simp only [flowPot, hD, fderiv_conj_comp X y, trOf_conj]
  simp only [trOf, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp

theorem flowVel_conj_any (X : E → ℂ) (y : E) :
    flowVel Jc (fun z => conj (X z)) y = conj (flowVel Jc X y) := by
  simp only [flowVel, flowPot_conj, fderiv_fun_neg, fderiv_conj_comp X y, ginvOf_conj, map_neg,
    ContinuousLinearMap.comp_apply, cL_apply, Jc_apply, map_mul, Complex.conj_I]
  ring
end Conj

/-- `flowVel (σ ∘ X) = σ (flowVel X)` for `σ = conj` (J appears twice, so the sign cancels:
`H' = σH`, `JH' = −σJH`, `flowPot' = −flowPot`, `flowVel' = J σ L g⁻¹ D(−pot) = σ flowVel`).
NOT PROVED; U1_sym.py S2 checks it at n = 1. -/
theorem flowVel_conj {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
    {X : E → ℂ} {U : Set E} (hU : IsOpen U) (hX : ContDiffOn ℝ ∞ X U)
    (hinj : ∀ y ∈ U, Function.Injective (fderiv ℝ X y)) {y : E} (hy : y ∈ U) :
    flowVel Jc (fun z => conj (X z)) y = conj (flowVel Jc X y) :=
  -- (S2) closed: the identity holds with no hypotheses at all (`flowVel_conj_any`)
  flowVel_conj_any X y

/-- info: 'CWFlow.Curve.flowVel_eq_neg_kss' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowVel_eq_neg_kss

end CWFlow.Curve
