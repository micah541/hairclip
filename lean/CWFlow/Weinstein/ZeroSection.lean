import CWFlow.Moser.Field

/-!
# The Weinstein primitive near the zero section

For `β = beta J F A` and `λ = lamA F A` (`Weinstein/Coord.lean`):

* `contDiffAt_lamA`, `contDiffAt_beta`: smooth where `F, A` are;
* `extDeriv_beta_zero`: `dβ = Φ₀*ω₀ − dλ_A` vanishes on the zero section;
* `curryLeft_lamA_injective`: `ι_· dλ_A` is injective on the zero section when `dλ_A` is nondegenerate there;
  in the cotangent model (`A` the metric lift, `λ_A = λ_can`) this always holds (`nondeg_of_lift`);
* `iotaT_beta_invertible`: hence `ι_· Ω_t` is invertible at `(x, 0)` for every `t`;
* `moserField_zero`: the Moser field vanishes on the zero section.
-/

open ContinuousAlternatingMap Filter

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {G : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G] [FiniteDimensional ℝ G]
variable {J : V →L[ℝ] V} {F : E → V} {A : E → G →L[ℝ] V}

theorem contDiffAt_pb1 {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {Φ : W → V} {p : W}
    (hΦ : ContDiffAt ℝ ∞ Φ p) : ContDiffAt ℝ ∞ (pb1 J Φ) p := by
  have h1 : ContDiffAt ℝ ∞ (fun q => lam0L J (Φ q)) p := (lam0L J).contDiff.contDiffAt.comp p hΦ
  exact h1.clm_comp (hΦ.fderiv_right (m := ∞) le_rfl)

theorem contDiffAt_lamA {x : E} {ξ : G} (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x) :
    ContDiffAt ℝ ∞ (lamA F A) (x, ξ) := by
  set p : E × G := (x, ξ)
  have hlift : ContDiffAt ℝ ∞ (fun q : E × G => A q.1 q.2) p :=
    (hA.comp p contDiffAt_fst).clm_apply contDiffAt_snd
  have hc : ContDiffAt ℝ ∞ (fun q : E × G => (innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) (A q.1 q.2)) p :=
    (innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).contDiff.contDiffAt.comp p hlift
  have hd : ContDiffAt ℝ ∞ (fun q : E × G => (fderiv ℝ F q.1).comp (ContinuousLinearMap.fst ℝ E G))
      p :=
    ((ContinuousLinearMap.compL ℝ (E × G) E V).flip (ContinuousLinearMap.fst ℝ E G)).contDiff.contDiffAt.comp
      p ((hF.fderiv_right (m := ∞) le_rfl).comp p contDiffAt_fst)
  exact hc.clm_comp hd

theorem contDiffAt_beta {x : E} {ξ : G} (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x) :
    ContDiffAt ℝ ∞ (beta J F A) (x, ξ) := by
  set p : E × G := (x, ξ)
  have hFst : ContDiffAt ℝ ∞ (fun q : E × G => F q.1) p := hF.comp p contDiffAt_fst
  have hlift : ContDiffAt ℝ ∞ (fun q : E × G => A q.1 q.2) p :=
    (hA.comp p contDiffAt_fst).clm_apply contDiffAt_snd
  have hPhi : ContDiffAt ℝ ∞ (Phi0 J F A) p := hFst.sub (J.contDiff.contDiffAt.comp p hlift)
  have hell : ContDiffAt ℝ ∞ (ell F A) p := contDiffAt_const.mul (hFst.inner ℝ hlift)
  exact (((contDiffAt_pb1 hPhi).sub (contDiffAt_lamA hF hA)).sub (contDiffAt_pb1 hFst)).sub
    (hell.fderiv_right (m := ∞) le_rfl)

theorem two_le_infty : (2 : WithTop ℕ∞) ≤ ∞ := WithTop.coe_le_coe.2 le_top

theorem one_add_one_le_infty : (1 : WithTop ℕ∞) + 1 ≤ ∞ := by
  rw [one_add_one_eq_two]; exact two_le_infty

/-- `dβ = 0` on the zero section. -/
theorem extDeriv_beta_zero (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) {x : E}
    (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x)
    (hL : ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0)
    (hA2 : ∀ ξ η, ⟪J (A x ξ), A x η⟫ = 0) :
    extDeriv (toAlt1 (beta J F A)) (x, 0) = 0 := by
  have hb := (contDiffAt_beta (J := J) (ξ := (0 : G)) hF hA).differentiableAt (by simp)
  have hF2 := (hF.fderiv_right (m := 1) one_add_one_le_infty).differentiableAt (by simp)
  ext v
  have hv : v = ![v 0, v 1] := by funext i; fin_cases i <;> rfl
  rw [hv, extDeriv_toAlt1 hb, dOne_beta hJi hJJ (hF.of_le two_le_infty) (hA.of_le two_le_infty) hL,
    pb2_Phi0_zero hJi hJJ (hF.differentiableAt (by simp)) (hA.differentiableAt (by simp)) hL hA2,
    dOne_lamA_zero hF2 (hA.differentiableAt (by simp))]
  simp

/-- `ι_· dλ_A` is injective on the zero section when `dλ_A(x, 0)` is nondegenerate. -/
theorem curryLeft_lamA_injective {x : E} (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x)
    (hnd : ∀ v : E × G, (∀ b : E × G, ⟪A x v.2, fderiv ℝ F x b.1⟫ - ⟪A x b.2, fderiv ℝ F x v.1⟫ = 0) →
      v = 0) :
    Function.Injective (extDeriv (toAlt1 (lamA F A)) (x, 0)).curryLeft := by
  have hl := (contDiffAt_lamA (ξ := (0 : G)) hF hA).differentiableAt (by simp)
  have hF2 := (hF.fderiv_right (m := 1) one_add_one_le_infty).differentiableAt (by simp)
  rw [injective_iff_map_eq_zero]
  intro v hv
  refine hnd v fun b => ?_
  have h := congrArg (fun f => f ![b]) hv
  simp only [curryLeft_apply_apply] at h
  have e : (Matrix.vecCons v ![b] : Fin 2 → E × G) = ![v, b] := rfl
  rw [← dOne_lamA_zero hF2 (hA.differentiableAt (by simp)), ← extDeriv_toAlt1 hl]
  simpa [e] using h

/-- Cotangent model (`G = E →L ℝ`): if `⟨A x ξ, DF u⟩ = ξ u` (`A x` is the metric lift), `dλ_A = ω_can` is
nondegenerate on the zero section. -/
theorem nondeg_of_lift {A : E → (E →L[ℝ] ℝ) →L[ℝ] V} {x : E}
    (hA1 : ∀ (ξ : E →L[ℝ] ℝ) u, ⟪A x ξ, fderiv ℝ F x u⟫ = ξ u) (v : E × (E →L[ℝ] ℝ))
    (h : ∀ b : E × (E →L[ℝ] ℝ), ⟪A x v.2, fderiv ℝ F x b.1⟫ - ⟪A x b.2, fderiv ℝ F x v.1⟫ = 0) :
    v = 0 := by
  have h2 : v.2 = 0 := by
    ext u
    have := h (u, 0)
    simpa [hA1] using this
  have h1 : v.1 = 0 := by
    refine SeparatingDual.eq_zero_of_forall_dual_eq_zero (R := ℝ) fun f => ?_
    have := h (0, f)
    simpa [hA1, h2] using this
  exact Prod.ext h1 h2

/-- At the zero section, `ι_· Ω_t` is invertible for every `t`. -/
theorem iotaT_beta_invertible (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) {x : E}
    (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x)
    (hL : ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0) (hA2 : ∀ ξ η, ⟪J (A x ξ), A x η⟫ = 0)
    (hnd : ∀ v : E × G, (∀ b : E × G, ⟪A x v.2, fderiv ℝ F x b.1⟫ - ⟪A x b.2, fderiv ℝ F x v.1⟫ = 0) →
      v = 0) (t : ℝ) :
    ∃ e : (E × G) ≃L[ℝ] ((E × G) [⋀^Fin 1]→L[ℝ] ℝ),
      (e : _ →L[ℝ] _) = iotaT (lamA F A) (beta J F A) t (x, 0) :=
  iotaT_invertible t (extDeriv_beta_zero hJi hJJ hF hA hL hA2) (curryLeft_lamA_injective hF hA hnd)

/-- The Moser field vanishes on the zero section. -/
theorem moserField_zero (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) {x : E} (hF : DifferentiableAt ℝ F x)
    (hA : DifferentiableAt ℝ A x) (t : ℝ) :
    moserField (lamA F A) (beta J F A) t (x, 0) = 0 := by
  have h0 : toAlt1 (beta J F A) (x, 0) = 0 := by
    ext v
    simp [beta_zero hJi hF hA]
  rw [moserField, h0, _root_.map_zero, neg_zero]

end CWFlow
