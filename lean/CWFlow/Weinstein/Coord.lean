import CWFlow.Moser.Chart

/-!
# The immersed Weinstein data in coordinates

`V` is a real inner product space with a complex structure `J` (`J² = −1`, `J` orthogonal); `λ₀(z) = ½⟨Jz, ·⟩`,
`ω₀(u, v) = ⟨Ju, v⟩ = dλ₀`. In a chart `x ∈ E` of `L`, with bundle coordinates `(x, ξ) ∈ E × G` (`G = E →L ℝ` for `T*L`).

Let `F : E → V` be the immersion and `A x : E →L V` the lift `ξ ↦ gⁱʲξᵢ∂ⱼF` (abstractly: `⟨A x ξ, DF u⟩ = ⟨ξ, u⟩`
and `⟨J A x ξ, A x η⟩ = 0`). Put
`Φ₀(x, ξ) = F x − J (A x ξ)`, `ℓ(x, ξ) = −½⟨F x, A x ξ⟩`, `β = Φ₀*λ₀ − λ_can − π*(F*λ₀) − dℓ`.

* `dOne_pb1`: `d(Φ*λ₀) = Φ*ω₀` for any `C²` map `Φ`;
* `dOne_beta`: `dβ = Φ₀*ω₀ − dλ_A` when `F` is Lagrangian;
* `beta_zero`: `β = 0` at every point `(x, 0)` of the zero section;
* `pb2_Phi0_zero`, `dOne_lamA_zero`: `Φ₀*ω₀ = dλ_A` at `(x, 0)`.

Here `λ_A(x, w)(a) = ⟨A x w, DF a₁⟩`. For the metric lift `A` (cotangent coordinates) `λ_A = λ_can`; for `A = DF`
(tangent coordinates, `Φ₀(x, v) = F x − J DF v`, used for the global construction on `TL`) `λ_A = g(w, dx)`.
-/

open scoped InnerProductSpace RealInnerProductSpace Topology

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
variable (J : V →L[ℝ] V)

/-- `z ↦ λ₀(z) = ½⟨Jz, ·⟩`, a linear map into `V →L ℝ`. -/
def lam0L : V →L[ℝ] V →L[ℝ] ℝ := (1 / 2 : ℝ) • (innerSL ℝ).comp J

/-- The pullback `Φ*λ₀`, a 1-form field on `W`. -/
def pb1 (Φ : W → V) (p : W) : W →L[ℝ] ℝ := (lam0L J (Φ p)).comp (fderiv ℝ Φ p)

variable {J}

theorem inner_J_swap (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (u v : V) :
    ⟪J u, v⟫ = -⟪J v, u⟫ := by
  have h := hJi u (J v)
  rw [hJJ, inner_neg_right] at h
  rw [← real_inner_comm (J v) u]
  linarith

/-- `d(Φ*λ₀) = Φ*ω₀`. -/
theorem dOne_pb1 (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) {Φ : W → V} {p : W}
    (hΦ : ContDiffAt ℝ 2 Φ p) (a b : W) :
    dOne (pb1 J Φ) p a b = ⟪J (fderiv ℝ Φ p a), fderiv ℝ Φ p b⟫ := by
  have hd : HasFDerivAt (fderiv ℝ Φ) (fderiv ℝ (fderiv ℝ Φ) p) p :=
    ((hΦ.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)).hasFDerivAt
  have hc : HasFDerivAt (fun q => lam0L J (Φ q)) ((lam0L J).comp (fderiv ℝ Φ p)) p :=
    (lam0L J).hasFDerivAt.comp p (hΦ.differentiableAt (by norm_num)).hasFDerivAt
  have h := hc.clm_comp hd
  have hsymm := hΦ.isSymmSndFDerivAt (by simp) a b
  rw [dOne, show pb1 J Φ = fun q => (lam0L J (Φ q)).comp (fderiv ℝ Φ q) from rfl, h.fderiv]
  simp only [ContinuousLinearMap.coe_comp', Function.comp_apply, lam0L,
    ContinuousLinearMap.compL_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_smul', Pi.smul_apply,
    innerSL_apply_apply, smul_eq_mul, ContinuousLinearMap.flip_apply]
  rw [hsymm, inner_J_swap hJi hJJ (fderiv ℝ Φ p b)]
  ring

/-- `d(df) = 0` for a `C²` function. -/
theorem dOne_fderiv {f : W → ℝ} {p : W} (hf : ContDiffAt ℝ 2 f p) (a b : W) :
    dOne (fderiv ℝ f) p a b = 0 := by
  rw [dOne, hf.isSymmSndFDerivAt (by simp) a b, sub_self]

theorem dOne_add {β₁ β₂ : W → W →L[ℝ] ℝ} {p : W} (h₁ : DifferentiableAt ℝ β₁ p)
    (h₂ : DifferentiableAt ℝ β₂ p) (a b : W) :
    dOne (fun q => β₁ q + β₂ q) p a b = dOne β₁ p a b + dOne β₂ p a b := by
  simp only [dOne, fderiv_fun_add h₁ h₂, ContinuousLinearMap.add_apply]
  ring

theorem dOne_sub {β₁ β₂ : W → W →L[ℝ] ℝ} {p : W} (h₁ : DifferentiableAt ℝ β₁ p)
    (h₂ : DifferentiableAt ℝ β₂ p) (a b : W) :
    dOne (fun q => β₁ q - β₂ q) p a b = dOne β₁ p a b - dOne β₂ p a b := by
  simp only [dOne, fderiv_fun_sub h₁ h₂, ContinuousLinearMap.sub_apply]
  ring

/-! ### Cotangent coordinates -/

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {G : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G]

/-! ### `Φ₀`, `ℓ` and the primitive `β` -/

variable (J) (F : E → V) (A : E → G →L[ℝ] V)

/-- `Φ₀(x, ξ) = F x − J (A x ξ)`. -/
def Phi0 (p : E × G) : V := F p.1 - J (A p.1 p.2)

/-- `ℓ(x, ξ) = −½⟨F x, A x ξ⟩`. -/
def ell (p : E × G) : ℝ := -(1 / 2) * ⟪F p.1, A p.1 p.2⟫

/-- `λ_A(x, w)(a) = ⟨A x w, DF a₁⟩`: for `A = ` the metric lift this is `λ_can = ⟨ξ, dx⟩`, and for `A = DF`
it is `g(w, dx)` (the tangent-bundle form). -/
def lamA (p : E × G) : E × G →L[ℝ] ℝ :=
  (innerSL ℝ (A p.1 p.2)).comp ((fderiv ℝ F p.1).comp (ContinuousLinearMap.fst ℝ E G))

/-- `β = Φ₀*λ₀ − λ_A − π*(F*λ₀) − dℓ`. -/
def beta (p : E × G) : E × G →L[ℝ] ℝ :=
  pb1 J (Phi0 J F A) p - lamA F A p - pb1 J (fun q : E × G => F q.1) p - fderiv ℝ (ell F A) p

variable {J F A}

theorem hasFDerivAt_lift {A' : E →L[ℝ] G →L[ℝ] V} {x : E} {ξ : G} (hA : HasFDerivAt A A' x) :
    HasFDerivAt (fun p : E × G => A p.1 p.2)
      ((A x).comp (ContinuousLinearMap.snd ℝ E G) + (A'.comp (ContinuousLinearMap.fst ℝ E G)).flip ξ)
      (x, ξ) :=
  (hA.comp (x, ξ) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt).clm_apply
    (ContinuousLinearMap.snd ℝ E G).hasFDerivAt

theorem fderiv_Phi0_zero {x : E} (hF : DifferentiableAt ℝ F x) (hA : DifferentiableAt ℝ A x)
    (v : E × G) : fderiv ℝ (Phi0 J F A) (x, 0) v = fderiv ℝ F x v.1 - J (A x v.2) := by
  have h1 : HasFDerivAt (fun p : E × G => F p.1) ((fderiv ℝ F x).comp (ContinuousLinearMap.fst ℝ E G))
      (x, (0 : G)) := hF.hasFDerivAt.comp (x, (0 : G)) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt
  have h2 := J.hasFDerivAt.comp (x, (0 : G)) (hasFDerivAt_lift (ξ := 0) hA.hasFDerivAt)
  have h : HasFDerivAt (fun p : E × G => F p.1 - J (A p.1 p.2)) _ (x, (0 : G)) := h1.sub h2
  rw [show Phi0 J F A = fun p => F p.1 - J (A p.1 p.2) from rfl, h.fderiv]
  simp

theorem fderiv_ell_zero {x : E} (hF : DifferentiableAt ℝ F x) (hA : DifferentiableAt ℝ A x)
    (v : E × G) : fderiv ℝ (ell F A) (x, 0) v = -(1 / 2) * ⟪F x, A x v.2⟫ := by
  have h1 : HasFDerivAt (fun p : E × G => F p.1) ((fderiv ℝ F x).comp (ContinuousLinearMap.fst ℝ E G))
      (x, (0 : G)) := hF.hasFDerivAt.comp (x, (0 : G)) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt
  have h2 := hasFDerivAt_lift (ξ := (0 : G)) hA.hasFDerivAt
  have h := (h1.inner ℝ h2).const_mul (-(1 / 2 : ℝ))
  rw [show ell F A = fun p => -(1 / 2) * ⟪F p.1, A p.1 p.2⟫ from rfl, h.fderiv]
  simp

/-- `β` vanishes at every point of the zero section. -/
theorem beta_zero (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) {x : E} (hF : DifferentiableAt ℝ F x)
    (hA : DifferentiableAt ℝ A x) (v : E × G) : beta J F A (x, 0) v = 0 := by
  have hFst : fderiv ℝ (fun q : E × G => F q.1) (x, 0) v = fderiv ℝ F x v.1 := by
    have h1 : HasFDerivAt (fun p : E × G => F p.1)
        ((fderiv ℝ F x).comp (ContinuousLinearMap.fst ℝ E G)) (x, (0 : G)) :=
      hF.hasFDerivAt.comp (x, (0 : G)) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt
    rw [h1.fderiv]
    rfl
  simp only [beta, pb1, lamA, ContinuousLinearMap.sub_apply, ContinuousLinearMap.coe_comp',
    Function.comp_apply, fderiv_Phi0_zero hF hA, fderiv_ell_zero hF hA, hFst, Phi0, map_zero,
    sub_zero, lam0L, ContinuousLinearMap.coe_smul', Pi.smul_apply, innerSL_apply_apply,
    smul_eq_mul, inner_zero_left, ContinuousLinearMap.coe_fst', ContinuousLinearMap.zero_apply,
    neg_zero, inner_sub_right, hJi]
  ring

/-- At the zero section, `Φ₀*ω₀(a, b) = ⟨A a₂, DF b₁⟩ − ⟨A b₂, DF a₁⟩` (`= dλ_A`, `dlamA_zero`), using only
that `F` is Lagrangian and `J A ⊥ A`. -/
theorem pb2_Phi0_zero (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) {x : E}
    (hF : DifferentiableAt ℝ F x) (hA : DifferentiableAt ℝ A x)
    (hL : ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0)
    (hA2 : ∀ ξ η, ⟪J (A x ξ), A x η⟫ = 0) (a b : E × G) :
    ⟪J (fderiv ℝ (Phi0 J F A) (x, 0) a), fderiv ℝ (Phi0 J F A) (x, 0) b⟫ =
      ⟪A x a.2, fderiv ℝ F x b.1⟫ - ⟪A x b.2, fderiv ℝ F x a.1⟫ := by
  rw [fderiv_Phi0_zero hF hA, fderiv_Phi0_zero hF hA, map_sub, hJJ, inner_sub_left, inner_sub_right,
    inner_sub_right, hL, hJi, inner_neg_left, inner_neg_left,
    real_inner_comm (A x b.2) (fderiv ℝ F x a.1),
    real_inner_comm (J (A x b.2)) (A x a.2), hA2]
  ring

theorem differentiableAt_pb1 {Φ : W → V} {p : W} (hΦ : ContDiffAt ℝ 2 Φ p) :
    DifferentiableAt ℝ (pb1 J Φ) p :=
  ((lam0L J).differentiableAt.comp p (hΦ.differentiableAt (by norm_num))).clm_comp
    ((hΦ.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num))

theorem fderiv_fst_apply {Gf : E → V} {x : E} (ξ : G) (hG : DifferentiableAt ℝ Gf x) (a : E × G) :
    fderiv ℝ (fun q : E × G => Gf q.1) (x, ξ) a = fderiv ℝ Gf x a.1 := by
  have h1 : HasFDerivAt (fun q : E × G => Gf q.1)
      ((fderiv ℝ Gf x).comp (ContinuousLinearMap.fst ℝ E G)) (x, ξ) :=
    hG.hasFDerivAt.comp (x, ξ) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt
  rw [h1.fderiv]
  rfl

theorem hasFDerivAt_lamA_aux {x : E} {ξ : G} (hF2 : DifferentiableAt ℝ (fderiv ℝ F) x)
    (hA : DifferentiableAt ℝ A x) : HasFDerivAt (lamA F A) (fderiv ℝ (lamA F A) (x, ξ)) (x, ξ) ∧
      ∀ a b, fderiv ℝ (lamA F A) (x, ξ) a b =
        ⟪A x ξ, fderiv ℝ (fderiv ℝ F) x a.1 b.1⟫ + ⟪A x a.2 + fderiv ℝ A x a.1 ξ, fderiv ℝ F x b.1⟫ := by
  have hc := (innerSL ℝ : V →L⋆[ℝ] V →L[ℝ] ℝ).hasFDerivAt.comp (x, ξ)
    (hasFDerivAt_lift (ξ := ξ) hA.hasFDerivAt)
  have hd := ((ContinuousLinearMap.compL ℝ (E × G) E V).flip
    (ContinuousLinearMap.fst ℝ E G)).hasFDerivAt.comp (x, ξ)
    (hF2.hasFDerivAt.comp (x, ξ) (ContinuousLinearMap.fst ℝ E G).hasFDerivAt)
  have h : HasFDerivAt (lamA F A) _ (x, ξ) := hc.clm_comp hd
  refine ⟨h.differentiableAt.hasFDerivAt, fun a b => ?_⟩
  rw [h.fderiv]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_comp', Function.comp_apply,
    ContinuousLinearMap.compL_apply, ContinuousLinearMap.flip_apply]
  simp [inner_add_left]
  rfl

theorem differentiableAt_lamA {x : E} {ξ : G} (hF2 : DifferentiableAt ℝ (fderiv ℝ F) x)
    (hA : DifferentiableAt ℝ A x) : DifferentiableAt ℝ (lamA F A) (x, ξ) :=
  (hasFDerivAt_lamA_aux hF2 hA).1.differentiableAt

/-- `dλ_A(a, b) = ⟨A a₂, DF b₁⟩ − ⟨A b₂, DF a₁⟩` on the zero section. -/
theorem dOne_lamA_zero {x : E} (hF2 : DifferentiableAt ℝ (fderiv ℝ F) x)
    (hA : DifferentiableAt ℝ A x) (a b : E × G) :
    dOne (lamA F A) (x, 0) a b = ⟪A x a.2, fderiv ℝ F x b.1⟫ - ⟪A x b.2, fderiv ℝ F x a.1⟫ := by
  rw [dOne, (hasFDerivAt_lamA_aux (ξ := 0) hF2 hA).2, (hasFDerivAt_lamA_aux (ξ := 0) hF2 hA).2]
  simp

/-- `dβ = Φ₀*ω₀ − dλ_A` wherever `F` is Lagrangian and `F, A` are `C²`. -/
theorem dOne_beta (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) {x : E} {ξ : G}
    (hF : ContDiffAt ℝ 2 F x) (hA : ContDiffAt ℝ 2 A x)
    (hL : ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0) (a b : E × G) :
    dOne (beta J F A) (x, ξ) a b =
      ⟪J (fderiv ℝ (Phi0 J F A) (x, ξ) a), fderiv ℝ (Phi0 J F A) (x, ξ) b⟫ -
        dOne (lamA F A) (x, ξ) a b := by
  set p : E × G := (x, ξ)
  have hFst : ContDiffAt ℝ 2 (fun q : E × G => F q.1) p := hF.comp p contDiffAt_fst
  have hlift : ContDiffAt ℝ 2 (fun q : E × G => A q.1 q.2) p :=
    (hA.comp p contDiffAt_fst).clm_apply contDiffAt_snd
  have hPhi : ContDiffAt ℝ 2 (Phi0 J F A) p := hFst.sub (J.contDiff.contDiffAt.comp p hlift)
  have hell : ContDiffAt ℝ 2 (ell F A) p := contDiffAt_const.mul (hFst.inner ℝ hlift)
  have hdell : DifferentiableAt ℝ (fderiv ℝ (ell F A)) p :=
    (hell.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have hlam : DifferentiableAt ℝ (lamA F A) p := differentiableAt_lamA
    ((hF.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num))
    (hA.differentiableAt (by norm_num))
  have hfd : fderiv ℝ (beta J F A) p = fderiv ℝ (pb1 J (Phi0 J F A)) p - fderiv ℝ (lamA F A) p -
      fderiv ℝ (pb1 J (fun q : E × G => F q.1)) p - fderiv ℝ (fderiv ℝ (ell F A)) p := by
    have h := (((differentiableAt_pb1 (J := J) hPhi).hasFDerivAt.sub
      hlam.hasFDerivAt).sub (differentiableAt_pb1 (J := J) hFst).hasFDerivAt).sub
      hdell.hasFDerivAt
    exact h.fderiv
  have e1 := dOne_pb1 hJi hJJ hPhi a b
  have e3 := dOne_pb1 hJi hJJ hFst a b
  have e4 := dOne_fderiv hell a b
  rw [fderiv_fst_apply ξ (hF.differentiableAt (by norm_num)),
    fderiv_fst_apply ξ (hF.differentiableAt (by norm_num)), hL] at e3
  simp only [dOne] at e1 e3 e4 ⊢
  rw [hfd]
  simp only [ContinuousLinearMap.sub_apply]
  linarith

end CWFlow
