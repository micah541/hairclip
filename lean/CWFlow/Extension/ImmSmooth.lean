import CWFlow.Extension.ImmCoord
import CWFlow.Extension.GraphCoeff
import CWFlow.Extension.GNorm
import CWFlow.Extension.TensorOps

/-!
# Basic properties of the coordinate tensors of an immersion

For `X : ℝⁿ → V`:

* the metric `gX` is symmetric, and positive definite where `fderiv ℝ X x` is injective, so
  `det g > 0`, `volX > 0` and `g⁻¹` is a two-sided inverse and symmetric;
* all tensors `dX, ddX, gX, ginvZ, volX, gamX, cX, SX, eDens, aNormSq` at `x` depend only on the
  germ of `X` at `x`;
* they are smooth on an open set where `X` is smooth (and the metric is nondegenerate);
* the recursion `SX (k+1)` in coordinates.
-/

open scoped ContDiff Matrix RealInnerProductSpace Topology

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X X' : EuclideanSpace ℝ (Fin n) → V} {x : EuclideanSpace ℝ (Fin n)}
  {Ω : Set (EuclideanSpace ℝ (Fin n))}

/-! ### Part 1: the metric -/

theorem gX_symm (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) (a b : Fin n) :
    gX X x a b = gX X x b a := by
  unfold gX
  exact real_inner_comm _ _

theorem ims_quad (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n))
    (v : Fin n → ℝ) :
    v ⬝ᵥ (gX X x *ᵥ v) = ‖fderiv ℝ X x (∑ a, v a • eb a)‖ ^ 2 := by
  have h1 : fderiv ℝ X x (∑ a, v a • eb a) = ∑ a, v a • dX X x a := by
    simp [dX, map_sum, map_smul]
  rw [h1, ← real_inner_self_eq_norm_sq, sum_inner]
  simp only [dotProduct, Matrix.mulVec, gX, inner_sum, inner_smul_left, inner_smul_right,
    Finset.mul_sum, conj_trivial]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
  ring

theorem ims_sum_ne_zero {v : Fin n → ℝ} (hv : v ≠ 0) : (∑ a, v a • eb a) ≠ 0 := by
  intro h
  apply hv
  funext b
  have := congrArg (fun w : EuclideanSpace ℝ (Fin n) => w b) h
  simpa [Pi.single_apply] using this

theorem gX_posDef (hinj : Function.Injective (fderiv ℝ X x)) : (gX X x).PosDef := by
  rw [Matrix.posDef_iff_dotProduct_mulVec]
  refine ⟨?_, fun v hv => ?_⟩
  · ext a b
    simp [Matrix.conjTranspose_apply, gX_symm X x a b]
  · rw [star_trivial, ims_quad]
    have hne : fderiv ℝ X x (∑ a, v a • eb a) ≠ 0 := by
      intro h
      exact ims_sum_ne_zero hv (hinj (by rw [h, map_zero]))
    positivity

theorem gX_det_pos (hinj : Function.Injective (fderiv ℝ X x)) : 0 < (gX X x).det :=
  (gX_posDef hinj).det_pos

theorem volX_pos (hinj : Function.Injective (fderiv ℝ X x)) : 0 < volX X x :=
  Real.sqrt_pos.2 (gX_det_pos hinj)

theorem gX_mul_inv (hinj : Function.Injective (fderiv ℝ X x)) : gX X x * (gX X x)⁻¹ = 1 :=
  Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.2 (gX_det_pos hinj).ne')

theorem gX_inv_mul (hinj : Function.Injective (fderiv ℝ X x)) : (gX X x)⁻¹ * gX X x = 1 :=
  Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.2 (gX_det_pos hinj).ne')

theorem gX_inv_symm (hinj : Function.Injective (fderiv ℝ X x)) (a b : Fin n) :
    (gX X x)⁻¹ a b = (gX X x)⁻¹ b a := by
  have h := (gX_posDef hinj).inv.isHermitian
  have := congrFun (congrFun h.eq b) a
  simpa [Matrix.conjTranspose_apply] using this

/-! ### Part 2: dependence on the germ -/

theorem dX_congr (h : X =ᶠ[𝓝 x] X') : dX X x = dX X' x := by
  funext a
  unfold dX
  rw [h.fderiv_eq]

theorem ddX_congr (h : X =ᶠ[𝓝 x] X') : ddX X x = ddX X' x := by
  funext a b
  unfold ddX
  rw [h.fderiv.fderiv_eq]

theorem gX_congr (h : X =ᶠ[𝓝 x] X') : gX X x = gX X' x := by
  funext a b
  unfold gX
  rw [dX_congr h]

theorem ginvZ_congr (h : X =ᶠ[𝓝 x] X') : ginvZ X x = ginvZ X' x := by
  unfold ginvZ
  rw [gX_congr h]

theorem volX_congr (h : X =ᶠ[𝓝 x] X') : volX X x = volX X' x := by
  unfold volX
  rw [gX_congr h]

theorem gamX_congr (h : X =ᶠ[𝓝 x] X') : gamX X x = gamX X' x := by
  unfold gamX
  rw [gX_congr h, ddX_congr h, dX_congr h]

theorem cX_congr (h : X =ᶠ[𝓝 x] X') : cX J X x = cX J X' x := by
  unfold cX
  rw [ddX_congr h, dX_congr h]

theorem SX_eventuallyEq (h : X =ᶠ[𝓝 x] X') (k : ℕ) : SX J X k =ᶠ[𝓝 x] SX J X' k := by
  induction k with
  | zero =>
    rw [SX_zero, SX_zero]
    filter_upwards [Filter.Eventually.eventually_nhds h] with y hy
    exact cX_congr hy
  | succ k ih =>
    have h0 : fderiv ℝ (SX J X k) =ᶠ[𝓝 x] fderiv ℝ (SX J X' k) := ih.fderiv
    filter_upwards [h0, ih, Filter.Eventually.eventually_nhds h] with y h1 h2 h3
    rw [SX_succ, SX_succ, h1, h2, gamX_congr h3]

theorem SX_congr (h : X =ᶠ[𝓝 x] X') (k : ℕ) : SX J X k x = SX J X' k x :=
  (SX_eventuallyEq h k).eq_of_nhds

theorem eDens_congr (h : X =ᶠ[𝓝 x] X') (k : ℕ) : eDens J X k x = eDens J X' k x := by
  unfold eDens
  rw [ginvZ_congr h, SX_congr h k, volX_congr h]

theorem aNormSq_congr (h : X =ᶠ[𝓝 x] X') : aNormSq J X x = aNormSq J X' x := by
  unfold aNormSq
  rw [ginvZ_congr h, cX_congr h]

/-! ### Part 3: smoothness -/

theorem ims_contDiffOn_fderiv {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    {Y : EuclideanSpace ℝ (Fin n) → F} (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ Y Ω) :
    ContDiffOn ℝ ∞ (fderiv ℝ Y) Ω :=
  hX.fderiv_of_isOpen hΩ (by simp)

theorem contDiffOn_dX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (a : Fin n) :
    ContDiffOn ℝ ∞ (fun x => dX X x a) Ω :=
  (ims_contDiffOn_fderiv hΩ hX).clm_apply (contDiffOn_const (c := eb a))

theorem contDiffOn_ddX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun x => ddX X x a b) Ω :=
  ((ims_contDiffOn_fderiv hΩ (ims_contDiffOn_fderiv hΩ hX)).clm_apply
    (contDiffOn_const (c := eb a))).clm_apply (contDiffOn_const (c := eb b))

theorem contDiffOn_gX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun x => gX X x a b) Ω :=
  ContDiffOn.inner ℝ (contDiffOn_dX hΩ hX a) (contDiffOn_dX hΩ hX b)

theorem ims_contDiffOn_det {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] {s : Set E}
    {f : E → Matrix (Fin n) (Fin n) ℝ} (hf : ∀ i j, ContDiffOn ℝ ∞ (fun x => f x i j) s) :
    ContDiffOn ℝ ∞ (fun x => (f x).det) s := by
  have hg : ContDiff ℝ ∞ (fun M : Fin n → Fin n → ℝ => (Matrix.of M).det) :=
    contDiff_det_of (f := fun M : Fin n → Fin n → ℝ => Matrix.of M) fun i j =>
      (contDiff_apply ℝ ℝ j).comp (contDiff_apply ℝ (Fin n → ℝ) i)
  exact hg.comp_contDiffOn (f := fun x i j => f x i j)
    (contDiffOn_pi.2 fun i => contDiffOn_pi.2 fun j => hf i j)

theorem ims_contDiffOn_adjugate {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] {s : Set E}
    {f : E → Matrix (Fin n) (Fin n) ℝ} (hf : ∀ i j, ContDiffOn ℝ ∞ (fun x => f x i j) s)
    (i j : Fin n) : ContDiffOn ℝ ∞ (fun x => (f x).adjugate i j) s := by
  have hg : ContDiff ℝ ∞ (fun M : Fin n → Fin n → ℝ => (Matrix.of M).adjugate i j) :=
    contDiff_adjugate_of (f := fun M : Fin n → Fin n → ℝ => Matrix.of M) (fun i j =>
      (contDiff_apply ℝ ℝ j).comp (contDiff_apply ℝ (Fin n → ℝ) i)) i j
  exact hg.comp_contDiffOn (f := fun x i j => f x i j)
    (contDiffOn_pi.2 fun i => contDiffOn_pi.2 fun j => hf i j)

theorem contDiffOn_gX_inv (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun x => (gX X x)⁻¹ a b) Ω := by
  have : ∀ x, (gX X x)⁻¹ a b = ((gX X x).det)⁻¹ * (gX X x).adjugate a b := by
    intro x
    rw [Matrix.inv_def, Ring.inverse_eq_inv']
    rfl
  simp only [this]
  exact ((ims_contDiffOn_det (contDiffOn_gX hΩ hX)).inv fun x hx => (gX_det_pos (hinj x hx)).ne').mul
    (ims_contDiffOn_adjugate (contDiffOn_gX hΩ hX) a b)

theorem contDiffOn_ginvZ (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) : ContDiffOn ℝ ∞ (ginvZ X) Ω :=
  contDiffOn_piLp' (𝕜 := ℝ) (p := 2) fun I => contDiffOn_gX_inv hΩ hX hinj (I 0) (I 1)

theorem contDiffOn_volX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) : ContDiffOn ℝ ∞ (volX X) Ω :=
  (ims_contDiffOn_det (contDiffOn_gX hΩ hX)).sqrt fun x hx => (gX_det_pos (hinj x hx)).ne'

theorem contDiffOn_gamX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) : ContDiffOn ℝ ∞ (gamX X) Ω := by
  refine contDiffOn_piLp' (𝕜 := ℝ) (p := 2) fun I => ?_
  show ContDiffOn ℝ ∞ (fun x => ∑ d, (gX X x)⁻¹ (I 2) d * ⟪ddX X x (I 0) (I 1), dX X x d⟫) Ω
  refine ContDiffOn.sum (s := (Finset.univ : Finset (Fin n))) fun d _ => ?_
  exact (contDiffOn_gX_inv hΩ hX hinj (I 2) d).mul
    (ContDiffOn.inner ℝ (contDiffOn_ddX hΩ hX (I 0) (I 1)) (contDiffOn_dX hΩ hX d))

theorem contDiffOn_cX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) : ContDiffOn ℝ ∞ (cX J X) Ω := by
  refine contDiffOn_piLp' (𝕜 := ℝ) (p := 2) fun I => ?_
  exact ContDiffOn.inner ℝ (contDiffOn_ddX hΩ hX (I 0) (I 1))
    (J.contDiff.comp_contDiffOn (contDiffOn_dX hΩ hX (I 2)))

theorem contDiffOn_SX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (k : ℕ) :
    ContDiffOn ℝ ∞ (SX J X k) Ω := by
  induction k with
  | zero => exact contDiffOn_cX hΩ hX
  | succ k ih =>
    have h1 : ContDiffOn ℝ ∞ (fun x => curryZ (fderiv ℝ (SX J X k) x)) Ω :=
      curryZ.contDiff.comp_contDiffOn (ih.fderiv_of_isOpen hΩ (by simp))
    have h2 : ContDiffOn ℝ ∞ (fun x => actG (gamX X x) (SX J X k x)) Ω :=
      (actG.contDiff.comp_contDiffOn (contDiffOn_gamX hΩ hX hinj)).clm_apply ih
    exact h1.sub h2

theorem contDiffOn_eDens (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (k : ℕ) :
    ContDiffOn ℝ ∞ (eDens J X k) Ω :=
  by
  have h := (contDiff_gNormSq (n := n) (m := k + 3)).comp_contDiffOn
    ((contDiffOn_ginvZ hΩ hX hinj).prodMk (contDiffOn_SX (J := J) hΩ hX hinj k))
  exact h.mul (contDiffOn_volX hΩ hX hinj)

theorem contDiffOn_aNormSq (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) : ContDiffOn ℝ ∞ (aNormSq J X) Ω :=
  by
  have h := (contDiff_gNormSq (n := n) (m := 3)).comp_contDiffOn
    ((contDiffOn_ginvZ hΩ hX hinj).prodMk (contDiffOn_cX (J := J) hΩ hX))
  exact h

/-! ### Part 4: the recursion in coordinates -/

theorem ims_fderiv_piLp_apply {E ι : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [Fintype ι]
    {f : E → PiLp 2 (fun _ : ι => ℝ)} {x : E} (hf : DifferentiableAt ℝ f x) (v : E) (I : ι) :
    fderiv ℝ f x v I = fderiv ℝ (fun y => f y I) x v := by
  have h := (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => ℝ) I).hasFDerivAt.comp x hf.hasFDerivAt
  have h2 : fderiv ℝ (fun y => f y I) x =
      (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => ℝ) I).comp (fderiv ℝ f x) := h.fderiv
  rw [h2]
  rfl

theorem SX_succ_apply (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (k : ℕ) (hx : x ∈ Ω)
    (I : Fin (k + 3 + 1) → Fin n) :
    SX J X (k + 1) x I =
      fderiv ℝ (fun y => SX J X k y (Fin.tail I)) x (eb (I 0))
      - ∑ i : Fin (k + 3), ∑ c : Fin n,
          gamX X x ![I 0, Fin.tail I i, c] * SX J X k x (Function.update (Fin.tail I) i c) := by
  have hd : DifferentiableAt ℝ (SX J X k) x :=
    ((contDiffOn_SX (J := J) hΩ hX hinj k).contDiffAt (hΩ.mem_nhds hx)).differentiableAt (by simp)
  rw [SX_succ, PiLp.sub_apply, curryZ_apply, actG_apply, ims_fderiv_piLp_apply hd]

/-- info: 'CWFlow.Extension.gX_posDef' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gX_posDef

/-- info: 'CWFlow.Extension.SX_congr' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms SX_congr

/-- info: 'CWFlow.Extension.contDiffOn_SX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_SX

/-- info: 'CWFlow.Extension.contDiffOn_eDens' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_eDens

/-- info: 'CWFlow.Extension.SX_succ_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms SX_succ_apply

end CWFlow.Extension
