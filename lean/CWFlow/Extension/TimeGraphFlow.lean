import CWFlow.Extension.TimeGraph
import CWFlow.Extension.PotReparam
import CWFlow.Extension.ImmSmooth
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# The time graph function solves `∂ₜφ = div(JH)`

Part 1: for a smooth family `Y_t` of immersions, the coordinate tensors of the slices (`potX`, `velX`)
are jointly smooth in `(t, y)`. The key device is the slice derivative `tf_sd F q = D_y F(q.1, ·) (q.2)`,
which for `F` smooth on an open set is again smooth (it is `D F ∘ inr`); iterating gives `dX`, `ddX`,
`SX` and everything built from them.

Part 2 (`exists_time_graph_flow`): the graph functions of `exists_time_graph` can be normalized (by adding
a function of `t` alone, the primitive of `potX Y_t (ξ(t, 0))`) so that `∂ₜφ = potX Y_t ∘ ξ_t` exactly. The
key step is `tf_step1`: differentiating `Y(t, ξ(t, x)) = y₀ + U x + J U ∇φ₀` in `t` and pairing with the
tangent vectors gives `∂ₜ ∇φ₀ = ∇(potX ∘ ξ)`; `tf_dt_phi` then integrates in `x`.
-/

open Set Metric
open CWFlow.Interp
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-! ### The slice derivative -/

/-- The derivative in the spatial variable of the slice `F (q.1, ·)` at `q.2`. -/
noncomputable def tf_sd {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] (F : ℝ × 𝔼 → W)
    (q : ℝ × 𝔼) : 𝔼 →L[ℝ] W :=
  fderiv ℝ (fun y => F (q.1, y)) q.2

theorem tf_sd_eq {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {F : ℝ × 𝔼 → W}
    {q : ℝ × 𝔼} (hF : DifferentiableAt ℝ F q) :
    tf_sd F q = (fderiv ℝ F q).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼) := by
  have h1 : HasFDerivAt (fun y : 𝔼 => ((q.1, y) : ℝ × 𝔼)) (ContinuousLinearMap.inr ℝ ℝ 𝔼) q.2 :=
    hasFDerivAt_prodMk_right q.1 q.2
  have h2 : HasFDerivAt F (fderiv ℝ F q) (q.1, q.2) := hF.hasFDerivAt
  exact (h2.comp q.2 h1).fderiv

theorem tf_sd_smooth {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {F : ℝ × 𝔼 → W}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hF : ContDiffOn ℝ ∞ F O) :
    ContDiffOn ℝ ∞ (tf_sd F) O := by
  have h1 : ContDiffOn ℝ ∞ (fun q => (fderiv ℝ F q).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼)) O :=
    (hF.fderiv_of_isOpen hO (by simp)).clm_comp
      (contDiffOn_const (c := ContinuousLinearMap.inr ℝ ℝ 𝔼))
  exact h1.congr fun q hq =>
    tf_sd_eq ((hF.contDiffAt (hO.mem_nhds hq)).differentiableAt (by simp))

/-! ### The tensors of the slices -/

theorem tf_dX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {Y : ℝ × 𝔼 → V}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O) (a : Fin n) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => dX (fun y => Y (q.1, y)) q.2 a) O :=
  (tf_sd_smooth hO hY).clm_apply (contDiffOn_const (c := eb a))

theorem tf_ddX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {Y : ℝ × 𝔼 → V}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ddX (fun y => Y (q.1, y)) q.2 a b) O :=
  ((tf_sd_smooth hO (tf_sd_smooth hO hY)).clm_apply (contDiffOn_const (c := eb a))).clm_apply
    (contDiffOn_const (c := eb b))

theorem tf_gX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {Y : ℝ × 𝔼 → V}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => gX (fun y => Y (q.1, y)) q.2 a b) O :=
  ContDiffOn.inner ℝ (tf_dX hO hY a) (tf_dX hO hY b)

theorem tf_gX_inv {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {Y : ℝ × 𝔼 → V}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) (a b : Fin n) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => (gX (fun y => Y (q.1, y)) q.2)⁻¹ a b) O := by
  have : ∀ q : ℝ × 𝔼, (gX (fun y => Y (q.1, y)) q.2)⁻¹ a b =
      ((gX (fun y => Y (q.1, y)) q.2).det)⁻¹ * (gX (fun y => Y (q.1, y)) q.2).adjugate a b := by
    intro q
    rw [Matrix.inv_def, Ring.inverse_eq_inv']
    rfl
  simp only [this]
  exact ((ims_contDiffOn_det (fun i j => tf_gX hO hY i j)).inv fun q hq =>
    (gX_det_pos (hinj q hq)).ne').mul
    (ims_contDiffOn_adjugate (fun i j => tf_gX hO hY i j) a b)

theorem tf_gamX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {Y : ℝ × 𝔼 → V}
    {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => gamX (fun y => Y (q.1, y)) q.2) O := by
  refine contDiffOn_piLp' (𝕜 := ℝ) (p := 2) fun I => ?_
  show ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ∑ d, (gX (fun y => Y (q.1, y)) q.2)⁻¹ (I 2) d *
    ⟪ddX (fun y => Y (q.1, y)) q.2 (I 0) (I 1), dX (fun y => Y (q.1, y)) q.2 d⟫) O
  refine ContDiffOn.sum (s := (Finset.univ : Finset (Fin n))) fun d _ => ?_
  exact (tf_gX_inv hO hY hinj (I 2) d).mul
    (ContDiffOn.inner ℝ (tf_ddX hO hY (I 0) (I 1)) (tf_dX hO hY d))

theorem tf_cX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    {Y : ℝ × 𝔼 → V} {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => cX J (fun y => Y (q.1, y)) q.2) O := by
  refine contDiffOn_piLp' (𝕜 := ℝ) (p := 2) fun I => ?_
  exact ContDiffOn.inner ℝ (tf_ddX hO hY (I 0) (I 1))
    (J.contDiff.comp_contDiffOn (tf_dX hO hY (I 2)))

theorem tf_SX {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    {Y : ℝ × 𝔼 → V} {O : Set (ℝ × 𝔼)} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) (k : ℕ) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => SX J (fun y => Y (q.1, y)) k q.2) O := by
  induction k with
  | zero => exact tf_cX hO hY
  | succ k ih =>
    have h1 : ContDiffOn ℝ ∞
        (fun q : ℝ × 𝔼 => curryZ (fderiv ℝ (SX J (fun y => Y (q.1, y)) k) q.2)) O :=
      curryZ.contDiff.comp_contDiffOn (tf_sd_smooth hO ih)
    have h2 : ContDiffOn ℝ ∞
        (fun q : ℝ × 𝔼 => actG (gamX (fun y => Y (q.1, y)) q.2)
          (SX J (fun y => Y (q.1, y)) k q.2)) O :=
      (actG.contDiff.comp_contDiffOn (tf_gamX hO hY hinj)).clm_apply ih
    exact h1.sub h2

theorem contDiffOn_potX_family {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))}
    (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => potX J (fun y => Y (q.1, y)) q.2) O := by
  unfold potX
  refine ContDiffOn.neg ?_
  refine ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => ContDiffOn.sum fun c _ =>
    ContDiffOn.sum fun d _ => ?_
  have h1 : ∀ I : Fin 2 → Fin n,
      ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ginvZ (fun y => Y (q.1, y)) q.2 I) O := fun I =>
    tf_gX_inv hO hY hinj (I 0) (I 1)
  have h2 : ContDiffOn ℝ ∞
      (fun q : ℝ × 𝔼 => SX J (fun y => Y (q.1, y)) 1 q.2 ![d, a, b, c]) O :=
    (contDiff_coordG (![d, a, b, c] : Fin 4 → Fin n)).comp_contDiffOn (tf_SX hO hY hinj 1)
  exact ((h1 ![c, d]).mul (h1 ![a, b])).mul h2

theorem contDiffOn_velX_family {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))}
    (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => velX J (fun y => Y (q.1, y)) q.2) O := by
  have hp := contDiffOn_potX_family (J := J) hO hY hinj
  have hdp : ∀ b : Fin n, ContDiffOn ℝ ∞
      (fun q : ℝ × 𝔼 => fderiv ℝ (potX J (fun y => Y (q.1, y))) q.2 (eb b)) O := fun b =>
    (tf_sd_smooth hO hp).clm_apply (contDiffOn_const (c := eb b))
  unfold velX
  refine ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => ?_
  exact ((tf_gX_inv hO hY hinj a b).mul (hdp b)).smul
    (J.contDiff.comp_contDiffOn (tf_dX hO hY a))

/-! ### Part 2: the equation -/

theorem tf_hasDerivAt_slice_t {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {F : ℝ × 𝔼 → W}
    {D : ℝ × 𝔼 →L[ℝ] W} {t : ℝ} {x : 𝔼} (h : HasFDerivAt F D (t, x)) :
    HasDerivAt (fun s => F (s, x)) (D (1, 0)) t :=
  h.comp_hasDerivAt t ((hasDerivAt_id t).prodMk (hasDerivAt_const t x))

theorem tf_Jvel_inner {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    (hJJ : ∀ v : V, J (J v) = -v) {X : 𝔼 → V} {y : 𝔼}
    (hinj : Function.Injective (fderiv ℝ X y)) (w : 𝔼) :
    ⟪J (velX J X y), fderiv ℝ X y w⟫ = - fderiv ℝ (potX J X) y w := by
  have hw : w = ∑ c, w c • eb c := by rw [gb_sum_eb]
  have hDX : fderiv ℝ X y w = ∑ c, w c • dX X y c := by
    conv_lhs => rw [hw]
    simp [map_sum, map_smul, dX]
  have hDp : fderiv ℝ (potX J X) y w = ∑ c, w c * fderiv ℝ (potX J X) y (eb c) := by
    conv_lhs => rw [hw]
    simp [map_sum, map_smul]
  have hJv : J (velX J X y) =
      -∑ a, ∑ b, ((gX X y)⁻¹ a b * fderiv ℝ (potX J X) y (eb b)) • dX X y a := by
    unfold velX
    simp [map_sum, map_smul, hJJ, Finset.sum_neg_distrib]
  have key : ∀ b c, ∑ a, (gX X y)⁻¹ a b * gX X y a c = if b = c then 1 else 0 := by
    intro b c
    have h := congrFun (congrFun (gX_inv_mul hinj) b) c
    rw [Matrix.mul_apply, Matrix.one_apply] at h
    rw [← h]
    exact Finset.sum_congr rfl fun a _ => by rw [gX_inv_symm hinj a b]
  rw [hJv, hDX, hDp, inner_neg_left, neg_inj]
  simp only [sum_inner, inner_sum, inner_smul_left, inner_smul_right, conj_trivial]
  have hg : ∀ a c, ⟪dX X y a, dX X y c⟫ = gX X y a c := fun a c => rfl
  simp only [hg]
  refine Finset.sum_congr rfl fun c _ => ?_
  congr 1
  calc ∑ a, ∑ b, (gX X y)⁻¹ a b * fderiv ℝ (potX J X) y (eb b) * gX X y a c
      = ∑ b, fderiv ℝ (potX J X) y (eb b) * ∑ a, (gX X y)⁻¹ a b * gX X y a c := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun b _ => ?_
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun a _ => by ring
    _ = fderiv ℝ (potX J X) y (eb c) := by
        simp only [key, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]

theorem tf_step1 {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : 𝔼 →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y₀ : V}
    {a b R : ℝ} {Y : ℝ × 𝔼 → V}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a b ×ˢ ball (0 : 𝔼) R))
    (hinj : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R,
      Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    {I : Set ℝ} (hI : IsOpen I) (hIab : I ⊆ Ioo a b) {r : ℝ}
    {ξ : ℝ × 𝔼 → 𝔼} {φ₀ : ℝ × 𝔼 → ℝ}
    (hξ : ContDiffOn ℝ ∞ ξ (I ×ˢ ball (0 : 𝔼) r))
    (hξR : ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r, ξ (t, x) ∈ ball (0 : 𝔼) R)
    (hφ : ContDiffOn ℝ ∞ φ₀ (I ×ˢ ball (0 : 𝔼) r))
    (hgr : ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r,
      Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ₀ (t, x')) x)
    {t : ℝ} (ht : t ∈ I) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) r) (i : Fin n) :
    fderiv ℝ (fderiv ℝ φ₀) (t, x) (1, 0) (0, eb i) =
      fderiv ℝ (fun z => potX J (fun y => Y (t, y)) (ξ (t, z))) x (eb i) := by
  classical
  have hSo : IsOpen (I ×ˢ ball (0 : 𝔼) r) := hI.prod isOpen_ball
  have hSn : I ×ˢ ball (0 : 𝔼) r ∈ 𝓝 (t, x) := hSo.mem_nhds ⟨ht, hx⟩
  have hyR : ξ (t, x) ∈ ball (0 : 𝔼) R := hξR t ht x hx
  have htab : t ∈ Ioo a b := hIab ht
  have hO : IsOpen (Ioo a b ×ˢ ball (0 : 𝔼) R) := isOpen_Ioo.prod isOpen_ball
  have hmem : (t, ξ (t, x)) ∈ Ioo a b ×ˢ ball (0 : 𝔼) R := ⟨htab, hyR⟩
  have hξd : DifferentiableAt ℝ ξ (t, x) := (hξ.contDiffAt hSn).differentiableAt (by simp)
  have hYd : DifferentiableAt ℝ Y (t, ξ (t, x)) :=
    (hY.contDiffAt (hO.mem_nhds hmem)).differentiableAt (by simp)
  have hφC : ContDiffOn ℝ ∞ (fderiv ℝ φ₀) (I ×ˢ ball (0 : 𝔼) r) :=
    hφ.fderiv_of_isOpen hSo (by simp)
  have hHd : DifferentiableAt ℝ (fderiv ℝ φ₀) (t, x) :=
    (hφC.contDiffAt hSn).differentiableAt (by simp)
  let Λ : (ℝ × 𝔼 →L[ℝ] ℝ) →L[ℝ] 𝔼 :=
    (EuclideanSpace.equiv (Fin n) ℝ).symm.toContinuousLinearMap.comp
      (ContinuousLinearMap.pi fun c => ContinuousLinearMap.apply ℝ ℝ ((0 : ℝ), eb c))
  have hΛ : ∀ h c, Λ h c = h (0, eb c) := fun h c => rfl
  let w : ℝ × 𝔼 → 𝔼 := fun q => Λ (fderiv ℝ φ₀ q)
  have hwgrad : ∀ q ∈ I ×ˢ ball (0 : 𝔼) r, gradE (fun x' => φ₀ (q.1, x')) q.2 = w q := by
    intro q hq
    have hd : DifferentiableAt ℝ φ₀ q :=
      (hφ.contDiffAt (hSo.mem_nhds hq)).differentiableAt (by simp)
    have e := tf_sd_eq hd
    ext c
    have := congrArg (fun L => L (eb c)) e
    simpa [tf_sd, gradE, hΛ, w] using this
  have hG : ∀ q ∈ I ×ˢ ball (0 : 𝔼) r,
      Y (q.1, ξ q) = y₀ + U q.2 + J (U (w q)) := by
    intro q hq
    have := hgr q.1 hq.1 q.2 hq.2
    rw [show ξ q = ξ (q.1, q.2) from rfl, this]
    simp [graphX, hwgrad q hq]
  -- two expressions for the derivative of `q ↦ Y (q.1, ξ q)`
  have hP : HasFDerivAt (fun q : ℝ × 𝔼 => ((q.1, ξ q) : ℝ × 𝔼))
      ((ContinuousLinearMap.fst ℝ ℝ 𝔼).prod (fderiv ℝ ξ (t, x))) (t, x) :=
    (hasFDerivAt_fst).prodMk hξd.hasFDerivAt
  have h1 : HasFDerivAt (fun q : ℝ × 𝔼 => Y (q.1, ξ q))
      ((fderiv ℝ Y (t, ξ (t, x))).comp
        ((ContinuousLinearMap.fst ℝ ℝ 𝔼).prod (fderiv ℝ ξ (t, x)))) (t, x) :=
    hYd.hasFDerivAt.comp (t, x) hP
  have hwd : HasFDerivAt w (Λ.comp (fderiv ℝ (fderiv ℝ φ₀) (t, x))) (t, x) :=
    Λ.hasFDerivAt.comp (t, x) hHd.hasFDerivAt
  have h2 : HasFDerivAt (fun q : ℝ × 𝔼 => y₀ + U q.2 + J (U (w q)))
      ((U.toContinuousLinearMap.comp (ContinuousLinearMap.snd ℝ ℝ 𝔼)) +
        ((J.comp U.toContinuousLinearMap).comp
          (Λ.comp (fderiv ℝ (fderiv ℝ φ₀) (t, x))))) (t, x) :=
    ((U.toContinuousLinearMap.comp (ContinuousLinearMap.snd ℝ ℝ 𝔼)).hasFDerivAt.const_add y₀).add
      ((J.comp U.toContinuousLinearMap).hasFDerivAt.comp (t, x) hwd)
  have h3 := h2.congr_of_eventuallyEq (f₁ := fun q : ℝ × 𝔼 => Y (q.1, ξ q))
    (by filter_upwards [hSn] with q hq using hG q hq)
  have hD := h1.unique h3
  have hv : ∀ v : ℝ × 𝔼, fderiv ℝ Y (t, ξ (t, x)) (v.1, fderiv ℝ ξ (t, x) v) =
      U v.2 + J (U (Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) v))) := by
    intro v
    have := congrArg (fun L => L v) hD
    simpa using this
  have hDy : ∀ z : 𝔼, fderiv ℝ Y (t, ξ (t, x)) (0, z) =
      fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) z := by
    intro z
    have := congrArg (fun L => L z) (tf_sd_eq hYd)
    simpa [tf_sd] using this.symm
  have hvY : fderiv ℝ Y (t, ξ (t, x)) (1, 0) = velX J (fun y' => Y (t, y')) (ξ (t, x)) :=
    (tf_hasDerivAt_slice_t hYd.hasFDerivAt).unique (hvel t htab _ hyR)
  have hE1 : J (U (Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) (1, 0)))) =
      velX J (fun y' => Y (t, y')) (ξ (t, x)) +
        fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (1, 0)) := by
    have := hv (1, 0)
    simp only [map_zero, zero_add] at this
    have e : ((1 : ℝ), fderiv ℝ ξ (t, x) (1, 0)) = ((1 : ℝ), (0 : 𝔼)) + (0, fderiv ℝ ξ (t, x) (1, 0)) := by
      simp
    rw [e, map_add, hvY, hDy] at this
    exact this.symm
  have hE2 : ∀ i : Fin n, fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i)) =
      U (eb i) + J (U (Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) (0, eb i)))) := by
    intro i
    have := hv (0, eb i)
    rw [← hDy]
    exact this
  have hpd : DifferentiableAt ℝ (potX J (fun y => Y (t, y))) (ξ (t, x)) :=
    ((vv_contDiffOn_potX (J := J) isOpen_ball (tg_slice hY htab) (hinj t htab)).contDiffAt
      (isOpen_ball.mem_nhds hyR)).differentiableAt (by simp)
  have hξtd : DifferentiableAt ℝ (fun z => ξ (t, z)) x :=
    ((tg_slice hξ ht).contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
  have hξi : fderiv ℝ (fun z => ξ (t, z)) x (eb i) = fderiv ℝ ξ (t, x) (0, eb i) := by
    have := congrArg (fun L => L (eb i)) (tf_sd_eq hξd)
    simpa [tf_sd] using this
  have hL : ⟪J (J (U (Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) (1, 0))))),
      fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i))⟫ =
      -(Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) (1, 0)) i) := by
    rw [hJJ, inner_neg_left, hE2 i, inner_add_right, gb_UU, gb_UJU hU]
    simp [eb, EuclideanSpace.inner_single_right]
  have hR : ⟪J (velX J (fun y' => Y (t, y')) (ξ (t, x)) +
        fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (1, 0))),
      fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i))⟫ =
      - fderiv ℝ (potX J (fun y => Y (t, y))) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i)) := by
    rw [map_add, inner_add_left, hlag t htab _ hyR, add_zero]
    exact tf_Jvel_inner hJJ (hinj t htab _ hyR) _
  have h6 : ⟪J (J (U (Λ (fderiv ℝ (fderiv ℝ φ₀) (t, x) (1, 0))))),
      fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i))⟫ =
      ⟪J (velX J (fun y' => Y (t, y')) (ξ (t, x)) +
        fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (1, 0))),
      fderiv ℝ (fun y' => Y (t, y')) (ξ (t, x)) (fderiv ℝ ξ (t, x) (0, eb i))⟫ := by
    rw [hE1]
  rw [hL, hR, neg_inj, ← hξi] at h6
  have hc : fderiv ℝ (fun z => potX J (fun y => Y (t, y)) (ξ (t, z))) x (eb i) =
      fderiv ℝ (potX J (fun y => Y (t, y))) (ξ (t, x)) (fderiv ℝ (fun z => ξ (t, z)) x (eb i)) := by
    rw [show (fun z => potX J (fun y => Y (t, y)) (ξ (t, z))) =
      potX J (fun y => Y (t, y)) ∘ (fun z => ξ (t, z)) from rfl, fderiv_comp x hpd hξtd]
    rfl
  rw [hc, ← h6]
  rfl

theorem tf_ext_eb {L M : 𝔼 →L[ℝ] ℝ} (h : ∀ i : Fin n, L (eb i) = M (eb i)) : L = M := by
  ext v
  have hv : v = ∑ c, v c • eb c := by rw [gb_sum_eb]
  have e1 : L v = ∑ c, v c * L (eb c) := by
    conv_lhs => rw [hv]
    simp [map_sum, map_smul]
  have e2 : M v = ∑ c, v c * M (eb c) := by
    conv_lhs => rw [hv]
    simp [map_sum, map_smul]
  rw [e1, e2]
  exact Finset.sum_congr rfl fun c _ => by rw [h c]

theorem tf_dt_phi {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : 𝔼 →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y₀ : V}
    {a b R : ℝ} {Y : ℝ × 𝔼 → V}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a b ×ˢ ball (0 : 𝔼) R))
    (hinj : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R,
      Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) R,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    {I : Set ℝ} (hI : IsOpen I) (hIab : I ⊆ Ioo a b) {r : ℝ}
    {ξ : ℝ × 𝔼 → 𝔼} {φ₀ : ℝ × 𝔼 → ℝ}
    (hξ : ContDiffOn ℝ ∞ ξ (I ×ˢ ball (0 : 𝔼) r))
    (hξR : ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r, ξ (t, x) ∈ ball (0 : 𝔼) R)
    (hφ : ContDiffOn ℝ ∞ φ₀ (I ×ˢ ball (0 : 𝔼) r))
    (hgr : ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r,
      Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ₀ (t, x')) x)
    (h0 : ∀ t ∈ I, φ₀ (t, 0) = 0) (hr : 0 < r) :
    ∀ t ∈ I, ∀ x ∈ ball (0 : 𝔼) r,
      HasDerivAt (fun s => φ₀ (s, x))
        (potX J (fun y => Y (t, y)) (ξ (t, x)) - potX J (fun y => Y (t, y)) (ξ (t, 0))) t := by
  classical
  intro t ht x hx
  have hSo : IsOpen (I ×ˢ ball (0 : 𝔼) r) := hI.prod isOpen_ball
  have hSn : ∀ z ∈ ball (0 : 𝔼) r, I ×ˢ ball (0 : 𝔼) r ∈ 𝓝 (t, z) := fun z hz =>
    hSo.mem_nhds ⟨ht, hz⟩
  have htab : t ∈ Ioo a b := hIab ht
  have hH : ContDiffOn ℝ ∞ (fderiv ℝ φ₀) (I ×ˢ ball (0 : 𝔼) r) :=
    hφ.fderiv_of_isOpen hSo (by simp)
  have hHd : ∀ z ∈ ball (0 : 𝔼) r, DifferentiableAt ℝ (fderiv ℝ φ₀) (t, z) := fun z hz =>
    ((hH.contDiffAt (hSn z hz)).differentiableAt (by simp))
  have hφd : ∀ z ∈ ball (0 : 𝔼) r, DifferentiableAt ℝ φ₀ (t, z) := fun z hz =>
    ((hφ.contDiffAt (hSn z hz)).differentiableAt (by simp))
  have hq : ∀ z ∈ ball (0 : 𝔼) r,
      HasDerivAt (fun s => φ₀ (s, z)) (fderiv ℝ φ₀ (t, z) (1, 0)) t := fun z hz =>
    tf_hasDerivAt_slice_t (hφd z hz).hasFDerivAt
  have hq0 : fderiv ℝ φ₀ (t, 0) (1, 0) = 0 := by
    have h := hq 0 (mem_ball_self hr)
    have h' : HasDerivAt (fun s => φ₀ (s, 0)) 0 t :=
      (hasDerivAt_const t (0 : ℝ)).congr_of_eventuallyEq
        (by filter_upwards [hI.mem_nhds ht] with s hs using h0 s hs)
    exact h.unique h'
  have hslice : ContDiffOn ℝ ∞ (fun z => (t, z)) (ball (0 : 𝔼) r) :=
    contDiffOn_const.prodMk contDiffOn_id
  have hmaps : MapsTo (fun z : 𝔼 => (t, z)) (ball (0 : 𝔼) r) (I ×ˢ ball (0 : 𝔼) r) :=
    fun z hz => ⟨ht, hz⟩
  have hqt : ContDiffOn ℝ ∞ (fun z => fderiv ℝ φ₀ (t, z) (1, 0)) (ball (0 : 𝔼) r) :=
    (hH.comp hslice hmaps).clm_apply contDiffOn_const
  have hpt : ContDiffOn ℝ ∞ (fun z => potX J (fun y => Y (t, y)) (ξ (t, z))) (ball (0 : 𝔼) r) := by
    refine (vv_contDiffOn_potX (J := J) isOpen_ball (tg_slice hY htab) (hinj t htab)).comp
      (tg_slice hξ ht) ?_
    intro z hz
    exact hξR t ht z hz
  have hder : ∀ z ∈ ball (0 : 𝔼) r,
      fderiv ℝ (fun z => fderiv ℝ φ₀ (t, z) (1, 0)) z =
        fderiv ℝ (fun z => potX J (fun y => Y (t, y)) (ξ (t, z))) z := by
    intro z hz
    refine tf_ext_eb fun i => ?_
    have hd1 : HasFDerivAt (fun z => fderiv ℝ φ₀ (t, z) (1, 0))
        ((ContinuousLinearMap.apply ℝ ℝ ((1 : ℝ), (0 : 𝔼))).comp
          ((fderiv ℝ (fderiv ℝ φ₀) (t, z)).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼))) z :=
      (ContinuousLinearMap.apply ℝ ℝ ((1 : ℝ), (0 : 𝔼))).hasFDerivAt.comp z
        ((hHd z hz).hasFDerivAt.comp z (hasFDerivAt_prodMk_right t z))
    rw [hd1.fderiv, ← tf_step1 hJi hJJ hU hY hinj hlag hvel hI hIab hξ hξR hφ hgr ht hz i]
    have hsymm : IsSymmSndFDerivAt ℝ φ₀ (t, z) :=
      (hφ.contDiffAt (hSn z hz)).isSymmSndFDerivAt (n := ∞)
        (by simp; exact WithTop.coe_le_coe.2 le_top)
    have := hsymm (1, 0) (0, eb i)
    simpa using this.symm
  have hqd : DifferentiableOn ℝ (fun z => fderiv ℝ φ₀ (t, z) (1, 0)) (ball (0 : 𝔼) r) :=
    hqt.differentiableOn (by simp)
  have hpd : DifferentiableOn ℝ (fun z => potX J (fun y => Y (t, y)) (ξ (t, z)) +
      (fderiv ℝ φ₀ (t, 0) (1, 0) - potX J (fun y => Y (t, y)) (ξ (t, 0)))) (ball (0 : 𝔼) r) :=
    (hpt.differentiableOn (by simp)).add_const _
  have heq := eqOn_ball_of_fderiv_eq (x₀ := (0 : 𝔼)) (R := r)
    (φ := fun z => fderiv ℝ φ₀ (t, z) (1, 0))
    (ψ := fun z => potX J (fun y => Y (t, y)) (ξ (t, z)) +
      (fderiv ℝ φ₀ (t, 0) (1, 0) - potX J (fun y => Y (t, y)) (ξ (t, 0))))
    hqd hpd (fun z hz => by
      rw [fderiv_add_const]; exact hder z hz) (by simp)
  have h2 : fderiv ℝ φ₀ (t, x) (1, 0) = potX J (fun y => Y (t, y)) (ξ (t, x)) +
      (fderiv ℝ φ₀ (t, 0) (1, 0) - potX J (fun y => Y (t, y)) (ξ (t, 0))) := heq hx
  have := hq x hx
  rw [h2, hq0] at this
  convert this using 1
  ring

theorem tf_primitive {c d t₀ : ℝ} (ht₀ : t₀ ∈ Ioo c d) {f : ℝ → ℝ}
    (hf : ContDiffOn ℝ ∞ f (Ioo c d)) :
    (∀ t ∈ Ioo c d, HasDerivAt (fun u => ∫ s in t₀..u, f s) (f t) t) ∧
    ContDiffOn ℝ ∞ (fun u => ∫ s in t₀..u, f s) (Ioo c d) := by
  have hcont : ContinuousOn f (Ioo c d) := hf.continuousOn
  have hder : ∀ t ∈ Ioo c d, HasDerivAt (fun u => ∫ s in t₀..u, f s) (f t) t := by
    intro t ht
    have hint : IntervalIntegrable f MeasureTheory.volume t₀ t :=
      (hcont.mono (Set.OrdConnected.uIcc_subset inferInstance ht₀ ht)).intervalIntegrable
    exact intervalIntegral.integral_hasDerivAt_right hint
      (hcont.stronglyMeasurableAtFilter isOpen_Ioo t ht)
      (hcont.continuousAt (isOpen_Ioo.mem_nhds ht))
  refine ⟨hder, (contDiffOn_infty_iff_deriv_of_isOpen isOpen_Ioo).2 ⟨fun t ht =>
    (hder t ht).differentiableAt.differentiableWithinAt, hf.congr fun t ht => (hder t ht).deriv⟩⟩

theorem exists_time_graph_flow {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) (y₀ : V)
    {a b t₀ : ℝ} (ht₀ : t₀ ∈ Ioo a b) {R : ℝ} {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) R))
    (hinj : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R,
      Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    (hgraph : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, horC U y₀ (Y (t₀, y)) = y)
    {r R' : ℝ} (hr : 0 < r) (hrR' : r < R') (hR'R : R' < R) :
    ∃ δ > 0, Ioo (t₀ - δ) (t₀ + δ) ⊆ Ioo a b ∧
      ∃ (ξ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
        (φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
        ContDiffOn ℝ ∞ ξ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) ∧
        ContDiffOn ℝ ∞ φ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) ∧
        ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ),
          InjOn (fun y => horC U y₀ (Y (t, y))) (ball (0 : EuclideanSpace ℝ (Fin n)) R') ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R',
            Function.Injective (fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y) ∧
            ‖fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4) ∧
          (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r,
            ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' ∧ horC U y₀ (Y (t, ξ (t, x))) = x ∧
            Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ (t, x')) x ∧
            HasDerivAt (fun s => φ (s, x)) (potX J (fun y => Y (t, y)) (ξ (t, x))) t) ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R',
            horC U y₀ (Y (t, y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r →
              ξ (t, horC U y₀ (Y (t, y))) = y) := by
  classical
  obtain ⟨δ, hδ, hsub, ξ, φ₀, hξ, hφ, hper⟩ :=
    exists_time_graph hJi hJJ hU hdim y₀ ht₀ hY hlag hgraph hr hrR' hR'R
  have ht₀I : t₀ ∈ Ioo (t₀ - δ) (t₀ + δ) := ⟨by linarith, by linarith⟩
  have hξR : ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ), ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r, ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R :=
    fun t ht x hx => ball_subset_ball hR'R.le ((hper t ht).2.2.1 x hx).1
  have hgr : ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ), ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => φ₀ (t, x')) x :=
    fun t ht x hx => ((hper t ht).2.2.1 x hx).2.2
  have h0 : ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ), φ₀ (t, 0) = 0 := fun t ht => (hper t ht).2.2.2.2
  have hdt := tf_dt_phi hJi hJJ hU hY hinj hlag hvel (isOpen_Ioo (a := t₀ - δ) (b := t₀ + δ)) hsub
    hξ hξR hφ hgr h0 hr
  -- the time-dependent constant
  have hO : IsOpen (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) R) := isOpen_Ioo.prod isOpen_ball
  have hP := contDiffOn_potX_family (J := J) hO hY (fun q hq => hinj q.1 hq.1 q.2 hq.2)
  have hp0 : ContDiffOn ℝ ∞ (fun s => potX J (fun y => Y (s, y)) (ξ (s, 0)))
      (Ioo (t₀ - δ) (t₀ + δ)) := by
    have h1a : ContDiffOn ℝ ∞ (fun s : ℝ => ((s, (0 : EuclideanSpace ℝ (Fin n))) : ℝ × EuclideanSpace ℝ (Fin n))) (Ioo (t₀ - δ) (t₀ + δ)) :=
      contDiffOn_id.prodMk (contDiffOn_const (c := (0 : EuclideanSpace ℝ (Fin n))))
    have h1b : ContDiffOn ℝ ∞ (fun s : ℝ => ξ (s, 0)) (Ioo (t₀ - δ) (t₀ + δ)) :=
      hξ.comp h1a (fun s hs => ⟨hs, mem_ball_self hr⟩)
    have h1 : ContDiffOn ℝ ∞ (fun s : ℝ => ((s, ξ (s, 0)) : ℝ × EuclideanSpace ℝ (Fin n))) (Ioo (t₀ - δ) (t₀ + δ)) :=
      contDiffOn_id.prodMk h1b
    have hmaps : MapsTo (fun s : ℝ => ((s, ξ (s, 0)) : ℝ × EuclideanSpace ℝ (Fin n)))
        (Ioo (t₀ - δ) (t₀ + δ)) (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) R) :=
      fun s hs => ⟨hsub hs, hξR s hs 0 (mem_ball_self hr)⟩
    have hc := hP.comp h1 hmaps
    refine hc.congr fun s hs => ?_
    simp only [Function.comp_apply]
  obtain ⟨hcd, hcC⟩ := tf_primitive ht₀I hp0
  refine ⟨δ, hδ, hsub, ξ, fun q => φ₀ q + ∫ s in t₀..q.1, potX J (fun y => Y (s, y)) (ξ (s, 0)),
    hξ, ?_, ?_⟩
  · exact hφ.add (hcC.comp contDiffOn_fst (fun q hq => hq.1))
  intro t ht
  obtain ⟨hinj', hder, hpos, hleft, -⟩ := hper t ht
  refine ⟨hinj', hder, fun x hx => ⟨(hpos x hx).1, (hpos x hx).2.1, ?_, ?_⟩, hleft⟩
  · rw [(hpos x hx).2.2]
    simp only [graphX, gradE, fderiv_add_const]
  · have h1 := hdt t ht x hx
    have h2 := h1.add (hcd t ht)
    exact h2.congr_deriv (by ring)

/-- info: 'CWFlow.Extension.contDiffOn_potX_family' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_potX_family

/-- info: 'CWFlow.Extension.contDiffOn_velX_family' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_velX_family

/-- info: 'CWFlow.Extension.exists_time_graph_flow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_time_graph_flow

end CWFlow.Extension
