import CWFlow.Extension.ImmDiv
import Mathlib.Topology.Instances.Matrix
import Mathlib.Analysis.Calculus.FDeriv.Symmetric

/-!
# First variation of the volume density, pointwise

For immersions `X : ℝⁿ → V` in coordinates:

* Jacobi's formula (`vv_hasFDerivAt_det`) and `∂_d √det g = √det g · Γ^e_{de}` (`hasFDerivAt_volX`);
* total symmetry of `C` in its last two indices for Lagrangian `X` and constant `J` (`cX_swap23`);
* the time derivative of `√det g` for a family `Y` (`hasDerivAt_volX_slice`);
* the variation along `velX` (`sum_gX_inv_inner_fderiv_velX`);
* the divergence identity `Σ_d ∂_d (√det g g^{dc} μ_c) = √det g · potX`
  (`sum_fderiv_volX_muC`).
-/

open scoped ContDiff Matrix Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X : EuclideanSpace ℝ (Fin n) → V} {Ω : Set (EuclideanSpace ℝ (Fin n))}
  {x : EuclideanSpace ℝ (Fin n)}

/-! ### Jacobi's formula -/

section Jacobi

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `det` as a continuous multilinear map of the rows. -/
noncomputable def vv_detCML : ContinuousMultilinearMap ℝ (fun _ : ι => ι → ℝ) ℝ :=
  { (Matrix.detRowAlternating : (ι → ℝ) [⋀^ι]→ₗ[ℝ] ℝ).toMultilinearMap with
    cont := (continuous_id.matrix_det : Continuous fun m : Matrix ι ι ℝ => m.det) }

theorem vv_detCML_apply (m : ι → ι → ℝ) : vv_detCML m = (Matrix.of m).det := rfl

theorem vv_det_updateRow_expand (A : Matrix ι ι ℝ) (i : ι) (u : ι → ℝ) :
    (A.updateRow i u).det = ∑ j, u j * A.adjugate j i := by
  let L := (Matrix.detRowAlternating : (ι → ℝ) [⋀^ι]→ₗ[ℝ] ℝ).toMultilinearMap.toLinearMap A i
  have hL : ∀ w, (A.updateRow i w).det = L w := fun w => rfl
  have hu : u = ∑ j, u j • (Pi.single j 1 : ι → ℝ) := by
    ext k; simp [Finset.sum_apply, Pi.single_apply]
  calc (A.updateRow i u).det = L (∑ j, u j • (Pi.single j 1 : ι → ℝ)) := by rw [← hu]; rfl
    _ = ∑ j, u j * A.adjugate j i := by
      rw [map_sum]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [map_smul, ← hL, Matrix.adjugate_apply, smul_eq_mul]

/-- Jacobi's formula: `d det M = Σᵢⱼ adj(M)ⱼᵢ dMᵢⱼ`. -/
theorem vv_hasFDerivAt_det {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    {M : F → Matrix ι ι ℝ} {M' : ι → ι → F →L[ℝ] ℝ} {x : F}
    (h : ∀ i j, HasFDerivAt (fun y => M y i j) (M' i j) x) :
    HasFDerivAt (fun y => (M y).det) (∑ i, ∑ j, (M x).adjugate j i • M' i j) x := by
  have hg : ∀ i, HasFDerivAt (fun y => (M y i : ι → ℝ)) (ContinuousLinearMap.pi (M' i)) x :=
    fun i => hasFDerivAt_pi.2 (h i)
  have key := HasFDerivAt.multilinear_comp (f := vv_detCML (ι := ι)) hg
  have key' : HasFDerivAt (fun y => (M y).det) _ x := key
  refine key'.congr_fderiv ?_
  ext v
  simp only [_root_.sum_apply, _root_.smul_apply,
    ContinuousLinearMap.comp_apply, ContinuousMultilinearMap.toContinuousLinearMap_apply, smul_eq_mul]
  refine Finset.sum_congr rfl fun i _ => ?_
  have := vv_det_updateRow_expand (M x) i fun j => M' i j v
  rw [Finset.sum_congr rfl fun j _ => mul_comm ..] at this
  rw [← this]
  rfl

/-- Jacobi for the square root of the determinant. -/
theorem vv_hasFDerivAt_sqrt_det {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    {M : F → Matrix ι ι ℝ} {M' : ι → ι → F →L[ℝ] ℝ} {x : F}
    (h : ∀ i j, HasFDerivAt (fun y => M y i j) (M' i j) x) (hpos : 0 < (M x).det) :
    HasFDerivAt (fun y => Real.sqrt (M y).det)
      ((1 / (2 * Real.sqrt (M x).det)) • ∑ i, ∑ j, (M x).adjugate j i • M' i j) x :=
  (vv_hasFDerivAt_det h).sqrt hpos.ne'

end Jacobi

theorem vv_adj_eq (A : Matrix (Fin n) (Fin n) ℝ) (h : A.det ≠ 0) (a b : Fin n) :
    A.adjugate a b = A.det * A⁻¹ a b := by
  rw [Matrix.inv_def, Ring.inverse_eq_inv']
  show _ = A.det * (A.det⁻¹ * A.adjugate a b)
  field_simp

theorem vv_ddX_symm (hXa : ContDiffAt ℝ ∞ X x) (a b : Fin n) : ddX X x a b = ddX X x b a := by
  have hs : IsSymmSndFDerivAt ℝ X x :=
    hXa.isSymmSndFDerivAt (by
      rw [minSmoothness_of_isRCLikeNormedField]; exact WithTop.coe_le_coe.2 le_top)
  exact hs (eb a) (eb b)

theorem vv_gamX_symm (hXa : ContDiffAt ℝ ∞ X x) (a b c : Fin n) :
    gamX X x ![a, b, c] = gamX X x ![b, a, c] := by
  simp only [ip_gamX_eq, vv_ddX_symm hXa a b]

theorem vv_trace (hXa : ContDiffAt ℝ ∞ X x) (d : Fin n) :
    ∑ i, ∑ j, (gX X x)⁻¹ j i * fderiv ℝ (fun y => gX X y i j) x (eb d) =
      2 * ∑ e, gamX X x ![d, e, e] := by
  simp only [ip_fderiv_gX hXa, mul_add, Finset.sum_add_distrib]
  have hA : ∑ i, ∑ j, (gX X x)⁻¹ j i * ⟪ddX X x d i, dX X x j⟫ = ∑ e, gamX X x ![d, e, e] := by
    simp only [ip_gamX_eq]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [ip_ginv_symm j i]
  have hB : ∑ i, ∑ j, (gX X x)⁻¹ j i * ⟪dX X x i, ddX X x d j⟫ = ∑ e, gamX X x ![d, e, e] := by
    rw [Finset.sum_comm]
    simp only [ip_gamX_eq]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [real_inner_comm]
  rw [hA, hB]; ring

theorem hasFDerivAt_volX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) :
    HasFDerivAt (volX X)
      (volX X x • ∑ e, ∑ d, gamX X x ![d, e, e] •
        (EuclideanSpace.proj d : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ)) x := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hpos := gX_det_pos (hinj x hx)
  have hJ := vv_hasFDerivAt_sqrt_det (M := gX X)
    (M' := fun i j => fderiv ℝ (fun y => gX X y i j) x)
    (fun i j => (ip_diffAt_gX hXa i j).hasFDerivAt) hpos
  refine hJ.congr_fderiv ?_
  refine ip_ext_eb fun d => ?_
  have hs : Real.sqrt (gX X x).det * Real.sqrt (gX X x).det = (gX X x).det :=
    Real.mul_self_sqrt hpos.le
  have hs0 : 0 < Real.sqrt (gX X x).det := Real.sqrt_pos.2 hpos
  have hT := vv_trace hXa d
  simp only [_root_.smul_apply, _root_.sum_apply, smul_eq_mul,
    EuclideanSpace.coe_proj, PiLp.single_apply]
  simp only [vv_adj_eq _ hpos.ne', mul_assoc, ← Finset.mul_sum] at *
  rw [hT]
  simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  unfold volX
  generalize Real.sqrt (gX X x).det = s at hs hs0 ⊢
  field_simp
  rw [← hs]
  ring

/-! ### Total symmetry of `C` -/

theorem cX_swap23 (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hx : x ∈ Ω)
    (a b c : Fin n) :
    cX J X x ![a, b, c] = cX J X x ![a, c, b] := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hev : (fun y => ⟪dX X y b, J (dX X y c)⟫) =ᶠ[𝓝 x] fun _ => (0 : ℝ) := by
    filter_upwards [hΩ.mem_nhds hx] with y hy
    exact ip_inner_L_J_L hJi hJJ (hlag y hy) (eb b) (eb c)
  have hg : HasFDerivAt (fun y => J (dX X y c)) (J.comp (fderiv ℝ (fun y => dX X y c) x)) x :=
    J.hasFDerivAt.comp x (ip_diff_dX hXa c).hasFDerivAt
  have h := fderiv_inner_apply ℝ (ip_diff_dX hXa b) hg.differentiableAt (eb a)
  rw [hg.fderiv] at h
  simp only [ContinuousLinearMap.comp_apply] at h
  rw [ip_fderiv_dX hXa a c] at h
  have h0 : fderiv ℝ (fun t => ⟪dX X t b, J (dX X t c)⟫) x = 0 := by
    rw [hev.fderiv_eq]; simp
  rw [h0, ip_fderiv_dX hXa a b] at h
  simp only [zero_apply] at h
  have h2 := ip_inner_J_left hJi hJJ (dX X x b) (ddX X x a c)
  have e1 : cX J X x ![a, b, c] = ⟪ddX X x a b, J (dX X x c)⟫ := rfl
  have e2 : cX J X x ![a, c, b] = ⟪ddX X x a c, J (dX X x b)⟫ := rfl
  rw [e1, e2, real_inner_comm (J (dX X x b)) (ddX X x a c)] at *
  linarith

/-! ### Time derivative of the volume density -/

section Slice

variable {G W : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G] [NormedAddCommGroup W]
  [NormedSpace ℝ W]

theorem vv_hasFDerivAt_apply {Y : G → W} {q : G} (hY : ContDiffAt ℝ ∞ Y q) (w : G) :
    HasFDerivAt (fun q' => fderiv ℝ Y q' w)
      ((ContinuousLinearMap.apply ℝ W w).comp (fderiv ℝ (fderiv ℝ Y) q)) q :=
  (ContinuousLinearMap.apply ℝ W w).hasFDerivAt.comp q
    ((hY.fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)).hasFDerivAt

theorem vv_hasDerivAt_slice_t {Y : ℝ × G → W} {t : ℝ} {y : G} (hY : DifferentiableAt ℝ Y (t, y)) :
    HasDerivAt (fun s => Y (s, y)) (fderiv ℝ Y (t, y) (1, 0)) t := by
  have h : HasDerivAt (fun s : ℝ => (s, y)) ((1 : ℝ), (0 : G)) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t y)
  exact hY.hasFDerivAt.comp_hasDerivAt t h

theorem vv_hasFDerivAt_slice_y {Y : ℝ × G → W} {t : ℝ} {y : G} (hY : DifferentiableAt ℝ Y (t, y)) :
    HasFDerivAt (fun y' => Y (t, y')) ((fderiv ℝ Y (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ G)) y :=
  hY.hasFDerivAt.comp y (hasFDerivAt_prodMk_right t y)

end Slice

theorem hasDerivAt_volX_slice {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))}
    (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2))
    {t : ℝ} {y : EuclideanSpace ℝ (Fin n)} (hq : (t, y) ∈ O) :
    HasDerivAt (fun s => volX (fun y' => Y (s, y')) y)
      (volX (fun y' => Y (t, y')) y * ∑ a, ∑ b, (gX (fun y' => Y (t, y')) y)⁻¹ a b *
        ⟪fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb a),
          dX (fun y' => Y (t, y')) y b⟫) t := by
  have hCA : ∀ q ∈ O, ContDiffAt ℝ ∞ Y q := fun q hq' => hY.contDiffAt (hO.mem_nhds hq')
  have hDY : ∀ q ∈ O, DifferentiableAt ℝ Y q := fun q hq' => (hCA q hq').differentiableAt (by simp)
  have hYa := hCA _ hq
  -- slices of `dX`
  have hdX : ∀ s, (s, y) ∈ O → ∀ p : Fin n,
      dX (fun y' => Y (s, y')) y p = fderiv ℝ Y (s, y) (0, eb p) := by
    intro s hs p
    unfold dX
    rw [(vv_hasFDerivAt_slice_y (hDY _ hs)).fderiv]
    simp
  -- the time derivative of the slice
  have hdt : ∀ y', (t, y') ∈ O → deriv (fun s => Y (s, y')) t = fderiv ℝ Y (t, y') (1, 0) :=
    fun y' hy' => (vv_hasDerivAt_slice_t (hDY _ hy')).deriv
  have hnhd : {y' | (t, y') ∈ O} ∈ 𝓝 y :=
    (continuous_const.prodMk continuous_id).continuousAt.preimage_mem_nhds
      (hO.mem_nhds hq)
  have hev : (fun y' => deriv (fun s => Y (s, y')) t) =ᶠ[𝓝 y]
      fun y' => fderiv ℝ Y (t, y') (1, 0) := by
    filter_upwards [hnhd] with y' hy' using hdt y' hy'
  have hVv : ∀ p : Fin n, fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb p) =
      fderiv ℝ (fderiv ℝ Y) (t, y) (0, eb p) (1, 0) := by
    intro p
    rw [hev.fderiv_eq]
    have h : HasFDerivAt (fun y' => fderiv ℝ Y (t, y') (1, 0)) _ y :=
      (vv_hasFDerivAt_apply hYa ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n)))).comp y
        (hasFDerivAt_prodMk_right t y)
    rw [h.fderiv]
    simp
  have hsymm : IsSymmSndFDerivAt ℝ Y (t, y) :=
    hYa.isSymmSndFDerivAt (by
      rw [minSmoothness_of_isRCLikeNormedField]; exact WithTop.coe_le_coe.2 le_top)
  -- derivative in `s` of `dX`
  have hnhd2 : {s | (s, y) ∈ O} ∈ 𝓝 t :=
    (continuous_id.prodMk continuous_const).continuousAt.preimage_mem_nhds (hO.mem_nhds hq)
  have hds : ∀ p : Fin n, HasDerivAt (fun s => dX (fun y' => Y (s, y')) y p)
      (fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb p)) t := by
    intro p
    have h : HasDerivAt (fun s : ℝ => (s, y)) ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) t :=
      (hasDerivAt_id t).prodMk (hasDerivAt_const t y)
    have h2 := (vv_hasFDerivAt_apply hYa ((0 : ℝ), eb p)).comp_hasDerivAt t h
    have h3 : HasDerivAt (fun s => fderiv ℝ Y (s, y) (0, eb p))
        (fderiv ℝ (fderiv ℝ Y) (t, y) (1, 0) (0, eb p)) t := by
      exact h2
    rw [hVv p]
    rw [hsymm (1, 0) (0, eb p)] at h3
    refine h3.congr_of_eventuallyEq ?_
    filter_upwards [hnhd2] with s hs using hdX s hs p
  have hpos : 0 < (gX (fun y' => Y (t, y')) y).det := gX_det_pos (hinj _ hq)
  have hM : ∀ p r : Fin n, HasDerivAt (fun s => gX (fun y' => Y (s, y')) y p r)
      (⟪dX (fun y' => Y (t, y')) y p,
          fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb r)⟫ +
        ⟪fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb p),
          dX (fun y' => Y (t, y')) y r⟫) t := fun p r =>
    HasDerivAt.inner ℝ (hds p) (hds r)
  have hF := vv_hasFDerivAt_sqrt_det (M := fun s => gX (fun y' => Y (s, y')) y)
    (M' := fun p r => (1 : ℝ →L[ℝ] ℝ).smulRight
      (⟪dX (fun y' => Y (t, y')) y p,
          fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb r)⟫ +
        ⟪fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb p),
          dX (fun y' => Y (t, y')) y r⟫))
    (fun p r => (hM p r).hasFDerivAt) hpos
  refine HasDerivAt.congr_deriv hF.hasDerivAt ?_
  simp only [_root_.smul_apply, _root_.sum_apply, smul_eq_mul, ContinuousLinearMap.smulRight_apply,
    ContinuousLinearMap.one_apply, one_mul]
  simp only [vv_adj_eq _ hpos.ne']
  set Xt : EuclideanSpace ℝ (Fin n) → V := fun y' => Y (t, y') with hXt
  set Vv : Fin n → V := fun p =>
    fderiv ℝ (fun y' => deriv (fun s => Y (s, y')) t) y (eb p) with hVvdef
  have hinjt : Function.Injective (fderiv ℝ Xt y) := hinj _ hq
  have hA : ∑ i, ∑ j, (gX Xt y)⁻¹ j i * ⟪dX Xt y i, Vv j⟫ =
      ∑ a, ∑ b, (gX Xt y)⁻¹ a b * ⟪Vv a, dX Xt y b⟫ := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    rw [real_inner_comm]
  have hB : ∑ i, ∑ j, (gX Xt y)⁻¹ j i * ⟪Vv i, dX Xt y j⟫ =
      ∑ a, ∑ b, (gX Xt y)⁻¹ a b * ⟪Vv a, dX Xt y b⟫ :=
    Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by
      rw [gX_inv_symm hinjt]
  have hs : Real.sqrt (gX Xt y).det * Real.sqrt (gX Xt y).det = (gX Xt y).det :=
    Real.mul_self_sqrt hpos.le
  have hs0 : 0 < Real.sqrt (gX Xt y).det := Real.sqrt_pos.2 hpos
  have hsum : ∑ i, ∑ j, (gX Xt y).det * (gX Xt y)⁻¹ j i * (⟪dX Xt y i, Vv j⟫ + ⟪Vv i, dX Xt y j⟫) =
      (gX Xt y).det * (2 * ∑ a, ∑ b, (gX Xt y)⁻¹ a b * ⟪Vv a, dX Xt y b⟫) := by
    have e : ∀ i j, (gX Xt y).det * (gX Xt y)⁻¹ j i * (⟪dX Xt y i, Vv j⟫ + ⟪Vv i, dX Xt y j⟫) =
        (gX Xt y).det * ((gX Xt y)⁻¹ j i * ⟪dX Xt y i, Vv j⟫ +
          (gX Xt y)⁻¹ j i * ⟪Vv i, dX Xt y j⟫) := fun i j => by ring
    simp only [e, ← Finset.mul_sum, Finset.sum_add_distrib]
    rw [hA, hB]
    ring
  rw [hsum]
  unfold volX
  generalize Real.sqrt (gX Xt y).det = s at hs hs0 ⊢
  field_simp
  rw [← hs]
  ring

/-! ### The variation along `velX` -/

theorem vv_contDiffOn_potX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) : ContDiffOn ℝ ∞ (potX J X) Ω := by
  unfold potX
  refine ContDiffOn.neg ?_
  refine ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => ContDiffOn.sum fun c _ =>
    ContDiffOn.sum fun d _ => ?_
  have h1 : ∀ I : Fin 2 → Fin n, ContDiffOn ℝ ∞ (fun y => ginvZ X y I) Ω := fun I =>
    (contDiff_coordG I).comp_contDiffOn (contDiffOn_ginvZ hΩ hX hinj)
  have h2 : ContDiffOn ℝ ∞ (fun y => SX J X 1 y ![d, a, b, c]) Ω :=
    (contDiff_coordG (![d, a, b, c] : Fin 4 → Fin n)).comp_contDiffOn
      (contDiffOn_SX (J := J) hΩ hX hinj 1)
  exact ((h1 ![c, d]).mul (h1 ![a, b])).mul h2

theorem sum_gX_inv_inner_fderiv_velX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hx : x ∈ Ω) :
    ∑ a, ∑ b, (gX X x)⁻¹ a b * ⟪fderiv ℝ (velX J X) x (eb a), dX X x b⟫ =
      ∑ c, ∑ d, (gX X x)⁻¹ c d * muC J X x c * fderiv ℝ (potX J X) x (eb d) := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hp : ContDiffOn ℝ ∞ (potX J X) Ω := vv_contDiffOn_potX hΩ hX hinj
  have hpa : ContDiffAt ℝ ∞ (potX J X) x := hp.contDiffAt (hΩ.mem_nhds hx)
  have hdp : ∀ b, DifferentiableAt ℝ (fun y => fderiv ℝ (potX J X) y (eb b)) x := fun b =>
    ((hpa.fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)).clm_apply
      (differentiableAt_const _)
  set f : Fin n → EuclideanSpace ℝ (Fin n) → ℝ :=
    fun a' y => ∑ b, (gX X y)⁻¹ a' b * fderiv ℝ (potX J X) y (eb b) with hfdef
  have hf : ∀ a', DifferentiableAt ℝ (f a') x := fun a' =>
    DifferentiableAt.fun_sum fun b _ => (ip_diffAt_ginv hΩ hX hinj hx a' b).mul (hdp b)
  have hvel : velX J X = fun y => ∑ a', f a' y • J (dX X y a') := by
    funext y
    unfold velX
    refine Finset.sum_congr rfl fun a' _ => ?_
    rw [Finset.sum_smul]
  have hH : HasFDerivAt (velX J X)
      (∑ a', (f a' x • J.comp ((ContinuousLinearMap.apply ℝ V (eb a')).comp
          (fderiv ℝ (fderiv ℝ X) x)) +
        (fderiv ℝ (f a') x).smulRight (J (dX X x a')))) x := by
    rw [hvel]
    exact HasFDerivAt.fun_sum fun a' _ =>
      (hf a').hasFDerivAt.smul (J.hasFDerivAt.comp x (ip_hasFDerivAt_dX hXa a'))
  have hvA : ∀ a, fderiv ℝ (velX J X) x (eb a) =
      ∑ a', (f a' x • J (ddX X x a a') + fderiv ℝ (f a') x (eb a) • J (dX X x a')) := by
    intro a
    rw [hH.fderiv]
    simp only [_root_.sum_apply, _root_.add_apply, _root_.smul_apply,
      ContinuousLinearMap.smulRight_apply, ContinuousLinearMap.comp_apply]
    rfl
  have hpair : ∀ a b, ⟪fderiv ℝ (velX J X) x (eb a), dX X x b⟫ =
      -∑ a', f a' x * cX J X x ![a, b, a'] := by
    intro a b
    rw [hvA a, sum_inner, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun a' _ => ?_
    rw [inner_add_left, real_inner_smul_left, real_inner_smul_left]
    have h0 : ⟪J (dX X x a'), dX X x b⟫ = 0 := hlag x hx _ _
    have h1 : ⟪J (ddX X x a a'), dX X x b⟫ = -cX J X x ![a, a', b] := by
      rw [ip_inner_J_left hJi hJJ]; rfl
    rw [h0, h1, cX_swap23 hJi hJJ hΩ hX hlag hx a a' b]
    ring
  have hmu : ∀ a', muC J X x a' = -∑ a, ∑ b, (gX X x)⁻¹ a b * cX J X x ![a, b, a'] :=
    fun a' => rfl
  calc ∑ a, ∑ b, (gX X x)⁻¹ a b * ⟪fderiv ℝ (velX J X) x (eb a), dX X x b⟫
      = ∑ a, ∑ b, ∑ a', -((gX X x)⁻¹ a b * (f a' x * cX J X x ![a, b, a'])) := by
        refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
        rw [hpair a b, mul_neg, Finset.mul_sum, ← Finset.sum_neg_distrib]
    _ = ∑ a', f a' x * muC J X x a' := by
        have key : ∀ a', f a' x * muC J X x a' =
            ∑ a, ∑ b, -((gX X x)⁻¹ a b * (f a' x * cX J X x ![a, b, a'])) := by
          intro a'
          rw [hmu, mul_neg, Finset.mul_sum, ← Finset.sum_neg_distrib]
          refine Finset.sum_congr rfl fun a _ => ?_
          rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
          refine Finset.sum_congr rfl fun b _ => ?_
          ring
        rw [Finset.sum_congr rfl fun a' _ => key a']
        symm
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = _ := by
        refine Finset.sum_congr rfl fun c _ => ?_
        simp only [hfdef, Finset.sum_mul]
        refine Finset.sum_congr rfl fun d _ => ?_
        ring

/-! ### The divergence identity -/

theorem vv_sum3_rev (f : Fin n → Fin n → Fin n → ℝ) :
    ∑ a, ∑ b, ∑ c, f a b c = ∑ c, ∑ b, ∑ a, f a b c :=
  calc ∑ a, ∑ b, ∑ c, f a b c = ∑ b, ∑ a, ∑ c, f a b c := Finset.sum_comm
    _ = ∑ b, ∑ c, ∑ a, f a b c := Finset.sum_congr rfl fun b _ => Finset.sum_comm
    _ = ∑ c, ∑ b, ∑ a, f a b c := Finset.sum_comm

theorem vv_alg5 (gi : Fin n → Fin n → ℝ) (hs : ∀ a b, gi a b = gi b a)
    (Γ : Fin n → Fin n → Fin n → ℝ) (hΓ : ∀ a b c, Γ a b c = Γ b a c) (μ : Fin n → ℝ)
    (Dm : Fin n → Fin n → ℝ) (v : ℝ) :
    ∑ d, (v * ∑ c, ((-∑ e, (gi d e * Γ d e c + gi c e * Γ d e d)) * μ c + gi d c * Dm d c) +
        (∑ c, gi d c * μ c) * (v * ∑ e, Γ d e e)) =
      v * ∑ a, ∑ c, gi a c * (Dm a c - ∑ e, Γ a c e * μ e) := by
  have hA : ∑ d, ∑ c, ∑ e, gi d e * Γ d e c * μ c = ∑ a, ∑ c, gi a c * ∑ e, Γ a c e * μ e := by
    refine Finset.sum_congr rfl fun d _ => ?_
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun e _ => by ring
  have hB : ∑ d, ∑ c, ∑ e, gi c e * Γ d e d * μ c = ∑ d, (∑ c, gi d c * μ c) * ∑ e, Γ d e e := by
    have h1 : ∑ d, (∑ c, gi d c * μ c) * ∑ e, Γ d e e =
        ∑ d, ∑ c, ∑ e, gi d c * μ c * Γ d e e := by
      refine Finset.sum_congr rfl fun d _ => ?_
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun c _ => ?_
      rw [Finset.mul_sum]
    rw [h1, vv_sum3_rev (fun d c e => gi d c * μ c * Γ d e e)]
    exact Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun c _ =>
      Finset.sum_congr rfl fun e _ => by rw [hs c e, hΓ d e d]; ring
  have e1 : ∀ d, v * ∑ c, ((-∑ e, (gi d e * Γ d e c + gi c e * Γ d e d)) * μ c + gi d c * Dm d c) +
        (∑ c, gi d c * μ c) * (v * ∑ e, Γ d e e) =
      v * (-(∑ c, ∑ e, gi d e * Γ d e c * μ c) - (∑ c, ∑ e, gi c e * Γ d e d * μ c) +
        ∑ c, gi d c * Dm d c + (∑ c, gi d c * μ c) * ∑ e, Γ d e e) := by
    intro d
    have hW : (∑ c, gi d c * μ c) * (v * ∑ e, Γ d e e) =
        v * ((∑ c, gi d c * μ c) * ∑ e, Γ d e e) := by ring
    rw [hW]
    simp only [add_mul, Finset.sum_add_distrib, Finset.sum_neg_distrib, neg_mul, Finset.sum_mul]
    ring
  rw [Finset.sum_congr rfl fun d _ => e1 d, ← Finset.mul_sum]
  congr 1
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_neg_distrib, hA, hB,
    mul_sub, Finset.sum_sub_distrib]
  ring

theorem vv_fderiv_volX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) (d : Fin n) :
    fderiv ℝ (volX X) x (eb d) = volX X x * ∑ e, gamX X x ![d, e, e] := by
  rw [(hasFDerivAt_volX hΩ hX hinj hx).fderiv]
  simp only [_root_.smul_apply, _root_.sum_apply, smul_eq_mul,
    EuclideanSpace.coe_proj, PiLp.single_apply]
  simp only [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]

theorem sum_fderiv_volX_muC (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) :
    ∑ d, fderiv ℝ (fun y => volX X y * ∑ c, (gX X y)⁻¹ d c * muC J X y c) x (eb d) =
      volX X x * potX J X x := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hvol : HasFDerivAt (volX X) (fderiv ℝ (volX X) x) x :=
    (hasFDerivAt_volX hΩ hX hinj hx).differentiableAt.hasFDerivAt
  have hmuH : ∀ c, HasFDerivAt (fun y => muC J X y c)
      (-∑ a, ∑ b, ((gX X x)⁻¹ a b • fderiv ℝ (fun y => cX J X y ![a, b, c]) x +
        cX J X x ![a, b, c] • fderiv ℝ (fun y => (gX X y)⁻¹ a b) x)) x := fun c =>
    (HasFDerivAt.fun_sum fun a _ => HasFDerivAt.fun_sum fun b _ =>
      ((ip_diffAt_ginv hΩ hX hinj hx a b).hasFDerivAt).fun_mul
        (ip_diff_cX hXa ![a, b, c]).hasFDerivAt).fun_neg
  have hmu : ∀ c d, fderiv ℝ (fun y => muC J X y c) x (eb d) =
      -∑ a, ∑ b, ((gX X x)⁻¹ a b * fderiv ℝ (fun y => cX J X y ![a, b, c]) x (eb d) +
        cX J X x ![a, b, c] * fderiv ℝ (fun y => (gX X y)⁻¹ a b) x (eb d)) := by
    intro c d
    rw [(hmuH c).fderiv]
    simp only [_root_.neg_apply, _root_.sum_apply, _root_.add_apply, _root_.smul_apply,
      smul_eq_mul]
  have hw : ∀ d, HasFDerivAt (fun y => ∑ c, (gX X y)⁻¹ d c * muC J X y c)
      (∑ c, ((gX X x)⁻¹ d c • fderiv ℝ (fun y => muC J X y c) x +
        muC J X x c • fderiv ℝ (fun y => (gX X y)⁻¹ d c) x)) x := fun d =>
    HasFDerivAt.fun_sum fun c _ =>
      ((ip_diffAt_ginv hΩ hX hinj hx d c).hasFDerivAt).fun_mul (hmuH c).differentiableAt.hasFDerivAt
  have hterm : ∀ d, fderiv ℝ (fun y => volX X y * ∑ c, (gX X y)⁻¹ d c * muC J X y c) x (eb d) =
      volX X x * ∑ c, ((-∑ e, ((gX X x)⁻¹ d e * gamX X x ![d, e, c] +
          (gX X x)⁻¹ c e * gamX X x ![d, e, d])) * muC J X x c +
        (gX X x)⁻¹ d c * fderiv ℝ (fun y => muC J X y c) x (eb d)) +
      (∑ c, (gX X x)⁻¹ d c * muC J X x c) * (volX X x * ∑ e, gamX X x ![d, e, e]) := by
    intro d
    rw [(hvol.fun_mul (hw d)).fderiv]
    simp only [_root_.add_apply, _root_.smul_apply, _root_.sum_apply, smul_eq_mul]
    rw [vv_fderiv_volX hΩ hX hinj hx d]
    congr 2
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [fderiv_gX_inv hΩ hX hinj hx d c d]
    ring
  rw [Finset.sum_congr rfl fun d _ => hterm d]
  refine (vv_alg5 (fun a b => (gX X x)⁻¹ a b) (fun a b => ip_ginv_symm a b)
    (fun a b c => gamX X x ![a, b, c]) (fun a b c => vv_gamX_symm hXa a b c)
    (fun c => muC J X x c) (fun a c => fderiv ℝ (fun y => muC J X y c) x (eb a))
    (volX X x)).trans ?_
  rw [potX_eq_div hΩ hX hinj hx]

/-- info: 'CWFlow.Extension.vv_hasFDerivAt_det' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms vv_hasFDerivAt_det

/-- info: 'CWFlow.Extension.hasFDerivAt_volX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasFDerivAt_volX

/-- info: 'CWFlow.Extension.cX_swap23' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cX_swap23

/-- info: 'CWFlow.Extension.hasDerivAt_volX_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_volX_slice

/-- info: 'CWFlow.Extension.sum_gX_inv_inner_fderiv_velX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sum_gX_inv_inner_fderiv_velX

/-- info: 'CWFlow.Extension.sum_fderiv_volX_muC' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sum_fderiv_volX_muC

end CWFlow.Extension
