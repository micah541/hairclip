import CWFlow.Scalar.Graph

/-!
# Time-dependent graphs on a closed time interval

Versions of the global-in-time statements of `CWFlow.Scalar.Graph` for `w : ℝ × E → ℝ` that is only
smooth on `Icc 0 T ×ˢ U` (`U` open), with one-sided time derivatives `dtW`.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
variable {w : ℝ × E → W} {T : ℝ} {U : Set E}

/-- The one-sided time derivative of `w` on `[0, T]`. -/
def dtW (w : ℝ × E → W) (T : ℝ) (p : ℝ × E) : W := derivWithin (fun s => w (s, p.2)) (Icc 0 T) p.1

private theorem uniqueDiffOn_S (hT : 0 < T) (hU : IsOpen U) :
    UniqueDiffOn ℝ (Icc 0 T ×ˢ U) := (uniqueDiffOn_Icc hT).prod hU.uniqueDiffOn

private theorem contDiffOn_D (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) :
    ContDiffOn ℝ ∞ (fderivWithin ℝ w (Icc 0 T ×ˢ U)) (Icc 0 T ×ˢ U) :=
  hw.fderivWithin (uniqueDiffOn_S hT hU) (m := ∞) (by simp)

private theorem hasDerivWithinAt_time_aux (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U) :
    HasDerivWithinAt (fun s => w (s, y)) (fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y) (1, 0))
      (Icc 0 T) t := by
  have hmem : (t, y) ∈ Icc 0 T ×ˢ U := ⟨ht, hy⟩
  have hwd : HasFDerivWithinAt w (fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y))
      (Icc 0 T ×ˢ U) (t, y) :=
    (hw.differentiableOn (by simp) _ hmem).hasFDerivWithinAt
  exact hwd.comp_hasDerivWithinAt (f := fun s => (s, y)) t
    ((hasDerivWithinAt_id t _).prodMk (hasDerivAt_const t y).hasDerivWithinAt)
    (fun s hs => ⟨hs, hy⟩)

private theorem hasFDerivAt_slice_aux (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U) :
    HasFDerivAt (fun z => w (t, z))
      ((fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y := by
  have hmem : (t, y) ∈ Icc 0 T ×ˢ U := ⟨ht, hy⟩
  have hwd : HasFDerivWithinAt w (fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y))
      (Icc 0 T ×ˢ U) (t, y) :=
    (hw.differentiableOn (by simp) _ hmem).hasFDerivWithinAt
  exact (hwd.comp y (hasFDerivAt_prodMk_right t y).hasFDerivWithinAt
    (fun z hz => ⟨ht, hz⟩)).hasFDerivAt (hU.mem_nhds hy)

theorem hasDerivWithinAt_time (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U) :
    HasDerivWithinAt (fun s => w (s, y)) (dtW w T (t, y)) (Icc 0 T) t := by
  have h := hasDerivWithinAt_time_aux hT hU hw ht hy
  have : dtW w T (t, y) = fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y) (1, 0) :=
    h.derivWithin (uniqueDiffOn_Icc hT t ht)
  rw [this]; exact h

private theorem dtW_eq (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {p : ℝ × E} (hp : p ∈ Icc 0 T ×ˢ U) :
    dtW w T p = fderivWithin ℝ w (Icc 0 T ×ˢ U) p (1, 0) :=
  (hasDerivWithinAt_time_aux hT hU hw hp.1 hp.2).derivWithin (uniqueDiffOn_Icc hT p.1 hp.1)

private theorem slice_eq (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {p : ℝ × E} (hp : p ∈ Icc 0 T ×ˢ U) :
    fderiv ℝ (fun z => w (p.1, z)) p.2 =
      (fderivWithin ℝ w (Icc 0 T ×ˢ U) p).comp (ContinuousLinearMap.inr ℝ ℝ E) :=
  (hasFDerivAt_slice_aux hT hU hw hp.1 hp.2).fderiv

theorem contDiffOn_dtW (hT : 0 < T) (hU : IsOpen U) (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) :
    ContDiffOn ℝ ∞ (dtW w T) (Icc 0 T ×ˢ U) :=
  ((contDiffOn_D hT hU hw).clm_apply contDiffOn_const).congr
    (fun p hp => dtW_eq hT hU hw hp)

theorem contDiffOn_slice_fderiv (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) :
    ContDiffOn ℝ ∞ (fun p : ℝ × E => fderiv ℝ (fun z => w (p.1, z)) p.2) (Icc 0 T ×ˢ U) :=
  (((ContinuousLinearMap.compL ℝ E (ℝ × E) W).flip
    (ContinuousLinearMap.inr ℝ ℝ E)).contDiff.comp_contDiffOn (contDiffOn_D hT hU hw)).congr
    (fun p hp => slice_eq hT hU hw hp)

theorem hasDerivWithinAt_slice_fderiv (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U) :
    HasDerivWithinAt (fun s => fderiv ℝ (fun z => w (s, z)) y)
      (fderiv ℝ (fun z => dtW w T (t, z)) y) (Icc 0 T) t := by
  set S : Set (ℝ × E) := Icc 0 T ×ˢ U with hSdef
  have hS : UniqueDiffOn ℝ S := uniqueDiffOn_S hT hU
  have hD := contDiffOn_D hT hU hw
  have hmem : (t, y) ∈ S := ⟨ht, hy⟩
  -- symmetry of the second derivative
  have hclos : (t, y) ∈ closure (interior S) := by
    rw [hSdef, interior_prod_eq, interior_Icc, hU.interior_eq, closure_prod_eq,
      closure_Ioo hT.ne]
    exact ⟨ht, subset_closure hy⟩
  have hsym : IsSymmSndFDerivWithinAt ℝ w S (t, y) :=
    (hw (t, y) hmem).isSymmSndFDerivWithinAt (by simpa using two_le_infty) hS hclos hmem
  -- LHS
  have hDd : HasFDerivWithinAt (fderivWithin ℝ w S) (fderivWithin ℝ (fderivWithin ℝ w S) S (t, y))
      S (t, y) := (hD.differentiableOn (by simp) _ hmem).hasFDerivWithinAt
  have hcurve := hDd.comp_hasDerivWithinAt (f := fun s => (s, y)) t
    ((hasDerivWithinAt_id t (Icc 0 T)).prodMk (hasDerivAt_const t y).hasDerivWithinAt)
    (fun s hs => ⟨hs, hy⟩)
  have hpre := ((ContinuousLinearMap.compL ℝ E (ℝ × E) W).flip
    (ContinuousLinearMap.inr ℝ ℝ E)).hasFDerivAt.comp_hasDerivWithinAt t hcurve
  have hL : HasDerivWithinAt (fun s => fderiv ℝ (fun z => w (s, z)) y)
      (((ContinuousLinearMap.compL ℝ E (ℝ × E) W).flip (ContinuousLinearMap.inr ℝ ℝ E))
        (fderivWithin ℝ (fderivWithin ℝ w S) S (t, y) (1, 0))) (Icc 0 T) t := by
    refine hpre.congr (fun s hs => ?_) ?_
    · exact slice_eq hT hU hw (p := (s, y)) ⟨hs, hy⟩
    · exact slice_eq hT hU hw (p := (t, y)) hmem
  -- RHS
  have hDy : HasFDerivAt (fun z => fderivWithin ℝ w S (t, z))
      ((fderivWithin ℝ (fderivWithin ℝ w S) S (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y :=
    (hDd.comp y (hasFDerivAt_prodMk_right t y).hasFDerivWithinAt
      (fun z hz => ⟨ht, hz⟩)).hasFDerivAt (hU.mem_nhds hy)
  have h2 := hDy.clm_apply (hasFDerivAt_const ((1 : ℝ), (0 : E)) y)
  have hev : (fun z => dtW w T (t, z)) =ᶠ[𝓝 y] fun z => fderivWithin ℝ w S (t, z) (1, 0) := by
    filter_upwards [hU.mem_nhds hy] with z hz
    exact dtW_eq hT hU hw (p := (t, z)) ⟨ht, hz⟩
  rw [hev.fderiv_eq, h2.fderiv]
  refine hL.congr_deriv ?_
  ext a
  simp only [ContinuousLinearMap.flip_apply, ContinuousLinearMap.compL_apply,
    ContinuousLinearMap.coe_comp', Function.comp_apply, ContinuousLinearMap.inr_apply,
    ContinuousLinearMap.add_apply, ContinuousLinearMap.zero_apply, map_zero,
    ContinuousLinearMap.comp_apply, zero_add]
  exact (hsym.eq (1, 0) (0, a))

theorem fderivWithin_prod_apply (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U)
    (τ : ℝ) (a : E) :
    fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y) (τ, a) =
      τ • dtW w T (t, y) + fderiv ℝ (fun z => w (t, z)) y a := by
  have h1 := dtW_eq hT hU hw (p := (t, y)) ⟨ht, hy⟩
  have h2 := slice_eq hT hU hw (p := (t, y)) ⟨ht, hy⟩
  simp only at h2
  rw [h1, h2]
  have : ((τ, a) : ℝ × E) = τ • ((1 : ℝ), (0 : E)) + (0, a) := by simp
  rw [this, map_add, map_smul]
  simp

theorem hasFDerivWithinAt_prod (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U) :
    HasFDerivWithinAt w (fderivWithin ℝ w (Icc 0 T ×ˢ U) (t, y)) (Icc 0 T ×ˢ U) (t, y) :=
  (hw.differentiableOn (by simp) _ ⟨ht, hy⟩).hasFDerivWithinAt

section Scalar

variable {w : ℝ × E → ℝ}

theorem hasDerivWithinAt_graphFu_time (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U)
    {Ψ : E × (E →L[ℝ] ℝ) → V}
    (hΨ : DifferentiableAt ℝ Ψ (graphMap (fun z => w (t, z)) y)) :
    HasDerivWithinAt (fun s => graphFu Ψ (fun z => w (s, z)) y)
      (fderiv ℝ Ψ (graphMap (fun z => w (t, z)) y)
        (0, fderiv ℝ (fun z => dtW w T (t, z)) y)) (Icc 0 T) t := by
  have hg : HasDerivWithinAt (fun s => graphMap (fun z => w (s, z)) y)
      ((0 : E), fderiv ℝ (fun z => dtW w T (t, z)) y) (Icc 0 T) t :=
    (hasDerivAt_const t y).hasDerivWithinAt.prodMk (hasDerivWithinAt_slice_fderiv hT hU hw ht hy)
  exact hΨ.hasFDerivAt.comp_hasDerivWithinAt t hg

theorem contDiffOn_graphFu_time (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {Ψ : E × (E →L[ℝ] ℝ) → V}
    {Ω : Set (E × (E →L[ℝ] ℝ))} (hΩ : IsOpen Ω) (hΨ : ContDiffOn ℝ ∞ Ψ Ω)
    (hmem : ∀ p ∈ Icc 0 T ×ˢ U, graphMap (fun z => w (p.1, z)) p.2 ∈ Ω) :
    ContDiffOn ℝ ∞ (fun p : ℝ × E => graphFu Ψ (fun z => w (p.1, z)) p.2) (Icc 0 T ×ˢ U) := by
  have h3 := contDiffOn_slice_fderiv hT hU hw
  exact hΨ.comp (contDiffOn_snd.prodMk h3) hmem

theorem contDiffWithinAt_graphFu_time (hT : 0 < T) (hU : IsOpen U)
    (hw : ContDiffOn ℝ ∞ w (Icc 0 T ×ˢ U)) {t : ℝ} {y : E} (ht : t ∈ Icc 0 T) (hy : y ∈ U)
    {Ψ : E × (E →L[ℝ] ℝ) → V} (hΨ : ContDiffAt ℝ ∞ Ψ (graphMap (fun z => w (t, z)) y)) :
    ContDiffWithinAt ℝ ∞ (fun p : ℝ × E => graphFu Ψ (fun z => w (p.1, z)) p.2)
      (Icc 0 T ×ˢ U) (t, y) := by
  have h3 := (contDiffOn_slice_fderiv hT hU hw) (t, y) ⟨ht, hy⟩
  exact ContDiffAt.comp_contDiffWithinAt (g := Ψ) (t, y) hΨ (contDiffWithinAt_snd.prodMk h3)

end Scalar

end CWFlow
