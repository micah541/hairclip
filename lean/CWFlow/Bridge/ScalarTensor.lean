import CWFlow.Bridge.Realizable
import ShortTime.Nonlinear.NemytskiiParamFibre

/-!
# A smooth function as a `(0,0)`-tensor field

`scalarTensor g f hf` is the tensor field `f • 1` (`unit00`), with `fnOf g (scalarTensor g f hf) = f`.
`scalarFam g f t` does this for a family `f : ℝ → M → ℝ` (zero at the times where `f t` is not smooth), and a
family that is jointly smooth on `S × M` gives a `JointlySmoothOn` family (`jointlySmoothOn_scalarFam`).
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped Manifold Topology ContDiff

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]

/-- `p ↦ c(p) • 1` over a smooth map `pr : P → M` is a smooth map into the bundle of `(0,0)`-tensors. -/
theorem contMDiffOn_smul_unit00 {EP : Type*} [NormedAddCommGroup EP] [NormedSpace ℝ EP]
    {HP : Type*} [TopologicalSpace HP] {IP : ModelWithCorners ℝ EP HP}
    {P : Type*} [TopologicalSpace P] [ChartedSpace HP P]
    {pr : P → M} {c : P → ℝ} {s : Set P} (hpr : ContMDiffOn IP 𝓘(ℝ, E) ∞ pr s)
    (hc : ContMDiffOn IP 𝓘(ℝ, ℝ) ∞ c s) :
    ContMDiffOn IP ((𝓘(ℝ, E)).prod 𝓘(ℝ, TensorRSModel 0 0 ℝ E)) ∞
      (fun p : P => TotalSpace.mk' (TensorRSModel 0 0 ℝ E)
        (E := fun x : M => TensorRSSpace 0 0 𝓘(ℝ, E) x) (pr p) (c p • unit00 (pr p))) s := by
  intro p₀ hp₀
  have hbase : pr p₀ ∈ (trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).baseSet :=
    mem_baseSet_trivializationAt _ _ (pr p₀)
  have hon : ContMDiffOn IP ((𝓘(ℝ, E)).prod 𝓘(ℝ, TensorRSModel 0 0 ℝ E)) ∞
      (fun p : P => TotalSpace.mk' (TensorRSModel 0 0 ℝ E)
        (E := fun x : M => TensorRSSpace 0 0 𝓘(ℝ, E) x) (pr p) (c p • unit00 (pr p)))
      (s ∩ pr ⁻¹' (trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).baseSet) := by
    rw [(trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).contMDiffOn_iff fun p hp =>
      (trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).mem_source.mpr hp.2]
    refine ⟨hpr.mono inter_subset_left,
      ((hc.mono inter_subset_left).smul (contMDiffOn_const (c := unit00M E))).congr fun p hp => ?_⟩
    have hlin := ((trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).linear ℝ hp.2).map_smul (c p) (unit00 (pr p))
    have hu := trivAt_unit00 (I := 𝓘(ℝ, E)) (M := M) (x₀ := pr p₀) (y := pr p)
      (show pr p ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) (pr p₀)).baseSet ∩
        (trivializationAt E (TangentSpace 𝓘(ℝ, E)) (pr p₀)).baseSet from hp.2).1
    exact hlin.trans (congrArg (fun z => c p • z) hu)
  exact (hon p₀ ⟨hp₀, hbase⟩).mono_of_mem_nhdsWithin
    (inter_mem self_mem_nhdsWithin ((hpr.continuousOn p₀ hp₀).preimage_mem_nhdsWithin
      ((trivAt 𝓘(ℝ, E) M 0 0 (pr p₀)).open_baseSet.mem_nhds hbase)))

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

/-- The `(0,0)`-tensor field `f • 1` of a smooth function `f`. -/
def scalarTensor (f : M → ℝ) (hf : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ f) : SmoothCcTensor g 0 0 where
  toSection :=
    { toFun := fun x => f x • unit00 x
      contMDiff_toFun := contMDiffOn_univ.mp
        (contMDiffOn_smul_unit00 (pr := fun x : M => x) contMDiffOn_id hf.contMDiffOn) }
  hasCompactSupport := HasCompactSupport.of_compactSpace _

theorem scalarTensor_toSection (f : M → ℝ) (hf : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ f) (x : M) :
    (scalarTensor g f hf).toSection x = f x • unit00 x := rfl

theorem lam_unit00 (x : M) : lam (I := 𝓘(ℝ, E)) x (unit00 x) = 1 := by
  show Tensor0SSpace.toModel (unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x) (fun i => i.elim0) = 1
  rw [unitZeroSec_apply (I := 𝓘(ℝ, E)) (M := M) x, Tensor0SSpace.toModel_ofModel,
    ContinuousMultilinearMap.constOfIsEmpty_apply]

theorem fnOf_scalarTensor (f : M → ℝ) (hf : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ f) :
    fnOf g (scalarTensor g f hf) = f := by
  funext x
  show lam (I := 𝓘(ℝ, E)) x (f x • unit00 x) = f x
  rw [map_smul, lam_unit00, smul_eq_mul, mul_one]

open Classical in
/-- The family of tensor fields of a family of functions (zero where the function is not smooth). -/
def scalarFam (f : ℝ → M → ℝ) (t : ℝ) : SmoothCcTensor g 0 0 :=
  if h : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t) then scalarTensor g (f t) h else 0

theorem scalarFam_eq {f : ℝ → M → ℝ} {t : ℝ} (h : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t)) :
    scalarFam g f t = scalarTensor g (f t) h := by
  unfold scalarFam
  rw [dif_pos h]

theorem fnOf_scalarFam {f : ℝ → M → ℝ} {t : ℝ} (h : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t)) :
    fnOf g (scalarFam g f t) = f t := by
  rw [scalarFam_eq g h, fnOf_scalarTensor]

/-- The slices of a jointly smooth family of functions are smooth. -/
theorem contMDiff_slice_fn {f : ℝ → M → ℝ} {S : Set ℝ}
    (hf : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => f q.1 q.2) (S ×ˢ univ))
    {t : ℝ} (ht : t ∈ S) : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t) := by
  intro x
  have hc : ContMDiffAt 𝓘(ℝ, E) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞ (fun y : M => ((t, y) : ℝ × M)) x :=
    contMDiffAt_const.prodMk contMDiffAt_id
  have h := ContMDiffWithinAt.comp (g := fun q : ℝ × M => f q.1 q.2)
    (f := fun y : M => ((t, y) : ℝ × M)) (s := univ) x (hf (t, x) ⟨ht, mem_univ _⟩)
    hc.contMDiffWithinAt (fun y _ => ⟨ht, mem_univ _⟩)
  exact contMDiffWithinAt_univ.mp h

/-- **A jointly smooth family of functions gives a jointly smooth family of tensor fields.** -/
theorem jointlySmoothOn_scalarFam {f : ℝ → M → ℝ} {S : Set ℝ}
    (hf : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => f q.1 q.2) (S ×ˢ univ)) :
    JointlySmoothOn g 0 0 S (scalarFam g f) := by
  have hc : ContMDiffOn ((𝓘(ℝ, E)).prod 𝓘(ℝ, ℝ)) 𝓘(ℝ, ℝ) ∞ (fun q : M × ℝ => f q.2 q.1)
      (univ ×ˢ S) :=
    ContMDiffOn.comp (g := fun q : ℝ × M => f q.1 q.2) (f := fun q : M × ℝ => ((q.2, q.1) : ℝ × M))
      hf (contMDiffOn_snd.prodMk contMDiffOn_fst) (fun q hq => ⟨hq.2, hq.1⟩)
  have h := contMDiffOn_smul_unit00 (M := M) (pr := fun q : M × ℝ => q.1)
    (c := fun q : M × ℝ => f q.2 q.1) (s := (univ : Set M) ×ˢ S) contMDiffOn_fst hc
  unfold JointlySmoothOn
  refine h.congr fun q hq => ?_
  rw [scalarFam_eq g (contMDiff_slice_fn hf hq.2)]
  rfl

end Bridge

end CWFlow
