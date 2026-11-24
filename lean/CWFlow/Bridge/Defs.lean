import CWFlow.Bridge.ChartJet
import CWFlow.Bridge.Contract
import CWFlow.Bridge.GraphChart

/-!
# The Mantegazza–Martinazzi data of the scalar equation: definitions

The scalar equation `∂ₜu = −div(JH)[F_u]` is written as MM data (`ShortTime.MMData _ M 2`) by gluing the chart
operators with a partition of unity.

* `patchU α ⊆ M`: the patch of MM's chart at `α` (`patchΩ α`, where the exact jet recursion holds), pulled back
  to `M`; `patchPts` finitely many patches covering `M`; `patchρ` a smooth partition of unity subordinate to them.
* For a covariant 4-jet `ζ` at `x` (MM's `CovJetFiber`): `rdJet α x ζ` its chart readings, `cJet α x ζ` the
  coordinate 4-jet (`coordJet`), `permTop` the slot permutation of the top coordinate jet which turns the pairing
  `(0,3),(1,2)` of `Pjet_add_top` into MM's pairing `(0,1),(2,3)` (harmless on symmetric jets), and
  `Fi α x ζ = Pjet(chartPsi α)(φ_α x, permTop (cJet α x ζ))` the chart operator.
* `mmDom t x`: the jets whose 1-jet lies in the Weinstein domain of every chart whose bump is nonzero near `x`;
  `mmEf t x ζ = Σ ρ_i(x) cometricAt α_i x (DF)`; `mmb t x ζ = Σ ρ_i(x) Fi α_i x (zeroTopC ζ)`.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow.Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-! ### Patches and the partition of unity -/

variable (E) in
/-- The patch of the chart at `α`, in `M`. -/
def patchU (α : M) : Set M := (chart E α).source ∩ chart E α ⁻¹' patchΩ (I := 𝓘(ℝ, E)) α

theorem isOpen_patchU (α : M) : IsOpen (patchU E α) :=
  (continuousOn_extChartAt (I := 𝓘(ℝ, E)) α).isOpen_inter_preimage
    (isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α) (isOpen_patchΩ (I := 𝓘(ℝ, E)) α)

theorem mem_patchU_self (α : M) : α ∈ patchU E α :=
  ⟨mem_extChartAt_source (I := 𝓘(ℝ, E)) α, mem_patchΩ_self (I := 𝓘(ℝ, E)) α⟩

theorem patchU_subset_source (α : M) : patchU E α ⊆ (chart E α).source := inter_subset_left

variable (E M) in
theorem exists_finite_patches : ∃ s : Finset M, ∀ x : M, ∃ α ∈ s, x ∈ patchU E α := by
  obtain ⟨s, hs⟩ := isCompact_univ.elim_finite_subcover (fun α : M => patchU E α)
    (fun α => isOpen_patchU α) (fun x _ => mem_iUnion.2 ⟨x, mem_patchU_self x⟩)
  refine ⟨s, fun x => ?_⟩
  have hx := hs (mem_univ x)
  simp only [mem_iUnion] at hx
  obtain ⟨α, hα, hx⟩ := hx
  exact ⟨α, hα, hx⟩

variable (E M) in
/-- Finitely many patch centres whose patches cover `M`. -/
def patchPts : Finset M := (exists_finite_patches E M).choose

variable (E M) in
/-- The index type of the patches. -/
abbrev PIdx : Type _ := ↥(patchPts E M)

variable (E M) in
theorem exists_patchPOU : ∃ ρ : SmoothPartitionOfUnity (PIdx E M) 𝓘(ℝ, E) M univ,
    ρ.IsSubordinate (fun i => patchU E i.1) := by
  refine SmoothPartitionOfUnity.exists_isSubordinate 𝓘(ℝ, E) isClosed_univ _ (fun i => isOpen_patchU _)
    fun x _ => ?_
  obtain ⟨α, hα, hx⟩ := (exists_finite_patches E M).choose_spec x
  exact mem_iUnion.2 ⟨⟨α, hα⟩, hx⟩

variable (E M) in
/-- A smooth partition of unity subordinate to the patches. -/
def patchρ : SmoothPartitionOfUnity (PIdx E M) 𝓘(ℝ, E) M univ := (exists_patchPOU E M).choose

theorem tsupport_patchρ_subset (i : PIdx E M) : tsupport (patchρ E M i) ⊆ patchU E i.1 :=
  (exists_patchPOU E M).choose_spec i

/-! ### Jets in the charts -/

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) (J : V →L[ℝ] V) (F : M → V)

/-- The chart readings of a covariant 4-jet at `x`, in the chart at `α`. -/
def rdJet (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : JetModel E 0 0 4 :=
  EquationData.trivJet (I := 𝓘(ℝ, E)) α x ζ

/-- The coordinate 4-jet, in the chart at `α`, of a covariant 4-jet at `x`. -/
def cJet (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : SJet E :=
  coordJet g α (chart E α x) (rdJet α x ζ)

/-- `(0, 3, 1, 2)`: the slot permutation of the top jet. -/
def σ4 : Fin 4 ≃ Fin 4 where
  toFun := ![0, 3, 1, 2]
  invFun := ![0, 2, 3, 1]
  left_inv := by decide
  right_inv := by decide

/-- The slot permutation of the top coordinate jet: `T ↦ (m ↦ T ![m 0, m 3, m 1, m 2])`. -/
def permTop (z : SJet E) : SJet E :=
  (z.1, z.2.1, z.2.2.1, z.2.2.2.domDomCongr σ4)

/-- The chart operator: `−div(JH)` of the graph, as a function of the covariant jet. -/
def Fi (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : ℝ :=
  Pjet J (chartPsi E J F α) (chart E α x, permTop (cJet g α x ζ))

/-- The differential of the graph immersion at `x`, in the chart at `α`. -/
def Li (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : E →L[ℝ] V :=
  LOf (chartPsi E J F α) (chart E α x) (cJet g α x ζ)

/-- The 1-jet lies in the Weinstein domain of the chart at `α`. -/
def good (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : Prop :=
  (chart E α x, (cJet g α x ζ).1) ∈ chartDom E J F α

/-- The jet with its top component replaced by `0`. -/
def zeroTopC {x : M} (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x :=
  Function.update ζ (Fin.last (2 * 2)) 0

/-- The unit `(0,0)`-tensor. -/
def unitT (x : M) : TensorRSSpace 0 0 𝓘(ℝ, E) x :=
  TensorRSSpace.ofCLM (ContinuousLinearMap.id ℝ (Tensor0SSpace 0 𝓘(ℝ, E) x))

/-! ### The MM data -/

/-- The domain: the 1-jet lies in the Weinstein domain of each chart whose bump function is nonzero near `x`. -/
def mmDom (_t : ℝ) (x : M) : Set (CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :=
  {ζ | ∀ i : PIdx E M, x ∈ tsupport (patchρ E M i) → good g J F i.1 x ζ}

/-- The principal coefficient: the glued cometric of the graph. -/
def mmEf (_t : ℝ) (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :
    (TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) →ₗ[ℝ] (TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) →ₗ[ℝ] ℝ :=
  ∑ i : PIdx E M, patchρ E M i x • cometricAt i.1 x (Li g J F i.1 x ζ)

/-- The scalar lower-order term. -/
def mmbS (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : ℝ :=
  ∑ i : PIdx E M, patchρ E M i x * Fi g J F i.1 x (zeroTopC ζ)

/-- The lower-order term. -/
def mmb (_t : ℝ) (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) : TensorRSSpace 0 0 𝓘(ℝ, E) x :=
  mmbS g J F x ζ • unitT x

/-! ### Basic facts -/

theorem permTop_apply_top (z : SJet E) (m : Fin 4 → E) :
    (permTop z).2.2.2 m = z.2.2.2 ![m 0, m 3, m 1, m 2] := by
  simp only [permTop, ContinuousMultilinearMap.domDomCongr_apply]
  congr 1
  ext i
  fin_cases i <;> rfl

theorem lam_unitT (x : M) : lam (I := 𝓘(ℝ, E)) x (unitT x) = 1 := by
  show Tensor0SSpace.toModel (unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x) (fun i => i.elim0) = 1
  rw [toModel_unitZeroSec]
  rfl

theorem lam_mmb (t : ℝ) (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :
    lam (I := 𝓘(ℝ, E)) x (mmb g J F t x ζ) = mmbS g J F x ζ := by
  rw [mmb, map_smul, lam_unitT, smul_eq_mul, mul_one]

end CWFlow.Bridge
