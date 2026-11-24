import CWFlow.Bridge.ScalarTensor
import CWFlow.Bridge.Existence

/-!
# A family of functions solving the scalar equation is a Mantegazza–Martinazzi solution

`isSolution_scalarFam`: let `f : ℝ → M → ℝ` be jointly smooth on `[0, T] × M`, `T ≤ 1`, such that the graph of
`df_t` lies in the Weinstein domain of every patch chart on the support of its bump function, and
`∂ₜf_t = −div(JH)` of the graph. Then `scalarFam g f` is a solution of the Mantegazza–Martinazzi equation
`cwData` with initial datum `scalarFam g f 0`.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)

theorem toModel_smul00 (x : M) (c : ℝ) (ξ : TensorRSSpace 0 0 𝓘(ℝ, E) x) :
    TensorRSSpace.toModel (c • ξ) = c • TensorRSSpace.toModel ξ := rfl

theorem eq_lam_smul_unit00 (x : M) (Ξ : TensorRSSpace 0 0 𝓘(ℝ, E) x) :
    Ξ = lam (I := 𝓘(ℝ, E)) x Ξ • unit00 x := by
  apply lam_injective (I := 𝓘(ℝ, E)) x
  rw [map_smul, lam_unit00, smul_eq_mul, mul_one]

theorem isSolution_scalarFam {f : ℝ → M → ℝ} {T : ℝ} (hT : 0 < T) (hT1 : T ≤ 1)
    (hf : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => f q.1 q.2)
      (Icc 0 T ×ˢ univ))
    (hdom : ∀ t ∈ Icc (0 : ℝ) T, ∀ i : PIdx E M, ∀ x ∈ tsupport (patchρ E M i),
      chart E i.1 x ∈ graphDom E J F (f t) i.1)
    (heq : ∀ t ∈ Icc (0 : ℝ) T, ∀ x α : M, x ∈ patchU E α →
      chart E α x ∈ graphDom E J F (f t) α →
      HasDerivWithinAt (fun s => f s x)
        (-flowPot J (chartGraph E J F (f t) α) (chart E α x)) (Icc 0 T) t) :
    (cwData g hJi hJJ hF).toEquationData.IsSolution g (scalarFam g f 0) T (scalarFam g f) := by
  have hsl : ∀ t ∈ Icc (0 : ℝ) T, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t) := fun t ht =>
    contMDiff_slice_fn hf ht
  have hjets : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M,
      covJet g 0 0 (2 * 2) (scalarFam g f t) x ∈ mmDom g J F t x := by
    intro t ht x i hi
    unfold good
    rw [cJet_covJet g _ (tsupport_patchρ_subset i hi), fnOf_scalarFam g (hsl t ht)]
    exact (hdom t ht i x hi).2
  refine ⟨hT, hT1, jointlySmoothOn_scalarFam g hf, rfl, hjets, fun t ht x => ?_⟩
  obtain ⟨i, hi⟩ := (patchρ E M).exists_pos_of_mem (mem_univ x)
  have hx : x ∈ tsupport (patchρ E M i) := subset_tsupport _ (Function.mem_support.2 hi.ne')
  have hα := tsupport_patchρ_subset i hx
  have hgα := hjets t ht x i hx
  have hlam := lam_cwData_F hJi hJJ hF (scalarFam g f t) t (hjets t ht x) hα hgα
  rw [fnOf_scalarFam g (hsl t ht)] at hlam
  have hd := (heq t ht x i.1 hα (hdom t ht i x hx)).smul_const
    (TensorRSSpace.toModel (unit00 (I := 𝓘(ℝ, E)) x))
  have hval : TensorRSSpace.toModel ((cwData g hJi hJJ hF).toEquationData.F t x
      (covJet g 0 0 (2 * 2) (scalarFam g f t) x)) =
      (-flowPot J (chartGraph E J F (f t) i.1) (chart E i.1 x)) •
        TensorRSSpace.toModel (unit00 (I := 𝓘(ℝ, E)) x) := by
    have h1 := eq_lam_smul_unit00 x
      ((cwData g hJi hJJ hF).F t x (covJet g 0 0 (2 * 2) (scalarFam g f t) x))
    rw [hlam] at h1
    show TensorRSSpace.toModel ((cwData g hJi hJJ hF).F t x
      (covJet g 0 0 (2 * 2) (scalarFam g f t) x)) = _
    rw [h1, toModel_smul00]
  rw [hval]
  refine hd.congr (fun s hs => ?_) ?_
  · rw [scalarFam_eq g (hsl s hs), scalarTensor_toSection, toModel_smul00]
  · rw [scalarFam_eq g (hsl t ht), scalarTensor_toSection, toModel_smul00]

end Bridge

end CWFlow
