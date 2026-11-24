import CWFlow.Bridge.Defs

/-!
# The Mantegazza–Martinazzi data of the scalar equation: the algebraic fields

For the data of `Bridge/Defs.lean`:

* `injective_LOf`: the differential of the graph, as a function of the jet, is injective at every jet whose
  1-jet lies in the Weinstein domain (the second-order part need not be symmetric);
* `mmEf_lower`, `mmb_lower`: `E` and `b` do not depend on the top component of the jet;
* `mmEf_symm`, `mmEf_pos`: `E` is symmetric, and positive definite on the domain.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

section LOf

variable {E V : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {Ψ : E × (E →L[ℝ] ℝ) → V} {O : Set (E × (E →L[ℝ] ℝ))}

/-- **`DF` of the jet is injective** wherever the 1-jet lies in the domain of a symplectic `Ψ`; no symmetry of
the second-order part is needed. -/
theorem injective_LOf
    (hs : ∀ q ∈ O, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {y : E} {z : SJet E} (hq : (y, z.1) ∈ O) : Function.Injective (LOf Ψ y z) := by
  rw [injective_iff_map_eq_zero]
  intro a ha
  refine SeparatingDual.eq_zero_of_forall_dual_eq_zero (R := ℝ) fun ζ => ?_
  have h := inner_J_fderiv_vertical (hs _ hq) ζ (a, slot1L z.2.1 a)
  rw [← LOf_apply, ha, inner_zero_right] at h
  exact h.symm

/-- `LOf` depends only on the first two components of the jet. -/
theorem LOf_congr {y : E} {z z' : SJet E} (h1 : z.1 = z'.1) (h2 : z.2.1 = z'.2.1) :
    LOf Ψ y z = LOf Ψ y z' := by
  ext a
  rw [LOf_apply, LOf_apply, h1, h2]

end LOf

namespace Bridge


variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) (J : V →L[ℝ] V) (F : M → V)

/-! ### Dependence on the lower part of the jet -/

theorem rdJet_lower (α x : M) {ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x}
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) :
    ∀ k, k ≠ Fin.last 4 → rdJet α x ζ k = rdJet α x ζ' k := by
  intro k hk
  show (trivAt 𝓘(ℝ, E) M 0 (0 + k) α).continuousLinearMapAt ℝ x (ζ k) =
    (trivAt 𝓘(ℝ, E) M 0 (0 + k) α).continuousLinearMapAt ℝ x (ζ' k)
  rw [h k hk]

theorem cJet_lower (α x : M) {ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x}
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) :
    (cJet g α x ζ).1 = (cJet g α x ζ').1 ∧ (cJet g α x ζ).2.1 = (cJet g α x ζ').2.1 ∧
      (cJet g α x ζ).2.2.1 = (cJet g α x ζ').2.2.1 :=
  coordJet_lower g α _ (rdJet_lower α x h)

theorem Li_lower (α x : M) {ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x}
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) : Li g J F α x ζ = Li g J F α x ζ' :=
  LOf_congr (cJet_lower g α x h).1 (cJet_lower g α x h).2.1

theorem good_lower (α x : M) {ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x}
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) : good g J F α x ζ ↔ good g J F α x ζ' := by
  unfold good
  rw [(cJet_lower g α x h).1]

theorem zeroTopC_apply_of_ne {x : M} (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) {k : Fin (2 * 2 + 1)}
    (hk : k ≠ Fin.last (2 * 2)) : zeroTopC ζ k = ζ k :=
  Function.update_of_ne hk _ _

theorem zeroTopC_last {x : M} (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :
    zeroTopC ζ (Fin.last (2 * 2)) = 0 :=
  Function.update_self _ _ _

theorem zeroTopC_congr {x : M} {ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x}
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) : zeroTopC ζ = zeroTopC ζ' := by
  funext k
  by_cases hk : k = Fin.last (2 * 2)
  · subst hk
    rw [zeroTopC_last, zeroTopC_last]
  · rw [zeroTopC_apply_of_ne ζ hk, zeroTopC_apply_of_ne ζ' hk, h k hk]

theorem good_zeroTopC (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :
    good g J F α x (zeroTopC ζ) ↔ good g J F α x ζ :=
  good_lower g J F α x fun _ hk => zeroTopC_apply_of_ne ζ hk

theorem mmEf_lower (t : ℝ) (x : M) (ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x)
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) : mmEf g J F t x ζ = mmEf g J F t x ζ' := by
  unfold mmEf
  exact Finset.sum_congr rfl fun i _ => by rw [Li_lower g J F i.1 x h]

theorem mmb_lower (t : ℝ) (x : M) (ζ ζ' : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x)
    (h : ∀ k, k ≠ Fin.last (2 * 2) → ζ k = ζ' k) : mmb g J F t x ζ = mmb g J F t x ζ' := by
  unfold mmb mmbS
  rw [zeroTopC_congr h]

/-! ### Symmetry and positivity -/

theorem mmEf_apply (t : ℝ) (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x)
    (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    mmEf g J F t x ζ ξ η =
      ∑ i : PIdx E M, patchρ E M i x * cometricAt i.1 x (Li g J F i.1 x ζ) ξ η := by
  simp only [mmEf, LinearMap.sum_apply, LinearMap.smul_apply, smul_eq_mul]

theorem mmEf_symm (t : ℝ) (x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x)
    (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) : mmEf g J F t x ζ ξ η = mmEf g J F t x ζ η ξ := by
  rw [mmEf_apply, mmEf_apply]
  exact Finset.sum_congr rfl fun i _ => by rw [cometricAt_symm]

theorem mem_source_of_mem_tsupport {i : PIdx E M} {x : M} (hx : x ∈ tsupport (patchρ E M i)) :
    x ∈ (chartAt E i.1).source := by
  have h := patchU_subset_source i.1 (tsupport_patchρ_subset i hx)
  rwa [extChartAt_source] at h

variable {g J F}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

/-- On good jets the differential of the graph is injective. -/
theorem injective_Li {α x : M} {ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x} (h : good g J F α x ζ) :
    Function.Injective (Li g J F α x ζ) :=
  injective_LOf (chartPsi_symplectic hJi hJJ hF α) h

theorem mmEf_pos {t : ℝ} {x : M} {ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x} (hζ : ζ ∈ mmDom g J F t x)
    {ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ} (hξ : ξ ≠ 0) : 0 < mmEf g J F t x ζ ξ ξ := by
  rw [mmEf_apply]
  have hterm : ∀ i : PIdx E M, 0 < patchρ E M i x →
      0 < patchρ E M i x * cometricAt i.1 x (Li g J F i.1 x ζ) ξ ξ := by
    intro i hi
    have hx : x ∈ tsupport (patchρ E M i) := subset_tsupport _ (Function.mem_support.2 hi.ne')
    exact mul_pos hi (cometricAt_pos (mem_source_of_mem_tsupport hx)
      (injective_Li hJi hJJ hF (hζ i hx)) hξ)
  obtain ⟨i₀, hi₀⟩ := (patchρ E M).exists_pos_of_mem (mem_univ x)
  refine Finset.sum_pos' (fun i _ => ?_) ⟨i₀, Finset.mem_univ _, hterm i₀ hi₀⟩
  rcases ((patchρ E M).nonneg i x).eq_or_lt with h0 | hpos
  · rw [← h0, zero_mul]
  · exact (hterm i hpos).le

end Bridge

end CWFlow
