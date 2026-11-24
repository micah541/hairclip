import CWFlow.Flow.ShortUnique

/-!
# Uniqueness of the flow `∂ₜF = J∇div(JH)` of Lagrangian immersions

`flow_unique`: two flows (`IsFlow`) on `[0, T]` of a compact manifold with the same initial immersion coincide
on `[0, T]`.

Proof: the set of times where the two flows agree is closed, contains `0`, and is open to the right by
short-time uniqueness (`flow_unique_short`) applied to the flows restarted at that time (`IsFlow.shift`).
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

section Shift

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- A flow restarted at the time `a`. -/
theorem IsFlow.shift {J : V →L[ℝ] V} {T a : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) (ha : 0 ≤ a) :
    IsFlow E J (T - a) (fun s => Φ (a + s)) where
  smooth := by
    have h1 : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => a + q.1) :=
      ContMDiff.comp (g := fun s : ℝ => a + s) (f := fun q : ℝ × M => q.1)
        (contDiff_const.add contDiff_id).contMDiff contMDiff_fst
    have h2 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((a + q.1, q.2) : ℝ × M)) (Icc 0 (T - a) ×ˢ univ) :=
      h1.contMDiffOn.prodMk contMDiffOn_snd
    exact ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2)
      (f := fun q : ℝ × M => ((a + q.1, q.2) : ℝ × M)) h.smooth h2
      (fun q hq => ⟨⟨by linarith [hq.1.1], by linarith [hq.1.2]⟩, mem_univ _⟩)
  imm := fun t ht => h.imm (a + t) ⟨by linarith [ht.1], by linarith [ht.2]⟩
  vel := fun t ht x z hx => by
    have h1 : HasDerivWithinAt (fun s : ℝ => a + s) 1 (Icc 0 (T - a)) t :=
      (hasDerivWithinAt_id t _).const_add a
    have h2 := HasDerivWithinAt.scomp t
      (h.vel (a + t) ⟨by linarith [ht.1], by linarith [ht.2]⟩ x z hx) h1
      (fun s hs => ⟨by linarith [hs.1], by linarith [hs.2]⟩)
    rw [one_smul] at h2
    exact h2

end Shift

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- **Uniqueness of the flow `∂ₜF = J∇div(JH)`.** Two solutions on `[0, T]` (jointly smooth families of
Lagrangian immersions of a compact manifold moving by `J∇div(JH)`) with the same initial immersion
coincide on `[0, T]`. -/
theorem flow_unique (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {T : ℝ} {Φ₁ Φ₂ : ℝ → M → V}
    (hΦ₁ : IsFlow E J T Φ₁) (hΦ₂ : IsFlow E J T Φ₂) (h0 : Φ₁ 0 = Φ₂ 0) :
    ∀ t ∈ Icc (0 : ℝ) T, Φ₁ t = Φ₂ t := by
  have hcont : ∀ {Φ : ℝ → M → V}, IsFlow E J T Φ → ∀ x : M,
      ContinuousOn (fun t => Φ t x) (Icc 0 T) := fun hΦ x t ht =>
    (hΦ.vel t ht x x (mem_extChartAt_source (I := 𝓘(ℝ, E)) x)).continuousWithinAt
  have hsub : Icc (0 : ℝ) T ⊆ {t | Φ₁ t = Φ₂ t} := by
    refine IsClosed.Icc_subset_of_forall_mem_nhdsWithin ?_ h0 ?_
    · have hcl : ∀ x : M, IsClosed (Icc (0 : ℝ) T ∩ (fun t => (Φ₁ t x, Φ₂ t x)) ⁻¹' diagonal V) :=
        fun x => ((hcont hΦ₁ x).prodMk (hcont hΦ₂ x)).preimage_isClosed_of_isClosed isClosed_Icc
          isClosed_diagonal
      have he : {t | Φ₁ t = Φ₂ t} ∩ Icc (0 : ℝ) T =
          Icc (0 : ℝ) T ∩ ⋂ x : M, (Icc (0 : ℝ) T ∩ (fun t => (Φ₁ t x, Φ₂ t x)) ⁻¹' diagonal V) := by
        ext t
        constructor
        · rintro ⟨ht, hI⟩
          exact ⟨hI, mem_iInter.2 fun x => ⟨hI, congrFun ht x⟩⟩
        · rintro ⟨hI, hx⟩
          exact ⟨funext fun x => (mem_iInter.1 hx x).2, hI⟩
      rw [he]
      exact isClosed_Icc.inter (isClosed_iInter hcl)
    · rintro a ⟨ha, haI⟩
      have ha' : Φ₁ a = Φ₂ a := ha
      have haT : a ∈ Icc (0 : ℝ) T := ⟨haI.1, haI.2.le⟩
      have hF := hΦ₁.imm a haT
      have i1 : (fun s => Φ₁ (a + s)) 0 = Φ₁ a := by
        show Φ₁ (a + 0) = Φ₁ a
        rw [add_zero]
      have i2 : (fun s => Φ₂ (a + s)) 0 = Φ₁ a := by
        show Φ₂ (a + 0) = Φ₁ a
        rw [add_zero, ha']
      obtain ⟨δ, hδ, -, hag⟩ := flow_unique_short g hJi hJJ hF hdim (sub_pos.2 haI.2)
        (hΦ₁.shift haI.1) (hΦ₂.shift haI.1) i1 i2
      rw [mem_nhdsGT_iff_exists_Ioc_subset]
      refine ⟨a + δ, lt_add_of_pos_right a hδ, fun s hs => ?_⟩
      have h : Φ₁ (a + (s - a)) = Φ₂ (a + (s - a)) :=
        hag (s - a) ⟨by linarith [hs.1], by linarith [hs.2]⟩
      have e : a + (s - a) = s := by ring
      rw [e] at h
      exact h
  exact fun t ht => hsub ht

end Bridge

end CWFlow
