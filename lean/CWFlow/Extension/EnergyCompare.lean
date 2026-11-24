import CWFlow.Extension.FlowM
import CWFlow.Extension.GraphBridge
import CWFlow.Interp.BallInterp
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.ContDiff
import Mathlib.LinearAlgebra.Determinant

/-!
# Comparing the global energy with chart integrals of the energy symbol

`energyM` is defined through a fixed atlas. Here it is compared with integrals in parametrizations:
`energyM_eq_sum` (partition of unity subordinate to parametrizations), `eDens_eq_eSym_of_graph`
(along a graph parametrization the energy density is the jet-space symbol `eSym`), and
`integral_eSym_le_energyM` (a chart integral of the symbol is bounded by the global energy).
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace Matrix

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {F : M → V}

/-- Monotonicity of the integral in the weight, for a nonnegative density. -/
theorem intM_densOf_mono [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) (hF : IsLagImm ι J F)
    {d} (hd : IsNatDens J (densOf d)) (hd0 : ∀ X y, 0 ≤ d X y) {h₁ h₂ : M → ℝ}
    (hh₁ : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h₁) (hh₂ : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h₂)
    (hle : ∀ p, h₁ p ≤ h₂ p) :
    intM ι A F h₁ (densOf d) ≤ intM ι A F h₂ (densOf d) := by
  have key := intM_densOf_sum ι A hF hd (Finset.univ : Finset (Fin 2))
    ![h₁, fun p => h₂ p - h₁ p] (by
      intro a _
      fin_cases a
      · exact hh₁
      · exact hh₂.sub hh₁)
  have hs : (fun p => ∑ a ∈ (Finset.univ : Finset (Fin 2)), (![h₁, fun p => h₂ p - h₁ p] : Fin 2 → M → ℝ) a p)
      = h₂ := funext fun p => by simp [Fin.sum_univ_two]
  rw [hs, Fin.sum_univ_two] at key
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one] at key
  have := intM_densOf_nonneg ι A (F := F) (h := fun p => h₂ p - h₁ p)
    (fun p => sub_nonneg.2 (hle p)) hd0
  linarith

/-- Along a graph parametrization the energy density is the jet-space energy symbol. -/
theorem eDens_eq_eSym_of_graph (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) {σ : EuclideanSpace ℝ (Fin n) → M} {B : Set (EuclideanSpace ℝ (Fin n))}
    (hB : IsOpen B) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ B)
    (hgraph : ∀ x ∈ B, F (σ x) = graphX J U y₀ φ x) {N l : ℕ} (hN : l + 3 ≤ N) {x} (hx : x ∈ B) :
    eDens J (F ∘ σ) l x = eSym l (jp N φ x) := by
  -- a global smooth function agreeing with `φ` near `x`
  obtain ⟨R, hR, hRB⟩ : ∃ R : ℝ, 0 < R ∧
      closedBall (0 : EuclideanSpace ℝ (Fin n)) R ⊆ (fun v => x + v) ⁻¹' B := by
    have hO : IsOpen ((fun v : EuclideanSpace ℝ (Fin n) => x + v) ⁻¹' B) :=
      hB.preimage (continuous_const.add continuous_id)
    obtain ⟨ε, hε, hεB⟩ := Metric.isOpen_iff.1 hO 0 (by simpa using hx)
    exact ⟨ε / 2, by positivity, (closedBall_subset_ball (by linarith)).trans hεB⟩
  have hO : IsOpen ((fun v : EuclideanSpace ℝ (Fin n) => x + v) ⁻¹' B) :=
    hB.preimage (continuous_const.add continuous_id)
  have hφ' : ContDiffOn ℝ ∞ (fun v => φ (x + v)) ((fun v => x + v) ⁻¹' B) :=
    hφ.comp (contDiff_const.add contDiff_id).contDiffOn fun v hv => hv
  obtain ⟨R', hR', g, hg, -, hgφ⟩ := exists_global_of_contDiffOn hO hR hRB hφ'
  set ψ : EuclideanSpace ℝ (Fin n) → ℝ := fun y => g (y - x) with hψ
  have hψc : ContDiff ℝ ∞ ψ := hg.comp (contDiff_id.sub contDiff_const)
  have hball : ball x R' ∈ 𝓝 x := ball_mem_nhds x (by linarith)
  have hψφ : ∀ y ∈ ball x R', ψ y = φ y := by
    intro y hy
    have : y - x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R' := by
      rw [mem_ball, dist_eq_norm] at hy ⊢
      simpa using hy
    simp only [hψ]
    rw [hgφ _ this]; simp
  have hev : φ =ᶠ[𝓝 x] ψ := by
    filter_upwards [hball] with y hy using (hψφ y hy).symm
  have hgrad : graphX J U y₀ φ =ᶠ[𝓝 x] graphX J U y₀ ψ := by
    have hev' : ∀ᶠ y in 𝓝 x, φ =ᶠ[𝓝 y] ψ := hev.eventually_nhds
    filter_upwards [hev'] with y hy
    unfold graphX gradE
    rw [hy.fderiv_eq]
  have hFσ : (F ∘ σ) =ᶠ[𝓝 x] graphX J U y₀ ψ := by
    filter_upwards [hB.mem_nhds hx, hgrad] with y hy hy2
    rw [Function.comp_apply, hgraph y hy, hy2]
  have hjp : jp N ψ x = jp N φ x := by
    unfold jp
    congr 1
    funext m
    exact jet_congr_of_eventuallyEq hev.symm _
  rw [eDens_congr hFσ l, eDens_graphX hJi hJJ hU hψc hN x, hjp]

/-- The energy as a sum over a partition subordinate to parametrizations. -/
theorem energyM_eq_sum [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) (hF : IsLagImm ι J F)
    {α : Type*} (s : Finset α) (h : α → M → ℝ) (hh : ∀ a ∈ s, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (h a))
    (hsum : ∀ p, ∑ a ∈ s, h a p = 1) (B : α → Set (EuclideanSpace ℝ (Fin n)))
    (σ : α → EuclideanSpace ℝ (Fin n) → M) (hσ : ∀ a ∈ s, IsParam ι (B a) (σ a))
    (hsupp : ∀ a ∈ s, tsupport (h a) ⊆ σ a '' B a) (l : ℕ) :
    energyM ι A J l F = ∑ a ∈ s, ∫ x in B a, h a (σ a x) * eDens J (F ∘ σ a) l x := by
  have h1 : (fun _ : M => (1 : ℝ)) = fun p => ∑ a ∈ s, h a p := funext fun p => (hsum p).symm
  have hnd := isNatDens_densOf_eDens (n := n) J l
  unfold energyM
  rw [h1, intM_densOf_sum ι A hF hnd s h hh]
  refine Finset.sum_congr rfl fun a ha => ?_
  rw [intM_eq_param ι A hF (hh a ha) hnd (hσ a ha) (hsupp a ha)]
  rfl

/-! ### Nonnegativity of the energy density -/

theorem ec_gX_posSemidef (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) :
    (gX X x).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun v => ?_
  · ext a b
    simp [Matrix.conjTranspose_apply, gX, real_inner_comm]
  · set w : V := ∑ a, v a • dX X x a with hw
    have key : star v ⬝ᵥ (gX X x *ᵥ v) = ⟪w, w⟫ := by
      simp only [hw, sum_inner, inner_sum, real_inner_smul_left, real_inner_smul_right,
        dotProduct, Matrix.mulVec, gX, Pi.star_apply, star_trivial, Finset.mul_sum]
      refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
      rw [real_inner_comm (dX X x a) (dX X x b)]
      ring
    rw [key]
    exact real_inner_self_nonneg

theorem ec_eDens_nonneg (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) (k : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) : 0 ≤ eDens J X k x := by
  unfold eDens
  have hP := (ec_gX_posSemidef X x).inv
  have hsymm : ∀ a b, ginvZ X x ![a, b] = ginvZ X x ![b, a] := by
    intro a b
    have := hP.1.apply a b
    simpa [ginvZ] using this.symm
  have hnn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, ginvZ X x ![a, b] * v a * v b := by
    intro v
    have := hP.dotProduct_mulVec_nonneg v
    simp only [dotProduct, Matrix.mulVec, Pi.star_apply, star_trivial, Finset.mul_sum] at this
    refine le_of_le_of_eq this ?_
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    simp [ginvZ]
    ring
  exact mul_nonneg (gNormSq_nonneg hsymm hnn _) (Real.sqrt_nonneg _)

/-! ### Extending a bump function from a parametrization -/

open Classical in
/-- `η ∘ κ` on the image of `σ`, zero elsewhere. -/
noncomputable def ec_ext (B : Set (EuclideanSpace ℝ (Fin n))) (σ : EuclideanSpace ℝ (Fin n) → M)
    (η : EuclideanSpace ℝ (Fin n) → ℝ) (κ : M → EuclideanSpace ℝ (Fin n)) (p : M) : ℝ :=
  if p ∈ σ '' B then η (κ p) else 0

theorem ec_ext_of_mem {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M}
    {η : EuclideanSpace ℝ (Fin n) → ℝ} {κ : M → EuclideanSpace ℝ (Fin n)} {p : M}
    (hp : p ∈ σ '' B) : ec_ext B σ η κ p = η (κ p) := by
  unfold ec_ext; rw [if_pos hp]

theorem ec_ext_of_notMem {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M}
    {η : EuclideanSpace ℝ (Fin n) → ℝ} {κ : M → EuclideanSpace ℝ (Fin n)} {p : M}
    (hp : p ∉ σ '' B) : ec_ext B σ η κ p = 0 := by
  unfold ec_ext; rw [if_neg hp]

theorem ec_contMDiff_ext [T2Space M] {B : Set (EuclideanSpace ℝ (Fin n))}
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    {K : Set (EuclideanSpace ℝ (Fin n))} (hK : IsCompact K) (hKB : K ⊆ B)
    {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη : ContDiff ℝ ∞ η) (hη0 : ∀ x, x ∉ K → η x = 0)
    {κ : M → EuclideanSpace ℝ (Fin n)} (hκc : ContinuousOn κ (σ '' B))
    (hκ : ∀ y ∈ B, κ (σ y) = y) :
    ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (ec_ext B σ η κ) := by
  intro p
  by_cases hp : p ∈ σ '' K
  · obtain ⟨x₀, hx₀K, rfl⟩ := hp
    have hx₀B : x₀ ∈ B := hKB hx₀K
    set p := σ x₀ with hpdef
    set g : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n) := chartFwd ι p ∘ σ with hg
    have hW := im_isOpen_W ι hσ p
    have hx₀W : x₀ ∈ B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) p).source := ⟨hx₀B, mem_extChartAt_source p⟩
    have hgat : ContDiffAt ℝ ∞ g x₀ := (hσ.trans_smooth p).contDiffAt (hW.mem_nhds hx₀W)
    have hdet : (fderiv ℝ g x₀).det ≠ 0 := by
      rw [im_det_eq_jacM]; exact hσ.trans_det p x₀ hx₀W
    set f' : EuclideanSpace ℝ (Fin n) ≃L[ℝ] EuclideanSpace ℝ (Fin n) :=
      (LinearMap.equivOfDetNeZero (fderiv ℝ g x₀ : EuclideanSpace ℝ (Fin n) →ₗ[ℝ]
        EuclideanSpace ℝ (Fin n)) hdet).toContinuousLinearEquiv with hf'
    have hf'c : (f' : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n)) = fderiv ℝ g x₀ := by
      ext w; rfl
    have hd : HasFDerivAt g (f' : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n)) x₀ := by
      rw [hf'c]; exact (hgat.differentiableAt (by simp)).hasFDerivAt
    have hn : (∞ : WithTop ℕ∞) ≠ 0 := by simp
    set ψ := ContDiffAt.localInverse hgat hd hn with hψdef
    have hψ : ContDiffAt ℝ ∞ ψ (g x₀) := ContDiffAt.to_localInverse hgat hd hn
    have hleft : ∀ᶠ x in 𝓝 x₀, ψ (g x) = x := (ContDiffAt.hasStrictFDerivAt' hgat hd hn).eventually_left_inverse
    have hκp : ContinuousAt κ p :=
      hκc.continuousAt (hσ.isOpen_image.mem_nhds ⟨x₀, hx₀B, rfl⟩)
    have hκx : κ p = x₀ := hκ x₀ hx₀B
    have h1 : ∀ᶠ q in 𝓝 p, ψ (g (κ q)) = κ q := by
      exact hκp.eventually (show ∀ᶠ x in 𝓝 (κ p), ψ (g x) = x by rw [hκx]; exact hleft)
    have h2 : ∀ᶠ q in 𝓝 p, q ∈ σ '' B := hσ.isOpen_image.mem_nhds ⟨x₀, hx₀B, rfl⟩
    have hev : ec_ext B σ η κ =ᶠ[𝓝 p] fun q => (η ∘ ψ) (chartFwd ι p q) := by
      filter_upwards [h1, h2] with q hq1 hq2
      rw [ec_ext_of_mem hq2]
      obtain ⟨y, hy, rfl⟩ := hq2
      have : g y = chartFwd ι p (σ y) := rfl
      rw [hκ y hy] at hq1 ⊢
      simp only [Function.comp_apply]
      rw [← this, hq1]
    have hchart : ContMDiffAt 𝓘(ℝ, E) 𝓘(ℝ, EuclideanSpace ℝ (Fin n)) ∞ (chartFwd ι p) p :=
      ι.symm.contDiff.contMDiff.contMDiffAt.comp p (contMDiffAt_extChartAt (I := 𝓘(ℝ, E)) (x := p))
    have hcd : ContDiffAt ℝ ∞ (η ∘ ψ) (chartFwd ι p p) := by
      have : g x₀ = chartFwd ι p p := rfl
      rw [← this]
      exact hη.contDiffAt.comp _ hψ
    exact (hcd.contMDiffAt.comp p hchart).congr_of_eventuallyEq hev
  · have hKc : IsClosed (σ '' K) :=
      (hK.image_of_continuousOn (hσ.cont.mono hKB)).isClosed
    have hev : ec_ext B σ η κ =ᶠ[𝓝 p] fun _ => (0 : ℝ) := by
      filter_upwards [hKc.isOpen_compl.mem_nhds hp] with q hq
      by_cases hqB : q ∈ σ '' B
      · rw [ec_ext_of_mem hqB]
        obtain ⟨y, hy, rfl⟩ := hqB
        rw [hκ y hy]
        exact hη0 y fun hyK => hq ⟨y, hyK, rfl⟩
      · exact ec_ext_of_notMem hqB
    exact contMDiffAt_const.congr_of_eventuallyEq hev

/-- A chart integral of the energy symbol is bounded by the global energy. -/
theorem integral_eSym_le_energyM [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M)
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (hF : IsLagImm ι J F)
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι (ball 0 6) σ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball 0 6))
    (hgraph : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, F (σ x) = graphX J U y₀ φ x)
    {R : ℝ} (hR : R < 6) (l N : ℕ) (hN : l + 3 ≤ N) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x) ≤ energyM ι A J l F := by
  have hnd := isNatDens_densOf_eDens (n := n) J l
  have hd0 : ∀ (X : EuclideanSpace ℝ (Fin n) → V) y, 0 ≤ eDens J X l y :=
    fun X y => ec_eDens_nonneg J X l y
  -- the cutoff
  set R₀ : ℝ := max R 1 with hR₀
  have hR₀6 : R₀ < 6 := max_lt hR (by norm_num)
  have hR₀pos : 0 < R₀ := lt_max_of_lt_right one_pos
  set R₁ : ℝ := (R₀ + 6) / 2 with hR₁
  have h1 : R₀ < R₁ := by rw [hR₁]; linarith
  have h2 : R₁ < 6 := by rw [hR₁]; linarith
  let bump : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨R₀, R₁, hR₀pos, h1⟩
  set η : EuclideanSpace ℝ (Fin n) → ℝ := fun x => bump x with hηdef
  have hη : ContDiff ℝ ∞ η := bump.contDiff
  have hη0 : ∀ x, x ∉ closedBall (0 : EuclideanSpace ℝ (Fin n)) R₁ → η x = 0 := by
    intro x hx
    have : x ∉ Function.support η := by
      rw [hηdef, bump.support_eq]
      exact fun hm => hx (ball_subset_closedBall hm)
    simpa using this
  have hη1 : ∀ x ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) R₀, η x = 1 :=
    fun x hx => bump.one_of_mem_closedBall hx
  have hηn : ∀ x, 0 ≤ η x := fun x => bump.nonneg
  have hηl : ∀ x, η x ≤ 1 := fun x => bump.le_one
  have hKB : closedBall (0 : EuclideanSpace ℝ (Fin n)) R₁ ⊆ ball 0 6 := closedBall_subset_ball h2
  obtain ⟨κ, hκc, hκ⟩ := hσ.inv
  set h : M → ℝ := ec_ext (ball 0 6) σ η κ with hh
  have hhs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h :=
    ec_contMDiff_ext ι hσ (isCompact_closedBall _ _) hKB hη hη0 hκc hκ
  have hσx : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, h (σ x) = η x := by
    intro x hx
    rw [hh, ec_ext_of_mem ⟨x, hx, rfl⟩, hκ x hx]
  have hh0 : ∀ p, 0 ≤ h p := by
    intro p
    by_cases hp : p ∈ σ '' ball 0 6
    · rw [hh, ec_ext_of_mem hp]; exact hηn _
    · rw [hh, ec_ext_of_notMem hp]
  have hh1 : ∀ p, h p ≤ 1 := by
    intro p
    by_cases hp : p ∈ σ '' ball 0 6
    · rw [hh, ec_ext_of_mem hp]; exact hηl _
    · rw [hh, ec_ext_of_notMem hp]; exact zero_le_one
  have hsupp : tsupport h ⊆ σ '' ball 0 6 := by
    have hKc : IsClosed (σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) R₁) :=
      ((isCompact_closedBall _ _).image_of_continuousOn (hσ.cont.mono hKB)).isClosed
    have hsub : Function.support h ⊆ σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) R₁ := by
      intro p hp
      by_cases hpB : p ∈ σ '' ball 0 6
      · obtain ⟨x, hx, rfl⟩ := hpB
        refine ⟨x, ?_, rfl⟩
        by_contra hxK
        apply hp
        rw [hσx x hx]; exact hη0 x hxK
      · exact absurd (by rw [hh, ec_ext_of_notMem hpB]) (Function.mem_support.1 hp)
    exact (closure_minimal hsub hKc).trans (Set.image_mono hKB)
  -- the integrals
  have hFp := hF.comp_param ι hσ
  have hcont : ContinuousOn (eDens J (F ∘ σ) l) (ball 0 6) :=
    (contDiffOn_eDens isOpen_ball hFp.1 hFp.2.1 l).continuousOn
  have hfc : ContinuousOn (fun x => η x * eDens J (F ∘ σ) l x) (ball 0 6) :=
    hη.continuous.continuousOn.mul hcont
  have hfi : IntegrableOn (fun x => η x * eDens J (F ∘ σ) l x) (ball 0 6) :=
    im_integrableOn_of_vanish isOpen_ball hfc (isCompact_closedBall _ _) hKB
      (fun y _ hy => by rw [hη0 y hy, zero_mul])
  have hint : intM ι A F h (densOf fun X => eDens J X l) =
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, η x * eDens J (F ∘ σ) l x := by
    rw [intM_eq_param ι A hF hhs hnd hσ hsupp]
    refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
    simp only [densOf, Function.comp_apply]
    rw [hσx x hx]
  have hmono := intM_densOf_mono ι A hF hnd hd0 hhs (h₂ := fun _ => 1) contMDiff_const hh1
  have hballR : ball (0 : EuclideanSpace ℝ (Fin n)) R ⊆ ball 0 6 := ball_subset_ball hR.le
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x)
      = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, η x * eDens J (F ∘ σ) l x := by
        refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
        have hxR : x ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) R₀ :=
          ball_subset_closedBall (ball_subset_ball (le_max_left R 1) hx)
        rw [hη1 x hxR, one_mul,
          eDens_eq_eSym_of_graph hJi hJJ hU isOpen_ball hφ hgraph hN (hballR hx)]
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 6, η x * eDens J (F ∘ σ) l x :=
        setIntegral_mono_set hfi (Filter.Eventually.of_forall fun x => mul_nonneg (hηn x) (hd0 _ _))
          hballR.eventuallyLE
    _ = intM ι A F h (densOf fun X => eDens J X l) := hint.symm
    _ ≤ intM ι A F (fun _ => 1) (densOf fun X => eDens J X l) := hmono
    _ = energyM ι A J l F := rfl

/-- info: 'CWFlow.Extension.integral_eSym_le_energyM' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_eSym_le_energyM

/-- info: 'CWFlow.Extension.intM_densOf_mono' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms intM_densOf_mono

/-- info: 'CWFlow.Extension.eDens_eq_eSym_of_graph' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms eDens_eq_eSym_of_graph

/-- info: 'CWFlow.Extension.energyM_eq_sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms energyM_eq_sum

end CWFlow.Extension
