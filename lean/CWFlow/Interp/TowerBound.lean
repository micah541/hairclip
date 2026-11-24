import Mathlib.Analysis.Calculus.ContDiff.Bounds
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.Calculus.ContDiff.FTaylorSeries
import Mathlib.Analysis.Normed.Module.FiniteDimension

/-!
# Bounds for iterated derivatives along a tower

If `G` is smooth with bounded derivatives on an open set `S N` and `(Z N)` is a tower of maps whose
derivatives are expressed through higher levels of the tower, `D (Z N) = Λ N (·, Z (N + d))`, with
`Λ N` smooth with bounded derivatives, then every iterated derivative of `q ↦ G (q.1, Z N q)` is
bounded on `U` (`tower_bounded`).

Universe deviation: `E` and `W` live in the same universe `u` (needed to run the induction over the
target space `W`, since the derivative of `G` takes values in `ℝ × E →L[ℝ] W`).
-/

open Set
open scoped ContDiff Topology

namespace CWFlow.Interp

/-- `f` is smooth on the open set `S` and all its iterated derivatives are bounded on `S`. -/
def SmoothBdd {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [NormedAddCommGroup Y]
    [NormedSpace ℝ Y] (S : Set X) (f : X → Y) : Prop :=
  ContDiffOn ℝ ∞ f S ∧ ∀ m : ℕ, ∃ C : ℝ, ∀ x ∈ S, ‖iteratedFDeriv ℝ m f x‖ ≤ C

theorem tw_le_infty (m : ℕ) : ((m : ℕ) : WithTop ℕ∞) ≤ ∞ := by
  exact_mod_cast (le_top : (m : ℕ∞) ≤ ⊤)

theorem tw_smoothBdd_const {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] (S : Set X) (c : Y) :
    SmoothBdd S (fun _ : X => c) := by
  refine ⟨contDiffOn_const, fun m => ⟨‖c‖, fun x _ => ?_⟩⟩
  by_cases hm : m = 0
  · subst hm; simp [norm_iteratedFDeriv_zero]
  · simp [iteratedFDeriv_const_of_ne hm]

theorem tw_smoothBdd_add {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (hS : IsOpen S) {f g : X → Y}
    (hf : SmoothBdd S f) (hg : SmoothBdd S g) : SmoothBdd S (fun x => f x + g x) := by
  refine ⟨hf.1.add hg.1, fun m => ?_⟩
  obtain ⟨C1, h1⟩ := hf.2 m
  obtain ⟨C2, h2⟩ := hg.2 m
  refine ⟨C1 + C2, fun x hx => ?_⟩
  rw [← iteratedFDerivWithin_of_isOpen m hS hx]
  have := iteratedFDerivWithin_add_apply (𝕜 := ℝ) (i := m) ((hf.1 x hx).of_le (tw_le_infty m))
    ((hg.1 x hx).of_le (tw_le_infty m)) hS.uniqueDiffOn hx
  have e : (fun x => f x + g x) = f + g := rfl
  rw [e, this]
  refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
  · rw [iteratedFDerivWithin_of_isOpen m hS hx]; exact h1 x hx
  · rw [iteratedFDerivWithin_of_isOpen m hS hx]; exact h2 x hx

theorem tw_smoothBdd_clm_left {X Y Z : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] [NormedAddCommGroup Z] [NormedSpace ℝ Z]
    {S : Set X} (hS : IsOpen S) (L : Y →L[ℝ] Z) {f : X → Y} (hf : SmoothBdd S f) :
    SmoothBdd S (fun x => L (f x)) := by
  refine ⟨L.contDiff.comp_contDiffOn hf.1, fun m => ?_⟩
  obtain ⟨C, hC⟩ := hf.2 m
  refine ⟨‖L‖ * C, fun x hx => ?_⟩
  rw [← iteratedFDerivWithin_of_isOpen m hS hx]
  have := L.norm_iteratedFDerivWithin_comp_left (hf.1 x hx) hS.uniqueDiffOn hx (tw_le_infty m)
  refine this.trans (mul_le_mul_of_nonneg_left ?_ (norm_nonneg _))
  rw [iteratedFDerivWithin_of_isOpen m hS hx]; exact hC x hx

theorem tw_smoothBdd_comp_clm {X Y Z : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] [NormedAddCommGroup Z] [NormedSpace ℝ Z]
    {S : Set X} {S' : Set Y} (hS' : IsOpen S') (L : X →L[ℝ] Y)
    (hL : MapsTo L S S') {f : Y → Z} (hf : SmoothBdd S' f) :
    SmoothBdd S (fun x => f (L x)) := by
  have hopen : IsOpen (L ⁻¹' S') := hS'.preimage L.continuous
  have hsub : S ⊆ L ⁻¹' S' := fun x hx => hL hx
  refine ⟨(hf.1.comp_continuousLinearMap L).mono hsub, fun m => ?_⟩
  obtain ⟨C, hC⟩ := hf.2 m
  refine ⟨C * ‖L‖ ^ m, fun x hx => ?_⟩
  have h1 := L.iteratedFDerivWithin_comp_right hf.1 hS'.uniqueDiffOn hopen.uniqueDiffOn
    (hL hx) (tw_le_infty m)
  have h2 : iteratedFDerivWithin ℝ m (f ∘ L) (L ⁻¹' S') x = iteratedFDeriv ℝ m (fun x => f (L x)) x :=
    iteratedFDerivWithin_of_isOpen m hopen (hsub hx)
  have h3 : iteratedFDerivWithin ℝ m f S' (L x) = iteratedFDeriv ℝ m f (L x) :=
    iteratedFDerivWithin_of_isOpen m hS' (hL hx)
  rw [h1, h3] at h2
  rw [← h2]
  refine (ContinuousMultilinearMap.norm_compContinuousLinearMap_le _ _).trans ?_
  simp only [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  exact mul_le_mul_of_nonneg_right (hC _ (hL hx)) (by positivity)

theorem tw_smoothBdd_fderiv {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (hS : IsOpen S) {f : X → Y}
    (hf : SmoothBdd S f) : SmoothBdd S (fderiv ℝ f) := by
  refine ⟨hf.1.fderiv_of_isOpen hS (by simp), fun m => ?_⟩
  obtain ⟨C, hC⟩ := hf.2 (m + 1)
  exact ⟨C, fun x hx => by rw [norm_iteratedFDeriv_fderiv]; exact hC x hx⟩

theorem tw_smoothBdd_bilinear {X Y Z V : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] [NormedAddCommGroup Z] [NormedSpace ℝ Z]
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    {S : Set X} (hS : IsOpen S) (B : Y →L[ℝ] Z →L[ℝ] V) {f : X → Y} {g : X → Z}
    (hf : SmoothBdd S f) (hg : SmoothBdd S g) : SmoothBdd S (fun x => B (f x) (g x)) := by
  refine ⟨B.isBoundedBilinearMap.contDiff.comp_contDiffOn (hf.1.prodMk hg.1), fun m => ?_⟩
  choose Cf hCf using hf.2
  choose Cg hCg using hg.2
  refine ⟨‖B‖ * ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * Cf i * Cg (m - i), fun x hx => ?_⟩
  rw [← iteratedFDerivWithin_of_isOpen m hS hx]
  refine (B.norm_iteratedFDerivWithin_le_of_bilinear hf.1 hg.1 hS.uniqueDiffOn hx
    (tw_le_infty m)).trans (mul_le_mul_of_nonneg_left ?_ (norm_nonneg B))
  refine Finset.sum_le_sum fun i _ => ?_
  rw [iteratedFDerivWithin_of_isOpen i hS hx, iteratedFDerivWithin_of_isOpen (m - i) hS hx]
  have a := hCf i x hx
  have b := hCg (m - i) x hx
  have a0 : 0 ≤ ‖iteratedFDeriv ℝ i f x‖ := norm_nonneg _
  have b0 : 0 ≤ ‖iteratedFDeriv ℝ (m - i) g x‖ := norm_nonneg _
  have c0 : (0 : ℝ) ≤ (m.choose i : ℝ) := by positivity
  exact mul_le_mul (mul_le_mul_of_nonneg_left a c0) b b0 (mul_nonneg c0 (a0.trans a))


/-- The derivative of `G` along the tower is again of the same form, one level up. -/
theorem tw_smoothBdd_step {E J J' W : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [NormedAddCommGroup J] [NormedSpace ℝ J] [NormedAddCommGroup J'] [NormedSpace ℝ J']
    [NormedAddCommGroup W] [NormedSpace ℝ W] {S : Set (ℝ × J)} {S' : Set (ℝ × J')}
    (hS : IsOpen S) (hS' : IsOpen S') (π : J' →L[ℝ] J)
    (hSπ : ∀ p ∈ S', (p.1, π p.2) ∈ S) (Λ' : ℝ × J' → (ℝ × E →L[ℝ] J))
    (hΛ : SmoothBdd S' Λ') (G : ℝ × J → W) (hG : SmoothBdd S G) :
    SmoothBdd S' (fun p : ℝ × J' =>
      (fderiv ℝ G (p.1, π p.2)).comp ((ContinuousLinearMap.fst ℝ ℝ E).prod (Λ' p))) := by
  let L1 : ℝ × J' →L[ℝ] ℝ × J :=
    (ContinuousLinearMap.fst ℝ ℝ J').prod (π.comp (ContinuousLinearMap.snd ℝ ℝ J'))
  have hA : SmoothBdd S' (fun p => fderiv ℝ G (L1 p)) :=
    tw_smoothBdd_comp_clm hS L1 (fun p hp => hSπ p hp) (tw_smoothBdd_fderiv hS hG)
  have hB0 : SmoothBdd S' (fun p : ℝ × J' =>
      ((ContinuousLinearMap.inl ℝ ℝ J).comp (ContinuousLinearMap.fst ℝ ℝ E)) +
        (ContinuousLinearMap.compL ℝ (ℝ × E) J (ℝ × J) (ContinuousLinearMap.inr ℝ ℝ J)) (Λ' p)) :=
    tw_smoothBdd_add hS' (tw_smoothBdd_const S' _)
      (tw_smoothBdd_clm_left hS' _ hΛ)
  have hB : SmoothBdd S' (fun p : ℝ × J' => (ContinuousLinearMap.fst ℝ ℝ E).prod (Λ' p)) := by
    convert hB0 using 2 with p
    ext v <;> simp
  exact tw_smoothBdd_bilinear hS' (ContinuousLinearMap.compL ℝ (ℝ × E) (ℝ × J) W) hA hB

universe u v

theorem tw_tower_aux {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {𝕁 : ℕ → Type v} [∀ N, NormedAddCommGroup (𝕁 N)] [∀ N, NormedSpace ℝ (𝕁 N)]
    {U : Set (ℝ × E)} (hU : IsOpen U)
    (Z : ∀ N, ℝ × E → 𝕁 N) (d : ℕ)
    (π : ∀ N, 𝕁 (N + d) →L[ℝ] 𝕁 N) (hπ : ∀ N, ∀ q ∈ U, Z N q = π N (Z (N + d) q))
    (Λ : ∀ N, ℝ × 𝕁 (N + d) → (ℝ × E →L[ℝ] 𝕁 N))
    (S : ∀ N, Set (ℝ × 𝕁 N)) (hS : ∀ N, IsOpen (S N))
    (hmem : ∀ N, ∀ q ∈ U, (q.1, Z N q) ∈ S N)
    (hSπ : ∀ N, ∀ p ∈ S (N + d), (p.1, π N p.2) ∈ S N)
    (hΛ : ∀ N, SmoothBdd (S (N + d)) (Λ N))
    (hZ : ∀ N, ∀ q ∈ U, HasFDerivAt (Z N) (Λ N (q.1, Z (N + d) q)) q) :
    ∀ m : ℕ, ∀ (N : ℕ) (W : Type u) [NormedAddCommGroup W] [NormedSpace ℝ W]
      (G : ℝ × 𝕁 N → W), SmoothBdd (S N) G →
      ∃ C : ℝ, ∀ q ∈ U, ‖iteratedFDeriv ℝ m (fun q => G (q.1, Z N q)) q‖ ≤ C := by
  intro m
  induction m with
  | zero =>
    intro N W _ _ G hG
    obtain ⟨C, hC⟩ := hG.2 0
    refine ⟨C, fun q hq => ?_⟩
    have := hC _ (hmem N q hq)
    rw [norm_iteratedFDeriv_zero] at this ⊢
    exact this
  | succ m ih =>
    intro N W _ _ G hG
    set G' : ℝ × 𝕁 (N + d) → (ℝ × E →L[ℝ] W) := fun p =>
      (fderiv ℝ G (p.1, π N p.2)).comp ((ContinuousLinearMap.fst ℝ ℝ E).prod (Λ N p)) with hG'def
    have hG' : SmoothBdd (S (N + d)) G' :=
      tw_smoothBdd_step (hS N) (hS (N + d)) (π N) (hSπ N) (Λ N) (hΛ N) G hG
    obtain ⟨C, hC⟩ := ih (N + d) (ℝ × E →L[ℝ] W) G' hG'
    refine ⟨C, fun q hq => ?_⟩
    rw [← norm_iteratedFDeriv_fderiv]
    have hev : fderiv ℝ (fun q : ℝ × E => G (q.1, Z N q)) =ᶠ[𝓝 q]
        fun q => G' (q.1, Z (N + d) q) := by
      filter_upwards [hU.mem_nhds hq] with q' hq'
      have hGd : DifferentiableAt ℝ G (q'.1, Z N q') :=
        (hG.1.contDiffAt ((hS N).mem_nhds (hmem N q' hq'))).differentiableAt (by simp)
      have h1 : HasFDerivAt (fun x : ℝ × E => (x.1, Z N x))
          ((ContinuousLinearMap.fst ℝ ℝ E).prod (Λ N (q'.1, Z (N + d) q'))) q' :=
        (hasFDerivAt_fst (𝕜 := ℝ)).prodMk (hZ N q' hq')
      have h2 := hGd.hasFDerivAt.comp q' h1
      change fderiv ℝ (G ∘ fun x : ℝ × E => (x.1, Z N x)) q' = _
      rw [h2.fderiv]
      simp only [hG'def]
      rw [← hπ N q' hq']
    rw [(hev.iteratedFDeriv ℝ m).self_of_nhds]
    exact hC q hq

theorem tower_bounded {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
    {𝕁 : ℕ → Type v} [∀ N, NormedAddCommGroup (𝕁 N)] [∀ N, NormedSpace ℝ (𝕁 N)]
    [∀ N, FiniteDimensional ℝ (𝕁 N)]
    {U : Set (ℝ × E)} (hU : IsOpen U)
    (Z : ∀ N, ℝ × E → 𝕁 N) (d : ℕ)
    (π : ∀ N, 𝕁 (N + d) →L[ℝ] 𝕁 N) (hπ : ∀ N, ∀ q ∈ U, Z N q = π N (Z (N + d) q))
    (Λ : ∀ N, ℝ × 𝕁 (N + d) → (ℝ × E →L[ℝ] 𝕁 N))
    (S : ∀ N, Set (ℝ × 𝕁 N)) (hS : ∀ N, IsOpen (S N))
    (hmem : ∀ N, ∀ q ∈ U, (q.1, Z N q) ∈ S N)
    (hSπ : ∀ N, ∀ p ∈ S (N + d), (p.1, π N p.2) ∈ S N)
    (hΛ : ∀ N, SmoothBdd (S (N + d)) (Λ N))
    (hZ : ∀ N, ∀ q ∈ U, HasFDerivAt (Z N) (Λ N (q.1, Z (N + d) q)) q)
    {W : Type u} [NormedAddCommGroup W] [NormedSpace ℝ W] (N : ℕ) (G : ℝ × 𝕁 N → W)
    (hG : SmoothBdd (S N) G) :
    ∀ m : ℕ, ∃ C : ℝ, ∀ q ∈ U, ‖iteratedFDeriv ℝ m (fun q => G (q.1, Z N q)) q‖ ≤ C :=
  fun m => tw_tower_aux hU Z d π hπ Λ S hS hmem hSπ hΛ hZ m N W G hG

/-- info: 'CWFlow.Interp.tower_bounded' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms tower_bounded

end CWFlow.Interp
