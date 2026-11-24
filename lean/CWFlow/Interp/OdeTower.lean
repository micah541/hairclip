import CWFlow.Interp.OdeJet
import CWFlow.Interp.TowerBound
import Mathlib.Analysis.ODE.Gronwall

/-!
# All joint derivatives of a solution of an ODE are bounded

For a solution `χ(t, y)` of `∂ₜχ = V(t, χ)` on `(t₁, b) × B` (smooth on the larger `(a, b) × B`), with `V`
smooth with bounded derivatives on `S` and the spatial jets of `χ(t₁, ·)` bounded, all joint derivatives
of `χ` are bounded on `(t₁, b) × B` (`ode_iteratedFDeriv_bounded`).

1. (`ot_smoothBdd_ojΘ`) the jet vector field `ojΘ V J` is smooth with bounded derivatives on the set
   `OjSet S ρ J` of jets whose `k`-th component has norm `< ρ k`;
2. (`ot_jet_bounds`) Grönwall bounds for the spatial jets, by induction on the order, using the top-linear
   structure `oj_top_clause`;
3. the tower lemma `tower_bounded` applied to `Z J = ojet n J χ`.
-/

set_option linter.unusedSectionVars false

open Set
open scoped ContDiff Topology

namespace CWFlow.Interp

/-! ### Generic lemmas on `SmoothBdd` -/

theorem ot_smoothBdd_clm_bdd {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (L : X →L[ℝ] Y) (M : ℝ)
    (hM : ∀ x ∈ S, ‖L x‖ ≤ M) : SmoothBdd S (fun x => L x) := by
  refine ⟨L.contDiff.contDiffOn, fun m => ?_⟩
  cases m with
  | zero =>
    refine ⟨M, fun x hx => ?_⟩
    rw [norm_iteratedFDeriv_zero]
    exact hM x hx
  | succ m =>
    obtain ⟨C, hC⟩ := (tw_smoothBdd_const S L).2 m
    refine ⟨C, fun x hx => ?_⟩
    rw [← norm_iteratedFDeriv_fderiv]
    have : fderiv ℝ (fun x => L x) = fun _ => L := by
      funext x; exact L.fderiv
    rw [this]
    exact hC x hx

theorem ot_smoothBdd_pi {X ι : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [Fintype ι]
    [DecidableEq ι] {φ : ι → Type*} [∀ i, NormedAddCommGroup (φ i)] [∀ i, NormedSpace ℝ (φ i)]
    {S : Set X} (hS : IsOpen S) {f : X → ∀ i, φ i}
    (hf : ∀ i, SmoothBdd S (fun x => f x i)) : SmoothBdd S f := by
  have key : ∀ s : Finset ι,
      SmoothBdd S (fun x => ∑ i ∈ s, ContinuousLinearMap.single ℝ φ i (f x i)) := by
    intro s
    induction s using Finset.induction_on with
    | empty => simpa using tw_smoothBdd_const S (0 : ∀ i, φ i)
    | insert j s hj ih =>
      simp only [Finset.sum_insert hj]
      exact tw_smoothBdd_add hS (tw_smoothBdd_clm_left hS _ (hf j)) ih
  have e : f = fun x => ∑ i, ContinuousLinearMap.single ℝ φ i (f x i) := by
    funext x
    ext j
    simp [Finset.sum_apply]
  rw [e]
  exact key _

theorem ot_smoothBdd_piLp {X ι W : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [Fintype ι]
    [NormedAddCommGroup W] [NormedSpace ℝ W] {S : Set X} (hS : IsOpen S)
    {f : X → PiLp 2 (fun _ : ι => W)} (hf : ∀ I, SmoothBdd S (fun x => f x I)) :
    SmoothBdd S f := by
  classical
  have h := tw_smoothBdd_clm_left hS
    (PiLp.continuousLinearEquiv 2 ℝ (fun _ : ι => W)).symm.toContinuousLinearMap
    (ot_smoothBdd_pi hS (f := fun x I => f x I) hf)
  exact h

/-! ### Linear maps on jets -/

section Jets

variable {n : ℕ}

local notation "E" => EuclideanSpace ℝ (Fin n)

/-- Truncation as a continuous linear map. -/
noncomputable def ojTruncLeL (n : ℕ) {i J : ℕ} (h : i ≤ J) : OJ n J →L[ℝ] OJ n i :=
  LinearMap.toContinuousLinearMap
    { toFun := ojTruncLe h
      map_add' := fun _ _ => rfl
      map_smul' := fun _ _ => rfl }

theorem ot_ojTruncLeL_apply {i J : ℕ} (h : i ≤ J) (z : OJ n J) :
    ojTruncLeL n h z = ojTruncLe h z := rfl

/-- The `a`-shift as a continuous linear map. -/
noncomputable def ojShiftL (n : ℕ) {J : ℕ} (a : Fin n) : OJ n (J + 1) →L[ℝ] OJ n J :=
  LinearMap.toContinuousLinearMap
    { toFun := ojShift a
      map_add' := fun x y => by funext k; ext K; rfl
      map_smul' := fun c x => by funext k; ext K; rfl }

theorem ot_ojShiftL_apply {J : ℕ} (a : Fin n) (z : OJ n (J + 1)) :
    ojShiftL n a z = ojShift a z := rfl

/-- Jets (with time) lying in the domain `S` whose `k`-th component has norm `< ρ k`. -/
def OjSet (S : Set (ℝ × E)) (ρ : ℕ → ℝ) (J : ℕ) : Set (ℝ × OJ n J) :=
  {p | (p.1, oj0 p.2) ∈ S ∧ ∀ k : Fin (J + 1), ‖p.2 k‖ < ρ (k : ℕ)}

theorem ot_isOpen_ojSet {S : Set (ℝ × E)} (hS : IsOpen S) (ρ : ℕ → ℝ) (J : ℕ) :
    IsOpen (OjSet S ρ J) := by
  have h1 : IsOpen {p : ℝ × OJ n J | (p.1, oj0 p.2) ∈ S} := oj_isOpen_dom hS J
  have h2 : IsOpen {p : ℝ × OJ n J | ∀ k : Fin (J + 1), ‖p.2 k‖ < ρ (k : ℕ)} := by
    have : {p : ℝ × OJ n J | ∀ k : Fin (J + 1), ‖p.2 k‖ < ρ (k : ℕ)} =
        ⋂ k : Fin (J + 1), {p : ℝ × OJ n J | ‖p.2 k‖ < ρ (k : ℕ)} := by
      ext p; simp
    rw [this]
    refine isOpen_iInter_of_finite fun k => ?_
    exact isOpen_lt (((continuous_apply k).comp continuous_snd).norm) continuous_const
  exact h1.inter h2

theorem ot_ojSet_norm_le (S : Set (ℝ × E)) (ρ : ℕ → ℝ) (J : ℕ) :
    ∃ M : ℝ, 0 ≤ M ∧ ∀ p ∈ OjSet S ρ J, ‖p.2‖ ≤ M := by
  refine ⟨∑ k ∈ Finset.range (J + 1), |ρ k|, Finset.sum_nonneg fun _ _ => abs_nonneg _, ?_⟩
  intro p hp
  refine (pi_norm_le_iff_of_nonneg (Finset.sum_nonneg fun _ _ => abs_nonneg _)).2 fun k => ?_
  have h1 : ‖p.2 k‖ ≤ |ρ (k : ℕ)| := (hp.2 k).le.trans (le_abs_self _)
  refine h1.trans ?_
  exact Finset.single_le_sum (f := fun k => |ρ k|) (fun _ _ => abs_nonneg _)
    (Finset.mem_range.2 k.2)

theorem ot_ojSet_trunc {S : Set (ℝ × E)} (ρ : ℕ → ℝ) {i J : ℕ} (h : i ≤ J) :
    ∀ p ∈ OjSet S ρ J, (p.1, ojTruncLeL n h p.2) ∈ OjSet S ρ i := by
  intro p hp
  refine ⟨hp.1, fun k => ?_⟩
  exact hp.2 ⟨k.1, Nat.lt_of_lt_of_le k.2 (Nat.succ_le_succ h)⟩

end Jets

/-! ### Step 1: the jet vector field is smooth with bounded derivatives -/

section Step1

variable {n : ℕ}

local notation "E" => EuclideanSpace ℝ (Fin n)

theorem ot_smoothBdd_ojθ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E} (hV : SmoothBdd S V)
    (ρ : ℕ → ℝ) : ∀ J : ℕ, SmoothBdd (OjSet S ρ J) (ojθ V J) := by
  intro J
  induction J with
  | zero =>
    classical
    refine ot_smoothBdd_piLp (ot_isOpen_ojSet hS ρ 0) fun I => ?_
    let L : ℝ × OJ n 0 →L[ℝ] ℝ × E := (ContinuousLinearMap.fst ℝ ℝ (OJ n 0)).prod
      ((oj0L n 0).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n 0)))
    exact tw_smoothBdd_comp_clm hS L (fun p hp => hp.1) hV
  | succ i ih =>
    classical
    refine ot_smoothBdd_piLp (ot_isOpen_ojSet hS ρ (i + 1)) fun A => ?_
    let L1 : ℝ × OJ n (i + 1) →L[ℝ] ℝ × OJ n i := (ContinuousLinearMap.fst ℝ ℝ (OJ n (i + 1))).prod
      ((ojTruncLeL n (Nat.le_succ i)).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n (i + 1))))
    have hD : SmoothBdd (OjSet S ρ (i + 1)) (fun p => fderiv ℝ (ojθ V i) (L1 p)) :=
      tw_smoothBdd_comp_clm (ot_isOpen_ojSet hS ρ i) L1 (fun p hp => ot_ojSet_trunc ρ (Nat.le_succ i) p hp)
        (tw_smoothBdd_fderiv (ot_isOpen_ojSet hS ρ i) ih)
    obtain ⟨M, hM0, hM⟩ := ot_ojSet_norm_le S ρ (i + 1)
    let g : ℝ × OJ n (i + 1) →L[ℝ] ℝ × OJ n i :=
      (0 : ℝ × OJ n (i + 1) →L[ℝ] ℝ).prod
        ((ojShiftL n (J := i) (A 0)).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n (i + 1))))
    have hg : SmoothBdd (OjSet S ρ (i + 1)) (fun p => g p) := by
      refine ot_smoothBdd_clm_bdd g (‖ojShiftL n (J := i) (A 0)‖ * M) fun p hp => ?_
      show ‖((0 : ℝ), ojShiftL n (J := i) (A 0) p.2)‖ ≤ _
      rw [Prod.norm_def]
      refine max_le ?_ ?_
      · simp only [norm_zero]
        exact mul_nonneg (norm_nonneg _) hM0
      · exact ((ojShiftL n (J := i) (A 0)).le_opNorm p.2).trans
          (mul_le_mul_of_nonneg_left (hM p hp) (norm_nonneg _))
    have hbil := tw_smoothBdd_bilinear (ot_isOpen_ojSet hS ρ (i + 1))
      (ContinuousLinearMap.id ℝ (ℝ × OJ n i →L[ℝ] PiLp 2 (fun _ : Fin i → Fin n => E))) hD hg
    exact tw_smoothBdd_clm_left (ot_isOpen_ojSet hS ρ (i + 1))
      (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin i → Fin n => E) (Fin.tail A)) hbil

theorem ot_smoothBdd_ojΘ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E} (hV : SmoothBdd S V)
    (ρ : ℕ → ℝ) (J : ℕ) : SmoothBdd (OjSet S ρ J) (ojΘ V J) := by
  classical
  refine ot_smoothBdd_pi (ot_isOpen_ojSet hS ρ J) fun k => ?_
  let L : ℝ × OJ n J →L[ℝ] ℝ × OJ n (k : ℕ) := (ContinuousLinearMap.fst ℝ ℝ (OJ n J)).prod
    ((ojTruncLeL n (Nat.lt_succ_iff.1 k.2)).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n J)))
  exact tw_smoothBdd_comp_clm (ot_isOpen_ojSet hS ρ k) L (fun p hp => ot_ojSet_trunc ρ (Nat.lt_succ_iff.1 k.2) p hp)
    (ot_smoothBdd_ojθ hS hV ρ k)

end Step1

/-! ### Step 2: Grönwall bounds for the spatial jets -/

theorem ot_gronwall {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] {a t₁ b : ℝ}
    (hat : a < t₁) (htb : t₁ < b) {f f' : ℝ → F} (hder : ∀ t ∈ Ioo a b, HasDerivAt f (f' t) t)
    {δ K ε : ℝ} (hδ : 0 ≤ δ) (hK : 0 ≤ K) (hε : 0 ≤ ε) (h0 : ‖f t₁‖ ≤ δ)
    (hbd : ∀ t ∈ Ico t₁ b, ‖f' t‖ ≤ K * ‖f t‖ + ε) :
    ∀ t ∈ Ico t₁ b, ‖f t‖ ≤ gronwallBound δ K ε (b - t₁) := by
  intro t ht
  have hsub : Icc t₁ t ⊆ Ioo a b := fun x hx => ⟨by linarith [hx.1], lt_of_le_of_lt hx.2 ht.2⟩
  have hc : ContinuousOn f (Icc t₁ t) := fun x hx =>
    (hder x (hsub hx)).continuousAt.continuousWithinAt
  have := norm_le_gronwallBound_of_norm_deriv_right_le (f := f) (f' := f') (a := t₁) (b := t) hc
    (fun x hx => (hder x (hsub ⟨hx.1, hx.2.le⟩)).hasDerivWithinAt) h0
    (fun x hx => hbd x ⟨hx.1, lt_of_lt_of_le hx.2 ht.2.le⟩) t ⟨ht.1, le_rfl⟩
  exact this.trans (gronwallBound_mono hδ hε hK (by linarith [ht.2]))

theorem ot_norm_comp_piLp {ι W : Type*} [Fintype ι] [NormedAddCommGroup W] [NormedSpace ℝ W]
    (L : W →L[ℝ] W) (u : PiLp 2 (fun _ : ι => W)) :
    ‖(WithLp.toLp 2 (fun I => L (u I)) : PiLp 2 (fun _ : ι => W))‖ ≤ ‖L‖ * ‖u‖ := by
  have h1 : ‖(WithLp.toLp 2 (fun I => L (u I)) : PiLp 2 (fun _ : ι => W))‖ ^ 2 ≤
      (‖L‖ * ‖u‖) ^ 2 := by
    rw [PiLp.norm_sq_eq_of_L2, mul_pow, PiLp.norm_sq_eq_of_L2, Finset.mul_sum]
    refine Finset.sum_le_sum fun I _ => ?_
    show ‖L (u I)‖ ^ 2 ≤ _
    rw [← mul_pow]
    exact pow_le_pow_left₀ (norm_nonneg _) (L.le_opNorm _) 2
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).1 h1

section Step2

variable {n : ℕ}

local notation "E" => EuclideanSpace ℝ (Fin n)

theorem ot_hasDerivAt_ojet {a b : ℝ} {B : Set E} (hB : IsOpen B) {S : Set (ℝ × E)}
    (hS : IsOpen S) {V : ℝ × E → E} (hV : ContDiffOn ℝ ∞ V S) {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1) (J : ℕ) :
    ∀ q ∈ Ioo a b ×ˢ B,
      HasDerivAt (fun s => ojet n J χ (s, q.2)) (ojΘ V J (q.1, ojet n J χ q)) q.1 := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  intro q hq
  obtain ⟨t, y⟩ := q
  obtain ⟨ht, hy⟩ := hq
  refine hasDerivAt_pi.2 fun k => ?_
  have h1 := CWFlow.Extension.hasDerivAt_jet_slice hΩ hχ (k : ℕ) (t := t) (x := y) ⟨ht, hy⟩
  rw [oj_jet_dt hS hV hB hχ hmaps (oj_dt_eq hB hχ hode) k t ht y hy] at h1
  exact h1

theorem ot_jet_step_zero {t₁ b : ℝ} (htb : t₁ < b) {B : Set E} {S : Set (ℝ × E)}
    {V : ℝ × E → E} (hV : SmoothBdd S V) {χ : ℝ × E → E} {a : ℝ} (hat : a < t₁)
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1)
    (hinit : ∀ j : ℕ, ∃ C : ℝ, ∀ y ∈ B, ‖jet j (fun y' => χ (t₁, y')) y‖ ≤ C) :
    ∃ r : ℝ, 0 ≤ r ∧ ∀ q ∈ Ico t₁ b ×ˢ B, ‖jet 0 (fun y => χ (q.1, y)) q.2‖ ≤ r := by
  obtain ⟨C0, hC0⟩ := hV.2 0
  obtain ⟨Ci, hCi⟩ := hinit 0
  refine ⟨max (gronwallBound (max Ci 0) 0 (max C0 0) (b - t₁)) 0, le_max_right _ _, ?_⟩
  rintro ⟨t, y⟩ ⟨ht, hy⟩
  refine le_trans ?_ (le_max_left _ _)
  rw [norm_jet_zero]
  have hmem : ∀ s ∈ Ico t₁ b, s ∈ Ioo a b := fun s hs => ⟨by linarith [hs.1], hs.2⟩
  refine ot_gronwall (f := fun s => χ (s, y)) (f' := fun s => V (s, χ (s, y))) hat htb
    (fun s hs => hode (s, y) ⟨hs, hy⟩) (le_max_right _ _) le_rfl (le_max_right _ _) ?_ ?_ t ht
  · have := hCi y hy
    rw [norm_jet_zero] at this
    exact this.trans (le_max_left _ _)
  · intro s hs
    have := hC0 _ (hmaps (s, y) ⟨hmem s hs, hy⟩)
    rw [norm_iteratedFDeriv_zero] at this
    simp only [zero_mul, zero_add]
    exact this.trans (le_max_left _ _)


theorem ot_norm_fderiv_slice_le {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : SmoothBdd S V) {K : ℝ} (hK : ∀ p ∈ S, ‖fderiv ℝ V p‖ ≤ K) {t : ℝ} {x : E}
    (hx : (t, x) ∈ S) : ‖fderiv ℝ (fun x' => V (t, x')) x‖ ≤ K := by
  have hVd : DifferentiableAt ℝ V (t, x) :=
    (hV.1.differentiableOn (by simp)).differentiableAt (hS.mem_nhds hx)
  have h2 : HasFDerivAt (fun x' : E => V (t, x'))
      ((fderiv ℝ V (t, x)).comp ((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E))) x :=
    hVd.hasFDerivAt.comp x ((hasFDerivAt_const t x).prodMk (hasFDerivAt_id x))
  rw [h2.fderiv]
  refine ContinuousLinearMap.opNorm_le_bound _ ((norm_nonneg _).trans (hK _ hx)) fun v => ?_
  have h3 : ‖((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E)) v‖ ≤ ‖v‖ := by
    simp [Prod.norm_def]
  calc ‖((fderiv ℝ V (t, x)).comp ((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E))) v‖
      ≤ ‖fderiv ℝ V (t, x)‖ * ‖((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E)) v‖ :=
        (fderiv ℝ V (t, x)).le_opNorm _
    _ ≤ K * ‖v‖ := mul_le_mul (hK _ hx) h3 (norm_nonneg _) ((norm_nonneg _).trans (hK _ hx))

theorem ot_jet_step_succ {t₁ b : ℝ} (htb : t₁ < b) {B : Set E} {S : Set (ℝ × E)}
    (hS : IsOpen S) (hB : IsOpen B)
    {V : ℝ × E → E} (hV : SmoothBdd S V) {χ : ℝ × E → E} {a : ℝ} (hat : a < t₁)
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1)
    (hinit : ∀ j : ℕ, ∃ C : ℝ, ∀ y ∈ B, ‖jet j (fun y' => χ (t₁, y')) y‖ ≤ C)
    (i : ℕ) (R : ℕ → ℝ)
    (hR : ∀ k < i + 1, 0 ≤ R k ∧ ∀ q ∈ Ico t₁ b ×ˢ B, ‖jet k (fun y => χ (q.1, y)) q.2‖ ≤ R k) :
    ∃ r : ℝ, 0 ≤ r ∧ ∀ q ∈ Ico t₁ b ×ˢ B, ‖jet (i + 1) (fun y => χ (q.1, y)) q.2‖ ≤ r := by
  set J := i + 1 with hJ
  set ρ : ℕ → ℝ := fun k => if k < J then R k + 1 else 1 with hρ
  obtain ⟨CΘ, hCΘ⟩ := (ot_smoothBdd_ojΘ hS hV ρ J).2 0
  obtain ⟨K1, hK1⟩ := (tw_smoothBdd_fderiv hS hV).2 0
  obtain ⟨Ci, hCi⟩ := hinit J
  have hK1' : ∀ p ∈ S, ‖fderiv ℝ V p‖ ≤ K1 := by
    intro p hp
    have := hK1 p hp
    rwa [norm_iteratedFDeriv_zero] at this
  set K := max K1 0 with hKdef
  set δ := max Ci 0 with hδdef
  set ε := max CΘ 0 with hεdef
  refine ⟨max (gronwallBound δ K ε (b - t₁)) 0, le_max_right _ _, ?_⟩
  rintro ⟨t, y⟩ ⟨ht, hy⟩
  refine le_trans ?_ (le_max_left _ _)
  have hmem : ∀ s ∈ Ico t₁ b, s ∈ Ioo a b := fun s hs => ⟨by linarith [hs.1], hs.2⟩
  have hder : ∀ s ∈ Ioo a b, HasDerivAt (fun s' => ojet n J χ (s', y) (Fin.last J))
      (ojΘ V J (s, ojet n J χ (s, y)) (Fin.last J)) s := fun s hs =>
    hasDerivAt_pi.1 (ot_hasDerivAt_ojet hB hS hV.1 hχ hmaps hode J (s, y) ⟨hs, hy⟩) (Fin.last J)
  refine ot_gronwall (f := fun s => ojet n J χ (s, y) (Fin.last J))
    (f' := fun s => ojΘ V J (s, ojet n J χ (s, y)) (Fin.last J)) hat htb hder
    (le_max_right _ _) (le_max_right _ _) (le_max_right _ _) ?_ ?_ t ht
  · exact (hCi y hy).trans (le_max_left _ _)
  · intro s hs
    have hsI : s ∈ Ioo a b := hmem s hs
    set z : OJ n J := ojet n J χ (s, y) with hz
    set z' : OJ n J := fun k => if (k : ℕ) < J then z k else 0 with hz'
    have hz0 : oj0 z = χ (s, y) := oj_oj0_ojet J χ (s, y)
    have hzS : (s, oj0 z) ∈ S := by rw [hz0]; exact hmaps (s, y) ⟨hsI, hy⟩
    have h0 : oj0 z = oj0 z' := by
      have : z' ⟨0, Nat.succ_pos J⟩ = z ⟨0, Nat.succ_pos J⟩ := if_pos (by simp [hJ])
      show z ⟨0, _⟩ _ = z' ⟨0, _⟩ _
      rw [this]
    have hlt : ∀ k : Fin (J + 1), (k : ℕ) < (Fin.last J : ℕ) → z k = z' k :=
      fun k hk => (if_pos hk).symm
    have hclause := fun I => oj_top_clause hS hV.1 J (Fin.last J) s z z' hzS h0 hlt I
    have hlast : z' (Fin.last J) = 0 := if_neg (by simp)
    set D := fderiv ℝ (fun x => V (s, x)) (oj0 z) with hD
    have hsplit : ojΘ V J (s, z) (Fin.last J) =
        ojΘ V J (s, z') (Fin.last J) +
          (WithLp.toLp 2 (fun I => D (z (Fin.last J) I)) :
            PiLp 2 (fun _ : Fin ((Fin.last J : Fin (J + 1)) : ℕ) → Fin n => E)) := by
      refine PiLp.ext fun I => ?_
      have := hclause I
      rw [hlast] at this
      simp only [PiLp.zero_apply, map_zero, sub_zero] at this
      rw [PiLp.add_apply]
      show _ = _ + D (z (Fin.last J) I)
      rw [← this]
      abel
    have hD1 : ‖D‖ ≤ K := by
      refine (ot_norm_fderiv_slice_le hS hV hK1' hzS).trans (le_max_left _ _)
    have hz'S : (s, z') ∈ OjSet S ρ J := by
      refine ⟨by show (s, oj0 z') ∈ S; rw [← h0]; exact hzS, fun k => ?_⟩
      by_cases hk : (k : ℕ) < J
      · have h1 : z' k = z k := if_pos hk
        show ‖z' k‖ < ρ (k : ℕ)
        rw [h1]
        have h2 : ρ (k : ℕ) = R k + 1 := if_pos hk
        rw [h2]
        have := (hR k hk).2 (s, y) ⟨hs, hy⟩
        have e : ‖z k‖ = ‖jet (k : ℕ) (fun y' => χ (s, y')) y‖ := rfl
        linarith
      · have h1 : z' k = 0 := if_neg hk
        have h2 : ρ (k : ℕ) = 1 := if_neg hk
        show ‖z' k‖ < ρ (k : ℕ)
        rw [h1, h2, norm_zero]
        exact one_pos
    have hΘ : ‖ojΘ V J (s, z') (Fin.last J)‖ ≤ ε := by
      have := hCΘ _ hz'S
      rw [norm_iteratedFDeriv_zero] at this
      exact ((norm_le_pi_norm (ojΘ V J (s, z')) (Fin.last J)).trans this).trans (le_max_left _ _)
    calc ‖ojΘ V J (s, z) (Fin.last J)‖
        = ‖ojΘ V J (s, z') (Fin.last J) +
          (WithLp.toLp 2 (fun I => D (z (Fin.last J) I)) :
            PiLp 2 (fun _ : Fin ((Fin.last J : Fin (J + 1)) : ℕ) → Fin n => E))‖ := by rw [hsplit]
      _ ≤ ‖ojΘ V J (s, z') (Fin.last J)‖ + ‖D‖ * ‖z (Fin.last J)‖ :=
        (norm_add_le _ _).trans (add_le_add le_rfl (ot_norm_comp_piLp D _))
      _ ≤ K * ‖z (Fin.last J)‖ + ε := by
        have := mul_le_mul_of_nonneg_right hD1 (norm_nonneg (z (Fin.last J)))
        linarith

theorem ot_jet_bounds {t₁ b : ℝ} (htb : t₁ < b) {B : Set E} {S : Set (ℝ × E)}
    (hS : IsOpen S) (hB : IsOpen B)
    {V : ℝ × E → E} (hV : SmoothBdd S V) {χ : ℝ × E → E} {a : ℝ} (hat : a < t₁)
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1)
    (hinit : ∀ j : ℕ, ∃ C : ℝ, ∀ y ∈ B, ‖jet j (fun y' => χ (t₁, y')) y‖ ≤ C) :
    ∃ R : ℕ → ℝ, ∀ j : ℕ, 0 ≤ R j ∧
      ∀ q ∈ Ico t₁ b ×ˢ B, ‖jet j (fun y => χ (q.1, y)) q.2‖ ≤ R j := by
  have key : ∀ j : ℕ, ∃ R : ℕ → ℝ, ∀ k < j, 0 ≤ R k ∧
      ∀ q ∈ Ico t₁ b ×ˢ B, ‖jet k (fun y => χ (q.1, y)) q.2‖ ≤ R k := by
    intro j
    induction j with
    | zero => exact ⟨fun _ => 0, fun k hk => absurd hk (Nat.not_lt_zero _)⟩
    | succ j ih =>
      obtain ⟨R, hR⟩ := ih
      have hnew : ∃ r : ℝ, 0 ≤ r ∧ ∀ q ∈ Ico t₁ b ×ˢ B,
          ‖jet j (fun y => χ (q.1, y)) q.2‖ ≤ r := by
        cases j with
        | zero => exact ot_jet_step_zero htb hV hat hmaps hode hinit
        | succ i => exact ot_jet_step_succ htb hS hB hV hat hχ hmaps hode hinit i R hR
      obtain ⟨r, hr0, hr⟩ := hnew
      refine ⟨fun k => if k = j then r else R k, fun k hk => ?_⟩
      by_cases h : k = j
      · subst h
        simp only [if_true]
        exact ⟨hr0, hr⟩
      · have := hR k (by omega)
        simpa [h] using this
  choose R hR using key
  exact ⟨fun j => R (j + 1) j, fun j => hR (j + 1) j (Nat.lt_succ_self j)⟩

end Step2

/-! ### Step 3: the tower of jets -/

section Step3

variable {n : ℕ}

local notation "E" => EuclideanSpace ℝ (Fin n)

/-- A continuous linear map out of `E` is determined by its values on the standard basis. -/
theorem ot_clm_ext_single {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F] (L M : E →L[ℝ] F)
    (h : ∀ i : Fin n, L (EuclideanSpace.single i 1) = M (EuclideanSpace.single i 1)) : L = M := by
  refine ContinuousLinearMap.coe_injective
    ((EuclideanSpace.basisFun (Fin n) ℝ).toBasis.ext fun i => ?_)
  simpa using h i

/-- The spatial derivative of a jet, as a linear function of the next jet. -/
noncomputable def ojSpL (n J : ℕ) :
    OJ n (J + 1) →L[ℝ] (EuclideanSpace ℝ (Fin n) →L[ℝ] OJ n J) :=
  ∑ a : Fin n, ((ContinuousLinearMap.smulRightL ℝ (EuclideanSpace ℝ (Fin n)) (OJ n J))
      (EuclideanSpace.proj a)).comp
    (ojShiftL n a)

theorem ot_ojSpL_apply {J : ℕ} (z : OJ n (J + 1)) (v : E) :
    ojSpL n J z v = ∑ a : Fin n, v a • ojShift a z := by
  simp [ojSpL, ot_ojShiftL_apply]

theorem ot_ojSpL_single {J : ℕ} (z : OJ n (J + 1)) (i : Fin n) :
    ojSpL n J z (EuclideanSpace.single i 1) = ojShift i z := by
  rw [ot_ojSpL_apply]
  simp

/-- The linear (spatial) part of the derivative of the jet tower, as a map on `ℝ × E`. -/
noncomputable def ojΛ₂ (n J : ℕ) :
    ℝ × OJ n (J + 1) →L[ℝ] (ℝ × EuclideanSpace ℝ (Fin n) →L[ℝ] OJ n J) :=
  (((ContinuousLinearMap.compL ℝ (ℝ × EuclideanSpace ℝ (Fin n)) (EuclideanSpace ℝ (Fin n))
      (OJ n J)).flip (ContinuousLinearMap.snd ℝ ℝ (EuclideanSpace ℝ (Fin n)))).comp
    ((ojSpL n J).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n (J + 1)))))

theorem ot_ojΛ₂_apply {J : ℕ} (p : ℝ × OJ n (J + 1)) (w : ℝ × E) :
    ojΛ₂ n J p w = ojSpL n J p.2 w.2 := rfl

theorem ot_norm_ojΛ₂_le {J : ℕ} (p : ℝ × OJ n (J + 1)) :
    ‖ojΛ₂ n J p‖ ≤ ‖ojSpL n J‖ * ‖p.2‖ := by
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun w => ?_
  rw [ot_ojΛ₂_apply]
  calc ‖ojSpL n J p.2 w.2‖ ≤ ‖ojSpL n J p.2‖ * ‖w.2‖ := (ojSpL n J p.2).le_opNorm _
    _ ≤ (‖ojSpL n J‖ * ‖p.2‖) * ‖w‖ :=
      mul_le_mul ((ojSpL n J).le_opNorm _) (norm_snd_le w) (norm_nonneg _) (by positivity)

/-- The derivative of the jet tower: `Λ J (t, z') (τ, v) = τ • Θ_J (t, z) + Σₐ vₐ • shiftₐ z'`. -/
noncomputable def ojΛ (V : ℝ × E → E) (J : ℕ) (p : ℝ × OJ n (J + 1)) : ℝ × E →L[ℝ] OJ n J :=
  ContinuousLinearMap.smulRightL ℝ (ℝ × E) (OJ n J) (ContinuousLinearMap.fst ℝ ℝ E)
      (ojΘ V J (p.1, ojTruncLeL n (Nat.le_succ J) p.2)) + ojΛ₂ n J p

theorem ot_ojΛ_apply (V : ℝ × E → E) {J : ℕ} (p : ℝ × OJ n (J + 1)) (τ : ℝ) (v : E) :
    ojΛ V J p (τ, v) = τ • ojΘ V J (p.1, ojTruncLeL n (Nat.le_succ J) p.2) + ojSpL n J p.2 v := by
  simp [ojΛ, ot_ojΛ₂_apply]

theorem ot_smoothBdd_ojΛ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E} (hV : SmoothBdd S V)
    (ρ : ℕ → ℝ) (J : ℕ) : SmoothBdd (OjSet S ρ (J + 1)) (ojΛ V J) := by
  have hopen := ot_isOpen_ojSet hS ρ (J + 1)
  let L : ℝ × OJ n (J + 1) →L[ℝ] ℝ × OJ n J := (ContinuousLinearMap.fst ℝ ℝ (OJ n (J + 1))).prod
    ((ojTruncLeL n (Nat.le_succ J)).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n (J + 1))))
  have h1 : SmoothBdd (OjSet S ρ (J + 1)) (fun p => ojΘ V J (L p)) :=
    tw_smoothBdd_comp_clm (ot_isOpen_ojSet hS ρ J) L
      (fun p hp => ot_ojSet_trunc ρ (Nat.le_succ J) p hp) (ot_smoothBdd_ojΘ hS hV ρ J)
  have h2 := tw_smoothBdd_clm_left hopen
    (ContinuousLinearMap.smulRightL ℝ (ℝ × E) (OJ n J) (ContinuousLinearMap.fst ℝ ℝ E)) h1
  obtain ⟨M, hM0, hM⟩ := ot_ojSet_norm_le S ρ (J + 1)
  have h3 : SmoothBdd (OjSet S ρ (J + 1)) (fun p => ojΛ₂ n J p) :=
    ot_smoothBdd_clm_bdd (ojΛ₂ n J) (‖ojSpL n J‖ * M) fun p hp =>
      (ot_norm_ojΛ₂_le p).trans (mul_le_mul_of_nonneg_left (hM p hp) (ContinuousLinearMap.opNorm_nonneg _))
  exact tw_smoothBdd_add hopen h2 h3

theorem ot_contDiffOn_ojet {a b : ℝ} {B : Set E} (hB : IsOpen B) {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B)) (J : ℕ) :
    ContDiffOn ℝ ∞ (ojet n J χ) (Ioo a b ×ˢ B) := by
  have hO : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  refine contDiffOn_pi.2 fun m => ?_
  refine contDiffOn_piLp' (𝕜 := ℝ) 2 (fun I => ?_)
  set e : Fin (m : ℕ) → ℝ × E := fun k => (0, EuclideanSpace.single (I k) 1) with he
  have hDm : ContDiffOn ℝ ∞ (iteratedFDeriv ℝ (m : ℕ) χ) (Ioo a b ×ˢ B) := by
    intro q hq
    exact ((hχ.contDiffAt (hO.mem_nhds hq)).iteratedFDeriv_right
      (infty_add_natCast_le (m : ℕ))).contDiffWithinAt
  have hev : ContDiffOn ℝ ∞ (fun q => iteratedFDeriv ℝ (m : ℕ) χ q e) (Ioo a b ×ˢ B) :=
    (ContinuousMultilinearMap.apply ℝ (fun _ : Fin (m : ℕ) => ℝ × E) E e).contDiff.comp_contDiffOn
      hDm
  refine hev.congr fun q hq => ?_
  show jet (m : ℕ) (fun y => χ (q.1, y)) q.2 I = _
  rw [jet_apply, CWFlow.Extension.iteratedFDeriv_slice hO hχ (m : ℕ) (show (q.1, q.2) ∈ _ from hq)]

theorem ot_hasFDerivAt_spatial {a b : ℝ} {B : Set E} (hB : IsOpen B) {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B)) (J : ℕ) {t : ℝ} (ht : t ∈ Ioo a b) {y : E}
    (hy : y ∈ B) :
    HasFDerivAt (fun y' => ojet n J χ (t, y')) (ojSpL n J (ojet n (J + 1) χ (t, y))) y := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  have hχt : ∀ y ∈ B, ContDiffAt ℝ ∞ (fun y => χ (t, y)) y := by
    intro y hy
    have h1 : ContDiffAt ℝ ∞ χ (t, y) := hχ.contDiffAt (hΩ.mem_nhds ⟨ht, hy⟩)
    exact ContDiffAt.comp y h1 (contDiffAt_const.prodMk contDiffAt_id)
  have hD : HasFDerivAt (fun y' => ojet n J χ (t, y'))
      (ContinuousLinearMap.pi fun k : Fin (J + 1) =>
        fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y) y := by
    refine hasFDerivAt_pi'.2 fun k => ?_
    have hdk : DifferentiableAt ℝ (jet (k : ℕ) (fun y => χ (t, y))) y := by
      have hd' : DifferentiableAt ℝ (iteratedFDeriv ℝ (k : ℕ) (fun y => χ (t, y))) y :=
        ((hχt y hy).iteratedFDeriv_right (infty_add_natCast_le k)).differentiableAt (by simp)
      exact (differentiableAt_piLp (𝕜 := ℝ) 2).2 fun I =>
        (hd'.hasFDerivAt.continuousMultilinear_apply_const
          (fun k => EuclideanSpace.single (I k) (1 : ℝ))).differentiableAt
    have := hdk.hasFDerivAt
    rw [ContinuousLinearMap.proj_pi]
    exact this
  refine hD.congr_fderiv (ot_clm_ext_single _ _ fun i => ?_)
  rw [ot_ojSpL_single]
  funext k
  refine PiLp.ext fun K => ?_
  show fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y (EuclideanSpace.single i 1) K = _
  rw [oj_fderiv_jet_apply (hχt y hy)]
  rfl

theorem ot_hasFDerivAt_ojet {a b : ℝ} {B : Set E} (hB : IsOpen B) {S : Set (ℝ × E)}
    (hS : IsOpen S) {V : ℝ × E → E} (hV : ContDiffOn ℝ ∞ V S) {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1) (J : ℕ) :
    ∀ q ∈ Ioo a b ×ˢ B, HasFDerivAt (ojet n J χ) (ojΛ V J (q.1, ojet n (J + 1) χ q)) q := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  rintro ⟨t, y⟩ hq
  have hd : DifferentiableAt ℝ (ojet n J χ) (t, y) :=
    ((ot_contDiffOn_ojet hB hχ J).contDiffAt (hΩ.mem_nhds hq)).differentiableAt (by simp)
  have HP := hd.hasFDerivAt
  have htime : fderiv ℝ (ojet n J χ) (t, y) (1, 0) = ojΘ V J (t, ojet n J χ (t, y)) := by
    have hs : HasDerivAt (fun s : ℝ => (s, y)) ((1 : ℝ), (0 : E)) t :=
      (hasDerivAt_id t).prodMk (hasDerivAt_const t y)
    have h1 := HP.comp_hasDerivAt t hs
    have h2 := ot_hasDerivAt_ojet hB hS hV hχ hmaps hode J (t, y) hq
    exact h1.unique h2
  have hspace : ∀ v, fderiv ℝ (ojet n J χ) (t, y) (0, v) = ojSpL n J (ojet n (J + 1) χ (t, y)) v := by
    have h2 : HasFDerivAt (fun y' => ojet n J χ (t, y'))
        ((fderiv ℝ (ojet n J χ) (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y :=
      HP.comp y (hasFDerivAt_prodMk_right t y)
    have h3 := ot_hasFDerivAt_spatial hB hχ J hq.1 hq.2
    intro v
    have := congrArg (fun L => L v) (h2.unique h3)
    simpa using this
  refine HP.congr_fderiv ?_
  refine ContinuousLinearMap.ext fun ⟨τ, v⟩ => ?_
  have hsplit : ((τ, v) : ℝ × E) = τ • ((1 : ℝ), (0 : E)) + ((0 : ℝ), v) := by
    ext <;> simp
  have e : fderiv ℝ (ojet n J χ) (t, y) (τ, v) =
      τ • ojΘ V J (t, ojet n J χ (t, y)) + ojSpL n J (ojet n (J + 1) χ (t, y)) v := by
    rw [hsplit, map_add, map_smul, htime, hspace]
  rw [e, ot_ojΛ_apply]
  rfl

end Step3

/-- For a solution `χ(t, y)` of `∂ₜχ = V(t, χ)` on `(a, b) × B`, with `V` smooth with bounded derivatives on
the open set `S` and the spatial jets of `χ(t₁, ·)` bounded on `B`, all joint derivatives of `χ` are bounded
on `(t₁, b) × B`. -/
theorem ode_iteratedFDeriv_bounded {n : ℕ} {t₁ b : ℝ} (htb : t₁ < b)
    {B : Set (EuclideanSpace ℝ (Fin n))} (hB : IsOpen B)
    {S : Set (ℝ × EuclideanSpace ℝ (Fin n))} (hS : IsOpen S)
    {V : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)} (hV : SmoothBdd S V)
    {χ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)} {a : ℝ} (hat : a < t₁)
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1)
    (hinit : ∀ j : ℕ, ∃ C : ℝ, ∀ y ∈ B, ‖jet j (fun y' => χ (t₁, y')) y‖ ≤ C) :
    ∀ k : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo t₁ b ×ˢ B, ‖iteratedFDeriv ℝ k χ q‖ ≤ C := by
  classical
  obtain ⟨R, hR⟩ := ot_jet_bounds htb hS hB hV hat hχ hmaps hode hinit
  set ρ : ℕ → ℝ := fun k => R k + 1 with hρ
  have hUo : IsOpen (Ioo t₁ b ×ˢ B) := isOpen_Ioo.prod hB
  have hUsub : ∀ q ∈ Ioo t₁ b ×ˢ B, q ∈ Ioo a b ×ˢ B :=
    fun q hq => ⟨⟨by linarith [hq.1.1], hq.1.2⟩, hq.2⟩
  have hUico : ∀ q ∈ Ioo t₁ b ×ˢ B, q ∈ Ico t₁ b ×ˢ B :=
    fun q hq => ⟨⟨hq.1.1.le, hq.1.2⟩, hq.2⟩
  have hmem : ∀ (N : ℕ), ∀ q ∈ Ioo t₁ b ×ˢ B, (q.1, ojet n N χ q) ∈ OjSet S ρ N := by
    intro N q hq
    refine ⟨?_, fun k => ?_⟩
    · show (q.1, oj0 (ojet n N χ q)) ∈ S
      rw [oj_oj0_ojet]
      exact hmaps q (hUsub q hq)
    · have := (hR k).2 q (hUico q hq)
      show ‖jet (k : ℕ) (fun y => χ (q.1, y)) q.2‖ < R k + 1
      linarith
  have hG : SmoothBdd (OjSet S ρ 0) (fun p : ℝ × OJ n 0 => oj0 p.2) := by
    obtain ⟨M, hM0, hM⟩ := ot_ojSet_norm_le S ρ 0
    refine ot_smoothBdd_clm_bdd ((oj0L n 0).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n 0)))
      (‖oj0L n 0‖ * M) fun p hp => ?_
    exact ((oj0L n 0).le_opNorm p.2).trans (mul_le_mul_of_nonneg_left (hM p hp) (norm_nonneg _))
  have key := tower_bounded (E := EuclideanSpace ℝ (Fin n)) (𝕁 := fun J => OJ n J)
    (U := Ioo t₁ b ×ˢ B) hUo (fun J => ojet n J χ) 1
    (fun J => ojTruncLeL n (Nat.le_succ J)) (fun N q _ => (oj_ojTrunc_ojet N χ q).symm)
    (fun J => ojΛ V J) (fun J => OjSet S ρ J) (fun J => ot_isOpen_ojSet hS ρ J) hmem
    (fun N p hp => ot_ojSet_trunc ρ (Nat.le_succ N) p hp)
    (fun J => ot_smoothBdd_ojΛ hS hV ρ J)
    (fun J q hq => ot_hasFDerivAt_ojet hB hS hV.1 hχ hmaps hode J q (hUsub q hq))
    (W := EuclideanSpace ℝ (Fin n)) 0 (fun p : ℝ × OJ n 0 => oj0 p.2) hG
  intro k
  obtain ⟨C, hC⟩ := key k
  refine ⟨C, fun q hq => ?_⟩
  have := hC q hq
  have e : (fun q : ℝ × EuclideanSpace ℝ (Fin n) => oj0 (ojet n 0 χ q)) = χ := by
    funext q; exact oj_oj0_ojet 0 χ q
  simpa only [e] using this

/-- info: 'CWFlow.Interp.ode_iteratedFDeriv_bounded' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms ode_iteratedFDeriv_bounded

end CWFlow.Interp
