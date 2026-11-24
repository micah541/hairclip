import CWFlow.Interp.TowerBound
import CWFlow.Extension.GaugeTerm
import CWFlow.Extension.SliceSmooth
import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct

/-!
# Bounded derivatives of the graph function in space and time

If `F (t, x)` solves `∂ₜF = potSym (jp 4 F_t x)` on `(a, b) × B_r` and all spatial jets are bounded,
then all joint `(t, x)`-derivatives of `F` and of `vSym (jp 5 F_t x) c` are bounded
(`graph_iteratedFDeriv_bounded`).  This is an instance of the abstract tower lemma `tower_bounded`
with `Z N q = jp N F_{q.1} q.2`.
-/

open Set Metric
open CWFlow.Interp
open scoped ContDiff Topology

namespace CWFlow.Extension

/-! ### Part 1: `SmoothBdd` for continuous linear maps, sums and symbols -/

theorem gtw_smoothBdd_clm {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (hS : Bornology.IsBounded S)
    (L : X →L[ℝ] Y) : SmoothBdd S (fun x => L x) := by
  refine ⟨L.contDiff.contDiffOn, fun m => ?_⟩
  obtain ⟨R, hR⟩ := isBounded_iff_forall_norm_le.1 hS
  rcases Nat.lt_or_ge m 2 with h | h
  · interval_cases m
    · refine ⟨‖L‖ * R, fun x hx => ?_⟩
      rw [norm_iteratedFDeriv_zero]
      exact (L.le_opNorm x).trans (mul_le_mul_of_nonneg_left (hR x hx) (norm_nonneg _))
    · refine ⟨‖L‖, fun x hx => ?_⟩
      refine ContinuousMultilinearMap.opNorm_le_bound (norm_nonneg _) fun w => ?_
      rw [iteratedFDeriv_one_apply, ContinuousLinearMap.fderiv]
      simpa using L.le_opNorm (w 0)
  · obtain ⟨k, rfl⟩ : ∃ k, m = k + 2 := ⟨m - 2, by omega⟩
    refine ⟨0, fun x hx => ?_⟩
    rw [iteratedFDeriv_clm_succ_succ]
    simp

theorem gtw_smoothBdd_sum {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (hS : IsOpen S) {ι : Type*}
    [DecidableEq ι] (s : Finset ι) {f : ι → X → Y} (hf : ∀ i ∈ s, SmoothBdd S (f i)) :
    SmoothBdd S (fun x => ∑ i ∈ s, f i x) := by
  induction s using Finset.induction_on with
  | empty =>
    simpa using tw_smoothBdd_const S (0 : Y)
  | insert i s hi ih =>
    have h1 := hf i (Finset.mem_insert_self i s)
    have h2 := ih fun j hj => hf j (Finset.mem_insert_of_mem hj)
    have := tw_smoothBdd_add hS h1 h2
    simpa [Finset.sum_insert hi] using this

theorem gtw_comp_le_norm {n N : ℕ} (s : Slot N) (δ : P n N) : comp s δ ≤ ‖δ‖ := by
  cases s with
  | none => exact norm_fst_le δ
  | some m => exact (norm_le_pi_norm δ.2 m).trans (norm_snd_le δ)

theorem gtw_aOf_le {n N : ℕ} (p : P n N) (j : ℕ) : aOf p j ≤ ‖p‖ := by
  unfold aOf
  split_ifs with h
  · exact (norm_le_pi_norm p.2 _).trans (norm_snd_le p)
  · exact norm_nonneg _

/-- A constant bounding `Om v ℓ a` when `a j ≤ R`. -/
noncomputable def gtwK (v ℓ : ℕ) (R : ℝ) : ℝ :=
  1 + ∑ j ∈ Finset.Icc 1 (min v ℓ), (max R 0) ^ ((v : ℝ) / j)

theorem gtwK_nonneg (v ℓ : ℕ) (R : ℝ) : 0 ≤ gtwK v ℓ R := by
  unfold gtwK
  have : 0 ≤ ∑ j ∈ Finset.Icc 1 (min v ℓ), (max R 0) ^ ((v : ℝ) / j) :=
    Finset.sum_nonneg fun j _ => Real.rpow_nonneg (le_max_right _ _) _
  linarith

theorem gtw_Om_le (v ℓ : ℕ) {a : ℕ → ℝ} {R : ℝ} (ha0 : ∀ j, 0 ≤ a j) (haR : ∀ j, a j ≤ R) :
    Om v ℓ a ≤ gtwK v ℓ R := by
  unfold Om gtwK
  refine add_le_add_right (Finset.sum_le_sum fun j _ => ?_) 1
  exact Real.rpow_le_rpow (ha0 j) ((haR j).trans (le_max_left _ _)) (by positivity)

theorem gtw_OmZ_le (v : ℤ) (ℓ : ℕ) {a : ℕ → ℝ} {R : ℝ} (ha0 : ∀ j, 0 ≤ a j)
    (haR : ∀ j, a j ≤ R) : OmZ v ℓ a ≤ gtwK v.toNat ℓ R := by
  unfold OmZ
  split_ifs
  · exact gtwK_nonneg _ _ _
  · exact gtw_Om_le _ _ ha0 haR

theorem gtw_smoothBdd_symbol {n N : ℕ} {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A B : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    (hBA : B ⊆ A) (hbd : Bornology.IsBounded B) : SmoothBdd B 𝔭 := by
  refine ⟨h.smooth.mono hBA, fun r => ?_⟩
  obtain ⟨R, hR⟩ := isBounded_iff_forall_norm_le.1 hbd
  obtain ⟨C, hC⟩ := h.bound r
  set Kt : ℝ := ∑ s : Fin r → Slot N,
    gtwK (((w : ℤ) - ∑ i, (wt (s i) : ℤ))).toNat ℓ R with hKt
  have hKt0 : 0 ≤ Kt := Finset.sum_nonneg fun s _ => gtwK_nonneg _ _ _
  refine ⟨|C| * Kt, fun p hp => ?_⟩
  refine ContinuousMultilinearMap.opNorm_le_bound (by positivity) fun δ => ?_
  refine (hC p (hBA hp) δ).trans ?_
  have hS := bdSum_nonneg (n := n) w ℓ r p δ
  have ha0 : ∀ j, 0 ≤ aOf p j := aOf_nonneg p
  have haR : ∀ j, aOf p j ≤ R := fun j => (gtw_aOf_le p j).trans (hR p hp)
  calc C * ∑ s : Fin r → Slot N,
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i)
      ≤ |C| * ∑ s : Fin r → Slot N,
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i) :=
        mul_le_mul_of_nonneg_right (le_abs_self C) hS
    _ ≤ |C| * (Kt * ∏ i, ‖δ i‖) := by
        refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg C)
        rw [hKt, Finset.sum_mul]
        refine Finset.sum_le_sum fun s _ => ?_
        refine mul_le_mul (gtw_OmZ_le _ _ ha0 haR)
          (Finset.prod_le_prod (fun i _ => comp_nonneg _ _) fun i _ => gtw_comp_le_norm _ _)
          (Finset.prod_nonneg fun i _ => comp_nonneg _ _) (gtwK_nonneg _ _ _)
    _ = |C| * Kt * ∏ i, ‖δ i‖ := by ring

/-- info: 'CWFlow.Extension.gtw_smoothBdd_symbol' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gtw_smoothBdd_symbol

/-! ### Part 2: localization by a cutoff -/

theorem gtw_loc {n : ℕ} {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    {x0 : EuclideanSpace ℝ (Fin n)} (hx0 : x0 ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r) :
    ∃ (G : ℝ × EuclideanSpace ℝ (Fin n) → ℝ) (ρ : ℝ), 0 < ρ ∧
      ball x0 ρ ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) r ∧
      ContDiffOn ℝ ∞ G (Ioo a b ×ˢ (univ : Set (EuclideanSpace ℝ (Fin n)))) ∧
      ∀ t ∈ Ioo a b, ∀ y ∈ ball x0 ρ, G (t, y) = F (t, y) := by
  obtain ⟨ε', hε', hsub'⟩ := Metric.isOpen_iff.1 (isOpen_ball (x := (0 : EuclideanSpace ℝ (Fin n)))
    (ε := r)) x0 hx0
  have hε : 0 < ε' / 2 := by positivity
  have hsub : closedBall x0 (ε' / 2) ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) r :=
    (closedBall_subset_ball (by linarith)).trans hsub'
  let f : ContDiffBump x0 := ⟨ε' / 4, ε' / 2, by positivity, by linarith⟩
  have hU : IsOpen (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) :=
    isOpen_Ioo.prod isOpen_ball
  refine ⟨fun q => f q.2 * F q, ε' / 4, by positivity, ?_, ?_, ?_⟩
  · exact (ball_subset_ball (by linarith)).trans hsub'
  · intro q hq
    obtain ⟨t, x⟩ := q
    obtain ⟨ht, -⟩ := hq
    by_cases hx : x ∈ closedBall x0 (ε' / 2)
    · have h1 : ContDiffAt ℝ ∞ F (t, x) := hF.contDiffAt (hU.mem_nhds ⟨ht, hsub hx⟩)
      have h2 : ContDiffAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => f q.2) (t, x) :=
        (f.contDiff.comp contDiff_snd).contDiffAt
      exact (h2.mul h1).contDiffWithinAt
    · have hev : (fun q : ℝ × EuclideanSpace ℝ (Fin n) => f q.2 * F q) =ᶠ[𝓝 (t, x)]
          fun _ => (0 : ℝ) := by
        filter_upwards [((isClosed_closedBall (x := x0) (ε := ε' / 2)).isOpen_compl.preimage
          continuous_snd).mem_nhds hx] with q hq
        have hq' : q.2 ∉ tsupport f := by
          rw [f.tsupport_eq]; exact hq
        rw [image_eq_zero_of_notMem_tsupport hq', zero_mul]
      exact (contDiffAt_const.congr_of_eventuallyEq hev).contDiffWithinAt
  · intro t ht y hy
    have : f y = 1 := f.one_of_mem_closedBall (ball_subset_closedBall hy)
    simp [this]

/-- Truncation of jets. -/
noncomputable def gtwProj (n N : ℕ) : P n (N + 4) →L[ℝ] P n N :=
  (ContinuousLinearMap.fst ℝ _ _).prod
    (ContinuousLinearMap.pi fun m : Fin (N + 1) =>
      slotProjL (n := n) (N := N + 4) (⟨(m : ℕ), by omega⟩ : Fin (N + 4 + 1)))

theorem gtwProj_jp {n : ℕ} (N : ℕ) (φ : EuclideanSpace ℝ (Fin n) → ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : gtwProj n N (jp (N + 4) φ x) = jp N φ x := rfl

theorem gtwProj_apply {n : ℕ} (N : ℕ) (p : P n (N + 4)) :
    gtwProj n N p = (p.1, fun m : Fin (N + 1) => p.2 ⟨(m : ℕ), by omega⟩) := rfl

/-- The jets agree for functions agreeing on an open set. -/
theorem gtw_jp_congr {n : ℕ} (N : ℕ) {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ}
    {s : Set (EuclideanSpace ℝ (Fin n))} (hs : IsOpen s) (h : ∀ y ∈ s, φ y = ψ y)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ s) : jp N φ x = jp N ψ x := by
  have hev : φ =ᶠ[𝓝 x] ψ := by
    filter_upwards [hs.mem_nhds hx] with y hy using h y hy
  refine Prod.ext rfl (funext fun m => ?_)
  exact jet_congr_of_eventuallyEq hev _

/-! ### Part 3: `potSym` along jets does not depend on the order of the jet space -/

theorem gtw_tSym_one_jp {n : ℕ} {N : ℕ} (hN : 4 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)}
    (hlow : ∀ m : ℕ, m ≤ 3 → ‖jet m φ x‖ < 2) (I : Fin 4 → Fin n) :
    tSym 1 (jp N φ x) I =
      fderiv ℝ (fun y => jet 3 φ y (Fin.tail I)) x (EuclideanSpace.single (I 0) 1)
      - ∑ i : Fin 3, ∑ c : Fin n,
          gamFun (jet 2 φ x, jet 3 φ x) ![I 0, Fin.tail I i, c] *
            jet 3 φ x (Function.update (Fin.tail I) i c) := by
  have hx : jp N φ x ∈ admSet n N 2 := fun m hm => hlow m hm
  have h := tSym_succ_jp_apply (isOpen_admSet (n := n) (N := N) 2) (admSet_low 2) 0 hN hφ hx I
  have e1 : ∀ y, tSym 0 (jp N φ y) = jet 3 φ y := fun y => tSym_zero_jp (by omega) φ y
  have e2 : gamSym (jp N φ x) = gamFun (jet 2 φ x, jet 3 φ x) := by
    unfold gamSym; rw [low23_jp (by omega)]
  simp only [e1, e2] at h
  exact h

theorem gtw_potSym_jp_eq {n : ℕ} {N : ℕ} (hN : 4 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)}
    (hlow : ∀ m : ℕ, m ≤ 3 → ‖jet m φ x‖ < 2) :
    potSym (jp N φ x) = potSym (jp 4 φ x) := by
  unfold potSym
  rw [low23_jp (by omega), low23_jp (by omega)]
  simp only [gtw_tSym_one_jp hN hφ hlow, gtw_tSym_one_jp le_rfl hφ hlow]

/-! ### Part 4: the tower operator -/

/-- `z ↦ ((τ, v) ↦ (0, Pi.single m (contract v z)))`. -/
noncomputable def gtwQ (n N : ℕ) (m : Fin (N + 1)) :
    Z n ((m : ℕ) + 1) →L[ℝ] (ℝ × EuclideanSpace ℝ (Fin n) →L[ℝ] P n N) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun z => (injSlot n N m).comp
        ((contractFlip n m z).comp (ContinuousLinearMap.snd ℝ ℝ (EuclideanSpace ℝ (Fin n))))
      map_add' := fun z z' => by
        ext w <;> simp
      map_smul' := fun c z => by
        ext w <;> simp }

/-- `z ↦ ((τ, v) ↦ τ • (0, Pi.single m z))`. -/
noncomputable def gtwT (n N : ℕ) (m : Fin (N + 1)) :
    Z n (m : ℕ) →L[ℝ] (ℝ × EuclideanSpace ℝ (Fin n) →L[ℝ] P n N) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun z => (ContinuousLinearMap.fst ℝ ℝ (EuclideanSpace ℝ (Fin n))).smulRight
        (injSlot n N m z)
      map_add' := fun z z' => by
        ext w <;> simp
      map_smul' := fun c z => by
        ext w <;> simp }

/-- The spatial part of the derivative of the jets: the shift. -/
noncomputable def gtwM (n N : ℕ) :
    P n (N + 4) →L[ℝ] (ℝ × EuclideanSpace ℝ (Fin n) →L[ℝ] P n N) :=
  ∑ m : Fin (N + 1), (gtwQ n N m).comp
    (slotProjL (n := n) (N := N + 4) (⟨(m : ℕ) + 1, by omega⟩ : Fin (N + 4 + 1)))

/-- The operator `Λ N` of the tower: the derivative of `q ↦ jp N F_{q.1} q.2` in terms of the
jets of order `N + 4`. -/
noncomputable def gtwΛ (n N : ℕ) (x : ℝ × P n (N + 4)) :
    ℝ × EuclideanSpace ℝ (Fin n) →L[ℝ] P n N :=
  (ContinuousLinearMap.inl ℝ (EuclideanSpace ℝ (Fin n)) (J n N)).comp
      (ContinuousLinearMap.snd ℝ ℝ (EuclideanSpace ℝ (Fin n))) + gtwM n N x.2 +
    ∑ m : Fin (N + 1), gtwT n N m (jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) x.2)

theorem gtw_sum_injSlot (n N : ℕ) (z : J n N) : ∑ m : Fin (N + 1), injSlot n N m (z m) = (0, z) := by
  refine Prod.ext (by simp only [injSlot_apply, Prod.fst_sum]; simp) ?_
  simp only [injSlot_apply, Prod.snd_sum]
  exact Finset.univ_sum_single z

theorem gtwΛ_apply (n N : ℕ) (t : ℝ) (p : P n (N + 4)) (τ : ℝ) (v : EuclideanSpace ℝ (Fin n)) :
    gtwΛ n N (t, p) (τ, v) =
      (v, fun m : Fin (N + 1) => contract v (p.2 ⟨(m : ℕ) + 1, by omega⟩ : Z n ((m : ℕ) + 1)) +
        τ • jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) p) := by
  have h1 : gtwM n N p (τ, v) = (0, fun m : Fin (N + 1) =>
      contract v (p.2 ⟨(m : ℕ) + 1, by omega⟩ : Z n ((m : ℕ) + 1))) := by
    have := gtw_sum_injSlot n N (fun m : Fin (N + 1) =>
      contract v (p.2 ⟨(m : ℕ) + 1, by omega⟩ : Z n ((m : ℕ) + 1)))
    rw [← this]
    simp only [gtwM, sum_apply, ContinuousLinearMap.comp_apply]
    exact Finset.sum_congr rfl fun m _ => rfl
  have h2 : (∑ m : Fin (N + 1), gtwT n N m (jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) p)) (τ, v) =
      τ • ((0 : EuclideanSpace ℝ (Fin n)), fun m : Fin (N + 1) =>
        jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) p) := by
    rw [← gtw_sum_injSlot n N (fun m : Fin (N + 1) =>
      jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) p), Finset.smul_sum]
    simp [gtwT]
  simp only [gtwΛ, add_apply, h1, h2]
  refine Prod.ext (by simp) ?_
  funext m
  simp

/-! ### Part 5: the derivative of the tower map -/

theorem gtw_dt_eq {n : ℕ} {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (heq : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      HasDerivAt (fun s => F (s, q.2)) (potSym (jp 4 (fun y => F (q.1, y)) q.2)) q.1)
    {t : ℝ} {y : EuclideanSpace ℝ (Fin n)}
    (hq : (t, y) ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) :
    dt F (t, y) = potSym (jp 4 (fun y' => F (t, y')) y) := by
  have hU : IsOpen (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) :=
    isOpen_Ioo.prod isOpen_ball
  have hd : DifferentiableAt ℝ F (t, y) :=
    (hF.contDiffAt (hU.mem_nhds hq)).differentiableAt (by simp)
  have hc : HasDerivAt (fun s : ℝ => (s, y)) ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t y)
  exact (hd.hasFDerivAt.comp_hasDerivAt t hc).unique (heq (t, y) hq)

theorem gtw_jet_dt {n : ℕ} {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (heq : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      HasDerivAt (fun s => F (s, q.2)) (potSym (jp 4 (fun y => F (q.1, y)) q.2)) q.1)
    (hlow : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2)
    {t : ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (hq : (t, x) ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) (N m : ℕ) (hm : m ≤ N) :
    jet m (fun y => dt F (t, y)) x =
      jetSym (potSym : P n (N + 4) → ℝ) m (jp (N + 4) (fun y => F (t, y)) x) := by
  obtain ⟨G, ρ, hρ, hsub, hG, hGF⟩ := gtw_loc hF hq.2
  have ht : t ∈ Ioo a b := hq.1
  set φ : EuclideanSpace ℝ (Fin n) → ℝ := fun y => G (t, y) with hφdef
  have hφ : ContDiff ℝ ∞ φ := contDiff_slice isOpen_Ioo hG ht
  have hx : x ∈ ball x ρ := mem_ball_self hρ
  have hjp : ∀ M : ℕ, ∀ y ∈ ball x ρ, jp M (fun y' => F (t, y')) y = jp M φ y := fun M y hy =>
    gtw_jp_congr M isOpen_ball (fun y' hy' => (hGF t ht y' hy').symm) hy
  have hlowφ : ∀ y ∈ ball x ρ, ∀ m ≤ 3, ‖jet m φ y‖ < 2 := by
    intro y hy m hm3
    have := hlow (t, y) ⟨ht, hsub hy⟩ m hm3
    have e := congrArg (fun p : P n 4 => p.2 ⟨m, by omega⟩) (hjp 4 y hy)
    simp only [jp_snd] at e
    simp only at this
    rwa [e] at this
  have hev : (fun y => dt F (t, y)) =ᶠ[𝓝 x]
      (fun y => potSym (jp (N + 4) φ y)) := by
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy
    rw [gtw_dt_eq hF heq ⟨ht, hsub hy⟩, hjp 4 y hy,
      (gtw_potSym_jp_eq (N := N + 4) (by omega) hφ (hlowφ y hy)).symm]
  rw [jet_congr_of_eventuallyEq hev m, hjp (N + 4) x hx]
  exact jet_symbol_jp (isSymbol_potSym (isOpen_admSet (n := n) (N := N + 4) 2) (admSet_low 2)
    (by omega)) (admSet_low 2) m (by omega) hφ (fun k hk => hlowφ x hx k hk)

theorem gtw_hasFDerivAt {n : ℕ} {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (heq : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      HasDerivAt (fun s => F (s, q.2)) (potSym (jp 4 (fun y => F (q.1, y)) q.2)) q.1)
    (hlow : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2)
    {t : ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (hq : (t, x) ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) (N : ℕ) :
    HasFDerivAt (fun q : ℝ × EuclideanSpace ℝ (Fin n) => jp N (fun y => F (q.1, y)) q.2)
      (gtwΛ n N (t, jp (N + 4) (fun y => F (t, y)) x)) (t, x) := by
  have hU : IsOpen (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) :=
    isOpen_Ioo.prod isOpen_ball
  obtain ⟨G, ρ, hρ, hsub, hG, hGF⟩ := gtw_loc hF hq.2
  have ht : t ∈ Ioo a b := hq.1
  set φ : EuclideanSpace ℝ (Fin n) → ℝ := fun y => G (t, y) with hφdef
  have hφ : ContDiff ℝ ∞ φ := contDiff_slice isOpen_Ioo hG ht
  have hx : x ∈ ball x ρ := mem_ball_self hρ
  have hjp : ∀ M : ℕ, ∀ t' ∈ Ioo a b, ∀ y ∈ ball x ρ,
      jp M (fun y' => F (t', y')) y = jp M (fun y' => G (t', y')) y := fun M t' ht' y hy =>
    gtw_jp_congr M isOpen_ball (fun y' hy' => (hGF t' ht' y' hy').symm) hy
  -- the smooth version of the tower map
  set Zt : ℝ × EuclideanSpace ℝ (Fin n) → P n N :=
    fun q => jp N (fun y => G (q.1, y)) q.2 with hZt
  have hZ : ContDiffOn ℝ ∞ Zt (Ioo a b ×ˢ (univ : Set (EuclideanSpace ℝ (Fin n)))) :=
    contDiffOn_jp_slice isOpen_Ioo hG
  have hat : ContDiffAt ℝ ∞ Zt (t, x) :=
    hZ.contDiffAt ((isOpen_Ioo.prod isOpen_univ).mem_nhds ⟨ht, mem_univ _⟩)
  have hD : HasFDerivAt Zt (fderiv ℝ Zt (t, x)) (t, x) :=
    (hat.differentiableAt (by simp)).hasFDerivAt
  have hevZ : (fun q : ℝ × EuclideanSpace ℝ (Fin n) => jp N (fun y => F (q.1, y)) q.2) =ᶠ[𝓝 (t, x)]
      Zt := by
    filter_upwards [(isOpen_Ioo.prod (isOpen_ball (x := x) (ε := ρ))).mem_nhds
      (show (t, x) ∈ Ioo a b ×ˢ ball x ρ from ⟨ht, hx⟩)] with q hq'
    exact hjp N q.1 hq'.1 q.2 hq'.2
  refine (hD.congr_of_eventuallyEq hevZ).congr_fderiv ?_
  -- identify the derivative
  set D := fderiv ℝ Zt (t, x) with hDdef
  set p : P n (N + 4) := jp (N + 4) (fun y => F (t, y)) x with hp
  have hpp : p = jp (N + 4) φ x := hjp (N + 4) t ht x hx
  -- time derivative
  have hcurve : HasDerivAt (fun s : ℝ => (s, x)) ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t x)
  have h1 := hD.comp_hasDerivAt t hcurve
  have h2 : HasDerivAt (fun s : ℝ => jp N (fun y => F (s, y)) x)
      ((0 : EuclideanSpace ℝ (Fin n)),
        fun m : Fin (N + 1) => jet (m : ℕ) (fun y => dt F (t, y)) x) t := by
    refine (hasDerivAt_const t x).prodMk (hasDerivAt_pi.2 fun m => ?_)
    exact hasDerivAt_jet_slice hU hF (m : ℕ) hq
  have h2' := h2.congr_of_eventuallyEq (f₁ := fun s : ℝ => Zt (s, x)) (by
    filter_upwards [isOpen_Ioo.mem_nhds ht] with s hs
    exact (hjp N s hs x hx).symm)
  have hDt : D ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) =
      ((0 : EuclideanSpace ℝ (Fin n)),
        fun m : Fin (N + 1) => jet (m : ℕ) (fun y => dt F (t, y)) x) := h1.unique h2'
  -- space derivative
  have hcurve2 : HasFDerivAt (fun y : EuclideanSpace ℝ (Fin n) => (t, y))
      (ContinuousLinearMap.inr ℝ ℝ (EuclideanSpace ℝ (Fin n))) x :=
    hasFDerivAt_prodMk_right t x
  have h3 := hD.comp x hcurve2
  have h4 := hasFDerivAt_jp (N := N) hφ x
  have hDx : D.comp (ContinuousLinearMap.inr ℝ ℝ (EuclideanSpace ℝ (Fin n))) =
      (ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin n))).prod
        (ContinuousLinearMap.pi fun m : Fin (N + 1) => fderiv ℝ (jet (m : ℕ) φ) x) :=
    h3.unique h4
  refine ContinuousLinearMap.ext fun w => ?_
  obtain ⟨τ, v⟩ := w
  have hsplit : ((τ, v) : ℝ × EuclideanSpace ℝ (Fin n)) =
      τ • ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) + ((0 : ℝ), v) := by
    ext <;> simp
  have hDv : D ((0 : ℝ), v) =
      (v, fun m : Fin (N + 1) => fderiv ℝ (jet (m : ℕ) φ) x v) := by
    simpa using congrArg (fun L => L v) hDx
  have hθ : ∀ m : Fin (N + 1), jet (m : ℕ) (fun y => dt F (t, y)) x =
      jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) p := fun m =>
    gtw_jet_dt hF heq hlow hq N m (by omega)
  have hc : ∀ m : Fin (N + 1), fderiv ℝ (jet (m : ℕ) φ) x v =
      contract v (p.2 ⟨(m : ℕ) + 1, by omega⟩ : Z n ((m : ℕ) + 1)) := fun m => by
    rw [fderiv_jet_eq_contract hφ, hpp]
    rfl
  conv_lhs => rw [hsplit]
  rw [map_add, map_smul, hDt, hDv, gtwΛ_apply]
  refine Prod.ext (by simp) ?_
  funext m
  simp [hθ, hc, add_comm]

/-! ### Part 6: the admissible bounded sets and `SmoothBdd` of the tower operator -/

/-- The bounded admissible set in jet space. -/
def gtwSp (n N : ℕ) (r : ℝ) (C : ℕ → ℝ) : Set (P n N) :=
  {p | p ∈ admSet n N 2 ∧ ‖p.1‖ < r ∧ ∀ m : Fin (N + 1), ‖p.2 m‖ < C m + 1}

theorem gtw_isOpen_Sp (n N : ℕ) (r : ℝ) (C : ℕ → ℝ) : IsOpen (gtwSp n N r C) := by
  have h3 : IsOpen {p : P n N | ∀ m : Fin (N + 1), ‖p.2 m‖ < C m + 1} := by
    have : {p : P n N | ∀ m : Fin (N + 1), ‖p.2 m‖ < C m + 1} =
        ⋂ m : Fin (N + 1), {p : P n N | ‖p.2 m‖ < C m + 1} := by
      ext p; simp
    rw [this]
    exact isOpen_iInter_of_finite fun m =>
      isOpen_lt (continuous_norm.comp ((continuous_apply m).comp continuous_snd)) continuous_const
  have h2 : IsOpen {p : P n N | ‖p.1‖ < r} :=
    isOpen_lt (continuous_norm.comp continuous_fst) continuous_const
  exact (isOpen_admSet (n := n) (N := N) 2).inter (h2.inter h3)

theorem gtw_bdd_Sp (n N : ℕ) (r : ℝ) (C : ℕ → ℝ) : Bornology.IsBounded (gtwSp n N r C) := by
  refine isBounded_iff_forall_norm_le.2 ⟨max r (∑ m : Fin (N + 1), (|C m| + 1)), fun p hp => ?_⟩
  obtain ⟨-, h1, h2⟩ := hp
  have hM : 0 ≤ ∑ m : Fin (N + 1), (|C m| + 1) :=
    Finset.sum_nonneg fun m _ => by positivity
  rw [Prod.norm_def]
  refine max_le_max h1.le ((pi_norm_le_iff_of_nonneg hM).2 fun m => ?_)
  have h3 : |C m| + 1 ≤ ∑ m : Fin (N + 1), (|C m| + 1) :=
    Finset.single_le_sum (f := fun m : Fin (N + 1) => |C m| + 1) (fun m _ => by positivity)
      (Finset.mem_univ m)
  have := h2 m
  have h4 : C m + 1 ≤ |C m| + 1 := by linarith [le_abs_self (C m)]
  linarith

theorem gtw_smoothBdd_Λ (n N : ℕ) (a b r : ℝ) (C : ℕ → ℝ) :
    SmoothBdd (Ioo a b ×ˢ gtwSp n (N + 4) r C) (gtwΛ n N) := by
  have hSp := gtw_isOpen_Sp n (N + 4) r C
  have hbd := gtw_bdd_Sp n (N + 4) r C
  have hS : IsOpen (Ioo a b ×ˢ gtwSp n (N + 4) r C) := isOpen_Ioo.prod hSp
  have hmaps : MapsTo (ContinuousLinearMap.snd ℝ ℝ (P n (N + 4)))
      (Ioo a b ×ˢ gtwSp n (N + 4) r C) (gtwSp n (N + 4) r C) := fun x hx => hx.2
  have f1 : SmoothBdd (Ioo a b ×ˢ gtwSp n (N + 4) r C)
      (fun _ => (ContinuousLinearMap.inl ℝ (EuclideanSpace ℝ (Fin n)) (J n N)).comp
        (ContinuousLinearMap.snd ℝ ℝ (EuclideanSpace ℝ (Fin n)))) := tw_smoothBdd_const _ _
  have f2 : SmoothBdd (Ioo a b ×ˢ gtwSp n (N + 4) r C) (fun x => gtwM n N x.2) :=
    tw_smoothBdd_comp_clm hSp (ContinuousLinearMap.snd ℝ ℝ (P n (N + 4))) hmaps
      (gtw_smoothBdd_clm hbd (gtwM n N))
  have f3 : SmoothBdd (Ioo a b ×ˢ gtwSp n (N + 4) r C) (fun x =>
      ∑ m : Fin (N + 1), gtwT n N m (jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) x.2)) := by
    refine gtw_smoothBdd_sum hS Finset.univ fun m _ => ?_
    have hsym := isSymbol_jetSym (isSymbol_potSym (isOpen_admSet (n := n) (N := N + 4) 2)
      (admSet_low 2) (by omega)) (admSet_low 2) (m : ℕ) (by have := m.2; omega)
    have hf : SmoothBdd (Ioo a b ×ˢ gtwSp n (N + 4) r C) (fun x =>
        jetSym (potSym : P n (N + 4) → ℝ) (m : ℕ) x.2) :=
      tw_smoothBdd_comp_clm hSp (ContinuousLinearMap.snd ℝ ℝ (P n (N + 4))) hmaps
        (gtw_smoothBdd_symbol hsym (fun p hp => hp.1) hbd)
    exact tw_smoothBdd_clm_left hS (gtwT n N m) hf
  exact tw_smoothBdd_add hS (tw_smoothBdd_add hS f1 f2) f3

/-! ### Part 7: the main theorem -/

theorem graph_iteratedFDeriv_bounded {n : ℕ} {a b r : ℝ}
    {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (heq : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      HasDerivAt (fun s => F (s, q.2)) (potSym (jp 4 (fun y => F (q.1, y)) q.2)) q.1)
    (hlow : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2)
    (hbd : ∀ m : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      ‖jet m (fun y => F (q.1, y)) q.2‖ ≤ C) :
    (∀ k : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      ‖iteratedFDeriv ℝ k F q‖ ≤ C) ∧
    (∀ (k : ℕ) (c : Fin n), ∃ C : ℝ, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      ‖iteratedFDeriv ℝ k (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        vSym (jp 5 (fun y => F (q.1, y)) q.2) c) q‖ ≤ C) := by
  have hU : IsOpen (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) :=
    isOpen_Ioo.prod isOpen_ball
  choose C hC using hbd
  let Zf : ∀ N : ℕ, ℝ × EuclideanSpace ℝ (Fin n) → P n N :=
    fun N q => jp N (fun y => F (q.1, y)) q.2
  let S : ∀ N : ℕ, Set (ℝ × P n N) := fun N => Ioo a b ×ˢ gtwSp n N r C
  have hS : ∀ N, IsOpen (S N) := fun N => isOpen_Ioo.prod (gtw_isOpen_Sp n N r C)
  have hSb : ∀ N, Bornology.IsBounded (S N) := fun N =>
    (Metric.isBounded_Ioo a b).prod (gtw_bdd_Sp n N r C)
  have hmem : ∀ N, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      (q.1, Zf N q) ∈ S N := by
    intro N q hq
    refine ⟨hq.1, fun m hm => hlow q hq (m : ℕ) hm, ?_, fun m => ?_⟩
    · have := hq.2
      rwa [mem_ball_zero_iff] at this
    · exact lt_of_le_of_lt (hC (m : ℕ) q hq) (lt_add_one _)
  have hSπ : ∀ N, ∀ p ∈ S (N + 4), (p.1, gtwProj n N p.2) ∈ S N := by
    intro N p hp
    obtain ⟨h1, h2, h3, h4⟩ := hp
    refine ⟨h1, fun m hm => h2 ⟨(m : ℕ), by omega⟩ hm, h3, fun m => h4 ⟨(m : ℕ), by omega⟩⟩
  have hΛ : ∀ N, SmoothBdd (S (N + 4)) (gtwΛ n N) := fun N => gtw_smoothBdd_Λ n N a b r C
  have hZ : ∀ N, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      HasFDerivAt (Zf N) (gtwΛ n N (q.1, Zf (N + 4) q)) q := fun N q hq =>
    gtw_hasFDerivAt hF heq hlow (t := q.1) (x := q.2) hq N
  have key : ∀ (N : ℕ) (G : ℝ × P n N → ℝ), SmoothBdd (S N) G → ∀ k : ℕ, ∃ C : ℝ,
      ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
        ‖iteratedFDeriv ℝ k (fun q => G (q.1, Zf N q)) q‖ ≤ C := fun N G hG =>
    tower_bounded (𝕁 := fun N => P n N) hU Zf 4 (fun N => gtwProj n N)
      (fun N q _ => (gtwProj_jp N _ _).symm) (gtwΛ n) S hS hmem hSπ hΛ hZ N G hG
  refine ⟨?_, ?_⟩
  · -- the function itself: the `0`-jet
    let L0 : ℝ × P n 0 →L[ℝ] ℝ :=
      ((PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 0 → Fin n => ℝ) (fun i => i.elim0)).comp
        (slotN (n := n) (N := 0) 0)).comp (ContinuousLinearMap.snd ℝ ℝ (P n 0))
    have hG : SmoothBdd (S 0) (fun p => L0 p) := gtw_smoothBdd_clm (hSb 0) L0
    have hE : (fun q : ℝ × EuclideanSpace ℝ (Fin n) => L0 (q.1, Zf 0 q)) = F := by
      funext q
      show (slotN (n := n) (N := 0) 0 (jp 0 (fun y => F (q.1, y)) q.2)) (fun i => i.elim0) = F q
      rw [slotN_jp 0 le_rfl, jet_apply, iteratedFDeriv_zero_apply]
    intro k
    simpa [hE] using key 0 _ hG k
  · intro k c
    have hsym := isSymbol_vSym (isOpen_admSet (n := n) (N := 5) 2) (admSet_low 2)
      (le_refl 5) c
    have hmaps : MapsTo (ContinuousLinearMap.snd ℝ ℝ (P n 5)) (S 5) (gtwSp n 5 r C) :=
      fun x hx => hx.2
    have hG : SmoothBdd (S 5) (fun p => vSym p.2 c) :=
      tw_smoothBdd_comp_clm (gtw_isOpen_Sp n 5 r C) (ContinuousLinearMap.snd ℝ ℝ (P n 5)) hmaps
        (gtw_smoothBdd_symbol hsym (fun p hp => hp.1) (gtw_bdd_Sp n 5 r C))
    exact key 5 _ hG k

/-- info: 'CWFlow.Extension.graph_iteratedFDeriv_bounded' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms graph_iteratedFDeriv_bounded

end CWFlow.Extension
