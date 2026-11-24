import CWFlow.Extension.EnergySymbol
import CWFlow.Extension.JetOfSymbol
import CWFlow.Extension.SliceSmooth
import CWFlow.Extension.WeightIntegral

/-!
# Pointwise time derivative of the energy density

For a time-dependent graph function `Φ` solving `∂ₜφ = potSym` with admissible jets, the pointwise time
derivative of the energy density `eSym k` is, up to a remainder bounded by the weight function
`Om (2k+3) (k+2)`, the sum of two top order terms `2 q(D^{k+3}∂ₜφ + cTop(D^{k+2}∂ₜφ), ∇^k C) vol`.
-/

open CWFlow.Interp Metric
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: the setting -/

/-- The admissible set: the jets of order `≤ 3` are bounded by `K`. -/
def admSet (n N : ℕ) (K : ℝ) : Set (P n N) := {p | ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ‖p.2 m‖ < K}

theorem isOpen_admSet (K : ℝ) : IsOpen (admSet n N K) := by
  have : admSet n N K = ⋂ m : Fin (N + 1), {p : P n N | (m : ℕ) ≤ 3 → ‖p.2 m‖ < K} := by
    ext p; simp [admSet]
  rw [this]
  refine isOpen_iInter_of_finite fun m => ?_
  by_cases h : (m : ℕ) ≤ 3
  · simp only [h, true_imp_iff]
    exact isOpen_lt (continuous_norm.comp ((continuous_apply m).comp continuous_snd))
      continuous_const
  · simp [h]

theorem admSet_low (K : ℝ) :
    ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K' : ℝ, ∀ p ∈ admSet n N K, ‖p.2 m‖ ≤ K' :=
  fun m hm => ⟨K, fun p hp => (hp m hm).le⟩

/-- A time-dependent graph function solving the scalar equation `∂ₜφ = pot` on the ball of radius `R`,
with admissible jets there. -/
structure GraphEvol (N : ℕ) (K R : ℝ) (T : Set ℝ)
    (Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ) : Prop where
  isOpen : IsOpen T
  smooth : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))
  adm : ∀ t ∈ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R,
    jp N (fun y => Φ (t, y)) x ∈ admSet n N K
  eqn : ∀ t ∈ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R,
    dt Φ (t, x) = potSym (jp N (fun y => Φ (t, y)) x)

/-! ### Part 2: the jets of the time derivative -/

theorem GraphEvol.jet_dt_eq {N : ℕ} {K R : ℝ} {T : Set ℝ}
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (h : GraphEvol N K R T Φ) {t : ℝ} (ht : t ∈ T) {x : EuclideanSpace ℝ (Fin n)}
    (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R) (m : ℕ) (hN : m + 4 ≤ N) :
    jet m (fun y => dt Φ (t, y)) x = jetSym potSym m (jp N (fun y => Φ (t, y)) x) := by
  have hφ := contDiff_slice h.isOpen h.smooth ht
  have hev : (fun y => dt Φ (t, y)) =ᶠ[nhds x]
      (fun y => potSym (jp N (fun y => Φ (t, y)) y)) := by
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy using h.eqn t ht y hy
  rw [jet_congr_of_eventuallyEq hev m]
  exact jet_symbol_jp (isSymbol_potSym (isOpen_admSet (n := n) (N := N) K) (admSet_low K)
    (by omega)) (admSet_low K) m (by omega) hφ (h.adm t ht x hx)

/-! ### Part 3: the derivative of the density, split -/

theorem evp_gInner_comm {m : ℕ} (h : Z n 2) (hs : ∀ a b, h ![a, b] = h ![b, a]) (S T : Z n m) :
    gInner h S T = gInner h T S := by
  unfold gInner
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
  have : ∏ i, h ![y i, x i] = ∏ i, h ![x i, y i] := Finset.prod_congr rfl fun i _ => hs _ _
  rw [this]; ring

theorem evp_qSym_comm (m : ℕ) (p : P n N) (S T : Z n m) : qSym m p S T = qSym m p T S := by
  unfold qSym
  rw [gForm_apply, gForm_apply]
  exact evp_gInner_comm _ (ginvFun_symm _) S T

/-- The part of the derivative of `sSym k` in the `(k+2)`-jet slot (the top slot on which it depends). -/
noncomputable def cTop (k : ℕ) (p : P n N) (z : Z n (k + 2)) : Z n (k + 3) :=
  fderiv ℝ (sSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)),
    fun m : Fin (N + 1) => if h : (m : ℕ) = k + 2 then h ▸ z else 0)

theorem fderiv_eSym_split {K : ℝ} (k : ℕ) (hN : k + 3 ≤ N) {p : P n N} (hp : p ∈ admSet n N K)
    (ζ' : J n N) :
    fderiv ℝ (eSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ') =
      2 * qSym (k + 3) p (slotN (k + 3) ((0 : EuclideanSpace ℝ (Fin n)), ζ') +
            fderiv ℝ (sSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ')) (tSym k p)
          * volFun (low23 p).1
      + (fderiv ℝ (fun r : P n N => qSym (k + 3) r) p ((0 : EuclideanSpace ℝ (Fin n)), ζ'))
          (tSym k p) (tSym k p) * volFun (low23 p).1
      + qSym (k + 3) p (tSym k p) (tSym k p)
          * fderiv ℝ (fun r : P n N => volFun (low23 r).1) p
              ((0 : EuclideanSpace ℝ (Fin n)), ζ') := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hq := isSymbol_qSym hA hlow (n := n) (k + 3)
  have ht := isSymbol_tSym hA hlow k hN (n := n)
  have hv := isSymbol_vol hA (bdd_low23 hlow) (n := n)
  have hqd : DifferentiableAt ℝ (fun r : P n N => qSym (k + 3) r) p :=
    (hq.smooth.contDiffAt (hA.mem_nhds hp)).differentiableAt (by simp)
  have htd : DifferentiableAt ℝ (tSym (N := N) k) p :=
    (ht.smooth.contDiffAt (hA.mem_nhds hp)).differentiableAt (by simp)
  have hvd : DifferentiableAt ℝ (fun r : P n N => volFun (low23 r).1) p :=
    (hv.smooth.contDiffAt (hA.mem_nhds hp)).differentiableAt (by simp)
  have h1 := hqd.hasFDerivAt.clm_apply htd.hasFDerivAt
  have h2 := h1.clm_apply htd.hasFDerivAt
  have h3 := h2.fun_mul hvd.hasFDerivAt
  have hfe : (eSym (n := n) (N := N) k) =
      fun r => qSym (k + 3) r (tSym k r) (tSym k r) * volFun (low23 r).1 :=
    funext (eSym_eq (n := n) (N := N) k)
  have hs : fderiv ℝ (sSym (N := N) k) p = fderiv ℝ (tSym (N := N) k) p - slotN (k + 3) := by
    have e : (sSym (n := n) (N := N) k) = fun r => tSym k r - slotN (k + 3) r := rfl
    rw [e, fderiv_fun_sub htd (slotN (k + 3)).differentiableAt, ContinuousLinearMap.fderiv]
  rw [hfe, h3.fderiv, hs]
  have e : ∀ a b : Z n (k + 3), a + (b - a) = b := fun a b => by abel
  simp only [add_apply, sub_apply, smul_apply, ContinuousLinearMap.comp_apply,
    ContinuousLinearMap.flip_apply, smul_eq_mul]
  rw [e]
  have hsym := evp_qSym_comm (N := N) (k + 3) p (tSym k p)
    (fderiv ℝ (tSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ'))
  rw [hsym]
  ring

/-! ### Part 4: the remainder is of odd weight -/

theorem evp_Om_nonneg (v ℓ : ℕ) {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j) : 0 ≤ Om v ℓ a :=
  zero_le_one.trans (one_le_Om v ℓ ha)

/-- A symbol bound with a nonnegative constant. -/
theorem evp_norm_le' {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p ∈ A, ‖𝔭 p‖ ≤ C * Om w ℓ (aOf p) := by
  obtain ⟨C, hC⟩ := h.norm_le
  refine ⟨|C|, abs_nonneg _, fun p hp => (hC p hp).trans ?_⟩
  exact mul_le_mul_of_nonneg_right (le_abs_self C) (evp_Om_nonneg _ _ (aOf_nonneg p))

/-- Product of three weight functions. -/
theorem evp_Om_mul3 (k v1 v2 v3 l1 l2 l3 : ℕ) (hv : v1 + v2 + v3 ≤ 2 * k + 3)
    (hl : max (max l1 l2) l3 ≤ k + 2) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ a : ℕ → ℝ, (∀ j, 0 ≤ a j) →
      Om v1 l1 a * Om v2 l2 a * Om v3 l3 a ≤ C * Om (2 * k + 3) (k + 2) a := by
  obtain ⟨C1, hC1, h1⟩ := exists_Om_mul v1 v2
  obtain ⟨C2, hC2, h2⟩ := exists_Om_mul (v1 + v2) v3
  refine ⟨C1 * C2 * (((v1 + v2 + v3 : ℕ) : ℝ) + 1), by positivity, fun a ha => ?_⟩
  have hm := Om_mono (v := v1 + v2 + v3) (v' := 2 * k + 3) (ℓ := max (max l1 l2) l3) (ℓ' := k + 2)
    hv hl ha
  calc Om v1 l1 a * Om v2 l2 a * Om v3 l3 a
      ≤ (C1 * Om (v1 + v2) (max l1 l2) a) * Om v3 l3 a :=
        mul_le_mul_of_nonneg_right (h1 _ _ a ha) (evp_Om_nonneg _ _ ha)
    _ = C1 * (Om (v1 + v2) (max l1 l2) a * Om v3 l3 a) := by ring
    _ ≤ C1 * (C2 * Om (v1 + v2 + v3) (max (max l1 l2) l3) a) :=
        mul_le_mul_of_nonneg_left (h2 _ _ a ha) hC1
    _ ≤ C1 * (C2 * ((((v1 + v2 + v3 : ℕ) : ℝ) + 1) * Om (2 * k + 3) (k + 2) a)) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hm hC2) hC1
    _ = C1 * C2 * (((v1 + v2 + v3 : ℕ) : ℝ) + 1) * Om (2 * k + 3) (k + 2) a := by ring

/-- The weight bookkeeping for the slot `m ≤ k + 1` in the derivative of `sSym k`. -/
theorem evp_slot (k : ℕ) (hk : 1 ≤ k) (m : ℕ) (hm : m ≤ k + 1) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ a : ℕ → ℝ, (∀ j, 0 ≤ a j) →
      OmZ (((k - 1 : ℕ) : ℤ) - ((m - 3 : ℕ) : ℤ)) (k - 1) a * Om (1 + m) (1 + m) a * Om k k a
        ≤ C * Om (2 * k + 3) (k + 2) a := by
  obtain ⟨C, hC, h⟩ := evp_Om_mul3 k (k - 1 - (m - 3)) (1 + m) k (k - 1) (1 + m) k
    (by omega) (by omega)
  refine ⟨C, hC, fun a ha => ?_⟩
  have e : ((k - 1 : ℕ) : ℤ) - ((m - 3 : ℕ) : ℤ) = ((k - 1 - (m - 3) : ℕ) : ℤ) := by omega
  rw [e, OmZ_natCast]
  exact h a ha

theorem evp_Dlow (k : ℕ) (p : P n N) (ζ' : J n N) (h2 : k + 2 < N + 1) :
    fderiv ℝ (sSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ') - cTop k p (ζ' ⟨k + 2, h2⟩) =
      fderiv ℝ (sSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)),
        fun m : Fin (N + 1) => if (m : ℕ) = k + 2 then 0 else ζ' m) := by
  have cast : ∀ (m : Fin (N + 1)) (h : (m : ℕ) = k + 2),
      (h ▸ ζ' ⟨k + 2, h2⟩ : Z n m) = ζ' m := by
    rintro ⟨m, hm⟩ h
    simp only at h
    subst h
    rfl
  have e : ((0 : EuclideanSpace ℝ (Fin n)), ζ') =
      ((0 : EuclideanSpace ℝ (Fin n)), fun m : Fin (N + 1) => if (m : ℕ) = k + 2 then 0 else ζ' m)
      + ((0 : EuclideanSpace ℝ (Fin n)), fun m : Fin (N + 1) =>
          if h : (m : ℕ) = k + 2 then h ▸ ζ' ⟨k + 2, h2⟩ else 0) := by
    refine Prod.ext (by simp) ?_
    funext m
    by_cases h : (m : ℕ) = k + 2
    · simp only [Prod.snd_add, Pi.add_apply, if_pos h, dif_pos h, cast m h, zero_add]
    · simp [h]
  unfold cTop
  rw [e, map_add]
  abel

/-- The derivative of a smooth function of `ζ₂` only, in the jet directions, is controlled by the
`2`-jet component of the variation. -/
theorem evp_fderiv_comp_L2 {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] (g : Z n 2 → W)
    (hg : ContDiff ℝ ∞ g) {A : Set (P n N)}
    (hb : Bornology.IsBounded (low23 '' A : Set (Z n 2 × Z n 3))) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p ∈ A, ∀ δ : P n N,
      ‖fderiv ℝ (fun r : P n N => g (low23 r).1) p δ‖ ≤ C * ‖(low23 δ).1‖ := by
  have hc : Continuous (fun q : Z n 2 × Z n 3 => fderiv ℝ g q.1) :=
    (hg.continuous_fderiv (by simp)).comp continuous_fst
  obtain ⟨C0, hC0⟩ := (Metric.isCompact_of_isClosed_isBounded isClosed_closure
    hb.closure).exists_bound_of_continuousOn hc.continuousOn
  refine ⟨max C0 0, le_max_right _ _, fun p hp δ => ?_⟩
  have hL : HasFDerivAt (fun r : P n N => (low23 r).1)
      ((ContinuousLinearMap.fst ℝ (Z n 2) (Z n 3)).comp (low23 (n := n) (N := N))) p :=
    ((ContinuousLinearMap.fst ℝ (Z n 2) (Z n 3)).comp (low23 (n := n) (N := N))).hasFDerivAt
  have hd : HasFDerivAt (fun r : P n N => g (low23 r).1)
      ((fderiv ℝ g (low23 p).1).comp
        ((ContinuousLinearMap.fst ℝ (Z n 2) (Z n 3)).comp (low23 (n := n) (N := N)))) p :=
    ((hg.differentiable (by simp)) (low23 p).1).hasFDerivAt.comp p hL
  rw [hd.fderiv]
  have hmem : low23 p ∈ closure (low23 '' A : Set (Z n 2 × Z n 3)) :=
    subset_closure ⟨p, hp, rfl⟩
  have h1 := hC0 _ hmem
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.coe_fst'] at h1 ⊢
  refine (ContinuousLinearMap.le_opNorm _ _).trans ?_
  exact mul_le_mul_of_nonneg_right (h1.trans (le_max_left _ _)) (norm_nonneg _)

/-- The jet-space estimate behind Part 4. -/
theorem evp_core {K : ℝ} (k : ℕ) (hk : 1 ≤ k) (hN : k + 7 ≤ N) (B : ℝ) (hB : 0 ≤ B) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p ∈ admSet n N K, ∀ ζ' : J n N,
      (∀ m : Fin (N + 1), (m : ℕ) ≤ k + 1 → ‖ζ' m‖ ≤ B * Om (1 + (m : ℕ)) (1 + (m : ℕ)) (aOf p)) →
      |fderiv ℝ (eSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ')
          - 2 * qSym (k + 3) p (ζ' ⟨k + 3, by omega⟩ + cTop k p (ζ' ⟨k + 2, by omega⟩))
              (tSym k p) * volFun (low23 p).1|
        ≤ C * Om (2 * k + 3) (k + 2) (aOf p) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hb := bdd_low23 hlow
  obtain ⟨Cq, hCq0, hCq⟩ := evp_norm_le' (isSymbol_qSym hA hlow (n := n) (k + 3))
  obtain ⟨CT, hCT0, hCT⟩ := evp_norm_le' (isSymbol_tSym hA hlow k (by omega) (n := n))
  obtain ⟨Cv, hCv0, hCv⟩ := evp_norm_le' (isSymbol_vol hA hb (n := n))
  obtain ⟨Cs, hCs0, hCs⟩ := (isSymbol_sSym hA hlow k hk (by omega) (n := n)).fderiv_jets_le
  obtain ⟨Cdq, hCdq0, hCdq⟩ := evp_fderiv_comp_L2 (N := N)
    (fun u : Z n 2 => (gForm (ginvFun u) : Z n (k + 3) →L[ℝ] Z n (k + 3) →L[ℝ] ℝ))
    (contDiff_gForm.comp contDiff_ginvFun) hb
  obtain ⟨Cdv, hCdv0, hCdv⟩ := evp_fderiv_comp_L2 (N := N) volFun contDiff_volFun hb
  obtain ⟨C3, hC30, hC3⟩ := evp_Om_mul3 k 3 k k 3 k k (by omega) (by omega)
  have hslot : ∀ m : Fin (N + 1), ∃ C : ℝ, 0 ≤ C ∧ ((m : ℕ) ≤ k + 1 → ∀ a : ℕ → ℝ,
      (∀ j, 0 ≤ a j) →
      OmZ (((k - 1 : ℕ) : ℤ) - (((m : ℕ) - 3 : ℕ) : ℤ)) (k - 1) a * Om (1 + (m : ℕ)) (1 + (m : ℕ)) a
        * Om k k a ≤ C * Om (2 * k + 3) (k + 2) a) := by
    intro m
    by_cases hm : (m : ℕ) ≤ k + 1
    · obtain ⟨C, hC0, hC⟩ := evp_slot k hk m hm
      exact ⟨C, hC0, fun _ => hC⟩
    · exact ⟨0, le_rfl, fun h => absurd h hm⟩
  choose Cm hCm0 hCm using hslot
  refine ⟨2 * Cq * Cs * CT * Cv * (B * ∑ m, Cm m) + Cdq * B * CT * CT * Cv * C3
    + Cq * CT * CT * (Cdv * B) * C3, ?_, fun p hp ζ' hζ => ?_⟩
  · have : 0 ≤ ∑ m, Cm m := Finset.sum_nonneg fun m _ => hCm0 m
    positivity
  have ha : ∀ j, 0 ≤ aOf p j := aOf_nonneg p
  set a := aOf p with ha_def
  set Ω := Om (2 * k + 3) (k + 2) a with hΩ
  set wk := Om k k a with hwk
  have hΩ0 : 0 ≤ Ω := evp_Om_nonneg _ _ ha
  have hwk0 : 0 ≤ wk := evp_Om_nonneg _ _ ha
  have h2N : k + 2 < N + 1 := by omega
  have h3N : k + 3 < N + 1 := by omega
  have h2N' : 2 < N + 1 := by omega
  -- the split
  have hsplit := fderiv_eSym_split (K := K) k (by omega) hp ζ'
  have hs3 : slotN (k + 3) ((0 : EuclideanSpace ℝ (Fin n)), ζ') = ζ' ⟨k + 3, h3N⟩ :=
    slotN_apply (k + 3) h3N _
  have hD := evp_Dlow k p ζ' h2N
  set D := fderiv ℝ (sSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)),
    fun m : Fin (N + 1) => if (m : ℕ) = k + 2 then 0 else ζ' m) with hDdef
  set Dq := (fderiv ℝ (fun r : P n N => qSym (k + 3) r) p ((0 : EuclideanSpace ℝ (Fin n)), ζ'))
    with hDq
  set DV := fderiv ℝ (fun r : P n N => volFun (low23 r).1) p
    ((0 : EuclideanSpace ℝ (Fin n)), ζ') with hDV
  set V := volFun (low23 p).1 with hV
  set T := tSym k p with hT
  have hE : fderiv ℝ (eSym (N := N) k) p ((0 : EuclideanSpace ℝ (Fin n)), ζ')
          - 2 * qSym (k + 3) p (ζ' ⟨k + 3, by omega⟩ + cTop k p (ζ' ⟨k + 2, by omega⟩)) T * V
      = 2 * qSym (k + 3) p D T * V + Dq T T * V + qSym (k + 3) p T T * DV := by
    rw [hsplit, hs3, ← hD]
    simp only [map_add, map_sub, add_apply, sub_apply]
    ring
  rw [hE]
  -- bounds on the factors
  have hq : ‖qSym (k + 3) p‖ ≤ Cq := by
    have := hCq p hp
    rwa [Om_zero, mul_one] at this
  have hVb : |V| ≤ Cv := by
    have := hCv p hp
    rwa [Om_zero, mul_one, Real.norm_eq_abs] at this
  have hTb : ‖T‖ ≤ CT * wk := hCT p hp
  have hζ2 : ‖ζ' ⟨2, h2N'⟩‖ ≤ B * Om 3 3 a := hζ ⟨2, h2N'⟩ (by simp; omega)
  have hlow2 : ‖(low23 ((0 : EuclideanSpace ℝ (Fin n)), ζ')).1‖ = ‖ζ' ⟨2, h2N'⟩‖ := by
    rw [low23_fst _ h2N']
  have hDqb : ‖Dq‖ ≤ Cdq * (B * Om 3 3 a) := by
    have := hCdq p hp ((0 : EuclideanSpace ℝ (Fin n)), ζ')
    rw [hlow2] at this
    exact this.trans (mul_le_mul_of_nonneg_left hζ2 hCdq0)
  have hDVb : |DV| ≤ Cdv * (B * Om 3 3 a) := by
    have := hCdv p hp ((0 : EuclideanSpace ℝ (Fin n)), ζ')
    rw [hlow2] at this
    rw [← Real.norm_eq_abs]
    exact this.trans (mul_le_mul_of_nonneg_left hζ2 hCdv0)
  -- the slot sum
  set X := ∑ m : Fin (N + 1), OmZ (((k - 1 : ℕ) : ℤ) - (wt (some m : Slot N) : ℤ)) (k - 1) a *
    ‖(fun m : Fin (N + 1) => if (m : ℕ) = k + 2 then (0 : Z n m) else ζ' m) m‖ with hXdef
  have hDb : ‖D‖ ≤ Cs * X := hCs p hp _
  have hX : X * wk ≤ (B * ∑ m, Cm m) * Ω := by
    rw [hXdef, Finset.sum_mul, Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_le_sum fun m _ => ?_
    have hO0 : 0 ≤ OmZ (((k - 1 : ℕ) : ℤ) - (wt (some m : Slot N) : ℤ)) (k - 1) a :=
      OmZ_nonneg ha
    have hwt : (wt (some m : Slot N) : ℕ) = (m : ℕ) - 3 := rfl
    by_cases h1 : (m : ℕ) = k + 2
    · simp only [h1, if_true, norm_zero, mul_zero, zero_mul]
      exact mul_nonneg (mul_nonneg hB (hCm0 m)) hΩ0
    by_cases h2 : (m : ℕ) ≤ k + 1
    · simp only [h1, if_false]
      have hz := hζ m h2
      have hc := hCm m h2 a ha
      rw [hwt]
      calc OmZ (((k - 1 : ℕ) : ℤ) - (((m : ℕ) - 3 : ℕ) : ℤ)) (k - 1) a * ‖ζ' m‖ * wk
          ≤ OmZ (((k - 1 : ℕ) : ℤ) - (((m : ℕ) - 3 : ℕ) : ℤ)) (k - 1) a *
              (B * Om (1 + (m : ℕ)) (1 + (m : ℕ)) a) * wk :=
            mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hz hO0) hwk0
        _ = B * (OmZ (((k - 1 : ℕ) : ℤ) - (((m : ℕ) - 3 : ℕ) : ℤ)) (k - 1) a *
              Om (1 + (m : ℕ)) (1 + (m : ℕ)) a * wk) := by ring
        _ ≤ B * (Cm m * Ω) := mul_le_mul_of_nonneg_left hc hB
        _ = B * Cm m * Ω := by ring
    · have hneg : OmZ (((k - 1 : ℕ) : ℤ) - (wt (some m : Slot N) : ℤ)) (k - 1) a = 0 := by
        apply OmZ_of_neg
        rw [hwt]
        omega
      rw [hneg]
      simp only [zero_mul]
      exact mul_nonneg (mul_nonneg hB (hCm0 m)) hΩ0
  have hX0 : 0 ≤ X := Finset.sum_nonneg fun m _ =>
    mul_nonneg (OmZ_nonneg ha) (norm_nonneg _)
  have hO3 : 0 ≤ Om 3 3 a := evp_Om_nonneg _ _ ha
  -- the three terms
  have hqD : |qSym (k + 3) p D T| ≤ Cq * (Cs * X) * (CT * wk) := by
    have h := ContinuousLinearMap.le_opNorm₂ (qSym (k + 3) p) D T
    rw [← Real.norm_eq_abs]
    refine h.trans ?_
    refine mul_le_mul (mul_le_mul hq hDb (norm_nonneg _) hCq0) hTb (norm_nonneg _) ?_
    positivity
  have hDqT : |Dq T T| ≤ Cdq * (B * Om 3 3 a) * (CT * wk) * (CT * wk) := by
    have h' : ‖Dq T T‖ ≤ ‖Dq‖ * ‖T‖ * ‖T‖ := ContinuousLinearMap.le_opNorm₂ Dq T T
    rw [← Real.norm_eq_abs]
    refine h'.trans ?_
    refine mul_le_mul (mul_le_mul hDqb hTb (norm_nonneg _) (by positivity)) hTb (norm_nonneg _)
      (by positivity)
  have hqTT : |qSym (k + 3) p T T| ≤ Cq * (CT * wk) * (CT * wk) := by
    have h' : ‖qSym (k + 3) p T T‖ ≤ ‖qSym (k + 3) p‖ * ‖T‖ * ‖T‖ :=
      ContinuousLinearMap.le_opNorm₂ _ T T
    rw [← Real.norm_eq_abs]
    refine h'.trans ?_
    exact mul_le_mul (mul_le_mul hq hTb (norm_nonneg _) hCq0) hTb (norm_nonneg _) (by positivity)
  have hC3a := hC3 a ha
  have hSum : 0 ≤ ∑ m, Cm m := Finset.sum_nonneg fun m _ => hCm0 m
  -- term 1
  have t1 : |2 * qSym (k + 3) p D T * V| ≤
      2 * Cq * Cs * CT * Cv * (B * ∑ m, Cm m) * Ω := by
    calc |2 * qSym (k + 3) p D T * V| = 2 * |qSym (k + 3) p D T| * |V| := by
          rw [abs_mul, abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 2)]
      _ ≤ 2 * (Cq * (Cs * X) * (CT * wk)) * Cv :=
          mul_le_mul (mul_le_mul_of_nonneg_left hqD (by norm_num)) hVb (abs_nonneg _)
            (by positivity)
      _ = 2 * Cq * Cs * CT * Cv * (X * wk) := by ring
      _ ≤ 2 * Cq * Cs * CT * Cv * ((B * ∑ m, Cm m) * Ω) :=
          mul_le_mul_of_nonneg_left hX (by positivity)
      _ = _ := by ring
  have t2 : |Dq T T * V| ≤ Cdq * B * CT * CT * Cv * C3 * Ω := by
    calc |Dq T T * V| = |Dq T T| * |V| := abs_mul _ _
      _ ≤ (Cdq * (B * Om 3 3 a) * (CT * wk) * (CT * wk)) * Cv :=
          mul_le_mul hDqT hVb (abs_nonneg _) (by positivity)
      _ = Cdq * B * CT * CT * Cv * (Om 3 3 a * wk * wk) := by ring
      _ ≤ Cdq * B * CT * CT * Cv * (C3 * Ω) :=
          mul_le_mul_of_nonneg_left hC3a (by positivity)
      _ = _ := by ring
  have t3 : |qSym (k + 3) p T T * DV| ≤ Cq * CT * CT * (Cdv * B) * C3 * Ω := by
    calc |qSym (k + 3) p T T * DV| = |qSym (k + 3) p T T| * |DV| := abs_mul _ _
      _ ≤ (Cq * (CT * wk) * (CT * wk)) * (Cdv * (B * Om 3 3 a)) :=
          mul_le_mul hqTT hDVb (abs_nonneg _) (by positivity)
      _ = Cq * CT * CT * (Cdv * B) * (Om 3 3 a * wk * wk) := by ring
      _ ≤ Cq * CT * CT * (Cdv * B) * (C3 * Ω) :=
          mul_le_mul_of_nonneg_left hC3a (by positivity)
      _ = _ := by ring
  calc |2 * qSym (k + 3) p D T * V + Dq T T * V + qSym (k + 3) p T T * DV|
      ≤ |2 * qSym (k + 3) p D T * V| + |Dq T T * V| + |qSym (k + 3) p T T * DV| :=
        (abs_add_three _ _ _)
    _ ≤ _ := by linarith

/-- `Om v ℓ a` only depends on `a j` for `1 ≤ j ≤ min v ℓ`. -/
theorem evp_Om_congr {v ℓ : ℕ} {a b : ℕ → ℝ} (h : ∀ j, 1 ≤ j → j ≤ min v ℓ → a j = b j) :
    Om v ℓ a = Om v ℓ b := by
  unfold Om
  congr 1
  refine Finset.sum_congr rfl fun j hj => ?_
  rw [h j (Finset.mem_Icc.1 hj).1 (Finset.mem_Icc.1 hj).2]

/-- The jets of the time derivative are bounded by weight functions, uniformly on the admissible set. -/
theorem evp_jet_const (K : ℝ) (k : ℕ) (hN : k + 7 ≤ N) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ p ∈ admSet n N K, ∀ m : ℕ, m ≤ k + 1 →
      ‖jetSym potSym m p‖ ≤ B * Om (1 + m) (1 + m) (aOf p) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hjet : ∀ m : ℕ, ∃ C : ℝ, 0 ≤ C ∧ (m ≤ k + 1 → ∀ p ∈ admSet n N K,
      ‖jetSym potSym m p‖ ≤ C * Om (1 + m) (1 + m) (aOf p)) := by
    intro m
    by_cases hm : m ≤ k + 1
    · obtain ⟨C, hC0, hC⟩ := evp_norm_le'
        (isSymbol_jetSym (isSymbol_potSym hA hlow (by omega : 4 ≤ N)) hlow m (by omega))
      have e : max 1 1 + m = 1 + m := by omega
      rw [e] at hC
      exact ⟨C, hC0, fun _ => hC⟩
    · exact ⟨0, le_rfl, fun h => absurd h hm⟩
  choose Cj hCj0 hCj using hjet
  refine ⟨∑ m ∈ Finset.range (k + 2), Cj m, Finset.sum_nonneg fun m _ => hCj0 m,
    fun p hp m hm => ?_⟩
  refine (hCj m hm p hp).trans (mul_le_mul_of_nonneg_right ?_ (evp_Om_nonneg _ _ (aOf_nonneg p)))
  exact Finset.single_le_sum (f := Cj) (fun i _ => hCj0 i) (Finset.mem_range.2 (by omega))

theorem exists_const_dt_eSym_sub_le (N k : ℕ) (hk : 1 ≤ k) (hN : k + 7 ≤ N) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (R : ℝ) (T : Set ℝ) (Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
      GraphEvol N K R T Φ → ∀ t ∈ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R,
        HasDerivAt (fun s : ℝ => eSym k (jp N (fun y => Φ (s, y)) x))
          (fderiv ℝ (eSym (N := N) k) (jp N (fun y => Φ (t, y)) x)
            ((0 : EuclideanSpace ℝ (Fin n)),
              fun m : Fin (N + 1) => jet (m : ℕ) (fun y => dt Φ (t, y)) x)) t ∧
        |fderiv ℝ (eSym (N := N) k) (jp N (fun y => Φ (t, y)) x)
            ((0 : EuclideanSpace ℝ (Fin n)),
              fun m : Fin (N + 1) => jet (m : ℕ) (fun y => dt Φ (t, y)) x)
          - 2 * qSym (k + 3) (jp N (fun y => Φ (t, y)) x)
              (jet (k + 3) (fun y => dt Φ (t, y)) x
                + cTop k (jp N (fun y => Φ (t, y)) x) (jet (k + 2) (fun y => dt Φ (t, y)) x))
              (tSym k (jp N (fun y => Φ (t, y)) x))
            * volFun (low23 (jp N (fun y => Φ (t, y)) x)).1|
          ≤ C * Om (2 * k + 3) (k + 2) (aJet (fun y => Φ (t, y)) x) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  obtain ⟨B, hB0, hB⟩ := evp_jet_const (n := n) K k hN
  obtain ⟨C, hC0, hC⟩ := evp_core (n := n) (N := N) (K := K) k hk hN B hB0
  refine ⟨C, hC0, fun R T Φ h t ht x hx => ⟨?_, ?_⟩⟩
  · exact hasDerivAt_symbol_slice (isSymbol_eSym hA hlow k (by omega))
      (h.isOpen.prod isOpen_univ) h.smooth ⟨ht, Set.mem_univ _⟩ (h.adm t ht x hx)
  · have hp := h.adm t ht x hx
    have hζ : ∀ m : Fin (N + 1), (m : ℕ) ≤ k + 1 →
        ‖(fun m : Fin (N + 1) => jet (m : ℕ) (fun y => dt Φ (t, y)) x) m‖ ≤
          B * Om (1 + (m : ℕ)) (1 + (m : ℕ)) (aOf (jp N (fun y => Φ (t, y)) x)) := by
      intro m hm
      show ‖jet (m : ℕ) (fun y => dt Φ (t, y)) x‖ ≤ _
      rw [h.jet_dt_eq ht hx (m : ℕ) (by omega)]
      exact hB _ hp (m : ℕ) hm
    have key := hC _ hp _ hζ
    have e : Om (2 * k + 3) (k + 2) (aOf (jp N (fun y => Φ (t, y)) x)) =
        Om (2 * k + 3) (k + 2) (aJet (fun y => Φ (t, y)) x) := by
      refine evp_Om_congr fun j hj1 hj2 => ?_
      have hj : j ≤ k + 2 := hj2.trans ((min_le_right _ _))
      rw [aOf_jp, if_pos (by omega)]
      rfl
    rw [e] at key
    exact key

/-- info: 'CWFlow.Extension.GraphEvol.jet_dt_eq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms GraphEvol.jet_dt_eq

/-- info: 'CWFlow.Extension.fderiv_eSym_split' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_eSym_split

/-- info: 'CWFlow.Extension.exists_const_dt_eSym_sub_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_dt_eSym_sub_le

end CWFlow.Extension
