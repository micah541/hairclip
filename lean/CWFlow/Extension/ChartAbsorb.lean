import CWFlow.Extension.ChartEnergy
import CWFlow.Extension.IBPBounds
import CWFlow.Interp.ThreeOrders

/-!
# Absorbing the top-order terms: chart inequalities (G2), (G3)

For a graph function `φ` with bounded low jets on `ball 0 3` and a cutoff `ρ` supported in `ball 0 2`:

* (G2) `∫ ρ² eSym (k+2) ≤ C ∫ ρ² |D^{k+5} φ|² + C (1 + energy3 (k+1) 3 φ)`;
* (G3) `∫ ρ² eSym (k+1) ≤ ε ∫_{B_3} |D^{k+5} φ|² + C (1 + energy3 k 3 φ)`.

Both follow from the pointwise bound `eSym (m+1) ≤ C |D^{m+4} φ|² + C Ω_m^{≤ m}(aJet φ)²`
together with the interpolation `exists_const_jet_three_orders`.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-- `volFun` is bounded on a ball. -/
theorem ca_vol_le (K : ℝ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ u : Z n 2, ‖u‖ ≤ K → volFun u ≤ C := by
  obtain ⟨C, hC⟩ := (isCompact_closedBall (0 : Z n 2) K).exists_bound_of_continuousOn
    (contDiff_volFun.continuous.continuousOn)
  refine ⟨max C 0, le_max_right _ _, fun u hu => ?_⟩
  have := hC u (by simpa using hu)
  exact (le_abs_self _).trans ((Real.norm_eq_abs _ ▸ this).trans (le_max_left _ _))

/-- Pointwise bound of the energy density by the top jet and the weight function. -/
theorem ca_pt (m : ℕ) (hN : m + 4 ≤ N) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ admSet n N K →
      eSym (m + 1) (jp N φ x) ≤ C * ‖jet (m + 1 + 3) φ x‖ ^ 2 + C * (Om m m (aJet φ x)) ^ 2 := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hN' : m + 1 + 3 ≤ N := by omega
  obtain ⟨Cs, hCs⟩ := norm_symbol_jp_le (isSymbol_sSym hA hlow (m + 1) (by omega) hN')
  obtain ⟨CK, hCK0, hCK⟩ := ca_vol_le (n := n) K
  refine ⟨2 * CK * (1 + Cs ^ 2), by positivity, fun φ x hx => ?_⟩
  set p := jp N φ x with hp
  have hu : ‖(low23 p).1‖ ≤ K := by
    have h2 : 2 < N + 1 := by omega
    rw [low23_fst p h2]
    exact (hx ⟨2, h2⟩ (by norm_num)).le
  have hvol : volFun (low23 p).1 ≤ CK := hCK _ hu
  have hjet : tSym (m + 1) p = jet (m + 1 + 3) φ x + sSym (m + 1) p := by
    have := slotN_jp (n := n) (m + 1 + 3) hN' φ x
    unfold sSym
    rw [← hp] at this
    rw [this]; abel
  have hs := hCs φ x hx
  simp only [Nat.add_sub_cancel] at hs
  have hOm : Om m m (fun j => if 3 + j ≤ N then ‖jet (3 + j) φ x‖ else 0) =
      Om m m (aJet φ x) := by
    refine evp_Om_congr fun j _ hj2 => ?_
    have : 3 + j ≤ N := by omega
    simp [aJet, this]
  rw [hOm] at hs
  have hs' : ‖sSym (m + 1) p‖ ≤ |Cs| * Om m m (aJet φ x) :=
    hs.trans (mul_le_mul_of_nonneg_right (le_abs_self _)
      (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x))))
  have hs2 : ‖sSym (m + 1) p‖ ^ 2 ≤ Cs ^ 2 * (Om m m (aJet φ x)) ^ 2 := by
    have := pow_le_pow_left₀ (norm_nonneg _) hs' 2
    rwa [mul_pow, sq_abs] at this
  have hsub : ‖tSym (m + 1) p‖ ≤ ‖jet (m + 1 + 3) φ x‖ + ‖sSym (m + 1) p‖ := by
    rw [hjet]; exact norm_add_le _ _
  have hsq : ‖tSym (m + 1) p‖ ^ 2 ≤
      2 * ‖jet (m + 1 + 3) φ x‖ ^ 2 + 2 * ‖sSym (m + 1) p‖ ^ 2 := by
    have h0 := norm_nonneg (tSym (m + 1) p)
    nlinarith [sq_nonneg (‖jet (m + 1 + 3) φ x‖ - ‖sSym (m + 1) p‖),
      norm_nonneg (jet (m + 1 + 3) φ x),
      norm_nonneg (sSym (m + 1) p), mul_self_le_mul_self h0 hsub]
  have hle := eSym_le (n := n) (N := N) (m + 1) p
  have h1 : eSym (m + 1) p ≤ ‖tSym (m + 1) p‖ ^ 2 * CK :=
    hle.trans (mul_le_mul_of_nonneg_left hvol (by positivity))
  have h3 : ‖tSym (m + 1) p‖ ^ 2 * CK ≤
      (2 * ‖jet (m + 1 + 3) φ x‖ ^ 2 + 2 * (Cs ^ 2 * (Om m m (aJet φ x)) ^ 2)) * CK :=
    mul_le_mul_of_nonneg_right (by linarith) hCK0
  have e1 := mul_nonneg (mul_nonneg hCK0 (sq_nonneg Cs)) (sq_nonneg ‖jet (m + 1 + 3) φ x‖)
  have e2 := mul_nonneg hCK0 (sq_nonneg (Om m m (aJet φ x)))
  nlinarith [h1.trans h3]

theorem ca_rho_zero {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2) : ρ x = 0 :=
  image_eq_zero_of_notMem_tsupport fun h => hx (hρ.supp h)

theorem ca_rho_sq_le {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    (x : EuclideanSpace ℝ (Fin n)) : ρ x ^ 2 ≤ 1 := by
  have h := hρ.le_one x
  have h0 := abs_nonneg (ρ x)
  nlinarith [sq_abs (ρ x)]

/-- Integral over the whole space of a function supported in `ball 0 2`. -/
theorem ca_integral_eq {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    (f : EuclideanSpace ℝ (Fin n) → ℝ) :
    ∫ x, ρ x ^ 2 * f x = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ρ x ^ 2 * f x :=
  (setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => by
    rw [ca_rho_zero hρ hx]; ring).symm

theorem ca_main {m : ℕ} (hN : m + 4 ≤ N) {K : ℝ} :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (Λ : ℝ) (ρ : EuclideanSpace ℝ (Fin n) → ℝ), IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∫ x, ρ x ^ 2 * eSym (m + 1) (jp N φ x) ≤
        C * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ρ x ^ 2 * ‖jet (m + 4) φ x‖ ^ 2) +
        C * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, (Om m m (aJet φ x)) ^ 2 := by
  obtain ⟨C, hC0, hC⟩ := ca_pt (n := n) (N := N) m hN K
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  refine ⟨C, hC0, fun Λ ρ hρ φ hφ hadm => ?_⟩
  rw [ca_integral_eq hρ]
  have hc : ContinuousOn (fun x => eSym (m + 1) (jp N φ x))
      (ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    (contDiffOn_symbol_jp (isSymbol_eSym hA hlow (m + 1) (by omega)) hφ).continuousOn.mono hadm
  have hF : ContinuousOn (fun x => ρ x ^ 2 * eSym (m + 1) (jp N φ x))
      (ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    ((hρ.smooth.continuous.pow 2).continuousOn).mul hc
  have hFint : IntegrableOn (fun x => ρ x ^ 2 * eSym (m + 1) (jp N φ x))
      (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    ((hF.mono (closedBall_subset_ball (by norm_num))).integrableOn_compact
      (isCompact_closedBall 0 2)).mono_set ball_subset_closedBall
  have hJ : IntegrableOn (fun x => ρ x ^ 2 * ‖jet (m + 4) φ x‖ ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    integrableOn_ball_of_continuous ((hρ.smooth.continuous.pow 2).mul
      ((contDiff_jet hφ (m + 4)).continuous.norm.pow 2)) 2
  have hO : IntegrableOn (fun x => (Om m m (aJet φ x)) ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    integrableOn_ball_of_continuous ((continuous_Om hφ m m).pow 2) 2
  have hsub : ball (0 : EuclideanSpace ℝ (Fin n)) 2 ⊆ ball 0 3 := ball_subset_ball (by norm_num)
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ρ x ^ 2 * eSym (m + 1) (jp N φ x)
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2,
          (C * (ρ x ^ 2 * ‖jet (m + 4) φ x‖ ^ 2) + C * (Om m m (aJet φ x)) ^ 2) := by
        refine setIntegral_mono_on hFint ((hJ.const_mul C).add (hO.const_mul C))
          measurableSet_ball fun x hx => ?_
        have h1 := hC φ x (hadm x (hsub hx))
        have h2 := ca_rho_sq_le hρ x
        have h3 : 0 ≤ C * (Om m m (aJet φ x)) ^ 2 := by positivity
        have h4 : (0 : ℝ) ≤ ρ x ^ 2 := sq_nonneg _
        have h5 := mul_le_mul_of_nonneg_left h1 h4
        have h6 : ρ x ^ 2 * (C * (Om m m (aJet φ x)) ^ 2) ≤ C * (Om m m (aJet φ x)) ^ 2 := by
          nlinarith
        show ρ x ^ 2 * eSym (m + 1) (jp N φ x) ≤ _
        have : jet (m + 1 + 3) φ x = jet (m + 4) φ x := rfl
        rw [this] at h5
        nlinarith
    _ = _ := by
        rw [integral_add (hJ.const_mul C) (hO.const_mul C), integral_const_mul,
          integral_const_mul]

/-- (G2) The energy two orders up is controlled by the dissipation and lower energies. -/
theorem exists_const_eSym_top_le {n N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (hK : 0 ≤ K) (Λ : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∫ x, ρ x ^ 2 * eSym (k + 2) (jp N φ x) ≤
        C * (∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2) + C * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨C, hC0, hC⟩ := ca_main (n := n) (N := N) (m := k + 1) (by omega) (K := K)
  obtain ⟨C2, hC2, h2⟩ := exists_const_integral_Om_sq (n := n) (k + 1) (k + 1)
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K
  refine ⟨C + C * C2, by positivity, fun ρ hρ φ hφ hadm => ?_⟩
  have h1 := hC Λ ρ hρ φ hφ hadm
  have hbd : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K := fun x hx =>
    ce_jet3_le (by omega) (hadm x hx)
  have h3 := h2 φ hφ hbd
  have hI : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ρ x ^ 2 * ‖jet (k + 1 + 4) φ x‖ ^ 2 =
      ∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2 :=
    (ca_integral_eq hρ (fun x => ‖jet (k + 5) φ x‖ ^ 2)).symm
  have hI0 : 0 ≤ ∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2 :=
    integral_nonneg fun x => mul_nonneg (sq_nonneg _) (sq_nonneg _)
  have hE0 := energy3_nonneg (k + 1) 3 φ
  rw [hI] at h1
  have h4 := mul_le_mul_of_nonneg_left h3 hC0
  nlinarith [mul_nonneg (mul_nonneg hC0 hC2) hI0, mul_nonneg (mul_nonneg hC0 hC2) hE0]

/-- (G3) The energy one order up is controlled by a small multiple of the top derivatives and lower energies. -/
theorem exists_const_eSym_next_le {n N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (hK : 0 ≤ K) (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∫ x, ρ x ^ 2 * eSym (k + 1) (jp N φ x) ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) φ x‖ ^ 2) +
          C * (1 + energy3 k 3 φ) := by
  obtain ⟨C, hC0, hC⟩ := ca_main (n := n) (N := N) (m := k) (by omega) (K := K)
  obtain ⟨C2, hC2, h2⟩ := exists_const_integral_Om_sq (n := n) k k
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K
  obtain ⟨C3, hC3, h3⟩ := exists_const_jet_three_orders (V := ℝ) (n := n) (ρ := 2) (R := 3)
    (ε := ε / (C + 1)) (by norm_num) (by norm_num) (by positivity)
  refine ⟨C * C3 + C * C2, by positivity, fun ρ hρ φ hφ hadm => ?_⟩
  have h1 := hC Λ ρ hρ φ hφ hadm
  have hbd : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K := fun x hx =>
    ce_jet3_le (by omega) (hadm x hx)
  have h4 := h2 φ hφ hbd
  have h5 := h3 (k + 3) φ hφ
  have hE0 := energy3_nonneg k 3 φ
  -- the cut-off integral is at most the unweighted one
  have h6 : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ρ x ^ 2 * ‖jet (k + 4) φ x‖ ^ 2 ≤
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 4) φ x‖ ^ 2 := by
    refine setIntegral_mono_on
      (integrableOn_ball_of_continuous ((hρ.smooth.continuous.pow 2).mul
        ((contDiff_jet hφ (k + 4)).continuous.norm.pow 2)) 2)
      (integrableOn_ball_of_continuous ((contDiff_jet hφ (k + 4)).continuous.norm.pow 2) 2)
      measurableSet_ball fun x _ => ?_
    have := ca_rho_sq_le hρ x
    nlinarith [sq_nonneg ‖jet (k + 4) φ x‖]
  -- the lowest term is a term of `energy3 k 3 φ`
  have h7 : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 3) φ x‖ ^ 2 ≤
      energy3 k 3 φ := by
    unfold energy3
    have := Finset.single_le_sum (f := fun i => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      ‖jet (3 + i) φ x‖ ^ 2) (s := Finset.range (k + 1))
      (fun i _ => setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _)
      (Finset.self_mem_range_succ k)
    rw [add_comm k 3]
    exact this
  have h8 : 0 ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) φ x‖ ^ 2 :=
    setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _
  have h5' : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 4) φ x‖ ^ 2 ≤
      ε / (C + 1) * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 5) φ x‖ ^ 2) +
        C3 * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet (k + 3) φ x‖ ^ 2) := h5
  have hq : C * (ε / (C + 1)) ≤ ε := by
    rw [mul_div_assoc', div_le_iff₀ (by positivity)]
    nlinarith
  have e1 := mul_le_mul_of_nonneg_left (h6.trans h5') hC0
  have e2 := mul_le_mul_of_nonneg_left h4 hC0
  have e3 := mul_le_mul_of_nonneg_left h7 (mul_nonneg hC0 hC3)
  have e4 := mul_le_mul_of_nonneg_right hq h8
  nlinarith [mul_nonneg (mul_nonneg hC0 hC3) hE0, mul_nonneg (mul_nonneg hC0 hC2) hE0]

/-- info: 'CWFlow.Extension.exists_const_eSym_top_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_eSym_top_le

/-- info: 'CWFlow.Extension.exists_const_eSym_next_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_eSym_next_le

end CWFlow.Extension
