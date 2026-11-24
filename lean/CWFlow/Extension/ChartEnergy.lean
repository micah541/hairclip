import CWFlow.Extension.EvolPointwise

/-!
# Chart energy: derivatives of order `3, …, m+3` controlled by the geometric energies

In a graph chart where the jets of order `≤ 3` of `φ` are bounded (admissibility on a ball), the
`L²` norms on a smaller ball of `D^{3}φ, …, D^{m+3}φ` are bounded by `C (1 + Σ_{l ≤ m} ∫ eSym l)`.
This is the chart part of Lemma 4.13 of Chen–Warren.  Admissibility is assumed on `ball 0 (R+1)` so
that the energy densities are continuous, hence integrable, on `ball 0 R`.
-/

open CWFlow.Interp Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-- The energy densities along `jp N φ` are integrable on `ball 0 R` if admissibility holds on
`ball 0 (R+1)`. -/
theorem ce_integrableOn_eSym {K : ℝ} {l : ℕ} (hN : l + 3 ≤ N) {R : ℝ}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hadm : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (R + 1), jp N φ x ∈ admSet n N K) :
    IntegrableOn (fun x => eSym l (jp N φ x)) (ball (0 : EuclideanSpace ℝ (Fin n)) R) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hc : ContinuousOn (fun x => eSym l (jp N φ x)) (ball (0 : EuclideanSpace ℝ (Fin n)) (R + 1)) :=
    (contDiffOn_symbol_jp (isSymbol_eSym hA hlow l hN) hφ).continuousOn.mono hadm
  exact (hc.mono (closedBall_subset_ball (by linarith))).integrableOn_compact
    (isCompact_closedBall 0 R) |>.mono_set ball_subset_closedBall

/-- The jets of order `3` are bounded by `K` on admissible points. -/
theorem ce_jet3_le {K : ℝ} (hN : 3 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    {x : EuclideanSpace ℝ (Fin n)} (hx : jp N φ x ∈ admSet n N K) : ‖jet 3 φ x‖ ≤ K := by
  have := hx ⟨3, by omega⟩ le_rfl
  exact this.le

/-- Pointwise comparison for the top derivative. -/
theorem ce_pointwise {K : ℝ} {k : ℕ} (hk : 1 ≤ k) (hN : k + 3 ≤ N) {Cs : ℝ}
    (hCs : ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ admSet n N K →
      ‖sSym k (jp N φ x)‖ ≤ Cs * Om (k - 1) (k - 1)
        (fun j => if 3 + j ≤ N then ‖jet (3 + j) φ x‖ else 0))
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n))
    (hx : jp N φ x ∈ admSet n N K) :
    ‖jet (k + 3) φ x‖ ^ 2 ≤ 2 * (1 + K ^ 2) ^ (k + 3) * eSym k (jp N φ x) +
      2 * Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2 := by
  set p := jp N φ x with hp
  have hjet : jet (k + 3) φ x = tSym k p - sSym k p := by
    have := slotN_jp (n := n) (k + 3) hN φ x
    unfold sSym
    rw [← hp] at this
    rw [this]; abel
  -- tSym bound
  have hu : ‖(low23 p).1‖ ≤ K := by
    have h2 : 2 < N + 1 := by omega
    rw [low23_fst p h2]
    exact (hx ⟨2, h2⟩ (by norm_num)).le
  have hK0 : 0 ≤ K := (norm_nonneg _).trans hu
  set B : ℝ := 1 + K ^ 2 with hB
  have hBpos : 0 < B := by positivity
  have hle := le_eSym (n := n) (N := N) k p
  have hB1 : 1 + ‖(low23 p).1‖ ^ 2 ≤ B := by
    have : ‖(low23 p).1‖ ^ 2 ≤ K ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hu 2
    linarith
  have hinv : B⁻¹ ≤ (1 + ‖(low23 p).1‖ ^ 2)⁻¹ := inv_anti₀ (by positivity) hB1
  have hpow : (B⁻¹) ^ (k + 3) ≤ ((1 + ‖(low23 p).1‖ ^ 2)⁻¹) ^ (k + 3) :=
    pow_le_pow_left₀ (by positivity) hinv _
  have h1 : (B⁻¹) ^ (k + 3) * ‖tSym k p‖ ^ 2 ≤ eSym k p :=
    (mul_le_mul_of_nonneg_right hpow (by positivity)).trans hle
  have hBB : B ^ (k + 3) * (B⁻¹) ^ (k + 3) = 1 := by
    rw [← mul_pow, mul_inv_cancel₀ hBpos.ne', one_pow]
  have ht : ‖tSym k p‖ ^ 2 ≤ B ^ (k + 3) * eSym k p := by
    calc ‖tSym k p‖ ^ 2 = B ^ (k + 3) * ((B⁻¹) ^ (k + 3) * ‖tSym k p‖ ^ 2) := by
          rw [← mul_assoc, hBB, one_mul]
      _ ≤ B ^ (k + 3) * eSym k p := mul_le_mul_of_nonneg_left h1 (by positivity)
  -- sSym bound
  have hs := hCs φ x hx
  have hOm : Om (k - 1) (k - 1) (fun j => if 3 + j ≤ N then ‖jet (3 + j) φ x‖ else 0) =
      Om (k - 1) (k - 1) (aJet φ x) := by
    refine evp_Om_congr fun j _ hj2 => ?_
    have : 3 + j ≤ N := by omega
    simp [aJet, this]
  rw [hOm] at hs
  have hs' : ‖sSym k p‖ ≤ |Cs| * Om (k - 1) (k - 1) (aJet φ x) :=
    hs.trans (mul_le_mul_of_nonneg_right (le_abs_self _)
      (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x))))
  have hs2 : ‖sSym k p‖ ^ 2 ≤ Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2 := by
    have := pow_le_pow_left₀ (norm_nonneg _) hs' 2
    rwa [mul_pow, sq_abs] at this
  have hsub : ‖jet (k + 3) φ x‖ ≤ ‖tSym k p‖ + ‖sSym k p‖ := by
    rw [hjet]; exact norm_sub_le _ _
  have hsq : ‖jet (k + 3) φ x‖ ^ 2 ≤ 2 * ‖tSym k p‖ ^ 2 + 2 * ‖sSym k p‖ ^ 2 := by
    have h0 := norm_nonneg (jet (k + 3) φ x)
    nlinarith [sq_nonneg (‖tSym k p‖ - ‖sSym k p‖), norm_nonneg (tSym k p),
      norm_nonneg (sSym k p), mul_self_le_mul_self h0 hsub]
  calc ‖jet (k + 3) φ x‖ ^ 2 ≤ 2 * ‖tSym k p‖ ^ 2 + 2 * ‖sSym k p‖ ^ 2 := hsq
    _ ≤ 2 * (B ^ (k + 3) * eSym k p) + 2 * (Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2) := by
        linarith
    _ = _ := by ring

/-- One step: the top derivative on `ball ρ` is bounded by `eSym k` and the energy of lower order. -/
theorem ce_step {k : ℕ} (hk : 1 ≤ k) (hN : k + 3 ≤ N) {ρ r' R : ℝ} (hρ : 0 < ρ) (hρr : ρ < r')
    (hr'R : r' ≤ R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, jp N φ x ∈ admSet n N K) →
      IntegrableOn (fun x => eSym k (jp N φ x)) (ball (0 : EuclideanSpace ℝ (Fin n)) R) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet (k + 3) φ x‖ ^ 2 ≤
        C * ((∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym k (jp N φ x)) + 1 +
          energy3 (k - 1) r' φ) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  obtain ⟨Cs, hCs⟩ := norm_symbol_jp_le (isSymbol_sSym hA hlow k hk hN)
  obtain ⟨C2, hC2, h2⟩ := exists_const_integral_Om_sq (n := n) (k - 1) (k - 1) hρ hρr K
  set D : ℝ := (1 + K ^ 2) ^ (k + 3) with hD
  have hD0 : 0 ≤ D := by positivity
  refine ⟨2 * D + 2 * Cs ^ 2 * C2, by positivity, fun φ hφ hadm hint => ?_⟩
  have hρR : ρ ≤ R := (hρr.trans_le hr'R).le
  have hsub : ball (0 : EuclideanSpace ℝ (Fin n)) ρ ⊆ ball 0 R := ball_subset_ball hρR
  have hint' : IntegrableOn (fun x => eSym k (jp N φ x)) (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
    hint.mono_set hsub
  have hOmint : IntegrableOn (fun x => (Om (k - 1) (k - 1) (aJet φ x)) ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
    integrableOn_ball_of_continuous ((continuous_Om hφ _ _).pow 2) ρ
  have hpt : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
      ‖jet (k + 3) φ x‖ ^ 2 ≤ 2 * D * eSym k (jp N φ x) +
        2 * Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2 := fun x hx =>
    ce_pointwise hk hN hCs φ x (hadm x (hsub hx))
  have hbd : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) r', ‖jet 3 φ x‖ ≤ K := fun x hx =>
    ce_jet3_le (by omega) (hadm x (ball_subset_ball hr'R hx))
  have hOm2 := h2 φ hφ hbd
  have hI0 : 0 ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym k (jp N φ x) :=
    setIntegral_nonneg measurableSet_ball fun x _ => eSym_nonneg _ _
  have hE0 := energy3_nonneg (k - 1) r' φ
  have hmono : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, eSym k (jp N φ x) ≤
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym k (jp N φ x) :=
    setIntegral_mono_set hint (Filter.Eventually.of_forall fun x => eSym_nonneg _ _)
      (ball_subset_ball hρR).eventuallyLE
  have hlhs : IntegrableOn (fun x => ‖jet (k + 3) φ x‖ ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
    integrableOn_ball_of_continuous ((contDiff_jet hφ (k + 3)).continuous.norm.pow 2) ρ
  have hrhs : IntegrableOn (fun x => 2 * D * eSym k (jp N φ x) +
      2 * Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
    (hint'.const_mul (2 * D)).add (hOmint.const_mul (2 * Cs ^ 2))
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet (k + 3) φ x‖ ^ 2
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, (2 * D * eSym k (jp N φ x) +
          2 * Cs ^ 2 * (Om (k - 1) (k - 1) (aJet φ x)) ^ 2) :=
        setIntegral_mono_on hlhs hrhs measurableSet_ball hpt
    _ = 2 * D * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, eSym k (jp N φ x)) +
          2 * Cs ^ 2 * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
            (Om (k - 1) (k - 1) (aJet φ x)) ^ 2 := by
        rw [integral_add (hint'.const_mul _) (hOmint.const_mul _), integral_const_mul,
          integral_const_mul]
    _ ≤ _ := by
        have h3 : 2 * D * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, eSym k (jp N φ x)) ≤
            2 * D * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym k (jp N φ x) :=
          mul_le_mul_of_nonneg_left hmono (by positivity)
        have h4 : 2 * Cs ^ 2 * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
            (Om (k - 1) (k - 1) (aJet φ x)) ^ 2 ≤ 2 * Cs ^ 2 * (C2 * (1 + energy3 (k - 1) r' φ)) :=
          mul_le_mul_of_nonneg_left hOm2 (by positivity)
        have h5 : 0 ≤ 2 * Cs ^ 2 * C2 := by positivity
        nlinarith [mul_nonneg hD0 hE0, mul_nonneg h5 hE0, mul_nonneg h5 hI0, mul_nonneg hD0 hI0]

/-- Induction on the order. -/
theorem ce_induct {K : ℝ} (hK : 0 ≤ K) (m : ℕ) (hN : m + 3 ≤ N) {ρ R : ℝ} (hρ : 0 < ρ)
    (hρR : ρ < R) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (R + 1), jp N φ x ∈ admSet n N K) →
      energy3 m ρ φ ≤
        C * (1 + ∑ l ∈ Finset.range (m + 1),
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x)) := by
  induction m generalizing ρ R with
  | zero =>
    refine ⟨K ^ 2 * volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ), by positivity,
      fun φ hφ hadm => ?_⟩
    have hsum : 0 ≤ ∑ l ∈ Finset.range (0 + 1),
        ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x) :=
      Finset.sum_nonneg fun l _ => setIntegral_nonneg measurableSet_ball fun x _ =>
        eSym_nonneg _ _
    have hsub : ball (0 : EuclideanSpace ℝ (Fin n)) ρ ⊆ ball 0 (R + 1) :=
      ball_subset_ball (by linarith)
    have hb : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet 3 φ x‖ ^ 2 ≤
        K ^ 2 * volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) := by
      calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet 3 φ x‖ ^ 2
          ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, K ^ 2 :=
            setIntegral_mono_on
              (integrableOn_ball_of_continuous ((contDiff_jet hφ 3).continuous.norm.pow 2) ρ)
              (integrableOn_const_ball _ ρ) measurableSet_ball fun x hx =>
              pow_le_pow_left₀ (norm_nonneg _) (ce_jet3_le (by omega) (hadm x (hsub hx))) 2
        _ = _ := by
            rw [setIntegral_const, smul_eq_mul, mul_comm]
    have he : energy3 0 ρ φ = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet 3 φ x‖ ^ 2 := by
      unfold energy3
      simp
    rw [he]
    have hV : 0 ≤ K ^ 2 * volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) := by positivity
    nlinarith
  | succ m ih =>
    have hr0 : 0 < (ρ + R) / 2 := by linarith
    obtain ⟨C1, hC1, h1⟩ := ih (by omega) hρ hρR
    obtain ⟨C2, hC2, h2⟩ := ih (by omega) hr0 (show (ρ + R) / 2 < R by linarith)
    obtain ⟨Cs, hCs, hs⟩ := ce_step (N := N) (n := n) (k := m + 1) (by omega) hN (ρ := ρ)
      (r' := (ρ + R) / 2) (R := R) hρ (by linarith) (by linarith) K
    refine ⟨C1 + Cs + Cs * C2, by positivity, fun φ hφ hadm => ?_⟩
    have hadmR : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, jp N φ x ∈ admSet n N K :=
      fun x hx => hadm x (ball_subset_ball (by linarith) hx)
    have hint := ce_integrableOn_eSym (K := K) (l := m + 1) hN hφ hadm
    have hs' := hs φ hφ hadmR hint
    rw [Nat.add_sub_cancel] at hs'
    have hA := h1 φ hφ hadm
    have hB := h2 φ hφ hadm
    have he : energy3 (m + 1) ρ φ = energy3 m ρ φ +
        ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet (m + 1 + 3) φ x‖ ^ 2 := by
      unfold energy3
      rw [Finset.sum_range_succ _ (m + 1), show 3 + (m + 1) = m + 1 + 3 from by omega]
    rw [he, Finset.sum_range_succ _ (m + 1)]
    set S := ∑ l ∈ Finset.range (m + 1),
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x) with hS
    set I := ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym (m + 1) (jp N φ x) with hI
    have hS0 : 0 ≤ S := Finset.sum_nonneg fun l _ => setIntegral_nonneg measurableSet_ball
      fun x _ => eSym_nonneg _ _
    have hI0 : 0 ≤ I := setIntegral_nonneg measurableSet_ball fun x _ => eSym_nonneg _ _
    have h3 : Cs * energy3 m ((ρ + R) / 2) φ ≤ Cs * (C2 * (1 + S)) :=
      mul_le_mul_of_nonneg_left hB hCs
    nlinarith [mul_nonneg hC1 hI0, mul_nonneg hCs hS0, mul_nonneg (mul_nonneg hCs hC2) hI0]

/-- Chart part of CW Lemma 4.13: the `L²` norms of `D^3 φ, …, D^{m+3} φ` on `ball 0 ρ` are bounded by
`C (1 + Σ_{l ≤ m} ∫_{ball 0 R} eSym l)` when the jets of order `≤ 3` are bounded by `K` (admissibility
on `ball 0 (R + 1)`). -/
theorem exists_const_energy3_le_eSym {N : ℕ} (m : ℕ) (hN : m + 3 ≤ N) {ρ R : ℝ} (hρ : 0 < ρ)
    (hρR : ρ < R) {K : ℝ} (hK : 0 ≤ K) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (R + 1), jp N φ x ∈ admSet n N K) →
      energy3 m ρ φ ≤
        C * (1 + ∑ l ∈ Finset.range (m + 1),
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, eSym l (jp N φ x)) :=
  ce_induct hK m hN hρ hρR

/-- info: 'CWFlow.Extension.exists_const_energy3_le_eSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_energy3_le_eSym

end CWFlow.Extension

