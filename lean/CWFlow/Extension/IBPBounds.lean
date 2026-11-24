import CWFlow.Extension.EvolPointwise
import CWFlow.Extension.LocalWeight
import CWFlow.Extension.JetIBP

/-!
# Integral bounds for jet coordinates against localized weights

For an admissible smooth `φ` and a cutoff `ρ`, three estimates for `∫ (jet coordinate) · ρ² 𝔣(jets of φ)`:
no integration by parts (odd weight), one integration by parts, and two integrations by parts with
the top term `ρ² ∂ₐ(∂_c 𝔣)` kept.  They combine integration by parts (`JetIBP`), pointwise bounds of
symbols by weight functions, Cauchy–Schwarz and the integral bounds of `WeightIntegral`.
-/

open CWFlow.Interp Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-- A cutoff supported in the ball of radius 2 with bounds on its first two derivatives. -/
structure IsCutoff (ρ : EuclideanSpace ℝ (Fin n) → ℝ) (Λ : ℝ) : Prop where
  smooth : ContDiff ℝ ∞ ρ
  supp : tsupport ρ ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 2
  le_one : ∀ x, |ρ x| ≤ 1
  d1 : ∀ x, ‖fderiv ℝ ρ x‖ ≤ Λ
  d2 : ∀ x, ‖fderiv ℝ (fderiv ℝ ρ) x‖ ≤ Λ

theorem ibb_hasCompactSupport {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (h : IsCutoff ρ Λ) :
    HasCompactSupport ρ :=
  IsCompact.of_isClosed_subset (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) 2)
    (isClosed_tsupport ρ) (h.supp.trans ball_subset_closedBall)

theorem ibb_nonneg {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (h : IsCutoff ρ Λ) : 0 ≤ Λ :=
  (norm_nonneg _).trans (h.d1 0)

/-- Pointwise bound of a symbol along jets by a weight function of `aJet`. -/
theorem ibb_pt {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)} {w ℓ : ℕ}
    {𝔤 : P n N → W} (h : IsSymbol A w ℓ 𝔤) (hN : 3 + min w ℓ ≤ N) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ A → ‖𝔤 (jp N φ x)‖ ≤ C * Om w ℓ (aJet φ x) := by
  obtain ⟨C, hC⟩ := h.norm_le
  refine ⟨|C|, abs_nonneg _, fun φ x hx => ?_⟩
  have h1 := hC _ hx
  have e : Om w ℓ (aOf (jp N φ x)) = Om w ℓ (aJet φ x) := by
    refine evp_Om_congr fun j _ h2 => ?_
    rw [aOf_jp, if_pos (by omega)]
    rfl
  rw [e] at h1
  exact h1.trans (mul_le_mul_of_nonneg_right (le_abs_self C)
    (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x))))

theorem ibb_norm3 {K : ℝ} (hN : 3 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    {x : EuclideanSpace ℝ (Fin n)} (hx : jp N φ x ∈ admSet n N K) : ‖jet 3 φ x‖ ≤ K := by
  have := hx ⟨3, by omega⟩ (by simp)
  exact this.le

theorem ibb_ball_sub (x : EuclideanSpace ℝ (Fin n)) (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2) :
    x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3 :=
  ball_subset_ball (by norm_num) hx

/-- (1) No integration by parts, odd weight. -/
theorem exists_const_integral_locWeight_odd {N k : ℕ} (hN : k + 5 ≤ N) {K : ℝ} {𝔣 : P n N → ℝ}
    (h𝔣 : IsSymbol (admSet n N K) (2 * k + 3) (k + 2) 𝔣) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (Λ : ℝ) (ρ : EuclideanSpace ℝ (Fin n) → ℝ), IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, locWeight ρ 𝔣 φ x| ≤ ε * energy3 (k + 2) 3 φ + C := by
  obtain ⟨C₀, hC₀, hpt⟩ := ibb_pt h𝔣 (by omega : 3 + min (2 * k + 3) (k + 2) ≤ N)
  obtain ⟨C', hC', hint⟩ := exists_const_integral_Om_odd (n := n) (m := k + 2) (by omega)
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K (ε := ε / (C₀ + 1)) (by positivity)
  refine ⟨C₀ * C', mul_nonneg hC₀ hC', fun Λ ρ hρ φ hφ hadm => ?_⟩
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hw := contDiff_locWeight h𝔣 hφ hρ.smooth hsupp
  have hom : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) ≤
      ε / (C₀ + 1) * energy3 (k + 2) 3 φ + C' := by
    have := hint φ hφ hK
    have e : 2 * (k + 2) - 1 = 2 * k + 3 := by omega
    rwa [e] at this
  have hb : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      |locWeight ρ 𝔣 φ x| ≤ C₀ * Om (2 * k + 3) (k + 2) (aJet φ x) := by
    intro x _
    by_cases hxs : x ∈ tsupport ρ
    · have h1 := hpt φ x (hsupp x hxs)
      rw [Real.norm_eq_abs] at h1
      have h2 : |ρ x ^ 2| ≤ 1 := by
        rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) (hρ.le_one x)
      unfold locWeight
      rw [abs_mul]
      calc |ρ x ^ 2| * |𝔣 (jp N φ x)| ≤ 1 * (C₀ * Om (2 * k + 3) (k + 2) (aJet φ x)) :=
            mul_le_mul h2 h1 (abs_nonneg _) zero_le_one
        _ = _ := one_mul _
    · rw [locWeight_eq_zero 𝔣 φ hxs, abs_zero]
      exact mul_nonneg hC₀ (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x)))
  have hz : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → locWeight ρ 𝔣 φ x = 0 :=
    fun x hx => locWeight_eq_zero 𝔣 φ (fun h => hx (hρ.supp h))
  have hq : C₀ * (ε / (C₀ + 1)) ≤ ε := by
    calc C₀ * (ε / (C₀ + 1)) = ε * (C₀ / (C₀ + 1)) := by ring
      _ ≤ ε * 1 := mul_le_mul_of_nonneg_left (div_le_one_of_le₀ (by linarith) (by positivity)) hε.le
      _ = ε := mul_one _
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero hz]
  calc |∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, locWeight ρ 𝔣 φ x|
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, |locWeight ρ 𝔣 φ x| := by
        simpa [Real.norm_eq_abs] using norm_integral_le_integral_norm (μ := volume.restrict
          (ball (0 : EuclideanSpace ℝ (Fin n)) 2)) (locWeight ρ 𝔣 φ)
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, C₀ * Om (2 * k + 3) (k + 2) (aJet φ x) :=
        setIntegral_mono_on (integrableOn_ball_of_continuous hw.continuous 2).abs
          ((integrableOn_ball_of_continuous (continuous_Om hφ _ _) 2).const_mul _)
          measurableSet_ball hb
    _ = C₀ * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) :=
        integral_const_mul _ _
    _ ≤ C₀ * (ε / (C₀ + 1) * energy3 (k + 2) 3 φ + C') := mul_le_mul_of_nonneg_left hom hC₀
    _ ≤ ε * energy3 (k + 2) 3 φ + C₀ * C' := by
        nlinarith [mul_le_mul_of_nonneg_right hq (energy3_nonneg (n := n) (k + 2) 3 φ)]

/-- Cauchy–Schwarz / Young step shared by the two integration-by-parts estimates. -/
theorem ibb_young {k : ℕ} {ε C₁ C₂ : ℝ} (hε : 0 < ε) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ)
    (hom : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, (Om (k + 1) (k + 1) (aJet φ x)) ^ 2 ≤
      C₂ * (1 + energy3 (k + 1) 3 φ))
    {h : EuclideanSpace ℝ (Fin n) → ℝ} (hcont : Continuous h)
    (hz : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → h x = 0)
    (hb : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      |h x| ≤ ‖jet (k + 5) φ x‖ * (C₁ * Om (k + 1) (k + 1) (aJet φ x))) :
    |∫ x, h x| ≤ ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
      C₁ ^ 2 / (4 * ε) * (C₂ * (1 + energy3 (k + 1) 3 φ)) := by
  have hJ : Continuous fun x => ‖jet (k + 5) φ x‖ ^ 2 :=
    (contDiff_jet hφ (k + 5)).continuous.norm.pow 2
  have hO : Continuous fun x => (Om (k + 1) (k + 1) (aJet φ x)) ^ 2 :=
    (continuous_Om hφ (k + 1) (k + 1)).pow 2
  have hpt : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, |h x| ≤
      ε * ‖jet (k + 5) φ x‖ ^ 2 + C₁ ^ 2 / (4 * ε) * (Om (k + 1) (k + 1) (aJet φ x)) ^ 2 := by
    intro x hx
    refine (hb x hx).trans ?_
    set J := ‖jet (k + 5) φ x‖
    set Ω := Om (k + 1) (k + 1) (aJet φ x)
    rw [← sub_nonneg]
    have : ε * J ^ 2 + C₁ ^ 2 / (4 * ε) * Ω ^ 2 - J * (C₁ * Ω) = (2 * ε * J - C₁ * Ω) ^ 2 / (4 * ε) := by
      field_simp
      ring
    rw [this]
    positivity
  have hint2 : IntegrableOn (fun x => ε * ‖jet (k + 5) φ x‖ ^ 2 +
      C₁ ^ 2 / (4 * ε) * (Om (k + 1) (k + 1) (aJet φ x)) ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    (continuous_const.mul hJ |>.add (continuous_const.mul hO)).continuousOn.integrableOn_compact
      (isCompact_closedBall 0 2) |>.mono_set ball_subset_closedBall
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero hz]
  calc |∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, h x|
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, |h x| := by
        simpa [Real.norm_eq_abs] using norm_integral_le_integral_norm (μ := volume.restrict
          (ball (0 : EuclideanSpace ℝ (Fin n)) 2)) h
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2,
          (ε * ‖jet (k + 5) φ x‖ ^ 2 + C₁ ^ 2 / (4 * ε) * (Om (k + 1) (k + 1) (aJet φ x)) ^ 2) :=
        setIntegral_mono_on (integrableOn_ball_of_continuous hcont 2).abs hint2
          measurableSet_ball hpt
    _ = ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
          C₁ ^ 2 / (4 * ε) * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2,
            (Om (k + 1) (k + 1) (aJet φ x)) ^ 2 := by
        rw [integral_add ((integrableOn_ball_of_continuous hJ 2).const_mul ε)
          ((integrableOn_ball_of_continuous hO 2).const_mul (C₁ ^ 2 / (4 * ε))), integral_const_mul,
          integral_const_mul]
    _ ≤ _ := by
        have : 0 ≤ C₁ ^ 2 / (4 * ε) := by positivity
        nlinarith [mul_le_mul_of_nonneg_left hom this]

/-- A coordinate of a jet is bounded by its norm. -/
theorem ibb_coord_le {φ : EuclideanSpace ℝ (Fin n) → ℝ} (m : ℕ) (x : EuclideanSpace ℝ (Fin n))
    (I : Fin m → Fin n) : |jet m φ x I| ≤ ‖jet m φ x‖ := by
  have := PiLp.norm_apply_le (jet m φ x) I
  rwa [Real.norm_eq_abs] at this

theorem ibb_coord_contDiff {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (m : ℕ) (I : Fin m → Fin n) : ContDiff ℝ ∞ (fun y => jet m φ y I) :=
  (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin m → Fin n => ℝ) I).contDiff.comp (contDiff_jet hφ m)

/-- A symbol and its total derivative are bounded by the weight `Ω_{k+1}^{≤k+1}`. -/
theorem ibb_f_td {k : ℕ} (hN : k + 4 ≤ N) {K : ℝ} {𝔣 : P n N → ℝ}
    (h𝔣 : IsSymbol (admSet n N K) k k 𝔣) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ admSet n N K →
      |𝔣 (jp N φ x)| ≤ C * Om (k + 1) (k + 1) (aJet φ x) ∧
      ‖totalD 𝔣 (jp N φ x)‖ ≤ C * Om (k + 1) (k + 1) (aJet φ x) := by
  have ht : IsSymbol (admSet n N K) (k + 1) (k + 1) (totalD 𝔣) := by
    have := h𝔣.totalD (by omega) (admSet_low K)
    rwa [max_self] at this
  obtain ⟨C₀, hC₀, h0⟩ := ibb_pt h𝔣 (by omega : 3 + min k k ≤ N)
  obtain ⟨C₁, hC₁, h1⟩ := ibb_pt ht (by omega : 3 + min (k + 1) (k + 1) ≤ N)
  refine ⟨C₀ * (k + 1) + C₁, by positivity, fun φ x hx => ?_⟩
  have hO0 : 0 ≤ Om (k + 1) (k + 1) (aJet φ x) :=
    zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x))
  have hm := Om_mono (v := k) (v' := k + 1) (ℓ := k) (ℓ' := k + 1) (by omega) (by omega)
    (aJet_nonneg φ x)
  have hm' : (0 : ℝ) ≤ (k : ℝ) + 1 := by positivity
  have a0 := h0 φ x hx
  have a1 := h1 φ x hx
  rw [Real.norm_eq_abs] at a0
  constructor
  · calc |𝔣 (jp N φ x)| ≤ C₀ * Om k k (aJet φ x) := a0
      _ ≤ C₀ * (((k : ℝ) + 1) * Om (k + 1) (k + 1) (aJet φ x)) :=
        mul_le_mul_of_nonneg_left hm hC₀
      _ ≤ _ := by nlinarith [mul_nonneg hC₁ hO0]
  · calc ‖totalD 𝔣 (jp N φ x)‖ ≤ C₁ * Om (k + 1) (k + 1) (aJet φ x) := a1
      _ ≤ _ := by nlinarith [mul_nonneg (mul_nonneg hC₀ hm') hO0]

theorem ibb_d1 {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    (x : EuclideanSpace ℝ (Fin n)) (a : Fin n) :
    |fderiv ℝ ρ x (EuclideanSpace.single a 1)| ≤ |Λ| := by
  calc |fderiv ℝ ρ x (EuclideanSpace.single a 1)| = ‖fderiv ℝ ρ x (EuclideanSpace.single a 1)‖ :=
        (Real.norm_eq_abs _).symm
    _ ≤ ‖fderiv ℝ ρ x‖ * ‖(EuclideanSpace.single a (1 : ℝ) : EuclideanSpace ℝ (Fin n))‖ :=
        (fderiv ℝ ρ x).le_opNorm _
    _ = ‖fderiv ℝ ρ x‖ := by simp
    _ ≤ |Λ| := (hρ.d1 x).trans (le_abs_self Λ)

theorem ibb_d2 {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    (x : EuclideanSpace ℝ (Fin n)) (a c : Fin n) :
    |fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1) (EuclideanSpace.single c 1)| ≤ |Λ| := by
  calc |fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1) (EuclideanSpace.single c 1)|
      = ‖fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1) (EuclideanSpace.single c 1)‖ :=
        (Real.norm_eq_abs _).symm
    _ ≤ ‖fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1)‖ *
          ‖(EuclideanSpace.single c (1 : ℝ) : EuclideanSpace ℝ (Fin n))‖ :=
        (fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1)).le_opNorm _
    _ = ‖fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1)‖ := by simp
    _ ≤ ‖fderiv ℝ (fderiv ℝ ρ) x‖ * ‖(EuclideanSpace.single a (1 : ℝ) : EuclideanSpace ℝ (Fin n))‖ :=
        (fderiv ℝ (fderiv ℝ ρ) x).le_opNorm _
    _ = ‖fderiv ℝ (fderiv ℝ ρ) x‖ := by simp
    _ ≤ |Λ| := (hρ.d2 x).trans (le_abs_self Λ)

theorem ibb_norm_td (a : Fin n)
    (T : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ) {B : ℝ} (hT : ‖T‖ ≤ B) :
    |T (EuclideanSpace.single a 1)| ≤ B := by
  calc |T (EuclideanSpace.single a 1)| = ‖T (EuclideanSpace.single a 1)‖ := (Real.norm_eq_abs _).symm
    _ ≤ ‖T‖ * ‖(EuclideanSpace.single a (1 : ℝ) : EuclideanSpace ℝ (Fin n))‖ := T.le_opNorm _
    _ = ‖T‖ := by simp
    _ ≤ B := hT

/-- First derivative of the localized weight in a coordinate direction. -/
theorem ibb_wfd {k : ℕ} (hN : k + 4 ≤ N) {K : ℝ} {𝔣 : P n N → ℝ}
    (h𝔣 : IsSymbol (admSet n N K) k k 𝔣) {Cf : ℝ}
    (hf : ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ admSet n N K →
      |𝔣 (jp N φ x)| ≤ Cf * Om (k + 1) (k + 1) (aJet φ x) ∧
      ‖totalD 𝔣 (jp N φ x)‖ ≤ Cf * Om (k + 1) (k + 1) (aJet φ x))
    {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)}
    (hx : jp N φ x ∈ admSet n N K) (a : Fin n) :
    |fderiv ℝ (locWeight ρ 𝔣 φ) x (EuclideanSpace.single a 1)| ≤
      (2 * |Λ| + 1) * Cf * Om (k + 1) (k + 1) (aJet φ x) := by
  obtain ⟨hf1, hf2⟩ := hf φ x hx
  rw [fderiv_locWeight h𝔣 hN hφ hρ.smooth hx]
  have hρ1 := hρ.le_one x
  have hd := ibb_d1 hρ x a
  have ht := ibb_norm_td a _ hf2
  have h1 : |2 * ρ x * fderiv ℝ ρ x (EuclideanSpace.single a 1) * 𝔣 (jp N φ x)| ≤
      2 * |Λ| * (Cf * Om (k + 1) (k + 1) (aJet φ x)) := by
    have := mul_le_mul (mul_le_mul hρ1 hd (abs_nonneg _) zero_le_one) hf1 (abs_nonneg _)
      (by positivity : (0 : ℝ) ≤ 1 * |Λ|)
    rw [abs_mul, abs_mul, abs_mul, abs_two]
    linarith
  have h2 : |ρ x ^ 2 * totalD 𝔣 (jp N φ x) (EuclideanSpace.single a 1)| ≤
      Cf * Om (k + 1) (k + 1) (aJet φ x) := by
    have h3 : |ρ x ^ 2| ≤ 1 := by rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) hρ1
    rw [abs_mul]
    calc |ρ x ^ 2| * |totalD 𝔣 (jp N φ x) (EuclideanSpace.single a 1)|
        ≤ 1 * (Cf * Om (k + 1) (k + 1) (aJet φ x)) := mul_le_mul h3 ht (abs_nonneg _) zero_le_one
      _ = _ := one_mul _
  calc _ ≤ _ := abs_add_le _ _
    _ ≤ _ := add_le_add h1 h2
    _ = _ := by ring

theorem ibb_term {r d t Λ B : ℝ} (hr : |r| ≤ 1) (hd : |d| ≤ Λ) (ht : |t| ≤ B) :
    |2 * r * d * t| ≤ 2 * Λ * B := by
  have := mul_le_mul (mul_le_mul hr hd (abs_nonneg _) zero_le_one) ht (abs_nonneg _)
    (by have := (abs_nonneg d).trans hd; positivity : (0 : ℝ) ≤ 1 * Λ)
  rw [abs_mul, abs_mul, abs_mul, abs_two]
  linarith

/-- The derivative of `2 ρ ∂_c ρ` in a coordinate direction. -/
theorem ibb_u {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    (x : EuclideanSpace ℝ (Fin n)) (a c : Fin n) :
    |fderiv ℝ (fun y => 2 * ρ y * fderiv ℝ ρ y (EuclideanSpace.single c 1)) x
      (EuclideanSpace.single a 1)| ≤ 2 * |Λ| ^ 2 + 2 * |Λ| := by
  have hρd : HasFDerivAt ρ (fderiv ℝ ρ x) x :=
    (hρ.smooth.differentiable (by simp)).differentiableAt.hasFDerivAt
  have hdρ : ContDiff ℝ ∞ (fderiv ℝ ρ) := hρ.smooth.fderiv_right (by simp)
  have hd2 : HasFDerivAt (fderiv ℝ ρ) (fderiv ℝ (fderiv ℝ ρ) x) x :=
    (hdρ.differentiable (by simp)).differentiableAt.hasFDerivAt
  have hg : HasFDerivAt (fun y => fderiv ℝ ρ y (EuclideanSpace.single c 1))
      ((ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single c (1 : ℝ) :
        EuclideanSpace ℝ (Fin n))).comp (fderiv ℝ (fderiv ℝ ρ) x)) x :=
    (ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single c (1 : ℝ) :
      EuclideanSpace ℝ (Fin n))).hasFDerivAt.comp x hd2
  have hu : HasFDerivAt (fun y => 2 * ρ y * fderiv ℝ ρ y (EuclideanSpace.single c 1)) _ x :=
    (hρd.const_mul 2).mul hg
  have key : fderiv ℝ (fun y => 2 * ρ y * fderiv ℝ ρ y (EuclideanSpace.single c 1)) x
      (EuclideanSpace.single a 1) =
      2 * ρ x * fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1) (EuclideanSpace.single c 1) +
        fderiv ℝ ρ x (EuclideanSpace.single c 1) * (2 * fderiv ℝ ρ x (EuclideanSpace.single a 1)) := by
    rw [hu.fderiv]
    simp
  rw [key]
  have h1 : |2 * ρ x * fderiv ℝ (fderiv ℝ ρ) x (EuclideanSpace.single a 1)
      (EuclideanSpace.single c 1)| ≤ 2 * |Λ| := by
    rw [abs_mul, abs_mul, abs_two]
    have := mul_le_mul (hρ.le_one x) (ibb_d2 hρ x a c) (abs_nonneg _) zero_le_one
    linarith
  have hd1a := ibb_d1 hρ x a
  have hd1c := ibb_d1 hρ x c
  have h2 : |fderiv ℝ ρ x (EuclideanSpace.single c 1) * (2 * fderiv ℝ ρ x (EuclideanSpace.single a 1))|
      ≤ |Λ| * (2 * |Λ|) := by
    rw [abs_mul, abs_mul, abs_two]
    exact mul_le_mul hd1c (mul_le_mul_of_nonneg_left hd1a (by norm_num)) (by positivity)
      (abs_nonneg _)
  calc _ ≤ _ := abs_add_le _ _
    _ ≤ 2 * |Λ| + |Λ| * (2 * |Λ|) := add_le_add h1 h2
    _ = _ := by ring

theorem ibb_D2_zero {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (𝔣 : P n N → ℝ)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport ρ)
    (e e' : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y e') x e = 0 := by
  have hev : (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y e') =ᶠ[nhds x] fun _ => 0 := by
    filter_upwards [(isClosed_tsupport ρ).isOpen_compl.mem_nhds hx] with y hy
    simp [fderiv_locWeight_eq_zero 𝔣 φ hy]
  rw [hev.fderiv_eq]
  simp

/-- One integration by parts. -/
theorem exists_const_integral_jet_locWeight {N k : ℕ} (hN : k + 5 ≤ N) {K : ℝ} (Λ : ℝ)
    {𝔣 : P n N → ℝ} (h𝔣 : IsSymbol (admSet n N K) k k 𝔣) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∀ L : Fin (k + 5 + 1) → Fin n,
        |∫ x, jet (k + 5 + 1) φ x L * locWeight ρ 𝔣 φ x| ≤
          ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
          C * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨Cf, hCf, hf⟩ := ibb_f_td (by omega : k + 4 ≤ N) h𝔣
  obtain ⟨C₂, hC₂, hom⟩ := exists_const_integral_Om_sq (n := n) (k + 1) (k + 1) (ρ := 2) (R := 3)
    (by norm_num) (by norm_num) K
  refine ⟨((2 * |Λ| + 1) * Cf) ^ 2 / (4 * ε) * C₂, by positivity,
    fun ρ hρ φ hφ hadm L => ?_⟩
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hw := contDiff_locWeight h𝔣 hφ hρ.smooth hsupp
  have hc := hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) 𝔣 φ
  obtain ⟨a, I, rfl⟩ : ∃ (a : Fin n) (I : Fin (k + 5) → Fin n), L = Fin.cons a I :=
    ⟨L 0, Fin.tail L, (Fin.cons_self_tail L).symm⟩
  rw [integral_jet_cons_mul hφ (k + 5) a I (hw.of_le (by simp)) hc, abs_neg]
  have := ibb_young (k := k) (ε := ε) (C₁ := (2 * |Λ| + 1) * Cf) (C₂ := C₂) hε hφ
    (hom φ hφ hK)
    (h := fun x => jet (k + 5) φ x I * fderiv ℝ (locWeight ρ 𝔣 φ) x (EuclideanSpace.single a 1))
    (((ibb_coord_contDiff hφ _ I).continuous).mul
      ((hw.continuous_fderiv (by simp)).clm_apply continuous_const))
    (fun x hx => by
      have : x ∉ tsupport ρ := fun h => hx (hρ.supp h)
      simp [fderiv_locWeight_eq_zero 𝔣 φ this])
    (fun x hx => by
      by_cases hxs : x ∈ tsupport ρ
      · rw [abs_mul]
        have := ibb_wfd (by omega) h𝔣 hf hρ hφ (hsupp x hxs) a
        exact mul_le_mul (ibb_coord_le _ x I) (by simpa [mul_assoc] using this) (abs_nonneg _)
          (norm_nonneg _)
      · simp only [fderiv_locWeight_eq_zero 𝔣 φ hxs, zero_apply, mul_zero,
          abs_zero]
        exact mul_nonneg (norm_nonneg _) (mul_nonneg (by positivity)
          (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x)))))
  refine this.trans (le_of_eq ?_)
  ring

/-- Second derivative of the localized weight minus the top term. -/
theorem ibb_w2 {k : ℕ} (hN : k + 5 ≤ N) {K : ℝ} {𝔣 : P n N → ℝ}
    (h𝔣 : IsSymbol (admSet n N K) k k 𝔣) {Cf : ℝ}
    (hf : ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)),
      jp N φ x ∈ admSet n N K →
      |𝔣 (jp N φ x)| ≤ Cf * Om (k + 1) (k + 1) (aJet φ x) ∧
      ‖totalD 𝔣 (jp N φ x)‖ ≤ Cf * Om (k + 1) (k + 1) (aJet φ x))
    {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)}
    (hx : jp N φ x ∈ admSet n N K) (a c : Fin n) :
    |fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y (EuclideanSpace.single c 1)) x
        (EuclideanSpace.single a 1) -
      ρ x ^ 2 * totalD (fun q => totalD 𝔣 q (EuclideanSpace.single c 1)) (jp N φ x)
        (EuclideanSpace.single a 1)| ≤
      (2 * |Λ| ^ 2 + 6 * |Λ|) * Cf * Om (k + 1) (k + 1) (aJet φ x) := by
  obtain ⟨hf1, hf2⟩ := hf φ x hx
  rw [fderiv_fderiv_locWeight h𝔣 (admSet_low K) hN hφ hρ.smooth hx]
  have e1 : ∀ s t u v : ℝ, s + t + u + v - v = s + t + u := fun _ _ _ _ => by ring
  rw [e1]
  have hρ1 := hρ.le_one x
  have ta := ibb_norm_td a _ hf2
  have tc := ibb_norm_td c _ hf2
  have t1 : |fderiv ℝ (fun y => 2 * ρ y * fderiv ℝ ρ y (EuclideanSpace.single c 1)) x
      (EuclideanSpace.single a 1) * 𝔣 (jp N φ x)| ≤
      (2 * |Λ| ^ 2 + 2 * |Λ|) * (Cf * Om (k + 1) (k + 1) (aJet φ x)) := by
    rw [abs_mul]
    exact mul_le_mul (ibb_u hρ x a c) hf1 (abs_nonneg _) (by positivity)
  have t2 := ibb_term hρ1 (ibb_d1 hρ x c) ta
  have t3 := ibb_term hρ1 (ibb_d1 hρ x a) tc
  refine (abs_add_three _ _ _).trans ?_
  linarith

/-- Two integrations by parts, keeping the top term. -/
theorem exists_const_integral_jet2_locWeight {N k : ℕ} (hN : k + 5 ≤ N) {K : ℝ} (Λ : ℝ)
    {𝔣 : P n N → ℝ} (h𝔣 : IsSymbol (admSet n N K) k k 𝔣) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∀ (a c : Fin n) (L : Fin (k + 5) → Fin n),
        |(∫ x, jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c L)) * locWeight ρ 𝔣 φ x)
          - ∫ x, jet (k + 5) φ x L *
              (ρ x ^ 2 * totalD (fun q => totalD 𝔣 q (EuclideanSpace.single c 1)) (jp N φ x)
                (EuclideanSpace.single a 1))| ≤
          ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
          C * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨Cf, hCf, hf⟩ := ibb_f_td (by omega : k + 4 ≤ N) h𝔣
  obtain ⟨C₂, hC₂, hom⟩ := exists_const_integral_Om_sq (n := n) (k + 1) (k + 1) (ρ := 2) (R := 3)
    (by norm_num) (by norm_num) K
  refine ⟨((2 * |Λ| ^ 2 + 6 * |Λ|) * Cf) ^ 2 / (4 * ε) * C₂, by positivity,
    fun ρ hρ φ hφ hadm a c L => ?_⟩
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hlow := admSet_low (n := n) (N := N) K
  have hw := contDiff_locWeight h𝔣 hφ hρ.smooth hsupp
  have hc := hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) 𝔣 φ
  have hG : IsSymbol (admSet n N K) (k + 1) (max k k + 1)
      (fun q => totalD 𝔣 q (EuclideanSpace.single c 1)) :=
    (h𝔣.totalD (by omega) hlow).clm_comp
      (ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single c (1 : ℝ) : EuclideanSpace ℝ (Fin n)))
  have hG2 : IsSymbol (admSet n N K) (k + 1 + 1) (max (max k k + 1) (k + 1) + 1)
      (fun q => totalD (fun q => totalD 𝔣 q (EuclideanSpace.single c 1)) q
        (EuclideanSpace.single a 1)) :=
    (hG.totalD (by omega) hlow).clm_comp
      (ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single a (1 : ℝ) : EuclideanSpace ℝ (Fin n)))
  set G : P n N → ℝ := fun q => totalD (fun q => totalD 𝔣 q (EuclideanSpace.single c 1)) q
    (EuclideanSpace.single a 1) with hGdef
  have hT := contDiff_locWeight hG2 hφ hρ.smooth hsupp
  have hTc := hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) G φ
  have hg : ContDiff ℝ ∞ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y (EuclideanSpace.single c (1 : ℝ))) :=
    (hw.fderiv_right (m := ∞) (by simp)).clm_apply contDiff_const
  have hD2c : Continuous (fun x => fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y
      (EuclideanSpace.single c (1 : ℝ))) x (EuclideanSpace.single a (1 : ℝ))) :=
    (hg.continuous_fderiv (by simp)).clm_apply continuous_const
  have hD2s : HasCompactSupport (fun x => fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y
      (EuclideanSpace.single c (1 : ℝ))) x (EuclideanSpace.single a (1 : ℝ))) :=
    (hc.fderiv_apply ℝ _).fderiv_apply ℝ _
  have i1 : Integrable (fun x => jet (k + 5) φ x L * fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y
      (EuclideanSpace.single c (1 : ℝ))) x (EuclideanSpace.single a (1 : ℝ))) := ((ibb_coord_contDiff hφ (k + 5) L).continuous.mul hD2c).integrable_of_hasCompactSupport
    hD2s.mul_left
  have i2 : Integrable (fun x => jet (k + 5) φ x L * locWeight ρ G φ x) := ((ibb_coord_contDiff hφ (k + 5) L).continuous.mul hT.continuous).integrable_of_hasCompactSupport
    hTc.mul_left
  rw [integral_jet_cons_cons_mul hφ (k + 5) a c L (hw.of_le (by first | exact WithTop.coe_le_coe.2 le_top | norm_cast; exact le_top | simp | exact_mod_cast le_top)) hc]
  change |(∫ x, jet (k + 5) φ x L * fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y
      (EuclideanSpace.single c (1 : ℝ))) x (EuclideanSpace.single a (1 : ℝ))) -
    ∫ x, jet (k + 5) φ x L * locWeight ρ G φ x| ≤ _
  have := ibb_young (k := k) (ε := ε) (C₁ := (2 * |Λ| ^ 2 + 6 * |Λ|) * Cf) (C₂ := C₂) hε hφ
    (hom φ hφ hK)
    (h := fun x => jet (k + 5) φ x L * (fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y
      (EuclideanSpace.single c (1 : ℝ))) x (EuclideanSpace.single a (1 : ℝ)) - locWeight ρ G φ x))
    ((ibb_coord_contDiff hφ (k + 5) L).continuous.mul (hD2c.sub hT.continuous))
    (fun x hx => by
      have hxs : x ∉ tsupport ρ := fun h => hx (hρ.supp h)
      simp [ibb_D2_zero 𝔣 φ hxs, locWeight_eq_zero G φ hxs])
    (fun x hx => by
      by_cases hxs : x ∈ tsupport ρ
      · rw [abs_mul]
        have := ibb_w2 (by omega) h𝔣 hf hρ hφ (hsupp x hxs) a c
        exact mul_le_mul (ibb_coord_le _ x L) (by simpa [mul_assoc, locWeight] using this)
          (abs_nonneg _) (norm_nonneg _)
      · simp only [ibb_D2_zero 𝔣 φ hxs, locWeight_eq_zero G φ hxs, sub_zero, mul_zero, abs_zero]
        exact mul_nonneg (norm_nonneg _) (mul_nonneg (by positivity)
          (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x)))))
  refine le_trans (le_of_eq ?_) (this.trans (le_of_eq ?_))
  · rw [← integral_sub i1 i2]
    congr 2
    funext x
    ring
  · ring

/-- info: 'CWFlow.Extension.exists_const_integral_locWeight_odd' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_integral_locWeight_odd

/-- info: 'CWFlow.Extension.exists_const_integral_jet_locWeight' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_integral_jet_locWeight

/-- info: 'CWFlow.Extension.exists_const_integral_jet2_locWeight' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_integral_jet2_locWeight

end CWFlow.Extension
