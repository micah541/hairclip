import CWFlow.Interp.Hamilton
import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
import Mathlib.Analysis.Calculus.BumpFunction.Basic
import Mathlib.Analysis.Calculus.ContDiff.Bounds
import Mathlib.Analysis.Normed.Module.Ball.Pointwise

/-!
# Localized interpolation on balls

Versions on balls of `ℝⁿ` of the interpolation inequality `integral_prod_norm_jet_le`, obtained by
multiplying with a cutoff function.

* `jet_congr_of_eventuallyEq`: jets are local;
* `exists_global_of_contDiffOn`: a function smooth near a closed ball agrees on a larger ball with a
  smooth compactly supported function;
* `exists_cutoff_energy`: the energy of `χ • f` is controlled by `ballEnergy m R f`;
* `exists_const_prod_jet_ball`, `exists_const_prod_jet_ball_lt`: localized product estimates.
-/

open Metric MeasureTheory
open scoped ContDiff

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace CWFlow.Interp

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- `Σ_{i ≤ m} ∫_{B_R} |Dⁱ f|²`. -/
noncomputable def ballEnergy (m : ℕ) (R : ℝ) (f : EuclideanSpace ℝ (Fin n) → V) : ℝ :=
  ∑ i ∈ Finset.range (m + 1), ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet i f x‖ ^ 2

/-! ### Part 1: jets are local; globalization -/

theorem jet_congr_of_eventuallyEq {f g : EuclideanSpace ℝ (Fin n) → V}
    {x : EuclideanSpace ℝ (Fin n)} (h : f =ᶠ[nhds x] g) (j : ℕ) : jet j f x = jet j g x := by
  have := (Filter.EventuallyEq.iteratedFDeriv ℝ h j).eq_of_nhds
  unfold jet
  rw [this]

theorem exists_global_of_contDiffOn {U : Set (EuclideanSpace ℝ (Fin n))} (hU : IsOpen U) {R : ℝ}
    (hR : 0 < R) (hRU : closedBall (0 : EuclideanSpace ℝ (Fin n)) R ⊆ U)
    {f : EuclideanSpace ℝ (Fin n) → V} (hf : ContDiffOn ℝ ∞ f U) :
    ∃ R' > R, ∃ g : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ g ∧ HasCompactSupport g ∧
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R', g x = f x := by
  obtain ⟨δ, hδ, hth⟩ := (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) R).exists_cthickening_subset_open
    hU hRU
  rw [cthickening_closedBall hδ.le hR.le] at hth
  have hR' : 0 < R + δ / 3 := by linarith
  have hlt : R + δ / 3 < R + δ / 3 * 2 := by linarith
  let χ : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨R + δ / 3, R + δ / 3 * 2, hR', hlt⟩
  have hsub : closedBall (0 : EuclideanSpace ℝ (Fin n)) (R + δ / 3 * 2) ⊆ U :=
    (closedBall_subset_closedBall (by linarith)).trans hth
  refine ⟨R + δ / 3, by linarith, fun x => χ x • f x, ?_, ?_, ?_⟩
  · rw [contDiff_iff_contDiffAt]
    intro x
    by_cases hx : x ∈ U
    · exact (χ.contDiffAt).smul ((hf.contDiffAt (hU.mem_nhds hx)))
    · have hx' : x ∈ (closedBall (0 : EuclideanSpace ℝ (Fin n)) (R + δ / 3 * 2))ᶜ := fun h => hx (hsub h)
      have hev : (fun y => χ y • f y) =ᶠ[nhds x] fun _ => (0 : V) := by
        filter_upwards [isClosed_closedBall.isOpen_compl.mem_nhds hx'] with y hy
        have : χ y = 0 := by
          apply χ.zero_of_le_dist
          rw [Set.mem_compl_iff, mem_closedBall, dist_zero_right, not_le] at hy
          rw [dist_zero_right]
          exact hy.le
        simp [this]
      exact contDiffAt_const.congr_of_eventuallyEq hev
  · refine HasCompactSupport.intro (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) (R + δ / 3 * 2)) ?_
    intro x hx
    have : χ x = 0 := by
      apply χ.zero_of_le_dist
      rw [mem_closedBall, dist_zero_right, not_le] at hx
      rw [dist_zero_right]
      exact hx.le
    simp [this]
  · intro x hx
    rw [mem_ball, dist_zero_right] at hx
    have : χ x = 1 := χ.one_of_mem_closedBall (by rw [mem_closedBall, dist_zero_right]; exact hx.le)
    simp [this]

/-- Real-variable estimate: `A y^θ ≤ ε y + C` for `0 < θ < 1`. -/
theorem exists_const_rpow_le_add {A θ ε : ℝ} (hA : 0 ≤ A) (hθ0 : 0 < θ) (hθ1 : θ < 1)
    (hε : 0 < ε) : ∃ C : ℝ, 0 ≤ C ∧ ∀ y : ℝ, 0 ≤ y → A * y ^ θ ≤ ε * y + C := by
  have hb : 0 < (A + 1) / ε := by positivity
  set T : ℝ := ((A + 1) / ε) ^ (1 / (1 - θ)) with hT
  have hT0 : 0 < T := Real.rpow_pos_of_pos hb _
  have hTe : T ^ (θ - 1) = ε / (A + 1) := by
    rw [hT, ← Real.rpow_mul hb.le]
    have : 1 / (1 - θ) * (θ - 1) = -1 := by
      have : 1 - θ ≠ 0 := by linarith
      field_simp; ring
    rw [this, Real.rpow_neg_one, inv_div]
  refine ⟨A * T ^ θ, by positivity, fun y hy => ?_⟩
  by_cases hyT : y ≤ T
  · have : y ^ θ ≤ T ^ θ := Real.rpow_le_rpow hy hyT hθ0.le
    nlinarith [mul_nonneg hε.le hy, mul_le_mul_of_nonneg_left this hA]
  · rw [not_le] at hyT
    have hy0 : 0 < y := hT0.trans hyT
    have h1 : y ^ θ = y * y ^ (θ - 1) := by
      rw [Real.rpow_sub hy0, Real.rpow_one]; field_simp
    have h2 : y ^ (θ - 1) ≤ T ^ (θ - 1) := Real.rpow_le_rpow_of_nonpos hT0 hyT.le (by linarith)
    have h3 : A * T ^ (θ - 1) ≤ ε := by
      rw [hTe, mul_div_assoc']
      rw [div_le_iff₀ (by positivity)]
      nlinarith
    have h4 : 0 ≤ A * T ^ θ := by positivity
    calc A * y ^ θ = A * y * y ^ (θ - 1) := by rw [h1]; ring
      _ ≤ A * y * T ^ (θ - 1) := mul_le_mul_of_nonneg_left h2 (by positivity)
      _ = y * (A * T ^ (θ - 1)) := by ring
      _ ≤ y * ε := mul_le_mul_of_nonneg_left h3 hy
      _ ≤ ε * y + A * T ^ θ := by nlinarith

/-! ### Part 2: energy of a cut-off function -/

theorem exists_cutoff_energy (m : ℕ) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) :
    ∃ C : ℝ, 0 ≤ C ∧ ∃ χ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ χ ∧ HasCompactSupport χ ∧
      (∀ x, 0 ≤ χ x ∧ χ x ≤ 1) ∧ (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ, χ x = 1) ∧
      tsupport χ ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) R ∧
      ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
        ∫ x, ‖jet m (fun y => χ y • f y) x‖ ^ 2 ≤ C * ballEnergy m R f := by
  let b : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨ρ, (ρ + R) / 2, hρ, by linarith⟩
  have hχ : ContDiff ℝ ∞ (b : EuclideanSpace ℝ (Fin n) → ℝ) := b.contDiff
  have hcs : HasCompactSupport (b : EuclideanSpace ℝ (Fin n) → ℝ) := b.hasCompactSupport
  have hts : tsupport (b : EuclideanSpace ℝ (Fin n) → ℝ) ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) R := by
    rw [b.tsupport_eq]
    intro x hx
    rw [mem_closedBall, dist_zero_right] at hx
    rw [mem_ball, dist_zero_right]
    have : ‖x‖ ≤ (ρ + R) / 2 := hx
    linarith
  have hbd : ∀ i, ∃ L, ∀ x, ‖iteratedFDeriv ℝ i (b : EuclideanSpace ℝ (Fin n) → ℝ) x‖ ≤ L :=
    fun i => (contDiff_iteratedFDeriv_infty hχ i).continuous.bounded_above_of_compact_support
      (hcs.iteratedFDeriv (𝕜 := ℝ) i)
  choose Lf hLf using hbd
  obtain ⟨L, hL0, hL⟩ : ∃ L : ℝ, 0 ≤ L ∧ ∀ i ≤ m, ∀ x,
      ‖iteratedFDeriv ℝ i (b : EuclideanSpace ℝ (Fin n) → ℝ) x‖ ≤ L := by
    refine ⟨∑ i ∈ Finset.range (m + 1), |Lf i|, Finset.sum_nonneg fun _ _ => abs_nonneg _,
      fun i hi x => ?_⟩
    calc _ ≤ Lf i := hLf i x
      _ ≤ |Lf i| := le_abs_self _
      _ ≤ _ := Finset.single_le_sum (f := fun i => |Lf i|) (fun _ _ => abs_nonneg _)
          (Finset.mem_range.mpr (Nat.lt_succ_of_le hi))
  obtain ⟨Cc, hCc0, hCC⟩ : ∃ Cc : ℝ, 0 ≤ Cc ∧ ∀ a : ℕ → ℝ,
      (∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * a i) ^ 2 ≤
        Cc * ∑ i ∈ Finset.range (m + 1), (a i) ^ 2 :=
    ⟨∑ i ∈ Finset.range (m + 1), ((m.choose i : ℝ)) ^ 2,
      Finset.sum_nonneg fun _ _ => sq_nonneg _, fun a => Finset.sum_mul_sq_le_sq_mul_sq _ _ _⟩
  refine ⟨(n : ℝ) ^ m * L ^ 2 * Cc, mul_nonneg (mul_nonneg (by positivity) (sq_nonneg _)) hCc0,
    (b : EuclideanSpace ℝ (Fin n) → ℝ), hχ, hcs, fun x => ⟨b.nonneg, b.le_one⟩,
    fun x hx => b.one_of_mem_closedBall (ball_subset_closedBall hx), hts, ?_⟩
  intro f hf
  have hξs : ContDiff ℝ ∞ (fun y => b y • f y) := hχ.smul hf
  have hpt : ∀ x, ‖jet m (fun y => b y • f y) x‖ ^ 2 ≤
      ((n : ℝ) ^ m * L ^ 2 * Cc) * ∑ k ∈ Finset.range (m + 1), ‖jet k f x‖ ^ 2 := by
    intro x
    have h1 := norm_iteratedFDeriv_smul_le (𝕜 := ℝ) hχ hf x (n := m) (by exact_mod_cast le_top)
    have h2 : ‖iteratedFDeriv ℝ m (fun y => b y • f y) x‖ ≤
        L * ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * ‖jet (m - i) f x‖ := by
      refine h1.trans ?_
      rw [Finset.mul_sum]
      refine Finset.sum_le_sum fun i hi => ?_
      have hi' : i ≤ m := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
      have e1 := hL i hi' x
      have e2 := norm_iteratedFDeriv_le_norm_jet (m - i) f x
      calc (m.choose i : ℝ) * ‖iteratedFDeriv ℝ i (b : EuclideanSpace ℝ (Fin n) → ℝ) x‖ *
            ‖iteratedFDeriv ℝ (m - i) f x‖
          ≤ (m.choose i : ℝ) * L * ‖jet (m - i) f x‖ := by gcongr
        _ = _ := by ring
    have h3 : ‖jet m (fun y => b y • f y) x‖ ≤
        (Real.sqrt n) ^ m * (L * ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * ‖jet (m - i) f x‖) :=
      (norm_jet_le m _ x).trans (mul_le_mul_of_nonneg_left h2 (by positivity))
    have h4 := hCC fun i => ‖jet (m - i) f x‖
    have h5 : ∑ i ∈ Finset.range (m + 1), ‖jet (m - i) f x‖ ^ 2 =
        ∑ k ∈ Finset.range (m + 1), ‖jet k f x‖ ^ 2 := by
      have := Finset.sum_range_reflect (fun k => ‖jet k f x‖ ^ 2) (m + 1)
      simpa using this
    rw [h5] at h4
    calc ‖jet m (fun y => b y • f y) x‖ ^ 2
        ≤ ((Real.sqrt n) ^ m * (L * ∑ i ∈ Finset.range (m + 1),
            (m.choose i : ℝ) * ‖jet (m - i) f x‖)) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h3 2
      _ = (n : ℝ) ^ m * L ^ 2 * (∑ i ∈ Finset.range (m + 1),
            (m.choose i : ℝ) * ‖jet (m - i) f x‖) ^ 2 := by
          rw [mul_pow, mul_pow, ← pow_mul, mul_comm m 2, pow_mul, Real.sq_sqrt (Nat.cast_nonneg _)]
          ring
      _ ≤ (n : ℝ) ^ m * L ^ 2 * (Cc * ∑ k ∈ Finset.range (m + 1), ‖jet k f x‖ ^ 2) :=
          mul_le_mul_of_nonneg_left h4 (by positivity)
      _ = _ := by ring
  have hjc : ∀ k, Continuous (fun x => ‖jet k f x‖ ^ 2) := fun k =>
    ((contDiff_jet hf k).continuous.norm).pow 2
  have hint : ∀ k, IntegrableOn (fun x => ‖jet k f x‖ ^ 2) (ball (0 : EuclideanSpace ℝ (Fin n)) R) :=
    fun k => ((hjc k).continuousOn.integrableOn_compact
      (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) R)).mono_set ball_subset_closedBall
  have hξint : IntegrableOn (fun x => ‖jet m (fun y => b y • f y) x‖ ^ 2)
      (ball (0 : EuclideanSpace ℝ (Fin n)) R) :=
    ((((contDiff_jet hξs m).continuous.norm).pow 2).continuousOn.integrableOn_compact
      (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) R)).mono_set ball_subset_closedBall
  have hsum : IntegrableOn (fun x => ((n : ℝ) ^ m * L ^ 2 * Cc) *
      ∑ k ∈ Finset.range (m + 1), ‖jet k f x‖ ^ 2) (ball (0 : EuclideanSpace ℝ (Fin n)) R) :=
    (integrable_finsetSum _ fun k _ => hint k).const_mul _
  calc ∫ x, ‖jet m (fun y => b y • f y) x‖ ^ 2
      = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet m (fun y => b y • f y) x‖ ^ 2 := by
        symm
        refine setIntegral_eq_integral_of_forall_compl_eq_zero fun x hx => ?_
        have hxt : x ∉ tsupport (fun y => b y • f y) := fun h =>
          hx (hts (tsupport_smul_subset_left (b : EuclideanSpace ℝ (Fin n) → ℝ) f h))
        rw [jet_eq_zero_of_notMem _ hxt]; simp
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ((n : ℝ) ^ m * L ^ 2 * Cc) *
          ∑ k ∈ Finset.range (m + 1), ‖jet k f x‖ ^ 2 :=
        setIntegral_mono hξint hsum fun x => hpt x
    _ = _ := by
        rw [integral_const_mul, integral_finsetSum _ fun k _ => hint k]
        rfl

/-! ### Part 3: products on balls -/

theorem ballEnergy_nonneg (m : ℕ) (R : ℝ) (f : EuclideanSpace ℝ (Fin n) → V) :
    0 ≤ ballEnergy m R f :=
  Finset.sum_nonneg fun i _ => setIntegral_nonneg measurableSet_ball fun x _ => sq_nonneg _

/-- Common core of the two product estimates: a bound `A (C₁ E)^θ` with `θ = (∑ jᵢ)/(2m)`. -/
theorem exists_prod_jet_ball_rpow {m q : ℕ} (hm : 1 ≤ m) (hq : 1 ≤ q) (j : Fin q → ℕ)
    (hj1 : ∀ i, 1 ≤ j i) (hjm : ∀ i, j i ≤ m) (hsum : ∑ i, j i ≤ 2 * m)
    {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) :
    ∃ A C₁ : ℝ, 0 ≤ A ∧ 0 ≤ C₁ ∧ ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖f x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ∏ i, ‖jet (j i) f x‖ ≤
        A * (C₁ * ballEnergy m R f) ^ ((∑ i, (j i : ℝ)) / (2 * m)) := by
  obtain ⟨C₁, hC₁, χ, hχ, hcs, h01, h1, hts, hE⟩ :=
    exists_cutoff_energy (n := n) (V := V) m hρ hρR
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have hS2 : (∑ i, (j i : ℝ)) ≤ 2 * m := by exact_mod_cast hsum
  have hθ0 : 0 ≤ (∑ i, (j i : ℝ)) / (2 * m) :=
    div_nonneg (Finset.sum_nonneg fun i _ => Nat.cast_nonneg _) (by positivity)
  have hθ1 : (∑ i, (j i : ℝ)) / (2 * m) ≤ 1 := by
    rw [div_le_one (by positivity)]; exact hS2
  refine ⟨((n : ℝ) + 2 * m) ^ ((q : ℝ) * m ^ 2) *
      (max K 0) ^ ((q : ℝ) - (∑ i, (j i : ℝ)) / m) *
      (volume (ball (0 : EuclideanSpace ℝ (Fin n)) R)).toReal ^ (1 - (∑ i, (j i : ℝ)) / (2 * m)),
    C₁, ?_, hC₁, fun f hf hK => ?_⟩
  · have : (0 : ℝ) ≤ (n : ℝ) + 2 * m := by positivity
    exact mul_nonneg (mul_nonneg (Real.rpow_nonneg this _) (Real.rpow_nonneg (le_max_right _ _) _))
      (Real.rpow_nonneg ENNReal.toReal_nonneg _)
  have hξs : ContDiff ℝ ∞ (fun y => χ y • f y) := hχ.smul hf
  have hξc : HasCompactSupport (fun y => χ y • f y) :=
    hcs.mono (Function.support_smul_subset_left _ _)
  have hξt : tsupport (fun y => χ y • f y) ⊆ tsupport χ := tsupport_smul_subset_left χ f
  have hM : ∀ x, ‖χ x • f x‖ ≤ max K 0 := by
    intro x
    by_cases hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R
    · rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (h01 x).1]
      calc χ x * ‖f x‖ ≤ 1 * ‖f x‖ := mul_le_mul_of_nonneg_right (h01 x).2 (norm_nonneg _)
        _ = ‖f x‖ := one_mul _
        _ ≤ K := hK x hx
        _ ≤ _ := le_max_left _ _
    · have : χ x = 0 := image_eq_zero_of_notMem_tsupport fun h => hx (hts h)
      simp [this]
  have hev : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
      (fun y => χ y • f y) =ᶠ[nhds x] f := by
    intro x hx
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy
    simp [h1 y hy]
  have hcont : Continuous (fun x => ∏ i, ‖jet (j i) (fun y => χ y • f y) x‖) :=
    continuous_finsetProd _ fun i _ => (contDiff_jet hξs (j i)).continuous.norm
  have hcomp : HasCompactSupport (fun x => ∏ i, ‖jet (j i) (fun y => χ y • f y) x‖) :=
    HasCompactSupport.intro (hξc : IsCompact (tsupport _)) fun x hx =>
      Finset.prod_eq_zero (Finset.mem_univ (⟨0, hq⟩ : Fin q))
        (by rw [jet_eq_zero_of_notMem _ hx, norm_zero])
  have hint := Continuous.integrable_of_hasCompactSupport (μ := volume) hcont hcomp
  have hQ0 : 0 ≤ ∫ x, ‖jet m (fun y => χ y • f y) x‖ ^ 2 := integral_nonneg fun x => sq_nonneg _
  have hmain := integral_prod_norm_jet_le hξs hξc hm hq j hj1 hjm hsum hM
  have hvol : (volume (tsupport (fun y => χ y • f y))).toReal ≤
      (volume (ball (0 : EuclideanSpace ℝ (Fin n)) R)).toReal :=
    ENNReal.toReal_mono measure_ball_lt_top.ne (measure_mono (hξt.trans hts))
  have hP : 0 ≤ ((n : ℝ) + 2 * m) ^ ((q : ℝ) * m ^ 2) * (max K 0) ^ ((q : ℝ) - (∑ i, (j i : ℝ)) / m) :=
    mul_nonneg (Real.rpow_nonneg (by positivity) _) (Real.rpow_nonneg (le_max_right _ _) _)
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ∏ i, ‖jet (j i) f x‖
      = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
          ∏ i, ‖jet (j i) (fun y => χ y • f y) x‖ := by
        refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
        refine Finset.prod_congr rfl fun i _ => ?_
        rw [jet_congr_of_eventuallyEq (hev x hx).symm (j i)]
    _ ≤ ∫ x, ∏ i, ‖jet (j i) (fun y => χ y • f y) x‖ :=
        setIntegral_le_integral hint
          (Filter.Eventually.of_forall fun x => Finset.prod_nonneg fun i _ => norm_nonneg _)
    _ ≤ _ := hmain
    _ ≤ _ := by
        refine mul_le_mul (mul_le_mul_of_nonneg_left
          (Real.rpow_le_rpow ENNReal.toReal_nonneg hvol (by linarith)) hP)
          (Real.rpow_le_rpow hQ0 (hE f hf) hθ0) (Real.rpow_nonneg hQ0 _) ?_
        exact mul_nonneg hP (Real.rpow_nonneg ENNReal.toReal_nonneg _)

theorem exists_const_prod_jet_ball {m q : ℕ} (hm : 1 ≤ m) (hq : 1 ≤ q) (j : Fin q → ℕ)
    (hj1 : ∀ i, 1 ≤ j i) (hjm : ∀ i, j i ≤ m) (hsum : ∑ i, j i = 2 * m)
    {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖f x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ∏ i, ‖jet (j i) f x‖ ≤ C * ballEnergy m R f := by
  obtain ⟨A, C₁, hA, hC₁, h⟩ := exists_prod_jet_ball_rpow (n := n) (V := V) hm hq j hj1 hjm hsum.le
    hρ hρR K
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have hθ : (∑ i, (j i : ℝ)) / (2 * m) = 1 := by
    have : (∑ i, (j i : ℝ)) = 2 * m := by exact_mod_cast hsum
    rw [this]; field_simp
  refine ⟨A * C₁, mul_nonneg hA hC₁, fun f hf hK => ?_⟩
  have := h f hf hK
  rw [hθ, Real.rpow_one] at this
  linarith

theorem exists_const_prod_jet_ball_lt {m q : ℕ} (hm : 1 ≤ m) (hq : 1 ≤ q) (j : Fin q → ℕ)
    (hj1 : ∀ i, 1 ≤ j i) (hjm : ∀ i, j i ≤ m) (hsum : ∑ i, j i < 2 * m)
    {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖f x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ∏ i, ‖jet (j i) f x‖ ≤
        ε * ballEnergy m R f + C := by
  obtain ⟨A, C₁, hA, hC₁, h⟩ := exists_prod_jet_ball_rpow (n := n) (V := V) hm hq j hj1 hjm hsum.le
    hρ hρR K
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have hS0 : 0 < (∑ i, (j i : ℝ)) :=
    Finset.sum_pos (fun i _ => by exact_mod_cast hj1 i) ⟨⟨0, hq⟩, Finset.mem_univ _⟩
  have hS2 : (∑ i, (j i : ℝ)) < 2 * m := by exact_mod_cast hsum
  have hθ0 : 0 < (∑ i, (j i : ℝ)) / (2 * m) := by positivity
  have hθ1 : (∑ i, (j i : ℝ)) / (2 * m) < 1 := by
    rw [div_lt_one (by positivity)]; exact hS2
  obtain ⟨C, hC0, hC⟩ := exists_const_rpow_le_add
    (mul_nonneg hA (Real.rpow_nonneg hC₁ ((∑ i, (j i : ℝ)) / (2 * m)))) hθ0 hθ1 hε
  refine ⟨C, hC0, fun f hf hK => ?_⟩
  have h1 := h f hf hK
  have h2 := hC _ (ballEnergy_nonneg m R f)
  rw [Real.mul_rpow hC₁ (ballEnergy_nonneg m R f), ← mul_assoc] at h1
  linarith

/-- info: 'CWFlow.Interp.exists_cutoff_energy' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_cutoff_energy

/-- info: 'CWFlow.Interp.exists_const_prod_jet_ball' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_prod_jet_ball

/-- info: 'CWFlow.Interp.exists_const_prod_jet_ball_lt' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_prod_jet_ball_lt

/-- info: 'CWFlow.Interp.exists_global_of_contDiffOn' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_global_of_contDiffOn

end CWFlow.Interp
