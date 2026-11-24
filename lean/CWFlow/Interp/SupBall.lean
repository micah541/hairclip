import CWFlow.Interp.SupBound
import CWFlow.Interp.BallHamilton

/-!
# A local sup bound for jets

For smooth `φ` on `ℝⁿ`, `sup_{B_ρ} |Dᵐ φ|` is bounded by a constant times
`1 + Σ_{i ≤ n} ∫_{B_R} |D^{m+i} φ|²`, by the Sobolev-type inequality for compactly supported
functions applied to a cutoff of `jet m φ`.
-/

open Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Interp

theorem exists_const_sup_jet {n : ℕ} (m : ℕ) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
        ‖jet m φ x‖ ≤ C * (1 + ∑ i ∈ Finset.range (n + 1),
          ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet (m + i) φ y‖ ^ 2) := by
  obtain ⟨C₀, hC₀, χ, hχ, hcs, h01, h1, hts, hE⟩ :=
    exists_cutoff_energy (n := n) (V := PiLp 2 (fun _ : Fin m → Fin n => ℝ)) n hρ hρR
  refine ⟨(volume (ball (0 : EuclideanSpace ℝ (Fin n)) R)).toReal + C₀,
    add_nonneg ENNReal.toReal_nonneg hC₀, fun φ hφ x hx => ?_⟩
  set B := ball (0 : EuclideanSpace ℝ (Fin n)) R with hB
  set f := jet m φ with hf
  have hfs : ContDiff ℝ ∞ f := contDiff_jet hφ m
  set u : EuclideanSpace ℝ (Fin n) → PiLp 2 (fun _ : Fin m → Fin n => ℝ) :=
    fun y => χ y • f y with hu
  have hus : ContDiff ℝ ∞ u := hχ.smul hfs
  have huc : HasCompactSupport u := hcs.smul_right
  set g : EuclideanSpace ℝ (Fin n) → ℝ := fun y => ‖jet n u y‖ with hg
  have hgc : Continuous g := (contDiff_jet hus n).continuous.norm
  have hgcs : HasCompactSupport g := (hasCompactSupport_jet huc n).norm
  have hgi : Integrable g := hgc.integrable_of_hasCompactSupport hgcs
  have hg2c : HasCompactSupport (fun y => g y ^ 2) := hgcs.comp_left (g := fun t : ℝ => t ^ 2) (by simp)
  have hg2i : Integrable (fun y => g y ^ 2) :=
    (hgc.pow 2).integrable_of_hasCompactSupport hg2c
  -- step 1
  have h1x : ‖f x‖ ≤ ∫ y, g y := by
    have e : ‖f x‖ = ‖u x‖ := by simp [hu, h1 x hx]
    rw [e]
    refine (norm_le_integral_norm_iteratedFDeriv hus huc x).trans ?_
    exact integral_mono_of_nonneg (Filter.Eventually.of_forall fun y => norm_nonneg _) hgi
      (Filter.Eventually.of_forall fun y => norm_iteratedFDeriv_le_norm_jet n u y)
  -- support
  have hsupp : ∀ y, y ∉ B → g y = 0 := by
    intro y hy
    have : y ∉ tsupport χ := fun h => hy (hts h)
    have hy0 : χ y = 0 := image_eq_zero_of_notMem_tsupport this
    have hev : u =ᶠ[nhds y] fun _ => 0 := by
      have : (tsupport χ)ᶜ ∈ nhds y := (isClosed_tsupport χ).isOpen_compl.mem_nhds this
      filter_upwards [this] with z hz
      simp [hu, image_eq_zero_of_notMem_tsupport hz]
    have hj : jet n u y = jet n (fun _ => (0 : PiLp 2 (fun _ : Fin m → Fin n => ℝ))) y :=
      jet_congr_of_eventuallyEq hev n
    simp only [hg, hj]
    have : jet n (fun _ : EuclideanSpace ℝ (Fin n) => (0 : PiLp 2 (fun _ : Fin m → Fin n => ℝ))) y = 0 := by
      ext I
      simp [jet_apply]
    simp [this]
  have hint : ∫ y, g y = ∫ y in B, g y :=
    (setIntegral_eq_integral_of_forall_compl_eq_zero hsupp).symm
  have hone : IntegrableOn (fun _ : EuclideanSpace ℝ (Fin n) => (1 : ℝ)) B :=
    integrableOn_const measure_ball_lt_top.ne
  have h2 : ∫ y in B, g y ≤ ∫ y in B, (1 + g y ^ 2) :=
    setIntegral_mono hgi.integrableOn (hone.add hg2i.integrableOn)
      fun y => by nlinarith [sq_nonneg (g y - 1)]
  have h3 : ∫ y in B, (1 + g y ^ 2) = (volume B).toReal + ∫ y in B, g y ^ 2 := by
    rw [integral_add hone hg2i.integrableOn, setIntegral_const]
    simp [measureReal_def]
  have h4 : ∫ y in B, g y ^ 2 ≤ ∫ y, g y ^ 2 :=
    setIntegral_le_integral hg2i (Filter.Eventually.of_forall fun y => by positivity)
  have h5 : ∫ y, g y ^ 2 ≤ C₀ * ballEnergy n R f := hE f hfs
  have hcomm : ∀ i y, ‖jet i f y‖ = ‖jet (m + i) φ y‖ := by
    intro i y; rw [hf, norm_jet_jet hφ, add_comm]
  have hEn : ballEnergy n R f = ∑ i ∈ Finset.range (n + 1),
      ∫ y in B, ‖jet (m + i) φ y‖ ^ 2 := by
    unfold ballEnergy
    simp only [hcomm]
    rfl
  rw [hEn] at h5
  have hS : 0 ≤ ∑ i ∈ Finset.range (n + 1), ∫ y in B, ‖jet (m + i) φ y‖ ^ 2 :=
    Finset.sum_nonneg fun i _ => setIntegral_nonneg measurableSet_ball fun y _ => by positivity
  have hv : 0 ≤ (volume B).toReal := ENNReal.toReal_nonneg
  nlinarith

/-- info: 'CWFlow.Interp.exists_const_sup_jet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_sup_jet

end CWFlow.Interp
