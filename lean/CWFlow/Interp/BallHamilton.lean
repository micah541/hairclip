import CWFlow.Interp.BallInterp

/-!
# Hamilton's inequality on balls

* `jet_jet_apply`, `norm_jet_jet`: the jet of a jet is a jet;
* `exists_const_jet_rpow_ball`: `∫_{B_ρ} |D^j f|^{2m/j} ≤ C ballEnergy m R f` for `f` bounded on `B_R`;
* `exists_const_jet_rpow_ball_lt`: for lower exponents, with an arbitrarily small constant in front of
  the energy.
-/

open Metric MeasureTheory
open scoped ContDiff

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace CWFlow.Interp

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-! ### Part 1: the jet of a jet -/

theorem jet_cast (ξ : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) {a b : ℕ}
    (h : a = b) (I : Fin a → Fin n) :
    jet a ξ x I = jet b ξ x (fun k => I (Fin.cast h.symm k)) := by
  subst h
  rfl

theorem jet_jet_apply {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ) (i j : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (I : Fin j → Fin n) (I' : Fin i → Fin n) :
    (jet j (jet i ξ) x I) I' = jet (j + i) ξ x (Fin.append I I') := by
  induction j generalizing x with
  | zero =>
    have hI : I = Fin.elim0 := Subsingleton.elim _ _
    subst hI
    have h1 : jet 0 (jet i ξ) x Fin.elim0 = jet i ξ x := by
      simp [jet_apply]
    rw [h1, jet_cast ξ x (Nat.zero_add i) (Fin.append Fin.elim0 I')]
    rw [Fin.elim0_append]
    rfl
  | succ j ih =>
    obtain ⟨k, I₀, rfl⟩ : ∃ (k : Fin n) (I₀ : Fin j → Fin n), I = Fin.cons k I₀ :=
      ⟨I 0, Fin.tail I, (Fin.cons_self_tail I).symm⟩
    have hη : ContDiff ℝ ∞ (jet i ξ) := contDiff_jet hξ i
    have hdj : Differentiable ℝ (jet j (jet i ξ)) :=
      (contDiff_jet hη j).differentiable (by simp)
    have hg : Differentiable ℝ (fun y => jet j (jet i ξ) y I₀) :=
      ((PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin j → Fin n => PiLp 2 (fun _ : Fin i → Fin n => V))
        I₀).differentiable).comp hdj
    have e1 : jet (j + 1) (jet i ξ) x (Fin.cons k I₀) I' =
        fderiv ℝ (jet j (jet i ξ)) x (EuclideanSpace.single k 1) I₀ I' := by
      rw [fderiv_jet_apply hη]
    rw [e1, fderiv_piLp_apply (hdj x), fderiv_piLp_apply (hg x)]
    have e2 : (fun y => jet j (jet i ξ) y I₀ I') =
        fun y => jet (j + i) ξ y (Fin.append I₀ I') := by
      funext y; exact ih y I₀
    rw [e2, ← fderiv_piLp_apply ((contDiff_jet hξ (j + i)).differentiable (by simp) x),
      fderiv_jet_apply hξ]
    rw [jet_cast ξ x (Nat.add_right_comm j i 1) _]
    rw [Fin.append_cons]
    rfl

theorem norm_jet_jet {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ) (i j : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) : ‖jet j (jet i ξ) x‖ = ‖jet (j + i) ξ x‖ := by
  refine (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp ?_
  rw [PiLp.norm_sq_eq_of_L2, norm_jet_sq, ← (Fin.appendEquiv (α := Fin n) j i).sum_comp,
    Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [PiLp.norm_sq_eq_of_L2]
  refine Finset.sum_congr rfl fun I' _ => ?_
  have := jet_jet_apply hξ i j x I I'
  rw [this]
  rfl

/-! ### Part 2: Hamilton's inequality on balls -/

theorem exists_const_jet_rpow_ball {m j : ℕ} (hj : 1 ≤ j) (hjm : j ≤ m) {ρ R : ℝ} (hρ : 0 < ρ)
    (hρR : ρ < R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖f x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ (2 * (m : ℝ) / j) ≤
        C * ballEnergy m R f := by
  obtain ⟨C₁, hC₁, χ, hχ, hcs, h01, h1, hts, hE⟩ :=
    exists_cutoff_energy (n := n) (V := V) m hρ hρR
  have hj0 : (0 : ℝ) < j := by exact_mod_cast hj
  have hP : 0 < 2 * (m : ℝ) / j := by
    have : (0 : ℝ) < m := by exact_mod_cast hj.trans hjm
    positivity
  set P : ℝ := 2 * (m : ℝ) / j with hPdef
  have hA : 0 ≤ ((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * (max K 0) ^ (P - 2) :=
    mul_nonneg (Real.rpow_nonneg (by positivity) _) (Real.rpow_nonneg (le_max_right _ _) _)
  refine ⟨((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * (max K 0) ^ (P - 2) * C₁,
    mul_nonneg hA hC₁, fun f hf hK => ?_⟩
  have hξs : ContDiff ℝ ∞ (fun y => χ y • f y) := hχ.smul hf
  have hξc : HasCompactSupport (fun y => χ y • f y) :=
    hcs.mono (Function.support_smul_subset_left _ _)
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
  have hcont : Continuous (fun x => ‖jet j (fun y => χ y • f y) x‖ ^ P) :=
    ((contDiff_jet hξs j).continuous.norm).rpow_const fun _ => Or.inr hP.le
  have hcomp : HasCompactSupport (fun x => ‖jet j (fun y => χ y • f y) x‖ ^ P) := by
    have := (hasCompactSupport_jet hξc j).comp_left (g := fun v => ‖v‖ ^ P)
      (by simp [Real.zero_rpow hP.ne'])
    exact this
  have hint := Continuous.integrable_of_hasCompactSupport (μ := volume) hcont hcomp
  have hmain := hamilton_interp hξs hξc hj hjm hM
  rw [← hPdef] at hmain
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ P
      = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j (fun y => χ y • f y) x‖ ^ P := by
        refine setIntegral_congr_fun measurableSet_ball fun x hx => ?_
        rw [jet_congr_of_eventuallyEq (hev x hx).symm j]
    _ ≤ ∫ x, ‖jet j (fun y => χ y • f y) x‖ ^ P :=
        setIntegral_le_integral hint
          (Filter.Eventually.of_forall fun x => Real.rpow_nonneg (norm_nonneg _) _)
    _ ≤ _ := hmain
    _ ≤ _ := by
        calc _ = (((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * (max K 0) ^ (P - 2)) *
              ∫ x, ‖jet m (fun y => χ y • f y) x‖ ^ 2 := by ring
          _ ≤ (((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * (max K 0) ^ (P - 2)) *
              (C₁ * ballEnergy m R f) := mul_le_mul_of_nonneg_left (hE f hf) hA
          _ = _ := by ring

/-! ### Part 3: lower exponents -/

theorem exists_const_jet_rpow_ball_lt {m j : ℕ} (hj : 1 ≤ j) (hjm : j ≤ m) {p : ℝ} (hp0 : 0 < p)
    (hp : p < 2 * (m : ℝ) / j) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ f : EuclideanSpace ℝ (Fin n) → V, ContDiff ℝ ∞ f →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖f x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ p ≤ ε * ballEnergy m R f + C := by
  obtain ⟨C₂, hC₂, h2⟩ := exists_const_jet_rpow_ball (n := n) (V := V) hj hjm hρ hρR K
  set P : ℝ := 2 * (m : ℝ) / j with hPdef
  have hP : 0 < P := hp0.trans hp
  have hθ0 : 0 < p / P := div_pos hp0 hP
  have hθ1 : p / P < 1 := (div_lt_one hP).mpr hp
  set δ : ℝ := ε / (C₂ + 1) with hδdef
  have hδ : 0 < δ := by positivity
  obtain ⟨Cδ, hCδ, hpt⟩ := exists_const_rpow_le_add (A := 1) (θ := p / P) (ε := δ) zero_le_one
    hθ0 hθ1 hδ
  have hδC : δ * C₂ ≤ ε := by
    rw [hδdef, div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith
  refine ⟨Cδ * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) ρ)).toReal,
    mul_nonneg hCδ ENNReal.toReal_nonneg, fun f hf hK => ?_⟩
  have hcp : ∀ q : ℝ, 0 < q → IntegrableOn (fun x => ‖jet j f x‖ ^ q)
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) := fun q hq =>
    ((((contDiff_jet hf j).continuous.norm).rpow_const fun _ => Or.inr hq.le).continuousOn.integrableOn_compact
      (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) ρ)).mono_set ball_subset_closedBall
  have hpt' : ∀ x, ‖jet j f x‖ ^ p ≤ δ * ‖jet j f x‖ ^ P + Cδ := by
    intro x
    have := hpt (‖jet j f x‖ ^ P) (Real.rpow_nonneg (norm_nonneg _) _)
    rw [one_mul, ← Real.rpow_mul (norm_nonneg _), mul_div_cancel₀ _ hP.ne'] at this
    exact this
  have hconst : IntegrableOn (fun _ : EuclideanSpace ℝ (Fin n) => Cδ)
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) := integrableOn_const measure_ball_lt_top.ne
  have hE0 := ballEnergy_nonneg m R f
  have hbound := h2 f hf hK
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ p
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, (δ * ‖jet j f x‖ ^ P + Cδ) :=
        setIntegral_mono (hcp p hp0) (((hcp P hP).const_mul δ).add hconst) hpt'
    _ = δ * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ P) +
          Cδ * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) ρ)).toReal := by
        rw [integral_add ((hcp P hP).const_mul δ) hconst, integral_const_mul, setIntegral_const]
        simp only [Measure.real, smul_eq_mul]
        ring
    _ ≤ ε * ballEnergy m R f + Cδ * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) ρ)).toReal := by
        have key : δ * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet j f x‖ ^ P
            ≤ ε * ballEnergy m R f := by
          calc _ ≤ δ * (C₂ * ballEnergy m R f) := mul_le_mul_of_nonneg_left hbound hδ.le
            _ = (δ * C₂) * ballEnergy m R f := by ring
            _ ≤ ε * ballEnergy m R f := mul_le_mul_of_nonneg_right hδC hE0
        linarith

/-- info: 'CWFlow.Interp.norm_jet_jet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_jet_jet

/-- info: 'CWFlow.Interp.exists_const_jet_rpow_ball' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_jet_rpow_ball

/-- info: 'CWFlow.Interp.exists_const_jet_rpow_ball_lt' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_jet_rpow_ball_lt

end CWFlow.Interp
