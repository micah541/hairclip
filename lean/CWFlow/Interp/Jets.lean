import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.Calculus.ContDiff.WithLp
import Mathlib.Analysis.Calculus.ContDiff.FTaylorSeries
import Mathlib.Analysis.Calculus.FDeriv.CompCLM
import Mathlib.Analysis.Calculus.ContDiff.Comp
import Mathlib.Analysis.Calculus.ContDiff.FiniteDimension
import Mathlib.Topology.Algebra.Module.Multilinear.Topology

/-!
# Coordinate jets as Hilbert-space valued functions

For `ξ : ℝⁿ → V` (`V` a finite dimensional real inner product space) the `j`-jet
`jet j ξ x ∈ ℓ²(Fin j → Fin n; V)` collects all coordinate partial derivatives
`∂_{I 0} ⋯ ∂_{I (j-1)} ξ (x)`.  We prove:

* the recursion `∂_i (jet j ξ)_I = (jet (j+1) ξ)_{i :: I}` and its consequences for sums of squares
  (`sum_norm_fderiv_jet_sq`, `sum_norm_fderiv_fderiv_jet_sq`);
* smoothness and compact support of `jet j ξ`;
* comparison with the operator norm of `iteratedFDeriv`:
  `‖iteratedFDeriv ℝ j ξ x‖ ≤ ‖jet j ξ x‖ ≤ n ^ (j/2) ‖iteratedFDeriv ℝ j ξ x‖`.
-/

open scoped ContDiff

set_option linter.unusedSectionVars false

namespace CWFlow.Interp

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

local notation "E" => EuclideanSpace ℝ (Fin n)

/-- The `j`-jet in coordinates: the family of all iterated partial derivatives `∂_{I 0} ⋯ ∂_{I (j-1)} ξ`,
`I : Fin j → Fin n`, as an element of the Hilbert space `ℓ²(Fin j → Fin n; V)`. -/
noncomputable def jet (j : ℕ) (ξ : E → V) (x : E) : PiLp 2 (fun _ : Fin j → Fin n => V) :=
  WithLp.toLp 2 (fun I => iteratedFDeriv ℝ j ξ x (fun k => EuclideanSpace.single (I k) (1 : ℝ)))

theorem jet_apply (j : ℕ) (ξ : E → V) (x : E) (I : Fin j → Fin n) :
    jet j ξ x I = iteratedFDeriv ℝ j ξ x (fun k => EuclideanSpace.single (I k) 1) := rfl

theorem norm_jet_sq (j : ℕ) (ξ : E → V) (x : E) :
    ‖jet j ξ x‖ ^ 2 =
      ∑ I : Fin j → Fin n, ‖iteratedFDeriv ℝ j ξ x (fun k => EuclideanSpace.single (I k) 1)‖ ^ 2 := by
  rw [PiLp.norm_sq_eq_of_L2]
  rfl

theorem norm_jet_zero (ξ : E → V) (x : E) : ‖jet 0 ξ x‖ = ‖ξ x‖ := by
  have h : ‖jet 0 ξ x‖ ^ 2 = ‖ξ x‖ ^ 2 := by
    rw [norm_jet_sq]
    simp
  exact (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp h

theorem infty_add_natCast_le (j : ℕ) : (∞ : WithTop ℕ∞) + (j : WithTop ℕ∞) ≤ ∞ := by
  induction j with
  | zero => simp
  | succ j ih =>
    rw [Nat.cast_succ, ← add_assoc]
    calc (∞ : WithTop ℕ∞) + (j : WithTop ℕ∞) + 1 ≤ ∞ + 1 := add_le_add_left ih 1
      _ = ∞ := rfl

theorem contDiff_iteratedFDeriv_infty {ξ : E → V} (hξ : ContDiff ℝ ∞ ξ) (j : ℕ) :
    ContDiff ℝ ∞ (iteratedFDeriv ℝ j ξ) :=
  hξ.iteratedFDeriv_right (infty_add_natCast_le j)

theorem contDiff_jet {ξ : E → V} (hξ : ContDiff ℝ ∞ ξ) (j : ℕ) : ContDiff ℝ ∞ (jet j ξ) := by
  refine contDiff_piLp' (𝕜 := ℝ) 2 (fun I => ?_)
  exact (ContinuousMultilinearMap.apply ℝ (fun _ : Fin j => E) V
    (fun k => EuclideanSpace.single (I k) (1 : ℝ))).contDiff.comp (contDiff_iteratedFDeriv_infty hξ j)

theorem hasCompactSupport_jet {ξ : E → V} (hc : HasCompactSupport ξ) (j : ℕ) :
    HasCompactSupport (jet j ξ) := by
  refine (hc.iteratedFDeriv (𝕜 := ℝ) j).mono ?_
  intro x hx
  rw [Function.mem_support] at hx ⊢
  intro h0
  apply hx
  ext I
  simp [jet_apply, h0]

/-- Partial derivatives of a `PiLp`-valued map are computed coordinatewise. -/
theorem fderiv_piLp_apply {ι : Type*} [Fintype ι] {f : E → PiLp 2 (fun _ : ι => V)} {x : E}
    (hf : DifferentiableAt ℝ f x) (v : E) (I : ι) :
    fderiv ℝ f x v I = fderiv ℝ (fun y => f y I) x v := by
  have h := (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => V) I).hasFDerivAt.comp x hf.hasFDerivAt
  have h2 : fderiv ℝ (fun y => f y I) x = (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => V) I).comp (fderiv ℝ f x) :=
    h.fderiv
  rw [h2]
  rfl

theorem fderiv_jet_apply {ξ : E → V} (hξ : ContDiff ℝ ∞ ξ) (j : ℕ) (x : E) (i : Fin n)
    (I : Fin j → Fin n) :
    (fderiv ℝ (jet j ξ) x (EuclideanSpace.single i 1)) I = jet (j + 1) ξ x (Fin.cons i I) := by
  have hd : DifferentiableAt ℝ (jet j ξ) x := (contDiff_jet hξ j).differentiable (by simp) x
  have hd' : DifferentiableAt ℝ (iteratedFDeriv ℝ j ξ) x :=
    (contDiff_iteratedFDeriv_infty hξ j).differentiable (by simp) x
  rw [fderiv_piLp_apply hd]
  show fderiv ℝ (fun y => iteratedFDeriv ℝ j ξ y (fun k => EuclideanSpace.single (I k) (1 : ℝ))) x
    (EuclideanSpace.single i 1) = _
  rw [fderiv_continuousMultilinear_apply_const hd', jet_apply, iteratedFDeriv_succ_apply_left]
  rfl

/-- Reindexing a sum over `Fin (j+1) → Fin n` by the first coordinate. -/
theorem sum_fin_cons (j : ℕ) (g : (Fin (j + 1) → Fin n) → ℝ) :
    ∑ J, g J = ∑ i : Fin n, ∑ I : Fin j → Fin n, g (Fin.cons i I) := by
  rw [← (Fin.consEquiv (fun _ : Fin (j + 1) => Fin n)).sum_comp, Fintype.sum_prod_type]
  rfl

theorem sum_norm_fderiv_jet_sq {ξ : E → V} (hξ : ContDiff ℝ ∞ ξ) (j : ℕ) (x : E) :
    ∑ i, ‖fderiv ℝ (jet j ξ) x (EuclideanSpace.single i 1)‖ ^ 2 = ‖jet (j + 1) ξ x‖ ^ 2 := by
  rw [norm_jet_sq, sum_fin_cons]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [PiLp.norm_sq_eq_of_L2]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [fderiv_jet_apply hξ]
  rfl

theorem sum_norm_fderiv_fderiv_jet_sq {ξ : E → V} (hξ : ContDiff ℝ ∞ ξ) (j : ℕ) (x : E) :
    ∑ i, ∑ k, ‖fderiv ℝ (fun y => fderiv ℝ (jet j ξ) y (EuclideanSpace.single k 1)) x
      (EuclideanSpace.single i 1)‖ ^ 2 = ‖jet (j + 2) ξ x‖ ^ 2 := by
  have hd : ∀ k : Fin n, Differentiable ℝ (fun y => fderiv ℝ (jet j ξ) y (EuclideanSpace.single k 1)) := by
    intro k
    have : ContDiff ℝ ∞ (fderiv ℝ (jet j ξ)) :=
      (contDiff_jet hξ j).fderiv_right (by simp)
    exact (this.clm_apply contDiff_const).differentiable (by simp)
  have hcoord : ∀ (k : Fin n) (i : Fin n) (I : Fin j → Fin n),
      fderiv ℝ (fun y => fderiv ℝ (jet j ξ) y (EuclideanSpace.single k 1)) x
        (EuclideanSpace.single i 1) I = jet (j + 2) ξ x (Fin.cons i (Fin.cons k I)) := by
    intro k i I
    rw [fderiv_piLp_apply (hd k x)]
    have : (fun y => fderiv ℝ (jet j ξ) y (EuclideanSpace.single k 1) I) =
        fun y => jet (j + 1) ξ y (Fin.cons k I) := by
      funext y; exact fderiv_jet_apply hξ j y k I
    rw [this, ← fderiv_piLp_apply ((contDiff_jet hξ (j + 1)).differentiable (by simp) x)]
    exact fderiv_jet_apply hξ (j + 1) x i (Fin.cons k I)
  rw [norm_jet_sq, sum_fin_cons]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [sum_fin_cons j]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [PiLp.norm_sq_eq_of_L2]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [hcoord k i I]
  rfl

omit [FiniteDimensional ℝ V] in
/-- Cauchy–Schwarz expansion of a multilinear map in the coordinate basis. -/
theorem norm_multilinear_sq_le {j : ℕ} (M : ContinuousMultilinearMap ℝ (fun _ : Fin j => E) V)
    (m : Fin j → E) :
    ‖M m‖ ^ 2 ≤ (∑ I : Fin j → Fin n, ‖M (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖ ^ 2) *
      ∏ k, ‖m k‖ ^ 2 := by
  have hm : m = fun k => ∑ i, m k i • EuclideanSpace.single i (1 : ℝ) := by
    funext k
    ext l
    simp [Pi.single_apply]
  have hexp : M m = ∑ I : Fin j → Fin n,
      (∏ k, m k (I k)) • M (fun k => EuclideanSpace.single (I k) (1 : ℝ)) := by
    conv_lhs => rw [hm]
    rw [ContinuousMultilinearMap.map_sum]
    refine Finset.sum_congr rfl fun I _ => ?_
    rw [ContinuousMultilinearMap.map_smul_univ]
  have hp : ∑ I : Fin j → Fin n, ∏ k, (m k (I k)) ^ 2 = ∏ k, ‖m k‖ ^ 2 := by
    have := Finset.prod_univ_sum (fun _ : Fin j => (Finset.univ : Finset (Fin n)))
      (fun k i => (m k i) ^ 2)
    rw [Fintype.piFinset_univ] at this
    rw [← this]
    refine Finset.prod_congr rfl fun k _ => ?_
    rw [PiLp.norm_sq_eq_of_L2]
    simp
  rw [hexp]
  calc ‖∑ I : Fin j → Fin n, (∏ k, m k (I k)) • M (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖ ^ 2
      ≤ (∑ I : Fin j → Fin n, |∏ k, m k (I k)| *
          ‖M (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖) ^ 2 := by
        refine pow_le_pow_left₀ (norm_nonneg _) ((norm_sum_le _ _).trans ?_) 2
        refine Finset.sum_le_sum fun I _ => ?_
        rw [norm_smul, Real.norm_eq_abs]
    _ ≤ (∑ I : Fin j → Fin n, |∏ k, m k (I k)| ^ 2) *
          ∑ I : Fin j → Fin n, ‖M (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖ ^ 2 :=
        Finset.sum_mul_sq_le_sq_mul_sq _ _ _
    _ = _ := by
        rw [mul_comm]
        congr 1
        rw [← hp]
        refine Finset.sum_congr rfl fun I _ => ?_
        rw [sq_abs, Finset.prod_pow]

theorem norm_iteratedFDeriv_le_norm_jet (j : ℕ) (ξ : E → V) (x : E) :
    ‖iteratedFDeriv ℝ j ξ x‖ ≤ ‖jet j ξ x‖ := by
  refine ContinuousMultilinearMap.opNorm_le_bound (norm_nonneg _) fun m => ?_
  have h := norm_multilinear_sq_le (iteratedFDeriv ℝ j ξ x) m
  rw [← norm_jet_sq] at h
  refine (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).mp ?_
  rw [mul_pow, ← Finset.prod_pow]
  exact h

theorem norm_jet_le (j : ℕ) (ξ : E → V) (x : E) :
    ‖jet j ξ x‖ ≤ (Real.sqrt n) ^ j * ‖iteratedFDeriv ℝ j ξ x‖ := by
  have hb : ∀ I : Fin j → Fin n,
      ‖iteratedFDeriv ℝ j ξ x (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖ ≤
        ‖iteratedFDeriv ℝ j ξ x‖ := by
    intro I
    have := (iteratedFDeriv ℝ j ξ x).le_opNorm (fun k => EuclideanSpace.single (I k) (1 : ℝ))
    simpa using this
  have h2 : ‖jet j ξ x‖ ^ 2 ≤ ((Real.sqrt n) ^ j * ‖iteratedFDeriv ℝ j ξ x‖) ^ 2 := by
    rw [norm_jet_sq, mul_pow, ← pow_mul, mul_comm j 2, pow_mul, Real.sq_sqrt (Nat.cast_nonneg _)]
    calc ∑ I : Fin j → Fin n,
          ‖iteratedFDeriv ℝ j ξ x (fun k => EuclideanSpace.single (I k) (1 : ℝ))‖ ^ 2
        ≤ ∑ _I : Fin j → Fin n, ‖iteratedFDeriv ℝ j ξ x‖ ^ 2 :=
          Finset.sum_le_sum fun I _ => pow_le_pow_left₀ (norm_nonneg _) (hb I) 2
      _ = (n : ℝ) ^ j * ‖iteratedFDeriv ℝ j ξ x‖ ^ 2 := by simp
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).mp h2

/-- info: 'CWFlow.Interp.sum_norm_fderiv_jet_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sum_norm_fderiv_jet_sq

/-- info: 'CWFlow.Interp.sum_norm_fderiv_fderiv_jet_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sum_norm_fderiv_fderiv_jet_sq

/-- info: 'CWFlow.Interp.norm_iteratedFDeriv_le_norm_jet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_iteratedFDeriv_le_norm_jet

/-- info: 'CWFlow.Interp.norm_jet_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_jet_le

end CWFlow.Interp
