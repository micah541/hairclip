import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.ContDiff.Comp
import Mathlib.Analysis.Calculus.FDeriv.CompCLM
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.Normed.Module.Multilinear.Curry
import Mathlib.Algebra.BigOperators.Pi

/-!
# The exact Euclidean jet recursion

Let `G i : E → F i` be smooth on an open set `Ω`, related by the first-order recursion

  `D(G i) y = P i (G (i+1) y) + Q i y (G i y)`      (`y ∈ Ω`),

with `P i : F (i+1) →L (E →L F i)` fixed linear maps and `Q i : E → F i →L (E →L F i)` smooth
fields of linear maps. (In a chart, `G i` is the chart reading of the covariant derivative `∇ⁱu`,
`P i` is the insertion into the first slot and `Q i` is the Christoffel term; see
`ShortTime.fderiv_chartRd_eq` in `../MM`, where the recursion holds on a closed ball, hence on its
interior.) `ShortTime.JetRec.exists_iter_le_rec` only bounds `D^n (G 0)`; here we compute it exactly:

* `jetLam P Q n y : Vals F n →L[ℝ] JetT E F n` is a continuous linear map of the values
  `v = (v 0, …, v n)` into the `n`-multilinear maps `JetT E F n = E [×n]→L[ℝ] F 0`, smooth in `y`
  (`contDiffOn_jetLam`), jointly smooth in `(y, v)` (`contDiffOn_jetLam_uncurry`), defined by the
  recursion `jetLam_succ_apply` (product rule);
* `iteratedFDeriv_eq_jetLam`: `D^n (G 0) y = jetLam P Q n y (fun k => G k y)` for `y ∈ Ω`;
* `jetTop P n : F n →L[ℝ] JetT E F n` is the iterated `P`-insertion
  `jetTop P (n+1) t m = jetTop P n (P n t (m 0)) (Fin.tail m)`, and `jetLam P Q n y` depends on the
  top value `v n` exactly through `jetTop P n` (`jetLam_snoc`), so that
  `D^n (G 0) y = jetTop P n (G n y) + (terms in G 0 y, …, G (n-1) y)` (`iteratedFDeriv_eq_jetTop_add`).

With `F i := TensorRSModel 0 i ℝ E` and `P i := slotIns E i`, `jetTop P n t m β = t β m`, i.e. the
`n`-th coordinate derivative is the uncurried `∇ⁿu` plus lower-order terms.

## Implementation note

`JetT E F n` is a (non-reducible) wrapper around `E [×n]→L[ℝ] F 0`, definitionally equal to it, with
the normed-space and `FunLike` structures transported. It is needed because instance resolution fails
on doubly nested spaces `X →L[ℝ] (Y →L[ℝ] E [×n]→L[ℝ] F 0)` written with the raw multilinear type
(the topological instances on `ContinuousMultilinearMap` found by elaboration do not unify with the
normed ones inside nested continuous linear maps), while `X →L[ℝ] (Y →L[ℝ] JetT E F n)` is fine.
Equalities between an `iteratedFDeriv` and a `JetT`-valued term typecheck by definitional unfolding.
-/

noncomputable section

open Set Filter
open scoped ContDiff Topology

namespace CWFlow.JetRec

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {F : ℕ → Type*} [∀ i, NormedAddCommGroup (F i)] [∀ i, NormedSpace ℝ (F i)]

/-- The values `(v 0, …, v n)` of a jet of order `n`. -/
abbrev Vals (F : ℕ → Type*) (n : ℕ) := ∀ k : Fin (n + 1), F k

variable (E F) in
/-- The target of the `n`-th coordinate derivative, `E [×n]→L[ℝ] F 0`, as a wrapper type (see the
implementation note in the module docstring). -/
def JetT (n : ℕ) : Type _ := E [×n]→L[ℝ] F 0

instance (n : ℕ) : NormedAddCommGroup (JetT E F n) :=
  inferInstanceAs (NormedAddCommGroup (E [×n]→L[ℝ] F 0))

instance (n : ℕ) : NormedSpace ℝ (JetT E F n) := inferInstanceAs (NormedSpace ℝ (E [×n]→L[ℝ] F 0))

instance (n : ℕ) : FunLike (JetT E F n) (Fin n → E) (F 0) :=
  inferInstanceAs (FunLike (E [×n]→L[ℝ] F 0) (Fin n → E) (F 0))

theorem JetT.add_apply (n : ℕ) (x y : JetT E F n) (m : Fin n → E) : (x + y) m = x m + y m := rfl

theorem JetT.ext (n : ℕ) {x y : JetT E F n} (h : ∀ m, x m = y m) : x = y :=
  ContinuousMultilinearMap.ext h

variable (E F) in
/-- Uncurrying the first slot, as a continuous linear map. -/
def uncurryL (n : ℕ) : (E →L[ℝ] JetT E F n) →L[ℝ] JetT E F (n + 1) :=
  (continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin (n + 1) => E)
    (F 0)).symm.toContinuousLinearEquiv.toContinuousLinearMap

theorem uncurryL_apply (n : ℕ) (L : E →L[ℝ] JetT E F n) (m : Fin (n + 1) → E) :
    uncurryL E F n L m = L (m 0) (Fin.tail m) := rfl

variable (F) in
/-- Forgetting the top value. -/
def restr (n : ℕ) : Vals F (n + 1) →L[ℝ] Vals F n :=
  ContinuousLinearMap.pi fun k : Fin (n + 1) =>
    ContinuousLinearMap.proj (R := ℝ) (φ := fun k : Fin (n + 2) => F k) k.castSucc

theorem restr_apply (n : ℕ) (v : Vals F (n + 1)) (k : Fin (n + 1)) :
    restr F n v k = v k.castSucc := rfl

variable (P : ∀ i, F (i + 1) →L[ℝ] (E →L[ℝ] F i)) (Q : ∀ i, E → F i →L[ℝ] (E →L[ℝ] F i))

/-- The derivative of the values `(G 0, …, G n)` at `y`, as a linear function of the values
`(G 0 y, …, G (n+1) y)`: `dV P Q n y v w = fun k => P k (v (k+1)) w + Q k y (v k) w`. -/
def dV (n : ℕ) : E → (Vals F (n + 1) →L[ℝ] (E →L[ℝ] Vals F n)) := fun y =>
  ∑ k : Fin (n + 1),
    (ContinuousLinearMap.compL ℝ E (F k) (Vals F n)
      (ContinuousLinearMap.single ℝ (fun k : Fin (n + 1) => F k) k)).comp
      ((P k).comp (ContinuousLinearMap.proj (R := ℝ) (φ := fun k : Fin (n + 2) => F k) k.succ) +
        (Q k y).comp (ContinuousLinearMap.proj (R := ℝ) (φ := fun k : Fin (n + 2) => F k) k.castSucc))

theorem dV_apply (n : ℕ) (y : E) (v : Vals F (n + 1)) (w : E) :
    dV P Q n y v w = fun k : Fin (n + 1) => P k (v k.succ) w + Q k y (v k.castSucc) w := by
  have hs : ∀ (k : Fin (n + 1)) (x : F k),
      ContinuousLinearMap.single ℝ (fun k : Fin (n + 1) => F k) k x = Pi.single k x := fun _ _ => rfl
  simp only [dV, sum_apply, ContinuousLinearMap.compL_apply, ContinuousLinearMap.comp_apply,
    add_apply, ContinuousLinearMap.proj_apply, hs]
  exact Finset.univ_sum_single _

theorem contDiffOn_dV {Ω : Set E} (hQ : ∀ i, ContDiffOn ℝ ∞ (Q i) Ω) (n : ℕ) :
    ContDiffOn ℝ ∞ (dV P Q n) Ω := by
  unfold dV
  refine ContDiffOn.sum fun k _ => ?_
  exact contDiffOn_const.clm_comp (contDiffOn_const.add ((hQ k).clm_comp contDiffOn_const))

/-- **The jet transfer map.** `jetLam P Q n y v` is the `n`-th coordinate derivative at `y` of a
function whose jet of values at `y` is `v = (v 0, …, v n)`; see `iteratedFDeriv_eq_jetLam`. It is
defined by the product rule: `D^(n+1) G 0 = D (y ↦ jetLam n y (V y))`, `V y = (G 0 y, …, G n y)`. -/
def jetLam : ∀ n : ℕ, E → (Vals F n →L[ℝ] JetT E F n)
  | 0 => fun _ =>
    ((continuousMultilinearCurryFin0 ℝ E (F 0)).symm.toContinuousLinearEquiv.toContinuousLinearMap :
      F 0 →L[ℝ] JetT E F 0).comp
      (ContinuousLinearMap.proj (R := ℝ) (φ := fun k : Fin 1 => F k) (Fin.last 0))
  | n + 1 => fun y =>
    (uncurryL E F n).comp
      ((ContinuousLinearMap.compL ℝ E (Vals F n) (JetT E F n) (jetLam n y)).comp (dV P Q n y) +
        (ContinuousLinearMap.flipₗᵢ ℝ E (Vals F n) (JetT E F n) (fderiv ℝ (jetLam n) y)).comp
          (restr F n))

theorem jetLam_zero_apply (y : E) (v : Vals F 0) :
    jetLam P Q 0 y v = (continuousMultilinearCurryFin0 ℝ E (F 0)).symm (v (Fin.last 0)) := rfl

theorem jetLam_succ_apply (n : ℕ) (y : E) (v : Vals F (n + 1)) :
    jetLam P Q (n + 1) y v =
      uncurryL E F n ((jetLam P Q n y).comp (dV P Q n y v) +
        (fderiv ℝ (jetLam P Q n) y).flip (restr F n v)) := rfl

/-- Linearity in the values (`jetLam P Q n y` is a continuous linear map). -/
theorem jetLam_add (n : ℕ) (y : E) (v w : Vals F n) :
    jetLam P Q n y (v + w) = jetLam P Q n y v + jetLam P Q n y w := map_add _ _ _

theorem jetLam_smul (n : ℕ) (y : E) (c : ℝ) (v : Vals F n) :
    jetLam P Q n y (c • v) = c • jetLam P Q n y v := map_smul _ _ _

/-- Smoothness of the jet transfer map in `y`. -/
theorem contDiffOn_jetLam {Ω : Set E} (hΩ : IsOpen Ω) (hQ : ∀ i, ContDiffOn ℝ ∞ (Q i) Ω) :
    ∀ n, ContDiffOn ℝ ∞ (jetLam P Q n) Ω
  | 0 => contDiffOn_const
  | n + 1 => by
    have ih := contDiffOn_jetLam hΩ hQ n
    have hd : ContDiffOn ℝ ∞ (fderiv ℝ (jetLam P Q n)) Ω :=
      ((contDiffOn_infty_iff_fderiv_of_isOpen hΩ).1 ih).2
    have h1 : ContDiff ℝ ∞ (ContinuousLinearMap.compL ℝ E (Vals F n) (JetT E F n)) :=
      ContinuousLinearMap.contDiff _
    have h3 : ContDiff ℝ ∞ (ContinuousLinearMap.flipₗᵢ ℝ E (Vals F n) (JetT E F n)) :=
      LinearIsometryEquiv.contDiff _
    exact contDiffOn_const.clm_comp
      (((h1.comp_contDiffOn ih).clm_comp (contDiffOn_dV P Q hQ n)).add
        ((h3.comp_contDiffOn hd).clm_comp contDiffOn_const))

/-- Joint smoothness of the jet transfer map in `(y, v)`. -/
theorem contDiffOn_jetLam_uncurry {Ω : Set E} (hΩ : IsOpen Ω) (hQ : ∀ i, ContDiffOn ℝ ∞ (Q i) Ω)
    (n : ℕ) :
    ContDiffOn ℝ ∞ (fun p : E × Vals F n => jetLam P Q n p.1 p.2) (Ω ×ˢ univ) :=
  ((contDiffOn_jetLam P Q hΩ hQ n).comp contDiff_fst.contDiffOn fun _ hp => hp.1).clm_apply
    contDiff_snd.contDiffOn

/-- **The exact jet recursion**: under `D(G i) = P i (G (i+1)) + Q i (G i)` on the open set `Ω`,
the `n`-th coordinate derivative of `G 0` at `y ∈ Ω` is `jetLam P Q n y` applied to the values
`(G 0 y, …, G n y)`. -/
theorem iteratedFDeriv_eq_jetLam {Ω : Set E} (hΩ : IsOpen Ω) (hQ : ∀ i, ContDiffOn ℝ ∞ (Q i) Ω)
    {G : ∀ i, E → F i} (hG : ∀ i, ContDiffOn ℝ ∞ (G i) Ω)
    (hrec : ∀ i, ∀ y ∈ Ω, fderiv ℝ (G i) y = P i (G (i + 1) y) + Q i y (G i y)) (n : ℕ) :
    ∀ y ∈ Ω, iteratedFDeriv ℝ n (G 0) y = jetLam P Q n y (fun k => G k y) := by
  induction n with
  | zero =>
    intro y _
    rw [iteratedFDeriv_zero_eq_comp, jetLam_zero_apply]
    rfl
  | succ n ih =>
    intro y hy
    have hnhds : Ω ∈ 𝓝 y := hΩ.mem_nhds hy
    have hΛ : DifferentiableAt ℝ (jetLam P Q n) y :=
      ((contDiffOn_jetLam P Q hΩ hQ n).differentiableOn (by simp)).differentiableAt hnhds
    have hGk : ∀ k : Fin (n + 1), DifferentiableAt ℝ (G k) y := fun k =>
      ((hG k).differentiableOn (by simp)).differentiableAt hnhds
    have hV : DifferentiableAt ℝ (fun z => fun k : Fin (n + 1) => G k z) y :=
      differentiableAt_pi.2 hGk
    -- the induction hypothesis, in the `JetT` world
    have heq : Filter.EventuallyEq (β := JetT E F n) (𝓝 y) (iteratedFDeriv ℝ n (G 0))
        (fun z => jetLam P Q n z (fun k => G k z)) :=
      Filter.eventually_of_mem hnhds fun z hz => ih z hz
    have h1 : (ContinuousLinearMap.pi fun k : Fin (n + 1) => fderiv ℝ (G k) y) =
        dV P Q n y (fun k => G k y) := by
      refine ContinuousLinearMap.ext fun w => ?_
      rw [dV_apply]
      funext k
      rw [ContinuousLinearMap.pi_apply, hrec k y hy]
      rfl
    have h2 : (fun k : Fin (n + 1) => G k y) = restr F n (fun k : Fin (n + 2) => G k y) := rfl
    -- the derivative of `D^n (G 0)`, computed in the `JetT` world by the product rule
    have key := heq.fderiv_eq (𝕜 := ℝ)
    rw [fderiv_clm_apply hΛ hV, fderiv_pi hGk, h1, h2] at key
    rw [iteratedFDeriv_succ_eq_comp_left, Function.comp_apply]
    exact congrArg (uncurryL E F n) key

/-! ### The top-order term -/

/-- The iterated insertion `jetTop P n t m = P 0 (P 1 (⋯ (P (n-1) t (m 0)) ⋯) (m (n-2))) (m (n-1))`:
the dependence of `jetLam P Q n y` on the top value. -/
def jetTop : ∀ n : ℕ, F n →L[ℝ] JetT E F n
  | 0 => (continuousMultilinearCurryFin0 ℝ E (F 0)).symm.toContinuousLinearEquiv.toContinuousLinearMap
  | n + 1 =>
    (uncurryL E F n).comp
      ((ContinuousLinearMap.compL ℝ E (F n) (JetT E F n) (jetTop n)).comp (P n))

theorem jetTop_zero_apply (t : F 0) (m : Fin 0 → E) : jetTop P 0 t m = t := rfl

theorem jetTop_succ_apply (n : ℕ) (t : F (n + 1)) (m : Fin (n + 1) → E) :
    jetTop P (n + 1) t m = jetTop P n (P n t (m 0)) (Fin.tail m) := rfl

theorem jetTop_one_apply (t : F 1) (m : Fin 1 → E) : jetTop P 1 t m = P 0 t (m 0) := rfl

theorem jetTop_two_apply (t : F 2) (m : Fin 2 → E) :
    jetTop P 2 t m = P 0 (P 1 t (m 0)) (m 1) := rfl

theorem jetTop_three_apply (t : F 3) (m : Fin 3 → E) :
    jetTop P 3 t m = P 0 (P 1 (P 2 t (m 0)) (m 1)) (m 2) := rfl

theorem jetTop_four_apply (t : F 4) (m : Fin 4 → E) :
    jetTop P 4 t m = P 0 (P 1 (P 2 (P 3 t (m 0)) (m 1)) (m 2)) (m 3) := rfl

omit [∀ i, NormedSpace ℝ (F i)] in
theorem snoc_eq_add (n : ℕ) (v : ∀ k : Fin n, F k) (t : F n) :
    (Fin.snoc (α := fun k : Fin (n + 1) => F k) v t) =
      Fin.snoc (α := fun k : Fin (n + 1) => F k) v 0 +
        Fin.snoc (α := fun k : Fin (n + 1) => F k) 0 t := by
  funext k
  cases k using Fin.lastCases with
  | last => simp only [Fin.snoc_last, Pi.add_apply, zero_add]
  | cast k => simp only [Fin.snoc_castSucc, Pi.add_apply, Pi.zero_apply, add_zero]

theorem restr_snoc_zero (n : ℕ) (t : F (n + 1)) :
    restr F n (Fin.snoc (α := fun k : Fin (n + 2) => F k) 0 t) = 0 := by
  funext k
  exact Fin.snoc_castSucc (α := fun k : Fin (n + 2) => F k) (p := 0) (x := t) (i := k)

theorem dV_snoc_zero (n : ℕ) (y : E) (t : F (n + 1)) (w : E) :
    dV P Q n y (Fin.snoc (α := fun k : Fin (n + 2) => F k) 0 t) w =
      Fin.snoc (α := fun k : Fin (n + 1) => F k) 0 (P n t w) := by
  rw [dV_apply]
  funext k
  cases k using Fin.lastCases with
  | last =>
    have h1 : (Fin.snoc (α := fun k : Fin (n + 2) => F k) 0 t) (Fin.last n).succ = t :=
      Fin.snoc_last (α := fun k : Fin (n + 2) => F k) (p := 0) (x := t)
    rw [h1, Fin.snoc_castSucc, Fin.snoc_last, Pi.zero_apply, map_zero, zero_apply, add_zero]
    rfl
  | cast k =>
    have h2 : (Fin.snoc (α := fun k : Fin (n + 2) => F k) 0 t) k.castSucc.succ = 0 :=
      Fin.snoc_castSucc (α := fun k : Fin (n + 2) => F k) (p := 0) (x := t) (i := k.succ)
    rw [h2, Fin.snoc_castSucc, Fin.snoc_castSucc]
    simp

/-- `jetLam` on a jet whose only nonzero value is the top one is the iterated insertion. -/
theorem jetLam_snoc_zero : ∀ (n : ℕ) (y : E) (t : F n),
    jetLam P Q n y (Fin.snoc (α := fun k : Fin (n + 1) => F k) 0 t) = jetTop P n t := by
  intro n
  induction n with
  | zero =>
    intro y t
    rw [jetLam_zero_apply, Fin.snoc_last]
    rfl
  | succ n ih =>
    intro y t
    rw [jetLam_succ_apply, restr_snoc_zero, map_zero, add_zero]
    have hc : (jetLam P Q n y).comp (dV P Q n y (Fin.snoc (α := fun k : Fin (n + 2) => F k) 0 t)) =
        (jetTop P n).comp (P n t) := by
      refine ContinuousLinearMap.ext fun w => ?_
      rw [ContinuousLinearMap.comp_apply, ContinuousLinearMap.comp_apply, dV_snoc_zero, ih]
    rw [hc]
    rfl

/-- **Top-term lemma**: `jetLam P Q n y` depends on the top value exactly through `jetTop P n`. -/
theorem jetLam_snoc (n : ℕ) (y : E) (v : ∀ k : Fin n, F k) (t : F n) :
    jetLam P Q n y (Fin.snoc (α := fun k : Fin (n + 1) => F k) v t) =
      jetTop P n t + jetLam P Q n y (Fin.snoc (α := fun k : Fin (n + 1) => F k) v 0) := by
  rw [snoc_eq_add n v t, map_add, jetLam_snoc_zero, add_comm]

/-- **Coordinate derivatives = iterated insertion of the top covariant derivative + lower-order
terms**: `D^n (G 0) y = jetTop P n (G n y) + jetLam P Q n y (G 0 y, …, G (n-1) y, 0)`. -/
theorem iteratedFDeriv_eq_jetTop_add {Ω : Set E} (hΩ : IsOpen Ω)
    (hQ : ∀ i, ContDiffOn ℝ ∞ (Q i) Ω) {G : ∀ i, E → F i} (hG : ∀ i, ContDiffOn ℝ ∞ (G i) Ω)
    (hrec : ∀ i, ∀ y ∈ Ω, fderiv ℝ (G i) y = P i (G (i + 1) y) + Q i y (G i y)) (n : ℕ)
    {y : E} (hy : y ∈ Ω) :
    iteratedFDeriv ℝ n (G 0) y =
      jetTop P n (G n y) +
        jetLam P Q n y (Fin.snoc (α := fun k : Fin (n + 1) => F k) (fun k : Fin n => G k y) 0) := by
  rw [iteratedFDeriv_eq_jetLam P Q hΩ hQ hG hrec n y hy, ← jetLam_snoc]
  congr 1
  exact (Fin.snoc_init_self (fun k : Fin (n + 1) => G k y)).symm

end CWFlow.JetRec
