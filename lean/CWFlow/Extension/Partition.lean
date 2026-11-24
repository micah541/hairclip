import CWFlow.Extension.GraphChart
import CWFlow.Extension.IBPBounds

/-!
# A partition of unity `Σ ρ_j² = 1` subordinate to a cover by graph charts

Given a finite cover of a compact manifold by graph charts of a Lagrangian immersion, with the second
and third derivatives of the graph functions bounded, we build `ρ_j` with `Σ_j ρ_j² = 1`,
`tsupport ρ_j ⊆ σ_j '' ball 0 2`, and such that each `ρ_j`, read in graph coordinates, is a cutoff
(`IsCutoff`) with a constant `Λ` depending only on `n`-independent data `N₀`, `c₀`, `K₃`.

Structure:
* Part A: pointwise `C²`-bounds `PcBd` for sums, products, scalar compositions.
* Part B: the abstract cutoff `f_j / √(Σ f_l²)` on `ℝⁿ`.
* Part C: derivative bounds for the graph coordinate change `g_l`.
* Part D: the global construction on `M` and the main theorem.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

/-! ### Part A: pointwise second-order bounds -/

section PartA

variable {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]

/-- `f` is `C²` at `x` with `|f x| ≤ a`, `‖Df x‖ ≤ b`, `‖D²f x‖ ≤ c`. -/
def PcBd (f : X → ℝ) (x : X) (a b c : ℝ) : Prop :=
  ContDiffAt ℝ 2 f x ∧ |f x| ≤ a ∧ ‖fderiv ℝ f x‖ ≤ b ∧ ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ c

theorem pc_hasD {Z : Type*} [NormedAddCommGroup Z] [NormedSpace ℝ Z] {f : X → Z} {x : X}
    (hf : ContDiffAt ℝ 2 f x) (w : X) :
    HasFDerivAt (fun y => fderiv ℝ f y w) ((fderiv ℝ (fderiv ℝ f) x).flip w) x := by
  have hd : DifferentiableAt ℝ (fderiv ℝ f) x :=
    (hf.fderiv_right (m := 1) (by norm_num)).differentiableAt (by simp)
  have := hd.hasFDerivAt.clm_apply (hasFDerivAt_const w x)
  simpa using this

theorem pc_d2_le {f : X → ℝ} {x : X} {C : ℝ} (hC : 0 ≤ C)
    (h : ∀ v w, |fderiv ℝ (fderiv ℝ f) x v w| ≤ C * ‖v‖ * ‖w‖) :
    ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ C := by
  refine ContinuousLinearMap.opNorm_le_bound _ hC fun v => ?_
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun w => ?_
  simpa [Real.norm_eq_abs] using h v w

theorem pc_d2_apply_le (f : X → ℝ) (x : X) (v w : X) :
    |fderiv ℝ (fderiv ℝ f) x v w| ≤ ‖fderiv ℝ (fderiv ℝ f) x‖ * ‖v‖ * ‖w‖ := by
  have := ContinuousLinearMap.le_opNorm₂ (fderiv ℝ (fderiv ℝ f) x) v w
  simpa [Real.norm_eq_abs] using this

theorem pc_d1_apply_le (f : X → ℝ) (x : X) (v : X) :
    |fderiv ℝ f x v| ≤ ‖fderiv ℝ f x‖ * ‖v‖ := by
  simpa [Real.norm_eq_abs] using (fderiv ℝ f x).le_opNorm v

theorem PcBd.mono {f : X → ℝ} {x : X} {a b c a' b' c' : ℝ} (h : PcBd f x a b c) (ha : a ≤ a')
    (hb : b ≤ b') (hc : c ≤ c') : PcBd f x a' b' c' :=
  ⟨h.1, h.2.1.trans ha, h.2.2.1.trans hb, h.2.2.2.trans hc⟩

theorem PcBd.congr {f g : X → ℝ} {x : X} {a b c : ℝ} (h : PcBd f x a b c) (hfg : f =ᶠ[𝓝 x] g) :
    PcBd g x a b c := by
  refine ⟨h.1.congr_of_eventuallyEq hfg.symm, ?_, ?_, ?_⟩
  · rw [← hfg.eq_of_nhds]; exact h.2.1
  · rw [← hfg.fderiv_eq]; exact h.2.2.1
  · rw [← hfg.fderiv.fderiv_eq]; exact h.2.2.2

theorem PcBd.congr_fun {f g : X → ℝ} {x : X} {a b c : ℝ} (h : PcBd f x a b c) (hfg : ∀ y, f y = g y) :
    PcBd g x a b c := by
  have : f = g := funext hfg
  subst this; exact h

theorem pc_bd_const (c₀ : ℝ) (x : X) : PcBd (fun _ => c₀) x |c₀| 0 0 := by
  have h : fderiv ℝ (fun _ : X => c₀) = fun _ => 0 := funext fun y => (hasFDerivAt_const c₀ y).fderiv
  have h' : fderiv ℝ (fun _ : X => (0 : X →L[ℝ] ℝ)) x = 0 := (hasFDerivAt_const _ x).fderiv
  refine ⟨contDiffAt_const, le_rfl, ?_, ?_⟩
  · rw [h]; exact norm_zero.le
  · rw [h, h']; exact ContinuousLinearMap.opNorm_zero.le

theorem pc_bd_zero (x : X) : PcBd (fun _ => (0 : ℝ)) x 0 0 0 := by
  simpa using pc_bd_const (0 : ℝ) x

theorem PcBd.add {f g : X → ℝ} {x : X} {a b c a' b' c' : ℝ} (hf : PcBd f x a b c)
    (hg : PcBd g x a' b' c') : PcBd (fun y => f y + g y) x (a + a') (b + b') (c + c') := by
  have hfd : DifferentiableAt ℝ f x := hf.1.differentiableAt (by norm_num)
  have hgd : DifferentiableAt ℝ g x := hg.1.differentiableAt (by norm_num)
  refine ⟨hf.1.add hg.1, ?_, ?_, ?_⟩
  · exact (abs_add_le _ _).trans (add_le_add hf.2.1 hg.2.1)
  · rw [fderiv_fun_add hfd hgd]
    exact (norm_add_le _ _).trans (add_le_add hf.2.2.1 hg.2.2.1)
  · have hev : fderiv ℝ (fun y => f y + g y) =ᶠ[𝓝 x] fun y => fderiv ℝ f y + fderiv ℝ g y := by
      filter_upwards [hf.1.eventually (by simp), hg.1.eventually (by simp)] with y h1 h2
      exact fderiv_fun_add (h1.differentiableAt (by norm_num)) (h2.differentiableAt (by norm_num))
    have hdf := hf.1.fderiv_right (m := 1) (by norm_num)
    have hdg := hg.1.fderiv_right (m := 1) (by norm_num)
    rw [hev.fderiv_eq, fderiv_fun_add (hdf.differentiableAt (by simp)) (hdg.differentiableAt (by simp))]
    exact (norm_add_le (fderiv ℝ (fderiv ℝ f) x) (fderiv ℝ (fderiv ℝ g) x)).trans
      (add_le_add hf.2.2.2 hg.2.2.2)

theorem PcBd.mul {f g : X → ℝ} {x : X} {a b c a' b' c' : ℝ} (hf : PcBd f x a b c)
    (hg : PcBd g x a' b' c') :
    PcBd (fun y => f y * g y) x (a * a') (a * b' + a' * b) (a * c' + 2 * b * b' + a' * c) := by
  have hfd : DifferentiableAt ℝ f x := hf.1.differentiableAt (by norm_num)
  have hgd : DifferentiableAt ℝ g x := hg.1.differentiableAt (by norm_num)
  have ha : 0 ≤ a := (abs_nonneg _).trans hf.2.1
  have ha' : 0 ≤ a' := (abs_nonneg _).trans hg.2.1
  have hb : 0 ≤ b := (norm_nonneg _).trans hf.2.2.1
  have hb' : 0 ≤ b' := (norm_nonneg _).trans hg.2.2.1
  have hc : 0 ≤ c := le_trans (ContinuousLinearMap.opNorm_nonneg _) hf.2.2.2
  have hc' : 0 ≤ c' := le_trans (ContinuousLinearMap.opNorm_nonneg _) hg.2.2.2
  refine ⟨hf.1.mul hg.1, ?_, ?_, ?_⟩
  · rw [abs_mul]; exact mul_le_mul hf.2.1 hg.2.1 (abs_nonneg _) ha
  · rw [fderiv_fun_mul hfd hgd]
    refine (norm_add_le _ _).trans ?_
    rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs]
    exact add_le_add (mul_le_mul hf.2.1 hg.2.2.1 (norm_nonneg _) ha)
      (mul_le_mul hg.2.1 hf.2.2.1 (norm_nonneg _) ha')
  · refine pc_d2_le (by positivity) fun v w => ?_
    have hev : (fun y => fderiv ℝ (fun y => f y * g y) y w) =ᶠ[𝓝 x]
        fun y => f y * fderiv ℝ g y w + g y * fderiv ℝ f y w := by
      filter_upwards [hf.1.eventually (by simp), hg.1.eventually (by simp)] with y h1 h2
      rw [fderiv_fun_mul (h1.differentiableAt (by norm_num)) (h2.differentiableAt (by norm_num))]
      simp
    have hD := (pc_hasD (f := fun y => f y * g y) (hf.1.mul hg.1) w)
    have hE : HasFDerivAt (fun y => f y * fderiv ℝ g y w + g y * fderiv ℝ f y w)
        ((f x • (fderiv ℝ (fderiv ℝ g) x).flip w + fderiv ℝ g x w • fderiv ℝ f x) +
          (g x • (fderiv ℝ (fderiv ℝ f) x).flip w + fderiv ℝ f x w • fderiv ℝ g x)) x :=
      (hfd.hasFDerivAt.mul (pc_hasD hg.1 w)).add (hgd.hasFDerivAt.mul (pc_hasD hf.1 w))
    have hEq := (hD.congr_of_eventuallyEq hev.symm |>.unique hE)
    have key : fderiv ℝ (fderiv ℝ fun y => f y * g y) x v w =
        f x * fderiv ℝ (fderiv ℝ g) x v w + fderiv ℝ g x w * fderiv ℝ f x v +
        (g x * fderiv ℝ (fderiv ℝ f) x v w + fderiv ℝ f x w * fderiv ℝ g x v) := by
      have := congrArg (fun L => L v) hEq
      simpa using this
    rw [key]
    have e1 := pc_d2_apply_le g x v w
    have e2 := pc_d2_apply_le f x v w
    have e3 := pc_d1_apply_le g x w
    have e4 := pc_d1_apply_le f x v
    have e5 := pc_d1_apply_le f x w
    have e6 := pc_d1_apply_le g x v
    have hx1 : |f x| ≤ a := hf.2.1
    have hx2 : |g x| ≤ a' := hg.2.1
    have hnv := norm_nonneg v
    have hnw := norm_nonneg w
    calc _ ≤ |f x| * |fderiv ℝ (fderiv ℝ g) x v w| + |fderiv ℝ g x w| * |fderiv ℝ f x v| +
          (|g x| * |fderiv ℝ (fderiv ℝ f) x v w| + |fderiv ℝ f x w| * |fderiv ℝ g x v|) := by
          refine (abs_add_le _ _).trans (add_le_add ((abs_add_le _ _).trans ?_) ((abs_add_le _ _).trans ?_))
          · rw [abs_mul, abs_mul]
          · rw [abs_mul, abs_mul]
      _ ≤ a * (c' * ‖v‖ * ‖w‖) + (b' * ‖w‖) * (b * ‖v‖) +
          (a' * (c * ‖v‖ * ‖w‖) + (b * ‖w‖) * (b' * ‖v‖)) := by
          have h1 : |fderiv ℝ (fderiv ℝ g) x v w| ≤ c' * ‖v‖ * ‖w‖ :=
            e1.trans (by gcongr; exact hg.2.2.2)
          have h2 : |fderiv ℝ (fderiv ℝ f) x v w| ≤ c * ‖v‖ * ‖w‖ :=
            e2.trans (by gcongr; exact hf.2.2.2)
          have h3 : |fderiv ℝ g x w| ≤ b' * ‖w‖ := e3.trans (by gcongr; exact hg.2.2.1)
          have h4 : |fderiv ℝ f x v| ≤ b * ‖v‖ := e4.trans (by gcongr; exact hf.2.2.1)
          have h5 : |fderiv ℝ f x w| ≤ b * ‖w‖ := e5.trans (by gcongr; exact hf.2.2.1)
          have h6 : |fderiv ℝ g x v| ≤ b' * ‖v‖ := e6.trans (by gcongr; exact hg.2.2.1)
          gcongr
      _ = _ := by ring

theorem PcBd.sum {ι : Type*} (s : Finset ι) (h : ι → X → ℝ) {x : X} {a b c : ℝ}
    (hh : ∀ l ∈ s, PcBd (h l) x a b c) :
    PcBd (fun y => ∑ l ∈ s, h l y) x (s.card * a) (s.card * b) (s.card * c) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using pc_bd_zero x
  | insert l s hl ih =>
    have := (hh l (Finset.mem_insert_self l s)).add
      (ih fun m hm => hh m (Finset.mem_insert_of_mem hm))
    simp only [Finset.sum_insert hl, Finset.card_insert_of_notMem hl]
    push_cast
    exact this.mono (le_of_eq (by ring)) (le_of_eq (by ring)) (le_of_eq (by ring))

theorem PcBd.comp_real {S : X → ℝ} {u : ℝ → ℝ} {x : X} {a b c ua ub uc : ℝ} (hS : PcBd S x a b c)
    (hu : ContDiffAt ℝ 2 u (S x)) (h0 : |u (S x)| ≤ ua) (h1 : |deriv u (S x)| ≤ ub)
    (h2 : |deriv (deriv u) (S x)| ≤ uc) :
    PcBd (fun y => u (S y)) x ua (ub * b) (uc * b ^ 2 + ub * c) := by
  have hSd : DifferentiableAt ℝ S x := hS.1.differentiableAt (by norm_num)
  have hud : ∀ s, ContDiffAt ℝ 2 u s → HasDerivAt u (deriv u s) s := fun s h =>
    (h.differentiableAt (by norm_num)).hasDerivAt
  have hub : 0 ≤ ub := (abs_nonneg _).trans h1
  have huc : 0 ≤ uc := (abs_nonneg _).trans h2
  have hb : 0 ≤ b := (norm_nonneg _).trans hS.2.2.1
  have hc : 0 ≤ c := le_trans (ContinuousLinearMap.opNorm_nonneg _) hS.2.2.2
  have hcomp : ContDiffAt ℝ 2 (fun y => u (S y)) x := hu.comp x hS.1
  refine ⟨hcomp, h0, ?_, ?_⟩
  · have := (hud _ hu).comp_hasFDerivAt x hSd.hasFDerivAt
    rw [show (fun y => u (S y)) = u ∘ S from rfl, this.fderiv, norm_smul, Real.norm_eq_abs]
    exact mul_le_mul h1 hS.2.2.1 (norm_nonneg _) hub
  · refine pc_d2_le (by positivity) fun v w => ?_
    have hev : (fun y => fderiv ℝ (fun y => u (S y)) y w) =ᶠ[𝓝 x]
        fun y => deriv u (S y) * fderiv ℝ S y w := by
      have hev1 : ∀ᶠ y in 𝓝 x, ContDiffAt ℝ 2 S y := hS.1.eventually (by simp)
      have hev2 : ∀ᶠ y in 𝓝 x, ContDiffAt ℝ 2 u (S y) :=
        hS.1.continuousAt.tendsto.eventually (hu.eventually (by simp))
      filter_upwards [hev1, hev2] with y h1 h2
      have := (hud _ h2).comp_hasFDerivAt y (h1.differentiableAt (by norm_num)).hasFDerivAt
      rw [show (fun y => u (S y)) = u ∘ S from rfl, this.fderiv]
      simp
    have hD := pc_hasD hcomp w
    have hDu : HasDerivAt (deriv u) (deriv (deriv u) (S x)) (S x) :=
      ((hu.derivWithin (m := 1) (by norm_num)).differentiableAt (by norm_num)).hasDerivAt
    have hE : HasFDerivAt (fun y => deriv u (S y) * fderiv ℝ S y w)
        (deriv u (S x) • (fderiv ℝ (fderiv ℝ S) x).flip w +
          fderiv ℝ S x w • (deriv (deriv u) (S x) • fderiv ℝ S x)) x :=
      (hDu.comp_hasFDerivAt x hSd.hasFDerivAt).mul (pc_hasD hS.1 w)
    have hEq := (hD.congr_of_eventuallyEq hev.symm).unique hE
    have key : fderiv ℝ (fderiv ℝ fun y => u (S y)) x v w =
        deriv u (S x) * fderiv ℝ (fderiv ℝ S) x v w +
          fderiv ℝ S x w * (deriv (deriv u) (S x) * fderiv ℝ S x v) := by
      have := congrArg (fun L => L v) hEq
      simpa using this
    rw [key]
    have e1 := pc_d2_apply_le S x v w
    have e3 := pc_d1_apply_le S x w
    have e6 := pc_d1_apply_le S x v
    have hnv := norm_nonneg v
    have hnw := norm_nonneg w
    have h3 : |fderiv ℝ (fderiv ℝ S) x v w| ≤ c * ‖v‖ * ‖w‖ := e1.trans (by gcongr; exact hS.2.2.2)
    have h4 : |fderiv ℝ S x w| ≤ b * ‖w‖ := e3.trans (by gcongr; exact hS.2.2.1)
    have h5 : |fderiv ℝ S x v| ≤ b * ‖v‖ := e6.trans (by gcongr; exact hS.2.2.1)
    calc _ ≤ |deriv u (S x)| * |fderiv ℝ (fderiv ℝ S) x v w| +
          |fderiv ℝ S x w| * (|deriv (deriv u) (S x)| * |fderiv ℝ S x v|) := by
          refine (abs_add_le _ _).trans ?_
          rw [abs_mul, abs_mul, abs_mul]
      _ ≤ ub * (c * ‖v‖ * ‖w‖) + (b * ‖w‖) * (uc * (b * ‖v‖)) := by gcongr
      _ = _ := by ring

theorem PcBd.comp_vec {Y : Type*} [NormedAddCommGroup Y] [NormedSpace ℝ Y] {η : Y → ℝ} {g : X → Y}
    {x : X} {a B₁ B₂ L K : ℝ} (hη : PcBd η (g x) a B₁ B₂) (hg : ContDiffAt ℝ 2 g x)
    (hL : ‖fderiv ℝ g x‖ ≤ L) (hK : ‖fderiv ℝ (fderiv ℝ g) x‖ ≤ K) :
    PcBd (fun y => η (g y)) x a (B₁ * L) (B₂ * L ^ 2 + B₁ * K) := by
  have hgd : DifferentiableAt ℝ g x := hg.differentiableAt (by norm_num)
  have hηd : DifferentiableAt ℝ η (g x) := hη.1.differentiableAt (by norm_num)
  have hB1 : 0 ≤ B₁ := (norm_nonneg _).trans hη.2.2.1
  have hB2 : 0 ≤ B₂ := le_trans (ContinuousLinearMap.opNorm_nonneg _) hη.2.2.2
  have hLn : 0 ≤ L := (norm_nonneg _).trans hL
  have hKn : 0 ≤ K := le_trans (ContinuousLinearMap.opNorm_nonneg _) hK
  have hcomp : ContDiffAt ℝ 2 (fun y => η (g y)) x := hη.1.comp x hg
  refine ⟨hcomp, hη.2.1, ?_, ?_⟩
  · have := hηd.hasFDerivAt.comp x hgd.hasFDerivAt
    rw [show (fun y => η (g y)) = η ∘ g from rfl, this.fderiv]
    exact (ContinuousLinearMap.opNorm_comp_le _ _).trans
      (mul_le_mul hη.2.2.1 hL (norm_nonneg _) hB1)
  · refine pc_d2_le (by positivity) fun v w => ?_
    have hev : (fun y => fderiv ℝ (fun y => η (g y)) y w) =ᶠ[𝓝 x]
        fun y => fderiv ℝ η (g y) (fderiv ℝ g y w) := by
      have hev1 : ∀ᶠ y in 𝓝 x, ContDiffAt ℝ 2 g y := hg.eventually (by simp)
      have hev2 : ∀ᶠ y in 𝓝 x, ContDiffAt ℝ 2 η (g y) :=
        hg.continuousAt.tendsto.eventually (hη.1.eventually (by simp))
      filter_upwards [hev1, hev2] with y h1 h2
      have := (h2.differentiableAt (by norm_num)).hasFDerivAt.comp y
        (h1.differentiableAt (by norm_num)).hasFDerivAt
      rw [show (fun y => η (g y)) = η ∘ g from rfl, this.fderiv]
      rfl
    have hD := pc_hasD hcomp w
    have hdη : HasFDerivAt (fderiv ℝ η) (fderiv ℝ (fderiv ℝ η) (g x)) (g x) :=
      ((hη.1.fderiv_right (m := 1) (by norm_num)).differentiableAt (by simp)).hasFDerivAt
    have hE : HasFDerivAt (fun y => fderiv ℝ η (g y) (fderiv ℝ g y w))
        ((fderiv ℝ η (g x)).comp ((fderiv ℝ (fderiv ℝ g) x).flip w) +
          ((fderiv ℝ (fderiv ℝ η) (g x)).comp (fderiv ℝ g x)).flip (fderiv ℝ g x w)) x :=
      (hdη.comp x hgd.hasFDerivAt).clm_apply (pc_hasD hg w)
    have hEq := (hD.congr_of_eventuallyEq hev.symm).unique hE
    have key : fderiv ℝ (fderiv ℝ fun y => η (g y)) x v w =
        fderiv ℝ η (g x) (fderiv ℝ (fderiv ℝ g) x v w) +
          fderiv ℝ (fderiv ℝ η) (g x) (fderiv ℝ g x v) (fderiv ℝ g x w) := by
      have := congrArg (fun L => L v) hEq
      simpa using this
    rw [key]
    have e1 := pc_d2_apply_le η (g x) (fderiv ℝ g x v) (fderiv ℝ g x w)
    have e2 := pc_d1_apply_le η (g x) (fderiv ℝ (fderiv ℝ g) x v w)
    have e3 : ‖fderiv ℝ (fderiv ℝ g) x v w‖ ≤ ‖fderiv ℝ (fderiv ℝ g) x‖ * ‖v‖ * ‖w‖ :=
      ContinuousLinearMap.le_opNorm₂ _ v w
    have e4 : ‖fderiv ℝ g x v‖ ≤ ‖fderiv ℝ g x‖ * ‖v‖ := (fderiv ℝ g x).le_opNorm v
    have e5 : ‖fderiv ℝ g x w‖ ≤ ‖fderiv ℝ g x‖ * ‖w‖ := (fderiv ℝ g x).le_opNorm w
    have hnv := norm_nonneg v
    have hnw := norm_nonneg w
    have f1 : ‖fderiv ℝ g x v‖ ≤ L * ‖v‖ := e4.trans (by gcongr)
    have f2 : ‖fderiv ℝ g x w‖ ≤ L * ‖w‖ := e5.trans (by gcongr)
    have f3 : ‖fderiv ℝ (fderiv ℝ g) x v w‖ ≤ K * ‖v‖ * ‖w‖ := e3.trans (by gcongr)
    calc _ ≤ |fderiv ℝ η (g x) (fderiv ℝ (fderiv ℝ g) x v w)| +
          |fderiv ℝ (fderiv ℝ η) (g x) (fderiv ℝ g x v) (fderiv ℝ g x w)| := abs_add_le _ _
      _ ≤ B₁ * (K * ‖v‖ * ‖w‖) + B₂ * (L * ‖v‖) * (L * ‖w‖) := by
          refine add_le_add ?_ ?_
          · refine e2.trans ?_
            exact mul_le_mul hη.2.2.1 f3 (norm_nonneg _) hB1
          · refine e1.trans ?_
            exact mul_le_mul (mul_le_mul hη.2.2.2 f1 (norm_nonneg _) hB2) f2 (norm_nonneg _)
              (by positivity)
      _ = _ := by ring

/-- The function `s ↦ s ^ (-1/2)`. -/
noncomputable def pcU : ℝ → ℝ := fun s => s ^ (-(1 / 2 : ℝ))

theorem pc_u_facts {s : ℝ} (hs : 1 ≤ s) :
    ContDiffAt ℝ 2 pcU s ∧ |pcU s| ≤ 1 ∧ |deriv pcU s| ≤ 1 / 2 ∧ |deriv (deriv pcU) s| ≤ 3 / 4 := by
  have hs0 : 0 < s := zero_lt_one.trans_le hs
  have hd1 : deriv pcU = fun t => (-(1 / 2 : ℝ)) * t ^ (-(1 / 2 : ℝ) - 1) := by
    funext t
    exact Real.deriv_rpow_const t _
  have hd2 : HasDerivAt (fun t : ℝ => (-(1 / 2 : ℝ)) * t ^ (-(1 / 2 : ℝ) - 1))
      ((-(1 / 2 : ℝ)) * ((-(1 / 2 : ℝ) - 1) * s ^ (-(1 / 2 : ℝ) - 1 - 1))) s :=
    (Real.hasDerivAt_rpow_const (Or.inl hs0.ne')).const_mul _
  have hp : ∀ q : ℝ, q ≤ 0 → s ^ q ≤ 1 := fun q hq => Real.rpow_le_one_of_one_le_of_nonpos hs hq
  have hp0 : ∀ q : ℝ, 0 ≤ s ^ q := fun q => Real.rpow_nonneg hs0.le q
  refine ⟨Real.contDiffAt_rpow_const_of_ne hs0.ne', ?_, ?_, ?_⟩
  · show |s ^ (-(1 / 2 : ℝ))| ≤ 1
    rw [abs_of_nonneg (hp0 _)]; exact hp _ (by norm_num)
  · rw [hd1, abs_mul, abs_of_nonneg (hp0 _)]
    have := hp (-(1 / 2 : ℝ) - 1) (by norm_num)
    have h0 := hp0 (-(1 / 2 : ℝ) - 1)
    rw [abs_neg, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 2)]
    nlinarith
  · rw [hd1, hd2.deriv, abs_mul, abs_mul, abs_of_nonneg (hp0 _)]
    have := hp (-(1 / 2 : ℝ) - 1 - 1) (by norm_num)
    have h0 := hp0 (-(1 / 2 : ℝ) - 1 - 1)
    rw [abs_neg, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 2),
      abs_of_neg (by norm_num : (-(1 / 2 : ℝ) - 1) < 0)]
    nlinarith

end PartA

/-! ### Part B: the abstract cutoff `f_j / √(Σ f_l²)` -/

theorem pc_cd2 {X F : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [NormedAddCommGroup F]
    [NormedSpace ℝ F] {f : X → F} {x : X}
    (h : ContDiffAt ℝ ∞ f x) : ContDiffAt ℝ 2 f x :=
  h.of_le (by first | exact WithTop.coe_le_coe.2 le_top | norm_cast; exact le_top | simp; exact_mod_cast le_top)

/-- The cutoff constant of `pc_abstract`. -/
noncomputable def pcLam (N₀ : ℕ) (B₁ B₂ C₁ C₂ : ℝ) : ℝ :=
  max 1 (max (1 / 2 * (N₀ * (2 * C₁)) + B₁)
    (3 / 4 * (N₀ * (2 * C₁)) ^ 2 + 1 / 2 * (N₀ * (2 * C₂ + 2 * C₁ ^ 2)) +
      B₁ * (N₀ * (2 * C₁)) + B₂))

section PartB

variable {n : ℕ}

open Classical in
theorem pc_abstract {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη : ContDiff ℝ ∞ η) (hη0 : ∀ x, 0 ≤ η x)
    (hη1 : ∀ x, η x ≤ 1) (hηs : ∀ x : EuclideanSpace ℝ (Fin n), 3 / 2 ≤ ‖x‖ → η x = 0)
    {B₁ B₂ : ℝ} (hB1 : ∀ x, ‖fderiv ℝ η x‖ ≤ B₁) (hB2 : ∀ x, ‖fderiv ℝ (fderiv ℝ η) x‖ ≤ B₂)
    {N N₀ : ℕ} (hN : N ≤ N₀) {C₁ C₂ : ℝ} (f : Fin N → EuclideanSpace ℝ (Fin n) → ℝ)
    (hfs : ∀ l, ContDiffOn ℝ ∞ (f l) (ball 0 6))
    (hf : ∀ l, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, PcBd (f l) x 1 C₁ C₂)
    (j : Fin N) (hfj : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, f j x = η x)
    (hS : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, 1 ≤ ∑ l, f l x ^ 2) :
    IsCutoff (fun x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then
      f j x * (∑ l, f l x ^ 2) ^ (-(1 / 2 : ℝ)) else 0) (pcLam N₀ B₁ B₂ C₁ C₂) := by
  set r : EuclideanSpace ℝ (Fin n) → ℝ := fun x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then
      f j x * (∑ l, f l x ^ 2) ^ (-(1 / 2 : ℝ)) else 0 with hr
  have hrz : ∀ x : EuclideanSpace ℝ (Fin n), 3 / 2 ≤ ‖x‖ → r x = 0 := by
    intro x hx
    by_cases hm : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6
    · simp only [r, if_pos hm, hfj x hm, hηs x hx, zero_mul]
    · simp only [r, if_neg hm]
  have hform : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      r =ᶠ[𝓝 x] fun y => η y * pcU (∑ l, f l y ^ 2) := by
    intro x hx
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy
    simp only [r, if_pos hy, hfj y hy, pcU]
  have hC1 : 0 ≤ C₁ := by
    have := (hf j 0 (mem_ball_self (by norm_num))).2.2.1
    exact (norm_nonneg _).trans this
  have hC2 : 0 ≤ C₂ := by
    have := (hf j 0 (mem_ball_self (by norm_num))).2.2.2
    exact le_trans (ContinuousLinearMap.opNorm_nonneg _) this
  have hNN : (N : ℝ) ≤ N₀ := by exact_mod_cast hN
  have hNn : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have hΛ1 : 1 ≤ pcLam N₀ B₁ B₂ C₁ C₂ := le_max_left _ _
  -- bounds on the sum of squares
  have hSbd : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      PcBd (fun y => ∑ l, f l y ^ 2) x N₀ (N₀ * (2 * C₁)) (N₀ * (2 * C₂ + 2 * C₁ ^ 2)) := by
    intro x hx
    have h1 : ∀ l ∈ (Finset.univ : Finset (Fin N)), PcBd (fun y => f l y ^ 2) x (1 * 1)
        (1 * C₁ + 1 * C₁) (1 * C₂ + 2 * C₁ * C₁ + 1 * C₂) := by
      intro l _
      exact ((hf l x hx).mul (hf l x hx)).congr_fun fun y => (sq _).symm
    have := PcBd.sum (Finset.univ : Finset (Fin N)) (fun l y => f l y ^ 2) h1
    rw [Finset.card_univ, Fintype.card_fin] at this
    refine this.mono ?_ ?_ ?_
    · nlinarith
    · nlinarith
    · nlinarith
  have hrbd : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      PcBd r x 1 (pcLam N₀ B₁ B₂ C₁ C₂) (pcLam N₀ B₁ B₂ C₁ C₂) := by
    intro x hx
    have hSx := hS x hx
    obtain ⟨hu1, hu2, hu3, hu4⟩ := pc_u_facts hSx
    have hT := (hSbd x hx).comp_real hu1 hu2 hu3 hu4
    have hηb : PcBd η x 1 B₁ B₂ :=
      ⟨pc_cd2 hη.contDiffAt, by rw [abs_of_nonneg (hη0 x)]; exact hη1 x, hB1 x, hB2 x⟩
    have hm := hηb.mul hT
    refine (hm.congr (hform x hx).symm).mono (by norm_num) ?_ ?_
    · exact le_trans (le_of_eq (by ring)) (le_trans (le_max_left _ _) (le_max_right _ _))
    · exact le_trans (le_of_eq (by ring)) (le_trans (le_max_right _ _) (le_max_right _ _))
  have hrz' : ∀ x : EuclideanSpace ℝ (Fin n), 6 ≤ ‖x‖ → r =ᶠ[𝓝 x] fun _ => (0 : ℝ) := by
    intro x hx
    filter_upwards [(isOpen_lt continuous_const continuous_norm).mem_nhds
      (show (3 / 2 : ℝ) < ‖x‖ by linarith)] with y hy using hrz y hy.le
  have hall : ∀ x, PcBd r x 1 (pcLam N₀ B₁ B₂ C₁ C₂) (pcLam N₀ B₁ B₂ C₁ C₂) := by
    intro x
    by_cases hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6
    · exact hrbd x hx
    · have hx6 : 6 ≤ ‖x‖ := by simpa using hx
      exact ((pc_bd_zero x).congr (hrz' x hx6).symm).mono (by norm_num) (by linarith) (by linarith)
  have hsm : ∀ x, ContDiffAt ℝ ∞ r x := by
    intro x
    by_cases hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6
    · have hat : ∀ l, ContDiffAt ℝ ∞ (f l) x := fun l =>
        (hfs l).contDiffAt (isOpen_ball.mem_nhds hx)
      have hSs : ContDiffAt ℝ ∞ (fun y => ∑ l, f l y ^ 2) x :=
        ContDiffAt.sum fun l _ => (hat l).pow 2
      have hne : (∑ l, f l x ^ 2) ≠ 0 := by
        have := hS x hx; linarith
      have hT : ContDiffAt ℝ ∞ (fun y => pcU (∑ l, f l y ^ 2)) x :=
        ContDiffAt.comp (g := pcU) x (Real.contDiffAt_rpow_const_of_ne hne) hSs
      exact (hη.contDiffAt.mul hT).congr_of_eventuallyEq (hform x hx)
    · have hx6 : 6 ≤ ‖x‖ := by simpa using hx
      exact contDiffAt_const.congr_of_eventuallyEq (hrz' x hx6)
  refine ⟨contDiff_iff_contDiffAt.2 hsm, ?_, fun x => (hall x).2.1, fun x => (hall x).2.2.1,
    fun x => (hall x).2.2.2⟩
  have hsup : tsupport r ⊆ closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2) := by
    refine closure_minimal (fun y hy => ?_) isClosed_closedBall
    by_contra hn
    rw [mem_closedBall, dist_zero_right, not_le] at hn
    exact hy (hrz y hn.le)
  exact hsup.trans (closedBall_subset_ball (by norm_num))

end PartB

/-! ### Part C: derivative bounds for the graph coordinate change -/

section PartC

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem pc_gradE_eq (φ : EuclideanSpace ℝ (Fin n) → ℝ) (y : EuclideanSpace ℝ (Fin n)) :
    gradE φ y = (InnerProductSpace.toDual ℝ (EuclideanSpace ℝ (Fin n))).symm (fderiv ℝ φ y) := by
  rw [eq_comm, LinearIsometryEquiv.symm_apply_eq]
  ext v
  rw [InnerProductSpace.toDual_apply_apply]
  have hv : v = ∑ c, v c • eb c := by
    ext l
    simp [eb, Pi.single_apply]
  conv_lhs => rw [hv]
  rw [map_sum]
  simp [gradE, PiLp.inner_apply, mul_comm]

theorem pc_g_affine {X Z W : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [NormedAddCommGroup Z]
    [NormedSpace ℝ Z] [NormedAddCommGroup W] [NormedSpace ℝ W] (f : X → Z) (A : X →L[ℝ] W)
    (B : Z →L[ℝ] W) (c : W) {x : X} (hf : ContDiffAt ℝ 2 f x) :
    ContDiffAt ℝ 2 (fun y => c + A y + B (f y)) x ∧
      fderiv ℝ (fun y => c + A y + B (f y)) x = A + B.comp (fderiv ℝ f x) ∧
      ∀ v w, fderiv ℝ (fderiv ℝ (fun y => c + A y + B (f y))) x v w =
        B (fderiv ℝ (fderiv ℝ f) x v w) := by
  have hd : ∀ y, DifferentiableAt ℝ f y → HasFDerivAt (fun y => c + A y + B (f y))
      (A + B.comp (fderiv ℝ f y)) y := by
    intro y hy
    have := ((hasFDerivAt_const c y).add A.hasFDerivAt).add (B.hasFDerivAt.comp y hy.hasFDerivAt)
    exact this.congr_fderiv (by simp)
  refine ⟨(contDiffAt_const.add A.contDiff.contDiffAt).add (B.contDiff.contDiffAt.comp x hf),
    (hd x (hf.differentiableAt (by norm_num))).fderiv, fun v w => ?_⟩
  have hcd : ContDiffAt ℝ 2 (fun y => c + A y + B (f y)) x :=
    (contDiffAt_const.add A.contDiff.contDiffAt).add (B.contDiff.contDiffAt.comp x hf)
  have hev : (fun y => fderiv ℝ (fun y => c + A y + B (f y)) y w) =ᶠ[𝓝 x]
      fun y => A w + B (fderiv ℝ f y w) := by
    filter_upwards [hf.eventually (by simp)] with y hy
    rw [(hd y (hy.differentiableAt (by norm_num))).fderiv]
    simp
  have hD := pc_hasD hcd w
  have hE : HasFDerivAt (fun y => A w + B (fderiv ℝ f y w))
      (B.comp ((fderiv ℝ (fderiv ℝ f) x).flip w)) x := by
    have := (hasFDerivAt_const (A w) x).add (B.hasFDerivAt.comp x (pc_hasD hf w))
    exact this.congr_fderiv (by simp)
  have hEq := (hD.congr_of_eventuallyEq hev.symm).unique hE
  have := congrArg (fun L => L v) hEq
  simpa using this

theorem pc_horCL_le (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (w : V) : ‖horCL U w‖ ≤ ‖w‖ := by
  have h := cv_horC_lip U (0 : V) w 0
  have e1 : horC U 0 w = horCL U w := by
    show horCL U (w - 0) = _
    rw [sub_zero]
  have e2 : horC U 0 0 = 0 := by
    show horCL U (0 - 0) = 0
    simp
  rw [e1, e2, sub_zero, sub_zero] at h
  exact h

theorem pc_norm_D2_eq {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] (φ : X → ℝ) (x : X) :
    ‖fderiv ℝ (fderiv ℝ φ) x‖ = ‖iteratedFDeriv ℝ 2 φ x‖ := by
  have h1 : ‖iteratedFDeriv ℝ 1 (fderiv ℝ φ) x‖ = ‖iteratedFDeriv ℝ 2 φ x‖ :=
    norm_iteratedFDeriv_fderiv
  rw [norm_iteratedFDeriv_one] at h1
  exact h1

theorem pc_D3_apply_le {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] (φ : X → ℝ) (x v w : X) :
    ‖fderiv ℝ (fderiv ℝ (fderiv ℝ φ)) x v w‖ ≤ ‖iteratedFDeriv ℝ 3 φ x‖ * ‖v‖ * ‖w‖ := by
  have h := iteratedFDeriv_two_apply (𝕜 := ℝ) (fderiv ℝ φ) x ![v, w]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one] at h
  rw [← h]
  have h2 : ‖iteratedFDeriv ℝ 2 (fderiv ℝ φ) x‖ = ‖iteratedFDeriv ℝ 3 φ x‖ :=
    norm_iteratedFDeriv_fderiv
  have := (iteratedFDeriv ℝ 2 (fderiv ℝ φ) x).le_opNorm ![v, w]
  simpa [Fin.prod_univ_two, h2, mul_assoc] using this

/-- Derivative bounds for `c + A y + B (Dφ y)` from bounds on `D²φ`, `D³φ`. -/
theorem pc_g_aff_bounds {X W : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup W] [NormedSpace ℝ W] (φ : X → ℝ) (A : X →L[ℝ] W)
    (B : StrongDual ℝ X →L[ℝ] W) (c : W) {x : X} (hf2 : ContDiffAt ℝ 2 (fderiv ℝ φ) x)
    (hAn : ∀ v, ‖A v‖ ≤ ‖v‖) (hBn : ∀ z, ‖B z‖ ≤ ‖z‖) {c₀ K : ℝ}
    (h2 : ‖iteratedFDeriv ℝ 2 φ x‖ ≤ c₀) (h3 : ‖iteratedFDeriv ℝ 3 φ x‖ ≤ K) :
    ContDiffAt ℝ 2 (fun y => c + A y + B (fderiv ℝ φ y)) x ∧
      ‖fderiv ℝ (fun y => c + A y + B (fderiv ℝ φ y)) x‖ ≤ 1 + c₀ ∧
      ‖fderiv ℝ (fderiv ℝ (fun y => c + A y + B (fderiv ℝ φ y))) x‖ ≤ K := by
  obtain ⟨hc, hd1, hd2⟩ := pc_g_affine (fderiv ℝ φ) A B c hf2
  have hc₀ : 0 ≤ c₀ := (norm_nonneg _).trans h2
  have hK : 0 ≤ K := (norm_nonneg _).trans h3
  refine ⟨hc, ?_, ?_⟩
  · rw [hd1]
    refine ContinuousLinearMap.opNorm_le_bound _ (by linarith) fun v => ?_
    have e1 : ‖fderiv ℝ (fderiv ℝ φ) x v‖ ≤ ‖fderiv ℝ (fderiv ℝ φ) x‖ * ‖v‖ :=
      (fderiv ℝ (fderiv ℝ φ) x).le_opNorm v
    have e2 : ‖fderiv ℝ (fderiv ℝ φ) x‖ ≤ c₀ := by rw [pc_norm_D2_eq]; exact h2
    have e3 := hAn v
    have e4 := hBn (fderiv ℝ (fderiv ℝ φ) x v)
    have := norm_nonneg v
    calc ‖(A + B.comp (fderiv ℝ (fderiv ℝ φ) x)) v‖ ≤ ‖A v‖ + ‖B (fderiv ℝ (fderiv ℝ φ) x v)‖ := by
          simpa using norm_add_le _ _
      _ ≤ ‖v‖ + ‖fderiv ℝ (fderiv ℝ φ) x‖ * ‖v‖ := add_le_add e3 (e4.trans e1)
      _ ≤ (1 + c₀) * ‖v‖ := by nlinarith
  · refine ContinuousLinearMap.opNorm_le_bound _ hK fun v => ?_
    refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun w => ?_
    rw [hd2 v w]
    have e1 := pc_D3_apply_le φ x v w
    have := norm_nonneg v
    have := norm_nonneg w
    refine (hBn _).trans (e1.trans ?_)
    calc _ ≤ K * ‖v‖ * ‖w‖ := by gcongr
      _ = _ := by ring

theorem pc_g_bounds {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (Ul Uj : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (yl yj : V) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    {x : EuclideanSpace ℝ (Fin n)} (hφ : ContDiffAt ℝ ∞ φ x) {c₀ K : ℝ} (h2 : ‖jet 2 φ x‖ ≤ c₀)
    (h3 : ‖jet 3 φ x‖ ≤ K) :
    ContDiffAt ℝ 2 (fun y => horC Ul yl (graphX J Uj yj φ y)) x ∧
      ‖fderiv ℝ (fun y => horC Ul yl (graphX J Uj yj φ y)) x‖ ≤ 1 + c₀ ∧
      ‖fderiv ℝ (fderiv ℝ (fun y => horC Ul yl (graphX J Uj yj φ y))) x‖ ≤ K := by
  set Θ : StrongDual ℝ (EuclideanSpace ℝ (Fin n)) →L[ℝ] EuclideanSpace ℝ (Fin n) :=
    (InnerProductSpace.toDual ℝ (EuclideanSpace ℝ (Fin n))).symm.toLinearIsometry.toContinuousLinearMap
    with hΘ
  have hΘn : ∀ z, ‖Θ z‖ = ‖z‖ := fun z => LinearIsometryEquiv.norm_map _ _
  have hJn : ∀ v : V, ‖J v‖ = ‖v‖ := by
    intro v
    have : ‖J v‖ ^ 2 = ‖v‖ ^ 2 := by
      rw [← real_inner_self_eq_norm_sq, ← real_inner_self_eq_norm_sq, hJi]
    exact (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp this
  set A : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n) :=
    (horCL Ul).comp Uj.toContinuousLinearMap with hA
  set B : StrongDual ℝ (EuclideanSpace ℝ (Fin n)) →L[ℝ] EuclideanSpace ℝ (Fin n) :=
    (horCL Ul).comp (J.comp (Uj.toContinuousLinearMap.comp Θ)) with hB
  have hAn : ∀ v, ‖A v‖ ≤ ‖v‖ := fun v => by
    simp only [hA, ContinuousLinearMap.comp_apply]
    exact (pc_horCL_le Ul _).trans (by simp)
  have hBn : ∀ z, ‖B z‖ ≤ ‖z‖ := fun z => by
    simp only [hB, ContinuousLinearMap.comp_apply]
    refine (pc_horCL_le Ul _).trans ?_
    rw [hJn]
    simp [hΘn]
  have hg : (fun y => horC Ul yl (graphX J Uj yj φ y)) =
      fun y => horCL Ul (yj - yl) + A y + B (fderiv ℝ φ y) := by
    funext y
    show horCL Ul (yj + Uj y + J (Uj (gradE φ y)) - yl) = _
    rw [pc_gradE_eq]
    simp only [hA, hB, hΘ, ContinuousLinearMap.comp_apply, map_add, map_sub]
    simp only [LinearIsometry.coe_toContinuousLinearMap, LinearIsometryEquiv.coe_toLinearIsometry]
    abel
  have hf2 : ContDiffAt ℝ 2 (fderiv ℝ φ) x := pc_cd2 (hφ.fderiv_right (m := ∞) (by simp))
  rw [hg]
  exact pc_g_aff_bounds φ A B (horCL Ul (yj - yl)) hf2 hAn hBn
    (((norm_iteratedFDeriv_le_norm_jet 2 φ x)).trans h2)
    (((norm_iteratedFDeriv_le_norm_jet 3 φ x)).trans h3)

end PartC

/-! ### Part D: the cutoffs `η_l` and their local forms -/

section PartD

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
  {J : V →L[ℝ] V} {F : M → V}

theorem pc_bump_bd1 {X F : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] [NormedAddCommGroup F]
    [NormedSpace ℝ F] {f : X → F}
    (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) : ∃ B₁ : ℝ, ∀ x, ‖fderiv ℝ f x‖ ≤ B₁ := by
  have h1 : Continuous (fderiv ℝ f) := hf.continuous_fderiv (by simp)
  have hc1 : HasCompactSupport (fderiv ℝ f) := hc.fderiv (𝕜 := ℝ)
  exact h1.bounded_above_of_compact_support hc1

theorem pc_bump_bd2 {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] {f : X → ℝ}
    (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) :
    ∃ B₂ : ℝ, ∀ x, ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ B₂ := by
  have hd : ContDiff ℝ ∞ (fderiv ℝ f) := hf.fderiv_right (m := ∞) (by simp)
  have hc1 : HasCompactSupport (fderiv ℝ f) := hc.fderiv (𝕜 := ℝ)
  exact pc_bump_bd1 hd hc1

theorem pc_bump_bd {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] {f : X → ℝ}
    (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) :
    ∃ B₁ B₂ : ℝ, (∀ x, ‖fderiv ℝ f x‖ ≤ B₁) ∧ (∀ x, ‖fderiv ℝ (fderiv ℝ f) x‖ ≤ B₂) := by
  obtain ⟨B₁, hB₁⟩ := pc_bump_bd1 hf hc
  obtain ⟨B₂, hB₂⟩ := pc_bump_bd2 hf hc
  exact ⟨B₁, B₂, hB₁, hB₂⟩

theorem pc_eta (n : ℕ) : ∃ η : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ η ∧ (∀ x, 0 ≤ η x) ∧
    (∀ x, η x ≤ 1) ∧ (∀ x, ‖x‖ ≤ 1 → η x = 1) ∧ (∀ x, 3 / 2 ≤ ‖x‖ → η x = 0) ∧
    ∃ B₁ B₂ : ℝ, (∀ x, ‖fderiv ℝ η x‖ ≤ B₁) ∧ (∀ x, ‖fderiv ℝ (fderiv ℝ η) x‖ ≤ B₂) := by
  let b : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨1, 3 / 2, by norm_num, by norm_num⟩
  have hη : ContDiff ℝ ∞ (b : EuclideanSpace ℝ (Fin n) → ℝ) := b.contDiff
  have hc : HasCompactSupport (b : EuclideanSpace ℝ (Fin n) → ℝ) := b.hasCompactSupport
  obtain ⟨B₁, B₂, hB₁, hB₂⟩ := pc_bump_bd hη hc
  refine ⟨b, hη, b.nonneg', fun _ => b.le_one, fun x hx => ?_, fun x hx => ?_, B₁, B₂, hB₁, hB₂⟩
  · exact b.one_of_mem_closedBall (by simpa [b] using hx)
  · exact b.zero_of_le_dist (by simpa [b] using hx)

/-! The cutoff `cvChi` with a bump supported in the ball of radius `3/2`. -/

theorem pc_chi_vanish {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη0 : ∀ x, 3 / 2 ≤ ‖x‖ → η x = 0)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    ∀ p, p ∉ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) (3 / 2) → cvChi η U y₀ σ F p = 0 := by
  intro p hp
  unfold cvChi
  split_ifs with h
  · obtain ⟨x, hx, rfl⟩ := h
    rw [hhor x hx]
    refine hη0 x ?_
    by_contra hlt
    exact hp ⟨x, by simpa using not_le.1 hlt, rfl⟩
  · rfl

theorem pc_chi_tsupport [T2Space M] {η : EuclideanSpace ℝ (Fin n) → ℝ}
    (hη0 : ∀ x, 3 / 2 ≤ ‖x‖ → η x = 0)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball 0 6) σ)
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    tsupport (cvChi η U y₀ σ F) ⊆ σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2) := by
  have hK : IsCompact (σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2)) :=
    (isCompact_closedBall _ _).image_of_continuousOn
      (hσ.cont.mono (closedBall_subset_ball (by norm_num)))
  refine closure_minimal (fun p hp => ?_) hK.isClosed
  by_contra hn
  exact hp (pc_chi_vanish hη0 hhor p (fun h => hn (image_mono ball_subset_closedBall h)))

theorem pc_chi_contMDiff [T2Space M] {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη : ContDiff ℝ ∞ η)
    (hη0 : ∀ x, 3 / 2 ≤ ‖x‖ → η x = 0) (hFs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, V) ∞ F)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball 0 6) σ)
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (cvChi η U y₀ σ F) := by
  have hg : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun p => η (horC U y₀ (F p))) :=
    (hη.comp (cv_contDiff_horC U y₀)).contMDiff.comp hFs
  have hK : IsCompact (σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2)) :=
    (isCompact_closedBall _ _).image_of_continuousOn
      (hσ.cont.mono (closedBall_subset_ball (by norm_num)))
  have hvan := pc_chi_vanish (η := η) hη0 (U := U) (y₀ := y₀) (F := F) hhor
  intro p
  by_cases hp : p ∈ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) 6
  · refine (hg p).congr_of_eventuallyEq ?_
    filter_upwards [hσ.isOpen_image.mem_nhds hp] with q hq
    simp [cvChi, hq]
  · have : cvChi η U y₀ σ F =ᶠ[𝓝 p] fun _ => 0 := by
      filter_upwards [hK.isClosed.isOpen_compl.mem_nhds (fun hm => hp
        (image_mono (closedBall_subset_ball (by norm_num : (3 / 2 : ℝ) < 6)) hm))]
        with q hq using hvan q (fun h => hq (image_mono ball_subset_closedBall h))
    exact (contMDiffAt_const).congr_of_eventuallyEq this

theorem pc_gradE_smooth {φ : EuclideanSpace ℝ (Fin n) → ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (hφ : ContDiffAt ℝ ∞ φ x) : ContDiffAt ℝ ∞ (gradE φ) x := by
  have : gradE φ = fun y => (InnerProductSpace.toDual ℝ (EuclideanSpace ℝ (Fin n))).symm
      (fderiv ℝ φ y) := funext (pc_gradE_eq φ)
  rw [this]
  exact (LinearIsometryEquiv.contDiff _).contDiffAt.comp x (hφ.fderiv_right (m := ∞) (by simp))

theorem pc_g_smooth (Ul Uj : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (yl yj : V)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} {x : EuclideanSpace ℝ (Fin n)} (hφ : ContDiffAt ℝ ∞ φ x) :
    ContDiffAt ℝ ∞ (fun y => horC Ul yl (graphX J Uj yj φ y)) x := by
  have h1 : ContDiffAt ℝ ∞ (fun y => graphX J Uj yj φ y) x := by
    unfold graphX
    exact (contDiffAt_const.add Uj.contDiff.contDiffAt).add
      (J.contDiff.contDiffAt.comp x (Uj.contDiff.contDiffAt.comp x (pc_gradE_smooth hφ)))
  exact (cv_contDiff_horC Ul yl).contDiffAt.comp x h1

/-- Near a point of `ball 0 6`, `η_l ∘ σ_j` is either `η ∘ g_l` or identically `0`. -/
theorem pc_local [T2Space M] {η : EuclideanSpace ℝ (Fin n) → ℝ} (hηs : ∀ x, 3 / 2 ≤ ‖x‖ → η x = 0)
    {c₀ : ℝ} {pj pl : M} {Uj Ul : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}
    {σj σl : EuclideanSpace ℝ (Fin n) → M} {φj φl : EuclideanSpace ℝ (Fin n) → ℝ}
    (hj : IsGraphChart ι J F c₀ pj Uj σj φj) (hl : IsGraphChart ι J F c₀ pl Ul σl φl)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :
    (fun y => cvChi η Ul (F pl) σl F (σj y)) =ᶠ[𝓝 x]
        (fun y => η (horC Ul (F pl) (graphX J Uj (F pj) φj y))) ∨
      (fun y => cvChi η Ul (F pl) σl F (σj y)) =ᶠ[𝓝 x] fun _ => 0 := by
  have hK : IsCompact (σl '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2)) :=
    (isCompact_closedBall _ _).image_of_continuousOn
      (hl.param.cont.mono (closedBall_subset_ball (by norm_num)))
  have hhor := fun z hz => (hl.graph z hz).2.1
  by_cases h : σj x ∈ σl '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (3 / 2)
  · left
    have hW : IsOpen (ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩ σj ⁻¹' (σl '' ball 0 6)) :=
      hj.param.cont.isOpen_inter_preimage hj.param.isOpen hl.param.isOpen_image
    have hxW : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩ σj ⁻¹' (σl '' ball 0 6) :=
      ⟨hx, image_mono (closedBall_subset_ball (by norm_num)) h⟩
    filter_upwards [hW.mem_nhds hxW] with y hy
    simp only [cvChi]
    rw [if_pos (show σj y ∈ σl '' ball (0 : EuclideanSpace ℝ (Fin n)) 6 from hy.2),
      (hj.graph y hy.1).1]
  · right
    have hW : IsOpen (ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩
        σj ⁻¹' (σl '' closedBall 0 (3 / 2))ᶜ) :=
      hj.param.cont.isOpen_inter_preimage hj.param.isOpen hK.isClosed.isOpen_compl
    filter_upwards [hW.mem_nhds ⟨hx, h⟩] with y hy
    exact pc_chi_vanish hηs hhor (σj y)
      (fun hm => hy.2 (image_mono ball_subset_closedBall hm))

theorem pc_f_bd [T2Space M] (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) {η : EuclideanSpace ℝ (Fin n) → ℝ}
    (hη : ContDiff ℝ ∞ η) (hη0 : ∀ x, 0 ≤ η x) (hη1 : ∀ x, η x ≤ 1)
    (hηs : ∀ x, 3 / 2 ≤ ‖x‖ → η x = 0) {B₁ B₂ : ℝ} (hB1 : ∀ x, ‖fderiv ℝ η x‖ ≤ B₁)
    (hB2 : ∀ x, ‖fderiv ℝ (fderiv ℝ η) x‖ ≤ B₂) {c₀ K : ℝ} (hc₀ : 0 ≤ c₀) (hK : 0 ≤ K)
    {pj pl : M} {Uj Ul : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}
    {σj σl : EuclideanSpace ℝ (Fin n) → M} {φj φl : EuclideanSpace ℝ (Fin n) → ℝ}
    (hj : IsGraphChart ι J F c₀ pj Uj σj φj) (hl : IsGraphChart ι J F c₀ pl Ul σl φl)
    (h3 : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φj x‖ ≤ K)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :
    PcBd (fun y => cvChi η Ul (F pl) σl F (σj y)) x 1 (B₁ * (1 + c₀))
        (B₂ * (1 + c₀) ^ 2 + B₁ * K) ∧
      ContDiffAt ℝ ∞ (fun y => cvChi η Ul (F pl) σl F (σj y)) x := by
  have hφat : ContDiffAt ℝ ∞ φj x := hj.smooth.contDiffAt (isOpen_ball.mem_nhds hx)
  have hB1' : 0 ≤ B₁ := (norm_nonneg _).trans (hB1 0)
  have hB2' : 0 ≤ B₂ := le_trans (ContinuousLinearMap.opNorm_nonneg _) (hB2 0)
  rcases pc_local ι hηs hj hl hx with h | h
  · obtain ⟨hgc, hg1, hg2⟩ := pc_g_bounds hJi Ul Uj (F pl) (F pj) hφat (hj.graph x hx).2.2.1
      (h3 x hx)
    have hbη : PcBd η (horC Ul (F pl) (graphX J Uj (F pj) φj x)) 1 B₁ B₂ :=
      ⟨pc_cd2 hη.contDiffAt, by rw [abs_of_nonneg (hη0 _)]; exact hη1 _, hB1 _, hB2 _⟩
    refine ⟨(hbη.comp_vec (g := fun y => horC Ul (F pl) (graphX J Uj (F pj) φj y)) hgc hg1 hg2).congr h.symm, ?_⟩
    exact ((hη.contDiffAt).comp x (pc_g_smooth Ul Uj (F pl) (F pj) hφat)).congr_of_eventuallyEq h
  · refine ⟨(pc_bd_zero x).congr h.symm |>.mono (by norm_num) (by positivity) (by positivity), ?_⟩
    exact contDiffAt_const.congr_of_eventuallyEq h

end PartD

/-! ### Part E: the partition of unity -/

section PartE

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

open Classical in
theorem exists_partition_cutoff [T2Space M] [CompactSpace M] {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (N₀ : ℕ) {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1) :
    ∃ Λ : ℝ, ∀ {F : M → V}, IsLagImm ι J F → ∀ (N : ℕ), N ≤ N₀ →
      ∀ (p : Fin N → M) (U : Fin N → (EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V))
        (σ : Fin N → EuclideanSpace ℝ (Fin n) → M) (φ : Fin N → EuclideanSpace ℝ (Fin n) → ℝ),
      (∀ j, IsGraphChart ι J F c₀ (p j) (U j) (σ j) (φ j)) →
      (∀ j, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 (φ j) x‖ ≤ K₃) →
      (∀ q : M, ∃ j, q ∈ σ j '' ball 0 1) →
      ∃ ρ : Fin N → M → ℝ, (∀ j, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (ρ j)) ∧ (∀ q, ∑ j, ρ j q ^ 2 = 1) ∧
        (∀ j, tsupport (ρ j) ⊆ σ j '' ball 0 2) ∧
        ∀ j, IsCutoff (fun x => if x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 then ρ j (σ j x) else 0)
          Λ := by
  obtain ⟨η, hη, hη0, hη1, hηone, hηs, B₁, B₂, hB1, hB2⟩ := pc_eta n
  refine ⟨pcLam N₀ B₁ B₂ (B₁ * (1 + c₀)) (B₂ * (1 + c₀) ^ 2 + B₁ * max K₃ 0), ?_⟩
  intro F hF N hN p U σ φ hch hjet hcov
  have hK : 0 ≤ max K₃ 0 := le_max_right _ _
  have h3 : ∀ j, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 (φ j) x‖ ≤ max K₃ 0 :=
    fun j x hx => (hjet j x hx).trans (le_max_left _ _)
  have hhor : ∀ l, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      horC (U l) (F (p l)) (F (σ l x)) = x := fun l x hx => ((hch l).graph x hx).2.1
  set χ : Fin N → M → ℝ := fun l => cvChi η (U l) (F (p l)) (σ l) F with hχ
  have hχs : ∀ l, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (χ l) := fun l =>
    pc_chi_contMDiff ι hη hηs (cv_contMDiff ι hF) (hch l).param (hhor l)
  have hχσ : ∀ l, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, χ l (σ l x) = η x := by
    intro l x hx
    simp only [hχ, cvChi]
    rw [if_pos (mem_image_of_mem _ hx), hhor l x hx]
  have hS1 : ∀ q, 1 ≤ ∑ l, χ l q ^ 2 := by
    intro q
    obtain ⟨l, x, hx, rfl⟩ := hcov q
    have hx6 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := ball_subset_ball (by norm_num) hx
    have h1 : χ l (σ l x) = 1 := by
      rw [hχσ l x hx6]
      exact hηone x (by simpa using hx.out.le)
    calc (1 : ℝ) = χ l (σ l x) ^ 2 := by rw [h1]; norm_num
      _ ≤ ∑ l', χ l' (σ l x) ^ 2 :=
        Finset.single_le_sum (f := fun l' => χ l' (σ l x) ^ 2) (fun i _ => sq_nonneg _)
          (Finset.mem_univ l)
  have hSsm : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun q => ∑ l, χ l q ^ 2) := by
    refine cv_contMDiff_sum Finset.univ (fun l q => χ l q ^ 2) fun l _ => ?_
    have : (fun q => χ l q ^ 2) = fun q => χ l q * χ l q := funext fun q => sq _
    rw [this]
    exact (hχs l).mul (hχs l)
  have hTsm : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun q => (∑ l, χ l q ^ 2) ^ (-(1 / 2 : ℝ))) := by
    intro q
    have hne : (∑ l, χ l q ^ 2) ≠ 0 := by have := hS1 q; linarith
    have hc : ContDiffAt ℝ ∞ (fun s : ℝ => s ^ (-(1 / 2 : ℝ))) (∑ l, χ l q ^ 2) :=
      Real.contDiffAt_rpow_const_of_ne hne
    exact hc.contMDiffAt.comp q (hSsm q)
  have hT2 : ∀ q, ((∑ l, χ l q ^ 2) ^ (-(1 / 2 : ℝ))) ^ 2 = (∑ l, χ l q ^ 2)⁻¹ := by
    intro q
    have h0 : 0 ≤ ∑ l, χ l q ^ 2 := by have := hS1 q; linarith
    rw [← Real.rpow_natCast, ← Real.rpow_mul h0]
    norm_num [Real.rpow_neg_one]
  refine ⟨fun l q => χ l q * (∑ l', χ l' q ^ 2) ^ (-(1 / 2 : ℝ)), fun l => (hχs l).mul hTsm,
    fun q => ?_, fun l => ?_, fun j => ?_⟩
  · have hne : (∑ l, χ l q ^ 2) ≠ 0 := by have := hS1 q; linarith
    calc ∑ l, (χ l q * (∑ l', χ l' q ^ 2) ^ (-(1 / 2 : ℝ))) ^ 2
        = (∑ l, χ l q ^ 2) * ((∑ l', χ l' q ^ 2) ^ (-(1 / 2 : ℝ))) ^ 2 := by
          rw [Finset.sum_mul]
          exact Finset.sum_congr rfl fun l _ => mul_pow _ _ _
      _ = 1 := by rw [hT2 q]; exact mul_inv_cancel₀ hne
  · refine (closure_mono (fun q hq => ?_)).trans
      ((pc_chi_tsupport ι hηs (hch l).param (hhor l)).trans
        (image_mono (closedBall_subset_ball (by norm_num))))
    exact left_ne_zero_of_mul (Function.mem_support.1 hq)
  · exact pc_abstract hη hη0 hη1 hηs hB1 hB2 hN (fun l x => χ l (σ j x))
      (fun l x hx => (pc_f_bd ι hJi hη hη0 hη1 hηs hB1 hB2 hc₀ hK (hch j) (hch l) (h3 j) hx).2.contDiffWithinAt)
      (fun l x hx => (pc_f_bd ι hJi hη hη0 hη1 hηs hB1 hB2 hc₀ hK (hch j) (hch l) (h3 j) hx).1)
      j (hχσ j) (fun x _ => hS1 (σ j x))

/-- info: 'CWFlow.Extension.exists_partition_cutoff' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_partition_cutoff

end PartE

end CWFlow.Extension
