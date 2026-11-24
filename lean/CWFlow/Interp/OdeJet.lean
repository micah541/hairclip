import CWFlow.Interp.Jets
import CWFlow.Extension.JetTime
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.Calculus.Deriv.Prod
import Mathlib.Analysis.Calculus.FDeriv.WithLp
import Mathlib.Analysis.Calculus.ContDiff.WithLp

/-!
# Spatial jets of a solution of a time-dependent ODE satisfy a closed ODE system

If `∂ₜ χ = V (t, χ)` with `V` smooth, then the jets `(χ, D_yχ, …, D_y^Jχ)` satisfy `∂ₜ Z_J = Θ_J (t, Z_J)`
with `Θ_J` smooth (`exists_ode_jet`), and `Θ_J` is affine in the top jet with coefficient `D_xV (t, χ)`.

`Θ` is built order by order (`ojθ`): `θ 0 (t, z) = V (t, z₀)` and
`θ (i+1) (t, z) (a :: I) = D_zθ_i (t, z_{≤ i}) [∂ₐ-shift of z] I`.

Encoding note: `z 0 : PiLp 2 (fun _ : Fin 0 → Fin n => E)` is not literally an element of `E`; we use
the coordinate `oj0 z : E` (the unique value of the order-zero jet) instead.
-/

set_option linter.unusedSectionVars false

open Set
open scoped ContDiff Topology

namespace CWFlow.Interp

/-- Jets of order `0, …, J` of an `E`-valued map. -/
abbrev OJ (n J : ℕ) : Type :=
  (i : Fin (J + 1)) → PiLp 2 (fun _ : Fin i → Fin n => EuclideanSpace ℝ (Fin n))

/-- The jets of `y ↦ χ (q.1, y)` at `q.2`. -/
noncomputable def ojet (n J : ℕ) (χ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (q : ℝ × EuclideanSpace ℝ (Fin n)) : OJ n J :=
  fun i => jet (i : ℕ) (fun y => χ (q.1, y)) q.2

/-- The order-zero component, as an element of `E`. -/
def oj0 {n J : ℕ} (z : OJ n J) : EuclideanSpace ℝ (Fin n) :=
  z ⟨0, Nat.succ_pos J⟩ (fun k => Fin.elim0 k)

/-- Truncation of a jet to lower order. -/
def ojTruncLe {n i J : ℕ} (h : i ≤ J) (z : OJ n J) : OJ n i :=
  fun k => z ⟨k.1, Nat.lt_of_lt_of_le k.2 (Nat.succ_le_succ h)⟩

/-- One-step truncation. -/
def ojTrunc {n J : ℕ} (z : OJ n (J + 1)) : OJ n J := ojTruncLe (Nat.le_succ J) z

/-- The `a`-slice of the next jet: `(shift_a z)_k K = z_{k+1} (a :: K)`. -/
noncomputable def ojShift {n J : ℕ} (a : Fin n) (z : OJ n (J + 1)) : OJ n J :=
  fun k => WithLp.toLp 2 (fun K : Fin k.1 → Fin n =>
    z ⟨k.1 + 1, Nat.succ_lt_succ k.2⟩ (Fin.cons a K : Fin (k.1 + 1) → Fin n))

/-- The top component. -/
def ojTop {n i : ℕ} (u : OJ n i) : PiLp 2 (fun _ : Fin i → Fin n => EuclideanSpace ℝ (Fin n)) :=
  u (Fin.last i)

section Basic

variable {n : ℕ}

local notation "E" => EuclideanSpace ℝ (Fin n)

theorem oj_oj0_ojet (J : ℕ) (χ : ℝ × E → E) (q : ℝ × E) : oj0 (ojet n J χ q) = χ q := by
  simp [oj0, ojet, jet_apply]

theorem oj_ojTruncLe_ojet {i J : ℕ} (h : i ≤ J) (χ : ℝ × E → E) (q : ℝ × E) :
    ojTruncLe h (ojet n J χ q) = ojet n i χ q := rfl

theorem oj_ojTrunc_ojet (J : ℕ) (χ : ℝ × E → E) (q : ℝ × E) :
    ojTrunc (ojet n (J + 1) χ q) = ojet n J χ q := rfl

theorem oj_contDiff_oj0 (J : ℕ) : ContDiff ℝ ∞ (fun z : OJ n J => oj0 z) :=
  (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (0 : ℕ) → Fin n => E) (fun k => Fin.elim0 k)).contDiff.comp
    (ContinuousLinearMap.proj (R := ℝ)
      (φ := fun i : Fin (J + 1) => PiLp 2 (fun _ : Fin i → Fin n => E))
      ⟨0, Nat.succ_pos J⟩).contDiff

theorem oj_contDiff_ojTruncLe {i J : ℕ} (h : i ≤ J) :
    ContDiff ℝ ∞ (fun z : OJ n J => ojTruncLe h z) :=
  contDiff_pi.2 fun k =>
    (ContinuousLinearMap.proj (R := ℝ)
      (φ := fun i : Fin (J + 1) => PiLp 2 (fun _ : Fin i → Fin n => E))
      ⟨k.1, Nat.lt_of_lt_of_le k.2 (Nat.succ_le_succ h)⟩).contDiff

theorem oj_contDiff_ojShift {J : ℕ} (a : Fin n) :
    ContDiff ℝ ∞ (fun z : OJ n (J + 1) => ojShift a z) :=
  contDiff_pi.2 fun k => contDiff_piLp' (𝕜 := ℝ) 2 fun K =>
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k.1 + 1) → Fin n => E)
      (Fin.cons a K : Fin (k.1 + 1) → Fin n)).contDiff.comp
    (ContinuousLinearMap.proj (R := ℝ)
      (φ := fun i : Fin (J + 2) => PiLp 2 (fun _ : Fin i → Fin n => E))
      ⟨k.1 + 1, Nat.succ_lt_succ k.2⟩).contDiff

/-- The order-by-order construction of the jet vector field. -/
noncomputable def ojθ (V : ℝ × E → E) :
    (J : ℕ) → ℝ × OJ n J → PiLp 2 (fun _ : Fin J → Fin n => E)
  | 0, p => WithLp.toLp 2 (fun _ => V (p.1, oj0 p.2))
  | J + 1, p => WithLp.toLp 2 (fun A =>
      fderiv ℝ (ojθ V J) (p.1, ojTrunc p.2) (0, ojShift (A 0) p.2) (Fin.tail A))

theorem oj_isOpen_dom {S : Set (ℝ × E)} (hS : IsOpen S) (J : ℕ) :
    IsOpen {p : ℝ × OJ n J | (p.1, oj0 p.2) ∈ S} :=
  hS.preimage (continuous_fst.prodMk ((oj_contDiff_oj0 J).continuous.comp continuous_snd))

theorem oj_contDiff_dom (J : ℕ) :
    ContDiff ℝ ∞ (fun p : ℝ × OJ n J => ((p.1, oj0 p.2) : ℝ × E)) :=
  contDiff_fst.prodMk ((oj_contDiff_oj0 J).comp contDiff_snd)

theorem oj_contDiffOn_θ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) (J : ℕ) :
    ContDiffOn ℝ ∞ (ojθ V J) {p : ℝ × OJ n J | (p.1, oj0 p.2) ∈ S} := by
  induction J with
  | zero =>
    refine contDiffOn_piLp' (𝕜 := ℝ) 2 fun I => ?_
    show ContDiffOn ℝ ∞ (fun p : ℝ × OJ n 0 => V (p.1, oj0 p.2)) _
    exact hV.comp (oj_contDiff_dom 0).contDiffOn (fun p hp => hp)
  | succ J ih =>
    refine contDiffOn_piLp' (𝕜 := ℝ) 2 fun A => ?_
    have hd : ContDiffOn ℝ ∞ (fun y => fderiv ℝ (ojθ V J) y)
        {p : ℝ × OJ n J | (p.1, oj0 p.2) ∈ S} :=
      ih.fderiv_of_isOpen (oj_isOpen_dom hS J) (by simp)
    have hq : ContDiff ℝ ∞ (fun p : ℝ × OJ n (J + 1) => ((p.1, ojTrunc p.2) : ℝ × OJ n J)) :=
      contDiff_fst.prodMk ((oj_contDiff_ojTruncLe (Nat.le_succ J)).comp contDiff_snd)
    have hcomp := hd.comp hq.contDiffOn (fun p hp => hp)
    have hs : ContDiff ℝ ∞ (fun p : ℝ × OJ n (J + 1) => (((0 : ℝ), ojShift (A 0) p.2) : ℝ × OJ n J)) :=
      contDiff_const.prodMk ((oj_contDiff_ojShift (A 0)).comp contDiff_snd)
    have hcl := hcomp.clm_apply hs.contDiffOn
    exact (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin J → Fin n => E) (Fin.tail A)).contDiff.comp_contDiffOn hcl

/-- The continuous linear order-zero coordinate. -/
noncomputable def oj0L (n J : ℕ) : OJ n J →L[ℝ] EuclideanSpace ℝ (Fin n) :=
  (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (0 : ℕ) → Fin n => EuclideanSpace ℝ (Fin n))
    (fun k => Fin.elim0 k)).comp
    (ContinuousLinearMap.proj (R := ℝ)
      (φ := fun i : Fin (J + 1) => PiLp 2 (fun _ : Fin i → Fin n => EuclideanSpace ℝ (Fin n)))
      ⟨0, Nat.succ_pos J⟩)

theorem oj_oj0L_apply (J : ℕ) (z : OJ n J) : oj0L n J z = oj0 z := rfl

theorem oj_fderiv_piLp_apply {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    {ι : Type*} [Fintype ι] {f : X → PiLp 2 (fun _ : ι => E)} {x : X}
    (hf : DifferentiableAt ℝ f x) (v : X) (I : ι) :
    fderiv ℝ f x v I = fderiv ℝ (fun y => f y I) x v := by
  have h := (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => E) I).hasFDerivAt.comp x hf.hasFDerivAt
  have h2 : fderiv ℝ (fun y => f y I) x =
      (PiLp.proj (𝕜 := ℝ) 2 (fun _ : ι => E) I).comp (fderiv ℝ f x) := h.fderiv
  rw [h2]
  rfl

theorem oj_differentiableAt_θ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) (J : ℕ) {t : ℝ} {w : OJ n J} (hw : (t, oj0 w) ∈ S) :
    DifferentiableAt ℝ (ojθ V J) (t, w) :=
  ((oj_contDiffOn_θ hS hV J).differentiableOn (by simp)).differentiableAt
    ((oj_isOpen_dom hS J).mem_nhds hw)

/-- The statement that `θ_i` has `D_xV` as derivative in the top jet direction. -/
def OjT (S : Set (ℝ × E)) (V : ℝ × E → E) (i : ℕ) : Prop :=
  ∀ (t : ℝ) (w : OJ n i), (t, oj0 w) ∈ S → ∀ u : OJ n i, (∀ k : Fin (i + 1), (k : ℕ) < i → u k = 0) →
    ∀ I : Fin i → Fin n,
      fderiv ℝ (ojθ V i) (t, w) (0, u) I = fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u I)

theorem oj_ojTop_zero (u : OJ n 0) (I : Fin 0 → Fin n) : ojTop u I = oj0 u := by
  have : I = fun k => Fin.elim0 k := Subsingleton.elim _ _
  subst this
  rfl

theorem oj_T_zero {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) : OjT S V 0 := by
  intro t w hw u _ I
  have hθ := oj_differentiableAt_θ hS hV 0 hw
  rw [oj_fderiv_piLp_apply hθ]
  have hVd : DifferentiableAt ℝ V (t, oj0 w) :=
    (hV.differentiableOn (by simp)).differentiableAt (hS.mem_nhds hw)
  set L : ℝ × OJ n 0 →L[ℝ] ℝ × E :=
    (ContinuousLinearMap.fst ℝ ℝ (OJ n 0)).prod ((oj0L n 0).comp (ContinuousLinearMap.snd ℝ ℝ (OJ n 0)))
    with hL
  have h1 : HasFDerivAt (fun p : ℝ × OJ n 0 => ojθ V 0 p I)
      ((fderiv ℝ V (t, oj0 w)).comp L) (t, w) :=
    hVd.hasFDerivAt.comp (t, w) L.hasFDerivAt
  have h2 : HasFDerivAt (fun x : E => V (t, x))
      ((fderiv ℝ V (t, oj0 w)).comp ((0 : E →L[ℝ] ℝ).prod (ContinuousLinearMap.id ℝ E))) (oj0 w) :=
    hVd.hasFDerivAt.comp (oj0 w) ((hasFDerivAt_const t (oj0 w)).prodMk (hasFDerivAt_id (oj0 w)))
  rw [h1.fderiv, h2.fderiv, oj_ojTop_zero]
  rfl

theorem oj_A_succ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) {i : ℕ} (hT : OjT S V i) (t : ℝ) (w u : OJ n (i + 1))
    (hw : (t, oj0 w) ∈ S) (hu : ∀ k : Fin (i + 2), (k : ℕ) < i + 1 → u k = 0)
    (A : Fin (i + 1) → Fin n) :
    ojθ V (i + 1) (t, w + u) A - ojθ V (i + 1) (t, w) A =
      fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u A) := by
  have htr : ojTrunc (w + u) = ojTrunc w := by
    funext k
    show w ⟨k.1, _⟩ + u ⟨k.1, _⟩ = w ⟨k.1, _⟩
    rw [hu _ (by have := k.2; simp only; omega), add_zero]
  have hsh : ∀ a : Fin n, ojShift a (w + u) = ojShift a w + ojShift a u := by
    intro a; funext k; ext K; rfl
  have hu' : ∀ a : Fin n, ∀ k : Fin (i + 1), (k : ℕ) < i → ojShift a u k = 0 := by
    intro a k hk
    refine PiLp.ext fun K => ?_
    show u ⟨k.1 + 1, _⟩ (Fin.cons a K : Fin (k.1 + 1) → Fin n) = 0
    rw [hu _ (by simp only; omega)]
    rfl
  show fderiv ℝ (ojθ V i) (t, ojTrunc (w + u)) (0, ojShift (A 0) (w + u)) (Fin.tail A) -
      fderiv ℝ (ojθ V i) (t, ojTrunc w) (0, ojShift (A 0) w) (Fin.tail A) = _
  have hp : (((0 : ℝ), ojShift (A 0) w + ojShift (A 0) u) : ℝ × OJ n i) =
      (0, ojShift (A 0) w) + (0, ojShift (A 0) u) := by simp
  rw [htr, hsh, hp, map_add, PiLp.add_apply, add_sub_cancel_left,
    hT t (ojTrunc w) hw (ojShift (A 0) u) (hu' _)]
  show fderiv ℝ (fun x => V (t, x)) (oj0 w) (u (Fin.last (i + 1)) (Fin.cons (A 0) (Fin.tail A) : Fin (i + 1) → Fin n)) = _
  rw [Fin.cons_self_tail]
  rfl

theorem oj_T_succ {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) {i : ℕ} (hT : OjT S V i) : OjT S V (i + 1) := by
  intro t w hw u hu A
  have hθ : DifferentiableAt ℝ (ojθ V (i + 1)) (t, w + (0 : ℝ) • u) := by
    simpa using oj_differentiableAt_θ hS hV (i + 1) hw
  have hpath : HasDerivAt (fun s : ℝ => ((t, w + s • u) : ℝ × OJ n (i + 1))) ((0 : ℝ), u) 0 := by
    have := (hasDerivAt_const (0 : ℝ) t).prodMk (((hasDerivAt_id (0 : ℝ)).smul_const u).const_add w)
    simpa using this
  have h2 := (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (i + 1) → Fin n => E) A).hasFDerivAt.comp_hasDerivAt
    (0 : ℝ) (hθ.hasFDerivAt.comp_hasDerivAt (0 : ℝ) hpath)
  simp only [zero_smul, add_zero] at h2
  have h3 : HasDerivAt (fun s : ℝ => ojθ V (i + 1) (t, w + s • u) A)
      (fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u A)) 0 := by
    have h4 : HasDerivAt (fun s : ℝ => ojθ V (i + 1) (t, w) A +
        s • fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u A))
        (fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u A)) 0 := by
      simpa using ((hasDerivAt_id (0 : ℝ)).smul_const
        (fderiv ℝ (fun x => V (t, x)) (oj0 w) (ojTop u A))).const_add (ojθ V (i + 1) (t, w) A)
    refine h4.congr_of_eventuallyEq (Filter.Eventually.of_forall fun s => ?_)
    have := oj_A_succ hS hV hT t w (s • u) hw (fun k hk => by simp [hu k hk]) A
    rw [sub_eq_iff_eq_add'] at this
    simp only
    rw [this]
    congr 1
    show _ = s • _
    rw [← map_smul]
    rfl
  exact h2.unique h3

theorem oj_T_all {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) (i : ℕ) : OjT S V i := by
  induction i with
  | zero => exact oj_T_zero hS hV
  | succ i ih => exact oj_T_succ hS hV ih

/-- Local version of `fderiv_jet_apply`. -/
theorem oj_fderiv_jet_apply {ξ : E → E} {x : E} (hξ : ContDiffAt ℝ ∞ ξ x) (j : ℕ) (i : Fin n)
    (I : Fin j → Fin n) :
    (fderiv ℝ (jet j ξ) x (EuclideanSpace.single i 1)) I = jet (j + 1) ξ x (Fin.cons i I) := by
  have hd' : DifferentiableAt ℝ (iteratedFDeriv ℝ j ξ) x :=
    (hξ.iteratedFDeriv_right (infty_add_natCast_le j)).differentiableAt (by simp)
  have hd : DifferentiableAt ℝ (jet j ξ) x :=
    (differentiableAt_piLp (𝕜 := ℝ) 2).2 fun I =>
      (hd'.hasFDerivAt.continuousMultilinear_apply_const
        (fun k => EuclideanSpace.single (I k) (1 : ℝ))).differentiableAt
  rw [fderiv_piLp_apply hd]
  show fderiv ℝ (fun y => iteratedFDeriv ℝ j ξ y (fun k => EuclideanSpace.single (I k) (1 : ℝ))) x
    (EuclideanSpace.single i 1) = _
  rw [fderiv_continuousMultilinear_apply_const hd', jet_apply, iteratedFDeriv_succ_apply_left]
  rfl

theorem oj_dt_eq {a b : ℝ} {B : Set E} (hB : IsOpen B) {V : ℝ × E → E} {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1) :
    ∀ q ∈ Ioo a b ×ˢ B, CWFlow.Extension.dt χ q = V (q.1, χ q) := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  rintro ⟨t, y⟩ hq
  have hd : DifferentiableAt ℝ χ (t, y) :=
    (hχ.contDiffAt (hΩ.mem_nhds hq)).differentiableAt (by simp)
  have h1 : HasDerivAt (fun s => χ (s, y)) (fderiv ℝ χ (t, y) (1, 0)) t := by
    have hs : HasDerivAt (fun s : ℝ => (s, y)) ((1 : ℝ), (0 : E)) t :=
      (hasDerivAt_id t).prodMk (hasDerivAt_const t y)
    exact hd.hasFDerivAt.comp_hasDerivAt t hs
  exact h1.unique (hode (t, y) hq)

theorem oj_jet_dt {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) {a b : ℝ} {B : Set E} (hB : IsOpen B) {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hdt : ∀ q ∈ Ioo a b ×ˢ B, CWFlow.Extension.dt χ q = V (q.1, χ q)) (i : ℕ) :
    ∀ t ∈ Ioo a b, ∀ y ∈ B,
      jet i (fun y => CWFlow.Extension.dt χ (t, y)) y = ojθ V i (t, ojet n i χ (t, y)) := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  have hψ : ∀ t ∈ Ioo a b, ∀ y ∈ B,
      ContDiffAt ℝ ∞ (fun y => CWFlow.Extension.dt χ (t, y)) y := by
    intro t ht y hy
    have h1 : ContDiffAt ℝ ∞ (CWFlow.Extension.dt χ) (t, y) :=
      (CWFlow.Extension.contDiffOn_dt hΩ hχ).contDiffAt (hΩ.mem_nhds ⟨ht, hy⟩)
    exact ContDiffAt.comp y h1 (contDiffAt_const.prodMk contDiffAt_id)
  have hχt : ∀ t ∈ Ioo a b, ∀ y ∈ B, ContDiffAt ℝ ∞ (fun y => χ (t, y)) y := by
    intro t ht y hy
    have h1 : ContDiffAt ℝ ∞ χ (t, y) := hχ.contDiffAt (hΩ.mem_nhds ⟨ht, hy⟩)
    exact ContDiffAt.comp y h1 (contDiffAt_const.prodMk contDiffAt_id)
  induction i with
  | zero =>
    intro t ht y hy
    refine PiLp.ext fun I => ?_
    show iteratedFDeriv ℝ 0 _ y _ = V (t, oj0 (ojet n 0 χ (t, y)))
    rw [oj_oj0_ojet, iteratedFDeriv_zero_apply]
    exact hdt (t, y) ⟨ht, hy⟩
  | succ i ih =>
    intro t ht y hy
    refine PiLp.ext fun A => ?_
    obtain ⟨a', I, rfl⟩ : ∃ a' I, A = Fin.cons a' I :=
      ⟨A 0, Fin.tail A, (Fin.cons_self_tail A).symm⟩
    rw [← oj_fderiv_jet_apply (hψ t ht y hy) i a' I]
    have hev : jet i (fun y => CWFlow.Extension.dt χ (t, y)) =ᶠ[𝓝 y]
        fun y' => ojθ V i (t, ojet n i χ (t, y')) :=
      Filter.eventuallyEq_of_mem (hB.mem_nhds hy) fun y' hy' => ih t ht y' hy'
    rw [hev.fderiv_eq]
    have hθ : DifferentiableAt ℝ (ojθ V i) (t, ojet n i χ (t, y)) := by
      refine oj_differentiableAt_θ hS hV i ?_
      rw [oj_oj0_ojet]
      exact hmaps (t, y) ⟨ht, hy⟩
    have hD : HasFDerivAt (fun y' => ojet n i χ (t, y'))
        (ContinuousLinearMap.pi fun k : Fin (i + 1) =>
          fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y) y := by
      refine hasFDerivAt_pi'.2 fun k => ?_
      have hdk : DifferentiableAt ℝ (jet (k : ℕ) (fun y => χ (t, y))) y := by
        have hd' : DifferentiableAt ℝ (iteratedFDeriv ℝ (k : ℕ) (fun y => χ (t, y))) y :=
          ((hχt t ht y hy).iteratedFDeriv_right (infty_add_natCast_le k)).differentiableAt (by simp)
        exact (differentiableAt_piLp (𝕜 := ℝ) 2).2 fun I =>
          (hd'.hasFDerivAt.continuousMultilinear_apply_const
            (fun k => EuclideanSpace.single (I k) (1 : ℝ))).differentiableAt
      have := hdk.hasFDerivAt
      rw [ContinuousLinearMap.proj_pi]
      exact this
    have hg := (hasFDerivAt_const t y).prodMk hD
    have hc := hθ.hasFDerivAt.comp y hg
    have hc' : HasFDerivAt (fun y' => ojθ V i (t, ojet n i χ (t, y'))) _ y := hc
    rw [hc'.fderiv]
    have hDe : (ContinuousLinearMap.pi fun k : Fin (i + 1) =>
          fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y) (EuclideanSpace.single a' 1) =
        ojShift a' (ojet n (i + 1) χ (t, y)) := by
      funext k
      refine PiLp.ext fun K => ?_
      show fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y (EuclideanSpace.single a' 1) K = _
      rw [oj_fderiv_jet_apply (hχt t ht y hy)]
      rfl
    show (fderiv ℝ (ojθ V i) (t, ojet n i χ (t, y))
        (0, (ContinuousLinearMap.pi fun k : Fin (i + 1) =>
          fderiv ℝ (jet (k : ℕ) (fun y => χ (t, y))) y) (EuclideanSpace.single a' 1))) I = _
    rw [hDe]
    show _ = fderiv ℝ (ojθ V i) (t, ojTrunc (ojet n (i + 1) χ (t, y)))
      (0, ojShift ((Fin.cons a' I : Fin (i + 1) → Fin n) 0) (ojet n (i + 1) χ (t, y)))
      (Fin.tail (Fin.cons a' I : Fin (i + 1) → Fin n))
    rw [oj_ojTrunc_ojet]
    simp only [Fin.cons_zero, Fin.tail_cons]

/-- The jet vector field. -/
noncomputable def ojΘ (V : ℝ × E → E) (J : ℕ) (p : ℝ × OJ n J) : OJ n J :=
  fun k => ojθ V k (p.1, ojTruncLe (Nat.lt_succ_iff.1 k.2) p.2)

theorem oj_top_clause {S : Set (ℝ × E)} (hS : IsOpen S) {V : ℝ × E → E}
    (hV : ContDiffOn ℝ ∞ V S) (J : ℕ) (j : Fin (J + 1)) (t : ℝ) (z z' : OJ n J)
    (hz : (t, oj0 z) ∈ S) (h0 : oj0 z = oj0 z')
    (hlt : ∀ i : Fin (J + 1), (i : ℕ) < j → z i = z' i) (I : Fin j → Fin n) :
    ojΘ V J (t, z) j I - fderiv ℝ (fun x => V (t, x)) (oj0 z) (z j I) =
      ojΘ V J (t, z') j I - fderiv ℝ (fun x => V (t, x)) (oj0 z') (z' j I) := by
  obtain ⟨jv, hj⟩ := j
  cases jv with
  | zero =>
    have e1 : ∀ w : OJ n J, w ⟨0, hj⟩ I = oj0 w := by
      intro w
      have : I = fun k => Fin.elim0 k := Subsingleton.elim _ _
      subst this
      rfl
    show V (t, oj0 z) - _ = V (t, oj0 z') - _
    rw [e1 z, e1 z', ← h0]
  | succ i =>
    set w : OJ n (i + 1) := ojTruncLe (Nat.lt_succ_iff.1 hj) z with hw
    set w' : OJ n (i + 1) := ojTruncLe (Nat.lt_succ_iff.1 hj) z' with hw'
    have hu : ∀ k : Fin (i + 2), (k : ℕ) < i + 1 → (w - w') k = 0 := by
      intro k hk
      have := hlt ⟨k.1, Nat.lt_of_lt_of_le k.2 (Nat.succ_le_succ (Nat.lt_succ_iff.1 hj))⟩ hk
      show w k - w' k = 0
      rw [sub_eq_zero]
      exact this
    have hA := oj_A_succ hS hV (oj_T_all hS hV i) t w' (w - w')
      (by show (t, oj0 z') ∈ S; rw [← h0]; exact hz) hu I
    rw [show w' + (w - w') = w from add_sub_cancel _ _] at hA
    have hop : ojTop (w - w') I = z ⟨i + 1, hj⟩ I - z' ⟨i + 1, hj⟩ I := rfl
    rw [hop, map_sub] at hA
    show ojθ V (i + 1) (t, w) I - _ = ojθ V (i + 1) (t, w') I - _
    have e : oj0 w' = oj0 z := h0.symm
    rw [sub_eq_sub_iff_sub_eq_sub, ← h0]
    rw [e] at hA
    exact hA

/-- For a solution `χ` of `∂ₜχ = V(t, χ)`, the spatial jets satisfy a closed smooth ODE system, linear in
the top jet with coefficient `D_xV(t, χ)`. (Encoding: `oj0 z : E` is the order-zero jet.) -/
theorem exists_ode_jet {a b : ℝ} {B : Set E} (hB : IsOpen B)
    {S : Set (ℝ × E)} (hS : IsOpen S)
    {V : ℝ × E → E} (hV : ContDiffOn ℝ ∞ V S)
    {χ : ℝ × E → E}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ B))
    (hmaps : ∀ q ∈ Ioo a b ×ˢ B, (q.1, χ q) ∈ S)
    (hode : ∀ q ∈ Ioo a b ×ˢ B, HasDerivAt (fun s => χ (s, q.2)) (V (q.1, χ q)) q.1) :
    ∃ Θ : ∀ J : ℕ, ℝ × OJ n J → OJ n J,
      (∀ J, ContDiffOn ℝ ∞ (Θ J) {p : ℝ × OJ n J | (p.1, oj0 p.2) ∈ S}) ∧
      (∀ J, ∀ q ∈ Ioo a b ×ˢ B,
        HasDerivAt (fun s => ojet n J χ (s, q.2)) (Θ J (q.1, ojet n J χ q)) q.1) ∧
      (∀ J (j : Fin (J + 1)) (t : ℝ) (z z' : OJ n J), (t, oj0 z) ∈ S → oj0 z = oj0 z' →
        (∀ i : Fin (J + 1), (i : ℕ) < j → z i = z' i) →
        ∀ I : Fin j → Fin n,
          Θ J (t, z) j I - fderiv ℝ (fun x => V (t, x)) (oj0 z) (z j I) =
          Θ J (t, z') j I - fderiv ℝ (fun x => V (t, x)) (oj0 z') (z' j I)) := by
  have hΩ : IsOpen (Ioo a b ×ˢ B) := isOpen_Ioo.prod hB
  refine ⟨fun J => ojΘ V J, fun J => ?_, fun J q hq => ?_, fun J j t z z' hz h0 hlt I =>
    oj_top_clause hS hV J j t z z' hz h0 hlt I⟩
  · refine contDiffOn_pi.2 fun k => ?_
    have hq : ContDiff ℝ ∞ (fun p : ℝ × OJ n J =>
        ((p.1, ojTruncLe (Nat.lt_succ_iff.1 k.2) p.2) : ℝ × OJ n (k : ℕ))) :=
      contDiff_fst.prodMk ((oj_contDiff_ojTruncLe (Nat.lt_succ_iff.1 k.2)).comp contDiff_snd)
    exact (oj_contDiffOn_θ hS hV (k : ℕ)).comp hq.contDiffOn (fun p hp => hp)
  · obtain ⟨t, y⟩ := q
    obtain ⟨ht, hy⟩ := hq
    refine hasDerivAt_pi.2 fun k => ?_
    have h1 := CWFlow.Extension.hasDerivAt_jet_slice hΩ hχ (k : ℕ) (t := t) (x := y) ⟨ht, hy⟩
    rw [show (fun y' => CWFlow.Extension.dt χ (t, y')) = fun y' => CWFlow.Extension.dt χ (t, y') from rfl,
      oj_jet_dt hS hV hB hχ hmaps (oj_dt_eq hB hχ hode) k t ht y hy] at h1
    exact h1

end Basic

/-- info: 'CWFlow.Interp.exists_ode_jet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_ode_jet

end CWFlow.Interp
