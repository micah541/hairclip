import CWFlow.Scalar.JetRec
import CWFlow.Scalar.Principal
import ShortTime.Nonlinear.FlowsNaturalChart
import ShortTime.MM.PairingBasis

/-!
# The exact jet recursion in a chart of `(M, g)`

Piece B of the Chen–Warren ↔ Mantegazza–Martinazzi bridge. In the chart and trivialization at
`α ∈ M` (MM's `chartRd`), the chart readings `G i = chartRd g α (∇ⁱu)` of the covariant
derivatives of a function `u : M → ℝ` (a `(0,0)`-tensor) satisfy the first-order recursion
`D(G i) = slotIns (G (i+1)) + kerL (G i)` on the open patch `patchΩ α` (`fderiv_chartRd_eq`). Feeding
this into the exact Euclidean recursion `CWFlow.JetRec.iteratedFDeriv_eq_jetLam` gives the
coordinate 4-jet of `u ∘ φ_α⁻¹` as an explicit continuous linear function `coordJet g α y` of the
chart readings of the covariant 4-jet:

* `scal : TensorRSModel 0 0 ℝ E →L[ℝ] ℝ` is the scalar value of a model `(0,0)`-tensor; it reads
  MM's `lam` through the trivialization (`scal_trivAt`, `scal_chartRd`);
* `coordJet g α y v : SJet E` is the coordinate 4-jet built from the covariant jet values `v`,
  jointly smooth in `(y, v)` (`contDiffOn_coordJet`); its components of order `≤ 3` depend only
  on `v 0, …, v 3` (`coordJet_lower`), those of order `≤ 2` only on `v 0, …, v 2` (`coordJet_lower2`);
* `coordJet_eq_addTop`: the top value enters affinely, through `topRead`, with
  `topRead t m = t unitM m` (`topRead_apply`);
* `jetOf_scalRd`: for `u : SmoothCcTensor g 0 0`, the coordinate 4-jet of `y ↦ scal (chartRd g α u y)`
  at `y ∈ patchΩ α` is `coordJet g α y` of the chart readings of `(u, ∇u, …, ∇⁴u)`;
* `trivJet_covJet_eq`: these chart readings are MM's `trivJet α x (covJet g 0 0 4 u x)`.

Conventions: as in MM, the model fibre of `∇ᵏu` is `TensorRSModel 0 (0 + k) ℝ E` and `0 + k` is kept in
this form (it is not definitionally `k` for variable `k`).
-/

noncomputable section

open Set Filter
open Bundle Manifold
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open ShortTime ShortTime.HFlow
open CWFlow.JetRec
open scoped Manifold ContDiff Topology

set_option linter.unusedSectionVars false

namespace CWFlow.JetRec

section Generic

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {F : ℕ → Type*} [∀ i, NormedAddCommGroup (F i)] [∀ i, NormedSpace ℝ (F i)]
variable (P : ∀ i, F (i + 1) →L[ℝ] (E →L[ℝ] F i)) (Q : ∀ i, E → F i →L[ℝ] (E →L[ℝ] F i))

/-- The jet of values with the top value set to zero. -/
def zeroLast (n : ℕ) (v : Vals F n) : Vals F n := Function.update v (Fin.last n) 0

/-- A jet of values is the jet with top value `0` plus the jet whose only nonzero value is the top one. -/
theorem eq_zeroLast_add_snoc (n : ℕ) (v : Vals F n) :
    v = zeroLast n v + Fin.snoc (α := fun k : Fin (n + 1) => F k) 0 (v (Fin.last n)) := by
  funext k
  unfold zeroLast
  cases k using Fin.lastCases with
  | last =>
    rw [Pi.add_apply, Function.update_self, Fin.snoc_last]
    exact (zero_add _).symm
  | cast k =>
    rw [Pi.add_apply, Function.update_of_ne (Fin.castSucc_ne_last k), Fin.snoc_castSucc,
      Pi.zero_apply, add_zero]

/-- **The top value enters `jetLam` through `jetTop`**, for every order `n`. (Stated for a variable `n`: at a
numeral the recursion `jetLam _ _ 4` unfolds, which is very expensive for the unifier and the kernel.) -/
theorem jetLam_zeroLast (n : ℕ) (y : E) (v : Vals F n) :
    jetLam P Q n y v = jetLam P Q n y (zeroLast n v) + jetTop P n (v (Fin.last n)) := by
  rw [← jetLam_snoc_zero P Q n y, ← map_add]
  exact congrArg (jetLam P Q n y) (eq_zeroLast_add_snoc n v)

end Generic

end CWFlow.JetRec

namespace CWFlow.Bridge

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]

/-! ### The scalar value of a `(0,0)`-tensor -/

variable (E) in
/-- The unit model `0`-form. -/
def unitM : Tensor0SModel 0 ℝ E := ContinuousMultilinearMap.constOfIsEmpty ℝ (fun _ : Fin 0 => E) 1

theorem unitM_apply (m : Fin 0 → E) : unitM E m = 1 := rfl

variable (E) in
/-- The scalar value of a model `(0,0)`-tensor: its value at the unit `0`-form, evaluated at the
empty family of vectors. (`E` is explicit because `TensorRSModel 0 0 ℝ E →L[ℝ] ℝ` cannot be coerced
to a function before `E` is known.) -/
def scal : TensorRSModel 0 0 ℝ E →L[ℝ] ℝ :=
  (ContinuousMultilinearMap.apply ℝ (fun _ : Fin 0 => E) ℝ (fun i => i.elim0)).comp
    (ContinuousLinearMap.apply ℝ (Tensor0SModel 0 ℝ E) (unitM E))

theorem scal_apply (T : TensorRSModel 0 0 ℝ E) : scal E T = T (unitM E) (fun i => i.elim0) := rfl

theorem unitM_compContinuousLinearMap (f : Fin 0 → E →L[ℝ] E) :
    (unitM E).compContinuousLinearMap f = unitM E :=
  ContinuousMultilinearMap.ext fun _ => rfl

/-- The chart twist of a `(0,0)`-tensor does not change its scalar value. -/
theorem scal_chartRSTwistInv (α b : M) (T : TensorRSModel 0 0 ℝ E) :
    scal E (chartRSTwistInv (I := 𝓘(ℝ, E)) (M := M) α b 0 0 T) = scal E T := by
  rw [scal_apply, scal_apply, chartRSTwistInv_apply,
    ContinuousMultilinearMap.compContinuousLinearMap_apply, unitM_compContinuousLinearMap]
  exact congrArg _ (Subsingleton.elim _ _)

theorem toModel_unitZeroSec (x : M) :
    Tensor0SSpace.toModel (unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x) = unitM E := by
  rw [unitZeroSec_apply, Tensor0SSpace.toModel_ofModel]
  rfl

/-- **The trivialization reading of a `(0,0)`-tensor has scalar value `lam`.** -/
theorem scal_trivAt {α x : M} (hx : x ∈ (extChartAt 𝓘(ℝ, E) α).source)
    (T : TensorRSSpace 0 0 𝓘(ℝ, E) x) :
    scal E ((trivAt 𝓘(ℝ, E) M 0 0 α).continuousLinearMapAt ℝ x T) = lam (I := 𝓘(ℝ, E)) x T := by
  have hs : x ∈ (chartAt E α).source := by rwa [extChartAt_source] at hx
  rw [triv_continuousLinearMapAt_eq_chartRSTwistInv_toModel 0 0 α hs T, scal_chartRSTwistInv,
    scal_apply, ← toModel_unitZeroSec x, TensorRSSpace.toModel_apply_toModel]
  rfl

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

/-- The chart reading of a `(0,0)`-tensor field has scalar value `lam` at the corresponding point. -/
theorem scal_chartRd (α : M) (u : SmoothCcTensor g 0 0) {y : E}
    (hy : y ∈ (extChartAt 𝓘(ℝ, E) α).target) :
    scal E (chartRd g α u y) =
      lam (I := 𝓘(ℝ, E)) ((extChartAt 𝓘(ℝ, E) α).symm y)
        (u.toSection ((extChartAt 𝓘(ℝ, E) α).symm y)) := by
  rw [chartRd_eq_trivAt]
  exact scal_trivAt ((extChartAt 𝓘(ℝ, E) α).map_target hy) _

/-! ### The data of the recursion -/

variable (E) in
/-- The model fibre of `∇ᵏu` for a function `u`. `Vals (TF E) 4` is `JetModel E 0 0 4`. -/
abbrev TF (k : ℕ) : Type _ := TensorRSModel 0 (0 + k) ℝ E

variable (E) in
/-- The insertion into the first slot, `P i` of the recursion. -/
def cP (i : ℕ) : TF E (i + 1) →L[ℝ] (E →L[ℝ] TF E i) := slotIns E (0 + i)

/-- The Christoffel term `Q i` of the recursion in the chart at `α`. -/
def cQ (α : M) (i : ℕ) : E → TF E i →L[ℝ] (E →L[ℝ] TF E i) :=
  kerL g α (0 + i) (patchB (I := 𝓘(ℝ, E)) α)

theorem contDiffOn_cQ (α : M) (i : ℕ) : ContDiffOn ℝ ∞ (cQ g α i) (patchΩ (I := 𝓘(ℝ, E)) α) :=
  (contDiffOn_kerL α (0 + i) (patchB (I := 𝓘(ℝ, E)) α)).mono (patchΩ_subset_target α)

variable (E) in
/-- Postcomposition with `scal`: the scalar `n`-linear map of a `TF 0`-valued one. -/
def flatJ (n : ℕ) : JetT E (TF E) n →L[ℝ] (E [×n]→L[ℝ] ℝ) :=
  ContinuousLinearMap.compContinuousMultilinearMapL ℝ (fun _ : Fin n => E) (TF E 0) ℝ (scal E)

theorem flatJ_apply (n : ℕ) (J : JetT E (TF E) n) (m : Fin n → E) :
    flatJ E n J m = scal E (J m) := rfl

variable (E) in
/-- The values `(v 0, …, v n)` of a 4-jet of values. -/
def trunc (n : ℕ) (h : n ≤ 4) : Vals (TF E) 4 →L[ℝ] Vals (TF E) n :=
  ContinuousLinearMap.pi fun k : Fin (n + 1) =>
    ContinuousLinearMap.proj (R := ℝ) (φ := fun k : Fin 5 => TF E k) (Fin.castLE (by omega) k)

theorem trunc_apply (n : ℕ) (h : n ≤ 4) (v : Vals (TF E) 4) (k : Fin (n + 1)) :
    trunc E n h v k = v (Fin.castLE (by omega) k) := rfl

theorem trunc_self (v : Vals (TF E) 4) : trunc E 4 le_rfl v = v := rfl

theorem trunc_congr (n : ℕ) (h : n ≤ 4) {v v' : Vals (TF E) 4}
    (hv : ∀ k : Fin 5, (k : ℕ) ≤ n → v k = v' k) : trunc E n h v = trunc E n h v' := by
  funext k
  rw [trunc_apply, trunc_apply]
  exact hv _ (by simp only [Fin.val_castLE]; omega)

/-! ### Two abstract lemmas on jets

They are stated for abstract jets and components: applied to `coordJet`, they keep the kernel from comparing
projections of `coordJet` with components of explicit tuples. -/

section Abstract

variable {E' : Type*} [NormedAddCommGroup E'] [NormedSpace ℝ E']

theorem eq_addTop_of {z z' : SJet E'} {ζ : E' [×4]→L[ℝ] ℝ} (h1 : z.1 = z'.1) (h2 : z.2.1 = z'.2.1)
    (h3 : z.2.2.1 = z'.2.2.1) (h4 : z.2.2.2 = z'.2.2.2 + ζ) : z = addTop z' ζ := by
  obtain ⟨a, b, c, d⟩ := z
  obtain ⟨a', b', c', d'⟩ := z'
  simp only at h1 h2 h3 h4
  subst h1 h2 h3 h4
  rfl

theorem contDiffOn_sjet {X : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X] {s : Set X}
    {a : X → E' [×1]→L[ℝ] ℝ} {b : X → E' [×2]→L[ℝ] ℝ} {c : X → E' [×3]→L[ℝ] ℝ}
    {d : X → E' [×4]→L[ℝ] ℝ} (ha : ContDiffOn ℝ ∞ a s) (hb : ContDiffOn ℝ ∞ b s)
    (hc : ContDiffOn ℝ ∞ c s) (hd : ContDiffOn ℝ ∞ d s) :
    ContDiffOn ℝ ∞ (fun p => ((continuousMultilinearCurryFin1 ℝ E' ℝ (a p), b p, c p, d p) : SJet E')) s := by
  have he : ContDiff ℝ ∞ (continuousMultilinearCurryFin1 ℝ E' ℝ) :=
    LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
  exact (he.comp_contDiffOn ha).prodMk (hb.prodMk (hc.prodMk hd))

end Abstract

/-! ### The jet maps as plain definitions

`lamJ g α n y v` is `flatJ E n (jetLam (cP E) (cQ g α) n y v)` and `topJ n t` is `flatJ E n (jetTop (cP E) n t)`.

They are definitions, and not just notation, for a reason of cost. An application `⇑(jetLam P Q 4 y) v` of a
coerced map is a projection application, and both the unifier and the kernel unfold the map before they compare
the arguments. So comparing two such terms whose arguments differ slightly (`zeroTop v` and `Function.update …`,
`p.2` at two definitionally equal types, …) unfolds the recursion `jetLam _ _ 4` on both sides: this exhausted
the recursion depth in the elaborator and more than 7 GB in the kernel. With a constant at the head the arguments
are compared first. The same happens under `+` (also a projection), so the lemmas below are proved for a
variable order `n`, where nothing unfolds, and stated so that at `n = 4` they are *syntactically* the statements
needed (`zeroTop v` unfolds to `zeroLast 4 v`, `topRead t` to `topJ E 4 t`). -/

/-- The `n`-th coordinate derivative, as a function of the chart readings `v = (u, ∇u, …, ∇ⁿu)`. -/
def lamJ (α : M) (n : ℕ) (y : E) (v : Vals (TF E) n) : E [×n]→L[ℝ] ℝ :=
  flatJ E n (jetLam (cP E) (cQ g α) n y v)

variable (E) in
/-- The scalar `n`-linear map read off the chart reading of `∇ⁿu`. -/
def topJ (n : ℕ) (t : TF E n) : E [×n]→L[ℝ] ℝ := flatJ E n (jetTop (cP E) n t)

/-- **The top value enters `lamJ` affinely**, through `topJ`. -/
theorem lamJ_top (α : M) (n : ℕ) (y : E) (v : Vals (TF E) n) :
    lamJ g α n y v = lamJ g α n y (zeroLast n v) + topJ E n (v (Fin.last n)) := by
  unfold lamJ topJ
  rw [jetLam_zeroLast (cP E) (cQ g α) n y v, map_add]

theorem contDiffOn_lamJ (α : M) (n : ℕ) :
    ContDiffOn ℝ ∞ (fun p : E × Vals (TF E) n => lamJ g α n p.1 p.2)
      (patchΩ (I := 𝓘(ℝ, E)) α ×ˢ univ) :=
  (flatJ E n).contDiff.comp_contDiffOn
    (contDiffOn_jetLam_uncurry (cP E) (cQ g α) (isOpen_patchΩ (I := 𝓘(ℝ, E)) α) (contDiffOn_cQ g α) n)

/-! ### The coordinate jet -/

/-- **The coordinate 4-jet** `(Du, D²u, D³u, D⁴u)` at the chart point `y`, as a function of the chart
readings `v = (u, ∇u, …, ∇⁴u)` of the covariant 4-jet. -/
def coordJet (α : M) (y : E) (v : JetModel E 0 0 4) : SJet E :=
  (continuousMultilinearCurryFin1 ℝ E ℝ (lamJ g α 1 y (trunc E 1 (by norm_num) v)),
    lamJ g α 2 y (trunc E 2 (by norm_num) v),
    lamJ g α 3 y (trunc E 3 (by norm_num) v),
    lamJ g α 4 y v)

theorem coordJet_four (α : M) (y : E) (v : JetModel E 0 0 4) :
    (coordJet g α y v).2.2.2 = lamJ g α 4 y v := rfl

/-- The components of order `≤ 3` depend only on `v 0, …, v 3`. -/
theorem coordJet_lower (α : M) (y : E) {v v' : JetModel E 0 0 4}
    (h : ∀ k, k ≠ Fin.last 4 → v k = v' k) :
    (coordJet g α y v).1 = (coordJet g α y v').1 ∧
      (coordJet g α y v).2.1 = (coordJet g α y v').2.1 ∧
      (coordJet g α y v).2.2.1 = (coordJet g α y v').2.2.1 := by
  have h' : ∀ n (hn : n ≤ 3) (hn' : n ≤ 4), trunc E n hn' v = trunc E n hn' v' := fun n hn hn' =>
    trunc_congr n hn' fun k hk => h k (by
      intro hk'
      subst hk'
      simp at hk
      omega)
  refine ⟨?_, ?_, ?_⟩
  · simp only [coordJet, h' 1 (by norm_num)]
  · simp only [coordJet, h' 2 (by norm_num)]
  · simp only [coordJet, h' 3 (by norm_num)]

/-- The components of order `≤ 2` depend only on `v 0, v 1, v 2`. -/
theorem coordJet_lower2 (α : M) (y : E) {v v' : JetModel E 0 0 4}
    (h : ∀ k : Fin 5, (k : ℕ) ≤ 2 → v k = v' k) :
    (coordJet g α y v).1 = (coordJet g α y v').1 ∧
      (coordJet g α y v).2.1 = (coordJet g α y v').2.1 := by
  have h' : ∀ n (hn : n ≤ 2) (hn' : n ≤ 4), trunc E n hn' v = trunc E n hn' v' := fun n hn hn' =>
    trunc_congr n hn' fun k hk => h k (hk.trans hn)
  refine ⟨?_, ?_⟩
  · simp only [coordJet, h' 1 (by norm_num)]
  · simp only [coordJet, h' 2 (by norm_num)]

/-! ### The top-order term -/

/-- The scalar `4`-linear map read off the chart reading of `∇⁴u`: `topRead t m = t unitM m`. -/
def topRead (t : TensorRSModel 0 (0 + 4) ℝ E) : E [×4]→L[ℝ] ℝ := topJ E 4 t

theorem topRead_apply (t : TensorRSModel 0 (0 + 4) ℝ E) (m : Fin 4 → E) :
    topRead t m = t (unitM E) m := by
  change t (unitM E) (Fin.cons (m 0) (Fin.cons (m 1) (Fin.cons (m 2) (Fin.cons (m 3)
    (fun i : Fin 0 => i.elim0))))) = t (unitM E) m
  congr 1
  ext i
  fin_cases i <;> rfl

/-- The 4-jet of values with the top value set to zero. -/
def zeroTop (v : JetModel E 0 0 4) : JetModel E 0 0 4 := zeroLast (F := TF E) 4 v

theorem zeroTop_apply_of_ne (v : JetModel E 0 0 4) {k : Fin 5} (hk : k ≠ Fin.last 4) :
    zeroTop v k = v k := Function.update_of_ne hk _ _

/-- **The top value enters affinely**, through `topRead`. -/
theorem coordJet_eq_addTop (α : M) (y : E) (v : JetModel E 0 0 4) :
    coordJet g α y v = addTop (coordJet g α y (zeroTop v)) (topRead (v (Fin.last 4))) := by
  obtain ⟨h1, h2, h3⟩ := coordJet_lower g α y (v := v) (v' := zeroTop v)
    fun k hk => (zeroTop_apply_of_ne v hk).symm
  have h4 : (coordJet g α y v).2.2.2 =
      (coordJet g α y (zeroTop v)).2.2.2 + topRead (v (Fin.last 4)) := by
    rw [coordJet_four, coordJet_four, topRead, zeroTop]
    exact lamJ_top g α 4 y v
  -- through the abstract `eq_addTop_of`, so that the kernel never compares `(coordJet ..).i` with the
  -- components of an explicit tuple
  exact eq_addTop_of h1 h2 h3 h4

/-! ### Smoothness -/

/-- **Joint smoothness of the coordinate jet** in `(y, v)` over the patch. -/
theorem contDiffOn_coordJet (α : M) :
    ContDiffOn ℝ ∞ (fun p : E × JetModel E 0 0 4 => coordJet g α p.1 p.2)
      (patchΩ (I := 𝓘(ℝ, E)) α ×ˢ univ) := by
  have hj : ∀ n (h : n ≤ 4), ContDiffOn ℝ ∞ (fun p : E × JetModel E 0 0 4 =>
      lamJ g α n p.1 (trunc E n h p.2)) (patchΩ (I := 𝓘(ℝ, E)) α ×ˢ univ) := by
    intro n h
    have h2 : ContDiffOn ℝ ∞ (fun p : E × JetModel E 0 0 4 => (p.1, trunc E n h p.2))
        (patchΩ (I := 𝓘(ℝ, E)) α ×ˢ univ) :=
      contDiffOn_fst.prodMk ((trunc E n h).contDiff.comp_contDiffOn contDiffOn_snd)
    -- elaborated without the expected type: otherwise the unifier has to split the goal as `?g ∘ ?f`
    have h3 := (contDiffOn_lamJ g α n).comp h2 fun p hp => ⟨hp.1, mem_univ _⟩
    exact h3
  have h4 : ContDiffOn ℝ ∞ (fun p : E × JetModel E 0 0 4 => lamJ g α 4 p.1 p.2)
      (patchΩ (I := 𝓘(ℝ, E)) α ×ˢ univ) := contDiffOn_lamJ g α 4
  exact contDiffOn_sjet (hj 1 (by norm_num)) (hj 2 (by norm_num)) (hj 3 (by norm_num)) h4

/-! ### The chart reading of a function and its coordinate jet -/

/-- The chart reading of a function `u` (a `(0,0)`-tensor field): `y ↦ scal (chartRd g α u y)`, i.e.
`u ∘ φ_α⁻¹` on the chart target. -/
def scalRd (α : M) (u : SmoothCcTensor g 0 0) : E → ℝ := fun y => scal E (chartRd g α u y)

theorem contDiffOn_scalRd (α : M) (u : SmoothCcTensor g 0 0) :
    ContDiffOn ℝ ∞ (scalRd g α u) (extChartAt 𝓘(ℝ, E) α).target :=
  (scal E).contDiff.comp_contDiffOn (contDiffOn_chartRd α u)

/-- The chart readings of the covariant derivatives satisfy the recursion on the patch. -/
theorem fderiv_chartRd_iteratedCovGrad (α : M) (u : SmoothCcTensor g 0 0) (i : ℕ) {y : E}
    (hy : y ∈ patchΩ (I := 𝓘(ℝ, E)) α) :
    fderiv ℝ (chartRd g α (iteratedCovGrad g 0 0 i u)) y =
      cP E i (chartRd g α (iteratedCovGrad g 0 0 (i + 1) u) y) +
        cQ g α i y (chartRd g α (iteratedCovGrad g 0 0 i u) y) :=
  fderiv_chartRd_eq (patch_subset α) (patchB_spec α) _ (patchΩ_subset_closedBall α hy)

/-- **Coordinate derivatives of the chart reading of `u`** from the chart readings of
`(u, ∇u, …, ∇ⁿu)`. -/
theorem iteratedFDeriv_scalRd (α : M) (u : SmoothCcTensor g 0 0) (n : ℕ) {y : E}
    (hy : y ∈ patchΩ (I := 𝓘(ℝ, E)) α) :
    iteratedFDeriv ℝ n (scalRd g α u) y =
      lamJ g α n y (fun k : Fin (n + 1) => chartRd g α (iteratedCovGrad g 0 0 k u) y) := by
  unfold lamJ
  have hΩ := isOpen_patchΩ (I := 𝓘(ℝ, E)) α
  have hG : ∀ i, ContDiffOn ℝ ∞
      (fun y => chartRd g α (iteratedCovGrad g 0 0 i u) y : E → TF E i) (patchΩ (I := 𝓘(ℝ, E)) α) :=
    fun i => (contDiffOn_chartRd α _).mono (patchΩ_subset_target α)
  have h := iteratedFDeriv_eq_jetLam (cP E) (cQ g α) hΩ (contDiffOn_cQ g α) hG
    (fun i y hy => fderiv_chartRd_iteratedCovGrad g α u i hy) n y hy
  have hc : ContDiffAt ℝ ∞ (fun y => chartRd g α (iteratedCovGrad g 0 0 0 u) y : E → TF E 0) y :=
    (hG 0).contDiffAt (hΩ.mem_nhds hy)
  have h2 := ContinuousLinearMap.iteratedFDeriv_comp_left (scal E) hc (i := n)
    (by exact_mod_cast le_top)
  rw [h] at h2
  exact h2

/-- **Main theorem: the coordinate 4-jet of the chart reading of `u` is `coordJet` of the chart
readings of its covariant 4-jet.** -/
theorem jetOf_scalRd (α : M) (u : SmoothCcTensor g 0 0) {y : E}
    (hy : y ∈ patchΩ (I := 𝓘(ℝ, E)) α) :
    jetOf (scalRd g α u) y =
      coordJet g α y (fun k => chartRd g α (iteratedCovGrad g 0 0 k u) y) := by
  have h := fun n => iteratedFDeriv_scalRd g α u n hy
  refine Prod.ext ?_ (Prod.ext (h 2) (Prod.ext (h 3) (h 4)))
  have h1 : lamJ g α 1 y (trunc E 1 (by norm_num)
      (fun k => chartRd g α (iteratedCovGrad g 0 0 k u) y)) = iteratedFDeriv ℝ 1 (scalRd g α u) y :=
    (h 1).symm
  change fderiv ℝ (scalRd g α u) y = continuousMultilinearCurryFin1 ℝ E ℝ
    (lamJ g α 1 y (trunc E 1 (by norm_num)
      (fun k => chartRd g α (iteratedCovGrad g 0 0 k u) y)))
  rw [h1]
  ext a
  rw [continuousMultilinearCurryFin1_apply, iteratedFDeriv_one_apply]
  simp

/-- The chart readings of the covariant 4-jet are MM's `trivJet` of `covJet`. -/
theorem trivJet_covJet_eq (α : M) (u : SmoothCcTensor g 0 0) {x : M}
    (hx : x ∈ (extChartAt 𝓘(ℝ, E) α).source) :
    EquationData.trivJet (I := 𝓘(ℝ, E)) α x (covJet g 0 0 4 u x) =
      fun k : Fin (4 + 1) => chartRd g α (iteratedCovGrad g 0 0 k u) (extChartAt 𝓘(ℝ, E) α x) := by
  funext k
  rw [chartRd_eq_trivAt, (extChartAt 𝓘(ℝ, E) α).left_inv hx]
  rfl

end CWFlow.Bridge
