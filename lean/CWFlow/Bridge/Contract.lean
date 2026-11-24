import ShortTime.MM.Main
import ShortTime.MM.SmoothF
import CWFlow.Scalar.Principal
import CWFlow.Scalar.Invariance

/-!
# The MM contraction versus the double trace (piece C of the Chen–Warren bridge)

MM's principal term is the frame sum `frameSum 2 x Ef b b' ζ = −Σ E₀^{ij} E₁^{kl} ζ(e_i, e_j, e_k, e_l)`
(`ShortTime.frameSum`, slot pairs `(0,1)`, `(2,3)`). Our principal term is the double `g`-trace
`−tr_g tr_g ζ` of the coordinate 4-jet (`Pjet_add_top`). Here the two are identified, for the data
`Ef = g⁻¹` (`cometricAt`, the cometric of the pulled-back metric read in the chart at `α`):

* `cometricAt α x L` is the symmetrized cometric on covectors at `x`, symmetric for every `L`
  (`cometricAt_symm`, as `MMData.Ef_symm` requires for all jets), equal to the cometric
  `ξ ↦ ξ(g⁻¹η)` for injective `L` (`cometricAt_eq`) and positive definite there (`cometricAt_pos`);
* `rd4 α x ζ` is the `4`-tensor `ζ` read in the chart at `α`, as an `E [×4]→L[ℝ] ℝ`;
* `lam_frameSum_cometricAt`: `lam (frameSum 2 x (fun _ => cometricAt α x L) …  ζ) = −tr_g tr_g (rd4 α x ζ)`
  with the slot pattern `(a, c, b, d)` of `Pjet_add_top`;
* `rd4_eq_trivAt`, `lam_eq_trivAt`: the chart readings in terms of MM's trivialisations `trivAt`.
-/

noncomputable section

set_option linter.unusedSectionVars false

open Bundle Manifold Set
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle
open DifferentialGeometry.Geometry.Connection
open ShortTime
open scoped Manifold ContDiff BigOperators InnerProductSpace RealInnerProductSpace

namespace CWFlow.Bridge

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-! ### Covectors in a chart and the cometric -/

/-- A covector at `x`, read in the chart at `α`: `ξ ∘ (Dφ_α)⁻¹`. -/
def chartCovec (α x : M) (ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) : E →L[ℝ] ℝ :=
  ξ.comp ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL ℝ x)

theorem chartCovec_add (α x : M) (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    chartCovec α x (ξ + η) = chartCovec α x ξ + chartCovec α x η :=
  ContinuousLinearMap.add_comp _ _ _

theorem chartCovec_smul (α x : M) (c : ℝ) (ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    chartCovec α x (c • ξ) = c • chartCovec α x ξ :=
  ContinuousLinearMap.smul_comp _ _ _

theorem mem_tangent_baseSet {α x : M} (hx : x ∈ (chartAt E α).source) :
    x ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).baseSet := by
  rw [TangentBundle.trivializationAt_baseSet]; exact hx

theorem chartCovec_eq_zero_iff {α x : M} (hx : x ∈ (chartAt E α).source)
    (ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) : chartCovec α x ξ = 0 ↔ ξ = 0 := by
  constructor
  · intro h
    ext v
    have h1 := congrArg (fun f : E →L[ℝ] ℝ =>
      f ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).continuousLinearMapAt ℝ x v)) h
    have h2 := (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL_continuousLinearMapAt (R := ℝ)
      (mem_tangent_baseSet hx) v
    simp only [chartCovec, ContinuousLinearMap.comp_apply, zero_apply] at h1
    rw [h2] at h1
    exact h1
  · rintro rfl
    ext v
    simp [chartCovec]

/-- The symmetrized cometric of `L` on covectors at `x`, read in the chart at `α`:
`(ξ, η) ↦ ½ (ξ_α(g⁻¹ η_α) + η_α(g⁻¹ ξ_α))`, `ξ_α = chartCovec α x ξ`. -/
def cometricAt (α x : M) (L : E →L[ℝ] V) :
    (TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) →ₗ[ℝ] (TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) →ₗ[ℝ] ℝ :=
  LinearMap.mk₂ ℝ (fun ξ η => (1 / 2 : ℝ) * (chartCovec α x ξ (ginvOf L (chartCovec α x η)) +
      chartCovec α x η (ginvOf L (chartCovec α x ξ))))
    (fun ξ ξ' η => by simp only [chartCovec_add, map_add, add_apply]; ring)
    (fun c ξ η => by simp only [chartCovec_smul, map_smul, smul_apply, smul_eq_mul]; ring)
    (fun ξ η η' => by simp only [chartCovec_add, map_add, add_apply]; ring)
    (fun c ξ η => by simp only [chartCovec_smul, map_smul, smul_apply, smul_eq_mul]; ring)

theorem cometricAt_apply (α x : M) (L : E →L[ℝ] V) (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    cometricAt α x L ξ η = (1 / 2 : ℝ) * (chartCovec α x ξ (ginvOf L (chartCovec α x η)) +
      chartCovec α x η (ginvOf L (chartCovec α x ξ))) := rfl

theorem cometricAt_symm (α x : M) (L : E →L[ℝ] V) (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    cometricAt α x L ξ η = cometricAt α x L η ξ := by
  simp only [cometricAt_apply]; ring

/-- For injective `L`, the symmetrized cometric is the cometric `ξ_α(g⁻¹ η_α)`. -/
theorem cometricAt_eq (α x : M) {L : E →L[ℝ] V} (hL : Function.Injective L)
    (ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ) :
    cometricAt α x L ξ η = chartCovec α x ξ (ginvOf L (chartCovec α x η)) := by
  rw [cometricAt_apply, cometric_symm hL (chartCovec α x η)]; ring

theorem cometricAt_pos {α x : M} (hx : x ∈ (chartAt E α).source) {L : E →L[ℝ] V}
    (hL : Function.Injective L) {ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ} (hξ : ξ ≠ 0) :
    0 < cometricAt α x L ξ ξ := by
  rw [cometricAt_eq α x hL]
  exact cometric_pos hL (mt (chartCovec_eq_zero_iff hx ξ).1 hξ)

/-- The coordinate coframe of MM, read in its own chart, is the dual basis of `E`. (The index is
typed as in `frameSum`, over `finrank ℝ (TangentSpace 𝓘(ℝ, E) x)`, so that `simp` can use it there.) -/
theorem chartCovec_coordCoframe {α x : M} (hx : x ∈ (chartAt E α).source)
    (i : Fin (Module.finrank ℝ (TangentSpace 𝓘(ℝ, E) x))) :
    chartCovec α x (coordCoframe (I := 𝓘(ℝ, E)) α x i) = dualBas (E := E) i := by
  ext v
  simp only [chartCovec, coordCoframe, dualBas, ContinuousLinearMap.comp_apply]
  rw [Trivialization.continuousLinearMapAt_symmL _ (mem_tangent_baseSet hx)]

/-! ### The 4-tensor in a chart -/

/-- The `4`-tensor `ζ` at `x`, read in the chart at `α`: `m ↦ ζ((Dφ_α)⁻¹ m₀, …, (Dφ_α)⁻¹ m₃)`. -/
def rd4 (α x : M) (ζ : TensorRSSpace 0 (0 + 2 * 2) 𝓘(ℝ, E) x) : E [×4]→L[ℝ] ℝ :=
  (Tensor0SSpace.toModel
    ((show Tensor0SSpace 0 𝓘(ℝ, E) x →L[ℝ] Tensor0SSpace (0 + 2 * 2) 𝓘(ℝ, E) x from ζ)
      (unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x)) :
      ContinuousMultilinearMap ℝ (fun _ : Fin 4 => E) ℝ).compContinuousLinearMap
    (fun _ => ((tangentSpaceModelContinuousLinearEquiv (I := 𝓘(ℝ, E)) x :
      TangentSpace 𝓘(ℝ, E) x →L[ℝ] E).comp
        ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL ℝ x)))

/-- `rd4 α x ζ m = lam (evalFirst 4 ((Dφ_α)⁻¹ m) ζ)`: the scalar value of `ζ` on the chart vectors. -/
theorem rd4_apply (α x : M) (ζ : TensorRSSpace 0 (0 + 2 * 2) 𝓘(ℝ, E) x) (m : Fin 4 → E) :
    rd4 α x ζ m = lam x (evalFirst (r := 0) (s := 0) 4
      (fun k => (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL ℝ x (m k)) ζ) := by
  rw [lam_evalFirst]
  simp only [rd4, ContinuousMultilinearMap.compContinuousLinearMap_apply]
  congr 1
  funext k
  simp only [appendIdx, dif_pos k.isLt]
  rfl

/-- `interleave` of two coordinate frames, as a `4`-vector of basis vectors. -/
theorem interleave_coordFrame (α x : M)
    (K K' : Fin 2 → Fin (Module.finrank ℝ (TangentSpace 𝓘(ℝ, E) x))) :
    interleave (fun i => coordFrame (I := 𝓘(ℝ, E)) α x (K i))
        (fun i => coordFrame (I := 𝓘(ℝ, E)) α x (K' i)) =
      fun k => (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL ℝ x
        (![bas E (K 0), bas E (K' 0), bas E (K 1), bas E (K' 1)] k) := by
  funext k
  obtain ⟨k, hk⟩ := k
  match k, hk with
  | 0, _ => simp [interleave, coordFrame]
  | 1, _ => simp [interleave, coordFrame]
  | 2, _ => simp [interleave, coordFrame]
  | 3, _ => simp [interleave, coordFrame]
  | n + 4, h => exact absurd h (by omega)

theorem lam_evalFirst_interleave (α x : M) (ζ : TensorRSSpace 0 (0 + 2 * 2) 𝓘(ℝ, E) x)
    (K K' : Fin 2 → Fin (Module.finrank ℝ (TangentSpace 𝓘(ℝ, E) x))) :
    lam x (evalFirst (r := 0) (s := 0) (2 * 2)
      (interleave (fun i => coordFrame (I := 𝓘(ℝ, E)) α x (K i))
        (fun i => coordFrame (I := 𝓘(ℝ, E)) α x (K' i))) ζ) =
      rd4 α x ζ ![bas E (K 0), bas E (K' 0), bas E (K 1), bas E (K' 1)] := by
  rw [interleave_coordFrame, rd4_apply]

/-! ### Multilinear expansions -/

theorem vec4_update_one (a y b w v : E) :
    Function.update ![a, y, b, w] 1 v = ![a, v, b, w] := by
  funext j
  fin_cases j <;> rfl

theorem vec4_update_three (a v b y w : E) :
    Function.update ![a, v, b, y] 3 w = ![a, v, b, w] := by
  funext j
  fin_cases j <;> rfl

theorem cml4_slot1 (T : E [×4]→L[ℝ] ℝ) (a b w v : E) :
    T ![a, v, b, w] = ∑ k, dualBas k v * T ![a, bas E k, b, w] := by
  have h : ∀ y, T ![a, y, b, w] = T.toContinuousLinearMap ![a, 0, b, w] 1 y := fun y => by
    rw [← vec4_update_one a 0 b w y]; rfl
  conv_lhs => rw [← sum_dualBas_smul v]
  simp only [h, map_sum, map_smul, smul_eq_mul]

theorem cml4_slot3 (T : E [×4]→L[ℝ] ℝ) (a v b w : E) :
    T ![a, v, b, w] = ∑ l, dualBas l w * T ![a, v, b, bas E l] := by
  have h : ∀ y, T ![a, v, b, y] = T.toContinuousLinearMap ![a, v, b, 0] 3 y := fun y => by
    rw [← vec4_update_three a v b 0 y]; rfl
  conv_lhs => rw [← sum_dualBas_smul w]
  simp only [h, map_sum, map_smul, smul_eq_mul]

/-- Expansion of slots `1` and `3` in the basis. -/
theorem cml4_expand (T : E [×4]→L[ℝ] ℝ) (a v b w : E) :
    T ![a, v, b, w] = ∑ k, ∑ l, (dualBas k v * dualBas l w) * T ![a, bas E k, b, bas E l] := by
  rw [cml4_slot1]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [cml4_slot3, Finset.mul_sum]
  refine Finset.sum_congr rfl fun l _ => ?_
  ring

/-- A sum over `Fin 2 → ι` as a double sum. -/
theorem sum_fin_two_arrow {ι : Type*} [Fintype ι] {W : Type*} [AddCommMonoid W]
    (f : (Fin 2 → ι) → W) : ∑ K, f K = ∑ i, ∑ j, f ![i, j] := by
  rw [Fintype.sum_equiv (finTwoArrowEquiv ι) f (fun p => f ![p.1, p.2]) (fun K => ?_),
    Fintype.sum_prod_type]
  congr 1
  funext j
  fin_cases j <;> rfl

/-! ### The main identification -/

/-- **MM's contraction is the double trace.** For `Ef = g⁻¹` (the cometric of `L` in the chart at
`α`), `lam (frameSum 2 x Ef (stdFrame x) (stdCoframe x) ζ) = −tr_g^{(a,c)} tr_g^{(b,d)} ζ_α(a, c, b, d)`,
where `ζ_α = rd4 α x ζ` is `ζ` read in the chart. -/
theorem lam_frameSum_cometricAt {α x : M} (hx : x ∈ (chartAt E α).source) {L : E →L[ℝ] V}
    (hL : Function.Injective L) (ζ : TensorRSSpace 0 (0 + 2 * 2) 𝓘(ℝ, E) x) :
    lam x (frameSum 2 x (fun _ => cometricAt α x L) (stdFrame x) (stdCoframe x) ζ) =
      - trOf L (fun a c => trOf L fun b d => rd4 α x ζ ![a, c, b, d]) := by
  have hx' := mem_tangent_baseSet hx
  rw [frameSum_eq_frameSum 2 x _ (stdFrame x) (coordFrame (I := 𝓘(ℝ, E)) α x) (stdCoframe x)
    (coordCoframe (I := 𝓘(ℝ, E)) α x) (stdCoframe_stdFrame x)
    (coordCoframe_coordFrame (I := 𝓘(ℝ, E)) α x hx') ζ]
  have h3 : ((-1 : ℝ) ^ (2 + 1)) = -1 := by norm_num
  simp only [frameSum, map_smul, map_sum, smul_eq_mul, Fin.prod_univ_two, cometricAt_eq α x hL,
    chartCovec_coordCoframe hx, lam_evalFirst_interleave, sum_fin_two_arrow, Matrix.cons_val_zero,
    Matrix.cons_val_one, h3, neg_one_mul, neg_inj]
  simp only [trOf]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  rw [cml4_expand]
  refine Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun l _ => ?_
  congr 2 <;> exact cometric_symm hL _ _

/-! ### The chart readings in MM's trivialisations -/

variable (E) in
/-- The unit `0`-form of the model: the constant `1`. -/
def unitModel0 : Tensor0SModel 0 ℝ E :=
  ContinuousMultilinearMap.constOfIsEmpty ℝ (fun _ : Fin 0 => E) (1 : ℝ)

/-- The unit model `0`-form, transported to the fibre at `x` by the trivialisation at `α`, is the
unit section. -/
theorem constInChart_unitModel0 {α x : M} (hx : x ∈ (chartAt E α).source) :
    Tensor0SSpace.constInChart (𝕜 := ℝ) (I := 𝓘(ℝ, E)) 0 α (unitModel0 E) x =
      unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x := by
  apply Tensor0SSpace.toModel_injective
  ext v
  simp only [Tensor0SSpace.toModel_apply_model_vector]
  rw [Tensor0SSpace.constInChart_apply (𝕜 := ℝ) (I := 𝓘(ℝ, E)) 0 (mem_tangent_baseSet hx),
    unitZeroSec_apply, ← Tensor0SSpace.eval_eq, Tensor0SSpace.eval_ofModel]
  rfl

/-- The scalar `lam x X` of a `(0,0)`-tensor is its reading in `trivAt 0 0`, at the unit. -/
theorem lam_eq_trivAt {α x : M} (hx : x ∈ (chartAt E α).source) (X : TensorRSSpace 0 0 𝓘(ℝ, E) x) :
    lam x X = ((trivAt 𝓘(ℝ, E) M 0 0 α).continuousLinearMapAt ℝ x X) (unitModel0 E)
      (fun i => i.elim0) := by
  have hx' := mem_tangent_baseSet hx
  rw [Trivialization.continuousLinearMapAt_apply,
    Trivialization.coe_linearMapAt_of_mem _ (trivAt_mem_of_tangent α x 0 0 hx'),
    TensorRSSpace.trivializationAt_apply (𝕜 := ℝ) (I := 𝓘(ℝ, E)) (x₀ := α) (x := x) 0 0 hx']
  have hc := constInChart_unitModel0 (E := E) hx
  unfold Tensor0SSpace.constInChart at hc
  rw [hc]
  show Tensor0SSpace.toModel ((show Tensor0SSpace 0 𝓘(ℝ, E) x →L[ℝ] Tensor0SSpace 0 𝓘(ℝ, E) x from X)
    (unitZeroSec (I := 𝓘(ℝ, E)) (M := M) x)) (fun i => i.elim0) = _
  rw [Tensor0SSpace.toModel_apply_model_vector]
  congr 1
  funext i
  exact i.elim0

/-- The chart reading `rd4` is the reading of `ζ` in `trivAt 0 4`, at the unit `0`-form. -/
theorem rd4_eq_trivAt {α x : M} (hx : x ∈ (chartAt E α).source)
    (ζ : TensorRSSpace 0 (0 + 2 * 2) 𝓘(ℝ, E) x) (m : Fin 4 → E) :
    rd4 α x ζ m =
      ((trivAt 𝓘(ℝ, E) M 0 (0 + 4) α).continuousLinearMapAt ℝ x ζ) (unitModel0 E) m := by
  have hx' := mem_tangent_baseSet hx
  rw [Trivialization.continuousLinearMapAt_apply,
    Trivialization.coe_linearMapAt_of_mem _ (trivAt_mem_of_tangent α x 0 (0 + 4) hx'),
    TensorRSSpace.trivializationAt_apply (𝕜 := ℝ) (I := 𝓘(ℝ, E)) (x₀ := α) (x := x) 0 (0 + 4) hx']
  have hc := constInChart_unitModel0 (E := E) hx
  unfold Tensor0SSpace.constInChart at hc
  rw [hc]
  simp only [rd4, ContinuousMultilinearMap.compContinuousLinearMap_apply,
    Tensor0SSpace.toModel_apply_model_vector]
  rfl

/-- info: 'CWFlow.Bridge.lam_frameSum_cometricAt' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lam_frameSum_cometricAt

/-- info: 'CWFlow.Bridge.rd4_eq_trivAt' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms rd4_eq_trivAt

/-- info: 'CWFlow.Bridge.cometricAt_pos' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cometricAt_pos

end CWFlow.Bridge
