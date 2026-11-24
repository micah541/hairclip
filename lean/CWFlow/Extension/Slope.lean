import CWFlow.Extension.GraphLocal
import CWFlow.Extension.ImmReparam
import CWFlow.Interp.BallInterp
import CWFlow.Extension.EnergySymbol

/-!
# The slope of a Lagrangian plane over a fixed Lagrangian plane

For a Lagrangian frame `U` of `V` and a linear map `L : ℝⁿ → V`, the slope `slopeZ J U L` is the matrix `S`
with `verCL ∘ L = S ∘ (horCL ∘ L)` (when `horCL ∘ L` is invertible), i.e. the plane `L(ℝⁿ)` is
`{U w + J U (S w)}`. The slope is invariant under reparametrization (`slopeZ_comp`) and for the graph of a
gradient it is the Hessian (`slopeZ_graphX`). Further:

* `norm_jet3_sq_le`: `|D³φ|² ≤ (1 + |D²φ|²)³ |A|²_g`;
* `norm_lt_of_deriv_le_Ico`: a one-dimensional comparison lemma keeping the slope small;
* `exists_local_graph`: local graph representations around points where `horCL ∘ DX` is injective.
-/

open Metric
open CWFlow.Interp
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-! ### Part 2: the third derivative is controlled by `|A|` -/

theorem norm_jet3_sq_le (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    ‖jet 3 φ x‖ ^ 2 ≤ (1 + ‖jet 2 φ x‖ ^ 2) ^ 3 * gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x) := by
  have h1 : ((1 + ‖jet 2 φ x‖ ^ 2)⁻¹) ^ 3 * ‖jet 3 φ x‖ ^ 2 ≤
      gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x) := by
    refine le_gNormSq (by positivity) (ginv_hsymm _) (ginv_hnn _) ?_ _
    intro v
    simp only [ginvFun_apply]
    exact ginvFun_quadratic_ge _ v
  have hp : 0 < (1 + ‖jet 2 φ x‖ ^ 2) ^ 3 := by positivity
  calc ‖jet 3 φ x‖ ^ 2
      = (1 + ‖jet 2 φ x‖ ^ 2) ^ 3 * (((1 + ‖jet 2 φ x‖ ^ 2)⁻¹) ^ 3 * ‖jet 3 φ x‖ ^ 2) := by
        field_simp
    _ ≤ _ := by gcongr

/-! ### Part 3: a one-dimensional comparison lemma -/

theorem norm_lt_of_deriv_le_Ico {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {b a c : ℝ}
    (ha : 0 ≤ a) (hc : 0 < c) {u : ℝ → W} (hcont : ContinuousOn u (Set.Ico 0 b))
    (hder : ∀ s ∈ Set.Ioo 0 b, ∃ u' : W, HasDerivAt u u' s ∧
      (‖u s‖ ≤ c → ‖u'‖ ≤ a * (1 + ‖u s‖ ^ 2) ^ ((3 : ℝ) / 2)))
    (h0 : u 0 = 0) (hsmall : b * (a * (1 + c ^ 2) ^ ((3 : ℝ) / 2)) < c) :
    ∀ s ∈ Set.Ico 0 b, ‖u s‖ < c := by
  intro s0 hs0
  by_contra hcon
  push Not at hcon
  choose! u' hu' using hder
  have hb0 : 0 < b := lt_of_le_of_lt hs0.1 hs0.2
  set S : Set ℝ := Set.Icc 0 s0 ∩ (fun t => ‖u t‖) ⁻¹' Set.Ici c with hS
  have hsub : Set.Icc 0 s0 ⊆ Set.Ico 0 b := fun t ht => ⟨ht.1, lt_of_le_of_lt ht.2 hs0.2⟩
  have hSc : IsClosed S :=
    (hcont.mono hsub).norm.preimage_isClosed_of_isClosed isClosed_Icc isClosed_Ici
  have hsS : s0 ∈ S := ⟨⟨hs0.1, le_rfl⟩, hcon⟩
  have hSb : BddBelow S := ⟨0, fun t ht => ht.1.1⟩
  have hts : sInf S ∈ S := hSc.csInf_mem ⟨s0, hsS⟩ hSb
  set ts := sInf S with hts_def
  have hts01 : ts ∈ Set.Icc 0 s0 := hts.1
  have hts_c : c ≤ ‖u ts‖ := hts.2
  have hlt : ∀ t ∈ Set.Ico 0 ts, ‖u t‖ < c := by
    intro t ht
    by_contra h
    push Not at h
    have : ts ≤ t := csInf_le hSb ⟨⟨ht.1, ht.2.le.trans hts01.2⟩, h⟩
    exact absurd ht.2 (not_lt.2 this)
  have hts_pos : 0 < ts := by
    rcases hts01.1.eq_or_lt with h | h
    · exfalso
      have : ‖u ts‖ = 0 := by rw [← h, h0]; simp
      linarith
    · exact h
  have htsb : ts < b := lt_of_le_of_lt hts01.2 hs0.2
  set K : ℝ := a * (1 + c ^ 2) ^ ((3 : ℝ) / 2) with hK_def
  have hK : 0 ≤ K := by positivity
  have hmain : ∀ ε > 0, c ≤ K * ts + ε := by
    intro ε hε
    obtain ⟨δ, hδ, hδε⟩ := Metric.continuousWithinAt_iff.1 (hcont 0 ⟨le_rfl, hb0⟩) ε hε
    set d := min (δ / 2) (ts / 2) with hd
    have hd0 : 0 < d := by positivity
    have hdts : d < ts := lt_of_le_of_lt (min_le_right _ _) (by linarith)
    have hdδ : d < δ := lt_of_le_of_lt (min_le_left _ _) (by linarith)
    have hdb : d ∈ Set.Ico 0 b := ⟨hd0.le, by linarith⟩
    have h1 := hδε hdb (by simpa [Real.dist_eq, abs_of_pos hd0] using hdδ)
    rw [h0, dist_zero_right] at h1
    have hI : Set.Icc d ts ⊆ Set.Ico 0 b := fun t ht => ⟨hd0.le.trans ht.1, by linarith [ht.2]⟩
    have hmv := norm_image_sub_le_of_norm_deriv_right_le_segment (f := u) (a := d) (b := ts)
      (f' := u') (C := K) (hcont.mono hI)
      (fun t ht => (hu' t ⟨by linarith [ht.1], by linarith [ht.2]⟩).1.hasDerivWithinAt)
      (fun t ht => by
        have ht0 : t ∈ Set.Ico 0 ts := ⟨by linarith [ht.1], ht.2⟩
        have h2 := (hu' t ⟨by linarith [ht.1], by linarith [ht.2]⟩).2 (hlt t ht0).le
        refine h2.trans ?_
        rw [hK_def]
        have : ‖u t‖ ≤ c := (hlt t ht0).le
        gcongr) ts ⟨hdts.le, le_rfl⟩
    have h3 : ‖u ts‖ ≤ ‖u d‖ + ‖u ts - u d‖ := by
      have := norm_add_le (u d) (u ts - u d)
      simpa using this
    have h4 : K * (ts - d) ≤ K * ts := by nlinarith
    linarith
  have hc' : c ≤ K * ts := le_of_forall_pos_le_add hmain
  have : K * ts ≤ K * b := by gcongr
  nlinarith

/-! ### Part 1: the slope matrix -/

/-- The slope of the plane `L(ℝⁿ)` over `U(ℝⁿ)`: `S_{ab} = (verCL (L ((horCL ∘ L)⁻¹ e_a)))_b`, i.e. the plane is
`{U w + J U (S w)}` with `S = verCL ∘ L ∘ (horCL ∘ L)⁻¹`. -/
noncomputable def slopeZ (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V)
    (L : EuclideanSpace ℝ (Fin n) →L[ℝ] V) : Z n 2 :=
  WithLp.toLp 2 (fun I => verCL J U (L ((horCL U ∘L L).inverse (eb (I 0)))) (I 1))

/-- The closed slope condition. -/
def slopeLE (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (c : ℝ)
    (L : EuclideanSpace ℝ (Fin n) →L[ℝ] V) : Prop :=
  ∀ a, ‖verCL J U (L a)‖ ≤ c * ‖horCL U (L a)‖

theorem sl_slope_apply (L : 𝔼 →L[ℝ] V) (I : Fin 2 → Fin n) :
    slopeZ J U L I = verCL J U (L ((horCL U ∘L L).inverse (eb (I 0)))) (I 1) := rfl

theorem sl_inverse {L : 𝔼 →L[ℝ] V} (hL : Function.Injective (horCL U ∘L L)) :
    ∃ E : 𝔼 ≃L[ℝ] 𝔼, (∀ x, E x = horCL U (L x)) ∧ horCL U ∘L L = (E : 𝔼 →L[ℝ] 𝔼) ∧
      (horCL U ∘L L).inverse = (E.symm : 𝔼 →L[ℝ] 𝔼) := by
  obtain ⟨E, hE⟩ := gl_exists_equiv hL
  have h2 : horCL U ∘L L = (E : 𝔼 →L[ℝ] 𝔼) := hE.symm
  refine ⟨E, fun x => ?_, h2, ?_⟩
  · have := congrArg (fun f : 𝔼 →L[ℝ] 𝔼 => f x) hE
    simpa using this
  · rw [h2, ContinuousLinearMap.inverse_equiv]

theorem slopeZ_comp {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V}
    (T : EuclideanSpace ℝ (Fin n) ≃L[ℝ] EuclideanSpace ℝ (Fin n))
    (hL : Function.Injective (horCL U ∘L L)) :
    slopeZ J U (L ∘L (T : 𝔼 →L[ℝ] 𝔼)) = slopeZ J U L := by
  obtain ⟨E, hEx, h2, hinv⟩ := sl_inverse hL
  have h1 : horCL U ∘L (L ∘L (T : 𝔼 →L[ℝ] 𝔼)) = ((T.trans E : 𝔼 ≃L[ℝ] 𝔼) : 𝔼 →L[ℝ] 𝔼) := by
    ext w; simp [hEx]
  unfold slopeZ
  congr 1
  funext I
  rw [h1, hinv, ContinuousLinearMap.inverse_equiv]
  simp

theorem horCL_comp_fderiv_graphX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) (y₀ : V) (x : EuclideanSpace ℝ (Fin n)) :
    horCL U ∘L fderiv ℝ (graphX J U y₀ φ) x = ContinuousLinearMap.id ℝ _ := by
  ext u a
  rw [ContinuousLinearMap.comp_apply, gb_fderiv_graphX_apply hφ, horCL_apply]
  simp only [inner_add_right, gb_UU, gb_UJU hU, add_zero]
  simp [eb, EuclideanSpace.inner_single_left]

theorem slopeZ_graphX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) (y₀ : V) (x : EuclideanSpace ℝ (Fin n)) :
    slopeZ J U (fderiv ℝ (graphX J U y₀ φ) x) = jet 2 φ x := by
  have hid := horCL_comp_fderiv_graphX hJi hJJ hU hφ y₀ x
  ext I
  have hI : I = ![I 0, I 1] := by
    funext i; fin_cases i <;> rfl
  rw [sl_slope_apply, hid, ContinuousLinearMap.inverse_id]
  rw [gb_fderiv_graphX_apply hφ]
  simp only [ContinuousLinearMap.id_apply, verCL_apply, inner_add_right, gb_JUU hU, gb_JUJU hJi,
    zero_add]
  rw [gb_gradE_comp hφ]
  conv_rhs => rw [hI]
  simp [eb, EuclideanSpace.inner_single_left]

theorem isClosed_slopeLE (c : ℝ) :
    IsClosed {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V | slopeLE J U c L} := by
  have : {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V | slopeLE J U c L} =
      ⋂ a : 𝔼, {L : 𝔼 →L[ℝ] V | ‖verCL J U (L a)‖ ≤ c * ‖horCL U (L a)‖} := by
    ext L; simp [slopeLE]
  rw [this]
  refine isClosed_iInter fun a => isClosed_le ?_ ?_
  · exact continuous_norm.comp ((verCL J U).continuous.comp
      (ContinuousLinearMap.apply ℝ V a).continuous)
  · exact continuous_const.mul (continuous_norm.comp ((horCL U).continuous.comp
      (ContinuousLinearMap.apply ℝ V a).continuous))

theorem sl_decomp_lin (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) (w : V) :
    w = U (horCL U w) + J (U (verCL J U w)) := by
  have := decomp_horC_verC (y₀ := (0 : V)) hJi hJJ hU hdim w
  simpa [horC_eq, verC_eq] using this

theorem injective_horCL_of_slopeLE (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) {c : ℝ} {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V}
    (hL : Function.Injective L) (hs : slopeLE J U c L) :
    Function.Injective (horCL U ∘L L) := by
  rw [injective_iff_map_eq_zero]
  intro a ha
  have h1 : horCL U (L a) = 0 := ha
  have h2 : verCL J U (L a) = 0 := by
    have := hs a
    rw [h1, norm_zero, mul_zero] at this
    exact norm_le_zero_iff.1 this
  have h3 := sl_decomp_lin hJi hJJ hU hdim (L a)
  rw [h1, h2] at h3
  simp only [map_zero, add_zero] at h3
  exact (injective_iff_map_eq_zero L).1 hL a h3

theorem sl_ver_eq {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V} (hL : Function.Injective (horCL U ∘L L))
    (w : 𝔼) : verCL J U (L w) =
      ∑ a, horCL U (L w) a • verCL J U (L ((horCL U ∘L L).inverse (eb a))) := by
  obtain ⟨E, hEx, h2, hinv⟩ := sl_inverse hL
  have hw : w = ∑ a, (E w) a • E.symm (eb a) := by
    calc w = E.symm (E w) := (E.symm_apply_apply w).symm
      _ = E.symm (∑ j, (E w) j • eb j) := by rw [← ir_expand (E w)]
      _ = _ := by simp [map_sum, map_smul]
  have h3 : L w = ∑ a, (E w) a • L (E.symm (eb a)) := by
    conv_lhs => rw [hw]
    simp [map_sum, map_smul]
  rw [hinv]
  conv_lhs => rw [h3]
  simp [map_sum, map_smul, ← hEx]

theorem sl_normsq2 (s : Z n 2) : ‖s‖ ^ 2 = ∑ a, ∑ b, s ![a, b] ^ 2 := by
  rw [normSq_eq_sum]
  rw [← Equiv.sum_comp (piFinTwoEquiv fun _ : Fin 2 => Fin n).symm (fun I => s I ^ 2)]
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
  have : (piFinTwoEquiv fun _ : Fin 2 => Fin n).symm (a, b) = ![a, b] := by
    funext i; fin_cases i <;> rfl
  rw [this]

theorem sl_frob (s : Z n 2) (w : 𝔼) :
    ‖(WithLp.toLp 2 (fun b => ∑ a, w a * s ![a, b]) : 𝔼)‖ ≤ ‖s‖ * ‖w‖ := by
  refine (sq_le_sq₀ (norm_nonneg _) (by positivity)).1 ?_
  rw [mul_pow, sl_normsq2, EuclideanSpace.norm_sq_eq, EuclideanSpace.norm_sq_eq]
  simp only [Real.norm_eq_abs, sq_abs]
  calc ∑ b, (∑ a, w a * s ![a, b]) ^ 2
      ≤ ∑ b, (∑ a, w a ^ 2) * (∑ a, s ![a, b] ^ 2) :=
        Finset.sum_le_sum fun b _ => Finset.sum_mul_sq_le_sq_mul_sq _ _ _
    _ = _ := by
      rw [← Finset.mul_sum, mul_comm]
      congr 1
      exact Finset.sum_comm

theorem slopeLE_of_norm_slopeZ_le (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) {c : ℝ} {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V}
    (hL : Function.Injective (horCL U ∘L L)) (hs : ‖slopeZ J U L‖ ≤ c) : slopeLE J U c L := by
  intro a
  have h1 : verCL J U (L a) = (WithLp.toLp 2 (fun b => ∑ a', horCL U (L a) a' *
      slopeZ J U L ![a', b]) : 𝔼) := by
    ext b
    rw [sl_ver_eq hL a]
    simp [sl_slope_apply]
  rw [h1]
  exact (sl_frob _ _).trans (mul_le_mul_of_nonneg_right hs (norm_nonneg _))

theorem continuousAt_slopeZ {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V}
    (hL : Function.Injective (horCL U ∘L L)) : ContinuousAt (slopeZ J U) L := by
  obtain ⟨E, hEx, h2, hinv⟩ := sl_inverse hL
  have hc : Continuous (fun L' : 𝔼 →L[ℝ] V => horCL U ∘L L') :=
    (ContinuousLinearMap.compL ℝ 𝔼 V 𝔼 (horCL U)).continuous
  have hT : ContinuousAt (fun L' : 𝔼 →L[ℝ] V => (horCL U ∘L L').inverse) L := by
    simp only [← ContinuousLinearMap.ringInverse_eq_inverse]
    have hu := NormedRing.inverse_continuousAt ((ContinuousLinearEquiv.unitsEquiv ℝ 𝔼).symm E)
    exact hu.comp_of_eq hc.continuousAt h2
  unfold slopeZ
  refine (PiLp.continuous_toLp 2 _).continuousAt.comp (continuousAt_pi.2 fun I => ?_)
  have h3 : ContinuousAt (fun L' : 𝔼 →L[ℝ] V => L' ((horCL U ∘L L').inverse (eb (I 0)))) L :=
    continuousAt_id.clm_apply (hT.clm_apply continuousAt_const)
  have h4 := ((verCL J U).continuous.continuousAt.comp h3)
  exact ((EuclideanSpace.proj (I 1) : 𝔼 →L[ℝ] ℝ).continuous.continuousAt.comp h4)

/-! ### Part 4: local graph representation -/

theorem sl_det_ne {τ : 𝔼 → 𝔼} {z : 𝔼} (h : Function.Injective (fderiv ℝ τ z)) :
    (jacM τ z).det ≠ 0 := by
  have hu : IsUnit (jacM τ z) := by
    rw [← Matrix.mulVec_injective_iff_isUnit]
    intro v w hvw
    have key : ∀ v : Fin n → ℝ, (jacM τ z).mulVec v =
        fun j => fderiv ℝ τ z (WithLp.toLp 2 v) j := by
      intro v
      funext j
      have h1 : fderiv ℝ τ z (WithLp.toLp 2 v) j =
          ∑ i, (WithLp.toLp 2 v : 𝔼) i * jacM τ z j i := by
        rw [ir_lin_expand (fderiv ℝ τ z) (WithLp.toLp 2 v)]
        simp [jacM]
      rw [h1]
      simp only [Matrix.mulVec, dotProduct]
      simp [mul_comm]
    rw [key, key] at hvw
    have h2 : fderiv ℝ τ z (WithLp.toLp 2 v) = fderiv ℝ τ z (WithLp.toLp 2 w) := by
      ext j; exact congrFun hvw j
    have h3 := h h2
    simpa using h3
  exact isUnit_iff_ne_zero.1 ((Matrix.isUnit_iff_isUnit_det _).1 hu)

theorem sl_graphX_congr {φ φ₀ : 𝔼 → ℝ} {x : 𝔼} (h : φ =ᶠ[𝓝 x] φ₀) (y₀ : V) :
    graphX J U y₀ φ x = graphX J U y₀ φ₀ x := by
  have : gradE φ x = gradE φ₀ x := by
    unfold gradE
    rw [h.fderiv_eq]
  unfold graphX
  rw [this]

theorem exists_local_graph (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) (y₀ : V)
    {X : EuclideanSpace ℝ (Fin n) → V} {W : Set (EuclideanSpace ℝ (Fin n))} (hW : IsOpen W)
    (hX : ContDiffOn ℝ ∞ X W) (hinj : ∀ y ∈ W, Function.Injective (fderiv ℝ X y))
    (hlag : ∀ y ∈ W, ∀ u w, ⟪J (fderiv ℝ X y u), fderiv ℝ X y w⟫ = 0)
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ W)
    (hpy : Function.Injective (horCL U ∘L fderiv ℝ X y)) :
    ∃ r > 0, ∃ (s : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
      (φ : EuclideanSpace ℝ (Fin n) → ℝ),
      ContDiff ℝ ∞ φ ∧ ContDiffOn ℝ ∞ s (ball (horC U y₀ (X y)) r) ∧
      Set.MapsTo s (ball (horC U y₀ (X y)) r) W ∧ s (horC U y₀ (X y)) = y ∧
      (∀ x ∈ ball (horC U y₀ (X y)) r, horC U y₀ (X (s x)) = x ∧ X (s x) = graphX J U y₀ φ x ∧
        slopeZ J U (fderiv ℝ X (s x)) = jet 2 φ x ∧
        aNormSq J X (s x) = gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x)) ∧
      ∀ᶠ y' in 𝓝 y, y' ∈ W ∧ horC U y₀ (X y') ∈ ball (horC U y₀ (X y)) r ∧
        s (horC U y₀ (X y')) = y' := by
  classical
  set χ : 𝔼 → 𝔼 := fun y => horC U y₀ (X y) with hχdef
  have hhor : ContDiff ℝ ∞ (horC U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
      (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hXat : ∀ z ∈ W, ContDiffAt ℝ ∞ X z := fun z hz => hX.contDiffAt (hW.mem_nhds hz)
  have hχat : ContDiffAt ℝ ∞ χ y := hhor.contDiffAt.comp y (hXat y hy)
  have hχder : HasFDerivAt χ (horCL U ∘L fderiv ℝ X y) y :=
    (hasFDerivAt_horC (X y)).comp y ((hXat y hy).differentiableAt (by simp)).hasFDerivAt
  obtain ⟨E, hE⟩ := gl_exists_equiv hpy
  have hn : (∞ : WithTop ℕ∞) ≠ 0 := by simp
  have hd : HasFDerivAt χ (E : 𝔼 →L[ℝ] 𝔼) y := by rw [hE]; exact hχder
  set e := hχat.toOpenPartialHomeomorph χ hd hn with he
  have he_coe : (e : 𝔼 → 𝔼) = χ := ContDiffAt.toOpenPartialHomeomorph_coe hχat hd hn
  have hy_src : y ∈ e.source := ContDiffAt.mem_toOpenPartialHomeomorph_source hχat hd hn
  -- openness of the invertibility condition near `y`
  have hDXc : ContinuousAt (fderiv ℝ X) y :=
    ((hXat y hy).fderiv_right (m := ∞) (by simp)).continuousAt
  have hMc : ContinuousAt (fun z => horCL U ∘L fderiv ℝ X z) y :=
    ((ContinuousLinearMap.compL ℝ 𝔼 V 𝔼 (horCL U)).continuous.continuousAt).comp hDXc
  have hMev : ∀ᶠ z in 𝓝 y, horCL U ∘L fderiv ℝ X z ∈
      Set.range ((↑) : (𝔼 ≃L[ℝ] 𝔼) → 𝔼 →L[ℝ] 𝔼) :=
    hMc.eventually (ContinuousLinearEquiv.isOpen.mem_nhds ⟨E, hE⟩)
  have hP : ∀ᶠ z in 𝓝 y, z ∈ W ∧ z ∈ e.source ∧ Function.Injective (horCL U ∘L fderiv ℝ X z) := by
    filter_upwards [hW.mem_nhds hy, e.open_source.mem_nhds hy_src, hMev] with z hz1 hz2 hz3
    obtain ⟨E', hE'⟩ := hz3
    refine ⟨hz1, hz2, ?_⟩
    rw [← hE']
    exact E'.injective
  obtain ⟨W', hW'P, hW'o, hyW'⟩ := _root_.eventually_nhds_iff.1 hP
  have hW'W : W' ⊆ W := fun z hz => (hW'P z hz).1
  have hW'src : W' ⊆ e.source := fun z hz => (hW'P z hz).2.1
  have hopen : IsOpen (χ '' W') := by
    have := e.isOpen_image_of_subset_source hW'o hW'src
    rwa [he_coe] at this
  obtain ⟨r', hr', hball'⟩ := Metric.isOpen_iff.1 hopen (χ y) ⟨y, hyW', rfl⟩
  have hinjOn : Set.InjOn χ W' := by
    have := e.injOn.mono hW'src
    rwa [he_coe] at this
  obtain ⟨s, φ₀, hsC, hmaps, hright, hleft, hφ₀C, -, hgraph⟩ :=
    exists_graph_of_proj (y₀ := y₀) (c := χ y) (R := r') hJi hJJ hU hdim hW'o
      (hX.mono hW'W) (fun z hz => hlag z (hW'W hz))
      (fun z hz => (hW'P z hz).2.2) hinjOn hball'
  -- a globally smooth potential
  set c := χ y with hc
  set r := r' / 2 with hr
  have hr0 : 0 < r := by positivity
  have hf : ContDiffOn ℝ ∞ (fun z : 𝔼 => φ₀ (c + z)) (ball 0 r') := by
    refine hφ₀C.comp (contDiff_const.add contDiff_id).contDiffOn ?_
    intro z hz
    simpa [mem_ball_iff_norm] using hz
  obtain ⟨R', hR', g, hgC, -, hg⟩ := exists_global_of_contDiffOn (U := ball (0 : 𝔼) r')
    isOpen_ball hr0 (closedBall_subset_ball (by rw [hr]; linarith)) hf
  set φ : 𝔼 → ℝ := fun x => g (x - c) with hφ
  have hφC : ContDiff ℝ ∞ φ := hgC.comp (contDiff_id.sub contDiff_const)
  have hBsub : ball c r ⊆ ball c r' := ball_subset_ball (by rw [hr]; linarith)
  have hφeq : ∀ x ∈ ball c r, φ x = φ₀ x := by
    intro x hx
    have h1 : x - c ∈ ball (0 : 𝔼) R' := by
      have := mem_ball_iff_norm.1 hx
      rw [mem_ball_zero_iff]
      linarith
    have := hg _ h1
    simp only [hφ]
    rw [this]
    simp
  have hφev : ∀ x ∈ ball c r, φ =ᶠ[𝓝 x] φ₀ := fun x hx =>
    Filter.eventuallyEq_of_mem (isOpen_ball.mem_nhds hx) hφeq
  have hcB : c ∈ ball c r := mem_ball_self hr0
  have hgr : ∀ x ∈ ball c r, X (s x) = graphX J U y₀ φ x := fun x hx => by
    rw [hgraph x (hBsub hx), sl_graphX_congr (hφev x hx)]
  have hsB : ContDiffOn ℝ ∞ s (ball c r) := hsC.mono hBsub
  have hmapsB : Set.MapsTo s (ball c r) W := fun x hx => hW'W (hmaps (hBsub hx))
  have hsW' : ∀ x ∈ ball c r, s x ∈ W' := fun x hx => hmaps (hBsub hx)
  have hcB' : c ∈ ball c r' := mem_ball_self hr'
  have hsc : s (horC U y₀ (X y)) = y := hleft y hyW' hcB'
  have key : ∀ z ∈ ball c r, DifferentiableAt ℝ s z ∧ DifferentiableAt ℝ X (s z) ∧
      fderiv ℝ (graphX J U y₀ φ) z = (fderiv ℝ X (s z)).comp (fderiv ℝ s z) ∧
      Function.Injective (fderiv ℝ s z) := by
    intro z hz
    have hDs : DifferentiableAt ℝ s z := (hsB.contDiffAt (isOpen_ball.mem_nhds hz)).differentiableAt
      (by simp)
    have hXd : DifferentiableAt ℝ X (s z) := (hXat _ (hmapsB hz)).differentiableAt (by simp)
    have hev : (X ∘ s) =ᶠ[𝓝 z] graphX J U y₀ φ :=
      Filter.eventuallyEq_of_mem (isOpen_ball.mem_nhds hz) (fun w hw => hgr w hw)
    have hDG : fderiv ℝ (graphX J U y₀ φ) z = (fderiv ℝ X (s z)).comp (fderiv ℝ s z) := by
      rw [← hev.fderiv_eq]
      exact fderiv_comp z hXd hDs
    have hid := horCL_comp_fderiv_graphX hJi hJJ hU hφC y₀ z
    have h1 : ∀ u, horCL U (fderiv ℝ X (s z) (fderiv ℝ s z u)) = u := fun u => by
      have := congrArg (fun L : 𝔼 →L[ℝ] 𝔼 => L u) hid
      simpa [hDG] using this
    exact ⟨hDs, hXd, hDG, fun u v h => by rw [← h1 u, ← h1 v, h]⟩
  refine ⟨r, hr0, s, φ, hφC, hsB, hmapsB, hsc, ?_, ?_⟩
  · intro x hx
    obtain ⟨hDs, hXd, hDG, hinjDs⟩ := key x hx
    have hev : (X ∘ s) =ᶠ[𝓝 x] graphX J U y₀ φ :=
      Filter.eventuallyEq_of_mem (isOpen_ball.mem_nhds hx) (fun w hw => hgr w hw)
    refine ⟨hright x (hBsub hx), hgr x hx, ?_, ?_⟩
    · obtain ⟨T, hT⟩ := gl_exists_equiv hinjDs
      have h1 := slopeZ_comp (J := J) (U := U) (L := fderiv ℝ X (s x)) T
        (hW'P (s x) (hsW' x hx)).2.2
      rw [hT] at h1
      rw [← h1, ← hDG, slopeZ_graphX hJi hJJ hU hφC y₀ x]
    · have h1 := aNormSq_comp (J := J) hW hX hinj hlag isOpen_ball hsB (fun z hz => hmapsB hz)
        (fun z hz => sl_det_ne (key z hz).2.2.2) hx
      rw [← h1, aNormSq_congr hev, aNormSq_graphX hJi hJJ hU hφC]
  · filter_upwards [hW'o.mem_nhds hyW', hχat.continuousAt.eventually_mem (isOpen_ball.mem_nhds hcB)]
      with z hz1 hz2
    exact ⟨hW'W hz1, hz2, hleft z hz1 (hBsub hz2)⟩

/-- info: 'CWFlow.Extension.slopeZ_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms slopeZ_comp

/-- info: 'CWFlow.Extension.slopeZ_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms slopeZ_graphX

/-- info: 'CWFlow.Extension.injective_horCL_of_slopeLE' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms injective_horCL_of_slopeLE

/-- info: 'CWFlow.Extension.norm_jet3_sq_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_jet3_sq_le

/-- info: 'CWFlow.Extension.norm_lt_of_deriv_le_Ico' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_lt_of_deriv_le_Ico

/-- info: 'CWFlow.Extension.exists_local_graph' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_local_graph

end CWFlow.Extension
