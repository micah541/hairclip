import CWFlow.Scalar.Flat
import CWFlow.Scalar.TaylorComp
import CWFlow.Scalar.Graph

/-!
# The scalar operator as a smooth function of the 4-jet

For a Weinstein-type map `Ψ : E × E* → V`, the operator `u ↦ −div(JH)[F_u]`, `F_u = Ψ(y, Du)`, is given at `y` by
an explicit function `Pjet J Ψ` of `(y, Du, D²u, D³u, D⁴u)(y)` (`neg_flowPot_graphFu_eq_Pjet`), smooth where
its first-order part is injective (`contDiffAt_Pjet`):

* `gSeries y z` is the Taylor series of `w ↦ (w, Du(w))` at `y` written in the jet `z` (`iteratedFDeriv_graphMap`);
* by Faà di Bruno (`iteratedFDeriv_comp`), `D^k F_u(y) = (D•Ψ(y, Du)).taylorComp (gSeries y z) k`;
* `Pjet = −PotM` of these (`flowPot_eq_potM`).
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

universe u

variable {E V : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) in
/-- The scalar 4-jet `(Du, D²u, D³u, D⁴u)`. -/
abbrev SJet := (E →L[ℝ] ℝ) × (E [×2]→L[ℝ] ℝ) × (E [×3]→L[ℝ] ℝ) × (E [×4]→L[ℝ] ℝ)

/-- The 4-jet of `u` at `y`. -/
def jetOf (u : E → ℝ) (y : E) : SJet E :=
  (fderiv ℝ u y, iteratedFDeriv ℝ 2 u y, iteratedFDeriv ℝ 3 u y, iteratedFDeriv ℝ 4 u y)

/-- `T₂(a, ·)` for a scalar bilinear `T₂`, as a linear map of `a`. -/
def slot1L (T₂ : E [×2]→L[ℝ] ℝ) : E →L[ℝ] (E →L[ℝ] ℝ) :=
  (continuousMultilinearCurryFin1 ℝ E ℝ : (E [×1]→L[ℝ] ℝ) →L[ℝ] (E →L[ℝ] ℝ)).comp
    (continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin 2 => E) ℝ T₂)

theorem slot1L_apply (T₂ : E [×2]→L[ℝ] ℝ) (a b : E) : slot1L T₂ a b = T₂ ![a, b] := by
  change T₂ (Fin.cons a fun _ => b) = T₂ ![a, b]
  congr 1
  ext i
  fin_cases i <;> rfl


/-- The Taylor series of `w ↦ (w, Du(w))` at `y`, from the jet `z`. -/
def gSeries (y : E) (z : SJet E) : FormalMultilinearSeries ℝ E (E × (E →L[ℝ] ℝ))
  | 0 => (continuousMultilinearCurryFin0 ℝ E (E × (E →L[ℝ] ℝ))).symm (y, z.1)
  | 1 => (continuousMultilinearCurryFin1 ℝ E (E × (E →L[ℝ] ℝ))).symm
      ((ContinuousLinearMap.id ℝ E).prod (slot1L z.2.1))
  | 2 => (0 : E [×2]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 2 E ℝ z.2.2.1)
  | 3 => (0 : E [×3]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ z.2.2.2)
  | _ => 0

theorem iteratedFDeriv_id_eq_zero {j : ℕ} (hj : 2 ≤ j) (y : E) :
    iteratedFDeriv ℝ j (fun w : E => w) y = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, j = k + 1 := ⟨j - 1, by omega⟩
  rw [iteratedFDeriv_succ_eq_comp_right]
  simp only [Function.comp_apply, fderiv_id']
  rw [iteratedFDeriv_const_of_ne (by omega)]
  simp

theorem iteratedFDeriv_fderiv_eq_curryRight (u : E → ℝ) (j : ℕ) (y : E) :
    iteratedFDeriv ℝ j (fderiv ℝ u) y =
      continuousMultilinearCurryRightEquiv' ℝ j E ℝ (iteratedFDeriv ℝ (j + 1) u y) := by
  rw [iteratedFDeriv_succ_eq_comp_right]
  simp

/-- `gSeries` is the Taylor series of the graph map. -/
theorem iteratedFDeriv_graphMap {u : E → ℝ} {U : Set E} (hU : IsOpen U) (hu : ContDiffOn ℝ ∞ u U)
    {y : E} (hy : y ∈ U) {j : ℕ} (hj0 : 0 < j) (hj3 : j ≤ 3) :
    iteratedFDeriv ℝ j (graphMap u) y = gSeries y (jetOf u y) j := by
  have huy := hu.contDiffAt (hU.mem_nhds hy)
  have hDu := huy.fderiv_right (m := ∞) le_rfl
  have hprod : iteratedFDeriv ℝ j (graphMap u) y =
      (iteratedFDeriv ℝ j (fun w : E => w) y).prod (iteratedFDeriv ℝ j (fderiv ℝ u) y) :=
    iteratedFDeriv_prodMk contDiffAt_id hDu (by exact_mod_cast le_top)
  interval_cases j
  · refine ContinuousMultilinearMap.ext fun m => ?_
    rw [iteratedFDeriv_one_apply, (hasFDerivAt_graphMap (huy.of_le two_le_infty)).fderiv]
    simp only [gSeries, continuousMultilinearCurryFin1_symm_apply, ContinuousLinearMap.prod_apply,
      ContinuousLinearMap.id_apply]
    refine Prod.ext rfl ?_
    ext b
    simp only [jetOf]
    rw [slot1L_apply, iteratedFDeriv_two_apply]
    rfl
  · rw [hprod, iteratedFDeriv_id_eq_zero le_rfl, iteratedFDeriv_fderiv_eq_curryRight]
    rfl
  · rw [hprod, iteratedFDeriv_id_eq_zero (by norm_num), iteratedFDeriv_fderiv_eq_curryRight]
    rfl

/-! ### The jet function -/

variable (J : V →L[ℝ] V) (Ψ : E × (E →L[ℝ] ℝ) → V)

/-- `D^k F_u(y)` from the jet, by Faà di Bruno. -/
def TkOf (k : ℕ) (y : E) (z : SJet E) : E [×k]→L[ℝ] V :=
  (ftaylorSeries ℝ Ψ (y, z.1)).taylorComp (gSeries y z) k

/-- `DF_u(y)` from the jet. -/
def LOf (y : E) (z : SJet E) : E →L[ℝ] V := continuousMultilinearCurryFin1 ℝ E V (TkOf Ψ 1 y z)

/-- **The scalar operator as a function of the 4-jet:** `Pjet(y, z) = −PotM(L, T₂, T₃)`. -/
def Pjet (p : E × SJet E) : ℝ := -PotMP J (LOf Ψ p.1 p.2, TkOf Ψ 2 p.1 p.2, TkOf Ψ 3 p.1 p.2)

variable {J Ψ}

theorem iteratedFDeriv_graphFu_eq_TkOf {O : Set (E × (E →L[ℝ] ℝ))} (hO : IsOpen O)
    (hΨ : ContDiffOn ℝ ∞ Ψ O) {u : E → ℝ} {U : Set E} (hU : IsOpen U) (hu : ContDiffOn ℝ ∞ u U)
    (hgr : ∀ w ∈ U, graphMap u w ∈ O) {y : E} (hy : y ∈ U) {k : ℕ} (hk : k ≤ 3) :
    iteratedFDeriv ℝ k (graphFu Ψ u) y = TkOf Ψ k y (jetOf u y) := by
  have hg : ContDiffAt ℝ ∞ (graphMap u) y :=
    (contDiffOn_id.prodMk (hu.fderiv_of_isOpen hU (m := ∞) le_rfl)).contDiffAt (hU.mem_nhds hy)
  rw [graphFu, iteratedFDeriv_comp (hΨ.contDiffAt (hO.mem_nhds (hgr y hy))) hg
    (by exact_mod_cast le_top), TkOf]
  refine taylorComp_congr (fun j _ => rfl) fun j hj0 hjk => ?_
  exact iteratedFDeriv_graphMap hU hu hy hj0 (hjk.trans hk)

/-- **`−div(JH)[F_u](y) = Pjet(y, J⁴u(y))`.** -/
theorem neg_flowPot_graphFu_eq_Pjet {O : Set (E × (E →L[ℝ] ℝ))} (hO : IsOpen O)
    (hΨ : ContDiffOn ℝ ∞ Ψ O)
    (hs : ∀ q ∈ O, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {u : E → ℝ} {U : Set E} (hU : IsOpen U) (hu : ContDiffOn ℝ ∞ u U)
    (hgr : ∀ w ∈ U, graphMap u w ∈ O) {y : E} (hy : y ∈ U) :
    -flowPot J (graphFu Ψ u) y = Pjet J Ψ (y, jetOf u y) := by
  have hX : ContDiffOn ℝ ∞ (graphFu Ψ u) U :=
    hΨ.comp (contDiffOn_id.prodMk (hu.fderiv_of_isOpen hU (m := ∞) le_rfl)) hgr
  have hinj : ∀ w ∈ U, Function.Injective (fderiv ℝ (graphFu Ψ u) w) := fun w hw =>
    injective_fderiv_graphFu ((hΨ.contDiffAt (hO.mem_nhds (hgr w hw))).differentiableAt (by simp))
      ((hu.contDiffAt (hU.mem_nhds hw)).of_le two_le_infty) (hs _ (hgr w hw))
  have hL : fderiv ℝ (graphFu Ψ u) y = LOf Ψ y (jetOf u y) := by
    ext a
    rw [LOf, ← iteratedFDeriv_graphFu_eq_TkOf hO hΨ hU hu hgr hy (by norm_num),
      continuousMultilinearCurryFin1_apply, iteratedFDeriv_one_apply]
    simp [Fin.snoc]
  rw [flowPot_eq_potM J hU hX hinj hy, hL,
    iteratedFDeriv_graphFu_eq_TkOf hO hΨ hU hu hgr hy (by norm_num),
    iteratedFDeriv_graphFu_eq_TkOf hO hΨ hU hu hgr hy (by norm_num)]
  rfl

/-! ### Smoothness of the jet function -/

omit [FiniteDimensional ℝ E] in
theorem contDiff_slot1L : ContDiff ℝ ∞ (slot1L : (E [×2]→L[ℝ] ℝ) → E →L[ℝ] (E →L[ℝ] ℝ)) := by
  have h1 : ContDiff ℝ ∞ (continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin 2 => E) ℝ) :=
    LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
  exact contDiff_const.clm_comp h1

omit [FiniteDimensional ℝ E] in
theorem contDiff_gSeries (j : ℕ) : ContDiff ℝ ∞ (fun p : E × SJet E => gSeries p.1 p.2 j) := by
  match j with
  | 0 =>
    have he : ContDiff ℝ ∞ (continuousMultilinearCurryFin0 ℝ E (E × (E →L[ℝ] ℝ))).symm :=
      LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
    exact he.comp (contDiff_fst.prodMk (contDiff_fst.comp contDiff_snd))
  | 1 =>
    have he : ContDiff ℝ ∞ (continuousMultilinearCurryFin1 ℝ E (E × (E →L[ℝ] ℝ))).symm :=
      LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
    have hS : ContDiff ℝ ∞ (fun p : E × SJet E =>
        (ContinuousLinearMap.inl ℝ E (E →L[ℝ] ℝ)) +
          (ContinuousLinearMap.inr ℝ E (E →L[ℝ] ℝ)).comp (slot1L p.2.2.1)) :=
      contDiff_const.add (contDiff_const.clm_comp
        (contDiff_slot1L.comp (contDiff_fst.comp (contDiff_snd.comp contDiff_snd))))
    have h := he.comp hS
    convert h using 1
    funext p
    simp only [Function.comp_apply, gSeries]
    congr 1
    ext a <;> simp
  | 2 =>
    have he : ContDiff ℝ ∞ (ContinuousLinearMap.compContinuousMultilinearMapL ℝ (fun _ : Fin 2 => E)
        (E →L[ℝ] ℝ) (E × (E →L[ℝ] ℝ)) (ContinuousLinearMap.inr ℝ E (E →L[ℝ] ℝ))) :=
      ContinuousLinearMap.contDiff (𝕜 := ℝ) _
    have hc : ContDiff ℝ ∞ (continuousMultilinearCurryRightEquiv' ℝ 2 E ℝ) :=
      LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
    have hz : ContDiff ℝ ∞ (fun p : E × SJet E => p.2.2.2.1) :=
      contDiff_fst.comp (contDiff_snd.comp (contDiff_snd.comp contDiff_snd))
    have h := he.comp (hc.comp hz)
    convert h using 1
    funext p
    ext m <;> simp [gSeries]
  | 3 =>
    have he : ContDiff ℝ ∞ (ContinuousLinearMap.compContinuousMultilinearMapL ℝ (fun _ : Fin 3 => E)
        (E →L[ℝ] ℝ) (E × (E →L[ℝ] ℝ)) (ContinuousLinearMap.inr ℝ E (E →L[ℝ] ℝ))) :=
      ContinuousLinearMap.contDiff (𝕜 := ℝ) _
    have hc : ContDiff ℝ ∞ (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ) :=
      LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
    have hz : ContDiff ℝ ∞ (fun p : E × SJet E => p.2.2.2.2) :=
      contDiff_snd.comp (contDiff_snd.comp (contDiff_snd.comp contDiff_snd))
    have h := he.comp (hc.comp hz)
    convert h using 1
    funext p
    ext m <;> simp [gSeries]
  | n + 4 => exact contDiff_const

theorem contDiffAt_TkOf {O : Set (E × (E →L[ℝ] ℝ))} (hO : IsOpen O) (hΨ : ContDiffOn ℝ ∞ Ψ O)
    {p : E × SJet E} (hp : (p.1, p.2.1) ∈ O) (k : ℕ) :
    ContDiffAt ℝ ∞ (fun p : E × SJet E => TkOf Ψ k p.1 p.2) p := by
  refine contDiffAt_taylorComp (fun j _ => ?_) (fun j _ => (contDiff_gSeries j).contDiffAt)
  exact ((hΨ.contDiffAt (hO.mem_nhds hp)).iteratedFDeriv_right (m := ∞) (i := j)
    (by norm_cast)).comp p (contDiffAt_fst.prodMk (contDiffAt_fst.comp p contDiffAt_snd))

/-- **The jet function is smooth** where `DF` is injective. -/
theorem contDiffAt_Pjet {O : Set (E × (E →L[ℝ] ℝ))} (hO : IsOpen O) (hΨ : ContDiffOn ℝ ∞ Ψ O)
    {p : E × SJet E} (hp : (p.1, p.2.1) ∈ O) (hinj : Function.Injective (LOf Ψ p.1 p.2)) :
    ContDiffAt ℝ ∞ (Pjet J Ψ) p := by
  have he : ContDiff ℝ ∞ (continuousMultilinearCurryFin1 ℝ E V) :=
    LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
  have hL : ContDiffAt ℝ ∞ (fun p : E × SJet E => LOf Ψ p.1 p.2) p :=
    he.contDiffAt.comp p (contDiffAt_TkOf hO hΨ hp 1)
  have hjet := hL.prodMk ((contDiffAt_TkOf hO hΨ hp 2).prodMk (contDiffAt_TkOf hO hΨ hp 3))
  have h := ContDiffAt.comp (g := PotMP J)
    (f := fun p : E × SJet E => (LOf Ψ p.1 p.2, TkOf Ψ 2 p.1 p.2, TkOf Ψ 3 p.1 p.2)) p
    (contDiffAt_PotMP J hinj) hjet
  exact h.neg

end CWFlow
