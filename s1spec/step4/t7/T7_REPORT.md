# T7 report: certificates and constants for P2 (production, base `t1/polished_v2`)

Everything P2 needs from T7 is in `t7/constants.json`. It holds 26 certificate intervals, plus one entry per slab for slabs 0–180 of `t1/polished_v2` (τ ∈ [0, 27.15], k = 26.25 → 2.37e4; τ = 0 is P1 slab 95's centre). All values are certified bounds unless an entry is marked otherwise.
- **Base: polished_v2 with the (SE) correction** of STEP4_RESULT §2(ii). The correction changes the base only on x ≥ 0.6ℓ̃, so every tip-region constant is that of polished_v2. The far constants are recomputed for the corrected base (§9).
- The work directory is `run_v2se/` (far7, nl7 and the farg increments; the tip-region units are copied from `run_v2/`, which holds the polished-base run). `constants.json` is published from `run_v2se/`.
- Gap 9 (tip box |e″| ≤ 1e-2) is closed by Lemma 4a: the tip nonlinear constants no longer assume any bound on e″.
- Review A3 and A2 (Nov 11): **Lemma 4b** certifies the tip box from Λ_e (|e|, |e′| ≤ 1.636e4·Λ_e^{1/2} on [0, 34.25]). **Lemma 4c** puts the nonlinear terms of (Y), with 𝒬 and the nonlinear/E_ij part of (M), into the certified closure. ϱ_c is now certified in nonlinear mode: min 0.04711 on slabs 2–180, losses ≤ 8e-5. The reserve and its split are unchanged (κ_D = ½), and DN = D4CUB = 0.
- The v1 results (store `t1/polished`, slabs 5–187) are archived as `constants_v1.json` and `t8req_v1.json`. §10 summarises them.

**One command** (restartable; units whose JSON exists are skipped):
`python run_all.py --base ../t1/polished_v2 --work run_v2 --jobs 2 --publish --t8` (polished base), then
`python run_all.py --base ../t1/polished_v2 --work run_v2se --jobs 2 --se ../p2 --prev run_v2 --steps sewidths,far,fargse,nl,merge --publish --nlclosure ../p2/t8se/constants/n3.json` (corrected base; the tip-region units `out`, `out_c2r3`, `out_tip` and `out_farg` are copied from `run_v2`; `--nlclosure` adds Lemma 4c's nonlinear-mode closure, §4c).
- Slabs default to the whole store (its `manifest.json`); `--slabs a-b` restricts them (the v1 store used 5–187).
- The steps are tip, far, cert, c2r3, farg, fargref, nl, merge (with `closure_cert.py`) and t8. Cost and timing are in §7.

**Base and fixed functions.**
- The base follows T6's definition. `t6base.py`, `tp6.py` and `t6engine.py` are copied from `t6/`; only the import names were changed.
- Every quantity is formed as a polynomial in t_c first, and only then evaluated at a t_c-ball. Evaluating the DOFs at the ball first loses a factor of about 600 through cancellations in the Hermite basis.
- **φ₁, φ₂ are fixed and stored; they do not depend on the base.**
  - They are χ₃₂p_L and χ₃₂n_L, where κ_L is T1's polynomial on the v1 slab-187 mesh, with DOFs equal to the exact midpoints of slab 187, stage 8.
  - The DOFs and the mesh are in `phi_store.json`, read by `base7.late_base`. κ_L is exactly C⁴.
- φ₃ = χ₃₂ψ̂₃ (`c4fun.PSI`). T9 must use `base7.late_base` and `c4fun`.
- polished_v2 has the same 41 tip knots as v1, so the lossless coarse-mesh cache `fixed_coarse.json` stays valid. It is keyed by the tip knots.

## 1. (C1) + (C2): 26 certificates covering τ ∈ [0, 27.15]

For every interval I and every τ ∈ I, with θ_I = 0.025:

  𝔅_τ(u) ≤ −c_I‖u‖²_ρ − θ_I‖u''‖²_ρ  for every even u ∈ H⁴_ρ with supp u ⊂ [0, 34.25] and ⟨u, φ_j⟩_ρ = 0.

**Plan** (`base7.schedule`, on the store's own τ axis):
- Intervals are groups of consecutive slabs. The group size is 2, 4, 10 or 20, chosen by k at the group start, with thresholds 40, 185 and 2500.
- τ_c is the stored start of the middle slab.
- On v2 the groups are: slabs 0–11 in pairs; 12–51 in fours; 52–121 in tens; then 122–141, 142–161 and 162–180.
- On the v1 store with slabs 5–187, the same rule reproduces the v1 plan exactly.

**Method.** One certificate is run at τ_c (the prototype's: Hermite block, Haynsworth inertia, tail, leak). It certifies Q_c − Σ_e(μ0_e‖w‖²_e + μ1_e‖w'‖²_e) − c_R‖w‖² ≥ 0, which covers every τ ∈ I as follows.
- **Local part.** C₂₂, C₂₁ and C₁₁ do not depend on τ. Integrating by parts (the boundary terms vanish by symmetry), the change of the local form is ∫[−2δC₂₀w'² + (δC₂₀'' − δC₁₀' + δC₀₀)w²]. μ1_e and μ0_e are the sups of these coefficients over each element and all of I, using 2 t_c-balls per slab.
- **Rank-3 part.** c_R is reserved as 1.6·n_side·(F_c e^{0.3H} L + 10⁻³).
  - F_c is the exact operator norm of dΠ/dτ on the constraint space Z at τ_c. It uses invariance under g → Tg: dΠ/dτ = Σδ_i m_i + Σâ_i l_{τ,i}, with δ = a_τ − M_τM⁻¹a.
  - The norm is taken against projected 6×6 Gram matrices, with the directions normalised.
  - `c2r3.py` certifies, per slab, V = ∫‖dΠ/dτ‖_Z dτ (exact rates on a t_c-grid plus a crude second-derivative term).
  - `merge.py` checks Σ V ≤ c_R on each side of τ_c.

| ids | slabs per interval | k_c | c_I | c(τ_c) | c_R | max μ0 |
|---|---|---|---|---|---|---|
| 0–5 | 2 | 27–40 | 0.0629–0.0660 | 0.0713–0.0730 | 0.0052–0.0064 | 7.6e-3 |
| 6–15 | 4 | 45–174 | 0.0599–0.0654 | 0.0700–0.0724 | 0.0050–0.0101 | 9.5e-3 |
| 16–22 | 10 | 226–2151 | 0.0576–0.0614 | 0.0694–0.0700 | 0.0086–0.0118 | 5.3e-3 |
| 23–25 | 20 / 20 / 19 | 3.8e3–1.6e4 | 0.0528–0.0530 | 0.0700 | 0.0170–0.0172 | 6.6e-4 |

- **min c_I = 0.05284** (id 23, k_c = 3776), against ω_eff = 0.045 needed by the closure.
- **Rank-3:** all 26 reserves are verified, with V ≤ 0.66 c_R on every side.
- **T6 residual corrections** (Ṽ = V_s + r in g₃; b₁ += ∫κ̃r) are now evaluated by `merge.py` from `t6/residuals.json` (v2). The cost per slab is 2m̄C_{x,full}‖χr‖_ρ/n_min + ½ sup_{[0,34.25]}(|κ̃r| + |φ'∫κ̃r|), and the cost per interval is the max over its slabs.
  - The cost is at most 9.7e-9; the smallest c_I minus its cost is 0.052840.

## 2. λ < 2μ1: certified from k = 26.25, margin 19.0%
`far7.py` computes the Gårding function of the (I0) operator on the whole far curve [10.25, ℓ̃], for every slab, with 4 t_c-balls and adaptive cells.
- The weight is W = (1+X²)⁻², where X is T8's similarity flattening. The term +½ψ_τ (the τ-derivative of the flattened weight) is essential at the far end.
- The raw infimum of μ − Hardy is 0.1359, at k = 26.3 and x = 10.25. The Hardy term is ≤ 1.8e-4.
- **Split at μ_f = 0.15:** the zone below the floor is x ≤ 10.52, where χ = 1, and a_μ ≤ 2.5e-13 is added to a_f.
- So μ1 = 0.15 on all 181 slabs, and μ2 ≥ 0.643.
- λ_max = 0.2520 is the max of T8's certified rates over slabs and levels. It covers both T8 runs: v2 (0.25201) and the assembly's `p2/t8se` on the corrected base (0.25194). So 2μ1/λ_max − 1 = **0.190**.
- These values hold for the (SE)-corrected base (§9). The raw infimum changes only in the 14th digit.

## 3. (C3)
- c_φ = (1.895, 1.534, 0.318) for normalised φ_j.
- At τ_c, from k = 27 to 1.6e4: m̄ 5.06 → 5.39, ℓ̄ 9.56 → 9.41, C_x on Z 7.07 → 6.33, n_V = ‖χV_s‖ from 0.167 down to 3.4e-4.
- Per slab (`c2r3`): m̄ ≤ 5.40, ℓ̄ ≤ 9.57, C_x ≤ 7.63, C_{x,full} ≤ 7.82, and n_min ≥ 2.33e-4 (at slab 180).

## 4. (C4): max over the 181 slabs
**E2 level (ρ, [0, 34.25]):** G₂₂ 15.63, G₁₂ 19.37, G₀₂ 12.74, all maximal at k = 26.3. Modulation is in AB form: G_x replaces G_m A² by G_x AB.

**Strip, σ ∈ [29, 32]:**
- s_B 3.47e-5, s_A 2.18e-4, s_A0 9.31e-7.
- The ε-level constants d₃, d₂, d₁, d₀ are 7.08e-3, 1.80e-2, 4.64e-2, 5.28e-2.

**ζ-ramp:**
- a_f 0.0696 (split included), b_f 0.0394;
- c₃..c₀ = 0.448, 0.539, 0.449, 0.425;
- Ψ₁ = 0.3885, c_{1W} = 0.0943.

**Modulation.** `farg7.py` computes the W-Gram matrices of ζĝ' and (ζĝ')'' and the ρ-Gram of (χĝ)''. merge.py then sets C_γ0, C_γ2, G_x = min(C_x√tr Γ, Σ_i C_{x,i}‖f_i‖).
- Final values: **C_γ0 ≤ 0.207** (slab 152), C_γ2 ≤ 0.0159, G_x ≤ 14.06.
- **The sawtooth in √Γ̂₀ is resolved.** With one t_c-ball per slab, √Γ̂₀ (the clock direction) rises over the 5–8 slabs before each outer-element insertion, up to 0.328 at slab 152, and then drops at the insertion. The cause is the t_c-hull of the moving outer knots. v1's "spikes" at slabs 132 and 160 were the same effect.
- The new step `fargref` (`farg7.py --refine 0.1 --nsub 4`) recomputes, with 4 t_c-balls, the 34 slabs whose √Γ̂₀ exceeds 0.1. For each entry it keeps the sharper of the two valid bounds.
- After refinement the max √Γ̂₀ is 0.147, and C_γ0 drops from 0.460 to 0.207. Most refined slabs fall to 0.05–0.07.

**Nonlinear terms** (`nl7.py`, all 181 slabs): n₂–n₄ are Hessian sups of the exact graph operator F.
- **Tip, x ∈ [0, 34.25]: box on (e, e′) only, no condition on e″** (Lemma 4a). The constants are n₄ 6.17, n₃* 9.26, n₂* 16.9 and m₂₃ 3.0003, all maximal at k ≈ 26.
  - These replace the earlier box values |e″| ≤ 1e-2 (n₃ 9.36, n₂ 18.7). On every slab, n₃* and n₂* are at most the old values.
- **Far boxes**, given per dyadic x-range up to ℓ̃. They are recomputed on the (SE)-corrected base (§9). ℓ̃ grows like k, reaching 6.7e5 at slab 180.
  - n₄ ≤ 7.8 in the T5 box and ≤ 0.42 in the sharp box (|e| ≤ 1e-3); n₃ ≤ 0.65.
  - On the outermost ranges (the bend, where κ̃e is not small), n₂ grows roughly linearly in x. In the T5 box it is 6.4 at slab 0, 193 at slab 90 and 5.6e3 at slab 180, giving x n₃ ≤ 1.0e5 and x²n₂ ≤ 6.2e14.
  - The sharp box gives 4.0e3 and 5.0e14. T8's Lemma N3 must take these per x-range.

## 4a. Lemma (tip nonlinearity without a pointwise bound on e″; Gap 9)
Here Ω = [0, 34.25] = supp χ₃₂ and w = |e| + |e′|. The graph operator is F = A₄e⁗ + A₃e‴ + A₂ (`nl7.derive`), and N = F(e) − F(0) − DF(0)e. Norms are those of T5 §1.1: B = ‖u″‖_ρ, W₃ = ‖u‴‖_ρ, D₄ = ‖u⁗‖_ρ, V_j = ‖v^{(j)}‖_W, with u = χe and v = ζe′.

**(i) Structure.** α = 1 − κ̃e and J = (α² + e′²)^{1/2} contain no e″. So F is exactly a polynomial in e″. All of the following are checked symbolically in `nl7.derive` (`poly_check`):
- A₄ = A₄(e, e′);
- A₃ = A₃₀(e, e′) + A₃₁(e, e′)e″, with A₃₁ = 10e′/J⁶;
- A₂ = A₂₀ + A₂₁e″ + A₂₂e″² + A₂₃e″³, with A₂₂ = −3κ̃ at e = e′ = 0 and A₂₃ = 3(α² − 5e′²)/J⁸.

**(ii) Bound with no e″ box.** Suppose |e|, |e′| ≤ η₁ = 1e-5 on Ω. Then, for all e″, e‴, e⁗,

  |N| ≤ n₄w|e⁗| + n₃*(w + |e″|)|e‴| + n₂*(w + |e″|)² + m₂₃|e″|³,

where
- n₄ = sup|∇A₄|;
- n₃* = max(sup|∇A₃₀|, sup|A₃₁|);
- n₂* = max(½sup|∇²A₂₀|, ½sup|∇A₂₁|, sup|A₂₂|);
- m₂₃ = sup|A₂₃|.

The gradients are in (e, e′). The sups are arb bounds over the cells of Ω, the slab, and the (e, e′)-box only.

*Proof.* Write N = [A₄ − A₄(0)]e⁗ + [A₃₀ − A₃₀(0) + A₃₁e″]e‴ + [A₂₀ − A₂₀(0) − ∇A₂₀(0)·(e, e′)] + [A₂₁ − A₂₁(0)]e″ + A₂₂e″² + A₂₃e″³. Apply the mean-value and Taylor theorems in (e, e′) to each bracket.

**(iii) Control of S := sup_Ω|e″|.**

  S² ≤ S_W² + 2K_reg·B·W₃,  with S_W² := (V₁²/18 + 2V₁V₂)/W(34.25), K_reg = 1/ρ(16.25) ≤ 3.021e4, W(34.25) ≥ 7.25e-7.

*Proof.*
- On [16.25, 34.25] we have ζ = 1, so e″ = v′. The 1-D Sobolev bound f(x)² ≤ L⁻¹∫f² + 2∫|ff′| on this interval (L = 18), with W ≥ W(34.25), gives |e″| ≤ S_W.
- On [0, 16.25] we have χ = 1 and ρ ≥ 1/K_reg. Hence e″(x)² ≤ e″(16.25)² + 2∫_x^{16.25}|e″e‴| ≤ S_W² + 2K_reg B W₃.

S_W is controlled by Λ_e alone, at each time. The D₄ hidden in W₃ ≤ (½ + √c₁)B + D₄/2 is dissipation, so S has no pointwise-in-time bound. None is needed.

**Consequence for T5 §2.5.**
- T5's estimates use the box only through the constants n_j; their proofs already bound e″ by S (Agmon) and by L² norms. With (ii) they hold verbatim, with n₃ → n₃* and n₂ → n₂*.
- The one new term m₂₃|e″|³ ≤ m₂₃S e″² adds:
  - to (E0): at most η₁m₂₃S(B² + σ²V₁²), with σ² = sup_strip ρ/W ≤ 2.8e-5 (on the strip, e″ = v′);
  - to (E2): at most εm₂₃S²(B + σV₁)(D₄ + 2Φ₁W₃ + Φ₂B).
- With (iii) and Young (weights 1e-9; `merge.tip_lemma`), both are at most δΛ_e + βεD₄².
- **At Ȳ = 1e-12, with Lemma 4b's box** (|e|, |e′| ≤ c_e1Ȳ = 1.64e-8, Λ_e ≤ Ȳ²): S_W ≤ 1.14e-6, δ ≤ 6.1e-9 and β ≤ 3.0e-9 (`merge.tip_lemma`, meta.nonlinear_tip).
  - The first version used T5's 2e3Ȳ/√μ_e = 1.05e-8 (review A3).
  - **Superseded by Lemma 4c.** The cubic term is now one of the closure's nonlinear terms, so DN and D4CUB are 0.
- **Far box.** S_W ≤ 1.14e-6 at every time also gives the far box's e″ component (x|e″| ≤ η_F) on [16.25, 34.25] from Λ_e alone. On [10.25, 16.25], bound (ii) replaces the far box.
- (C1)–(C4) other than n₂–n₄ are linear constants of the base; they contain no e.

## 4b. Lemma (pointwise tip box from Λ_e; review A3)
Let Λ = Λ_e = A² + εB² + μ_e(V₀² + ε_IV₂²) and Ω = [0, X_E], X_E = 34.249 (the end of supp χ₃₂, with χ₃₂ = 1 up to X₁ = 31.249). At every time,

  |e(x)| ≤ c_e(x)Λ^{1/2}, |e′(x)| ≤ c_e1(x)Λ^{1/2} on Ω, with **sup_Ω c_e ≤ 1.6031e4 and sup_Ω c_e1 ≤ 1.6352e4.**

All constants are arb upper bounds (`nlterms7.lemma4b`). They are published in meta.tip_box, together with the table.

| x | 0 | 10.25 | 16.25 | 24.25 | 28 | 30.25 | 31.25 | 34.249 |
|---|---|---|---|---|---|---|---|---|
| c_e | 5.81 | 84.4 | 890.95 | 6425 | 9510 | 11647 | 12670 | 16031 |
| c_e1 | 92.0 | 1336 | 14105 | 8205 | 10934 | 12759 | 13615 | 16352 |

*Proof.*
- **(i) |e|, x ≤ X₁ (u = e).** ρu²(x) = ∫_x^{X_E}(φ′ρu² − 2ρuu′) ≤ Φ₁A² + 2AW₁. Here Φ₁ = sup_Ω φ′ = 1.0822, and W₁² ≤ AB + c₁A² ≤ s²Λ with s² = 1/(2√ε) + c₁. So c_e(x)² = (Φ₁ + 2s)/ρ(x); for example K_e(16.25) = 890.95.
- **(ii) |e′|, x ≤ X₁.** ρu′²(x) ≤ Φ₁W₁² + 2W₁B ≤ 2C_AΛ. Here C_A = (m₃₁ε^{−3/2} + c₁/(4ε))^{1/2} + (Φ₁/2)(1/(2√ε) + c₁) = 3292.8, where m₃₁ = max cos t sin³t = 3√3/16.
- **(iii) |e′|, x₀ ∈ [16.25, X_E].** Here ζ = 1 and v = e′.
  - Let f = ψ_c(1+x²)^{−1}v, where ψ_c = 1 on [10.25, x₀], ψ_c = cos²(π(x−x₀)/(2ℓ_c)) on [x₀, x₀+ℓ_c], and ψ_c = 0 beyond. Take ℓ_c = 100. Then ψ_c ∈ W^{2,∞}, |ψ_c′| ≤ π/(2ℓ_c) and |ψ_c″| ≤ π²/(2ℓ_c²).
  - The support lies in [10.25, ℓ̃/2], because min ℓ̃ = 411.8 ≥ 268.5. So W = (1+x²)^{−2} exactly there.
  - The sharp line inequality is sup f² ≤ K_line(‖f‖² + ε_I‖f″‖²), with K_line = ε_I^{−1/4}/(2√2) = 4.7773.
  - Also ‖f″‖ ≤ aV₂ + bV₀, with a = 1 + σ₁ + π/(2ℓ_c) = 1.2090 and b = 2σ₁k_{V1} + σ₂ + (π/ℓ_c)(k_{V1} + σ₁) + π²/(2ℓ_c²) = 0.3998. This uses σ₁ = 2x/(1+x²) and σ₂ = (6x²−2)/(1+x²)² at x = 10.25, and V₁ ≤ k_{V1}V₀ + V₂/2.
  - Hence |e′(x₀)| ≤ (1+x₀²)(K_line c_ω/μ_e)^{1/2}Λ^{1/2}, where c_ω = max(1 + ε_Ib²(1+1/δ), a²(1+δ)) = 1.4619 with δ = 1e-4.
  - On [16.25, X₁], c_e1 is the smaller of the two bounds (ii) and (iii).
- **(iv) |e|, x ≥ 16.25.** |e(x)| ≤ |e(16.25)| + H(x)V₀, with H(x)² = ∫_{16.25}^x(1+s²)²ds. This uses W ≥ (1+x²)^{−2}, which holds on the flattened part too, since Φ′ ∈ [0, 1] gives X ≤ x. So c_e(x)² = K_e(16.25)² + H(x)²/μ_e. The same bound holds on all of [16.25, ℓ̃]; Lemma 4c uses it for the far e-terms.

*Consequences.*
- η₁ := c_e1Ȳ ≤ 1e-5 (T5's and T8's tip box) holds for every Ȳ ≤ 6.11e-10. This is asserted. At the bootstrap box Ȳ = 1.7e-12, η₁ = 2.78e-8.
- The old constant 2e3Ȳ/√μ_e = 1.054e4·Ȳ was the value of (iii) at x = 30.25; it is too small on [30.25, 34.25].
- `merge.tip_lemma` (Lemma 4a's numbers) now uses η₁ = c_e1Ȳ. T8's N3′ uses the 1e-5 box, which this lemma covers. `n3far --eta1` may be set to c_e1·Y.

## 4c. Lemma (the nonlinear terms of (Y); review A2, T8's note on (M))
**Setting.** The nonlinear contributions to ½Λ′ ≤ E0 + E2 + μ_e(I0 + I2) are:
- N₀ = ⟨χ(𝒩+𝒬), u⟩_ρ;
- N₂ = ε⟨χ(𝒩+𝒬), D_ρu⟩_ρ;
- N_{I0} = ∫Wζ²e′∂_x(𝒩+𝒬);
- N_{I2} = ε_I∫(Wv″)″ζ∂_x(𝒩+𝒬);
- the part δx of (M) beyond its linear and residual parts. δx enters through the modulation couplings √3|δx̂|A (E0), εg_{2E}|δx̂|B (E2), μ_eg₀|δx̂|V₀ (I0) and μ_eε_Ig₂|δx̂|(V₄ + 2Ψ₁V₃ + 2c_{1W}V₂) (I2).

T8's exact split is used throughout: 𝒩 + 𝒬 = (1−ϑ)𝒩 + ϑ(e⁗ + κ̃²e″ − c₀e + P) + ge′ − (ϑc + x₁)e, with 𝒩 = Ãe⁗ + B̃e‴ + C̃ (LEMMA_F §6.1).

**Box 𝔹_Y** (bootstrap; closed by (Y) itself, by Lemma 4b and by T8's n3boxes):
- Λ_e ≤ Ȳ_box², with Ȳ_box = 1.7e-12, the Y of T8's n3.json run;
- Lemma 4b's tip box;
- T8's far box 𝔹′ on x ≥ 30.249 (T8's core);
- |x̂| ≤ x̄ = C_xȲ_box + x_r + x_NL and |ϑ| ≤ x̄/n_V, per interval.

**Statement.** NL ≤ Σ_k c_k a_kb_k over 55 explicit terms per interval (closure atoms, `nl_terms.json`). With them, closure_cert certifies Q(z) − D_s − NL(z) ≥ ϱ_c^{NL}Λ_e. The reserve D_s is the same as in the linear mode.

**Estimates.** w_m := e^{(m+1)}, T := Wζ², and σ = (ρ(X₁)/W(X_E))^{1/2} = 5.28e-3.
- **Tip (E0, E2).**
  - Use nl7's bound |𝒩| ≤ n₄w|e⁗| + n₃₀w|e‴| + |A₃₁||e″||e‴| + n₂₀w² + n₂₁w|e″| + m₂e″² + m₂₃|e″|³ (Lemma 4a), with A₃₁ = 10e′/J⁶ ≤ 10η₁/α_min⁶.
  - Use |u| ≤ η₁ and ‖D_ρu‖ ≤ D₄ + 2Φ₁W₃ + Φ₂B.
  - The ρ-norms on Ω are N_k ≤ U_k + σV_{k−1} (on the strip, e^{(k)} = v^{(k−1)}), and N₀ ≤ A + c_e(3ρ(X₁))^{1/2}Λ^{1/2}.
  - sup_Ω|e″| ≤ S ≤ S_W + (2K_reg)^{1/2}(B^{3/4}D₄^{1/4} + c₁^{1/4}B), by Lemma 4a(iii).
- **(I0): one global weak identity on [10.25, ℓ̃].** There are no boundary terms: ζ vanishes to order 4 at 10.25, and (SE) holds at ℓ̃.

    N_{I0} = ∫TÃw₂² + ∫(2T′Ã + TÃ′)w₁w₂ + ∫(T″Ã + T′Ã′)w₀w₂ − ∫(T′w₀ + Tw₁)(B̃w₂ + C̃) + ½∫(g′T − gT′)w₀² − (ϑc + x₁)∫Tw₀².

  Its integrands are bounded region by region:
  - **far core:** T8's per-range Ã⁺, |Ã|, |Ã′|, |B̃|, γ and G (𝒬-effective). Here |W′/W| ≤ min(Ψ₁, 8/x) and |W″/W| ≤ min(2c_{1W}, 2264/x²); on the flattened part these use X ≥ x/2 and |Φ″| ≤ 152.25.
  - **[16.25, 30.25]:** analytic Ã, Ã′ (J ≥ α ≥ 1 − K₀η₁), nl7's constants, and |e″| ≤ S_W.
  - **ramp [10.25, 16.25]:** the same with |e″| ≤ S, and jets in ρ-norms (κ_R = W(10.25)/ρ(16.25) = 2.685).
  - **e-terms:** |e(x)| ≤ c_e(x)Λ^{1/2}, by 4b(iv).
- **(I2): strong form for (1−ϑ)𝒩 + ϑ(…) − (ϑc+x₁)e.**
  - The test function is (Wv″)″/W = v⁗ + 2ω₁v‴ + ω₂v″, paired per range with the coefficient. On ζ = 1 the diagonal term ∫WÃ_eff(v⁗)² ≤ Ã⁺V₄² keeps its sign.
  - The transport term uses the global identity ∫Wζ²e‴∂³(ge′) = ½∫[5g′T − gT′]w₂² + 3∫Tw₂g″w₁ + ∫Tw₂g‴w₀, plus the ramp commutators (strong form, ramp only).
  - On [10.25, 30.25], ∂_x𝒩 is a polynomial in e″, …, e⁽⁵⁾ with 12 monomials (`nlsym/rampforms.json`). Its coefficients are enclosed over the (e, e′) box by second-order AD (`nlsym/td2.py`, `nlsym/nlforms_gen.py`): the value for order-0 coefficients, the gradient for order 1, the Hessian for the constant term.
  - On the ramp, the single ζ goes to the highest jet (ζe^{(k)} = v^{(k−1)} minus ζ-commutators), and e″ ≤ S. The e‴² monomial (coefficient A₃₁) uses sup_R|ζe‴|² ≤ (2/W(16.25))‖√Wζe‴‖‖√W(ζe‴)′‖.
  - On [16.25, 34.25], sup|v″|² ≤ (V₂²/18 + 2V₂V₃)/W(34.25).
  - Far core: T8's B̃′, dA3e, dC̃ and β, with T8's e‴ box, which is certified on the core. **No e‴ box is used below x = 30.25.**
- **Reduction.** Products of more than two jets are expanded with W₃ ≤ (BD₄)^{1/2} + c₁^{1/2}B and its analogues. The excess degree is boxed on A, Λ^{1/2}, V₀, B, V₂ at Ȳ_box, and a weighted AM–GM split finishes. **The dissipation degree is at most 7/4 < 2 for every term (asserted).**
- **(M).**
  - Write ψ = ρφ̂_j. The weak form ∫ψ(−(J/α)κ_ss) = −∫((ψ/α)′J^{−1})′κ is exact (`nlgen.py` checks it symbolically) and quadratic in e″. Hence |b̂_NL| ≤ √3 Σ_m S_φm[½H0_m(N₀+N₁)² + L1_m(N₀+N₁)N₂ + V2_mN₂²], plus the P-term.
  - S_φm = sup|s_m/ρ| ≤ (81.0, 124.6, 207.9), from cert7's fixed functions. H0, L1 and V2 are TD2 enclosures on the tip cells.
  - |Ê| comes from Q₁ = q̃e′/α − e, Q₂ = t̃₁e′/α, and Q₃ = coefficient of ϑ over n_V; the ϑe⁗ part goes via c_φ.
  - Then x_NL ≤ m̄|b̂_NL| + m̄‖Ê‖/(1 − m̄‖Ê‖)·(C_xA + x_r + m̄|b̂_NL|).
  - **At the box: x_NL ≤ 1.95e-15 and m̄‖Ê‖ ≤ 1.04e-4.** T8's slack is C_x(1e-12 − Y_cert) ≥ 2.44e-12 with each interval's C_x (2.9e-12 with T8's 7.82). So x_NL is at most 8.0e-4 of the slack.

**Results** (`nlterms7.py`, 35 s; `closure_cert.py --nl`, 2 min):
- **All 25 intervals on slabs 2–180 are certified. min ϱ_c^{NL} = 0.04711** (linear mode: 0.04713); losses are ≤ 8e-5.
- The top terms are μ_eÃ⁺V₂² ≤ 5.8e-4 and μ_eε_IÃ⁺V₄² ≤ 1.7e-8, with T8's sup Ã⁺ ≤ 0.0161 (interval 25). Every other coefficient is ≤ 7.1e-5 (an A² term, interval 1).
- Why the loss is far below the review's 1e-2 estimate: the μ_eV₂² dissipation that remains after the reserve (0.0135) is not what binds the closure, and the top term takes ≤ 4.5% of it. Interval 25 is still certified with all its NL terms ×10 (ϱ = 0.04934) and fails at ×20.
- Interval 0 (slabs 0–1) is outside P2 (SLAB0 = 2; n3.json starts at slab 2). Its rho_c is null in nonlinear mode.
- **Reserve and shares are unchanged.** κ_D = ½ for Lemma J (as in T8). ycomp's shares are θ₀/8, ε/16, μ_e/16, μ_eε_I/16, and **DN = D4CUB = 0**.
- No dissipation slice is used, so the hand-off window (one slab) does not interact with this lemma.
- The (Y) residual level with ϱ_c^{NL} and DN = D4CUB = 0 is 1.3295e-13 (was 1.3281e-13). Y_cert = max(Y(τ_a) ≤ 6.20e-13, level) is unchanged, and so is T8's slack.
- **Published** (constants.json): meta.closure (nonlinear mode, with rho_c per interval, shares, DN, D4CUB and Ȳ_box), meta.closure_linear, meta.nonlinear_closure (per interval: loss, top terms, box, parts and the n-dependence), meta.tip_box and meta.modulation_nl.
- **Rerun:** `run_all.py ... --steps merge --nlclosure ../p2/t8se/constants/n3.json`.

## 5. Closure, ϱ_c and K_e
- `closure_cert.py` at μ_e = 0.036 (AM–GM relaxation plus arb LDLᵀ on 6×6), linear mode: all 26 intervals are certified.
  - min ϱ_c = 0.0471 (id 24, k_c = 7994); the range is 0.0471–0.0634.
  - Before `fargref` it was 0.0338 at the same interval, where C_γ0 was 0.46.
- **Nonlinear mode** (Lemma 4c, the published ϱ_c): all 25 intervals of slabs 2–180 are certified, with **min ϱ_c = 0.04711**.
- **K_e ≤ 890.6**, such that |e(x)| ≤ K_e Ȳ for x ≤ x_a = 14 + X₀. The formula does not depend on the base.

## 6. T8 on the new base (`run_all.py --t8`)
**Changes to T8's code** (cconst/merge/camp):
- The slab range is T6's (0–180).
- camp's window is the trailing run of slabs with k ≥ 2.0e4 at their start: slabs 177–180, k 2.04e4 → 2.37e4. The same rule gives 184–187 on v1.
- camp checks that slabs are contiguous with length 0.15; F₀ is taken at the start of slab 0.
- cconst now takes μ_e, ϱ_c and K_e from `t7/constants.json` when present. K_e = 890.6 replaces the 1.4e4 placeholder.
- The t8 step always recomputes the per-slab constants (about 2 min). Earlier sets are moved aside, and v1's are in `t8/constants/slabs_v1`, `constants_v1.json` and `amp_v1.json`.
- On v1 data, merge and camp reproduce v1's Ȳ to all digits.

**Results:**
- **Ȳ requirement 2.092e-11** (binding item 3, E₃). It was 1.950e-11 before `fargref`; v1 gave 1.964e-11.
- **ϱ_c has no effect on Ȳ.** It enters T8 only through (ϱ_c − λ*/2)₊ in Lemma J's B, and λ*/2 ≥ 0.1196. camp with ϱ_c = 0.0338, 0.0381 and 0.10 gives identical Ȳ.
- The change in Ȳ comes from C_γ0 (about 7% between 0.46 and 0.207).
- Inputs still from v1: BASE_T1 (T1's window energies at k = 2.3e4) and KCH = 1.02 are the v1 T1-report values.

## 7. Batches and costs (v2 run, 2 jobs, nice 10)
| step | CPU per unit | CPU total | wall |
|---|---|---|---|
| tip7, far7 | 1.4 s per slab each | 0.14 h | 4 min |
| cert7 (26 intervals) | 2.1 min | 0.92 h | 30 min |
| c2r3 | 35 s | 1.76 h | 53 min |
| farg7 + fargref (34 slabs × 4 balls) | 28 s + 80 s | 1.43 + 0.76 h | 43 + 24 min |
| nl7 | 39 s | 1.97 h | 60 min |
| merge, closure_cert, T8 (cconst, merge, camp) | — | 0.1 h | 6 min |

The total is about 7.1 CPU-h, or 3.7 h wall.

The corrected-base run (`run_v2se`) took 51 min wall in all:
- far7 --se: 3 min;
- se7 widths and farg increments: seconds;
- nl7 --se: about 40 s per slab, 47 min. The new tip constants dominate the cost. The far cells below x_c are reused from `run_v2` via `--prev`;
- merge and closure_cert: 2 min.

Portability is unchanged: relative paths, `os.path`, and spawn plus `subprocess` via `sys.executable`. Only cert7 needs step 3's `c4fun`; the cache rebuild also needs it, but only for new tip knots.

## 8. Flags
1. The certified claim uses directions V_s. The T6 corrections are covered by the computed sensitivity in §1, which is at most 9.7e-9.
2. T8's BASE_T1 window energies and KCH are still the v1 values; T1 should confirm them for polished_v2.
3. d_j and s_A (§4) are only used at the ε-level.
4. The far-box n₂ grows with x on the bend (§4). T8's Lemma N3 must use the values per range, preferably from the sharp box.
5. The modulation constants are basis-free per component.
6. **`p2/t6/jets` (the assembly's widened T6 jets) are not strictly rigorous.**
   - `se_corr.py` widens the float64 hull ends by adding the width. Where the width is below half an ulp of the end, the addition is lost.
   - Lost widths were found for Q and θ on all far cells, and for k₀, p, n, t₁ and V on many (checked on every 10th slab).
   - Fix: after adding, move each end one ulp outward with `np.nextafter`.
   - This affects T8's `t8se` run, not T7: T7 widens with exact arb balls (§9).
   - Lemma 4c's far-core part (x ≥ 30.25) reads T8's n3.json, which is computed on these jets. After T8's `nextafter` fix and an n3far rerun, repeat `nlterms7.py` and `closure_cert.py --nl` (about 3 min). The effect is at most an ulp of each hull end.
7. **T8's residual part of (M).** n3far/camp take x_r = m̄·max_j c_φ,j·ℛ₂. The bound needs the 2-norm ‖(c_φ,j)_j‖₂, which is ×1.30 larger: x_r ≤ 2.6e-16 instead of 2.0e-16. This is inside the slack; ycomp and nlterms7 use the 2-norm.
8. **Interval 0** (slabs 0–1) has no nonlinear-mode ϱ_c, since n3.json starts at slab 2. P2 starts at slab 2 (SLAB0 = 2), so it is not needed.
9. Lemma 4c uses no e‴ box below x = 30.25. T8's N3′ still assumes x²|e‴| ≤ η_F for y < 28; that is T8's flag, and T7's results do not depend on it.

## 9. Far constants on the (SE)-corrected base (Gap 10)
The (SE) correction changes κ̃ only on u = x/ℓ̃ ∈ [0.6, 0.8]. Beyond 0.8ℓ̃, the integrated quantities shift. `se7.py widths` rebuilds the certified widths per slab with se_corr's formulas, from `p2/out/se_corr.json` (a₁, a₂, a₁′, a₂′), T6's v2 window sups and the bump norms; the result is `se_widths.json`. Max widths over the slabs:

| quantity | max width |
|---|---|
| κ̃ | 3.1e-17 |
| Q | 1.1e-19 |
| θ | 5.7e-16 |
| p, q | 4.9e-13 |
| n, t₁ | 5.7e-16 |
| V | 5.1e-12 |
| ϖ | 5.1e-12 |
| θ_τ (from \|∂_τδθ\| ≤ \|a₁′\| + … with δθ = a₁I_B(u) + (a₂/ℓ)B(u)) | 1.2e-13 |

How each constant is recomputed:
- **far7 (μ₁, μ₂, a_μ, Hardy) and nl7 far boxes:** recomputed with the jets widened by these balls on every cell meeting x ≥ x_c = 0.6ℓ_lo (`--se`).
- **farg Γ₀, Γ₂:** given certified increments (`se7.farg_se`): ∫_{x≥x_c}W̄(2|g^{(m)}|Dg^{(m)} + Dg^{(m)2})/‖χg‖². Here sup|g^{(m)}| comes from T6's hulls by the Frenet recursion, and Dg^{(m)} from the widths. The increments are ≤ 2.6e-19 on a normalised diagonal.
- **Not recomputed:** Γ₂E, the certificates, C3, E2, strip, ramp, the tip n_j and K_e. These live on x ≤ 35 < 0.6ℓ̃, where the corrected base is the polished one.

Results (published):
- μ₁ = 0.15, μ₂ ≥ 0.643, a_μ ≤ 2.5e-13 and Hardy ≤ 1.83e-4. These are unchanged except in the 14th digit.
- C_γ0 ≤ 0.2068, C_γ2 ≤ 0.0159, G_x ≤ 14.06.
- The far n_j change by at most 8.7e-7 relative. They are taken as max(new cells, unchanged cells of `run_v2`), so they can only increase.
- The closure was rerun: **min ϱ_c = 0.0471**, unchanged.
- `constants.json` now also carries, per slab, the raw farg values (`slabs[s].farg`) that the P2 (Y) majorant reads, and a `meta.se_correction` block.

**9a. Direct inputs (review B1, Nov 12).** T8 recomputed the corrected base's residual with T6's engine (`p2/t6c`: slabc.py, residuals.json, jets). se_corr's majorants are retired, so T7 no longer reads `p2/out/se_corr.json` for slabs 2–180.
- **se7 widths.**
  - a₁, a₂, a₁′, a₂′ are now t6c's certified enclosures (`se_direct.sharp`, the max over sub-slabs).
  - The residual widths are direct: |r̂ − r̃| ≤ sup|r̂| + sup|r̃| over the far cells, and the same for A. Here r̂ comes from p2/t6c/jets and r̃ from t6/jets.
  - Slabs 0–1 are outside P2 and outside t6c; they keep the se_corr widths, flagged.
  - `se_widths.json` carries an args block with the paths and sha256 of every input and of the code. The old file is kept as `se_widths_secorr.json`.
- **Comparison.** The direct widths of κ̃^{(m)}, Q, θ, p, q, n, t₁, V* and θ_τ are ≤ the se_corr ones on every slab (max ratio 0.814; 0.623 for θ_τ). The r and A widths are larger, up to ×148, which brings the V and ϖ widths to ≤ 9.9e-11.
- **Effect on the published constants:**
  - **far7** widens only κ̃^{(0..3)} and Q, so `out_far` (computed with the se_corr widths) contains the direct case and stays valid. merge checks this (meta.se_correction.far7_check, ratio 0.814).
  - **farg** increments were recomputed (unchanged, ≤ 2.6e-19).
  - **nl7** far boxes were rerun with the direct widths (`out_nl`; the se_corr run is kept as `out_nl_secorr`). No certified chain uses them: Lemma 4c reads T8's n3, and the tip constants are identical.
- **Lemma 4c with T8's B1 far box** (`p2/t8se/constants/n3_b1.json`) and the direct residuals/jets (`--res p2/t6c/residuals.json --jets p2/t6c/jets`; the jets equal t6/jets for x ≤ 40):
  - **min ϱ_c^{NL} = 0.04711**, all 25 intervals certified, losses ≤ 8e-5;
  - top terms μ_eÃ⁺V₂² ≤ 2.5e-4 (sup Ã⁺ ≤ 6.95e-3);
  - x_NL ≤ 1.95e-15 against a slack ≥ 2.44e-12 (ratio ≤ 8.0e-4).
- **(Y) level** (ycomp with T6RES = p2/t6c/residuals.json, SLAB0 = 2, ϱ_c^{NL}, DN = D4CUB = 0): residual level 3.21e-14, and Y(τ_a) = 6.21e-13, so Y_cert = 6.21e-13 ≤ 6.3e-13.
- **Provenance.** Every output now has an args block with its input paths and sha256 (input files, the jets directories as one fingerprint of their per-file hashes, and the code): nl_terms.json, closure_cert(_nl).json, se_widths.json, and meta.closure / nonlinear_closure / tip_box / modulation_nl / se_correction. The helpers are in `provenance.py`.

## 10. v1 (superseded)
- The store `t1/polished` (slabs 5–187, τ ∈ [0.75, 28.2]) does not pass through P1's centre: it is off by 2.8e-9. Everything on it was certified:
  - 27 intervals, with min c_I 0.0528;
  - min ϱ_c 0.0381, C_γ0 0.389;
  - T8's Ȳ 1.964e-11.
- These results are archived in `constants_v1.json` and `t8/constants/*_v1`.
- Carry-over by perturbation would have cost as much as re-certifying, so the pipeline was rerun on v2.
