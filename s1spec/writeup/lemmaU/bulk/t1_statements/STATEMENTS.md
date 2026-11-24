# Bulk task 1: statements for Lemma U (summary)

Fable, Nov 13 2026. Full text: `STATEMENTS.tex` (not compiled; braces, environments, `$`, `\[ \]`, labels checked by script). Contracts for tasks 2–5: `CONTRACTS.md`. Hypothesis (O) wording: `HYPOTHESIS_O.md`.
Labels: **[cert]** arb project script; **[est]** float; **[rev]** kill-check 2 scratch; **[paper]** proved by hand here, not scripted.

## Decision

**Global box (0.15, 0.52)**, U_max 0.09, sup|v| ≤ 0.16, design data56: Ĝ₂₂ = 338258.8, g = 3.503e5, rate margin ×29.07 [cert]. Fallback (0.20, 0.60): rate ×16.65 [cert], about 1% of M1 [rev], needs K rerun for g = 6.29e5. Task 5 decides by certifying b₀ = sup|f|, b₁ = sup|r ∂_r f|, b_z = sup|f − U₁| on 190λ₁ ≤ r ≤ 379.5λ₁ from D5a.

## 1. Lemma D in window form (§3)

- **Statement.** For k = 0, 1 and every window [τ_a, τ_a + T] in the bootstrap interval: ∬ over the window and y ∈ [28, 350] of εω (z^A_{k+2})² ≤ C0_k + C1_k·T, with C0_k = (Ê^A_k + Dsrc_k)/(r_k a_c) and C1_k = (2F_k²/γ′_k + coup_k + ô_k + ρ_J·Dsrc_k)/(r_k a_c).
- **Hypotheses.** (D1) (A)'s level-k inequality with retained dissipation r_k = 1.249; (D2) the impulse sources obey the window bound Dsrc_k(1 + ρ_J T); (D3) E^A_k ≤ Ê^A_k.
- **Summation by parts.** Layer-cake lemma: for a nonincreasing m, ∫ m·x ≤ m(0)·C0 + C1·∫m. Applied with m̄(x) = sup over the ramp of W(ln(y/y_c) − x)·D·φ·y/28, which is nonincreasing because W increases on σ ≤ ½ ln 1.5.
- **Data term.** 𝓓_k = S_X·C0_k + C1_k·Σ_X replaces S_X·D_k in Prop. U, with S_X = m̄(0) and Σ_X = (1/c⁻)∫₀^∞ m̄. At y_c = 56: S_X = 1.522, Σ_X = 2.99 [est]. C0 = (2.925e-7, 1.5104e-4), C1 = (5.06e-9, 5.85e-9) [cert, at the old placeholders].
- **H-J-W.** ∫ over every window of V_m² ≤ J_m + s_m T. Sufficient: (W1) L_up_k ≤ Ly1 + 2P_k on every march cell; (W2) the tail bound; (W3) the floor term is accounted on march cells. If (W1) fails, J_m(1 + κ_J) with κ_J ≤ 0.13 [est].
- `U_ramp.tex` l.35, 148, 158 are withdrawn.

## 2. R5 in L² form (§4, Proposition R5′)

- **Structure.** No cut of the identities at y = 350, so no pointwise boundary terms. [350, 700] is a core segment with box (½, 1); its base parts sit inside the core rate (segment condition certified today: `weights_rep.log`). Only terms carrying a derivative of the outer cutoff remain.
- **Result.** Outer source ô_j = Σ_{k ≤ 2} κ^o_jk·N̄_k² per unit τ; enters as ô_j/γ′_j in Ê^A and as ô_k/(r_k a_c) in C1.
- **Consumes.** H-O1: box (½, 1) on [350, 700], from B3U where r < 190λ₁ and from D5a (Z2) beyond. H-O2′: N̄_k² ≥ sup_τ ∫ over [350, 700] of (z^A_k)², k = 0, 1, 2, with N̄_k² = 2(√(2B̂_k/w_min) + p^k·1.41e-5·√ln 1.5)² + n_III_k: Lemma U's level-0–2 energies at the hypothesis level on two windows, plus D5a (Z3). The outer dissipation absorbs z_{j+2}.
- **Tools** (Lemma 11): box for z₀, z₁; energy data; Young against the outer dissipation; one integration by parts for m ≤ 2; interpolation with the C^{1,1} weight εω·t²(1−t)³; an L³ identity ∫ϑ|z₂|³ ≤ Z₁(∫|ϑ′|z₂² + 2∫ϑ|z₂||z₃|); the (2,2,3) handler. No pointwise z₂, z₃; no z₃ or z₄ datum.
- **Closure rule.** The enumerator must assert that every class closes. Fallback (a design change, to be approved): Lemma U's level-1 dissipation as an extra B3U hypothesis.
- **Size.** 1e-6 to 1e-4 on E₂ᴬ [est]; trigger 0.01.

## 3. The bootstrap (§5)

- **Set.** B1; B2′ (core box, E^A ≤ 2Ê^A on [14, 700]); B3U: (a) 𝔈^ρ_j ≤ 2B̂_j, (b) box (0.15, 0.52) and (c) |v| ≤ 0.16 on {y ≥ 56, r < 190λ₁}; B4.
- **Order.** verifier (B1, B2′ → J's in window form, ā, b̄, β̄) → D5a (B4, hand-off, D1 → zone-III box, n_III, Δ^R) → (A) (D1, (H2), (H3), boxes, H-O from Lemma N at 2B̂ → Ê^A, γ′, core box, C0, C1) → U-crude (Ê^A, C0, C1, boxes, neighbours at 2B̂, Δ^R → B̂, B3U strictly) → U-fine → D5b/H6 → D6.
- **H-O is derived** (Lemma 10), not assumed. (A) uses U only at the hypothesis level; U uses the derived (A).
- **Closure** (Proposition 16): any pair (Ê^A, B̂) with 𝒜(2) ≤ Ê^A and 𝒰(2) ≤ B̂ plus the box and rate conditions closes the argument; 𝒜 is affine in N̄².
- **Loop gain.** G_j = (𝒰_j(2) − 𝒰_j(0))/B̂_j. Neighbour launch data: 2q = 0.037; through (A): ≤ 1e-3. Total ≤ 0.05 [est]; acceptance ≤ 0.5.

## 4. Box (§6)

See "Decision". Zone III must deliver Z1 (the deciding box on [190, 379.5]λ₁), Z2 (box (½, 1) out to 700λ₁), Z3 (n_III_k), Z4 (Δ^R, cone family). The same box becomes (A)'s coefficient box on (60, 350].

## 5. Hypothesis (O) (§7, `HYPOTHESIS_O.md`)

(O1) |f| ≤ 0.11, |y f_y| ≤ 0.52 on 60 < y ≤ 350; (O2) |f| ≤ ½, |y f_y| ≤ 1 on [350, 700]; (O3) |z_k| ≤ 1 there, k = 0..3. The (F) macro already postulates the conclusion of Lemma S(A) in (F2), so the accurate conditionality is "(F1) and (O)".

## New findings of this task [paper]

| # | finding | effect [est] | who checks |
|---|---|---|---|
| C5 | The per-centre frozen profile U_ρ is unbounded on the far tail of late centres (r/ρ > 3e3), where the certified U_max = 0.09 fails. Use U₁ for all centres. | none on margins; removes Δ_j | tasks 3, 4 |
| C6 | The cutoffs are C²: (η²)‴ jumps by 360.3 where η = 1. Integrating the m = 3 monomial of pair (j, j+1) by parts leaves +g_c·z_j(28)², g_c ≤ 5.86e-4. Level 2 of the repair is correct (net −g_c); levels 0, 1 miss it; v3.1 misses it at all levels. | about 2% of E₀ᴬ, E₁ᴬ | task 2 (A3) |
| C7 | (A)'s initial datum is taken as (H2), which covers [14, 60]; the energy lives on [14, 700]. Honest: H2 + 5·H3 + zone-III part. | E₀ᴬ ×6, E₁ᴬ ×1.5, E₂ᴬ unchanged; late M1 ×2.28 → about ×1.8–2.0 | task 2 (A4) |
| C8 | Global identities with region-dependent boxes need Lemma L. At y = 350 it is free; at y = 28 it may add about 1e4 to the core Ĝ₂₂. | core Ĝ₂₂ margin ×1.94 → about ×1.7 | task 2 (A5) |
| C9 | (A)'s segments on (60, 350] were certified with M₀ = 0.11; the far-field box is 0.15 (0.236 in v3.2's per-window boxes). | segment checks to redo; large room | task 2 (A5) |
