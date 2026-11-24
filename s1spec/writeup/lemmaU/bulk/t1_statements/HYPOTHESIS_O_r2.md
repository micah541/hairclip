# Hypothesis (O), revision r2 (supersedes `HYPOTHESIS_O.md`)

Bulk task 1 (Fable), Nov 13 2026. For `writeup/main/main_v2.tex`. Reasons: `R2_DECISIONS.md`.
Labels: **[cert]** arb, project script. **[est]** fit or float model on certified rows. **[probe]** arb cells on step-4 data, base curve only, not a pinned stage.

## What changed against r1
- **(O1)–(O3): no change.** (O1) stays at |f| ≤ 0.11: the certified segments of the repaired chain use 0.11. It can be relaxed to 0.15 after bulk task 2 certifies the segments at (0.15, 0.52); not before.
- **New: two hand-off conditions, (O-H).** The repaired Lemma S(A) needs its initial datum on its whole support [14, 700]. v3.1/v3.2 supply it on [14, 60] only.
- **The "facts" of r1 at the v3.1 thresholds are withdrawn.** They used the datum (H2) alone.

## The paragraph to quote

> **Hypothesis (O) (far part of the shell energies of Lemma S(A)).**
> Let γ satisfy the hand-off conditions (H1)–(H6) of step 3 v3.2 at t₁, and in addition the two conditions (O-H) below. Let τ′ ≥ τ₁ be such that the bootstrap hypotheses B1–B4 of step 3 hold on [τ₁, τ′]. Let f be the slope of the curve over the S₁ ray in the rescaled variable y, let v = f − U(α(τ), ·) with U(α, y) = 1.3056 α y^1.0246, and let z_k = (y∂_y)^k v. Let the shell energies of Lemma S(A) be cut off by a fixed smooth outer cutoff that equals 1 for y ≤ 350 and 0 for y ≥ 700. Then for all τ ∈ [τ₁, τ′]:
>
> - **(O1)** on 60 < y ≤ 350: |f| ≤ 0.11 and |y f_y| ≤ 0.52;
> - **(O2)** on 350 ≤ y ≤ 700: |f| ≤ 1/2 and |y f_y| ≤ 1;
> - **(O3)** on 350 ≤ y ≤ 700: |z_k| ≤ 1 for k = 0, 1, 2, 3.
>
> **(O-H) (hand-off conditions at t₁, in addition to (H1)–(H6)).**
> - **(H2, tightened)** the shell energies on [14, 60] satisfy E_j ≤ (10⁻⁸, 10⁻⁶, 10⁻³), j = 0, 1, 2;
> - **(H2f)** ∫ from y = 60 to 700 of (28/y)·(z_j)² dy/y ≤ (5·10⁻⁷, 5·10⁻⁶, 10⁻⁵), j = 0, 1, 2.

LaTeX form (balanced; same content):

```latex
Let $\gamma$ satisfy the hand-off conditions \textup{(H1)--(H6)} of step~3 \textup{v3.2} at $t_1$, and in addition the two conditions \textup{(O-H)} below. Let $\tau'\ge\tau_1$ be such that the bootstrap hypotheses \textup{B1--B4} of step~3 hold on $[\tau_1,\tau']$. Let $f$ be the slope of the curve over the $S_1$ ray in the rescaled variable $y$, let $v=f-U(\alpha(\tau),\cdot)$ with $U(\alpha,y)=1.3056\,\alpha\,y^{1.0246}$, and let $z_k=(y\partial_y)^kv$. Let the shell energies of Lemma~S\textup{(A)} be cut off by a fixed smooth outer cutoff that equals $1$ for $y\le350$ and $0$ for $y\ge700$. Then for all $\tau\in[\tau_1,\tau']$:
\begin{enumerate}[label=\textup{(O\arabic*)},leftmargin=2.6em,itemsep=1pt]
\item on $60<y\le350$: $|f|\le0.11$ and $|yf_y|\le0.52$;
\item on $350\le y\le700$: $|f|\le\frac12$ and $|yf_y|\le1$;
\item on $350\le y\le700$: $|z_k|\le1$ for $k=0,1,2,3$.
\end{enumerate}
\textup{(O-H)} Hand-off conditions at $t_1$, in addition to \textup{(H1)--(H6)}: \textup{(H2, tightened)} the shell energies on $[14,60]$ satisfy $E_j\le(10^{-8},10^{-6},10^{-3})$, $j=0,1,2$; \textup{(H2f)} $\int_{60}^{700}\frac{28}{y}\,z_j^2\,\frac{dy}{y}\le(5\cdot10^{-7},\,5\cdot10^{-6},\,10^{-5})$, $j=0,1,2$.
```

## Status of each part (state this next to the quotation)

| part | kind | status |
|---|---|---|
| (O1)–(O3) | bootstrap hypotheses, far field | assumed; to be derived in v3.3 |
| (H2, tightened) | hand-off, step 4 | **certified**: step 4 meets it by ×25.7 or more (`gonogo/t7_step4_H2_tight`) |
| (H2f) | hand-off, step 4 | **not certified.** Feasibility probe on step 4's base curve: (3.13e-7, 1.00e-6, 2.43e-6) [probe], so about ×1.5, ×4.9, ×4.0 before the error term. Certification is a new step-4 item (bulk task 6; needs Micah's approval). |

So the conditionality is "(F1), (O), and the hand-off condition (H2f)", until (H2f) is certified. After that it is "(F1) and (O)", with the hand-off list extended.

## What it is for (quotable)

Lemma S(A) of v3.1 was written with its coefficients beyond the core "bounded by |V| ≤ ½", with no estimate of its outer cutoff, and with an initial bound on [14, 60] for an energy that integrates far beyond. The repaired Lemma S(A) (corrected constants; erratum in preparation) is being certified **given (O) and (O-H)**: (O1) is the coefficient box its rate conditions use on (60, 350], (O2)–(O3) bound the fourth-order terms of the outer cutoff on [350, 700], and (O-H) is its initial datum.

## Facts the main text may state

- **Size of the outer cutoff.** Given (O), it adds at most (2.5e-8, 2.9e-8, 1.0e-6) to the shell energies [cert]. Unchanged from r1.
- **At the v3.1/v3.2 hand-off thresholds the repaired Lemma S(A) does not close**, even given (O) [cert, `bulk/t2_A/probe/asm_0_h3_v32.json`, setting v31]. With the initial datum taken honestly from (H2) and the (H3) windows, and nothing counted beyond y = 379.5, the flux constant is K0 = 2.80e-21 at the derived energies. At the bootstrap level 2Ê that is 5.6e-21, above the verifier's input 4.1e-21 (×0.73). The energy at level 0 is the cause: 1.39e-6 against v3.1's quoted 3.12e-7.
- **With (O-H) it is expected to close** [est; certification is bulk task 2]: Ê ≈ (7.0e-7, 1.0e-4, 0.073); K0 ≈ 1.40e-21, which is ×2.9 below the verifier's input at Ê and ×1.47 at 2Ê; core sup|v| ≈ 6.2e-3 against 0.0153 (×2.5); core M1 about ×4.
- **Do not quote** r1's "E = (3.86e-7, 1.92e-4, 0.269), S₀ᵛ ×3.25, M₀ ×5.74, M₁ ×2.73, flux ×4.9 and ×5.3": those are for the datum (H2) alone.
- **Where it sits.** All of (O) concerns y > 60, the far field, which is the region of (F).

## Relation to (F): as r1, with one addition

1. (O3) with k = 3 follows from nothing in v3.1/v3.2 (a fourth derivative of the curve; Lemma O has the third in L² only).
2. (O1) is stronger than what B3 / (F1) state (per-window boxes up to 0.236).
3. (O2) and (O3) for k ≤ 2 would follow from (F1) where r is small and from Lemma O beyond, after a chart conversion that has not been scripted.
4. The (F) macro's clause (F2) postulates the conclusion of Lemma S(A). The accurate structure is: (F1), (O) and (O-H) are the hypotheses; (F2) is then a theorem (the repaired Lemma S(A)).
5. **New.** That theorem is not yet certified end to end. The r1 sentence "(F2) is implied by (F1) + (O) through the repaired Lemma S(A)" should read: "(F2) is expected to follow from (F1), (O) and the hand-off conditions (O-H) through the repaired Lemma S(A); its certification with the full initial datum is in progress." Wording decision for Micah.

## Caveats on the numbers (not part of (O))

- The repaired Lemma S(A) is not yet in a versioned write-up.
- Verifier: the J_m constants are used as J_m(1 + κ_J), κ_J = 0.042348 [cert]. This covers window uniformity and a term omitted by the verifier script (the flux floor on its march cells, relative size 4.8e-7). The rows above include it.
- Point terms of the C² cutoffs at y = 28 (about 2% of E₀, E₁ [est]) and the y = 28 interface convention: bulk task 2, as in r1.
- v3.3 works with the bootstrap level (5/4)·Ê for the shell energies instead of 2Ê (`STATEMENTS_r2.tex` §6). For the closed versions, whose B2 has the factor 2, the flux margin above is the ×1.47 one.

## After v3.3

(O) disappears together with (F): the box on (60, 700] is derived (Lemma U's box where r < 253λ₁, Lemma O beyond), and the outer cutoff is bounded through L² energies and the outer dissipation. (H2f) and the tightened (H2) remain, as hand-off conditions certified by step 4.
