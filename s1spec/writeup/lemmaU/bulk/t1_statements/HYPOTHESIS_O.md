# Hypothesis (O), as it stands for closed step 3 (v3.1 / v3.2) until v3.3 exists

Bulk task 1 (Fable), Nov 13 2026. For quotation in `writeup/main/main_v2.tex`.
Labels: **[cert]** arb, project script (`lemmaU/phase2a/A_repair/`, `./run_chain.sh`). **[rev]** reviewer float run (`KILLCHECK2.md`). **[paper]** derived by hand in `STATEMENTS.tex`, not scripted.

## The paragraph to quote

> **Hypothesis (O) (far part of the shell energies of Lemma S(A)).**
> Let γ satisfy the hand-off conditions (H1)–(H6) of step 3 v3.2 at t₁, and let τ′ ≥ τ₁ be such that the bootstrap hypotheses B1–B4 of step 3 hold on [τ₁, τ′]. Let f be the slope of the curve over the S₁ ray in the rescaled variable y, let v = f − U(α(τ), ·) with U(α, y) = 1.3056 α y^1.0246, and let z_k = (y∂_y)^k v. Let the shell energies of Lemma S(A) be cut off by a fixed smooth outer cutoff that equals 1 for y ≤ 350 and 0 for y ≥ 700. Then for all τ ∈ [τ₁, τ′]:
>
> - **(O1)** on 60 < y ≤ 350: |f| ≤ 0.11 and |y f_y| ≤ 0.52;
> - **(O2)** on 350 ≤ y ≤ 700: |f| ≤ 1/2 and |y f_y| ≤ 1;
> - **(O3)** on 350 ≤ y ≤ 700: |z_k| ≤ 1 for k = 0, 1, 2, 3.

LaTeX form (balanced; same content):

```latex
Let $\gamma$ satisfy the hand-off conditions \textup{(H1)--(H6)} of step~3 \textup{v3.2} at $t_1$, and let $\tau'\ge\tau_1$ be such that the bootstrap hypotheses \textup{B1--B4} of step~3 hold on $[\tau_1,\tau']$. Let $f$ be the slope of the curve over the $S_1$ ray in the rescaled variable $y$, let $v=f-U(\alpha(\tau),\cdot)$ with $U(\alpha,y)=1.3056\,\alpha\,y^{1.0246}$, and let $z_k=(y\partial_y)^kv$. Let the shell energies of Lemma~S\textup{(A)} be cut off by a fixed smooth outer cutoff that equals $1$ for $y\le350$ and $0$ for $y\ge700$. Then for all $\tau\in[\tau_1,\tau']$:
\begin{enumerate}[label=\textup{(O\arabic*)},leftmargin=2.6em,itemsep=1pt]
\item on $60<y\le350$: $|f|\le0.11$ and $|yf_y|\le0.52$;
\item on $350\le y\le700$: $|f|\le\frac12$ and $|yf_y|\le1$;
\item on $350\le y\le700$: $|z_k|\le1$ for $k=0,1,2,3$.
\end{enumerate}
```

## What it is for (two sentences, quotable)

Lemma S(A) of v3.1 was written with its coefficients beyond the core "bounded by |V| ≤ ½" and with no estimate of its outer cutoff: the certified energy integrates to y = 350 and stops there. The repaired Lemma S(A) (corrected constants, erratum in preparation) is certified **given (O)**: (O1) is the coefficient box its rate conditions use on (60, 350], and (O2)–(O3) bound the fourth-order terms of the outer cutoff on [350, 700].

## Facts the main text may state

- **Size.** Given (O), the outer cutoff adds at most (2.5e-8, 2.9e-8, 1.0e-6) to the shell energies (E₀, E₁, E₂) [cert]. That is 13% of E₀ and negligible at levels 1 and 2.
- **Room.** With all three thresholds of (O3) multiplied by 10, the repaired E₂ rises from 0.0696 to 0.077 at the tightened hand-off, and the far-field M0 margin stays ×3.3 [rev]. (O2) is certified as a valid coefficient box for the outer segment (`weights_rep.log`, segment [350, 700]) [cert].
- **Repaired Lemma S(A) at the v3.1 hand-off thresholds, given (O)** [cert]: E = (3.86e-7, 1.92e-4, 0.269), below the bootstrap level 2Ê = (6.24e-7, 3.62e-4, 0.496); core margins S₀ᵛ ×3.25, M₀ ×5.74, M₁ ×2.73; flux constants ×4.9 and ×5.3 below the verifier inputs.
- **Where it sits.** All of (O) concerns y > 60, the far field, which is the region of (F).

## Relation to (F): read this before quoting

1. **(O3) with k = 3 follows from nothing in v3.1/v3.2.** z₃ is a fourth derivative of the curve. Lemma O bounds the third derivative of the graph only in L² (`step3_proof_v32.tex` l.265), and the window energies stop at level 3 on the strip and level 2 elsewhere, which gives z₃ in L² at best. So (O) is a genuine extra hypothesis, not a corollary of (F).
2. **(O1) is stronger than what B3 / (F1) state.** The per-window boxes of v3.2 §3.7 allow |f| ≤ U(1.5 y_w) + 0.15, up to 0.236 on the last window. The repaired chain is certified with |f| ≤ 0.11 on all of [28, 350]. The *derived* window bound of D4 (|f| ≤ 0.093) would give (O1), but D4 comes after D3 in the derivation order and is itself conditional on (F), so this is not available as written.
3. **(O2) and (O3) for k ≤ 2** would follow from (F1) where r < 190λ₁ and from Lemma O (D5a) beyond, after a chart conversion that has not been scripted.
4. **Overlap with the (F) macro of `main.tex`.** Its clause (F2) postulates the *conclusion* of Lemma S(A) (E^(A) < 2Ê and the core bounds), "cut off at an outer ramp in y > 60 whose fourth-order terms are kept". So if (F) is kept verbatim, (O) is logically redundant. The accurate structure is:
   - (F1) (windows) and (O) are the hypotheses;
   - (F2) is then a theorem: the repaired Lemma S(A), certified.
   Recommendation for main_v2: keep the (F) macro text for now if Micah prefers, but say in §"The gap" that (F2) is implied by (F1) + (O) through the repaired Lemma S(A), so that the conditionality is really "(F1) and (O)". This is a wording decision for Micah.

## Caveats on the certified numbers (not part of (O))

- **Lemma S(A) repaired is not yet in a versioned write-up** (constants E, γ′, E_B, K0, K2, I₃ differ from v3.1; `A_REPAIR.md`).
- **One more omitted term, found in this task [paper], not yet in the certified chain.** The cutoffs are C² only: the third derivative of η² jumps at the end where η = 1 (y = 28 for the inner ramp). Integrating the m = 3 monomial of the pair (j, j+1) by parts therefore leaves a point term +g_c·z_j(28, τ)², g_c = 360.3·28⁻⁴·a ≤ 5.9e-4 per unit τ. At level 2 the repaired chain also integrates the top pair by parts, which contributes −2g_c·z₂(28)², so the net term is favourable and may be dropped. At levels 0 and 1 the top pair goes by route J, so the positive term is missing. Estimated size: about 2% of E₀ and of E₁ [est., crude]. Task 2 of the bulk certifies or removes it. v3.1 itself omits the term at all three levels. Details: `STATEMENTS.tex`, Lemma P.
- The ramp/core interface conventions at y = 28 are inherited from v3.1 and not yet audited (bulk task 2).

## After v3.3

(O) disappears together with (F). In v3.3 the box on (60, 700] is *derived* (Lemma U's box where r < 190λ₁, Lemma O beyond), and the outer cutoff is bounded through L² energies of z₀, z₁, z₂ and the outer dissipation, with no pointwise bound on z₂ or z₃ (`STATEMENTS.tex`, Proposition R5′).
