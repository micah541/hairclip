# (P1') literature check: DKS 2002 and Chou 2003

**Verdict: CONFIRMED WITH CHANGES.**
The L^2 blow-up criterion that `continuation_v33.tex` needs is stated for closed immersed curves in R^n (so for plane curves too) as **Dziuk–Kuwert–Schätzle 2002, Theorem 3.1** (SIAM J. Math. Anal. 33(5), p. 1234). It is also stated for closed immersed plane curves as **Chou 2003, Proposition A** (Hokkaido Math. J. 32, p. 2; proof pp. 9–10).

The changes needed are to wording and citations, not to mathematics:
- DKS do not prove **uniqueness**, and they only *sketch* short-time existence. The "unique" in (P1') needs a separate citation.
- Chou's result is called "Proposition A", not a theorem. Its rate holds only "when t is close to ω", and its constant depends on the area and the turning number.
- The DKS form is the clean one: it holds for all t < T.

Sources read: the text layers of both PDFs, extracted with PyMuPDF (read only, nothing executed). Page numbers below are journal page numbers.

---

## 1. Which results give LWP, smoothness, T_max and the L^2 criterion; their exact form

### DKS 2002, Theorem 3.1 (journal p. 1234; PDF p. 7)

Verbatim:
> "Theorem 3.1. Let f : [0, T) × I → R^n be a maximal solution of the curve diffusion equation ∂_t f = −∇_s^2 κ. If T < ∞, then ∫_I |κ|^2 ds ≥ c (T − t)^{−1/4}."

Here I = R/Z (p. 1228, "a closed curve f : I = R/Z → R^n"), κ = ∂_s^2 f is the curvature vector, and ∇_s is the normal part of ∂_s.

- **Same flow as CDF.** For n = 2, write κ = kN. Then ∇_s^2 κ = k_ss N, so the equation is ∂_t f = −k_ss N, i.e. normal velocity −k_ss. DKS p. 1229 says so: "for n = 2, this is the one-dimensional analogue of surface diffusion".
- **Form of the criterion.** The statement gives a lower bound **for every t ∈ [0, T)**. It therefore implies lim_{t↑T} ∫κ^2 ds = ∞, a limit and not just a limsup.
  - The proof first shows lim sup_{t↗T} ‖κ‖_{L^2} = ∞ by contradiction (p. 1234: "Now we prove by contradiction that lim sup_{t↗T} ‖κ‖_{L^2} = ∞").
  - It then integrates d/dt ∫|κ|^2 ≤ c(∫|κ|^2)^5 (p. 1235) "on [t, t_l], where t_l ↗ T with ‖κ‖_{L^2} → ∞".
- **Rate.** ∫κ^2 ds ≥ c (T−t)^{−1/4}. Integrating the ODE inequality gives c = (4c_0)^{−1/4}, where c_0 is the constant in (3.2) with m = 0.
- **Dependence of the constant (my reading of the proof; the statement says only "c").**
  - c_0 comes from (2.16) of Prop. 2.5 with k = 2, μ = 2, ν = 4. Both exponents there equal 5, and the constant is c(n, k, μ, ν).
  - So **c depends only on n**: it is universal for plane curves. It does not depend on length, area, turning number or the initial data.
- **Short-time existence (p. 1233).** It is only sketched: "short-time existence is a standard matter, and we only briefly sketch the argument for the case of (1.6)". The sketch is DeTurck-type (tangential term plus ODE reparametrisation) with "standard theory" in parabolic Hölder spaces.
  - **Uniqueness is neither stated nor proved.** "Maximal solution" presupposes it.
- **Smoothness up to and beyond T if ∫κ^2 stays bounded.** This is inside the proof (p. 1234–1235): "f extends smoothly to [0, T]×I, and in fact even beyond T by short-time existence, which contradicts the maximality of T".
- **Length decrease.** DKS (2.10), p. 1230: d/dt L(f) = −∫|∇_s κ|^2 ds ≤ 0. Cite this for the line "CDF decreases length" in continuation_v33.

### Chou 2003, Proposition A (p. 2; proof pp. 9–10)

Verbatim (p. 2):
> "Proposition A. Let γ(·, t) be a solution of (1) where γ_0 is a closed, immersed, smooth curve with nonzero area. Then the solution exists as long as the L^2-norm of the curvature of γ(·, t) is finite. Furthermore, when ω is finite, ∫_{γ(·,t)} k^2(s, t) ds ≥ C(ω − t)^{−1/4} for some constant C."

- **The flow.** Chou's equation (1) is ∂γ/∂t = −k_ss n, i.e. CDF.
- **Maximal solution.** p. 2: "the unique solvability of (1) for small time enables us to define a unique maximal solution in a maximal interval [0, ω), ω ≤ ∞."
- **Local existence and uniqueness are taken from Giga–Ito [12]** (Comm. Appl. Anal. 2 (1998) 383–405). p. 9: "From [12] we know that for any closed, immersed H^4-curve, there exists a unique solution of (1) in some [0, t_1), t_1 > 0 where t_1 depends on the H^4-norm of the initial curve." I have **not** seen Giga–Ito; this is Chou's description of it.
- **Rate.** It is proved only asymptotically. p. 10: integrating (8) gives "E^4(t) ≥ (1/(8C))(ω − t)^{−1}, when t is close to ω".
  - Here E = ∫(k − k̄)^2 ds (p. 7). Since ∫k^2 = E + L k̄^2 ≥ E, this gives the stated bound.
  - C_3 in (8) "depends only on the initial area, total curvature, and the best constants in these interpolation inequalities" (p. 8).
  - Area enters through Radó's isoperimetric inequality, used to bound |k̄| (p. 7).
- **Form of the criterion.** It is again effectively a limit, because the lower bound holds for t near ω. The existence statement "exists as long as the L^2-norm … is finite" is the limsup form.
- **Weaker data.** Chou remarks (pp. 11–12) that the results extend to initial curves with only L^2 curvature. This is not needed here.

**Bottom line for Q1.** Cite DKS Thm 3.1 for the criterion and the rate. It is cleaner: the bound holds for all t, the constant is universal, and there are no area or turning-number hypotheses. Cite Chou Prop. A as the plane-curve statement, with uniqueness via Giga–Ito. Cite an LWP-with-uniqueness reference separately (§4).

## 2. Hypotheses

- **Immersed, not embedded.** Both results cover immersed curves. Neither assumes embeddedness.
  - DKS: f : I = R/Z → R^n is a parametrised closed curve. Regularity of the parametrisation, γ = |∂_x f| > 0, is used throughout, for example in (3.5). Embeddedness is never mentioned.
  - Chou: "γ_0 is a closed, immersed, smooth curve" (Prop. A). p. 2 notes that embeddedness is not preserved under (1).
- **Dimension.** DKS: any n ≥ 2, including n = 2. Chou: plane curves only.
- **Initial regularity.**
  - DKS: smooth. This is implicit: Lemma 2.4 says "smooth closed curve", and the LWP sketch is in Hölder spaces.
  - Chou: smooth in Prop. A; LWP for H^4 via Giga–Ito.
  - Our γ(t_1) is smooth, since the solution is smooth for t > 0 and Γ0 is an exact trigonometric polynomial.
- **Length.** No hypothesis in either. L ≤ L(f_0) is used inside DKS's proof, via (2.10).
- **Turning number and area.**
  - DKS: none.
  - Chou: Prop. A assumes **nonzero (algebraic) area**, and its constant depends on the area and the total curvature 2nπ.
  - For the hair clip, Chou's area hypothesis holds. Γ0 is embedded and closed, so its algebraic area is ± its enclosed area ≠ 0. CDF preserves algebraic area (Chou (4)–(5); DKS p. 1231).
  - DKS needs no such hypothesis, which is one more reason to make DKS the primary citation.

## 3. Is DKS Theorem 3.1 about CDF? Yes.

DKS 2002 is mostly about the elastic flows (Theorems 3.2 and 3.3), but **Theorem 3.1 is specifically about the curve diffusion flow (1.6), ∂_t f = −∇_s^2 κ**, not the elastic flow.
- The abstract says "a lower bound for the lifespan of solutions to the curve diffusion flow is observed".
- p. 1229: "For the curve diffusion flow (1.6), we prove a lower bound on the existence time in terms of the L^2-norm of κ".
- p. 1233: "We now start with an estimate for the curve diffusion flow (1.6)."

Citing DKS Thm 3.1 for CDF is **correct**. Citing Thm 3.2 or 3.3, or "DKS" generically for the elastic flow, would be wrong.

## 4. Is it enough for continuation_v33.tex?

**How (P1') is used** (Proposition C:prop:global):
- Assume τ_max < ∞. By B1 and (P2'), λ ≥ λ_* > 0, so T_max = t_1 + ∫λ^4 dτ < ∞.
- Show sup_{[t_1, T_max)} ∫_γ κ^2 ds < ∞:
  - tip: L^2_ρ control of w̃'' via B1, divided by λ_*;
  - elsewhere: pointwise κ ≤ C/(14λ_*) or 0.3/(50λ_1), times L(t) ≤ L(t_1).
- Contradict (P1'). The only input is: **T_max < ∞ ⇒ ∫κ^2 is unbounded on [t_1, T_max)**, for the maximal smooth solution started from the smooth immersed curve γ(t_1).

**Mathematically, DKS Thm 3.1 suffices.**
- It gives a lower bound for every t < T, which is stronger than unboundedness.
- It is stated for an arbitrary initial time. Apply it with time shifted by t_1 and f_0 = γ(t_1).
- It holds for immersed curves.
- No rate is needed for the contradiction. Neither the rate nor the constant enters the argument.

**Gaps (all citation or wording, none quantitative):**
1. **Uniqueness and LWP.** DKS gives no uniqueness proof and only sketches existence. "Unique smooth solution on a maximal interval" needs a citation. Candidates:
   - Giga–Ito 1998 (as Chou cites it, for closed immersed H^4 curves);
   - Escher–Mayer–Simonett 1998 (immersed hypersurfaces; cited by both DKS and Chou).

   I have **not** read either. Whether EMS 1998 covers curves in R^2 (hypersurfaces of dimension 1) must be checked by the authors, and likewise Giga–Ito's exact statement.

   Uniqueness matters because step 4's validated solution, the modulated solution of step 3, and "the maximal solution" of (P1') must be the same object.
2. **Identifying τ_max with T_max.** The proof sets T_max = t_1 + ∫_{τ_1}^{τ_max} λ^4. That presupposes τ_max = "end of the rescaled interval of existence" corresponds to the end of the physical maximal interval: the modulation cannot break down before T_max while B1 holds.

   That is a (P2') matter (implicit function theorem for the modulation parameters), not (P1'). The text is fine if (P2') says so. One sentence would close it: "by (P2') the modulation extends as long as the physical solution exists and B1 holds, so t(τ_max) = T_max".
3. **No extra estimate is needed.** No scaling or energy argument has to be added. The text should simply use DKS's form (a lower bound for all t, hence → ∞) or the weaker "unbounded" form, both of which are proved.

## 5. Inaccuracies in continuation_v33.tex

Lines refer to `continuation_v33.tex`.

- **l.51, "Theorem 3.1 (closed immersed curves in R^n)".** Correct theorem number and setting.
- **l.51, "and/or of Chou: as we recall it, a lower bound ∫κ^2 ds ≥ c (T_max − t)^{−1/4}".** Misattributed.
  - This exact lower bound is the *statement of DKS Theorem 3.1*.
  - Chou has the same form as **Proposition A**, which is not a theorem number. Its rate is proved only "when t is close to ω", with a constant depending on the area and the turning number, under a nonzero-area hypothesis.
- **l.51, "has a unique smooth solution".** Not supplied by either cited paper.
  - DKS sketches existence and is silent on uniqueness.
  - Chou cites Giga–Ito [12] for existence and uniqueness.
  - A uniqueness citation must be added (§4, item 1).
- **l.51, "∫κ^2 → ∞ as t ↑ T_max".** Correct as a consequence of either result: a lower bound for all t (DKS) or near ω (Chou).
- **l.16 and l.107, "[to be confirmed by the authors]".** These can now say: confirmed, DKS Thm 3.1 and Chou Prop. A. The LWP and uniqueness reference is still to be confirmed.
- **l.75, "which is what (P1') needs to identify T as the maximal time".** Harmless but slightly misattributed. That T is maximal (the solution cannot be continued smoothly past T) follows from ∫_tip κ^2 → ∞ alone, because a smooth solution on [t_1, T] has bounded ∫κ^2. (P1') is not needed there. Optional rewording.
- **"CDF decreases length" (l.55).** Can cite DKS (2.10) or Chou (5).
- **Wheeler remark ("Wheeler's lifespan theorem is not used").** Neither PDF bears on Wheeler, so I cannot check the claim about the earlier draft. As a statement of what is used, it is accurate: the argument uses only DKS Thm 3.1 or Chou Prop. A. DKS Thm 3.1 is itself a lifespan lower bound (T − t ≥ (c/∫κ^2)^4), so if anything the "lifespan" role is filled by DKS.
- **Bibliography.**
  - DKS: title, SIAM J. Math. Anal. 33(5), 1228–1245, 2002. Matches the PDF (p. 1 header and metadata).
  - Chou: title, Hokkaido Math. J. 32 (2003) 1–19. Matches the PDF (p. 1).
  - **The issue number (1) is not shown in the PDF.** It is plausible (pages 1–19) but unverified here.
  - main.bib's `GigaIto1998` omits pages. Chou's reference list gives "Comm. Appl. Anal. 2 (1998), 383-405" (unverified beyond Chou).

---

## Suggested replacement text (LaTeX) for (P1')

```latex
\item \emph{Local well-posedness and blow-up criterion.} For a smooth closed
immersed (not necessarily embedded) plane curve, CDF has a unique smooth solution
on a maximal interval $[t_1,T_{\max})$ (local existence and uniqueness:
\cite{EscherMayerSimonett1998,GigaIto1998}; see also the sketch in
\cite[\S3]{DKS2002}) \textbf{[authors: confirm the LWP/uniqueness statement in
EMS 1998 or Giga--Ito 1998 for closed immersed plane curves]}. If
$T_{\max}<\infty$, then
\[
  \int_{\gamma(t)}\kappa^2\,ds\;\ge\;c\,(T_{\max}-t)^{-1/4}\qquad(t_1\le t<T_{\max}),
\]
in particular $\int_{\gamma(t)}\kappa^2ds\to\infty$ as $t\uparrow T_{\max}$, with
$c>0$ depending only on the dimension: Dziuk--Kuwert--Sch\"atzle
\cite[Theorem~3.1]{DKS2002}, stated for maximal solutions $f:[0,T)\times\mathbb
R/\mathbb Z\to\mathbb R^n$ of $\partial_tf=-\nabla_s^2\kappa$ (for $n=2$ this is
CDF), with no embeddedness assumption. For closed immersed plane curves with
nonzero algebraic area (which holds here, since $\Gamma_0$ is embedded and CDF
preserves algebraic area) the same criterion is Chou's
\cite[Proposition~A]{Chou2003}, with the rate proved for $t$ near $T_{\max}$. Only
the qualitative consequence ``$T_{\max}<\infty\Rightarrow\int\kappa^2ds$ is
unbounded on $[t_1,T_{\max})$'' is used below.
```

Also: in the gap list (l.107), replace "[to be confirmed by the authors]" with "confirmed (DKS Thm 3.1; Chou Prop. A); LWP/uniqueness reference to be confirmed". At l.55 add "\cite[(2.10)]{DKS2002}" after the length identity.

## Bibliography entries

```latex
\bibitem{DKS2002} G. Dziuk, E. Kuwert, R. Sch\"atzle, Evolution of elastic curves
in $\mathbb R^n$: existence and computation, SIAM J. Math. Anal. 33(5) (2002)
1228--1245.  % verified against PDF
\bibitem{Chou2003} K.-S. Chou, A blow-up criterion for the curve shortening flow by
surface diffusion, Hokkaido Math.\ J.\ 32 (2003) 1--19.  % verified against PDF;
% issue no. (1) not shown in PDF
\bibitem{GigaIto1998} Y. Giga, K. Ito, On pinching of curves moved by surface
diffusion, Comm. Appl. Anal. 2 (1998) 383--405.  % as cited by Chou; not verified
\bibitem{EscherMayerSimonett1998} J. Escher, U. F. Mayer, G. Simonett, The surface
diffusion flow for immersed hypersurfaces, SIAM J. Math. Anal. 29(6) (1998)
1419--1433.  % as cited by DKS [6] and Chou [11]; not read
```

### Uncertainty, stated plainly
- I read only DKS 2002 and Chou 2003.
- I have not checked the statements of Giga–Ito 1998 or EMS 1998, so the "unique" part of (P1') is not confirmed.
- That DKS's c depends only on n is my reading of the proof, from (2.16) and (3.2); the theorem statement itself says only "c".
- The constant plays no role in the continuation argument.
