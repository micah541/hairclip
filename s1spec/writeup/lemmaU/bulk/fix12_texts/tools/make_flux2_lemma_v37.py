# make_flux2_lemma_v37.py -- builds fix12_texts/flux2_lemma_v37.tex (flux lemma v2, final text) from
# fix8_flux2_impl/texts/flux2_lemma_v2.tex (unchanged) by asserted-unique replacements, applying the text-fix list T1-T4 of
# fix9_review/FLUX2_REVIEW.md Section 7: status, "to be asserted", the march first-touch argument with Y = Lambda - Z as a proof,
# the bootstrap set with strict improvements, (I3') as the only flux hypothesis, B1 -> Lambda <= L + 2C', tail-floor factor.
import sys
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
src = open(B + 'fix8_flux2_impl/texts/flux2_lemma_v2.tex').read()
R = []
def rep(o, n): R.append((o, n))
rep(src[:src.index('\\documentclass')],
    "% flux2_lemma_v37.tex -- flux lemma v2 for the outer-ramp functional K^OR (levels 0 and 2), the data lemma, and the time-integrated\n"
    "% absorption, FINAL TEXT (text v37, fix12_texts, 20 Nov 2026).  Built by fix12_texts/tools/make_flux2_lemma_v37.py from\n"
    "% fix8_flux2_impl/texts/flux2_lemma_v2.tex (unchanged), with the text fixes T1-T4 of fix9_review/FLUX2_REVIEW.md.\n"
    "% Numbers [cert] from fix8_flux2_impl/out/flux2_cert.json (and the stage JSONs there).  Not compiled (no LaTeX on this machine).\n"
    "% Standalone, or \\input by step3_v37_main.tex (which defines \\stepmain); in that case the preamble below is skipped.\n"
    "\\ifdefined\\stepmain\\else\n")
rep("\\begin{document}\n\\title{Flux lemma v2 for the outer ramp of $E_I$ (step 3): statements and proof sketches}\n\\date{20 November 2026 (fix round 8, implementation). \\emph{Not compiled.} Numbers [cert] from \\fn{fix8_flux2_impl/out/}.}\n\\maketitle\n",
    "\\begin{document}\n\\title{Flux lemma v2 for the outer ramp of $E_I$ (step 3): statements and proofs}\n\\date{20 November 2026 (text v37: final text of flux lemma v2). \\emph{Not compiled.} Numbers [cert] from \\fn{fix8_flux2_impl/out/}.}\n\\maketitle\n\\fi\n"
    "\\providecommand{\\cst}[1]{\\ensuremath{\\langle\\texttt{\\detokenize{#1}}\\rangle}}\n"
    "\\section{Flux lemma v2 (outer ramp of $E_I$)}\\label{F:sec:top}\n"
    "\\noindent\\textbf{Status (text v37).} \\emph{Implemented and certified; independently reviewed; not closed.} The certified computation is \\fn{fix8_flux2_impl} (stages X0--X7, two byte-identical production passes; \\fn{FLUX2_IMPL_REPORT.md}); every scripted condition (C0)--(C7) passes at both production points of \\fn{fin13v4b}. "
    "The independent review \\fn{fix9_review/FLUX2_REVIEW.md} (Fable, Nov 20) found no error in the mathematics or in the certified numbers (verdict ``sound modulo listed items'') and recomputed $W(57)$, $\\bar N$, $K_2^{\\rm pt}$, $C'$, $S'$ and (C2)--(C6) independently. "
    "Remaining conditions: the regularity lemma U-1 (C8; written in \\fn{fix10_proofs/U1_regularity.tex}, \\tHd, not yet reviewed) and the written level-2 export (G1; \\fn{fix10_proofs/G1_export_k2.tex}, not yet reviewed). The proofs of Lemmas~\\ref{F:lem:data}, \\ref{F:lem:chart4}, \\ref{F:lem:flux} are hand proofs whose arithmetic is certified; the proofs of Lemma~\\ref{F:lem:Z} and Proposition~\\ref{F:prop:abs} are complete hand proofs. "
    "\\textbf{The only flux hypothesis of the verifier is the combined (I3$'$)} of Proposition~\\ref{F:prop:abs}; a separate inequality $K_2^{\\OR}\\le\\fn{FLUX2}$ is not claimed and is false ($K_2^{\\rm pt}=1.74\\times\\fn{FLUX2}$).\n\n")
rep("\\newtheorem{lemma}{Lemma}\\newtheorem{proposition}{Proposition}\\newtheorem{remark}{Remark}\n",
    "\\newtheorem{lemma}{Lemma}\\newtheorem{proposition}{Proposition}\\newtheorem{remark}{Remark}\\newtheorem{definition}{Definition}\n")
# sections -> subsections (so the file nests under \section in the main document)
for t in ["\\section{Setting}\\label{F:sec:set}", "\\section{Data}\\label{F:sec:data}", "\\section{The flux lemma}\\label{F:sec:flux}",
          "\\section{Absorption of the time-integrated part}\\label{F:sec:abs}", "\\section{Implementation notes (fix round 8)}\\label{F:sec:impl}"]:
    rep(t, t.replace('\\section', '\\subsection'))
# T1: data lemma proof with the implemented Minkowski form
rep("so\n\\[M_3^2\\le\\frac3{\\sigma_{\\min}^5}\\int_S\\bigl[(D^3g)^2+9(D^2g)^2+4(Dg)^2\\bigr]ds,\\qquad M_4^2\\le\\frac4{\\sigma_{\\min}^7}\\int_S\\bigl[(D^4g)^2+36(D^3g)^2+121(D^2g)^2+36(Dg)^2\\bigr]ds .\\]",
    "so, by Minkowski in $L^2(S\\times I;ds\\,d\\tau)$ (the form implemented, \\S\\ref{F:sec:impl}(i); $ds=d\\sigma/\\sigma$),\n"
    "\\[\\|M_3\\|_{L^2(I)}\\le\\sigma_{\\min}^{-5/2}\\bigl(\\|D^3g\\|+3\\|D^2g\\|+2\\|Dg\\|\\bigr),\\qquad\\|M_4\\|_{L^2(I)}\\le\\sigma_{\\min}^{-7/2}\\bigl(\\|D^4g\\|+6\\|D^3g\\|+11\\|D^2g\\|+6\\|Dg\\|\\bigr),\\]\n"
    "and $M_2\\le\\sigma_{\\min}^{-3/2}(\\|D^2g\\|+\\|Dg\\|)$; sums of window norms $\\sum_iw_i(a_i+b_iT)^{1/2}$ are squared with fixed rational weights $\\theta_i$, $\\sum\\theta_i\\le1$ (Cauchy--Schwarz), which keeps every bound affine in $T$. (The squared-sum form $M_3^2\\le3\\sigma_{\\min}^{-5}\\int_S[(D^3g)^2+9(D^2g)^2+4(Dg)^2]$ of the design is also valid and is weaker.)")
rep("Products are split by AM--GM so that every bound is affine in $T$.", "The weights are those of (B$\\Phi$) at the bootstrap level; $\\int\\varepsilon\\omega(z^A_{k+2})^2$ is the export of Lemma~S(A) (\\fn{lemmaS_A_v37.tex}, Lemma~D, $k=1,2$).")
# T1: chart lemma "To be asserted"
rep("\\textbf{To be asserted in sympy} (as \\fn{n1_chart.py} does for $k\\le3$).",
    "This structure is asserted in sympy in stage X3 (\\fn{chart4_shell.py}; the assertion passed) and was confirmed independently with a construction sharing no code (\\fn{fix9_review/chi4_check.py}: $\\chi_4$ has no $b_5,b_6$, no $k^{(5)}$; $\\partial^2\\chi_4/\\partial b_i\\partial b_j=0$ for $(2,3),(2,4),(3,3),(3,4),(4,4)$; $\\partial^3\\chi_4/\\partial b_2^3=0$; $\\chi_k=0$ on a ray) \\tS.")
rep("\\item The level-4 chart structure of Lemma~\\ref{F:lem:chart4} is asserted in sympy (\\fn{chart4_shell.py}): $b_2^2=(u^{(3)})^2$ is the only quadratic monomial in $(b_2,b_3,b_4)$, no $b_5$, $\\chi_4=0$ for a ray.",
    "\\item The level-4 chart structure of Lemma~\\ref{F:lem:chart4} is asserted in sympy (\\fn{chart4_shell.py}; passed) and independently confirmed (\\fn{fix9_review/chi4_check.py}): $b_2^2=(u^{(3)})^2$ is the only quadratic monomial in $(b_2,b_3,b_4)$, no $b_5$, $\\chi_4=0$ for a ray.")
rep("Only $4W\\zeta\\zeta'j\\,g\\,g'''$ qualifies; its coefficient is differentiated once, so the reduced coefficients are affine in $(j,j')$ and contain no $j''$ \\textbf{(asserted)}.",
    "Only $4W\\zeta\\zeta'j\\,g\\,g'''$ qualifies; its coefficient is differentiated once, so the reduced coefficients are affine in $(j,j')$ and contain no $j''$ (asserted in stage X1, \\fn{flux2_forms.py}, 35 checks PASS).")
# T3/T4: absorption proposition: full proof of the march first touch (FLUX2_REVIEW 1.1-1.4)
old = src[src.index("\\begin{proof}[Sketch]\\tH\\ Apply the verifier's first-touch arguments"):src.index("\\begin{remark}[bootstrap order]")]
new = ("\\begin{proof}\\tH\\ Fix one of the two production points and let $Y=\\Lambda-Z$ be as in Lemma~\\ref{F:lem:Z}; by (C0), (C1) the lemma applies, so $Y\\le\\Lambda\\le Y+2C'$ and $\\frac12Y'\\le\\mathrm{tot}(\\alpha(\\tau),x(\\tau))$ a.e. "
 "Throughout, ``the state is in cell $k$'' means $\\Lambda\\in[\\mathrm{Ly}^{lo}_k,\\mathrm{Ly}^{hi}_k]$ and $|\\alpha|\\le\\alpha^{\\rm up}_k$; the alpha-barrier itself is checked by the verifier at $\\mathrm{Ly}^{hi}_k$, and $\\beta$-bounds are evaluated there.\n\n"
 "\\emph{What a cell certifies} \\tR. \\fn{closure.py}:129 calls \\fn{orc.rate}, i.e.\\ \\fn{bb.bound_cell}$(0,\\alpha_k,\\mathrm{Ly}^{lo}_k,\\mathrm{Ly}^{hi}_k,\\text{floors}=\\text{True})$, which returns $U_k\\ge\\sup\\mathrm{tot}/\\Lambda$ over all shapes, $|\\alpha|\\le\\alpha_k$ and $\\Lambda\\in[\\mathrm{Ly}^{lo}_k,\\mathrm{Ly}^{hi}_k]$ (\\fn{bb.py}:142--185; the $\\Lambda$-range by evaluating $H-\\mathrm{NEG}$ at $s_{hi}$ and $L$ at $s_{lo}$; the flux constants enter $L$ as the degree-0 constant $F_0$, taken at $\\mathrm{Ly}^{lo}$; \\fn{TAIL_LEMMA.tex}, ``s-monotone split'' and ``leaf certificate''). "
 "Cell ends: $\\mathrm{Ly}^{hi}_k=\\max(L_k,g)(1+10^{-3})$, $\\mathrm{Ly}^{lo}_k=\\min(L_k,g)(1-10^{-3})$, and the step is accepted only if $L_{k+1}\\in[\\mathrm{Ly}^{lo}_k,\\mathrm{Ly}^{hi}_k]$ (\\fn{closure.py}:112--115, 131--132). The barrier on step $k$ is $L'=2(U_k+\\mathrm{pad}_k)L$, $\\mathrm{pad}_k=10^{-9}\\max(1,|U_k|)$ (\\fn{closure.py}:130).\n\n"
 "\\emph{March (first touch with $Y$).} Claim: $Y\\le L$ on every step. At $\\tau_1$, $Z=0$ and $Y=\\Lambda(\\tau_1)\\le\\mathrm{Ly}_1<L(0)$. Suppose the claim fails on step $k$ and put $G:=Y-L$, absolutely continuous. Let $t_0$ be the last time before a point with $G>0$ at which $G\\le0$; then $G(t_0)=0$ and $G>0$ on $(t_0,t_0+\\delta)$. "
 "For $t\\in(t_0,t_0+\\delta)$ and $\\delta$ small: $\\Lambda\\ge Y>L\\ge\\mathrm{Ly}^{lo}_k$, and $\\Lambda\\le Y+2C'=L+G+2C'<\\mathrm{Ly}^{hi}_k$ because $L+2C'\\le\\mathrm{Ly}^{hi}_k$ holds with room on every step by (C2) (checked against $\\max(L_k,L_{k+1})$, margin $\\ge\\times5.0\\cdot10^{7}$) and $G$ is small. The $\\alpha$-barrier is unchanged (it is evaluated at $\\mathrm{Ly}^{hi}_k$ and $\\Lambda\\le\\mathrm{Ly}^{hi}_k$). So the state is in cell $k$ and $\\frac12Y'\\le\\mathrm{tot}\\le U_k\\Lambda$ a.e. on $(t_0,t_0+\\delta)$.\n"
 "If $U_k\\ge0$: $\\Lambda\\le Y+2C'$ gives $Y'\\le2U_k(Y+2C')$, and $L'=2(U_k+\\mathrm{pad}_k)L$, so $G'\\le2U_kG+(4U_kC'-2\\,\\mathrm{pad}_kL)$; the bracket is $<0$ by (C3) (with $L\\ge\\mathrm{Ly}^{lo}_k$, margin $\\ge\\times1547$). Gronwall from $G(t_0)=0$ gives $G\\le0$ on $(t_0,t_0+\\delta)$, a contradiction. "
 "If $U_k<0$: $\\Lambda\\ge Y$ gives $Y'\\le2U_k\\Lambda\\le2U_kY$, so $G'\\le2U_kG-2\\,\\mathrm{pad}_kL<0$ wherever $G\\ge0$ is small; again $G\\le0$. Hence $Y\\le L$ and $\\Lambda\\le L+2C'$ on the march. No march step of either table used the slack path of \\fn{closure.py}:126--128 (\\fn{recompute.json}, $n_{\\rm slack}=0$).\n\n"
 "\\emph{Tail.} The anchor certificate gives $\\mathrm{rate}_s\\le-\\rho$ on $\\Lambda\\in[L^*,(1+\\eta)L^*]$, $|\\alpha|\\le a^*$, floors off, $\\eta=0.01$ (\\fn{closure.py}:165). Let $\\ell$ be the tail barrier, $\\ell'=-\\rho\\ell$ above the floor and constant at $L_F=2F_0/\\rho$. At a first touch $G=Y-\\ell=0$, $G>0$ just after: $\\Lambda\\in(\\ell,\\ell+2C'+G]$; with $\\lambda^2=\\ell/L^*$ the scaled state has $\\Lambda/\\lambda^2\\in[L^*,(1+\\eta)L^*]$ as soon as $2C'<\\eta L_F$ (C4, margin $\\ge\\times1.7\\cdot10^5$), and $|\\alpha|/\\lambda\\le a^*$ by the $\\alpha$-barrier. Joint scaling gives $\\mathrm{tot}_s\\le-\\rho\\Lambda$, so $Y'\\le-2\\rho\\Lambda+2F_0\\le-2\\rho Y+2F_0$. This is the differential inequality of the floor paragraph of \\fn{TAIL_LEMMA.tex} with $Y$ in place of $\\Lambda$: above the floor $G'\\le-2\\rho G$ (using $\\rho\\ell\\ge2F_0$ there), on the flat part $G'\\le-2\\rho G-\\rho L_F$; so $G\\le0$, with the same $L_F$. "
 "The tail anchor boxes are asserted at $L^*$ without cushion (\\fn{closure.py}:162); they cover $\\Lambda\\le L^*+2C'$ because $L^*+2C'\\le\\mathrm{Ly}^{hi}$ of the last march cell, whose boxes are asserted at $\\mathrm{Ly}^{hi}$, and $a^*\\le\\alpha^{\\rm up}$ of that cell (margins $\\times1.97\\cdot10^8$, $\\times5.1\\cdot10^7$; $a^*=3.68\\cdot10^{-5}$, $3.58\\cdot10^{-5}$ against $5\\cdot10^{-5}$; \\fn{fix9_review/recompute.json}).\n\n"
 "\\emph{$\\alpha$-equation.} $\\beta$ is nondecreasing in $\\Lambda$ and each $\\Lambda$-dependent term has degree $\\le1$ in $\\Lambda$ ($\\sqrt\\Lambda$ and $\\mathrm{QBL}\\cdot\\Lambda$, \\fn{closure.py}:22); with $\\Lambda\\le(1+2C'/L_F)\\ell$ the right side of the whole-tail condition (\\fn{closure.py}:185) grows by at most the factor $1+2C'/L_F$, which (C5) allows (margin $\\ge\\times1512$). Floor quotes: $\\sqrt{\\Lambda}\\le(1+2C'/L_F)^{1/2}\\sqrt{\\ell}$, relative change $\\le2.7\\cdot10^{-8}$.\n\n"
 "\\emph{Integrals.} The march uses $r_1$ at the point cell $\\mathrm{Ly}^{hi}$ with $L$ dropped (\\fn{closure.py}:222; \\fn{TAIL_LEMMA.tex} \\S V1), valid for every $\\Lambda\\le\\mathrm{Ly}^{hi}$, which (C2) keeps. Then $2\\theta\\int_a^b\\mathrm{NEGV}\\le Y(a)-Y(b)+2\\int_a^b(r_1^+\\Lambda+\\Phi+F_0)$, $Y(a)\\le\\Lambda(a)$, $-Y(b)\\le2C'$: the constant gains $C'$ (or $2C'$ if march and tail are bounded separately), the slope nothing. $\\kappa_J$ enters as $J_m(1+\\kappa_J)$ with $\\kappa_J=\\kappa_J^{\\rm cert}+10^{-9}$ (\\fn{fix6_v36/t2_A/rigor/hjw_rep.py}:38--39), so the pad covers $C'\\le10^{-9}X_{\\rm const}$: (C6), margin $\\ge\\times141$ ($\\times70$ for $2C'$).\n\n"
 "\\emph{Other consumers} (C7) \\tR: \\fn{mod_sups.py}:33--37 reads the march rows at $\\mathrm{Ly}^{hi}$ and one tail row at $L^*$ (scripted reserve $\\ge\\times2.5\\cdot10^8$); \\fn{fix6_v36/v36cfg.py}:66--84 rounds the $J$'s, slopes and peak up to 5 digits and $\\rho$ down; the raw $\\sqrt{L_F}$ in \\fn{assembly_arb.f4b} feeds only $I_{3c}$, $I_{3s}$, which v37 does not use. Hence all conclusions of the verifier hold with $\\Lambda\\le L+2C'$ (march), $\\Lambda\\le\\ell+2C'$ (tail), the $\\alpha$-barriers unchanged, and $\\int_IV_m^2\\le J_m(1+\\kappa_J)+s_mT$.\n"
 "\\end{proof}\n\n"
 "\\begin{definition}[bootstrap set of the flux lemma]\\label{F:def:boot}\n"
 "On $[\\tau_1,\\tau']$ the following are assumed; each is \\emph{strictly improved} by the conclusions obtained under it (the improvement factor, certified, is given):\n"
 "\\begin{enumerate}[label=(F\\arabic*),leftmargin=3em]\n"
 "\\item (B$\\Phi$): for $k=1,2$ and every window, $\\int_I\\!\\int_{\\ln28}^{\\ln350}\\varepsilon\\omega(z^A_{k+2})^2\\le2(C0_k+C1_kT)$, $C0_1\\le1.6637\\cdot10^{-4}$, $C1_1\\le7.6843\\cdot10^{-8}$, $C0_2\\le116.3782$, $C1_2\\le0.2606004$ (setting $4N_2^{\\rm ref}$). Improved by Lemma~S(A), Lemma~D: the bound holds without the factor 2 (factor $\\times2$).\n"
 "\\item $|(a,b)|\\le\\bar a:=10^{-6}$ (typed, used in $C_{\\rm eff}$ and $K_0^w$). Improved by the verifier: $\\sup|(a,b)|\\le2.14041\\cdot10^{-7}$ (\\fn{mod_sups}, $\\times4.67$).\n"
 "\\item $|\\alpha|\\le2.5\\cdot10^{-4}$ (coefficient tables, chart boxes; (B-a), (E0), (E2) of \\fn{ei_lemma_v4.tex}). Improved: $|\\alpha|\\le1.6\\cdot10^{-4}$ on every march cell and $\\le a^*$ in the tail (worst ratio $0.64$).\n"
 "\\item $\\max(|w|,|w'|)\\le10^{-4}$ on $[0,8]$ (E2 tip-jet box (E2)). Improved: asserted at every cell's $\\mathrm{Ly}^{hi}$ with worst ratio $0.3245$ (\\key{E2_tip_jet}).\n"
 "\\item the standing boxes on OR and the rest of (B-b)--(B-e) of \\fn{ei_lemma_v4.tex}: $|g|\\le S^V_0=0.0153$, $|g'|\\le S^V_1/\\sigma_{\\min}$, $|w|\\le\\frac12$ (ratios $0.7475$, $0.765$, $0.789$ for the asserted ones; $S^V_0,S^V_1$ improved by Lemma~S(A), $\\times2.116$ and $\\times5.745$ at $4N_2^{\\rm ref}$).\n"
 "\\item B1 in the form $\\Lambda\\le L+2C'$ on the march, $\\le\\ell+2C'$ in the tail. Improved: Proposition~\\ref{F:prop:abs} gives $Y\\le L$ strictly (pad and cushion, (C2), (C3)).\n"
 "\\end{enumerate}\n"
 "\\end{definition}\n\n")
assert src.count(old) == 1; R.append((old, new))
rep("(B$\\Phi$) is added to the bootstrap set. Under the bootstrap set on $[\\tau_1,\\tau']$:",
    "(B$\\Phi$) and the other items of Definition~\\ref{F:def:boot} are added to the bootstrap set. Under the bootstrap set on $[\\tau_1,\\tau']$:")
rep("The window integrals are continuous in the window end, so (B$\\Phi$) is closed and open in the sense of the closure proposition. Lemma~S(A) itself does not use (B$\\Phi$).",
    "The window integrals are continuous in the window end, so (B$\\Phi$) is closed (monotone convergence) and open (the extra piece $\\int_{\\tau'}^{\\tau''}$ is independent of the window start and tends to $0$, given local integrability of $\\varepsilon\\omega(z^A_{k+2})^2$ past $\\tau'$, Theorem U-1); (F2)--(F6) are closed conditions on continuous quantities and are improved strictly, so they are open by continuity (Theorem U-1). Lemma~S(A) itself does not use (B$\\Phi$).")
# status paragraph at the end
old_st = src[src.index("\\paragraph{Status.} Design."):src.index("\\end{document}")]
new_st = ("\\paragraph{Status (text v37).} Certified computation complete (\\fn{flux2_cert.json}, verdict PASS; all of (C0)--(C7) pass at both production points, (C9) by construction). "
 "Independently reviewed (\\fn{FLUX2_REVIEW.md}): sound modulo the listed items; the review notes that $\\bar K_0$ is conservative by $\\times3.4$ (a direct pair-form reduction gives $1.003\\cdot10^{-21}$ and a combined margin $\\times16.29$; not needed). "
 "Open: (C8) the regularity lemma U-1 (written: \\fn{fix10_proofs/U1_regularity.tex}, \\tHd, not reviewed); the written $k=2$ export (written: \\fn{fix10_proofs/G1_export_k2.tex}, not reviewed); the outputs \\fn{fix8_flux2_impl/out/*} are in no manifest and no kit; LaTeX compile.\n")
assert src.count(old_st) == 1; R.append((old_st, new_st))
rep("\\end{document}", "\\ifdefined\\stepmain\\else\n\\end{document}\n\\fi")
R.append(('\\providecommand{\\cst}[1]', '\\providecommand{\\fD}{\\mathfrak D}\\providecommand{\\fX}{\\mathfrak X}\\providecommand{\\fa}{\\mathfrak a}\\providecommand{\\OR}{\\mathrm{OR}}\\providecommand{\\eI}{\\varepsilon_I}\\providecommand{\\tHd}{\\textup{[hand]}}\\providecommand{\\tS}{\\textup{[sym]}}\\providecommand{\\tC}{\\textup{[cert]}}\\providecommand{\\tE}{\\textup{[est]}}\\providecommand{\\tR}{\\textup{[read]}}\n\\providecommand{\\cst}[1]'))
out = src
for o, n in R:
    c = out.count(o)
    if c != 1: sys.exit('replacement count %d for: %s' % (c, o[:100]))
    out = out.replace(o, n)
import re
out = re.sub(r'\\tH(?![a-zA-Z])', lambda m: '\\tHd', out)
open(B + 'fix12_texts/flux2_lemma_v37.tex', 'w').write(out)
print('ok', len(R))
