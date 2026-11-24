# make_continuation_v37.py -- builds fix12_texts/continuation_v37.tex from fix1/b12_text/continuation_v33.tex (unchanged) by exact,
# asserted-unique replacements.  Changes: (GS) added to the bootstrap; (P1')/(P2') cited as in fix10_proofs/U1_regularity.tex
# (Lean CWFlow for existence and uniqueness, DKS Thm 3.1 proof step (D1) for the criterion); fin13v4b / v36 numbers.
import sys
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
src = open(B + 'fix1/b12_text/continuation_v33.tex').read()
R = []
def rep(o, n): R.append((o, n))
rep("% continuation_v33.tex -- step 3 v3.3: from the bootstrap to the theorem (fix round 1, review finding B2).",
    "% continuation_v37.tex -- step 3, text v37: from the bootstrap to the theorem.  Built by fix12_texts/tools/make_continuation_v37.py\n"
    "% from fix1/b12_text/continuation_v33.tex (unchanged).  New: (GS) in the bootstrap; (P1'), (P2') from fix10_proofs/U1_regularity.tex; v36/fin13v4b numbers.")
rep("% Standalone, or \\input by step3_v33_main.tex", "% Standalone, or \\input by step3_v37_main.tex")
rep("\\input{numbers_A_v33}\n\\title{Step 3, v3.3: continuation", "\\input{numbers_A_v37}\n\\title{Step 3, text v37: continuation")
rep("\\author{}\\date{14 November 2026 (fix round 1, review finding B2). \\emph{Not compiled.}}",
    "\\author{}\\date{20 November 2026 (text v37). \\emph{Not compiled.}}")
rep("\\noindent\\emph{Labels} as in \\fn{lemmaS_A_v33.tex}", "\\noindent\\emph{Labels} as in \\fn{lemmaS_A_v37.tex}")
rep("B2$'$ and the repaired Lemma~S(A) (\\fn{lemmaS_A_v33.tex})", "B2$'$ and the repaired Lemma~S(A) (\\fn{lemmaS_A_v37.tex})")
rep("provided the continuation criterion (P1) is used in its $L^2$ form (the blow-up criterion of Dziuk--Kuwert--Sch\\\"atzle \\cite{DKS2002}, Theorem~3.1, and of Chou \\cite{Chou2003}; \\textbf{[to be confirmed by the authors]}; Proposition~\\ref{C:prop:global}).",
    "provided the continuation criterion (P1) is used in its $L^2$ form (Proposition~\\ref{C:prop:global}). In text v37 (P1$'$) and (P2$'$) are no longer bare literature items: they are the regularity lemma U-1 and Proposition~(P2$'$) of \\fn{fix10_proofs/U1_regularity.tex} \\tH, whose inputs are Micah Warren's Lean theorems \\fn{exists_flow_of_immersion}, \\fn{flow_unique_of_immersion} \\cite{CWFlow} (existence and uniqueness) and the smooth-extension step (D1) inside the proof of Theorem~3.1 of Dziuk--Kuwert--Sch\\\"atzle \\cite{DKS2002} (the criterion); see \\S\\ref{C:sec:cont}.")
rep("\\SU, \\SZ, \\SV\\ denote \\fn{t4_cert/out/U_cert_v35.json} (sha256 \\fn{960809fd...}), \\fn{t5_zone3/out/zone3_interface.json} (\\fn{45e8d4dc...}, the N1 rerun) and the verifier \\fn{step3nl_v33} (\\fn{V33_VERIFIER.md}, tag \\fn{fin12f4}).",
    "\\SU, \\SZ, \\SV\\ denote \\fn{fix6_v36/t4_cert/out/U_cert_v36_pass2.json} (sha256 \\fn{d512596a...}), \\fn{fix6_v36/t5_zone3/out/zone3_interface.json} (\\fn{d6813ddb...}) and the verifier \\fn{step3nl_v4b} (\\fn{V4B_REPORT.md}, tag \\fn{fin13v4b}, \\fn{runs/P}). The v35 files and \\fn{fin12f4} are superseded (erratum section F).")
# (GS) definition and lemma, inserted after Lemma G's hand-off paragraph
anchor = "the wording of (H4) must be extended to $y\\in[48,700]$ (\\S\\ref{C:sec:gaps}, item 2).\n"
gs = (anchor +
 "\n\\begin{definition}[GS: graph over $S_1$ near the tip]\\label{C:def:GS}\n"
 "$\\mathrm{GS}(\\tau)$: there is $u\\in C^1([-X_S,X_S])$, $X_S=60+\\frac14$, with $\\sup|k\\,u|\\le\\frac12$ and $\\sup|u'|\\le1$, such that $x\\mapsto\\Psi(x,u(x))$ (normal coordinates of $S_1$) reparametrises a sub-arc of $\\Gamma_{\\rm tip}(\\tau)$ containing the tip. "
 "This is the chart in which the verifier's $w$ and $E_I$ (on $x\\in[10.25,60]$) are defined (\\fn{U1_regularity.tex}, Definition ``charts'', item (GS)).\n"
 "\\end{definition}\n"
 "\\begin{lemma}[GS is part of the bootstrap]\\label{C:lem:GS}\n"
 "Add $\\mathrm{GS}$ to $\\mathrm{BS}_U$ (next to G). The extended set is nonempty, closed and open in the sense of the closure proposition of \\fn{lemmaU.tex}.\n"
 "\\end{lemma}\n"
 "\\begin{proof}\\tH\\ \\emph{At $\\tau_1$}: (GS) on $|x|\\le60$ is implicit in step~4's certification of $\\Lambda(\\tau_1)$ (the functional $E_I$ is defined through it); the hand-off list must state it (\\S\\ref{C:sec:gaps}, item 2). "
 "\\emph{Strictness} (sketch added in text v37; it is not in \\fn{U1_regularity.tex} or any report, and must be reviewed): under the bootstrap, $|w|\\le0.5$ on $[34.25,60]$ and $\\le5.91\\cdot10^{-4}$ on $[16.25,34.25]$, $|w'|\\le0.0153$ on $[34.25,60]$ and $\\le1.1606\\cdot10^{-3}$ on $x\\ge16.25$, $\\max(|w|,|w'|)\\le7.48\\cdot10^{-4}$ on $[10.25,16.25]$ (verifier box hypotheses B-b to B-e, asserted in \\fn{closure.box_assert} of \\fn{fin13v4b}), and $|w|,|w'|\\le10^{-4}$ on $[0,8]$ (hypothesis H-E2b, asserted); with $u=u_0+w$, $|u_0|\\le|\\alpha|\\sup|\\hat\\psi_3|+\\alpha^2\\sup|\\hat h_2|$, these should give $|ku|<\\frac12$ and $|u'|<1$; the sup bounds of $\\hat\\psi_3$, $\\hat h_2$ and of $w$ on $[8,10.25]$ needed for this are not collected anywhere yet. "
 "\\emph{Open}: Lemma ``smooth dependence of the graph functions'' of \\fn{U1_regularity.tex} (the set where (GS) holds strictly is open and $u$ is $C^\\infty$ in $(t,\\lambda,\\xi,x)$). "
 "\\emph{Closed}: the compactness step (c) of Proposition~(P2$'$) of \\fn{U1_regularity.tex} (limits of graph functions converge in $C^1$ with the same bounds). \\emph{Note:} no script checks the strictness implication.\n"
 "\\end{proof}\n")
rep(anchor, gs)
# bounds table numbers (fin13v4b)
rep("B1: $\\Lambda\\le\\Lambda_{\\rm bar}$, so $\\|\\tilde w\\|_\\rho\\le5.31\\cdot10^{-8}$, $\\|\\tilde w''\\|_\\rho\\le\\sqrt{\\Lambda/\\varepsilon}\\le1.68\\cdot10^{-5}$ ($\\varepsilon=10^{-5}$); $\\max(|w|,|w'|)\\le7.48\\cdot10^{-4}$ on $[10.25,16.25]$ (\\SV\\ \\key{S_ramp})",
    "B1: $\\Lambda\\le\\Lambda_{\\rm bar}$, so $\\|\\tilde w\\|_\\rho\\le5.31\\cdot10^{-8}$, $\\|\\tilde w''\\|_\\rho\\le\\sqrt{\\Lambda/\\varepsilon}\\le1.68\\cdot10^{-5}$ ($\\varepsilon=10^{-5}$); $\\max(|w|,|w'|)\\le7.48\\cdot10^{-4}$ on $[10.25,16.25]$ (\\SV\\ \\key{S_ramp_max_along_trajectory}$\\le\\SrampV$; peak $\\sqrt\\Lambda\\le\\peakV$)")
rep("D1: $|w'|\\le1.1604\\cdot10^{-3}$, $|w''|\\le1.5457\\cdot10^{-2}$ (\\SV\\ \\key{sI1}, \\key{sI2})",
    "D1: $|w'|\\le\\sIoneV$, $|w''|\\le\\sItwoV$ (\\SV\\ \\key{far_max.sI1}, \\key{sI2})")
# continuation standing facts
oldP = src[src.index("\\noindent\\textbf{Standing facts} \\tH\\ (literature;"):src.index("CDF decreases length:")]
newP = ("\\noindent\\textbf{Standing facts.} (P1$'$) and (P2$'$) are proved in \\fn{fix10_proofs/U1_regularity.tex} \\tH\\ from the cited inputs below; (P3$'$) remains a literature item \\tH.\n"
 "\\begin{enumerate}[label=(P\\arabic*$'$),leftmargin=3em]\n"
 "\\item \\emph{Local well-posedness, uniqueness, and the blow-up criterion.} For a smooth closed immersed plane curve, CDF has a unique smooth solution on a maximal interval $[t_1,T_{\\max})$ (Lemma ``maximal smooth solution'' of \\fn{U1_regularity.tex}), and if $T_{\\max}<\\infty$ then $\\sup_{t<T_{\\max}}\\int_{\\gamma(t)}\\kappa^2ds=\\infty$ (Lemma ``blow-up criterion, the form needed'' there). Inputs:\n"
 "  \\begin{itemize}[leftmargin=*]\n"
 "  \\item \\emph{existence} (L1) \\fn{exists_flow_of_immersion} and \\emph{uniqueness} (L2) \\fn{flow_unique_of_immersion} of the Lean development \\fn{HSLflow-Whitney/CWFlow} \\cite{CWFlow} (commit \\fn{d72b7447edb0c830331a76152559902c23a4092e}, files \\fn{Flow/Main.lean}, \\fn{Flow/UniqueMain.lean}), instantiated at $E=\\mathbb R^1$, $M=\\mathbb S^1$, $V=\\mathbb C$, $J=i$; at $n=1$ the Lean velocity is $-k_{ss}N$ with no tangential part \\tS\\ (\\fn{U1_sym.py} S1--S4). The statements were read from source; the build was not run, and the $n=1$ instance has not been elaborated in Lean;\n"
 "  \\item \\emph{criterion} (D1): ``bounded $\\int\\kappa^2ds$ on $[0,T)$, $T<\\infty$, implies smooth extension to $[0,T]$'', a step inside the proof of Theorem~3.1 of Dziuk--Kuwert--Sch\\\"atzle \\cite{DKS2002} (\\fn{writeup/lit_P1/P1_CONFIRMATION.md}); the rate of Theorem~3.1, $\\int\\kappa^2ds\\ge c(T-t)^{-1/4}$, and Chou \\cite{Chou2003} are not needed;\n"
 "  \\item \\emph{smoothness at $t_1$} (H6/S4.1): step~4's statement ``smooth for $t>0$ by parabolic regularity'' (\\fn{step4_proof.tex} l.672, labelled [A] there), not proved in step~3.\n"
 "  \\end{itemize}\n"
 "\\item \\emph{Modulation} (P2$'$): Proposition ``(P2$'$): modulation'' of \\fn{U1_regularity.tex}: on the closed set where (GS) and the modulation validity box hold, the modulation ODE has a $C^1$ solution, the constraints are conserved, and it continues until $t\\to T_{\\max}$ or the box fails \\tH\\ (matrix invertibility from \\fn{MOD_IBP.md}, Neumann bound $\\nu\\le0.0466$, \\tR). Theorem U-1 there also gives that $\\Lambda$ and every energy used by the verifier, by Lemma~S(A) and by the flux lemma are $C^1$ in $\\tau$, with their identities at every $\\tau$.\n"
 "\\item \\emph{Continuity and interior parabolic regularity}: on a space--time region where the curve is a graph with bounded slope, bounded curvature and $\\kappa_s$ bounded in $L^2$, all derivatives are bounded on interior subregions \\tH\\ (literature; reference to be supplied by the authors).\n"
 "\\end{enumerate}\n")
assert src.count(oldP) == 1; R.append((oldP, newP))
rep("\\begin{proof}\\tP\\ Suppose $\\tau_{\\max}<\\infty$. By B1 and (P2$'$) the modulation exists",
    "\\begin{proof}\\tP\\ Suppose $\\tau_{\\max}<\\infty$. By B1, GS (Lemma~\\ref{C:lem:GS}) and (P2$'$) the modulation exists")
rep("the closure proposition of \\fn{lemmaU.tex} with Lemma~\\ref{C:lem:G} then gives $\\mathrm{BS}_U$ and G for all $\\tau\\ge\\tau_1$.",
    "the closure proposition of \\fn{lemmaU.tex} with Lemmas~\\ref{C:lem:G} and \\ref{C:lem:GS} then gives $\\mathrm{BS}_U$, G and GS for all $\\tau\\ge\\tau_1$.")
rep("$\\tau_{\\max}=\\infty$, and $\\mathrm{BS}_U$ together with $\\mathrm{G}$ holds for all $\\tau\\ge\\tau_1$.",
    "$\\tau_{\\max}=\\infty$, and $\\mathrm{BS}_U$ together with $\\mathrm{G}$ and $\\mathrm{GS}$ holds for all $\\tau\\ge\\tau_1$.")
# conclusion numbers
rep("Unchanged from v3.2 \\S5.1--5.2, with the \\fn{fin12f4} verifier values \\tC\\ (\\SV): $\\bar a,\\bar b\\le2.19103\\cdot10^{-7}$, $\\bar\\beta\\le3.33811\\cdot10^{-9}$.",
    "Unchanged from v3.2 \\S5.1--5.2, with the \\fn{fin13v4b} verifier values \\tC\\ (\\SV\\ \\fn{modsups_fin13v4b}): $\\bar a,\\bar b\\le\\abarV$, $\\bar\\beta\\le\\betabarV$.")
rep("\\item[(W1) trapping] $|\\alpha|\\le1.5\\cdot10^{-4}$, $\\sqrt\\Lambda\\le5.31\\cdot10^{-8}$; after the tail start the barrier decays at rate $\\rho\\ge2.6930\\cdot10^{-3}$ to the floor $\\sqrt{\\Lambda_F}\\le9.872\\cdot10^{-10}$, and eventually $|\\alpha|\\le\\alpha_F=6.8666\\cdot10^{-7}$ \\tC\\ (\\SV, fin12f4 table; review B3: the fin11 values $4.94\\cdot10^{-10}$, $3.43\\cdot10^{-7}$ do not apply to v3.3); $A\\le5.31\\cdot10^{-8}$, $W_1\\le6.682\\cdot10^{-7}$; $\\bar a,\\bar b\\le2.19103\\cdot10^{-7}$, $\\bar\\beta\\le3.33811\\cdot10^{-9}$.",
    "\\item[(W1) trapping] $|\\alpha|\\le1.5\\cdot10^{-4}$, $\\sqrt\\Lambda\\le5.31\\cdot10^{-8}$ (peak $\\peakV$); after the tail start the barrier decays at rate $\\rho\\ge\\rhoTailV$ to the floor $\\sqrt{\\Lambda_F}\\le\\sqrtLFV$, and eventually $|\\alpha|\\le\\alpha_F=\\alphaFV$ \\tC\\ (\\SV, \\fn{rate_table_fin13v4b_a0.00015.json}; the whole-tail $\\alpha$-condition holds with margin $\\times\\tailAlphaMarV$, review R-3; the \\fn{fin12f4} values $2.6930\\cdot10^{-3}$, $9.872\\cdot10^{-10}$, $6.8666\\cdot10^{-7}$ are superseded). "
    "With the time-integrated outer-ramp flux (flux lemma v2, \\fn{flux2_lemma_v37.tex}) these hold for $Y=\\Lambda-Z$, and $\\Lambda\\le L(\\tau)+2C'$, $2C'\\le5.305\\cdot10^{-26}$ (relative $2.7\\cdot10^{-8}$ at the floor); $A\\le5.31\\cdot10^{-8}$, $W_1\\le6.683\\cdot10^{-7}$; $\\bar a,\\bar b\\le\\abarV$, $\\bar\\beta\\le\\betabarV$.")
rep("\\item[(W2) tip] $\\max(|w|,|w'|)\\le7.48\\cdot10^{-4}$ on $[10.25,16.25]$; on $x\\ge16.25$: $|w'|\\le1.1604\\cdot10^{-3}$, $|w''|\\le1.5457\\cdot10^{-2}$ \\tC\\ (\\SV).",
    "\\item[(W2) tip] $\\max(|w|,|w'|)\\le7.48\\cdot10^{-4}$ on $[10.25,16.25]$; on $x\\ge16.25$: $|w'|\\le\\sIoneV$, $|w''|\\le\\sItwoV$ \\tC\\ (\\SV).")
# gaps
rep("\\item \\textbf{(P1) in $L^2$ form} \\tH. The argument needs the $L^2$ blow-up criterion for closed immersed curves, Dziuk--Kuwert--Sch\\\"atzle \\cite{DKS2002} Theorem~3.1 and/or Chou \\cite{Chou2003} \\textbf{[to be confirmed by the authors]}. This is a literature item, as (P1) was in v3.2; the change is that v3.2's wording (curvature and first derivatives) must be replaced. No quantitative gap.",
    "\\item \\textbf{(P1$'$), (P2$'$)} \\tH: proved in \\fn{U1_regularity.tex} from (L1), (L2) \\cite{CWFlow} and DKS (D1) \\cite{DKS2002}. Open there: the Lean build at commit \\fn{d72b744} was not run and the $n=1$ instance ($E=\\mathbb R^1$, $M=\\mathbb S^1$, $V=\\mathbb C$) was not elaborated (suggested check for Micah: a short Lean file instantiating both theorems); the gluing lemma and the closure of \\fn{IsFlow} under the reflection are written mathematics, not Lean; (D1) is a step inside the proof of DKS Theorem~3.1, not its statement; (H6) smoothness at $t_1$ is step~4's [A] item. No independent review of U-1 exists yet.")
rep("\\item \\textbf{(H4) range.} G at $\\tau_1$ needs the graph property on $y\\in[14,700]$",
    "\\item \\textbf{(GS) at $\\tau_1$} must be added to the hand-off list (it is implicit in step~4's certification of $\\Lambda(\\tau_1)$; \\fn{PROOFS_REPORT.md} gap~4). The strictness step of Lemma~\\ref{C:lem:GS} is a hand reading of asserted boxes.\n"
    "\\item \\textbf{(H4) range.} G at $\\tau_1$ needs the graph property on $y\\in[14,700]$")
rep("\\bibitem{Chou2003} K.-S. Chou, A blow-up criterion for the curve shortening flow by surface diffusion, Hokkaido Math.\\ J.\\ 32(1) (2003) 1--19. (Entry \\fn{Chou2003} of \\fn{writeup/main/main.bib}.)\n",
    "\\bibitem{Chou2003} K.-S. Chou, A blow-up criterion for the curve shortening flow by surface diffusion, Hokkaido Math.\\ J.\\ 32(1) (2003) 1--19. (Entry \\fn{Chou2003} of \\fn{writeup/main/main.bib}.)\n"
    "\\bibitem{CWFlow} M. Warren, Lean development \\fn{HSLflow-Whitney}, directory \\fn{CWFlow/}, \\fn{github.com/micah541/HSLflow-Whitney}, commit \\fn{d72b7447edb0c830331a76152559902c23a4092e} (19 Nov 2026); theorems \\fn{exists_flow_of_immersion}, \\fn{flow_unique_of_immersion}.\n")
R.append(('\\input{b12_preamble}\n', '\\input{b12_preamble_v37}\n'))
out = src
for o, n in R:
    c = out.count(o)
    if c != 1: sys.exit('replacement count %d for: %s' % (c, o[:100]))
    out = out.replace(o, n)
open(B + 'fix12_texts/continuation_v37.tex', 'w').write(out)
print('ok', len(R))
