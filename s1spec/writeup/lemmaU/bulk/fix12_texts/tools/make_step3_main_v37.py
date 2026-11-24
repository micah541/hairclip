# make_step3_main_v37.py -- builds fix12_texts/step3_v37_main.tex from fix1/b12_text/step3_v33_main.tex (unchanged).
import sys
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
src = open(B + 'fix1/b12_text/step3_v33_main.tex').read()
R = []
def rep(o, n): R.append((o, n))
rep(src[:src.index('\\def\\stepmain{1}')],
    "% step3_v37_main.tex -- assembled step-3 document, text v37 (fix12_texts, 20 Nov 2026).  Built by tools/make_step3_main_v37.py from\n"
    "% fix1/b12_text/step3_v33_main.tex (unchanged).  \\input's: numbers_A_v37, lemmaS_A_v37, ei_lemma_v4 (with ei_constants_v37), flux2_lemma_v37,\n"
    "% continuation_v37.  Cited, not input: lemmaU_v37.tex, U1_regularity.tex, G1_export_k2.tex, E0_LEMMA_v4b.tex, erratum_v37.tex.  Not compiled.\n")
rep("\\input{b12_preamble}\n\\input{numbers_A_v33}", "\\input{b12_preamble_v37}\n\\input{numbers_A_v37}")
rep("\\title{Step 3 of the hair-clip proof, version 3.3:", "\\title{Step 3 of the hair-clip proof, text v37 (version 3.3 line, chain v36, verifier fin13v4b):")
rep("\\date{14 November 2026 (fix round 1). \\emph{Not compiled.} v3.1 and v3.2 files are unchanged; this is a new file.}",
    "\\date{20 November 2026 (text v37). \\emph{Not compiled.} v3.1, v3.2 and the v3.3 texts are unchanged; this is a new file.}")
old_box = src[src.index("\\begin{center}\\fbox{\\parbox{0.94\\linewidth}{\\small\\textbf{Status.}"):src.index("\\tableofcontents")]
new_box = ("\\begin{center}\\fbox{\\parbox{0.94\\linewidth}{\\small\\textbf{Status (text v37). Step~3 is NOT closed.} "
 "This document assembles step~3 after the fix rounds of Nov 16--20: the $E_I$ inequality was re-derived (\\fn{fix5/ei_spec/EI_SPEC.md}) and re-implemented as verifier \\fn{fin13v4b} (\\fn{step3nl_v4b}, \\S\\ref{EI:sec}); "
 "the E0 and E2 blocks were audited and fixed; Lemma~S(A) was rerun as chain v36 with the R3 forcing constant and the level-2 export $r_2=\\frac1{1000}$ (\\S\\ref{A:sec:status}); "
 "the outer-ramp flux is flux lemma v2 with time-integrated absorption (\\S\\ref{F:sec:top}); (P1$'$), (P2$'$) and the regularity of the energies are the lemma U-1 (\\fn{fix10_proofs/U1_regularity.tex}); (GS) is added to the bootstrap (\\S\\ref{C:sec:graph}). "
 "All defects found Nov 16--20 are in \\fn{fix12_texts/erratum_v37.tex}, section F. "
 "\\textbf{Open:} independent reviews of the \\fn{fin13v4b} delta, of the v36 chain, of U-1 and of G1; the Lean build and the $n=1$ instance; (H6) smoothness at $t_1$ (step~4 [A]); (GS) in the hand-off list; (P3$'$); the [read] and not-audited verifier inputs (\\S\\ref{EI:sec:open}); the step-4 re-pin (stages 8, 9, 18, 19, 20 of manifest\\_v33) and new kits; LaTeX compile. "
 "Lemma~U is cited (\\fn{fix12_texts/lemmaU_v37.tex}).}}\\end{center}\n")
assert src.count(old_box) == 1; R.append((old_box, new_box))
old_ver = src[src.index("\\section{The verifier (D1, D6)}\\label{M:sec:verifier}"):src.index("%======================================================================\n% ---------------- repaired Lemma S(A)")]
new_ver = ("\\section{The verifier (D1, D6)}\\label{M:sec:verifier}\n"
 "The inner and intermediate closure is verifier \\fn{fin13v4b} (\\fn{step3nl_v4b}, \\fn{V4B_REPORT.md}; production \\fn{runs/P}, double run byte-identical apart from timing). It replaces \\fn{fin12f4}, whose $E_I$ inequality failed audit (G0; erratum section F). "
 "Its $E_I$ block is Lemma~\\ref{EI:lem} below; the flux inputs are $\\fn{FLUX0}=1.64\\cdot10^{-20}$, $\\fn{FLUX2}=8.0\\cdot10^{-20}$, consumed only through the combined hypothesis (I3$'$), which flux lemma v2 proves (\\S\\ref{F:sec:top}). "
 "Outputs used downstream \\tC\\ (\\fn{rate_table_fin13v4b_a0.00015.json}, \\fn{modsups_fin13v4b}): peak $\\sqrt\\Lambda\\le\\peakV$, $\\rho\\ge\\rhoTailV$, $\\sqrt{L_F}\\le\\sqrtLFV$, $\\alpha_F\\le\\alphaFV$, whole-tail $\\alpha$-condition $\\times\\tailAlphaMarV$ ($\\times\\tailAlphaMarBV$ at $\\alpha_1=1.25\\cdot10^{-4}$); "
 "$\\bar a=\\bar b\\le\\abarV$, $\\bar\\beta\\le\\betabarV$; $(J_m,s_m)$ for $m=0,1,2,4$: $(\\JzeroV)$, $(\\JoneV)$, $(\\JtwoV)$, $(\\JfourV)$; $s_{I1}\\le\\sIoneV$, $s_{I2}\\le\\sItwoV$, $M_1^{\\rm ramp}\\le\\MoneRampV$, $S_{\\rm ramp}\\le\\SrampV$. "
 "Inputs from (A) (interface (I2)): $S^v_0\\le0.0153$, $S^v_1\\le0.6$ on $[34.25,60]$ \\tB. H-J-W in window form ($J_m(1+\\kappa_J)$, $\\kappa_J\\le\\kappaJV$, including the pad $10^{-9}$ for $C'$) is Lemma~\\ref{A:lem:HJ}. "
 "With the time-integrated flux part, B1 reads $\\Lambda\\le L(\\tau)+2C'$ (Proposition~\\ref{F:prop:abs}).\n\n"
 "\\input{ei_lemma_v4}\n\n")
assert src.count(old_ver) == 1; R.append((old_ver, new_ver))
rep("\\input{lemmaS_A_v33}", "\\input{lemmaS_A_v37}\n\n%======================================================================\n% ---------------- flux lemma v2 (final text) ----------------\n\\input{flux2_lemma_v37}")
rep("\\section{Lemma U (summary; proof in \\fn{lemmaU.tex})}\\label{M:sec:U}\n\\fn{bulk/t3_proof/lemmaU.tex} proves,",
    "\\section{Lemma U (summary; proof in \\fn{lemmaU.tex}, revision \\fn{fix12_texts/lemmaU_v37.tex})}\\label{M:sec:U}\n\\fn{bulk/t3_proof/lemmaU.tex} (revised in \\fn{lemmaU_v37.tex}: Lemma $\\delta U$ with $\\alpha_\\tau=\\tilde\\lambda\\alpha+\\beta$, $c_{\\delta U}\\le\\cdUV$, $\\times\\cdUMarV$ below the typed $2\\cdot10^{-8}$; Lemma~D for $k\\in\\{0,1,2\\}$) proves,")
rep("Certified values \\tC\\ (\\fn{t4_cert/out/U_cert_v35.json}, sha256 \\fn{960809fdc18dd67c1770fa5b77b2658081d0c81cb24f1b246f48ca4a583d5f54}):",
    "Certified values \\tC\\ (\\fn{fix6_v36/t4_cert/out/U_cert_v36_pass2.json}, sha256 \\fn{d512596a...}; pass 3 equal apart from args):")
rep("\\fn{bulk/t5_zone3/Z5_REPORT.md}, \\fn{out/zone3_interface.json} (N1 rerun, v35, sha256 \\fn{45e8d4dcaff2bd75fcea725f80e98625af70a40709fc2cf5d3df78b5becb3ea3}; $\\Delta^R$, $n^{\\rm III}$ and the boxes identical to v33 and v34, only the H6 entries differ in the sixth digit)",
    "\\fn{bulk/t5_zone3/Z5_REPORT.md}, \\fn{fix6_v36/t5_zone3/out/zone3_interface.json} (v36, sha256 \\fn{d6813ddb...}; $\\Delta^R$ and the boxes identical to v35, $n^{\\rm III}$ and the H6 entries differ by less than $2\\cdot10^{-5}$ relative and stay below the quoted values, \\fn{V36_REPORT.md} \\S4)")
rep("\\item[B1] $|\\alpha|$, $\\sqrt\\Lambda$ below the barrier of the verifier (\\fn{fin12f4}).",
    "\\item[B1] $|\\alpha|$ below the barrier and $\\Lambda\\le L(\\tau)+2C'$ (march), $\\le\\ell(\\tau)+2C'$ (tail) for the verifier \\fn{fin13v4b}, with the box hypotheses (B-a)--(B-f), (E0), (E2) of \\S\\ref{EI:sec:hyp}.\n"
    "\\item[BF] the flux bootstrap items (F1) (B$\\Phi$) and (F2) $\\bar a\\le10^{-6}$ of Definition~\\ref{F:def:boot}, each strictly improved.\n"
    "\\item[GS] the tip arc is a normal graph over $S_1$ on $|x|\\le60$ (Definition~\\ref{C:def:GS}).")
rep("D1 verifier & B1; B2$'$ ($S^v_0,S^v_1$; fluxes at $\\beta_A\\Eh$) & $J_m,s_m$ (H-J-W); $\\bar a,\\bar b,\\bar\\beta$; ramp values\\\\",
    "D1 verifier & B1, GS; B2$'$ ($S^v_0,S^v_1$); (I3$'$) from flux lemma v2 under BF & $J_m,s_m$ (H-J-W); $\\bar a,\\bar b,\\bar\\beta$; ramp values; $\\Lambda\\le L+2C'$\\\\\n"
    "flux v2 & (A)'s $\\bar N_l$, $\\bar M_l$ at $\\beta_A\\Eh$; BF & (I3$'$); BF strictly\\\\")
rep("G & all boxes above; Lemma~O graph persistence & G strictly (Lemma~\\ref{C:lem:G})\\\\\\bottomrule",
    "G & all boxes above; Lemma~O graph persistence & G strictly (Lemma~\\ref{C:lem:G})\\\\\nGS & B1 boxes; U-1 & GS (Lemma~\\ref{C:lem:GS}; strictness is an unreviewed sketch)\\\\\\bottomrule")
rep("fluxes at $\\beta_A=2$ $\\times\\MarFluxZeroV$, $\\times\\MarFluxTwoV$ \\tB\\ (G-A1 applied);",
    "flux lemma v2: combined pointwise margin $\\times\\FluxCombMarV$, conditions (C0)--(C7) of Proposition~\\ref{F:prop:abs} \\tC\\ (the old flux-block margins $\\times\\MarFluxZeroV$, $\\times\\MarFluxTwoV$ are records only);")
rep("\\item $\\mathcal A(2)\\le\\Eh=\\EhV$;", "\\item $\\mathcal A(2)\\le\\Eh=\\EhV$; level-2 export $C0_2\\le\\CzeroTwoRefFourV$, $C1_2\\le\\ConeTwoRefFourV$;")
rep("Under (K1)--(K6) $\\mathrm{BS}_U$ holds on $[\\tau_1,\\tau_{\\max})$: implication step, nonempty, closed, open, as proved in \\fn{lemmaU.tex} (Proposition ``closure''), with G added by Lemma~\\ref{C:lem:G}.",
    "Under (K1)--(K6) $\\mathrm{BS}_U$ holds on $[\\tau_1,\\tau_{\\max})$: implication step, nonempty, closed, open, as proved in \\fn{lemmaU.tex} (Proposition ``closure''), with G added by Lemma~\\ref{C:lem:G}, GS by Lemma~\\ref{C:lem:GS}, and BF by Remark~\\ref{F:rem:boot}.")
rep("\\input{continuation_v33}", "\\input{continuation_v37}")
old_rem = src[src.index("\\section{What remains, and the B0 slots}\\label{M:sec:remains}"):src.index("\\begin{thebibliography}")]
new_rem = ("\\section{What remains}\\label{M:sec:remains}\n\\begin{enumerate}[leftmargin=*]\n"
 "\\item \\textbf{Hypotheses of the theorem} \\tH: (H1)--(H6) and (H2f) (step 4, values unchanged; step~4 not re-pinned on v36); (GS) at $\\tau_1$ to be added to the hand-off list; the verifier's tail lemma (W2).\n"
 "\\item \\textbf{Reviews outstanding}: the \\fn{fin13v4b} delta (E0 F-B, E2 fix, R-2, R-3); the v36 chain as a whole; U-1 and G1 (\\fn{fix10_proofs}); a final review of the assembled step~3.\n"
 "\\item \\textbf{Verifier inputs} not audited or [read] (\\S\\ref{EI:sec:open}); the float check of the $\\beta$ decomposition.\n"
 "\\item \\textbf{Literature / Lean}: the Lean build at commit \\fn{d72b744} and the $n=1$ instance; DKS (D1) is a step inside a proof; (P3$'$).\n"
 "\\item \\textbf{Records}: \\fn{fix8_flux2_impl/out/*} in no manifest; step-4 stages 8, 9, 18, 19, 20 and a new driver/manifest; kits.\n"
 "\\item \\textbf{Not claimed}: convergence to $S_1$ (upgrade C).\n\\end{enumerate}\n\n")
assert src.count(old_rem) == 1; R.append((old_rem, new_rem))
rep("\\bibitem{DKS2002}", "\\bibitem{CWFlow} M. Warren, Lean development \\fn{HSLflow-Whitney}, \\fn{CWFlow/}, commit \\fn{d72b7447edb0c830331a76152559902c23a4092e}.\n\\bibitem{DKS2002}")
rep("$\\Delta^R\\le(4.4986\\cdot10^{-8},6.7955\\cdot10^{-6},0.011480)$", "$\\Delta^R\\le(4.4986\\cdot10^{-8},6.7972\\cdot10^{-6},0.011480)$ (v36 \\key{DeltaR}; the v33--v35 texts printed $6.7955\\cdot10^{-6}$, below the certified $6.79715\\cdot10^{-6}$ of both v35 and v36, a text slip)")
out = src
for o, n in R:
    c = out.count(o)
    if c != 1: sys.exit('replacement count %d for: %s' % (c, o[:100]))
    out = out.replace(o, n)
open(B + 'fix12_texts/step3_v37_main.tex', 'w').write(out)
print('ok', len(R))
