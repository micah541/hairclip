#!/usr/bin/env python3
"""Stage X7 (FLUX2_DESIGN.md sec. 10): texts generated from the certified outputs (no hand-typed numbers).
Reads fix7_flux2/flux2_lemma.tex (read-only) and out/{flux2_coef, chart4_shell, flux2_data, flux2_cert}.json; writes, in texts/ (new files):
  flux2_lemma_v2.tex          the design statements with every \\cst{...} placeholder filled, plus an implementation section (deviations);
  absorption_lemma_v2.tex     the reflected-integral lemma and the proposition with the certified C0-C7 table;
  erratum_flux2.tex           erratum entries (cutoff shape, jet box, drift double use, trace route, W(57) blend);
  BPhi_continuation_v36.tex   the (B Phi) paragraph for the continuation text (v36 suffix);
  numbers_flux2.json          every number used, with its source key.
Balance check: braces and \\begin/\\end pairs per file (no LaTeX on this machine; not compiled).  Output out/gen_texts.json."""
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
T_ = FC.Timer()
TX = os.environ.get('F2_TX', os.path.join(FC.ROOT, 'texts')); os.makedirs(TX, exist_ok=True)
J = {n: json.load(open(os.path.join(FC.OUT, n + '.json'))) for n in ('flux2_coef', 'chart4_shell', 'flux2_data', 'flux2_cert')}
c2 = J['flux2_coef']['res']; c3 = J['chart4_shell']['res']; c4 = J['flux2_data']['res']; c5 = J['flux2_cert']['res']; cond = J['flux2_cert']['conditions']
NUM = {}


def g(x, key, nd=4, dirn='up'):
    """decimal with nd digits after the point, rounded UP (upper bounds) or DOWN (lower bounds: margins, s)"""
    NUM[key] = dict(value=x, rounded=dirn)
    if isinstance(x, str): return x
    from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR
    d = Decimal(repr(x)); e = d.adjusted(); m = d.scaleb(-e).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_CEILING if dirn == 'up' else ROUND_FLOOR)
    if abs(m) >= 10: m = (m / 10).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_CEILING if dirn == 'up' else ROUND_FLOOR); e += 1
    return str(m) if e == 0 else '%s\\cdot10^{%d}' % (m, e)


fill = {
    'kappa0': '(' + ', '.join(g(r['up'], 'kappa0_%d' % l, 5) for l, r in enumerate(c2['kappa0'])) + ')',
    'coef.OR': '$\\bar c=%s$, $\\bar c_1=%s$, $d_{010}=%s$, $d_{110}=%s$, $d_{300}=%s$, $A_1=%s$, $B_1=%s$' % (
        g(c2['CBAR_OR']['up'], 'CBAR_OR'), g(c2['C1_OR']['up'], 'C1_OR'), g(c2['coef_OR']['d010']['up'], 'd010'), g(c2['coef_OR']['d110']['up'], 'd110'),
        g(c2['coef_OR']['d300']['up'], 'd300'), g(c2['coef_OR']['A1']['up'], 'A1'), g(c2['coef_OR']['B1']['up'], 'B1')),
    'data.Nbar': '$\\bar N=(%s, %s, %s)$, $\\bar M=(%s, %s, %s)$' % tuple([g(c4['Nbar'][str(l)]['up'], 'Nbar_%d' % l) for l in range(3)] + [g(c4['Mbar'][str(l)]['up'], 'Mbar_%d' % l) for l in range(3)]),
    'data.m3': '$m_{3c}=%s$, $m_{3s}=%s$' % (g(c4['m3c']['up'], 'm3c'), g(c4['m3s']['up'], 'm3s')),
    'data.m4': '$m_{4c}=%s$, $m_{4s}=%s$' % (g(c4['m4c']['up'], 'm4c'), g(c4['m4s']['up'], 'm4s')),
    'chart4.sup': '$\\bar A_4=%s$, $\\max_i\\bar B_{4i}=%s$, $\\bar B_{422}=%s$' % (g(c3['Abar'][4]['up'], 'Abar_4'), g(max(c3['Bbar'][k]['up'] for k in ('42', '43', '44')), 'Bbar_4max'), g(c3['B422']['up'], 'B422')),
    'chart4.beta44': g(c3['Bbar']['44']['up'], 'Bbar_44'),
    'phi2.C': g(c5['Cpp']['up'], 'Cpp'), 'phi2.S': g(c5['Spp']['up'], 'Spp'),
    'slack': '$s=%s$, $S\'=%s$, $C\'=%s$' % (g(c5['slack_s']['dn'], 'slack_s', 4, 'dn'), g(c5['Sp']['up'], 'Sp'), g(c5['Cp']['up'], 'Cp')),
}
src = open(os.path.join(FC.DESIGN, 'flux2_lemma.tex')).read()
out = src
for k, v in fill.items():
    out = out.replace('\\cst{%s}' % k, '\\ensuremath{%s}' % v.replace('$', '') if k in ('kappa0', 'chart4.beta44', 'phi2.C', 'phi2.S') else v)
assert '\\cst{' not in out.replace('\\cst{...}', '').replace('\\newcommand{\\cst}', ''), 'unfilled placeholder'
out = out.replace('DESIGN STAGE (fix round 7, 20 Nov 2026): statements and proof sketches.  Nothing is certified.',
                  'IMPLEMENTED (fix round 8, fix8_flux2_impl, 20 Nov 2026): numbers filled by code/gen_texts.py from out/*.json.')
out = out.replace('\\date{20 November 2026 (fix round 7, design). \\emph{Not compiled. Not certified.}}',
                  '\\date{20 November 2026 (fix round 8, implementation). \\emph{Not compiled.} Numbers [cert] from \\fn{fix8_flux2_impl/out/}.}')
impl = r'''
\section{Implementation notes (fix round 8)}\label{F:sec:impl}
Certified values (python-flint arb, outward rounding; \fn{fix8_flux2_impl/out/flux2_cert.json}):
$\bar K_0=%(K0)s$ against $\fn{FLUX0}=1.64\cdot10^{-20}$; $K_2^{\rm pt}=%(K2pt)s$; combined pointwise margin $(\fn{FLUX0}+\eI\fn{FLUX2})/(\bar K_0+\eI K_2^{\rm pt})\ge%(marg)s$;
$C''=%(Cpp)s$, $S''=%(Spp)s$, $C'=%(Cp)s$, $S'=%(Sp)s$, $s=%(s)s$; $W(57)\le%(W57)s$.
Deviations from the design, each a valid alternative:
\begin{enumerate}[label=(\roman*)]
\item The conversions of Lemma~\ref{F:lem:data}(b) are used in Minkowski form, $M_3\le\sigma_{\min}^{-5/2}(\|D^3g\|+3\|D^2g\|+2\|Dg\|)$, $M_4\le\sigma_{\min}^{-7/2}(\|D^4g\|+6\|D^3g\|+11\|D^2g\|+6\|Dg\|)$ (and $M_2\le\sigma_{\min}^{-3/2}(\|D^2g\|+\|Dg\|)$), and sums of window norms $\sum_iw_i(a_i+b_iT)^{1/2}$ are squared with fixed rational weights $\theta_i$, $\sum\theta_i\le1$.
\item $\|g''^2\|^2_{L^2(S\times I;ds\,d\tau)}\le\sigma_{\min}^{-1}\bar M_2^2\int_Is_2^2$ (the factor $\sigma_{\min}^{-1}$ from $ds=d\sigma/\sigma$).
\item $W$ is the verifier's weight ($S_7$ blend on $[30.25,38.25]$, $W=(1+x^2)^{-2}$ below): $W(57)\le%(W57)s$ by a validated integral. The v3.1 value $1.43716\cdot10^{-18}$ (quintic blend) is not used.
\item The coefficient sups on $\OR$ are the certified C2 tables of the verifier restricted to the boxes meeting $(57,60)$ (\fn{c2_t3d.json}, sub-region $[57,60]$), each asserted $\le$ the zone value.
\item $|(a,b)|\le\bar a:=10^{-6}$ is taken as a typed hypothesis ($\ge2.19103\cdot10^{-7}$, every level in use).
\item The level-4 chart structure of Lemma~\ref{F:lem:chart4} is asserted in sympy (\fn{chart4_shell.py}): $b_2^2=(u^{(3)})^2$ is the only quadratic monomial in $(b_2,b_3,b_4)$, no $b_5$, $\chi_4=0$ for a ray.
\end{enumerate}
''' % dict(K0=g(c5['K0bar']['up'], 'K0bar'), K2pt=g(c5['K2pt']['up'], 'K2pt'), marg=g(c5['margin_combined']['dn'], 'margin_combined', 3, 'dn'), Cpp=g(c5['Cpp']['up'], 'Cpp'),
           Spp=g(c5['Spp']['up'], 'Spp'), Cp=g(c5['Cp']['up'], 'Cp'), Sp=g(c5['Sp']['up'], 'Sp'), s=g(c5['slack_s']['dn'], 'slack_s', 4, 'dn'), W57=g(c2['W57']['up'], 'W57'))
out = out.replace('\\paragraph{Status.}', impl + '\n\\paragraph{Status.}')
open(os.path.join(TX, 'flux2_lemma_v2.tex'), 'w').write(out)

rows = []
for k in sorted(cond):
    v = cond[k]
    m = v.get('margin')
    rows.append('%s & %s & %s\\\\' % (k.replace('_', '\\_'), {True: 'pass', False: 'FAIL', None: 'open'}[v['pass_']], ('$\\ge%s$' % g(m, 'margin_' + k, 3, 'dn')) if isinstance(m, float) else '--'))
absorb = r'''%% absorption_lemma_v2.tex -- generated by fix8_flux2_impl/code/gen_texts.py from out/flux2_cert.json. Not compiled.
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}\usepackage{amsmath,amssymb,amsthm,booktabs,enumitem}
\newtheorem{lemma}{Lemma}\newtheorem{proposition}{Proposition}
\newcommand{\fn}[1]{\texttt{\detokenize{#1}}}
\begin{document}
\section*{Absorption of the time-integrated outer-ramp flux (flux lemma v2)}
Setting as in \fn{flux2_lemma_v2.tex}, \S\ref{F:sec:abs}. With $\mu=0.08$, $\varepsilon_I=3\cdot10^{-5}$, $F_0=\mu(\fn{FLUX0}+\varepsilon_I\fn{FLUX2})=%(F0)s$, for a.e.\ $\tau$
\[\tfrac12\Lambda'\le\mathrm{tot}(\alpha(\tau),x(\tau))-s+\varphi(\tau),\qquad\int_\sigma^\tau\varphi\le C'+S'(\tau-\sigma),\]
with $s\ge%(s)s$, $C'\le%(Cp)s$, $S'\le%(Sp)s$ [cert].
\begin{lemma}[reflected integral] If $S'\le s$, then $Z(\tau):=2\sup_{\sigma\le\tau}\int_\sigma^\tau(\varphi-s)$ satisfies $0\le Z\le2C'$, $Z(\tau_1)=0$, $Z$ is absolutely continuous with $Z'\ge2(\varphi-s)$ a.e., and $Y:=\Lambda-Z$ satisfies $\frac12Y'\le\mathrm{tot}$ a.e.\ and $Y\le\Lambda\le Y+2C'$.\end{lemma}
\begin{proof} As in \fn{flux2_lemma_v2.tex}, Lemma~\ref{F:lem:Z}.\end{proof}
\begin{proposition} Conditions (C0)--(C7) of \fn{flux2_lemma_v2.tex} hold for both production points of \fn{fin13v4b} [cert, \fn{flux2_assemble.py}]:
\begin{center}\begin{tabular}{@{}lll@{}}\toprule condition & status & margin\\\midrule
%(rows)s
\bottomrule\end{tabular}\end{center}
(C8) the regularity lemma U-1 is open; (C9) all constants are at the (B$\Phi$) level. Hence the verifier's conclusions hold with $\Lambda\le L(\tau)+2C'$ (march) and $\Lambda\le\ell(\tau)+2C'$ (tail), conditionally on U-1.\end{proposition}
\end{document}
''' % dict(F0=g(c5['F0']['up'], 'F0'), s=g(c5['slack_s']['dn'], 'slack_s', 4, 'dn'), Cp=g(c5['Cp']['up'], 'Cp'), Sp=g(c5['Sp']['up'], 'Sp'), rows='\n'.join(rows))
open(os.path.join(TX, 'absorption_lemma_v2.tex'), 'w').write(absorb)

errat = r'''%% erratum_flux2.tex -- erratum entries for the outer-ramp flux (generated; not compiled).
\begin{itemize}
\item[(E-F1)] \emph{Cutoff shape.} The v3.1 flux tables used a quintic cutoff; the verifier's cutoff is the septic $S_7((60-x)/3)$, $C^3$, with $\zeta^{(4)}(57^+)=-280/27$. All v2 tables are generated from $S_7$ (\fn{flux2_coef.py}).
\item[(E-F2)] \emph{Jet box.} v3.1 modelled the commutator with quasilinear coefficients and an a-priori jet box (M2--M5). v2 uses the commutator of $L_1$ with the quasilinear factor $j=J^{-4}$ inside it, $j\in[1-\bar c,1+\bar c]$, $|j'|\le\bar c_1$, and no jet box.
\item[(E-F3)] \emph{Drift double use.} The drift sign ($C_{\rm eff}>0$, $\zeta\zeta'\le0$) is used once, in the one-sided diagonal, with the modulation transport inside the same ball $C_{\rm eff}\in C\pm\bar a(|Q_u|+|C_u|)$.
\item[(E-F4)] \emph{Trace route.} The endpoint route through $g^{(3)}(57)$ does not close $K_2^{\rm OR}$ (residue $\frac52 g''(g^{(3)})^2$, \fn{flux2_forms.py} (f)); v2 uses the level-2 export of (A) as a space-time $L^2$ bound of $g^{(4)}$ on $[57,60]$, with no integration by parts at level 2.
\item[(E-F5)] \emph{$W(57)$.} The value $1.43716\cdot10^{-18}$ of \fn{lemmaS_cflux.json} belongs to the quintic weight blend; for the verifier's $S_7$ blend $W(57)\le%(W57)s$. The (A)-chain flux block (records only) still uses the old value; flux lemma v2 does not.
\end{itemize}
''' % dict(W57=g(c2['W57']['up'], 'W57'))
open(os.path.join(TX, 'erratum_flux2.tex'), 'w').write(errat)
bphi = r'''%% BPhi_continuation_v36.tex -- paragraph for the continuation text (v36 suffix; generated; not compiled).
\paragraph{Bootstrap hypothesis (B$\Phi$).} For $k=1,2$ and every window $I=[\tau_a,\tau_a+T]\subset[\tau_1,\tau']$:
$\int_I\int_{\ln28}^{\ln350}y^{-4}(28/y)(z^A_{k+2})^2\,ds\,d\tau\le2(C0_k+C1_kT)$, with $C0_1=%(C01)s$, $C1_1=%(C11)s$, $C0_2=%(C02)s$, $C1_2=%(C12)s$ (Lemma~S(A), v36, setting $4N_2^{\rm ref}$).
Under the bootstrap set including (B$\Phi$), flux lemma v2 holds with $C'\le%(Cp)s$, $S'\le%(Sp)s$; the absorption proposition gives the verifier's conclusions with $\Lambda\le L+2C'$; Lemma~S(A) then gives the window bounds with $(C0_k,C1_k)$, a strict improvement. The window integrals are continuous in $\tau'$, so (B$\Phi$) is closed and open. Lemma~S(A) itself does not use (B$\Phi$).
''' % dict(C01=g(c4['export']['1']['C0']['up'], 'C0_1'), C11=g(c4['export']['1']['C1']['up'], 'C1_1'), C02=g(c4['export']['2']['C0']['up'], 'C0_2'),
           C12=g(c4['export']['2']['C1']['up'], 'C1_2'), Cp=g(c5['Cp']['up'], 'Cp'), Sp=g(c5['Sp']['up'], 'Sp'))
open(os.path.join(TX, 'BPhi_continuation_v36.tex'), 'w').write(bphi)
lemD = r"""%% lemmaD_k2_sketch.tex -- statement and proof SKETCH of the level-2 window export of Lemma S(A) (unit G1). [hand]; generated; not compiled.
%% Status: the export is CERTIFIED as a computation in v36 (assembly_rep.export, A_interface r2_export); this written proof is a sketch, NOT reviewed.
\begin{lemma}[Lemma D, $k=2$]\label{D:k2}
Under the bootstrap hypotheses of Lemma~S(A) on $[\tau_1,\tau']$, for every $\tau_a\ge\tau_1$ and $T\ge0$ with $\tau_a+T\le\tau'$,
\[X_2(\tau_a,T):=\int_{\tau_a}^{\tau_a+T}\!\int_{\ln28}^{\ln350}y^{-4}\tfrac{28}{y}\,(z^A_4)^2\,ds\,d\tau\le C0_2+C1_2\,T,\]
\[C0_2=\frac{\hat E_2+D_2}{r_2a_c},\qquad C1_2=\frac{2F_2^2/\gamma'_2+\mathrm{coup}_2+o_2+D_2^{\rm slope}}{r_2a_c},\qquad r_2=\tfrac1{1000},\ a_c=(1+0.11^2)^{-2},\]
with, at the setting $4N_2^{\rm ref}$ used by flux lemma v2, $C0_2\le%(C02)s$, $C1_2\le%(C12)s$ [cert, v36].
\end{lemma}
\begin{proof}[Sketch \textup{[hand]}]
(1) \emph{Retained dissipation.} In the level-2 energy inequality of Lemma~S(A) the dissipation $2\int\Theta(z^A_4)^2$ is split by the weight tables of \fn{weights_t2}: the parts $B_c$ used to absorb the level-2 sources and a remainder. With $B_c^{\rm core}=1688999/10^6$ and, on every zone of $[28,350]$, $B_c\le2-\mathrm{used}/a-r_2a_c/a_z$, the remainder is $\ge r_2a_c\,y^{-4}\frac{28}{y}(z^A_4)^2$ pointwise on $[\ln28,\ln350]$ (asserted in arb in \fn{weights_t2}). So
\[\tfrac{d}{d\tau}E_2+r_2a_c\int_{\ln28}^{\ln350}y^{-4}\tfrac{28}{y}(z^A_4)^2\,ds\le\textstyle\sum\mathrm{sources}.\]
(2) \emph{Sources are window-affine.} Every level-2 source of \fn{assembly_rep} is a \texttt{Tr} object (inflow, ramp commutators, IBP top term with its Young users, the class costs, foot, S33): over any window its integral is at most a window constant plus a slope times $T$, the damped constant being $\le D_2$ (asserted $D_2^c\le D_2$) and the slopes summing to $D_2^{\rm slope}$. The forcing is bounded by Young against the retained part of $\gamma'_2E_2$, giving $2F_2^2/\gamma'_2$ per unit $\tau$; the couplings to $k<2$ use the pointwise $E_k\le\hat E_k$ and give $\mathrm{coup}_2$ per unit $\tau$; the outer part gives $o_2$ per unit $\tau$.
(3) \emph{Integration.} Integrate (1) over $[\tau_a,\tau_a+T]$, drop $E_2(\tau_a+T)\ge0$ and use $E_2(\tau_a)\le\hat E_2$ (pointwise R-sup bound, every $\tau_a\ge\tau_1$).
\end{proof}
\noindent\emph{Open:} the written verification, source by source, that each level-2 \texttt{Tr} bound holds in window form (the code asserts it by construction; no text exists), and an independent review. Until then the export, and with it the integrated part of flux lemma v2, rests on the certified computation plus this sketch.
""" % dict(C02=g(c4['export']['2']['C0']['up'], 'C0_2'), C12=g(c4['export']['2']['C1']['up'], 'C1_2'))
open(os.path.join(TX, 'lemmaD_k2_sketch.tex'), 'w').write(lemD)
json.dump(NUM, open(os.path.join(TX, 'numbers_flux2.json'), 'w'), indent=1, sort_keys=True)


def balance(path):
    t = open(path).read()
    t2 = re.sub(r'(?<!\\)%.*', '', t)
    br = t2.count('{') - t2.count('\\{') - (t2.count('}') - t2.count('\\}'))
    b = re.findall(r'\\begin\{([^}]*)\}', t2); e = re.findall(r'\\end\{([^}]*)\}', t2)
    return dict(brace_balance=br, begin_end_equal=sorted(b) == sorted(e), ok=(br == 0 and sorted(b) == sorted(e)))


files = ['flux2_lemma_v2.tex', 'absorption_lemma_v2.tex', 'erratum_flux2.tex', 'BPhi_continuation_v36.tex', 'lemmaD_k2_sketch.tex', 'numbers_flux2.json']
bal = {f: balance(os.path.join(TX, f)) for f in files if f.endswith('.tex')}
res = dict(stage='X7 gen_texts', balance=bal, all_balanced=all(v['ok'] for v in bal.values()), outputs={f: FC.sha(os.path.join(TX, f)) for f in files},
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), os.path.join(FC.DESIGN, 'flux2_lemma.tex')] +
                              [os.path.join(FC.OUT, n + '.json') for n in J]))
res['timing'] = dict(runtime_s=T_.s())
h = FC.dump(os.path.join(FC.OUT, 'gen_texts.json'), res)
print('balance', bal); print('X7 done sha', h)
assert res['all_balanced']
