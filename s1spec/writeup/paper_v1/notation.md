# Notation table (paper v1) -- unit U1

Final symbols for the main paper and the supplement. "Where defined" uses the main-paper sections of OUTLINE.md (§1 intro, §2 objects, §3 rescaled formulation, §4 step 1, §5 step 2, §6 step 3, §7 step 4). "Replaces / clash removed" names the archive symbol. Parts B--K of the supplement use the same symbols; the archive symbols appear only in Part K. Macro names are in `preamble.tex`.

Rules: (i) a symbol has one meaning in the whole paper and supplement; (ii) a subscript or accent that is part of a symbol in this table is not optional; (iii) local auxiliary symbols may be reused inside one proof, if introduced there with "write ... for" and not used outside it; (iv) never use a renamed archive symbol in its old meaning.

| symbol | meaning | where defined | replaces / clash removed |
|---|---|---|---|
| **Flow and curves** | | | |
| $\gamma(t)$ | closed immersed plane curve evolving by curve diffusion flow (CDF) | §1 | — |
| $\mathfrak s$ | arclength of the flow curve; CDF is $\partial_t\gamma=-\kappa_{\mathfrak{ss}}\nu$ (macros \arclen, \kss) | §1 | $s$ in the flow equation ($s$ is now the log-radius, below) |
| $\kappa,\ k,\ \kappa_\Gamma$ | curvature of $\gamma$; curvature of the profile curve of $S_1$; curvature of the rescaled curve | §1, §2, §3 | keep (OUTLINE 5.2) |
| $T,\ \nu,\ \theta$ | unit tangent $T=e^{i\theta}$, normal $\nu=iT$, tangent angle | §1 | $\tau$ used for the tangent (step 3, v3.2: $q=\langle\Gamma,\tau\rangle$) |
| $\operatorname{Area},\ \mathcal L,\ \operatorname{rot}$ | signed enclosed area; length of a curve; turning number | §1 | $A$; $L(\gamma)$ (the symbol $L$ is the linearised operator); $\omega$ (turning number) |
| $T_*,\ T_{\max}$ | singular time; maximal time of existence | §1, §6.7 | $T$ (clashes with the tangent $T$) |
| $p_*$ | image point of the singularity (D13) | §1 | $x_T$ (step 1), $p(t)$ |
| $\Gamma_0,\ \vartheta\in[0,2\pi]$ | initial curve; its parameter (proportional to arclength) | §2, §7 | $\alpha$ in step 4 (clash with the slow amplitude $\alpha$) |
| $u_n,\ b_n,\ \varpi_n=n^4-n^2$ | Fourier coefficients of the tangent angle $\theta(\vartheta)=-\frac\pi2+\vartheta+\sum_nu_n\sin n\vartheta$ (evolving), their value $b_n$ for $\Gamma_0$; linear damping rates | §7 | $\lambda_n$ |
| $\delta_t$ | clock defect: $t_\tau=\lambda^4(1-\delta_t)$ in step 4 (macro \dclk) | §7 | $\vartheta$ in step 4 |
| **Rescaling, modulation** | | | |
| $\lambda(t),\ \lambda_1=\lambda(t_1)$ | scale of the singular profile; scale at the hand-off | §3 | spectral parameter (now $\mu$); slow eigenvalue $\tilde\lambda$ (now $\mu_3$) |
| $\tau,\ \tau_1,\ \tau_{\max}$ | rescaled time $d\tau=\lambda^{-4}dt$; value at $t_1$; end of the rescaled interval of existence | §3, §6.7 | $\tau_{\rm tip}$, $\tau=S-\sigma$ of step 1 (now $\sigma^{\rm rev}$); tangent $\tau$ |
| $c(t)e_1,\ \Gamma(\tau)$ | centre of the frame, $\gamma=c\,e_1+\lambda\Gamma$; rescaled curve $\Gamma(\tau)=(\gamma-ce_1)/\lambda$ | §3 | $p(t)$ (step 3 v3.2; $p$ is the support function); $c(\tau)=\frac14-a$ (now $\varsigma$) |
| $a,\ b$ | scale and translation modulation: $a=\frac14+\lambda_\tau/\lambda$, $b=-c_\tau/\lambda$; bounds $\bar a,\bar b$ | §3, §6 | keep |
| $\varsigma:=\tfrac14-a$ | scale rate, $\lambda_\tau/\lambda=-\varsigma$ (macro \rate); $\varsigma^\pm=\frac14\pm\bar a$ | §3, §6.4 | $c(\tau)$, $c^\pm$ (Lemma U), $c:=\frac14-\tilde a$ (step 4) |
| $\xi_a:=(1-4\bar a)^{-1}$ | constant of Lemma U | §6.4 | $q$ (Lemma U; $q$ is $\langle\Gamma,T\rangle$) |
| $\alpha,\ \beta$ | slow-mode amplitude; its forcing, $\alpha_\tau=\mu_3\alpha+\beta$; bound $\bar\beta$ | §3 | $\beta$ as weight exponent (step 2, now literal $\frac14$), $\beta_s$ (now $b_s=b\cos\theta\,e^{-s}$), $\beta_A$ keeps its subscript |
| $e_1,\ (X_1,X_2)$ | unit vector of the symmetry axis of the clip; ambient coordinates | §2 | keep |
| **The shrinker $S_1$ and step 1** | | | |
| $S_1,\ \Gamma_{\!*},\ \mathcal S$ | the shrinker; its arclength-parametrised profile curve; the certified set of profiles (macros \Sone, \Gsh, \cSet) | §2 | $\Gamma$ for the profile curve in steps 1 and 2 ($\Gamma$ is the rescaled curve) |
| $\sigma,\ \sigma_{\rm tip}$ | arclength of $\Gamma_*$; position of the tip, $\sigma_{\rm tip}\in[-2.249127,-2.249019]$ | §2 | keep |
| $x,\ y,\ X_0$ | arclength from the tip, $x=\sigma-\sigma_{\rm tip}$; $y=x-X_0$, $X_0=2.2490776515$; far field: $y=r/\lambda$ | §2, §3 | keep (OUTLINE 5.2) |
| $p,\ q,\ k,\ m$ | on $\Gamma_*$: $p=\langle\Gamma_*,\nu\rangle$, $q=\langle\Gamma_*,T\rangle$, $k=\theta'$, $m=k'$ | §2 | $p=1.0246$ (now $p_U$); $p(t)$; $q=(1-4\bar a)^{-1}$ (now $\xi_a$); slab index $m$ (now $i$) |
| $D,\ k_0$ | tip distance $D=2.500806277$; tip curvature $k_0=1.539820638$ | §2 | $D$ as Evans function, as $y\partial_y$, as a weight, as a derivation step |
| $\varepsilon,\ M(\varepsilon)$ | shooting parameter and shooting function of step 1 (unsubscripted $\varepsilon$ occurs only in §4 and Part B) | §4 | $\varepsilon=10^{-5}$ (now $\varepsilon_E$); $\varepsilon=y^{-4}$ (written out); chart quantity $\varepsilon_c$ (now $\epsilon_{\rm ch}$) |
| $\mathcal Y$ | certified parallelepiped of step 1 | §4 | $Y$ (step 1; $\bar Y$ is the step-4 majorant) |
| $\Sigma$ | symmetry axis of $S_1$ (reflection line) | §2 | $\ell$ (the symbol $\ell$ is the cylinder unit) |
| $\sigma_{\rm c},\ \sigma^{\rm rev}$ | start of the conical-end construction ($=8$); reversed variable $\sigma_{\rm c}-\sigma$ of the shooting | §4 | $S=8$, $\tau=S-\sigma$, $\tau_1=8.5$ of step 1 |
| **Linear theory (step 2 and the linearisation)** | | | |
| $L$ | linearised operator at $S_1$: $Lu=-u''''-k^2u''-(5km+\frac q4)u'+(\frac14-3m^2-kp)u$ | §3, §5 | $\mathcal L$ (step 4 base operator, now $\widetilde L$) |
| $\rho,\ \|\cdot\|_\rho,\ \mathcal H_\rho$ | weight $\rho=\exp(-\frac14(1+x^2)^{2/3})$; norm; Hilbert space $L^2(\rho\,dx)$ | §3, §5 | $H$ for the Hilbert space (step 2); $\rho$ as radius (Lemma U, now $r_c$); $\rho_\theta$ (step 4, now $r_\theta$) |
| $\mu,\ \mu_3$ | spectral parameter; slow even eigenvalue $\mu_3=-0.2562$ | §5, §3 | $\lambda$; $\tilde\lambda$ |
| $\hat\psi_3,\ \hat h_2,\ w,\ \tilde w$$ | slow eigenfunction, second-order correction, remainder: $u=\alpha\hat\psi_3+\alpha^2\hat h_2+w$; $\tilde w=\chi w$ ($\chi$ cutoff) | §3 | keep |
| $\mathfrak r_{\rm re},\ \mathfrak r_{\rm im}$ | constants $1.790819$, $1.768282$ of the trapezoid of Theorem B | §5 | $\Lambda_0,\Lambda_1$ (step 2; $\Lambda$ is the energy) |
| $\mathrm{Ev}_{\rm even},\ \mathrm{Ev}_{\rm odd}$ | Evans functions | §5 | $D_{\rm even}$, $D_{\rm odd}$ |
| $\mathbf g$ | translation vector in $\mathbb R^2$ ($\mathbf g_\parallel,\mathbf g_\perp$) | §5 | $a\in\mathbb R^2$ (step 2; $a$ is the modulation) |
| $\mathfrak f_j,\ j=1,2,3$ | constraint functions of the modulation | §3 | $f_j$ ($f$ is the slope) |
| $\varphi_\rho$ | exponent, $\rho=e^{-\varphi_\rho}$, $\varphi_\rho=\frac14(1+x^2)^{2/3}$ | §5 | $\varphi$ (also the Lemma U weight factor, which keeps $\phi$) |
| **Energies and the tip** | | | |
| $E,\ E_I,\ \Lambda,\ \mu_I$ | inner energy, intermediate energy, $\Lambda=E+\mu_IE_I$, $\mu_I=0.08$ | §3, §6.2 | $\mu=0.08$ (clash with spectral $\mu$); $\mu_e=0.036$ of step 4 (now $\mu_{\rm err}$) |
| $\varepsilon_E,\ \varepsilon_I$ | $10^{-5}$ and $3\cdot10^{-5}$: weights of the second-derivative terms | §3, §7 | $\varepsilon$ (clash with the shooting parameter) |
| $\mathsf w_j=\|\tilde w^{(j)}\|_\rho$ | shape variables $j=0,\dots,4$ (macro \wn{j}) | §3 | $A,W_1,B,W_3,D_4$ (clash with bootstrap $B$, tip $D$, weight/window $W$) |
| $v_I=\zeta w'$ | intermediate-zone variable on $x\in[10.25,60]$, cutoff $\zeta$ | §3 | $v$ (step 3 and step 4; $v$ is $f-U_1$) |
| $J_j,\ \Lambda_{\rm bar}$ | time integrals of $V_j^2$; barrier level of $\Lambda$ | §6.2 | keep |
| **Far field (Lemma S, Lemma U, zone III)** | | | |
| $s=\ln y,\ \partial_s$ | log-radius; $\partial_s=y\partial_y$ | §3 | $D_y=y\partial_y$ (clash with $D$) |
| $f$ | slope of the curve as a graph over the $S_1$ ray | §3 | keep |
| $p_U,\ U,\ U_1$ | $p_U=1.0246$; $U(\alpha,y)=1.3056\,\alpha\,y^{p_U}$; frozen profile $U_1(r)=1.3056\,\alpha(\tau_1)(r/\lambda_1)^{p_U}$ | §3 | $p=1.0246$ |
| $v,\ v_A,\ z_k$ | $v=f-U_1$ (Lemma U), $v_A=f-U$ (Lemma S), $z_k=\partial_s^kv$ | §3, §6 | keep |
| $d,\ \delta U$ | profile difference $d=v_A-w'$; mismatch $\delta U=v_A-v$ | §6.3, §6.4 | keep |
| $\omega=28/y,\ \operatorname{sm}(t)$ | shell weight; smoothstep $t^3(10-15t+6t^2)$ | §6.3 | $\omega$ (turning number, now $\operatorname{rot}$); $S(t)$ (clash with $S_1$, Lemma S) |
| $\eta_{\rm in},\ \eta_o,\ \eta_A,\ \Theta_A$ | ramp cutoffs, $\eta_A=\eta_{\rm in}\eta_o$, $\Theta_A=y^{-4}\omega\eta_A^2$ | §6.3 | $\Theta$ (Lemma S and chart lemma use $\Theta$ for different things; Lemma U $\Theta:=WD$ is written out) |
| $E^A_j,\ \hat E,\ \mathfrak E^{r_c}_j,\ \hat B,\ \hat B^{\rm fine}$ | shell energies; their certified bounds; uniformly local energies centred at $r_c$; their certified bounds | §6.3, §6.4 | $\mathfrak E^\rho_j$ (centre $\rho$ now $r_c$) |
| $r,\ r_c,\ r_\chi$ | physical radius; centre of a local energy; outer centre limit $r_\chi=253\lambda_1$ | §6.4 | $\rho$ for the centre |
| $\varrho,\ \varrho_c$ | tail rate; coercivity rate of the step-4 comparison system ($\min\varrho_c=0.04711$) | §6.4, §7.4 | $\rho$ as rate |
| $I_w=[y_w,1.5y_w]$ | hand-off window | §6, §7 | $W$ for a window (clash with the weight $W$, which keeps $W$) |
| $\mathfrak q_1(V'),\ \mathfrak q_2(V')$ | quasilinear coefficients of the slope equation | §6.4 | $\mathfrak a,\mathfrak b$ (clash with $a,b$) |
| $b_s=b\cos\theta\,e^{-s}$ | translation drift term of the slope equation | §6.4 | $\beta_s$ |
| $g^A_j$ | level-$j$ coefficients of the (A) energy inequality (certified; Part E) | §6.3 | $\gamma'_j$ (clash with the curve $\gamma$) |
| $K_4$ | fourth-derivative time-average bound on the cone cylinders at the hand-off (Lemma O data) | §6.5, §7 | keep; the $K$ of Lemma O is $K_O$ |
| **Cylinders (Lemma O, Lemma R)** | | | |
| $\ell$ | cylinder length unit; cone family $\ell=50\lambda_1$, far family $\ell=113\lambda_1$ | §6.5 | $L$ in \texttt{fixed-L}, $L=50\lambda_1$, $L=113\lambda_1$ (\texttt{fixed-}$\ell$ in the paper) |
| $X,\ h_0,\ t^*$ | half-length in cylinder units, graph function at the start time, certified time-ratio threshold | §6.5 | keep |
| $M_O,\ K_O$ | Lemma O slope/jet constants, $(M_O,K_O)=(0.12,0.18)$ cone, $(0.15,0.3)$ far | §6.5 | $(M,K)$ ($M$ is the shooting function) |
| $H_{\rm cyl},\ C^O_j,\ \Lambda^O,\ \Xi_O,\ f_O$ | height of the cylinder box; Lemma O constants $C_*,C_1,C_2,C_3,\Lambda_1$; $\kappa:=K_3+\min\{\sqrt{2X}M,K_3\}$; slack fraction $f_O=\frac12$ | §6.5 | $H$; $C_j,\Lambda_1$; $\kappa$ (Lemma O v3.2 §1(a)); $f$ |
| $M_0,\ M_1$ | core boxes of Lemma S: $\sup\|f\|$, $\sup y\|f_y\|$ | §6.3 | keep |
| **Step 4 (Phase I, Phase II, hand-off)** | | | |
| $\bar Y,\ Y$ | error majorant and its certified constant; $Y=\sqrt{\Lambda_e}$ | §7 | keep ($\mathcal Y$ is step 1) |
| $\Lambda_e,\ \mu_{\rm err}$ | error energy $\Lambda_e=A_e^2+\varepsilon_EB^2+\mu_{\rm err}(V_0^2+\varepsilon_IV_2^2)$, $\mu_{\rm err}=0.036$ | §7.4 | $\mu_e$ |
| $\widetilde L$ | linearisation about the corrected base | §7.4 | $\mathcal L$ (step 4) |
| $\Pi,\ \mathrm{res},\ \mathrm{Res}_1,\mathrm{Res}_2$ | transport integral $\int_0^x\tilde\kappa\tilde V$; residual of the base; its first and second integrals | §7.4 | $\varpi$; $r$; $S_1,S_2$ (clash with the shrinker) |
| $\tilde\nu_1,\ \tilde T_1$ | $\langle e_1,\tilde\nu\rangle$ and $\langle e_1,\tilde T\rangle$ | §7.4 | $\tilde n$, $\tilde t_1$ ($n$ is the Fourier index) |
| $\Box_{\rm tip},\ \Box_{\rm far}$ | validity boxes of the error in the tip and far regions | §7.4 | $\mathbb B,\mathbb B'$ (clash with the bootstrap $B$) |
| $\mathcal E_{\rm env}(x)$ | far envelope of the error | §7.4 | $E(x)$ |
| $\beta_{\rm base}$ | $\beta:=Q-cX_0-6\tilde\kappa\tilde\kappa'+A$ in the coefficient $b_1$ of $\widetilde L$ | §7.4 | $\beta$ (clash with modulation $\beta$) |
| $i$ | slab index, slab $i$ covers $\tau\in[0.15i,0.15(i+1)]$ | §7 | $m$ (clash with $m=k'$) |
| $\tau^{\rm I}$ | Phase I clock, $d\tau^{\rm I}=(2\pi/\mathcal L)^4dt$ | §7.2 | $\tau'$ |
| $\text{Phase I, Phase II}$ | the two phases of step 4 | §7 | P1, P2 |
| **Names of conditions** | | | |
| $\mathrm{(H1)}$--$\mathrm{(H7)}$ | hand-off hypotheses of Theorem C | §2 | see namemap.md §1 |
| $\mathrm{(O1)}$--$\mathrm{(O9)}$ | hypotheses of Lemma O | §6.5 | see namemap.md §2 |
| $\mathrm{(L1)}$--$\mathrm{(L3)}$ | literature facts | §2, §6.7 | (P1$'$)--(P3$'$) |
| $\mathrm{(B1)}$--$\mathrm{(B5)}$ | bootstrap conditions | §6.1 | B1, B2$'$, B3U, B4, G |
| $\mathrm{(K1)}$--$\mathrm{(K6)}$ | closure conditions | §6.6 | keep |
| $\mathrm{(W1)}$--$\mathrm{(W5)}$ | conclusions of Theorem C | §2 | keep |

## Clashes resolved beyond OUTLINE §5.2 (found by skimming the sources)

Each line: the clash, then the resolution (new symbol in bold). All are mechanical renamings; none changes a statement.

1. **Arclength $s$ vs log-radius $s=\ln y$.** The flow equation uses $s$ (arclength) in every source, and Lemma S/U use $s=\ln y$ (and §5.2 fixes $\partial_s=y\partial_y$). Resolution: log-radius keeps $s$ (as §5.2 requires); arclength of the flow curve is $\mathfrak s$ ($\kappa_{\mathfrak{ss}}$, macro `\kss`). See open question Q1.
2. **$c$:** frame point $c(t)e_1$ vs rate $c(\tau)=\frac14-a$ (Lemma U, with $c^\pm$) vs $c:=\frac14-\tilde a$ (step 4) vs constants $c_I,c_e,c_g,c_{\delta U}$. Resolution: rate is $\varsigma$ (so $\varsigma^\pm$); $c(t)e_1$ keeps $c$; subscripted constants stay.
3. **$q$:** $q=\langle\Gamma,T\rangle$ (steps 1--3) vs $q=(1-4\bar a)^{-1}$ (Lemma U). Resolution: $\xi_a$.
4. **$\kappa:=5$** in Lemma S(A) ($(\varepsilon\omega)'=-\kappa\varepsilon\omega$) and **$\kappa:=K_3+\min\{\sqrt{2X}M,K_3\}$** in Lemma O v3.2 vs curvature. Resolution: write the number $5$; Lemma O's bound is $\Xi_O$.
5. **$\varepsilon$:** shooting parameter (step 1: $\varepsilon,\varepsilon_a,\varepsilon_b,\varepsilon^*$) vs energy weight $10^{-5}$ (steps 3, 4) vs $\varepsilon=y^{-4}=e^{-4s}$ (Lemma S/U) vs chart quantity $\varepsilon_c$. Resolution: shooting keeps $\varepsilon$ (only §4/Part B); energy weight is $\varepsilon_E$ (and $\varepsilon_I$); $y^{-4}$ is written out; chart quantity $\epsilon_{\rm ch}$.
6. **$\omega$:** turning number (main_v3 §1) vs shell weight $\omega=28/y$. Resolution: $\operatorname{rot}$ for the turning number.
7. **$\vartheta$:** §5.2 gives $\vartheta$ to the parameter of $\Gamma_0$, but step 4 already uses $\vartheta$ for the clock defect ($t_\tau=\lambda^4(1-\vartheta)$, $|\vartheta|\le3.37\cdot10^{-8}$). Resolution: clock defect is $\delta_t$.
8. **$T$:** tangent vector vs singular time $T$ vs $T_{\max}$. Resolution: $T_*$ (singular), $T_{\max}$ (maximal); $T$ is the tangent.
9. **$S$:** shrinker $S_1$ / certified set $\mathcal S$ / Lemma S / smoothstep $S(t)$ / step-4 residual integrals $S_1,S_2$ / sup quantities $S^v_0,S^v_1,S^f_1$ and Lemma U's $S_1(2\hat B)$. Resolution: smoothstep is $\operatorname{sm}$; residual integrals are $\mathrm{Res}_1,\mathrm{Res}_2$; the sup quantities become $\operatorname{sp}^v_0,\operatorname{sp}^v_1,\operatorname{sp}^f_1,\operatorname{sp}_1(2\hat B)$ (macro `\spp`).
10. **$A$:** enclosed area vs $A=\|\tilde w\|_\rho$ vs $A_{\rm interface}$ (file) vs $A'=\tilde\kappa r$ (step 4). Resolution: area is $\operatorname{Area}(\gamma)$ (also because Lemma U's closure maps $\mathcal A(\beta),\mathcal U(\beta)$ keep their calligraphic letters); shape variables $\mathsf w_j$ (also removes $B,W_1,W_3,D_4$); the (A) of Lemma S(A) only as the lemma name.
11. **$B$:** bootstrap (B1)--(B5) vs $\hat B$ (Lemma U bound) vs $\mathfrak B_j$ vs $\mathbb B$, $\mathbb B'$ (step-4 boxes) vs $B=\|\tilde w''\|$. Resolution: $\Box_{\rm tip},\Box_{\rm far}$ for the step-4 boxes; $B$ no longer a norm; $\hat B,\mathfrak B_j$ keep (calligraphic/hat, certified numbers use them).
12. **$H$:** hypotheses (H$n$) vs Hilbert space $H$ (step 2) vs thresholds $H_{2,j}$, $H^f_{2,j}$ (Lemma S/U) vs height $H$ (Lemma O). Resolution: $\mathcal H_\rho$; thresholds are written $\mathsf H^{(2)}_j$, $\mathsf H^{(3)}_j$ (consistent with the renamed (H2), (H3)); $H_{\rm cyl}$.
13. **$G$:** far datum $G_j$ of (H2f) vs property "G" of the bootstrap vs $\hat G_{22}$ (weights). Resolution: the property G becomes (B5); $G_j$ keeps (now the datum of (H3)); $\hat G$ keeps, local to Parts E/F.
14. **$D$:** tip distance $D$ vs Evans functions $D_{\rm even}$ vs $D_y=y\partial_y$ vs the weight $D(s)$ of Lemma U vs derivation steps D1--D6 vs $D_4$. Resolution: $\mathrm{Ev}$; $\partial_s$; the Lemma U weight is written $D_U(s)$; D1--D6 only in Part H; $D_4\to\mathsf w_4$.
15. **$\Lambda$:** energy $\Lambda=E+\mu E_I$ vs $\Lambda_0,\Lambda_1$ (step 2 trapezoid) vs $\Lambda_1$ (Lemma O) vs $\Lambda_e$. Resolution: $\mathfrak r_{\rm re},\mathfrak r_{\rm im}$; $\Lambda^O$; $\Lambda_e$ keeps (subscript).
16. **$\lambda_n=n^4-n^2$** (step 4 Fourier rates) vs scale $\lambda$. Resolution: $\varpi_n$ (and the step-4 transport integral $\varpi$ becomes $\Pi$).
17. **$v$:** $\zeta w'$ (step 3 intermediate energy, step 4) vs $f-U_1$ (Lemma U) vs $v=\mathrm{th}_\alpha$ (Phase I). Resolution: $v_I$; $v=f-U_1$; Phase I writes $\mathrm{th}_\vartheta$, $\mathrm{th}_{\vartheta\vartheta}$ explicitly.
18. **$w$:** remainder $w$ vs $w=\mathrm{th}_{\alpha\alpha}$ (Phase I) vs weight/window $W$. Resolution: as in 17; hand-off windows are $I_w$; the weight $W$ keeps $W$.
19. **$\beta$:** modulation $\beta$ vs weight exponent $\beta=\frac14$ (step 2) vs $\beta_A=2$ (flux multiplier) vs $\beta_s$ vs $\beta:=Q-cX_0-\dots$ (step 4). Resolution: step 2 writes $\frac14$; $\beta_A$ keeps; $\beta_s\to b_s$; step 4 $\beta_{\rm base}$.
20. **$\mathfrak a,\mathfrak b$** (quasilinear coefficients of the slope equation) vs modulation $a,b$. Resolution: $\mathfrak q_1,\mathfrak q_2$.
21. **$a\in\mathbb R^2$** (step 2 translation vectors) vs modulation $a$. Resolution: $\mathbf g$.
22. **$\ell$:** symmetry line (steps 1, 2) vs cylinder unit vs $\tilde\ell$ (end of the step-4 domain) vs $\ell_j$ (functionals). Resolution: $\Sigma$ for the axis; $\ell$ cylinder; $\tilde\ell$ keeps; $\ell_j$ local.
23. **$L$:** operator $L$ vs length vs cylinder length $L$ ("fixed-$L$") vs Phase I length in $(2\pi/L)^4$. Resolution: §5.2 plus Lemma O's "fixed-$L$" becomes "fixed-$\ell$"; length is $\mathcal L$ everywhere (Phase I: $(2\pi/\mathcal L)^4$); the step-4 base operator $\mathcal L$ becomes $\widetilde L$.
24. **$M$:** shooting function vs $(M,K)$ of Lemma O vs $M_0,M_1$ vs $M_{1z}$ vs matrices $M_{ij}$ (modulation) . Resolution: §5.2, and $M_{ij}$ is the modulation matrix, local to Part D.
25. **$K$:** closure conditions (K1)--(K6) vs $K_4$ vs $(M,K)$ vs $K_2,K_3$ (Lemma O) vs $K_J$ ($\kappa_J$). Resolution: $K_O$; $K_4$ and (K$n$) keep.
26. **$\mu$:** energy weight $0.08$ vs spectral parameter vs $\mu_e=0.036$ (step 4) vs weight $\mu$ in AM--GM. Resolution: §5.2 plus $\mu_{\rm err}$.
27. **$\rho$:** weight $\rho=e^{-\varphi}$ vs centre $\rho$ (Lemma U) vs $\rho_\theta$ (step 4 residual density) vs $\varrho_c,\varrho_m$. Resolution: §5.2; $r_\theta$.
28. **$n$:** Fourier index vs $\tilde n=\langle e_1,\tilde\nu\rangle$ vs $n=\|\chi\tilde V\|_\rho$ vs $n^{\rm III}$ (zone-III data). Resolution: $\tilde\nu_1$, $\tilde T_1$; $\mathsf n_V:=\|\chi\tilde V\|_\rho$ (local, Part I); $n^{\rm III}$ keeps.
29. **$m$:** $m=k'$ vs slab index $m$ (step 4). Resolution: slab index $i$.
30. **$r$:** physical radius vs residual $r=\tilde V+\tilde\kappa_{xx}-\dots$ (step 4) vs $r_k$ (coefficients in Lemma O). Resolution: $\mathrm{res}$.
31. **$E$:** inner energy vs $E(x)$ (step-4 far envelope) vs $\mathfrak E$ vs $E^A_j$. Resolution: $\mathcal E_{\rm env}$.
32. **$\Gamma$:** rescaled curve $\Gamma(\tau)$ vs profile curve of $S_1$ ($\Gamma$ in steps 1, 2) vs $\Gamma_0$ vs $\Gamma_{\rm tip}$ (arc). Resolution: $\Gamma_{\!*}$ for the profile; $\Gamma_{\rm tip}$ keeps.
33. **$\tau$:** rescaled time vs tangent (v3.2) vs $\tau'$ (Phase I clock) vs $\tau_\rho,\tau_e,\tau_a$ (entry times). Resolution: Phase I clock $\tau^{\rm I}$; entry times keep (subscripts).
34. **$\sigma$ and $S=8$, $\tau=S-\sigma$ in step 1.** Resolution: $\sigma_{\rm c}$, $\sigma^{\rm rev}$.
35. **$\Theta$:** $\Theta=\varepsilon\omega\eta_A^2$ (Lemma S(A)) vs $\Theta=\varepsilon\omega\eta_o^2$ (Lemma tools) vs $\Theta:=WD$ (Lemma U) vs $\Theta$ = limit angle (chart lemma notes). Resolution: $\Theta_A$ for Lemma S(A) only; others local; the limit angle is $\theta_\infty$ (as the source already does).
36. **$f$:** slope vs constraint functions $f_j$ vs slack $f=\frac12$ (Lemma O) vs $f(u)$ (polynomial in Lemma O proof). Resolution: $\mathfrak f_j$; $f_O$; proof-local $f$ stays local.
37. **$\gamma'_j$** (coefficients of Lemma S(A) energy inequality, with $\gamma=1$ damping) vs the curve $\gamma$. Resolution: $g^A_j$; damping "$\gamma=1$" is dropped (the number $1$).
38. **$p_U$ and $p_*$:** $p_U=1.0246$; $p_*$ the singular point; $p=\langle\Gamma_*,\nu\rangle$. No further clash.
39. **Phase naming:** P1, P2 (step 4) vs (P1$'$)--(P3$'$) (literature) vs P1/P2 as files. Resolution: Phase I, Phase II; (L1)--(L3).
40. **(H$n$) numbering itself:** the archive's "H4" means three different things: (a) hand-off graph (v3.2/v3.3 (H4), now (H5)), (b) Lemma O H4 (fifth-derivative bound, in step-4 file names `h4loc`, `h4far`; now (O4)), (c) "H4 ×4.18 / ×5.68" in CLAUDE.md = (b). The paper uses only (H$n$) for hand-off and (O$n$) for Lemma O; see namemap.md §1--2.

## Open questions (mathematical or editorial judgement, not decided here)

- **Q1 (taste).** Following §5.2, log-radius is $s$, so the arclength of the flow curve is $\mathfrak s$ and the flow equation reads $\partial_t\gamma=-\kappa_{\mathfrak{ss}}\nu$. A reader of §1 may prefer plain $s$ there. The alternative (arclength $s$, log-radius renamed) contradicts §5.2 and touches every Lemma S/U formula. The authors should confirm; changing it later means changing only the two macros `\arclen`, `\kss` and the symbol `s` in Parts E, F.
- **Q2.** Lemma D and Lemma H-J-W occur twice in the archive (`lemmaS_A_v33.tex` `A:lem:D`, `A:lem:HJ` and `lemmaU_v34.tex` `lem:D`, `lem:HJ`), each time with "window form" or "sufficient check" in the statement. Whether the two copies are the same statement (then one lemma in the paper) or two different ones needs a statement-by-statement comparison; the namemap lists them as two entries until then.
- **Q3.** Lemma O's hypothesis H6 is no longer a hypothesis (it is the recovery proposition). Keep the numbering (O1)--(O9) with (O6) left as "proved, see Proposition", or renumber (O1)--(O8)? OUTLINE §5.3 says (O1)--(O9); I kept that and mark (O6) as unused.
- **Q4.** Names for the step-4 lemmas (Lemma 4a--4c, F, B, J, N3$'$, Corollary C) clash with step-3 names (Lemma B, J, N, R...). I propose Lemma E1--E8 in namemap.md §4, only as a proposal; U14/U15 own the final choice.
- **Q5.** Whether the step-3 symbol $\tilde w=\chi w$ (cutoff $\chi=\chi_{32}$) clashes with the chart terms $\chi_k$ of the chart lemma and the cutoff $\chi$ in the Lemma U weight $\Psi$. I left $\chi$ for cutoffs and renamed nothing; the chart terms could become $\mathsf c_k$. Needs a decision by the unit that writes Part E.
- **Q6.** The symbol $\theta$ also names a "$\theta$-form" of constrained coercivity (OUTLINE §6.2) and $\theta_0$ in the step-4 coercivity certificates. Not renamed; U7 and U14 should check that no formula mixes it with the tangent angle.
- **Q7.** $x$ is the arclength from the tip, but $x\in[0,\tilde\ell]$ is also the step-4 variable, and $x_c$ appears in "step-3 modulation $(x_c,\lambda,\alpha)$". Probably the same object ($x_c\leftrightarrow c$); I did not check. U14 should confirm and use $c$.
