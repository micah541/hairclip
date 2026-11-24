# T5 design: the shadowing theorem for P2

Scope: design for the paper proof and for T1, T6–T9. Read: STEP4_ARCH, STEP4_DESIGN §6–9, FEAS, FEAS2 (+ `feas2/p1.py`, `p2.py`), `step3_proof.tex`, `ASSEMBLY.md`, `verifier/constants.py`, Lemma O v3, Lemma S summary v3, M1 and T3 reports.
Tags: **[F]** float check done for this note (scripts in `t5/chk/`, 1 process); **[D]** derived here on paper, not machine-checked; **[S3]** taken from step 3.

## 0. Design decisions and the float facts behind them

**D1. All three non-coercive directions are exact symmetries; the third is the clock.**
The rescaled flow with clock τ is autonomous, so its symmetries are scaling (λ), axial translation (c) and the shift of τ along the base orbit. The mode of the τ-shift is the base's own normal velocity $\tilde V=\langle\partial_\tau\tilde\Gamma,\tilde\nu\rangle\approx\alpha_\tau\psi_3$, i.e. the slow direction. So the true solution is compared with $\tilde\Gamma(\tau)$ at a modulated label τ, and no residual is created in the third direction.
This replaces the scalar $a_3=\langle e,f_3\rangle$ of the brief. The two are the same object: $a_3=n(\tau)\,s$ with $s$ the accumulated clock shift and $n=\lVert\chi\tilde V\rVert_\rho$. The 2×2 comparison system survives in triangular form ($c_{21}=0$, §2.3), so the closure condition $c_{12}c_{21}<\mu_3c$ is empty. A scalar carried over fixed S1 functions would instead need $c_{12}c_{21}\sim(66.5\lvert\alpha\rvert)\cdot(10$–$60)\lvert\alpha\rvert<0.256\cdot0.05$, i.e. $k\gtrsim600$–1500 **[D]**.

**D2. Constraints are fixed step-3 functions; directions come from the base.**
Constraints: $\langle\chi e,\varphi_j\rangle_\rho=0$ with $\varphi_1=f_1=\chi p_{S1}$, $\varphi_2=f_2$, $\varphi_3=\chi\hat\psi_3$ (the direct slow eigenfunction, `c4fun.py`; not the adjoint $f_3$). They do not depend on τ and are arb-evaluable. Directions: $g=(\tilde p,\tilde n,\tilde V/n)$, each normalised in $L^2_\rho$ on $\operatorname{supp}\chi$ (the factors go into $x_i$).
What must be coercive is the *projected* form (§2.4), not $\mathrm{sym}_\rho L$ on a complement. Float values of its top eigenvalue (Hermite FEM, Σ = 40, h = 0.1, transport included; `chk/concheck.py`, `oblcheck.py`, `chk3.py`) **[F]**:

| constraints Φ / directions G | k=10 | 14 | 20 | 40 | 160 | 1810 | S1 |
|---|---|---|---|---|---|---|---|
| Φ = fixed $\{f_1,f_2,\chi\hat\psi_3\}$, G = base modes (**adopted**) | −0.061 | −0.072 | −0.079 | −0.083 | −0.083 | −0.082 | −0.082 |
| same, with $+\theta_0\lVert u''\rVert^2$, θ₀ = 0.05 | −0.038 | −0.050 | – | −0.061 | – | – | −0.061 |
| same, θ₀ = 1/9 (step-3 value) | −0.009 | −0.021 | – | −0.033 | – | – | −0.033 |
| Φ = G = base modes | −0.066 | −0.076 | −0.081 | −0.083 | −0.083 | −0.082 | −0.082 |
| Φ = 3 S1-sym eigenvectors (FEAS2 certificates), G = base modes | +0.063 | +0.071 | +0.083 | +0.106 | +0.132 | +0.141 | +0.142 |
| only $(a,b)$ modulated, no third direction | +0.469 | +0.325 | +0.178 | +0.038 | +0.044 | +0.064 | +0.067 |
| frozen-time adjoint as third constraint (sym form on complement) | +0.426 | +0.286 | +0.142 | −0.030 | −0.083 | −0.074 | −0.072 |

Consequences:
- The FEAS/FEAS2 statement "coercive on the complement of 3 S1-sym eigenvectors" (−0.17 to −0.24) does **not** give a contracting energy once the modulation term is included (row 5 is positive). The FEAS2 certificates certify the wrong form and must be redone for row 1.
- The usable coercivity is 0.06–0.083, not 0.17–0.24. Use θ₀ ≈ 0.02, not 1/9.
- $\lVert M^{-1}\rVert$ for the adopted pair is benign: cond ≤ 1.03 after normalisation, k ≥ 10 **[F]**.
- Frame rate of the base: $\lvert\tilde a\rvert\le3.4\cdot10^{-3}$ for k ≥ 7 **[F]** (`chk/atil.py`), so it does not move these numbers.

**D3. Contraction is not needed, only non-growth.** P1's ball (≈ 3e-14 physical) is far below the budget, so the error sits at the residual-driven level $Y\approx\Phi/\varrho$. A certificate interval with $\omega_c\le0$ is admissible (it costs $e^{\int\lvert\omega\rvert}$).

**D4. The intermediate and far energies of the error run over the whole far curve with pure power weights.** For the error there is no reason to cut at x = 60 (the residual is small everywhere). This removes the outer ramp, the flux term (c) and the $S_0^v,S_1^v$ inputs, and makes the chain tip → $E_I$ → far triangular. Float check of the pointwise level-0 rate with $W=(1+x^2)^{-2}$ on $[16.25,\ell]$: sup = −0.26 (k=10), −0.27 (14), −0.28 (40), −0.29 (640), against the asymptotic −3/8 **[F]**.

**D5. The accuracy driver is the far field at derivative levels 3–4, not the tip.** Items 3–4 need $Y=\sqrt{\Lambda_e}\lesssim2\cdot10^{-11}$ (§2.7), about 1000× below the tip budget of ARCH §3. This fixes the residual targets (§3) and forces an arb-polished base.

## 1. Theorem

### 1.1 Objects
- **Base (T1).** For $\tau\in[\tau_a-\delta,\tau_1]$: a half curve $\tilde\Gamma(x,\tau)$, $x\in[0,\tilde\ell(\tau)]$ = arclength from the tip, given by $\tilde\kappa$ (piecewise polynomial in $x$ and τ) and exact quadrature; both end points are symmetry points on the axis. Scalars $\tilde a(\tau),\tilde b(\tau)$. Derived: $\tilde T,\tilde\nu$, $\tilde p=\langle\tilde\Gamma,\tilde\nu\rangle$, $\tilde q=\langle\tilde\Gamma,\tilde T\rangle$, $\tilde n=\langle e_1,\tilde\nu\rangle$, $\tilde t_1=\langle e_1,\tilde T\rangle$, $\tilde V=\langle\partial_\tau\tilde\Gamma|_x,\tilde\nu\rangle$, transport $\varpi(x)=\int_0^x\tilde\kappa\tilde V$.
- **Residual.** $r:=\tilde V+\tilde\kappa_{xx}-(\tfrac14-\tilde a)\tilde p-\tilde b\,\tilde n$. Put $R_0:=r-\tilde\kappa_{xx}$ (no derivative of $\tilde\kappa$) and the primitives from the tip
$$S_1:=\int_0^xr=\tilde\kappa_x+\int_0^xR_0,\qquad S_2:=\int_0^x(x-x')r\,dx'=\tilde\kappa-\tilde\kappa(0)+\int_0^x(x-x')R_0\,dx'.$$
$S_2$ needs only $\tilde\kappa$; $S_1$ needs $\tilde\kappa_x$.
- **Linearisation** (FEAS F3 form with frame rates; signs to be re-derived symbolically and checked by finite differences as in `f3fd.py`):
$$\mathcal Le=-e''''+b_2e''+b_1e'+b_0e,\quad b_2=-\tilde\kappa^2,\quad b_1=-5\tilde\kappa\tilde\kappa'-(\tfrac14-\tilde a)\tilde q-\tilde b\tilde t_1+\varpi,\quad b_0=\tfrac14-\tilde a-3\tilde\kappa'^2-4\tilde\kappa\tilde\kappa''.$$
- **Weights, cutoffs [S3].** ρ, χ = χ₃₂, ε = 1e-5, ε_I = 3e-5 as in step 3. $\zeta_\infty$ = the step-3 inner ramp (0→1 on [10.25, 16.25]) and 1 beyond. $W_\infty=(1+x^2)^{-2}$ on all of $[10.25,\tilde\ell]$ (no blend). η = the Lemma S ramp on $y\in[14,28]$ and 1 beyond.
- **Energies of the error.** $u=\chi e$, $v=\zeta_\infty e'$; shape variables $A,W_1,B,W_3,D_4$ (ρ-norms of $u^{(j)}$) and $V_j=\lVert v^{(j)}\rVert_{W_\infty}$;
$$\Lambda_e=A^2+\varepsilon B^2+\mu_e\,(V_0^2+\varepsilon_IV_2^2),\qquad Y=\sqrt{\Lambda_e},\qquad \mu_e\in[0.01,0.08]\ \text{(set by T7, default 0.02)},$$
$$\mathfrak F_j=28\int_{14}^{\tilde\ell}\eta^2x^{2j-2}(\partial_x^je')^2dx,\qquad \mathfrak D_j=28\int\eta^2x^{2j-2}(\partial_x^{j+2}e')^2dx,\qquad j=0,\dots,4 .$$
Since $W_{\rm step3}\le W_\infty$ and $\zeta_{\rm step3}=\zeta_\infty\zeta_{\rm out}$, the step-3 functional satisfies $\Lambda^{S3}[e]\le(0.08/\mu_e)(1+10^{-12})\Lambda_e$. $\mathfrak F_j$ is the ω-normalised Lemma S energy with $x^j\partial^j$ in place of $(y\partial_y)^j$ (triangular equivalence).

### 1.2 Hypotheses (each is a certified output)
- **(I) Input (P1/T4).** At $t_a$ the true curve lies in a ball around P1's centre: ≤ 3e-14 in physical C², and ≤ 1e-9 in physical C⁵ away from the tip (expected from T3's weighted ℓ¹ bound and tail $Tn^{-6}$; T4 to confirm; it is used only for $\mathfrak F_j(\tau_0)$, where the tolerance is crude). The base is started from that centre: $\lVert\gamma_{P1}(t_a)-(\tilde c+\tilde\lambda\tilde\Gamma)(\tau_a)\rVert_{C^2}\le$ 1e-14.
- **(B) Base (T1).** (B1) exact symmetry, arclength parametrisation, $\tilde\kappa(\cdot,\tau)\in C^4$ piecewise polynomial, polynomial in τ on slabs and continuous across slabs. (B2) enclosures of the jets used below. (B3) $n(\tau)\ge n_{\min}>0$.
- **(R) Residuals (T6),** sup over each slab:
$\mathcal R_2=\lVert S_2\rVert_{L^2_\rho(0,34.25)}$, $\mathcal R_s=\lVert r\rVert_{L^2_\rho(0,34.25)}$, $\mathcal R_1=\lVert\zeta_\infty S_1\rVert_{W_\infty}$, $\mathcal R^I_s=\lVert\zeta_\infty r_x\rVert_{W_\infty}$, and $\mathcal R^F_j=\bigl(28\int\eta^2x^{2j-2}(\partial_x^{j-1}r)^2\bigr)^{1/2}$ for $j=0..4$ ($\partial^{-1}r:=S_1$).
- **(C) Tip constants (T7),** for each certificate interval $I_c\ni\tau_c$:
  - (C1) *Projected coercivity at $\tau_c$.* With $M_{ij}=\langle\chi g_i,\varphi_j\rangle_\rho$ and $\ell_j=\rho^{-1}\mathcal L^\dagger(\rho\varphi_j)$, $\mathcal L^\dagger\psi=-\psi''''+(b_2\psi)''-(b_1\psi)'+b_0\psi$: for all even $u\in H^4_\rho$ supported in $[0,34.25]$ with $\langle u,\varphi_j\rangle_\rho=0$,
  $$\mathfrak B_{\tau_c}(u):=\langle\mathcal Lu,u\rangle_\rho-\sum_{i,j}\langle\chi g_i,u\rangle_\rho\,(M^{-T})_{ij}\,\langle u,\ell_j\rangle_\rho\ \le\ -(\omega_c+m_c)\lVert u\rVert_\rho^2-\theta_c\lVert u''\rVert_\rho^2 .$$
  - (C2) *Spacing:* $s_c(\tau)\le m_c$ on $I_c$, with $s_c$ the explicit function of §2.6.
  - (C3) *Modulation:* $\bar m\ge\lVert M^{-1}\rVert$, $\bar\ell\ge(\sum_j\lVert\ell_j\rVert_\rho^2)^{1/2}$, $c_\varphi\ge\max_j\lVert\rho^{-1}(\rho\varphi_j)''\rVert_\rho$ (S1 reference **[S3]**: $\lVert M^{-1}\rVert\le2.0005$, $\lvert(a,b)\rvert\le3.708A$).
  - (C4) *Gårding, coupling and $E_I$ constants over the base* (sups of explicit expressions in the $b_j$; S1 reference values **[S3]**): $G_{22},G_{12},G_{02},G_m$ (14.34, 16.3, 10.7, 43.33); $a_f,b_f$ (0.0696, 0.0394); $c_{0..3}$ (0.425, 0.424, 0.526, 0.442); strip $s_A,s_B,d_{0..3}$; $\mu_1(\tau),\mu_2(\tau)$ (0.1552, 0.6555) on the *whole* far curve; modulation norms $\gamma_0=\max_i\lVert\zeta_\infty g_i'\rVert_{W_\infty}$, $\gamma_2$; nonlinear Hessian sups $n_4,n_3,n_2$ (§2.5).
- **(F) Far constants (T8),** per slab: rate defects $\sigma_j$, level couplings $c_{ji}$ ($i<j$), ramp inflow constants $k_{jm}$ ($m\le j$), modulation norms $\gamma^F_j$.

### 1.3 Statement
**Theorem T5.** Assume (I), (B), (R), (C), (F), and that the bars $\bar Y,\bar F_j$ defined by the comparison system of §2.8 stay inside the box 𝔹 of §2.5 on $[\tau_0,\tau_1]$. Then the CDF solution through the P1 ball exists up to a time $t_1$, and there are $C^1$ functions $\lambda(\tau)>0$, $c(\tau)$, $t(\tau)$ with $t_\tau=\lambda^4(1-\vartheta)$, $\lvert\vartheta\rvert\le10^{-3}$, such that for all $\tau\in[\tau_0,\tau_1]$ ($\lvert\tau_0-\tau_a\rvert\le\delta$)
$$\frac{\gamma(t(\tau))-c(\tau)e_1}{\lambda(\tau)}=\tilde\Gamma(\cdot,\tau)+e(\cdot,\tau)\,\tilde\nu(\cdot,\tau),\qquad\langle\chi e,\varphi_j\rangle_\rho=0\ (j=1,2,3),$$
$$\sqrt{\Lambda_e(\tau)}\le\bar Y(\tau),\qquad\mathfrak F_j(\tau)\le\bar F_j(\tau)\quad(j=0,\dots,4).$$
At $\tau_1$ this is the enclosure "$\tilde\Gamma(\tau_1)$ + ball" from which T9 checks hand-off items 1–5 (§2.9). The sizes of $\ln(\lambda/\tilde\lambda)$, $c-\tilde c$ and of the clock shift are not constrained and not needed.

## 2. Proof

### 2.1 Representation and error equation
For the true solution put $\Gamma=(\gamma-ce_1)/\lambda$, viewed as a function of the label τ, with $t_\tau=\lambda^4(1-\vartheta)$, $a:=\tfrac14+\lambda_\tau/\lambda$, $b:=-c_\tau/\lambda$. Then, exactly,
$$\langle\Gamma_\tau,\nu\rangle=-(1-\vartheta)\kappa_{ss}+(\tfrac14-a)\langle\Gamma,\nu\rangle+b\langle e_1,\nu\rangle .$$
Subtract the base identity $\tilde V=-\tilde\kappa_{xx}+(\tfrac14-\tilde a)\tilde p+\tilde b\tilde n+r$, write the true curve as a normal graph over $\tilde\Gamma(\tau)$, and use $\tilde\kappa_{xx}=-\tilde V+(\tfrac14-\tilde a)\tilde p+\tilde b\tilde n+r$ to turn the ϑ-direction into $\tilde V$:
$$e_\tau=\mathcal Le-\sum_{i=1}^3x_i\,g_i-(1-\vartheta)\,r+\mathcal N(e)+\mathcal Q(x;e),\qquad g=(\tilde p,\ \tilde n,\ \tilde V/n),$$
$$x_1=(a-\tilde a)-\vartheta(\tfrac14-\tilde a),\qquad x_2=-(b-\tilde b)-\vartheta\tilde b,\qquad x_3=n\,\vartheta .$$
$\mathcal N$ is the quadratic remainder of the graph operator. $\mathcal Q$ is bilinear in $(x,\text{jets of }e)$: $-x_1(e-\tilde qe')$, $x_2\tilde t_1e'$, and $\vartheta\times$(the linearised $\kappa_{ss}$), which only multiplies the principal part by $(1-\vartheta)$.
Regularity: $e''=\kappa-\tilde\kappa+\dots$ with κ smooth, so $e\in W^{m+3,\infty}$ if $\tilde\kappa\in C^m$ piecewise smooth. The tip identities need $m=1$, $E_I$ needs $m=2$, level $j$ of the far field needs $m=j$. This is why (B1) asks for $C^4$. In τ everything is piecewise $C^1$ with $e$ continuous; the inequalities hold a.e. and are integrated.

### 2.2 Modulation (two symmetries plus the clock)
Differentiating the three constraints and integrating by parts onto the smooth $\varphi_j$ (as in MOD_IBP **[S3]**; no $D_4$, $W_3$ appear):
$$\sum_ix_i\bigl(M_{ij}+E_{ij}(e)\bigr)=\langle e,\ell_j\rangle_\rho-(1-\vartheta)\,\bigl\langle S_2,\rho^{-1}(\rho\varphi_j)''\bigr\rangle_\rho+\langle\chi\mathcal N,\varphi_j\rangle_\rho+\text{strip},$$
with $\lvert E_{ij}\rvert\le\bar e\,(A+W_1)$ from $\mathcal Q$. Hence, in the box,
$$\lvert x\rvert\le\frac{\bar m}{1-\bar m\bar e(A+W_1)}\Bigl(\bar\ell A+c_\varphi\mathcal R_2+n_\varphi\,\eta_1(A+W_1+B)+\text{strip}\Bigr).\tag{M}$$
The residual enters through $S_2$ only. Existence and continuity of $(\lambda,c,t)$ along the flow is the implicit function theorem with matrix $M+E$, valid while $\bar m\bar e(A+W_1)<\tfrac12$ (margin ≥ 1e6).
At the start (T4), the same IFT gives $(\lambda_0,c_0,\tau_0)$ near $(\tilde\lambda,\tilde c,\tau_a)$ with $\lvert\tau_0-\tau_a\rvert\le\bar m\,\lVert\text{ball}\rVert/n(\tau_a)$. This is why the base starts at $\tau_a-\delta$.

### 2.3 The third direction: a triangular 2×2 system
From (M), $\lvert\vartheta\rvert=\lvert x_3\rvert/n$. The accumulated clock shift $s$ ($s_\tau=\vartheta$), or equivalently the scalar $a_3:=n\,s$ of the brief, obeys
$$\frac{d\lvert a_3\rvert}{d\tau}\le\nu_3(\tau)\lvert a_3\rvert+c_{12}Y+\rho_3,\qquad\frac{dY}{d\tau}\le-\varrho_cY+\Phi_c,\qquad \nu_3=\frac{n_\tau}{n}\approx\tilde\lambda=-0.256,\quad c_{12}=\bar m\bar\ell,\quad\rho_3=\bar mc_\varphi\mathcal R_2 .$$
$Y$ is not forced by $a_3$ because the shift is exact: $\tilde\Gamma(\tau+s)$ is again the base. Only $\lvert\vartheta\rvert\le10^{-3}$ is used downstream. It holds if $c_{12}Y+\rho_3\le10^{-3}n$; at the end $n\approx0.256\lvert\alpha\rvert\lVert\chi\hat\psi_3\rVert\approx2.7\cdot10^{-4}$, so the margin is ≥ 1e3 for $Y\le10^{-10}$. T7 certifies $n_{\min}$.
The same bookkeeping gives $\lvert\ln(\lambda/\tilde\lambda)\rvert_\tau\le\lvert x_1\rvert+\tfrac14\lvert\vartheta\rvert$ and the drift of $c$; none is used.

### 2.4 Energy identities
Insert (M) into the equation. The part linear in $e$ of $\tfrac12\frac d{d\tau}\lVert\chi e\rVert_\rho^2$ is exactly $\mathfrak B_\tau(u)$ plus strip commutators $[\mathcal L,\chi]$. With $D_\rho u:=\rho^{-1}(\rho u'')''$, $\Phi_1=1.083$ and $\Phi_2:=\sup_{[0,34.25]}\lvert\rho''/\rho\rvert$, the four inequalities are (step-3 structure **[S3]**, constants from (C)):

$$\tfrac12\tfrac d{d\tau}A^2\le-(\omega_c+m_c-s_c)A^2-\theta_0B^2+(s_A+s_B)V_0^2+\mathcal R_2\bigl(B+2\Phi_1W_1+(\Phi_2+\bar m\,c_\varphi)A\bigr)+N_0\tag{E0}$$
$$\tfrac\varepsilon2\tfrac d{d\tau}B^2\le\varepsilon\Bigl[-\tfrac12D_4^2+G_{22}B^2+G_{12}W_1B+G_{02}AB+G_mA^2+\mathcal R_s(D_4+2\Phi_1W_3+\Phi_2B)+\text{strip}_2\Bigr]+N_2\tag{E2}$$
$$\tfrac12\tfrac d{d\tau}V_0^2\le-\mu_1V_0^2-\tfrac12V_2^2+a_fW_1^2+b_fB^2+\mathcal R_1(V_2+2\Psi_1V_1+\Psi_2V_0)+\lvert x\rvert\gamma_0V_0+N_{I0}\tag{I0}$$
$$\tfrac{\varepsilon_I}2\tfrac d{d\tau}V_2^2\le\varepsilon_I\Bigl[-\tfrac12V_4^2-\mu_2V_2^2+(c_3D_4+c_2W_3+c_1B+c_0W_1+\mathcal R^I_s+\lvert x\rvert\gamma_2)(V_4+2\Psi_1V_3+2c_{1W}V_2)\Bigr]+N_{I2}\tag{I2}$$

Notes.
- **Weak pairing.** $\langle\chi r,u\rangle_\rho=\int S_2\,(\rho\chi u)''$, so (E0) needs only $\tilde\kappa$. In (I0), $\langle\zeta r_x,v\rangle_W=\int S_1(W\zeta v)''$ needs $\tilde\kappa_x$ and is absorbed by $\tfrac12V_2^2$. The ε-levels pair $r$ and $r_x$ strongly ($\tilde\kappa_{xx}$, $\tilde\kappa_{xxx}$) but carry the factors ε, ε_I.
- **(E2)** is the identity $\langle u_\tau,D_\rho u\rangle_\rho$; it needs $u''''\in L^2$ only, and coefficients $b_j$ undifferentiated.
- **(I0).** $v_\tau=\partial_x[-v'''+b_2v'+b_1v+(b_0-\tfrac14+\tilde a)e]+(\tfrac14-\tilde a)v$ is in divergence form. The term with $e$ itself is handled by Hardy ($\int x^{-6}e^2\le\tfrac4{25}\int x^{-4}e'^2$ + ramp term) and is inside $\mu_1$. $\mu_1,\mu_2$ are sups of explicit pointwise expressions; they need no spectral information.
- **Far end.** $b_1(\tilde\ell)=\varpi(\tilde\ell)=-\tilde\ell_\tau$ (there $\tilde q=\tilde t_1=\tilde\kappa'=0$), so the transport boundary term cancels the moving-domain term exactly. No boundary terms occur, here or in §2.7 **[D]**.
- **Modulation in (I0).** The coupling $\bar m\bar\ell\gamma_0AV_0$ is new relative to step 3 (there $\beta$ was small because $f_3$ was the adjoint). Expected $\bar m\bar\ell\gamma_0\approx0.15$–0.4 **[D]**. It is why $\mu_e$ may have to drop to ≈ 0.02.

**Closure on $I_c$.** Let $z=(A,B,D_4,V_0,V_2,V_4)$. The terms of (E0)+(E2)+$\mu_e$[(I0)+(I2)] that are quadratic in $z$ form $-\mathfrak Q_c(z)$; the step-3 verifier (`rhs_iv` + branch-and-bound over shapes), run in *linear mode* with the constants of (C), certifies $\varrho_c$ with
$$\mathfrak Q_c(z)\ \ge\ \varrho_c\Lambda_e(z)+\tfrac{\theta_0}4B^2+\tfrac\varepsilon8D_4^2+\tfrac{\mu_e}8V_2^2+\tfrac{\mu_e\varepsilon_I}8V_4^2 .$$
Young on the residual terms then gives $\tfrac12\Lambda_e'\le-\varrho_c\Lambda_e+a_cY+\Xi_c^2$, with
$$a_c=(\Phi_2+\Phi_1(2\sqrt{c_1}+1)+\bar mc_\varphi)\mathcal R_2+\mu_e^{1/2}(\Psi_2+\bar m c_\varphi\gamma_0)\mathcal R_1,\qquad\Xi_c^2=\frac{(1+\Phi_1)^2\mathcal R_2^2}{\theta_0}+2\varepsilon\mathcal R_s^2+2\mu_e\mathcal R_1^2+2\mu_e\varepsilon_I(\mathcal R^I_s)^2 ,$$
and therefore
$$Y(\tau)\le\bar Y(\tau),\qquad\bar Y'=-\varrho_c\bar Y+\Phi_c,\qquad\Phi_c:=a_c+\sqrt{\varrho_c^+}\,\Xi_c,\qquad\bar Y(\tau_0)=\max\bigl(Y(\tau_0),\ \sup\Xi_c/\sqrt{\varrho_c^+}\bigr).\tag{Y}$$
Indicative size: $\varrho\approx0.03$–0.04, $\bar Y\approx200\,\mathcal R_2+2\mathcal R_1+0.02\,\mathcal R_s+0.01\,\mathcal R^I_s$ **[D]**.

### 2.5 Nonlinear terms
Box 𝔹 (bootstrap, closed strictly at the end): $\sup_{[0,34.25]}(\lvert e\rvert,\lvert e'\rvert)\le\eta_1=10^{-5}$; $B\le\eta_B=10^{-4}$; $\sup_{x\ge10.25}x^{j-1}\lvert\partial^je\rvert\le\eta_F=10^{-2}$ ($j=1,2,3$); $\lvert\vartheta\rvert\le10^{-3}$.
The graph operator is quasilinear, linear in $e''''$ and $e'''$, with top coefficient depending on $(e,e')$ only. So pointwise, with Hessian sups over the box (arb automatic differentiation of the explicit rational expression; they depend on $\tilde\kappa,\tilde\kappa',\tilde\kappa''$ and are O(1–30)):
$$\lvert\mathcal N(e)\rvert\le n_4(\lvert e\rvert+\lvert e'\rvert)\lvert e''''\rvert+n_3(\lvert e\rvert+\lvert e'\rvert+\lvert e''\rvert)\lvert e'''\rvert+n_2(\lvert e\rvert+\lvert e'\rvert+\lvert e''\rvert)^2 .$$
- $L^2$ level: every term of $\langle\chi\mathcal N,u\rangle$ contains a factor $e$; put the two highest derivatives in $L^2$: $N_0\le\eta_1\bigl[n_4D_4A+n_3(2AW_3+BW_3)+n_2(A+W_1+B)^2\bigr]$.
- ε level: $\lVert\chi\mathcal N\rVert_\rho\le\eta_1(n_4D_4+n_3W_3)+n_3S_2^eW_3+\dots$ with $S_2^e=\sup\lvert e''\rvert$. Use the regional Agmon bounds **[S3]** $S_2^{e\,2}\le K(\Phi_1B^2+2BW_3)$ and $W_3^2\le BD_4+c_1B^2$, so $S_2^eW_3D_4\le CB^{5/4}D_4^{7/4}\le\delta D_4^2+C'\delta^{-7}\eta_B^8B^2$.
- Net effect: $\omega_c\to\omega_c-\delta_N$ with $\delta_N\le C(n_4^2/\varepsilon)\eta_1^2+Cn_2\eta_1\le10^{-3}$ for $n\le30$. $\mathcal Q$ is treated the same way with $\lvert x\rvert$ from (M).
- The box closes from $\bar Y$ by the step-3 pointwise constants ($K_{\rm reg}\le3.02\cdot10^4$, $K_W\le2.3\cdot10^6$): $\sup\lvert e\rvert,\lvert e'\rvert\lesssim2\cdot10^3\,\bar Y/\sqrt{\mu_e}\le10^{-7}$ at $\bar Y=10^{-11}$. Margins ≥ 50.

### 2.6 Coercivity between certificate times
For even $u$ supported in $[0,34.25]$: $\langle\mathcal Lu,u\rangle_\rho=-\lVert u''\rVert^2+\int\rho P_1u'^2+\int\rho P_0u^2$ with
$$P_1=2\rho''/\rho-b_2,\qquad P_0=b_0-\tfrac12(\rho b_1)'/\rho+\tfrac12(\rho b_2)''/\rho-\tfrac12\rho''''/\rho .$$
$P_0\to-\infty$ like $-x^{4/3}/24$; the far part of the ρ-zone is confined by the weight. Fix $C_0$ with $\Pi_c:=C_0-P_0(\tau_c)\ge1$. For $\tau\in I_c$ let
$$d_1=\sup\lvert P_1(\tau)-P_1(\tau_c)\rvert,\quad\vartheta_c=\sup\frac{(P_0(\tau)-P_0(\tau_c))_+}{\Pi_c},\quad\bar p_1=\sup(P_1(\tau_c))_+,\quad\lVert\Pi^{(3)}_{\tau_c}\rVert\le\bar\pi_c,\quad\lVert\Pi^{(3)}_\tau-\Pi^{(3)}_{\tau_c}\rVert\le\Delta_\Pi ,$$
where $\Pi^{(3)}$ is the rank-3 form in $\mathfrak B$ (norms from 3×3 Gram matrices). Then, with $e_1:=\vartheta_c\bar p_1+d_1$ and any $\theta_1>0$,
$$\mathfrak B_\tau(u)\le-(\omega_c+m_c-s_c)\lVert u\rVert^2-\theta_0\lVert u''\rVert^2,\qquad\theta_0=(1-\vartheta_c)\theta_c-\theta_1,$$
$$s_c(\tau)=\vartheta_c(\omega_c+m_c+C_0+\bar\pi_c)+\Delta_\Pi+c_1e_1+\frac{e_1^2}{4\theta_1},\qquad c_1=0.5804 .$$
Proof: $\int\rho\,\delta P_0u^2\le\vartheta_c\int\rho\Pi_cu^2=\vartheta_c(C_0\lVert u\rVert^2-\langle\mathcal L_cu,u\rangle-\lVert u''\rVert^2+\int\rho P_1^cu'^2)$, then $W_1^2\le AB+c_1A^2$ and (C1). The relative form of $\vartheta_c$ is essential: the unweighted sup of $\delta b_1$ at $x\approx30$ is O(1) per unit τ at small k.
Placement: FEAS2 measured a relative form variation ≈ $1.2k^{-1.05}$ per unit of $s$. Certificates are cheap ($N\le60$), so do not economise: one every Δτ ≈ 0.6 (≈ 45–50 in all), and ask $s_c\le m_c/2$. One S1-based certificate with the tail reference serves all $k\ge1600$; there (C2) is checked per slab.

### 2.7 Far field (Lemma F)
In the far field the equation is $v_\tau\approx-v''''-\tfrac x4v'$ plus terms of relative size $x^{-4}$. With weight $x^{-m}$ the level-$j$ rate is $(1-m)/8-j/4$, so the scale-covariant choice $m=2-2j$ damps every level at rate 1/8 (1/4 in the energy). This is the Lemma S damping $\gamma'\approx0.246$. In physical terms $\mathfrak F_j/(28\lambda)$ is conserved up to inflow: the far error is frozen.

**Lemma F [D].** In 𝔹, for $j=0,\dots,4$:
$$\mathfrak F_j'\le-(\tfrac14-\sigma_j)\mathfrak F_j-\mathfrak D_j+\sum_{i<j}c_{ji}\sqrt{\mathfrak F_i\mathfrak F_j}+\sum_{m\le j}k_{jm}\bigl(V^{\rm ramp}_m\bigr)^2+4(\mathcal R^F_j)^2+4\lvert x\rvert^2(\gamma^F_j)^2 ,$$
with $V^{\rm ramp}_m=\lVert\partial^mv\rVert_{W_\infty,[14,28]}$.
- Proof: test $\partial^jv$ in divergence form and move two derivatives to the test function (Lemma S Prop. 1 **[S3]** for $j\le2$, the (B) code for $j=3$; $j=4$ is the same computation). Level $j$ then needs ramp data of order ≤ $j$ and $\partial^{j-1}$ of the coefficients, i.e. $\tilde\kappa^{(j+1)}\in L^\infty$. The system is linear and triangular in $j$. There are no $L^4$ terms: quadratic terms are bounded with $\eta_F$ and go into $\sigma_j$.
- $\sigma_j$ collects $\sup x^{-4}\times$(scale-invariant curvature products) and $\lvert b_1+x/4\rvert/x$. Float, all $k\ge10$: $\lvert b_1+x/4\rvert/x\le0.05$; the level-0 pointwise rate with weight $x^{-2}$ is ≤ −0.074 on [16.25, 60] and ≤ −0.11 beyond, against −1/8 **[F]**. So $\sigma_0\le0.1$ and the energy rate stays ≥ 0.15.
- $c_{ji}$ can be large on the macroscopic bend (there $x\tilde\kappa\approx19$ **[F]**), but $\int c_{ji}\,d\tau\lesssim s^{j-i+2}\lvert\partial_s^{j-i}\kappa^2\rvert_{\rm phys}\,(T-t_a)$ is finite and only feeds higher levels, whose tolerance out there is crude.
- **Inflow.** $\int k_{jm}(V_m^{\rm ramp})^2\le k_{jm}J_m$, $J_m=\int V_m^2d\tau$ over the transit window Δ ≈ 3. From (I0), (I2), (E2) integrated in time (the V1 lemma of the verifier **[S3]**) and a mean-value choice of the initial time in each step: $J_0\le\Delta\bar Y^2/\mu_e$; $J_2\le2V_0^2(\tau')+\dots$; $J_4\le2V_2^2(\tau^*)+\dots$ with $V_2^2(\tau^*)\le J_2/\Delta$; $J_1,J_3$ by interpolation. This avoids the factor $1/(\mu_e\varepsilon_I)$ that $\Lambda_e$ alone would give for $V_2$.
- **Expected amplification** $\sqrt{\bar F_j}/\bar Y$, from the step-3 constants $K_j=(1.72\cdot10^5,8.56\cdot10^7,4.29\cdot10^{10})$ and a factor ≈ 500 per further level **[D, to be certified by T8]**: ≈ 3e3, 1e5, 2e6, 8e7, 3e9 for $j=0..4$. Step 3's own numbers give the same per-level ratios (1e4, 2.5e5, 9e6 for $j\le2$).
- **Pointwise** on a window $[y,1.5y]$: $\sup\lvert x^j\partial^jv\rvert^2\le\frac{1.5y}{28}\bigl(2.47\,\mathfrak F_j+2\sqrt{\mathfrak F_j\mathfrak F_{j+1}}\bigr)$, $j\le3$.

**Resulting accuracy requirement.** Item 3 at level 3 needs $\sqrt{\mathfrak F_3}\le\sqrt{10^{-5}}-\sqrt{10^{-6}}=2.2\cdot10^{-3}$, i.e. $\bar Y\le2.5\cdot10^{-11}$. Item 4 ($K_4$, pointwise level 3 at $y=50$, tolerance 0.05) needs $\bar Y\le4.5\cdot10^{-11}$. Levels 0–2 need only $\bar Y\le10^{-9}$. **Target: $\bar Y\le10^{-11}$ on all of $[\tau_0,\tau_1]$** (there is no useful contraction between early and late times).
The material points of the windows $y\le380$ at $\tau_1$ crossed the ramp at $k\gtrsim850$. Points that were never in the tip (physical $s\gtrsim14\lambda_a$) keep their P1 error, plus $\int(\mathcal R^F)^2$, and need only the crude item-5 tolerances.

### 2.8 Closing the argument
Comparison system: (Y) for $\bar Y$; the linear triangular system of Lemma F for $\bar F_0..\bar F_4$ with inflow from $\bar Y$; the scalar bound of §2.3 for ϑ. Continuity argument on the largest interval where 𝔹 holds: (M), (Y), Lemma F, and the pointwise bounds of §2.5 and §2.7 give 𝔹 with strict inequalities; the bounds control curvature, so the solution continues (the continuation fact P1 of step 3 **[S3]**). No tail lemma and no first-touch barrier are needed: the system is linear with small coefficients, and the horizon is finite.

### 2.9 From the enclosure to items 1–5 at $\tau_1$ (T9)
All functionals are evaluated on $\tilde\Gamma(\tau_1)$ in arb; the error enters by the triangle inequality in the corresponding norm. Let $\kappa_{\rm ch}\le1.02$ bound the chart conversion (normal graph over the base → graph over S1; the base is within 2e-4 of S1 on the ρ-zone, and has slope ≤ 0.1 over the ray up to 380).
1. IFT for $(x_1,\lambda_1,\alpha_1)$ with the step-3 matrix **[S3]**: $\lvert\alpha_1-\tilde\alpha_1\rvert\le C_\alpha\bar Y$, and $\sqrt{\Lambda[w]}\le\sqrt{\Lambda[\tilde w]}+C_\Lambda\sqrt{0.08/\mu_e}\,\bar Y$, where $\tilde w=\tilde u-u_0(\tilde\alpha_1)$ and $C_\alpha,C_\Lambda\le50$ contain $\lVert M_{S3}^{-1}\rVert$ and the Λ-norms of $p,n_1,\hat\psi_3$. The error share is ≤ 1e-9 against 5e-8 − 4e-10.
2. Shells [14, 60]: $\sqrt{E_j}\le\sqrt{E_j[\tilde v]}+\kappa_j\sqrt{\bar F_j}$, $j=0,1,2$, $\tilde v$ = base deviation; $\kappa_j$ = chart factor × the triangular $(y\partial_y)^j\leftrightarrow x^j\partial^j$ factor.
3. Windows: the same with $j\le3$ on each $[y,1.5y]$ ($\mathfrak F_j$ restricted to a window is ≤ $\bar F_j$).
4. $K_4$: jet of the base on $y\in[50,350]$ plus the pointwise level-3 bound of §2.7, in cylinder units $L=50$.
5. Far cylinders: base jets plus pointwise levels ≤ 3 (for $h'',h''',h''''$), $\lVert h'''\rVert_{L^2}$ from $\bar F_2$, $\int(h^{(5)})^2$ from $\bar F_4$, slope from levels 0–1; graph property (H4) from the level-0 pointwise bound.

### 2.10 Where a number must be certified
Every symbol in (C), (R), (F); $n_{\min}$; $\Phi_2$, $c_\varphi$ (fixed functions, once); the conversion constants $\kappa_j,C_\alpha,C_\Lambda$ (T9); the J-cascade constants (verifier V1); the box margins. Nothing else is numerical.

## 3. Interfaces

| task | must enclose (set, norm) | target |
|---|---|---|
| **T1** base | $\tilde\kappa$: piecewise polynomial, $C^4$ in $x$ (B-splines of degree ≥ 6 on the geometric mesh, or Hermite-$C^4$ interfaces), exact symmetry at both ends, polynomial in τ per slab, exactly continuous across slabs; scalars $\tilde a,\tilde b$; normalisation $\tilde\kappa(0)=K_0$, tip at the S1 tip position. Weak (Legendre–Galerkin) discretisation, so that $S_1,S_2$ vanish at element ends up to quadrature. Polished against arb residuals at 128 bits. Starts at $\tau_a-0.02$ from P1's centre. | meets the T6 targets; ≈ 600 DOF for θ at 1e-15, ≈ 900 with the residual (F1 + 45 DOF per decade) **[D]** |
| **T6** residuals | $\mathcal R_2$ (ρ, $[0,34.25]$) | ≤ 2e-14 |
| | $\mathcal R_1$ ($W_\infty$, whole far curve) | ≤ 5e-13 |
| | $\mathcal R_s$, $\mathcal R^I_s$ (strong) | ≤ 5e-11, 1e-10 |
| | $\mathcal R^F_j$, $j=0..4$ | ≤ 1e-7 ($j\le3$), 1e-5 ($j=4$) |
| | jets of the base on each slab: $\tilde\kappa^{(i)}$ ($i\le5$), $\tilde p,\tilde q,\tilde n,\tilde V,\varpi$, as arb hulls per element | rel. 1e-6 is enough |
| **T7** tip | (C1) at ≈ 45–50 times $\tau_c$, spacing Δτ ≤ 0.6; port of `p1.py` to arb with (i) the projected form (rank-3 term exact in the Galerkin block, $\bar\pi$ in the tail), (ii) general constraints (null space of $\Phi_N^T$ in the block, leak $\lVert P_{>N}\varphi_j\rVert$), (iii) $\theta_c=0.03$ | $\omega_c=0.03$, $m_c\ge0.02$ (float truth 0.05–0.07 at θ = 0.02–0.03); $\omega_c\le0$ admissible for $k<40$ |
| | (C2) $s_c(\tau)$ per slab | ≤ $m_c/2$ |
| | (C3) $\bar m,\bar\ell,c_\varphi,n_{\min}$; (C4) Gårding, ramp, strip, $\mu_1,\mu_2$ on the whole far curve, $\gamma_0,\gamma_2$, $n_{2,3,4}$; then the verifier in linear mode → $\varrho_c$, $\mu_e$ | $\varrho_c\ge0.02$; $\int\varrho_c^-\le0.5$ |
| **T8** far | per slab: $\sigma_j$, $c_{ji}$, $k_{jm}$, $\gamma^F_j$; J-cascade; integration of the $\bar F_j$ system; pointwise constants | $\sqrt{\bar F_j}\le$ (3e-8, 1e-6, 2e-5, 8e-4, 3e-2) on the part ending in $y\le380$; ≤ 1e-2 elsewhere |
| **T9** hand-off | functionals of $\tilde\Gamma(\tau_1)$ for items 1–5; $\kappa_{\rm ch},\kappa_j,C_\alpha,C_\Lambda$; the box check | truth ≥ 10× inside (M1), error share ≤ 2× |
| **T4/P1** | ball at $t_a$: C² ≤ 3e-14, C⁵ ≤ 1e-9 (physical); gives $Y(\tau_0)\le2\cdot10^{-12}$ and $\mathfrak F_j(\tau_0)$ | uses ≈ 20% of $\bar Y$ |

Budget check: $\bar Y\approx200\cdot2\cdot10^{-14}+2\cdot5\cdot10^{-13}+0.02\cdot5\cdot10^{-11}+0.01\cdot10^{-10}\approx7\cdot10^{-12}$, against the target $10^{-11}$ and the need $2.5\cdot10^{-11}$. The factors 200, 2, 0.02, 0.01 are estimates **[D]**; T7 replaces them by certified $a_c,\Xi_c,\varrho_c$.

## 4. Risks

**Hardest step: the far field at levels 3–4 (T8 with T1/T6).** It sets $\bar Y\le10^{-11}$, hence weak residuals of 2e-14, hence an arb-polished $C^4$ base; and Lemma F at levels 3–4 is new (linear, but with five levels of commutators). The amplification factors 8e7 and 3e9 are estimates. Fallbacks, in order:
1. *Time-differentiated energies.* $e_\tau$ satisfies the same equation with forcing $\mathcal L_\tau e+r_\tau$; $\Lambda_e[e_\tau]$ with the same certificates gives $D_4$, $V_4$ pointwise without the $1/\sqrt\varepsilon$ losses. Expected gain 10–60× at levels 3–4, so $\bar Y\approx5\cdot10^{-10}$ would do. Cost: $r_\tau$ in T6, P1 data in C³ of θ.
2. *Smoother exit.* Levels 3–4 only with a $C^2$ base, by an interior smoothing lemma for the true flow on cylinders (Lemma-O type; it has four orders of slack). New nonlinear lemma; avoid if possible.
3. *Ask step 3 for less.* Only $E_3$ on the outer windows and $K_4$ bind. If Lemma S(B) can take level 3 from levels ≤ 2 plus a crude level-4 bound, the requirement drops to $\bar Y\le10^{-9}$.

**Second risk: the certificate (C1) at small k (T7).** The true margin is 0.06–0.08, a third of what FEAS2 assumed, and FEAS2's Galerkin+tail scheme lost 0.05–0.1 at its first closing $N$. Expect $N\approx40$–80 at $k\le40$ and $N\approx10$–20 for the S1 certificate (not 3). It also needs a rigorous lower bound for the reference form on the complement of the Galerkin space (step-3 linear-certificate technique). Fallbacks: (a) accept $\omega_c\le0$ for $k<40$: over Δτ = 4 at rate 0.05 the loss is a factor 1.2; (b) raise $k_a$; (c) a Lyapunov matrix on the first $N$ modes instead of the identity (removes the non-normality loss; gap 0.7). The 2-direction scheme (row 6 of §0) is *not* a fallback for accuracy: its growth is only $e^{2}$, but the error then accumulates instead of saturating.

**Reuse of the step-3 verifier.** Yes for the engine: `rhs_iv` with branch-and-bound over shapes, the comparison ODE, and the V1 time integrals, driven by a per-interval constants table. It runs in linear mode (no α, no tail lemma, no barrier march). No for the linear theorem and its constants: they are S1-specific, with the adjoint $f_3$, and the perturbation $\approx(60$–$100)\lvert\alpha\rvert$ fits the margin 1/32 only for $k\gtrsim10^4$, the last 3 of 28 units of τ. Lemma S Prop. 1 identities are reused for Lemma F; Lemma S itself is not (its boxes hold only for $k\gtrsim2200$, and it is nonlinear).

**Starting later.** Recommended $k_a\approx20$ (range 14–28): coercivity 0.079 against 0.072 at 14 and 0.061 at 10; the base is half as far from S1, so $N$ and the spacing losses halve. P1 then costs ≈ 5–8 h instead of ≈ 2 h (T3 fit: ≈ 23 slabs per doubling of k, time per slab × 2.5–3 per doubling), still with block resizing (T3 R-a) untested. Beyond $k\approx40$ nothing is gained in P2 (0.083 is reached) and P1 passes 15 h. Rule: take the smallest $k_a$ at which (C1) closes with $m_c\ge0.02$ at $N\le40$.

**Expected P2 cost.** Δτ = 4 ln(2.3e4/20) ≈ 28, i.e. ≈ 190 slabs of 0.15. Per slab at 128 bits: base polish and residual enclosure 2–8 min, tip and far constants ≈ 1 min; ≈ 50 certificates at ≈ 2 min. Total ≈ 10–30 CPU-hours on one core. The human cost is in T8 (Lemma F, 2–3 weeks) and T7 (rigorous tail for the projected form).

**Other caveats.**
- Everything in §0 is float, one trajectory, FEM on [0, 40]; the ε-level of (C1) is covered by Gårding, not by a spectral check.
- The sign conventions of $\mathcal L$ with $(\tilde a,\tilde b,\varpi)$, and $b_1(\tilde\ell)=-\tilde\ell_\tau$, must be confirmed symbolically.
- $\mu_e<0.08$ costs $\sqrt{0.08/\mu_e}$ in item 1 only.
- $n\to0$ like $\lvert\alpha\rvert$: harmless up to $k_1$, but the scheme is not uniform as $k\to\infty$ (it does not have to be).
