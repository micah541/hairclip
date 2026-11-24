# Part C: interface lemmas `fix14_interface/` (C-1 identification, (GS)) -- fork of the fix15 reviewer

Scope: `I` = `bulk/fix14_interface`. Proof tree read-only; everything written is in `fix15_review/c/`. `I/SHA256SUMS` verifies after my runs.
Labels: [cert] arb, outward; [sym] sympy; [hand] my derivation; [read] read from the cited file; [est] estimate.
Runs (flock, nice 19, PYTHONHASHSEED=0, -B; peak RSS 173 MB): copies of `C1_sym.py` (26/26 pass), `C1_ke.py`, `GS_tau1.py`, `C1_num.py` -> JSON equal to
the originals apart from script path / cpu time. `GS_boot_v4c.py` (repointed to fin13v4c) and `c_frame_shift.py` are new.

## Verdicts
- **C-1: sound modulo items.** The chain Lean flow -> gauge -> P1 (Gronwall) -> P2 -> t_1 is correct and not circular. One text item is owed (C1-1: (P2-a) is a
  reading of step 4's proof, not a proved lemma). Inherited, unchanged: Lean not built / n=1 not elaborated; DKS not re-read (PDF is on disk, not opened by me).
- **(GS): sound modulo items** (minor: records and wording). Tube lemma, tau_1 certificate, strictness, open/closed and the endpoint at 60 hold.

## GS_boot rerun on fin13v4c production (`runs/P/rate_table_fin13v4c_a0.00015.json`)
sup|ku| <= 5.3902e-4 (unchanged), sup|u'| <= 0.040935 (unchanged), |w'| on [8,10.25] <= 1.414421e-3 (was 1.414426e-3), sqrtL_max 5.308499e-8 [cert].
All 8 box ratios used are <= 1 (max 0.7891, Bd_w_34_60); the three new v4c asserts are <= 0.116. Only diffs against `I/GS_boot.json`: box ratios in the 5th-6th
digit, sqrtL_max, B_max, W1_on_8_10, region '8-10.25' up (4.16502e-3), paths. The 1.25e-4 table has smaller ratios and sqrtL (5.2142e-8) [read]. Lemma GS-2's numbers stand.

## C-1: what I checked
- **Gauge lemma (c) [hand].** Re-derived from theta_t|_p = V_s, d_t|g_p| = -kV|g_p|, V = -k_ss: s_t|_p = (2pi/L)^3 int_0^a (1+v)w', L_t = -(2pi/L)^3 int w^2,
  d_sigma a|_p = int_0^a(1+v)w' + a<w^2>, and th_sigma = -th'''' - (1+v)^2 w + (1+v) int_0^a (w^2 - <w^2>) = step 4's P1 equation (`step4_proof.tex:227-230`, same N,
  same clock, same slaved equations; X_0' sign and N(p_0) = 1 checked). Agrees with Y3, Y4, Y5, Y7.
- **Sympy map.** Y1 (kinematics), Y2, Y3, Y4, Y5 -> Lemma gauge (c); Y6, Y7 -> parity/Parseval/slaved; Y8 -> Lemma kP2; Y9 -> tube lemma; Y10 -> Lemma G0, gauge (b).
  26 = 2+4+1+2+2+5+3+2+3+2. Independently derived by hand: Y3/Y4/Y5, Y7a-b, Y8b, Y9a-c, Y10. The identities proved are the ones the tex uses.
- **Parity / Parseval [hand].** For smooth odd periodic th: w odd, w^2 - <w^2> even with mean 0, its antiderivative from 0 is odd and periodic, hence mean 0 = step 4's
  zero-mean antiderivative. <w^2> = (1/2) sum n^4 u_n^2 and th'''(0) = -sum n^3 u_n hold for the infinite series (smooth th; for the P1 solution |u_n| <= T n^-6). Correct.
- **Uniqueness lemma [hand].** Correct. Norm: L^2 of delta; needs only sup||v||_inf and sup||w||_{L^2} of both solutions (H^2 would do) and f in C^1 (from C^1L^2);
  all integrations by parts are periodic and legitimate in H^5; ||Pf||_inf <= 2||f||_{L^1} (the 4 used is safe); the nonlocal term is handled by
  P(delta''(w+w~)); no gauge term remains because both solutions are in the same gauge (the gauge is fixed by Lemma gauge before uniqueness is applied).
- **Class of P1 [read].** `step4_proof.tex:244-246, 823-831`: Lemma gal gives exactly C^0(I;H^5) n C^1(I;L^2), solving the untruncated PDE in L^2; existence is
  constructed. Its uniqueness statement is inside W; C1's lemma removes W correctly. P1 starts from Gamma_0 (b_2 interval contains b_2*).
- **Prop P1 (no circularity) [hand].** The a-priori bound int k^2 is used only where th^gamma = u is already proved (on [0,S']); Case B contradiction uses L <= L_0. Correct.
- **P2 [read + hand].** Step 4's Theorem P2 proof (Steps 1-7, `:339-370`) is a-priori for a given smooth solution; the only existence input is the last sentence of
  Step 7 ("(P1) literature [A]") and the theorem's own "exists on [t_96,t_1]". Replacing that by Lean + DKS is legitimate and not circular: Prop P2 argues by
  contradiction on T_*, with int k^2 bounded on [t_96, T_*) from the boxes (Lemma kP2) -- see C1-1 for what is owed. Hand-over at t_96 = switch lemma (`:282-288`),
  base slab 2 from tau_a; C1_ke covers slabs 2..180 (179 files).
- **C1_ke [cert, rerun].** sup|k~ e| = 1.74285e-3 (slab 180, x in [3.55e5, 3.70e5], |e| box 32.2, k~ 5.41e-5); boxes are step 4's (`:314`: |e| <= 1e-5 for y < 28,
  1.7E(x) for y >= 28, E = c_e Ybar + sqrtF0 (y^3/84)^(1/2), c_e = 1.6031e4), at Ybar_box 1.7e-12; jets widened outward. a >= 1/2 with margin x287.
- **DKS.** (D1) as quoted in `P1_CONFIRMATION.md:34` ("f extends smoothly to [0,T] x I ..."). Applied correctly via U-1's blow-up lemma (extend, restart with (L1), glue).
  The C-4 completion (immersion at T from d_t ln|Phi_x| = k k_ss bounded) is right given uniform bounds on k, k_ss; that DKS's proof supplies them is not re-read.

## (GS): what I checked
- **Tube lemma T [hand].** g' = -1 + k<c-S,nu> <= -1 + 3K delta; window argument, cell ends, X' > 0, u' = (1-ku) tan d, coverage: all correct.
- **GS_tau1.py [cert, rerun].** delta, eta, K are what Lemma T needs (mean-value form inside x_13 is right: |D'| <= |th~ - th_S| + |1/lam - 1| + (|k~e| + |e'|)/lam;
  K on [x_lo - 2 delta, x_hi + 2 delta], mirrored). 16 t_c balls of radius 1/16 cover all of slab 180; frame = P_hull +- 1e-9; S_1 chart with tip ball;
  coverage >= 61.27. sup|ku| <= 6.517e-3, sup|u'| <= 0.05418. Inputs hashed (2 JSON, 2 profile files, 3 modules).
- **GS-2 / GS_boot.py [cert, rerun on v4c].** Boxes used are asserted in `step3nl_v4c/verifier/closure.py:86-97`, called at `:145` (each step) and `:164` (anchor).
  The [8,10.25] bound is right (chi = 1, rho decreasing, Lambda >= eps B^2). No circularity: (GS) is used only to define u, w; the bounds come from BS_U.
- **Open / closed.** GS-3 (endpoint: the curve continues past theta(+-60), d_theta X > 0) is correct. GS-4: with (P2') the path P(tau) is C^1 on [tau_1, tau_mod),
  so P_n -> P_b is immediate; the Arzela-Ascoli step is correct (closed set, non-strict bounds c_1, c_2; C^0 limit, then smooth by the local inverse).
- **Range 60.** Verifier boxes end at 60; U-1(b) needs one-sided derivatives at 60 only. Posing (GS) on |x| <= 60 is consistent; tau_1 holds on 60.25 anyway.
- **Section 6.** The seven objections to the v37 sketch are right (|u| reaches 1.096 on [34.25,60] [cert], so sup k * sup u fails).

## Findings
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| C1-1 | `I/C1_identification.tex:36-40, 111` | **major (text owed; no numerical risk)** | (P2-a) -- "maximal tau_e; either tau_e = 27.15 with t < T', or t(tau) -> T'" -- is stated as a reading [hand]. Step 4's Step 7 (`step4_proof.tex:370`) is two sentences. Nowhere is it written why, at tau_e < 27.15 with t(tau_e^-) < T', the decomposition continues: (i) (lambda,c,t) have limits (Lipschitz from Step 2), (ii) the global normal-graph representation over the base on [0, l~] is closed (limit argument as in GS-4, with a >= 1/2 from C1_ke) and open (as in GS-3; both ends on the axis), (iii) IFT at the limit (m||E|| < 1/2), (iv) strict box improvement. All ingredients exist. | Write (P2-a) as a lemma with this proof (about 15 lines, same pattern as GS-3/GS-4); cite C1_ke for (ii). |
| C1-2 | `I/C1_sym.py:72-76` | minor | Y4's inputs d_t s and L_t are typed in, so Y4 checks only the last algebra step; the tex labels the drift formula [sym]. My hand derivation confirms it. | Label "[hand] + [sym]". |
| C1-3 | `I/C1_identification.tex:105` | text | Lemma kP2 mixes x and y ("on y <= 29 (chi = 1) ... rho(31.25)"), and "x >= 10.25, y < 28" leaves the regions to the reader. Qualitative finiteness is right. | Reword with one variable. |
| C1-4 | `I/C1_identification.tex:42, 44, 135-136` | inherited | Lean: read, not built, n = 1 not elaborated (C-3). DKS: (D1) and the uniform sup-bounds on d_s^m kappa not re-read. | Authors / Micah. |
| GS-1 | `I/GS_lemma.tex:57`, `I/GS_tau1.py:10, 48` | minor (records) | The frame widening 1e-9 cites a "true-curve shift <= 3.3e-11 (step4_proof)". No such number is in `step4_proof.tex` or in `t9_h1_clam_s180.json` (clam.py exports only the alpha row, C_alpha). From the stored M, f_norms, tail_t: rows (x_c, ln lam) <= 9.66, 5.18, shift <= 1.65e-11 at Ybar_box [cert on 5-digit ball strings, `c_frame_shift.json`]. The widening covers it x60. | Cite a script; or export the two rows from clam.py at the re-pin. |
| GS-2 | `I/GS_boot.py:61, 66-70` | minor (rigor, negligible) | For x in (X13, T0+5) (width 1e-6) K is the tail envelope at sigma = 13, valid only for sigma >= 13. k ~ 2e-7 there; margin x900. | Drop `max(13.0, ...)` or extend the chart branch to T0+5. |
| GS-3 | `I/GS_lemma.tex:71-83, 122` | minor (hypothesis list; matters for Part D) | Lemma GS-2 uses, besides B1 and the alpha barrier, the verifier input S0V = 0.0153 (sup|w'| on [34.25,60], `closure.py:88, 90`; a Lemma-U/S(A) bootstrap member, not a consequence of Lambda). The tex says "it uses D1's boxes only". | List it: GS strictness depends on B1, alpha <= 2.5e-4 and the B2' bound S0V. |
| GS-4 | `I/GS_lemma.tex:71, 126, 129` | minor (inherited) | That the asserted boxes hold beyond the tail anchor rests on VERIFIER.md's monotonicity statement [read, not audited]. | Part of the E_I lemma text. |
| GS-5 | `I/GS_lemma.tex:71, 83`; `I/GS_boot.py` | text | Refers to fin13v4b, closure.py l.143/162. In fin13v4c: l.145/164; W1 on [8,10.25] 1.4144e-3; results unchanged. | Update at the next text round; pin `GS_boot_v4c.json` or rerun in place. |
| GS-6 | `I/GS_tau1.py:155` | info | Contiguity assert tolerates 1e-9; coverage is structural (s-balls tile [-1,1]), so harmless. | None. |

## Uses / proves (for Part D)
**C-1 uses:** step 4 Prop Gamma_0 [cert]; step 4 Galerkin-limit lemma (existence, class C^0H^5 n C^1L^2) on the 96 P1 slabs; switch lemma; Theorem P2's a-priori
part (Steps 1-6, boxes B u B', constants [C], analytic steps [D]) in the form (P2-a) [C1-1]; C1_ke.json; Lean (L1), (L2) at n = 1; DKS (D1) + C-4 completion;
U-1's lemmas gluing / maximal flow / symmetry / blow-up (which use only (L1), (L2), (D1)). **Does not use:** step 4's (H6), step 3's (P1'), (GS), any step-3 bootstrap.
**C-1 proves:** T_* > t(27.15) >= t_1; gamma smooth, symmetric, immersed on [0, t(27.15)]; P1 and P2 enclosures are enclosures of gamma; (S4.1), (S4.2) of U-1;
U-1's T_max = T_* - t_1.

**(GS) uses:** Theorem C-1(i),(iii) and P2's boxes at slab 180, the hand-off frame P_hull (step 4 H1), S_1 charts (GS-1); BS_U members B1 (Lambda <= barrier,
alpha <= 2.5e-4) and S0V (B2'), the verifier's asserted boxes [cert, fin13v4c], c4fun psi3^, h2^ (GS-2); smoothness of the flow on [t_1, t_1+T_max) and the C^1
path P(tau) of U-1 (P2')(b) (GS-3, GS-4). **(GS) proves:** GS at tau_1 on |x| <= 60.25; under BS_U + GS, strict bounds 5.40e-4 / 0.0410; GS is an open and closed
member of the extended bootstrap. Note the mutual dependence with (P2'): (P2') needs GS in its set B (to define u, w, c_j); GS-4 needs P(tau) from (P2')(b).
This is a joint continuation on [tau_1, tau_mod), not a circle: (P2')(a),(b) need GS only on the closed set / its GS-neighbourhood, and GS-2 never uses (P2').

## Not checked
Lean build and source (not re-read here; fix11 read it at d72b7447, same HEAD now); DKS2002.pdf; step 4's P2 analytic steps [D] and its jets/base engine
(`basegeo`, `s1chart2`, `c4fun` taken as given); VERIFIER.md monotonicity beyond the anchor; `C1_num.py` beyond rerunning it (float sanity only).
