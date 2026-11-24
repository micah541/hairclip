# E2_AUDIT: the E2 energy block of the step-3 tip-zone verifier (fin12f4 / fin13v4)

Nov 20, 2026. Independent audit (Fable). Scope: `step3nl_v4/verifier/rhs_iv.py:243-251` (E2 block), its forcing line `rhs_iv.py:106`, its typed constants `rhs_iv.py:41-43`, `constants.py:14-15`, and the LAMT question of V4_REVIEW section 6 item 3. The E2 lines are byte-identical in `step3nl_v33/verifier` (fin12f4).
Proof tree read only (`step3nl_v4/SHA256SUMS` verifies after my runs; no file under `step3nl*`, `step3`, `s1_tight` is newer than this audit). Everything new is in this directory.
Labels: **[cert]** arb/interval arithmetic on certified inputs (inputs taken as given); **[read]** read off code or a certified JSON; **[hand]** hand derivation, not machine-checked; **[est]** estimate; **[sym]** exact sympy identity.
Paths: `V/` = `step3nl_v4/verifier/`, `S/` = `step3nl_v4/`.

## 0. Verdict

**FIX (small). Not "E2 OK as coded", not serious.**

1. **The modulation term `43.33 A^2` is a correct, correctly rounded bound of one piece only**: the (a,b)-lead coupling `-a<Df1,w~> - b<Df2,w~>` with `(a,b) = Minv (<w~,L*f1>, <w~,L*f2>)` on {f1,f2}^perp (certified 43.3271 [read `agent_mod/mod_forms.json`]). It is not the PI3 pathology: no mistyped norm, no hidden projection. But everything else in the E2 modulation is absent from the code:
   - the beta coupling `-beta <chi psi3 + 2 alpha chi h2, D_rho w~>` (all of it: lead, alpha-linear, quadratic, forcing, strip);
   - (a,b) beyond the 2x2 lead: the f3 row of H (0.0131 A [cert]), the alpha-linear part (622.4 alpha A), the quadratic part (51.34 Q), the forcing part (3.84 Fn), the strip part;
   - the couplings `a <chi A1(u), D_rho w~> + b <chi B1(u), D_rho w~>` (order alpha and order w).
   With these the forcing factor of the E2 level is **13.98** [cert], coded 1 (`V/rhs_iv.py:106`).
2. **The nonlinear constants nd1..nd8, nf1..nf7 are valid, but not for the stated reason.** 8 of the 15 are typed below their own certificate (alpha <= 1e-3, jet box 1e-4), and that jet box does not hold on [8,12] (worst case 1.35e-4 [cert]) nor in the far field (asserted box 0.02). Recertified with the boxes the verifier actually asserts, at **alpha <= 2.5e-4**, every coded constant is an upper bound [cert]. At the alpha-level 1e-3 four of them fail (up to x1.45). Two hypotheses are therefore used and not asserted: alpha <= 2.5e-4 in the E2 block, and |w|,|w'| <= 1e-4 on [0,8] (holds: 3.24e-5 [cert]).
3. **`r143 = 14.34` is below its certificate 14.343003** (`S/garding_cert.json`). Effect 5e-10 on the rate [cert probe].
4. **A small missing term**: on the cutoff strip `chi N(w)` is bounded as if `chi w^(j) = w~^(j)`. The strip part needs W-norms: `SRW (...)·Dq`, SRW = 4.43e-3 [cert].
5. **Sizes** [cert given the constants; production cells of `runs/P`, B&B tolerance 5e-4]: with every item added the first-cell rate goes +0.03035649 -> **+0.03036965** (+1.3e-5) and the tail anchor -0.00268437 -> **-0.00268121** (rho x0.99882). With all new terms multiplied by 3: +4.4e-5 and rho x0.9969. The degree conditions of `tail_cert.py` still pass.
6. **Not absorbed by existing slack in a way the code can claim.** The new terms are of type A·D4, Q·D4, Fn·D4, not A^2; the E0 slack (PI3) is A^2-type. Numerically the E0 slack is about 1600 times larger at the anchor (E0 audit: -5.2e-3 against +3.2e-6 here), so nothing is at risk; but the coded inequality is not an upper bound until the terms are added.
7. **LAMT is consistent, in the safe direction.** All residual norms use the double -0.2561532224577463; the verifier uses -0.2561532, which is 2.2458e-8 less negative [cert]. The comparison equation stays valid; no change needed (section 6).

## 1. What the E2 inequality must contain (derivation) [hand]

Setting as E0_AUDIT section 1. `D_rho g = rho^-1 (rho g'')'' = g'''' - 2 phi' g''' + (phi'^2 - phi'') g''`; B = ||w~''||, Dq := D4 + 2.2 W3 + 1.2 B.

    (1/2) d/dtau B^2 = < d_tau w~, D_rho w~ >_rho,
    d_tau w~ = L w~ - [L,chi] w + chi Phi_F + chi N(w) + chi M,
    M = -beta (psi3 + 2 alpha h2) - a psi1 - b psi2 + a A1(u) + b B1(u),   A1(u) = -(u - q u'/alpha_g),  B1(u) = u' cos(theta)/alpha_g.

No boundary term at x = 0 (even class: (d_tau w~)'(0) = 0 and (rho w~'')'(0) = 0). `||D_rho w~|| <= Dq` because sup 2|phi'| = 2.16445 <= 2.2 and sup|phi'^2 - phi''| = 1.16067 <= 1.2 on [0, 34.25] [cert, `e2_consts.json` b_Drho].

| # | term | what it needs | coded (`V/rhs_iv.py`) |
|---|---|---|---|
| T1 | `<L w~, D_rho w~> = <(L w~)'', w~''>`, (Lw)'' = Mt w'' + e1 w' + e0 w | Gaarding: -(1/2) D4^2 + sup(-mu2) B^2 + sup\|e1\| W1 B + sup\|e0\| A B = 14.343003, 16.276331, 10.668647 [read `garding_cert.json`] | :249 `r143 BB + r163 W1B + r107 AB`, :258 NEG. **r143 = 14.34 too small** |
| T2 | `-<[L,chi]w, D_rho w~>` | \|\|[L,chi]w\|\|_rho <= d3 V2 + d2 V1 + d1 V0 + d0 \|\|w\|\|_{W,strip}, d_j = sup_strip \|coef_j\|(rho/W)^(1/2) [read `agent_coup/run_coup.py:44,53`] | :247-249 DS·Dq, ok (constants from `constants_v4.json`, rounded up) |
| T3 | `<chi Phi_F, D_rho w~>` | Fn·Dq | :106 (L list), ok |
| T4 | `<chi N(w), D_rho w~>`, N = F(u0+w) - F(u0) - Lw = sum d_tjl w''^t w'''^j w''''^l + d_000 | \|\|chi N\|\|_rho·Dq; contains the alpha-linear part Lin(w), the quasilinear top order (1 - J^-4) w'''' (Cauchy-Schwarz, no derivative of the coefficient), the drift. Coefficient sups on jet boxes; on the strip chi w^(j) != w~^(j) | :244-246 ND·Dq. Constants: section 3. **Strip part missing** |
| T5a | `-a<f1, D w~> - b<f2, D w~> = -a<D f1, w~> - b<D f2, w~>`, lead part (a,b) = y := Minv·l, l_j = <w~, L* f_j> | sup over {f1,f2}^perp of the quadratic form: Q2 <= 43.3271 A^2 (inf -0.189) | :250 `modE2 AA`, 43.33, ok |
| T5b | same with (a,b) - y | \|(a,b) - y\| <= (dab + dlead + MAB alpha) A + QAB Q + ABs + PHIAB Fn, times NF12·Dq (Cauchy-Schwarz) or NDF·A (after integration by parts, w~ in {f1,f2}^perp) | **missing** |
| T5c | `-beta <chi psi3 + 2 alpha chi h2, D_rho w~>` | \|beta\| (n3 + 2 alpha nh2)·Dq, \|beta\| <= (BETA0 + dlead + BETA1 alpha) A + QB Q + BTs + BETAF Fn | **missing** |
| T5d | `a <chi A1(u), D w~> + b <chi B1(u), D w~>` | \|(a,b)\| [alpha NGAB + alpha^2 (RA2 + RB2) + (1 + CA0 + CB0) A + (CA1 + CB1)(W1 + \|\|chi' w\|\|)]·Dq | **missing** |

Notes on structure.
- **Projection.** Only T5a/T5b see it: `<D f_j, w~>` may use the {f1,f2}^perp-projected norm (NDF = 12.341 against 13.278 unprojected [cert]). The beta and A1/B1 terms pair with `D_rho w~`, which is not orthogonal to anything, so no projection and no rank-one symmetrisation is available there; the full norms n3 = ||chi psi3|| = 7.69698, ||h2|| <= 13.1084, ||(gA,gB)||_op = 8.0286 enter [cert].
- **Lead mismatch.** Q2 was certified with the 2x2 block Minv. The modulation lemma's lead is `sum_{j<=3} H_ij <w~, L* f_j>`. The difference is `(H_2x2 - Minv^T) l + H_i3 <w~, L* f3>`, bounded by dab·A with **dab = 0.0131137** [cert; of which the f3 column 78.65, 147.84 times eps1 gives 0.011128]. In the integration-by-parts form this is +0.1589 A^2: needed 43.4889 against 43.33.
- **Strip parts of (a,b), beta** enter in the state form ABs, BTs of v4 (`V/rhs_iv.py:141-142`), so no constant term without a factor alpha appears and the hypothesis sup|w'''| <= 0.02 is not used.
- **Why Cauchy-Schwarz against D_rho w~ is enough.** Every T5b-d term is eps·X·Dq with X of x-degree >= 1; the D4 dissipation absorbs it. The alternative (move D_rho to the other factor) needs new certified norms (||D_rho(chi psi3)||, ...) and gains nothing visible (section 4).

## 2. Audit of the whole E2 block (file:line in `V/`)

| # | finding | size | type |
|---|---|---|---|
| **E2-1** | `rhs_iv.py:250` `modE2·AA` covers T5a only. T5b, T5c, T5d missing. Coefficient of A·Dq they add: 0.1131 at the tail anchor (alpha = 3.68e-5), 0.3380 at alpha = 1.6e-4, 0.5025 at 2.5e-4 [cert], against the coded A·Dq coefficient 0.54 (nd4). Quadratic part 148.4 Q·Dq; cubic w-part 127.7 A·W1·Dq [cert] | rate +1.3e-5 / +3.2e-6 (with E2-2,3,5) | missing terms |
| **E2-2** | `rhs_iv.py:106` E2 forcing `eps·Fn·Dq`: the Fn parts of beta and (a,b) are missing. Needed factor 1 + BETAF (n3 + 2 alpha nh2) + PHIAB (NF12 + alpha(NGAB + ...)) = **13.984** at alpha <= 2.5e-4 [cert] | in the numbers of E2-1 | missing terms |
| **E2-3** | `rhs_iv.py:43` r143 = 14.34 < 14.343002687670673 (`S/garding_cert.json`, region x<=1). r163 = 16.3 >= 16.276331 and r107 = 10.7 >= 10.668647 are rounded up | rate +5.2e-10 / +3.8e-10 [cert probe] | rounded down |
| **E2-4** | `rhs_iv.py:41-42` typed below their own certificate `S/nlin_E2_0.001_0.0001_0.0001.json`: nd2 0.012 < 0.0120544; nd4 0.54 < 0.543880; nd7 0.136 < 0.1360149; nd8 3.05 < 3.052044; nf2 6.9e-4 < 6.910e-4; nf3 2.9e-5 < 2.9062e-5; nf5 0.46 < 0.464015; nf7 8.6e-6 against A1 = 5.0e-3 in that file (a q/4 dependency artefact; the artefact-free value at that box is 1.034e-5 [cert], still above 8.6e-6). Rounded up: nd1, nd3, nd5, nd6, nf1, nf4, nf6 | none at alpha <= 2.5e-4 (section 3) | rounded down; valid by a different box |
| **E2-5** | jet box of the certificate (\|w\|,\|w'\| <= 1e-4 on [0, 34.25], alpha <= 1e-3) is asserted nowhere (`closure.py:72-85`) and is false in the worst case: sup(\|w\|,\|w'\|) <= 3.244e-5 on [0,8], **1.355e-4 on [0,12]**, 7.474e-4 on [0,16.25] at the peak [cert, `runs/e2_probe_tol0.0005.json` box]; far field: asserted box 0.02 | none after section 3 | hypothesis not asserted |
| **E2-6** | the E2 block has no alpha-level switch. `K.v4_level` lets alpha run to 1e-3; there nd7, nd8, nf5, nf6 are too small: needed 0.1453, 3.0519, 0.6688, 3.1201 [cert `nl_regions_a1e-3.json`] | none in production (alpha <= 1.6e-4) | hypothesis not asserted |
| **E2-7** | `rhs_iv.py:244-246`: ND uses \|\|w~^(j)\|\| for `chi w^(j)` on the strip sigma in [29,32]. Correct split: core (w~ = w) + strip, the strip part = SRW (nf1 V3 + nf2 V2 + nf3 V1 + nf7 (V0 + wstrip) + supw2f (nf4 V1 + nf5 V2) + nf6 supw2f^2 V1), SRW = sup_strip (rho/W)^(1/2) <= 4.4292e-3 [cert] | 2e-5·V3·Dq at most; in the numbers of E2-1 | missing term |
| E2-8 | `rhs_iv.py:244-246`: the tip and far coefficients of the same L2 norm are added (nd1 + nf1 etc.); the maximum would do. Conservative; this is what makes E2-4 harmless | slack | info |
| E2-9 | `rhs_iv.py:43` t22 = 2.2, t12 = 1.2: valid (2.16445, 1.16067 [cert]) | none | ok |
| E2-10 | `rhs_iv.py:247` wstrip: sqrt(6 W(31.25)) sup\|w\| >= \|\|w\|\|_{W,strip} (3 W would do); `rhs_iv.py:243` supw2, supw2f cover [0,12] and [12,34.25] (sr uses K_rho(16.25), which covers [0,16.25]) [read] | none | ok |
| E2-11 | DS0..DS3, EPS-weighting, NEG `eps D4^2/2` (:258): consistent with T1, T2 [read]. Not re-derived: the Gaarding identity behind mu2, e1, e0 (`garding_cert.py:1-26`), the Gram behind Q2, the commutator coefficients | - | not audited |
| E2-12 | The coefficient functions are certified for the C^3 functions psi3^, h2^; v4 uses the C^4 ones. My recertification pads the jets by sup\|psi3~ - psi3^\|^(j) (`c4switch_S1.json`, <= 3.9e-6); the change is invisible at the printed digits | none | info |

What E2 does **not** suffer from, against the E0 list: no mistyped norm inside a product (43.33 is a certified eigenvalue bound); the alpha-linear part of N is inside nd/nf (N contains Lin(w), and the drift, unlike E0's T_c forms; `nlin_E2.py:22-24`); the quasilinear top order is there (nd1, nf1).

## 3. The nonlinear constants, recertified [cert: `e2_nlin_regions.py`, `nl_regions_a2.5e-4.json`, `nl_regions_a1e-3.json`]

Same machinery as `S/nlin_E2.py` (read-only import of `h2res`), with per-region jet boxes, the drift part of A0/A1 differentiated analytically ([sym], asserted in the script), and the C^4 pad. A run with the original setting reproduces the far column of the original JSON to all printed digits (`nl_repro_a1e-3_eta1e-4.json`).

alpha <= 2.5e-4:

| region, box (\|w\|, \|w'\|) | w'''' | w''' | w'' | A0 | A1 | w''^2 | w''w''' | w''^3 |
|---|---|---|---|---|---|---|---|---|
| [0,8], 1e-4 | 2.166e-3 | 3.785e-3 | 1.145e-2 | 0.2052 | 0.2295 | 4.6574 | 0.02064 | 3.0163 |
| [0,8], 1e-3 (if the new assertion is not wanted) | 7.788e-3 | 1.359e-2 | 2.783e-2 | 0.4979 | 0.5510 | 4.7433 | 0.02964 | 3.0590 |
| [8,12], 1e-3 | 8.23e-5 | 1.69e-4 | 1.62e-4 | 2.84e-5 | 1.74e-4 | 0.01052 | 0.04373 | 3.0007 |
| [12,16.25], 1e-3 | 6.94e-5 | 2.11e-5 | 3.99e-6 | 4.0e-7 | 2.99e-6 | 0.00324 | 0.05890 | 3.0009 |
| [16.25,34.25], 0.02 | 1.978e-3 | 1.15e-4 | 1.36e-6 | 5.6e-9 | 4.1e-8 | 0.00333 | 0.3150 | 3.0268 |
| **coded** (L2 terms: nd + nf; sup terms: tip / far) | 0.0112 | 0.01269 | 0.04203 | 0.54001 | 0.62001 | 4.75 / 0.014 | 0.136 / 0.46 | 3.05 / 3.06 |

Every coded entry is >= the maximum of its column over the regions it has to cover (L2 terms: all regions; sup terms: tip = first two regions, far = last two). With the box 1e-3 on [0,8] two fail (w''' 0.01359 > 0.01269, w''^3 3.0590 > 3.05), so the box 1e-4 on [0,8] is needed; it holds with ratio 0.324.

alpha <= 1e-3 (same boxes): [0,8] w''^3 3.0519 > 3.05; [8,12] w''w''' 0.1453 > 0.136; [16.25,34.25] w''w''' 0.6688 > 0.46, w''^3 3.1201 > 3.06. So the coded constants are **not** valid at the 1e-3 level.

Hypotheses of the E2 nonlinear constants, to be asserted (worst case over shapes, as `closure.box_assert`):
- (H-E2a) alpha <= 2.5e-4: **new**. Production maximum 1.6e-4.
- (H-E2b) sup(|w|,|w'|) <= 1e-4 on [0,8] (K_rho(8) Sobolev, `KREG[4,8]`): **new**. Peak value 3.244e-5.
- (H-E2c) <= 1e-3 on [8,16.25]: this is B-b (the code's S_ramp is the sup over [0,16.25]). Peak 7.474e-4.
- (H-E2d) <= 0.02 on [16.25,34.25]: B-c, B-d, asserted.

## 4. Sizes [cert given the constants: `e2_probe.py`, `runs/e2_probe_tol0.0005.json`]

Production cells of `S/runs/P/rate_table_fin13v4_a0.00015.json`; `bb.bound_cell` of fin13v4 unmodified; `rhs_iv.py` replaced by `vp/rhs_iv.py` (= fin13v4 plus the switchable E2 terms, built by `mk_vp.py`; with the switch off it reproduces the production U to the last bit).

| variant | first cell (alpha <= 1.6e-4) | tail anchor | rho ratio |
|---|---|---|---|
| coded fin13v4 | +0.0303564931 | -0.0026843659 | 1 |
| r143 = 14.3431 only | +0.0303564937 | -0.0026843655 | 1.0000 |
| **E2 fix, Cauchy-Schwarz form** (section 5) | **+0.0303696537** | **-0.0026812096** | 0.99882 |
| E2 fix, (a,b)-rest by parts (NDF·A; forcing 6.07·Fn·Dq + 47.4·eps·Fn·A) | +0.0303742910 | -0.0026786289 | 0.99786 |
| Cauchy-Schwarz form, all new terms x3 | +0.0304000045 | -0.0026759287 | 0.99686 |

- The Cauchy-Schwarz form is the cheaper one; use it.
- `tail_cert.run` on the patched rhs: radial H (0,2), L (0,1); s H (2,5), L (0,1); joint without floors H (2,5), L (2,5): **PASS** for both forms [cert-mechanical, `runs/e2_tailcert.json`].
- Not run: a full march, the unconditional pass (J constants), `mod_sups`. Expected [est]: J and peak move by about 1e-3 relative (rho x0.9988; the E0 fix F-B moved J0 by x1.0002 for a rate change of 1.9e-6 at the anchor, this one is 3.2e-6).
- Slack: the A·Dq coefficient actually needed at alpha <= 2.5e-4 is 0.2052 (A0) + 0.5025 = 0.708 > 0.54 coded; at alpha <= 1.6e-4 it is about 0.47 [est: A0 not recertified at 1.6e-4], i.e. inside the coded 0.54. The Q·Dq, A·W1·Dq, strip and Fn·Dq parts have no coded counterpart. So: partly absorbed by accident, not as a statement one can write down.

## 5. The fix, for fin13v4b (`step3nl_v4b/verifier/`)

All constants as outward entries of the v4b constants JSON (raw -> rounded up). Names of v4b objects as in its `rhs_iv.py` (`ABx`, `BTx`, `ABs`, `BTs`, `ABc = PHIAB·fn`, `Q`, `S0s`, `supv2`, `XW = XW_COEF·S0s`) and its `e0v4b_constants.json` (CA0, CA1, CB0, CB1, RA2, RB2, GABN, XW_COEF).

**F1. Constants** (new stage or an extension of `e0_lemma_consts.py`; `e2_consts.py` here computes all of them, section d and e of `e2_consts.json`):

| name | raw [cert] | meaning, source |
|---|---|---|
| R143 | 14.343002687670673 | max over regions of `-mu2`, `garding_cert.json`. Replaces 14.34 |
| DAB | 0.013113620982 (+ CAB0_dlead) | \|lemma lead - Minv·l\|/A: sum_j \|H_ij - Minv_ji\| c_j + \|H_i3\| eps1, Euclidean over i = 1,2 (`ibp_consts.py` H, Mi, c_j, eps1) |
| NF12 | 2.0619618907 | lambda_max^(1/2) of the Gram of (f1, f2) (`gram8_g1s4.json` entries 6,7) |
| N3 | 7.6969814468 | \|\|chi psi3\|\|: sqrt of the upper end of g0·g0 in `bloop/gram7_g3.json`, plus 1.94e-10 (C^4) |
| NH2R | 13.1084396532 | \|\|h2\|\|_rho (`h2res.json` h2_norm2_half), plus 7.3e-10 (C^4) |
| SRW | 4.4291959663e-3 | sup over the strip of (rho/W)^(1/2), at x = 31.24 (psi - phi decreasing on the strip [hand]) |

The nd/nf literals may stay (valid under H-E2a..d, section 3); cleaner is to read the column maxima of `nl_regions_a2.5e-4.json` rounded up, with nd := max over [0,8], [8,12] and nf := max over [12,16.25], [16.25,34.25].

**F2. E2 block** (after `V/rhs_iv.py:250`; inside `eps·[...]`, H group):

    RABx = (DAB + CAB0_dlead + MAB·al)·A + QAB·Q + ABs                 # |(a,b) - Minv l|, x-part
    n3e  = N3 + 2·al·NH2R
    wprt = (1 + CA0 + CB0)·A + (CA1 + CB1)·(W1 + XW)                    # ||chi (A1,B1)(u) - chi (A1,B1)(u0)||
    u0pr = al·(GABN + al·(RA2 + RB2))                                   # ||chi (A1,B1)(u0)||
    RM   = RABx·NF12 + BTx·n3e + ABx·(u0pr + wprt) + ABc·wprt
    NS   = SRW·(nf1·V3 + nf2·V2 + nf3·V1 + nf7·(V0 + wstrip) + supw2f·(nf4·V1 + nf5·V2) + nf6·supw2f^2·V1)
    R2p += (RM + NS)·Dq          # and r143 -> R143;  modE2 stays 43.33 (the lead; its mismatch is DAB inside RABx)

**F3. L list** (`V/rhs_iv.py:106`): replace `eps·fn` in the three E2 entries (D4, W3, B) by `eps·fn·FE2`,

    FE2 = 1 + BETAF·n3e + PHIAB·(NF12 + u0pr)        (= 13.984 at al = 2.5e-4; P-arithmetic in al, like np3e)

(the x-degree-0 parts ABc, BTc = BETAF·fn of |(a,b)|, |beta| times Dq have x-degree 1 and carry a factor alpha, so they belong to L; `ABc·wprt·Dq` has x-degree 2 and stays in H, as written in F2). `unc.py` picks the change up through the shared list.

**F4. Assertions** in `closure.box_assert`: `E2_alpha = (aq, 2.5e-4)`; `E2_tip_jet = (sup over shapes of max(s0, s1) with K = KREG[4,8], 1e-4)` (same formula as `sramp_sup`, one more K). Both hold in production with ratios 0.64 and 0.324.

**F5. Checks**: `tail_cert.py` on the new code (passes on my patched copy); regression switch (old E2 block) reproducing fin13v4; a float check of the E2 identity of section 1 on random even w (as `agent_mod/ibp/float_check.py`), because sections 1 and 5 are [hand].

**F6. Text**: one paragraph "E2 couplings" next to the E0 lemma: the table of section 1, the constants of F1, hypotheses H-E2a..d, and the statement that 43.33 bounds the Minv-lead form only.

Cost: 2-3 agent-hours (Opus) on top of the running v4b work, compute under 10 minutes; no new downstream consequence (fin13v4 is not pinned anywhere yet).

Not recommended: absorbing the E2 shortfall into the PI3 slack by a lemma. The terms are not A^2-type, and v4b already codes the analogous E0 pieces explicitly.

## 6. LAMT against the lamt of the forcing residuals [cert + read]

- Residual side: `step3/gen_hs3.py:19` `LAMT = -0.2561532224577463` (a double). It is the only source: `psi3fun.py:14-16` (LAM), `h2fun.py:16-18` (LAM = 2·LAMT), hence `psi3res3.py:16`, `h2res.py:17`, `agent_int/v3/int_res3.py:15` and `agent_int/c4switch/int_res3_base.py:15` (`LAMS = {psi: F3.LAM, h2: H3.LAM}`), which produce FN0, FN1 (`c4switch_S2_1.json`) and FNW*, FNWd* (`int_constants_v3.json`) [read `constants_v4.json` sources]. The cubic constant C3 (`n3_forcing`) contains no lamt. The modulation lemma uses the same double (`agent_mod/ibp/ibp_boxes.py:20`, `lemma/lyap_copy.py:102`, `step3/lyap_consts2.py:102`).
- Verifier side: `V/constants.py:14` `LAMT = -0.2561532`, used only in `V/closure.py:18, 96, 106, 123, 149, 156` (barrier and tail alpha-condition).
- Difference: LAMT_verifier - lamt_residual = **+2.2458e-8** (relative 8.8e-8) [cert, `e2_consts.json` a_LAMT]. The verifier's value is the less negative one.
- Consequence [hand]: beta is defined with the residual's lamt (beta = alpha_tau - lamt·alpha is what the w-equation and the modulation lemma use). Then |alpha|' <= lamt|alpha| + |beta| <= LAMT_verifier|alpha| + |beta|, so the barrier (`closure.py:96-125`) and the tail condition (`closure.py:149, 156`, which only get harder with a smaller |LAMT|) are valid. **No mismatch of the R3 kind.**
- For the record: read the other way (residual taken with respect to the verifier's LAMT), the forcing norm would grow by 2.2458e-8 x 7.697 = 1.73e-7 per unit alpha, 21.8% of FN0 = 7.93e-7. That reading is not the one the proof uses, but the write-up should say which lamt defines beta. A typed LAMT rounded the other way (-0.2561533) would have been invalid.

## 7. What this audit does not establish

- [hand]: the derivation of section 1 and the forms of section 5 (two-line Cauchy-Schwarz steps, but not machine-checked; F5 asks for a float check).
- [read] inputs taken as certified: `gram8_g1s4.json`, `gram7_g3.json`, `mod_forms.json` (Q2), `garding_cert.json`, `coupling_constants.json`, `final_5e-4.json`, `lyap_consts2.json`, `h2res.json`, `c4switch_S1.json`, the profile hulls.
- The v4b constants CA0, CA1, CB0, CB1, RA2, RB2, XW_COEF were read from `step3nl_v4b/v4bcert/out/e0v4b_constants.json` while that work was running (values of Nov 20 00:26); my probe uses its own cruder XW coefficient 5.71e-6. They hold for alpha <= 2.5e-4, |w| <= 0.02.
- Not re-derived: the Gaarding identity (mu2, e1, e0), the level-2 commutator constants DS, the Q2 Gram computation, the K_rho Sobolev constants.
- No full march; J constants and peak for the fixed E2 block are [est] (section 4).
- Step 3 is not closed by this audit.

## 8. Files (all in `writeup/lemmaU/bulk/fix5/e2_audit/`)

- `e2_consts.py` -> `e2_consts.json`, `e2_consts.log`: LAMT, D_rho constants, rounding tables, modulation constants, SRW [cert; arb 200 bits; inputs with sha256].
- `e2_nlin_regions.py`, `run_nl.sh` -> `nl_regions_a2.5e-4.json`, `nl_regions_a1e-3.json` (+ logs), `nl_repro_a1e-3_eta1e-4.json` [cert].
- `mk_vp.py` -> `vp/rhs_iv.py`, `vp/PROVENANCE.txt`: the patched copy of the fin13v4 rhs (switch `E2X`; off = fin13v4).
- `e2_probe.py`, `run_probe.sh` -> `runs/e2_probe_tol0.0005.json`, `runs/probe_5e-4.out`; `e2_tailcert.py` -> `runs/e2_tailcert.json`.
- `SHA256SUMS`.
- Runs: one process at a time, nice 19, PYTHONHASHSEED=0, `python3 -B`, PYTHONDONTWRITEBYTECODE=1, peak RSS 212 MB (nonlinear recertification); about 1 min wall per nonlinear run (3 runs), 2.5 min for the probe, seconds for the rest.
