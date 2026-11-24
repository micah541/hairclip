# E0_AUDIT: the E0 energy block of the step-3 tip-zone verifier (fin12f4 / fin13v4)

Nov 19, 2026. Independent audit (Fable). Scope: `step3nl_v4/verifier/rhs_iv.py:218-241` (E0 block), its constants, and V4_REPORT open item 1.
Proof tree read only. Everything new is in this directory. Labels: **[cert]** arb/interval arithmetic on certified inputs (inputs taken as given), **[read]** read off code or a certified JSON, **[hand]** hand derivation, not machine-checked, **[est]** float estimate.

## 0. Verdict

**FIX (small). Not "E0 OK as coded", not serious.**

1. **The two constants V4 flagged are valid upper bounds, but not for the reason written in ASSEMBLY.md, and V4's "new" values are not the right ones either.**
   - The typed `n3 = 7.52` is not `||chi psi3^||_rho`. The certified Gram gives **7.6970** [cert, `gram7_g3.json`]. With the full lemma bound that would mean Pi3 = 0.01878 and a forcing factor 6.05, worse than V4's probe.
   - But `w~` is orthogonal to f1, f2, f3, so only the projected norm enters: **n3_123 <= 6.5882** [cert]. And the j = 3 part of beta's lead is a rank-one form `<w,u><w,g0>` with `u` almost orthogonal to `g0`, so it costs half of Cauchy-Schwarz.
   - Result: **Pi3 needed = 0.010385 [cert + read]** against 0.01607 coded; **forcing factor needed = 5.334 [cert]** against 5.72 coded.
2. **A real term is missing from E0:** beta's alpha-linear part from `<chi Lin(w), f3>` (the IBP lemma's `m_new` = 4.48). The alpha-form needs **38.3** (37.2 [cert] + 1.1 [est]) where the run uses MODA = 25. The shortfall, 13.3 alpha A^2 <= 2.13e-3 A^2, is covered by the Pi3 slack 5.68e-3 A^2 (x2.67 at alpha <= 1.6e-4) [cert given 1].
3. **Small defects that the A^2 slack cannot absorb** (they are W1-, B-type): two alpha-linear constants typed below their certified values (59.3 < 59.327, 6.2 < 6.21990), and the alpha-linear form is not derived on the strip x in [31.25, 34.25]. Effect of the first: +8.4e-6 on the first-cell rate, +1.9e-6 at the anchor [cert probe].
4. **So fin13v4's E0 inequality is numerically right to about 1e-5 of rate, but it is not a certified upper bound as coded and its written justification is wrong.** The fix is a short E0 modulation lemma plus one verifier pass with corrected constants (section 6, F-B: J0 x1.0002). Nothing downstream is pinned to fin13v4 yet, so the fix costs no step-4 or kit work beyond the re-pin already planned for v4.
5. **The tail is not short.** With correct constants the tail anchor rate is the same or better than coded (section 5). V4's -2.85e-4 came from a bound that is too crude; and even in that scenario the march would simply have gone on to a smaller alpha (-3.97e-3 at alpha = 2e-5).

## 1. What the E0 inequality must contain (derivation)

Setting (`BOOTSTRAP.md:7-13`, `ASSEMBLY.md:36-55`), half line, L2(rho), `w~ = chi w`, `g0 = chi psi3^`, A = ||w~||:

    w_tau = L w + Phi_F - beta psi3^ - 2 alpha beta h2^ - a psi1 - b psi2 + N(w) + a A1(u) + b B1(u),   beta = alpha_tau - lamt alpha.

Pair with rho chi^2 w:

    (1/2) d/dtau A^2 = <L w~, w~> - <[L,chi] w, w~>                         (theta-certificate; strip)
                     - beta <g0, w~> - 2 alpha beta <chi h2^, w~>            (beta coupling)
                     + a <chi A1(u), w~> + b <chi B1(u), w~>                 ((a,b) coupling)
                     + <chi Phi_F, w~> + <chi N(w), w~>.                     (forcing; alin + cubic)

Exact facts [hand]:
- `<chi psi1, w~> = <f1, w~> = 0` and `<chi psi2, w~> = <f2, w~> = 0`. So (a, b) couple only through A1(u), B1(u): at order alpha through `gA = chi(psi - q psi')`, `gB = chi psi' cos(theta)`, and at order w through `-(w - q w'/alpha_g)`, `w' cos/alpha_g`.
- `<g0, w~> = <Pi_perp g0, w~>`, Pi_perp the projection onto {f1, f2, f3}^perp. So every "n3 |beta| A" term may use `n3_123 = ||Pi_perp g0||`.
- Modulation (`MOD_LEMMA.md`): x = (a, b, beta) = H R - H E^T x, R_j = <w~, L* f_j> - <[L,chi] w, f_j> + <chi Phi_F, f_j> + <chi N(w), f_j>. Hence

      beta = beta_3 + beta_12 + beta_Lin + beta_Q + beta_F + beta_S + beta_E,

  | piece | definition | bound |
  |---|---|---|
  | beta_3 | H33 <w~, L* f3> = H33 <w~, u>, u = (L* - lamt) f3 (uses <w~, f3> = 0) | H33 <= 32.149, \|\|u\|\| <= r3 = 5.264e-5 <= eps1 = 6.645e-5 |
  | beta_12 | H31 <w~, L* f1> + H32 <w~, L* f2> | 3.0345e-4 A [cert] |
  | beta_Lin | sum_j H3j <chi Lin(w), f_j>, Lin = (DF(u0) - L) w | m_new alpha A = 4.4819 alpha A [read `final_5e-4.json`] |
  | beta_Q | quadratic and cubic in w, and the w-part of E | QB (A + W1 + B)^2, QB = 5.5264 |
  | beta_F | forcing | BETAF Fn = 0.65653 Fn; without feedback \|\|h\|\| Fn, h = sum_j H3j f_j, \|\|h\|\| = 0.62614 |
  | beta_S | strip | state form (v4, U-3) |
  | beta_E | -sum H3j E_kj x_k, alpha-part | first order: -alpha H33 (e13 a + e23 b), exact in Qfull (BLOOP) |

Four structures decide the constants.

**S1, projection [cert].** `||g0||^2 in [59.2411, 59.2435]`, so `||g0|| = 7.6970` (`agent_mod/bloop/gram7_g3.json`, entry g0-g0). Projecting off f1, f2 removes 13.2887; projecting off f3 removes at least (G33_lo - |b M^-1 c|)^2 / ||f3||^2. Result: n3_12 <= 6.7790, **n3_123 <= 6.5882** (lemma's G33 >= 0.031107) or 6.5805 (G33 >= 0.031717, the independent enclosure in `bloop/pair_f3.json`).

**S2, rank-one symmetrisation [hand + cert numbers].** The E0 term from beta_3 is `-H33 <w~,u> <w~,g0>`. For any u, v: sup over unit w of |<w,u><w,v>| = (||u|| ||v|| + |<u,v>|)/2. Here u is the residual of the adjoint eigenfunction: it lives on the cutoff strip (norm r3_far = 5.135e-5) and on the tip bumps (r3_bump = 1.292e-6) (`step3/lyap_consts2.py:143-163, 223`). So

    |<u, g0>| <= r3_far ||g0||_{L2(rho), strip} + r3_bump ||g0|| <= 2.90e-7 + 9.95e-6 = 1.024e-5   [cert; sup|psi3^| <= 719.62 on the strip],

against ||u|| ||g0|| <= 5.115e-4. (Duality gives the same conclusion: <u, g0> = <f3, (L - lamt)(chi psi3~)> <= ||f3|| FN0 + strip = 5.2e-8 [hand].) The beta_3 term therefore costs (1/2) H33 (eps1 ||g0|| + 1.02e-5) = 0.008386 A^2 instead of 0.016443 A^2. The conclusion is robust: Pi3 <= 0.01607 holds as long as |<u,g0>| <= 3.64e-4, i.e. 71% of ||u|| ||g0||.

**S3, oblique projection [hand].** `<chi Phi_F, w~> - beta_F0 <g0, w~> = <chi Phi_F, (I - h (x) g0) w~>`, and `<g0, h> = (G G^-1)_33 = 1` exactly, so `h (x) g0` is a projection and `||I - h (x) g0|| = ||h|| ||g0||` (no "+1"). The forcing factor is 0.62614 x 7.6970 = 4.819, plus the feedback part n3_123 (BETAF - ||h||) = 0.200: **5.020**. The plain bound with S1 alone is 1 + n3_123 BETAF = **5.325**.

**S4, the Q1 + loop cancellation [cert, BLOOP].** `Qfull <= 6.1942 alpha A^2` on {f1,f2,f3}^perp covers a <chi A1(u0), w~> + b <chi B1(u0), w~> at first order together with beta_E at first order. It does **not** contain beta_Lin (`bloop_forms.py:1-5`: the loop is built from e13, e23 only; `BLOOP.md` caveats).

## 2. Which constants are correct [cert: `e0_cert.py`, `e0_cert.json`]

Pi3 (coefficient of A^2 from beta_3 + beta_12 + dlead), coded 0.01607 (`constants.py:14`, `rhs_iv.py:39`, `unc.py:15`, `rhs_iv.py:364`):

| treatment | value | <= 0.01607 ? |
|---|---|---|
| ASSEMBLY.md: n3 eps1 / G33 with n3 = 7.52 (diagonal only) | 0.016065 | this is what was typed; it drops beta_12 and uses a wrong n3 |
| T0, V4 probe: 7.52 x BETA0 | 0.018347 | no |
| T1, true norm: 7.6970 x BETA0 | 0.018779 | no |
| T2, S1 only: n3_123 x BETA0, lemma's G33 | 0.0160737 | **no, by 3.7e-6** |
| T2, S1 only, G33 >= 0.031717 (both enclosures intersected) | 0.015746 | yes (2.0%) |
| **T3, S1 + S2**: n3_123 (beta_12 + dlead) + (1/2) H33 (eps1 n3_true + X) | **0.010385** | yes (x1.55) |
| T3 with G33 intersected | 0.010183 | yes |
| T3 with \|\|u\|\| <= r3 instead of eps1 | 0.008508 | yes |

Forcing factor, coded 5.72 (`rhs_iv.py:40` literal, `unc.py:22, 37` literals):

| treatment | value |
|---|---|
| T0, V4 probe: 1 + 7.52 BETAF | 5.937 |
| T1, true norm | 6.053 |
| **T2, S1: 1 + n3_123 BETAF** | **5.325** |
| T3, S3 | 5.020 |
| uncoded extras at alpha <= 1.6e-4: (a,b)-forcing x alpha \|\|(gA,gB)\|\|, and -2 alpha beta_F <chi h2^, w~> | + 0.0077 |
| T2 + extras | 5.3331, quoted as 5.334 <= 5.72 |

So: **5.72 is a valid bound with 7% room; 0.01607 is a valid bound only through S2 (or, by 2%, through S1 with the tighter G33).** V4's 0.01835 / 5.937 are neither the documented nor the correct values: they mix the typed 7.52 with the full lemma bound.

The other "7.52" literals (`rhs_iv.py:50` qb752, `rhs_iv.py:56, 184` b752) are valid for the same reason (7.52 >= 6.5882), not because 7.52 is the norm.

## 3. Audit of the whole E0 block (file:line in `step3nl_v4/verifier/` unless noted)

| # | finding | size | type |
|---|---|---|---|
| E0-1 | `constants.py:14` PI3: number valid (section 2), derivation in `ASSEMBLY.md:49` wrong (n3 mistyped, beta_12 dropped, S1/S2 not stated) | none after lemma | doc + lemma |
| E0-2 | `rhs_iv.py:40` f572 = 5.72: valid, derivation wrong; hard-coded again in `unc.py:22, 37` | none | doc |
| **E0-3** | **beta_Lin missing.** `rhs_iv.py:221` mod = MOD_ALPHA alpha A^2, MODA = 25 (`run_v4.sh:7`). 25 was chosen as "Q1 + loop without cancellation" (`STEP3_REVIEW2.md:22-25`); it was never meant to hold beta_Lin. Needed: Qfull 6.194 + n3_123 m_new = 29.527 [cert] | +10.7 alpha A^2 over the coded 25 | missing term |
| E0-4 | second-order and small first-order pieces, none coded: E13/E23 remainders 0.040; (a,b) alpha-corrections in Q1 0.799 and in the loop 0.561; -2 alpha beta <chi h2^, w~> 0.064; E33 feedback 0.017 (coefficients of alpha A^2 at alpha <= 1.6e-4) [cert] | +1.48 alpha A^2 | missing terms |
| E0-5 | **nonlinear drift absent.** The E0 forms alin/cubic are built from the curvature part T_c only (`nlin_E0.py:1-4`, `nlin_E0_check.py:39-42`). The drift (1/4)(p + u - u' q/alpha_g) has the nonlinear part -(q/4) u' k u/(1 - k u). Its alpha-form is (alpha/8) int rho w^2 [((qk)' - phi' qk) psi - qk psi'], sup of the bracket/8 = 0.549 [est, `e0_drift_est.py`]; cubic part <= 0.37 sup\|w\| A^2 [est]. Same structure in a A1(u): \|a\| 15.7 A^2, \|b\| 0.81 A^2 [est], i.e. <= 3.4e-6 A^2 with sup\|(a,b)\| = 2.14e-7 | +1.1 alpha A^2 (allowance, x2 on the float sup) | missing terms |
| E0-6 | `rhs_iv.py:237` strip: 2 IA K_rho(16.25) (A^2 + 2 A W1). The pointwise bound (I2) of `agent_sob/SOB_NOTES.md:11-14` is K (Phi1 \|\|g\|\|^2 + 2 \|\|g\|\| \|\|g'\|\|) with Phi1 = 1.083 | +4.8e-8 A^2 [cert] | missing factor |
| **E0-7** | `constants.py:16` ALIN: sq = 59.3 < 59.32715 and BB = 6.2 < 6.21990 (`alin_E0_spec_0.001.json`, region x<=2). AB = 35.9 and W1B = 7.9 are rounded up | first cell +8.4e-6, anchor +1.9e-6 [cert probe] | rounded down |
| E0-8 | alin covers x <= 31.25 only (`alin_E0_spec.py:6-7`). On the strip [31.25, 34.25] the alpha-linear form and the chi'-terms of the nonlinear pairing are not derived; only the cubic part has a strip term (`constants.py:38`, 2.6e-3). If the far coefficient continues like psi', W1B would be about 8.7 > 7.9 [est] | W1 B-type, unquantified | not derived |
| E0-9 | `constants.py:15, 72` MOD_ALPHA literal 6.19 < certified 6.194179; superseded by MODA = 25 in every run | none | rounded down (dead) |
| E0-10 | `rhs_iv.py:50, 56, 184` literal 7.52 called "n3 = \|\|chi psi3^\|\|": valid as a bound of n3_123 only | none | doc |
| E0-11 | H31, H32 parts of beta_E (about 3e-3 alpha A^2), not in Qfull | negligible | missing |

Sum of the A^2-type shortfalls (E0-3, 4, 5) at alpha <= 1.6e-4: needed alpha-coefficient **38.30** (37.204 [cert] + 1.1 [est]) against 25, i.e. **2.13e-3 A^2**; plus E0-5/E0-6 constants about 3.5e-6 A^2.
Pi3 slack from section 2: 0.01607 - 0.010385 = **5.685e-3 A^2**. Ratio **x2.67** at alpha <= 1.6e-4 (x1.61 at 2.5e-4) [cert given the reads; the drift part is est]. The run never exceeds alpha = 1.6e-4 (`runs/P/rate_table_fin13v4_a0.00015.json`, step 0 alpha_up).

So the A^2 budget of E0 closes with room. E0-7 and E0-8 are not A^2 terms and need corrected constants.

Rounding of the rest [read]: lin1 = 1/32 - PI3 is exact arb (`rhs_iv.py:39`); 1/9 and the cubic coefficients are widened outward (`rhs_iv.py:37, 40`); SB, SA, IA, KR16, QB come from `constants_v4.json` (outward). No other down-rounding found in the E0 lines.

Not audited: the theta-certificate (1/32, 1/9) for <L w~, w~>; the coefficients in `cubic_E0fine_0.001_0.01_0.01.json`; ADD_STRIP_CUBIC = 2.6e-3; the E2 block (its modulation term 43.33 A^2 has the same history as Pi3 and deserves the same check: no beta_Lin, no forcing or strip part).

## 4. Is the V4 probe right?

- Mechanically yes [read `probe_E0/e0probe.py`]: it patches K.PI3 and cl['f572'], resets the constants cache, and calls the production `bb.bound_cell` on the anchor cell and the first cell. I reproduce its numbers exactly (`runs/e0_probe_tol0.0005.json`: -2.850049e-4, +3.340046e-2). They are certified cell bounds for the patched constants, not estimates.
- Its premise is wrong (section 2).
- Its reading is too pessimistic in two ways:
  1. "-2.85e-4 against the -0.002 the tail needs" is not a failure. `closure.py:153` only declines the anchor and marches on; -0.002 is the design threshold rho_min, not a lemma requirement. With the probe's constants the anchor rate is -3.97e-3 at alpha = 2e-5 and -6.13e-3 at 1e-5.
  2. U is the B&B upper bound; U - lb = 5e-4 exactly (the tolerance). At tolerance 1e-4 the coded anchor bound is -3.084e-3 instead of -2.684e-3 (`runs/probe_1e-4_partial_stopped.out`). A tighter anchor tolerance is a free lever worth 15% of rho.

## 5. How short is the tail? [cert given the constants; `e0_probe.py`, production tolerance 5e-4]

Certified upper bounds U of the rate. Cells: first cell (alpha <= 1.6e-4, floors on); tail anchor a* = 3.6846e-5, L* = 2.80672e-15 (floors off); the same anchor at smaller alpha.

| variant | Pi3 | f | MODA | first cell | anchor | a = 2e-5 | a = 1e-5 | a = 1e-6 |
|---|---|---|---|---|---|---|---|---|
| coded fin13v4 | 0.01607 | 5.72 | 25 | +0.030356 | -0.002684 | -0.006314 | -0.008437 | -0.010337 |
| V4 probe "E0fix" | 0.018347 | 5.937 | 25 | +0.033400 | -0.000285 | -0.003972 | -0.006127 | -0.008056 |
| coded + ALIN rounded up (E0-7) | 0.01607 | 5.72 | 25 | +0.030365 | -0.002683 | -0.006313 | -0.008437 | -0.010337 |
| T1: true norm, nothing else (worst case) | 0.018779 | 6.053 | 43.3 | +0.037179 | +0.000888 | -0.003138 | -0.005495 | -0.007604 |
| T2: S1 only, no symmetrisation, G33 intersected | 0.015746 | 5.72 | 38.4 | +0.032185 | -0.002512 | -0.006369 | -0.008627 | -0.010647 |
| **T3: correct constants** | 0.0104 | 5.72 | 38.4 | +0.026839 | **-0.007859** | -0.011715 | -0.013973 | -0.015994 |
| T3 with MODA = 40 | 0.0104 | 5.72 | 40 | +0.027095 | -0.007800 | -0.011683 | -0.013957 | -0.015992 |

(The last four rows also have the ALIN constants rounded up.)

Reading:
- With the correct constants the tail is **not short**: the anchor rate is -7.86e-3, i.e. rho/rho_min = 3.9 instead of 1.34, and the first-cell rate is lower than coded.
- Even without S2 (row T2) the same anchor still passes (-2.512e-3).
- Only if neither S1 nor S2 is accepted (row T1) does this anchor fail; the tail then starts near alpha = 2e-5.

Full marches at alpha1 = 1.5e-4 (the production driver `closure.py`, unmodified, constants patched in memory by `e0_closure.py`):

| run | Pi3 / MODA | steps | peak sqrt(Lambda) | tail a* | rho | sqrt(L_F) | J0 | J1 | J2 | J4 | J0 ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| fin13v4 production (`runs/P`) | 0.01607 / 25 | 10 | 5.307698e-08 | 3.6846e-05 | 2.684366e-03 | 9.8876e-10 | 6.020438e-12 | 8.461143e-12 | 1.872555e-12 | 6.242096e-08 | 1 |
| `c_conservative`: coded + ALIN rounded up (F-B) | 0.01607 / 25 | 10 | 5.307835e-08 | 3.6846e-05 | 2.682591e-03 | 9.8909e-10 | 6.021614e-12 | 8.462797e-12 | 1.872921e-12 | 6.243316e-08 | 1.0002 |
| `c_nosym`: S1 only, no symmetrisation (fallback) | 0.015746 / 38.4 | 10 | 5.333810e-08 | 3.6846e-05 | 2.527654e-03 | 1.0190e-09 | 6.117858e-12 | 8.598058e-12 | 1.902856e-12 | 6.343103e-08 | 1.0162 |
| `c_sharp`: correct constants, anchor rule unchanged | 0.0104 / 40 | 7 | 5.222690e-08 | 5.8425e-05 | 2.584374e-03 | 1.0077e-09 | 7.293446e-12 | 1.025023e-11 | 2.268503e-12 | 7.561973e-08 | 1.2114 |
| `c_sharp_rho75`: correct constants, anchor at rate <= -0.0075 (F-A) | 0.0104 / 40 | 10 | 5.222690e-08 | 3.6846e-05 | 7.724905e-03 | 5.8286e-10 | 4.744401e-12 | 6.667798e-12 | 1.475666e-12 | 4.919079e-08 | 0.7880 |

All five close: tail found, alpha-condition OK, V1 delta_t = 0.2, box hypotheses asserted at every step (worst ratios unchanged to 3 digits: |w| on [34.25,60] 0.789, strip S2 0.76-0.78, S0V 0.765, IR jet 0.74-0.75). f572 = 5.72 in all runs. These are certified marches for the stated constants [cert given the constants]; they are audit runs, not production records (alpha1 = 1.5e-4 only, single pass, no mod_sups).

Reading:
- F-B moves nothing beyond the 4th digit (J0 x1.0002, peak x1.00003, rho x0.9993).
- The fallback without S2 still closes with the same anchor: peak +0.49%, J0 +1.6%, rho 2.53e-3.
- With the correct constants the peak falls 1.6%. The tail rule "first anchor with rate <= -0.002" then fires three steps early (a* = 5.84e-5), and that makes the J's 21% **worse**. With the anchor kept at the production a* (threshold -0.0075) rho is 7.72e-3 (x2.88), the floor sqrt(L_F) falls to 5.83e-10 (x0.59) and the J's fall 21%.

## 6. Fix options and cost

**F-B (recommended): minimal. Keep PI3 = 0.01607, MODA = 25, f572 = 5.72; fix E0-6, E0-7, E0-8; add the E0 modulation lemma.**
- Lemma (text, about 2 pages): the decomposition of beta and S1-S4 of section 1, with the statement
  `Pi3_true + (MOD_true(alpha) - 25) alpha <= 0.01607 for alpha <= 1.6e-4` (now x2.67) and `f_true <= 5.72` (now 5.334).
- Certification: promote `e0_cert.py` to a v4cert stage (seconds). Replace the two [est] items by arb passes over the existing boxes, same machinery as `alin_E0_spec.py`: the drift alpha-form (E0-5) and the alin form on the strip (E0-8); minutes each. Add a float check of the beta decomposition on random w (as `agent_mod/ibp/float_check.py`).
- Verifier: ALIN sq 59.33, BB 6.22 (or the raw values rounded up), W1B from the extended alin, Phi1 in the strip term. New tag (fin13v4 is closed as a record; new files). One double run, about 2 x 12 min.
- Effect: J0 x1.0002, rho x0.9993, peak x1.00003 [cert, `c_conservative`; E0-8 not included, it will move W1B by about 10% [est], a term that is small at these alpha].
- Downstream: nothing extra. fin13v4 is not yet pinned by the (A) chain, step 4 or the kits; the re-pin already planned for v4 (A chain v36 with the R3 patches, manifest_v33, kits) takes the new tag's numbers instead.
- Cost: 3-5 agent-hours (Opus for text and code; Fable only to read the lemma), compute under 1 hour.

**F-A (optional, same lemma): use the sharp constants.** PI3 = 0.0104, MODA = certified need (38.3, typed 40), and the tail threshold raised so the anchor stays at the production a* (one default argument, `closure.py:143`).
- Effect: rho x2.88, floor x0.59, J's x0.79, peak x0.984 [cert, `c_sharp_rho75`].
- Without the threshold change the J's get 21% worse (`c_sharp`), so F-A is only worth doing together with it.
- Every coded constant is then the bound of its own term, which is cleaner to write up. The price is that the slack now held in reserve for the estimated items (drift, strip) is spent, so those must be certified first.
- This is a design choice for Micah; nothing forces it. The (A) chain uses kappa_J = 0.04235 against a limit of 0.13, so it does not need the smaller J's.

**If S2 is not wanted at all:** PI3 = 0.015746 (S1 with the intersected G33), MODA = 38.4. Closes with the same anchor, J0 +1.6%, rho 2.53e-3 [cert, `c_nosym`]. Thin (rho/rho_min 1.26) but sufficient.

**Not needed:** moving the tail anchor to smaller alpha (L5); the sharper state-form strip (already in v4, and the strip part is 1e-9); using v4's smaller sup|beta| (it is an output, and E0 uses the shape form, not the sup); changing the E0 weights.

**Free lever, any option:** anchor cell at B&B tolerance 1e-4 (worth +4e-4 of rho; three times the cell time).

## 7. What this audit does not establish

- [read] inputs taken as certified: `gram7_g3.json`, `pair_f3.json`, `lyap_consts2.json` (eps1, r3_far, r3_bump, and that Q2/Q4/Q5 are the norms of the strip and bump parts of (L* - lamt) f3), `final_5e-4.json`, `bloop.json`.
- [hand] the decomposition of beta and S1-S4 are hand derivations; S2 and S3 are two-line Hilbert-space facts, but the identification of each piece with a coded term should be machine-checked when the lemma is written (a float check on random w, as `ibp/float_check.py` does, is enough).
- [est] the drift terms (E0-5) and the strip alin (E0-8).
- The Gram is for psi3^ (C^3 version); the C^4 psi3~ changes it by about 1e-6 relative (`CERT_PLAN.md`, Nov 9 "C^4 switch"); not re-derived here.
- Step 3 is not closed by this audit. Open besides E0: V4_REPORT items 2-7, the E2 modulation term, the written E_I lemma, the flux lemma.

## 8. Files

- `e0_cert.py` -> `e0_cert.json`: certified constants (arb, 200 bits), with input paths and sha256.
- `e0_drift_est.py` -> `e0_drift_est.json`: [est] drift terms.
- `e0_probe.py` -> `runs/e0_probe_tol0.0005.json`, `runs/probe_5e-4.out`; `runs/probe_1e-4_partial_stopped.out` (3 variants at tolerance 1e-4; stopped by me to save time).
- `e0_closure.py`, `run_e0.sh` -> `runs/c_conservative/`, `runs/c_sharp/`, `runs/c_nosym/`; `e0_closure_rho.py`, `run_e0b.sh` -> `runs/c_sharp_rho75/` (logs and rate tables, alpha1 = 1.5e-4). The driver `closure.py` is executed unmodified; constants (and, for rho75, the default `rho_min`) are changed in memory only.
- Proof tree check after the runs: no file under `step3nl/`, `step3/`, `step3nl_v33/`, `step3nl_v4/`, `s1_tight/` is newer than the start of this audit; `step3nl_v4/SHA256SUMS` verifies; no `__pycache__` written.
- Runs: one process at a time, nice 19, PYTHONHASHSEED=0, `python3 -B` and PYTHONDONTWRITEBYTECODE=1, peak RSS under 100 MB. Another agent's closure run on a copy (`fix5/v4_review/copy/`) was active at the same time; it is not mine.
