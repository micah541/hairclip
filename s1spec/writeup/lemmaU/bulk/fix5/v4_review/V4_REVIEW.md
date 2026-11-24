# V4_REVIEW: independent review of verifier v4 (tag fin13v4)

Nov 19, 2026. Fable (review). Object: `/home/claude/s1spec/step3nl_v4/` against `fix5/ei_spec/EI_SPEC.md`.
The proof tree was read-only; everything I wrote is in this directory. The E0 beta constants are out of scope (`fix5/e0_audit/`).

Labels: **[cert]** arb/interval result (mine unless a file is named); **[read]** verified by reading the cited file:line; **[hand]** my derivation, not machine-checked; **[est]** estimate, not a bound.
Paths: `V/` = `step3nl_v4/verifier/`, `C/` = `step3nl_v4/v4cert/`, `M/` = `step3nl_v4/agent_mod/ibp/`.

## 0. Verdict

**v4's E_I part is sound, modulo the items listed in §5 and §6.**

- Code against spec: every term of EI_SPEC §5.1 and §5.2 is in `V/rhs_iv.py:110-193` with the right coefficient, norm and sign. I found no mismatch (§1).
- Constants: no wrong-side rounding and no typed certified literal in the active E_I path. My independent arb recomputation of about 60 constants found none below its true value (§2).
- Spec derivations marked [hand]: I re-derived them and found no error (§3).
- V4_REPORT claims: regression, double run, tail_cert, float check and box assertions all hold; my own full rerun from a copy reproduces the production outputs exactly (§4).

**Step 3 is not closed by this.** What v4 certifies is: "the inequality of EI_SPEC §5 closes at both points, given hypotheses B-a to B-f, the strip hypothesis on w''', the flux hypothesis (I3) for the redefined K^OR, the legacy E/E2 block, and the E0 constants under separate audit". The open items are in §6.

Top findings (details in §5):
1. **R-1 (open, medium).** FLUX0/FLUX2 were certified (in part) for the old outer-ramp functional. EI_SPEC §4.10 redefines K0^OR, K2^OR; nothing bounds the new ones, and the g''' parts are known not to be pointwise-bounded.
2. **R-2 (open, low in size).** sup|w'''| <= 0.02 on the strip is still a hypothesis in the alpha-equation and in `mod_sups`. Its weight is 4.9% of STRIP_BETA [cert], 0.17% of sup|beta|. One integration by parts in the commutator pairing would remove it.
3. **R-3 (formal, low).** The tail alpha-condition is checked only at the anchor (`V/closure.py:157`), but the constant STRIP_BETA/sqrt(l) grows as l falls. At the floor it is 0.081 against a margin of 177 [hand from cert log]: true, but no script checks it.
4. **R-4 (dependency).** rho, the peak and the J constants all depend on the E0 block (PI3 enters n_A in `V/unc.py:15` and NEG). They are final only after the E0 audit.
5. Records: three small inaccuracies in V4_REPORT (R-6, R-7, R-8). None changes a number.

## 1. Code against spec, term by term

All in `V/rhs_iv.py`, function `_ei_v4` (lines 110-193) and `L_list` (95-109). "ok" = coefficient, norm and sign agree with the spec; the bound itself is re-derived in §3 where the spec marks it [hand].

### 1.1 Derived sups (EI_SPEC §1.2)
| Quantity | Code line | Check |
|---|---|---|
| S1r, S2r | 122 | ok. Regions [8,12] and [12,16.25]; the second uses the arb 1/rho(16.25) (`V/constants.py:80, 101-102`). K_rho(12) also covers [0,12], so B-b covers the whole core [hand] |
| sZ = max(sI2r, S1V/32) | 123 | ok. sigma(34.25) = 32.0009 >= 32 [hand, TIP = -2.24908] |
| sOR = S1V/54.75 + (35/48) S0V | 124 | ok. sigma(57) = 54.7509 >= 54.75; sup|zeta'| on OR = 35/48 exact [cert, sympy in `indep_consts.py`] |
| sg1, sG | 125-126 | ok. sG contains sOR and the IR value S2r + (35/96) S1r |
| S1z, S2z | 127-128 | ok. Left-end Sobolev with KWIR = (1+16.25^2)^2 = 70258.129 [cert]; needs psi' >= 0: inf psi' >= 0.1303 on [10.25,60] [cert] |
| S0s, supv2, S_all34, DV | 130-134 | ok. The strip is x in [31.249, 34.249], inside [16.25, 34.25] [read `agent_sob/sob_constants_P2.json` x_strip_enclosure] |

### 1.2 Level 0 (EI_SPEC §5.1), lines 150-160
| Spec term | Line | Check |
|---|---|---|
| 4.2 AF W1^2 + BF B^2 | 150 | ok |
| 4.3 CBAR V2^2 + (C1 + 2 PSI1 CBAR) V1V2 + (PSI1 C1 + PSI2 CBAR) V0V2 | 151, 136 | ok. CBAR, CA, CB are the max over IR, ZA, ZB; C1 uses sg1 (not sG), as the spec says |
| 4.4 V0 RQ | 152, 137 | ok. IR coefficients; k4 with (4,6,4,1), k3 with (3,3,1) |
| 4.5 Z1 block, incl. hA0 sw16 + HA0 V0 | 153-154 | ok. Far coefficients = max(ZA, ZB); KAP_Z1 |
| 4.6 IR block, incl. the three nonlinear monomials | 155-156 | ok. Table row eta = 1e-3 as a constant (no chord) |
| 4.7 V0 C0W | 157, 138 | ok |
| 4.9 T8, T9, T10 | 158 | ok. T9 carries QUF = 1/(1 - sup|kU|) for Q_u, C_u. ABx part of T10 here, ABc part in L |
| 4.9 T11 beta and psi1, psi2 | 159 | ok. BTx (with the state-form strip BTs) here, BTc in L |
| 4.9 T12 | 160, 148 | ok |

### 1.3 Level 2 (EI_SPEC §5.2), lines 162-182
| Spec term | Line | Check |
|---|---|---|
| 4.1 E1 V1V2 + E0 V0V2 | 172 | ok |
| 4.2 CF block times DV | 173 | ok |
| 4.3 CBAR V4 DV + C1 V3 DV | 174 | ok |
| 4.4 RQ DV | 175 | ok |
| 4.5 table, 12 rows | 162-166 | ok, row by row. 2 d200 sZ V2 (index fixed), 3 d300 (factor fixed), d110 sG (3V3 + PSI1 V2), all d^a, d^b rows, A1 V1, B1 V0 + hB0 sw16 + HB0 V0 |
| 4.6 RWR R5 DV | 167-171, 177 | ok. d110 (S1z D4 + S2z W3); A1 B + B1 W1 + B0 A |
| 4.7 C0W DV | 178 | ok |
| 4.9 T8 (T2, Q2B) | 179 | ok. c62 is gone |
| 4.9 T9 (TC2), T10, T11, T12 | 180-182 | ok |

### 1.4 L group and floors
- `L_list` (95-109): FNW(al)(V1 + KAP V0); FNWd(al) DV; ABc al NAB, BTc (NP3 + 2 al NH2), ABc NP12 at both levels; the E0/E2 forcing in the fin12f4 form. ok. Every coefficient carries a factor alpha.
- FLUX0, FLUX2 are added at line 192 only when floors are on. ok.
- `V/unc.py:16-31` (`lk_v4`) linearises the same list; W3, V1, V3 by AM-GM with the right constants (c1, C1W). ok. SF-3 is fixed.
- NEG (`rhs_iv.py:258`) is unchanged in form; MU1, MU2 now come from the JSON, rounded down.

### 1.5 Things I looked for and did not find
- A term in H with x-valuation below 2, or in L without a factor alpha: none (confirmed by `tail_cert`, §4).
- A use of the max-with-S0V sups where the raw sups are needed (BTs, ABs): none; lines 141-142 use sI1r, sI2r.
- A subtraction or a constant that could be negative in POS: none.

## 2. Constants

### 2.1 Provenance and rounding [read]
- `V/constants.py:85-121`: with NEW=1 every E_I constant is read from `C/out/constants_v4.json`. The file has 269 entries, each with raw, value, dir, source and key; `C/c5_constants.py:93-95` asserts value >= raw (up) or <= raw (down).
- MU1, MU2 are the only "down" entries (`c5_constants.py:42-44`). ok.
- Typed numbers left in the active E_I path: 54.75, 32, 25.75 (geometry, checked above), 35/96 and 35/48 (exact), and the inputs S0V, S1V, FLUX0, FLUX2. No certified literal.
- Legacy typed literals remain in the E0/E2 blocks and the alpha-equation (LAMT, PI3, MODA, MOD_E2, ALIN, nd*, nf*, 5.72, 7.52, 2.2, 1.2, r143, r163, r107). They are listed honestly in `constants_v4.json` `legacy_typed`. See §6.
- "Raw" is not always a raw ball. The entries from `agent_sob/sob_constants_P2.json` and `agent_int/v3/int_constants_v3.json` are 4-digit summaries. Both producers round up (`sob_constants_P2.py:57-60`, `collate_v3.py:5`), so they are valid; I recomputed the ones I could (below).
- The modulation-lemma constants (BETA1, BETAF, MAB, PHIAB, dlead, gamma_ref) come from plain float arithmetic in `M/ibp_final.py:53-67`. The feedback parts all carry fac = 1/(1 - nu), and nu is padded by (1 + 1e-12) at line 43, which covers the float errors by a factor of about 50 [hand]. Acceptable; a note in the write-up would be enough.

### 2.2 Independent recomputation [cert, `indep_consts.py` -> `indep_consts.json`]
My own arb code (200 bits), written from the spec definitions. "lower"/"upper" are my rigorous enclosures of the true sup. Every v4 value is >= my lower bound. Where v4 is below my upper bound, my cover is the coarser one (ratio >= 0.999).

| Constant | v4 | mine (lower, upper) |
|---|---|---|
| PSI1 | 1.3048713 | 1.3048230, 1.3048834 |
| PSI2 | 1.6954377 | 1.6953098, 1.6954698 |
| RWR | 0.6557918 | 0.6555644, 0.6558321 |
| KWIR | 70258.12890625 | exact |
| KAP_IR / Z1 / OR | 0.533967 / 1.282740 / 1.588324 | 0.533066..0.533960 / 1.282690..1.282752 / 1.582225..1.589610 |
| k4_1..4 | 0.442904, 0.526704, 0.424151, 0.425050 | lower 0.441868, 0.525585, 0.423252, 0.424266 |
| k3_1..3 | 0.332178, 0.263352, 0.106038 | lower 0.331401, 0.262792, 0.105813 |
| Z2i | 0.211321 | 0.208700, 0.209077 |
| KW = 1/W(34.25) | 2.323e6 | 2322472.74 (both) |
| Phi(34.25) | 9854657.57 | 9850950.9, 9858367.4 |
| K_rho(1,2,4,8,12) | 1.488, 2.078, 5.222, 56.92, 992.4 | 1.48713, 2.07717, 5.22196, 56.9150, 992.3207 |
| K_rho(16.25) | 30208.006139451245 | 30208.00613945124 (the old typed 3.0208e4 was indeed too small) |
| T0, T2 | 37.184487, 39.184201 | 37.182004..37.184836, 39.181716..39.184550 (profile bounds at sigma = 8 taken from `C/out/c1_cut_sups.json`) |
| SF_beta_0..3, SF_beta_al | 1.6159e-9, 1.3703e-9, 8.1316e-10, 1.9778e-10, 4.5386e-8 | identical to the last digit (my own code from `mod_ibp.json`) |
| SF_ab_0..3, SF_ab_al | 8.7823e-9, 7.4396e-9, 4.4154e-9, 1.0795e-9, 2.4396e-7 | identical |
| STRIP_BETA | 8.0396e-11 | identical |

Top-order coefficients of C2, against the closed forms c = 1 - J^-4, c_U1 = 4 U1 J^-6, d110 = 10 U1 J^-6, d300 = 3 J^-6 (1 - 6 U1^2 J^-2) (I13, [sym] in `C/out/s1_identities.json`), with sup|U1| = sup|u0'| + box:
- level 2.5e-4: ZB c 3.3655e-3 against my 3.3430e-3 (x1.007); c_b 0.16512 against 0.16374; d110 0.41142 against 0.40935; d300 3.0456 against 3.0000. ZA and IR likewise within x1.01, v4 above mine.
- IR c at 2.5e-4: v4 6.92e-5 is below my crude 7.71e-5. Mine adds the worst U1 and the worst kU over the whole zone; v4's per-box value is the correct one (2 U1^2 = 6.92e-5 at |U1| = 5.88e-3) [hand].
- c_a equals c_b * sup|u0''| to 9 digits on ZA and ZB, as it must when k is negligible [cert arithmetic on the table].

MU1 = 0.155523: I read `agent_int/v2/int_garding.py` and `v2_T4_P2_0.25.json`. All five regions of [10.25, 60] are covered; the minimum is on [10.25, 16.25]. My float evaluation of the same form at x = 10.25 gives about 0.1559 [est], consistent with both 0.155523 and the older 0.155167. The switch is legitimate.

## 3. Spec soundness: the hand-derived items

### 3.1 Level-2 inner-ramp reformulation: correct [hand]
- zeta (c g''')' = (c v''')' - c K4 - c' K3: from v''' = zeta g''' + K3 and K3' + zeta' g''' = K4 (binomials (4,6,4,1)). Checked by hand; also I2 [sym].
- (c v''')' = c v'''' + c' v''' is bounded (v in C^3, v'''' in L^inf), and the level-2 pairing uses no integration by parts on it: <D2v, (c v''')'>_W <= DV (CBAR V4 + sup|c'| V3). Sound.
- Level 0, identity I1: I re-derived it; the coefficients are c, c' - 2 psi' c and -((psi' c)' - psi'^2 c). No boundary term (v = v' = 0 at 10.25 and 60).
- |c'| <= CA + CB |g'| with U1' = a2 + g': exact chain rule, and sup|g'| = sg1 covers IR (S2r), [16.25, 34.25] (sI2r) and [34.25, 60] (S1V/32).
- What is left on IR is of order <= w'''' (D4). w^(5) never appears. The answer to GI 2-10 is right.
- S1z, S2z: zeta g' = v' - zeta' g and zeta g'' = v'' - 2 zeta' g' - zeta'' g; v' and v'' vanish at 10.25; the left-end Sobolev needs only psi' >= 0. Sound.
- I8 (||v''^2||_W <= sup|v'| (3 V3 + PSI1 V2)): re-derived; the sup must be global, and the code uses sG.

### 3.2 N_low' and D(d000): correct [hand + read]
- `agent_int/planB/pb_nlin.py:16-24`: Flow depends on (U0, U1) only through the symbols U0, U1, so differentiating H = D G in U0, U1 (`C/c2_t3d.py:77-83`) gives exactly B0 = sup|H_U0 - c0'| and B1 = sup|H_U1 - c0 + (5km + q/4)'|. I checked c0' = -(5/2) m p + k^2 q and (5km + q/4)' = 5(m^2 + kp/4) + (1 + kp)/4 by hand.
- E1c is the same expression as A1, so |E1c| <= A1 holds with the certified A1.
- a5 enters only through d000, and C2 passes 6 jets. ok.
- T12: E_A', E_B', dQ', dC' contain no u'' (asserted at `c2_t3d.py:89`). The spec's remark that EA1 "contains u''" is wrong but harmless.

### 3.3 beta state-form strip and Neumann feedback: correct [read + cert]
- `linear_alpha_S0` does multiply |alpha| S0. `M/ibp_codegen.py:41` defines LamA = d_alpha Lam; `M/ibp_run.py:40` evaluates it on the ball |alpha| <= alpha1; lines 44-49 feed the same kap into m_j (non-strip, "coefficient of |alpha| A") and into L1k (strip). Lam vanishes at alpha = 0, so the mean value form is exact. The old STRIP_BETA did count it with alpha = 1.
- The feedback: `M/ibp_final.py:64-67` gives x_i <= s_i + gamma_i max_k s_k/d_k for the forcing source (cPhi_total); lines 76-81 store the strip sources without it. `C/u3_strip_form.py:47-48` adds it coefficient-wise, which is an upper bound because all coefficients are nonnegative. Correct, and my recomputation agrees to the last digit.
- Linearisation of the quadratic and cubic strip forms (`u3_strip_form.py:37-38`): uses S0, S1, S2 <= 0.02 on the second factor only. The key set of `quad_gamma_L1` is exactly {00, 01, 11, 02, 12, 22} (asserted in my script), so no S3 product exists. S0, S1, S2 <= 0.02 are asserted in `closure.box_assert`.
- In E_I and E0 the S3 coefficient multiplies supv2, a state sup. No hypothesis is needed there. Sound.
- The E-matrix entries n3, n4 use the strip box on |w|, |w'| only [read `ibp_final.py:39`], both asserted.

### 3.4 Radial-degree and valuation conditions: hold [hand + cert]
- By hand: CF0 D4 DV, CBAR V4 DV, d010i D4 DV, d110i S1z D4 DV have degree exactly 2; 3 d300i S2r^2 W3 DV and d110^b S2r^2 W3 DV exactly 2; d110i S2z W3 DV, C1 V3 DV, (CA + CB S2r) k3_1 W3 DV, d110 sG 3V3 DV are 1.75. Nothing exceeds 2.
- `tail_cert.py` runs the verbatim `rhs_pn` on degree objects. `P.__mul__`/`__add__` with a constant on the left dispatch to the `Deg` reflected methods (Deg subclasses P and overrides them), and `gm` falls through to `_dsqrt` [read `V/ival.py:101-107`, `V/tail_cert.py:16-27`]. So the check is mechanical and covers the new code.

### 3.5 Hardy constant for c0' w: correct [hand]
- |w| <= S_all34 on [10.25, 34.25]; on [34.25, 57], |w(x)| <= |w(34.25)| + (Phi(x) - Phi(34.25))^{1/2} V0 because g = v there. So ||zeta c0' w||_W <= C0P S_all34 + C0Hb V0, with C0P the norm over the whole support.
- `C/c4_hardy.py:53-61`: the difference of two upper table sums is itself an upper sum of the integral over [34.25, x], because 34.25 is a grid node (16.25 + 1152/64). ok. C0Hb = 1.9e-42.
- The same Hardy form for A0, B0, EA1 (`c4_hardy.py:31-41`): W at the left end, Phi at the right end, both the correct sides.

### 3.6 v-equation (T1..T12): correct [hand]
I re-derived the modulation part from M = -beta(psi3 + 2 alpha h2) - a psi1 - b psi2 + a A1(u) + b B1(u): the splits A1(u) = A1(u0) - w + g Q_u + E_A w and zeta (g Q_u)' = (Q_u v)' - zeta' Q_u g give T8, T9, T10, T12 as stated. T8 level 0 gives psi' q/2 + q'/2 - 1; level 2 gives psi' q/2 + 5q'/2 - 1 plus the Q'' terms of I6, which I checked by one integration by parts.

## 4. Claims in V4_REPORT

| Claim | My check | Result |
|---|---|---|
| Regression to fin12f4 with NEW=0 | `compare_runs.py`: `runs/reg_fin12f4` against `step3nl_v33/verifier`, JSON-equal apart from secs, runtime_s, args, box_hypotheses; both alpha points, final_alpha_points, tail_cert log | **holds** [cert] |
| `*_pre_v4.py` are byte copies of v33; bb, ival, tail_cert unchanged | `cmp` | holds |
| Byte-identical double run | P against P_A with my comparator | **holds** apart from timing |
| Reproducibility | **my own full rerun from a copy** (`rerun/run_copy.sh` -> `rerun/P_R`, 11 min, one process, nice 19): both rate tables, final_alpha_points, tail_cert log and both verify logs equal to `runs/P` apart from timing and args | **holds** [cert] |
| tail_cert PASS on the production path | in my rerun: radial H (0,2), L (0,1); s H (2,5); joint, floors off, H (2,5), L (2,5) | holds |
| 400-point float check | rerun in my copy against the current constants file (590102ee…): PASS | holds, with two caveats (R-7) |
| Box assertions | `V/closure.py:72-85, 129, 148`; worst case over all shapes with Lyap <= Ly_hi of each cell; 20 checks (1.5e-4), 18 (1.25e-4) | holds; complete except S3 and B-e (inputs) |
| Worst box ratios | 0.789, 0.773, 0.765, 0.747, 0.32 | confirmed from the rate table |
| step3nl_v33 unchanged | `sha256sum -c audit/v33_after.sha`: 555 OK | holds |
| A_interface v35, U_cert_v35, zone3, manifest_v33, closed manifest unchanged | recomputed the five hashes | holds |
| SHA256SUMS | `sha256sum -c`: all OK (134 lines, the report says 130) | holds |
| Old/new table | rho 2.684366e-3 (x0.99678), peak 5.307698e-8, J ratios 0.99887 / 0.99928 / 1.00095 / 1.00095, delta_t 0.2 | confirmed |

**Are the J constants now true bounds?** Of the coded inequality, yes [hand]:
- `lk_v4` majorises every L term (§1.4), so Phi >= sup(L - delta N).
- The march sum in `V/closure.py:199` omits F0 * tau_t (the docstring formula at line 187 has F0, the sum does not). This is 7.19e-21, i.e. 1.9e-6 of X_const [cert arithmetic]. It is covered because the code also adds L*/2 (line 213) without subtracting Ly(tau_t)/2 at the end of the march; that slack is 1.40e-15, larger by x1.95e5. So the totals are valid; the write-up should say so.
- They remain conditional on E0 (n_A = 1/32 - PI3 in `V/unc.py:15`), on (I3) and on R-2.

## 5. Findings

| # | Severity | Finding |
|---|---|---|
| R-1 | open, medium | **Outer ramp.** FLUX0 = 1.64e-20, FLUX2 = 8.0e-20 enter as inputs (`V/rhs_iv.py:191-192`). They were derived for the old decomposition; EI_SPEC §4.10 changes K^OR (commutator form, three added items). No script bounds the new functionals, and the g''' parts need the time-integrated lemma (G3), which the tail lemma does not yet consume. v4 is honest about this |
| R-2 | open, low | **Strip hypothesis sup|w'''| <= 0.02.** Used only in the constants STRIP_BETA, STRIP_AB (`V/closure.py:22, 157`, `V/mod_sups.py:24-25`). Share of the S3 term: 4.92% of STRIP_BETA [cert], 3.96e-12 against sup|beta| 2.357e-9 (0.17%). Not implied by Lambda. Cleanest fix [hand]: integrate the 4 chi' w''' piece of the commutator pairing by parts once, so that only S0..S2 appear (needs the kappa coefficients of `lem_integrals` recomputed with one more derivative of the test function); no verifier rerun beyond `mod_sups` and the alpha-equation |
| R-3 | formal, low | **Tail alpha-condition** (`V/closure.py:156-158`). The barrier needs (-LAMT - rho/2) kappa >= ... + STRIP_BETA/sqrt(l) for every l in [L_F, L*]. The code evaluates at l = L* only. At l = L_F: 8.04e-11 / 9.888e-10 = 0.0813, against the logged margin 177 [hand]. Inherited from fin12f4 (where it was 1.07 against 177). Fix: evaluate the strip term at sqrt(L_F), or use the state form there too |
| R-4 | dependency | rho, peak, L_F and J all move if PI3 and 5.72 change (V4_REPORT item 1: tail rate -2.85e-4 with the other constants [probe]). Nothing downstream should be pinned on fin13v4 before the E0 audit decides |
| R-5 | info | The E2 modulation term (43.33 A^2, `V/rhs_iv.py:249-250`) has no alpha-linear, quadratic, forcing or strip part, while E0 and E_I now have all four. I did not audit E2; it should be on the E0 auditor's or a follow-up list |
| R-6 | records | V4_REPORT/V4_LOG say the 7 non-identical C2 regression entries are "lower because v4 evaluates at 200 bits". All 7 are in [40, 60], where `C/c2_t3d.py:149` uses box width 0.025 and `agent_int/planB/r60_t3.py:11` used 0.05. The finer cover is the explanation [read]. The values are valid either way |
| R-7 | records | `V/t_float_check.json` pins constants_v4.json 1e0fb711…, not the production file 590102ee… (it was run before pass B). My rerun against the production file passes. The test mirrors the code's structure and only level 2.5e-4, so it checks transcription and arithmetic, not the derivation |
| R-8 | records | Timing fields (runtime_s) sit inside the pinned JSONs, so the sha256 pins differ between passes and `t_double_run.py:14` has to mask them. Moving timing to a side file would make the outputs byte-identical, as the conventions ask. `V/closure.py:282` still prints "box certified: 1e-2" |
| R-9 | nit | `V/constants.py:104` forms BETA0 + dlead in float, and `V/closure.py:18`, `V/mod_sups.py:24` convert doubles through `str` to arb. In the worst case this loses the one-ulp outward margin (1e-16 relative). Covered in the march by the (1 + 1e-9) factor at `closure.py:105`; not provably covered in `mod_sups` |
| R-10 | nit | EI_SPEC §4.9 says EA1 contains u''; it does not (R3.2). `c1_cut_sups.py:119` writes qmax = 60 - TIP.lower (it should be 60 + TIP.upper); the slip is on the safe side |

I found no item of the kind "term missing", "wrong coefficient", "wrong norm", "wrong sign" or "constant below its certified value" in the E_I block.

## 6. What is still unjustified

1. **Hypothesis (I3), outer ramp** (R-1): flux lemma v2 (G2) for the new K^OR, and the time-integrated consumer (G3).
2. **E0 block** (other auditor): PI3, 5.72, MODA = 25; and whether n3 = 7.52 in the new strip term `V/rhs_iv.py:184` is the right coupling constant.
3. **E2 block and other legacy typed constants** (`constants_v4.json` `legacy_typed`): LAMT, MOD_E2, ALIN, nd1..nd8, nf1..nf7, 2.2, 1.2, r143, r163, r107, ADD_STRIP_CUBIC, and `cubic_E0fine` (box 1e-2). Not re-derived by v4, not checked by me. In particular nobody has checked that LAMT (`V/constants.py:14`) is the same number as the lamt inside the residual norms FN, FNW, FNWd (the R3 finding was a mismatch of exactly this kind in (A)).
4. **Strip hypothesis on w'''** (R-2).
5. **U-1, regularity lemma** for the two energy identities: not written. The ingredients are there (v in C^3, v'''' bounded with jumps at 10.25, 16.25, 57, 60; u0 in C^4 with bounded piecewise u0^(5)); no step in §3 needed more.
6. **B-e**: S0V = 0.0153, S1V = 0.6 are inputs from (A); the (A) chain must be rerun on fin13v4 (v36, with the R3 patches) and must return them.
7. **Modulation lemma**: I checked the strip part and the feedback (§3.3). I did not re-derive the core of the IBP lemma (lead, m, q, cPhi) or the commutator coefficients kappa in `lem_integrals.json`.
8. **Upstream certificates taken as given**: coupling (AF, BF, CF), Garding (MU1, MU2, E0, E1), the Sobolev constants c1, C1W, PHX, the forcing norms, c4fun jets, hstail profile bounds.
9. **Write-up**: `ei_lemma_v4.tex` does not exist; `sanity.py`, `sanity_unc.py` were not run.

## 7. Suggested order of work (none started; each needs Micah's go)
1. E0 audit verdict (running elsewhere). It decides whether fin13v4's numbers stand.
2. R-3 (one line in `closure.tail`) and R-2 (integration by parts in the strip commutator), together with any E0 change, in one new tag.
3. The (A) chain v36 with the R3 patches, then G2/G3 on the new K^OR.
4. `ei_lemma_v4.tex` with U-1, citing §3 here for the hand checks.

## Files in this directory
- `V4_REVIEW.md`: this file.
- `indep_consts.py`, `indep_consts.json`: my independent arb recomputation (§2.2).
- `compare_runs.py`, `compare_runs.json`: regression, double run, rerun comparison, box ratios, J numbers (§4).
- `copy/step3nl_v4/`: the copy of the verifier I ran (hashes equal to the tree, recorded in `compare_runs.json`).
- `rerun/run_copy.sh`, `rerun/P_R/`: my production rerun; `rerun/tests/`: the float check rerun.
