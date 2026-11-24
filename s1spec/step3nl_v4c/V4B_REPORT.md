# V4B_REPORT: verifier fin13v4b, E0 fix F-B + E2 fix + V4_REVIEW R-2/R-3

Nov 20, 2026. Base: `step3nl_v4` (fin13v4), never edited. Everything new is in `/home/claude/s1spec/step3nl_v4b/`.
Work orders:
- `fix5/e0_audit/E0_AUDIT.md` §6 F-B;
- V4_REVIEW R-3, and R-2 if cheap;
- added mid-task by the coordinator: `fix5/e2_audit/E2_AUDIT.md` §5 (E2 fix).

Log: `V4B_LOG.md`. Labels: [cert] arb on certified inputs, [read] certified input taken as given, [sym] exact sympy, [hand] hand step, [num] float check, [est] estimate.

## Verdict
**fin13v4b closes at both points, alpha1 = 1.5e-4 and alpha1 = 1.25e-4.** Every certified item came in within its estimate. No margin fails, and no stop rule fired:
- J moves by at most +0.085%;
- rho moves by ×0.9975;
- the A² absorption holds with margin ×1.68.

Step 3 is **not** closed by this. See "Open" below.

Checks:
| check | result |
|---|---|
| Regression: V4B=0, final code, TAG fin13v4, against `step3nl_v4/runs/P` (`runs/reg_fin13v4`) | **PASS**. Both rate tables, final_alpha_points, tail_cert log and both verify logs are equal apart from timing and args |
| Double run: cert chain + production, twice (`v4bcert/out_A` + `runs/P_A`, `v4bcert/out` + `runs/P`) | **PASS**. Production JSON is equal apart from timing. All cert JSON is equal; it is byte-identical apart from `runtime_s` in tc_weighted.json and the sha pin of that file inside e0_identities.json. The pickle is byte-identical (`runs/t_compare_v4b.json`) |
| tail_cert | **PASS** (V4B=0 and V4B=1): radial H (0,2), L (0,1); s H (2,5); joint without floors H (2,5), L (2,5) |
| Box hypotheses (asserted at every step and at the anchor) | All hold. Worst ratios: E0/E2 alpha 0.64; E2 tip jet 0.324; E0 \|w\| on [0,34.25] 0.030; the rest unchanged (0.789, 0.773, 0.765, 0.747) |
| step3nl_v4 (589 files) and step3nl_v33 (560 files) | **unchanged**, every file by sha256 (`audit_v4b/*_before.sha` = `*_after.sha`) |
| A_interface v35, U_cert_v35, zone3, manifest_v33, closed FINAL_manifest | unchanged (sha256) |
| Proof-tree writes and caches | No file under step3nl, step3, s1_tight, step3nl_v4 or step3nl_v33 is newer than the start. No `__pycache__` was written anywhere. (step3nl_v33/verifier has 5 pre-existing .pyc files from Nov 13; not touched) |
| Resources | One process at a time, nice 19, PYTHONHASHSEED=0, `python3 -B`. Peak RSS 271 MB (tc_weighted), 212 MB (e2_nlin_regions), 79 MB (production). Production about 6 min per alpha point |

## Old vs new [cert; ratio = fin13v4b / fin13v4]
Columns:
- "E0 only" is a single complete pass with the E0 fix, R-2 and R-3 but without the E2 fix (`runs/E0only_P_A`, `v4bcert/out_E0only_A`).
- "fin13v4b" is the final run (E0 + E2 + R-2 + R-3), double-run checked.

alpha1 = 1.5e-4 (sqrt(Lambda_1) = 5e-8):
| quantity | fin13v4 | E0 only | fin13v4b | ratio |
|---|---|---|---|---|
| steps to the tail | 10 | 10 | 10 | |
| peak sqrt(Lambda) | 5.307698e-8 | 5.308856e-8 | 5.308562e-8 | 1.000163 |
| first cell rate U | +0.0303565 | +0.0304384 | +0.0304299 | 1.00242 |
| tail rho | 2.684366e-3 | 2.667904e-3 | 2.677666e-3 | 0.99750 |
| a* | 3.684593e-5 | same | same | 1 |
| sqrt(L_F) | 9.887644e-10 | 9.918102e-10 | 9.900007e-10 | 1.00125 |
| alpha_F | 6.876746e-7 | | 6.884239e-7 | 1.00109 |
| J0 const / slope | 6.020438e-12 / 3.452361e-18 | 6.031227e-12 | 6.025378e-12 / 3.455292e-18 | 1.000821 / 1.000849 |
| J1 const / slope | 8.461143e-12 / 4.85196e-18 | | 8.468087e-12 / 4.856079e-18 | 1.000821 / 1.000849 |
| J2 const / slope | 1.872555e-12 / 1.073798e-18 | | 1.874092e-12 / 1.07471e-18 | 1.000821 / 1.000849 |
| J4 const / slope | 6.242096e-8 / 3.579469e-14 | 6.253282e-8 | 6.247218e-8 / 3.582508e-14 | 1.000821 / 1.000849 |
| sup\|(a,b)\| | 2.14024e-7 | 2.14047e-7 | 2.14041e-7 | 1.00008 |
| sup\|beta\| | 2.35741e-9 | 2.35488e-9 | 2.35488e-9 | 0.99893 |
| V1 delta_t | 0.2 | 0.2 | 0.2 | |
| tail alpha-condition, whole tail | (anchor only: 177) | | margin 177.08, strip term at L_F 0.0787 | |

alpha1 = 1.25e-4:
- 9 steps;
- peak 5.213696e-8 → 5.214276e-8 (×1.000111);
- rho 2.855715e-3 → 2.849387e-3 (×0.99778; E0 only 2.838202e-3);
- J0..J4 const ×1.000657, slope ×1.000734;
- sup|(a,b)| 2.10289e-7 → 2.10295e-7;
- sup|beta| 2.22438e-9 → 2.22185e-9;
- whole-tail alpha margin 175.10.

The full tables are in `runs/v4b_compare.md` and `runs/v4b_compare.json`.

Reading the numbers:
- The E0 part alone costs 0.6% of rho and +0.18% of J0. Its largest new term is beta beyond its lead in the h2 and E33 terms (H6, ≈ alpha²·33·66.5·A²).
- The E2 part gives back +0.37% of rho and −0.10% of J. The reason is that the recertified nonlinear constants at the asserted boxes are much smaller than the typed ones (e.g. nd4 0.54 → 0.205); this outweighs the new RM, NS and FE2 terms.
- Net: rho ×0.9975, J ×1.0008.

## E0 part (fix F-B of E0_AUDIT) — status per item
| item | status |
|---|---|
| Keep PI3 = 0.01607, MODA = 25, f572 = 5.72 | **kept**. Justified by the lemma (next two rows) |
| A² budget (absorption by the PI3 slack) | **[cert given reads]**. Pi3_true (S1 + S2) = 0.0103854, slack 5.685e-3. Alpha-form need at alpha ≤ 2.5e-4 is 38.541: Qfull 6.194, n3_123·m_new 29.527, E0-4 pieces 2.27, **drift 0.549308 [cert]** (est was 0.549 with allowance 1.1). Excess 3.385e-3, **margin ×1.679** (×2.786 at alpha ≤ 1.6e-4). The lemma is stated for alpha ≤ AMAX_E0 = 2.5e-4, which is asserted |
| Forcing budget | **[cert]** 5.3439 ≤ 5.72 (×1.070) at alpha ≤ 2.5e-4. This includes every Fn piece: T2 5.3253, (a,b) in Q1 and in the loop, beta_F in the h2/E33 terms, H31/H32, alpha² remainders |
| Drift alpha-form (E0-5, was [est]) | **[cert]** C_DR_CORE = 0.549308 (`e0_strip_drift.py`). The integration-by-parts identity is [sym] exact (`tc_weighted.py`) |
| Strip alpha-linear form (E0-8, was [est] W1B ≈ 8.7) | **[cert], coded in state form instead of W1B.** The weighted divergence form T_c^{rho chi²} [sym]: X=1 reduces exactly to the closed nlin_E0.pkl, and the Euler residual is exactly 0. Strip coefficients per unit alpha: m00 4.2e-15, m01 1.60, m02 23.45, m11 1.15, m12 16.94, m22 3.21. On the strip w' = v, so the term is alpha·[m00·IR·S0s² + (m01·V0 + m02·V1)·IRW·S0s + RW·(m11 V0² + m12 V0V1 + m22 V1²)] with IR 1.84e-11, IRW 2.13e-8, RW 4.70e-5. W1B is now the raw core maximum 7.8967 |
| ALIN rounding (E0-7) | **done**: read from raw `alin_E0_spec_0.001.json` and rounded up: sq 59.327153, AB 35.868736, W1B 7.896699, BB 6.219904 |
| Phi1 in the strip commutator term (E0-6) | **done**: 2·IA·K_rho(16.25)·(PHX·A² + 2·A·W1) |
| n3 at rhs_iv.py:184 (beta_S coupling) | **n3_123 = 6.588186 [cert]**. -beta_S<g0, w~> = -beta_S<Pi g0, w~> exactly, because w~ ⊥ f1, f2, f3. So the projected norm applies; the full norm 7.6970 is not needed. The other typed 7.52 (`qb752`, beta_Q) is valid because 7.52 ≥ n3_123; it is left as typed |
| E0 modulation lemma | `e0_lemma/E0_LEMMA_v4b.tex`, about 2 pages, **not compiled** (no LaTeX here). Covers the beta decomposition, S1 projection, S2 rank-one symmetrisation, S3 oblique projection, the budgets, the coded terms and a constants table |
| Machine checks | `e0_identities.py` I1–I8 [sym], all exact: rank-one formula, oblique projection, modulation algebra, (a,b) w-parts, cubic drift, alpha² remainders, chi_32; plus random-float I1/I2. `e0_lemma_consts.py` is the promoted e0_cert.py; its shared sections equal the audit's JSON |
| **New finding, coded** | The audit's E0-4 list misses several higher-order modulation pieces. All tiny, all now **coded as explicit H terms** with certified constants (none absorbed): <br>• H3: (a,b) w-parts of A1(u), B1(u) (CA1 32.0, CB1 1.032, ...); <br>• H2: cubic drift (CD1 0.607); <br>• H4: QAB·Q and strip parts of (a,b) in Q1 and in the loop (GABN 8.03, LOOPC 5.65); <br>• H5: alpha² remainders of A1(u0), B1(u0) (RA2 26.9, RB2 5.50); <br>• H6: beta beyond its lead in the h2/E33 terms (HOB 33.37); <br>• H7: H31/H32 part of beta_E. <br>New hypothesis: \|w\| ≤ 0.02 on [0,34.25], asserted (ratio 0.030) |

## E2 part (E2_AUDIT §5, added at the coordinator's request) — status per item
| item | status |
|---|---|
| (1) RM + NS terms × Dq | **done**, exactly as in §5 F2. modE2 stays 43.33 (≥ 43.3271) |
| (2) r143 | **R143 = 14.343002687670673**, read from JSON and rounded up |
| (3) E2 forcing factor | **FE2 = 1 + BETAF·(n3 + 2α‖h2‖) + PHIAB·(NF12 + u0pr)** multiplies the three E2 entries of L_list. It is ≈ 13.98 at α = 2.5e-4. It is computed in P-arithmetic in alpha in rhs and in arb in `unc.lk_v4` (same shared list) |
| (4) Assertions | E2_alpha ≤ 2.5e-4 (ratio 0.64); E2_tip_jet sup(\|w\|,\|w'\|) ≤ 1e-4 on [0,8] with K = KREG[4,8] (ratio 0.324) |
| nd1..nd8, nf1..nf7 | **from the recertified JSON at the asserted boxes, rounded up**: nd = max over [0,8] (box 1e-4) and [8,12]; nf = max over [12,16.25] and [16.25,34.25]. Every typed value covers its new value |
| Constants | `e2_nlin_regions.py` and `e2_consts.py` were promoted into `v4bcert` and run in the chain; their outputs equal the audit's JSON in every section. `e2_v4b_constants.py` writes `e2v4b_constants.json` (raw → one ulp up; no typed literals except the two asserted hypotheses) |
| Float check of the E2 identity | `e2_identity_check.py`. [sym] D_rho expansion; rho f''g'' − rho f D_rho g = (B)' exactly; B(0) = 0 for even data. [num] quadrature on random even data, symmetry error 4e-13; ‖D_rho g‖ ≤ Dq-form on the test data (ratio 0.38). The Cauchy–Schwarz forms of RM themselves stay **[hand]** |
| Text | An "E2 couplings" paragraph has been added to `E0_LEMMA_v4b.tex` |
| LAMT | No change. The audit shows the typed value is on the safe side |

## V4_REVIEW items
- **R-3 (tail alpha-condition along the whole tail): done.** `closure.tail` now checks rhs_all = (BETA1, BETAF·Fn, QBL terms at L*) + STRIP_BETA/sqrt(L_F). This is rigorous by monotonicity: the first terms are nondecreasing in l and the strip term is decreasing. Margin 177.08 (1.5e-4) and 175.10 (1.25e-4). The result is stored as `tail.alpha_condition_whole_tail`. Independent cross-check `verifier/t_tail_alpha.py` (arb, endpoints + 4001 log-spaced l): min margin 177.11 / 175.13 [num-arb].
- **R-2 (strip hypothesis sup|w'''| ≤ 0.02): removed.** `v4bcert/r2_strip_ibp.py` [cert] integrates the 4χ'w''' commutator term by parts once (χ' vanishes at both strip ends). The new w'' coefficient is rho|2fχ'' − 4f'χ' + 4φ'fχ'|; f_j' comes from the outA5 strip rows. The old kappa match lem_integrals.json to < 1e-9 relative, and the mod_ibp commutator coefficients are reproduced. After the IBP the S3 coefficient is 0: STRIP_BETA 8.040e-11 → 7.787e-11 and STRIP_AB 4.368e-10 → 4.227e-10. These now need only S0, S1, S2 ≤ 0.02 (asserted). They are used in the alpha equation, the tail and mod_sups.

## Files
- Certification chain: `v4bcert/run_cert_v4b.sh`. Stages: tc_weighted, e0_strip_drift, e0_lemma_consts, e0_identities, r2_strip_ibp, e2_nlin_regions, e2_consts, e2_v4b_constants, e2_identity_check. Outputs go to `v4bcert/out/` (pass A: `out_A/`; E0-only pass: `out_E0only_A/`). Verifier inputs: `e0v4b_constants.json`, `e2v4b_constants.json`, `r2_strip_ibp.json`.
- Verifier: `verifier/{constants,rhs_iv,closure,unc,mod_sups}.py`. Changes are gated by V4B (default 1). Backups `*_pre_v4b.py` are byte copies of fin13v4. Also `run_v4b.sh`, `t_compare_v4b.py`, `t_tail_alpha.py`, `v4b_compare.py`.
- Runs: `runs/P` (production, pass B), `runs/P_A`, `runs/reg_fin13v4` (final regression), `runs/reg_fin13v4_prefinal` (regression before the R-2/E2 patches; also PASS), `runs/E0only_P_A`, `runs/t_compare_v4b.json`, `runs/v4b_compare.{md,json}`, `runs/P/t_tail_alpha_fin13v4b.json`.
- Driver: `run_prod_v4b.sh`. Hashes: `SHA256SUMS`; `audit_v4b/`.
- Copied v4 scripts that pointed at step3nl_v4 (`v4cert/v4common.py`, `run_cert_chain.sh`, `s1_identities.py`, `verifier/t_double_run.py`, `run_v4.sh`) were repointed to step3nl_v4b, with backups `*_pre_v4b.*`. This is so that nothing in v4b can write into v4. None of them was run.

## Open (not done or not justified here)
1. **[read] inputs** of both lemmas: gram7/gram8 Grams, pair_f3, bloop.json (Qfull), final_5e-4.json (m_new, BETA*, QB, CAB0, MAB, QAB, PHIAB), lyap_consts2, mod_forms Q2 (43.3271), garding_cert, coupling constants, the profile hulls.
2. **[hand] steps**:
   - the beta decomposition and its identification with coded terms (I3 checks only the algebra);
   - the Cauchy–Schwarz forms of H3–H7 and of the E2 RM/NS;
   - the audit suggested a float check of the beta decomposition on random w, as in `agent_mod/ibp/float_check.py`; this was not done.
3. **Not audited**: the theta-certificate (1/32, 1/9), `cubic_E0fine` and ADD_STRIP_CUBIC = 2.6e-3 (the cubic strip term is still the scaling hack), the Q2 Gram, the Gårding identity behind mu2/e1/e0, the level-2 commutator constants DS.
4. **Unchanged from V4_REPORT/V4_REVIEW**:
   - outer-ramp FLUX0/FLUX2 for the new K^OR (R-1, G2/G3);
   - the U-1 regularity lemma;
   - B-e (S0V, S1V are inputs from (A));
   - `ei_lemma_v4.tex`;
   - legacy typed constants outside E0/E2: MOD_E2 43.33 (valid), t22/t12 (valid), qb752 (valid), LAMT (valid by the E2 audit), nd/nf (now from JSON).
5. `E0_LEMMA_v4b.tex` has not been compiled.
6. **Downstream not touched**: the (A) chain v36 with the R3 patches, step 4, kits, manifests. Re-pinning on fin13v4b needs Micah's go. The J, peak, alpha_F, L_F, abar and betabar to use are those of `runs/P`.
