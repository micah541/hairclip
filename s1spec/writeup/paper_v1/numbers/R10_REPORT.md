# R10 report (unit U0): number disagreements and unsourced numbers

Ledger: 352 entries (285 cert, 2 est, 46 hand, 19 doc); `check_numbers.py` PASS (0 FAIL, 0 WARN).
Rule used: the JSON wins; upper bounds round up, lower bounds and margins round down at the printed digit (exact rational arithmetic). Labels are the new ones of OUTLINE 5.3: new (H3) = archive H2f (far datum), new (H4) = archive H3 (windows, K4). R10's "(H3) margins" are therefore the windows, ledger prefix `NHFour`.

## 1. The items listed in R10

| item | what the documents say | JSON value that settles it (source: key) | printed (macro) |
|---|---|---|---|
| windows margins (archive H3, new H4) | main_v3 Table 3: only "x15.7 in Ybar"; V33_INDEX/addendum: x5.700, x79.26, x26.50, x7.799; lemmaU_v34 Sec. 10 and `zone3_interface.json` handoff_thresholds: x5.7, x78.7, x26.6 | `FINAL_manifest_v33.json` results.M9_H3_H2v32.rows.win[j].margin_arb = 5.70027, 79.2696, 26.5089, 7.79962 (float rows 5.70029, 79.2698, 26.5089, 7.79963) | x5.700, x79.26, x26.50, x7.799 (`NHFourMargin{Zero..Three}`). The x78.7 and x26.6 are stale typed "step4_margin_if_known" values (pre-M9); x78.7 is wrong in the 3rd digit (arb 79.27). |
| windows Ybar limit | main_v3 and STEP4_RESULT: "x15.7"; addendum: x15.75 | same file: binding_Ybar_limit = 1.5756791e-11, so limit/1e-12 = 15.7568 | x15.75 (`NHFourYbarMargin`); limit 1.575e-11 (`NHFourYbarLimit`). "x15.7" is a truncation, not the rule's rounding. |
| far-family t* margin | outline x15.5; STEP4_RESULT Sec. 3 still shows the v3.1 pair 6.13e-9 vs 1.23e-8 (M = 1) | `lemmaO_v32_consts.json` far_v32.case.margin_factor_lower = 15.50518 (t* >= 9.51088e-8 against needed 6.134e-9) | x15.5 (`NFarTstarMargin`). The 1.23e-8 in STEP4_RESULT is the v3.1 (M=1) value, superseded. |
| cone fixed-L t* margin | x2.30 (main_v3, addendum) | `zone3_interface.json` acceptance.tstar_margin = 2.30194; step-4 `results.H7_fixedL.margin_tstar_over_ratio` = 2.30209 | x2.30 from the zone-III file (`NConeTstarMargin`); x2.302 from step 4 (`NHSixConeTstarMargin`). Same to 3 digits. |
| far-family Lemma O H4 | main_v3: x15.58; M9/T6 reports print x15.59 | `results.M9_H4_far.min_margin` = 15.587791 | x15.58 (`NHSixFarHFiveMargin`); x15.59 is rounded up, forbidden by the rule. |

## 2. Further disagreements found (JSON value settles each)

| item | quoted | JSON value | printed |
|---|---|---|---|
| half-width of the b_2 interval (upper bound) | 8.137e-24 (main_v3, step4_proof, outline) | `gamma0_final.json` b_j_halfwidth = 8.137112909545183e-24 | 8.138e-24 (`NGammaHalfwidth`); quoted value rounds down |
| min c_I over the 26 certificates (lower bound) | 0.05284 (outline) | `t7/constants.json` meta.summary.min_c_I = 0.0528399023 | 0.05283 (`NMinCI`) |
| smallest relative margin of the S1 containment | 1.3e-9 (main_v3) | `check_outA5.json`: min over all three comparisons = 1.28793e-9 (prof_tight, St2) | 1.2e-9 (`NSOneMinMargin`). Source file is the R7 exception file (hash differs from the closed manifest pin; ledger records `pin_status` = PIN-MISMATCH, the checker prints a NOTE). |
| cone H7 needed ratio (upper bound) | 1.6000014e-7 (main_v3, addendum) | `results.H7_fixedL.ratio_upper` = 1.6000014022604e-7 | 1.6000015e-7 (`NHSixConeTstarNeeded`); the thresholds `needed` 1.6001e-7 of the zone-III file still pass |
| C_alpha * Ybar, the error term of |alpha| | 7.3e-13 (main_v3 Table 3, STEP4_RESULT) | `t9_h1_clam_s180.json` C_alpha = 0.72061, times Ybar = 1e-12 | 7.21e-13 (`NHOneAlphaErr`); 7.3e-13 is valid but looser. Margin 1.5e-4/(|alpha|+err) = 1.2778 -> x1.27 unchanged. |
| (H3) (archive H2f) level-0 threshold | 5e-7 in `zone3_interface.json` (probe label), r2 statements, R2_DECISIONS | canonical 1e-6 (beta_A = 2, addendum Sec. 5, `NHThreeThrZero`); step-4 margin x3.1885 | the certified value 3.1363e-7 passes both |
| cone/far H4 margins x4.18 / x5.68 | CLAUDE.md "Key numbers", STEP4_RESULT | these are the v3.2 unit-L cone and far (centre 376) values: ratio 0.23923 -> x4.18005, 0.17604 -> x5.68044 (`NHSixConeHFiveMarginOld`, `NHSixFarHFiveMarginOld`). STEP4_RESULT prints the ratios as 0.2393 and 0.1761 (rounded up); 1/0.2393 would be x4.17 | In v3.3 the final hypotheses carry x1119 (fixed-L cone, `NHSixConeHFiveMargin`) and x15.58 (far, centres >= 490). Per D15 the old values go to Part K only. |
| far-family data margins | STEP4_RESULT Sec. 2: x48, x2698, x89, x5.0, x1509 (against the v3.1 thresholds) | `t9_far_v32.json` margins (v3.2 thresholds 3.53e-4, 1.17e-5, 1.03e-7, 0.08, 3.32e-6): 1.964, 450.98, 14.79, 1.926, 6.019 | use the JSON (`NHSixFar*Margin`); step4_proof prints x1.96, x450.9, x14.79 for the same rows |
| lambda_1 interval, T - t1 | outline: [6.49e-5, 6.75e-5], 2.1e-17 | `t1clock.json` lnlam, exp evaluated to 50 digits: [6.4920e-5, 6.7401e-5]; lambda_1^4/(1-4 abar) = 2.0638e-17 | agree: 6.49e-5, 6.75e-5, 2.1e-17 (derived entries) |
| (H3) K4 y in [50,350] | 5.60e-3, x10.7 | only the base value 3.1633e-4 (`amp_A1_w180.json` inputs.BASE_T1.K4) is in a result file; K4 "with error" matches `t9_cone_kch.json` cone.0.K4 = 0.0055898 by value only | kept as doc entry (flagged), see Sec. 3 |
| far boxes x6.42 / x4.69 / x1.99 / x1.69 | STEP4_RESULT: 6.4287, 4.6917, 1.9940, 1.69999995 (unrounded) | no JSON; `n3boxes_final_B.log` prints x6.43, x4.69, x1.99, x1.70 (nearest) | floors 6.42, 4.69, 1.99, 1.69 (doc entries, flagged) |

No disagreement was found in the 184 quoted numbers of `constants_needed_v34.json` (Lemma S, Lemma U, zone III, step-4 M9 rows): each equals the minimal outward rounding of its JSON value at its digits (where a number is quoted at two precisions, the ledger uses the longer).

## 3. Numbers that could not be sourced from a JSON (kind `doc`, flagged; or absent from the ledger)

In the ledger as `doc` (the checker verifies the snippet is still in the document, and the rounding of the raw document text):
- S1 enclosures of D, k(0), tau_tip (tight and conditioned) and the sigma = 8 state uncertainty: READMEs of `s1_tight` and `ext_tip` only. `NSOneStateUncert` is 4.4e-13; the outline's "|d| <= 4.37e-13" has no source.
- far boxes x6.42/4.69/1.99/1.69 (log and STEP4_RESULT only); (H4) K4 with error 5.60e-3 and x10.7.
- semigroup bound 14.8 e^{-0.0213 tau} (step3_proof_v32.tex; no certificate JSON found among the step3nl JSONs).

Not in the ledger (no machine-readable source found, or out of the numbers needed):
- M(1.45e-4), M(1.62e-4) enclosures (README text); supplement-level D, k(0) enclosures; contour vertices of step 2; SDP-side numbers.
- Theta = 1.8617364190148 +- 5.2e-11 (STEP4_RESULT; not located in a JSON); "181 base slabs"; k running from 26 to 2.4e4 (F7); wall times and sizes (T5).
- Clip parameters (arm 24, chamber 0.4, slot 0.05, W 0.85), S1 numerical spectrum (1, 1/4, -0.2562, -0.6999) (CLAUDE.md; [est]).
- Lambda_0 = 1.790819 and Lambda_1 = 1.768282 are in the ledger as `est`: they are the float SDP values in `apr/cert_b0.25_th*_tau0.03.json` (status optimal_inaccurate for theta = 0.48); the Arb verification is in `apr/ver_t03_*.log`, which is a log.
- eta <= 0.1163 and contraction <= 0.272 (Step 2) come from float values in `cap/eta_cache.json` (Arb provenance is in the script, not in the file); ledger kind `cert` with that caveat in the note.
- No `lit` entries: no literature statement used in the paper is a number.

## 4. Notes for the coordinator
- Ledger sources not pinned by either manifest: step-1 and step-2 files, `gamma0_final.json`, `P1S` (status.json), the v33 manifest itself (`pin_status` = unpinned). One source is pinned with a mismatch: `check_outA5.json` (R7/D17).
- Auto-generated macro names for the Lemma S/U/zone-III chain follow `N` + source tag + key words with digits spelled out (e.g. `\NUSupersolutionN`...); the headline ones have curated names (`NUEhatZero`, `NUBhatTwo`, `NUMarginLateMOne`, `NSAMarginCoreMOne`, `NZoneIIIBoxOne`, `NRecoveryMinMargin`, ...).
- Margins are macros of the bare number: write `$\times\NHFourMarginOne$`.
- Run: `nice -n 19 python3 -I build_ledger.py` (regenerates ledger.json and numbers_paper.tex), then `nice -n 19 python3 -I check_numbers.py [paper .tex files]`.
