# Documentation sync, v3.3 of step 3 (Nov 14, 2026)

Only wording and numbers were changed. No mathematics, no closed file, nothing in `step4/`. Every edited file was backed up first as `*_pre_sync.<ext>`. Numbers come from the final certified sources (A_interface 34832f43..., asm_design_N2ref4 6ad175c1..., U_cert_v34 271430d7..., zone3_interface a381b1f2..., ramp_U 55317ce5..., FINAL_manifest_v33 8f6f7f95... with 19 stages, V33_VERIFIER fin12f4). Bounds rounded up, margins rounded down. All paths are under `/home/claude/s1spec/writeup/lemmaU/bulk/` unless absolute.

## Checks
| check | result |
|---|---|
| `make_constants_v34.py` (PIN = 8f6f7f95...) | wrote 104 rows, 184 numbers |
| `check_keys_v34.py` | **PASS**: 7 sources pinned (S4 = 8f6f7f95...), 132 cited keys covered, 184 quoted numbers outward, side conditions hold (K6 1.19013, B2' 2.11604, sqrt2*S1 = 0.147486) |
| balance `lemmaU_v34.tex` | PASS (36 labels, 26 refs) |
| balance `lemmaS_A_v33.tex`, `numbers_A_v33.tex`, `erratum.tex` | PASS |
| balance `main_v3.tex` (`t3_proof/check_balance.py` and `check_balance_main.py`) | PASS (no undefined refs) |
| balance `step3_v33_main.tex` | FAIL, **pre-existing and unchanged**: the same 9 undefined refs (`A:*`, `C:*`) that come from the other files (`lemmaS_A_v33.tex`, `continuation_v33.tex`) when it is checked alone; the pre_sync copy gives the identical result |

## Edit 1. `fix1/lemmaU_rev/lemmaU_v34.tex` (backup `lemmaU_v34_pre_sync.tex`; new sha256 dc178ba0...)
- l.27: S4 is now 19 stages, `8f6f7f95...`, final after M9; the M9 arb stages are listed (m9_balls/win_arb, h3_arb, h2t_v33_arb, h2u_v33_arb, h2f_camp_arb, t9_lemmaO_v33_arb, h4far_v33) and the arb rows are said to supersede the float rows.
- l.505 ((H2) item): added the M9 arb row `results.M9_H2_tightened_fU1`, margins x27.34, x562.8, x600.9.
- l.507 ((H2f) item): added the M9 arb row `results.M9_H2f` (x3.1885, x4.9834, x4.1067). The float cell half-width caveat (inherited step-4 caveat) is replaced by the statement that it is closed by `win_arb.py` (arb element length and sums; arb values within 7.5e-5 relative).
- l.510 ((H5) item): the **far-family coverage caveat is removed** and replaced by the certified `h4far_v33` result: every far centre at radius >= 490, 37 rows, minimum margin x15.58 (`results.M9_H4_far`), plus the centre-map row `results.M9_H1_centre_map`.
- l.531 (O2): B2' wording now says it is aligned in `lemmaS_A_v33.tex` and `step3_v33_main.tex`.
- l.535: "S4 in flux" item replaced by "S4 and review M9: final", citing the M9 rows and the closed caveat.
- l.536: the "far-family coverage item" is no longer a condition of the statement.
- `make_constants_v34.py` (backup `_pre_sync.py`, regenerated `constants_needed_v34.json` with backup `_pre_sync.json`): l.96 the M9 row is now final and also lists `results.M9_H1_centre_map`; l.185 PIN = 8f6f7f95...; l.186 on: new `num()` lines for the M9 margins (27.34 / 562.8 / 600.9; 3.1885 / 4.9834 / 4.1067; 15.58). No other script line was changed.

## Edit 2. `fix1/b12_text/numbers_A_v33.tex` (backup `numbers_A_v33_pre_sync.tex`)
- l.4, 7, 8: header hashes (A 34832f43..., U 271430d7..., S4 19 stages 8f6f7f95...).
- l.16-17: `\HtwoCertV` (3.658e-10, 1.777e-9, 1.665e-6) (the old 3.657e-10 and 1.664e-6 were not outward); `\HtwoMarV` x27.34, x562.8, x600.9 (old x563 and x601 were rounded up).
- l.26-27: `\EhV` (1.2021e-6, 1.0440e-4, 0.056686); `\EAzeroV` (1.1980e-6, 1.0434e-4, 0.056163).
- l.37, 39, 41: `\coupV` 1.1571e-12, `\CzeroV` (1.1345e-6, 1.6621e-4), `\NbarV` (0.071472, 0.33012, 5.0403).
- l.44-49: `\SzerovV` 7.1881e-3, `\SonevV` 0.10429, `\MarFluxZeroV` 3.363, `\MarFluxTwoV` 5.953, `\MarSzeroV` 2.116, `\MarSoneV` 5.748, `\MarMzeroV` 5.449.
- l.51, 55: `\MarRampMzeroV` 7.020; `\MarFluxRefV` 3.395, 6.010.
- l.57: `\PhiV` 3.3197e-5.
- l.63, 65-70: `\BhV`, `\UmarV` (x2.353, x4.892, x2.770, x2.225), `\BhFineV`, `\KsixV` 1.190, `\GAoneV` (6.7e-9, 6.5e-8, 3.0e-6) [est, measured in M14]; the comment says G-A1 is applied.
- Checked and left: `\MarMoneV` 4.409, `\MarGV`, `\MarRampMoneV`, `\MarRampZoneV`, the d tables, the verifier macros.

## Edit 3. B2' wording
- `fix1/b12_text/lemmaS_A_v33.tex` (backup `_pre_sync`): l.86 (Lemma z1) now reads sup|v_A| + sup_core|d| <= 0.0153, margin x\MarSzeroV (= 2.116) for the derived value, with |z0| <= sup|v_A| and sup|w'| <= 0.0153.
- `fix1/b12_text/step3_v33_main.tex` (backup `_pre_sync`): l.66 (B2') same wording.
- Matches lemmaU_v34.tex l.383.
- Also (the "script values until rerun" wording, since G-A1 is fixed): `lemmaS_A_v33.tex` l.32, 119, 280, 312 and `step3_v33_main.tex` l.12, 84, 102 now say G-A1 is applied (M14 rerun); l.312 "Fix:" became "Fixed (M14)" and no longer claims Ê^true <= Ê^scr + GA1 or a margin 3.363 -> 3.36 jump; the macro `\GAoneV` is the measured increase. `numbers_A_v33.tex` l.69 comment.

## Edit 4. `t7_assembly/main_v3.tex` (backup `main_v3_pre_sync.tex`; log appended to `main_v3_CHANGES.md`, backup `_pre_sync.md`)
- l.111 (B3): floor sqrt(Lam_F) = 9.872e-10, |alpha| <= 6.867e-7 (fin12f4: 9.871704e-10, 6.866519e-7), replacing 4.94e-10 and 3.43e-7. Rounded up: the 3-digit values 9.87e-10 and 6.87e-7 in the task text are 9.87e-10 < 9.871704e-10 (not outward), so 4 digits are used.
- l.290 (H4): [48, 700].
- l.294: (W1) 2.192e-7 (was 2.19e-7, not outward); (W2) 1.161e-3 (was 1.16e-3, not outward); (W3) and (W4) rewritten (M7) to the lemmaU_v34 and continuation_v33 deliveries (E-hat, B-hat, boxes, zone III, Lemma O bounds); the v3.2 box |f| <= 1/2 on [48, R(tau)] is stated as replaced.
- l.301-303, 305: forcing bullet (peak <= 5.308e-8, floor <= 9.872e-10); Lemma S(A) bullet with the final margins (the "phase-2a/v3.3 margins" at the old l.316 area: E2 0.056686, M1 x4.409, S0v x2.116, G22 x1.917, flux x3.363 / x5.953); Lemma U bullet (x2.353, x4.892, x2.770, x2.225, loop gain 0.03904, K6 x1.190, c_dU 1.7187e-8, zone III box (0.10032, 0.41974, 0.014291) x1.238); H6 bullet (U-fine B-hat).
- l.312: status paragraph margins x3.188, x4.983, x4.106; x124.7; x27.34, x562.8, x600.9.
- l.344: step-4 v3.3 stages: 19 stages, 281 paths, 35 scripts, 85 reused pins, sha256 8f6f7f95...
- l.354-356, 361, 363, 367: hand-off table, with non-outward roundings fixed (H2 base 8.02e-11; H2f base and with-error values; x3.03 for the cone h0'' margin, was x3.04), new f - U_1 row and far-family H4 row (x15.58).
- l.384, 411: manifest description (19 stages, 281 paths, 35 scripts, 85 pins); open item 3 (H4 on [48, 700]; far-family caveat closed).
- The boxed "unconditional subject to the pending independent review" note (l.64, l.90) is untouched.

## Edit 5. `t7_assembly/V33_INDEX.md`, `STEP4_RESULT_v33_addendum.md`, `KITS_TODO.md` (backups `*_pre_sync.md`)
- Each got a banner at the top (l.5-6) saying the old sections are the Nov 13 text, and a new final section: `V33_INDEX.md` section 4 (l.166-249), `STEP4_RESULT_v33_addendum.md` section 5 (l.51-83), `KITS_TODO.md` section D (l.63-97).
- Contents: the new components (t2_A_v34, B0 stages 10-11, M9 stages 12-18, M14 fixes, lemmaU_v34 and its checkers, fix1/b12_text); final sha256 table; the 19-stage manifest and the arb margins; the rerun order including assemble pass 2 and pass 3 (V33_INDEX 4.4); kit file lists including the 8 M9 scripts, `b0_dprof.py` and `h2u_v33.py` with byte sizes, the `b0_dprof` inputs, and the reference directory (281 paths).
- The old sections were not rewritten; they are marked superseded.

## Edit 6. `t7_assembly/erratum/ERRATUM.md` and `erratum.tex` (backups `*_pre_sync`)
- ERRATUM.md l.35-36, 179-181, 191, 194-195, 313-314, 356-366; erratum.tex l.38, 139, 147, 152-153, 211, 213, 238-246.
- x1.1245 is stated as certifying the flux line only, with the A5 ramp box and (H2f) at 5e-7, and the closed-version claim is reworded per review M8; M9 now certifies G0 = 3.1363e-7 (x1.59 against 5e-7).
- x0.7317 is labelled [est from cert components] (hand arithmetic from a probe), not [cert].
- (H2f) margins x3.188 / x4.983 / x4.106; tightened (H2) x124.7 / x2124 / x641; B0 f - U_1 row x27.34 / x562.8 / x600.9; far-family H4 minimum x15.58; new "v3.3 final" block with the t2_A_v34 values at N2_ref and at 4 x N2_ref, B2' and G-A1.

## Rounding changes against the source texts (for the PI)
- Far-family H4 minimum: M9_REPORT prints x15.59; the JSON value is 15.5878, so x15.58 is used.
- (H2) for f - U_1, level 2: the task text and M9_REPORT give x601.0; the JSON value is 600.959, so x600.9 is used. Level 1: 562.81 gives x562.8; level 0: 27.341 gives x27.34.
- (H2f): x3.1885 / x4.9834 / x4.1067; at 3 digits x3.188 / x4.983 / x4.106 (the old x3.19, x4.11 were rounded up).
- Tightened (H2) of w': 124.7895, 2124.96, 641.54 give x124.7, x2124, x641 (the old x124.8, x2125 were rounded up).
- Thinnest zone III margin 1.2389: x1.238 (the old x1.24 was rounded up); cone h0'' margin 3.0378: x3.03 (old x3.04); loop gain 0.039035: 0.03904 (old 0.0390); c_dU 1.7186043e-8: 1.7187e-8 (old 1.7186e-8).

## TODO markers and items not edited
- `[TODO: source]` x1 in `V33_INDEX.md` 4.4: the position of `ramp_U.py` in the M14 rerun order is not recorded in M14_REPORT.
- `[TODO: source]` x1 in `KITS_TODO.md` D.2: total run time of the M9 stages is not in M9_REPORT.
- `[TODO: source]` in ERRATUM.md and erratum.tex: the (O3) outer-data values for the v3.3 column on t2_A_v34 (the existing v3.3 column comes from t2_A_v33; `o3_outer_check.py` was not rerun).
- **Not edited (outside the list), stale and the owner should align them:** `fix1/b12_text/continuation_v33.tex`: its header cites U_cert 4c3f23bc... and zone3 8152c650..., (W3) at l.~112 says "subject to G-A1" and "<= GAoneV more", and its gap item 3 still says far-family H4 is certified only at centres 376 and 602 and conditions (W4) beyond 602 on it. Also `zone3_interface.json handoff_thresholds` (5e-7), `STATEMENTS_r2.tex`, `R2_DECISIONS.md` (M9 notes they are being fixed separately), and `ramp_U.json`'s pre-M4 u1 pin.
- Not rerun: `o3_outer_check.py`, `u5_enum`, any kit. Nothing was compiled.
