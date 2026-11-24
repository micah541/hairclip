# LEMMAU_MAIN_v37_NOTES: change log for lemmaU_v37.tex and main_v3_v37.tex

Nov 20, 2026. Files (all new; the old ones are unchanged):
- `fix12_texts/lemmaU_v37.tex` (from `fix1/lemmaU_rev/lemmaU_v34.tex`, copy + targeted edits)
- `fix12_texts/main_v3_v37.tex` (from `t7_assembly/main_v3.tex`, copy + targeted edits)
- this file.
Number source: `fix12_texts/numbers_A_v37.tex` / `.json` (macro value, source file, key). Other sources as listed per item. Labels: all numbers kept as [cert] exactly as before; nothing was newly labelled. Not compiled (no LaTeX). Helper scripts were in the scratchpad (not part of the proof tree). No proof-tree script was run; only JSON/MD/tex were read.
`B` = `/home/claude/s1spec/writeup/lemmaU/bulk/`.

## 0. Counts
- lemmaU_v37: 59 occurrences replaced by macro mapping (v35 macro value -> v37 macro value, 29 distinct pairs, table in section 2) plus about 55 manual occurrences (section 3), plus new text (Lemma dU, Lemma D k=2, (I3'), status box, hashes/paths).
- main_v3_v37: 14 occurrences by macro mapping plus 12 manual number edits (section 4), plus the status box and open item.

## 1. Sources and hashes cited (full sha256, computed with sha256sum)
| file | sha256 |
|---|---|
| fix6_v36/t2_A/out/A_interface.json | 799b087cf192f41df54986feba7a0ec715b407e2b53f0cc89db245d7c7415c4a |
| fix6_v36/t2_A/out/asm_design_N2ref4.json | a9267ad701088a781656d65562cdae2469c7c8fbd21127cfedac7c8a556689a1 |
| fix6_v36/t4_cert/out/U_cert_v36_pass2.json | d512596ab0868da5219ef41acd1d62e7d09266503407ea4b43b2502ed1d72bb5 |
| fix6_v36/t4_cert/out/U_cert_v36_pass3.json | 6d52881a5607b08167086d5e48886c3ef7b71e16a2c822bd7f46ef4f6b84a27c |
| fix6_v36/t4_cert/out/u34_primary.json | 0967bed548d5abf2c46d4020e8c91f98c6dc256da795bcd6cd9160c19590b01b |
| fix6_v36/t4_cert/out/u1_ghat.json | 03d367c941a07ca7cbf583db7657619415fc5e19e5187747ef78a5f9387ef424 |
| fix6_v36/t5_zone3/out/zone3_interface.json | d6813ddb17051d18905e721434d27d99e2e29c31998e991287673f844a6dd4a2 |
| fix6_v36/t4_cert/out/U_cert_for_zone3_v36.json | 1344d5bf0399f0df96cd95fcc5e7e7138e0538ab15eac373eefdb07e5968f1c1 |
| fix6_v36/t5_zone3/out/ramp_U.json (byte copy of the v35 file, equal sha to t5_zone3/out/ramp_U.json) | 55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de |
| fix8_flux2_impl/out/flux2_cert.json | 300bf7ff35d41da6cf90e86b824cc9520771aaa45e326729a8fbdd04e778fe44 |
| fix6_v36/r3/r3_forcing.json | f763f73d169aa7cbcd308347254a2f1c2696520f3572306fd5c49696f232cc5a |
| fix10_proofs/G1_export_k2.tex | 4004a774f5ef41764c9feb32327b8d42c972f70340ec4cdca9cf79b98cf45e05 |
| fix10_proofs/U1_regularity.tex | 827a5f51ed7af4d7a21c663183b6cdb06a63fffb64c420fa795d49f78ff2e9e7 |
| step3nl_v4b/runs/P/rate_table_fin13v4b_a0.00015.json | d1dcd10f639c83a5de15075dfa31702535fb93f1bd7c3fa3e6fb5af465a05007 |
| step3nl_v4b/runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json | 3c4b17228dc8aa11f603a42a808a4637432213679738eb078a35fcf720f0ce65 |
v35 hashes that were replaced: A_interface 78fe2f8e... -> 799b087c...; asm_design_N2ref4 46bae149... -> a9267ad7...; U_cert_v35 960809fd... -> d512596a... (pass 2); u1_ghat 5a06d24c... -> 03d367c9...; zone3_interface 45e8d4dc... -> d6813ddb...; U_cert_for_zone3 22c13e44... -> 1344d5bf.... (`ramp_U.json` 55317ce5... unchanged.) Not cited in the texts: A_interface of R3-on-fin12f4 (d8d41fbe...), and the n1_chart v36 file (7f64040f..., `fix6_v36/n1/out/n1_chart_fin13v4b.json`); V36_REPORT says (A) consumes it, but the texts cite the step-4 `n1_chart.json` (341ed1f4..., still in the text of lemmaU) -- see conflict C6.

## 2. lemmaU_v37 / main_v3_v37: changes by v35->v37 macro mapping (bounds up, margins down, same digits as before)
Columns: old token -> new token | occurrences (lemmaU / main) | macro (source: JSON key in numbers_A_v37.json).
| change | occurrences | macro (source: key) |
|---|---|---|
| `4.8723\cdot10^{-6}` -> `4.8722\cdot10^{-6}` | lemmaU 2 / main 1 | macro `BhFineV` (U: ufine.Bhat) |
| `2.6266\cdot10^{-4}` -> `2.6267\cdot10^{-4}` | lemmaU 2 / main 1 | macro `BhFineV` (U: ufine.Bhat) |
| `0.025842` -> `0.025850` | lemmaU 2 / main 1 | macro `BhFineV` (U: ufine.Bhat) |
| `5.5933\cdot10^{-4}` -> `5.5949\cdot10^{-4}` | lemmaU 2 / main 1 | macro `BhV` (U: supersolution.Bhat_num) |
| `0.36675` -> `0.36718` | lemmaU 2 / main 1 | macro `BhV` (U: supersolution.Bhat_num) |
| `1.5241\cdot10^{-9}` -> `1.5227\cdot10^{-9}` | lemmaU 2 / main 0 | macro `ConeV` (U: supersolution.C1) |
| `5.5255\cdot10^{-8}` -> `5.5261\cdot10^{-8}` | lemmaU 2 / main 0 | macro `ConeV` (U: supersolution.C1) |
| `7.2060\cdot10^{-10}` -> `7.1914\cdot10^{-10}` | lemmaU 2 / main 0 | macro `ConeZeroV` (A: res.tight.C1_0[0:2]) |
| `4.3143\cdot10^{-8}` -> `4.3148\cdot10^{-8}` | lemmaU 2 / main 0 | macro `ConeZeroV` (A: res.tight.C1_0[0:2]) |
| `1.1340\cdot10^{-6}` -> `1.1328\cdot10^{-6}` | lemmaU 2 / main 0 | macro `CzeroV` (U: supersolution.C0) |
| `1.6623\cdot10^{-4}` -> `1.6628\cdot10^{-4}` | lemmaU 2 / main 0 | macro `CzeroV` (U: supersolution.C0) |
| `1.8092\cdot10^{-7}` -> `1.8090\cdot10^{-7}` | lemmaU 2 / main 0 | macro `DsrcV` (A: res.tight.Dsrc[0:2]) |
| `9.8276\cdot10^{-5}` -> `9.8315\cdot10^{-5}` | lemmaU 2 / main 0 | macro `DsrcV` (A: res.tight.Dsrc[0:2]) |
| `1.2018\cdot10^{-6}` -> `1.2004\cdot10^{-6}` | lemmaU 2 / main 1 | macro `EhV` (U: supersolution.EA) |
| `1.0441\cdot10^{-4}` -> `1.0443\cdot10^{-4}` | lemmaU 2 / main 1 | macro `EhV` (U: supersolution.EA) |
| `0.056694` -> `0.056798` | lemmaU 2 / main 2 | macro `EhV` (U: supersolution.EA) |
| `1.189` -> `1.188` | lemmaU 2 / main 1 | macro `KsixV` (U: supersolution.EA[0] / 1.01e-6) |
| `1.917` -> `1.916` | lemmaU 2 / main 1 | macro `MarGV` (A: res.tight.core.margins.Ghat22) |
| `4.409` -> `4.407` | lemmaU 2 / main 1 | macro `MarMoneV` (A4: res.tight.margins.M1core) |
| `5.748` -> `5.745` | lemmaU 2 / main 0 | macro `MarSoneV` (A4: res.tight.margins.S1v_vs_0.6) |
| `0.071472` -> `0.071471` | lemmaU 2 / main 0 | macro `NbarV` (U: supersolution.N2) |
| `0.33012` -> `0.33011` | lemmaU 2 / main 0 | macro `NbarV` (U: supersolution.N2) |
| `5.0406` -> `5.0429` | lemmaU 2 / main 0 | macro `NbarV` (U: supersolution.N2) |
| `0.10430` -> `0.10435` | lemmaU 3 / main 0 | macro `SonevV` (A4: res.tight.S1) |
| `7.1878\cdot10^{-3}` -> `7.1860\cdot10^{-3}` | lemmaU 2 / main 1 | macro `SzerovV` (A4: res.tight.S0) |
| `4.892` -> `4.893` | lemmaU 2 / main 1 | macro `UmarV` (U: margins) |
| `0.11961` -> `0.11954` | lemmaU 2 / main 0 | macro `gpV` (A: res.tight.gp) |
| `0.042351` -> `0.042432` | lemmaU 2 / main 0 | macro `kappaJV` (HJW: kappa_J (certified + pad 1e-9)) |
| `5.4995\cdot10^{-7}` -> `5.5016\cdot10^{-7}` | lemmaU 2 / main 0 | macro `rhoJV` (A: res.tight.rhoJ) |

## 3. lemmaU_v37: manual changes (old -> new, source file + JSON key)
| item | old -> new | source |
|---|---|---|
| abar = bbar (sec. 1 text and table) | 2.19103e-7 -> 2.14041e-7 | `step3nl_v4b/runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json` key `1.5e-4.sup_ab` (macro abarV) |
| betabar | 3.33811e-9 -> 2.35488e-9 | same file, `1.5e-4.sup_beta` (betabarV) |
| c_dU (Lemma dU, table row) | 1.7187e-8 -> 1.46343e-8; x1.366 under the typed 2e-8 (new remark) | `fix6_v36/r3/r3_forcing.json` `c_dU.new`, `2e-8/c_dU.new`; U key `launch.c_dU` |
| D_L (new) | 3.2224577463e-6 (exact), betabar' = 2.83825e-9 | r3_forcing.json `DL`, `c_dU.betabar_prime` |
| Lemma dU hypothesis | alpha_tau = -(p/4) alpha + beta -> alpha_tau = lambda~ alpha + beta, lambda~ = -0.2561532224577463; new formula and proof identity (as specified) | R3_REPORT.md section 2 (iv) |
| kappa_J (W1 part, W3 part) | 0.0423485 -> 0.0424300; 1.93e-6 -> 1.92e-6 (pad 1e-9 noted) | A_interface `HJW.kappa_J_W1`, `HJW.kappa_J_W3` |
| (W2) tail text | fin12f4 closure statement -> fin13v4b: margins x177.07 / x175.09, alpha_F <= 6.8843e-7, rho >= 2.6776e-3, sqrt(L_F) <= 9.9001e-10, peak sqrt(Lambda) <= 5.30857e-8 | RT15 `tail.alpha_condition_whole_tail.margin_whole_tail` (RT125 for x175.09), `tail.alpha_F`, `tail.rho`, `tail.sqrt_LF`, `sqrtL_max` |
| F_k^2 (k=0,1; k=2 new) | (5.6506, 5.9320)e-13 -> (3.5932, 3.8198, 4.1573)e-13 | A `res.tight.F2` (macro FtwoV) |
| coup_k | (0, 1.1568e-12) -> (0, 1.1555e-12, 2.0144e-8) | A `res.tight.coup` (coupV) |
| D^s | (8.7399e-10, 5.2599e-8) -> (8.7390e-10, 5.2607e-8, 6.05091e-6) | A `res.tight.B0.D_window_slope`; k=2: `r2_export.window_data.D_window_slope` |
| D^src (k=2 new) | 0.0552611 | A `res.tight.r2_export.window_data.D_2` (DtwoV) |
| D^c | (1.7732e-7, 9.8058e-5) -> (1.7730e-7, 9.8097e-5, 0.0552107) | `asm_design_N2ref4.json` `res.tight.info.{0,1}.D_window_const` (arb mantissa decoded by hand: 1.772957e-7, 9.809675e-5); k=2 value 0.0552107 from V36_REPORT section 6, NOT re-read from JSON (see TODO T3) |
| C1^0 (k=2 new) | 6.21886e-3 | A `res.tight.r2_export.C1_0` (ConeTwoZeroV) |
| Lemma D, k = 0,1 -> k in {0,1,2}; r_k | r_export = (1249/1000, 1249/1000, 1/1000) | A `res.tight.r_export` |
| C0_2, C1_2 | 116.3782, 0.2606004 at 4 N2_ref; 114.783, 0.0698143 at N2_ref | A `res.tight.r2_export.C0_ref4/C1_ref4/C0_ref/C1_ref` |
| Delta^in (3 values, text + table) | (3.131,3.287,3.451)e-12 -> (2.271,2.384,2.502)e-12 (table 3.451e-12 -> 2.502e-12) | U `launch.Delta_in` (u34_primary) |
| data terms D (late centre) | (1.7330e-6, 2.5353e-4) -> (1.7313e-6, 2.5361e-4) | U `data_terms_late.D` |
| Phi_j (late), table | (6.434e-9, 1.4217e-8, 3.8404e-8) -> (6.387e-9, 1.4169e-8, 3.8357e-8) | U `launch.Phi` |
| S_1(2 Bhat) zones I/II (3 places) | 0.28585 -> 0.28595 | U `supersolution.M1z_hyp_check.zones_I_II.S1_2Bhat` |
| n_III (text, table) | (0.071445, ...) -> (0.071444, 0.32707, 3.0445) | U `inputs.n_III` (zone III n_III 0.0714437) |
| Delta^{R,b} (text, table) | (2.66e-12, 3.26e-9, 4.35e-8) -> (2.60e-12, 3.18e-9, 4.25e-8); table 4.35e-8 -> 4.25e-8 | key table v36 rows "zone III DeltaR_b k=0..2" (2.59184e-12, 3.17669e-9, 4.24675e-8), 3 digits rounded up |
| loop gain level 2 through A / total | 0.0020172 -> 0.0020183; 0.039035 -> 0.039036 (2 places) | U `loop_gain.2.through_A`, `loop_gain.2.total` |
| U-fine K2, K3 (text, table) | 0.046956 -> 0.046958; 0.054140 -> 0.054148 | Z `cone_family.H6.at_Bhat_ref.K2_ufine`, `K3_ufine_strip` |
| verifier name | fin12f4 -> fin13v4b (B1, provenance item G2, (W2), table, sources) | V4B_REPORT.md |
| flux hypothesis | "K0 <= FLUX0, K2 <= FLUX2" -> (I3'): pointwise mu(K0bar + eps_I K2pt) <= mu(FLUX0 + eps_I FLUX2) plus time-integrated part absorbed, C' <= 2.6524e-26; K0bar 3.4099e-21, K2pt 1.3901e-19, FLUX0/K0bar x4.8096 (not in text digits beyond 4.8096), combined x4.804; K2pt alone = 1.74 x FLUX2. Old (A)-chain margins 3.365/5.955 removed (record only, v36 values 3.368/5.961 mentioned once) | `fix8_flux2_impl/out/flux2_cert.json` keys `res.K0bar`, `res.K2pt`, `res.Cp`, `res.Cpp`, `(FLUX0+eps_I FLUX2)/(K0bar+eps_I K2pt)` |
| (K1) and table row for fluxes | rewritten to (I3') | same |
| sources/hashes block, S4 note, header, title/date, status box | v36 paths and hashes; S4 re-pin pending (stages 8, 9, 18, 19, 20) | V36_REPORT.md section 7 |
| Lemma D proof for k=2 | cites G1_export_k2.tex Theorem "Lemma D, k = 2", conditional on U-1 (U1_regularity.tex); no independent review | PROOFS_REPORT.md |
Numbers found equal after rounding (left unchanged on purpose): kappa_J total 0.042432, C^D_A, q^D, Sigma sweeps (5.0892, 0.30846, 1.5243, 3.0081, 4.335e-3), rate margins 28.231 / 4.432, 6.697, theta etc., K table, G-hat table, Gamma/kappa of Theta phi, M0z margins 12.99/14.37 (new 12.9969/14.3787, floor same), late M1 2.353, present M1 2.770, M0 2.225, N_k^2 outer ramp, Delta^R (4.4986e-8, 6.7972e-6, 0.011480), H6 margins 1.860/2.46/2.54/2.51/2.30, (H2)/(H2f) values (step 4, unchanged).

## 4. main_v3_v37: manual changes
| item | old -> new | source |
|---|---|---|
| Main Thm (i): abar | <= 2.19103e-7 -> <= 2.14041e-7 | MS `1.5e-4.sup_ab` |
| Main Thm (iv): tail rate | "= 2.69e-3" -> ">= 2.6776e-3" (equality sign changed to a bound, because the source is a lower bound) | RT15 `tail.rho` |
| (iv): floor | "= 9.872e-10" -> "<= 9.9001e-10" | RT15 `tail.sqrt_LF` |
| (iv): eventual alpha | "<= 6.867e-7" -> "<= 6.8843e-7"; tag fin12f4 -> fin13v4b | RT15 `tail.alpha_F` |
| Thm W (W1): sup abar, bbar | 2.192e-7 -> 2.14041e-7 | MS |
| (W1): sup beta | 3.34e-9 -> 2.35488e-9 | MS `1.5e-4.sup_beta` |
| Components, verifier: peak | 5.308e-8 -> 5.30857e-8 | RT15 `sqrtL_max` |
| tail rate / floor | 2.693e-3 -> 2.6776e-3; 9.872e-10 -> 9.9001e-10 (+ alpha_F 6.8843e-7, margins x177.07 / x175.09 added) | RT15, RT125 |
| J_4 | 6.24e-8 -> 6.25e-8 (J_{0,1,2} (6.03, 8.47, 1.88)e-12 unchanged after rounding) | RT15 `unc.V4sq` (JfourV 6.2473e-8) |
| Lemma S(A) item: sources, interface hash | fix2/n12/t2_A_v35, 78fe2f8e... -> fix6_v36/t2_A, 799b087c...; flux margins x3.365 / x5.955 -> (I3') values; thinnest margin now ramp M1 x1.184 | A key table; flux2_cert.json |
| Lemma U item: c_dU | 1.7187e-8 -> 1.46343e-8; proof-text pointer lemmaU_v34.tex -> lemmaU_v37.tex; certificates fix6_v36 | r3_forcing.json |
| macro-mapped numbers (section 2) | 14 occurrences: Bhat, Ehat_2, M1/M0 margins, G22 margin, ramp, (K6) 1.188, U-fine Bhat, S0 7.1860e-3 | numbers_A_v37 |
| New: status box (v37), open item 1 (status v37), appendix rows V4B/FX2, kit sentence, abstract of rel. work sentence | text only | |
Unchanged on purpose: all step-4 hand-off values (H1)-(H6), the Lemma O t* numbers, the theorem structure, Lemma U closure margins that did not change after rounding (loop gain 0.03904, rate x28.2, b1 thinnest x1.238, Delta^R_2 = 0.01148).

## 5. Conflicts between sources (listed, not silently resolved)
- C1. V36_REPORT section 6 prints C0_2 = 114.782, C1_2 = 0.0698142 (N2_ref) and 116.378, 0.260600 (4 N2_ref), rounded DOWN; numbers_A_v37 / G1_export_k2.tex give 114.783, 0.0698143 and 116.3782, 0.2606004 (rounded UP). The texts use the up-rounded values; the superseded quote is mentioned in lemmaU_v37 (Lemma D inputs).
- C2. FLUX2_IMPL_REPORT section 1 prints K2^pt = 1.3900e-19 and C' = 2.6523e-26; numbers_A_v37 (from flux2_cert.json, rounded up) gives 1.3901e-19 and 2.6524e-26. The texts use the numbers_A_v37 values (also the task's C' <= 2.6524e-26). D_2: V36_REPORT 0.0552610 vs numbers_A_v37 0.0552611 (up); text uses 0.0552611.
- C3. The D^c_2 window constant 0.0552107 exists only in V36_REPORT section 6 (not a key I could read in the JSON); quoted from the report, flagged in the text.
- C4. hjw_rep JSON key `kappa_J` includes the 1e-9 pad (0.0424318707645664); V36_REPORT quotes the certified kappa_J 0.042431869764566 without it. The text bound 0.042432 holds either way; the pad is noted.
- C5. Old text in main_v3 and lemmaU_v34 (not rewritten, flagged): main_v3 status boxes say "Steps 1, 2, 3 (v3.3) and 4 are complete as written" and "unconditional subject to ...", Theorem W says "The bootstrap closes (... unconditional subject to the pending review)", Section 1 of related work says "unconditional only subject to ...", and lemmaU says "Nothing in this item is open", "No mathematical gap ... was found", "(F1) is closed by this lemma". The sources do not certify step 3 as closed (open items in the v37 status box). I added the v37 status box and softened the three lemmaU sentences with explicit qualifiers; the main_v3 sentences were left and are covered by the new box ("read them subject to this box") -- decide whether to rewrite them.
- C6. lemmaU cites the step-4 `n1_chart.json` (341ed1f4...) for the chart term; V36_REPORT section 3 says chain v36 used a rerun on fin13v4b (`fix6_v36/n1/out/n1_chart_fin13v4b.json`, 7f64040f...). The "chart term moves every margin by at most 0.06%" sentence is from v35 and was not re-verified for v36.
- C7. lemmaU "Open items (O1)" cites `fix2/n35_asserts.py` and `fix1/lemmaU_rev/check_keys_v34.py` as checking the v35 rows; for v36 these have not been rerun (the v36 N35-type assertions are listed as values in V36_REPORT section 2).

## 6. TODOs
- T1. Regenerate `constants_needed_v34.json` and `check_keys_v34.py` for v36/v37 (the text now says the v34 list carries v35 values).
- T2. "Increase against v33 (3.89e-12, 4.91e-11)" in Lemma D inputs was not recomputed for v36 (flagged in text). Also "slope parts of D_k at most 2.7e-3 and 6.6e-4" (data_terms_late) not rechecked.
- T3. Read D^c_2 (k=2 window constant) and, for the key table, a few 3-digit quotes (e.g. "Delta^R,b") directly from JSON; I rounded Delta^{R,b} from the key table and did a JSON-diff scan for 4-6 digit tokens only. Tokens with <= 3 significant digits were scanned against `key_table_v35_r3_v36.json` only.
- T4. main_v3 TODO "source for 0.2686, x4.9, x5.3" (line 315 of main_v3, now ~330) untouched; erratum column (O3) not touched.
- T5. In lemmaU_v37 Lemma D the weight for k=2 is written as eps*omega = y^{-4}(28/y) on [ln 28, ln 350] as in G1; the sentence "the range ... is the range of the lower levels" is my gloss of the V36 report section 6 (b); check against G1_export_k2.tex before compiling.
- T6. Kits (`kit_step3_v3.zip`, `kit_step4b_v3.zip`) predate chain v36; sentence added to the main text; rebuild is Micah's decision.
- T7. Table row for U-fine H6 minimum "x1.86" etc. unchanged by the scan; thinnest (A) margin is ramp M1 x1.18446 (floor x1.184).

## 7. Checks (script, LaTeX not installed; ignoring `\{ \}` and %-comments)
| file | braces {/} | \begin / \end | mismatches | \ref without \label |
|---|---|---|---|---|
| lemmaU_v34.tex (old, for comparison) | 1771 / 1771 | 52 / 52 | none | none |
| lemmaU_v37.tex | 1947 / 1947 | 52 / 52 | none | none (36 labels) |
| main_v3.tex (old) | 909 / 909 | 45 / 45 | none | none |
| main_v3_v37.tex | 988 / 988 | 46 / 46 | none | none (35 labels) |
External (\cite and \dref) targets are bibliography keys / appendix labels and are not checked by this script; the new appendix labels \dref{V4B}, \dref{FX2} are defined only in the appendix table text, as for the existing \dref labels.
