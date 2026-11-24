# Phase 1 checkpoint numbers (Nov 13, 2026)
Script `t1_checkpoint.py` (-> `.out`, `.json` with args block and input sha256). Runs in 1.5 s. Margins are [est] (float model) on [cert] inputs.
Self-check (t1 inputs): t1 sweep (4,4,0.3) KS=1 M1 x1.7248 = t1 json, rel dev 0; launch LATE M1 x2.1300, conservative x1.979, PRESENT rows: all rel dev 0; t1 with (H3) x1.666. PASSED (assert in script).

## Headline: weight (4,4,0.3), certified-input row (worst Ghat22 = 276917.6, Gw = 1.8199), with (2,2,3) placeholder [est, unbounded in proof]
| centre / variant | KS=1 M0 | KS=1 M1 | KS=2 M0 | KS=2 M1 |
|---|---|---|---|---|
| LATE, launch-lemma (I_W^in) | x3.83 | **x2.02** | x3.70 | x1.93 |
| LATE, conservative (full I_W ramp+foot) | x3.60 | x1.86 | x3.35 | x1.68 |
| PRESENT y_c<=126.7, launch | x1.96 (y_c=126.6) | x2.41 (y_c=56) | x1.96 | x2.27 |
| PRESENT y_c<=126.7, conservative | x1.76 | x1.42 | x1.72 | x1.34 |
| t1 sweep (full I_W, data + (H3)) | | x1.576 | | x1.468 |

Binding (smallest) number: PRESENT conservative, M1 x1.34 (KS=2); LATE launch M1 x1.93-2.02. M0 is never below x1.72. Excluding (2,2,3): identical to 3 digits (the S33-type term is ~1e-4 relative). Level-2 late launch KS=1: B2 = 0.485 = data 0.211 + inflow 0.215 + ramp 0.042 + foot 0.001 + outer 0.016.

## (4,4,0.0), worst Ghat22, with (2,2,3): LATE launch M1 x2.17 / x2.09 (KS 1/2), M0 x3.94 / x3.83; LATE cons M1 x2.00 / x1.81; PRESENT launch M1 x2.67 / x2.53, M0 x1.98; PRESENT cons M1 x1.52 / x1.43, M0 x1.79 / x1.75.

## Rows by Ghat22, (4,4,0.3), KS=1/2, LATE launch M1 / LATE cons M1 / PRESENT cons M1
| Ghat22 source | value | Gw | LATE launch M1 | LATE cons M1 | PRESENT cons M1 |
|---|---|---|---|---|---|
| P1 Umax=0.05 (info only) [cert] | 239750.3 | 1.639 | x2.12 / x2.03 | x1.96 / x1.77 | x1.49 / x1.40 |
| P1 Umax=0.09 [cert] | 250917.8 | 1.690 | x2.09 / x2.00 | x1.93 / x1.75 | x1.47 / x1.38 |
| P2 conventions [cert] | 276917.6 | 1.820 | x2.02 / x1.93 | x1.86 / x1.68 | x1.42 / x1.34 |
| WORST = max(P1 0.09, P2) | 276917.6 | 1.820 | same as P2 | | |

## Flip points (worst Ghat22, with (2,2,3), min(M0,M1) = 1, LATE; Ghat22 flip held with Gw recomputed and the foot recomputed)
| weight | variant | Ghat22 flip | Gw at flip | KS flip (absolute) |
|---|---|---|---|---|
| (4,4,0.3) | launch KS=1 / KS=2 | 5.30e5 / 5.23e5 (x1.91 / x1.89) | 7.29 / 6.70 | 35.2 |
| (4,4,0.3) | conservative KS=1 / KS=2 | 5.16e5 / 4.95e5 (x1.86 / x1.79) | 6.24 / 5.14 | 12.3 |
| (4,4,0.0) | launch KS=1 / KS=2 | 5.46e5 / 5.41e5 | 8.98 / 8.30 | 49.3 |
| (4,4,0.0) | conservative KS=1 / KS=2 | 5.34e5 / 5.16e5 | 7.58 / 6.21 | 15.3 |
Note: Ghat22 flip is far above 28^4 x 0.45; the rate criterion (Ghat22 < 28^4 = 614656) is the real limit, margin x2.22 for P2's value.

## Inputs
| input | label | source |
|---|---|---|
| E^A = (1.211005e-7, 8.174175e-5, 0.0480448585), gamma' = (0.2463, 0.2462, 0.12355) | [cert] | gonogo/lemmaS_A_ibp/rigor/assembly_arb.json |
| Ghat22 P1 (239750.3 / 250917.8), P2 (276917.6), worst | [cert] | phase1/ghat/ghat_cert.json (cases, sensitivity_Umax_0.09); phase1/ramp/ramp_K_cert.json foot['2'] |
| K_jk incl. new K[1][3] = 0.5716 (used with X_3 = gamma'_1 E_1/(2 a_min) + eps28 p^6 dU_0, as level 2 uses) | [cert] | ramp_K_cert.json ["4.0"]["K"] |
| K_jj rescale | none needed | P2's K_jj use P2's Y (the largest, Ghat 2.77e5); every row's Ghat is <=, so conservative |
| footw 1.9407, 2068.74 (no extra rr factor) | [cert] | ramp_K_cert.json foot |
| S33 coefficient 99.216 = 26.917 + 72.299 (2,2,2); times 28^2 for t1's normalisation (t1: 1.17e4*rr = 41*28^2 within 0.1%) | [cert] x [est] normalisation | ramp_K_cert.json foot['2'], cubic_remaining |
| (2,2,3) placeholder: 7.328 * C_eta/3 = 32.6 added to the S33 coefficient (C_eta^Xi = 13.349) | **[est, unbounded in proof]** | ramp_K_cert.json cubic_remaining |
| H2 = (1e-8, 1e-6, 1e-3), H3 = (2e-7, 1e-5, 1e-5) | [paper] | LAUNCH.md / launch_model.py (typed) |
| weights, IW, IW_in, w_min, C_A, q | [est] | float quadrature on the weight, as launch_model |

## Not wired in (still typed or missing)
- H2, H3, bbar, abar, alpha_1 = 1.5e-4, DT6 (Lemma R outer ramp, from t5_t6_fine.out): typed as in t1/launch_model; DT6 not recomputed with the Phase-1 Gw/weights. Also dU uses |alpha| <= alpha_1 = 1.5e-4, not step 4's 1.1739e-4.
- (2,2,3) is only a placeholder (no proof); (2,2,2) enters as the S33 term only. Possible (A) v3.1 omission (cutoff x jet terms on [14,28], P2 problem 2): not modelled.
- K was certified for g in [3e5, 6e5] (note text says 4.5e5), Gw uses the Ghat fixed point formula 1/(1 - Ghat/28^4) without re-deriving g; P2's Ghat is the worst of two inconsistent convention sets, not a reconciled certified value.
- The new Ĝ is not used in K_jj Y (see above), M0z/Umax of P1 vs P2 not reconciled; the transport term of the moving weight (P1 open) is absent from the model.
- Shell-model neighbour term uses q, Dd from [est] grid sups (+0.1%); (A)'s ramp/outer terms for the late centre unchanged from launch_model.
- t3 Young-route addition to E2^A (3e-5) not re-run.
