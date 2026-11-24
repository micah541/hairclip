# v3.3 verifier rerun (tag fin12f4): FLUX0 = 1.64e-20, FLUX2 = 8.0e-20 (x4)

Nov 13, 2026. Copy `step3nl_v33/` of `step3nl/` (`cp -a`); the original is untouched (`diff -rq` shows only added files in the copy). Driver `verifier/run_v33.sh` (environment of `run_final.sh` + FLUX0/FLUX2, TAG=fin12f4, `nice -n 19`, sequential, logs in `verifier/`).

## Status
- tail_cert: **PASS** (rc 0; `tail_cert_log.txt` is byte-identical to the fin11 log, kept as `tail_cert_log_fin11.txt`; the tail_cert structural conditions do not depend on the floors).
- closure alpha = 1.5e-4: **success** (rc 0, tail=True, 10 steps, runtime 215 s).
- closure alpha = 1.25e-4 (remark point): **success** (rc 0).
- `mod_sups.py fin12f4 modsups_fin12f4` run (needs both alpha tables): sup|(a,b)|, sup|beta| below.
- Nothing moved by more than 1e-3 relative among the quantities that should not depend on FLUX (largest: slope-independent J constants 1e-6, rho 2.3e-6, peak sqrt(Lambda) 2.0e-6, M1_ramp 2.0e-6). The quantities that move are the flux-driven ones, as expected: F0 x4, slopes x4, sqrt_LF x2, alpha_F x2.

## fin11 vs fin12f4 at alpha1 = 1.5e-4 (certified, doubles as stored)
| quantity | fin11 | fin12f4 | ratio |
|---|---|---|---|
| J0 (V0sq const) | 6.0272312e-12 | 6.0272369e-12 | 1.0000009 |
| J1 | 8.4672025e-12 | 8.4672105e-12 | 1.0000009 |
| J2 | 1.8707790e-12 | 1.8707808e-12 | 1.0000009 |
| J4 (V4sq const) | 6.2361752e-08 | 6.2361811e-08 | 1.0000009 |
| slope0 | 8.636879e-19 | 3.454743e-18 | 4.0000 |
| slope1 | 1.213330e-18 | 4.853308e-18 | 4.0000 |
| slope2 | 2.680782e-19 | 1.072310e-18 | 4.0000 |
| sJ4 (V4sq slope) | 8.936291e-15 | 3.574508e-14 | 4.0000 |
| peak sqrt(Lambda) | 5.307447915e-08 | 5.307458343e-08 | 1.000002 |
| tail rho | 2.6930356e-03 | 2.6930419e-03 | 1.0000023 |
| tail a_star | 3.68488194e-05 | 3.68488194e-05 | 1.0000000 |
| tail L_star | 2.80644427e-15 | 2.80645530e-15 | 1.000004 |
| tail sqrt_LF | 4.935858e-10 | 9.871704e-10 | 2.0000 |
| tail alpha_F | 3.433270e-07 | 6.866519e-07 | 2.0000 |
| alpha-condition margin a_star/alpha_F | 107.33 | 53.66 | 0.5 |
| floor F0 | 3.28048e-22 | 1.312192e-21 | 4.0000 |
| X_slope (unc) | 5.361775e-22 | 2.144705e-21 | 4.0000 |
| sup\|(a,b)\| (abar), 6 digits | 2.19103e-07 | 2.19103e-07 | 1 |
| sup\|beta\|, 6 digits | 3.33811e-09 | 3.33811e-09 | 1 |
| A_max | 5.307447915e-08 | 5.307458343e-08 | 1.000002 |
| W1_max | 6.681316e-07 | 6.681329e-07 | 1.000002 |
| sI1, sI2 | 1.160335e-03, 1.545615e-02 | 1.160338e-03, 1.545618e-02 | 1.000002 |
| M1_ramp | 0.4327722 | 0.4327730 | 1.000002 |
| S_ramp | 7.473577e-04 | 7.473592e-04 | 1.000002 |
| floor terms | f1 5e-324, F0 3.28048e-22 | f1 5e-324, F0 1.312192e-21 | F0 x4 |

Remark point alpha1 = 1.25e-4: J0..J4 and all non-flux quantities agree to 1e-6; slopes x4 (V0sq 8.364763e-19 to 3.345899e-18, V1sq 1.175103e-18 to 4.700401e-18, V2sq 2.596320e-19 to 1.038526e-18, V4sq 8.654742e-15 to 3.461890e-14); sup|(a,b)| 2.15369e-07 to 2.15370e-07 (6-digit round-up), sup|beta| 3.20505e-09 unchanged; peak 5.21351020e-08 to 5.21351845e-08.

Note: fin11 reference values are read from `rate_table_fin11_*.json` and the committed `final_alpha_points_5e-8.json`. The "J" values above are the raw `unc` constants; the proof inputs are these rounded up (e.g. 6.0273e-12).

## sha256 of new outputs
```
b20e87153440144b0537afa4cd4b9dff2d788ffd2ce2163626faf973b13d2ad1  rate_table_fin12f4_a0.00015.json
79b60067f0f935eec8807a7641ba0af0b514755f7ac0ae553488f1a0e36fef65  rate_table_fin12f4_a0.000125.json
f028ed4a1275563de73dee145934f5ce49dac856da1582c5f53f8592246b9a60  verify_log_fin12f4_a0.00015.txt
303f0b38d06315a5c3ff9ecb2b736b09f751bf4094227370260989783ee627f2  verify_log_fin12f4_a0.000125.txt
d2237124555e0f64802797ce28a8ad5457427839eb933513d1e7a3aaaaaa945f  fin12f4_tail.out
f028ed4a1275563de73dee145934f5ce49dac856da1582c5f53f8592246b9a60  fin12f4_1.5e-4.out
303f0b38d06315a5c3ff9ecb2b736b09f751bf4094227370260989783ee627f2  fin12f4_1.25e-4.out
ee382a38550210f9ed8d56f54f8588c5d1193cd07d186162243a582fa6a75561  fin12f4_status.txt
d5b14797e932bb01c3ce2c4cd89c56fd994a42647ac6e4c0133fbb808d3cfde6  modsups_fin12f4/final_alpha_points_5e-8.json
6c7c158ef94df1552eed01f0202b77dcfe2584b81373ec6f912d5a0709671638  modsups_fin12f4/mod_sups.out
d2237124555e0f64802797ce28a8ad5457427839eb933513d1e7a3aaaaaa945f  tail_cert_log.txt
6dae054b69525a4cb6a6a176123277449e2953ce871cdd00cfae7fe4724f71b0  run_v33.sh
```
