# M14: fixes M1-M4 (REVIEW_v33) and G-A1 (B12_REPORT), chain rerun

Nov 14, 2026. Labels: [cert] arb, outward rounding. No stop trigger hit: the largest margin move is 0.55% (flux), every check passes (K1-K3, covering asserts, rule_reproduces_EA_ref).

## Fixes (backups `*_pre_m14.py`, and `*_pre_ga1.py` for the three G-A1 files)
- **M1** `t5_zone3/ramp_U.py`: K[j][j] now uses the primary Y_next read from `t4_cert/out/u1_ghat.json` (457.0916, 356.9396, 1320.2191), asserted to dominate the contract Y. foot_out is unchanged, so u1 does not change (no second u1 pass needed). Only K changed: K[0][0] 13435.567 -> 13436.024, K[1][1] 10456.94 -> 10526.18, K[2][2] 38901.06 -> 39482.29. Old file kept as `ramp_U.v34pre.json`.
- **M2** `t4_cert/assemble_U.py`: coupling term of the Ê rule restored (was NOT yet in the file; the scratch copy was only a diff target). Done together with G-A1.
- **M3** typed constants:
  - `t2_A_v34/rigor/assembly_rep.py`: ramp sups M1 and sI1 read from fin12f4 (alpha_1 = 1.5e-4 row), rounded up to 4 significant digits: 0.4328 (unchanged), sI1 0.00116 -> 0.001161 (verifier 0.0011603375).
  - `zone3.py`: abar = bbar read from the verifier (sup_ab, string 2.19103e-7, asserted equal to the old typed value); c_dU read from U_cert_v34.json launch.c_dU (certified 1.7186e-8, asserted < 2e-8) instead of typed 2e-8 [hyp]. New option `--ucert_full`.
  - a_W = 3/10: `u2_K.py` uses arb(3)/10 (degenerate ball, contains 3/10). `ramp_U.py` uses the box upper end `math.nextafter(0.3, inf)` (outward, contains 3/10); I did not use an arb ball there because bb_sup splits float boxes.
- **M4** u1_ghat.json and u34_primary.json regenerated with the current t4lib. All values byte-identical strings; only the args block changes (inputs now pin fin12f4 and the current t4lib). The abar record now points at the fin12f4 file.
- **G-A1** (R-sup of lemmaS_A_v33.tex), in `assembly_rep.py` (assemble), `export_t2.py` (rule `EA_of_N2`, plus new fields `A0_Rsup`, `Fg`) and `assemble_U.py`:
  A_j = E(tau1) + Dsrc(gamma') + (chat + oh)/gamma' (coupling now with coefficient 1, was 2), E_j = A_j + 2Fg_j^2 + 2 Fg_j sqrt(A_j + Fg_j^2), Fg = F/gamma'. The rule is nonlinear in o_j now, so the interface exports A at N2 = 0 and Fg, and assemble_U recomputes A_j(N2) then E_j. Note: the coupling coefficient 2 -> 1 follows lemmaS_A_v33.tex (A_j contains chat/gamma'). Its effect is below 1e-4 of the o/gamma' part.

## Chain order run
t2_A_v34/run_t2.sh (38 s) -> u1_ghat, u34 (identical) -> u2_K (primary, nophi, fallback; K, footw, R, L_max identical; c33 moves by 6e-13 relative) -> assemble pass 1 (zone3 v33 interface, as the previous pin) -> make_ucert -> zone3 -> assemble pass 2 (new zone3 interface) -> zone3 again -> assemble pass 3. **Pass 3 == pass 2 apart from args (fixed point).** Pass 1 differs from pass 2 (n_III, DeltaR moved slightly through M3/M1, expected). Old U_cert kept as `U_cert_v34_pre_m14.json`, pass files as `U_cert_v34_pass1/pass2.json`.

## Old -> new [cert]
- Ê supersolution (level 0, 1, 2): 1.195288e-6 -> 1.202035e-6 (+0.56%); 1.043274e-4 -> 1.043926e-4 (+0.06%); 0.05668273 -> 0.05668573 (+0.0053%). The G-A1 increase matches the expected 6.8e-9 / 6.6e-8 / 3.2e-6 (measured 6.7e-9, 6.5e-8, 3.0e-6). The M2 effect is below this.
- Ê_ref at N2_ref (A_interface): 1.19506e-6 -> 1.20181e-6; 1.04309e-4 -> 1.04374e-4; 0.0566785 -> 0.0566815.
- Flux margins (A): FLUX0 x3.4139 -> x3.3952, FLUX2 x6.0410 -> x6.0105 (beta_A = 2); beta = 5/4: x5.4332 -> x5.4035, x9.6164 -> x9.5678.
- S0v: x2.7783 -> x2.7738 (S0v vs 0.0153: 2.1254 -> 2.1220); S1v vs 0.6: 5.7893 -> 5.7883.
- Core M1: x4.4369 -> x4.4362; core M0 5.4586 -> 5.4554; Ghat22 x1.9174 (same); ramp M0 7.0213 -> 7.0203 (sI1 0.001161); ramp M1 1.18474 and ramp z1 1.20143 unchanged.
- Late M1 x2.3543 -> x2.3539; late M0 x4.8955 -> x4.8924; late M0z 13.013 -> 12.992. Present M1 min x2.7712 -> x2.7709; present M0 min x2.2251 -> x2.2251 (2.225143 -> 2.225136); M0z min 14.387 -> 14.377.
- Loop gain (total, levels 0, 1, 2): 1.2212e-7, 0.037169, 0.039034 -> 1.2251e-7, 0.037169, 0.039034.
- N2 (beta = 2): (0.0714724, 0.330112, 5.039987) -> (0.0714713, 0.330111, 5.040298).
- Delta^R (levels 0, 1, 2): 4.498520e-8, 6.795432e-6, 0.0114797984 -> 4.498521e-8, 6.797148e-6, 0.0114799909. Delta^R_2 moved by +1.9e-7 (as the reviewer predicted).
- n_III: (0.0714457, 0.327068, 3.044420) -> (0.0714447, 0.327066, 3.044413) (certified c_dU 1.7186e-8 instead of 2e-8).
- H6 at U-fine B-hat: min margin 1.8609 (no visible change); K3 cone margin 2.51080 -> 2.51071; slope margin 2.46300 -> 2.46299; t* margin 2.3019; box margin 1.2389.
- U-fine B-hat: (4.871949e-6, 2.626495e-4, 0.025835735) -> (4.872222e-6, 2.626579e-4, 0.025839794). Main B-hat: (4.873310e-6, 5.589784e-4, 0.36663939) -> (4.873482e-6, 5.592692e-4, 0.36669782). K1, K2, K3 true.

## Surprises
1. The coordinator's M2 premise ("the reviewer's scratch version") was only a copy with the fix, not applied to assemble_U.py; I applied it.
2. The G-A1 flux estimate (x3.382 -> x3.36) was relative to a different baseline; here x3.4139 -> x3.3952 (-0.55%).
3. The old zone3 pin used `zone3_interface.v33.json` for pass 1; I kept that for pass 1 and used the new interface from pass 2.

## New hashes (sha256)
- t4_cert/out/U_cert_v34.json 271430d7ffc0306fdb8b7d4942d1e98262f02ff1abde6706fd46334eb99b2507
- t4_cert/out/U_cert_for_zone3_v34.json 9feebbb2b672879f332858caf2ac209c29978a3b1df913e7bd780a4127e4c1cd
- t5_zone3/out/zone3_interface.json a381b1f22248a3adc62fc085086e170d75d1388dd086ea6fc1aff9de7ec704c0
- t5_zone3/out/ramp_U.json 55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de
- fix1/b0/t2_A_v34/out/A_interface.json 34832f434b7abd7e2dde14be0ee217d99d64236d07983f03863b4239a3ac4226
- t4_cert/out/u1_ghat.json 5a06d24cae8591f508f9f7499c7d2b4355e06cd5b9d8d2419a9a0108c5638466
- t4_cert/out/u34_primary.json abaf73724b84a0845b549d74a32050e5efafc8f246360ea2d8c388dd0dd70333
- t4_cert/out/u2_K.json 10ae3c383f7098284e5ae16f4355db9a2645f64531f62c6565c5824213f73261
- scripts: assemble_U.py d14baeea..., u2_K.py 059a2734..., zone3.py b9f23489..., ramp_U.py 0df7ce69..., assembly_rep.py 91dbd2be..., export_t2.py a395f048... (full lists in `t4_cert/out/SHA256SUMS_scripts_m14`, `t5_zone3/SHA256SUMS_scripts`, `t2_A_v34/out/SHA256SUMS_scripts_m14`).
- Regenerated: t4_cert/out/SHA256SUMS, t5_zone3/out/SHA256SUMS, t2_A_v34/out/SHA256SUMS (previous copies `*.pre_m14`).

Other kept backups: t2_A_v34/out_pre_m14 (outputs before M3/G-A1), out_m14_run1 (M3 only, no G-A1), u1/u2/u34 `*.pre_m14.json`, zone3_interface.pre_m14.json.
Not rerun: fallback-case assemble, u5_enum. Not touched: fix1/b12_text, step4.
