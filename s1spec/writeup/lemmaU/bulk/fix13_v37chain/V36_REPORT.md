# V36_REPORT: Lemma S(A) chain (t2 -> t4 -> zone III) on verifier fin13v4b + R3 + r_2 export (tag v36)

Nov 20, 2026. Approved by Micah (Nov 20). `N` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36`, `B` = `.../bulk`. Resumable log: `N/V36_LOG.md`.
Labels: [cert] = python-flint arb, outward rounding, inputs from raw JSON rounded outward; [float of cert] = a float quotient of certified values (the inequality itself is asserted in arb); [est] = estimate (as in v35).

## Verdict
**The v36 chain closes. No margin fails. Nothing threatened its bound, and no stop rule fired.**
- c_dU = 1.46343e-8 [cert], x1.367 under the typed bound 2e-8. It was x1.018 with R3 on fin12f4 and x1.164 in v35. The gain comes from sup|beta|, which fell from 3.338e-9 to 2.355e-9.
- Thinnest margins (v36): ramp M1 x1.18446, ramp z1 x1.20115, zone III box b1 x1.2389, H6 far K3 x1.8609, Ghat22 x1.9164.
- Largest relative loss on any margin against v35: x0.99948 (Ghat22, from r_2).
- r_2 export: C0_2 = 114.78, C1_2 = 0.0698 [cert]. The flux-lemma kill value is about 4e4 [est], so this has about x350 of room.
- Step 3 is **not** closed by this. The open items of V4B_REPORT and the flux lemma v2 (fix7_flux2) remain open.

## 1. What was built (all in N; nothing outside N was written)
| item | what |
|---|---|
| `t2_A/`, `t4_cert/`, `t5_zone3/` | Copies of the v35 scripts: `fix2/n12/t2_A_v35` (rigor, base_ibp, base_v31H2, prov, ref_phase2a, run_t2.sh), `B/t4_cert/*.py` and `B/t5_zone3/*.py`. `t5_zone3/out/ramp_U.json` is a byte copy of the v35 input; it was not rerun (see §6). |
| fix4_r3 patches | `paths_*.diff` and `r3_*.diff` were applied with `patch`. Each patched file was byte-identical to its fix4_r3 version before the v36 edits. Copies are in `N/patches/`. |
| `v36cfg.py` | The single source of the verifier inputs, plus three switches: `V36_VER` (fin13v4b, or fin12f4 for regression), `V36_R3` and `V36_R2`. It reads the raw verifier JSON and rounds outward: J_m, slopes, J4, sJ4 and the sqrt(Lambda) peak are rounded up to 5 significant digits, and tail rho is rounded down to 4. abar and betabar are the mod_sups strings, which mod_sups already rounds up to 6 digits. With fin12f4 it reproduces the v35 typed argv exactly; this is asserted. It also checks the rate table's own env: FLUX0 1.64e-20, FLUX2 8.0e-20, S0V 0.0153, S1V 0.6. These are the thresholds (A) is checked against. |
| `patches/apply_v36.py` | The v36 edits, each matched exactly once. Two later edits are logged in V36_LOG items 8 and 12 (the weights_t2 guards). The pad edit is item 11. |
| `n1/n1_chart.py` | A copy of step-4 stage 19, `VER = CFG.RT`. It writes only to `N/n1/out/`; the default destination into step 4 was removed. Rerun on the fin13v4b sups, because (A) consumes n1_chart.json and that file depends on the verifier's far_w, sI1, sI2 and S_ramp. |
| `run_chain_v36.sh` | Driver: n1 (fin13v4b) -> base_ibp reference chain -> r3 -> run_t2.sh -> u1, u34, u2_K -> assemble/zone3 passes 1-3 (same order as v35/fix4_r3) -> snapshot. The base_ibp chain (run_F10p_ibp, capture_F2, weights_exact, assembly_arb) is rerun every time, because `lemmaS_F10p_M0.52.json` args are the verifier inputs (A) consumes. |
| `cmp_runs.py`, `key_table_v36.py` | Read-only comparison tools. |
| Backups | `*_pre_v36.py` (14 scripts) and `run_t2_pre_v36.sh`, all byte copies of v35. `diff -u X_pre_v36.py X.py` gives the full change, R3 included. |

Gating: R3 is gated by `V36_R3` in assembly_rep, assembly_arb and u34_sweeps. r_2 is gated by `V36_R2` in weights_t2, assembly_rep.export and export_t2. The kappa_J pad is active only with fin13v4b. With the gates off the code paths are the v35 ones; regression 1 checks this.

## 2. Regression, double run, resources
| check | result |
|---|---|
| **Reg 1:** fin12f4, R3=0, R2=0 vs pinned v35 (`runs/reg_v35`, `cmp_vs_v35.json`) | **PASS.** Identical (args and timing ignored): A_interface, hjw, ic×4, a3, weights_t2 (+fb), ramp_rep×2, outer_L2, u2_K, u34, and all 5 base_ibp outputs (lemmaS_F10p, cap, rigcap, weights_exact, assembly_arb). The asm_×7, u1_ghat, U_cert passes, U_cert_for_zone3 and zone3 differ **only in path strings**. The base_ibp rerun with the generated argv reproduces the v35 files. |
| **Reg 2:** fin12f4, R3=1, R2=0 vs fix4_r3 `runs/r3run2` (`cmp_vs_r3run2.json`) | **PASS.** Every difference is a path string. r3_forcing.json is identical. |
| **Final-code regressions** (`runs/reg_v35_final`, `runs/reg_r3_final`) | See §2a. Reason for them: hjw_rep (kappa_J pad) and weights_t2 (r_2 grid) were edited after Reg 1 and 2. |
| **Production double run:** `runs/prod1` vs `runs/prod2` (fin13v4b, R3=1, R2=1) | **All 39 chain JSON outputs byte-identical, args included.** Logs differ only in timing lines. |
| Lemma U fixed point | U_cert pass 3 == pass 2 apart from args, as in v35. |
| ramp_U input | u1_ghat's primary Y_next (457.0916, 356.9396, 1320.2191) equals ramp_U.json's `Y_next_primary_used` exactly. ramp_U depends on nothing else that changed, because Y_next is abar-independent. So the v35 ramp_U.json is valid unchanged. |
| N35 assertions (not inside the scripts; values) | Ramp sups <= box: ramp M1 x1.18446, z1 x1.20115, M0 x7.02 (all > 1). Nbar²/(4 N2ref) <= 0.41264 (<= 1). |
| Resources | One process at a time, nice 19, PYTHONHASHSEED=0, python3 -B, no .pyc. Max RSS 93 MB. About 7 min per chain (base_ibp F2 + capture: 2 × 115 s). |
| Failed attempts (kept) | `runs/failed_reg_v35_assert`: the new r_2 assertion fired with R2=0. `runs/failed_prod1_r2grid`: on the 1/1000 Bc grid, the retained ball equalled r_2 up to rounding. Both fixed in the R2 branch only (V36_LOG items 8, 12). |

### 2a. Final-code regressions
Final code (with the kappa_J pad and the r_2 grid fix):
- `runs/reg_v35_final` (fin12f4, R3=0, R2=0) vs pinned v35: **PASS**, path strings only.
- `runs/reg_r3_final` (fin12f4, R3=1, R2=0) vs fix4_r3 r3run2: **PASS**, path strings only.
- The n1_chart copy run with fin12f4 (`n1/regtest/`) reproduces the step-4 `n1_chart.json` exactly (args ignored).
- The in-place outputs were then restored from `runs/prod2` by overwrite and checked byte-equal to the snapshot.

## 3. Inputs: fin12f4 vs fin13v4b (as consumed)
| input | fin12f4 (v35, R3) | fin13v4b (v36) | rounding |
|---|---|---|---|
| J0, J1, J2 | 6.0273e-12, 8.4673e-12, 1.8708e-12 | 6.0254e-12, 8.4681e-12, 1.8741e-12 | up, 5 digits |
| slopes | 3.4548e-18, 4.8534e-18, 1.0724e-18 | 3.4553e-18, 4.8561e-18, 1.0748e-18 | up |
| J4, sJ4 | 6.2362e-8, 3.5746e-14 | 6.2473e-8, 3.5826e-14 | up |
| sqrt(Lambda) peak | 5.3075e-8 | 5.3086e-8 | up |
| tail rho | 2.693e-3 | 2.677e-3 | down, 4 digits |
| abar = sup\|(a,b)\| | 2.19103e-7 | 2.14041e-7 | mod_sups up6 |
| betabar = sup\|beta\| | 3.33811e-9 | 2.35488e-9 | mod_sups up6 |
| ramp M1, sI1 (assembly_rep, up 4 digits) | 0.4328, 0.001161 | 0.4329, 0.001161 | up |
| kappa_J (hjw_rep) | 0.0423504 | 0.0424319 (certified 0.042431869764566 + pad 1e-9) | arb, fup |
| n1_chart (chart terms) | step-4 stage 19 output 341ed1f4... | `n1/out/n1_chart_fin13v4b.json` 7f64040f...: sup_add up by at most 1.0e-6 relative, beta22/32 and rho_meas by at most 6.6e-5 relative | [cert] |

## 4. Old (v35) vs R3 on fin12f4 vs v36: every named margin and interface quantity
The full generated table is `N/key_table_v35_r3_v36.md`, with JSON in `.json` (133 margin rows plus the interface rows). It comes from `key_table_v36.py fix4_r3/runs/base fix4_r3/runs/r3run2 runs/prod1`; `runs/base` has the v35 numbers (fix4_r3 R3_REPORT §4). Extract:

| quantity | v35 | R3 on fin12f4 | v36 | label |
|---|---|---|---|---|
| **c_dU** | 1.71860e-8 | 1.96498e-8 | **1.46343e-8** | [cert] |
| 2e-8 / c_dU (typed bound in u1_ghat M1z, zone3) | x1.16373 | x1.01782 | **x1.36665** | [float of cert]; `c_dU < 2e-8` asserted in arb |
| ramp M1 (tight, N2ref; all settings equal) | 1.18473 | 1.18473 | **1.18446** | [cert] |
| ramp z1 | 1.20143 | 1.20143 | 1.20115 | [cert] |
| ramp M0 | 7.02045 | 7.02045 | 7.02045 | [cert] |
| Ghat22 margin (153664/Ghat_22) | 1.91741 | 1.91741 | **1.91640** | [cert] |
| Ghat_22 core | 80141.6 | 80141.6 | 80183.8 | [cert] |
| FLUX0/K0 at beta=2 (tight, N2ref) | 3.39666 | 3.39535 | **3.40048** | [cert] |
| FLUX2/K2 at beta=2 (tight, N2ref) | 6.01283 | 6.01070 | **6.01860** | [cert] |
| FLUX0/K0 at beta=2 (tight, 4 N2ref) / (tight, hi) | 3.36511 / 2.92158 | 3.36382 / 2.92054 | 3.36888 / 2.92463 | [cert] |
| FLUX0/K0 at beta=2 (v31 H2, N2ref) | 2.88939 | 2.88838 | 2.89237 | [cert] |
| K0, K2 at E (tight, N2ref) | 2.43409e-21, 6.70519e-21 | 2.43502e-21, 6.70756e-21 | 2.43136e-21, 6.69879e-21 | [cert] |
| core M1 (tight, N2ref) | 4.43598 | 4.43592 | 4.43398 | [cert] |
| S0v vs 0.0153 / S1v vs 0.6 (tight, N2ref) | 2.12212 / 5.78795 | 2.12187 / 5.78788 | 2.12265 / 5.78501 | [cert] |
| windows: M1, M0, z2, \|f\|<=1/2, K2, K3 | 2.0516, 2.5438, 4.2684, 5.3881, 5.0804, 4.8245 | same | same (to 1e-6) | [cert] |
| **late M1** (Lemma U) | 2.35381 | 2.35377 | **2.35305** | [cert] |
| late M0 / late M0z | 4.89246 / 12.9929 | 4.89209 / 12.9904 | 4.89307 / 12.9969 | [cert] |
| **present M1 min** | 2.77080 | 2.77077 | **2.77017** | [cert] |
| present M0 min / M0z min | 2.22514 / 14.3774 | 2.22514 / 14.3763 | 2.22514 / 14.3787 | [cert] |
| rate margin (R2) / R<2 j=0,1 | 28.2313 / 4.432, 6.698 | same | same | [cert] |
| loop gain total j=1, 2 | 0.0371686, 0.0390342 | same | 0.0371688, 0.0390353 | [cert] |
| **H6 min** (= far K3, at Bhat_ref) | 1.86087 | 1.86087 | **1.86087** | [cert] |
| H6 cone K2 / K3 / slope (at Bhat_ref) | 2.54818 / 2.51066 / 2.46299 | same / 2.51065 / same | 2.54818 / 2.51046 / 2.46299 | [cert] |
| H6 far K2 / slope | 2.28566 / 2.85617 | same | same | [cert] |
| H6 "at fine est" rows | (as v35) | same | same | [est] |
| **zone III box** [253, 379.5]: b0, b1, bz | 0.100314, 0.419731, 0.0142902 | same | same | [cert] |
| zone III box margins b0, b1, bz (0.15/b0, 0.52/b1, 0.16/bz) | 1.4953, **1.23889**, 11.1965 | same | same | [float of cert] |
| zone III t* margin / DeltaR_2 | 2.30194 / 0.01148 | same | same | [cert] |
| zone III n_III k=0,1,2 | 0.0714447, 0.327066, 3.04441 | 0.0714456, 0.327068, 3.04442 | 0.0714437, 0.327064, 3.04441 | [cert] |
| M1z checks: zone III b1 + pU1 (vs 0.534) / S1(2 Bhat) | 0.50786 / 0.285844 | 0.50786 / 0.285849 | 0.50786 / 0.285945 | [cert] |
| Nbar²/(4 N2ref) max | 0.412639 | 0.412642 | 0.412637 | [cert] |

**Interface quantities pinned by step 4** (A_interface, U_cert, zone3 interface):

| quantity | v35 | R3 on fin12f4 | v36 | label |
|---|---|---|---|---|
| Ehat = EA_ref j=0,1,2 | 1.20150e-6, 1.04384e-4, 0.0566892 | 1.20198e-6, 1.04388e-4, 0.0566895 | 1.20011e-6, 1.04410e-4, 0.0567931 | [cert] |
| EA_ref4 j=0,1,2 | 1.21296e-6, 1.04511e-4, 0.0582464 | 1.21343e-6, 1.04516e-4, 0.0582466 | 1.21155e-6, 1.04538e-4, 0.0583511 | [cert] |
| EA0 j=0,1,2 | 1.19768e-6, 1.04341e-4, 0.0561702 | 1.19816e-6, 1.04346e-4, 0.0561705 | 1.19629e-6, 1.04367e-4, 0.0562737 | [cert] |
| F_j² j=0,1,2 | 5.65052e-13, 5.93194e-13, 6.22738e-13 | 6.46462e-13, 6.85051e-13, 7.38836e-13 | 3.59312e-13, 3.81978e-13, 4.15721e-13 | [cert] |
| C0_ref j=0,1 | 1.13376e-6, 1.66208e-4 | 1.13415e-6, 1.66212e-4 | 1.13260e-6, 1.66261e-4 | [cert] |
| C1_ref j=0,1 / C1_0 j=0,1 | 1.47907e-9, 5.15665e-8 / 7.20596e-10, 4.31427e-8 | 1.47962e-9, 5.15671e-8 / 7.21146e-10, 4.31433e-8 | 1.47761e-9, 5.15714e-8 / 7.19134e-10, 4.31476e-8 | [cert] |
| r_export | (1.249, 1.249, 0) | same | **(1.249, 1.249, 1/1000)** | exact |
| S0v, S1v | 0.00716755, 0.103578 | 0.00716838, 0.103579 | 0.00716575, 0.103631 | [cert] |
| U-fine Bhat, Bhat, EA (U) | see key table | | x0.9988 to x1.0018 | [cert] |
| Delta_in j=0..2 | 3.131e-12, 3.287e-12, 3.451e-12 | 4.093e-12, 4.297e-12, 4.511e-12 | 2.270e-12, 2.383e-12, 2.502e-12 | [cert] |
| Nbar² (beta=2) k=0..2 | 0.0714713, 0.330111, 5.04054 | 0.0714723, 0.330114, 5.04058 | 0.0714704, 0.330110, 5.04286 | [cert] |
| sha256 A_interface.json | 78fe2f8e... | d8d41fbe... | **799b087c...** | |
| sha256 U_cert pass 2 (Lemma U result) | 960809fd... (U_cert_v35) | 466ac45c... | **d512596a...** (`U_cert_v36_pass2.json`; pass 3 6d52881a..., equal apart from args) | |
| sha256 U_cert_for_zone3 | 22c13e44... | 6fced8ac... | **1344d5bf...** | |
| sha256 zone3_interface.json | 45e8d4dc... | 58bc8925... | **d6813ddb...** | |
| sha256 n1_chart (input) | 341ed1f4... (step 4) | same | 7f64040f... (N/n1/out) | |

Why the forcing fell although R3 is on: r0' is dominated by c_U·betabar, and betabar fell by x0.705. F_j² is x0.64 to x0.67 of v35, and the forcing share of E_j is x0.80 to x0.82.

Ehat_2 rose by x1.0018. The cause is the n1/rate-table data (J2, slopes) and r_2 (Ghat_22 up x1.00053). It is harmless: every margin that contains Ehat_2 moved by less than 0.06%.

## 5. kappa_J reserve and the raw rate-table consumers (flux-lemma v2 request, fix7_flux2/FLUX2_DESIGN.md C6/C7)
**(1) Reserve on the J's.**
- (A) uses J_m^used = J_m^5dig·(1 + kappa_J^used). Here J_m^5dig is the raw verifier constant rounded up to 5 digits, and kappa_J^used = kappa_J^cert + 1e-9. The 1e-9 is an absolute pad, added in hjw_rep (`v36cfg.KJ_PAD = 1/10^9`, production only).
- hjw_rep.json records `kappa_J_certified` = 0.042431869764566 [cert] and `kappa_J` = 0.042431870764566.
- The reserve over the certified J_m(1 + kappa_J^cert) is therefore at least 1e-9/1.0425 ≈ **9.6e-10 relative from the pad alone** [cert]. This covers C′/X_const ≤ 3e-10.
- The 5-digit rounding adds more, for the constants: J0 +3.6e-6, J1 +1.5e-6, J2 +4.4e-6, J4 +1.3e-5 relative [cert].
- The slopes are unchanged by C′ (FLUX2_DESIGN C6); they carry rounding reserves of +2.3e-6 or more anyway.
- Cost of the pad: no printed margin moved at 6 digits.

**(2) Every consumer in the chain that reads rate-table / mod_sups values.** Labels: "cushioned" = a Ly_hi-based quantity of the verifier, per ROUTE §5 C7 for fin12f4 (not re-read on fin13v4b by me); "raw" = a floor or tail double with only a one-ulp inflation.

| consumer | values read | treatment | reserve |
|---|---|---|---|
| `v36cfg.verifier_inputs` -> base_ibp argv -> `assembly_arb` (AS.J, SL, J4, SJ4, SQ, ABAR, BETA) | unc V0/1/2/4 const and slope, sqrtL_max, tail rho, mod_sups sup_ab, sup_beta | up 5 digits (rho down 4); mod_sups up6 | ≥ 1.5e-6 relative (J's); sqrtL_max +2.6e-6 |
| `hjw_rep` | Ly1, unc X_const, F0, steps Ly_hi and D, tail sqrt_LF, L_star, rho, rate_s_U | arb of the stored doubles. kappa_J = (max(Ly_hi − Ly1)/2 + F0·ΣD)/X_const, then fup, then + pad. W2: L_F <= L* | pad 1e-9 on kappa_J; Ly_hi cushioned |
| `assembly_arb.f4b` (module RT) | steps Ly_hi, sqrtL_next, D (UP), tail L_star (UP), rho (DN), **sqrt_LF (UP, one ulp)** | PHI_C, PHI_S -> I3c/I3s, the v3.1 page claims, and (B)-exit/v31-only rows | **PHI_S = raw sqrt_LF floor, no (1 + 1.2e-7) cushion.** Not consumed by Lemma U or zone III. It is exported only as records (`A_interface ... B0.flux_parts.I3*`, `v31_only`). If the flux lemma v2 consumes I3s, it must apply the (1 + 2C′/L_F)^{1/2} factor itself. |
| `assembly_rep.margins` | mod_sups M1, sI1 | up 4 digits | cushioned (mod_sups) |
| `t4lib.verifier_sups` -> u1_ghat (rate margins), u34 (c_dU), assemble_U | sup_ab, sup_beta | strings (mod_sups up6) | cushioned |
| `zone3.py` | sup_ab | string; assert <= 2.19103e-7 | cushioned |
| `r3_forcing` (record) | sup_ab, sup_beta; LAMT from constants.py (asserted equal to LAMT_Q) | strings | cushioned |
| `n1_chart` (copy) | far_max far_w, sI1, sI2, S_ramp_max | one ulp up | cushioned (far sups along the trajectory) |

Nothing else in t2/t4/t5 reads a verifier file. I checked by grep for step3nl, fin12f4, rate_table and modsups in the v36 scripts.

## 6. The r_2 export (ROUTE_DECISION.md unit G1 / F2′)
**Spec status.**
- ROUTE asks for "level-2 retained r_2 = 1/1000 on the core" and an export for k = 2 "with C0, C1, r", in the Lemma D window form.
- G1 (the written check that Lemma D extends to k = 2) has **not** been done; no `fix3/export2/` exists.
- So I implemented the natural certified quantity, the exact k = 2 analogue of the k = 0, 1 export. It is flagged for the flux-lemma designer below.

**Definition.** For every tau_a >= tau_1 and T >= 0,

  X_2(T) := ∫_{tau_a}^{tau_a+T} ∫_{ln 28}^{ln 350} y^{-4}(28/y) (z_4)² ds dtau <= C0_2 + C1_2·T,  with z_4 = d_s^4(f − U), s = ln y,

  C0_2 = (EA_2 + D_2)/(r_2 a_c),  C1_2 = (2F_2²/gp_2 + coup_2 + o_2 + D_2^slope)/(r_2 a_c),  r_2 = 1/1000,  a_c = (1 + 0.11²)^{-2}.

**Values [cert]** (`A_interface.json res.tight.r2_export`; doubles rounded up):

| setting | C0_2 | C1_2 |
|---|---|---|
| Nbar² = N2_ref (design) | **114.782** | **0.0698142** |
| Nbar² = 4 N2_ref | 116.378 | 0.260600 |
| Nbar = 0 (C1_0; then C1_2(Nbar) = C1_0 + o_2/(r_2 a_c), affine in Nbar²) | — | 0.00621885 |

The designer's guess was about 119 [est] and the kill value about 4e4 [est], so the room is about x350.

Window data:
- EA_2 = 0.0567931;
- D_2 = 0.0552610 (damped), with window constant 0.0552107 <= D_2 (asserted in `export`);
- D_2^slope = 6.0509e-6;
- F_2² = 4.157e-13, gp_2 >= 0.119547;
- coup_2 = 2.014e-8;
- rho_w = 5.50e-7.

**How it is certified.**
1. *Retained dissipation, weights_t2 at level 2.*
   - Core reference: Bc_core = floor((2 − r_2 − used/a_c − 1e-6)·10^6)/10^6 = 1688999/10^6, so the retained fraction is ≥ 0.001001 > r_2 [cert, asserted].
   - Zones [40, 60] (rise and S33 partition), segments (60, 175], transition [175, 350]: each cell's Bc is taken from avail = 2 − used/a − r_2·a_c/a_z, the same absolute-retain rule as at levels 0 and 1.
   - So r_2·a_c·∫Θ z_4² stays on the left of the level-2 energy inequality on all of [28, 350] [cert].
   - The ramp [14, 28] (Bc 137/100) and O [350, 700] are unchanged and not part of the export, as at k = 0, 1.
2. *Window integration (Lemma D form), assembly_rep.export.*
   - Every level-2 source in assembly_rep is a `Tr` object (damped value, window constant, window slope). This holds for the inflow, the ramp commutators, the IBP top term with its Young users, the class costs (pair, phi, ibp223), foot and S33. So the integrated inequality over any window gives r_2 a_c X_2(T) <= E_2(tau_a) + D_2^c + (2F_2²/gp_2 + coup_2 + o_2 + D_2^s) T.
   - E_2(tau_a) <= EA_2 (pointwise R-sup bound, valid at every tau_a >= tau_1).
   - The forcing enters as 2F²/gp per unit time, and the coupling as Σ_{k<2} Ghat_2k 28^{-4} E_k (pointwise E_k). Both are as at k = 0, 1.
   - Everything is arb. C0 uses D_2 >= D_2^c (asserted).
3. Two consecutive production runs are byte-identical.

**Cost.** Ghat_22 rises from 80141.6 to 80183.8 (x1.00053), so the Ghat22 margin goes from 1.91741 to 1.91640 [cert]. This is the only margin r_2 moves visibly. ramp M1 and K6 are unaffected, as ROUTE §7 predicted.

**FLAG for the flux-lemma designer.**
- (a) This is the k = 2 copy of the k = 0, 1 export. The G1 check that every level-2 source is window-affine was done here by construction (every source is a `Tr`, and the code asserts D^c <= D). There is no written Lemma D proof for k = 2; that text is still owed.
- (b) The export covers [ln 28, ln 350] with weight y^{-4}(28/y). The Agmon shell [54.75, 57.75] is inside it.
- (c) The variable is (A)'s z^A_4 = d_s^4(f − U), not the verifier's w′. Moving to g or w′ on the shell needs the chart/d terms at k = 4. n1_chart and b0_dprof give k <= 3 only (FLUX2_DESIGN stage X3).
- (d) Use C1_2 at the Nbar level of the bootstrap: 4 N2_ref gives 0.2606, and the window bound is affine in Nbar².

## 7. Step-4 stages and kit files that will need re-pinning (not done)
Step 4, manifest_v33 (`S4/p2/out_v33/FINAL_manifest_v33.json`, 21 stages, 8bc47e83...):

| stage | why |
|---|---|
| 8 `ratio_v33` (H7 time ratio) | reads fin12f4 modsups (abar) |
| 9 `notes_v33` | pins `t5_zone3/out/zone3_interface.json` |
| 18 `h4far_v33` | pins zone3_interface.json |
| 19 `n1_chart` | reads the fin12f4 rate table. Its v36 equivalent is `N/n1/out/n1_chart_fin13v4b.json` (from the copied script; the step-4 script needs a VER update) |
| 20 `f12_sliver` | re-derives the 7 asm runs, A_interface, U_cert and zone3 with the F1/F2-corrected constants (`out_v33/f12/{repro,new}`). It reads fin12f4 RT, modsups, step3nl_v33 constants.py and agent_sob, t2_A_v35/t4_cert/t5_zone3 and n1_chart. It must be repointed to fin13v4b, the v36 chain and R3/r_2 |
| driver | a new `final_tail_v36.py` (final_tail_v33.py is pinned). Per CLAUDE.md it must consult `EXCEPTION_v32_*.json` before any proof-tree rerun |
| manifest | a new manifest version. Stages 0-7 and 10-17 do not read the chain or the verifier, so their pins carry over |

That is 5 stages to re-pin (8, 9, 18, 19, 20), plus a new driver version and manifest. The closed FINAL_manifest (e01cb4da...) is untouched and needs nothing.

Kits:
- `kit_step3_v3.zip`:
  - `reference/lemmaSA/out/*` (A_interface, asm×7, hjw, weights_t2×2);
  - `reference/lemmaU/*` (U_cert_v35 + passes, U_cert_for_zone3, u1_ghat, u34, zone3_interface; ramp_U is unchanged);
  - `reference/verifier/*fin12f4*` (rate tables, verify logs);
  - `step3nl_v33/verifier/{rate_table_fin12f4, modsups_fin12f4}`;
  - `step4/p2/out_v33/n1_chart.json`;
  - the chain scripts (`t2_A_v35/rigor/*`, `t4_cert/*`, `t5_zone3/*`).
  - All are to be replaced by the fin13v4b verifier outputs and the v36 chain: `v36cfg.py`, r3 module, n1 copy, base_ibp argv.
- `kit_step4b_v3.zip`:
  - `reference/v33_out/p2/out_v33/f12/{new,repro}/*` (17 files) and `reference/v33_out/p2/out_v33/n1_chart.json`;
  - `step3nl_v33/verifier/{rate_table, modsups, constants.py}`, `step3nl_v33/agent_sob/sob_constants_P2.json` (if f12 moves to v4b's);
  - `step4/p2/v33/n1_chart.py`;
  - `writeup/lemmaU/bulk/fix2/n12/t2_A_v35/out/*`.
  - It also needs the new manifest and driver.

Texts (not edited; documents for later): lemmaS_A_v33.tex step 6 (R3, r_2 / Lemma D k = 2), numbers_A_v33.tex, lemmaU_v34.tex Lemma dU (c_dU 1.4634e-8), main_v3.tex, erratum, V33_INDEX successor.

## 8. Protected trees unchanged (sha256)
- step3nl_v4b, step3nl_v4, step3nl_v33 and the whole of `B` outside fix6_v36: no file present at the start changed or disappeared (`audit/tree_before.sha` vs `tree_after.sha`).
  - The only differences are 40 **added** files in `B/fix7_flux2` and `B/fix8_flux2_impl`. They were written by other agents (the flux-lemma v2 work) during this run, not by me.
  - `/home/claude/CLAUDE.md` was also modified by someone else.
- step4/p2/out, out_v33, v33: identical (1812 files).
- Manifests: closed FINAL_manifest e01cb4da..., manifest_v33 8bc47e83... (unchanged).
- `sha256sum -c` OK: v35 t2_A_v35/out, B/t4_cert/out, B/t5_zone3/out, fix4_r3/SHA256SUMS, step3nl_v4b/SHA256SUMS.
- No `__pycache__` or `.pyc` in N.
- Hash lists of the v36 outputs: `N/SHA256SUMS` (38 files; sha256 of the list d942c0ca...) and `N/SHA256SUMS_scripts` (100 files; 6b749a34...).

## 9. Files
- Outputs (production, in place = `runs/prod2`): `N/t2_A/out/` (A_interface.json etc.), `N/t2_A/base_ibp/{lemmaS_F10p_M0.52.json, rigor/*.json}`, `N/t4_cert/out/` (U_cert_v36_pass{1,2,3}, U_cert_for_zone3_v36[_pass1], u1_ghat, u2_K, u34_primary), `N/t5_zone3/out/` (zone3_interface[_pass1], ramp_U input), `N/r3/r3_forcing.json`, `N/n1/out/n1_chart_fin13v4b.json`.
- Hash lists: `N/SHA256SUMS` (outputs) and `N/SHA256SUMS_scripts`.
- Runs: `runs/{reg_v35, reg_r3, prod1, prod2, reg_v35_final, reg_r3_final, failed_*}`; `stages.txt` per run gives wall time and RSS.
