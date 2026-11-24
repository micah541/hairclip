# v3.3 rerun of the task-2 (A) chain with the fin12f4 verifier outputs

Nov 13, 2026. New files only (`t2_A_v33/`, copy of `t2_A/`; originals untouched, original `out/SHA256SUMS` still verifies). Verifier: `step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json` (sha256 b20e8715...), see `step3nl_v33/V33_VERIFIER.md`. Labels: cert = arb with outward rounding; hyp = hypothesis. Driver: `base_ibp/run_chain_v33.sh` (4 steps, 2.5 min), then `run_t2.sh` (40 s). **Two consecutive `run_t2.sh` runs are byte-identical** (19 JSON, `out/SHA256SUMS`). Old outputs kept in `out_pre_v33/`.

## Changed inputs (every scripts backed up as `*_pre_v33.*`)
| item | file | old | new |
|---|---|---|---|
| tree root | `rigor/arlib.py` ROOT, `run_t2.sh` cd | `bulk/t2_A` | `bulk/t2_A_v33` (writes only into the copy) |
| verifier table (H-J-W, F0, L_F, L*, rho, steps) | `rigor/hjw_rep.py` RT_DEF, `rigor/assembly_rep.py` provenance, `base_ibp/rigor/assembly_arb.py:161` | step3nl/.../rate_table_fin11_a0.00015.json | step3nl_v33/.../rate_table_fin12f4_a0.00015.json |
| slopes (3) | `base_ibp/run_F10p_ibp.py`, `base_ibp/rigor/capture_F2.py` argv | 8.6369e-19, 1.2134e-18, 2.6808e-19 | 3.4548e-18, 4.8534e-18, 1.0724e-18 (raw 3.45474e-18, 4.85331e-18, 1.07231e-18, rounded up) |
| sJ4 | same | 8.9363e-15 | 3.5746e-14 (raw 3.57451e-14) |
| J0,J1,J2,J4 | same | 6.0273e-12, 8.4673e-12, 1.8708e-12, 6.2362e-8 | unchanged (still >= fin12f4 raw values; hjw_rep asserts it) |
| sqLam, rho, abar, beta, alpha1 | same | 5.3075e-8, 2.693e-3, 2.19103e-7, 3.33811e-9, 1.5e-4 | unchanged (fin12f4: peak 5.30746e-8, rho 2.693042e-3, sup_ab/sup_beta equal at 6 digits) |
| FLUX0, FLUX2 in the margin/claims | `rigor/assembly_rep.py` (claims l.170, FLX l.220, margins l.225), `rigor/export_t2.py:62` | 4.1e-21, 2.0e-20 | 1.64e-20, 8.0e-20 |
| (H2f) design | `run_t2.sh` (`ic_rep.py --h2f a,b,c --name design`) | (5e-7, 5e-6, 1e-5) | **(1e-6, 5e-6, 1e-5)** |
| (H2f) lo, hi (+-20% of first entry, as before) | `run_t2.sh` | 4e-7 / 6e-7 | 8e-7 / 1.2e-6 (probe unchanged) |
| beta_A | assembly_rep `--betas` | 5/4 and 2 reported | unchanged; both reported, beta_A = 2 below |

Regenerated in `base_ibp/` (backups `*_pre_v33.*`): `lemmaS_F10p_M0.52.json`, `rigor/cap_F10p.json`, `rigor/weights_exact.json`, `rigor/assembly_arb.json/.out`; logs `F2_v33.log`, `rigor_*_v33.log`. `base_v31H2/` is read only for comparison columns and was not regenerated (it still carries fin11 inputs).

## Headline (design (H2f) = (1e-6, 5e-6, 1e-5), N2_ref = (0.1, 0.2, 5.0) [est]) [cert unless noted]
| | tight (H2) = (1e-8, 1e-6, 1e-3) | v3.1 (H2) = (2e-7, 1e-4, 0.2) |
|---|---|---|
| Ê (N2_ref) | 1.18766e-6, 1.03883e-4, 0.056585 | 1.37766e-6, 2.02883e-4, 0.255585 |
| EA0 (no outer) | 1.18385e-6, 1.03840e-4, 0.056066 | 1.37385e-6, 2.02840e-4, 0.255066 |
| Ê at 4 x N2_ref | 1.19908e-6, 1.04011e-4, 0.058142 | 1.38908e-6, 2.03011e-4, 0.257142 |
| gamma' (lower) | 0.242926, 0.241065, 0.119615 | same |
| C0 (z2, z3) | 1.11659e-6, 1.65440e-4 | 1.27241e-6, 2.46633e-4 |
| C1_0 | 3.894e-12, 4.910e-11 | 3.894e-12, 4.925e-11 |
| C1_ref | 7.624e-10, 8.473e-9 | 7.624e-10, 8.473e-9 |
| K0(Ê), K2(Ê) | 2.339e-21, 6.456e-21 | 2.758e-21, 8.605e-21 |
| FLUX0/K0, FLUX2/K2 at beta_A = 2 (target >= 2.0) | **x3.505, x6.196** | **x2.973, x4.649** |
| same at beta = 5/4 | x5.608, x9.914 | x4.757, x7.438 |
| FLUX0_needed, FLUX2_needed (4K) | 9.36e-21, 2.58e-20 | 1.10e-20, 3.44e-20 |
| S0v vs 0.0153 (target >= 2.0) | **x2.144** | x1.761 (note only) |
| core M1 (vs 0.52), M0 (vs 0.11) | **x4.443**, x5.464 | x2.723, x5.072 |
| core Ghat_22 (153664) | 80141.6, x1.917 | same |

The stop level (flux margin < 1.5 at beta_A = 2) is not approached. beta_A = 2 holds for Lemma S(A) at both (H2) settings.

(H2f) settings at the tight (H2), beta_A = 2 (FLUX0/K0, FLUX2/K2) and S0v margin: lo (8e-7): see `out/A_interface.json` H2f_settings; hi (1.2e-6): S0v x2.056 (still >= 2.0); design x2.144; lo x2.251; probe (3.5e-7, 1.1e-6, 2.7e-6) x2.666. Ê_0 = 0.99e-6 (lo), 1.19e-6 (design), 1.39e-6 (hi), 0.54e-6 (probe). FLUX0_needed: 7.8e-21 / 9.4e-21 / 1.09e-20 / 4.3e-21, all below 1.64e-20 (so beta = 2 margins > 1.5 at every setting).

## Acceptance (CONTRACTS_r2 rows), design setting, tight (H2)
| quantity | threshold | value | status |
|---|---|---|---|
| FLUX0/(beta_A K0), FLUX2/(beta_A K2), beta_A = 2 | >= 2.0 (stop < 1.5) | x3.505, x6.196 | pass |
| S0v against 0.0153 | >= 2.0 | x2.144 | pass |
| core M1 | >= 2.0 | x4.443 | pass |
| E2A at N2_ref | <= 0.1 (report > 0.08) | 0.05658 | pass |
| o0/gamma'0 at 4 x N2_ref | <= 1e-7 | 1.52e-8 (3.8e-9 at N2_ref) | pass |
| R5' addition to E2A at 4 x N2_ref | <= 0.01 | 2.08e-3 | pass |
| every class on O closed by T1-T7, P | yes | yes (asserted in outer_L2) | pass |
| H-J-W kappa_J | <= 0.13 | 0.0423504 (W1 0.0423485 + W3 1.92e-6) | pass |
| core Ghat_22 margin | >= 1.5 | x1.917 | pass |
| retained r0, r1 on [28,350] | >= 1.249 a_c | 1249/1000 | pass |
| segment conditions, both boxes | all true | true for (0.15,0.52) and (0.20,0.60) | pass |
| reruns | byte-identical | yes (two runs, 19 JSON) | pass |

## Movements vs the old t2_A run (for quantities that should not depend on FLUX)
gamma', Ghat, ramp tables, weights, outer kappa, segments, r, core M0/Ghat_22: identical. kappa_J 0.042348 to 0.0423504 (+5e-5 relative; W3 part grows with F0 x4, as it should). Ê, C0, S0v, M1core changed because (H2f) changed (design 5e-7 to 1e-6), not because of FLUX (old design Ê_0 6.877e-7 vs new 1.1877e-6, difference 5.0e-7 = the (H2f) shift plus 3.8e-9 outer). C1_0 level 1: 1.56e-11 to 4.91e-11 and level 0: 3.83e-12 to 3.89e-12, driven by the x4 slope ratio rho_J (expected); C1_ref moved <1% (7.62e-10 unchanged, 8.44e-9 to 8.47e-9). No unexpected movement; nothing stopped.

## Caveats
- (H2f) remains [hyp] (task 6 must certify it on step-4 data); the outer data N2 are [est].
- `base_ibp` rigor_assembly_arb_v33.log lists 5 "QUOTE TOO TIGHT" legacy page quotes versus 3 before (gamma'_2, gamma'_B, Ghat_22 core as before; new: I3 slope page 1.39e-27 -> 2.15e-27 and Phi slope page 4.94e-10 -> 9.87e-10). Both new ones are slope terms and scale with the x4 flux slope; they are v3.1 page quotes, not task-2 acceptance rows. They need updating in any v3.3 text that quotes them.
- `v31_only.claims_failing` in A_interface.json is the legacy v3.1 quote list (informational, as in r2); `base_v31H2` comparison columns use fin11 inputs.
- The verifier closure and the chain used 1.5e-4 only for task 2; the 1.25e-4 remark point was run in the verifier but not propagated.

## Files
`t2_A_v33/`: `V33_A_REPORT.md`, `out/A_interface.json` (+ 18 stage JSONs, logs, `SHA256SUMS`), `out_pre_v33/`, `run_t2.sh`, `run_v33.log`, `run_v33_2.log`, `base_ibp/run_chain_v33.sh`, `base_ibp/chain_v33.out`, regenerated base_ibp outputs listed above, `rigor/{arlib,hjw_rep,assembly_rep,export_t2}.py` (+ `*_pre_v33.py`).
