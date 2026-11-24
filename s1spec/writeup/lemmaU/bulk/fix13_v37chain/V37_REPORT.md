# V37_REPORT: verifier fin13v4c → (A) chain v37 → flux lemma v2 reassembly (fix round 13)

Nov 20, 2026. Trigger: `fix11_review/REVIEW_v4b_v36_proofs.md` (D-3 blocking; A-1, D-2a, A-2). Approved by Micah.
- `N` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain`; `C` = `/home/claude/s1spec/step3nl_v4c`.
- Logs: `N/V37_LOG.md`, `C/V4C_LOG.md`. Verifier report: `C/V4C_REPORT.md`.
- Labels: [cert] arb with outward rounding; [float of cert] float quotient of certified values.

## Verdict
**Everything closes.** No margin fails and no stop rule fired:
- J moves by −0.016% (limit 1%);
- C0_2 = 116.359 at 4·N2_ref (limit 1e3);
- c_dU = 1.463427e-8, ×1.36665 under 2e-8.

The blocking D-3 item is replaced by a certified term (H8). Every downstream margin is equal or slightly larger, because the certified strip cubic term is smaller than the old 2.6e-3 scaling in the states that matter.

**Step 3 is still not closed.** The other items of the review's §5 remain: (GS), C-1, the ei_lemma_v4 review, the DS/Sobolev audits, records B1–B7 and D-1a, texts, re-pin and kits, the LaTeX compile, and the final review.

## Checks
| check | result |
|---|---|
| v4c regression (V4C=0, TAG fin13v4b) vs `step3nl_v4b/runs/P` | **PASS**: rate tables, verify logs, final_alpha_points and tail_cert log are equal, ignoring timing and args |
| v4c double run (`C/runs/P_A` vs `C/runs/P`) and cert chain (`out_A` vs `out`) | **PASS**. The cert JSON shows the same runtime-only pattern as v4b; record note in V4C_REPORT |
| v4c tail_cert, box asserts | PASS; all box ratios ≤ 1. New D-2a: CF_w1_16_31 0.116 (×8.6) |
| v37 regression: chain with fin13v4b inputs (`N/runs/reg_v36`) vs `fix6_v36/runs/prod2` | **PASS**: path strings only (`runs/cmp_reg_vs_v36prod2.txt`); all 215 key-table rows are equal |
| v37 double run `runs/prod1` vs `runs/prod2` (fin13v4c) | **all 39 JSON byte-identical, args included** |
| Flux regression (`F2_VTAG=fin13v4b`, inputs fix6_v36 prod2) vs `fix8_flux2_impl/out` | **PASS**: all 8 JSON equal apart from args/timing |
| Flux double run (`flux/snap/v37_pass1` vs `v37_pass2`) | all 8 JSON and the pickle are byte-identical |
| R3, r_2 = 1/1000 and κ_J pad 1e-9 | kept (V36_R3 = V36_R2 = 1; pad active for fin13v4c) |
| Resources | one process at a time, nice 19, PYTHONHASHSEED=0, python3 -B; peak RSS 271 MB |

## Old vs new
### Verifier (fin13v4b → fin13v4c), α1 = 1.5e-4 [cert]
| quantity | v4b | v4c | ratio |
|---|---|---|---|
| peak √Λ | 5.308562e-8 | 5.308499e-8 | 0.999988 |
| tail ρ | 2.677666e-3 | 2.679225e-3 | 1.00058 |
| J0/J1/J2/J4 const | 6.025378e-12 / 8.468087e-12 / 1.874092e-12 / 6.247218e-8 | 6.024419e-12 / 8.466739e-12 / 1.873793e-12 / 6.246224e-8 | 0.999841 |
| slopes | 3.455292e-18 / … | 3.454540e-18 / … | 0.999782 |
| sup\|(a,b)\|, sup\|β\| | 2.14041e-7, 2.35488e-9 | 2.14039e-7, 2.35488e-9 | |
| whole-tail α-condition (ratio, A-2) | ×1554.75 (1.25e-4: ×1512.43) | ×1554.46 (×1512.16) | |
| E0 A² absorption (A-1 added) | ×1.679 as written (×1.657 true) | **×1.65695** [cert] | |
| strip cubic | 2.6e-3 × far-region term (not a bound) | H8 with certified C3S_abc | |

At α1 = 1.25e-4, ρ goes 2.849387e-3 → 2.850927e-3 and J changes by ×0.999861.

### Inputs consumed by (A) (`v36cfg.verifier_inputs`, outward-rounded strings)
| | v36 | v37 |
|---|---|---|
| J0, J1, J2 | 6.0254e-12, 8.4681e-12, 1.8741e-12 | 6.0245e-12, 8.4668e-12, 1.8738e-12 |
| slopes | 3.4553e-18, 4.8561e-18, 1.0748e-18 | 3.4546e-18, 4.8551e-18, 1.0745e-18 |
| J4, sJ4 | 6.2473e-8, 3.5826e-14 | 6.2463e-8, 3.5818e-14 |
| √Λ peak, ρ | 5.3086e-8, 2.677e-3 | 5.3085e-8, 2.679e-3 |
| abar, betabar | 2.14041e-7, 2.35488e-9 | 2.14039e-7, 2.35488e-9 |
| κ_J (certified + pad) | 0.0424319 | 0.0424297 |
| n1_chart sup_add / e | (v36) | equal to 1e-7 relative |

### (A) chain, Lemma U, zone III (v36 → v37) [cert]
The full 215-row table is `N/key_table_v36_v37.md` and `.json`. Its columns are v36 prod2 | v37 regression (= v36) | v37 prod2; read the "v35" header as v36 and the "R3" header as reg_v36.

| quantity | v36 | v37 | ratio |
|---|---|---|---|
| **c_dU** | 1.4634275e-8 | **1.4634274e-8** | 1.000000 |
| 2e-8 / c_dU | ×1.36665 | ×1.36665 | |
| ramp M1 (thinnest) | 1.18446 | 1.18446 | 1.000000 |
| ramp z1 | 1.20115 | 1.20115 | 1 |
| zone III box b1 | 1.23889 | 1.23889 | 1 |
| H6 far K3 | 1.86087 | 1.86087 | 1 |
| Ghat22 | 1.91640 | 1.91640 | 1 |
| core M1 (tight, N2ref) | 4.43398 | 4.43429 | 1.000068 |
| S0v vs 0.0153 / S1v vs 0.6 | 2.12265 / 5.78501 | 2.12273 / 5.78546 | |
| FLUX0/K0, FLUX2/K2 at β=2 (tight, N2ref) | 3.40048 / 6.01860 | 3.40057 / 6.01881 | 1.00002 |
| late / present M1 (Lemma U) | 2.35305 / 2.77017 | 2.35321 / 2.77032 | 1.00007 |
| present M0 min | 2.22514 | 2.22514 | 1 |
| rate_margin (Lemma U) | 28.2313 | 28.2313 | 1 |
| Nbar²/(4 N2ref), max | 0.412637 | 0.412637 | ≤ 1 |
| **C0_2 / C1_2 (N2ref)** | 114.782 / 0.0698142 | **114.763 / 0.0698137** | 0.99984 |
| **C0_2 / C1_2 (4 N2ref, used by the flux lemma)** | 116.378 / 0.260600 | **116.359 / 0.260600** | 0.99984 |
| E_j (A) j = 0, 1, 2 | 1.20011e-6, 1.04410e-4, 5.67931e-2 | 1.20008e-6, 1.04394e-4, 5.67838e-2 | ≤ 1 |
| K0, K2 at β=2 | 4.82285e-21, 1.32921e-20 | 4.82273e-21, 1.32917e-20 | |

The smallest ratio over all 133 margin rows is 1.000000, so no margin decreased.

### Flux lemma v2 (fix8 → v37 reassembly, `N/flux/out/flux2_cert.json`, verdict PASS) [cert]
| quantity | fix8 (v36) | v37 |
|---|---|---|
| K̄0 | 3.40985e-21 | 3.40977e-21 |
| K2^pt | 1.39001e-19 | 1.38995e-19 |
| FLUX0/K̄0 | ×4.80959 | ×4.80971 |
| **combined (FLUX0 + ε_I FLUX2)/(K̄0 + ε_I K2^pt)** | **×4.80442** | **×4.80453** |
| C′ / S′ | 2.65232e-26 / 1.05510e-28 | 2.65189e-26 / 1.05502e-28 |
| slack s (C0) | 1.03907e-21 (0.7919 F0) | 1.03908e-21 (0.7919 F0) |
| C1 | ×9.85e6 | ×9.85e6 |
| C2 (1.25e-4 / 1.5e-4) | ×5.00e7 / ×5.34e7 | ×5.00e7 / ×5.34e7 |
| **strictness pad C3** (1.5e-4 / 1.25e-4) | ×1547.2 / ×2235.8 | **×1547.6 / ×2236.4** |
| C4 | ×1.85e5 / ×1.74e5 | ×1.85e5 / ×1.74e5 |
| C5 (whole-tail α with 2C′) | ×1554.75 / ×1512.43 | ×1554.46 / ×1512.16 |
| C6 (κ_J reserve) | ×141.32 | ×141.32 |
| C7 mod_sups / tail | ×2.5e8 / ×5.1e7 | ×2.5e8 / ×5.1e7 |
| C8 (U-1) | proved [hand] in fix10 (not scripted) | unchanged |

FLUX0 = 1.64e-20 and FLUX2 = 8.0e-20 are unchanged. The verifier was run with them, so no further verifier pass is needed.

## Step-4 / kit re-pin list (not done; needs Micah's go)
**Confirmed: the same stages as V36_REPORT §7, with fin13v4c and v37 replacing fin13v4b and v36.** The verifier change touches only the rate tables, mod_sups and the cert constants. The chain outputs that step 4 reads (zone3_interface, A_interface, U_cert, n1_chart) are the same files.

Step-4 stages:
- **8 ratio_v33:** modsups, abar 2.14039e-7.
- **9 notes_v33:** zone3_interface.
- **18 h4far_v33:** zone3_interface.
- **19 n1_chart:** the rate table; v37 equivalent `N/n1/out/n1_chart_fin13v4c.json`.
- **20 f12_sliver:** the rate table, modsups and constants.py, now step3nl_v4c, which also needs `v4bcert/out/strip_cubic_v4c.json`; plus the v37 chain and n1.

Also needed:
- the new driver `final_tail_v37` (must consult `EXCEPTION_v32_*.json`) and a new manifest;
- flux outputs `N/flux/out/*` (8 JSON + pickle) and texts `N/flux/texts/*`;
- review items B1 (base_ibp outputs, no args) and B3 (ramp_U hash). ramp_U is still the v35 byte copy, and Y_next is unchanged.

The closed FINAL_manifest needs nothing.

Kits: the same files as in V36_REPORT §7. In kit_step3_v3 the reference verifier becomes `step3nl_v4c/runs/P` plus `v4bcert/out/{strip_cubic_v4c,e0_A2_budget_v4c}.json`, and the chain becomes `fix13_v37chain`. In kit_step4b_v3, f12 is repointed to v4c/v37.

## Nothing pre-existing changed (sha256, snapshots in `C/audit_v4c/`, taken before any write)
| tree | files | result |
|---|---|---|
| step3nl_v4b | 652 | all unchanged, no new file |
| step3nl_v4 | 589 | unchanged |
| step3nl_v33 | 560 | unchanged |
| fix6_v36 | 1039 | unchanged, no new file |
| fix8_flux2_impl | 243 | unchanged, no new file |
| fix10_proofs | 10 | unchanged |
| fix12_texts | 29 | unchanged |
| fix11_review | 20 | unchanged |
| FINAL_manifest.json (closed) | e01cb4da… | OK |
| FINAL_manifest_v33.json | 8bc47e83… | OK |

- No file under step3, step3nl, s1_tight, step4 or the kits is newer than the start marker.
- No `__pycache__` was written anywhere under s1spec.
- Step 4, the kits and the manifests were not touched.

## Files
- `N/run_v37_all.sh` and `run_chain_v37.sh` (drivers); `v36cfg.py` (fin13v4c added; backup `_pre_v37`); `cmp_runs.py`.
- Runs: `runs/{reg_v36, prod1, prod2}`; comparisons `runs/cmp_*.{txt,json}`; `key_table_v36_v37.{md,json}`.
- Flux: `N/flux/code/` (copies, with `F2_VTAG`; changed scripts backed up as `*_pre_v37.*`); `N/flux/run_flux_v37.sh`; `flux/out/`, `flux/texts/` (regenerated); `flux/snap/{reg_v36, v37_pass1, v37_pass2}`.
- The texts in `flux/texts` still carry the fix-8 wording and labels (gen_texts strings). The numbers are v37. Text updates belong to the next text round, together with `fix12_texts` (which predates D-3).
