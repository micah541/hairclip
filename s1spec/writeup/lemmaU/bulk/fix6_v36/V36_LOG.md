# V36_LOG: (A) chain v36 on verifier fin13v4b + R3 + r_2 export (resumable log)

Nov 20, 2026. Approved by Micah (Nov 20). Dir `N` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36`, `B` = `.../bulk`.
Rules: new files only; one process, nice 19, PYTHONHASHSEED=0, python3 -B, < 800 MB; outputs only into N.

## Resume instructions
- State is on disk. Each entry below says what is done. Runs are snapshotted in `N/runs/<name>/`.
- Driver: `bash N/run_chain_v36.sh <name>` with env V36_VER (fin13v4b|fin12f4), V36_R3 (1|0), V36_R2 (1|0). The driver refuses an existing snapshot name.
- Before-snapshots of protected trees: `N/audit/tree_before.sha` (step3nl_v4b, step3nl_v4, step3nl_v33, bulk excluding fix6_v36: 5391 files), `N/audit/step4_before.sha` (step4/p2/out, out_v33, v33: 1812 files). Start marker `N/audit/START_MARKER`.

## Log
1. Snapshots taken; manifests at start: closed FINAL_manifest e01cb4da..., manifest_v33 8bc47e83... (both as pinned).
2. Copied v35 sources (identical to fix2/n12/t2_A_v35 rigor/base_ibp; B/t4_cert/*.py; B/t5_zone3/*.py; ramp_U.json input), r3 module from fix4_r3, step-4 `n1_chart.py` (to rerun on fin13v4b sups, output inside N only).
3. Backups `*_pre_v36.py` (14 scripts) + `run_t2_pre_v36.sh` made before any edit (byte copies of v35).
4. fix4_r3 patches (paths_*, r3_*) applied with `patch`: every patched file is byte-identical to the fix4_r3 version. Patches copied to `N/patches/`.
5. `v36cfg.py`: verifier-source switch + outward rounding of the base_ibp argv from raw JSON. In fin12f4 mode it reproduces the v35 typed argv exactly (asserted at import). fin13v4b argv: J (6.0254e-12, 8.4681e-12, 1.8741e-12), slopes (3.4553e-18, 4.8561e-18, 1.0748e-18), J4 6.2473e-8, sJ4 3.5826e-14, sqLam 5.3086e-8, rho 2.677e-3 (down), abar 2.14041e-7, beta 2.35488e-9. The rate table's own args env is checked: FLUX0 1.64e-20, FLUX2 8.0e-20, S0V 0.0153, S1V 0.6.
6. `patches/apply_v36.py` applied (all replacements unique): v36cfg wiring in arlib/t4lib/zone3/n1_chart/r3_forcing/base_ibp argv; R3 gated by V36_R3 (assembly_rep, assembly_arb, u34_sweeps); r_2 in weights_t2 (level-2 retain r_2 = 1/1000 on [28,350], same absolute rule as levels 0,1), assembly_rep.export k = 2, export_t2 `r2_export` block; zone3 abar assert relaxed to <= 2.19103e-7; n1_chart copy writes only into N/n1/out.
7. Driver `run_chain_v36.sh` written (base_ibp reference chain is rerun every time, because its lemmaS_F10p args are the verifier inputs of (A)).
8. reg_v35 attempt 1 stopped in weights_t2: the new r_2 assertion also fired with R2=0 (v35 retained is exactly 0, ball slightly below). Assertion now guarded by V36_R2 (edit after apply_v36.py; in the final diff). Failed snapshot renamed runs/failed_reg_v35_assert.
9. reg_v35 (fin12f4, R3=0, R2=0): REGRESSION PASS vs pinned v35 (cmp_runs.py; runs/reg_v35/cmp_vs_v35.json). A_interface, hjw, ic, weights, ramp, outer, u2_K, u34 and all 5 base_ibp outputs identical (args ignored); asm_*, u1_ghat, U_cert passes, zone3 differ only in path strings. Base_ibp rerun with the generated argv reproduces the v35 base_ibp files.
10. reg_r3 (fin12f4, R3=1, R2=0): REGRESSION PASS vs fix4_r3/runs/r3run2 (runs/reg_r3/cmp_vs_r3run2.json): all differences are path strings; r3_forcing.json identical.
11. Coordinator request (flux-lemma v2, fix7_flux2/FLUX2_DESIGN.md C6/C7) before production: kappa_J pad KJ_PAD = 1e-9 absolute added in hjw_rep (production only; v36cfg.KJ_PAD); hjw_rep.json records kappa_J_certified and the padded kappa_J. Raw rate-table consumers listed in V36_REPORT.
12. prod attempt 1 stopped in weights_t2: with Bc on the 1/1000 grid the core retained ball equals r_2 = 1e-3 up to rounding and the assertion retained >= r_2 failed. Fix (R2 branch only): Bc_core(level 2) = floor((2 - r_2 - used - 1e-6)*1e6)/1e6. Snapshot kept as runs/failed_prod1_r2grid.
13. prod1, prod2 (fin13v4b, R3=1, R2=1, KJ_PAD 1e-9): all 39 chain JSON outputs byte-identical (args included); logs differ only in timing lines. Key table: key_table_v35_r3_v36.md/.json (key_table_v36.py runs v35=fix4_r3/runs/base, R3=fix4_r3/runs/r3run2, v36=runs/prod1). No margin < 1.
14. Final-code regressions (hjw_rep KJ pad and weights_t2 grid edited after reg_v35/reg_r3): reg_v35_final, reg_r3_final launched; afterwards the in-place outputs are restored from runs/prod2.
15. reg_v35_final PASS (vs v35), reg_r3_final PASS (vs fix4_r3 r3run2), n1_chart fin12f4 regression identical to step-4 n1_chart.json. In-place outputs restored from runs/prod2 (overwrite; byte-checked).
16. SHA256SUMS (38 outputs, d942c0ca...), SHA256SUMS_scripts (100, 6b749a34...). Protected trees: no pre-existing file changed (added only: fix7_flux2, fix8_flux2_impl by other agents); step4 identical; manifests e01cb4da / 8bc47e83 unchanged. V36_REPORT.md written. DONE.
