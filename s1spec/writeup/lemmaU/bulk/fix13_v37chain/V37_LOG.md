# V37_LOG: (A) chain v37 on verifier fin13v4c (resumable log)

Trigger: REVIEW_v4b_v36_proofs.md D-3 (new verifier tag fin13v4c, `s1spec/step3nl_v4c/`). Approved by Micah (Nov 20).
Final production snapshot (named here for the flux stage X0, `F2_V36RUN`): `runs/prod2` (double run `runs/prod1`).

1. `rsync -a fix6_v36/ fix13_v37chain/` without `runs/` and `run.pid`. 17 scripts with the absolute fix6_v36 path repointed, backups `*_pre_v37.*`; `t4_cert/t4lib.py:12` (relative ROOT + '/bulk/fix6_v36/...') fixed by hand.
2. `v36cfg.py`: new VER fin13v4c (step3nl_v4c/runs/P), now the default; KJ_PAD 1e-9 for fin13v4b and fin13v4c; env names V36_* kept so args stay comparable with v36.
3. `run_chain_v37.sh` = run_chain_v36.sh with default fin13v4c and the n1 chart step for every VER != fin12f4 (file n1_chart_<VER>.json).
4. `cmp_runs.py`: both chain roots normalised (PATH class), n1 fin13v4c pair added.
5. Driver `run_v37_all.sh`: reg_v36 (V36_VER=fin13v4b, R3=1, R2=1; must equal fix6_v36/runs/prod2 up to paths), prod1, prod2 (fin13v4c).
6. Flux copy in `flux/` (code from fix8_flux2_impl/code; changed scripts backed up `*_pre_v37.*`): F2_VTAG (default fin13v4c) selects verifier dir, (A)-chain dir and the report that must name the run (this file for v37, V36_REPORT.md for fin13v4b).
7. 10:08-10:28 run_v37_all.sh: reg_v36 PASS (path-only vs fix6_v36/runs/prod2), prod1 = prod2 (39 JSON byte-identical). Key table key_table_v36_v37.md: no margin decreased; c_dU 1.4634274e-8; C0_2 116.359 (4 N2ref). ramp_U Y_next_primary equals v37 u1_ghat primary Y_next.
8. 10:28-10:39 flux/run_flux_v37.sh: regression (F2_VTAG=fin13v4b) equals fix8_flux2_impl/out apart from args/timing; v37_pass1 = v37_pass2 byte for byte; verdict PASS, combined x4.80453, C3 pad x1547.6 / x2236.4. (The driver's exit code 1 came only from a trailing `tail` on several files after ALLDONE.)
9. Final checks: protected trees and manifests unchanged (sha256), no proof-tree writes, no __pycache__. Report: V37_REPORT.md. Checksums: SHA256SUMS_v37, step3nl_v4c/SHA256SUMS_v4c.
