# V4C_LOG: verifier fin13v4c (resumable log)

Trigger: `writeup/lemmaU/bulk/fix11_review/REVIEW_v4b_v36_proofs.md` (D-3 blocking, A-1, D-2a, A-2). Micah approved the verifier and (A)-chain pass (Nov 20).

1. 09:18 Snapshot sha256 of protected trees -> `audit_v4c/*.sha` (step3nl_v4b 652 files, v4 589, v33 560, fix6_v36 1039, fix8 243, fix10 10, fix12 29, fix11 20; manifests e01cb4da / 8bc47e83).
2. 09:18 `rsync -a step3nl_v4b/ step3nl_v4c/` without `runs/` and the pid file. 14 scripts with the absolute path step3nl_v4b repointed to step3nl_v4c, backups `*_pre_v4c.*` (so nothing in v4c can write into v4b).
3. New cert stages in `v4bcert/`: `strip_cubic_v4c.py` (D-3; promoted from fix11_review/d_cubic_strip.py, all sums in arb, values nextafter-up; matches the review to <1e-15 rel.), `e0_A2_budget_v4c.py` (A-1; margin x1.65695). Chain driver `v4bcert/run_cert_v4c.sh` = run_cert_v4b.sh + the two stages. The copied v4b cert outputs were moved to `v4bcert/out_from_v4b/` (reference only); the chain is rerun.
4. Verifier patches (backups `verifier/{constants,rhs_iv,closure}_pre_v4c.py`): constants.py reads `strip_cubic_v4c.json` (flag SC_V4C, env V4C; first named V4C, which clashed with the fin13v4 constants dict K.V4C -> renamed before any run); rhs_iv.py term H8 replaces `astrip*cubfar` when SC_V4C; closure.py box asserts CF_w_16_31, CF_w1_16_31, CF_alpha (D-2a; boxes read from cubic_E0fine JSON) and args record. Smoke test: tail_cert PASS with V4C=0 and 1, identical degree table.
5. 09:21-09:25 cert pass A (`runs/certA.out`), all rc=0.
6. 09:3x `run_prod_v4c.sh` (resumable): regression V4C=0 TAG fin13v4b -> runs/reg_fin13v4b; prod A -> runs/P_A; cert pass B; prod B -> runs/P.
7. Finding (records): `v4bcert/out_A` already existed as a copy of v4b's out_A (rsync), so run_prod_v4c.sh skipped `mv out out_A` and the pass-A cert outputs were deleted before pass B. Production P_A and P pin the same sha256 for every cert input (both used identical files). Fix: stale copy renamed `out_A_from_v4b/`; cert chain rerun as an independent pass into `out_A/` (09:5x, all 11 stages rc=0), compared with `out/` (pass B, the one production P consumed).
8. Production numbers, comparisons and final checks: see V4C_REPORT.md and fix13_v37chain/V37_REPORT.md.
