# FLUX2_LOG: implementation of flux lemma v2 (fix round 8). Resumable log.

Nov 20, 2026. Approved by Micah (Nov 20). Spec: `B/fix7_flux2/FLUX2_DESIGN.md` (sha 4405dce0...), `flux2_lemma.tex`, `probe/`.
Dir `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix8_flux2_impl` (code/, out/, snap/). `B` = `.../lemmaU/bulk`.
Rules: new files only in D; proof tree read-only (imports via sys.path, python3 -B, PYTHONDONTWRITEBYTECODE=1); one process, nice 19,
PYTHONHASHSEED=0, < 800 MB (maxrss logged in out/<stage>.time). Runner: `bash code/run_stage.sh <script.py>`.

## Resume instructions
- Each stage writes out/<stage>.json with an args block (argv, env, inputs with sha256) and a separate `timing` key.
- Before-snapshots: snap/before_*.sha (step3nl_v4b, fix3, fix4_r3, fix5, fix7_flux2, step4 manifests), snap/pycache_before.txt.
- X0/X4/X5 need `B/fix6_v36/V36_REPORT.md` (the (A) chain v36 agent writes it; do not write in fix6_v36).

## Log
1. Snapshots taken (before any run). fix6_v36 was still in its R3 regression run (pid 29695) at 07:3x.
2. X1 `code/flux2_forms.py` -> out/flux2_forms.{pkl,json}: all 35 [sym] checks PASS (identity of sec. 2.1 from the v-equation; level-0
   reduction with generic j: pairs (0,0),(1,0),(2,0),(2,1), only 4WZZ'J g g''' moved, no j'', no C_eff/k^2 derivative, affine in (j,j'),
   Euler PASS, F1 regression exact; D2v coefficients; N_low structure from graph_struct.pkl and the 12-monomial N_low'; conversions;
   trace-route residue 5/2). 1.6 s, 76 MB.
3. X2 `code/flux2_coef.py` -> out/flux2_coef.json, all checks PASS (1.5 s, 91 MB). [cert] kappa^(0) = (8.46429, 4.17394, 2.26745)
   (F1 regression j=1, abar=0 reproduces (8.42439, 4.15920, 2.25952) to 1e-12); z = (1.00002*, 0.72917, 0.83482, 1.94446, 10.37038)
   (*cell enclosure; exactly <= 1); frak a_0..3 = (10.3974, 10.4181, 8.9734, 4.4332); KAP2 1.74288; x_3 2.93116; C_eff >= 13.
   Deviations from the design, recorded: (i) OR coefficient sups read from the pinned C2 output c2_t3d.json tables/0.00025/ZB/sub/[57,60]
   (same code path, boxes meeting (57,60)) instead of a re-run; each asserted <= the ZB value. (ii) W(57) for the verifier's weight
   (S7 blend on [30.25,38.25], W=(1+x^2)^-2 below) by acb.integral: W(57) <= 1.340305e-18 [cert]; the v3.1 W57_upper 1.4371630e-18 was
   for the quintic blend (reproduced bitwise by a mirror of its code at prec 100). The (A)-chain flux block still uses the quintic value
   (not consumed by v2). (iii) abar typed as hypothesis 1e-6 (>= 2.19103e-7).
4. X3 `code/chart4_shell.py` -> out/chart4_shell.json, all checks PASS (3 min 47 s, 135 MB). [sym] chi_0..chi_4 built as in n1_chart;
   ASSERTED for chi_4: no b_5, degree <= 2 in (b2, b3, b4) with b2^2 = (u''')^2 the ONLY quadratic monomial (all other second derivatives
   and d^3/db2^3 vanish exactly), chi_4 = 0 for a ray; k <= 3 structure of n1_chart re-asserted. [cert] on the shell [54.75, 57.76]
   (boxes B-a alpha <= 2.5e-4, B-d |w| <= 1/2, B-e): Abar_0..4 = (1.7e-44, 2.3e-42, 3.0e-40, 4.1e-38, 5.5e-36), Bbar_4i <= 2.2e-35,
   B422 <= 1.4e-36; sup|D^k d| = (1.3e-4, 2.4e-4, 4.7e-4, 9.4e-4, 1.9e-3); rho_meas 2.2e-16; k <= 3 <= n1_chart core values.
   Float spot check (geometric route Gamma' = T1((1-ku)+ib), synthetic curve): T^4 f - D^4 b == chi_4 to all 16 digits (diagnostic).
   Note: first attempt was stopped in a sympy simplify on the nonzero B's (replaced by the asserted zero pattern; no result changed).
5. While waiting for V36_REPORT: v36 prod1 exists (runs/prod1), prod2 running. prod1 values (not yet final): C0_2 = 116.38, C1_2 = 0.2606,
   C0_1 = 1.6637e-4, C1_1 = 7.684e-8; kappa_J pad 1e-9 absolute (for C6). X4 `code/flux2_data.py` written.
6. X0 `inputs_manifest.py`, X4 `flux2_data.py`, X5 `flux2_assemble.py` written; X6 `flux2_sanity.py` run (float, not proof): reduced
   level-0 form == defining integral (3 tests, 30 digits), bounds hold, reflected-integral simulation Z <= 2C', Z' >= 2(phi - s).
   DRY RUN of X0/X4/X5 against fix6_v36/runs/prod1 (scratch output only; F2_OUT/F2_DRYRUN): verdict PASS. Kbar0 3.4099e-21,
   K2pt 1.390e-19, C' 2.65e-26, S' 1.06e-28, s 1.039e-21, combined margin x4.80; C2 x5e7, C3 x1547, C4 x1.8e5, C5 x1512, C6 x141,
   C7 mod_sups tail rows x2.5e8. Waiting (background loop, 600 s) for V36_REPORT.md before the production run.
7. V36_REPORT.md appeared 07:56 (sec. 2a, 8 still to be filled; final-code regressions running). Production outputs per the report =
   runs/prod2 (byte-identical to prod1). Report values: C0_2 = 114.78 (N2ref) / 116.378 (4 N2ref, used here), C1_2 0.0698 / 0.2606 (used);
   kappa_J pad 1e-9; raw consumers listed (sec. 5): PHI_S/PHI_C raw in assembly_arb.f4b feed only I3 and v3.1-only rows (unused by v2).
   X5 extended: C6 uses the certified reserve kappa_J(dn) - kappa_J_cert(up); C7 table + scripted checks (mod_sups tail row reserve,
   L* + 2C' <= last march Ly_hi). X7 `gen_texts.py` (numbers rounded in the safe direction; balance check) and X8 `run_all.sh`,
   `cmp_runs.py` written; dry runs in scratch pass. No new __pycache__. Waiting for the v36 agent to finish before the production double run
   (one heavy job at a time).
8. Coordinator (08:1x): V36 final. Added texts/lemmaD_k2_sketch.tex (statement + [hand] proof sketch of the k=2 export; written proof
   listed OPEN). I3 not used by v2 (no cushion needed). Production: run1 (old texts), run2/run3 identical apart from sha of upstream timing
   fields -> timing moved to *.timing.json sidecars (f2common.dump); run4/run5: ALL outputs byte-identical (snap/cmp_run4_run5.json OK).
   Inputs: fix6_v36/runs/prod2. Verdict PASS: Kbar0 3.40985e-21, C' 2.65232e-26, S' 1.05510e-28, s >= 1.03907e-21, combined x4.8044.
9. After-snapshots: step3nl_v4b, fix3, fix4_r3, fix5, fix7_flux2, step4 manifests unchanged (sha); fix6_v36 stable across two snapshots and
   its own SHA256SUMS 38/38 OK; no new __pycache__. Report: FLUX2_IMPL_REPORT.md. DONE.
