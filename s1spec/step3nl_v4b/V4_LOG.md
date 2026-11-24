# V4_LOG: verifier v4 (tag fin13v4), running log

Work order: /home/claude/s1spec/writeup/lemmaU/bulk/fix5/ei_spec/EI_SPEC.md §10 (13 stages). Base: step3nl_v33 (fin12f4), never edited.
Rules: arb outward, constants from raw JSON rounded up, args blocks, nice 19, PYTHONHASHSEED=0, <800 MB, one process.

## Baseline hashes (taken at start, Nov 19)
- step3nl_v33: 555 files hashed -> scratch copy, will be rechecked at end (also stored as step3nl_v4/audit/v33_before.sha).
- pinned: A_interface v35 78fe2f8e..., U_cert_v35 960809fd..., zone3 45e8d4dc..., manifest_v33 8bc47e83..., closed FINAL_manifest e01cb4da...

## Stage log
### S0 setup (done, Nov 19)
- `cp -a step3nl_v33/. step3nl_v4/`; removed old run outputs from `verifier/` (rate_table_*, verify_log_*, fin12f4_*, unc/integ logs, modsups_fin12f4, tail_cert_log*, E_I_sups*, final_alpha_points*, rate_table.json, verify_log.txt, __pycache__).
- External read-only inputs pinned by sha256 in every args block: `v4cert/v4common.py: EXT` (step3nl psi3fun/h2fun/h2res + data, step3 hstail/hs_uforms/gen_hs3, s1_tight/bracket2.json), graph_struct.pkl (v4 copy == step3nl original, 646b4444...).
- New code lives in `step3nl_v4/v4cert/` (certification stages, outputs in `v4cert/out/`) and `step3nl_v4/verifier/` (rhs etc., backups *_pre_v4.py).
### Scripts written (not all run yet)
- `v4cert/s1_identities.py` (S1), `c1_cut_sups.py` (C1), `c2_t3d.py` (C2, checkpoints in out/c2_ckpt/), `c3_modx.py` (C3), `c4_hardy.py` (C4, needs C2), `u3_strip_form.py` (U-3).
### Findings so far (read)
- U-3: `linear_alpha_S0` multiplies |alpha| S0 (ibp_run.py: kap = d_alpha Lam, same object as m_j "coefficient of |alpha| A"). The strip bound is a sum of L1 coefficients x strip sups => linear form goes through. The old STRIP_BETA counted it with alpha = 1 (overcount).
- NEW: the lemma omits the Neumann feedback (x = HR + HEx) on the strip source terms (ibp_final.py:75-81 vs cPhi_total). v4 adds it (gamma_i max_k c_k/d_k).
- NEW: typed KREG[12<x<=16.25] = 3.0208e4 < raw K_rho_16.25 = 30208.0061 (coupling_constants.json); v4 uses the arb 1/rho(16.25).
- Raw sources found: FN[0], FN[1] = c4switch_S2_1.json r_half (psi, h2); C3 = c4switch_S3_1e-3.json; KW, PSI1W, C1W, K's in agent_sob/sob_constants_P2.json; AF, BF, CF, DS, SB, SA, IA, PHI34 in agent_coup/coupling_constants.json; MU1/MU2/E0/E1 in agent_int/v2/v2_T4_P2_0.25.json; FNW/FNW2/NAB in agent_int/v3/int_constants_v3.json (R=60); modulation in agent_mod/ibp/final_5e-4.json.
### S1 identities (done, Nov 19) -- `v4cert/s1_identities.py` -> `out/s1_identities.json`
- I12 (chain rule, D X = X^a + X^b g', X^b = X_U1, all 6 coefficients), I13 (c = 1 - J^-4, d110, d300 exact), I14 (D(d000) = E1c g' + E0c with E1c = G_U1 + 5km + q/4,
  E0c = H(U) - H(u0) - c0' w - c0 g + (5km+q/4)' g, E0c(0,0) = 0; generic G), I15 (N_low' table of EI_SPEC 4.5), I16 (the full v-equation T1..T12 = zeta g_tau): all residuals exactly 0 (sympy expand/cancel).
- FD check (mpmath odefun profile, random u0, w; 6 points): table vs numerical derivative of N_low built from F itself, rel. agreement <= 2.9e-10 (limited by odefun accuracy; threshold 1e-7). First run had threshold 1e-20 -> recorded FAIL; rerun with 1e-7 in batch1.
- I14 was restructured to the generic-G form (the concrete-G expansion did not finish in 5 min); C2 then bounds B0 = sup|H_U0(U) - c0'|, B1 = sup|H_U1(U) - c0 + (5km+q/4)'|.
### Verifier code (written, Nov 19): `verifier/constants.py` (NEW block reading constants_v4.json), `rhs_iv.py` (`L_list`, `_ei_v4`, NEW branch), `unc.py` (`lk_v4` from L_list), `closure.py` (box_assert, args), `mod_sups.py` (args), `run_v4.sh`, `t_float_check.py`. Backups *_pre_v4.py (== v33, cmp checked).
### Side finding (read, not fixed; outside EI_SPEC): E0 beta-loop constants
- E0 uses PI3 = 0.01607 = n3 eps1/G33 (old lemma lead 0.002137) and f572 = 5.72 = 1 + n3 * 0.6261 (pre-feedback cPhi), while the IBP lemma used everywhere else gives BETA0 = 0.0024398 (n3 BETA0 = 0.01835) and BETAF = cPhi_total = 0.6565 (1 + n3 BETAF = 5.937). If E0's beta coupling is |beta| n3 A, both E0 constants are low (PI3 by 14%, f572 by 3.8%). Needs an E0 audit; to be quantified by a probe cell.
### C1 cut_sups (done) -- `out/c1_cut_sups.json`: PSI1 1.304871, PSI2 1.695438, RWR 0.655792, KWIR 70258.13, KAP IR/Z1/OR 0.53397/1.28274/1.58832,
  k4 = (0.442904, 0.526704, 0.424151, 0.425050), k3 = (0.33218, 0.26335, 0.10604), Z2i 0.21132, T0 37.18449, T2 39.18420, TC0 0.35992, TC2 1.27089, Q2B 2.10e-4 [cert].
  sup|zeta'| = 35/96 (IR), 35/48 (OR) exact (S7' = 140 t^3 (1-t)^3, sympy-checked); the piece sups exceed them by 1e-3 (first run asserted piece sup <= exact and failed; replaced by the exact check).
### U-3 strip form (done) -- `out/u3_strip_form.json`: goes through as a linear form (no fallback needed). With Neumann feedback:
  beta: (1.616e-9, 1.370e-9, 8.132e-10, 1.978e-10) on (S0..S3) + |alpha| 4.539e-8 S0;  (a,b): (8.782e-9, 7.440e-9, 4.415e-9, 1.079e-9) + |alpha| 2.440e-7 S0 [cert].
  Constants (S_r = 0.02, |alpha| = 5e-4): STRIP_BETA 8.04e-11 (fin12f4 typed 1.06e-9: counted the alpha-linear coefficient with alpha = 1), STRIP_AB 4.37e-10 (was 5.52e-9). Feedback adds +4.9% / +7.6%.
### C3 modx (done) -- `out/c3_modx.json` (C^4 functions): NP3 0.2595063, NH2 1.795356, NAB 0.273467, NP12 1.7168e-7, C0P 9.875e-12 [cert]; C^4 vs C^3 rel. diff <= 1e-8.
### C2 t3d: timing test 0.28 s/box/level -> ~28 min for 2950 boxes x 2 levels; full run started in background (checkpoints out/c2_ckpt/), then C4, C5.
### Cert chain pass A (done) -- `v4cert/run_cert_chain.sh` -> out/chain_status.txt: all rc=0, max RSS 571 MB (S1), 214 MB (C2), others < 100 MB.
- C2 restructured for memory (first run peaked at 1.1 GB: sp.expand of D d100; now unexpanded, reuse pb_nlin's Fpert; 209 MB, 0.02 s/box).
- C2 regression: same code path as r60_t3/pb_ramp; 30 of 37 table entries bitwise equal, the other 7 ([40,60]) are lower by <= 1.2e-3 rel because v4 evaluates at 200 bits (fin12f4 tables: 53 bits). Criterion: v4 <= fin12f4 table (holds).
- C2 values at level 2.5e-4 (sup over [10.25,60]): CBAR 3.365e-3 (fin12f4 typed d001 2.27e-2), CA 6.16e-5, CB 0.165, d010 1.54e-4, d110 0.411, d300 3.05; IR: CBAR_IR 6.92e-5, d110i 0.0588, d300i 3.00, d110^b 10.0-10.2 [cert].
- C4: hA0 1.21e-11, HA0 1.94e-9, hB0 1.96e-11, HB0 3.14e-9, hEA1 7.98e-12, HEA1 1.30e-9 (level 2.5e-4); C0Hb 1.9e-42; Phi(34.25) upper sum 9.872e6 (coarse; constants.py uses min with the coupling raw 9.8547e6).
- C5: constants_v4.json 269 entries, all outward.
- Smoke (NEW=1): tail_cert PASS (radial H (0,2), s H (2,5), joint floors-off H (2,5) L (2,5)); t_float_check PASS at 400 points.
- Running: regression NEW=0 (runs/reg_fin12f4), then production NEW=1 (runs/P).
### T(i) regression (done, PASS) -- runs/reg_fin12f4 (NEW=0): rate tables (both alpha1), final_alpha_points, tail_cert log identical to fin12f4 apart from timing/new keys; spec cells reproduced exactly (first cell U 0.03034405538749633, tail anchor U -0.0026930418625762047).
### Production pass A (done) -- runs/P_A (NEW=1, 11 min, 79 MB): tail_cert PASS, closure 1.5e-4 and 1.25e-4 close (10 / 9 steps, tail, delta_t 0.2), mod_sups ok. rho 2.684366e-3 (fin12f4 2.693042e-3), peak 5.307698e-8, J x0.9989..x1.0009, sup|(a,b)| 2.14024e-7, sup|beta| 2.35741e-9. Box hypotheses all hold (worst: |w| on [34.25,60] 0.789 of 0.5-box, strip S2 0.773, S0V 0.765, IR jet 0.747).
- First production attempt crashed in the args block (cwd-relative paths); fixed (paths from __file__).
### Pass B running (cert chain into v4cert/out, production into runs/P); A kept as v4cert/out_A, runs/P_A. Compare: verifier/t_double_run.py.
### Pass B + double run (done, PASS) -- runs/t_double_run.json (8 identical, 17 equal apart from timing / timing-carrying sha pins).
### E0 side-finding probe [probe] -- probe_E0/e0probe.json: PI3 0.01835, f572 5.937 -> tail anchor U -2.85e-4 (needs <= -2e-3), first cell +3.340e-2. Reported, not fixed.
### Final checks: step3nl_v33 (555 files) and pinned v35/manifest hashes unchanged (audit/*_after.sha). Removed __pycache__ files my imports wrote into step3nl/, step3/, s1_tight/.
### DONE: V4_REPORT.md written; SHA256SUMS regenerated. Not done: W (ei_lemma_v4.tex, U-1), sanity.py/sanity_unc.py, Fable review. Next (needs Micah): E0 audit (item 1 of V4_REPORT), then (A) chain v36 etc.
