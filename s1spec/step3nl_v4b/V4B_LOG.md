# V4B_LOG: verifier v4b (tag fin13v4b), fix F-B of the E0 audit — running log

Work order: E0_AUDIT.md (fix5/e0_audit) section 6 F-B, plus V4_REVIEW R-3 (tail alpha-condition along the tail), R-2 if cheap.
Base: step3nl_v4 (fin13v4), never edited. Rules: arb outward, constants raw -> rounded up, args blocks, nice 19, PYTHONHASHSEED=0, python3 -B, < 800 MB, one process.

## S0 setup (Nov 19-20)
- Baseline sha256 of every file of step3nl_v4 (589) and step3nl_v33 (560; 555 + 5 pre-existing Nov 13 __pycache__ files in v33/verifier) -> audit_v4b/*_before.sha.
- `cp -a step3nl_v4 step3nl_v4b`; removed from the copy: all __pycache__, runs/, v4cert/out_A, probe_E0/, audit/ (v4 records stay in step3nl_v4).
- New stage dir `v4bcert/` (common module v4bcommon.py, outputs v4bcert/out/).
## Certification stages (v4bcert)
- tc_weighted.py [sym]: T_c with weight rho*X (X = chi^2): X=1 reduces exactly to the closed nlin_E0.pkl; exact Euler-operator residual 0 for the curvature
  divergence form and for the drift IBP (c_dr); drift split and linearisation exact.
- e0_strip_drift.py [cert]: C_DR_CORE 0.549308 (audit est 0.549, allowance 1.1); strip forms M_STRIP (00 4.2e-15, 01 1.60, 02 23.45, 11 1.15, 12 16.94, 22 3.21) per unit alpha;
  IR 1.84e-11, IRW 2.13e-8, RW 4.70e-5; cubic drift CD1 0.607, CD0 2.5e-4; (a,b) w-parts CA1 32.0, CA0 7.1e-4, CB1 1.032, CB0 4.2e-4; AG_MIN 0.9688;
  alpha^2 remainders RA2 26.90, RB2 5.50.  AMAX 2.5e-4, |w| <= 0.02.  6 s, 144 MB.
- e0_lemma_consts.py (promotion of the audit's e0_cert.py; shared sections byte-equal to the audit's e0_cert.json): A^2 budget x1.679 at alpha <= 2.5e-4
  (x2.786 at 1.6e-4) with the certified drift; forcing budget 5.3439 <= 5.72 (x1.070) incl. all Fn pieces; n3_123 6.588185; constants -> out/e0v4b_constants.json.
- e0_identities.py: I1-I8 [sym] all exact (rank-one symmetrisation, oblique projection, modulation algebra, (a,b) w-parts, cubic drift, alpha^2 remainders, chi_32).
- Finding while writing the lemma: the audit's E0-4 list misses some higher-order modulation pieces (QAB Q and strip parts of (a,b) in Q1 and in the loop;
  beta - beta_lead in the h2 and E33 terms; the alpha^2 remainders of A1(u0), B1(u0) in Q1; the H31/H32 part of beta_E). All tiny; v4b CODES them as explicit E0 terms
  (H4-H7) instead of absorbing them; their Fn parts are in the forcing budget.
## Verifier (fin13v4b), gated by V4B (default 1); backups constants/rhs_iv/closure/mod_sups _pre_v4b.py (byte copies of fin13v4)
- constants.py: ALIN from e0v4b_constants.json (raw, rounded up); AMAX_E0 = 2.5e-4, WB_E0 = 0.02; STRIP_BETA/STRIP_AB from r2_strip_ibp.json (R-2).
- rhs_iv.py: E0 strip commutator term with Phi1 (E0-6); beta_S coupling n3_123 = 6.588186 instead of 7.52 (line 184);
  coded terms H1 (strip alpha-linear form), H2 (cubic drift), H3 ((a,b) w-parts), H4-H7 (higher-order modulation pieces).
- closure.py: box_assert adds E0_alpha (alpha <= 2.5e-4) and E0_w_all34 (sup|w| on [0,34.25] <= 0.02, worst case over shapes);
  tail alpha-condition checked on the whole tail [L_F, L*] by the monotone split (V4_REVIEW R-3); args pin e0v4b_constants.json, r2_strip_ibp.json.
- Smoke: tail_cert PASS with V4B=0 and 1; first cell U +0.030356 -> +0.030438, anchor -0.002684 -> -0.002667 (single cells).
## R-2 (v4bcert/r2_strip_ibp.py) [cert]: one IBP of the 4 chi' w''' commutator term; kappa_{j,0..3} regression vs lem_integrals.json rel < 1e-9;
  new S3 coefficient 0; STRIP_BETA 8.0396e-11 -> 7.7866e-11, STRIP_AB 4.3678e-10 -> 4.2271e-10; hypothesis sup|w'''| <= 0.02 removed.
## Regression (runs/reg_fin13v4, V4B=0, TAG fin13v4) vs step3nl_v4/runs/P: PASS (rate tables both alpha, verify logs w/o timing, final_alpha_points, tail_cert log).
  Run with the code before the R-2 patch; the R-2 change is gated by V4B, so V4B=0 is unaffected.
## Production: run_prod_v4b.sh (cert chain + production, twice: v4bcert/out_A + runs/P_A, v4bcert/out + runs/P).
## Interruption (Nov 20): coordinator asked to add the E2 audit fix (fix5/e2_audit/E2_AUDIT.md section 5) before the final runs.
- The interrupted sequence had finished one full E0-only pass: kept as runs/E0only_P_A + v4bcert/out_E0only_A (reference for the E0/E2 split; single pass).
- E2 stages promoted into v4bcert: e2_nlin_regions.py, e2_consts.py (outputs equal to the audit's JSON, all sections); e2_v4b_constants.py -> out/e2v4b_constants.json
  (R143 14.343003, DAB 0.0131136, NF12 2.06196, N3E2 7.69698, NH2R 13.10844, SRW 4.4292e-3, ND1..8/NF1..7 recertified at alpha <= 2.5e-4 and the asserted boxes);
  e2_identity_check.py: D_rho expansion and rho f''g'' - rho f D_rho g = (B)' exact [sym], B(0) = 0 for even data [sym], quadrature 4e-13 [num].
- Verifier (V4B): E2 block + (RM + NS) Dq; r143 -> R143; nd/nf from JSON; E2 forcing in L_list times FE2; unc.lk_v4 gets the FE2 constants;
  box_assert E2_alpha (<= 2.5e-4) and E2_tip_jet (sup(|w|,|w'|) <= 1e-4 on [0,8], K = KREG[4,8]). tail_cert PASS (V4B=0 and 1).
  Single cells: first +0.030430 (E0-only +0.030438, fin13v4 +0.030356), anchor -0.002679 (E0-only -0.002667, fin13v4 -0.002684).
- Relaunched run_prod_v4b.sh (two full passes).
## Final runs (Nov 20, 00:51-01:20): two complete passes (cert chain incl. E2 stages + production): all rc=0.
- t_compare_v4b: regression (rerun with the FINAL code, V4B=0) PASS; production double run equal apart from timing; cert chain equal apart from
  runtime_s in tc_weighted.json and its sha pin in e0_identities.json (pickle byte-identical).
- fin13v4b at 1.5e-4: rho 2.677666e-3 (x0.9975), peak 5.308562e-8, J x1.00082, sup|beta| 2.35488e-9, sup|(a,b)| 2.14041e-7; whole-tail alpha margin 177.08.
- step3nl_v4 and step3nl_v33 unchanged (sha256 of every file); pinned (A)/step-4 hashes OK; no caches; no proof-tree writes.
- Copied v4 scripts with step3nl_v4 paths repointed to v4b (backups *_pre_v4b.*), not run.
## DONE: V4B_REPORT.md, SHA256SUMS.
