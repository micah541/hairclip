# V4C_REPORT: verifier fin13v4c (REVIEW_v4b_v36_proofs D-3, A-1, D-2a, A-2)

Nov 20, 2026. Base: `step3nl_v4b` (fin13v4b), never edited. Everything new is in `/home/claude/s1spec/step3nl_v4c/`. Log: `V4C_LOG.md`.
Labels: [cert] arb with outward rounding on certified inputs; [read] input taken as given; [hand] hand step.

## Verdict
**fin13v4c closes at both points, α1 = 1.5e-4 and α1 = 1.25e-4.** No margin fails; no stop rule fired (J moves by −0.016%, ρ by +0.058%).
The blocking D-3 item is replaced by a certified term. **Step 3 is still not closed** (other items of the review's §5 closing list remain).

| check | result |
|---|---|
| Regression: V4C=0, TAG fin13v4b vs `step3nl_v4b/runs/P` (`runs/reg_fin13v4b`) | **PASS**: both rate tables, final_alpha_points, tail_cert log and verify logs are equal, ignoring timing and args |
| Double run, production `runs/P_A` vs `runs/P` | **PASS**, equal apart from timing; both pin the same sha256 of every cert input |
| Double run, cert chain `v4bcert/out_A` vs `out` | **PASS**: every JSON is equal; all are byte-identical except tc_weighted.json (runtime_s) and e0_identities.json (its pin of tc_weighted.json), the same pattern as v4b; the pickle is byte-identical. See the record note below |
| Cert chain vs fin13v4b's `v4bcert/out` | equal apart from path strings (e0v4b_constants `sources` path) and timing |
| tail_cert | **PASS** with V4C=1: H (2,5), L (0,1); joint (2,5). The degree table is identical to fin13v4b, so H8 is of x-degree 3 |
| Box hypotheses (asserted at every step and at the anchor) | all ≤ 1. New D-2a asserts: CF_w1_16_31 0.116 (×8.6), CF_w_16_31 0.059, CF_alpha 0.16 (α1 = 1.5e-4); 0.114 / 0.058 / 0.13 (1.25e-4) |
| Resources | one process at a time, nice 19, PYTHONHASHSEED=0, python3 -B; peak RSS 271 MB (tc_weighted), 79 MB (production); about 6 min per α point |

Record note: `v4bcert/out_A` already existed as a copy of v4b's directory, so the driver skipped the move and pass A's cert outputs were overwritten. The stale copy is now `out_A_from_v4b/`. The cert chain was then rerun as an independent pass into `out_A/` and compared against `out/` (pass B, which production P consumes). `out_from_v4b/` holds the v4b outputs as copied (reference only).

## What changed
1. **D-3, term H8** (`v4bcert/strip_cubic_v4c.py` → `out/strip_cubic_v4c.json`; `verifier/rhs_iv.py`; `verifier/constants.py`).
   - The script is promoted from `fix11_review/d_cubic_strip.py`. Every sup and every sum is kept in arb. Each value is `nextafter(float(upper), +inf)`, and the raw arb string is stored.
   - The constants read the JSON, and an assertion checks the set of nine keys. There is no typed literal. The review values are reproduced to < 1e-15 relative.
   - The term replaces `ADD_STRIP_CUBIC·cubfar` when the env V4C=1 (module flag `SC_V4C`):
     H8 = RW[(c030 sI1 + c120 S0s)V0² + (c021 sI1 + c111 S0s)V0V1 + (c012 sI1 + c102 S0s)V1²] + IRW S0s²(c201 V1 + c210 V0) + c300 IR S0s³,
     with sI1 the raw sup|w′| on [16.25, 34.25] (x-degree 1) and IR, IRW, RW from `e0v4b_constants.json`.
   - The script checks that its own strip integrals are ≤ those of e0v4b.
   - Its hypotheses are |α| ≤ 2.5e-4 and S0, S1, S2 ≤ 0.02 on the strip. All are already asserted (E0_alpha, U3_strip_S0..S2), and constants.py asserts that the JSON box matches.
   - Coefficients [cert]: c012 0.567028, c021 2.571937, c030 0.765674, c102 0.283185, c111 3.541830, c120 1.067966, c201 2.45e-17, c210 5.59e-18, c300 4.0e-37.
2. **D-2a.** `closure.box_assert` asserts the cubic_E0fine box on [16.25, 31.25]: far_w ≤ η0 = 0.01, far_w1 ≤ η1 = 0.01 and α ≤ 1e-3, all read from that JSON.
3. **A-1** (`v4bcert/e0_A2_budget_v4c.py` → `out/e0_A2_budget_v4c.json`; a new chain stage after the E2 stages). It adds the (a,b) lead mismatch, 0.179382 per αA², and the third-order loop piece, 0.002641 at α_E.
   - Need: 38.7232. Excess: 3.4308e-3. Slack: 5.6846e-3.
   - **Margin ×1.65695** at α ≤ 2.5e-4, and ×2.747 at 1.6e-4 [cert]. The review expected ×1.657.
   - The verifier constants are unchanged.
4. **Lemma text.** `e0_lemma/E0_LEMMA_v4c.tex` (not compiled) has item (iv) with A-1, the support sentence for u₃ (A-4), H8 with its constants table, the D-2a asserts, and an updated "Not covered" paragraph (A-3: the w‴ hypothesis is removed, with R-2/R-3 lines).
5. **A-2.** The whole-tail α-condition is now reported as a ratio lhs/rhs (below).

## Old (fin13v4b) vs new (fin13v4c) [cert; ratio = new/old]
α1 = 1.5e-4:
| quantity | fin13v4b | fin13v4c | ratio |
|---|---|---|---|
| steps to the tail | 10 | 10 | |
| peak √Λ | 5.308562e-8 | 5.308499e-8 | 0.999988 |
| first-cell rate U | +0.0304299 | +0.0304270 | 0.99990 |
| tail ρ | 2.677666e-3 | 2.679225e-3 | 1.00058 |
| √L_F | 9.900007e-10 | 9.897126e-10 | 0.99971 |
| J0 const / slope | 6.025378e-12 / 3.455292e-18 | 6.024419e-12 / 3.454540e-18 | 0.999841 / 0.999782 |
| J1 const / slope | 8.468087e-12 / 4.856079e-18 | 8.466739e-12 / 4.855023e-18 | same ratios |
| J2 const / slope | 1.874092e-12 / 1.074710e-18 | 1.873793e-12 / 1.074476e-18 | same |
| J4 const / slope | 6.247218e-8 / 3.582508e-14 | 6.246224e-8 / 3.581728e-14 | same |
| sup\|(a,b)\| / sup\|β\| (mod_sups, 6 digits up) | 2.14041e-7 / 2.35488e-9 | 2.14039e-7 / 2.35488e-9 | |
| whole-tail α-condition lhs/rhs | ×1554.75 (difference 177.08) | **×1554.46** | |

α1 = 1.25e-4:
- 9 steps (unchanged);
- peak 5.214276e-8 → 5.214218e-8;
- ρ 2.849387e-3 → 2.850927e-3 (×1.00054);
- J const ×0.999861, slope ×0.999810;
- sup|(a,b)| 2.10295e-7 → 2.10294e-7; sup|β| 2.22185e-9 (unchanged);
- whole-tail α ratio ×1512.43 → **×1512.16**.

How to read this: the certified H8 is smaller, in the states that matter, than the old uncertified scaling 2.6e-3 × (far-region cubic term). The numbers therefore move slightly in the favourable direction. The review's estimate of H8 was ≤ 1.3e-4 of the rate.

## Files
- Chain: `v4bcert/run_cert_v4c.sh` (the v4b chain + strip_cubic_v4c + e0_A2_budget_v4c). Outputs in `v4bcert/out` (pass B, consumed) and `out_A` (independent pass).
- Verifier: `verifier/{constants,rhs_iv,closure}.py`, with backups `*_pre_v4c.py`; `verifier/run_v4c.sh` (`OUTDIR V4CFLAG TAG`); `verifier/t_compare_v4c.py` → `runs/t_compare_v4c.json`.
- Driver: `run_prod_v4c.sh`.
- Runs: `runs/P` (production), `runs/P_A`, `runs/reg_fin13v4b`.
- Repointed copies: 14 scripts in which step3nl_v4b → step3nl_v4c, with backups `*_pre_v4c.*`. Snapshots of the protected trees before the work: `audit_v4c/`.

## Open (unchanged by this pass)
Review §5: (GS) bootstrap lemma, C-1 interface, ei_lemma_v4.tex review, DS / Sobolev constants unaudited, D-1a pins, B1–B7 records, texts, re-pin and kits, LaTeX compile, final review.
