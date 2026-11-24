# Lemma U bulk proof (approved by Micah, Nov 13 2026: "Start the bulk, staged")
Plan and hours: `../phase2a/KILLCHECK2.md` (50–62 h). Report to Micah at each major milestone.

| # | task | dir | model | status |
|---|---|---|---|---|
| 1 | Restate Lemma D (window form), R5 (L² form, D3′), the bootstrap order with H-O derived; decide the box | t1_statements/ | Fable | **done** (STATEMENTS.md, CONTRACTS.md, HYPOTHESIS_O.md) |
| 2 | (A): R5 redone and certified, C0/C1 interface, H-J window check, y = 28 interface audit | t2_A/ | Opus | **done** (A2_REPORT.md, out/A_interface.json) |
| 3 | Lemma U proof text (no Lemma Y) | t3_proof/ | Opus | **done** (lemmaU.tex, T3_REPORT.md; uncompiled) |
| 4 | Arb: sweeps, launch constants, c₂₂₃, K and Ĝ at the chosen box with φ, enumerator cross-check | t4_cert/ | Opus (+Sonnet reruns) | **done** (T4_REPORT.md, U_cert.json); v33 rerun pending |
| 5 | Zone III box from D5a, Lemma R ramp, cone family | t5_zone3/ | Opus | **done, final** (Z5_REPORT.md; zone3_interface.json rerun with the U-fine B̂) |
| 6 | Step-4 hand-off reruns, manifest, kits | step4/p2/v33/ (+ out_v33/) | Opus | **done** (T6_REPORT.md, FINAL_manifest_v33.json) |
| 7 | v3.3 assembly, v3.1 erratum, main / STEP4_RESULT | t7_assembly/ | Sonnet (+Opus erratum) | **done** (V33_INDEX.md, erratum/ERRATUM.md, main_v3.tex) |
| 8 | Independent review | t8_review/ | Fable | **done: DOES NOT YET CLOSE** (REVIEW_v33.md) |

**Stop triggers** (stop and ask Micah):
- certified late M1 < ×1.5 after task 4;
- D5a cannot give zone III ≤ (0.20, 0.60);
- R5 in L² form adds > 0.01 to E₂ᴬ;
- H-J fails on windows.

Side task (Micah): draft `writeup/main/main_v2.tex` with the new conditionality, (F) + (O); main.tex untouched.
Resources: one process per agent, nice 19, < 800 MB; no hurry.

## Log
- Task 1 done. Box: global (0.15, 0.52) (rate ×29 [cert]); fallback (0.20, 0.60).
  - New risks: the initial datum of (A) on [14, 700] (late M1 may fall to ~×1.8–2.0 [est]); point terms of the C² cutoffs (~2% of E₀, E₁; also missed in v3.1); a single profile U₁ for all centres; the y = 28 interface (core Ĝ₂₂ margin → ~×1.7).
  - (O) wording is in HYPOTHESIS_O.md. Note: (F)'s (F2) already postulates (A)'s conclusion, so the accurate conditionality is "(F1) and (O)". **Micah's wording decision** before the (O) text goes into main_v2.tex (placeholder kept).
- main_v2.tex draft done (main_v2_CHANGES.md); manifest figures corrected to 34 stages / 97 upstream.
- Tasks 2 and 5 started (two Opus agents). Tasks 3 and 4 next.
- Micah chose "(F1) and (O)" for main_v2 (Sonnet applying it).
- **Task 2 stopped at A4** (T2_REPORT.md). Step 4 has nothing on [379.5, 700] at τ₁.
  - With the crude zone-III datum (A) fails (S0v ×0.107); with an (H3)-level far datum it passes (M1 ×3.8, E₂ᴬ 0.0726).
  - H-J-W passes via κ_J = 0.0423.
  - **New verifier error:** flux floor F0 = 3.28e-22 missing on march cells (closure.py:182), absorbed by κ_J.
- **Task 5 stopped on design** (T5_REPORT.md). Zone III box [cert, prelim]: b₀ = 0.1003, b₁ = 0.328, b_z = 0.0143, primary (×1.5).
  - The χ ramp [190, 379.5] is not reached by the fixed-L family (r ≥ 251.44). Proposed: r_χ ≈ 252.
- Fable resolving both (t1_statements/*_r2). Tasks 2 and 5 paused until then.
- **Fable r2 decisions** (t1_statements/R2_DECISIONS.md, STATEMENTS_r2, CONTRACTS_r2, HYPOTHESIS_O_r2).
  - (1) A new hand-off threshold **(H2f)** on [60, 700] at t₁: (5e-7, 5e-6, 1e-5). Probe on step-4 data: ×1.5 / ×4.9 / ×4.0 [probe]. K0 (flux), not S0v, binds; β_A = 5/4 in B2′.
  - (2) **r_χ = 253**: ramp [253, 379.5]; Δ^R ≈ 0.0175 [est]; present M0 ×2.23.
  - Verifier error E-V1 recorded.
  - Tasks 2 and 5 resumed under r2.
  - **Needs Micah:** approve (H2f) as a new step-4 item (task 6); β_A = 5/4 vs a verifier rerun with FLUX ×4; the main_v2 (O-H) clause.
- **Micah:** (H2f) APPROVED as a task-6 item (extended window stage; existing step-4 outputs unchanged). β_A = 5/4 vs verifier rerun with FLUX ×4: **decide after task 2's numbers**.
- **Task 5 done** (Z5_REPORT.md, zone3_interface.json, provisional).
  - Box [253, 379.5] [cert]: (0.1003, 0.4197, 0.0143), primary (thinnest ×1.24 on b₁). HO1 ok.
  - Δ^R = 0.0115 [cert]. n_III = (0.0714, 0.327, 3.04).
  - Cone family t* ×2.30; H6 ≥ ×1.86 at the fine estimate, **open until task 4 certifies U-fine B̂₂ ≤ 0.106**.
  - Notes for task 3: write ramp_U's hand-derived cubic/foot steps into Prop. U; B4's far family uses centres ≥ 490; state the dissipation split (outer foot a_min/10 at level 2).
  - Notes for task 4: global Y adds (Y₂ +0.108, Y₃ +0.374 at level 2; Y₂ +0.0044 at level 1); certify c_δU.
  - Task 6 list: t9_cone_kch, h4loc, Lemma O H1 for fixed-L cylinders; the (H2f) window stage.
- **Task 2 done** [cert]. Tightened (H2) + H2f design: Ê₂ = 0.0566; S0v ×2.48; core M1 ×4.44; Ĝ₂₂ ×1.92.
  - FLUX0/K0 ×2.39 at (5/4)Ê, ×1.50 at 2Ê. R5′ adds 2.1e-3 at 4 N2_ref.
  - **A5 design choice (agent's own, to be checked in review):** ramp box = core box (0.11, 0.52); ramp M1 margin ×1.185. (O1) can go to 0.15.
  - More v3.1 errors: per-region convention at y = 28 (weights_exact.py:246,254; ramp_arb.py:29,57).
  - Interface sent to task 4.
- **Micah: verifier rerun, FLUX ×4 (v3.3)**, β_A = 2, (H2f)₀ = 1e-6. Sonnet agent: step3nl_v33/ (copy; tag fin12f4), then t2_A_v33/. Task 4 will need the v33 interface afterwards.
- **Task 4 done** [cert]: late M1 ×2.358, M0 ×5.17; present M1 ×2.78, M0 ×2.23; loop gain 0.039; U-fine B̂₂ = 0.0258; Ĝ₂₂ = 348341 with φ, rate ×28.2; c_δU = 1.7186e-8.
  - Enumerator cross-check OK. New hypothesis M1z = 0.534 (sup|z₁| ≠ M1).
  - Fallback box: late M1 ×2.64.
- **Task 5 final**: zone3 rerun with U_cert_for_zone3.json: H6 ×1.86 [cert]; all acceptance rows pass; provisional = false. Old output kept as zone3_interface.provisional.json.
- Task 3 (proof text, Opus) started.
- **v3.3 verifier rerun done** (step3nl_v33/, fin12f4): only the slopes/sJ4/F0 moved, ×4; tail_cert identical.
  - t2_A_v33 [cert]: β_A = 2, H2f₀ = 1e-6: FLUX0/K0 ×3.51, FLUX2/K2 ×6.20 (tight); S0v ×2.14; core M1 ×4.44; E₂ᴬ 0.0566.
  - Legacy page quotes (I3 slope, Φ slope) scale ×4; update them in v3.3.
- Lemma U on v33 (t4_cert/out/U_cert_v33.json) [cert]: late M1 ×2.358, M0 ×4.90; present ×2.775 / ×2.225; G 0.039.
- zone3 rerun on the v33 U-fine (U_cert_for_zone3_v33.json): H6 ×1.86; all pass (previous: zone3_interface.t2base.json).
- **v33 chain = the final inputs.** Task 3 told.
- Task 6 started (Opus): new v33 stages only, in step4/p2/v33 and out_v33, with FINAL_manifest_v33.json; FINAL_manifest.json untouched. Items: (H2f) (1e-6, 5e-6, 1e-5); the fixed-L cone cylinders (jets, h4loc, Lemma O H1, time ratio).
- **Task 3 done**: lemmaU.tex (495 lines; 99 cited keys verified). No math gaps.
  - G1–G3 bookkeeping fixes (Sonnet, running). G4 is a report row.
  - G5 (N₃² normalisation) → task 8 review. G6 (no IBP of (η²)‴/(χ²)‴ in ramp_U/Ĝ) and G7 (design choices) → task 8. G8 → task 7.
- G1/G2 fixed (Sonnet); G3 made strict by the coordinator (b-term constants = Γ₁ + 1, Γ₁ + 4). No margins changed.
  - New hashes: A_interface 5237a96a…, U_cert_v33 42cf9617…; lemmaU.tex and constants_needed.json updated (balance PASS).
- **Task 6 done** [cert]: (H2f) ×3.19 / ×4.98 / ×4.11; tightened (H2) ×125+; fixed-L cylinders: slope ×2.64, h″ ×3.04, jets ×12+, H4 ×1119, H1 ok (sup|h₀| 0.216), H7 t* ×2.30.
  - FINAL_manifest.json unchanged (523/523). New manifest: out_v33/FINAL_manifest_v33.json (52af64e6…).
  - **For review:** the v3.2 t9_lemmaO.py never checked |ξ| > XL at the arc ends; the v33 copy checks it, with the arc lengthened to 1.3·XL. Does closed v3.2's cone H1 need an erratum?
  - Lemma O copy agrees with v3.2 only to 1e-9 (python-flint version?).
- **Task 7b done** (Sonnet): t7_assembly/V33_INDEX.md, STEP4_RESULT_v33_addendum.md, main_v3.tex (+CHANGES; unconditional "subject to review"), KITS_TODO.md.
  - Coordinator fixes: zone3 SUMS regenerated; make_ucert_for_zone3.py replaces the hand copy; zone3 rerun (numbers identical, 0d6338d0…); circular-pin note; hash citations synced.
  - For review: u5_enum reads t2_A (not v33; class structure only); main_v3 (W3)/(W4) wording; arc-end question.
- **Task 7a erratum done** (erratum/ERRATUM.md, erratum.tex): 15 items (A1–A10, B1–B2, C1–C2, D1).
  - **v3.1/v3.2 conditional status needs (H2f) too:** without it, the flux fails at 2Ê (×0.7317 [cert]); with it, ×1.1245 [cert]. (H2f) is certified only in FINAL_manifest_v33.json. main_v2's (O) needs the (O-H) clause (HYPOTHESIS_O_r2.md).
  - D1: v3.2 cone H1 stands via H8 (min distance 0.136 L [cert]); wording fix at step4_proof.tex:715.
  - New: lemmaS.tex:256 uses J on every window unchecked (covered by κ_J's W1 part).
- Task 8 (Fable independent review) started.
- **Task 8 review (Fable): v3.3 does NOT yet close step 3.** No failing inequality; the margins re-derived hold.
  - **Blocking:**
    - B0: tightened (H2) is certified for the v3.2 shell variable w′ = f − αψ̂₃′ − α²ĥ₂′, but v3.3 uses f − U₁ on [14, 700]. The difference on [14, 60] is unbounded.
    - B1: no proof text for the repaired Lemma S(A).
    - B2: no v3.3 continuation / single-point argument and no assembled step-3 document.
    - B3: main_v3 quotes fin11 tail numbers (floor, |α|).
  - Must-fix M1–M9 (small). G5, G6, A5, far ≥ 490, arc 1.3·XL, the assemble/zone3 loop and u5_enum all pass.
  - **STOPPED for Micah.**

## Fix round 1 (Micah, Nov 14: "Fix all, staged")
- fix1/b0 (Opus): (H2) variable mismatch. Certify f − U₁ on [14, 60] and/or the w′ difference; rerun t2_v33 → t4 → zone3.
- fix1/b12_text (Opus): B1 proof text of the repaired Lemma S(A); B2 v3.3 continuation/single-point argument and the assembled step-3 v3.3 document.
- Then B3 + M1–M9 (Sonnet/coordinator), then a Fable re-review.
- **B0 fixed** (fix1/b0/B0_REPORT.md), route (b): (A) stays in f − U, with a certified difference bound d at every interface.
  - New step-4 stages b0_dprof and h2u_v33: (H2) for f − U₁ ×27.3 / ×563 / ×601 [cert]. Threshold unchanged.
  - v34 chain: t2_A_v34 → U_cert_v34 → zone3. Flux ×3.41 / ×6.04; S0v ×2.125; late M1 ×2.354; H6 ×1.861.
  - FINAL_manifest_v33 is now 12 stages (8202d2c6…); the closed manifest verifies.
  - Open for re-review: the chart identification (S1-chart w′ vs the ray-graph variable) and chart terms at levels ≥ 1, inherited from v3.1/v3.2.
  - Docs to update: the list in B0_REPORT.
- Fix round 1 running: M1–M4 code fixes + chain rerun (Sonnet, fix1/M14_REPORT.md); M9 step-4 arb replacements (Opus, step4/p2/v33/M9_REPORT.md); B1/B2 text (Opus, fix1/b12_text). Next: lemmaU.tex M5 + v34 numbers; B3/M6/M7/M8 doc updates; Fable re-review.
- **B1/B2 written** (fix1/b12_text/: lemmaS_A_v33.tex, continuation_v33.tex, step3_v33_main.tex, numbers_A_v33.tex; uncompiled).
  - Continuation: τ_max = ∞ with no level-3 or pointwise-z₂ quantity, provided (P1) is Chou's L² criterion. **The authors must confirm Chou's criterion covers immersed curves.**
  - (H4) must state the graph property on [48, 700].
  - New G-A1 (the (A) assembly mixes γ′/2 and γ′ integration routes; ~0.6% on Ê₀) handed to the M14 agent.
  - lemmaU.tex needs 8 changes (B12_REPORT §4).
- **M1–M4 + G-A1 done** (fix1/M14_REPORT.md). Max move −0.55% (flux).
  - Final: Ê₀ 1.20203e-6; flux ×3.395 / ×6.011 (β_A = 2); S0v ×2.122; late M1 ×2.354; present ×2.771 / ×2.225; H6 ×1.861.
  - Hashes: A_interface 34832f43…, U_cert_v34 271430d7…, zone3 a381b1f2…, ramp_U 55317ce5….
- lemmaU.tex revision (Opus, fix1/lemmaU_rev) started. M9 still running.
- **lemmaU_v34.tex done** (fix1/lemmaU_rev/LEMMAU_REV.md). 131 keys resolve; 177 numbers are outward; only the step-4 manifest pin is pending M9.
  - Fixed B2′ to sup|v_A| + sup_core|d| ≤ 0.0153 (margin ×2.116); lemmaS_A_v33.tex Lemma z1 and step3_v33_main.tex l.66 need the same change.
  - numbers_A_v33.tex is stale (list in LEMMAU_REV.md).
  - M9 interim: far H4 at all centres ≥ 490 (min ×15.58).
- **M9 done** (step4/p2/v33/M9_REPORT.md) [cert]: arb replacements, all thresholds hold (max move 7.5e-5); far H4 at all centres ≥ 490 (min ×15.59); manifest_v33 is 19 stages (8f6f7f95…); closed manifest verifies.
- b0_dprof.py checked by the coordinator: arb cells; float only as Taylor expansion centres. OK.
- Doc sync pass (Sonnet, fix1/SYNC_REPORT.md) running. Then the Fable re-review.
- **Doc sync done** (fix1/SYNC_REPORT.md). Keys PASS (132 keys, 184 numbers outward); balance PASS (step3_v33_main has cross-file refs only).
  - The coordinator synced continuation_v33.tex: hashes, G-A1 and the far-family item (now closed by h4far_v33).
  - Left as known minor items: 3 TODO markers (ramp_U position in the M14 order, M9 run time, v34 (O3) column in the erratum); historical design docs (STATEMENTS_r2, R2_DECISIONS) keep (H2f)₀ = 5e-7; ramp_U.json has a pre-M4 u1 pin (args only).
- **Fable re-review started** (t8_review/REREVIEW_v33.md).
- **Re-review (Fable): v3.3 closes step 3 with listed fixes; no blocking item** (t8_review/REREVIEW_v33.md).
  - Must-fix:
    - N1: certify the S1-chart vs ray-graph term D^k(θ₁ − Θ) on the ramp [14, 28] (est < 1%);
    - N2: env recording in final_tail_v33;
    - N3: stale citations (lemmaS_A_v33, step3_v33_main, erratum, ic_rep, make_constants);
    - N4: remaining M3/M6/M7 parts;
    - N5: cite the blow-up criterion precisely; Z2 for all τ.
  - Authors to confirm: (P1′) the DKS 2002 / Chou 2003 L² blow-up rate for immersed curves; (P3′); (P2′); Lemma O arc by arc; the single-point statement.
  - STOPPED for Micah.

## Fix round 2 (Micah, Nov 14: "Do N1–N5")
- N1 + N2 (Opus, fix2/n12): chart-term certificate + env records; chain copy v35. Then N3–N5 doc pass (Sonnet), then a short Fable check.
- **N1 + N2 done** (fix2/n12/N12_REPORT.md) [cert]. Chart lemma; new stage n1_chart; manifest_v33 is 20 stages (a0d348ec…). Max move 0.06%.
  - v35 chain: A_interface 78fe2f8e…, U_cert_v35 960809fd…, zone3 45e8d4dc…. Late M1 ×2.354; flux ×3.365 (4·N2ref).
  - Env now scratch-built and recorded. Note: stages run with -I, so PYTHONHASHSEED is ignored; byte-identity rests on the two passes.
  - The coordinator removed the reviewer's stray step4/p2/t9h/__pycache__ (created Nov 14 07:27, not pinned).
- Next: N3–N5 + the N1 doc updates (Sonnet), then a short Fable check.
- **N3–N5 done** (fix2/N35_REPORT.md): docs on v35; chart lemma merged; keys/balance PASS. K6 margin is ×1.189 (not 1.190).
  - The coordinator moved the M3/M6 asserts into fix2/n35_asserts.py (PASS) and restored the three scripts, keeping the pinned outputs reproducible.
  - Remaining TODOs: (O3) v35 column in the erratum (unused by v3.3); main_v3:315 source for the old repaired-v3.1 numbers; ramp_U rerun position; step-4 total run time.
  - Not updated (agent reports / closed): T6_REPORT, M9_REPORT, REREVIEW §6.
- Next: a short Fable check of N1–N5.
- **Fable check N1–N5: all closed** (t8_review/CHECK_N1_N5.md). Both manifests verify. Thinnest margins: ramp M1 ×1.1847, K6 ×1.1898, H6 ×1.860.
  - Remaining minor items:
    - F1: sliver_G domain does not cover σ down to 14 − 2⁻¹⁹ (effect ~2e-8 on Ê₂);
    - F2: sI2 is used on the sliver x ∈ [16.24907, 16.25] (inherited; est 1e-4 rel);
    - F3: 5 doc places still say the asserts are inside the scripts;
    - F4: cosmetics and TODOs.
  - v3.3 status: a closed certified chain, conditional on the author literature confirmations ((P1′), (P2′), (P3′)), a LaTeX compile, and F1/F2.
  - STOPPED for Micah.
- **Micah: F1/F2 by a standalone pinned check.** Opus is adding stage f12_sliver (step4/p2/v33/F12_REPORT.md); no chain rerun. Then F3/F4 docs (Sonnet).
- **F1/F2 closed** (step4/p2/v33/F12_REPORT.md) [cert]: stage f12_sliver; 326/326 claims hold; pinned numbers move ≤ 4.8e-8 relative; no quoted number is exceeded.
  - manifest_v33 is 21 stages (8bc47e83…); the closed manifest verifies.
  - Final doc pass (Sonnet, fix2/F34_REPORT.md) running.
- **Final doc pass done** (fix2/F34_REPORT.md): check_keys PASS on the 21-stage pin; balance PASS. (Earlier log lines in this file are history; current counts are 21 stages, 8bc47e83….)
- **STATUS (Nov 14): v3.3 step 3 is a closed certified chain.** Open: author literature confirmations ((P1′) L² blow-up criterion for immersed curves, (P2′), (P3′), Lemma O arc by arc, the single-point statement); a LaTeX compile; non-blocking TODOs (main_v3:315 historical numbers; erratum (O3) v35 column; kit input list); kits not yet built/run (KITS_TODO.md).
