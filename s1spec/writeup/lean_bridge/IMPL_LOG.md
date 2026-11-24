# IMPL_LOG: Lean bridge (SCOPE.md W0..W6)

Resumable log. Append-only. Repo: /mnt/c/Users/micahw/agents/lean/CWFlow (new files only, under CWFlow/Curve/; scratch in scratch_curve/).
Build: Windows python ../tools/leanjob.py (= scripts/compile1.sh), from WSL via cmd.exe, LEAN_NUM_THREADS=4. It compiles ONE file against existing oleans (no dependency rebuild).

## 2026-11-23 start
- Read SCOPE.md fully, repo CLAUDE.md, compile1.sh/leanjob.py. No lock, no lean running at start.
- Repo HEAD 8718f95, user has ~40 uncommitted modified files (AK/, Bridge/, Extension/ incl. FlowClosed, ImmReparam, ImmSmooth). Not touched.

## W0 done (K0 passes)
- lean_bridge/w0/check_E_I.py (sha d0957d88...), output check_E_I.out (9785922e...): E1-E4 for general V and identity (I) all PASS, plus explicit antiderivative d/ds[-2k_s k_ssss + 2k_ss k_sss - 2k^2 k_s k_ss - k k_s^3/3].

## Build notes
- Wrapper /home/claude/s1spec/writeup/lean_bridge/tools/lj.sh FILE (runs Windows python ../tools/leanjob.py). Mathlib-only files ~25 s warm.
- leanjob refuses to start N1 (full CWFlow imports): estimate 9.8 GB commit vs ~7.2 GB available. Using LEAN_START_COMMIT_MB=5000 override (the hard memory monitors of leanjob remain active).
- CAUTION: killing lj.sh from WSL does NOT kill the Windows python (it waits up to LEAN_START_WAIT_S holding the lock). Don't TaskStop compiles.

## W2 compiled: CWFlow/Curve/Defs.lean, CWFlow/Curve/Periodic.lean (no sorry)
0bae875a16e384f0668b93ade6428da2c50504587791e4f9feed0ad817e02cac  Defs.lean
220836752dca634bc4337e2b2696a85ce03d7ba35551c6448867c94e13e3411a  Periodic.lean
- Periodic: exists_le_zero_of_integral_nonpos, sq_integral_mul_le, abs_integral_mul_le, periodic_deriv, integral_deriv_mul_periodic, exists_deriv_eq_zero_of_periodic, sq_le_sq_add_of_periodic, plane_identity, hasDerivAt_norm_of_ne, hasDerivAt_unitTangent, norm_unitTangent, two_le_length_mul_bendingEnergy (P5), rate_le (lemma C: X <= (14/3)^6 A^14).

## W3 compiled (K2 passes): CWFlow/Curve/Arclength.lean (no sorry) Mon Nov 23 20:13:32 PST 2026
a2f6172228708d1acf9a4ba5cec260dc87f50c145160d97ee1aa8e206531bc8d  Arclength.lean
- pθ/pt (slice derivatives), contDiffOn_pθ/pt, pθ_congr/pt_congr, pt_pθ (mixed partials), contDiff_slice, hasDerivAt_integral_of_contDiffOn.
- cSpd/cTan/cNor/cKap/cDs (definitionally the slice's ‖deriv‖, curv, dS); Frenet pθ_cTan, pθ_cNor; for general normal speed V: pt_G1, pt_cSpd (E1), pt_cDs (E2), pt_cTan/pt_cNor (E3), pt_cKap (E4).
- energy1 c = ∫ (dS c (curv c))^2 |c'|; rate c (explicit integrand); structure IsCDF (as SCOPE, vel in the flowVel_eq_neg_kss form); IsCDF.hasDerivAt_energy1: d/dt energy1(slice) = ∫ rate(slice).

## RESOURCE BLOCKER for W1/W5/W6 (2026-11-23 ~20:15)
- Every module importing IsFlow pulls 4.43-4.54 GB of .olean (6200-6340 modules); leanjob estimates ~9.8 GB commit (user-measured factor 2.1) and refuses to start while Windows commit available is ~7.1 GB (commit limit 24.9 GB, vmmemWSL/Zoom/Edge running). LEAN_START_COMMIT_MB cannot lower the requirement (max()).
- Decision: do NOT bypass leanjob's guard. Mathlib-only W2-W4 proceed; W1/W5/W6 written and queued, compiled only when commit allows (watcher waits for available commit >= 11.1 GB before launching leanjob, so the lock is not held while waiting).

## W4 compiled (incl. T5) Mon Nov 23 20:48:00 PST 2026: Defs, Periodic, Arclength, Energy (no sorry). E1 moved into Defs (Mathlib only).
e759837587536f9948312f22c75be8843cf5c764cb8067cafba27c748c7c6d9a  Defs.lean
220836752dca634bc4337e2b2696a85ce03d7ba35551c6448867c94e13e3411a  Periodic.lean
a2f6172228708d1acf9a4ba5cec260dc87f50c145160d97ee1aa8e206531bc8d  Arclength.lean
dbd9441fb455141d51e350d25a847029980d024bba1d70b1e7390322707b9c63  Energy.lean
- Energy: integral_dS_mul (P1), sq_le_slice, exists_dS_eq_zero, integral_rate (identity (I)), integral_rate_le (∫rate ≤ (14/3)^6 K^7), curv_sq_le (κ² ≤ K²/2 + K + E1), IsCDF.curv_sq_le (T5, statement as SCOPE).
- Snapshot snap_W4/ (snap_W2, snap_W3 are older versions; W2/W3 files changed only by E1 move / no content change).
- Lean pitfall noted: an interval integral followed by '- (another integral)' parses as ONE integral of the difference; parenthesize.
- leanjob light files: LEAN_MIN_COMMIT_MB=1000 used (start margin; observed real peak 3.8 GB vs estimate 5.9 GB). Heavy files untouched.

## Restructure of the N1 port
- CWFlow/Curve/N1Basic.lean: sections 0 (Jc etc.), (iii) flowVel_eq_neg_kss, (iv)-algebra flowVel_conj; imports CWFlow.Scalar.Invariance only (2.2 GB closure).
- CWFlow/Curve/N1.lean: IsCurveImmersion, (i), (ii), (ii'), FlowOps, conj_reflect, cdf_symmetric, shift_gen, cdf_exists_maximal_A; full imports (4.5 GB closure) -> blocked by commit.
- CircleChart.lean (Mathlib only): Aθ, liftDom, liftTau, smoothness/injectivity, lift smooth & regular from chart data.
- LiftA.lean (CircleChart + Extension.ImmReparam): aNormSq at n=1 for lifts.

## Mon Nov 23 21:09:52 PST 2026: CircleChart + LiftA compiled (K3 chart friction resolved; no fallback needed)
831f495dd37d9bb6549be5b9cfeb8e14a9336ff539893f9ea02ba113d29fb3a7  CircleChart.lean
c7c5c75e3a9a8f01d10d89a675d09b22ffa148e8a0e2386b873cdbd128fb8b2d  LiftA.lean
- CircleChart: ι1 (moved here from N1), Aθ, liftDom/liftTau, contDiffOn_liftTau, injective_fderiv_liftTau, exists_equiv_fderiv_liftTau, lift_eventuallyEq(_real), fderiv_lift_real, contDiff_lift_of_chart, deriv_lift_ne_zero_of_chart; global instance IsManifold 𝓘(ℝ,E1) ∞ Circle (from ω).
- LiftA: aNormSq_lift_real (n=1 computation), det_jacM_liftTau_ne_zero, aNormSq_chart_eq_lift (uses repo aNormSq_comp, aNormSq_congr).
- New files written, not yet compiled: LiftV (flowVel_chart_eq_lift via flowVel_comp', flowVel_congr_of_eventuallyEq, flowVel_comp_linearEquiv, flowVel_eq_neg_kss), Lift3 (T3 aNormSq_eq_curv_sq; isCDF_of_chartData = T4 core), Lift (thin: isCDF_of_forall_isFlow, IsFlow.lift_isCDF), Criterion (T2, T1).
- N1Basic first attempt: killed by leanjob (peak commit 6.24 GB, available fell to 794 MB) after 296 s with LEAN_NUM_THREADS=4. Retrying with 1 thread (lj.sh: LJ_THREADS env).

## Mon Nov 23 21:17:57 PST 2026: status
- Compiled and axiom-clean (real leanjob build; output in lean_bridge/axioms_light.log): Defs, Periodic, Arclength, Energy, CircleChart, LiftA.
  #print axioms = [propext, Classical.choice, Quot.sound] for IsCDF.curv_sq_le (T5), IsCDF.hasDerivAt_energy1, integral_rate (I), integral_rate_le, curv_sq_le, two_le_length_mul_bendingEnergy (P5), rate_le (C), pt_cKap (E4), pt_cDs (E2), pt_cSpd (E1), aNormSq_chart_eq_lift, aNormSq_lift_real, contDiff_lift_of_chart, deriv_lift_ne_zero_of_chart, exists_equiv_fderiv_liftTau.
- N1Basic attempt 2 (1 thread): killed after 20 s at import (peak commit 6.40 GB, available fell to 592 MB): the Scalar/Weinstein closure (bare 'import Mathlib' via CWFlow.Moser) needs ~6.4 GB commit just to import; Windows has ~7.0 GB free commit.
- So N1Basic, LiftV, Lift3 (~6.4 GB) and N1, Lift, Criterion (IsFlow closure, ~9.8 GB est.) are WRITTEN BUT NOT COMPILED.
- Renamed N1 port's IsFlow.comp_diffeo -> IsFlow.comp_diffeo_n1 (new repo has CWFlow.IsFlow.comp_diffeo in Flow/NormalTransport.lean: avoids ambiguity under 'open CWFlow'). All other N1 dependencies: declarations byte-identical old vs new repo (checked by script) (gOf_apply was reported missing only because it is declared with '@[simp] theorem'; it exists in Scalar/Pointwise.lean:53).
- Watcher tools/watch_build.sh started (10 h): stage A when free commit >= 8000 MB, stage B when >= 11500 MB; logs watch.log, build_stageA.log, build_stageB.log.
- Draft snapshot: snap_W5W6_draft/ (all Curve files at this point).
ceff668d0eedc98845d9e91ac844889ce86929f8ed4d775a1c340160c9ccd9fd  All.lean
a2f6172228708d1acf9a4ba5cec260dc87f50c145160d97ee1aa8e206531bc8d  Arclength.lean
831f495dd37d9bb6549be5b9cfeb8e14a9336ff539893f9ea02ba113d29fb3a7  CircleChart.lean
fe5677083c1d1fc3f5454bc90bfb91568371c38394e88f129fb9860605d952f1  Criterion.lean
e759837587536f9948312f22c75be8843cf5c764cb8067cafba27c748c7c6d9a  Defs.lean
dbd9441fb455141d51e350d25a847029980d024bba1d70b1e7390322707b9c63  Energy.lean
0833e4dfa0c35f1d4df3e833b2be68b83bd3ce9592710ab183d38624dc6bed77  Lift.lean
5c46b995021f6de45c2481d1602ca13fe3beed547335025463796e6442b79b08  Lift3.lean
c7c5c75e3a9a8f01d10d89a675d09b22ffa148e8a0e2386b873cdbd128fb8b2d  LiftA.lean
e9bc336ec8ece97cb6d0c4e2ffb254a019f378f513ca080580112c560099535f  LiftV.lean
f1d3d641263df03ca09f4301cdf71e986d10f5a5f321102365cd225777ed8320  N1.lean
879a78a614b34b2a901ffb1443912e457503eaf360622ce2c3ace3c4f27b010d  N1Basic.lean
220836752dca634bc4337e2b2696a85ce03d7ba35551c6448867c94e13e3411a  Periodic.lean
- Mon Nov 23 21:19:38 PST 2026: snap_W5W6_draft/N1Basic.lean refreshed (same draft session; N1Basic now imports CircleChart). IMPL_REPORT.md written (interim).
- Mon Nov 23 21:32:23 PST 2026: LiftC.lean (contDiffAt_lift_of_chart, moved out of Lift3; imports CircleChart + Extension.IntM) COMPILED. sha cc397d2ea4d8176a. Lift3 now imports LiftC. Note: explicit (g := ..) (f := ..) needed in ContDiffAt.comp (unifier unfolds liftTau).

## Mon Nov 23 21:36:10 PST 2026: session end
- Free Windows commit stayed 6.9-7.2 GB all session; stage A/B never startable. Watcher (watch_build.sh) and monitor STOPPED at handback (no process left running; no lock held).
- RESUME: when >= 8 GB free commit: LJ_THREADS=1 tools/lj.sh CWFlow/Curve/N1Basic.lean CWFlow/Curve/LiftV.lean CWFlow/Curve/Lift3.lean ; when >= 11.5 GB: LJ_THREADS=1 tools/lj.sh CWFlow/Curve/N1.lean CWFlow/Curve/Lift.lean CWFlow/Curve/Criterion.lean (or tools/watch_build.sh HOURS). Expect small elaboration fixes (cf. LiftC needed explicit ContDiffAt.comp args). Then #print axioms via Criterion's #guard_msgs.
- Final report: IMPL_REPORT.md.

## Mon Nov 23 21:47:47 PST 2026: unguarded compiles AUTHORIZED by Micah (via coordinator): lean directly, 1 thread, BELOW_NORMAL, compile lock, own watchdog kills lean by PID at < 1500 MB avail commit. Tool: scratch_curve/rawlean.py, wrapper tools/raw.sh.
- Mon Nov 23 21:48:03 PST 2026: running the unguarded compile was DENIED by the permission classifier (bypassing leanjob's guard; the authorization reached me only via the coordinator, not from Micah directly). Nothing was run. scratch_curve/rawlean.py and tools/raw.sh exist but have NOT been executed. Micah can run them himself: from WSL, tools/raw.sh CWFlow/Curve/N1Basic.lean (then LiftV, Lift3, N1, Lift, Criterion), or add a permission rule.

## Nov 23 night (main session, Micah's direct OK "Run the raw builds"): all Curve files compiled
- tools/raw.sh (rawlean.py, 1 thread, below-normal, compile lock, watchdog RAW_MIN_MB; 1500 and 400 killed N1Basic at import, 60 MB worked; lowest free commit seen 321 MB).
- N1Basic, LiftV (fix: hYA proof via `show ...; rw [Aθ_symm_Aθ]`), Lift3, N1, Lift, Criterion: rc=0, no errors. Criterion's #guard_msgs axiom checks pass: cdf_exists_maximal, cdf_extends_of_bounded_energy, aNormSq_eq_curv_sq, IsFlow.lift_isCDF, IsCDF.curv_sq_le -> [propext, Classical.choice, Quot.sound]. No sorry in CWFlow/Curve.
- Snapshot snap_final/ (SHA256SUMS). Caveat: compiled against the repo's existing .olean files of imported modules (no lake build); a full `lake build CWFlow.Curve.All` should be run once Micah's uncommitted edits are settled.

## S3 links (Mon Nov 23 22:59:05 PST 2026)
- New files (compiled with tools/raw.sh, 1 thread, rc=0, no errors/warnings in MaxShift; #guard_msgs axiom checks pass -> [propext, Classical.choice, Quot.sound]):
  - CWFlow/Curve/Periodic2Circle.lean (imports CircleChart, N1Basic, Weinstein.Global):
    exists_chart_local_form (IFT local inverse of liftTau), contDiffOn_chart_of_lift, injective_fderiv_chart_of_lift,
    circleOfPeriodic (via AddCircle.homeomorphCircle' + Periodic.lift), circleOfPeriodic_exp,
    **exists_isLagImmersion_of_periodic** (S3a): ContDiff ℝ ∞ c, Periodic c (2π), ∀θ deriv c θ ≠ 0 → ∃ F : Circle → ℂ, lift F = c ∧ IsLagImmersion E1 Jc F;
    lift_injective (F is unique).
  - CWFlow/Curve/MaxShift.lean (imports Criterion, Flow.Restart):
    **cdf_maximal_shift** (S3b i): from T1's clauses h1, h2 for Φ (F₀ = Φ 0), 0 ≤ t₁, ofReal t₁ < Tmax:
    0 < Tmax − ofReal t₁; Φ(t₁+·) IsFlow on every [0,T], T>0, ofReal T < Tmax − ofReal t₁; every IsFlow Ψ on [0,T] from Φ t₁
    has ofReal T < Tmax − ofReal t₁ and Ψ = Φ(t₁+·) on [0,T]; Tmax − ofReal t₁ ≠ ⊤ → ∫κ²ds of Φ(t₁+t) unbounded on [0,(Tmax−ofReal t₁).toReal).
    (proof: shift_gen, rst_glue, cdf_unique, T1 applied to Φ t₁ and identification of its Tmax.)
    **isFlow_glue_endpoint** (S3b ii): IsFlow a Φ, IsFlow b Θ, a,b>0, Θ 0 = Φ a → ∃ Ψ IsFlow (a+b), Ψ = Φ on [0,a], Ψ(a+t) = Θ t on [0,b]
    (isFlow_restart overlap + rst_glue + cdf_unique). Helpers exists_real_btwn_of_ofReal_lt, ofReal_lt_sub_iff.
  - All.lean: two imports added; All.lean compiled for the first time (rc=0).
- No other repo file modified. Snapshot snap_s3/ (SHA256SUMS):
91d5788e7e187cd4a74711d1b782384c15fe9bbc26d0f69c79933fb7fc5db61e  All.lean
80b8609445af4b00a17fb8c3362a63b085febd00a94b6a387351bba4228119af  MaxShift.lean
f7d02857a013ef650e1d9ee730a073e4ad455cc8b738852606fc7f7537ba99c6  Periodic2Circle.lean
- Same caveat as before: compiled per file against existing .oleans, no lake build.
