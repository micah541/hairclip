# Lemma U, Phase 2a (approved by Micah, Nov 13 2026: "Start Phase 2a")
Source: `../phase1/CHECKPOINT.md` (Fable, CONTINUE WITH CONDITIONS) and `../phase1/PHASE1.md`.

| dir | task | owner |
|---|---|---|
| A_repair/ | certified repair of Lemma S(A), in new files (v3.1 untouched): top pair (j+1, j+1) on the ramp [14, 28]; cutoff × jet classes; U_j (j ≥ 2) in the quadratic coefficients; IBP {2,3,3} core; retained-dissipation export; outer cutoff [350, 700]; tightened (H2) | agent A (Opus) |
| U_ramp/ | Lemma U ramp inequality in time-integrated form; (2,2,3) in S33 form; transport identity for ∂_τΘ; Ĝ design (ramp share as data → rate criterion from y = 56?; localised Ĝ for zone III) | agent B (Opus) |
| checkpoint2/ | rerun of the checkpoint with the repaired E^A and honest X_{j+1}; Fable kill check | Sonnet + Fable |

Kill check 2 (from CHECKPOINT.md): continue only if the late-centre M1 is ≥ ×1.3 at the honest KS and the repaired E₂ᴬ is ≲ 0.1.
Resources (Micah): one process per agent, nice -n 19, < 800 MB; no hurry.

## Status
- **Agent B done** (U_ramp/U_RAMP.md; U_ramp.tex not compiled).
  - Prop. U, time-integrated, with the exported D_k kept symbolic.
  - X-term factor 4.14·c_D/κ: ×6.6 at κ = 1.25 (c_D = 2); ×28–41 with (A) as is.
  - Transport (Lemma T) is exact. The b-term does not close without a new factor φ = exp(K_b e^{−s}), K_b = 1e-3, which was not modelled before.
  - Ĝ design "data56", global box (0.15, 0.52) [cert]: Ĝ₂₂ = 338259, D(56) = 1.0356, rate margin ×29.
  - U_model2, defaults [est]: late M1 ×2.40 (KS = 1), ×2.09 (KS = 2); ×1.50 if (A) is not re-split.
  - Open: Lemma Y cubic classes; c₂₂₃ in arb; zone III box from D5a ≤ (0.15, 0.52), or (0.20, 0.60) as fallback; φ in Γ/K.
- Agent A (certified (A) repair + D_k export): running. Then rerun `U_model2.py --inputs <A json> --tag A1` → kill check 2.
- **Agent A done** (A_repair/A_REPAIR.md; `./run_chain.sh`, 11 s).
  - Repaired E^A [cert], tightened (H2): (1.96e-7, 9.26e-5, **0.0696**); v3.1 (H2): E₂ = 0.2686. Core M1 margin 4.36 (tightened) and 2.73 (v3.1 (H2)).
  - Costs at level 2: top pair +0.0098; cutoff × jet +0.0118 (the reviewer's ≲ 1.5e-3 is NOT confirmed); U_j; re-split r ≥ 1.249 a_min; outer cutoff on placeholder hypotheses H-O1/H-O2 (M̄ = 1, not certified).
  - Export (export_for_U.json): D = (1.40e-7, 7.59e-5, 71) + D_slope·T. **There is no horizon-free D**, because of persistent sources (forcing, coupling, outer placeholder rate). U_model2 ignores D_slope. Level 2 exports only the ε₄ reserve.
  - More v3.1 errors: ramp_arb.py:78,82–83 (jet-free IBP derivatives, z₅ jets); I₃ ×10 too small (assembly_arb.py:175); UQ = 1/20 used where U reaches 0.086 (harmless, now certified); outer cutoff never treated.
- Running: U_model2 --inputs A_repair/export_for_U.json --tag A1.
- U_model2 A1 run: late M1 ×2.57 / ×2.36. NOT valid: it used the from-τ₁ D and ignored D_slope (U_model2.py:98), and its SENS rows are dead (l.307, l.310).
- **Kill check 2 (Fable): GO WITH CONDITIONS** (KILLCHECK2.md; scratch in killcheck2_scratch/).
  - Honest late M1 ×2.28 at KS = 1.33 [rev]; worst joint stress ×1.81. E₂ᴬ = 0.0696 [cert].
  - Horizon: Lemma D as written is false. Fix: restate it in window form (C0 + C1·Σ, by summation by parts); costs 7–11% on M1. Owed: the verifier's J(T) ≤ J + sT check on every window.
  - H-O2 (pointwise z₂, z₃ on [350, 700]) has no source. Redo R5 in L² form (D3′).
  - **v3.1/v3.2 step 3 is conditional on (F) AND a new (O) (outer cutoff of (A)), and (A)'s proof as written is incomplete.** main.tex wording needs Micah's decision.
  - Bulk estimate 50–62 h (tasks 1–8 in KILLCHECK2.md). **STOPPED for Micah.**
