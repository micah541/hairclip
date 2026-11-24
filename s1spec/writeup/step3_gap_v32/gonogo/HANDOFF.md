# HANDOFF: clean-go step finished (Nov 13, 2026)

Labels: [cert] = arb; [est] = float.

## Done today
- **t8 fixed and rerun.**
  - The fix: an atomic JSON write with `default=float`. The backup is `t8_L4_cost_pre_jsonfix.py`. The run-3 partial is kept as `t8_L4_cost.run3_partial_INVALID.json`.
  - All 6 cases closed. Crude (4,4,0.3): Ĝ₂₂ = 2.398e5, rate margin ×2.56, Gw = 1.639 [est]. The damping factor costs nothing for μ = 4 (Γ_m = 4^m, κ̃ = ν+4).
  - Inner-ramp L⁴ foot terms at level 2: 0.0085 [est].
- **t1 and t5 rerun** with the clean-go inputs. The previous outputs are kept as `*.prev0849.*`.
  - T1 (crude), M1 box: ×1.72 / 1.61 / 1.44 at KS = 1 / 2 / 4 for (4,4,0.3); ×1.53–1.27 for (4,2,0.5); M0 ×3.0–3.6 [est].
  - T5 (fine): K2 ×3.0–3.1, K3 ×2.3–2.4 [est].
- **`L4_general_weight.md` §4–5 filled** from t8. The backup is `L4_general_weight_pre_t8fill.md`.
- **GONOGO.md rewritten** by the Fable review. The v1 file is kept as `GONOGO_v1.md`.
  - Verdict: **GO-WITH-RISKS**, start staged, with a kill checkpoint after the first 10–12 h.

## Errors Fable found (not fixed; feed into phase 1)
1. `../feas_P.py:11,22`: γ = −μ in |coef| tables gives spurious cancellations. The Ĝ values there, and ε(28)Ĝ₂₂ = 0.0044 in `gap_v32.tex`, are wrong. Effect on T1's ramp term ≤ 2%. The honest general-weight ramp constants have never been computed. (I verified the code.)
2. `t3_A_outer.py:40–42` reads t1's JSON as a list, so T3 crashes if rerun. T3 is stale; its effect is about 3e-5.
3. T1 omits the (H3) window data: M1 ×1.72 → ×1.67.
4. The launch accounting is inconsistent between t1 and t8/plan §1. Done the plan's way: ν = 4 gives M1 ×1.56–1.98; ν = 2 gives ×1.23–1.55.
5. The single global box |f| ≤ 0.11 out to 380λ₁ has margin only ×1.07 (×1.31 with step 4's certified |α|), so the box must be radius-dependent.
6. Labelling: E^A is certified arithmetic for a statement that is conditional on the box. The foot weights on [28,56] are an estimate.

## Recommended plan (if Micah approves)
- **Parameters:** weight (4,4,0.3); inner ramp [28,56]; (H2) = (1e-8, 1e-6, 1e-3) with (H3) included; (A)'s outer cutoff by the Young route; a box that depends on radius.
- **Phase 1 (10–12 h):**
  - certify Ĝ₂₂ for the general weight;
  - compute honest general-weight ramp constants (fixing feas_P);
  - write the launch lemma;
  - fix and rerun T3.
- **Checkpoint:** continue only if the M1 margin is about ×1.3 or better with the certified inputs.
- **Then** the rest of next-steps item 2 in CLAUDE.md.

Nothing outside `gonogo/` was modified.
