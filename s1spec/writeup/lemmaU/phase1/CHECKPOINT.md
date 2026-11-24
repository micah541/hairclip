# Phase-1 kill checkpoint, Lemma U — independent review (Nov 13, 2026)

**Verdict: CONTINUE WITH CONDITIONS.** The rule is met (M1 ≥ ×1.3 with certified Ĝ₂₂ and K), but the review found a real gap in closed v3.1 Lemma S(A) that feeds E^A into every Lemma-U number. Do a short Phase 2a (10–14 h) that repairs (A) and the ramp accounting, then a second kill check, before the bulk.

Labels: **[cert]** arb. **[est]** float model on certified inputs. **[rev]** this reviewer's scratch run (float, or arb arithmetic without a project script); redo in Phase 2.

## Margins, weight (4,4,0.3)

| quantity | value | label |
|---|---|---|
| Rate 28⁴/Ĝ₂₂, worst Ĝ₂₂ = 276918 (Gw 1.82) | ×2.22 | [cert] |
| Same, reconciled conventions: Ĝ₂₂ = 268126 (Gw 1.77) | ×2.29 | [rev, arb] |
| T1 as GONOGO defined it (full sweep + (H3)), KS = 1 / 2 | M1 ×1.58 / ×1.47 | [est] |
| Late centre, launch lemma, KS = 1 / 2 | M1 ×2.02 / ×1.93; M0 ×3.83 | [est] |
| Late centre, conservative | M1 ×1.86 / ×1.68 | [est] |
| Present centre, launch lemma | M1 ×2.41 / ×2.27; M0 ×1.96 (set by U at y_c = 126.6) | [est] |
| Present centre, "conservative" | M1 ×1.42 / ×1.34 — not a possible configuration | [est] |
| Ramp K: effective KS ≤ 1.33; flip at 12.3 (cons.) / 35 (launch) | ×9 or more | [cert] / [est] |
| (A) with core {2,3,3} by IBP, tightened (H2) | M1 core ×4.86, Ĝ₂₂ ×1.98 | [cert], but see finding 1 |
| **Late centre after the corrections below** (E^A repaired, honest ramp data) | **M1 ≈ ×1.5–1.7** | [rev] |
| Worst joint stress (E₂ᴬ ×1.65, KS = 8, Ĝ₂₂ 3.3e5) | M1 ×1.27 | [rev] |

I reproduced the checkpoint rows from a scratch copy of `t1_checkpoint.py` (self-check passes, rows identical), P1's Ĝ₂₂ (239750.3, 250917.8) and P2's Ĝ₂₂ (276561 with P1's exact Γ, against P2's 276918). Not rerun: `ramp_K_cert.py`, `weight_consts_arb.py`. The launch tex is uncompiled.

## Answers to the checkpoint questions

- **The ×1.34 does not govern.** It applies the full-weight sweep I_W = 4.14 to a centre already at y_c = 126.6, whose launch-lemma B₂ is 0.042, not 0.874 (`t1_checkpoint.py:97`). The governing rows are the late centres: ×2.02 (launch) or ×1.86 (conservative).
- **Taking the worst Ĝ₂₂ is correct in direction and about 3% conservative.** The P1–P2 gap splits into: U and sup|z₀| 0.05/0.02 → 0.11, +7.5%; dfrac 2/25 with foot adds, +4.2%; profile jets U_j (j ≥ 2), which P1 drops, +3.0%; looser Γ in P2, +0.13% [rev, arb].
- **Conventions the proof needs:**
  - dfrac 2/25 and the foot adds: yes, as long as one Ĝ must hold on the ramp too (level-2 budget slack is 0).
  - U = 0.09 (the support's value), not 0.11. With zones I/II only it is 0.05, but that needs the localised Ĝ.
  - sup|z₀| ≥ 0.152 unless B3U carries a bound on v. It moves Ĝ₂₂ by 1e-4.
  - U_j parts in the coefficient bounds: yes. P2 has them; P1 does not.

## Findings

1. **Real gap in closed v3.1 (A): the top-pair ramp commutators are dropped.** The cutoff-derivative part of the pair (j+1, j+1) cannot be reduced by the stated integration by parts, and the code sets its table to zero.
   - Where: `lemmaS/lemmaS_F2.py:26,56,81`; `rigor/ramp_arb.py:31,72`; `rigor/assembly_arb.py:111`; the claim in `lemmaS.tex:209`.
   - My port reproduces the certified k ≤ j tables within 2% and the stored ramp terms (4.0e-9, 3.4e-6, 0.027).
   - Size with the verifier's data [rev]: E₀ +4.9e-8 (J₁); E₁ +1.04e-5 (J₂); E₂ **+3.0** by the only scripted route (J₃).
   - Fix [rev]: at level 2 integrate by parts once and use Young against the ramp dissipation, which has 0.32 a_min spare. This gives E₂ +0.005 (δ = 0.2) to +0.010 (δ = 0.1), with z₂ data only.
   - Effect on v3.1 itself: E ≈ (3.6e-7, 1.91e-4, 0.255); S0v ×1.05; core M1 ×2.82 → about ×2.75. **Absorbed, but it needs a new argument, a rerun and an erratum.**
   - Effect on Lemma U: E^A ×(1.41, 1.13, 1.15–1.25) at the tightened (H2). Late M1 ×2.02 → about ×1.9 at KS = 1.
2. **P2's flag is confirmed: cutoff × jet terms are not bounded on (A)'s ramp.** The ramp commutators use jet-free coefficients (`lemmaS_F2.py:138`, `rigor/ramp_arb.py:76`), and the L⁴ treatment covers only the Θ-weighted part.
   - Size [rev]: foldable classes add 4e-5 to E₂; (2,2,2), weight 24 against the existing 14.9, adds about 8e-4; (2,2,3), weight 3.4, and (1,3,3), weight 0.9, are small after one integration by parts.
   - Total ≲ 1.5e-3 on E₂: 0.6% of 0.247, 3% of 0.048. Absorbed.
3. **Profile jets U_j (j ≥ 2) are dropped from the quadratic coefficients** (`lemmaS/rigor/weights_exact.py:81–83`, inherited by `phase1/ghat/ghat_arb.py:57–58`).
   - (A) core Ĝ₂₂: 77721 → 79271, margin ×1.98 → ×1.94.
   - P1's Ĝ₂₂ is low by 1.4% (U = 0.05) to 2.5% (U = 0.09). P2's is right.
4. **The ramp data X_{j+1} are a steady-state model** (`t1_crude_shell_model.py:7,31`; `t1_checkpoint.py:65`; `launch_lemma.tex:81`).
   - They assume (A) keeps 2·a_min of dissipation. It keeps 0.3 / 0.2 / 0 a_min at levels 0 / 1 / 2 (`lemmaS/rigor/weights_exact.out:2–4`).
   - Only time integrals of dissipation are bounded, so a one-off term sup W · 2E/(κ a_min) replaces the rate.
   - Honest factor on the X-terms: 8.3/κ in the launch variant. As (A) stands that is 28–41, beyond the flip at 35. With (A) levels 0–1 re-split to keep about 1.25 a_min (cheap: γ′₁ 0.246 → about 0.238), it is about 6.6.
   - v3.1's flux constant I₃ makes the same assumption (`rigor/assembly_arb.py:175`); it has ×10³ room.
5. **The IBP finding is confirmed.** With `--s33 ibp` the two scratch chains differ from v3.1 only in core Ĝ₂₂ (77720.9), γ′₂ (0.12355) and the quotes PHASE1.md lists; every box margin is unchanged. v3.1 (A) needs only an erratum on this point. Finding 1 is separate.

## Open items: can they eat the margin?

| item | size | verdict |
|---|---|---|
| (2,2,3) cutoff cubic | about 1e-3 of B₂ [rev]; reduces to the S33 form by one integration by parts | no |
| ∂_τΘ transport | W(σ) is co-moving with the scaling transport, so it cancels; residual 4·b̄·sup\|W′/W\|/28 ≈ 1e-7; the D part is the damping itself [rev, paper] | no |
| DT6 typed | 3% of B₂; +0.4% if it scales with Gw | no |
| Float grid sups | 0.5% or less | no |
| **E^A after the (A) repair** | ×1.15–1.65 at level 2 | −5 to −13% on M1 |
| **Honest X-terms** | KS 4–8 on 97% of the ramp term | −10 to −25% on M1 |
| **Box by radius** | needs a localised Ĝ, not derived; a global M0 box 0.11 → 0.15 costs Ĝ₂₂ ×1.26 [rev, arb] | open design |

## Recommendation

**Phase 2a (10–14 h), in order:**
1. Repair (A) in new files (v3.3 scratch; v3.1 untouched): the top-pair term, the cutoff × jet classes, U_j in the coefficients, levels 0–1 re-split to export retained dissipation, outer cutoff [350, 700]. Certified.
2. Write Lemma U's ramp inequality in time-integrated form; put (2,2,3) into S33 form; write the transport identity.
3. Decide the Ĝ design. Option to test: treat the ramp share of Ĝ as data, as (A) does. The rate criterion then starts at y = 56 and Gw falls to about 1.03 [rev, not run]. Derive the localised Ĝ for zone III.
4. Rerun the checkpoint. **Continue only if late-centre M1 ≥ ×1.3 at the honest KS and repaired E₂ᴬ ≤ about 0.1.**

**Needs Micah's decision:** finding 1 changes a closed step's statement (new files only, erratum in v3.3).

**Hours for the remainder:** 2a 10–14 h; then 45–55 h; plus 3–5 h of step-4 reruns. Total 55–70 agent-hours, against GONOGO's 40–50.

Scratch scripts behind the [rev] numbers (session scratchpad, not in the project): `topk_est.py`, `young.py`, `cutjet.py`, `recon.py`, `stress.py` in `/tmp/claude-1000/-home-claude/d6051bf5-c549-4610-9ba2-b9ab0aad2538/scratchpad/`.
