> **Superseded by `Z5_REPORT.md` (revision r2, r_χ = 253). This is the r1 stop report, kept for the record.**

# Task 5 report: STOPPED on a design question. The zone-III box is primary under every covering [cert, preliminary]

Opus, Nov 13 2026. The contract asks for `Z5_REPORT.md`; this file is that report. Only `t5_zone3/` was written. No step-4, lemmaS*, lemmaO*, step3* or step3nl file was touched, and nothing was run in step4/.

## 1. Why I stopped (contract stop rule: "a contract item needs a design change")

The contract places zone III and the χ ramp of Lemma U on **[190, 379.5]λ₁** (STATEMENTS.tex §2, §6; CONTRACTS Z1, Z4). It also fixes the cone family as **L = 50, X = 3, centres ≥ 300λ₁** (Z5, plan §4). These two cannot both hold, for three reasons:

- **The ramp is not covered.** A member of that family has conclusion region |r − r_c| ≤ 50(1 − e). With the frame tilt e = 0.0287 [cert], the family covers only **r ≥ 251.44λ₁**. Far members cover r ≥ 265.4. So **nothing in B4 covers [190, 251.4)**. D5a gives no box there, Lemma R has no unit-50 members there (Z4: "Lemma R with unit-50 members"), and n_III has no source there for τ − τ₁ ≤ ln(700/190)/c⁻.
- **Moving the centres down to 238 breaks H6.** The cylinders would then start at 88λ₁. H6 there needs 50|k| ≤ 0.18, but the box of B3U(b) gives only about (50/88)·0.52 ≈ 0.30 [est]. That fails unless U-fine is extended to [88, 150), which the plan does not include.
- **The v3.2 family has the old gap.** With the v3.2 members (L_m = r_m/4, centres 200–380), H6 needs the strip [50, 150) at hand-off level. That is exactly the gap (F). U-crude does not supply it.

**Where the inconsistency comes from.** The go/no-go design already used the ramp **[250, 379.5]**: `gonogo/t5_t6_fine.py` (T6, unit-50 members, `t6_ramp_out_250_379.5.json`). The typed Δ^R = (9.86e-8, 1.04e-6, 0.0161) is that ramp. Task 1 carried the v3.2 Lemma R ramp [190, 379.5] into STATEMENTS, and `U_ramp.tex` l.28 has it too. `U_model2.py` l.175 mixes the two: it uses the [250, 379.5] Δ with a factor (1 + g/190⁴).

**Proposed fix (option A; needs approval, Fable/Micah):** put the χ ramp on **[r_χ, 379.5]λ₁** with r_χ ≈ 252. The value must be at least 251.44; alternatively keep 250 and let the family centres start at 298.5. Then:
- zone III becomes [r_χ, 379.5]λ₁;
- B3U(b, c) applies on {y ≥ 56, r < r_χλ₁};
- Lemma N splits O at r_χ, so τ_III = ln(700/r_χ)/c⁻ ≈ 4.09;
- on the ramp, Ψ_ρ ≤ D(r_χ)φW;
- in U_ramp.tex l.60 the b-term gets 1/r_χ in place of 1/190.

What this changes for the other tasks:
- task 1 / STATEMENTS: §2 Setting (χ), B3U, Lemma N, §6.
- task 3: the same text.
- task 4: the outer term in the Prop. U assembly.
- task 2: the segment box on (60, 350] in the r ≥ r_χ part. This is the same box, so there is no numerical change.

## 2. What is done: Z1 and Z2, certified but preliminary (`z1_box.py` → `out/z1_box_prelim.json`)

**Inputs.**
- arb, 128 bits, cell covers of width ¼λ₁ (Z2: 1λ₁), monotone enclosures.
- D5a deviations at the needed ratio from `lemmaO_v32_consts.json`:
  - cone (C₀, C₁, C₂)√t = (0.00227, 0.00748, 0.0395), times (50/L)²;
  - far: (0.00092, 0.00437, 0.0381).
- Hand-off **thresholds** only:
  - (H3) windows: E₀ ≤ 2e-7, E₁, E₂ ≤ 1e-5, weight 28/y;
  - |α₁| ≤ 1.5e-4;
  - (H5): s₀ = 0.08, m₂ = 0.03 (cone) and 0.04 (far).

**Chart conversion.** f = tan ψ, where ψ is the angle to the straight ray. Then ψ(t) ≤ atan F0(r′) + d₁, where r′ is the t₁-point at the same member-frame x and |r − r′| ≤ d₀L. Next, r|f′| ≤ (1 + b₀²)^{3/2}[D1·r/r′ + (r/L)(d₂ + 1.5·m₂·d₁(2s₀ + d₁))]. Finally, |f − U₁| ≤ (1 + b₀²)d₁ + S₀ + p·U₁·d₀L/r.

**Window data.** From the (H3) thresholds by Sobolev: F0_max = 0.09277 and S₀ ≤ 0.0067 on the last window.

| covering | zone | b₀ | b₁ | b_z | decision | margins vs (0.15, 0.52, 0.16) |
|---|---|---|---|---|---|---|
| A: fixed-L, centres ≥ 300, + far (recommended) | [251.44, 379.5] | 0.10032 | 0.32807 | 0.01429 | **primary** | ×1.495, ×1.585, ×11.2 |
| B: fixed-L, centres ≥ 238, + far | [190, 379.5] | 0.10032 | 0.32812 | 0.01429 | primary | ×1.495, ×1.585, ×11.2 |
| C: v3.2 members, + far | [190, 379.5] | 0.09548 | 0.15882 | 0.01028 | primary | ×1.571, ×3.27, ×15.6 |

- **The box decision does not depend on how the design question is settled: primary, with the thinnest margin ×1.495 on b₀.** That is above the report level of ×1.1. b₀ is almost all the hand-off |f| ≤ 0.0928 at y ≈ 380 (U₁ with α = 1.5e-4) plus C₁√t.
- **Z2, beyond the windows, [379.5, 700]λ₁:** far members, with the angle chained by the far curvature threshold 0.04/113. This gives |f| ≤ 0.2276 and |r f_r| ≤ 0.5423, so `HO1_ok` = true (×2.20, ×1.84). On [190, 379.5] Z2 follows from Z1.
- **Reproducibility:** a rerun is byte-identical. `out/SHA256SUMS` is written. Peak memory 32 MB.
- **What these numbers depend on:** B4 for the members used, the horizon from D1 (ā), and the thresholds above. They do not depend on B2′ or B3U.
- **The new hand-off thresholds this computation needs (Z6):**
  - (H5) cone s₀ = 0.08 and m₂ = 0.03 in **L = 50 units** for every fixed-L cylinder (centres 300–380, whole cylinder);
  - the per-cylinder H1–H4 and H7 for the fixed-L family;
  - the far m₂ = 0.04 as now;
  - nothing beyond (H3) and (H5).

## 3. Not done (waiting for the design decision)

- **Z3, n_III_k:** the region depends on r_χ.
- **Z4, Δ^R certified:** under option A this is a rerun of `ramp_out.py` on [r_χ, 379.5] (copied), then `ramp_cert` with unit-50 members and Lemma R at L = 50 (R4 ≤ 0.3661 [cert]), plus the b-term. The float estimate for [250, 379.5] is Δ^R = (9.86e-8, 1.04e-6, 0.0161) [est]; Δ^R₂ is about half of the acceptance limit 0.03.
- **Z5:** the fixed-L family (covering, t* ×2.30 [cert, unchanged], H6 redo, U-fine K₂, K₃).
- **Z6:** the full list, including the step-4 reruns.

**Step-4 stages that will need rerunning (task 6), as far as is known now:**
- the per-cylinder hand-off checks of `t9_cone_kch` (slope, h″, N3, H1 graph) and `h4loc` (H4), for the **fixed-L geometry** (L = 50, X = 3, centres 300–380, cylinders [r_c − 150, r_c + 150]) in place of r_m/4 units;
- `t9_h6_v32` with the new covering;
- the Lemma O H1 C₀ check for the fixed-L cylinders.

The (H3) and (H5) thresholds themselves are unchanged. The final list will come with Z5.

## Files
- `z1_box.py`: preliminary Z1 and Z2 for coverings A, B and C.
- `out/z1_box_prelim.json`: carries an `args` block and `provisional: true`.
- `out/SHA256SUMS`.

`zone3_interface.json` is **deliberately not written**, so that tasks 2 and 4 keep using their provisional inputs until the design is settled.
