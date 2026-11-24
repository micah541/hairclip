# P3: the launch lemma and the T3 rerun (Nov 13, 2026)

Labels: [cert] = arb; [est] = float; [paper] = proved by hand in `launch_lemma.tex`.

- `launch_lemma.tex` has **not been compiled**: this machine has no LaTeX. Braces, `$` and environments balance. Compile it elsewhere before review.
- Model: `launch_model.py`, with `.out` and `.json` (float, about 50 MB RSS, under 1 s).

## Launch lemma, in short
- **Launch.** Present centres (ρ ≥ 56λ₁) launch at τ₁. Late centres launch when ρ/λ = 56.
  - The deviation is v_ρ = f − U_ρ, where U_ρ is the profile frozen at launch. All present centres share U₁.
- **Present centres' data** [paper]: D = G_w [c₂(y_c)·(H2) + c₃(y_c)·(H3)].
  - (H2) covers [28, 60]; the five (H3) windows cover [50, 379.7]; χ_out = 0 beyond 379.5λ₁.
- **Late centres' data** [paper]: D = G_w [C_A E^A + (√(q·B̂′) + √Δ)²].
  - At launch v_ρ = v_A exactly.
  - (A) is uncut on [28, 350], and C_A = sup W·y/28.
  - Beyond y = 350, the launched neighbour windows contribute q = Σ sup W / w_min.
  - Δ is the mismatch between the neighbours' frozen profiles. It is finite only if ν > 2p.
  - Neighbour factor [est]: 1/(1 − G_w q) = 1.031 for ν = 4, but 2.10 for ν = 2. So ν = 2 is excluded.
- **δU** (Lemma 3, [paper]): |U(α(τ)) − U_ρ| ≤ 1.3056 y^p sup|β − paα| / (p(¼ − ā)).
  - The constant is 1.719e-8; at y = 56 the bound is 1.06e-6.
  - Hypotheses: |α| ≤ α₁, |a| ≤ ā < ¼ and |β| ≤ β̄ on [τ_ρ, τ], and α_τ = −(p/4)α + β.
- **After launch, only the inner tail crosses the ramp**: I_W^in(56) = 1.44, against T1's I_W = 4.14.
- **The box by radius** (three zones):
  - core and zones I/II (y ≥ 28, r < 190λ₁): box (0.11, 0.52). Here |U_ρ| ≤ U(190) = 0.042.
  - zone III (r ≥ 190λ₁, where χ_out < 1): an O(1) box from Lemma O (D5a). Its rate costs a Gronwall factor exp(Ĝ·190⁻⁴/(1 − 4ā)) ≤ 1.008 for Ĝ ≤ 1e7.
  - This removes the ×1.07 of GONOGO error 5.
- **Order and open–closed argument:**
  - The order is D1 → D5a (needs only B4, (H5) and the horizon) → (A) with outer cutoff [350, 700] (Young route; consumes U only at the hypothesis level 2B̂) → U-crude (consumes the derived E^A) → U-fine → H6 (D5b).
  - The new hypothesis B3U (all launched energies ≤ 2B̂, plus the zone I/II box) replaces B3.
  - Openness uses joint continuity of (τ, ρ) ↦ 𝔈^ρ(τ) on a compact set.

## Numbers [est], weight (4,4,0.3), late centre, KS = 1 / 2 / 4
- B₂ = 0.440 / 0.481 / 0.561. Level-2 parts: launch 0.19, inflow 0.193, ramp 0.039·KS, foot 0.003, outer 0.016.
- **Launch-limited M1 margin: ×2.13 / 2.05 / 1.91.** M0: ×3.90 / 3.76 / 3.56.
- Reviewer variant (ramp and foot terms swept with the full I_W, conservative): B₂ = 0.521 / 0.636 / 0.867 and M1 ×1.98 / 1.81 / 1.57. This reproduces GONOGO's [rev] values (0.52; ×1.98 / 1.80 / 1.56).
- T1 is reproduced exactly. Adding (H3) to its data (error 3) gives M1 ×1.725 → 1.666, 1.612 → 1.563, 1.441 → 1.406.
- Present centres: the M1 margin is ≥ ×2.19 everywhere. The M0 margin is ×1.98 at the end of zone II (y_c = 126.7) and ×1.12 at y_c = 253, which is zone III and not used.
- Other weights: (4,4,0) gives M1 ×2.28–2.08; (4,2,0.5) gives ×1.42–1.26 (×1.01 in the conservative variant).

## T3 rerun (`gonogo/t3_A_outer.py`)
- Fix: the script now reads t1's dict. Backup: `t3_A_outer_pre_dictfix.py`. Previous outputs: `t3_A_outer.prev0824.{out,json}`.
- The script now writes an args block with input hashes. It runs from `writeup/lemmaS_v32` with `python3 -B`, and no bytecode is written there.
- CLEAN-GO rows: the Young route adds **3.20e-5 (KS = 2)** or **4.09e-5 (KS = 4)** to E₂ᴬ = 0.048 [est], about ×1500. Previously 2.67e-5, and that row reproduces exactly.
- With the window energies of the launch model, the addition is 1.9–2.2e-5 at level B̂ and 3.8–4.5e-5 at level 2B̂.

## To certify
See §9 of the tex for the full list:
1. Lemma U's energy inequality with general-weight Ĝ and K_jk, the inflow constant, the L⁴ foot, and I_W^in;
2. C_A, c₂, c₃, q, w_min and Δ in arb;
3. c_δU in arb;
4. window Sobolev constants for a non-flat W;
5. T3 in arb;
6. the zone III box from D5a, converted to the S₁ chart, with Ĝ(M^O);
7. the outer-ramp Δ^R with the G_w factor;
8. (A) rerun with its outer cutoff and the B3U/D5a boxes;
9. joint continuity, written out.
