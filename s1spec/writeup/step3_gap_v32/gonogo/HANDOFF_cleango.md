# Handoff: "clean go" step (Nov 12, ~13:20 PST), stopped on the coordinator's instruction

Labels: **[cert]** = certified (arb, same machinery as v3.1); **[est]** = float estimate; **[paper]** = derived by hand, not yet in the proof.

- Nothing was written outside `gonogo/`.
  - Imports from `writeup/lemmaS` ran with `PYTHONDONTWRITEBYTECODE=1`.
  - `step4/` was only read. A step-4 agent was active in `step4/p2/out` (files written at 12:13).
- No job is running.

## 0. Status per task

| task | state | result |
|---|---|---|
| 1. Lemma S(A) rerun at (H2) = (1e-8, 1e-6, 1e-3) | **done, verified** | E₂ = 0.048044818 [cert]; core margins 5.38 / 6.58 / 4.86 [cert] |
| 2. step 4 meets the tightened (H2) | **done** | Ȳ margins ×4.97e4 / ×97.5 / ×25.7 [cert] |
| 3. general-weight L⁴ identity | identities **done** [sympy]; derivation **written** (§4/§5 of the note are placeholders); cost script **2 of 6 cases finished, then crashed** | Ĝ₂₂ = 2.398e5, rate margin ×2.56 [est] |
| 4. re-evaluate T1, update GONOGO.md | scripts **updated, not run**; GONOGO.md **not updated** (stale) | see §5 |
| 5. HANDOFF.md | **not written**; this file stands in for it | — |

## 1. `lemmaS_A_rerun/` (task 1)

**Setup, done by the previous agent:**
- The folder is a full scratch copy of `writeup/lemmaS/`.
- It differs from v3.1 in exactly two places (I diffed every `.py`):
  - `rigor/capture_F2.py`: `--handoff 1e-8 1e-6 1e-3` (was `2e-7 1e-4 0.2`);
  - `rigor/assembly_arb.py` l.161: absolute path to `step3nl/verifier/rate_table_fin11_a0.00015.json`.
- A direct F2 run (`F2_tight.log`) overwrote the scratch `lemmaS_F10p_M0.52.json`.
  - That is intended: `assembly_arb.py` reads its args from this file (l.88). Its args differ from `ref/lemmaS_F10p_M0.52.json` only in `handoff`.
- `ref/` holds the v3.1 originals.

**What I did:**
1. Checked consistency. The handoff value is identical in the scratch F10p args, in `cap_F10p.json` (both assemble calls) and in `rigcap`.
2. Stashed the previous outputs in `prev_agent_outputs/`.
3. Re-ran `rigor/assembly_arb.py`. Its log and JSON are **bit-identical** to the previous run (`rigor_assembly_arb.recheck.log`).
4. Re-ran `rigor/capture_F2.py` with `PYTHONHASHSEED=0`.
   - The log is identical (`rigor_capture_F2.recheck.log`).
   - `cap_F10p.json` differs only by the documented hash-order spread of 1.4e-9 in the *double* level-0 ramp table (RIGOR_PASS §1). The arb assembly does not use that table.
5. Re-ran the assembly on the fresh capture. E and γ′ are identical (`rigor_assembly_arb.recheck2.log`).
6. `weights_exact.py` and `ramp_arb.py` were **not** re-run, because they do not depend on (H2). `ramp_tables.json` is identical to v3.1's, and the `res` of `weights_exact.json` is identical (only last-ulp diagnostics differ).

**Results [cert]** (from `rigor/assembly_arb.json` / `rigor_assembly_arb.log`; "claims failing: []"):
- E = (1.2110055e-7, 8.1741747e-5, **0.048044818**). These are arb upper bounds, sliver included. v3.1: (3.1110055e-7, 1.8074175e-4, 0.24704507).
- γ′ = (0.24628897, 0.24617403, 0.14306229) (lower bounds), unchanged.
- Level-2 parts:
  - inflow 0.018664931
  - ramp_comm 0.027577956
  - extra_outer (foot + S33) 0.00080171889
  - coupling 2.1097e-7
  - forcing 5.78e-11
  - handoff 0.001
- **Floor as (H2)₂ → 0: 0.047045** (sum of the certified parts).
- S0v = 0.0037181278 and S1v = 0.093435501; core M0 = 0.016713693 and core M1 = 0.10708084.
- **Margins:**
  - changed: S0v 5.3791 (v3.1: 3.4809), M0core 6.5814 (5.8694), **M1core 4.8561 (2.8200)**;
  - unchanged: ramp M1 2.0505, ramp M0 7.0304, Ĝ₂₂ 2.3378, windows 2.0516 / 2.5438 / 4.2684, |f| ≤ ½ 5.3881, K2/K3 5.0746 / 4.8140.
- The (B)/F4B lines in the same log (for example E_B2(y48) = 0.0486) use the edge-free computation, i.e. the gap itself. They are **not** valid bounds.

## 2. Step 4 against the tightened (H2) (task 2)

`t7_step4_H2_tight.py`, `.out`, `.json`, **[cert]** (arb), < 1 s.

- **Method.** The script reproduces camp's item 2 (`step4/p2/t8se/camp.py` l.236–242) from the same inputs:
  - `amp_A1_w180.json`;
  - `out/t9_cone_kch.json[shell]` (shell route: KCH = c0 = 0, so amp = c0tip = 1832.2036647, 1.0137605e7, 1.2313676e9 per unit Ȳ);
  - `out/h1slab/win_summary.json` (BASE shells 8.0096e-11, 1.3352e-10, 2.9357e-10).
- **Self-check.** At the v3.1 thresholds it gives camp's exact limits, bit for bit: **2.392004545435252e-07, 9.852864439976837e-10, 3.631705598847083e-10**.
- **Tightened (H2) = (1e-8, 1e-6, 1e-3):**
  - Ȳ limits: 4.9694450281e-8, 9.7502788732e-11, 2.5667105973e-11.
  - Against the certified Ȳ = 1e-12 (ycomp) these are margins **×4.97e4, ×97.5, ×25.7**. The v3.1 thresholds gave ×2.39e5, ×985, ×363.
  - Step 4's binding item stays item3_E3 = 1.575679703611599e-11 (×15.7).
  - In energy form, at Ȳ = 1e-12: E^sh ≤ 8.013e-11, 4.706e-10, 1.559e-6, i.e. ×125, ×2125, ×641.

## 3. L⁴ for a general weight (task 3)

- **`l4_identity_check.py`/`.out`** (previous agent's I1, I2, extended by me), 2.5 s. **All identities verified for an arbitrary weight:**
  - I1 (L⁴), I2 ({2,3,3} by IBP);
  - I3, I4 (cutoff foot split, and the exact foot ⅓(Ξηη′)′|z₂|³ for both signs of z₂);
  - I5 (interpolation);
  - I6 (log-derivatives of εW(1+g y⁻⁴));
  - I7 (two-sided W: −ν < W′/W < μ, log-concave);
  - I8 (Θ^(m)/Θ through the factor D).
- **`L4_general_weight.md`**: the derivation, written. Conclusion [paper]:
  - The power-weight bounds hold verbatim with γ^m → Γ_m = sup|Θ^(m)/Θ| (tables, m ≤ 3), κ → κ̃ = sup|(ln εΘ)′| (L⁴ classes), κ² → κ̃₂ = sup((εΘ)″/εΘ)₊ (interpolation) and C_η → C_η^Ξ (ramp foot).
  - No new classes and no new data. "The identity fails" is excluded.
  - **§4 table (`RESULTS_TABLE`) and §5 (`COST_TEXT`) are placeholders.** Fill them from t8.
- **`t8_L4_cost.py`** computes the honest Γ_m, κ̃ and κ̃₂ for Θ = W·(1+g y⁻⁴) on the geometry (support y ≥ 28; ρ/λ ≥ 56 crude, ≥ 150 fine), and iterates the fixed point θ_a = Ĝ₂₂/28⁴.
  - Run 1: killed by me (hang from `sp.simplify` plus a 2401×1201 grid).
  - Run 2: killed by `timeout 1100` after 2 cases.
  - Run 3: **crashed at the first `json.dump`** (`TypeError: sympy Float`). **`t8_L4_cost.json` is a 209-byte partial file and is INVALID.**
  - Partial numbers **[est]**, from the run-2 console output (in this conversation; the run-3 `.out` has the references only):

| case | Γ_m | κ̃ | κ̃₂ | Ĝ_jj no L⁴ | Ĝ_jj with L⁴ | L⁴ factor (22) | rate margin 28⁴/Ĝ₂₂ | Gw |
|---|---|---|---|---|---|---|---|---|
| power ref (γ=1, κ=5), crude box | 1, 1, 1 | 5 | 25 | 2278, 1360, 1.906e4 | 2281, 2351, 7.772e4 | ×4.08 | — | — |
| previous T2 v1 (Γ=4, κ=5) | — | 5 | 25 | 4.052e4, 1.751e4, 6.702e4 | 4.054e4, 2.198e4, 1.786e5 | ×2.67 | (3.44) | — |
| **crude, μ=ν=4, a=0.3**, θ_a = 0.3901, g = 3.931e5 | 4, 16, 64 | **8.000** | **64.00** | 4.93e4, 2.303e4, 7.845e4 | 4.932e4, 2.856e4, **2.398e5** | ×3.06 (power, same box ×4.08); extras ×2.75 those of the power weight | **2.564** | **1.639** |
| fine, μ=ν=4, a=0.3 (box 0.21), θ_a = 0.1776 | 4, 16, 64 | 8.000 | 64.00 | 5.053e4, 2.196e4, 5.466e4 | 5.054e4, 2.441e4, 1.092e5 | ×2.00 (power ×2.45) | 5.63 | 1.216 |

- **Reading [est]:**
  - The damping factor's derivatives cost nothing for μ = 4. On the inner tail, e^{4σ}(1+g y⁻⁴) ∝ y⁴ + g, and on the outer tail θ ≲ 4e-3.
  - So the previous agent's shortcut (Γ = 4, κ = ν+4) was exact here.
  - The GONOGO.md figure "T2 margin ×3.4" is stale (it used κ = 5). The correct figure is **×2.56**.
  - The fine family's damping must use the **crude** Ĝ: its weight reaches y = 28, where only the crude box holds. So T5 uses Gw = 1.639, not 1.216.
- **Not run:** the cases (μ,ν,a) = (4,2,0.5) and (4,4,0.0), crude and fine, and the `ramp_foot_est` block (inner-ramp L⁴ foot terms).
  - My hand estimate [est] of the foot terms is ~1e-2 at level 2 (~1–2% of B₂ ≈ 0.67); the S33 sup-free term is ~1e-4.

## 4. Finding: the profile mismatch δU is negligible [paper; inputs cert]

- T5 used |U(α(τ)) − U_ρ| ≤ 0.0245 on the inner ramp; T1 omitted δU.
- That crude bound would dominate levels 0 and 1, but it is far too pessimistic.
- Prop. B′ (`step3_proof_v32.tex` l.236) gives the exact identity d/dτ[αλ^{−p}] = (β − paα)λ^{−p}. Hence
  |δU(y,τ)| ≤ 1.3056 y^p sup|β − paα| / (p(¼ − ā)) ≈ **1.06e-6 at y = 56**,
  using sup|β| ≤ 3.33811e-9 and ā = 2.19103e-7 (verifier, certified) and |α| ≤ α₁ = 1.5e-4 (assumed; the paα term is only 3.4e-11).
- So T1's omission is harmless, and T5 is over-conservative.
- The full proof must state this bound. It needs |α| ≤ α₁ and a ≤ ā on [τ_ρ, τ].

## 5. T1, T5 and GONOGO (task 4): modified, not run, outputs stale

**`t1_crude_shell_model.py`, rewritten:**
- The v3.1 and "tightened (estimated floor)" rows are unchanged: identical formulas when the new terms are zero.
- New **CLEAN-GO** row:
  - reads E^A and γ′ from `lemmaS_A_rerun/rigor/assembly_arb.json` [cert];
  - reads Ĝ₂₂ from `t8_L4_cost.json`;
  - asserts the t7 tightened check;
  - adds δU (Prop. B′, combined with Minkowski (√E + √dU)²), the L⁴ inner-ramp foot terms (est.) and the outer ramp (T6's Δ_j = 9.86e-8, 1.04e-6, 0.0161, times W ≤ 1, which is conservative for every ρ ≤ 380λ₁);
  - bisects the flip points for E₂ᴬ, Ĝ₂₂ and the ramp-constant factor KS.

**`t5_t6_fine.py`:** a CLEAN-GO row was added, with E^A and γ′ from the JSON, Gw from the t8 crude (4,4,0.3) case, δU from Prop. B′, and Minkowski.

**Both scripts crash until `t8_L4_cost.json` is valid.** The `.out`/`.json` on disk (08:49) come from the **previous agent's** versions.

**Interim numbers [est]** (previous agent's "FINAL" rows, Ĝ₂₂ = 2.40e5, Gw = 1.64, consistent with t8):
- T1, M1 box (y|f′| ≤ 0.52):
  - μ=ν=4, a=0.3: ×1.74 (KS=1), ×1.45 (KS=4);
  - μ=4, ν=2: ×1.54 / ×1.28;
  - a=0: ×1.88 / ×1.57;
  - M0 box: ×3.0–3.6.
  - The new terms (foot ~1e-2, outer ≤ 0.016, δU ~0) raise B₂ by ≤ ~4%, which should lower the M1 margin to ≈ ×1.7 (not run).
- T5 (crude δU, so these are lower bounds): K2 ×2.85, K3 ×2.20 (KS=1); ×2.79, ×2.12 (KS=4).
- T3 (Young route): adds 2.67e-5 to E₂ᴬ = 0.048.

**Interim reading (not a final verdict):** every finished input supports a go. E₂ᴬ = 0.048 is ×2.1 below the ~0.1 no-go trigger. The rate margin is ×2.56. Step-4 (H2) has ×25.7 in Ȳ. The formal verdict needs t8 → t1/t5 → GONOGO.

## 6. How to resume (from `/home/claude/s1spec/writeup/step3_gap_v32/gonogo`)

```bash
# 3a. FIX t8 FIRST: add default=float to BOTH json.dump calls (sympy Floats in the L4 results), e.g.
#     json.dump(out, open(...,'w'), indent=1, default=float)        (lines ~112 and ~130)
PYTHONDONTWRITEBYTECODE=1 python3 t8_L4_cost.py > t8_L4_cost.out 2>&1    # background; ~45 s setup + ~3 min refs + ~1-6 min per case (6 cases): ~25-35 min as is
python3 -c "import json; d=json.load(open('t8_L4_cost.json')); print(list(d['cases']))"   # needs crude|mu=4,nu=4,a=0.3 (+ the other two crude cases for t1)
# 4. then (seconds; t1's flip bisection ~1 min)
python3 t1_crude_shell_model.py > t1_crude_shell_model.out
python3 t5_t6_fine.py > t5_t6_fine.out
# 5. fill L4_general_weight.md section 4 (RESULTS_TABLE) and section 5 (COST_TEXT) from t8_L4_cost.out; rewrite GONOGO.md (corrections in sections 3-4 above); write HANDOFF.md
# optional re-verification of task 1 [cert], ~3 min total:
cd lemmaS_A_rerun && PYTHONHASHSEED=0 python3 rigor/capture_F2.py && PYTHONHASHSEED=0 python3 rigor/assembly_arb.py > rigor_assembly_arb.log
#   (expect 'claims failing as quoted ...: []'; weights_exact.py ~15 s and ramp_arb.py ~3 min need no rerun: independent of (H2))
# task 2 / identities:  python3 -I t7_step4_H2_tight.py  (<1 s);  PYTHONDONTWRITEBYTECODE=1 python3 -I l4_identity_check.py  (2.5 s)
```

## 7. Pitfalls and speed-ups

1. **The `ghat` class bookkeeping is slow and returns sympy Floats.**
   - Measured warm cost: 0.7–1.7 s per call on the crude box with Γ = 4.
   - But the four reference calls took ~3 min, and each fixed-point case 1–6 min (the fine box is slower). I did not profile further.
   - Suspects: per-box recomputation inside `L4.classify` (`Aterms` sympy expansions, `F.A_bound_set` arb sups), and sympy-Float arithmetic in `l4_extra_YC`.
   - Fixes:
     - copy `l4_extra_YC` (~25 lines) into the gonogo script, call `L4.classify(..., At=At)` with `At = L4.Aterms(M0, M1)` precomputed once per box, and cast every weight to `float`;
     - keep the memo on (Γ_eff, κ̃, κ̃₂, box) (already in t8);
     - skip fixed-point iterations when the constants don't change (here Γ_m, κ̃, κ̃₂ were constant, so 2 iterations suffice);
     - start by profiling one fine-box `ghat` call with cProfile.
   - Expected after the fixes: < 5 min for all six cases.
2. **`F.constants` costs ~23 s per call** (the arb ramp tables are not needed here). t8 uses `core_YC`, a replica of F1 l.86–95, asserted equal to `F.constants` at Γ = 4.
3. **No `sp.simplify` on W^(k)/W** (minutes). Lambdify the raw ratio instead.
4. **Grid size.** 2401×1201 used > 30% of RAM and thrashed. Use 401×1601 plus the exact σ → ±∞ limits (already in t8).
5. **Foreground limit.** Foreground tool calls cap at 10 min, and a `timeout 1100` wrapper killed run 2. Run long jobs in the background and write the JSON after each case (already in t8).
6. **Rounding.** camp's `L`/`U` round outward by one ulp (`math.nextafter`). Replicate that to reproduce camp bit for bit (t7 does).
7. **capture_F2.** Run it with `PYTHONHASHSEED=0`, alone (it is CPU-heavy). Double-ramp diffs of ~1e-9 are the known hash-order spread and are harmless.
8. **(B)/F4B lines in the Lemma S(A) log** are edge-free (the gap). Never quote them as bounds.
9. **Never edit** `writeup/lemmaS*`, `step3_proof*` or `step4/`. Import with `PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode = True` so that no `__pycache__` lands outside `gonogo/`.

## 8. Files I created or changed (all under `gonogo/`)

- **New:**
  - `t7_step4_H2_tight.{py,out,json}`;
  - `t8_L4_cost.{py,out}` and `t8_L4_cost.json` (partial, invalid);
  - `L4_general_weight.md`;
  - `HANDOFF_cleango.md`;
  - in `lemmaS_A_rerun/`: `prev_agent_outputs/`, `rigor_assembly_arb.recheck{,2}.log`, `rigor_capture_F2.recheck.log`.
- **Changed:**
  - `l4_identity_check.{py,out}` (extended);
  - `t1_crude_shell_model.py` (rewritten; outputs stale);
  - `t5_t6_fine.py` (CLEAN-GO row; outputs stale);
  - `lemmaS_A_rerun/rigor/{cap_F10p,lemmaS_rigcap_M0.52,assembly_arb}.json`, regenerated (certified content identical).
- **Untouched:** `GONOGO.md` (stale: T2 margin ×3.4 → ×2.56; T4 is no longer "rough").
