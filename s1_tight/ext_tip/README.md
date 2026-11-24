# ext_tip: tight S1 profile at the tip end (σ ∈ [σ_tip, 8])

This is the follow-up request from the step-4 hand-off (√Λ is limited by the tip-end θ widths of prof_tight).

**Result.** `prof_tip.txt` is a drop-in replacement for `prof/prof_tight.txt`. It has the same 51,246 rows, the same grid and the same columns, and every box is nested inside the prof_tight box.
- On the tip segment x = τ_tip − τ ∈ [0, 2], the θ width drops from **1.76e-9** (radius 8.8e-10) to **3.2e-12**. This is pointwise, including the second-order Taylor form.
- Requested: ≲ 2e-10. Achieved: about 60× better.
- Every component tightens by 2–3 orders of magnitude (table below).
- The tip time is now τ_tip ∈ [10.249077651120778, 10.249077651121361], width 5.8e-13 (it was 9.4e-10).
- Cost: about 15 s of CPU.

## What drove the 8.8e-10 (diagnosis: `diag_float.py`, linearised, floating point)
θ half-width at the last prof_tight row, split by direction of the σ = 8 seed set `../in_br2.txt`:

| source | θ half-width |
|---|---|
| ε* bracket (half-width 3.17e-14) | 6.13e-10 (70%) |
| d = Z_nl,1(8) box (\|d\| ≤ 4.37e-13), Re and Im | 0.90e-10 + 1.79e-10 (30%) |
| q remainder | 3e-15 |
| integrator (CAPD / linear) | ×1.001, i.e. none |

- **The integrator.** CAPD adds essentially nothing for θ, q, k and m. For p it overestimates by ×5.5, but p is not the bottleneck.
- **τ_tip.** It does not enter the fixed-τ rows of prof_tight. It does enter any use through x = τ_tip − τ: 4.7e-10 × |θ′| ≈ 7e-10.
- **The ε* bracket.** It cannot be made much tighter on its own:
  - M = m at q = 0 has width 2.3e-10 over the d-box, partly because CAPD's Poincaré map loses the q–m correlation;
  - |∂M/∂ε| = 4850;
  - so ε* is pinned only to the d-box: even with that correlation used, ε* is still 2.8e-14 wide.
- **The real issue: correlation.** The box enclosure treats ε and d as independent. On the set {M = 0} that actually contains S1, ε and d are tied together, and the spread at the tip is only ~1e-12. So neither a finer step nor a tighter ε* would have helped much. Conditioning on M = 0 does.

## Method (rigorous)
**The constraint.** The S1 state at σ = 8 is a point y* of the certified set Y = x_m + C·r0 + rr (`../in_br2.txt`). It satisfies M(y*) = 0, where M is m at the first return to q = 0 after τ₁ = 8.5. This is the function whose zero defines ε* (`s1proof.cpp`, `bracket2.py`).

**The formula.** Write y* = x_m + Cδ + e, with δ ∈ r0 and e ∈ rr. Apply the mean value theorem on the segment [x_m, y*], which lies in the convex set Y. With A = DM·C and G = DΦ_τ·C:

    Φ_τ(y*) ∈ Φ_τ(x_m) − (G₀/A₀)·(M(x_m) + DM·e) + Σ_{j≥1} (G_j − G₀A_j/A₀)·δ_j + DΦ_τ·e

**`s1cond.cpp`.** This program evaluates the formula with CAPD, using the same CAPD 1.8.1 build as prof_tight:
- **DM over all of Y** comes from a C1 Poincaré map with `computeDP`. A = DM·C = (−4850.43, 51.7765, −103.045, −0.93749, 0), which agrees with the float diagnosis.
- **DΦ_τ over all of Y** comes from a C1 solution curve (`derivative(tc)`), at every row centre.
- **Φ_τ(x_m) and M(x_m)** come from point runs.
- **Intersection.** The result is intersected with the plain box run, which is s1prof's computation and is reproduced exactly.

**Scope.** The formula holds for **every** zero of M in Y. So no uniqueness of ε* is used, which matters because step 1 does not prove uniqueness.

The same formula is applied to the tip map P and to the return time T, using DP and dT from `computeDP`. This gives the conditioned tip state and the conditioned τ_tip.

**Hulls.** `ext_tip_post.py` (arb, outward rounding) tightens the hulls. The new hull is

    (c + f′(c)·h + f″(old hull)·[0, H²]/2) ∩ old hull

At the tip this brings the p hull from 9e-5 to 1.1e-7 and the k hull from 1.3e-5 to 1.6e-8, i.e. down to the true variation.

## Format (as prof_tight)
Each row is `ta tb tc` (τ = 8 − σ, grid τ_k = Tend·k/N with Tend = 10.249077651519443 and N = 51246, exactly as prof_tight), followed by 10 numbers for the centre at τ = tc and 10 numbers for the hull over [ta, tb]. The 10 numbers are lo/hi for each of p, q, k, m, θ.
- The normalisation is that of prof_tight: θ(σ = 8) = 0, the step-3 coordinate is x = τ_tip − τ, and the ODE is k′ = m, m′ = p/4, p′ = −kq, q′ = 1 + kp, θ′ = k (′ = d/dσ).
- All numbers are doubles printed with 17 significant digits, so they read back exactly.
- Pointwise use is the same as for prof_tight: f ∈ f(c) + f′(c)·h + f″(hull)·h²/2.

Widths of this pointwise form (maximum over rows), with prof_tight's centre widths for comparison:

| x = τ_tip − τ | p | q | k | m | θ |
|---|---|---|---|---|---|
| [0, 2) | 8.5e-12 (1.9e-9) | 9.9e-12 (2.7e-9) | 9.1e-13 (7.0e-10) | 1.3e-12 (3.0e-10) | **3.2e-12 (1.8e-9)** |
| [2, 5) | 2.3e-12 (5.3e-10) | 2.4e-12 (6.3e-10) | 4.4e-13 (4.4e-10) | 4.7e-13 (2.6e-10) | 2.2e-12 (5.7e-10) |
| [5, 10.25] | 1.5e-11 (1.3e-10) | 1.1e-13 (9.1e-12) | 1.9e-12 (3.8e-11) | 3.0e-12 (3.7e-11) | 2.5e-12 (3.4e-11) |

- **Centre widths alone** have these maxima over all rows: p 1.5e-11 (near σ = 8, where the d-box dominates and nothing changes), q 1.6e-12, k 1.9e-12, m 3.0e-12, θ 2.5e-12.
- **Near the tip, the q column is dominated by the h²-term.** The q hull is 5.7e-4 wide, which is its true variation over a row, against a centre width of 1.6e-12.

**Tip data** (`ext_tip_check.json`, conditioned; all inside the s1_tight intervals):
- τ_tip ∈ [10.249077651120778, 10.249077651121361], so σ_tip ∈ [−2.2490776511213615, −2.249077651120778].
- D = −p(tip) ∈ [2.500806276619695, 2.500806276620049].
- k(0) ∈ [1.539820637871836, 1.5398206378721015].
- θ_tip ∈ [−3.432419942099867, −3.4324199420976464].
- q(tip) and m(tip) contain 0, with widths ≤ 5e-13.
- The ε of all zeros of M in Y lies in [1.5346198894026843e-4, 1.5346198896815588e-4], width 2.8e-14.

**For the consumer:** use this τ_tip interval in x = τ_tip − τ. With the old interval, the τ_tip width alone costs ~7e-10 in θ(x).

## Checks (`ext_tip_post.py`, results in `ext_tip_check.json` and `.log`; all pass)
1. **Against prof_tight.** The grid is identical to prof_tight. Each centre lies inside prof_tight's centre, each new hull inside the old hull, and each centre inside its new hull. The box-run hulls equal prof_tight's byte for byte.
2. **Widths.** As in the table above.
3. **Tip consistency.**
   - The last row, Taylor-transported to the conditioned τ_tip, meets the conditioned tip state and contains q = m = 0.
   - τ_tip, D and k(0) lie inside the s1_tight intervals.
4. **Float sanity check (not part of the certificate).** A DOP853 trajectory from a zero of M in Y stays within 0.12 radii of the enclosure centres.

## Rerun
```
cd /home/claude/s1_tight/ext_tip && bash run_tip.sh      # about 15 s; one process, nice 10
```
If `capd-config` is not on PATH, set `CAPDCFG`; the default is `/home/claude/ext/CAPD/build/bin/capd-config`.

The seed is the stored `../in_br2.txt`. For its provenance check, see `../ext813/README.md`.
