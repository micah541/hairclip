# T3 report: production step (embedded Γ0, initial phase, error budget, resizing, approach, quadratures)

## P1 production run (runner)
```
mkdir -p /home/claude/s1spec/step4/t3/p1 && cd /home/claude/s1spec/step4/t3
nohup python3 -u -I run_p1.py --target-k 12 >> p1/run.log 2>&1 &   # start; the SAME command resumes after a kill or crash
cat p1/status.json                                                  # slab, t, k_max, M, block, rho, |e|_0..4, C2/C5 monitors, CPU
pkill -f run_p1.py                                                  # safe at any time; rerun the start command to resume
nohup python3 -u -I run_p1.py --target-k 20 >> p1/run.log 2>&1 &   # later: extend the same chain to k = 20
```
- **Process.** One process, nice 10, single-threaded BLAS.
- **Checkpoints.** p1/ckpt/slab_NNNNN.{json,npz}, one per slab. Each holds the exact 128-bit centre, G, w, box, tail T, M, block, h, and the t/L/X0 enclosures. Writes are atomic, and an unreadable last file is skipped on resume.
- **Retries.** A missed absorption cap is redone inside the slab. A failed slab is retried with h/2 and a ×10 cap.
- **Stop.** The run stops when the rigorous k_max enclosure reaches the target. The last checkpoint is the P2 switch state.

| phase (expected, 1 core) | range | slabs | CPU per slab | CPU |
|---|---|---|---|---|
| I graded (Δt = t/2) | t 0.001 → 0.03 | 9 | 8–10 min | 1.0 h |
| II (h = 2/θ_max⁴) | t 0.03 → 4 | ≈ 60 | 4 → 0.5 min | 1.5 h |
| III swing | θ_max 16 → 41 | 15–25 | 1–4 min | 0.5 h |
| IV approach | k 5 → 12 | ≈ 10 | 3 → 9 min | 0.7 h |
| extension | k 12 → 20 | ≈ 5 | 10 → 26 min | +1.5 h |

Totals: **≈ 3.6 h to k = 12** and **≈ 5.1 h to k = 20**, i.e. ≤ 5.5 / 8 h with ×1.5 for redos. The test of the first 3 slabs took 0.50 h CPU including one redo (§2).

Development commands:
```
(cd ../t2 && G0FINAL=gamma0_final_unslaved.json python3 -I gamma0.py final)   # T2 curve, tight b_2 (input of 'slave'; 5 min)
python3 -I run_proto.py slave             # option (c), two passes -> ../t2/gamma0s_src.npz (2 min)
(cd ../t2 && G0SNAP=gamma0s_src.npz python3 -I gamma0.py final)              # re-close b_2, re-certify embedding -> gamma0_final.json (2.5 min)
python3 -I run_p1.py --out p1test --max-slabs 3                              # driver test, 3 slabs (as run)
python3 -I run_proto.py budget p1test 1 12,20                                # expected error at k = 12, 20: float tangent run (25 min)
python3 -I run_proto.py swing 45          # block and M resizing, t = 4 -> theta_max = 41 (25 slabs, K=16)
python3 -I run_proto.py approach 2 1500   # hls regenerated to k=12, then 2 slabs at k ~ 13 (M=1500, K=24), k_max at the end
python3 -I run_proto.py estimate 12       # P1 slab count and 1-core time (estimate 20: to k = 20)
```
Code: `slab.py` (engine), `run_proto.py` (driver for the tests), `run_p1.py` (production), `../t2/gamma0.py` (certificate; source and output set by `G0SNAP`/`G0FINAL`). Logs are in `logs/` and `p1test/`. Method and lemma: `T3_DESIGN.md` §5–6. Other jobs (T1) shared the cores, so wall times are up to 2× the CPU times quoted.

## 0. Error budget at the P1 → P2 switch (budget: physical C² ≤ 3e-14, C⁵ ≤ 1e-9)
- **Where the set grew.** Almost all of the earlier overshoot came from slab 1.
  - With the one-pass slaved Γ0, the slab's start Γ0 − ū(0) was 5.5e-18 at n = 250, its first free-start mode.
  - The a priori J majorant couples that box mode into the block edge n = 185 as J_blk ≈ 1.6e7·max Z_rest (9e-11 in sub-intervals 0–2, `logs/diag_slab1.txt`).
  - It was absorbed as η = 2.5e-19, so ρ went 8.1e-24 → **1.27e-20 (×1560)** and the hull Σn|G|w reached 1.2e-16.
- **Fix: second slaving pass** (§1). Modes ≥ 250 of Γ0 are taken from slab 1's own base, so the start mismatch is the double rounding, 9.4e-23. Measured with `run_p1.py` (§2):
  - Slab 1: ρ 8.14e-24 → **8.26e-24** (η 6.9e-24), hull 2.9e-21.
  - Slabs 2–3: ρ 8.27e-24 → 8.27e-24, with η ≈ 4.5e-25 per slab; the hull falls to 1.1e-21.
  - Cost: 2.5 min of re-certification; slab cost is unchanged.
- **Growth along P1** (`run_proto.py budget`). This is a float tangent run (ETDRK4) of the slab-1 set G·diag(w). The SVD at t = 0.03 keeps 69 directions with S_2/S_1 = 0.045: the absorbed edge components have decayed.

| t / k | t 0.0015 | t 0.03 | t 1 | t 3.9 | k 3.9 | k 9.6 | **k 12** | **k 20.3** |
|---|---|---|---|---|---|---|---|---|
| \|e\|_1 | 2.9e-21 | 8.1e-22 | 1.4e-20 | 2.9e-20 | 1.0e-18 | 7.5e-17 | 2.0e-16 | 2.6e-15 |

  - Phase II grows ×36, slowly and roughly ∝ t. This is not the time-shift mode, whose |∂_τ′θ|_1 falls ×20 there.
  - The swing plus approach give ×7e3 to k = 12; beyond that |e|_1 grows ∝ k^4.9 and C⁵ ∝ k^7.9.
  - The validated swing grew 2.2× faster than the tangent run over the same stretch (×288 vs ×130). Hence an allowance of **×3**, which also covers later absorptions (≈ 5% of w per slab, at the decaying block edge).
- **Expected at the switch** (tangent × 3; C² = (2π/L)|e|_1 + κ_max|e|_0, C⁵ = (2π/L)⁴|e|_4 globally, which bounds the value away from the tip):

| k_a | 12 | 20 | 28 (extrapolated) |
|---|---|---|---|
| C² (budget 3e-14) | **2.6e-16** (115×) | **3.3e-15** (9×) | 1.6e-14 (2×) |
| C⁵ (budget 1e-9) | **8e-14** | **4.8e-12** | 6e-11 |

  - The one-pass Γ0 would have given C² ≈ 1.3e-13 at k = 12 and 1.7e-12 at k = 20: that was the 10× overshoot.
- **More margin if k_a ≥ 24 is chosen** (not run). The dominant remaining direction is the b_2 closure interval (S_1 = 8.3e-24 ≈ its halfwidth).
  - Tighten the closure quadrature in gamma0.py (TAU 1e-22 → 1e-26) and impose the initial value for all modes in slab 1, which the two-pass Γ0 allows.
  - Expected ×20 lower, for about 10 min of re-certification and one slab.

## 1. Final Γ0 (task 1; option (c) chosen)
- **Why the T2 curve fails as a datum.** It is a float snapshot at t = 0.001. Its fast modes are off the slow manifold by 1e-10 (n = 186–350), 3e-12, 9e-14, 1e-16 and 1e-20 up to n = 1000. A Δt = t/2 slab from it had residual 2.0e-3, ρ → 5.3e-17 and box |·|_3 → 5.3e-6.
- **Options.**
  - (a) Graded slabs from the raw curve: h₁ ≈ 1e-3·t and about 9 extra doublings at M = 1000, +1.3 h.
  - (b) A later start reads the same float snapshots.
  - (c) Slaving: about 5 min once. **Chosen.**
- **Slaving, pass 1.** Modes 186..1000 come from a free-start Lobatto collocation over [0, 1.8t] in τ′ (ODE at all nodes for n ≥ 186; arb-polished; λ_186h = 29; free start fails at λh = 11 and works from 26). |Δb_n| ≤ 1.6e-10, Σ|Δb_n| ≈ 9e-9.
- **Slaving, pass 2.** Modes ≥ 250 come from slab 1's own base (h₁ = t/2, free start λh₁ > 26). The largest change is 5.5e-18, at n = 250. Modes ≤ 185 and > 1000 keep their T2 values; all b_n are exact doubles except b_2.
- **Closure and embedding** (`t2/gamma0_final.json`, source `gamma0s_src.npz`).
  - b_2 by arb Newton plus Krawczyk: halfwidth **8.14e-24**, contraction 3.0e-12. The centre moved −6.1e-14 from the T2 curve and −5.0e-22 from pass 1.
  - Embedding: 701 arcs, no intersecting pair, certified slot half-gap ≥ **0.0248199**, max arc turning 0.050, X(π) on the axis to ±3.3e-22.
  - Earlier certificates are kept as `gamma0_final_unslaved.json` and `gamma0_final_slaved1.json`. Γ0 is at t = 0.001, before the crossing at t ≈ 0.0028.

## 2. Initial phase (graded slabs from Γ0) and the driver test
- **Rule.** Δt = ½t with K = 24 (Chebyshev ratio 12.5–13.3, 10 required) while tκ⁴ < 4, i.e. to t ≈ 0.03 (9 slabs). Then h = 2/θ⁴.
- **Slab 1.** Slaved start for λh > 26, with a preset absorption cap.
- **Rows 3–5** are `run_p1.py --out p1test` with kill and resume during slab 3; the rows from `logs/initial_run*.txt` are for comparison.

| slab (M=1000) | Δt | M_a | residual max | ρ (max w) | η intra | hull Σn\|G\|w | box \|·\|_3 | CPU s |
|---|---|---|---|---|---|---|---|---|
| 1, T2 curve | 5.0e-4 | 185 | 2.0e-3 | 8.1e-24 → 5.3e-17 | 2.5e-14 | 7.4e-13 | → 5.3e-6 | 525 (+redo) |
| 1, pass-1 Γ0 | 5.0e-4 | 185 → 167 | 6.0e-12 | 8.1e-24 → 1.27e-20 | 2.5e-19 | 1.2e-16 | → 4.5e-14 | 525 (+redo) |
| **1, final Γ0** | 5.0e-4 | 185 → 167 | 6.0e-12 | 8.14e-24 → **8.26e-24** | 6.9e-24 | 2.9e-21 | → 4.8e-15 | 534 |
| **2** | 7.5e-4 | 167 → 151 | 8.6e-14 | → **8.27e-24** | 4.9e-25 | 1.5e-21 | → 1.1e-16 | 463 |
| **3** | 1.1e-3 | 151 → 136 | 2.1e-15 | → **8.27e-24** | 4.4e-25 | 1.1e-21 | → 2.0e-19 | 561 (+redo) |

- **Readings.**
  - The 6.0e-12 residual is the truncation part n > 1000; it does not reach the block. The dropped modes n > 1000 of Γ0 (box |·|_3 = 2.4e-12) decay within the first sub-interval.
  - Slab 3's redo came from the cap left by the old normalisation, and it exercised the retry path. The cap now uses the max radius.
- **Quadratures** (t from Γ0) after slab 3: t = 0.00236725442407002033 ± 2.8e-23, L = 103.473103510274117 ± 5.7e-19, X0 = −0.876331537214834 ± 6.9e-19.

## 3. Block and M resizing across the swing (task 2)
25 validated slabs from t = 4 (θ_max 16.3) to θ_max = 41.2, with K = 16. The test datum's (L, t) are the snapshot values, assigned after a Δt = 0.135 relaxation, so its t runs ≈ 0.02 ahead of feas2. This does not matter for the test.

| slabs | θ_max | M | block | ρ growth/slab | box \|·\|_3 | tail T | residual | wall s |
|---|---|---|---|---|---|---|---|---|
| 1–7 | 16→21 | 480 | 23→30 | 1.06–1.09 | 2–3e-20 | ≤ 2e-59 | ≤ 2e-18 | 42–50 |
| 8–17 | 21→31 | 480→560 | 30→43 | 1.10–1.155 | 3–9e-20 | ≤ 3e-45 | ≤ 5e-18 | 54–119 |
| 18–25 | 31→41 | 640→720→800 | 43→58 | 1.15 | 1–2e-19 | ≤ 5e-43 | ≤ 1e-17 | 93–241 |

- **Growth.** It is the true linearised growth: |B|w equals |diag B|w to 1e-3. ρ grows ×17.6 overall, from 1.0e-20 to 1.8e-19.
- **Resizing.** The block grew by 1–2 passengers per slab (shifted seeds) with no feedback. M grew 480 → 800, by +80 whenever the top 40 modes of x̄ exceeded 1e-30. After each growth the tail drops by 1e-5 to 1e-7.
- **Drops.** Tested in `smoke` and in the initial phase (185 → 167 → 151).
- **Slab length.** With K = 16 the adaptive test halved every c = 2 slab. K = 24 accepts c = 2 (§4), which would halve the swing's 25 slabs.

## 4. Approach phase (task 3)
- **Datum.** hls from p2_t4.59 with tail_tol 1e-15, regenerated to k = 12.002 (t = 4.5985973, N = 2048, 2 s). The consistent preparation took 151 s at M = 1000 and 351 s at M = 1500. Afterwards the datum sits at θ_max = 100.7 (k ≈ 13, T − t ≈ 1.9e-4).
- **Slab length.** K = 24. The adaptive test accepts c = 2 (ratio 15.1 vs 10), so Δt = 0.36(T − t).
- **Block.** M_a ∝ θ_max: +16 to +21 passengers per slab, with no feedback.

| run | slab: k | M_a | residual max | box \|·\|_3 / \|·\|_1 | tail T | ρ | hull Σn\|G\|w | X0 width | CPU s |
|---|---|---|---|---|---|---|---|---|---|
| M=1000 | 1: 13.0→14.5 | 142→158 | 3.0e-6 | 9.1e-8 / 8.6e-14 | 1e-14 | 1.001e-20 | 3.9e-20 | 4.4e-14 | ≈450 |
| M=1000 | 2: 14.5→16.2 | 158→177 | 1.1e-4 | 3.6e-6 / 3.4e-12 | 1.5e-11 | 1.003e-20 | 2.0e-18 | 8.7e-13 | ≈470 (+redo) |
| M=1000 | 3: 16.2→18.1 | 177→198 | 2.6e-3 | 1.0e-4 / 9.5e-11 | 1.1e-8 | 5.5e-20 | 9.1e-16 | 7.5e-12 | ≈520 (+redo) |
| **M=1500** | 1: 13.0→14.5 | 142→158 | **2.3e-13** | 4.5e-15 / 2.2e-21 | 2.3e-29 | 1.001e-20 | 3.8e-20 | 8.9e-21 | **≈720** |
| **M=1500** | 2: 14.5→16.2 | 158→177 | 4.8e-11 | 1.1e-12 / 4.5e-19 | 1.2e-24 | 1.001e-20 | 9.5e-20 | 1.5e-19 | **≈800** |

- **M sets the cost here.** The residual is the truncation part n > M, and it grows ×40 per slab once M is short.
  - M = 1000 still validates, but its box degrades fast with k.
  - M = 1500 is clean to k ≈ 14.5 and adequate to k ≈ 16. Required M ≈ 1300 at k = 12, ≈ 1500 at k = 14, ≈ 2000 at k = 16.
- **Block.** With M = 1500, ρ does not grow at k ≈ 13–16; the hull grows ×2.5–3.6 per slab.
- **Cost.** CPU matches the model 0.645·M_a·(M/360)^1.35 + 10·(M/360)^1.63 within 5%: 730 and 800 s predicted at M = 1500. Wall time was 850–1030 s under load.

## 5. Quadratures (task 4)
- **t, L, X0.** All three are propagated every slab by exact Chebyshev algebra in arb (T3_DESIGN §5). Widths at the end of the M = 1500 test: t ± 2.6e-20, L ± 3.8e-19, X0 ± 1.5e-19.
  - After the Γ0 slab the widths are t ± 2.9e-23, L ± 3.3e-18 and X0 ± 5.5e-16; X0 then contracts to 2.1e-16.
  - The X0 width follows the box |e|_3: 4.4e-14 at M = 1000, k ≈ 13.
- **k_max at the last slab end.** Grid with Taylor bounds, branch and bound, plus |e(h)|_1.
  - M = 1500: **k_max ∈ [16.225515619703, 16.225515619704]** at α = 3.1408, with |e|_1 ≤ 5.5e-19.
  - M = 1000: k_max ∈ [18.139092593859, 18.139092593884], with |e|_1 ≤ 9.5e-11.
  - Each enclosure takes 1.6 s.
- **Switch state for P2.** Every item is produced: x̄ (exact), G, w, box, T, t, L, X0 and k_max.

## 6. P1 cost; go/no-go (task 6)
`run_proto.py estimate 12|20` uses the slab-count rules above, the feas2 κ_max(t) history and the measured cost model (1 core):

| phase | range | slabs | M | M_a | hours |
|---|---|---|---|---|---|
| I graded, Δt = t/2, K = 24 | t 0.001 → 0.03 | 9 | 1000 | 185 → 81 | 0.95 |
| II h = 2/θ⁴ (∫κ⁴dt = 118) | t 0.03 → 4.0 | 59 | 930 → 360 | 74 → 23 | 1.54 |
| III swing (measured, K = 16) | θ 16 → 41 | 25 | 480 → 800 | 23 → 58 | 0.46 |
| IV approach (K = 24, c = 2) | k 5 → 12 (→ 20) | 10 (15) | 800 → 1200 (2000) | 55 → 132 (218) | 0.66 (2.15) |
| **total** | | **103 (108)** | | | **3.6 h (5.1 h)** |

- **Allowance.** ×1.5 for redos and for the M uncertainty in phase II.
- **Error.** At the switch the budget is met for k_a ≤ 28 (§0). The new Γ0 is the only change that was needed.
- **GO for P1 production** with `run_p1.py`. Its settings:
  - Γ0 = gamma0_final.json (two-pass slaved); graded slabs; preset slab-1 cap.
  - K = 24 throughout.
  - M rule: +80 while the top 40 modes exceed 1e-26, except in the graded phase; −80 while the top 120 are below 1e-31.
  - Adaptive absorption cap.
  - Stop at the target k_max.

## 7. Risks
- **R1. M near the singularity.** M ≈ 1200 at k = 12 and ≈ 2000 at k = 20, with cost ∝ M^1.35·M_a. Slabs at k = 20 cost ≈ 26 min each.
- **R2. M in phase II** is assumed (930 → 360). The feas2 snapshots there are N = 1024 at tail_tol 1e-10, so they cannot resolve 1e-30. Phase II carries ±50% (±0.8 h).
- **R3. Error estimate.** It comes from a float tangent run with a ×3 allowance calibrated on the swing. `status.json` reports the validated |e|_j, C² and C⁵ every slab. If C² heads past 1e-14 before k_a, apply the b_2-tightening knob of §0.
- **R4. Redos.** A missed absorption cap costs a full slab. The cap is now adaptive and normalised by the max radius (the earlier normalisation caused slab 3's redo), but it is untested on a long chain.
- **R5. Test data taken as exact.** The swing and approach start from float states, and the swing datum's (L, t) are offset by ≈ 0.02 in t. Production chains from Γ0, so this does not affect the proof.
- **R6. Speed.** The arb defect is 75–85% of the time, about M_a·(2K_J + K_P + 1) products; it parallelises over columns.
