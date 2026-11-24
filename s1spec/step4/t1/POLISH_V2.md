# POLISH_V2: base re-seeded from P1's exact centre (slab 95), arb-128 polished store

**Commands** (in `t1/`; 1 process each, `nice -n 10`, OMP/OPENBLAS/MKL threads = 1):
```
python3 seed_p1.py 95 80                                   # exact seed -> seed_p1_95.pkl (4 s)
python3 run3.py P11v2 2.31e4 seed_p1_95.pkl                # float base -> slabs_P11v2.pkl (181 slabs, 483 s)
python3 -u polish_v2.py P11v2 180 0,1,2,3,4,5,180 --seed seed_p1_95.pkl --store polished_v2   # [--resume I]
python3 verify_v2.py 0 1 160 170 180                       # fresh-process reload check
python3 se_v2.py polished_v2 [i0:i1]                        # (SE) tilt/offset, arb 192 -> se_polished_v2.json
STORE=polished_v2 python3 -I f_err2_v2.py ../t3/p1/ckpt/slab_00095.json 0    # p2/f_err2.py, store path only
```
The polish was interrupted by a server restart after slab 164. It was resumed with `--resume 165`, which reloads stored slab 164 losslessly and builds a fresh float Jacobian. Newton paths after slab 164 therefore differ from an uninterrupted run, but only within the tolerance |E| < 1e-26.

**Seed** (`seed_p1.py`). The input is P1 ckpt slab 95: θ = −π/2 + α + Σ u_n sin nα with M = 2800 exact dyadic u_n, plus L, X0 and t. The seed is computed in arb at 192 bits.
- Batch evaluation: the sum is split as Σ_q e^{iqBα} Σ_r c_{qB+r} e^{irα} (B = 64), using arb_mat products. Horner in acb on |z| = 1 inflates radii by (|cos|+|sin|)^n, which reaches 1e152 at n = 2800.
- Values: k = 26.2537786005622, λ = K0/κ_tip, ℓ̃ = L/(2λ) = 411.97336, mesh N = 120 (n = 840 + 4).
- Projection: the nodal κ̃^{(0..4)} are exact at the knots; bubbles match Legendre moments 0, 1 of the remainder (80-point Gauss).
- |κ_h − κ_P1| ≤ 4.7e-17 (x ≤ 40) and ≤ 6.3e-20 (x > 40), at 300 random points.
- The far point Im Γ̃(ℓ̃) = 3.4e-21 is removed by a zero-mean bump ε w′ (physical s = 12, width 1.5, ε = 7.6e-23). The bump increment uses 0 for the fixed DOF; a version that used their values (K0) moved the turning.
- Stored as exact 128-bit dyadics (radius 0); scalars c = 0 and t = P1's t.
- Seed tilt = −3.0e-38 and offset = −1.5e-36 (192-bit evaluation of the stored seed).

**Store** `polished_v2/`: 181 slabs (0–180), τ ∈ [0, 27.15], k = 26.2538 → 23718.7, N = 120 → 171 (n + 4 = 844 → 1201), 53 MB. Format and `manifest.json` are as in POLISH_STORE.md, and `polstore.load_slab` works unchanged. Slab 0, row 0 is the exact seed. Max polish residual is 9.7e-27. The reload check on slabs 0, 1, 160, 170 and 180 gives bit-identical norms.

**Residual norms** (24 sampled slabs: 0–5 and every 10th to 180; sup over 6 times + 8 nodes). **No misses.**
| slab (k) | R2 | R1 | Rs | RIs | RF0 | RF1 | RF2 | RF3 | RF4 |
|---|---|---|---|---|---|---|---|---|---|
| target | 2e-14 | 5e-13 | 5e-11 | 1e-10 | 1e-7 | 1e-7 | 1e-7 | 1e-7 | 1e-5 |
| 0 (26.3) | 2.5e-15 | 6.9e-15 | 7.4e-14 | 5.7e-18 | 1.1e-11 | 2.6e-11 | 1.6e-10 | 3.1e-9 | 7.7e-8 |
| 1 (27.3) | 1.8e-16 | 6.3e-15 | 4.0e-14 | 2.7e-18 | 1.1e-11 | 2.5e-11 | 1.7e-10 | 2.9e-9 | 7.8e-8 |
| 10 (38.4) | 1.8e-17 | 4.8e-16 | 3.8e-14 | 2.3e-18 | 1.2e-12 | 2.8e-12 | 1.9e-11 | 3.1e-10 | 1.2e-8 |
| 90 (781) | 1.5e-17 | 6.3e-18 | 3.1e-14 | 2.1e-18 | 2.9e-13 | 7.1e-13 | 5.1e-12 | 8.0e-11 | 2.3e-9 |
| 180 (22846) | 1.5e-17 | 2.9e-19 | 3.1e-14 | 2.1e-18 | 1.4e-13 | 7.5e-13 | 1.2e-11 | 1.9e-10 | 5.7e-9 |
| max (all sampled) | 2.5e-15 ×8 | 6.9e-15 | 7.4e-14 | 5.7e-18 | 3.4e-11 | 8.4e-11 | 5.9e-10 | 9.6e-9 ×10 | 2.6e-7 ×38 |

- There is no start-up layer: slab 0 meets every target. Seeding at slab 93 was therefore not needed, and P2 can start exactly at the seed (base τ = 0 is P1 slab 95).
- The RF3/RF4 maxima are at slab 160, which starts with an outer element insertion.

**Switch mismatch** (float probe `f_err2_v2.py`, P1 slab 95 vs base slab 0):
- Y = 9.7e-15, with A = 9.7e-15, B = 1.3e-14, V0 = 4.5e-16, clock −2.9e-13.
- This is the probe's float floor: dksup = 2.8e-14 comes from float sums of the 2800-term series, while the arb value is |κ_h − κ_P1| ≤ 4.7e-17.
- The old base gave 2.8e-9.

**(SE)** (`se_v2.py`, arb 192 on the stored midpoints, all 9 rows of every slab):
- Slab ends: |tilt| ≤ 4.2e-28. It accumulates at the polish tolerance, about 2e-30 per slab.
- Slab ends: the offset h = Im Γ̃(ℓ̃) drifts to −1.7e-16 rescaled (−1.1e-20 physical, 6e-22 of |Re Γ̃(ℓ̃)|) by k = 2.37e4. This is residual-driven and amplified by ¼−ã.
- **Inside slabs** (Radau stage rows): tilt ≤ 3.2e-18 (slab 0; ≤ 2.8e-20 for k ≥ 100) and offset ≤ 1.05e-15 rescaled (6.1e-17 physical).
- Cause: ∫κ̃ and Im Γ̃(ℓ̃) are nonlinear in (DOF, ℓ̃) on the moving outer mesh. Leg₀ gives d∫κ̃/dτ = 0 only at the collocation nodes, so the stage values carry O(Δτ^{s+1}) and the slab ends are superconvergent. The stage tilt at slab 0 is the same with Newton tolerance 1e-31 (`t_tilt_tol.py`).
- Consequence: exact (SE) for all τ needs the static blended correction of STEP4_RESULT §2(ii). It must absorb both quantities: tilt ≤ 3.2e-18 and offset ≤ 1.05e-15 rescaled.
