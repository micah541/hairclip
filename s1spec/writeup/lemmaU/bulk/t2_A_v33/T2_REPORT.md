> **Superseded (r2, Nov 13): task 2 was resumed under CONTRACTS_r2.md; the final report is `A2_REPORT.md`, the interface `out/A_interface.json`.** The text below is the r1 stop report, kept as a record (its probe files moved to `probe/`).

# Task 2: STOPPED at A4 (initial datum of (A)). Status report, Nov 13 2026

**Stop:** with the contract's provisional zone-III input, the honest initial datum of (A) pushes the result past the acceptance thresholds. At the tightened (H2): E₂ᴬ = 0.1685 against the 0.1 rule; S0v margin ×0.107 and core M1 margin ×1.14, both against ≥ 2 [cert, given the provisional n_III1]. A2 (R5′), A3, A5 and A6 are **not started**. Waiting for a decision.

## What was done (all in `bulk/t2_A/`; phase2a untouched)
- The phase2a `A_repair` chain was copied here with its paths repointed. It reproduces phase2a's outputs (args blocks excepted). The original outputs are saved in `ref_phase2a/`.
- **A1 H-J-W** (`rigor/hjw_rep.py` → `out/hjw_rep.json`) [cert]
  - Verifier file: `step3nl/verifier/rate_table_fin11_a0.00015.json`, sha256 in the JSON. The script asserts that (A)'s J's are its `unc` block rounded up.
  - **(W1) fails at k = 0.** On every cell L_up_k = Ly_hi > Ly1, because it is the barrier's upper end. The per-cell r1_k and Φ_k are not stored, so the check used P_k ≥ 0. The worst cell is k = 9.
  - **(W2) holds:** L_F = 2.4e-19 ≤ L* = 2.81e-15. It relies on the verifier's tail lemma, which was not re-run.
  - **(W3) fails as written. This is a new verifier error.** `closure.unconditional()` marches with `set_mode(delta, dropL=True)`: `rhs_iv.py:198` multiplies L by 0, and `unc.phi_val` has no constant term (`unc.py:17–27`). So the degree-0 flux floor F0 = 3.28e-22 per unit τ is omitted on march cells (`closure.py:182`). It enters only the tail slope (`closure.py:197`). The τ_a = τ₁ J's of v3.1 carry the same omission, relative size 4.8e-7.
  - Fallback: **κ_J = 0.042348**, made of a W1 part of 0.042347 and a W3 part (F0·τ_march/X_c) of 4.8e-7. This is ≤ 0.13, so it is not a trigger. It enters `assembly_rep.py --hjw` as J_m(1 + κ_J) for m = 0, 1, 2, 4.
- **A4 audit** (`rigor/ic_rep.py` → `out/ic_rep.json`)
  - Step 4 certifies (H2) on the shells [14, 60] only (`STEP4_RESULT.md` l.33; camp item 2). Remark 15 therefore applies.
  - The step-4 (H3) window energies at Ȳ = 1e-12 are 3.509e-8, 1.2615e-7 and 3.7723e-7 [cert]. They reproduce camp's item-3 limits to 1 ulp, on the safe side.
  - Three H3 settings: v32 (2e-7, 1e-5, 1e-5); x3 (levels 0, 1 divided by 3); min = 1.5 × step 4 = (5.27e-8, 1.90e-7, 5.66e-7).
- `assembly_rep.py` (backed up as `assembly_rep_pre_t2.py`) gained `--ic`, `--h3key`, `--hjw` and `--out`.

## Certified results (κ_J applied; outer term still at the phase2a H-O placeholders)

| zone-III datum n_III1 at τ₁ | H3 | E₀ᴬ, E₁ᴬ, E₂ᴬ (tight) | S0v | M0 core | M1 core |
|---|---|---|---|---|---|
| provisional (0.05, 0.12, 1.3) | v32 | 3.69e-3, 9.00e-3, **0.1685** | **0.107** | 0.549 | **1.138** |
| same | min | 3.69e-3, 8.95e-3, 0.1685 | 0.107 | 0.550 | 1.140 |
| 0.1 × provisional | v32 | 3.70e-4, 1.03e-3, 0.0822 | 0.329 | 1.49 | 2.36 |
| 0.01 × provisional | v32 | 3.81e-5, 2.35e-4, 0.0735 | 0.883 | 3.09 | 3.47 |
| 0.001 × provisional | v32 / min | 4.9e-6, 1.55e-4, 0.0727 | 1.74 / 1.98 | 4.49 / 4.77 | 3.84 / 4.19 |
| as two more (H3) windows (2·H3) | v32 / min | 1.60e-6, 1.66e-4, 0.0726 | 2.32 / 3.45 | 5.09 / 5.85 | 3.78 / 4.27 |
| 0 (H2 + 5·H3 only) | v32 / min | 1.20e-6, 1.46e-4, 0.0726 | 2.58 / 3.64 | 5.30 / 5.95 | 3.89 / 4.27 |

At the v3.1 (H2) with the provisional datum: E₂ᴬ = 0.3675 and M1 ×0.96.

**Reading.**
- The H3 part of C7 is harmless. At H3 = v32: E₂ᴬ = 0.0726, S0v ×2.58, M1 ×3.89. At H3 = min: S0v ×3.64, M1 ×4.27. These versus phase2a's 0.0696, ×4.61 and ×4.36.
- The zone-III part decides everything, and its binding effect is through E₀ᴬ and S0v, not through E₂ᴬ.
- S0v ≥ 2 needs the ω-weighted energy of the hand-off on [379.5, 700] to be about ≤ 2e-6 at level 0 (H3 = v32) or about ≤ 3.5e-6 (H3 = min) [est, interpolated between the certified rows].
- The provisional n_III is a supremum over the 5.2 units of τ of zone III, not a τ₁ value. If task 5 or 6 certifies the τ₁ far energy at the (H3) window level, everything passes.

## Files
- `rigor/hjw_rep.py`, `rigor/ic_rep.py`, `rigor/assembly_rep.py` (and `_pre_t2`)
- `prov/zone3_provisional.json`
- `out/hjw_rep.json`, `out/ic_rep.json`, `out/assembly_A4probe_{h3_v32,h3_min}.json`
- `probe/` (sensitivity runs)
- `out/SHA256SUMS_t2_stop`

`hjw_rep` and `ic_rep` were run twice; the reruns are byte-identical. `run_chain.sh` has not been extended yet.
