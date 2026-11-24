# F12_REPORT: review items F1 and F2 (standalone pinned check, step-4 v33 stage 20 `f12_sliver`)

Nov 14, 2026. Source of the items: `writeup/lemmaU/bulk/t8_review/CHECK_N1_N5.md` (F1, F2). Script `p2/v33/f12_sliver.py` (sha256 124f13be…). Outputs go only to `p2/out_v33/f12/`. Summary: `f12/f12_sliver.json` (1fbc0f80…).

**Verdict: F1 and F2 are closed.** Both corrected constants are certified in arb. Every claim built on them was recomputed and still holds: **326 of 326 checks pass**. No acceptance threshold fails, and no quoted (rounded) number is exceeded. No existing pinned output changed. The pinned *numbers* move in the 8th to 12th significant digit (§4).

## 1. F1: sliver sups G_i on the extended domain [cert]

Domain: y ∈ [14 − 2⁻¹⁹, 16.25 − X0lo], with X0lo = 2.2490776505. X0lo is ≤ −σ_tip (from the s1chart ball) and ≤ the code's X0 = 2.2490776515, so the domain covers both. The formula and the 64 cells are those of `assembly_arb.sliver_G`. The old G is recomputed by the same code and asserted bitwise equal (arb midpoint and radius) to `sliver_G()` inside every assembly run.

| i | G_i old [14, 14.00092] | G_i new [14 − 2⁻¹⁹, …] | new/old | needed domain (14 − 4.388e-10) / old |
|---|---|---|---|---|
| 0 | 1.95383e-14 | 1.97004e-14 | 1.00830 | 1.0000062 |
| 1 | 1.18623e-9 | 1.19361e-9 | 1.00622 | 1.0000047 |
| 2 | 5.40125e-5 | 5.42362e-5 | 1.00414 | 1.0000031 |
| 3 | 1.63941 | 1.64280 | 1.00207 | 1.0000016 |

These agree with the reviewer's probe (0.21–0.83 % for the stated domain, 4.7e-7 to 1.9e-6 for the needed one). The recomputation uses the **stated** (larger) domain.

## 2. F2: sup|w″| on the sliver x ∈ [16.249, 16.25] [cert]

This sliver contains every x = σ − σ_tip with σ ≥ 14 − 2⁻¹⁹, so it is wider than the reviewer's [16.24907, 16.25].

**Argument.** No Taylor step was needed.
- On the sliver, w′ = g·v, where g = 1/ζ and v = ζw′ is the verifier's cut-off variable (RIGOR_PASS.md item 3). Hence |w″| ≤ sup|g|·|v′| + sup|g_x|·|v|.
- The verifier's far sups are the (I2_W) bound of `agent_sob/SOB_NOTES.md`, sI_k = √(K_W(Ψ1V² + 2VV′)), applied to the full weighted norms of v and v′ (`closure.py farsup`).
- (I2_W) holds at every x in [6, R], with the factors 1/W(x) and Ψ1(x) = sup over [x, R] of ψ′. On the sliver:
  - 1/W(x) ≤ (1 + 16.25²)² = 70258.1 ≤ K_W = 2.323e6. W = (1+x²)⁻² there and W is decreasing (ψ′ ≥ 0.130 on [6, 80]).
  - ψ′ ≤ 0.2453 on the sliver, and Ψ1 from 10.25 is 1.305 (`sob_constants_P2.json`).
- So sup|v′| ≤ sI2 and sup|v| ≤ sI1 on the sliver, with the verifier's own numbers.

**Numbers:**
- sup|g − 1| ≤ 2.70e-14 and sup|g_x| ≤ 1.08e-10 on the sliver.
- **sI2′ = (1 + 2.70e-14)·sI2 + 1.08e-10·sI1 = 0.015456178888449044**, against sI2 = 0.015456178888323341 (α₁ row; relative +8.13e-12). The 1.25e-4 row is also +8.13e-12.
- M1′ = 28·sI2′ = 0.43277300888 ≤ 0.4328, the 4-digit value that (A) reads. Asserted, so (A)'s ramp boxes are unchanged.

## 3. How the claims were re-checked (no rerun of the pinned chain)

Each step below is a scratch recomputation from pinned inputs only. Each harness was first validated by reproducing its pinned file exactly (all keys except `args`):

| step | what was recomputed | validation |
|---|---|---|
| χ bounds | `n1_chart.py` executed in memory with the single line `W2SUP = arb(UPd(FM['sI2']))` replaced by sI2′ | repro equals `n1_chart.json` |
| (A) | the 7 `run_t2.sh` runs of `assembly_rep.py` (as `__main__`), with `sliver_G` replaced by the new G and `--chart` set to the new n1 file | repro equals each pinned `asm_*.json` |
| export | `export_t2.main` with the asm paths and the output directory redirected | repro equals `A_interface.json` (78fe2f8e…) |
| Lemma U | `assemble_U.py --A <new> --out <abs f12 path>` | repro equals `U_cert_v35.json` (960809fd…) |
| zone III | `zone3.py` executed in memory with the output path replaced, using the new U-fine B̂ | repro equals `zone3_interface.json` (45e8d4dc…) |

Lemma U was then run a second time with the new zone3 file, and pass 2 equals pass 1. The fixed point holds because the zone3 quantities that U consumes (ΔR, ΔR_b, n_III, box_III_380/700, τ_III) do not depend on B̂ and are byte-equal.

The stage also asserts that the closed trees are unchanged: hashes of 1007 files under t2_A_v35, t4_cert, t5_zone3, lemmaU_rev, the verifier and agent_sob, with no new files and no `.pyc`.

## 4. Re-checked claims (certified lower bounds for margins, upper bounds for bounds; tight setting unless stated)

| claim | pinned | new | threshold / quote | status |
|---|---|---|---|---|
| ramp M1 | ×1.184734096 | ×1.184734096 (identical) | > 1; quoted ×1.1847 | holds |
| ramp z1 / ramp M0 | ×1.2014269 / ×7.0204547 | identical | > 1 | holds |
| S0v (0.02) / S0v vs 0.0153 | ×2.77401145064 / ×2.12211875974 | ×2.77401145050 / ×2.12211875963 | > 1 | holds |
| flux, β = 2, at N2_ref (FLUX0, FLUX2) | ×3.39665691442, ×6.01283394122 | ×3.39665691436, ×6.01283393983 | > 1 | holds |
| flux, β = 2, at 4·N2_ref | ×3.36510972230, ×5.95559013661 | ×3.36510972224, ×5.95559013525 | quoted ×3.365, ×5.955 | holds |
| K6 | ×1.18983006271 | ×1.18983006271 (+5e-14) | ≥ 1.189 | holds |
| H6 (zone III, min margin) | ×1.86086713742 | identical (the binding row does not use B̂) | quoted ×1.860 | holds |
| cone K3 (H6 row) | ×2.51066253229 | ×2.51066253188 | quoted ×2.510 | holds |
| U: late M1 / M0 / M0z | ×2.35381253 / ×4.89245701 / ×12.9928857 | ×2.35381253 / ×4.89245701 / ×12.9928857 (−2e-9 relative) | > 1 | holds |
| U: present M1 / M0 / M0z min | ×2.77079544 / ×2.22513585 / ×14.3773783 | ×2.77079544 / ×2.22513585 / ×14.3773783 | > 1 | holds |
| U: K1, K2, K3, N2 ≤ 4·N2_ref, M1z, U-fine requirement; loop gain | all true; 0.0390342 | all true; 0.0390342 | < 1 | holds |
| χ (n1): e_ramp < 2⁻¹⁹, sup χ₀^add ≤ 2e-7, β < 1e-6 | | sup χ^add changes only in the 16th digit (levels 2, 3) | | holds |
| check_keys_v34 side conditions: Ê_ref4 dominates Ê, D^c ≤ D^src, √2·S1v ≤ 0.1475, B2′ ≥ 2.116 | | all hold (B2′ ×2.11619) | | holds |
| M3 / M6 (`n35_asserts`) | | ramp margins > 1 in both settings; N̄²/(4·N2_ref) < 1 | | holds |
| every asm run (7 runs, both settings): margins > 1, flux, claims | | no new failing claim; every margin above 1 stays above 1 | | holds |

**Pinned numbers that are exceeded.** These are bound values, not thresholds:

| quantity | pinned | new | relative change |
|---|---|---|---|
| Ê_ref, A_interface `EA_ref` | (1.2015018217681e-6, 1.04383558809e-4, 0.0566892258905) | (…7681149e-6, 1.04383558831e-4, 0.0566892262276) | up to +5.9e-9 |
| Ê_ref4, 4·N2_ref (`asm_design_N2ref4` E) | 0.05824635883 | 0.05824635917 | +5.8e-9 |
| U supersolution Ê | 0.0566934746 | 0.0566934749 | +6.0e-9 |
| B̂₂ | 0.36674177770 | 0.36674177897 | +3.5e-9 |
| U-fine B̂₂ | 0.025841280428 | 0.025841280445 | +6.3e-10 |

- Largest relative change in any pinned file: 4.8e-8 (a few internal asm parts such as JG and `J_route`); 1.6e-8 in A_interface; 7.6e-9 in U_cert; 6.3e-10 in zone3.
- **Quoted numbers.** All 184 quoted numbers of `lemmaU_v34.tex` (constants_needed_v34 `numbers`) are still rounded outward on the new files. Every macro of `numbers_A_v33.tex` (regenerated with `n35_gen_numbers.py`) is identical. The generic leaf diff over n1, the 7 asm runs, A_interface, U_cert and zone3 finds no number whose outward rounding to 3, 4, 5 or 6 significant digits changes.
- **Consequence.** No document quotes a value that is now too small. The only places with more digits are the JSON files themselves, which stay pinned as they are.

## 5. Mechanics and manifests

- **Backups:** `p2/v33/final_tail_v33_pre_f12.py` (92308040…) and `p2/out_v33/FINAL_manifest_v33_pre_f12.json` (a0d348ec…, the 20-stage manifest).
- **Driver:** `final_tail_v33.py` (now 57c20a44…) gains one stage, index 20, `f12_sliver`. It has 67 inputs (all read-only, all existing pins) and 26 outputs (25 JSON under `out_v33/f12/` plus `logs/f12_sliver.log`). It also adds a `results.F12_sliver` block.
  - The environment is the N2 one: `env_child` = {PATH, LC_ALL=C.UTF-8, PYTHONHASHSEED=0, PYTHONDONTWRITEBYTECODE=1}. `env_read` is over the unchanged READ_VARS. Every variable that the stage's modules read is already in READ_VARS.
  - Children run with `python3 -I -B`, one at a time, under nice 19. The whole stage took 28 s.
  - Not pinned: the child logs in `out_v33/f12/logs/`, because some print wall times.
- **Runs:** two clean passes of `--from f12_sliver`, each with no partial file and `out_v33/f12/` deleted first.
  - Stages 0–19 were reused with hash, definition and environment checks; their records are identical to the pre-f12 manifest.
  - Pass 1 and pass 2 agree byte for byte on all 26 stage outputs and on every manifest field apart from `seconds`.
  - Pass 1 manifest: `FINAL_manifest_v33_f12_pass1.json` (f324adc0…). Pass 2 log: `out_v33/final_tail_v33_f12_pass2.log`.
- **New manifest: `p2/out_v33/FINAL_manifest_v33.json`, sha256 8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec.** It has 21 stages, 366 paths, 39 scripts and 85 reused v3.2 pins. All hashes were verified at the end, `inputs_with_two_versions` is empty, and `all_items_ok` is true.
- **Closed `FINAL_manifest.json`** (e01cb4da…): verified before (`verify_v32_manifest.py`, 0 mismatches), at the driver's start and end in both passes (523 paths, 97 upstream, 26 scripts), and again afterwards (0 mismatches).

## 6. Follow-ups (docs not edited)

- `check_keys_v34.py` now fails only its sha check of S4: `constants_needed_v34.json` pins the old manifest a0d348ec…. Everything else passes (K6 1.18983, B2′ 2.11619). Fix it by updating the PIN in `make_constants_v34.py` (l.189) to 8bc47e83… and regenerating.
- The texts that cite "20 stages, a0d348ec…" can cite 21 stages and 8bc47e83…, and can drop the F1/F2 caveats: `step3_v33_main.tex`, `continuation_v33.tex`, `lemmaS_A_v33.tex`, `numbers_A_v33.tex`, `main_v3.tex`, `lemmaU_v34.tex`, `V33_INDEX.md`, `BULK.md`, `STEP4_RESULT_v33_addendum.md`, `KITS_TODO.md`.
- `chart_lemma.tex` Remark (4) and Lemma C:lem:bd still describe the sliver as an "inherited convention". It can now cite §2: sI2′ ≥ sup|w″| on the sliver, certified.
- No pinned number quoted in any document is exceeded, so no number needs to change. If the full JSON values are ever quoted, the new ones are in `out_v33/f12/new/`.
