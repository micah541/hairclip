# N1 and N2 (re-review REREVIEW_v33.md): chart term certified, step-4 environment recorded

Nov 14, 2026. Labels: [cert] = python-flint arb with outward rounding; [hyp] = hypothesis; [est] = estimate.

**Result.**
- **N1:** the chart term between the verifier's S1-chart variable w′ and (A)'s ray-graph variable is now an exact lemma, certified by a new step-4 v33 stage, `n1_chart`, and consumed by (A) (`t2_A_v35`). The typed 2e-7 is no longer used.
- **Margin moves:** every margin moved by at most 0.06%. No threshold fails, and K1–K3 hold. None of the stop triggers fired.
- **N2:** the driver now builds each stage's environment from scratch and pins it in the manifest. Two fresh full passes, rerun with explicit environments, agree byte for byte apart from timing fields. They also agree with the outputs produced before this fix, apart from timing and the records of which runs were reused.
- **Closed step 4:** `FINAL_manifest.json` (e01cb4da…) verifies with 0 mismatches. It was checked before the work, at the start and end of every driver pass, and again at the end.

## 1. Chart lemma (`fix2/n12/chart_lemma.tex`; balance check PASS, not compiled)

**Setting.**
- S1 is Γ₁ = (q + ip)e^{iθ₁}, with θ₁′ = k and the profile ODE.
- δ := θ₁ − Θ = −∫_σ^∞ k.
- The curve is Γ = Γ₁ + u·iT₁ (this is graph_eq.py's convention), with b = u′ = w′ + P and P = αψ̂₃′ + α²ĥ₂′.
- The ray chart is y = Re(e^{−iΘ}Γ), with slope f.

**Lemma C1 (exact, [paper]).**
- (a) Y(σ) = q cos δ − (p + u) sin δ.
  - e := Y − σ = (q − σ)cos δ − 2σ sin²(δ/2) − (p + u) sin δ.
  - Y′ = 1 − η, with η = 2 sin²(δ/2) + b sin δ + ku cos δ.
- (b) f(Y(σ)) = b + Φ, with Φ = [(1 − ku + b²) sin δ + kub cos δ] / [(1 − ku)cos δ − b sin δ].
  - This contains the nonlinear terms k·u·u′ and the shift u·sin δ.
- (c) y∂_y = (1 + ε)D, with D = σ∂_σ and 1 + ε = (1 + e/σ)/(1 − η).
- (d) **z^A_k∘Y = D^k(w′ + d) + χ_k**, where χ_k = (T^k − D^k)b + T^kΦ + p^k(U(σ) − U(Y)) and T = (1 + ε)D.
- (e) χ_k = A_k + B_k2·u‴ + B_k3·u⁗.
  - This form is affine in u‴ and u⁗ and vanishes for a ray.
  - B_22 does not depend on u″.
  - Every one of these properties is checked symbolically by sympy inside `n1_chart.py`.

**Lemma C2 (bounds, [cert]).**
- Envelopes:
  - |δ| ≤ k_b/μ, by convexity of E in the step-1 tail theorem;
  - |k^{(i)}| ≤ k_b, m_b, p_b/4, and k_b(σ+q_s)/4.
  - k_b, m_b, p_b and q_s come from `hstail.profile_bounds`.
- Verifier sups from fin12f4 (α₁ row) and |α| ≤ α₁ = 1.5e-4 [hyp].
- Split: χ_k = χ_k^add + β_k2·D²w′ + β_k3·D³w′, with the D w′ parts moved into χ^add through σ·sI2.

What (A) uses (all implemented in `assembly_rep.py --chart`):
- Ramp w-part: (√C_k + β_k2√C_2 + β_k3√C_3)².
- Ramp d-part: (√conv_k(d) + √(a₀X_k))².
- Both parts are multiplied by 1 + ρ_meas, which accounts for transferring the weight from y to σ.
- Φ: Z_w = (1 + β22)D²w′, and ∂_sZ_w = (1 + ε)((1 + β22)D³w′ + (Dβ22)D²w′).
- Sups get + sup χ^add.
- The core gets the χ^add terms, a shift term √(e_c Ê_{k+1}/27), a factor (1 + ρ_c), and a factor 1/(1 − β_c)².
- d's W-norms gain the ramp-end slivers [14 − 2⁻¹⁹, 14] and [28, 28 + 2⁻¹⁹]. These slivers are needed because Y⁻¹([14, 28]) ⊂ [14 − 2⁻¹⁹, 28 + 2⁻¹⁹]: |e| ≤ 4.39e-10, which is asserted.

**Datum.** (H2) and (H2f) are computed by step 4 directly in the ray chart, and f − U₁ = (f − P) + d exactly, so the datum needs no chart term. The ±2e-7 at `win.py:70` only widens those enclosures. It is harmless, but its comment is misleading. No closed file was edited.

## 2. Certified chart bounds against the reviewer's probe

| σ | quantity | probe (float) | certified (`n1_chart.json` → `probe`) |
|---|---|---|---|
| 14 | \|δ\|, level 0 | 2.2e-8 | 2.086e-8; sup χ^add_0 = 2.087e-8 |
| 14 | level 1 | 4.4e-7 | σk_b 4.434e-7; χ^add_1 4.444e-7 |
| 14 | level 2 | 9.4e-6 | σ²m_b 9.425e-6; χ^add_2 9.910e-6 |
| 14 | level 3 | 2.0e-4 | σ³p_b/4 2.003e-4; χ^add_3 2.307e-4 |
| 16 | levels 0–3 | about ×20 smaller | 8.54e-10, 2.17e-8, 5.75e-7, 1.58e-5 |
| 28 | level 3 | below 1e-13 | 6.46e-14 |

Other ramp values [cert]:
- X_k = ∫W(χ^add_k)² = (1.876e-21, 8.989e-19, 4.725e-16, 2.708e-13), as rates.
- β22 = 2.41e-9, β32 = 1.82e-7, β33 = 3.20e-9, Dβ22 = 8.50e-8.
- |ε| ≤ 7.83e-10, ρ_meas ≤ 1.05e-9.

Core [28, 60] [cert]:
- sup χ^add ≤ (3.95e-19, 2.15e-17, 1.19e-15, 6.63e-14);
- energies ≤ 6.6e-29;
- β ≤ 2.8e-17, ρ_c ≤ 4.5e-16.

d on the slivers [cert]: W-norms ≤ 3.3e-22; the sups equal b0_dprof's to 4 digits.

## 3. New stage and manifest
- Stage 19: `step4/p2/v33/n1_chart.py` (sha256 b42cca46be465b3794148870941412cf39a1a941e4bb48771e4a394219a81f4f), 8 s, about 100 MB. It writes `p2/out_v33/n1_chart.json` (sha256 341ed1f442ac64114bea6943cf37466a1515361bfcf838b41a7f7d860bf1ece0, identical in passes A and B).
- Driver: `final_tail_v33.py` (sha256 92308040bd5892aee5d617347b131d3e0165b7ce5623779d8cc71288a40a93b8).
- **`FINAL_manifest_v33.json` is now 20 stages, 289 paths, 37 scripts, 85 reused v3.2 pins.** sha256 **a0d348ec2dfb84bd3c62482517a0aa26f5d29ab3e4230081384054f395f8e357**. It supersedes 8f6f7f95…
- Backups:
  - `v33/final_tail_v33_pre_n12.py`;
  - `out_v33/FINAL_manifest_v33_pre_n12.json` (= 8f6f7f95);
  - `out_v33/pre_n12/` (all old outputs and sub-run JSONs);
  - `out_v33/logs_pre_n12/`;
  - the pass manifests `FINAL_manifest_v33_n12_pass0_reuse.json`, `…_passA.json` and `…_passB.json`, and the copies in `passA_out/`, `passA_runs/` and `passB_runs/`.
- Pass history:
  - pass 0 reused the sub-runs;
  - **passes A and B were clean.** The 64 win_arb runs and 32 h2f_win runs were moved aside, so everything was recomputed under the explicit environment.
  - Pass C (`--from 9`, with m9 sub-runs fresh again) was needed because `notes_v33` and `h4far_v33` pin `zone3_interface.json`, which N1 changed.

## 4. Old → new (v34/M14 → v35; all [cert])

| quantity | old | new |
|---|---|---|
| Ê, supersolution | (1.2020347e-6, 1.0439264e-4, 0.05668573) | (1.2017284e-6, 1.0440222e-4, 0.05669347) |
| Ê_ref at N2ref | (1.20181e-6, 1.04374e-4, 0.0566815) | (1.20150e-6, 1.04384e-4, 0.0566892) |
| FLUX0/K0, FLUX2/K2, β = 2, at N2ref | ×3.395, ×6.010 | ×3.397, ×6.013 |
| same at 4·N2ref (cited macros) | ×3.364, ×5.953 | ×3.365, ×5.956 |
| β = 5/4, at N2ref | ×5.404, ×9.568 | ×5.406, ×9.572 |
| S0v against 0.0153, N2ref / 4·N2ref | ×2.1220 / ×2.1160 | ×2.1221 / ×2.1162 |
| S1v against 0.6 | ×5.7883 / ×5.7485 | ×5.788 / ×5.7482 |
| core M1, N2ref / 4·N2ref | ×4.4362 / ×4.4092 | ×4.436 / ×4.409 |
| core M0 | ×5.4554 | ×5.4555 |
| ramp M1 / z1 / M0 | ×1.1847 / ×1.2014 / ×7.0203 | ×1.1847 / ×1.2014 / ×7.0205 |
| late M1 / M0 / M0z | ×2.3539 / ×4.8924 / ×12.992 | ×2.3538 / ×4.8925 / ×12.993 |
| present M1 / M0 | ×2.7709 / ×2.2251 | ×2.7708 / ×2.2251 |
| loop gain, total, levels 0, 1, 2 | 1.2251e-7, 0.037169, 0.039034 | 1.2252e-7, 0.037169, 0.039034 |
| B̂ | (4.873482e-6, 5.592692e-4, 0.3666978) | (4.873474e-6, 5.593227e-4, 0.3667418) |
| U-fine B̂ | (4.872222e-6, 2.626579e-4, 0.0258398) | (4.872209e-6, 2.626595e-4, 0.0258413) |
| H6 min / cone K3 / slope / t* / box | 1.8609 / 2.51070 / 2.46299 / 2.3019 / 1.2389 | 1.8609 / 2.51066 / 2.46299 / 2.3019 / 1.2389 |
| ΔR, n_III, τ_III | — | unchanged (zone3 differs only in the H6/U-fine entries listed) |
| K1, K2, K3 | true | true |

Notes on the table:
- Ê₀ decreases because the certified 2.09e-8 replaces the typed 2e-7. Ê₂ rises by 1.4e-4 relative.
- K6 is not recomputed here. Ê₀ moved −0.025%, so K6 ≈ ×1.190 [est].

## 5. Hashes (sha256)

**Chain outputs**
- `fix2/n12/t2_A_v35/out/A_interface.json`: 78fe2f8e38a340ba515fbc33cbab60cdd85231095af6ccd192e5f6707f457e8c. Two runs of `run_t2.sh` gave identical `SHA256SUMS`.
- `t4_cert/out/U_cert_v35.json` (pass 2): 960809fdc18dd67c1770fa5b77b2658081d0c81cb24f1b246f48ca4a583d5f54.
  - Pass 3 equals pass 2 apart from args.
  - Pass 1, with the v34 zone3 interface, also equals pass 2 apart from args.
- `t4_cert/out/U_cert_for_zone3_v35.json`: 22c13e4403a1b41ab09ab15a9d8852f153832a97742895757fdb1ce3fd0ca19f.
- `t5_zone3/out/zone3_interface.json`: 45e8d4dcaff2bd75fcea725f80e98625af70a40709fc2cf5d3df78b5becb3ea3.
  - The second zone3 run equals the first apart from args.
  - The previous version is kept as `zone3_interface.v34.json` (a381b1f2…), with `zone3.v34.log`.

**Scripts in v35**
- `assembly_rep.py`: 77d4abe6878b82e1daae794e9e23077ed1de43d411f92c91476d84407d5dc1db (backup `_pre_n1`).
- `export_t2.py`: 067e8c4d315c45cc2d07d8869eeb682f1ec797b1f8aae623f705f15d52da16f7 (backup `_pre_n1`).
- `run_t2.sh`: eca67cbd20dc0b8d0964511853bfcb12df4540bd7756ad36ab61c66385656eb1 (backup `_pre_n1`).
- `arlib.py`: only ROOT is changed (backup `_pre_n1`).
- The previous outputs are in `out_pre_n1/`.

**Checksum files.** `t4_cert/out/SHA256SUMS` and `t5_zone3/out/SHA256SUMS` were regenerated; the old copies are `SHA256SUMS.v34`.

## 6. Environment records (N2)
- **Child environment.** Every stage now runs with `env_child` = {PATH: <venv bin>:/usr/bin:/bin, LC_ALL: C.UTF-8, PYTHONHASHSEED: 0, PYTHONDONTWRITEBYTECODE: 1} plus the stage's own `env` (T6DIR, N3JSON, LO_NS, H4OUT, T9JSON as before). Nothing is inherited from the parent.
- **`env_read`.** Each stage record also pins `env_read`: the values in the child of 55 variables, found by grepping os.environ and getenv in the v33 stages and every module they import. The list includes S1TT, S1PROF, FORM_MULT, COHERENT, HSTEP, OMEGA_Q, T9TC, T9IFT and V33_VALIDATE.
  - S1TT, S1PROF and FORM_MULT are unset (null) for every stage.
  - Sub-processes inherit this environment. win_arb and h2f_win set T9TC and T9IFT themselves, and record them.
- **Manifest block.** The manifest has an `environment` block containing base_env, read_vars, the interpreter and flags, and `parent_values_not_passed`, which was {} in these runs.
- **Reuse check.** A `--from` run now asserts that each reused stage ran with the same env_child and env_read.
- **Caveat.** Children run with `-I`, which implies `-E`, so PYTHONHASHSEED is set and recorded but ignored by the child interpreter. Str hashing is therefore randomized. Byte-identity rests on the two clean passes.
- **Byte-identity.** All 289 paths in pass B equal pass A apart from `cpu` and `seconds` and the "(Ns)" timings in the logs. Against the pre-N2 outputs, the only differences are `cpu`, the reused-run lists and the input hashes that follow from them.
- **Bytecode.** No `.pyc` file was written by these runs. The existing `step4/p2/t9h/__pycache__/s1chart.cpython-314.pyc` (07:27, before this work) and the Nov 13 caches in `step3/`, `lin/` and `s1_cap/` are untouched; remove them if wanted.

## 7. Plainly stated limits
- **Inherited sliver convention.** As before, the verifier's sI2 (and so M1 = 28·sI2) is used on the sliver x ∈ [16.2491, 16.25]. This convention is inherited from (A), and the chart bounds use it too. The other sups (|w| and |w′|) are covered by S_ramp there.
- **Coarse chart conversion.** χ^add enters through the global table entry a₀. It does not use the localisation of the ramp cutoff, so it is conservative. Even so, the effect is ≤ 0.06%.
- **I3.** I3 gets the core factors too; it is not consumed downstream.

## 8. Documents to update
- **`fix1/b12_text/lemmaS_A_v33.tex`** (and `numbers_A_v33.tex`):
  - merge `chart_lemma.tex` into Lemma A:lem:conv;
  - delete "A level-0 chart allowance 2e-7 is added" (l.53);
  - change the "Not addressed by B0" sentence (l.64) and item B0 (l.315);
  - Φ, the ramp boxes and the B2′ text get the χ terms;
  - the macros take the new values from §4;
  - change the header hashes (A 78fe2f8e, U 960809fd, S4 a0d348ec with 20 stages).
- **`fix1/lemmaU_rev/lemmaU_v34.tex`**, `make_constants_v34.py` (set PIN = a0d348ec…) and `check_keys_v34.py`: the new sources and their numbers.
- **`t7_assembly/main_v3.tex`** and `main_v3_CHANGES.md`:
  - the (A) and U numbers;
  - "20 stages, 289 paths, 37 scripts, a0d348ec";
  - remove N1 from the open conditions.
- **`fix1/b12_text/step3_v33_main.tex`**: the hashes and the chart condition.
- **`t7_assembly/V33_INDEX.md`**, **`STEP4_RESULT_v33_addendum.md`**, **`step4/p2/v33/T6_REPORT.md`** and `M9_REPORT.md`: the new stage n1_chart and the explicit environment.
- **`ERRATUM.md`** and `erratum.tex`: v3.1 identified the S1 chart with the ray graph without proof; this is now closed.
- **`KITS_TODO.md`**: kit B/v33 must ship `n1_chart.py`, the new driver, `hstail.py` with `hs_uforms.pkl`, `lintail2.py`, `s1_cap/tail.py`, the fin12f4 rate table, and t2_A_v35.
- **`t4_cert/run_all.sh`**: its default should point at the v35 interface.
- **`REREVIEW_v33.md` §6:** N1 and N2 are closed.
