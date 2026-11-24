**Verdict: N1–N5 are closed. No inequality fails and nothing blocking was introduced.** Four small items remain (F1–F4 below); F1 is a real but tiny rigor gap left by the N1 fix (effect 2e-8 relative on Ê₂ [probe]).

Check of the N1–N5 fixes (Fable), Nov 14 2026. Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4/p2`. Labels: **[checked]** = I read the code or text; **[rerun]** = I ran it (scratch copy or read-only checker, scratch in `B/t8_review/scratch3/`); [probe] = my scratch computation, arb but not a pinned stage.

## 1. Item table

| Item | Status | Evidence |
|---|---|---|
| **N1** chart term | **Closed**, with F1 and F2 below | Lemma C1 re-derived by hand: Y, Y′ = 1 − η, Φ, y∂_y = (1+ε)D, χ_k all correct. `n1_chart.py` builds exactly these with sympy and asserts affinity and vanishing for a ray. The arb boxes are right (see §2). Rerun three times under random hash seeds: all numbers equal the pinned `n1_chart.json` (341ed1f4…) [rerun]. `assembly_rep.py --chart` puts the term in every consumer (see §2). Scratch rerun of the 4·N2ref assembly reproduces the pinned output to double precision [rerun]. Merged text `chart_lemma_body.tex` differs from `chart_lemma.tex` only by the two renamings. |
| **N2** env records | **Closed** | `final_tail_v33.py:161–168` builds the child environment from scratch and records `env_child` and `env_read` per stage; `--from` reuse asserts they match. Manifest a0d348ec…: 20 stages, 289 paths, 37 scripts, 85 reused pins, all hashes verified by my own script [rerun]. Pass-A copies against current outputs: 115 JSON files compared, 23 differ, only in `cpu` and in pinned input hashes [rerun]. |
| **N3** stale citations | **Closed**, two leftovers (F4) | No stale hash or stage count in `lemmaS_A_v33.tex`, `step3_v33_main.tex`, `continuation_v33.tex`, `numbers_A_v33.tex`, `lemmaU_v34.tex`, `main_v3.tex`, `erratum.tex` [grep]. `numbers_A_v33.tex`: 40 of 42 macros equal the regenerated values; the two differences are the documented tweaks [rerun]. `check_keys_v34.py`: PASS, 132 keys, 184 numbers outward [rerun]. |
| **N4** M3, M6, M7 | **Closed**, with F3 | `n35_asserts.py` (scratch copy) PASS: ramp margins ×1.18473, ×1.20143, ×7.02046 (certified lower bounds, both settings); N̄²/(4·N2ref) = 0.179, 0.413, 0.252. The three restored scripts are byte-identical to their `_pre_n35` copies, and the `_n35assert` copies differ by the one assert line each [rerun]. `main_v3.tex:309` defines R(τ); `:312` says Lemma N gives (H-O1) and (H-O2′), not (O3); `:315` quotes the erratum values (0.25558, ×2.723, ×1.146, ×1.124, ×1.185, all found in `ERRATUM.md`) [checked]. |
| **N5** criterion, Z2 | **Closed**, pending the authors | `continuation_v33.tex:16, 51, 107`, `main_v3.tex:122`: DKS 2002 Thm 3.1 and/or Chou 2003, marked "to be confirmed by the authors"; Wheeler dropped. `step3_v33_main.tex:59`: Z2 "for all τ ≥ τ₁"; `zone3_interface.json` `inputs_used.box_III_700` lists only B4, D1, (H3), (H5) [checked]. |
| Closed step 4 | **Verifies** | `verify_v32_manifest.py`: `FINAL_manifest.json` e01cb4da…, 0 mismatches [rerun]. No `.pyc` written by my runs. |

## 2. N1 in detail

- **Exact relation** [checked by hand]: with Γ = (q + i(p+u))e^{iθ₁}, Y = q cos δ − (p+u) sin δ; Γ′ = T₁(1 − ku + ib) gives f∘Y = tan(δ + arctan(b/(1−ku))) = b + Φ with the stated Φ; χ_k = (T^k − D^k)b + T^kΦ + p^k(U(σ) − U(Y)). The split u‴ = (D²w′ − Dw′)/σ² + P″, u⁗ = (D³w′ − 3D²w′ + 2Dw′)/σ³ + P‴ is right, and so are β_k2 ≤ |B_k2|/σ² + 3|B_k3|/σ³, β_k3 = |B_k3|/σ³.
- **Arb bounds** [checked]: symmetric balls for δ, k, k′, k″, k‴, u, u′, u″; envelopes at the left cell end (decreasing); |δ| ≤ k_b/μ holds because E′ = μ exactly (A_E = ¾·4^(−1/3)); the tail theorem is based at σ = 8, so 14 − 2⁻¹⁹ is inside its range. The weight-transfer factor covers Y′, (σ/Y)^(2m) and W(Y+X₀)/W(σ+X₀); with it, a₀ ≤ a₁/196 is valid because Y ≥ 14.
- **Consumers** [checked in `assembly_rep.py`]: all ramp conversions (inflow, commutator, top, classes, Young users) through `cv`; Φ with Z_w = (1+β₂₂)D²w′ and its T-derivative; the cubic d-part with sup|D²d + χ₂^add|; foot and S33; ramp boxes; core sups with the end-shift term; fluxes K0, K2 (levels 0–2, so β₂₂ suffices) and I3 (level 3, uses max(β₂₂, β₃)). The datum needs no chart term. I found no consumer left out.

## 3. Questions asked

- **`-I` and PYTHONHASHSEED: acceptable.** Hash randomization can change only iteration order, never the validity of an enclosure. Identity across passes run with different random seeds is a stronger test than a fixed seed. One wart: `n1_chart.json` records `PYTHONHASHSEED: "0"` in its own `args.env` although the interpreter ignored it; the manifest note says so, the stage file does not.
- **zone3-pin relaxation in `check_keys_v34.py`: legitimate.** `zone3_interface.json` (45e8d4dc…) pins `U_cert_v35.json` (960809fd…), so U cannot also pin that zone3 file; the cycle has no hash fixed point. `U_cert_v35_pass3.json` pins 45e8d4dc… and equals pass 2 in every key except `args` [rerun]. The checker tests both conditions live. Cleaner would be to cite pass 3 as the U file, but it is not required.

## 4. New findings, ranked

- **F1 (minor, rigor; from the N1 fix). `sliver_G` does not cover the extended ramp domain.**
  - The ζ < 1 sups G_i are taken on y ∈ [14, 14.0009] (`base_ibp/rigor/assembly_arb.py:58–75`). After N1 the w-part is integrated over σ ∈ Y⁻¹([14, 28]), which reaches down to 14 − 4.39e-10 (the text and code use 14 − 2⁻¹⁹).
  - G_i must grow by 4.7e-7 to 1.9e-6 relative for the needed domain, and by 0.21% to 0.83% for the stated one [probe, `scratch3/gsliver.py`]. G₃ is 1.64, not small.
  - Effect: with all G_i × 1.0083, Ê₂ moves 0.058246359 → 0.058246360 at 4·N2ref and no printed margin changes [probe, `scratch3/asm_G1.0083.log`].
  - Fix: compute G on [14 − 2⁻¹⁹, 14.0009] in the v35 assembly. That is a rerun of t2 (40 s), U, zone3 and the step-4 tail from stage 9, with new hashes; or a pinned standalone stage showing the claims hold with the larger G. Coordinator's choice.
- **F2 (minor; inherited, slightly widened by N1). sI2 on the sliver.** sup|w″| ≤ sI2 is certified on x ≥ 16.25 but used on x ∈ [16.24907, 16.25] (ramp M1 box, and now the χ bounds). The authors state this in N12_REPORT §7; neither earlier review examined it. From the formula sI2² = K_W(…), I expect a correction of order 1e-4 relative [est] against a ramp margin of ×1.18. It needs one sentence of proof or a certified constant.
- **F3 (minor, docs; from the coordinator's restore).** Four places still say the M3/M6 asserts are inside the scripts: `lemmaS_A_v33.tex:19, 87, 317`, `lemmaU_v34.tex:530`, and `V33_INDEX.md:316`. On disk the asserts are in `fix2/n35_asserts.py` (and the `_n35assert` copies). Reword to cite `n35_asserts.py`.
- **F4 (cosmetic).**
  - `ERRATUM.md:36` still says "19 stages, sha256 8f6f7f95".
  - `ic_rep.py:20` still names the float outputs (restored on purpose; it is a label string inside the pinned A_interface, no number depends on it).
  - Visible TODOs: `main_v3.tex:315`, `ERRATUM.md:361`, `erratum.tex:240` (the (O3) values, not used by v3.3).
  - Still at 19 stages / ×15.59: `S4/v33/T6_REPORT.md`, `M9_REPORT.md`, `BULK.md:120`. `REREVIEW_v33.md` §6 still lists N1, N2 as open.
  - `check_keys_v34.py` prints 132 cited keys; N35_REPORT says 134.

## 5. Status of v3.3

Step 3 v3.3 is closed as a certified chain: N1–N5 are done, both manifests verify, and the thinnest margins are unchanged (ramp M1 ×1.1847, K6 ×1.1898, H6 ×1.860).

It is conditional on:
1. the authors' literature confirmations: (P1′) the L² blow-up criterion (DKS 2002 Thm 3.1 and/or Chou 2003, for immersed non-embedded curves), (P3′) interior regularity, and (P2′) with the frame conventions;
2. a LaTeX compile (nothing has been compiled; only brace balance is checked);
3. beyond those two: **F1 and F2**, two sliver formalities of width 4e-10 and 9e-4 with effects of about 2e-8 and 1e-4 against margins of ×1.18 or more. Neither can change a verdict, but by the project's rule ("no 'closed' without a script that checks it") they should be certified or fixed before the word "unconditional" is used.

Not re-examined here: everything listed as "not checked" in REREVIEW_v33 §8, the items the authors labelled [hyp] in the texts, and the step-4 hand-off beyond the hash checks.
