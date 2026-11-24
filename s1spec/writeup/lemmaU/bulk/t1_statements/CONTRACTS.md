# Task contracts for bulk tasks 2–5

Bulk task 1 (Fable), Nov 13 2026. Statements and notation: `STATEMENTS.tex` (same directory; read it first, it is the specification). Labels as there: [cert], [est], [rev], [paper].

## Rules for all four tasks

- **Write only in your own task directory** (`bulk/t2_A/`, `t3_proof/`, `t4_cert/`, `t5_zone3/`). Copy scripts you need to change; never edit `phase1/`, `phase2a/`, `lemmaS*`, `step3_proof*`, `step4/`.
- **Arithmetic:** python-flint `arb`, outward rounding; `math.nextafter` when a double is exported. No sampling-based claim. Suprema by cell covers with Taylor or monotonicity enclosures, as `ramp_rep.py` does.
- **Records:** every output JSON has an `args` block (argv, env, input paths with sha256). Two consecutive runs must be byte-identical apart from timing. Write `out/SHA256SUMS`.
- **Run:** `PYTHONHASHSEED=0 nice -n 19 python3 -B ...`, one process, under 800 MB.
- **Interface numbers** are exported twice: as an arb string and as an outward-rounded double. Every number has a label field: `cert`, `est` or `hyp`.
- **Inputs from another task are read from a path given on the command line**, so the stage can be rerun when that task's file changes. Until the other task has delivered, use the provisional file named below and mark the output `provisional: true`.
- **Stop and report to the parent** (do not work around) when a stop trigger or a "report" item below occurs.
- **Do not use:** the per-centre profile U_ρ (use U₁, Definition 1); a horizon-free D; pointwise bounds on z₂ or z₃ on [350, 700]; any Lemma Y number.

## Who needs what from whom

| producer | file | consumers |
|---|---|---|
| task 2 | `t2_A/out/A_interface.json` | 3 (names), 4 (numbers), 6, 7 |
| task 5 | `t5_zone3/out/zone3_interface.json` | 2 (n_III, box), 4 (Δ^R, box decision), 6 (hand-off thresholds) |
| task 4 | `t4_cert/out/U_cert.json` | 3 (names), 6, 7, 8 |
| task 3 | `t3_proof/lemmaU.tex`, `t3_proof/constants_needed.json` | 7, 8 |

No task blocks another: task 2 exports (A) as an affine function of the outer data, so it does not wait for tasks 4 or 5.

Provisional inputs: `phase2a/A_repair/export_for_U.json` (for task 4 until task 2 delivers); Δ^R = (9.9e-8, 1.0e-6, 0.0161) [est, typed] and n_III = (0.05, 0.12, 1.3) [est] (until task 5 delivers).

---

## Task 2 — (A): R5′, window export, H-J-W, audits (Opus; `t2_A/`)

**Start from** a copy of `phase2a/A_repair/` (scripts `rigor/*.py`, `run_chain.sh`, `base_ibp/`, `base_v31H2/`).

### Work items

| # | item | kind | spec |
|---|---|---|---|
| A1 | H-J-W | mechanical | `STATEMENTS.tex` Def. 3, Lemma 4. From the verifier result that feeds (A)'s J's (locate the JSON under `s1spec/step3nl/`; record its sha256), compute per march cell k: L_up_k, c_k, P_k; check (W1) L_up_k ≤ Ly1 + 2P_k for every k; check (W2); audit (W3): is the floor term F0 inside `orc.rate` on march cells? If (W1) fails, compute κ_J and rerun the chain with J_m(1 + κ_J). |
| A2 | R5′ | design, then mechanical | Proposition 14. New stage `rigor/outer_L2.py` replacing `outer_rep.py`: enumerate every cutoff term on [350, 700] (cutoff parts m ≥ 1 of all pairs with all jets, interpolation remainder, L⁴ and S33 foot terms) with the class enumerator of `ramp_rep.py` mirrored to the decreasing cutoff η_o = 1 − S(t); reduce each class with tools T1–T7 and Lemma P case (i) at level 2, case (iii) at levels 0, 1; assert that every class ends in N̄₀², N̄₁², N̄₂² or the reserved dissipation. Output κ^o_jk (3×3). In `weights_rep.py`: remove the `lam` and `eps4` additions (no cut at 350); outer segment with box (1/2, 1), U ≤ 0.162 and B_c lowered by δ^o_j; rerun the segment check. Write the localisation at y = 350 by Lemma L(ii) with the transition on [175, 350]. |
| A3 | point terms at y = 28 | mechanical | Lemma P. Confirm the signs of the table entries (2,2,0):4 on (j+1, j+1) and (3,2,0):2 on (j, j+1) on the stored forms. Levels 0, 1: bound the cutoff part of the m = 3 monomial of pair (j, j+1) directly by data (Cauchy–Schwarz in space-time with conv_j, conv_{j+1}), no integration by parts. Level 2: keep both integrations by parts and record that the net point term is −g_c z₂(28)² ≤ 0. Report the change of E₀, E₁. |
| A4 | initial datum | audit, then mechanical | Remark 15. Find on which y-range step 4 certifies (H2) (`step4/STEP4_RESULT.md`, the H1/Ȳ stages, `gonogo/t7_step4_H2_tight.*`). If it is [14, 60]: set E^A_j(τ₁) ≤ H2_j + 5·H3_j + (28/379.5)·n_III1_j, with H3 as a command-line parameter. Run three settings: H3 = (2e-7, 1e-5, 1e-5); H3 with levels 0, 1 tightened ×3; the smallest H3 that step 4 meets by ×1.5 (read step 4's values, do not rerun step 4). |
| A5 | interface at y = 28; segments | audit, then mechanical | Lemma L. State how v3.1 splits the interpolation and the L⁴ identity between the ramp box (0.05, 0.9) and the core box. If the convention is not covered by Lemma L(i), apply L(ii) (transition of Ỹ on [28, 56]) and L(iii), and recompute the core Ĝ, γ′. Redo the segments (60, 220] and [220, 350] with the box (0.15, 0.52), and with the fallback (0.20, 0.60); U ≤ 1/20 and 0.08 as now. The export keeps a_c = (1 + 0.11²)⁻² because [28, 56] lies in the core box. |
| A6 | export | mechanical | `out/A_interface.json`, below. |

### `A_interface.json` (per setting `tight`, `v31`; keys under `res.<setting>`)

| key | meaning |
|---|---|
| `EA0` [3] | Ê^A without the outer source |
| `kappa_o` [3][3] | κ^o_jk of Proposition 14; ô_j = Σ_k κ^o_jk N̄_k² |
| `gp` [3] | γ′_j (lower bounds) |
| `EA_of_N2` | the rule Ê^A_j = EA0_j + ô_j/γ′_j (string), and `EA_ref` [3] at `N2_ref` |
| `N2_ref` [3] | (0.1, 0.2, 5.0), the reference outer data [est]; also evaluate at 4 × `N2_ref` (`EA_ref4`) |
| `Dsrc` [2], `r_export` [2], `amin_c`, `rhoJ`, `F2` [2], `coup` [2] | ingredients of Lemma D |
| `C0_rule`, `C1_rule` (strings); `C0_ref` [2], `C1_0` [2], `C1_ref` [2] | C0_k = (Ê^A_k + Dsrc_k)/(r_k a_c); C1_k = C1_0_k + ô_k/(r_k a_c) |
| `IC` | `E_tau1` [3] used, its three parts, the H3 setting, the y-range of step 4's (H2) |
| `HJW` | `W1_ok`, `W2_ok`, `W3_note`, `kappa_J`, `worst_cell`, verifier file sha256 |
| `pointterm` | `gc_28`, routes per level, `dE` [3] |
| `interface28` | convention found, `Ghat_core` [3], `gp` before/after |
| `segments` | per segment and per box: `ok`, rows, retained |
| `core` | `S0v`, `S1v`, margins S0v, M0, M1, Ĝ₂₂; `K0`, `K2` at Ê^A and at 2Ê^A against FLUX0, FLUX2 |
| `classes_outer` | the class list on [350, 700] with the tool chain used for each |

### Acceptance

| quantity | threshold | now |
|---|---|---|
| E₂ᴬ at `N2_ref`, tightened (H2) | ≤ 0.1 (rule); report if > 0.08 | 0.0696 |
| R5′ addition to E₂ᴬ at 4 × `N2_ref` | ≤ 0.01 (**stop trigger** if larger) | 1e-6 at placeholders |
| every class on [350, 700] closed by T1–T7, P | yes (**report** if not; do not bound pointwise) | — |
| H-J-W | (W1)–(W3) hold, or κ_J ≤ 0.13 with all thresholds still met (**stop trigger** otherwise) | unchecked |
| core M1 margin; S0v margin | ≥ 2.0 each | 4.36; 4.61 |
| core Ĝ₂₂ margin after A5 | ≥ 1.5 (report if < 1.7) | 1.94 |
| K0, K2 at 2Ê^A | ≤ FLUX0, FLUX2 with margin ≥ 2 | ×4.9, ×5.3 at v3.1 (H2) |
| retained r₀, r₁ on [28, 56] | ≥ 1.249 a_c | 1.249 |
| segment conditions, both boxes | all true | true for (0.11, 0.52) |
| reruns | byte-identical | — |

Deliver `A2_REPORT.md`: results table at both (H2) settings, what each of A1–A5 changed, v3.1 file:line for every new error, the list of numbers still [est] or [hyp].

---

## Task 3 — Lemma U proof text (Opus; `t3_proof/`)

Not mechanical. Output `lemmaU.tex` (cannot be compiled here: run a balance script and say so) and `constants_needed.json`.

### Content, in this order

1. Setting with **v = f − U₁ for all centres** (Definition 1, Remark 2); the weight with φ; energies; zones.
2. Lemma T (transport identity with φ) as in `U_ramp.tex` §2, with full proof.
3. Launch data: present centres as in `launch_lemma.tex` Lemma 2(a) with pointwise Dφ; late centres in the form of Remark 2 (no Δ_j; the term Δ^in on [28, 350]).
4. Lemma D in window form, Corollary 6, Definition 7, Proposition 8 (copy the statements; expand the proofs).
5. Prop. U of `U_ramp.tex` §3 with 𝓓_k in place of S_X D_k, proof written out term by term: global quadratic terms and Ĝ; cutoff commutators and K (top pair included, no integration by parts on the ramp); cubic classes with Lemma C; foot terms and S33; couplings; b-terms; forcing Φ_j; rates under (R2), (R<2). No Lemma Y.
6. Window Sobolev: from 𝔈 ≤ B̂ to sup|v| and sup|y v_y| on a window with the non-flat W ≥ w_min.
7. Lemma N (H-O derived) and the bootstrap: B3U, the order table, Proposition 16 with the open–closed argument written against (P1)–(P3), including the joint continuity of (τ, ρ) ↦ 𝔈^ρ_j(τ) and of the launch set.
8. U-fine: statement only, with the constants it takes from task 5.

### Rules

- Every constant is a named symbol with a JSON key; no number is typed without a label and a source file. `constants_needed.json` lists `{symbol, key, producer_task, used_in}` for each.
- Mark each step [paper] or cite the script that certifies it.
- Where the proof needs something the statements do not give, write it as a numbered "gap for task 4/5/8" at the end; do not assume it.

### Acceptance

- Each inequality used in Prop. U and Proposition 16 is either proved in the text or tied to one key of `A_interface.json`, `U_cert.json` or `zone3_interface.json`.
- Balance script passes (braces, environments, `$`, `\[ \]`, labels).
- The key names agree with the tables in this file.

---

## Task 4 — arb certification of Lemma U (Opus; reruns by Sonnet; `t4_cert/`)

### Work items

| # | item | kind |
|---|---|---|
| U1 | Γ_m, κ̃, κ̃₂ of Θφ (φ = exp(K_b e^{−s}), K_b = 1e-3); Ĝ_jk at the box (0.15, 0.52), U_max 0.09, M0z 0.16, dfrac 2/25, foot adds, U_j included; g, D(56), D(28), θ(28); (R2) with y_R = 56 and (R<2). Same for the fallback (0.20, 0.60), M0z 0.20. | mechanical (rerun `phase1/ghat`, `phase2a/U_ramp/ghat_design_arb.py` with φ) |
| U2 | K_jk on [28, 56] with φ at the chosen g; for the fallback, extend the certified g-range to 6.3e5. footw, c₃₃, and c₂₂₃ in arb (Lemma C: sup\|c′_nj\|, L_max, sup\|c\|, the quartic c_J part). | mechanical (rerun P2), c₂₂₃ new |
| U3 | Sweeps as arb enclosures: Σ_in, Σ_E, Σ_rs^(j), S_X(y_c), Σ_X(y_c) (Definition 7) for y_c = 56 and on a cover of y_c ∈ [56, 253]; prove monotonicity in y_c or take the cover sup. Tail of the x-integral bounded analytically (W ≤ e^{4(σ + a_W)}). | mechanical |
| U4 | Launch constants with pointwise Dφ: C_A^D, q^D, w_min, c₂(y_c), c₃(y_c), Δ^in_j; c_δU in arb; Φ_j; the window Sobolev constants for W ≥ w_min. | mechanical |
| U5 | Class-enumerator cross-check: run the enumerator of `A_repair/rigor/ramp_rep.py` on [28, 56] with the weight Θφ and compare with P2's direct classes. | design |
| U6 | Assembly `assemble_U.py`: reads `A_interface.json` and `zone3_interface.json`; finds (Ê^A, B̂) satisfying (K1)–(K3) of Proposition 16 with β = 2 (iterate from β = 0, then verify the inequalities in arb); loop gain G_j (Definition 17); margins. Levels in the order 0, 1, 2; late centre and a cover of present centres y_c ∈ [56, 126.7]. | design |

### `U_cert.json`

| key | meaning |
|---|---|
| `box` | {M0, M1, Umax, M0z, design: "data56"} and `fallback` with the same fields |
| `Gamma`, `Ghat` [3][3], `g`, `D56`, `D28`, `theta28`, `rate_margin`, `R_lt2_margin` [2] | U1 |
| `K` (per level, per column), `footw`, `c33`, `c223`, `g_range` | U2 |
| `sweeps` | {Sigma_in, Sigma_E, Sigma_rs [3], S_X, Sigma_X}, each as {yc56, cover_sup} |
| `launch` | {CA_D, qD, wmin, c2, c3, Delta_in [3], c_dU, Phi [3], sobolev} |
| `crosscheck` | class lists of both enumerators, max relative difference of weights, unmatched classes |
| `supersolution` | {EA [3], Bhat [3], N2 [3], C0 [2], C1 [2], K1_ok, K2_ok, K3_ok, slack per inequality} |
| `loop_gain` [3] | G_j, with the split neighbour / through (A) |
| `margins` | {late_M1, late_M0, present_M1_min, present_M0_min, at_yc}, each with the box value it is measured against |
| `parts_B2_late` | launch, inflow, ramp E, X-terms, rs, foot+S33, couplings, b, outer, forcing |

### Acceptance

| quantity | threshold | now |
|---|---|---|
| certified late M1 | ≥ 1.5 (**stop trigger** below); report if < 1.8 | ×2.28 [rev], before (C6), (C7) |
| late M0, present M0, present M1 | ≥ 2.0 each against the chosen box | ×4.2, ×2.85, ×2.84 [rev] |
| rate margin (R2) | ≥ 10 | ×29.07 [cert] without φ |
| (R<2) margins | ≥ 5 | not certified |
| θ(28) | ≤ 3/5 | holds at both boxes |
| loop gain G_j | ≤ 0.5 | ≈ 0.04 [est] |
| φ effect on Γ_m, Ĝ, K | ≤ 1e-3 relative (report otherwise) | ≤ 1e-4 [est] |
| S33 + foot at level 2 | ≤ 1e-2 of B̂₂ | 2e-4 [est] |
| cross-check | no class of the A-enumerator missing from P2; weights within 5% or explained | — |
| slope part C1·Σ_X | report its share of 𝓓_k | 1e-4 [rev] |

If the certified late M1 is below 1.8, also report the row with H3 levels 0, 1 tightened (task 2, A4).

Deliver `U4_REPORT.md` with the table of margins at both boxes and the list of anything still [est].

---

## Task 5 — zone III, Lemma R ramp, cone family (Opus; `t5_zone3/`)

### Work items

| # | item | kind |
|---|---|---|
| Z1 | **The deciding quantity.** b₀ = sup \|f\|, b₁ = sup \|r ∂_r f\|, b_z = sup \|f − U₁\| over 190λ₁ ≤ r ≤ 379.5λ₁ and t ∈ [t₁, t₁ + λ₁⁴/(1 − 4ā)], f the slope over the S₁ ray. Derive from D5a (\|g_x\| ≤ C₁√t, \|g_xx\| ≤ C₂√t per conclusion region, scaled (50/L_m)²), the hand-off bounds on h₀ (thresholds of (H3), (H5), not step 4's actual values), and the chart conversion cylinder → S₁ ray (rotation by the cylinder axis; `lemmaS_F5_geom.json`). All in arb. | design |
| Z2 | The same on 190λ₁ ≤ r ≤ 700λ₁ against (1/2, 1): b₀⁷⁰⁰, b₁⁷⁰⁰. | mechanical after Z1 |
| Z3 | n_III_k ≥ sup_τ ∫ over {350 ≤ y ≤ 700, r ≥ 190λ₁} of (z_k^A)² ds, k = 0, 1, 2, with z^A = ∂_s^k (f − U(α(τ), y)); and n_III1_k at τ₁ on [379.5, 700]λ₁. z₂ needs h_xx pointwise (D5a gives it); no z₃. | mechanical after Z1 |
| Z4 | Lemma R ramp [190, 379.5]λ₁ for the weight Ψ_ρ ≤ D(190)φ·W: certified Δ^R_j replacing the typed (9.9e-8, 1.0e-6, 0.0161); add the b-term part of Lemma T(iv). | mechanical (rerun `lemmaS_v32/ramp_cert`, `lemmaR_interval.py`) |
| Z5 | Cone family with fixed L = 50λ₁, X = 3, centres ≥ 300λ₁ (plan §4): covering, Lemma R with unit-50 members, the H6 geometry, and the constants U-fine needs (K₂, K₃ on [150, 250)λ₁ as functions of B̂). | design |
| Z6 | List of hand-off thresholds that task 6 must certify on the step-4 data: existing ones in L = 50 units, and any new one introduced in Z1–Z3. | mechanical |

### `zone3_interface.json`

| key | meaning |
|---|---|
| `box_III_380` | {b0, b1, bz} |
| `decision` | "primary" if b0 ≤ 0.15, b1 ≤ 0.52, bz ≤ 0.16; "fallback" if b0 ≤ 0.20, b1 ≤ 0.60, bz ≤ 0.20; else "fail" |
| `box_III_700` | {b0, b1}, `HO1_ok` (≤ (1/2, 1)) |
| `n_III` [3], `n_III1` [3], `tau_III` | Z3 |
| `DeltaR` [3], `DeltaR_b` [3] | Z4 |
| `cone_family` | {L, X, centres, covering, tstar_margin, R4, H6 margins, slope check} |
| `ufine` | {K2, K3 rules, values at `Bhat_ref` = B̂ from `U_cert.json` or (5e-6, 1e-3, 0.53) provisional} |
| `handoff_thresholds` | list of {name, region, value, new: bool, step4_margin_if_known} |
| `inputs_used` | which hypotheses (B4, H1–H7, H9, thresholds) each number depends on |

### Acceptance

| quantity | threshold | expected |
|---|---|---|
| `decision` | "primary"; "fallback" is acceptable and is reported at once; "fail" is a **stop trigger** | b₀ ≈ 0.09, b₁ ≈ 0.26 [est] |
| margin of the chosen box in Z1 | report if < ×1.1 | ×1.6 on b₀ [rev] |
| `HO1_ok` | true | — |
| Δ^R₂ | ≤ 0.03 (report if larger: it is 4% of B̂₂ per 0.016) | 0.0161 [est] |
| fixed-L family: t* margin, H6 margins, cone slope | ≥ 1.3 each (report the thinnest) | ×2.30; v3.2 cone slope ×1.40 |
| K₂, K₃ of U-fine against 0.18 | ≥ 1.5 | ×3.0, ×2.3 [est] |
| every number | depends only on B4, the hand-off thresholds and D1, never on B2′ or B3U | — |

Deliver `Z5_REPORT.md` with the decision in its first line.

---

## Stop triggers (from `BULK.md`), with the exact test

| trigger | test | task |
|---|---|---|
| certified late M1 < ×1.5 | `U_cert.json: margins.late_M1` | 4 |
| D5a cannot give zone III ≤ (0.20, 0.60) | `zone3_interface.json: decision == "fail"` | 5 |
| R5 in L² form adds > 0.01 to E₂ᴬ | `EA_ref4[2] − EA0[2] > 0.01`, or a class not closed by the tools | 2 |
| H-J fails on windows | (W1)–(W3) fail and κ_J > 0.13 or a threshold breaks with J(1 + κ_J) | 2 |
