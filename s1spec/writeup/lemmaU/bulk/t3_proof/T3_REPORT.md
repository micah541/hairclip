# Task 3: Lemma U proof text. Done; no mathematical gap found that the certified inputs do not cover. 8 bookkeeping/review items (G1–G8), none changes a margin.

Opus, Nov 13 2026. Contract: `CONTRACTS.md` Task 3 as amended by `CONTRACTS_r2.md`, with the coordinator's v3.3 values (β_A = 2, (H2f) = (1e-6, 5e-6, 1e-5), verifier fin12f4). Only `t3_proof/` was written.

**LaTeX is not installed here.** `lemmaU.tex` was not compiled. `check_balance.py` passes: braces, environments, `$` (per paragraph), `\[ \]`, `\( \)`, labels and refs (36 labels, 25 refs). Every `\key{...}` cited in the text (99) was checked by script to resolve in the cited JSON.

## Cited sources (v3.3 only; the non-v33 files are superseded and not cited)
| tag | file | sha256 |
|---|---|---|
| A | `t2_A_v33/out/A_interface.json` (report `V33_A_REPORT.md`) | 5237a96a… |
| U | `t4_cert/out/U_cert_v33.json` | 42cf9617… |
| Z | `t5_zone3/out/zone3_interface.json` (rerun on the v33 U-fine B̂) | 0d6338d0… |

## Lemma list (`lemmaU.tex`)
| § | item | status |
|---|---|---|
| 1 | Def. 1 (one frozen profile U₁ for all centres, δU); **Lemma 2** (profile mismatch, c_δU) | paper + cert |
| 2 | **Lemma T** (transport identity with φ: co-moving cancellation, signs, b-term, outer b-part with r_χ = 253, sup\|χ′\| = 1.875/ln 1.5) | paper, full proof |
| 3 | **Lemma (launch data)**: (a) present centres from (H2), (H3) with pointwise Dφ; (b) late centres, Remark-2 form (no Δ_j; Δ^in on [28, 350]; neighbours via q^D) | paper + cert |
| 4 | **Lemma H-J-W** (cited, κ_J = 0.042351); **Lemma D** (window form); **Cor.** (layer cake); **Def.** sweeps; **Prop. X-term** (sharpened to the cutoff-free weight εΘφ, which S33/T223 need) | paper + cert |
| 5 | **Lemma (comparison)** for y′ + h ≤ σ + 2φ√y; **Lemma S** (sup-free cubic/quartic on an interval) | paper |
| 5 | **Prop. U** with full hypothesis list (U1)–(U6), proof in Steps 0–10: energy identity and base/cutoff split; point terms (Lemma P case (iii) throughout); base quadratic terms and Ĝ (general-weight L⁴, {2,3,3}, interpolation, dissipation budget with the a_min/10 reserve used locally on disjoint ramps); rates (R2), (R<2), Σ_rs; inner commutators K (no IBP, top pair included, K_S = 1); cubic classes ((2,2,2) in c₃₃; (2,2,3) by Cauchy–Schwarz + sup-free L⁴); foot terms and S33; conversion to (A)'s data; couplings; **outer ramp with task 5's hand-derived Sobolev/Jensen/Cauchy–Schwarz steps** (time measure ε_e dτ = (r_χλ₁)⁻⁴ dt); forcing Φ_j; summation | paper + cert |
| 5 | **§ M1z**: sup\|z₁\| ≤ 0.534 asserted per region (core, zones I/II via S₁(2B̂), zone III via b₁ + pU₁) | cert |
| 6 | **Lemma W** (window Sobolev, W ≥ w_min, windows with 1.5ρ ≤ r_χλ₁) and the covering of zones I/II | paper |
| 7 | **Lemma N** (H-O derived); bootstrap set BS_U (B2′ with β_A = 2, written parametrically; B4 with far centres ≥ 490); order table; **Prop. closure** (K1)–(K6) with the open–closed argument written out (joint continuity of (τ, ρ) ↦ 𝔈^ρ_j(τ), continuity of the launch set, newly launched centres handled by uniform continuity at the launch boundary, no hypothesis-level launch bound needed); **Def.** loop gain | paper + cert |
| 8 | **Prop. U-fine** (statement only; K₂, K₃, H6 from task 5) | cert |
| 9–11 | certified-inputs table ([cert]/[paper]/[hyp]); hypotheses; gaps | — |

Headline numbers quoted (all [cert]): B̂ = (4.8732e-6, 5.5667e-4, 0.36561); Ê = (1.1879e-6, 1.0391e-4, 0.056589); late M1 ×2.358, M0 ×4.900; present M1 ×2.775, M0 ×2.225; loop gain ≤ 0.0390; rate (R2) ×28.23; flux at β_A = 2 ×3.505 / ×6.196; Δ^R₂ = 0.01148; H6 min ×1.86; U-fine B̂ = (4.8717e-6, 2.6258e-4, 0.02578).

## Gaps found
**None at the design level, and none that the certified inputs leave uncovered.** Items for other tasks (§11 of the text):
- **G1 (task 4/2), label.** `U_cert_v33.json: supersolution.beta_A` = 1.25, because the A export still carries `beta_A: "5/4"` although v3.3 adopts β_A = 2. Lemma U does not use β_A; the text cites the flux rows at β_A = 2 (`core.flux_margin.2`). Relabel.
- **G2 (task 4), provenance.** `t4lib.verifier_sups` reads ā, β̄ from `step3nl/` (fin11), not `step3nl_v33/`. The six-digit values coincide (`V33_VERIFIER.md` l.32–33), so nothing changes; repoint for the v3.3 record.
- **G3 (task 4), rounding.** The b-term constants are used as C₁ = 5, C₂ = 8; the proof needs Γ₁(Θ) + …, with Γ₁(Θ) ≤ 4 + 5.4e-22. Effect ~1e-25 relative on a 5e-7 term.
- **G4 (task 4), acceptance row.** (R<2) level-0 margin ×4.43 < 5 (report line, not a trigger). Strict; enters only the coupling surplus.
- **G5 (task 5), to confirm.** Step 8 needs N₃² as an *integral* bound over [t₁, T_τ′] normalised by the full horizon λ₁⁴q, so the same prefactor q r_χ⁻⁴ applies to the actual (shorter) interval. Z describes N₃² as a "time average over the horizon"; confirm that it is this form (a sup-in-time bound would also do). Size if it were not: none expected (the horizon normalisation only overestimates).
- **G6 (task 8), code check.** Step 0 asserts that no Lemma-U step integrates by parts a term containing (η²)‴ or (χ²)‴, so no point term arises at y = 56 or at r = r_χλ₁ (Lemma P case (iii)). `u5_enum.json` supports this for the inner ramp; confirm for `ramp_U.py` and the global Ĝ pipeline.
- **G7 (task 8), design choices to review:** (A)'s ramp box set to the core box (ramp M1 ×1.184); far family restricted to centres ≥ 490 (zone-III b₁ ×1.238); M1z = 0.534 per region.
- **G8 (task 7), wording** of the main_v2 conditionality (see §10).

Small corrections made in the text against the r2 statements (no number changes):
- Prop. 8 is stated for the cutoff-free weight εΘφ (the proof gives it; S33 and T223 need it).
- The couplings use D(28)·B_k/s_k; the certified assembly also adds the ramp term J·E^k_k and a factor φ(28) (conservative, noted).
- P₂ carries a factor D(28)φ(28) that is not needed if R covers εΘφ (conservative, noted).
- Lemma-U's Ĝ uses one box on the whole support, so Lemma L of the statements is not needed inside Lemma U.

## Still hypotheses (§10 of the text)
- **(H2), tightened** (1e-8, 1e-6, 1e-3) on [14, 60]: step 4 meets it ×25.7, but it must be added to the step-4 manifest (task 6).
- **(H3)** (2e-7, 1e-5, 1e-5), unchanged from v3.2.
- **(H2f)** (1e-6, 5e-6, 1e-5) on [60, 700] at t₁: **not certified**; task 6 (approved). Probe, base curve only: about ×3.2, ×5, ×4 [est].
- **α₁** = 1.5e-4.
- **Fixed-L cone-cylinder data** ((H5), Lemma O H1/H4/H7, L = 50λ₁, centres [300, 380]): new geometry, step-4 stages `t9_cone_kch`, `h4loc` and the Lemma O H1 check to be rerun (task 6).
- **(P1)–(P3)** and the graph property (literature, as v3.2), and the verifier's tail lemma (W2), not rerun.

## Files (`/home/claude/s1spec/writeup/lemmaU/bulk/t3_proof/`)
- `lemmaU.tex`: the proof text (495 lines, not compiled).
- `check_balance.py`: the balance checker (`python3 -I check_balance.py lemmaU.tex` gives PASS).
- `constants_needed.json`: 76 rows of {symbol, file, key, producer_task, used_in, label}, with the sha256 of the three sources. Every key was checked to resolve.
- `T3_REPORT.md`: this report.
