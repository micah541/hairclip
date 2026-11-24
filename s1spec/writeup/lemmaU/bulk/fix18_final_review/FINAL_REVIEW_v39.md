# FINAL_REVIEW_v39: independent final review of step 3 (nonlinear trapping), text set v39

Nov 20-21, 2026. Fable, independent reviewer. The proof tree was read-only; everything written is in this directory (`fix18_final_review/`).

- Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`; `T` = `B/fix17_texts_v39`; `V` = `/home/claude/s1spec/step3nl_v4c`; `N` = `B/fix13_v37chain`; `M38` = `/home/claude/s1spec/step4/p2/out_v38/FINAL_manifest_v38.json`.
- Labels: **[cert]** arb, outward rounding, on pinned inputs taken as given (mine: `scratch/recompute.py`); **[read]** verified by reading the cited file:line; **[hand]** my derivation, not machine-checked; **[est]** estimate.
- Runs: one process, nice 19, `python3 -B`, a few seconds each, under 100 MB.

## 0. Verdict

**CLOSES WITH FIXES (none blocking), modulo the listed external and author items.**

- The assembled argument is coherent end to end: one bootstrap set, an order table without circularity, and a correct continuation to T < ∞. I found no failed inequality, no constant rounded the wrong way, and no dropped term in what I checked.
- I found **one gap in the written closure proof** (F-1: the member B1(a,b) has no gap, so the "Open" step as written does not cover it). It is repaired with material already certified; no computation is needed. I classify it non-blocking, but it must be fixed before the text is called closed.
- I found **one defect that makes the text set uncompilable and leaves two hypothesis thresholds of the theorem unprinted** (F-2: eight macros used but never defined). Trivial to fix; the values exist in pinned JSON.
- Step 3 should still **not be described as closed** until F-1 to F-4 are applied and the external items (Lean build with the n = 1 instance, the authors' reading of DKS, Windows kit runs, LaTeX compile) are done. The texts themselves say "NOT closed" everywhere; nothing over-claims.

**Not reviewed here:** the inherited mathematics of Lemma S(A) steps 1-7, Lemma U's sweeps, Lemma O/R and zone III beyond their interfaces (reviewed in `t8_review`, `fix11`, `fix15`); `L_list` and the E0/E2 blocks of `rhs_iv.py` line by line (covered by `fix5/v4_review`, `fix11` part A); step 4's P1/P2 analysis; Lean was not built; DKS was not re-read; the chains and kits were not rerun.

## 1. Q1: does the assembled argument prove the theorem?

### 1.1 Bootstrap members: where each is improved

| member (`T/bootstrap_v39.tex` Def. BS) | hypothesis level | improved to | by | my check |
|---|---|---|---|---|
| GS | sup\|ku\| ≤ ½, sup\|u′\| ≤ 1 on \|x\| ≤ 60 | 5.3903e-4, 0.040936 | GS-2, stage 24 | text and arithmetic of the [8, 10.25] step re-derived [hand]: 1e-3 + 1.5·ρ(10.25)^(-1/2)·√(Λ/ε) = 1.414e-3 |
| B1(a,b) | \|α\| ≤ barrier; Λ ≤ L + 2C′ (ℓ + 2C′ in the tail) | Y = Λ − Z ≤ L, i.e. again Λ ≤ L + 2C′ | first touch, `flux2_lemma_v39.tex` Prop. abs | **no gap: F-1** |
| B1(c) boxes | 20 box inequalities | worst ratio 0.7892 | `closure.box_assert` at Ly_hi of every cell and at the anchor | all 18 + D-2a ratios < 1 in both rate tables [cert]; Bd = 0.78913 recomputed |
| BF (F1) | window integrals ≤ 2(C0_k + C1_k T) | without the factor 2 | Lemma S(A), Lemma D, G1 | the flux constants equal (A)'s exports (116.359231…, 0.260599…) [read] |
| BF (F2) | \|(a,b)\| ≤ 1e-6 | 2.14039e-7 (×4.67) | `mod_sups` | equal to the H7 stage's ā [cert] |
| B2′ | S0V 0.0153, S1V 0.6, \|f\| ≤ 0.11, \|yf_y\| ≤ 0.52, E^A ≤ 2Ê | ×2.116, ×5.745, ×5.450, ×4.407; E^A ≤ Ê | (A) chain v37 at 4·N2_ref | N̄² ≤ 4·N2_ref confirmed in `U_cert_v36_pass2.json` [read] |
| B3U | E^ρ ≤ 2B̂; \|f\| ≤ 0.15, \|yf_y\| ≤ 0.52; \|f − U₁\| ≤ 0.16 | E^ρ ≤ B̂; M1 ×2.353 / ×2.770, M0 ×4.893 / ×2.225 | Lemma U, window Sobolev | read from the certificate |
| B4 | Lemma O a-priori bounds | H6 min ×1.860 | zone III / U-fine | read; H7 ratio 1.6000014e-7 = q/50⁴ re-derived [hand]; t* margin 2.3019 |
| G | slope ≤ ½ on {14 ≤ y, r ≤ 700λ₁} | 0.11, 0.15, 0.10032, 0.2276 | Lemma C:lem:G | read |

Thinnest margins of the whole set: ramp M1 ×1.18446, (K6) level 0 ×1.18842, zone III \|rf_r\| ×1.23889, (H1) \|α\| ×1.27, c_δU ×1.366, A-1 budget ×1.65694, H6 ×1.860. All confirmed (§2).

### 1.2 Hypotheses: which step-4 stage certifies each

| hypothesis | certified by | checked here |
|---|---|---|
| (H1) decomposition, \|α\| ≤ 1.5e-4, √Λ ≤ 5e-8 | closed v3.2 manifest (H1 stage `p2/t9h/h1.py`), verified by the v38 driver with the one B8 exception | **Λ of step 4 = Λ of fin13v4c**: `h1.py:4,22,103` uses ε = 1e-5, ε_I = 3e-5, μ = 0.08, ζ built from S₇ on [10.25, 16.25] and [57, 60], and W ≤ (1+x²)^-2; the blended W of the verifier satisfies this bound because φ′ ≥ 4x/(1+x²) for x ≥ 6.4 [hand]. `c4fun.py`, `ibp_boxes.py`, `sob_constants_P2.json` are byte-identical in `step3nl` and `step3nl_v4c` [cert, sha256]. The production run starts from Ly1 = 2.5e-15 = (5e-8)² [read]. |
| (H2) for f − U₁ on [14, 60] | stages 10, 11, 15 | bounds below (1e-8, 1e-6, 1e-3): ×27.34, ×562.8, ×600.9 [cert] |
| (H2f) on [60, 700] | stages 0-3, 16 | ×3.188, ×4.983, ×4.106 [cert] |
| (H3) windows, K₄ | stages 12, 13 | ×5.70, ×79.3, ×26.5, ×7.80 [cert]; variable is f − U with U = U₁ at τ₁ (`win_arb.py` header) |
| (H4) graph on [48, 700] | stage 12 asserts cos φ > 0 on every cell of [14, 60], [50, 380.25], [60, 700] (`win_arb.py` item (e)) | so the range is **certified**, not only "wording"; the texts may say so |
| (H5) cone, far | stages 5-8, 17, 18 | H7 and far H4 (×15.58) read |
| (H6) | Theorem C-1 | §4 |
| (GS) at τ₁ | stage 23 | read |

`scratch/verify_manifest.py`: M38 has sha256 be909bbf…, **992 pin entries, 0 mismatches** on disk now [cert]. (After path normalisation there are 646 distinct files, not "629 paths + 55 scripts = 684": 38 files are pinned under two spellings. Records only, R-1.)

### 1.3 Continuation and conclusion

- `continuation_v39.tex` Prop. C:prop:global is correct [hand]: if τ_max < ∞ then λ ≥ λ* > 0, T_max < ∞, and ∫κ²ds is bounded (tip: L² through B1; everything else pointwise with length decreasing), which contradicts the criterion. The case "T_max = ∞ and τ_mod < ∞" is excluded by t′ = λ⁴ with λ bounded.
- Finite time: T − t₁ ≤ λ₁⁴/(1 − 4ā). Blow-up of ∫_tip κ²ds: the deviation is ≤ ρ(16.25)^(-1/2)·1.68e-5 ≈ 2.9e-3 in L²(x ≤ 16.25) [hand], small against ∫κ²_{S₁}; fine.
- The whole evolution lies inside the Lemma O horizon: (T − t₁)/(50λ₁)⁴ ≤ 1.6000014e-7 against t* ≥ 3.683e-7 [cert]. So B4 holds for all τ without a time restriction.
- **(P3′) is used only for smoothness away from p\*** (part of (W5)). The finite-time singularity and the single singular point do not need it.

## 2. Q2: numbers

- `scratch/recompute.py` → `scratch/recompute.json`: **43 checks, 0 fail** [cert] (the JSON keeps 38 names, because five check names repeat for the two production points). It recomputes from the pinned JSONs: the tail floor (F0 = μ(FLUX0 + ε_I FLUX2) = 1.312192e-21, √L_F = √(2F0/ρ)), the whole-tail α-condition (×1554.456, ×1512.164), cell consistency of both marches (10 and 9 steps), L* ≤ Ly_hi and a* ≤ α_up of the last march cell, the box ratios, flux (C0)-(C4) (combined ×4.80453, C2 ×5.34e7, C3 ×1547.6, C4 ×1.85e5), ramp M1 (my value 1.18446424 against the chain's 1.18446424), (K6) (×1.18842, ×17.40, ×56.23), every hand-off threshold of §1.2, H7, zone III.
- Ramp M1 traced through the code [read]: `N/t2_A/rigor/assembly_rep.py:258-262`: 0.4329 (= 28·sI2 rounded up, `_up4` uses ROUND_CEILING and asserts) × (28 + 2⁻¹⁹)/28 + p·U(28) + sup\|Dd\|.
- `closure.py` march and tail re-read (`:100-205`): barrier padding, α fixed point, anchor certificate with floors off, whole-tail α-condition (the STRIP_BETA term is correctly taken at L_F). No direction error.
- **E_I lemma text against the code** (the line-by-line check that `ei_lemma_v4c_v39.tex` §open item 6 leaves to this review): I compared I₀, I₂, Σ_Z, R₅, R_Q, C_0W, the derived sups and the modulation forms AB, BT, E_At of Lemma EI:lem with `V/verifier/rhs_iv.py:118-200` term by term. **They agree exactly** [read]. Not compared: L₀/L₂ (`L_list`), the E0/E2 blocks.
- Printed numbers: every rounding I looked at goes outward. One inconsistency: n^III level 0 is printed 0.071446 in `step3_v39_main.tex:69` and 0.071444 in `lemmaU_v39.tex:387, 503`; the certified value is 0.0714437…, so both are valid bounds (text, F-6).
- The limits of the text set's own number check stand (value matching only, except 56 keyed rows). It did not detect F-2.

## 3. Q3: earlier closing lists

Master list: `fix15_review` §5 (which subsumes `fix11` §5 and `fix9` §7), plus the G0/E0/E2/v4 audits.

| item | source | status in v39 | evidence |
|---|---|---|---|
| G0 KILL: 18 terms absent, 17 constants low | AUDIT_RAMP, G0_INDEP | **closed** | fin13v4 (V4_REVIEW: code = spec); here: text = code for I₀, I₂, Σ_Z, R₅ |
| E0 F-B, E2 fix, R-2, R-3 | E0_AUDIT, E2_AUDIT, V4_REVIEW | **closed** | fin13v4b, fix11 part A |
| R3 forcing constant | paper stage R3 | **closed** | chain v36/v37, c_δU 1.46343e-8 |
| Flux lemma rigor gap, (I3′) | fix3 FLUX_REPORT | **closed** | flux v2, stage 21; (C0)-(C4) recomputed here |
| D-3 ADD_STRIP_CUBIC | fix11 | **closed** | term H8, fin13v4c |
| D-2a box assert; A-1 budget; A-2 ratios | fix11 | **closed** | asserted; ×1.65694; ratios |
| A-3, A-4, A-6 (E0 lemma text) | fix11, fix15 | **closed** | `E0_LEMMA_v4c_v39.tex` |
| U-1 written; (P2′) | fix9 item 1 | **closed** [hand] | `U1_regularity_v39.tex` |
| D-1, D-5, C-5 (X_S = 60, GS°, B1 form) | fix15 | **closed** | U1 Def. charts, Prop. (P2′) |
| C-2 / (GS) lemma with script | fix11 | **closed** | stages 23, 24 |
| C-1 interface; C1-1 (P2-a) written | fix11, fix15 | **closed** [hand] | `C1_identification_v39.tex` §P2a; read here, sound |
| C1-2, C1-3, C-4, GS-1, GS-3, GS-5 | fix15 | **closed** | texts |
| G1 written; (A)/Lemma D text for k = 2 | fix9 item 2 | **closed** | stage 25 rerun on v37 |
| T1-T4 flux text; B-v37-2 | fix9, fix15 | **closed** | `flux2_lemma_v39.tex` |
| D-2 one bootstrap set, one closure proof | fix15 | **closed with a defect** | `bootstrap_v39.tex`; **F-1** |
| D-3 texts on v4c/v37; D-4 density sentence; C-6 list of IBPs | fix15, fix11 | **closed** | `ei_lemma_v4c_v39.tex` |
| B5, B7 wording | fix11 | **closed** | lemmaS_A, lemmaU texts |
| Re-pin, flux outputs pinned, interface certs pinned | fix15 item 10 | **closed** | M38; 992 pins verified here |
| Step-4 wording ((H6) by C-1, GS in the list) | fix15 item 10 | **open, listed** | handoff file + patch note; `step4_proof.tex` not edited (needs the authors' go) |
| Kits rebuilt and run | item 11 | **WSL done; Windows runs open, listed** | zip hashes verified here |
| A-7 W monotone on [31.25, 34.25] | fix11 A-5, fix15 | listed open; **I close it** | ψ′ > 0 for x > 0 [hand]; certified key `W.psi1_min_lower_bound_on_[6,80]` = 0.1303 in `sob_constants_P2.json` |
| GS-4 boxes beyond the tail anchor | fix15 | listed open; **I close the box part** | L* + 2C′ ≤ Ly_hi and a* ≤ α_up of the last march cell [cert], and that cell's boxes are asserted for all shapes; no monotonicity statement is needed. The scaling step of the rate (`tail_cert`) stays [read] |
| K_ρ, K_W Sobolev constants | fix11 | listed open; **spot-audited: form valid, 5 values reproduced** | K = 1/ρ(right end): 992.4, 1.212e12; Φ₁ = φ′(34.25) = 1.083; Ψ₁ = φ′(60) = 1.305; c₁ = 0.5804 [hand] |
| DS0-DS3 level-2 commutator constants | fix11 | **open, listed** | not audited by anyone |
| D-1a pins (θ-certificate) | fix11 | **open, listed** | |
| [read] inputs of E0/E2 (Grams, pair_f3, bloop, lyap_consts2, …) | fix11 | **open, listed** | `ei_lemma` §open 2 |
| D-6a β float check with true functions | fix11 | open, optional, **not listed by name** | |
| A-1r | fix15 | exception record (stage 28); reports uncorrected; **listed** (R1) | |
| A-2r, A-3r, B-v37-1 | fix15 | **open, listed** (R2, R3) | |
| B1 (base_ibp args), B2, B3 (ramp_U hash), B4 (asserts) | fix11 | **open, listed** (R6); B1 now pinned as inputs | min margin 1.18446 > 1 confirmed here |
| C-7 per-cell arb assert in `weights_t2.py` | fix11 | **open, listed** (R6) | ≤ 1 ulp against r₂ = 1e-3 |
| GS-2, GS-6 | fix15 | **open, listed** | negligible |
| A-4 (float cell hulls in strip_cubic), C-8 (count in PROOFS_REPORT), `closure.py:130` pad formed in double | fix15, fix11, fix9 §8 | open, **A-4 and C-8 not listed in RECORDS_ERRATA_v39** | immaterial |
| M6 (no script asserts the 4·N2_ref rows with the domination) | REVIEW_v33 | **open, listed** | N̄² ≤ 4·N2_ref true [read] |
| Lean build, n = 1 instance | fix11 C-3 | **open, listed (author)** | |
| DKS re-read | fix11 C-4 | **open, listed (author)** | see §4 |
| (P3′), Lemma O arc by arc, single-point statement | REREVIEW_v33 | **open, listed (author)** | (P3′) not needed for blow-up |
| LaTeX compile | item 13 | **open, listed**; will fail until F-2 is fixed | |
| Remark C1:rem:C1 (base only continuous in τ at slab ends) | flagged for this review | **closed here** | Lipschitz in τ suffices: step 4 integrates a.e. inequalities, and (P2-a) uses only continuity and C¹ per slab [hand] |

## 4. Q4: interfaces

- **C-1** (`C1_identification_v39.tex`): the chain Lean flow → step-4 gauge → P1 (uniqueness by Grönwall in C⁰H⁵∩C¹L²) → P2 through (P2-a) → t₁ is correct and not circular [hand]. Cases A/B of Prop. P1 and the contradiction in Prop. P2 are right. (P2-a) (a)-(f) is a standard open/closed argument; its inputs (Steps 2-6 of step 4's P2 proof, the switch lemma) are labelled [cited] with step 4's own status [D]. Claims are stated with the right strength: [hand] for the argument, [cert] only for sup\|κ̃e\| ≤ 1.7429e-3. One inherited point I did not pursue: (P2-a) states e ∈ W^{7,∞} in x while the base is described as C⁴ in x; this is step 4's Step 1 and belongs to step 4.
- **(GS)**: GS-1 (tube lemma, certified on \|x\| ≤ 60.25), GS-2 (cellwise, necessary because \|u\| reaches about 1.1 near x = 60), GS-3/GS-4 (open/closed with the endpoint) are correct [hand/read]. Labels are right.
- **U-1 / (P2′)**: correct [hand]. (P2′)(c) is the standard "a maximal ODE solution leaves every compact subset of the domain" argument; the compact set K is built correctly. Remark for F-1: U-1 is stated under P(τ) ∈ 𝔅_mod, which contains B1(b); the regularity statements (a)-(d) need only (GS), (G) and τ′ < τ_mod.
- **Lean**: labelled everywhere as "read from source at d72b744, not built, n = 1 instance not elaborated". Correct and necessary. This is the main external item: existence, uniqueness, the symmetry lemma, gluing and C-1 all rest on (L1), (L2).
- **DKS**: labelled "(D1), a step inside the proof of Theorem 3.1, not re-read". This is more cautious than needed. The verbatim statement recorded in `writeup/lit_P1/P1_CONFIRMATION.md` §1 ("Let f be a maximal solution of the curve diffusion equation. If T < ∞, then ∫κ²ds ≥ c(T − t)^(-1/4)") gives Lemma "blow-up criterion" directly: by (L2) and gluing, the Lean-maximal flow is a maximal smooth solution in DKS's sense. The C-4 completion (immersion at T) is then not needed. Suggest citing the **statement** of Theorem 3.1 and keeping (D1) as a remark; the authors should confirm that DKS's "maximal solution" means maximal among smooth solutions (author item A-2).
- **Theorem statement**: symmetry of γ(t₁) is used (U-1 (S4.2), even class of the verifier) but is not in the hypothesis list of `step3_v39_main.tex` Theorem M:thm:W; it appears only as "even" in the Setting and in `step4_handoff_v39.tex` item (3). F-3.

## 5. Q5: new findings

| # | class | file:line | finding | fix |
|---|---|---|---|---|
| **F-1** | non-blocking, **proof text must be repaired** | `T/bootstrap_v39.tex:27` (Def. B1), `:43-44` (table), `:85` (Open); `T/lemmaU_v39.tex:417, 446` ("D6 … B1 strictly") | B1(b) is assumed as Λ ≤ L + 2C′ and "improved" to Y = Λ − Z ≤ L, which gives Λ ≤ L + 2C′ again (Z can reach 2C′). Likewise B1(a): the α-barrier is its own conclusion (in the tail the cell bound a\*λ is the barrier, with no cushion). So these members have **no uniform gap**, and the Open step's "scalar members: continuity and the gaps" does not extend them past τ\*. The underlying argument is sound: B1(a,b) are proved by first touch (`flux2_lemma_v39.tex` Prop. abs; `TAIL_LEMMA.tex`; cooperative in α, direct in Y), which needs the other members only. | Either (i) state in the Open step that B1(a,b) on [τ\*, τ\*+δ] follow from Prop. abs applied on [τ₁, τ\*+δ], where the remaining members hold by continuity and their gaps; say that U-1's regularity part needs only (GS), (G), τ′ < τ_mod; correct the table rows ("no gap; first touch"). Or (ii) define the hypothesis level as cell membership, Λ ≤ Ly_hi (gap certified by (C2), ≥ 5e7·2C′), with a cushion for α in the tail. (i) needs no computation. |
| **F-2** | non-blocking, **must fix before compile** | `T/step3_v39_main.tex:34, 36, 81`; `T/lemmaS_A_v39.tex:300, 302, 334`; `T/bootstrap_v39.tex:28`; `T/flux2_lemma_v39.tex:188`; `T/numbers_A_v39.tex:148-151` | Eight macros are used and never defined: `\HtwoV`, `\HtwofV`, `\EtauOneV`, `\GAoneV` (defined last in `fix1/b12_text/numbers_A_v33.tex`, dropped from the generated v37-v39 files) and `\BPhiCzerooneV`, `\BPhiConeoneV`, `\BPhiCzerotwoV`, `\BPhiConetwoV` (the generator emits `\BPhiC0oneV` …, names with digits, invalid in LaTeX). So **the thresholds of (H2) and (H2f) in the theorem, the datum of (K6), and the four (F1) constants of BF have no printed value** in the text set. Also `lemmaU_v39.tex` (`\cdUV`, `\DcTwoV`) and `main_v3_v39.tex` (`\cdUV`, `\CprimeV`, `\EzeroAoneMarV`, `\tailAlphaMarV`, `\tailAlphaMarBV`) never `\input{numbers_A_v39}`. `scratch/macro_check.py` lists them. | Add the four hand-off macros with values from `A_interface.json` `res.tight.IC` ((1e-8, 1e-6, 1e-3), (1e-6, 5e-6, 1e-5), (1.01e-6, 6e-6, 1.01e-3)); rename the BPhi macros; add the two `\input`s; add a "used ⊂ defined" test to the balance check. |
| F-3 | non-blocking (statement) | `T/step3_v39_main.tex:41` (H6) | Symmetry of γ(t₁) about the axis is used but not listed as a hypothesis of Theorem M:thm:W. | Add it to (H6), as in `step4_handoff_v39.tex` item (3). |
| F-4 | non-blocking (text) | `T/bootstrap_v39.tex:75, 87` against `T/continuation_v39.tex:64, 78` | T_max is a duration in one file ("t(τ) ↑ t₁ + T_max") and an absolute time in the other ("[t₁, T_max)", "T_max = t₁ + ∫λ⁴"). | One convention. |
| F-5 | non-blocking (labelling) | `T/U1_regularity_v39.tex` §1.3, `T/continuation_v39.tex:66` | DKS: the statement of Theorem 3.1 suffices (§4). | Cite the statement. |
| F-6 | text | `T/step3_v39_main.tex:69` | n^III level 0 printed 0.071446; `lemmaU_v39.tex` prints 0.071444; certified 0.0714437…. Both valid. | Unify. |
| F-7 | text | `T/step4_handoff_v39.tex` item (2); `T/continuation_v39.tex` gaps item 3 | (H4) on [14, 700] is certified by the coverage assert of stage 12 (`win_arb.py` (e)), not only implied. | Cite the stage. |
| R-1 | records | `T/step4_handoff_v39.tex` §1, M38 | "629 distinct paths and 55 scripts": after normalising relative/absolute spellings there are 646 distinct files; 992 pin entries is right. | Reword. |
| R-2 | records | M38 | The verifier code (`rhs_iv.py`, `unc.py`, `bb.py`, `ival.py`) and `cubic_E0fine_…json` are pinned only transitively, through the args block of the pinned rate tables (12 inputs, all matching the disk now [cert]); `tail_cert.py`, `mod_sups.py`, `integ.py` not at all. | Pin them directly at the next manifest touch. |
| R-3 | records | `T/RECORDS_ERRATA_v39.md` | fix15 A-4, fix11 C-8, fix11 D-6a and the fix9 §8 remark on `closure.py:130` are not carried in the open list. | Add four lines. |

No hidden circularity found. The loop verifier → (A) → G1 → flux → verifier is a simultaneous bootstrap with hypothesis-level inputs (0.0153, 0.6, 2Ê, doubled exports) and certified factors on return. (F2) and the ramp boxes of (A) are consequences of B1, not extra assumptions. The a-priori bound \|α\| ≤ α₁ used by (A) is the barrier value α₁(1 + 1e-12); the α-barrier is nonincreasing in both runs (max 1.4691e-4 after the first step) [cert]; immaterial.

## 6. Closing list

**Fixes (text only; no rerun):** F-1, F-2, F-3, F-4, then F-5 to F-7 and R-1 to R-3.

**External / author items (unchanged in kind):**
- A-1. Lean: `lake build` at d72b744 and a file instantiating (L1), (L2) at E = ℝ¹, M = S¹, V = ℂ. Everything about existence, uniqueness and symmetry rests on it.
- A-2. DKS: confirm the statement of Theorem 3.1 and the meaning of "maximal solution".
- A-3. Lemma O arc by arc; the single-point statement; (P3′) (only for smoothness away from p\*).
- A-4. Windows runs of `kit_step3_v4.zip` (0f899eb2…) and `kit_step4b_v4.zip` (ab47b541…); hashes verified here. The WSL runs are not independent.
- A-5. LaTeX compile after F-2.
- A-6. Step-4 text: a new version applying `STEP4_PATCH_NOTE_v39.md` (authors' go).

**Still unaudited, listed, not blocking:** DS0-DS3; the [read] inputs of the E0/E2 lemmas; the D-1a pins; the `tail_cert` degree table; M6; records R2, R3, R6 of `RECORDS_ERRATA_v39.md`.

## 7. Files in this directory
- `FINAL_REVIEW_v39.md`: this file. `REVIEW_LOG.md`: log.
- `scratch/verify_manifest.py`: re-hash of all 992 pins of M38.
- `scratch/pinned_cited.py`: which cited certificates are pinned, and by which stages.
- `scratch/recompute.py` → `scratch/recompute.json`: 43 arb checks (38 distinct names in the JSON).
- `scratch/macro_check.py`: macros used but not defined (F-2).

Nothing outside this directory was written.
