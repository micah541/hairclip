# Paper v1: outline for the final write-up (main paper + supplement)

Design document, Nov 16 2026 (Fable). Outline only; no text of the paper is fixed here except the abstract draft. Nothing existing was edited.
Abbreviations for paths: `W` = `/home/claude/s1spec/writeup`, `B` = `W/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`, `V` = `/home/claude/s1spec/step3nl_v33`.
Labels in this file: **[cert]** certified (arb, outward rounding, pinned output); **[est]** floating point or estimate; **[lit]** cited; **[hand]** hand proof.

**State of the inputs when this was written.**
- Step 3 is the v3.3 chain (v35 numbers): `A_interface` 78fe2f8e…, `U_cert_v35` 960809fd…, zone3 45e8d4dc…, step-4 `FINAL_manifest_v33.json` 8bc47e83… (21 stages), closed `FINAL_manifest.json` e01cb4da… (34 stages).
- **(P1′) is confirmed with changes** (`W/lit_P1/P1_CONFIRMATION.md`, read from the two PDFs): the $L^2$ blow-up criterion is Dziuk–Kuwert–Schätzle 2002, Theorem 3.1 (closed immersed curves in $\mathbb R^n$, curve diffusion flow, bound for every $t<T$, no area or embeddedness hypothesis; primary citation) and Chou 2003, Proposition A (plane curves, nonzero area, rate near the final time). **Still open: a citation for local existence with uniqueness** (candidates Giga–Ito 1998, Escher–Mayer–Simonett 1998; neither read), and one sentence identifying $\tau_{\max}$ with the physical maximal time through (L2). Drafters use the replacement text at the end of that file.
- No v3.3 LaTeX has been compiled. The Windows runs of the v3 kits are not in yet (WSL: kit3 v3 63 PASS / 1 FAIL provenance hash; kit4b v3 137 PASS / 0 FAIL).
- The Ogden–Warren notes `cdf_blowup_notes.pdf` are **not on disk** (CLAUDE.md lists them as an upload). Drafters cannot quote them; see decision D9.

---

## 1. Title, abstract, venue

### 1.1 Title options
1. **A finite-time singularity of curve diffusion flow from an embedded curve** (recommended: says the result, makes no priority claim, short).
2. A curve diffusion singularity from an embedded hair clip, modelled on a non-embedded shrinker (current draft title; accurate, long).
3. An embedded closed curve that becomes singular in finite time under curve diffusion flow: a computer-assisted proof.

Do not use "first" in the title. Running head: "CDF singularity from an embedded curve".

### 1.2 Abstract draft (about 150 words)
> Curve diffusion flow moves a plane curve with normal velocity $-\kappa_{ss}$. It has no maximum principle, and no soft criterion forces an embedded closed curve to become singular. We exhibit an explicit embedded closed curve $\Gamma_0$, a symmetric "hair clip" whose tangent angle is a finite sine series, and prove that its evolution becomes singular in finite time: the curvature blows up at one point, at the self-similar rate $(T-t)^{-1/4}$, and the rescaled curve stays trapped near a self-similar shrinker $S_1$ with one loop. Since $S_1$ has a transverse double point, the flow loses embeddedness before the singular time. The proof is computer-assisted in four steps: existence of $S_1$ with tight enclosures; its linear stability modulo symmetries; a nonlinear trapping theorem in rescaled variables; and a validated integration of the flow from $\Gamma_0$ into the trapping region. All numerical statements are certified in interval or ball arithmetic, and the computations are reproducible from a public archive.

Drafter notes: "one point" must match the wording chosen in D13; drop "public archive" if D8 is not settled.

### 1.3 Venue options
| option | fit | page implication |
|---|---|---|
| **A. Long-paper analysis/PDE journal** (candidates: Annals of PDE, Comm. Pure Appl. Math., Arch. Ration. Mech. Anal., Mem. AMS if merged) | singularity formation with rigorous numerics is in scope; referees expect a full proof architecture in the paper | main 45–55 pp + supplement 150–180 pp as electronic material |
| B. Computer-assisted-proof-friendly journal (candidates: SIAM J. Appl. Dyn. Syst., Found. Comput. Math., J. Differential Equations) | referees used to interval arithmetic, code review, supplements | same split; emphasis shifts to sections 7–8 |
| C. Geometric analysis journal (candidates: J. Differential Geom., Calc. Var. PDE, J. reine angew. Math.) | best audience for the statement; least used to 150 pp of certificates | main must shrink to about 35 pp, step 4 to 5 pp |
| D. Short announcement + long paper | fastest priority | 6–8 pp note; costs a second document |

**Recommendation: A, with arXiv first** (main PDF + supplement PDF as one submission; code and data as a DOI archive, not as arXiv ancillary files, because of size: kits are 2 MB and 41 MB, the proof tree is larger). Reasons: the main theorem needs sections 5–7 at full length to be checkable in outline; option C forces cuts exactly where a referee needs detail. I have not checked any journal's current policy on supplements, page limits or AI-assisted work; the authors must do that before choosing (D1). The plan below is sized for A and degrades to C by cutting section 7 to its statement and tables.

---

## 2. Main paper: section plan (target 50 pp, amsart 11pt)

"Self-contained" is defined as: every object, hypothesis and theorem is stated in full; every argument that is not a computation and is under two pages is proved in the paper (bootstrap logic, continuation, conclusions, the reduction of the Main Theorem to Theorems A, C, D); long estimates are stated precisely with a pointer "proved in Supplement X". No internal file path appears in the main paper; paths live in the supplement (Part J).

Theorem names in the paper: **Main Theorem** (1.1), **Theorem A** (shrinker), **Theorem B** (spectrum), **Theorem C** (trapping), **Theorem D** (validated integration). Section numbers below are the paper's.

### §1 Introduction (4 pp)
- Purpose: problem, why embedded data are hard, the scenario, the statement, the method in one page, related work.
- Content: CDF and its basic identities; no maximum principle, loss of embeddedness and convexity; the "no soft criterion" remark (turning number and area conserved, circle admissible), proved in three lines here; the clip scenario (slot + chamber = $S_1$ with its crossing resolved); Main Theorem; "what is not claimed"; related work with the Wheeler paragraph; organisation.
- Sources: `B/t7_assembly/main_v3.tex` §1, §8 (related work; reuse its Wheeler paragraph almost verbatim, minus the status clause), `W/main/main.bib`.
- Numbers: $t_1\in[4.5978574775,4.5978586144]$, $\lambda_1\in[6.49\cdot10^{-5},6.75\cdot10^{-5}]$, $T-t_1\le2.1\cdot10^{-17}$, $\bar a\le2.19103\cdot10^{-7}$ (sources: `S4/p2/out/FINAL_manifest.json` stage t1clock; `V/verifier/modsups_fin12f4/…`; step4_proof §10 for $\lambda_1$).
- Figures: **F1** (three panels: $\Gamma_0$, the bend, $S_1$ near the tip; exists as `W/main/fig_curves.pdf`, script `make_fig.py`); **F2** new: 4–6 floating-point snapshots of the evolution showing chamber → loop → crossing, captioned "illustration, not a certificate" (D14).
- Main Theorem wording changes against main_v3: remove every status clause and version tag; item (ii) per D13; item (iv) keeps peak, rate, floor; state the certified-set caveat once ("for every $S_1\in\mathcal S$").

### §2 The shrinker, the initial curve, and the four theorems (5 pp)
- Purpose: all data and all four component statements in one place, so the reader can see the logical skeleton before any method.
- Content: shrinker equation and the profile system (derive it; do not cite the notes); reversing symmetry; Theorem A; Definition of the certified set $\mathcal S$; Table T1; Definition of $\Gamma_0$ (exact data, $b_2$ by Krawczyk) and Proposition (closed, embedded); Theorem B; rescaled frame in five lines, then Theorem C with hypotheses (H1)–(H7) written out in full (D6 renumbering) and conclusions (W1)–(W5); Theorem D; proof of Main Theorem from A, C, D (half a page); Remark "what the proof uses beyond computation" (literature facts (L1)–(L3); the role of Theorem B, D4).
- Sources: main_v3.tex §2, Theorems `thm:S1`, `thm:step2`, `thm:W`, `thm:4`; `B/fix1/b12_text/step3_v33_main.tex` §2 (the authoritative (H) list) and `continuation_v33.tex` §5.4 ((W1)–(W5)); `/home/claude/s1_cap/s1_proof.tex` §1; `W/step2_proof.tex` §1; `W/step4_proof/step4_proof.tex` §1.3, §2.
- Tables: **T1** $S_1$ enclosures (main_v3 Table 1; sources `s1_tight/bracket2.json`, `s1_tight/ext_tip/README.md`).
- Figure: **F3** dependency diagram (tikz in main_v3, simplified: A→B, A→C, A→D, D→C, {A,C,D}→Main).

### §3 Strategy and the rescaled formulation (4 pp)
- Purpose: the analytic set-up of steps 3 and 4 and the map of regions; the one place where notation is fixed.
- Content: frame $\gamma=c\,e_1+\lambda\Gamma$, $d\tau=\lambda^{-4}dt$, modulation $(a,b)$; linearisation $L$ and weight; slow mode and decomposition $u=\alpha\hat\psi_3+\alpha^2\hat h_2+w$; energies $E$, $E_I$, $\Lambda$; far-field variable $f$ (slope over the ray), $U_1$, $z_k=(y\partial_y)^kv$; region table; why three different mechanisms are needed (tip: spectral gap; cone: transport outruns diffusion; far: curve does not move before $T$); how steps hand off.
- Sources: main_v3.tex §3; `W/step3_proof_v32/step3_proof_v32.tex` §1 (frame, chart, weights, regions); step3_v33_main.tex §1 (v3.3 regions).
- Tables: **T2** regions (main_v3 Table 2 updated to: tip, intermediate, shells $[14,60]$ + segments $(60,350]$ + annulus $[350,700]$, zones I/II, zone III $[253,379.5]\lambda_1$, cone family, far family).
- Figure: **F4** new schematic: one log-radius axis from the tip to the far curve, regions shaded, each labelled with its controlling lemma and its quantity; a second row showing how the regions move as $\lambda\to0$.

### §4 Step 1: the shrinker (4 pp)
- Purpose: proof sketch of Theorem A with every constant that matters.
- Content: conical ends by Lyapunov–Perron with explicit constants; the sharp end state by the exact decaying solution (contour integral, validated quadrature); shooting function $M(\varepsilon)$, two certified signs, bisection; geometry (double point, embedding of half-branches, turning); exports used downstream and the containment check (step 3 and step 4 use the same certified set); trust base of this step.
- Sources: `/home/claude/s1_cap/s1_proof.tex` (§§3–5), `/home/claude/s1_tight/s1_tight.tex`, READMEs of `s1_tight/{ext813,ext_tip,consistency}`; main_v3.tex §4.
- Numbers: $M(1.45\cdot10^{-4})$, $M(1.62\cdot10^{-4})$ enclosures; $|d|\le4.37\cdot10^{-13}$; 51,246 boxes, widths $3.2\cdot10^{-12}$, $5.7\cdot10^{-10}$; containment 1,025 + 501 rows, smallest relative margin $1.3\cdot10^{-9}$ (source `s1_tight/consistency/check_outA5.json`; see risk R7).
- Figure: none beyond F1(c); optional small plot of $M(\varepsilon)$ [est].

### §5 Step 2: linear stability modulo symmetries (4 pp)
- Purpose: Theorem B and its method; state precisely what later steps take from it (D4).
- Content: numerical-range certificates (SDP finds, Arb verifies) → trapezoid; admissible plane at the conical end; Evans functions on the contour, winding number 2 per class, analyticity cover; the multiplicity remark [lit] and that its hypotheses are not checked in the weighted half-line setting; numerical spectrum [est].
- Sources: `W/step2_proof.tex` (6 pp, all sections), `/home/claude/s1spec/README_step2.md`.
- Numbers: $\Lambda_0=1.790819$, $\Lambda_1=1.768282$, contour vertices, 202 boxes, $\eta\le0.1163$, contraction $\le0.272$, 5,058 segments, 174 boxes (sources: the step-2 logs `wind_even.log`, `wind_odd.log`, `apr/`; the drafter must record the exact file per number in the ledger).
- Figure: **F5** spectrum region and contour (`W/spectrum_region.png` exists; regenerate as vector).

### §6 Step 3: nonlinear trapping (14 pp)
The heart of the paper. Subsections:
- **6.1 Statement and plan (1 p).** Theorem C restated by reference; bootstrap set $\mathrm{BS}$ = (B1) barrier, (B2) core box and shell energies, (B3) uniformly local energies and far box, (B4) cylinder bounds, (B5) graph property (the "G" of continuation_v33). Source: step3_v33_main.tex §7.
- **6.2 Tip and intermediate zone (2.5 pp).** Constrained coercivity (Theorem: $\mathrm{Re}\langle Lv,v\rangle_\rho\le-\frac3{64}\|v\|^2$ and the $\theta$-form), energies and Lyapunov constants, forcing, modulation by integration by parts, the comparison system and its certified solution (barrier: peak, rate, floor). Sources: step3_proof_v32.tex §3.1–3.5; `W/step3_linear_certificate.tex`; `V/V33_VERIFIER.md`. Numbers: $14.8e^{-0.0213\tau}$; peak $5.308\cdot10^{-8}$, rate $2.693\cdot10^{-3}$, floor $9.872\cdot10^{-10}$, $\alpha_F\le6.867\cdot10^{-7}$, $\bar\beta\le3.33811\cdot10^{-9}$ (source `V/verifier/rate_table_fin12f4_a0.00015.json`, `tail_cert_log.txt`).
- **6.3 Shells: Lemma S (2.5 pp).** Variable $f-U$ on $[14,700]$; divergence-form energy identities (statement; symbolic check); the level-$j$ inequality; what the ramp, the interface at $y=28$ and the outer annulus cost; chart identification (verifier variable vs ray graph) as a lemma; statement of the theorem with $\hat E$ and the core box. Sources: `B/fix1/b12_text/lemmaS_A_v33.tex` §§2–5, `chart_lemma_body.tex`, `numbers_A_v33.tex` (the number macros with their JSON keys). Numbers: $\hat E=(1.2018\cdot10^{-6},1.0441\cdot10^{-4},0.056694)$; margins core $M_1$ ×4.409, $S_0^v$ ×2.116, $\hat G_{22}$ ×1.917, fluxes ×3.365 / ×5.955, ramp $M_1$ ×1.184 (source `B/fix2/n12/t2_A_v35/out/A_interface.json`, `asm_design_N2ref4.json`).
- **6.4 Far field: Lemma U (3 pp).** The idea in plain terms first: energies centred at every physical radius, a weight that never vanishes, one frozen profile, so there are no window edges to estimate; transport identity (Lemma T); launch data; the data term via the shells; Proposition U (comparison argument); window Sobolev → far box; Lemma N (outer annulus data derived, not assumed); the fine version for the strip. Sources: `B/fix1/lemmaU_rev/lemmaU_v34.tex` §§1–8, `B/t1_statements/STATEMENTS_r2.tex`. Numbers: $\hat B=(4.8735\cdot10^{-6},5.5933\cdot10^{-4},0.36675)$, $\hat B^{\rm fine}$, margins ×2.353 / ×4.892 / ×2.770 / ×2.225, loop gain $\le0.03904$, rate ×28.2, $M_{1z}=0.534$, $c_{\delta U}\le1.7187\cdot10^{-8}$ (source `B/t4_cert/out/U_cert_v35.json`).
- **6.5 Cylinders: Lemma O, Lemma R, zone III, the closure of the cylinder bounds (2 pp).** Graph persistence on cylinders up to a time ratio $t^*$; fixed-$L$ cone family and far family; time-averaged fourth-derivative bound; zone III boxes; recovery of the a-priori bounds (the "H6" theorem; rename, D6). Sources: `W/lemmaO_v32/lemmaO_v32.tex`, `lemmaR.tex`; step3_v33_main.tex §6; `B/t5_zone3/Z5_REPORT.md`. Numbers: $(M,K)$ = (0.12, 0.18), (0.15, 0.3); $t^*$ ×2.30, ×15.5; zone III box (0.10032, 0.41974, 0.014291), thinnest ×1.238; recovery minimum ×1.860 (source `B/t5_zone3/out/zone3_interface.json`, `W/lemmaO_v32/lemmaO_v32_consts.json`, `lemmaR_consts.json`).
- **6.6 The bootstrap closes (1.5 pp).** The order table (who consumes what; no step uses its own conclusion), the closure conditions (K1)–(K6) as one table with margins, the open–closed argument including the graph property. Full proof in the paper. Sources: step3_v33_main.tex §7; lemmaU_v34.tex §7.3–7.5; continuation_v33.tex §2. Table **T4**.
- **6.7 Global existence in rescaled time and the conclusions (1.5 pp).** Region-by-region bounds table; $\tau_{\max}=\infty$ by bounding total $\int\kappa^2ds$ and applying (L1) = DKS Theorem 3.1 (use the wording of `W/lit_P1/P1_CONFIRMATION.md`; add the sentence that the modulation extends as long as the solution exists and (B1) holds, so $t(\tau_{\max})=T_{\max}$; cite DKS (2.10) for length decrease); $\lambda\to0$, $T<\infty$, rate; single singular point; smooth convergence away from it. Full proof in the paper. Source: continuation_v33.tex §§3–5 (apply the re-review's wording fix: use $60\lambda$ not $56\lambda$ when converting distance to $y$).
- Figures: **F6** bootstrap dependency graph (new; from the order table). F4 is reused.

### §7 Step 4: validated integration from $\Gamma_0$ to $t_1$ (9 pp)
- Purpose: Theorem D with enough of the method that a rigorous-numerics reader can judge it.
- Subsections and sources (`W/step4_proof/step4_proof.tex` = S4doc; `S4/STEP4_RESULT.md`; `B/t7_assembly/STEP4_RESULT_v33_addendum.md` = ADD):
  - 7.1 $\Gamma_0$: closure by Krawczyk, embeddedness by arcs and boxes (S4doc §2). 0.5 p.
  - 7.2 Phase I, Fourier slabs: equation for the tangent angle, one slab, the first-exit/Galerkin-limit lemma (statement; proof in supplement), results (S4doc §3). 1.5 pp.
  - 7.3 The base solution and the switch to rescaled variables (S4doc §§4–5). 1 p.
  - 7.4 Phase II, energy shadowing: objects, the shadowing theorem, residuals on the corrected base, comparison system, coercivity certificates, far-field lemmas, nonlinear boxes (S4doc §6; statements only for the far-field lemmas). 3 pp.
  - 7.5 The hand-off on one slab: implicit-function/Krawczyk at 16 points and 32 derivative balls; shells, far datum, windows; cylinder data; the time $t_1$; the certified self-intersection at $t_1$ (S4doc §§7–8; ADD §§1, 5; `S4/p2/v33/{T6,M9,F12}_REPORT.md`, `B/fix2/n12/N12_REPORT.md`). 2.5 pp.
  - 7.6 Where step 4 is tight (S4doc §1.5). 0.5 p.
- Write the hand-off for the final hypothesis list only (D15): no v3.1/v3.2 intermediate thresholds, no "float row superseded by arb row" narrative; quote the arb rows.
- Tables: **T3** hand-off items: required / certified / margin, one row per item of (H1)–(H7). Sources: `S4/p2/out/FINAL_manifest.json` outputs for (H1), (H3) $K_4$, (H4); `S4/p2/out_v33/{h2u_v33_arb,h2f_result_arb,h3_arb,h2t_v33_arb,t9_cone_v33,h4loc_v33,t9_lemmaO_v33_arb,h4far_v33,ratio_v33}.json` for the rest. See risk R10 on the (H3) margins.
- Numbers: 3762 modes; half-width $8.137\cdot10^{-24}$, contraction $3.01\cdot10^{-12}$; 701 arcs, separation $8.89\cdot10^{-3}$, half-gap $\ge0.0248219$; 96 slabs, $t_{96}$, $k_{\max}$, 2800 modes; 181 base slabs; $Y(\tau_a)\le6.21\cdot10^{-13}$ vs $\bar Y=10^{-12}$; 26 certificates, $\min c_I\ge0.05284$, $\min\varrho_c\ge0.04711$; far boxes ×6.42 / ×4.69 / ×1.99 / ×1.69.
- Figure: **F7** new: timeline of step 4 on a log scale of maximal curvature against time: Phase I (slabs 0–96), switch, Phase II (181 slabs, $k$ from 26 to $2.4\cdot10^4$), hand-off slab; [est] illustration.

### §8 Computer assistance, trust base, reproducibility (4 pp)
- Content as specified in section 4 of this outline. Tables: **T5** inventory (step; what is computed; tool and arithmetic; wall time; where pinned); **T6** thinnest margins of the whole proof (one table, about ten rows).
- Sources: main_v3.tex §7; S4doc App. A; `B/t7_assembly/V33_INDEX.md` §4.6–4.7; kit READMEs (`/home/claude/s1spec/kit_step3_v3/README.md`, `kit_step4b_v3/README.md`); run notes (`kit3v3_run_wsl1/RUN_NOTES.md`, `kit4b_v3_run/RUN_NOTES.md`, `kit_step4a/user_runs/`).

### §9 Discussion (1.5 pp)
- Not claimed and open: convergence of the rescaled flow to $S_1$; uniqueness of $S_1$; a certified first crossing time; an open set of data (D10); other clips; continuation past $T$.
- Source: main_v3.tex remark "What is not claimed", `step3_proof_v32.tex` §5.5.

### Back matter (2.5 pp)
Acknowledgements; statement on the use of AI (D3; draft in main_v3.tex §10, to be rewritten by the authors); data and code availability; references (30 keys of `W/main/main.bib`, keys unchanged).

**Page sum:** 4 + 5 + 4 + 4 + 4 + 14 + 9 + 4 + 1.5 + 2.5 = 52.

---

## 3. Supplement: parts (target 160–180 pp)

Parts are lettered A–K, not "S1…Sn", because $S_1$ is the shrinker and "S1" is already a document label in the archive. Each part opens with: what it proves, which statement of the main paper it supports, and its trust base.

| part | content | pp | sources |
|---|---|---|---|
| **A** | Conventions: notation table, labels, rounding rule, the number ledger (how every quoted number maps to a key of a pinned file) | 5 | section 5 of this outline; `B/fix1/lemmaU_rev/constants_needed_v34.json` and `numbers_A_v33.tex` as models |
| **B** | Step 1 in full: conical ends, sharp end, shooting, geometry, exports, containment check | 14 | `/home/claude/s1_cap/s1_proof.tex`, `/home/claude/s1_tight/s1_tight.tex`, READMEs |
| **C** | Step 2 in full | 8 | `W/step2_proof.tex` |
| **D** | Step 3, tip and intermediate zone: coercivity certificates, slow mode and slow-manifold correction, energies, forcing, modulation, the comparison system and its verifier, the tail lemma | 22 | `W/step3_linear_certificate.tex`, `step3_slowmode_cert.tex`, `step3_slowmanifold2.tex`, `step3_proof_v32.tex` §§1–3.5, `/home/claude/s1spec/step3nl/CERT_PLAN.md`, `V/V33_VERIFIER.md` |
| **E** | Lemma S (shells) in full, with the chart lemma and the profile difference $d$ | 18 | `B/fix1/b12_text/lemmaS_A_v33.tex`, `chart_lemma_body.tex`, `numbers_A_v33.tex`; `W/lemmaS/lemmaS.tex` Prop. 1 (energy identities) and `lemmaS/rigor/RIGOR_PASS.md` |
| **F** | Lemma U in full; zone III | 24 | `B/fix1/lemmaU_rev/lemmaU_v34.tex`, `B/t1_statements/STATEMENTS_r2.tex`, `B/t5_zone3/Z5_REPORT.md`, `B/t4_cert/{T4,U4}_REPORT.md` |
| **G** | Lemma O, Lemma R, recovery of the cylinder bounds | 12 | `W/lemmaO_v32/lemmaO_v32.tex`, `lemmaR.tex`, `W/lemmaO/lemmaO_v3.tex` |
| **H** | Step 3 assembly: full closure proposition with every (K) condition tied to a key; table "which region supplies which box" | 6 | `step3_v33_main.tex` §§7–9, `lemmaU_v34.tex` §§7, 9–11, `continuation_v33.tex` §6 |
| **I** | Step 4 in full, written for the final hand-off list | 45 | S4doc (39 pp) merged with ADD and the `S4/p2/v33/*_REPORT.md` files; `S4/STEP4_DESIGN.md`, `S4/LEMMA_F.md` for proofs |
| **J** | Computation records | 15 | see 3.1 |
| **K** | Provenance: version map, review history, record exceptions | 4 | see section 4.3 |

Parts D–H must be one consolidated text of step 3 with no "as in v3.2 §x" pointers (D7). This is rewriting for exposition only: no statement, constant or script changes; new files.

### 3.1 Part J: how the computation records are presented
1. **Environment**: Python 3.13, python-flint 0.9.0 (Arb), sympy 1.14.0, CAPD 6.1.0; `PYTHONHASHSEED=0`; one process; the per-stage environment recorded in the manifests.
2. **Two manifests, printed as tables**: stage index, stage name, one-line content, wall time; then the SHA-256 of each manifest and of the headline outputs (the table of `V33_INDEX.md` §4.6–4.7). Full per-file hashes stay in the archive, not in the PDF.
3. **Step-3 chain**: dependency order verifier → shells → Lemma U ↔ zone III (state plainly that this pair is a fixed point computed in passes and that the two files cannot pin each other's hash) → text checks. Commands from `V33_INDEX.md` §4.6.
4. **What a reader can rerun, in three tiers**:
   - Tier 1, minutes, no compiler: kit "step 3" (`kit_step3_v3.zip`, e40b444a…, about 35–40 min): verifier, modulation, shells, Lemma U and zone III, Lemma O, text checks.
   - Tier 2, hours, no compiler: kit A (`kit_step4a_v3.zip`: $\Gamma_0$ and Phase I, about 6 h one core) and kit B (`kit_step4b_v3.zip`, 646f2552…, base, Phase II, hand-off, 39 stages, about 3–6.5 h).
   - Tier 3, needs CAPD and a C++ toolchain: steps 1 and 2 (`s1_cap/run_all.sh`, `s1_tight/run_tight.sh`, `run_step2.sh`). No kit exists; say so (R12).
5. **What "reproduced" means**: byte-identical apart from timing and provenance fields on the same platform; on another platform, overlapping enclosures and identical PASS tables; known platform effects (step sizes in Phase I differ, flint-version last digits). Quote the actual check summaries verbatim, with machine and date.
6. **Record exceptions** (pointer to Part K).

---

## 4. Presenting the computer assistance, the trust base and the history

### 4.1 Three classes, stated once in §8 and used everywhere
- **Certified**: an inequality between exact objects (dyadic numbers, stored balls, finite series) evaluated in ball or interval arithmetic with outward rounding by a pinned script. Floating point only proposes objects (SDP certificates, bases, step sizes, the base solution); the statement is about the exact object. No sampling.
- **Proved by hand**: every reduction from the PDE to those inequalities (energy identities, comparison arguments, the open–closed argument, the Galerkin-limit lemma). Symbolic identities were checked by computer algebra offline; say that this is a check, not part of the proof.
- **Cited**: (L1) local well-posedness with uniqueness (reference still to be confirmed) and the $L^2$ blow-up criterion for closed immersed curves (DKS 2002 Theorem 3.1; Chou 2003 Proposition A; confirmed); (L2) regularity of the modulation; (L3) interior regularity for graphs; in step 2 only, order of an Evans zero = algebraic multiplicity.
The Main Theorem carries **no unproved hypothesis**; it depends on the three classes above. Say exactly that, and do not write "unconditional" without the list next to it.

### 4.2 Trust base, as a short explicit list in §8
Arb through python-flint; CAPD (steps 1, 2); IEEE-754 double arithmetic with `nextafter` widening in the step-3 verifier (not Arb: say so); Python and sympy as carriers; the authors' and agents' scripts (about 39 + 34 stage scripts in step 4 and the step-3 chains), which are the largest item and are not formally verified. State plainly: the code and most of the proof text were written by AI agents under the authors' direction; reviews were by separate AI runs; give the list of what human authors checked by hand (D11). A referee will ask; the answer must be in the paper.

### 4.3 The v3.1 → v3.3 history
**The erratum does not belong in the paper.** v3.1 and v3.2 were never published; an erratum to unpublished drafts confuses a reader and adds risk of misquoted numbers (main_v3 already carries a withdrawn sentence with a TODO). The paper and Parts D–H present one proof, with no version tags.
What does belong:
- §8, one paragraph "Development and review": the far-field estimate was first written with window energies whose edge terms were dropped; independent review found this and several omitted ramp terms in the shell lemma; the final proof replaces windows by uniformly local energies. This is honest, short, and explains why Lemma U looks the way it does.
- Part K, 4 pp: (i) map from paper names to archive names and versions (the archive keeps closed v3.1/v3.2 files byte-identical, and file names such as `lemmaS_A_v33`, `U_cert_v35` carry tags); (ii) list of reviews and what each found; (iii) pointer to `B/t7_assembly/erratum/` in the archive for anyone reading the old files; (iv) record exceptions (R7).

---

## 5. Consistency rules for the drafters

### 5.1 Numbers
1. **No number without a ledger entry.** Before drafting, unit U0 builds `paper_v1/numbers/ledger.json` (name, value string as printed, kind cert/est, source file, key, sha256 of the source) and generates `numbers_paper.tex` with one macro per number. Drafters use macros only; a literal certified number in the text is an error. A script checks macro use and re-evaluates each value against the pinned source (model: `B/fix1/lemmaU_rev/check_keys_v34.py`).
2. **Rounding**: upper bounds up, lower bounds and margins down, at the printed digit. Margin = threshold / bound.
3. **Labels**: in the main paper certified is the default and unmarked; every non-certified number carries "(floating point)" or "≈" and never appears in a theorem. In the supplement keep the explicit tags [cert], [est], [hand], [lit].
4. A number may be quoted in two places only through the same macro.
5. When two archive documents disagree, the pinned JSON wins and the disagreement is reported to the coordinator, not resolved silently.

### 5.2 Notation (fix in §3 and Part A; these symbols clash in the sources)
| object | paper symbol | clash being removed |
|---|---|---|
| slow-mode amplitude | $\alpha$ | curve parameter of $\Gamma_0$ becomes $\vartheta$ |
| rescaled time | $\tau$ | $\tau_{\rm tip}$ of step 1 becomes $\sigma_{\rm tip}$-based; tangent vector is $T$ |
| scale | $\lambda$ | spectral parameter becomes $\mu$; slow eigenvalue $\mu_3$ (was $\tilde\lambda$); energy weight $\mu=0.08$ becomes $\mu_I$ |
| Gaussian-type weight | $\rho$ | centre radius of Lemma U becomes $r_c$; tail rate becomes $\varrho$ |
| shooting function | $M$ | Lemma O constants become $(M_O,K_O)$; core boxes stay $M_0,M_1$ |
| linearised operator | $L$ | length is $\mathcal L(\gamma)$; cylinder unit is $\ell$ |
| arclength of $S_1$ from tip | $x$; $y=x-2.2490776515$ | ambient coordinates are $(X_1,X_2)$ |
| tip distance of $S_1$ | $D$ | $y\partial_y$ becomes $\partial_s$, $s=\ln y$ |
| point of the frame | $c(t)e_1$ | do not use $p(t)$; $p$ is the support function; profile exponent is $p_U=1.0246$ |
| step-4 error majorant | $\bar Y$ | step-1 parallelepiped becomes $\mathcal Y$ |
| curvature | $\kappa$ (flow), $k$ ($S_1$) | keep |

### 5.3 Names and numbering
- Hand-off hypotheses (H1)–(H7): (H1) modulation and smallness, (H2) shells, (H3) far datum (archive: H2f), (H4) windows and fourth jet (archive: H3), (H5) graph (archive: H4), (H6) cylinder data (archive: H5), (H7) regularity (archive: H6). Part K carries the map. (D6)
- Lemma O's own hypotheses: (O1)–(O9), never "H". The recovery theorem (archive "H6") is "Proposition (recovery of the cylinder bounds)".
- Literature facts: (L1)–(L3) (archive: P1′–P3′). Step-4 phases: Phase I, Phase II (archive: P1, P2). Bootstrap: (B1)–(B5); derivation steps keep D1–D6 only in Part H.
- "Hypothesis (F)", "(O)", "B0", "M9", "N1", "G-A1", "fin12f4", "v35" and all review-item codes do not appear in the main paper. In the supplement they appear only in Parts J and K.
- Theorem numbering by section; equations numbered by section; labels prefixed per unit (`intro:`, `obj:`, `s3:` …) to avoid collisions when sections are merged.
- One preamble file `paper_v1/preamble.tex` owned by the coordinator; drafters add macros by request only.
- Citation keys: exactly those of `W/main/main.bib` (30 entries, 29 verified per CLAUDE.md; the drafter of §1 lists any key whose verification status is unknown). New references need a verified entry (publisher or arXiv page) before use. The paper cites DKS as `DziukKuwertSchatzle2002`; the keys `DKS2002` in `continuation_v33.tex` are local and must not be copied.

### 5.4 Style
Plain declarative sentences; no status boxes, dates or "pending" in the paper body; every theorem states its trust class in a closing sentence ("The proof is computer-assisted; see §8"); a drafter who finds a gap writes it into `paper_v1/ISSUES.md` and stops that paragraph with `\GAP{…}`, never papers over it.

---

## 6. Decisions only the authors can make

| # | decision | recommendation |
|---|---|---|
| D1 | Venue | Long-paper PDE journal after arXiv; check supplement and AI policies first. |
| D2 | Title | Option 1. |
| D3 | AI-use statement, acknowledgement wording, who carried out the referee check, and authorship statement of responsibility | Full disclosure in back matter plus the trust-base paragraph of §8; authors write it themselves. |
| D4 | Role of Theorem B. From what I read, step 3's linear input is its own coercivity certificate (Morse index and Gram matrix), not the Evans count; if that is confirmed, the Main Theorem does not depend on Theorem B or on the multiplicity citation | Have one agent trace it; if confirmed, present Theorem B as an independent result that explains the spectral picture, and say so. This shrinks the cited class. |
| D5 | Erratum and version history | Out of the paper; one paragraph in §8, Part K in the supplement. |
| D6 | Renaming (H2f → H3 etc., P1′ → L1, Lemma O hypotheses) | Rename, with the map in Part K. Cost: the paper no longer matches JSON keys one to one. |
| D7 | Consolidated step-3 text in Parts D–H instead of the layered v3.2 + v3.3 files | Consolidate; no mathematical change; verify by the number ledger and a statement-by-statement diff review. |
| D8 | Public archive: where, licence, and whether the tree is made relocatable (code has absolute paths under `/home/claude`) | DOI archive with the three kits plus a tarball of the proof tree as record; kits are the supported way to rerun. |
| D9 | The Ogden–Warren notes: cited in drafts for the shrinker equation, the soft-criterion remark and the scenario; not on disk | Absorb the two short derivations into §1–2 and do not cite unpublished notes; credit the numerical discovery of $S_1$ in the text. |
| D10 | An open-set corollary: all hand-off conditions have strict margins, so symmetric curves near $\Gamma_0$ in a strong topology should also become singular | Worth one agent-day to check whether every item is open under continuous dependence; include only if a clean proof exists. Not proved now. |
| D11 | What the authors verify by hand before submission | Minimum: §6.6–6.7 (bootstrap logic, continuation, conclusions), the statements of Theorems A–D, the uniqueness reference in (L1) and (L2)–(L3) against the literature, Lemma O read arc by arc for immersed curves. State the list in the paper. |
| D12 | Novelty sentence against Wheeler (arXiv 2608.21882) | Keep "complementary" paragraph; claim only "we are not aware of an earlier proof of a finite-time singularity of curve diffusion flow from an embedded closed curve"; recheck the arXiv listing for a new version on submission day. |
| D13 | Main Theorem (ii): the curve is immersed, so "the singular set is a single point" needs a precise form | State it in the parameter: curvature stays bounded on compact sets of parameters away from the tip parameter, and blows up there; the image point is $p_*$. Other arcs may pass near $p_*$. |
| D14 | Floating-point illustrations (snapshots, first crossing at $t\approx0.003$) | Include F2 and the crossing time, clearly marked as not certified. |
| D15 | Step 4 written for the final hand-off list only | Yes; intermediate thresholds go to Part K. |
| D16 | Wait for the independent Windows runs of both v3 kits, and fix the two checker items, before posting | Yes for the runs; the checker fixes are small and remove the one FAIL line a referee would see. |
| D17 | The record exception for `check_outA5.json` (R7) | Disclose in Part K as is; optionally build a release manifest that pins the current file. Do not rerun closed stages for cosmetics. |
| D18 | Certified-set caveat ("for every $S_1\in\mathcal S$") in the Main Theorem | Keep it once in the statement; it is the honest form, since uniqueness is not proved. |

---

## 7. Gaps, risks, and the drafting plan

### 7.1 Gaps and risks to flag in the paper (or resolve before it)
- **R1 Literature facts.** The blow-up criterion in (L1) is confirmed (DKS 2002 Theorem 3.1; Chou 2003 Proposition A). Not confirmed: the local existence **and uniqueness** statement for smooth closed immersed plane curves (Giga–Ito 1998 or Escher–Mayer–Simonett 1998, neither read; uniqueness matters because step 4's validated solution and step 3's modulated solution must be the same maximal solution), and precise references for (L2) and (L3). `main.bib`'s `GigaIto1998` lacks pages; Chou's issue number is unverified.
- **R2 No human line-by-line check of the step-3 text.** The re-review lists what it did not check: the shell proof Steps 2–5 against the scripts, verifier internals, inherited error-transfer constants, the window formulas of the hand-off code (`B/t8_review/REREVIEW_v33.md` §8).
- **R3 Inherited, not re-audited items of the shell lemma** (`lemmaS_A_v33.tex` G-A5): $\cos\theta$ taken constant along the ray in two coefficients; the forcing constant; the flux lemma. The consolidation unit must either audit them or list them in the paper.
- **R4 Evans multiplicity** hypotheses unchecked in the weighted half-line setting (matters only if D4 fails).
- **R5 Finite but uncertified constants** $c_g$, $C_{S_1}$, $\ell_{\rm tip}$ in §6.7: only finiteness is used; say so.
- **R6 Verifier arithmetic** is outward-rounded doubles, not Arb; kits load sympy pickles. Both go in the trust-base list.
- **R7 Record exception.** `/home/claude/s1_tight/consistency/check_outA5.json` was overwritten by a kit run on Nov 15 (same script and inputs, result unchanged, hash differs from the closed manifest's pin; record `S4/p2/v33/EXCEPTION_v32_B8_check_outA5.json`). Also: the Lemma U / zone III pair is a fixed point whose files cannot pin each other; a rerun of the shell chain changes one provenance string and hence downstream hashes (`V33_INDEX.md` §4.6 open items); stale pins listed there.
- **R8 Nothing in v3.3 has been compiled.** Expect LaTeX errors in the sources the drafters copy from.
- **R9 Reproduction is not yet independent**: WSL only for the v3 kits; kit 3 has one FAIL (provenance hash) and kit 4B one WARN.
- **R10 Number disagreements already visible**: (H3) is quoted as "×15.7 in $\bar Y$" in main_v3 Table 3, as ×5.700 / ×79.26 / ×26.50 / ×7.799 with $\bar Y$ limit ×15.75 in `V33_INDEX.md` §4.3, and as ×5.7 / ×78.7 / ×26.6 in `lemmaU_v34.tex` §10; the far $t^*$ margin is ×15.5; the far H4 margin is ×15.58 (reports print ×15.59). U0 must settle each from the JSON.
- **R11 Leftover wording in sources**: main_v3 line 315 has an open TODO (withdrawn numbers 0.2686, ×4.9, ×5.3); it disappears with D5. `zone3_interface.json` still lists the level-0 far-datum threshold as $5\cdot10^{-7}$ with a probe label (certified value passes both).
- **R12 Steps 1–2 have no kit**; they need CAPD. A reader's rerun path for them is weaker than for steps 3–4.
- **R13 Step-4 document predates the final clean-up** and the 21 later stages; Part I is a real merge, not a copy.
- **R14 Not claimed, to be said explicitly**: convergence to $S_1$; uniqueness of $S_1$; certified first crossing time; maximal curvature of the true curve at $t_1$.

### 7.2 Drafting plan
Units are sized at 2–5 hours for one Opus agent unless marked. All output goes to `W/paper_v1/` (new files only). Each unit delivers its `.tex`, a list of ledger macros used, and an `ISSUES` list. Mechanical units go to Sonnet.

| unit | content | model | depends on |
|---|---|---|---|
| U0 | Number ledger, `numbers_paper.tex`, checker script; resolve R10 from the JSON | Sonnet | — |
| U1 | `preamble.tex`, notation table (5.2), name map (5.3), Part A | Opus | D6 |
| U2 | Trace of D4 (does step 3 use Theorem B?) and of R3; a two-page memo | Opus | — |
| U3 | §2 objects and Theorems A–D, proof of Main Theorem | Opus | U0, U1, D13 |
| U4 | §3 strategy, T2, F4 schematic | Opus | U1 |
| U5 | §4 + Part B (step 1) | Opus | U0, U1 |
| U6 | §5 + Part C (step 2) | Opus | U0, U1, U2 |
| U7 | Part D (tip, verifier) | Opus | U0, U1 |
| U8 | Part E (shell lemma) | Opus | U0, U1, U2 |
| U9 | Part F (Lemma U, zone III); two agents in sequence: §§1–5, then §§6–11 | Opus ×2 | U0, U1 |
| U10 | Part G (Lemma O, R, recovery) | Opus | U0, U1 |
| U11 | Part H + §6.6 (bootstrap closure) | Opus | U7–U10 |
| U12 | §6.1–6.5 (main-paper summaries of Parts D–G) | Opus | U7–U10 |
| U13 | §6.7 continuation and conclusions | Opus | U11, `W/lit_P1/P1_CONFIRMATION.md`, D13 |
| U14 | Part I (step 4 in full); three agents in sequence: $\Gamma_0$ + Phase I; base, switch, Phase II; hand-off and cylinder data | Opus ×3 | U0, U1, D15 |
| U15 | §7 (main-paper step 4), T3, F7 | Opus | U14 |
| U16 | §8 + Part J + Part K, T5, T6 | Opus | U0, D5, D8, D17, kit run results |
| U17 | §1 introduction, abstract, §9, related work, bibliography check | Opus | U3, U12, U15, D12 |
| U18 | Figures F1–F7 (scripts kept with the paper) | Sonnet | U4, D14 |
| U19 | Compile, fix LaTeX, cross-reference and macro audit, run the ledger checker on the final PDF sources | Sonnet | all |
| U20 | **Independent review** (Fable, fresh context, given only the PDFs and the archive): (a) does each main-paper statement match its supplement statement and its pinned numbers; (b) does the Main Theorem follow from A, C, D and (L1)–(L3) as written; (c) a hunt for anything the consolidation changed; verdict and list | Fable | U19 |
| U21 | Fix round, then a second, narrower check of the fixes | Opus, then Fable | U20 |

Order: U0–U2 first (parallel). Then U3–U10 and U14 in any order, one or two at a time (RAM and usage limits). Then U11–U13, U15–U18. Then U19–U21. Rough total: 24 agent runs, about 80–100 agent-hours, plus the two Fable reviews.

Before U3 starts the authors should settle D4 (after U2), D5, D6, D7, D13, D15. D1 and D8 can wait until U16.
