# Name map: archive name -> paper name (source of supplement Part K) -- unit U1

Authors' decision D6 (rename, map in Part K), D5 (no erratum in the paper), D7 (consolidated step-3 text), D15 (step 4 written for the final hand-off list). Rules of OUTLINE §5.3 apply. Codes below that are **not** in the paper appear only in Parts J and K ("K only").

How to read: "archive" = the name in the files under `/home/claude` (v3.2/v3.3 sources named in OUTLINE.md); "paper" = the name the drafters must use in the main paper and in Parts A--I. A "proposed" name is my suggestion where OUTLINE fixes none; the unit that writes that part may change it and must then update this file (and tell the coordinator).

## 1. Hand-off hypotheses of Theorem C (archive: Theorem W)

The archive numbering is the same in v3.2 and v3.3 except that v3.3 inserts (H2f). `step4_proof.tex` still uses the v3.2 numbering (it has no (H2f)); the v3.3 step-4 stages use the v3.3 numbering.

| paper | archive (v3.3, `step3_v33_main.tex`, `main_v3.tex`) | archive in `step4_proof.tex` (v3.2) | content |
|---|---|---|---|
| (H1) | (H1) | (H1) | modulation and smallness: decomposition exists at $\tau_1$, $\abs\alpha\le\alpha_1=1.5\cdot10^{-4}$, $\sqrt\Lambda\le5\cdot10^{-8}$ |
| (H2) | (H2), "tightened (H2)" | (H2) (thresholds $(2\cdot10^{-7},10^{-4},0.2)$ are superseded; K only) | shell energies on $y\in[14,60]$ for $f-U_1$ |
| (H3) | **(H2f)** | -- | far datum of Lemma S: $G_j$, $j=0,1,2$, on $[60,700]$ |
| (H4) | (H3) | (H3) | windows $[y_w,1.5y_w]$, $y_w\in[50,253.125]$; $K_4\le0.06$ on cone cylinders |
| (H5) | (H4) | (H4) | graph over the $S_1$ ray on $y\in[48,700]$ (v3.2: $[48,350]$) |
| (H6) | (H5) | (H5) | Lemma O data, arc by arc: hypotheses (O1)--(O5), (O7), (O9) (archive: H1--H7, H9) for the fixed-$\ell$ cone family and the far family |
| (H7) | (H6) | (H6) | $\gamma(t_1)$ smooth; forcing norm determined by $\alpha(\tau_1)$ |

Trap: "H4" in file and stage names `h4loc`, `h4loc_v33`, `h4far_v33`, `h4loc_A1_w180.json`, and "H4 ×4.18 / ×5.68" in CLAUDE.md means **Lemma O H4 = (O4)**, not the hand-off graph hypothesis. Similarly `h2f_*` = (H3), `h2t_*`/`h2u_*` = (H2), `h3_arb` = (H4) with the v3.2 (H2) shells.

Other hand-off names: "(H1)--(H6) of step 3 v3.2" = (H1), (H2), (H4)--(H7) (K only); "hand-off conditions" = (H1)--(H7).

## 2. Lemma O: hypotheses and results

| paper | archive (`lemmaO_v32.tex` §"Hypotheses the step-4 hand-off must check", `lemmaO.tex`) | content |
|---|---|---|
| (O1) | H1 | graph over the cylinder, $\sup\abs{h_0}\le H_{\rm cyl}-C_0\sqrt{t_1}$ |
| (O2) | H2 | slope $\le\frac14-C_1\sqrt{t_1}$ (data: $\le s_0=0.08$) |
| (O3) | H3 | curvature jets $\abs{h_0^{(j)}}\le M_O$, $j=2,3,4$ |
| (O4) | H4 | $\int(h_0^{(5)})^2\le M_O^2$ ("good time") |
| (O5) | H5 | data thresholds (slope, $\sup\abs{h_0''}$, $\norm{h_0'''}_{L^2}$; the "six data thresholds") |
| (O6) | H6 | **no longer a hypothesis**: the a-priori bounds close; proved as the Proposition below. Label (O6) is left unused (open question Q3 in notation.md) |
| (O7) | H7 | time: time ratio $t_1/\ell^4\le$ needed (cone $50^{-4}$, far $6.134\cdot10^{-9}$) |
| (O8) | H8 | other arms (not a requirement; arc by arc) |
| (O9) | H9 | regularity of $h$ on $[-X,X]\times[0,t_1]$ |

| paper | archive |
|---|---|
| Lemma O | Lemma O (v3, v3.2; `lemmaO_v3.tex`, `lemmaO_v32.tex`; v1, v2 = `lemmaO_v2.tex` are K only) |
| Proposition (recovery of the cylinder bounds) | "H6", `thm:H6` of `lemmaO_v32.tex`, "the H6 theorem", "D5b/H6", "H6 under v3.1" (K only), "U-fine H6 entries" |
| Lemma R | Lemma R (`lemmaR.tex`, `thm:R`); its `cor:h4` = "Corollary (of Lemma R)" |
| fixed-$\ell$ cone family | "fixed-$L$ cone family" ($L=50\lambda_1$, $X=3$, centres $[300,380]\lambda_1$, $(M,K)=(0.12,0.18)$) |
| far family | "far family" (centres $\ge490\lambda_1$, $\ell=113\lambda_1$, $(M_O,K_O)=(0.15,0.3)$) |
| $(M_O,K_O)$ | $(M,K)$ of Lemma O |

## 3. Literature facts

| paper | archive | content (OUTLINE §4.1, `lit_P1/P1_CONFIRMATION.md`) |
|---|---|---|
| (L1) | (P1$'$); older (P1) of v3.2 (K only) | (a) local well-posedness with uniqueness for smooth closed immersed plane curves (reference not yet confirmed); (b) the $L^2$ blow-up criterion, DKS 2002 Thm 3.1, Chou 2003 Prop. A. OUTLINE treats (a) and (b) together as (L1); if the drafters split it, call the parts (L1a), (L1b) and tell the coordinator |
| (L2) | (P2$'$) | regularity of the modulation (reference to be confirmed) |
| (L3) | (P3$'$) | interior regularity for graphs (reference to be confirmed) |

Note: step 2 alone cites "order of an Evans zero = algebraic multiplicity" as in the step-2 text; it is not a numbered fact and matters only if D4 fails.

## 4. Theorems, propositions and lemmas

### 4.1 Top level

| paper | archive |
|---|---|
| Main Theorem | "Main Theorem" (`main_v3.tex` §1.3; `main.tex` "conditional on Hypothesis (F)" is K only) |
| Theorem A | step 1: `thm:S1` (`s1_proof.tex` Thm 1, `s1_tight.tex`) |
| Theorem B | step 2: `thm:step2` / `thm:main` of `step2_proof.tex` |
| Theorem C | "Theorem W", `thm:W`, "trapping near $S_1$ with a floor, version W/v3.1/v3.2/v3.3" |
| (W1)--(W5) | (W1)--(W5), `continuation_v33.tex` §5.4 / `C:sec:W` |
| Theorem D | step 4: "Theorem 4", `thm:4` |
| Proposition ($\Gamma_0$ is closed and embedded) | `prop:G0` (`main_v3.tex`, `step4_proof.tex`) |
| Phase I, Phase II | P1 (validated Fourier slabs), P2 (energy shadowing) of step 4 |
| Theorem (shadowing) | "Theorem P2 (shadowing; consolidated from T5 §1.3)" |
| (B1)--(B5) | B1, B2$'$, B3U (parts (a),(b),(c) kept as (B3)(a)--(c)), B4, G (Definition G of `continuation_v33.tex`) |
| (K1)--(K6) | (K1)--(K6) "closure conditions" of `lemmaU_v34.tex`, Prop. "closure" |
| $\mathrm{BS}$ | $\mathrm{BS}_U$ (the bootstrap set; K only for the subscript) |

### 4.2 Step 3: Lemma S (shells), Lemma U (far field), chart lemma, zone III, continuation

| paper | archive |
|---|---|
| Lemma S | "repaired Lemma S(A)", `A:thm:A` of `lemmaS_A_v33.tex`; v3.1 Lemma S(A) (`lemmaS.tex`) is K only |
| (withdrawn, K only) | Lemma S(B), Proposition B$'$ (v3.2 `prop:B`), Hypothesis (F), (F1), (F2), "window energies", "strip" energies of v3.2; replaced by Lemma U |
| (withdrawn, K only) | Hypothesis (O) (outer cutoff of (A) on $[350,700]$): now inside Lemma S via the outer annulus and Lemma N |
| Lemma U, Proposition U | Lemma U; `prop:U` of `lemmaU_v34.tex` |
| Proposition (closure) | `prop:close` "closure" |
| Proposition (U-fine) | `prop:fine` "U-fine" |
| Lemma T (transport) | `lem:T` of `lemmaU_v34.tex` (not Lemma T of `STATEMENTS_r2.tex`, see next) |
| Lemma (annulus tools) | `A:lem:tools` "tools on $O$" (= `STATEMENTS_r2.tex` Lemma T). Proposed name; avoids the second Lemma T |
| Lemma (sup-free bounds) | `lem:S` "sup-free cubic and quartic bounds on an interval" (= "Lemma S" inside `lemmaU.tex`; **not** the shell Lemma S). Proposed name |
| Lemma ($\Phi$) | `A:lem:Phi` ($\Phi$, "Lemma S (sup-free cubic)" with $\Xi=\varepsilon\omega$) -- the same lemma as the previous row specialised? Not checked: open question Q2 in notation.md |
| Lemma (profile mismatch) | `lem:dU` |
| Lemma (launch data) | `lem:launch` |
| Lemma D | `lem:D` and `A:lem:D` "D, window form" (one or two lemmas: Q2) |
| Lemma HJW | `lem:HJ` and `A:lem:HJ` "H-J-W" (expansion of the initials not in the text I skimmed; U9 must state it; one or two lemmas: Q2) |
| Corollary (layer cake) | `cor:abel` |
| Definition (sweeps) | `def:SX` |
| Proposition (data term) | `prop:Xterm`, `STATEMENTS_r2.tex` Prop. 8 |
| Lemma (comparison) | `lem:cmp` |
| Lemma (window Sobolev) | `lem:W` |
| Lemma N | `lem:N` (outer-annulus data); its hypotheses H-O1, H-O2$'$ become (Na), (Nb) (proposed; not "(N1)" because N1 is a review code) |
| Definition (loop gain) | `def:gain` |
| Definition (one frozen profile) | `def:v`, "change (C5)" |
| (U1)--(U6) | hypotheses (U1)--(U6) of Proposition U (box, rate conditions, (A)'s data, neighbour bound, zone-III data, modulation bounds); keep, Part F only |
| Lemma (core slope bound) | `A:lem:z1` "$\abs{z^A_1}$ on the core; review M5" |
| Proposition (level-$j$ inequality) | `A:prop:level` |
| Lemma (divergence form) | `A:lem:id` (`lemmaS.tex` Prop. 1) |
| Lemma P | `A:lem:P` (ramp-end jump $J=120/\ln^32$; `STATEMENTS_r2.tex` Lemma P) |
| Lemma (L-term) | `A:lem:L` (`STATEMENTS_r2.tex` Lemma L; proposed name to avoid "Lemma L" next to $L$ and (L$n$)) |
| Lemma (conversion) | `A:lem:conv` |
| Corollary (exponential weights) | `A:cor:exp` |
| Lemma (integration) | `A:lem:int` "integration of (D1)" |
| Lemma (chart identification), Lemma (chart bounds) | `C:lem:id`, `C:lem:bd` of `chart_lemma_body.tex`; "chart lemma", "N1 chart term" |
| Definition G, Lemma G | become item (B5) and "Lemma ((B5) is part of the bootstrap)" (`C:def:G`, `C:lem:G`) |
| Proposition (global existence in rescaled time) | `C:prop:global`; v3.2 `prop:global` (Prop. 4.1) is K only |
| Theorem (constrained coercivity) | `thm:lin` of `step3_proof_v32.tex`; "Lemma 3.1" in `step3_linear_certificate.tex` (not checked) |
| Zone III; Z1, Z2, Z3 | "zone III", Z1--Z3, "D5a" (Part F; region table T2) |
| region names | inner ramp $[14,28]$; core $[28,60]$; segments $(60,350]$; outer annulus $O=[350,700]$; zones I/II; zone III $[253,379.5]\lambda_1$ -- kept |

### 4.3 Step 4 (Part I)

Step 4's lemma names clash with step 3 (Lemma B, J, N, ...). Proposed paper names use "II" for Phase II; the units U14/U15 decide.

| paper (proposed) | archive (`step4_proof.tex`) |
|---|---|
| Lemma (Galerkin limit) | `lem:gal` "first exit and Galerkin limit" |
| Lemma (switch) | `lem:switch` |
| Lemma (corrected base and residual) | `lem:SE` (T6; "(SE)-corrected base"; LEMMA_F.md §6.1) |
| Lemma (constant majorant) | `lem:Y` |
| Lemma II-a, II-b, II-c | Lemma 4a (tip nonlinearity without a pointwise bound on $e''$), 4b (pointwise tip box from $\Lambda_e$), 4c (nonlinear closure); `4a'`, `4b'`, `O'` are variants (K only) |
| Lemma II-F, II-B, II-J, II-N | Lemma F (far energies of the error), B (transport defect), J (time-weighted ramp data), N3$'$; the paper's Lemma B is not used; note the archive Lemma F is not "Lemma F" of any step-3 file |
| Corollary II-C | Corollary C (comparison and hand-off) |
| Lemma (H1 on slab 180) | `lem:H1` |
| Lemma ($t_1$) | `lem:t1` (`t1clock.json`) |
| Lemma (Lemma O (O4) per cylinder) | `lem:H4` (`h4loc`) |

## 5. Derivation steps and quantities of the step-3 bootstrap (Part H only)

| paper | archive |
|---|---|
| D1, D3, D5a, D5b, D6 (Part H only) | D1 verifier; D3$'$ = Lemma S(A) repaired; D5a zone III and Lemma O; D5b/H6 recovery; D6 barrier (tail); "U-crude", "U-fine" are the two passes of Lemma U (keep as "Lemma U, crude pass / fine pass") |
| verifier | "D1", `step3nl_v33`, tag `fin12f4` (`fin11` superseded); K/J only |
| $\bar a=\bar b$, $\bar\beta$ | `modsups`, D1 outputs; same symbols |
| $\hat E$ | $\hat E^A$ (`\Eh`), $\hat E=(1.2018\cdot10^{-6},1.0441\cdot10^{-4},0.056694)$ |
| $\hat B$, $\hat B^{\rm fine}$ | $\hat B$ (`\Bh`), $\hat B^{\rm fine}$ |
| $\beta_A$ | $\beta_A=2$ (flux multiplier); keep, Part E |
| profile difference $d$ | "B0 route (b)" ($d=v_A-w'$), `b0_dprof` |
| (G-A1), integration rule (R-sup) | K only |

## 6. Step-4 stages (Part J/K, never in the main paper)

Stage names stay as file/stage names in Part J. For the paper they are described by the hypothesis they serve.

### 6.1 `FINAL_manifest_v33.json` (21 stages, sha256 8bc47e83...)

| archive stage | serves (paper) |
|---|---|
| h2f: widened parameter balls (`h2f_ift`) | (H3) |
| h2f: base sums on $[60,700]$, 32 $t_c$ balls (`h2f_balls`, `h2f_win`) | (H3) |
| h2f: error-transfer constants (`t9_h2f`) | (H3) |
| h2f: (H2f) for the true curve at $\bar Y=10^{-12}$ (`h2f_camp`) | (H3) |
| tightened (H2) shells (`h2t_v33`) | (H2) (float row; superseded by the arb row `h2t_v33_arb`) |
| fixed-$L$ cone: jets, (H5), H3 (`t9_cone_v33`) | (H6) data thresholds and (O3) |
| fixed-$L$ cone: H4 (`h4loc_v33`) | (O4) cone |
| fixed-$L$ cone: Lemma O H1 graph (`t9_lemmaO_v33`) | (O1) cone |
| fixed-$L$ cone: H7 time ratio (`ratio_v33`) | (O7) |
| covering notes (`notes_v33`) | -- (record) |
| B0: profile difference $d$ (`b0_dprof`) | (H2), Lemma S |
| B0: tightened (H2) for $f-U_1$ (`h2u_v33`) | (H2) |
| M9: arb window/shell/(H2f) base sums (`m9_balls`, `win_arb`) | (H2)--(H4) |
| M9: (H3) windows and v3.2 (H2) shells on arb bases (`h3_arb`) | (H4) |
| M9: tightened (H2) shells, arb (`h2t_v33_arb`) | (H2) |
| M9: tightened (H2) for $f-U_1$, arb (`h2u_v33_arb`) | (H2) |
| M9: (H2f) on arb bases (`h2f_camp_arb`) | (H3) |
| M9: Lemma O H1 with certified centre map (`t9_lemmaO_v33_arb`) | (O1) (H6) |
| M9: far-family H4 at every centre $\ge490$ (`h4far_v33`) | (O4) far |
| N1: chart terms (`n1_chart`) | chart lemma / (B2) |
| f12_sliver: F1, F2 sliver constants (`f12_sliver`) | Lemma S conversion |

### 6.2 `FINAL_manifest.json` (34 stages, sha256 e01cb4da...; closed, pinned)

Titles as printed in the manifest, grouped (paper: Part J table; the main paper says "the 34 stages of the closed chain"):

- Phase II base, switch and modulation: T8 cconst, T8 merge, switch (P1 to P2 at slab 2) T6 = t6c, ycomp ($Y_0=6.21\cdot10^{-13}$, slab 2), clam.
- (H1): H1 exact points (16 points), H1 slab summary.
- Hand-off passes A and B: T8 camp (no N3), T8 n3far it0, T8 camp --n3 it0, T8 n3far (final boxes), camp --n3 n3_b1 --window 180-180 (T9WIN), n3boxes; t9.py (cone jets, KCH, c0, c0tip, shells, frame); pass A and pass B differ by the ball set.
- Far: t9_far.py; v3.2: t9_far.py --v32; v3.2: t9_v32_checks.py.
- Lemma O: h4loc on the final amp, Lemma O H1/H8 (`t9_lemmaO.py`), Lemma O H6 (`t9_h6.py`), v3.2: h4loc FAR_M=0.15, v3.2: `t9_h6_v32.py`, v3.2: `t9_lemmaO.py --v32`, Lemma O H1 far family (`t9_lemmaO_far.py`).
- Others: axis crossing (informational), t1clock, S1 hull containment (B8; `check_outA5.py`; run last).

## 7. Review, version and fix codes (K only)

None of these appears in the main paper or Parts A--I (OUTLINE §5.3); each gets one line in Part K's review list (U16 completes meanings from the review files).

| code | where seen | what it names |
|---|---|---|
| v3.1, v3.2, v3.3, v34, v35, "step3 v3.3" | file names | versions; the paper presents one proof |
| (F), (F1), (F2), "Hypothesis (F)", "former Hypothesis" | `main.tex`, `main_v3.tex` §1.3 | far-field window estimates of v3.2 (replaced by Lemma U) |
| (O) | `main_v2.tex` | outer cutoff of Lemma S(A) on $[350,700]$ (replaced by the annulus treatment and Lemma N) |
| A1--A10 | `erratum/ERRATUM.md` | the v3.1 Lemma S(A) defects |
| B0 | `fix1/b0/` | variable mismatch between verifier $w'$ and (A)'s variable; resolved by $d$ (route (b)) |
| B1, B2, B3 (blocking, v3.3 review) | `t8_review/REVIEW_v33.md` | blocking review items; not the bootstrap (B1)--(B5) |
| M1--M9 | `fix1/M14_REPORT.md`, `step4/p2/v33/M9_REPORT.md` | review items of the v3.3 chain (M3: ramp sup assertion; M5: core slope; M6: (A) rows quoted at $\beta=2$; M7: (H4) range; M9: arb stages) |
| N1--N5, N35 | `t8_review/REREVIEW_v33.md`, `fix2/` | re-review items (N1: chart identification; N2: per-stage environment) |
| F1, F2, F3, F4 | `fix2/F34_REPORT.md` | sliver items and docs |
| G-A1 .. G-A5 | `lemmaS_A_v33.tex` §gaps | bookkeeping items of the shell lemma (R3 in OUTLINE) |
| step-4 review A1--A3, B1--B9 | `STEP4_REVIEW.md`, `step4_proof.tex` §review | A1 one-slab good-time window; A2 nonlinear terms allowance DN; B1 (SE)-corrected base; B4, B5 shells/windows, KCH and the chart; B7 hand-off time $t_1$; B8 $S_1$ chart and certified set |
| T1, T3, T5--T9 | step-4 task reports | task names of step 4; K only |
| R1--R14 | OUTLINE §7.1 | risk list of this plan, not archive names |
| D1--D18 | OUTLINE §6 | author decisions, not archive names |

## 8. Documents and kits (for Part J/K)

| archive | paper |
|---|---|
| `s1_cap/s1_proof.tex`, `s1_tight/s1_tight.tex` | Part B (step 1) |
| `step2_proof.tex` | Part C |
| `step3_proof_v32.tex` + `step3_v33_main.tex` + `lemmaU_v34.tex` + `lemmaS_A_v33.tex` + `continuation_v33.tex` + `lemmaO_v32.tex` + `lemmaR.tex` | Parts D--H (consolidated, D7) |
| `step4_proof.tex` (39 pp) + `STEP4_RESULT_v33_addendum.md` | Part I |
| `kit_step3_v3.zip` (e40b444a), `kit_step4a_v3.zip` ("kit A"), `kit_step4b_v3.zip` (646f2552, "kit B") | Part J Tiers 1 and 2 |
| `s1_cap/run_all.sh`, `s1_tight/run_tight.sh`, `run_step2.sh` | Part J Tier 3 (no kit) |
| `EXCEPTION_v32_B8_check_outA5.json` | Part K record exception (D17) |
