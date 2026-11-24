**Verdict: v3.3 closes step 3 with the listed fixes. No blocking item remains and no inequality fails.** "Unconditional" is not yet justified: it waits on one small missing certificate (N1, the chart term), two literature statements the authors must confirm, and compilation.

Re-review (Fable), Nov 14 2026. Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4/p2`. Labels: **[checked]** = I read the code or text, or probed it; **[agent]** = found by a read-only sub-reviewer (Opus for step 4, Sonnet for documents), not re-checked by me unless stated. Numbers are certified unless marked [est] or [probe].

## 1. Item table

| Item | Status | Evidence |
|---|---|---|
| B0 variable mismatch | **Fixed**, except the chart term (N1) | See §2. |
| B1 proof text of (A) | **Fixed** (text exists; not compiled; not checked line by line) | `B/fix1/b12_text/lemmaS_A_v33.tex`: Prop. level (D1), Thm A, tools, Steps 1–7, integration lemma, (D2) on every window. I checked the integration lemma, the conversion lemma with d, Lemma z1 and the (c)–(e) checks. Stale header, see N3. |
| B2 continuation | **Fixed**, subject to author confirmation of (P1′), (P3′) | `continuation_v33.tex`; see §3. |
| B3 fin11 numbers | **Fixed** | `main_v3.tex:111`: 9.872e-10, 6.867e-7 (sources 9.8717043e-10, 6.8665194e-7). Old values appear only as "do not apply". [agent] |
| M1 outer-ramp Y | **Fixed** | `t5_zone3/ramp_U.py:183–193` reads the primary `Y_next` from `u1_ghat.json` and asserts it dominates. [checked] |
| M2 coupling term | **Fixed** | `t4_cert/assemble_U.py:116–123`. [checked] |
| M3 typed constants | **Partly fixed** | Ramp sups are read from fin12f4 and rounded up (`assembly_rep.py:202–207`); `zone3.py:24, 134` and `u2_K.py:40` are fixed; `ramp_U.py:181` uses an outward box. Still missing: an assert "verifier ramp sups ≤ ramp box" (only a reported margin, ×1.185). [agent, lines seen by me] |
| M4 provenance | **Fixed** | u1 and u34 pin the fin12f4 file and the current `t4lib.py`. Left over: `zone3_interface.json` pins the pass-2 `U_cert_v34` (9e6bc84a, differs in args only); `ramp_U.json` holds a stale `u1_ghat_sha`. [agent] |
| M5 lemmaU text | **Fixed** | `lemmaU_v34.tex:277–282` (one Young step on √(a²+b²)·N₄), `:383` (B2′ with d), `:146` ((W2) rerun). [checked] |
| M6 margins at N2_ref | **Partly fixed** | Texts quote the 4·N2_ref run. `assemble_U.py:245` stores `N2_le_4N2ref` but does not assert it; only `check_keys_v34.py` does (ratios 0.18, 0.41, 0.25). [agent] |
| M7 main_v3 wording | **Partly fixed** | (H4), (W3), (W4) are fixed. Not fixed: `main_v3.tex:312` still says Lemma N derives (O1)–(O3); `:309` uses R(τ) undefined; `:315` keeps "all repaired margins at least ×2" with no correction. [agent] |
| M8 erratum scope | **Fixed** | `ERRATUM.md:179, 314`; ×0.7317 is now [est]. [agent] |
| M9 step-4 items | **Fixed** | `win_arb.py` (arb element length, remainders, sums; coverage asserted); arb base sums with Minkowski error transfer; centre map asserted; far H4 by a contiguous row cover of every centre ≥ 490 (min ×15.5878); (H2f)₀ is 1e-6 in code. Both manifests verify (closed: 0 mismatches; v33: 19 stages, 8f6f7f95…). [agent] |

## 2. B0 in depth

- **d is bounded correctly.** `S4/v33/b0_dprof.py`: exact dyadic cells covering [14, 60] (ends asserted), mean-value form with the derivative enclosed on the whole cell, constants from strings, same operator and weight as the shells (D = y∂_y, 28/y, ds). [agent]
- **d enters at every interface** [checked in `assembly_rep.py`]:
  - ramp J conversion: `trq` (quadratic), `trc` (cubic, and the convexity inequality is right), `trp` (products);
  - Φ, foot and S33: yes;
  - fluxes: (√K(E) + √K(E_d))², valid because √K is a seminorm;
  - S0V/S1V and the H-B2 claims: + sup_core|d|, + sup_core|Dd|;
  - ramp boxes: yes;
  - datum: `h2u_v33`, Minkowski, ×27.34 / ×562.8 / ×600.9;
  - B2′: consistent in the three texts.
- **Lemma U needs no d.** It takes (A)'s export and δU only. [checked]
- **Chart identification (B0_REPORT §7): the datum is fine, the ramp is not certified.**
  - Datum: `t9h/win.py:28–72` computes the ray-graph quantity itself at every level (f = tan(θ̃ − Θ), minus the profiles at y − σ_tip). So (H2) needs no chart term at levels ≥ 1.
  - (A)'s equation is exact for a graph over the ray; it does not need S₁ to be a ray.
  - What is missing: on the ramp [14, 28], (A) uses the verifier's S₁-chart w′ (J_m, sI1, M1) as ray-graph data. The two differ by δ(x) = θ₁(x) − Θ and its derivatives, to leading order. No lemma or script bounds this; `win.py:70` types 2e-7 with a comment.
  - Size [probe, `scratch2/chartprobe.py`, step-1 envelopes `hstail.profile_bounds`]: at σ = 14, |δ| ≤ 2.2e-8 (so the typed 2e-7 is valid), y·k ≤ 4.4e-7, y²·k′ ≤ 9.4e-6, y³·k″ ≤ 2.0e-4; at σ = 16 these are ×20 smaller; at σ = 28 below 1e-13.
  - In (A)'s W-norm bookkeeping the m = 3 chart term is about 1.4 times d's own m = 3 entry and the m = 2 term about 3% of d's [est]. d moved Ê by at most 0.6%, so the certified margins are not at risk.
  - In the other direction (S0V/S1V, fluxes, y ≥ 32) the term is below 1e-17.

## 3. Continuation

- **Complete and, as far as I can check, correct.** Lemma G (graph inside the bootstrap), the region table, τ_max = ∞ by a bound on total ∫κ² ds, T < ∞, tip blow-up in L², the single singular parameter from B3U(b) and B4. No level-3 or pointwise-z₂ quantity is used.
- **The criterion is the right one.** The proof bounds the total L² curvature of the whole closed curve, which is exactly what "T_max < ∞ ⇒ ∫κ² ds → ∞" needs. Two cautions:
  - I recall the statement as Dziuk–Kuwert–Schätzle 2002, Theorem 3.1 (closed immersed curves in ℝⁿ, ∫κ² ≥ c(T − t)^(−1/4)), and as Chou 2003. I could not check either here. `main.bib` has both entries; the text cites only Chou and Wheeler.
  - The Wheeler lifespan theorem is not a drop-in alternative as written; drop it or state it exactly.
- **Z2 for all τ.** Lemma G uses |f| ≤ 0.2276 on [253, 700]λ₁ at every τ, while `step3_v33_main.tex:59` states Z2 "for τ − τ₁ ≤ 4.0708". The bound comes from the Lemma O chains and holds over the whole horizon, so the mathematics is fine; the wording must say so.
- **Coverage of the whole closed curve** now rests on `h4far_v33` (every centre ≥ 490, far end and mirror arm). [agent]
- **Small wording gaps:** "r ≥ r₀/2 ≥ 56λ ⇒ y ≥ 56" needs a constant relating distance to y (use 60λ, say); c_g, C_S₁, ℓ_tip are finite constants, correctly not called certified.

## 4. G-A1

- **The rule is valid.** With M = sup E: E ≤ A + 2(F/γ′)√M gives √M ≤ F/γ′ + √(A + F²/γ′²), and squaring gives (R-sup). The coupling coefficient 1 is right for a single integration at rate γ′; levels close in order 0, 1, 2. [checked]
- **Applied consistently:** `assembly_rep.py:169–173`, `export_t2.py:36–41` (exports `A0_Rsup`, `Fg`), `assemble_U.py:116–123`. `rule_reproduces_EA_ref` is true. [checked]
- Minor: the coupling increment subtracts the upward-rounded `EA0`; the deficit is below 1e-22.

## 5. New findings

**Blocking:** none.

**Must-fix**
- **N1. Certify the chart term on the ramp.** Add a stage that bounds D^k(θ₁ − Θ), k ≤ 3, on σ ∈ [14, 28] from `profile_bounds`, with the nonlinear terms (k·u·u′, the shift u·sin δ), and add it to the d tables (W-norms and ramp sups). Replace the typed 2e-7 by the computed value. Write the identification as a lemma in `lemmaS_A_v33.tex`. Expected move well below 1% [est].
- **N2. Step-4 environment not recorded.** `final_tail_v33.py` passes the parent environment to the stages; `s1chart.py` reads `S1TT`, yet `b0_dprof` records `env={}`. Clear or record it. [agent]
- **N3. Stale citations presented as current.** [agent]
  - `lemmaS_A_v33.tex:19, 282` and `step3_v33_main.tex:41, 59`: 04947b7f, "12 stages", 8202d2c6, 8152c650, and the M9 caveats now closed.
  - `lemmaS_A_v33.tex` §boxes and G-A2 still say the ramp values are typed.
  - `ERRATUM.md:350–353`, `erratum.tex:236`: Lemma U values from `U_cert_v33.json`.
  - `ic_rep.py:20` and `make_constants_v34.py:181` cite the float outputs, not the `_arb` ones.
- **N4. Remaining M3, M6, M7 parts** (table above): the ramp-box assert, the N̄² ≤ 4·N2_ref assert in `assemble_U.py`, and the three `main_v3.tex` wordings.
- **N5. Cite the blow-up criterion precisely** in `continuation_v33.tex` and `main_v3.tex:130` (which still says the (P1) references are to be supplied), and state Z2 for all τ.

**Minor**
- Non-outward roundings in prose: `main_v3.tex` H4 1.286e-5 (source 1.2864e-5); M9_REPORT ×15.59 (15.5878) and one table row; B0_REPORT 3.657e-10, 1.664e-6. [agent]
- A `.pyc` was written into the closed tree on Nov 14 (`S4/t9h/__pycache__/`); it is not pinned. [agent]
- `M9_base_sums` hard-codes `ok=True`; some imported modules are missing from the stage input pins. [agent]
- The far-end flattened-weight step of `h4far_v33` (LEMMA_F §1.3) is applied, not asserted; margin there about ×8e6. [agent]
- Literal "[TODO: source …]" in the erratum; pin leftovers under M4; nothing is compiled.

## 6. Is "unconditional, subject to review" justified?

Not yet. After this review step 3 v3.3 is conditional on exactly:
1. N1, the chart term on the ramp (a missing certificate with large room, inherited from v3.1);
2. (P1′), the L² blow-up criterion for closed immersed curves, and (P3′), interior regularity (literature, authors);
3. the step-4 hand-off as pinned (both manifests verify);
4. the doc fixes N3–N5 and a compile.

## 7. For the authors to confirm by hand

- **(P1′):** DKS 2002 Thm 3.1 and/or Chou 2003: smooth closed *immersed* plane curves, T < ∞ ⇒ ∫κ² ds ≥ c(T − t)^(−1/4). Cite the theorem number.
- **(P3′):** interior regularity for graphs with bounded slope, bounded curvature and κ_s in L² (used only for smoothness of the far curve up to T).
- **(P2′) and the frame conventions** (v3.2 §6 item 5), unchanged.
- **Lemma O read arc by arc** for an immersed curve (inherited).
- **The single-point statement:** the singular set is one parameter value; other arcs pass through p_* at T.

## 8. Not checked

- `lemmaS_A_v33.tex` Steps 2–5 line by line against `weights_t2.py`, `ramp_rep.py`, `outer_L2.py`.
- The verifier internals, and whether it uses the same C⁴ ψ̂₃ as `c4fun.py`.
- The error-transfer constants inherited from v3.2; the `win.py` formulas.
- No certified stage was rerun; I relied on the hash checks and the reported byte-identical passes.
