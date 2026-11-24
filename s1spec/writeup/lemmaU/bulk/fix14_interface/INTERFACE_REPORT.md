# INTERFACE_REPORT: fix round 14, the step-4 to step-3 interface (review items C-1, C-2)

Nov 20, 2026. Source review: `fix11_review/REVIEW_v4b_v36_proofs.md` §3.3, items C-1 and C-2. Everything is in this new directory. The proof tree was only read: python was run with `-B`, and after the runs `find -newer` shows no change under `step4`, `step3nl_v4b`, `step3nl`, `step3` or `s1_tight`. Resource use: one process at a time, `nice -n 19`, peak RSS 174 MB, about 1 min of CPU in total.

Labels used below:
- [hand]: hand argument, not machine-checked.
- [sym]: sympy identity.
- [cert]: arb certificate.
- [cited]: taken from another document, with the label it has there. For the Lean theorems, [cited] means the statement was read from source, not compiled.

## Verdict
- **C-1: closed, modulo the inherited caveats** (Lean not built; DKS read secondhand). The step-4 solution is identified with the smooth Lean flow from Γ0 on [0, t(27.15)], and this interval contains t₁. γ(t₁) is a smooth, symmetric, closed immersion. Step 3 continues the same flow. Step 4's two [A] sentences are no longer needed: "smooth for t>0 by parabolic regularity" and "continuation by the literature fact (P1)".
- **C-2: closed.** (GS) is certified at τ₁. It is added to the bootstrap with a certified strict improvement and a large gap, plus hand-written open and closed steps. The range is now **|x| ≤ 60**, which matches the verifier boxes. At τ₁ it is certified on |x| ≤ 60.25 as well.
- **The v37 sketch is not valid as written.** Its global-sup argument fails because |u| ≈ 1.1 near x = 60. The replacement text is in `GS_lemma.tex` §6.

## C-1: `C1_identification.tex`

| step | content | label |
|---|---|---|
| L 1.1 (`C1:lem:G0`) | Γ0 is a smooth, real-analytic, symmetric immersion; it closes because C(b₂*)=0 | [cited] step 4 Prop. Γ0 [cert]; [hand]; [sym] Y10 |
| (I3) | Theorem P2 read as an a-priori theorem for a given smooth flow (the form P2-a). Only the last sentence of its Step 7 used literature. | [hand] (reading of step 4's proof); P2 itself [cited], with step 4's [C]/[D] labels |
| (I4) | Lean (L1), (L2) at n=1; the velocity is −k_ss N; U-1's lemmas on gluing, maximal flow, symmetry and blow-up carry over to F = Γ0 | [cited] (read, **not compiled**, n=1 instance **not elaborated**); [hand] |
| (I5) | DKS (D1). Added: the C-4 completion, i.e. the extension is an immersion at T because ∂_t ln\|Φ_x\| = k k_ss is bounded | [cited] (via P1_CONFIRMATION, not re-read); [hand]; [sym] Y1b |
| L 2.1 (`C1:lem:max`) | maximal smooth symmetric flow γ on [0, T*) with the L² criterion | [hand] |
| L 2.2 (`C1:lem:gauge`) | in step 4's gauge (normalised arclength from the axis point, θ(0) = −π/2, clock dσ = (2π/L)⁴dt), th^γ is in the sine class, solves the P1 equation exactly, and L, X₀, t satisfy the slaved ODEs | [sym] Y1–Y7, Y10; [hand] (parity and rotation index, Parseval for infinitely many modes) |
| (num) | one explicit normal step of a test curve, re-gauged, reproduces the P1 equation (relative defect 3e-7 at dt=1e-7, O(dt)) | float check only (`C1_num.json`), not a proof step |
| L 2.3 (`C1:lem:uniq`) | uniqueness in C⁰H⁵∩C¹L² with **no isolating box**. The Galerkin lemma states uniqueness only inside W; its Grönwall proof does not use W. | [hand] |
| L 3.1, Prop 3.2 | ∫k² is bounded along P1; T* > t₉₆; γ equals the P1 curve on [0, t₉₆] | [hand] |
| L 4.1 (`C1:lem:ke`) | on P2's boxes, slabs 2..180, every cell satisfies \|κ̃ e\| ≤ 1.743e-3, so a = 1 − κ̃e ≥ ½ and J ≥ ½ | **[cert]** `C1_ke.py` → `C1_ke.json` (180 inputs hashed) |
| L 4.2 | ∫k² is bounded on P2's interval (normal-graph curvature formula, boxes, B, a lower bound for λ) | [hand]; [sym] Y8 |
| Prop 4.3 | T* > t(27.15) ≥ t₁; P2's conclusions hold for γ | [hand] |
| Thm C-1 | (i) smoothness and symmetry on [0, t(27.15)], including t=0 and t=t₁; (ii) P1 identification; (iii) P2 identification; (iv) (H1)–(H6) concern γ(t₁), so U-1's (S4.1) and (S4.2) hold; (v) U-1's T_max = T* − t₁ | [hand] |

**Parametrisation, symmetry class and existence time** (the points the task asked to check):
- **Parametrisation.** Lean's velocity is purely normal; step 4 uses a constant-speed gauge. The two are related by an explicit tangential reparametrisation (Y4: d_σα|_p = ∫₀^α(1+v)w′ + α⟨w²⟩).
- **Symmetry class.** Uniqueness gives σ∘γ = γ∘ι. The fixed point of ι stays on the axis, so th is odd: this is the sine class.
- **Existence time.** It is not assumed to reach t₁; it is proved, in two continuation steps (P1, then P2), from DKS's L² criterion and step 4's bounds. For P2, Lean's class matches P2's directly, because P2 uses only normal velocity and the normal graph.

## C-2: `GS_lemma.tex`

| step | content | label |
|---|---|---|
| §1 | range 60 vs 60.25: nothing uses (60, 60.25], so (GS) is posed on \|x\| ≤ 60; GS° (relaxed, open) added for C-5 | [hand] |
| L T (tube lemma) | cellwise δ, η, K with 3Kδ < ½ give a unique C¹ foot X with X′ > 0, \|ku\| ≤ 3Kδ and \|u′\| ≤ (1+3Kδ) tan(η+2Kδ) | [hand]; [sym] Y9 |
| L GS-1 (τ₁) | for every τ* in slab 180 and every S1 in the certified set, on \|x\| ≤ 60.25: sup\|ku\| ≤ **6.52e-3**, sup\|u′\| ≤ **0.0542**; coverage to x ≥ 61.27 | **[cert]** `GS_tau1.py` → `GS_tau1.json`; inputs: step 4's base of slab 180, its S1 chart, its frame P_hull, its P2 boxes, √F₀ of slab 180 [cited] |
| L GS-2 (strict) | under BS_U, on \|x\| ≤ 60: sup\|ku\| ≤ **5.40e-4**, sup\|u′\| ≤ **0.0410**. Uses fin13v4b's asserted boxes (ratios read from the run file and asserted ≤ 1), B1 (√Λ ≤ 5.3086e-8), and on [8, 10.25], where no w′ box is asserted, \|w′\| ≤ 1.415e-3. ψ̂₃, ĥ₂ and k_S1 are bounded cell by cell. | **[cert]** `GS_boot.py` → `GS_boot.json`; derivation on [8, 10.25] [hand]; box assertions [cited] (closure.py l.143, 162) |
| L GS-3 (open, with the endpoint at 60), L GS-4 (closed; GS is in the bootstrap) | | [hand] |
| §6 | review of the v37 sketch (7 points); drop-in replacement text; list of edits owed | [hand] |

## Gaps, stated plainly
1. **Lean.** (L1), (L2) were read, not compiled; the n=1 instance is not elaborated (C-3). Gluing, the time shift and closure under the reflection are hand arguments.
2. **DKS.** (D1) and the C-4 completion rest on the reading in `P1_CONFIRMATION.md`. The paper was not re-read; specifically, its uniform bounds on ‖∂_s^m κ‖_∞ on [0, T) are my reading of the proof's structure.
3. **Theorem P2.** It is used as an a-priori theorem. That is a reading of step 4's proof, and its analytic steps keep step 4's [D] status. Lemma 4.2 uses P2's boxes and Step 2's bound on |x| qualitatively.
4. **Galerkin existence.** Steps 1–5 of step 4's appendix are cited. Only uniqueness is re-proved here.
5. **Parity and Parseval.** The parity and Parseval steps for infinitely many modes are [hand]; the sympy checks use degree-3 polynomials.
6. **GS-2 depends on fin13v4b.** It uses fin13v4b's box assertions, and assumes they persist past the tail anchor (the monotonicity statement in VERIFIER.md, read, not audited). If N1 (D-3) produces a new production run, rerun `GS_boot.py` on that run's rate table. It takes 3 s; change the `RUN` path.
7. **GS-1 inherits** step 4's P2 boxes, its hand-off frame, and Theorem C-1.
8. **Texts not edited.** No closed or draft text was changed. The owed edits are listed in `GS_lemma.tex` §6:
   - `U1_regularity.tex`: set X_S := 60; define 𝒪 by GS° (C-5); gap items 1 and 4 can be closed.
   - `continuation_v37.tex`: the replacement lemma; gap item 2 closed.
   - `lemmaU_v37.tex`: add GS to BS_U.
   - step 4: the (H6) and l.130/336/370 wording, and naming GS in the hand-off list.
9. **No LaTeX compile.** Only brace, environment and `$` balance were checked.

## Files
- `C1_identification.tex`, `GS_lemma.tex`: the two lemmas.
- `C1_sym.py` → `C1_sym.json`, `C1_sym.log`: 26 sympy identities Y1–Y10, all pass.
- `C1_num.py` → `C1_num.json`: float cross-check, ok.
- `C1_ke.py` → `C1_ke.json`: [cert], sup\|κ̃e\| = 1.7428e-3.
- `GS_tau1.py` → `GS_tau1.json`, `GS_tau1.log`: [cert], (GS) at τ₁.
- `GS_boot.py` → `GS_boot.json`, `GS_boot.log`: [cert], strictness.
- `SHA256SUMS`.

Every JSON has an args block: argv, versions, the script's own sha256, and the sha256 of every input. Rerun each script with `nice -n 19 python3 -B <script>`.
