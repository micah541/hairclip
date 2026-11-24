# Fix round 1, B1 and B2: proof text of the repaired Lemma S(A), the v3.3 continuation, and the assembled step-3 document

Nov 14, 2026 (Opus). I wrote only in `fix1/b12_text/`. The LaTeX is not compiled, because there is no LaTeX on this machine. The balance check (`t3_proof/check_balance.py`) passes for every file and for the expanded main document. The B0 slots are filled with the route-(b) values (v34), per the coordinator's message and `fix1/b0/B0_REPORT.md`.

## 1. What was written

| file | content |
|---|---|
| `lemmaS_A_v33.tex` | **B1.** A full proof text of the repaired Lemma S(A). Contents: setting; the energy variable f − U with the profile difference d (four interfaces); boxes by region; Proposition A, the level-j inequality (D1); Theorem A (Ê, Lemma D window export, core box, ramp box, fluxes). Tools: the divergence-form identity, Lemma P (point terms), Lemma L (localisation), conversion with d, H-J-W, the exponential-weight corollary, Φ. The proof runs step by step: transport and inflow; base terms and Ĝ with one C^{1,1} Ỹ and the segment/O rate conditions; the y = 28 interface (A5); the ramp terms R1/R2/R3 with the point-term routes (A3); foot and S33; R5′ on O (tools T1–T7, the class/tool table, the reserve); forcing. Then the integration lemma, the datum H2 + H2f, why (D2) holds on every window, the (c)/(d)/(e) checks, and a gaps section. Every number is a macro with its JSON key. |
| `continuation_v33.tex` | **B2.** Lemma G puts the graph property inside the bootstrap. Then: uniform bounds region by region (which norms); (P1′) in L² form (Chou's blow-up criterion); Proposition τ_max = ∞; λ → 0, T < ∞ and the rate; a single singular point with B3U(b) in place of the (B) windows; smoothness of the far curve; (W1)–(W5) in v3.3 form (fin12f4 numbers); a table mapping each v3.2 use of B3/(B)/level 3 to its v3.3 replacement; gaps. |
| `step3_v33_main.tex` | The assembled v3.3 step-3 driver: setting; Theorem W v3.3 with (H1)–(H6) and (H2f) (and (H2) for f − U₁); unchanged components; verifier; `\input{lemmaS_A_v33}`; Lemma U summary; zone III and the cone family; the bootstrap set, order table and (K1)–(K6); `\input{continuation_v33}`; what remains. |
| `numbers_A_v33.tex` | **All numerical values**, one macro each, with source file, key and label. These are the v34 values from `fix1/b0/t2_A_v34/out/A_interface.json` (04947b7f…), `asm_design_N2ref4.json`, `U_cert_v34.json` (4c3f23bc…), `b0_dprof.json` and `h2u_v33.json`. |
| `b12_preamble.tex` | The shared preamble. Labels: [cert] [paper] [hyp] [est] [B0]. |

The two sub-texts compile standalone or inside the driver (`\ifdefined\stepmain`). Labels are prefixed A:, C: and M:.

Read-only probe (scratchpad; it imports `ramp_rep.py` and writes nothing to the project):
- m = 3 occurs only on the pair (j, j+1);
- on (a, j+1) with a < j, m ≤ 2 (m ≤ 1 at level 2);
- the coefficient κ of every (2,2,3) class has (m, q) = (1, 0) only.

So there are no hidden point terms in R1, R2 or ibp223.

## 2. Gaps found

**G-A1 (new, bookkeeping, small; affects the certified Ê).** The (A) assembly rule mixes two different integrations of (D1):

E^scr = H + D(J + s/γ′) + 4F²/γ′² + 2ĉ/γ′ + ô/γ′ (+ t_d/γ′)

- The terms 4F²/γ′² and 2ĉ/γ′ come from absorbing 2√E·F with rate γ′/2.
- The terms ô/γ′, t_d/γ′ and J + s/γ′ use the full rate γ′.

Neither valid route gives E^scr. The same mixture is already in v3.1 (`lemmaS_F1.py`:148, `assembly_arb.py`:126).
- The valid sup-route is (R-sup): E ≤ A + 2F²/γ′² + 2(F/γ′)√(A + F²/γ′²).
- It gives Ê^true ≤ Ê^scr + 2(F/γ′)√Ê^scr ≤ (6.8e-9, 6.6e-8, 3.2e-6) [est, from certified components]. That is 0.57%, 0.063% and 0.0055% relative.
- Effect: FLUX0/K0 at β_A = 2 goes from ×3.382 to about ×3.36. (K6) is unaffected. Lemma U's inputs move by less than 0.6%.
- Fix: one change in `assembly_rep.py` (add the cross term, use 2F²/γ′²), the same in the rule `EA_of_N2` and in `assemble_U.py`, then rerun t2 → t4 → zone3.
- The text states that the quoted Ê are script values until that rerun.

**Continuation (B2): no quantitative gap.** No level-3 or pointwise-z₂ quantity is needed anywhere. The items:
1. **(P1) must be stated in L² form** [hyp, literature]. Use Chou 2003 (if T_max < ∞, then ∫κ² ds → ∞), or Wheeler 2013's lifespan theorem. v3.2's wording ("curvature and its first derivatives") is not supplied by v3.3 in zones I/II/III or on the Lemma O members, where κ_s is only in L². It is not supplied even by v3.2 near the tip, where Λ controls w̃″ only in L²_ρ. v3.3 supplies the total L² of κ uniformly whenever τ_max < ∞. The authors should confirm the exact form and hypotheses of Chou's criterion (immersed curves).
2. **(H4) range.** Lemma G needs the graph property on y ∈ [14, 700] at τ₁. (H4) says [48, 350] (review M7). The step-4 (H2f) stage implies it on x ∈ [45, 725], but the hand-off list must say so.
3. **Far-family coverage** (review M9, inherited). Step 4 certifies far H4 only at centres 376 and 602. Centres in (490, 602) rest on an unwritten monotonicity argument, and centres above 602 have no certificate. Proposition τ_max = ∞ (bound (ii)) and (W4) beyond 602λ₁ are conditional on this hand-off item.
4. **Constants c_g, C_{S₁} and ℓ_tip** (the curvature conversion; S₁'s conical end). Only their finiteness is used. They are not computed and not quoted as certified.

**Also found in (A), small:**
- `ramp_rep.py` does not assert the m-pattern above. It holds by the probe; add the assert.
- The A_REPAIR R1 claim "c, c′ vanish at t = 1" is false for m = 2. Level 2 is correct by Lemma P(i). This is stated in the text.
- M3: the typed 0.4328 and 0.00116 remain. The typed 2e-5 is now asserted to be at least sup|Dd|.
- Inherited and not re-audited: cos θ constant along the ray in c_s and R; r₀^F; the cflux lemma; the verifier tail lemma (W2).

## 3. B0: slots and how they were filled

Route (b): (A) stays in f − U(α(τ), y) on [14, 700], and d = α(ψ̂₃′ − 1.3056 y^p) + α²ĥ₂′ is inserted. In the texts these places carry [B0]:
- (I1) the ramp J conversion: an η-split per term, with t_d a rate per unit τ;
- (I2) S0V/S1V: + sup_core|d|, + sup_core|Dd|; ×2.119 against 0.0153 at 4·N2_ref;
- (I3) the fluxes: (√K(E) + √K(E_d))²; ×3.382 and ×5.983 at 4·N2_ref and β_A = 2;
- (I4) (H2) for f − U₁ (`h2u_v33`): ×27.3, ×563, ×601;
- the ramp boxes;
- the export: C1 now contains the d-rate; C1⁰ = (7.52e-10, 4.22e-8), C1 = (1.556e-9, 5.431e-8).

**Margins quoted at 4·N2_ref.** I quote the run at 4·N2_ref, which dominates the certified N̄² (review M6). The coordinator's headline numbers at N2_ref [est] (×3.414, ×6.041, S0v ×2.125, M1 ×4.437) appear only as a comparison macro (`\MarFluxRefV`).

**Still open, not part of B0** (B0_REPORT §7): the identification of the verifier's S₁-chart w′ with the ray-graph f − αψ̂₃′ − α²ĥ₂′ at levels ≥ 1.

## 4. Statements of `lemmaU.tex` that need to change

1. **Definition "one frozen profile" (l.60, (C5) remark).** The claim "(H2) … already stated for f − U₁" is now true via `h2u_v33`; cite it. Also say that (A)'s ramp data, fluxes and S0V/S1V carry d.
2. **Lemma D inputs.** (D2) takes the window form D^c + D^s·T; with d the slope is not ρ_J·D^src. Update C0/C1 to the v34 values. Ê is subject to G-A1.
3. **§7.1 standing facts.** The graph property becomes part of the bootstrap: Lemma G of `continuation_v33.tex`. (P1) must be stated in the L² form (Chou / Wheeler).
4. **Proposition "closure".** Its conclusion is on [τ₁, τ_max). Add a pointer to the continuation proposition for τ_max = ∞, and add G to BS_U and to the order table.
5. **B2′ / (A)'s core.** Add the sentence |z^A₁| ≤ √2·S₁ᵛ(Ê) ≤ 0.148 (Lemma A:lem:z1, review M5).
6. **(K6), §"What remains a hypothesis" and the Table.**
   - (H2f) and the fixed-L data are certified (12-stage manifest, 8202d2c6…).
   - (H2) is certified for f − U₁.
   - The K6 margin is ×1.183.
   - All v34 numbers: B̂, late ×2.354/×4.895, present ×2.771/×2.225, U-fine B̂, S₁(2B̂) = 0.28578, cone K3 ×2.511.
7. **Cited sources.** Replace `t2_A_v33/out/A_interface.json`, `U_cert_v33.json` and zone3 `0d6338d0` with the v34 files (04947b7f…, 4c3f23bc…, 8152c650…).
8. **Step 1 reserve argument (l.272)** as review M5: one Young step on √(a² + b²)·N₄.

## 5. Files

All files are in `/home/claude/s1spec/writeup/lemmaU/bulk/fix1/b12_text/`:
- `lemmaS_A_v33.tex`
- `continuation_v33.tex`
- `step3_v33_main.tex`
- `numbers_A_v33.tex`
- `b12_preamble.tex`
- `B12_REPORT.md`
