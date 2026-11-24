# Route decision: closing the level-2 flux gap (fix round 3)

Nov 16, 2026. Fable (executive decision). Read-only: no file other than this one was written, nothing was run except reading the verifier's stored tables.
Labels: **[read]** = verified by reading the cited file:line; **[cert]** = certified number quoted from an existing output; **[est]** = my hand estimate from certified inputs; **[inferred]** = concluded from code and notes, derivation not checked.
All numbers of mine below are [est] unless tagged otherwise. Nothing here is closed.

## 1. Verdict

**Route: (a) + (b)(i), with (i) used in a lighter "trace" form. No verifier rerun, no change to the definition of E_I, no step-4 rerun.**

1. **Pointwise pieces** of the outer-ramp terms (order ≤ 2 in g) are charged to the pointwise slack of the flux constants.
2. **g‴ pieces** (T3 and the (3,·) monomials) go through a **time-integrated absorption lemma** (a), fed by (A)'s existing level-1 window export.
3. **The endpoint (4,0) term** is integrated by parts once, accepting the point term W(57)·(280/27)·g(57)g‴(57). The trace g‴(57) is bounded, integrated in time, by Agmon on the shell from (A)'s z₃ and z₄. This needs a **level-2 dissipation export r₂ > 0 from (A)**, set to r₂ = 1/1000 in the already-approved F2 chain rerun. It needs **no** fifth jets and **no** level-4 chart identity: the existing level ≤ 3 tables on the core suffice [inferred].
4. **Before any of this, an audit (unit G0)**. I believe F1's lemma covers less than the full outer-ramp remainder (§3), and the side note's worry is real (§4). The audit decides the final scope of the flux lemma and whether a verifier-side term is missing.

Fallback if the level-2 window export fails at kill-check G1: route (ii′), the PDE substitution with a corrected energy (§6). Rejected: (iii) and three routes the report did not list (§6).

## 2. Why this route

### 2.1 The room is much larger than the report says
- The verifier consumes only the **sum** F = μ(FLUX0 + ε_I·FLUX2): both constants are added to the same degree-0 group L (`rhs_iv.py:171-172, 184-185`) and to the floor F0 (`closure.py:144`) [read].
- ε_I·FLUX2 = 2.4e-24 against FLUX0 = 1.64e-20. So the level-2 flux is 1.5e-4 of the floor.
- Consequence: the "×1.359" margin of P2 is not a real constraint. The real pointwise condition is K0 + ε_I·K2 ≤ FLUX0 + ε_I·FLUX2, which holds with ×4.2. In level-2 units the pointwise room is about 1.25e-20/3e-5 ≈ 4e-16, roughly 7000 × P2.
- Pointwise slack: s = μ[(FLUX0 − K0) + ε_I(FLUX2 − P2)] ≥ 1.0e-21 (from the F1 certified K0, P2).
- Time-integrated pieces enter with the factor με_I = 2.4e-6: the report's I3 ≈ 2.0e-21 + 7.1e-25·T becomes C′ ≈ 5e-27, S′ ≈ 2e-30.

### 2.2 Why (i) for the endpoint term, and why in trace form
- g⁗ on [57, 60] is tiny only because of the y⁻⁴ scaling that (A)'s variables carry. Every bound that goes through the verifier's own V4 loses that (§6, rejected routes). So the data must come from (A).
- (A)'s level-2 budget is fully spent: core B_c = 169/100, used 0.31, retained 0 (`lemmaS_A_v33.tex:199`) [read]. Reserving r₂ = 1/1000 changes B_c to about 1.689.
- The need is minute. With r₂ = 1/1000 the trace term is about 7e-27·√T in Λ-units, against a binding allowance of 8e-23 (§5). Even r₂ = 1e-6 would do.
- **Trace form avoids the two costly conversions.** A direct bound of ∫∫W(g⁗)² would need the chart identity at level 4 (which is nonlinear in u‴, u⁗) and D⁴d (fifth profile jets). Instead:
  - (4,0) = point term − ∫(c g)′g‴, with c = Wζζ⁗;
  - the integral is of (3,0), (3,1) type and is served by the level-1 export;
  - g‴(57) is expressed through z₃, z₂, z₁ at one point by the existing identity `chart_lemma_body.tex` Lemma (d), k ≤ 3, with the certified core sups of χ^add and β (item (iii) there) and D^k d, k ≤ 3 (`b0_dprof`);
  - sup over the shell of z₃² ≤ (1/ℓ)∫z₃² + 2‖z₃‖‖z₄‖ (ℓ = ln(57.75/54.75)), integrated in time by the exports at k = 1 and k = 2;
  - |g(57)| is pointwise from interface (I2).
- Side remark: fifth jets on [57, 60] are in fact cheap. For x ≥ 28 both profiles are the explicit far-field extension and `c4fun.ext_jets` gives jets to order 8 [read, `agent_int/v3/c4fun.py:6, 55`]. The "[F] fifth jets" note in ASSEMBLY.md §6 concerns the piecewise region. This matters only for the fallback.

### 2.3 Why not wait for a cleverer cutoff
T3 exists for every cutoff, so (a) is needed anyway. Once (a) exists, the endpoint term costs one export parameter and one Agmon step. Changing ζ instead would change the Lyapunov functional itself (§6).

## 3. New finding: F1's lemma is probably not the whole outer-ramp remainder

- F1 replaced v3.1's commutator of the **quasilinear** operator (H = ∂²X(ζg) − ζ∂²Xg with coefficients A_i(V), `lemmaS/lemmaS_cflux.py:3`) by the commutator of the **linear** L1, on the grounds that "the quasilinear and forcing terms enter multiplied by ζ and are bounded by the verifier itself" (FLUX_REPORT §1).
- But the verifier bounds ζ𝒩′ in the cut norms V_j, which requires commuting ζ through. On the inner ramp the remainder is the certified term (a) (`constants.py:35`, `rhs_iv.py:153-164`, S_ramp on [10.25, 16.25]). **I find no outer-ramp counterpart** in the verifier: ASSEMBLY.md §5 lists "inner-ramp couplings of 𝒩" only, and term (c) was the place for everything on [57, 60] [inferred, not verified by derivation].
- If so, the following are covered by neither the verifier nor F1:
  1. the quasilinear part of the commutator, −[(J⁻⁴ − 1)∂⁴ + lower, ζ]g, of relative size ≤ sup|1 − J⁻⁴| ≤ 2.27e-2 of the linear flux (far box), about 5e-4 with |g| ≤ 0.0153;
  2. the modulation transport commutators (§4);
  3. possibly ζ′-terms of the forcing (I checked one by hand: sup(|ζ′| + ψ′ζ) ≈ 1.59 on [57, 60] is below the verifier's 1.67, so that one looks covered [est]).
- Numerically none of this matters (W(57) ≤ 1.44e-18, margins ×4.2). Structurally it is fine too: at level 2 the differentiated coefficients give cubic terms such as g·g″·g‴, which are (sup-bounded) × N2 (pointwise) × N3 (L² in time). The same three channels of §1 take them.
- So the flux lemma must be restated for the **full** operator of the g-equation on [57, 60]. That is unit G2, after the audit G0 fixes the exact list.

## 4. Side note: does the level-2 modulation bound cover the outer ramp?

**No, as coded. It is a further gap. It is formal, and it closes inside the flux lemma at no numerical cost.**

- Level 0 (`rhs_iv.py:147-149`): CAB3985·Am·V0² + CAB25·Am·W1². The W1² term is the inner-ramp commutator (sup q ζ′W/ρ). There is no outer-ramp term [read].
- Level 2 (`rhs_iv.py:178-180`): CAB3985·Am·V2² + 62·Am·V2V3 + CAB03·Amx·α·D2v. **No ramp term at all, for either ramp** [read]. ASSEMBLY.md §6 tags these constants "[F constants]"; I found no written derivation of the form (searched ASSEMBLY, CERT_PLAN, VERIFIER.md, MOD_FORMS, step3_proof.tex).
- Outer ramp, level 0: the missing term is a·q·|ζζ′|W g² (sign of a unknown). F1's remark that the drift sign absorbs it is right in size (|a| ≤ 2.19e-7 against 1/4), **but F1's level-0 bound already uses the full drift sign** for its one-sided diagonal. The same good term cannot be used twice. Fix: put the transport terms into the commutator operator, i.e. replace C by the ball C − a·q − b·(translation coefficient) with |a|, |b| ≤ 2.19103e-7 in `flux_cert`. Effect on κ: about 1e-6 relative.
- Outer ramp, level 2: the commutator −a(qζ′g)″ paired with W v″ gives pairs of order ≤ (2,2) only, pointwise in τ, with coefficient ≤ 4ā ≈ 9e-7 of the drift part of κ⁽²⁾. No new derivative order, no L¹ channel needed.
- **Inner ramp, level 2: not verified, and outside the flux lemma.** If the modulation commutator on [10.25, 16.25] is truly absent, the missing term is about |a|·q·ζ′·(W/ρ)^{1/2}·(W3, B, W1)·‖D2v‖, i.e. ≲ 2e-6 of the certified inner-ramp coupling CF = (0.4421, 0.5259, 0.4235, 0.4250). It could be absorbed in the outward rounding of CF only if the raw certified values leave that reserve. If not, it is a verifier-side term and needs a rerun. G0 must settle this.

## 5. Is (a)'s "no verifier rerun" claim right? Yes, with one correction

The reflected-integral construction is sound. Let φ ≥ 0 be the time-integrated part of the forcing, with ∫_σ^τ φ ≤ C′ + S′(τ − σ) on every window, and S′ < s. Put Z(τ) = 2·sup_{σ≤τ}∫_σ^τ(φ − s) ∈ [0, 2C′]. Then Z′ ≥ 2(φ − s) a.e., and Y = Λ − Z satisfies the verifier's inequality evaluated at the true state, whose energy is Λ ∈ [Y, Y + 2C′].

**Correction: the binding condition is not the cell width (1e-3) but the strictness pad.** At a first touch Y = L(t) with a positive cell rate U, Y′ ≤ 2UΛ ≤ 2U(L + 2C′), while the barrier has L′ = 2(U + 1e-9)L (`closure.py:103`). So one needs U·2C′ < 1e-9·L. At step 0 (U = 0.0303, L ≥ 2.4975e-15, `rate_table_fin12f4_a0.00015.json`) this is **2C′ < 8.2e-23** [est from the stored table]. All ten march steps have U > 0. The estimated 2C′ ≈ 1e-26 leaves about ×8000.

What must be proved and checked (unit G3), for both α points (1.5e-4 and 1.25e-4):

| # | Condition | Current numbers |
|---|---|---|
| C0 | Pointwise slack s certified after **all** pointwise outer-ramp terms are charged (G2) | s ≥ 1.0e-21 before the new terms |
| C1 | S′ < s | 2e-30 against 1e-21 |
| C2 | March, state stays in the cell: L(t) + 2C′ ≤ Ly_hi_k | 2C′ ≤ 2.5e-18 |
| C3 | March, strictness: U_k⁺·2C′ < 1e-9·max(1, |U_k|)·L | **2C′ < 8.2e-23 (binding)** |
| C4 | Tail, anchor cell [L*, 1.01L*] (ETA = 0.01, `run_v33.sh:4`): 2C′ ≤ η·L_F | 2C′ ≤ 9.7e-21 |
| C5 | Tail α-condition with Λ ≤ (1 + 2C′/L_F)ℓ; the α-barrier κ√ℓ itself must stay unchanged, because the anchor rate is certified for α ≤ a* with no cushion | logged margin 177 (`verify_log_fin12f4_a0.00015.txt:22`), trivially fine |
| C6 | Unconditional J's: X_const gains ≤ C′ (relative 1e-12 of 3.74e-15). Must fit in the reserve between the certified κ_J and the rounded 0.042351 used by (A) | reserve about 6e-7 |
| C7 | Every consumer of "Λ ≤ barrier" uses Ly_hi (cushioned) and not L(t): `mod_sups.py:30` does [read]; `farsup`, `sramp_sup`, `sqrtL_max` do (`closure.py:112-113`) [read]; quoted tail floors get the factor (1 + 2C′/L_F) | to be listed exhaustively |
| C8 | Regularity: φ ∈ L¹_loc, Λ absolutely continuous, strict first-touch version | text |
| C9 | **Bootstrap order.** The window bound on φ comes from (A)'s exports, which use the verifier's J's (H-J-W), which now contain C′. Break the loop with a time-integrated hypothesis (BΦ): the window bound holds with (2C′, 2S′) on [τ1, τ′]. The set of such τ′ is closed and contains small τ′; the derived (C′, S′) is strictly smaller. So C2–C6 must hold with the **doubled** constants | still ×4000 at C3 |

The existing L¹_loc lemma of `TAIL_LEMMA.tex` (discounted Z, needs L_F ≥ Z_max/η) covers the tail only and replaces the floor instead of using the slack. Write the new lemma as a separate file; do not edit TAIL_LEMMA.tex.

## 6. Options considered and rejected

| Option | Verdict | Reason |
|---|---|---|
| (ii) / (ii′) PDE substitution g⁗ = J⁴(R̃ − g_τ), giving a corrected energy V2² + ∫Wζζ⁗Ψ(g) | **Fallback** | No change to (A), and fifth jets are cheap on [57, 60]. But it is one heavy bespoke derivation (full g-equation, w_τ through ∂_wJ, modulation and forcing on the ramp) and a less uniform statement for the paper |
| (iii) C⁵ outer cutoff | Reject | Does not remove T3. Changes the definition of E_I, hence the hand-off (H1) stage of step 4, every certified integral that contains ζ, the verifier tag and all texts. A cutoff that is C⁵ at 57 cannot stay below the septic one, so no monotonicity shortcut |
| Bound the (4,0) term by 10.37·N0·V4 and absorb in time through the verifier's J4 | Reject | With N0 ≤ 2.2e-11 and ∫V4² ≤ 6.5e-8 + 3.6e-14·T this gives 2C′ ≈ 1e-17, against 8e-23 (C3) and 2.5e-18 (C2) |
| Same term, Young against the V4 dissipation pointwise | Reject | Needs a fraction ≥ 6e-5 of ½με_I V4² in every rate cell; the cells have no certified reserve there, so it means a verifier rerun |
| Trace of g‴(57) from the left by V3, V4 | Reject | Same loss as the two rows above |
| A-priori jet box on [57, 60] from interior parabolic regularity | Reject | A new explicit-constant PDE lemma of Lemma-O weight, for a term of size 1e-25 |
| Level-4 chart identity plus D⁴d to bound ∫∫W(g⁗)² directly | Not needed | The trace form avoids it |

## 7. Risk to existing margins

- **Verifier (fin12f4): untouched.** No rerun, same tag, same J's up to the rounding reserve (C6).
- **(A) chain, from r₂ = 1/1000** [est; to be confirmed by the v36 run]:
  - core Ĝ₂₂ margin ×1.917 → about ×1.916; Ê₂ up by ≲ 0.1%;
  - ramp M1 ×1.1847: unaffected. It compares the verifier's M1 with the ramp box (`lemmaS_A_v33.tex:83`) and does not depend on (A)'s budgets [read];
  - K6 ×1.1898: unaffected. It is the level-0 datum check, and the levels are triangular 0 → 1 → 2 [inferred];
  - consumers of Ê₂ (S₁ᵛ ×4.4 core, flux) move by ≲ 0.1%.
- **Flux**: K0 ×4.20 and the combined pointwise margin ×4.2 lose at most a few percent to the quasilinear and modulation terms (≤ 2.3% and 1e-6). The separate "P2 ≤ FLUX2" check (×1.359) may drop toward ×1.3; if it ever fails, state the combined condition of §2.1 instead.
- **F2 itself** moves Ê by ≲ 0.1% and c_δU by ×1.143 (U2_TRACE); the flux numbers must be re-transferred from the v36 assembly, not from v35.

## 8. Work units

Models per CLAUDE.md: Opus for certified mathematics, Sonnet for mechanical work, Fable for the final review only. One process, nice 19, < 800 MB. All output in new files; nothing closed is edited.

| Unit | Type | Model, time | Content | Output |
|---|---|---|---|---|
| **G0** ramp audit | certified math, read-only | Opus, 3–4 h | Write out the E_I inequality at levels 0 and 2 on [34.25, 60] as the verifier codes it. List every term where ζ is commuted or integrated by parts: linear, quasilinear, lower nonlinear, forcing, c₀′w, modulation (a, b, β). For each ramp say: covered by which constant / belongs to the flux lemma / uncovered. Derive the form behind 38.5, 0.66, 62, 0.3. **Kill check:** an uncovered inner-ramp term that the rounding of CF, AF, BF cannot absorb means a verifier rerun; stop and report | `fix3/audit/AUDIT_RAMP.md` |
| **G1** level-2 export | certified math, read-only | Opus, 3–4 h | Check that Lemma D ((D2) window form) extends to k = 2: every level-2 source density in `assembly_rep.py` (class Tr) is window-affine, including S33, route IBP and the outer data. Write the statement with r₂. **Kill check:** if a level-2 source is only bounded in the exponentially weighted form, stop; go to fallback (ii′) | `fix3/export2/EXPORT2.md` (+ tex fragment) |
| **F2′** chain rerun | mechanical | Sonnet, F2 + 1–2 h | F2 as approved (forcing constant, c_δU hypothesis), plus: (1) `weights_t2` level-2 retained r₂ = 1/1000 on the core; (2) `assembly_rep`/`export_t2` export for k = 2 and dump of the level-2 window data. t2 → t4 → zone3 into v36 files; two consecutive runs byte-identical | `t2_A_v36/`, U cert v36, zone3 v36, `A_interface.json` v36 with C0, C1, r for k = 0, 1, 2 |
| **G2** flux lemma v2 | certified math | Opus, 5–6 h (split 2a/2b if needed) | 2a: exact reduction for the full g-operator on [57, 60] per G0's list (quasilinear coefficients with the pointwise box, transport with the (a, b) ball); Euler-operator check; no derivative moved across an endpoint except the one accepted point term. 2b: arb coefficient sups; three outputs: pointwise K0, K2^pt; window bounds (C′, S′) for the g‴ class from the level-1 export and for the trace term from the k = 1, 2 exports with the level ≤ 3 chart and d tables; slack s. Develop against v35 numbers, finalise on v36 | `fix3/flux/flux_reduce_v2.py`, `flux_cert_v2.py`, `.json`, `flux_lemma_v36.tex` revised (old one kept as `_pre_v2`) |
| **G3** L¹ consumer | certified math | Opus, 3–4 h | Lemma of §5 with proof (march, tail, α-equation, unconditional J's, bootstrap hypothesis (BΦ)). Checker reading both rate tables, the tail block, `unc`, `hjw_rep.json` and `flux_cert_v2.json`; asserts C0–C7 with doubled constants; args block | `fix3/l1/L1_LEMMA_v36.tex`, `l1_check.py`, `l1_check.json` |
| **G4** text sync | mechanical with care | Sonnet, 2–3 h | v36 copies: Lemma S(A) text (r₂, Lemma D for k ≤ 2, (e) replaced by the flux lemma), continuation text (hypothesis (BΦ), order D1 ⇄ D3), erratum entries (A5; new ones for the cutoff shape, the jet box, the drift double use), numbers file, V33_INDEX successor, U2 items in the paper ISSUES log | new `*_v36.tex` files |
| **G5** review | executive | Fable, one pass | Independent review of G0–G3 from the scripts and texts, not from the reports | `fix3/REVIEW_v36.md` |

Total about 20–25 agent-hours plus the review. Compute is minutes.

## 9. Order relative to F2

1. **G0 and G1 first** (parallel is fine by content, but run them one after the other for RAM). Both are read-only kill checks.
2. **Then F2′**, so that the chain is rerun once. F2 has not started (no v36 files exist) [read].
3. G2a and G3 can be written while F2′ runs; G2b and the G3 checker are finalised on the v36 interface.
4. G4, then G5.

If Micah prefers F2 now: run F2 as approved and accept a second chain rerun (minutes of compute, one more version tag and re-pin).

## 10. What is rerun, what appears

- **Rerun:** the (A) → U → zone3 chain, once, into v36 (already approved as F2; now with r₂). `flux_cert` on the v36 assembly.
- **Not rerun:** the verifier (`step3nl_v33`, tag fin12f4); step 4 (the closed manifest and manifest_v33 are untouched; `b0_dprof.json` and `n1_chart.json` are reused as inputs); Lemma O/R.
- **To check in F2′:** whether any step-4 v33 stage pins the v35 `A_interface.json` as an input. If one does, that is a records question for Micah, not a mathematical one [not verified].
- **New:** `fix3/audit/`, `fix3/export2/`, `fix3/l1/`, v2 files in `fix3/flux/`, the v36 chain directories and texts, later kit_step3 v4 (adds flux_cert_v2 and l1_check).

## 11. Questions for Micah

1. **Approve r₂ = 1/1000** (level-2 retained dissipation in (A)) as part of the F2 rerun? It is a design change to (A) beyond the approved F2 scope. Cost: core Ĝ₂₂ ×1.917 → about ×1.916 [est]; thinnest margins untouched [est].
2. **Hold F2** until G0 and G1 are done (about one day of agent time), or run it now and rerun once more later?
3. **If G0 finds an uncovered inner-ramp term** that rounding cannot absorb: a verifier rerun with a new tag and re-pinning of everything downstream would be needed. I will stop and ask; nothing is started on that path.
4. **Paper wording until G5 passes:** the flux lemma stays listed as open at level 2, and the list should now also name the quasilinear and modulation outer-ramp remainders (§3, §4) as "not yet covered, numerically ≤ 2.3% of a ×4.2 margin".

## 12. Honest limits of this note

- §3 and the inner-ramp part of §4 are inferences from the code and notes; I did not re-derive the verifier's E_I inequality. G0 exists to settle them.
- That Lemma D extends to k = 2, that the level ≤ 3 chart tables suffice for the trace, and that K6 does not move are judgements, not checks (G1, G2, F2′).
- No number in this file is certified by me.
