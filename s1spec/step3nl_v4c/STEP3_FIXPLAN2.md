# Far-field design, second pass (STEP3_FIXPLAN2.md)

All numbers are FLOAT/design estimates (scripts in `step3nl/fixplan/`), with J = 3× the conditional verifier values as before. Not certified.

## 1. Verdict

The energy-based far field is the right tool. Two things in the current implementation are wrong, and both are repairable without new machinery:

1. **Ramp commutators are taken at order j+1 (or by fourth-power Young).** They can be taken at order j with bounded coefficients. This is what makes E2 ≈ 150; with the correct algebra E2 ≈ 1–5.
2. **M1 is treated as one number.** It must be two: a *ramp* value read off E_I pointwise in time (derived in the closure, not assumed), and a *core* value derived from the shells. With E2 small, the core value is ≈ 0.3, level 2 closes on the whole core [28, 60], and the nested design and its [28, 45] circularity disappear.

The jet boxes M2–M4 are not needed at all (§3).

## 2. Recommended design

- **Hand-off:** α1 = 5e-4, **√Λ1 ≤ 1e-7** (idea (a), adopted as the margin provider, not as the fix). 2e-7 still closes but with margins 1.6–2× instead of ≥ 2×.
- **E_I:** unchanged (blend from 30.25, outer ramp [57, 60]). Idea (b) rejected: it moves the uncontrolled zone outward and costs a re-certification.
- **Shells (A):** one system, levels 0, 1, 2, deviation v = w′, core y ∈ [28, 60], ramp y ∈ [14, 28] (x ∈ [16.25, 30.25]), γ_w = 1. Fed by J0, J1, J2 only. No J3, J4; no nested level 2.
- **(B), Lemma O, F5 geometry:** as in STEP3_FIXPLAN.md (physically fixed windows; y_O = 200; R1 = 175).
- **E_I coefficients on x ≥ 34.25:** S0_v and S1_v/32 from the shells; sup|w‴| removed by the identity ‖v″²‖_W ≤ sup|v′|(3V3 + 1.31V2). The weighted GN route through sup|w′| alone is rejected: the hand-provable constant is 105 (chain 3·3·7·5/3), giving 105·S0_v·d110 ≈ 1.2 against a dissipation of 1/2.

## 3. The two lemmas that carry the design

**L1 (order-j ramp commutators).** For Θ = ωη²e^{−4s} and w = ∂_s^j ṽ:

  ∫Θ w ∂²(a w″) = ∫aΘ w″² + 2∫aΘ″ w w″ − ∫(aΘ″)′ w w′ − ∫a_sΘ′ w′²

(for constant a: ∫Θw″² + 2∫Θ″ww″ − ½∫Θ⁗w²; checked by hand). There is no Θ′w′w″ term. Every remaining w′ on the ramp is removed by one integration by parts against a smooth majorant m ≥ |Θ″|, |Θ′| with m²/Θ bounded (η a quintic smoothstep: η ~ t³, (η²)″ ~ t⁴, so Θ″²/Θ → 0 at the ramp foot). One Young step against δ∫Θw″² then leaves (bounded coefficient)·∫_ramp w². The lower-order parts of y⁴∂_y⁴ = D⁴ − 6D³ + 11D² − 6D reduce the same way. Result: level j needs ramp data of order ≤ j, i.e. J_m with m ≤ j.

**L2 (only M0, M1).** Test ∂_s^jṽ and move two derivatives to the test function. For j = 2, ∂²[a(V)ṽ_ss] = aṽ_ssss + 2a′V_sṽ_sss + (a″V_s² + a′V_ss)ṽ_ss, and V_ss = U_ss + ṽ_ss with U explicit. The only genuinely nonlinear product is ṽ_ss², bounded in L⁴ by ∫Θṽ_ss⁴ = −3∫Θṽ_sṽ_ss²ṽ_sss − ∫Θ′ṽ_sṽ_ss³, i.e. ‖ṽ_ss²‖_Θ ≤ M1(3‖ṽ_sss‖_Θ + ramp term of order 2). The b-terms give ṽ_sṽ_sss and ṽ_ss² only. So the coefficients are M0 = sup|V|, M1 = sup|V_s|, and jets of U up to order 4. **M2, M3, M4 do not appear.** This is the "GN-idealised" column of `level2_attack`, and it should be the only version.

**Ramp nonlinear terms as data ("ε-split", stated as a lemma).** On the ramp, ε(y)·G_j(M1^ramp)·Θ·(∂_s^jṽ)² is not absorbed in the rate; it is bounded by sup(εG_jη²y^{2m−1}/W)·J_m. On the core the rate condition is ε(28)·G_22(M1^core) < γ.

## 4. Dependency graph (no bound used before it is derived)

Inputs fixed before the bootstrap: certified constants; hand-off data H (Λ1, α1, shell energies at τ1, far-field data for (B) and Lemma O).

Bootstrap hypotheses on [τ1, τ*] — exactly these, nothing else:
- **B1** Λ ≤ Λ_cap, |α| ≤ α1.
- **B2** on the core y ∈ [28, 60]: |v| ≤ 0.02, y|v′| ≤ 0.75; shell energies E_j ≤ 2× their derived values.
- **B3** (B) windows on [48, R]: |f| ≤ 1/2, window energies ≤ 2× entry values.
- **B4** Lemma O K-bounds on r ≥ 150λ1.

Derivation order at any τ ≤ τ*:

| step | uses | gives |
|---|---|---|
| D1 closure | B1; B2 as coefficients on x ∈ [34.25, 60]; B2 energies for the absolute (c) constants | Λ ≤ peak; J0, J1, J2; **pointwise on x ∈ [16.25, 34.25]: \|w′\| ≤ 3.2e-3, \|w″\| ≤ 0.042** |
| D2 ramp M-values | D1 only | M0^ramp ≤ 0.03, **M1^ramp = y\|V′\| ≤ 1.2 on y ∈ [14, 28]** (derived, not a hypothesis) |
| D3 shells (A), levels 0→1→2 (triangular) | D1 (J's), D2 (ramp coefficients), B2 (core coefficients M0, M1), B3 (\|V\| ≤ 1/2 on the outer ramp y > 60), H | E0, E1, E2; **S0_v ≤ 0.007, S1_v ≤ 0.25** ⇒ B2 strictly |
| D4 (B) | D3 on [48, 56]; H; B4 at [150, 175] | B3 strictly; curvature ≤ S1_f/r |
| D5 Lemma O | H; D4 on [50, 150] | B4 strictly |
| D6 | D1 | B1 strictly (barrier) |

Every arrow into a hypothesis is strict, so τ* = ∞. M1 enters twice and is never circular: on the ramp it is D2 (from Λ alone); on the core it is hypothesis B2, re-derived in D3 with margin. No time-integrated quantity is used as a coefficient.

## 5. Margins (hand-off 5e-4, 1e-7)

Peak √Λ ≈ 1.45e-7 is an **estimate** (hand-off + ∫5.72·Fn(α) ≈ 2e-8, times the 1.2 growth seen at 2e-7); my verification run failed because the verifier code changed under my patched copy. J's scaled accordingly (design/3.5). Shell numbers from `fixplan/shellA_est2.py`.

| bootstrap quantity | derived | box / limit | margin |
|---|---|---|---|
| peak √Λ | 1.45e-7 | cap 7.5e-6 | 50× |
| sup\|w′\| on [16.25, 34.25] | 3.2e-3 | 1e-2 | 3.1× |
| M1^ramp (y\|V′\|, y ≤ 28) | 1.2 | 2.5 (used in G_ramp = 2.5e6) | 2.1× |
| S0_v (core) | 0.007 | 0.02 | 2.9× |
| M0^core = \|U\| + S0_v | 0.053 | 0.11 (re-evaluate G at this box) | 2.1× |
| M1^core = y\|U′\| + S1_v | 0.05 + 0.25 = 0.30 | 0.75 | 2.5× |
| core rate: G_22(M1 = 0.8) | 3.9e4 | 1.2e5 | 3.1× |
| E_I coefficient sup\|w″\|, x ≥ 34.25 | 0.75/32 = 0.023 (at the box) | 0.094 (tested) | 4× |
| level-1 rate: G_11(M1 = 1.4) | 8.0e3 | ~1e5 | >10× |

Shell energies at this point: E0 = 8.9e-7, E1 = 5.9e-4, E2 = 1.2–1.7 (inflow 0.026, commutators 0.003, ramp nonlinear 0.33–0.56, hand-off 0.2; ×2.14). The ramp nonlinear term dominates E2 and scales like (M1^ramp)⁴·J2.

At hand-off 2e-7 the same table gives M1^ramp 2.0/2.5 (1.25×), S1_v 0.46, M1^core 0.51/0.75 (1.5×): closes, but not at 2×.

**Hand-off requirements that changed:** E0(τ1) ≤ 2e-7, E1(τ1) ≤ 1e-4, **E2(τ1) ≤ 0.2** on y ∈ [14, 60] (was 20). Step 4 must confirm this is attainable; it is the one input I cannot estimate.

## 6. The four ideas

- **(a) √Λ1 = 1e-7:** adopt. Alone it does not fix anything (E2 would still be ~40 with the old commutators), but with L1 it turns 1.5× margins into ≥ 2×.
- **(b) shift the blend to 48:** reject (see §2).
- **(c) parabolic smoothing:** keep as the reserve. Local energies with a time cutoff on unit τ-windows in the polynomial zone give pointwise ‖w″‖, ‖w‴‖_{L²(I)} from J1, J2 (and J3), hence M1^ramp ≈ 0.15–0.25 instead of 1.2. That would cut the dominant E2 term by >100 and restore the 2e-7 hand-off. Cost: medium (a local quasilinear estimate at levels 1–2, plus J3). Not needed if §5 certifies.
- **(d) effective-time view:** it is the reason the design works rather than a new lemma. Scale-invariant jets are frozen along characteristics after y ≈ 20, so the far field only transports what the ramp delivers; (B) with physically fixed windows is exactly the statement "fourth-order flow for effective time y_in⁻⁴". It also explains L1: nothing above order j needs to cross the ramp.

## 7. Tasks

| # | task | size |
|---|---|---|
| T1 | L1: symbolic identity with variable a, majorant m, interval sups of the ramp coefficients for [14, 28] | small–medium |
| T2 | L2: Prop. 1 for levels 0–2 with (M0, M1) only; delete M2–M4 from `lemmaS_prop1_*`; G_jj at boxes (0.11, 0.75) core and (0.05, 2.5) ramp | medium |
| T3 | ε-split lemma and assembly (`run_F1.py`) with two-zone M1; recompute (B) entry E2 and the (c) constants | small |
| T4 | Verifier: output pointwise sI1, sI2 on [16.25, 34.25] as certified bounds (D1 → D2); S1_v box 0.75; production run at (5e-4, 1e-7) with unconditional J0–J2 (V1; J4 no longer needed) | small (+ V1 already dispatched) |
| T5 | Bootstrap write-up following the table in §4 | small–medium |
| R | Reserve: smoothing lemma (c) | medium, only if T2's G_ramp or step 4's E2(τ1) miss |

## 8. Risks

- G_22 on the ramp is extrapolated as M1⁴ from the `level2_attack` table (2.5e5 at 1.4 → 2.5e6 at 2.5). If the certified value is >3× larger, use reserve R.
- L1 is checked by hand only; the majorant step is standard but unwritten.
- The peak at 1e-7 and the J's there are scaled, not run.
- E2(τ1) ≤ 0.2 depends on step 4.
