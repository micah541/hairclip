# F1: flux lemma (outer cutoff of E_I on x ∈ [57, 60]) — report

Fix round 3, unit F1. Nov 16, 2026. All work is in this directory; no existing file was edited. Runs: nice 19, one process, < 150 MB, seconds.
Labels: [cert] = arb with outward rounding; [est] = estimate; [read] = verified by reading the cited file.

## Verdict: NOT CLOSED (level 2). Level 0 is CLOSED.

| input | what the verifier uses | new certified bound | margin |
|---|---|---|---|
| FLUX0 (level 0) | 1.64e-20 per unit τ, pointwise | K0 ≤ **3.9038e-21** at β_A = 2 (with d and chart) [cert] | **×4.201** (was ×3.365) |
| FLUX2 (level 2), part P2 that is pointwise in τ | 8.0e-20 per unit τ, pointwise | P2 ≤ **5.8864e-20** [cert] | **×1.359** (was ×5.956) |
| FLUX2, rest of K2: T3 + O4 | no consumer exists | T3 ≤ 12.488 N3² (N3 integrated in time only); O4 needs ∫∫W(g'''')² near x = 57, which v3.3 does not supply | — |

No verifier rerun is needed for FLUX0 or for P2: both new bounds are below the verifier inputs. The level-2 input FLUX2 is **not** proved, for two structural reasons. Neither is a numerical shortfall: the unbounded pieces are [est] about 1e-21 and 1e-25 against 8e-20.
1. The commutator contains a g'''² term, T3. Its coefficient c33 = W(14ζ'² + 6(W'/W)ζζ' − 8ζζ'') is nonzero for every cutoff. N3 is only bounded integrated in time, from (A)'s level-1 dissipation export. The verifier takes only constants per unit τ (closure.py:144, rhs_iv.py:171–187). The L¹ lemma of TAIL_LEMMA.tex is not implemented (STEP3_REVIEW2.md item 4). So the old "I3" piece had, and still has, no consumer [read].
2. **New finding:** with the actual septic C³ cutoff, ζ''''(57+) = −280/27 ≠ 0. The level-2 form then contains O4 = ∫ Wζζ''''(g g'''' − 4g'g''') + 2W'ζζ''''g g'''. This cannot be integrated by parts below order 4 without the point term W(57)·(280/27)·g(57)g'''(57), so K2 is not bounded on H³(57, 60). It needs local 4th-derivative data (z4) near y = 54.75, and nothing in v3.3 supplies that: (A) exports z2 and z3 only, with r_2 = 0.

## 1. What the verifier consumes (task item 1)
- `constants.py:65`: ABS_FLUX_I0 = FLUX0, ABS_FLUX_I2 = FLUX2. With FIX=1, CMODE = "abs".
  - `rhs_iv.py:171-172` adds μ·FLUX0, and `:184-185` adds μ ε_I·FLUX2, to the degree-0 part L of every rate cell.
  - `closure.py:144`: floor F0 = μ(FLUX0 + ε_I FLUX2).
  - So the verifier needs upper bounds on the outer-cutoff contributions to ½ d/dτ V0² and ½ d/dτ V2² **for every τ**: K0(τ) ≤ FLUX0 and K2(τ) ≤ FLUX2 [read].
  - `grep I3` in verifier/ finds nothing, so there is no time-integrated channel [read].
- The objects themselves (ASSEMBLY.md §5): v_τ = L1 v − [L1, ζ]g + ζ(c0'w + Φ_F' + 𝒩' + M'), with L1 = −∂⁴ − k²∂² − (7km + q/4)∂ + e0 (CERT_PLAN.md:192). The quasilinear and forcing terms enter **multiplied by ζ** and are bounded by the verifier itself. So the flux is the commutator of the *linear* L1:
  - Hc = −[L1, ζ]g;
  - K0 = ∫ W v Hc;
  - K2 = ∫ (W v'')'' Hc, in the H⁴-continuous form: v = ζg is W^{4,∞} but not C⁴, and ½ d/dτ V2² = ∫(Wv'')'' v_τ.
  - **No jet box is needed.** v3.1 modelled the commutator with quasilinear coefficients A_i(V), which is what created M2…M5.
  - **v3.1 also omitted the drift commutator (q/4)ζ'g.** At level 0 it has a good sign (C ≥ 13.68, ζζ' ≤ 0). At level 2 it matters: the drift-off diagnostic `flux_cert_nodrift_DIAG.json` gives κ⁽²⁾ 5–10% smaller.

## 2. Lemma and proof (task item 2)
The statement and proof are in `flux_lemma_v36.tex`. LaTeX is not installed, so the file has not been compiled.
- `flux_reduce.py` (sympy, exact):
  - expands both integrands into monomials;
  - integrates by parts only monomials whose coefficient vanishes identically at 57+ and at 60−, checked with the exact one-sided values of ζ^(i). So no derivative crosses an endpoint; this replaces v3.1's reduce();
  - runs an Euler-operator check that integrand − reduced is a total derivative: PASS at both levels.
  - Level 0 reduces fully to pairs ≤ (2,·).
  - Level 2 leaves (3,3) and exactly three monomials (4,0), (3,1), (3,0). All three carry ζζ'''' and do not vanish at 57+.
- `sanity_quad.py` (float, not part of the proof): a 30-digit quadrature of the defining integrals and the reduced forms agrees to all printed digits, for two test functions with k² ≠ 0 (`sanity_quad.json`).
- Hypotheses now used:
  - (i) the verifier definitions of ζ, W, L1, V_m;
  - (ii) the step-1 profile tail bounds on σ ≥ 54.75: |k| ≤ 3.7e-44, |m| ≤ 8.9e-44, |p| ≤ 8.5e-43, |q − σ| ≤ 2e-86 (`step3/hstail.profile_bounds`);
  - (iii) the pointwise data N0, N1, N2, through the certified v35 pipeline: E^A at β_A = 2, d (b0_dprof) and chart (n1_chart).
  - **No a-priori jet box (M2…M5) is used.**

## 3. Certified constants (task item 3; `flux_cert.py` → `flux_cert.json`)
- Run: arb at 200 bits, 3000 cells, mean-value form for ζ^(i), integer powers by products. Diagonal sups are one-sided, which is valid because only upper bounds are needed. Off-diagonal pairs use weighted AM–GM with fixed rationals: level 0 λ = (3/4, 3/4, 1/2), level 2 λ = (4, 16, 3/2), for the pairs (1,0), (2,0), (2,1).
- Coefficients:
  - κ⁽⁰⁾ = (8.4244, 4.1592, 2.2596); v3.1 had (10.90, 5.27, 2.82).
  - κ⁽²⁾ for levels ≤ 2 = (120.94, 281.14, 439.13); v3.1 had (27.6, 65.8, 113.7).
  - Level-3 diagonal: 12.488.
  - O4 coefficient sups: (4,0) 10.371, (3,1) 41.483, (3,0) 26.606.
- Transfer: the v35 pipeline value is monotone and degree-1 homogeneous in κ (assembly_rep.py:336–362), so K_new ≤ t·K_old.
  - Level 0: t0 = 0.80103, giving K0 ≤ 0.80103 × 4.8735e-21 = 3.9038e-21.
  - Level 2 (P2): t2 = 4.3821, giving P2 ≤ 4.3821 × 1.3433e-20 = 5.8864e-20.
- Rerun check: a second run is byte-identical apart from the args block.

## 4. What is left open
1. **L¹ consumer.** T3 is unavoidable for any cutoff, so the verifier needs a lemma that absorbs a time-integrated flux into the per-unit-τ slack s = μ(FLUX0 − K0) + με_I(FLUX2 − P2) ≥ 0.08 × 1.25e-20 = 1.0e-21 [cert inputs].
   - A reflected-integral (Skorokhod) argument works. The rate cells must have room 2C' ≈ 1e-26 [est] above the barrier, and they do:
     - march cushion Ly_hi/max(L_k, L_k+1) − 1 ≥ 1.0e-3 with L ≥ 2.5e-15 [read, rate_table_fin12f4_a0.00015.json];
     - tail anchor cell [L*, 1.01 L*] (ETA = 0.01 in run_v33.sh).
   - The lemma is sketched in the .tex remark. It still has to be written properly, and the H-J-W J shifts (relative ≈ 1e-12 [est]) still have to be checked. No rerun is needed.
2. **The (4,0) monomial of O4** needs one of the following (all need approval):
   - (i) an (A) export with r_2 > 0 at level 2, giving window bounds on ∫∫ y^{-4}(28/y) z4². This is a t2 → t4 → zone3 rerun with a tiny effect on Ĝ22;
   - (ii) the PDE identity g'''' = J⁴(R̃ − g_τ), giving a total τ-derivative plus order ≤ 3 terms. It needs S0V/S1V (available) and fifth jets of ψ̂3, ĥ2 on [54.75, 57.75] (not certified; the same item is listed as [F] in ASSEMBLY.md §6);
   - (iii) a verifier rerun with a C⁵ outer cutoff, which removes O4 but not T3.
3. **I3 sizes:** with the new coefficients the level-3 parts are about 1.04× the v3.1 "honest" I3 ≈ 2.0e-21 + 7.1e-25 T [est]. O4 (4,0) is about 1e-25 per unit τ if |z4| is of the size of |z2| [est].
4. **Side observation, not part of F1:** the verifier's level-0 modulation bound (ASSEMBLY.md §5) shows a ramp term for the inner ramp only. The outer-ramp part of the a/b drift (|a|q|ζζ'|) is absorbed by the sign of L1's drift (|a| ≤ 2.2e-7 ≪ 1/4). I did not check whether the verifier's level-2 modulation bound covers the outer ramp [not verified].

## Files
- `flux_reduce.py`, `flux_reduce.json`, `flux_forms.pkl`: exact reduction, Euler check.
- `flux_cert.py`, `flux_cert.json`: certified tables and transfer; args block with the sha256 of every input.
- `flux_cert_nodrift_DIAG.json`: diagnostic only (drift switched off; not a bound).
- `sanity_quad.py`, `sanity_quad.json`: float cross-check.
- `flux_lemma_v36.tex`: replacement text for lemmaS_A_v33 (e) and lemmaS.tex §CFLUX.
- `SHA256SUMS`.
