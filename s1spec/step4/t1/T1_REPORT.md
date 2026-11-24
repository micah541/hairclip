# T1 report: native float base for P2 (rescaled frame), arb polish, residual feasibility, hand-off

```
cd /home/claude/s1spec/step4/t1                                  # everything below: 1 process each, OMP_NUM_THREADS=1
python3 signs.py                                                  # sympy sign checks (~3 min)
python3 run2.py P11 2.31e4 11 40 0.35 0.08 80 0.12 0.02 1.5 8 0.15  # float base k=20 -> 2.36e4 (950 s) -> slabs_P11.pkl
python3 polish_all.py P11 187 3,5,30,61,90,122,150,183,187        # arb-128 polish of the whole base from the seed + T5 norms (40 min)
python3 resfeas.py P11 6 24.2 200; python3 resfeas.py P11 5 1984 20000   # polished chains started from the float base; mkres.py: table
python3 handoff.py P11 2020 6080.7 10131.6 20000 23000; python3 hotab.py P11   # hand-off quantities, M1 comparison
```
Code: `rep.py` (representation, mesh), `ops.py` (float operator), `radau.py`, `run2.py` (driver, knot insertion), `arbops.py` (arb residual, polish, arb insertion), `resid.py` (T5 norms), `handoff.py` (reuses m1/en.py, m1/jets.py), `signs.py` (output `signs.out`); checks `t_verify.py`, `t_rdirect4.py`, `t_eig3.py`; pre-fix artefacts in `old/`. Seed: even part of HLS `feas2/snaps/p2_k20.npz`.

## 1. Representation and discretisation
- **Unknown** κ̃(x,τ), x ∈ [0, ℓ̃(τ)] arclength from the tip, rescaled frame; θ̃(0) = −π/2, Γ̃(0) = −D_S1 (S1 conventions), κ̃(0) = K0 = 1.5398206378719341 built in. Per element a polynomial of degree P: Hermite-C⁴ nodal DOF (κ̃,…,κ̃'''' at knots, shared) + (P−9) bubbles u⁵(1−u)⁵Leg_i. Exactly C⁴ (node jumps 3e-37 in arb), κ̃' = κ̃''' = 0 at both ends (exact symmetry). Dyadic basis normalisations (identical at both ends of an element: a 1e-12 C⁴ defect was found and removed).
- **Equations.** The κ-form of the rescaled flow is local except one primitive: κ̃_τ|ₓ = ∂ₓΦ, Φ = θ̃_τ = −κ̃''' − κ̃²κ̃' + κ̃Q − (¼−ã)xκ̃, Q = ∫₀ˣκ̃'²; ã = ¼ + (κ̃''''(0) + K0²κ̃''(0))/K0 (keeps κ̃(0) = K0); ℓ̃_τ = −Q(ℓ̃) + (¼−ã)ℓ̃; d ln λ/dτ = ã − ¼; b̃ = κ̃''(0) + (¼−ã)D_S1 (tip at rest, r(0) = 0) enters only r and c̃, not κ̃.
- **Test space (Legendre–Galerkin).** Leg₀..Leg_{P−6} on each element (flux form, exact) + one Galerkin row per knot (the nodal-4th-derivative basis function). Pure Leg₀..Leg_{P−5} moments (the T5 suggestion) have a zig-zag near-null mode: mass-matrix cond 4e5 (uniform) → 8e12 (graded mesh, k = 2e4); the hybrid has cond 1e3–2e4. Leg₀ ∈ test space ⇒ ∫κ̃ = π exactly (turning drift 3e-13 over the float run); I¹..I^{P−5} of the κ-residual vanish at knots at collocation times.
- **Mesh.** 40 fixed elements on [0, 40] (h = 0.26 → 2.5). Outer knots xⱼ(ℓ̃) = 40(1 + Cⱼℓ̃³)^{1/3}, x_N = ℓ̃: Lagrangian far out, at rest at 40; elements split (exact refinement, arb-exact transfer) at slab boundaries when h > 1.5 h_target, h_target ≈ min(0.35 + 0.12x, 0.02ℓ̃). A non-Lagrangian outer mesh is **unstable** (Petrov–Galerkin transport on h ≈ 12 at x ≈ 80: eigenvalues +1.5 ± 20i).
- **Time.** Radau IIA, s = 8 per slab, Δτ = 0.15: κ̃ polynomial of degree 8 in τ per slab (in element coordinates), exactly continuous; knots are explicit functions of ℓ̃(τ). Simplified Newton, float FD Jacobian at mid-slab.
- **Sizes (P = 11).** N = 120 → 174 elements (54 insertions), **n = 840 → 1218** DOF + 4 scalars (ℓ̃, ln λ, c, t); 188 slabs, τ ∈ [0, 28.2], k = 20.0 → 2.36e4. Spectrum checked at k = 193, 1841, 8253, 22716: only the ℓ̃-mode ¼−ã is positive; next −0.24…−0.25; stiffness −2.9e8.

## 2. Runtime (shared 2-core machine)
Float: 5.1 s/slab, 950 s total (P = 9: 3.5 s/slab, 656 s); Newton 8.4 it/slab (max 23); float base vs arb-polished base: 1e-13 (early) to 5e-12 (late) in the normalised unknowns. Arb polish (128 bits, exact polynomial test integrals, float Jacobian inverse, |E| < 1e-26): 0.35 s per 8-stage residual (idle machine), 32 simplified-Newton iterations on average (max 43), 12 s/slab mean for the iterations; the whole base (188 slabs, incl. Jacobian refreshes and 9 norm evaluations) took 2407 s on the shared machine. T5 norms: 0.5 s per sample time.

## 3. Sign checks (sympy, `signs.py`): all confirmed
- e_τ = 𝓛e − r with T5's b₂ = −κ̃², b₁ = −5κ̃κ̃' − (¼−ã)q̃ − b̃t̃₁ + ϖ, b₀ = ¼ − ã − 3κ̃'² − 4κ̃κ̃'' (order ε⁰ vanishes; no κ̃Ṽe term at first order).
- ϖ = −κ̃κ̃' + Q + (¼−ã)(q̃−x) + b̃t̃₁ satisfies ϖ' = κ̃Ṽ; κ̃_τ = Ṽ'' + (κ̃ϖ)' equals the local form above; Φ = θ̃_τ.
- Far end: t̃₁(ℓ̃(τ)) ≡ q̃(ℓ̃(τ)) ≡ 0 give ℓ̃_τ = −ϖ(ℓ̃) and θ̃_τ(ℓ̃) = κ̃ϖ; hence b₁(ℓ̃) = ϖ(ℓ̃) = −ℓ̃_τ and the boundary terms cancel.
- Tip: b̃ and ã as above (S1: 0, 0). Residual identity used for T6: r = Re[e^{iθ̃}∫₀ˣe^{−iθ̃}ρ_θ], ρ_θ = θ̃_τ − Φ; checked against an independent acb reconstruction of Ṽ = ⟨∂_τΓ̃|ₓ, ν̃⟩ (e.g. x = 10: −1.340e-16 vs −1.375e-16 at the nearest sample).

## 4. Residual feasibility (arb-polished, sup over 6 times + 8 nodes per slab)
| slab (k) | start | n | R2 | R1 | Rs | Rᴵs | RF0 | RF1 | RF2 | RF3 | RF4 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **target** | | | 2e-14 | 5e-13 | 5e-11 | 1e-10 | 1e-7 | 1e-7 | 1e-7 | 1e-7 | 1e-5 |
| 0 (20) | seed | 840 | **9.6e-06** | **2.0e-06** | **1.5e-06** | **8.6e-09** | **1.0e-03** | **4.0e-03** | **5.3e-02** | **1.1e+00** | **2.8e+01** |
| 1 (21) | seed | 840 | **6.4e-08** | **1.3e-08** | **9.6e-09** | 5.7e-11 | **7.1e-06** | **2.8e-05** | **3.6e-04** | **7.7e-03** | **2.0e-01** |
| 3 (22) | seed | 840 | **2.0e-12** | 4.4e-13 | 3.1e-13 | 1.8e-15 | 2.5e-10 | 9.9e-10 | 1.3e-08 | **2.7e-07** | 7.0e-06 |
| 5 (24) | seed | 840 | 5.4e-16 | 6.3e-15 | 4.2e-14 | 2.9e-18 | 9.5e-12 | 2.3e-11 | 1.5e-10 | 2.6e-09 | 7.0e-08 |
| 30 (62) | seed | 896 | 1.7e-17 | 4.2e-16 | 3.5e-14 | 2.2e-18 | 1.6e-12 | 3.7e-12 | 2.5e-11 | 4.1e-10 | 1.2e-08 |
| 61 (201) | seed | 973 | 1.5e-17 | 8.2e-17 | 3.2e-14 | 2.1e-18 | 9.8e-13 | 2.4e-12 | 1.7e-11 | 2.7e-10 | 7.5e-09 |
| 122 (1984) | seed | 1092 | 1.5e-17 | 2.5e-18 | 3.1e-14 | 2.1e-18 | 3.0e-13 | 7.4e-13 | 5.5e-12 | 8.3e-11 | 2.4e-09 |
| 183 (19552) | seed | 1218 | 1.5e-17 | 4.4e-19 | 3.1e-14 | 2.1e-18 | 5.1e-13 | 1.4e-12 | 1.6e-11 | 1.9e-10 | 7.0e-09 |
| 187 (22716) | seed | 1218 | 1.5e-17 | 2.9e-19 | 3.1e-14 | 2.1e-18 | 1.3e-13 | 7.2e-13 | 1.1e-11 | 1.9e-10 | 5.5e-09 |
| 59 (186) | float+5 | 966 | 1.5e-17 | 5.9e-17 | 3.2e-14 | 2.1e-18 | 6.6e-13 | 1.8e-12 | 9.5e-10 | **1.4e-06** | **2.9e-03** |
| 122 (1984) | float+4 | 1092 | 1.5e-17 | 4.3e-18 | 3.1e-14 | 2.1e-18 | 3.0e-13 | 7.5e-12 | 1.7e-09 | **9.2e-07** | **1.2e-03** |
| 184 (20299) | float+4 | 1218 | 1.5e-17 | 5.7e-18 | 3.1e-14 | 2.1e-18 | 2.0e-12 | 5.4e-11 | 1.1e-08 | **5.8e-06** | **7.7e-03** |

"seed": the base arb-polished slab by slab from the seed (`polish_all.py`, exact continuity, arb knot insertion); "float+m": a chain of m+1 polished slabs started from the float base. Norms as in T5 (R): R2, Rs with ρ on [0, 34.25]; R1, Rᴵs with ζ∞, W∞ on the whole half curve; RFⱼ = (28∫η²x^{2j−2}(∂^{j−1}r)²)^{1/2} to x = ℓ̃ (η: Lemma-S ramp on y ∈ [14, 28]). Time sampling does not matter: R2, Rs are constant in θ to 2 digits (spatial error only; s = 8, Δτ = 0.15 is ample).
- **Reachable, P = 11 with 40 tip elements (n = 840 → 1218):** past the start-up layer (k ≥ 24) every target holds when the base is polished from the seed: from k ≈ 60 to 2.3e4 R2 ≈ 1.5e-17, R1 ≤ 4e-16, Rs ≈ 3e-14, Rᴵs ≈ 2e-18, RF0–RF2 ≤ 7e-11 (margins ≥ 1000) and RF3, RF4 ≤ 1.2e-9, 3.3e-8 (≥ 80, ≥ 300; at k = 2e4: 1.9e-10, 7.0e-9). Chains started from the float base meet everything except RF3/RF4 at k ≳ 150 (risk 2).
- **P = 9 (n = 600 → 870):** with 40 tip elements R2 = 4e-13 (fails ×20), Rs = 3.6e-11 (×1.4); with 60 tip elements (h = 0.18 → 1.9, n ≈ 715 → 1000) R2 = 5.6e-15 (×3.6), Rs = 9e-13, RF4 = 6.7e-6 (×1.5) at k ≈ 35. **Recommendation: P = 11 / 40 tip elements, ≈ 0.85–1.2k DOF** (T5 estimated 600–900).

## 5. Hand-off (v3 thresholds; U = 1.3056 α y^1.0246 unless marked; base P = 11)
| k | 2020 | 6081 | 10132 | **2.0e4** | **2.3e4** | thr. |
|---|---|---|---|---|---|---|
| α base / M1 | −1.40932e-3 / −1.40930e-3 | −4.55632e-4 / −4.55635e-4 | −2.70044e-4 / −2.70043e-4 | −1.3453e-4 ×1.12 | −1.1658e-4 ×1.29 | 1.5e-4 |
| √Λ base / M1 | 6.10e-7 / 6.10e-7 | 2.25e-8 / 2.26e-8 | 4.90e-9 / 7.74e-9 | 6.5e-10 ×77 | 4.3e-10 ×117 | 5e-8 |
| [14,60] E0,E1,E2 | 1.2e-6, 3.8e-6, 1.0e-5 | 1.8e-8, 6.9e-8, 2.6e-7 | 2.3e-9, 9.4e-9, 3.7e-8 | 1.4e-10 ×1400, 6.0e-10, 2.5e-9 | 8.0e-11 ×2500, 3.4e-10, 1.4e-9 | 2e-7, 1e-4, 0.2 |
| windows max E0 | 9.2e-5 | 1.2e-6 | 2.8e-7 | 2.45e-8 ×8.2 | 1.45e-8 ×13.8 | 2e-7 |
| windows max E1, E2, E3 | 3.2e-3, 7.8e-2, 2.35 | 1.9e-6, 2.2e-6, 8.3e-5 | 7.5e-7, 1.3e-6, 1.6e-6 | 8.6e-8, 2.7e-7, 6.6e-7 ×15 | 5.2e-8, 1.7e-7, 4.6e-7 ×22 | 1e-5 |
| cone K2, K3, K4 (r=200, L=50, X=3) | 0.112, 1.3e-2, 4.0e-3 | 3.4e-2, 3.6e-4, 7.3e-4 | 2.0e-2, 2.3e-4, 4.7e-4 | 1.0e-2, 1.9e-4, 2.4e-4 ×252 | 8.8e-3, 1.7e-4, 2.1e-4 ×291 | K4 0.06 |

- **All v3 items hold at k = 2e4 and 2.3e4.** Binding: |α| (×1.12, ×1.29), then window E0 (×8.2, ×13.8) and E3 (×15, ×22); window maxima sit at the last window [253, 380]. Far-field deviation still −v/(y/k)³ = 40.4–40.6 (U of M1). Chart extent (|f| ≤ ½) to y = 2348 / 2700.
- **Comparison with M1 (U = αψ̃₃' + α²h̃₂', M1's definition)** at k = 2020 / 6081 / 10132: [14,60] E0 9.96e-8 / 1.381e-10 / 6.49e-12 (M1 9.96e-8 / 1.382e-10 / 6.62e-12); sliding windows E0..E3 and K4 agree to ≤ 1% (mostly ≤ 1e-3; E3 at 10132: 3.679e-5 both). A agrees to 0.03% / 0.4% at 2020 / 6081; at 10132 the base gives A = 4.81e-9 = 244|α|³ against M1's 7.68e-9 (390|α|³), i.e. M1's "unexplained excess" was an HLS floor.
- With the v3 U the [14,60] energies are ~1000× those with M1's U (the α²h̃₂ and near-field ψ̃₃ corrections), still ≥ 1400× under threshold at 2e4. Base sensitivity (P = 9 vs P = 11): α ≤ 2e-8, windows and K4 ≤ 5e-7, A ≤ 0.8% relative.

## 6. Risks
1. **Start-up layer.** A seed is not on the discrete slow manifold; the first ~4 slabs (Δτ ≈ 0.6) carry a transient (table rows k = 20.0–23.3). Start the base ≈ 4 slabs before τ_a from a P1 state and use it from τ_a on (T5 assumed τ_a − 0.02), or relax the P1 centre first; T4 must budget base − centre at τ_a.
2. **Far-end strong levels.** RF3, RF4 are dominated by the last ~5% of the half curve (x ≈ ℓ̃, weights x⁴, x⁶ ≈ 1e20): element-scale noise of only 1e-14 there gives RF4 ≈ 3e-3 at k ≈ 180 and does not decay (no diffusion on h ≈ 40; rate ≈ ¼). Short polished chains started from the float base therefore fail RF3/RF4 for k ≳ 150; polished from the seed on (`polish_all.py`) the far end stays clean (RF3 ≤ 1.2e-9, RF4 ≤ 3.3e-8 to k = 2.27e4). T6 must enclose from the seed on, and T5/T8 should confirm that level-4 weights up to x = ℓ̃ are intended.
3. **Scheme stability is numerical, not proven.** Pure moments are ill-conditioned on graded meshes; Petrov–Galerkin is unstable for transport on large non-Lagrangian elements. Both are avoided (hybrid test space; Lagrangian outer mesh with insertions); spectra checked at 4 times only. T6 must handle knots that are explicit functions of ℓ̃ and element splits.
4. **α margin.** 1.12 at k = 2e4: hand off at k₁ = 2.3e4 (×1.29).
5. **Seed closure.** The HLS seed inherits HLS's closure defect (8e-7 physical); corrected by a zero-mean C^∞ bump at physical s = 12 (κ change ≤ 4e-9). P1/T2 data must close exactly. Float closure drift over the run: Im Γ̃(ℓ̃) = 1e-7 rescaled (6e-12 physical), turning 3e-13; T6 must enclose it.
6. **Cost.** Float alone is not enough for T6 (start-up and far-end items need arb). Arb polish of the whole base: 40 min; enclosures (T6) will dominate.
7. **κ-form vs T5's weak form.** I¹..I^{P−5} of the κ-residual vanish at knots; S1, S2 of r do so only up to the coupling (1+K)⁻¹, K = Iκ̃Iκ̃. Immaterial at the polished levels.
8. **U in the shells.** v3 (Lemma S finding 2) uses U = 1.3056 α y^1.0246 everywhere, M1 used αψ̃₃' + α²h̃₂'; both reported.
