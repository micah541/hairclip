# Lean n = 1 instance of CWFlow: report (Nov 21, 2026)

**Source.** github.com/micah541/HSLflow-Whitney at commit d72b7447edb0c830331a76152559902c23a4092e, fetched with curl into `src_d72b744/`. That directory holds the CWFlow/ package (76 files) plus the root README, and nothing in it was executed. Hashes are in `SRC_SHA256SUMS.txt`. Main.lean (f0996906…), UniqueMain.lean (12b29d3c…), GraphUnique.lean (2b9d734c…), Pointwise.lean (d4efaf4a…) and Global.lean (3013f2ff…) match the prefixes quoted in U1_regularity.tex §1.2.

**Status.**
- `CWFlowN1.lean` has been drafted, but neither it nor the repo has been compiled.
- Everything below marked [read] comes from reading the source.

## 1. What the repo theorems assume [read]

`exists_flow_of_immersion` (Flow/Main.lean l.29) and `flow_unique_of_immersion` (Flow/UniqueMain.lean l.72) have the same setting:

- **E:** finite-dimensional real normed space, with `[NeZero (finrank ℝ E)]`.
- **M:** charted on E, with `IsManifold 𝓘(ℝ,E) ∞`, `CompactSpace`, `T2Space`, `SigmaCompactSpace`. No boundary, and no orientability assumption.
- **V:** finite-dimensional real inner product space, in the same universe as E.
- **J:** a CLM with `⟪Ju,Jv⟫ = ⟪u,v⟫` and `J(Jv) = −v`.
- **Dimension:** `finrank V = 2·finrank E`.
- **Initial datum** (existence only): `IsLagImmersion E J F`. In every chart the map is C^∞ with injective derivative and `⟪J DF u, DF w⟫ = 0`.

**Hidden dim ≥ 2 assumption: none.**
- No hypothesis bounds `finrank E` from below other than `NeZero`.
- A text search of CWFlow/ finds no `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by` or `opaque`.
- Both theorems, and `exists_flow_of_normal_family`, end with `#guard_msgs in #print axioms … [propext, Classical.choice, Quot.sound]`.
- So if the repo builds, the theorems hold for every E satisfying the stated typeclasses, and in particular for E = ℝ¹. The guard also covers the dependency cone (ShortTime, DifferentialGeometry, Mathlib): a `sorry` anywhere in it would show up as `sorryAx` and break the build.
- At n = 1 the Lagrangian condition is automatic. This is proved in the file (`lag_dim_one`), from hJi and hJJ alone.

**Regularity class (`IsFlow`, GraphUnique.lean l.42-46):**
- Φ is **jointly C^∞ on the closed set [0,T] × M**, including t = 0. So the initial curve must be C^∞, and a solution that is smooth only for t > 0 is not in the class.
- Every Φ(t) is an immersion.
- For every x and every chart, `HasDerivWithinAt (s ↦ Φ s x) (flowVel J (chartRep (Φ t) z) (chart z x)) (Icc 0 T) t`, one-sided at 0 and T. That is, the **velocity is exactly flowVel, with no tangential term**.

**Uniqueness is equality of maps, not "up to reparametrization".** Because the velocity is pinned to the purely normal `flowVel`, there is no gauge freedom: `Φ₁ t = Φ₂ t` as maps M → V.
- Geometric uniqueness modulo diffeomorphisms of S¹ follows from the repo's `exists_flow_of_normal_family` (Flow/NormalFlow.lean l.282). That lemma turns any jointly smooth family whose velocity is flowVel plus a tangent vector into an IsFlow by composing with diffeomorphisms χ_t, where χ_0 = id.
- Combined with uniqueness, this is `cdf_unique_up_to_reparam` in the file, proved with no sorry.

## 2. What the file proves (modulo compiling)

| Item | Lean name | Status |
|---|---|---|
| setting: E1 = EuclideanSpace ℝ (Fin 1), M = `Circle`, V = ℂ, J = I·; hJi, hJJ, hdim; immersion ⇒ Lagrangian | `Jc_isom`, `Jc_JJ`, `hdim1`, `lag_dim_one`, `isLagImmersion_iff` | proved |
| (i) existence | `cdf_exists` | proved from (L1) |
| (ii) uniqueness, exact equality | `cdf_unique` | proved from (L2) |
| (ii′) uniqueness up to reparametrization, for smooth families with tangential term | `cdf_unique_up_to_reparam` | proved from `exists_flow_of_normal_family` + (L2) |
| (iii) flowVel = −k_ss N at n = 1 (stated for X : ℝ → ℂ) | `flowVel_eq_neg_kss` | **sorry (S1)**; exact sympy check exists (U1_sym.py S1) |
| conj-equivariance of flowVel | `flowVel_conj` | **sorry (S2)**; sympy U1_sym.py S2 |
| IsFlow closed under conj∘·∘ι, ι p = p⁻¹ | `IsFlow.conj_reflect` | **sorry (S3)**; plumbing + S2 + repo `flowVel_comp` |
| (iv) symmetry preserved | `cdf_symmetric` | proved from (ii) + S3 |
| restriction | `IsFlow.mono` (repo) | in the repo |
| time shift / gluing (U1 hand arguments) | `IsFlow.shift`, `IsFlow.glue` | **sorry (S4, S5)**; S4 should be easy, S5 is not |

**Choice of M.** `Circle` (the unit circle in ℂ), because Mathlib gives it charts on `EuclideanSpace ℝ (Fin 1)`, a Lie group structure, and the reflection ι = inverse = conj. U1 wrote "the unit sphere of EuclideanSpace ℝ (Fin 2)", which is equivalent mathematically; it makes no difference to the theorems.

**Not covered by the file.**
- The chart-level link between the `IsCurveImmersion` hypothesis and step 3's concrete curve: a 2π-periodic γ : ℝ → ℂ with γ′ ≠ 0, F(e^{iθ}) = γ(θ). The repo's `param_chart` (Flow/LocalParam.lean l.41) is the tool for that link; it was not stated here.
- The maximal-solution lemma.

## 3. Match with what U-1 / C-1 cite

| U-1 / C-1 use | Repo | Verdict |
|---|---|---|
| (L1) existence from a smooth immersion, T > 0, Φ(0) = F, three clauses of IsFlow | Main.lean l.29-36 | **matches** |
| (L2) Φ₁ = Φ₂ on [0,T] if same initial value | UniqueMain.lean l.72-76 | **matches** (exact equality, stronger than "up to reparametrization") |
| no dim ≥ 2 hypothesis | only `NeZero` | **matches** |
| class: jointly C^∞ on [0,T] × S¹, immersions, purely normal velocity | IsFlow | **matches**. U1 §1.2 states it correctly |
| "every immersion is Lagrangian at n = 1" (U1 S4) | — | correct; now a Lean proof (`lag_dim_one`) |
| flowVel = −k_ss N, commutes with conj and with reparametrisations (U1 Lemma "n = 1") | `flowVel_comp` (Scalar/Reparam.lean l.118) covers reparametrisation (any orientation; needs only invertible Dτ) | conj part and the −k_ss N formula are **not in the repo** (S1, S2 here); U1 labels them [sym], which is accurate |
| `IsFlow.mono` (U1 cites GaugeEvolution.lean l.34) | l.34-38 | **matches** |
| gluing, time shift, closure under σ∘·∘ι | not in repo | U1 and C1 already say these are hand arguments; accurate |

**Mismatches and caveats to flag plainly:**
1. **The time-shift wording in C-1.** C-1 Thm (v) says the IsFlow clauses are "invariant under a time shift". That is true mathematically but is not a repo lemma. It is S4 here.
2. **Smoothness up to t = 0 is part of the class.** For (L1) to apply at γ(t₁), the datum γ(t₁) must be C^∞ as a map S¹ → ℂ. That is C-1's job (identifying the step-4 solution with the Lean flow, P1 Gronwall). The point matters in exactly one place, applying (L2) to a solution that is only known to be C^∞ for t > 0. Nothing in C-1 does that: it applies (L2) only to IsFlow objects. So this is OK, but the reliance should be stated explicitly.
3. **(L2) does not identify step 4's object directly.** The P1 Galerkin solution in C⁰H⁵ ∩ C¹L², in a tangential gauge, is not in the IsFlow class. C-1 correctly uses its own Gronwall uniqueness there rather than (L2). The new `cdf_unique_up_to_reparam` would cover tangential gauges, but only for C^∞ families, so it does not replace C-1's argument.
4. **"Uniqueness up to reparametrization"** in CLAUDE.md (Nov 19) and in the review questions. The Lean statement is exact equality within IsFlow. The up-to-reparametrization form follows (ii′), but it is a consequence, not the repo's theorem. The texts should cite (L2) as exact equality; U1/C1 already do.
5. **Lean not built.** Everything rests on the repo compiling at d72b744 and on this file elaborating. Until Micah runs `lake build CWFlow.N1`, the [Lean] labels in U1/C1 should stay "statement read from source, not compiled".
