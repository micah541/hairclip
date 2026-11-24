# Phase 2a / A: certified repair of Lemma S(A) (Nov 13, 2026)

New files only. v3.1 (`writeup/lemmaS*`, `step3_proof*`, `step4/`) and the gonogo chains are untouched. `base_ibp/` and `base_v31H2/` are byte-identical copies of `gonogo/lemmaS_A_ibp{,_v31H2}`; the new scripts import their certified functions and do not edit them.

Labels: **[cert]** arb with outward rounding, computed by a script here. **[hyp]** an explicit hypothesis with a placeholder value, not certified. **[rev]** the reviewer's number (CHECKPOINT.md).

## Rerun
`./run_chain.sh` runs four stages in sequence, at `nice 19`, in about 11 s with 93 MB peak RSS:
`rigor/weights_rep.py` → `rigor/ramp_rep.py --N 1200` → `rigor/outer_rep.py` → `rigor/assembly_rep.py`.
- Outputs: `out/{weights_rep,ramp_rep,outer_rep,assembly_rep}.json`. Each has an `args` block (argv, env, input sha256). Logs are `out/*.log`; hashes are in `out/SHA256SUMS`.
- Two consecutive runs gave byte-identical outputs.
- Both (H2) settings are assembled in the same run, in `assembly_rep.json` under `res.tight` and `res.v31`.

## Results [cert; the outer term uses the H-O placeholders, see R5]

| | repaired, tightened (H2) (1e-8, 1e-6, 1e-3) | lemmaS_A_ibp | repaired, v3.1 (H2) (2e-7, 1e-4, 0.2) | lemmaS_A_ibp_v31H2 | v3.1 |
|---|---|---|---|---|---|
| E₀ | **1.9603e-7** | 1.2110e-7 | 3.8603e-7 | 3.1110e-7 | 3.1110e-7 |
| E₁ | **9.2597e-5** | 8.1742e-5 | 1.91597e-4 | 1.8074e-4 | 1.8074e-4 |
| E₂ | **0.069612** | 0.048045 | 0.26861 | 0.24705 | 0.24705 |
| γ′ (lower) | 0.24292, 0.24139, 0.12087 | 0.24629, 0.24617, 0.12355 | same as tightened | same as lemmaS_A_ibp | 0.24629, 0.24617, 0.14306 |
| S0_v / S1_v | 0.004337 / 0.10559 | 0.003718 / 0.09344 | 0.006160 / 0.1769 | — | 0.00575 / 0.171 (page) |
| margin S0v | 4.611 | 5.379 | 3.247 | 3.481 | 3.481 |
| margin M0 core | 6.346 | 6.581 | 5.743 | 5.869 | 5.869 |
| **margin M1 core** | **4.361** | 4.856 | **2.729** | 2.820 | 2.820 |
| margin Ĝ₂₂ core (153664/Ĝ₂₂) | 1.936 (Ĝ₂₂ = 79364.6) | 1.977 | 1.936 | 1.977 | 2.338 (assumed sup) |
| ramp M1 / M0; window M1 / M0 / z₂; \|f\| ≤ ½ | 2.050 / 7.030; 2.052 / 2.544 / 4.268; 5.388 | unchanged | same | unchanged | unchanged |
| K2 / K3 jets | 5.080 / 4.825 | 5.075 / 4.814 | same | same | same |
| flux K0 / K2 (verifier FLUX0 = 4.1e-21, FLUX2 = 2.0e-20) | 4.15e-22 / 1.60e-21 | — | 8.34e-22 / 3.75e-21 (margins ×4.9 / ×5.3) | — | 6.83e-22 / 3.26e-21 (page) |
| I3 const / slope (honest, see R4) | 1.59e-21 / 1.11e-25 | — | 3.19e-21 / 1.12e-25 | — | 1.84e-21 / 1.39e-27 (page) |
| (B) exit y_E = 50: E_B; γ′_B | (1.99e-7, 9.38e-5, 0.0705); 0.23730 | 0.0486; 0.23756 | (3.91e-7, 1.94e-4, 0.272); 0.23730 | 0.2501; 0.23756 | 0.250; 0.2394 (page) |

Parts of the repaired E₂ at the tightened (H2): hand-off 1e-3, inflow 0.01866, ramp k ≤ j 0.02754, **top pair 0.00984**, **cutoff × jet 0.01176**, foot + S33 8.0e-4, outer 1.0e-6, coupling 3.0e-7, forcing 1.7e-10.

**Claims failing.** These are the v3.1 page-1 quotes the repaired numbers exceed. Every structural threshold holds (boxes, S0V/S1V, |f| ≤ ½, Ĝ₂₂ < 153664, K2/K3 ≤ 0.12, FLUX0/FLUX2).
- Tightened (H2):
  - γ′₀ = 0.2429, γ′₁ = 0.2414, γ′₂ = 0.1209, against the quoted 0.246 / 0.246 / 0.143;
  - γ′_B = 0.2373, against the quoted 0.2394;
  - the I3 slope.
- v3.1 (H2), additionally:
  - E = (3.86e-7, 1.92e-4, 0.269), against the quoted (3.12e-7, 1.81e-4, 0.248);
  - E_B, K0, K2 and the I3 constant exceed their page quotes.

**Kill check 2 input.** The repaired E₂ᴬ = **0.0696** at the tightened (H2). That is ×1.44 below 0.1.
- With the outer placeholders ×10 it is 0.0771.
- The late-centre M1 rerun is checkpoint2's job.

## Repairs: the argument and what each one cost

Notation:
- Ramp: y = 14·2^t, η² = S(t)² with S the quintic smoothstep, ε = y⁻⁴, ω = 28/y, Θ = εωη².
- The cutoff part of a level-j pair coefficient is c = ε·co·ω·rem_m(t)·A_i^{(r)} with m ≥ 1, where rem_m·ω = (ωη²)^{(m)} − (−1)^m ωη².
- All suprema are maxima of arb upper ends over a 1200-cell cover of [0,1]. Each cell value is a second-order Taylor enclosure, which is rigorous.
- Data conversions use the v3.1 `convk`: verifier J's, J₃ = √(J₂J₄) + c₁W·J₂, including the ζ < 1 sliver.

### R1. The top pair (j+1, j+1) on the ramp [14, 28]
v3.1 sets this table to zero (`lemmaS_F2.py:81`, `rigor/ramp_arb.py:105`, `rigor/assembly_arb.py:111`). The top-pair monomials are (m,i,r) ∈ {(1,1,0), (1,2,0), (1,2,1), (2,2,0)}, so c contains only A₁, A₂ and A₂′. It is therefore **jet-free**; the script asserts this.

Two routes:
- **Route J:** ∬ c z_{j+1}² ≤ conv_{j+1}(|c|).
- **Route IBP (exact):**
  - Since c and c′ vanish at t = 0 and t = 1: ∫ c z_{j+1}² = −∫ c z_j z_{j+2} + ½∫ c″ z_j².
  - |c z_j z_{j+2}| ≤ δ a_min Θ z_{j+2}² + c²/(4δ a_min Θ) z_j². The weight c²/Θ = (c/εω)² εω/S² is evaluated with c/t³ and S/t³ (both polynomial), so it is finite at the foot.
  - c″ is expanded by Leibniz in A^{(r+u)} with u ≤ 2. Its jet-free part goes to data on z_j². Its jet parts z₂ and z₃ give the classes (j,j,2) and (j,j,3), treated as in R2.

Routes used and costs [cert]:
- Level 0, route J: **+4.96e-8** (reviewer +4.9e-8). Route IBP is impossible here, because c″ carries z₃ and that would leave z₃ linear.
- Level 1, route J: **+1.04e-5** (reviewer +1.04e-5). Route IBP would cost ≥ 3e-5.
- Level 2, route IBP with δ = 0.3178 a_min, inside the certified ramp slack 0.3200 a_min (the reviewer's "0.32 spare" is confirmed):
  - Young 0.00317, c″ jet-free 0.00129, jets of c″ 0.00539 (Φ-class 0.00483, two ibp223 classes 5.6e-4);
  - **total +0.00984**. Route J would cost 3.03 (reviewer 3.0).
- The reviewer's +0.005 to 0.010 leaves out the jets of c″.

### R2. Cutoff × jet classes on (A)'s ramp
Every jet monomial X of A_i^{(r+u)} defines a class: the indices of the pair, together with X, with an explicit weight function. This includes the jets created by the integration-by-parts derivatives, which neither the reviewer nor P2 counted. v3.1's `ramp_level` uses the jet-free Abr[(i, r+u)] only (`rigor/ramp_arb.py:78`).

**Pairs (a, j+1) with a < j now get one integration by parts**, c z_a z_{j+1} → −c′ z_a z_j − c z_{a+1} z_j, instead of v3.1's full reduction (`expand_ab`, `ramp_arb.py:83`, `lemmaS_F2.py:56`).
- With honest jets, the full reduction of (0,3) and (1,3) at level 2 creates z₅ (from A₀⁽⁵⁾ and A₁⁽⁴⁾) and z₄ against non-vanishing weights. No datum controls those.
- The jet-free level-2 tables get smaller: q₀ 511.7 → 2.29 and q₁ 636429 → 77031. q₂ is about the same.

Folding: low factors are folded by sup|z₀| ≤ 1/50 and sup|z₁| ≤ 9/10, keeping two or three factors.

Handlers:
- **pair (x, y), y ≤ 3:** space-time Cauchy–Schwarz against verifier data.
- **pair (x, j+2):** Young against the spare ramp dissipation, when w/S is bounded; otherwise data, whichever is cheaper.
- **(2,2,2):** sup(w/εω)·Φ, where Φ = (εω)_min^{−1/2} √F_max (Jf₂/√ln2 + √2 Jf₂^{3/4} JG^{1/4}) bounds ∬_ramp εω|z₂|³. This is v3.1's sup-free S33 data formula (pointwise Λ). Value: Φ = 3.18e-5.
- **(2,2,3):** the z₃ coefficient is always a *top* jet coefficient f(V) ∈ {a′(V), b(V)}; the script asserts this.
  - Then ∫ κ z₂² z₃ = −⅓∫ κ′ z₂³, with |κ′| ≤ Σ|coef|·εω·[sup|f|·(|D_{q+1}| + 4|D_q|) + sup|f′|·M1·|D_q|].
  - This reduces the class to Φ.

Classes and costs at the tightened (H2) [cert]:
- Level 0: (0,2) via Young, 1.9e-13.
- Level 1: (1,2) via data, 4.56e-7; (1,3) via Young, 1.7e-9.
- Level 2, total **0.01176**:

| class | route | cost | sup weight (per εω) and sources |
|---|---|---|---|
| (2,2,2) | Φ | 0.00777 | 245; from (2,2)×{z₂} 190 (direct pair (2,2), plus the IBP-derived jets of pairs (2,3) and (1,3)), and from (1,2)×{z₂,z₂} 56 |
| (2,3) | data with J₃ | 0.00316 | 82 |
| (2,2) | data | 1.05e-4 | — |
| (2,2,3) | ibp223 | 7.1e-4 | — |
| (2,4) | Young | 1.9e-5 | — |

- The reviewer's "≲ 1.5e-3" is **not confirmed**. The difference is the IBP-derived (2,2,2) jets (weight 190, against the reviewer's 24) and the (2,3) class, which needs z₃ data.

### R3. Profile jets U_j (j ≥ 2) in the quadratic coefficients
Ab is now the Z-free part of the symbolic Aterms. v3.1's `weights_exact.ab_exact` sets V_j = 0 for j ≥ 2 (`weights_exact.py:81–83`); the honest jet-free part of V^{(j)} + U^{(j)} is 2p^jU.
- Core Ĝ_jj: 2280.85 / 2351.49 / 77720.88 → **2281.53 / 2372.05 / 79270.71**. The reviewer's 79271 is confirmed.
- Ramp Ĝ₂₂: 235592 → 239176.
- The windows use the honest U(1.5 y_E) in place of 1/20. Example: window 2025/8 Ĝ₂₂ 33467 → 37268. The window margins are unchanged.
- Cost: γ′₂ −0.0025. E₂ is barely affected.

### R4. Retained dissipation and export
The core interpolation share Bc_j (v3.1: 3/2 a_min) is lowered to **Bc₀ = 11/20 and Bc₁ = 9/20**. The L⁴ usage is 0.2 a_min at level 0 and 0.3 a_min at level 1, so the core keeps **r_j ≥ 1.249 a_min** of 2a(V)εω z_{j+2}² at levels 0 and 1 [cert lower bound; exported as the rational 1249/1000].
- Cost: Ĝ₀₀ 2282 → 4348 and Ĝ₁₁ 2372 → 5288, so γ′₀ 0.24629 → 0.24292 and γ′₁ 0.24617 → 0.24139. The reviewer estimated about 0.238.
- E changes only through the 1/γ′ terms, which are negligible.

**Segment check (new):** on [220, 350], where U ≤ 0.08 > 1/20, with the core box from (B), the inequality ε(220)·Ĝ(U = 0.08) ≤ ε(28)·Ĝ_core holds for every jk, and the same retained fraction holds [cert]. So the core γ′ and the export are valid on all of [28, 350].

**Exported statement [cert]**, for j = 0, 1, every τ_a ≥ τ₁ and every T ≥ 0:

  X_j(T) := ∫_{τ_a}^{τ_a+T} ∫_{ln 28}^{ln 350} y⁻⁴ (28/y) z_{j+2}² ds dτ ≤ C0_j + C1_j·T,

where z_k = ∂_s^k v, v = f − U, s = ln y, y = x − 2.2490776515. Derivation: integrate dE_j + γ′E_j + r_j a_min^c ∫_core εω z_{j+2}² ≤ (data sources) + 2√E·F + coup + outer.
- C0_j = (E_j^A + D_j)/(r_j a_min^c).
- C1_j = (2F_j²/γ′_j + coup_j + outer_j + ρ_J·D_j)/(r_j a_min^c), where ρ_J = max(slope/J) = 1.43e-7 and a_min^c = (1 + 0.11²)⁻².

| | C0₀ (for z₂) | C1₀ | C0₁ (for z₃) | C1₁ |
|---|---|---|---|---|
| tightened (H2) | **2.925e-7** | 5.06e-9 | **1.5104e-4** | 5.85e-9 |
| v3.1 (H2) | 4.483e-7 | 5.06e-9 | 2.3223e-4 | 5.86e-9 |

- C1 is dominated by the outer placeholder rate (R5), so it scales with H-O2.
- Level 2 exports nothing; ε₄ = 1/1000 a_min is reserved for R5.
- **Use in Lemma U:** replace the steady-state X_{j+1} by these time integrals. Lemma U's [28, 56] lies inside [28, 350].
- **The flux constant I₃ is now honest:** it uses r₁ in place of the 2 a_min that v3.1 assumed (`assembly_arb.py:175`), and its level-1 budget includes every new data part.

### R5. Outer cutoff [350, 700]
η_out² = (1 − S)² on t = ln(y/350)/ln 2. The forcing integral now runs to 700; it was 350.

The core identities are taken over [14, 350], where η_out = 1. These are the z_{j+1} interpolation, the L⁴ identity and the S33 integration by parts. Their boundary terms at y = 350:
- **Pointwise boundary terms (R_b, R_e):** bounded by the hypotheses.
- **Inside the L⁴ square root (R_d):** they enter as w·√O_b·N_r with O_b = Θ(350)·M̄₁·M̄₂³. This is bounded by λN_r² + w²O_b/(4λ):
  - λ = 1/10 is added to YC[r] in the core;
  - the S33 N₄ part uses ε₄ = 1/1000 a_min of level-2 dissipation, so Bc₂ = 1499/1000;
  - cost: Ĝ₂₂ 79271 → 79365.

On [350, 700], every term is bounded term by term:
- **R_a:** the cutoff parts of all pairs, including the top pair, with their jets.
- **R_c:** the base parts (weight Θ) of all pairs except (n,n). A factor z_n is absorbed by Young against the outer dissipation (δ_o = ¼ a_min^o); every other factor by its pointwise bound. The outflow term c(η²)′ < 0 has the favourable sign.

**Hypotheses** (not certified; to be supplied by Lemma U zone III, Lemma O or (B)), for τ ≥ τ₁ and y ∈ [350, 700]:
- **H-O1:** |f| ≤ M0o, y|f′| ≤ M1o. Placeholders ½ and 1. |U| ≤ 0.162 is certified.
- **H-O2:** |z_k| ≤ M̄_k for k = 0..3. Placeholders 1, 1, 1, 1.

The rate is stored as an explicit polynomial in M̄ (`outer_rep.json`, `poly`).
- Rates per unit τ at the placeholders: (6.2e-9, 7.1e-9, 1.26e-7). At 10× M̄ they are (6.2e-7, 1.3e-6, 9.1e-4).
- Entering as rate/γ′, they add (2.5e-8, 2.9e-8, 1.0e-6) to E.
- This is visible only in E₀ (13%). The dominant level-0 term is YC_c[1]·Θ(350)·M̄₀M̄₁.

## What remains open or estimated
1. **The H-O1 / H-O2 values are placeholders [hyp].** Everything else in E^A is [cert] given the v3.1 hypotheses (H-B1, H-J, H-R, H-B2, H-B3 and the hand-off).
   - For the v3.1 erratum, some source must supply H-O on [350, 700]. v3.1 never treated the outer cutoff, so this is a new hypothesis for v3.1 as well.
2. **The ramp/core interface conventions of v3.1 are inherited, not re-audited.** This covers how the global L⁴ and interpolation identities split between the ramp box and the core box at y = 28, and the ramp-share term.
3. **Possible savings, not done:**
   - the (2,3) class (0.0032) by interpolation against dissipation instead of J₃ data;
   - (2,2,2) through a sharper sup-free bound than Φ;
   - weighted AM–GM in place of v3.1's ½–½ split for the jet-free mixed pairs.
4. **Not touched:** the step-4 hand-off and the (B) / Lemma O parts beyond what the assembly recomputes. The window loop, the F4 flux, the (B) exit and H6 are recomputed with the new E and Ĝ.

## Further errors found in v3.1 (file:line)
- `lemmaS/lemmaS_F2.py:56` and `rigor/ramp_arb.py:82–83` (`expand_ab`). The full reduction of (0,3) and (1,3) at level 2 requires A^{(r+u)} up to order 7. With honest jets this creates z₅, which is uncontrollable.
  - `rigor/ramp_arb.py:78` bounds all IBP derivatives with jet-free Abr. So the **IBP-derived cutoff × jet classes are an omitted class beyond CHECKPOINT finding 2.** They are the largest part of R2: (2,2,2) has weight 190.
- `rigor/weights_exact.py:24` (UQ = 1/20), used in:
  - `window()` (l.277–287), for windows up to y = 380, where U(380) = 0.086;
  - the core, for y ∈ [220, 350], where U reaches 0.079.
  - Effect: harmless (segment check certified; window margins unchanged), but unproved in v3.1.
- `rigor/assembly_arb.py:175` (I₃). It assumes the full 2 a_min of level-1 dissipation. v3.1 (A) keeps 0.2 a_min (`weights_exact.out`), so I₃ was ×10 too small. The repaired I₃ is honest.
- `rigor/assembly_arb.py:102–104` and the `lemmaS_F1.py` docstring ("outer cutoff far away"): the outer cutoff is never treated, and the forcing integral stops at Rout = 350.
- `lemmaS.tex:209`: the claim "no Young step, only J_k, k ≤ j" is false, because the top pair is dropped. (This is CHECKPOINT finding 1, confirmed.)
- `lemmaS.tex` page 1, v3.1 (H2), after repair: E, γ′, E_B, γ′_B (already known), K0, K2 and I3 exceed their quotes. FLUX0 and FLUX2 still hold at ×4.9 and ×5.3. The erratum needs the new numbers above.

## Files
- `rigor/arlib.py`: shared helpers (args block, honest Ab, U(y)).
- `rigor/weights_rep.py`: R3, R4, the R5 adds and the segment checks.
- `rigor/ramp_rep.py`: R1 and R2.
- `rigor/outer_rep.py`: R5.
- `rigor/assembly_rep.py`: both settings, claims, margins, export and the honest I₃.
- `run_chain.sh`.
- `base_ibp/` and `base_v31H2/`: reference copies.

## Interface for agent B: `export_for_U.json` and Lemma D's c_D
`rigor/export_for_U.py` is the last stage of `run_chain.sh`. It reads the tightened-(H2) assembly and writes `export_for_U.json` (with an args block) in the `U_model2.py --inputs` format. The keys are EA, gp, D, amin, plus kappa, cD_constant_part and D_slope.

**What is certified.** For any τ_ρ ≥ τ₁ and horizon T, measured from τ₁:

  ∫_{τ_ρ}^{τ₁+T} ∫_{[28,56]} εω z_{k+2}² ≤ D_k + D_slope_k·T.

- Start from the level-k inequality, integrated from τ₁, with E_k(τ₁) ≤ H2_k.
- **D_k** = (H2_k + Dsrc_k)/(r_k a_min^c), where Dsrc_k collects the undiscounted, time-integrated data sources: inflow, ramp k ≤ j, top pair, cutoff × jet classes, foot and S33.
- **D_slope_k** = (2F_k²/γ′_k + coup_k + outer_k + ρ_J·Dsrc_k)/(r_k a_min^c).

| | D (for z₂ / z₃ / z₄) | D_slope per unit τ |
|---|---|---|
| tightened (H2) | 1.399e-7 / 7.592e-5 / 71.4 | 5.06e-9 / 5.85e-9 / 1.6e-4 |

- r = (1.249, 1.249, 0.000999), a_min^c = 0.97623.
- **Level 2 is useless.** Only the ε₄ reserve is retained there; nothing is exported at level 2.

**Lemma D's c_D ≈ 1** (U_RAMP.md §2, U_ramp.tex Lemma D) is **confirmed for the horizon-free part** and **refuted as a horizon-free bound.**
- The horizon-free part: (H2_k + Dsrc_k)/E^A_k = 0.870 and 1.000 at levels 0 and 1, so c_D ≤ 1.
  - This holds because E^A = H2 + Dsrc + (sup-form terms), and the sup-form terms are tiny.
  - It needs the integral to start from τ₁, or any later τ_ρ, since it is dominated by the integral from τ₁.
  - If one instead restarts at τ_a with only E(τ_a) ≤ E^A, the factor becomes (E^A + Dsrc)/E^A ≈ 1.82 and 1.99 (`D_from_any_tau_a_with_E_tau_a_le_EA`).
- The horizon-free claim fails because ∫f_k⁺ contains persistent per-unit-τ sources: the β-forcing (2F²/γ′ ≈ 1e-11), the coupling (≈ 1e-12) and, dominating, the R5 outer rate at the H-O placeholders (≈ 6–7e-9).
  - These give D_slope ≠ 0, so ∬_{τ_ρ}^∞ is not finite with these inputs.
  - At level 0 the slope doubles D after about 28 units of τ; at level 1 after about 1.3e4.
  - A horizon-free D needs time-integrable versions of the H-O hypotheses and of the forcing; that is a design decision for Lemma U.
  - **U_model2 ignores D_slope.** Treat its output as conditional on that slope, or add D_slope·T for the horizon used.
