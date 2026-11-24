# H6 (Lemma O: the a-priori bounds close): step-3 amendment plan

This is analysis only. Step 3 and `writeup/lemmaO` are unchanged.
- Scripts: `h6_cases.py` imports `lemmaO_v3_interval.case` read-only and writes `h6_cases.json`; `h6_margins.py` writes `h6_margins.json`.
- Hand-off data: `p2/out/t9_cone_kch.json`, `t9_far.json` and `h4loc_A1_w180.json` (slab-180 window).

## Why the current constants fail
The current families are: cone X = 3, L = r/4, M = 0.12, K₂ = K₃ = 0.18; far X = 2, L = 113, M = 1, K = 1.5. Two checks fail:

| check | failure |
|---|---|
| cone, inner part of a cylinder | even with the |x| ≤ 1 covering below (ratio 3), 3·(0.0189 + 0.0395) = 0.175 against 0.18. That is ×1.03, too thin. K₃ fails too: √37·0.0239 against 0.18 is ×1.2 at best. |
| cone cylinders reaching r > 475 | their bounds can only come from far members: (95/113)·(0.0208 + 0.353) = 0.314 against 0.18 (×0.57) |

Raising the cone K alone does not close it, since C₂ grows with K. At K = 0.25 and ratio 4 the bound is 4·(0.0189 + 0.0545) = 0.294 > 0.25.

## Covering argument (verified)
- Lemma O's conclusion holds on |x| ≤ 1 of each member (lemmaO_v2 l.19).
- Take a point p in cone cylinder c, so r_p ∈ [r_c/4, 7r_c/4]. The member with the largest unit whose |x| ≤ 1 contains p is r′ = 4r_p/3, with L′ = r_p/3. That member exists if 200 ≤ r′ ≤ 380, i.e. r_p ∈ [150, 285]. For r_p ∈ (285, 475], take r′ = 380.
- So L_c/L′ = 3r_c/(4r_p) ≤ 3.
- The bound recovered at p in c-units is ratio·(data |h₀″| + C₂√t). C₂ is taken at the certified M, which bounds it for the smaller data by monotonicity.
- The two edge strips have no cone member:
  - **inner** r_p < 150: covered by step 3's time-uniform core bound;
  - **outer** r_p > 475: covered by the far family.

## Proposed amendment
- **Cone family:** X = 3, M = 0.12, **K₂ = K₃ = 0.25** (was 0.18), with an unchanged L = r/4 and centres in [200, 380]. Re-certified by `lemmaO_v3_interval.case`: t* ≥ 9.107e-7 against the needed 1.6e-7 (×5.69). At the needed time C₀√t ≤ 2.86e-3, C₁√t ≤ 9.51e-3, C₂√t ≤ 0.0545, √(C₃t) ≤ 0.0362.
- **Far family:** X = 2, L = 113, **M = 0.1, K = 0.2** (was M = 1, K = 1.5). t* ≥ 1.004e-7 against 6.13e-9 (×16.4). C₀√t ≤ 6.2e-4, C₁√t ≤ 2.94e-3, C₂√t ≤ 0.0248, √(C₃t) ≤ 0.0180.
- **H6 argument:** the |x| ≤ 1 covering above (ratio ≤ 3). The inner strip uses step 3's M₁^core ≤ 0.1844; the outer strip uses far members.

## Margins (`h6_margins.json`, arb, rounded down)
| item | value against bound | margin |
|---|---|---|
| cone K₂ by covering, ratio 3: 3·(0.0189 + 0.0545) | 0.220 ≤ 0.25 | ×1.137 |
| cone K₃ in L² over [−3, 3], pieces summed in quadrature (Σ ratio³ ≤ 37, conservative): √37·0.0362 + data 1.2e-3 | ≤ 0.25 | ×1.129 |
| cone inner strip from M₁^core = sup y\|f′\| ≤ 0.1844 (L_c ≤ r_p, y ≥ 50) | 0.194 ≤ 0.25 | ×1.289 |
| cone outer strip from far: (95/113)(0.0208 + 0.0248) | 0.038 ≤ 0.25 | ×6.54 |
| far by self-covering (ratio 1) | 0.0455 ≤ 0.2 | ×4.40 |
| far inner strip [150, 376) from the cone family, ratio 113/50 | 0.166 ≤ 0.2 | ×1.21 |
| far data against the new M = 0.1: slope, j = 2, j = 4 | 0.0428, 0.0208, 0.0121 | ×2.34, ×4.81, ×8.27 |
| far H5 (‖h₀‴‖_L² ≤ K − √(2C₃t) = 0.175) | data 7.2e-4 | ×243 |
| cone H5 (≤ K − √(3C₃t) = 0.187) | data 1.21e-3 | ×156 |
| far H4 (∫(h₀⁽⁵⁾)² ≤ M² = 0.01) | 0.0052 | ×1.93 |
| t* | cone ×5.69, far ×16.4 | |

## What would change in step 3 (user approval needed)
- **`writeup/lemmaO/lemmaO_v3.tex` and `lemmaO_v3_consts.json`:**
  - the cone case becomes K₂ = K₃ = 0.25, so t*, C_j and the H5 threshold become 0.25 − √(3C₃t) = 0.187;
  - the far "margin2" case becomes M = 0.1, K = 0.2;
  - the H6 text: the |x| ≤ 1 covering with ratio 3, and the two edge strips.
- **`step3_proof.tex`:**
  - B4: 0.18 → 0.25, and 1.5 → 0.2 with M = 1 → 0.1;
  - (H5) hand-off text: cone ‖h₀‴‖_L² ≤ 0.187, far M = 0.1 and ‖h₀‴‖ ≤ 0.175;
  - l.105: the far-curve a-priori bounds;
  - the Lemma O table: needed time against t*.
- **Unchanged in step 3:** the verifier (H1–H4 constants), Lemma S, and the step-3 core bound M₁^core ≤ 0.1844, which is used and not changed.
- **Step 4:** only the far-jet thresholds change, to h₂ ≤ M/L = 8.8e-4, h₃ ≤ 7.8e-6, h₄ ≤ 6.9e-8, ‖h‴‖_L² ≤ 0.175/113^1.5. The far data pass them with ×4.8, ×2.7e2, ×8.2 and ×243.

## Open items resolved (`h6_k3.py`, `h6_k3.json`)
- **Cone K₃ by covering: verified.** lemmaO_v2's conclusion gives ∫_{−1}^{1} g_xxx² ≤ C₃t per member (g = h − h₀, member units).
  - Partition [−3, 3] of cylinder c into pieces, each inside one member's |x| ≤ 1.
  - A piece in member m contributes (L_c/L_m)³·C₃,m·t, because ∫(h_c‴)²dX_c = L_c³∫(h‴)²dx.
  - Covering greedily with the largest unit, the worst cylinder is r_c = 380, with members 200, 333, 380 and one far member. Its deviation is ≤ 0.1115 in L².
  - So the cone K₃ margin is **×2.23**, excluding the inner strip. This replaces the conservative ×1.13 above.
- **Inner strip K₃ (r < 150): fails with step 3's present bounds.**
  - Step 3's time-uniform core energies E₁ ≤ 1.81e-4 and E₂ ≤ 0.248 give 28∫y²v″²dy ≤ 2(E₂ + E₁) per window.
  - In c-units that allows up to **2.11 in L²**, against the room 0.25 − 0.1115 − 0.0012 = 0.137.
  - Closing it needs a time-uniform **E₂ ≤ 1.0e-3 on the windows meeting y ∈ [50, 150)**, which is ×240 below step 3's 0.248; the hand-off value is E₂ ≤ 3.8e-7. The alternative is another mechanism for the strip. Either one is a step-3 item.
  - The inner strip's K₂ closes (×1.29, from M₁^core ≤ 0.1844).
- **Where step 3 uses the cone a-priori bound 0.18** (read-only grep):
  - `step3_proof.tex` l.105 (far curve: K₂ = K₃ = 0.18 cone, 1.5 far), l.240 (Lemma O table) and l.263 (B4);
  - `lemmaO_v3.tex` l.19 (table), plus the H5 threshold 0.14 derived from it.
  - Not in `lemmaS_consts.json`, the (B) window boxes or the verifier constants. These use S₀^v, S₁^v, M₀ and M₁ and do not use the Lemma O K.
  - The time ratio (l.58, l.311) does not depend on K.
- **H8 dropped.** Step 3 applies Lemma O arc by arc (l.97, l.105, l.111). The axis crossing (`p2/out/t9_axis_crossing.json`) is informational: the curve is immersed after the hand-off, which is expected.

## Open points inside this plan
- *(superseded by the exact covering above)* **Cone K₃ quadrature.** It assumes the per-member L² deviation on |x| ≤ 1 is √(C₃t) (as in H5's √(⌈X⌉C₃t)). Σ ratio³ ≤ 37 counts at most 3 members on the inner piece (ratios 3, 1.8, 1.08) plus at most 2 on the rest. A finer count only raises the margin.
- *(done above: it fails with E₂ ≤ 0.248 and needs E₂ ≤ 1.0e-3 on [50, 150))*
- *(withdrawn: the crossing is informational; H8 is not required)*
