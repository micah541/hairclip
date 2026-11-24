# Lemma O H6: review of the failure analysis and recommended amendment

Analysis only. No step-3 file and no step-4 production file was touched; everything here is under `step4/p2/h6/fable/`.

- `h6f_cases.py` -> `h6f_cases.json`: Lemma O constants, re-certified by a read-only import of `lemmaO_v3_interval.case` (arb, 200 bits).
- `h6f_margins.py [far_tag] [380] [cone_tag]` -> `h6f_margins_<coneK>_<farM>_<farK>_380.json`: the H6 closure, arb, outward rounding, no sampling (cone centres in cells of width 1). `MU=mu` reruns it with all deviations taken at mu times the horizon.

## 1. Summary

1. **The inner-strip K3 "failure" is an artefact.** The assembly agent compared against E2 <= 0.248, which is the bound for the core y in [28, 60] and for windows that enter later through y = 48. Those never lie in r >= 50 lambda_1. For the strip r in [50, 150) lambda_1, step 3 already asserts a much stronger statement: the (B) windows present at tau_1 are fixed in physical space and keep their hand-off energies up to the factor e^(C 50^-4) (FIXPLAN F5; step3_proof sec. 3.7 and D5; lemmaS H-H6; `lemmaS_Bbox_1e-05_a0.00015_dv0.15.json`, E_B about 1.07e-5). With those per-window sups the strip contributes 0.018 to the cone K3 and 0.050 to K2.
2. **The real H6 failure is the far family.** With far (M, K) = (1, 1.5), C2 sqrt t = 0.35, so the cone cylinders that reach r > 475 recover only x0.54 (K2) and x0.44 (K3). This needs a step-3 change.
3. **The cone family does not need to change.** With the strip taken from the windows, the largest unit ratio in the |x| <= 1 covering is 380/200 = 1.9, not 3, and K2 = K3 = 0.18 closes.
4. **Recommended amendment:** cone unchanged; far re-certified at (M, K) = (0.15, 0.3); H6 proved by strip + covering; six explicit data thresholds added to (H5). Thinnest closure margin x1.69; the binding time margin stays the existing cone t* x2.30.
5. **One apparent step-3 gap found on the way (sec. 6, F-U):** the strip statement is derived for the deviation from U(alpha(tau)), whose forcing grows at fixed physical r. It should be stated for the deviation from the profile frozen at tau_1. H6 depends on it.

## 2. Verdict on each claim of H6_PLAN.md

| claim | verdict |
|---|---|
| Lemma O conclusion holds on \|x\| <= 1 of each member, incl. int g_xxx^2 <= C3 t | Correct (lemmaO_v2, Theorem). |
| Cone K = 0.25: t* x5.69, C2 sqrt t = 0.0545, sqrt(C3 t) = 0.0362 | Reproduced (0.05449, 0.03621, x5.692). |
| Far (0.1, 0.2): t* x16.4, C2 sqrt t = 0.0248 | Reproduced (x16.37, 0.02471). |
| Covering ratio <= 3 | Too pessimistic. Ratio 3 needs r_p = r_c/4 >= 150, i.e. r_c >= 600. Once the strip r_p < 150 is handled separately the ratio is <= 3 r_c/(4*150) <= 1.9. |
| Recovered bound = ratio * (data + C2 sqrt t) | Over-counts the data: the data part is L_c \|k_0(p)\|, which is the cylinder's own datum with no ratio. Only the deviation carries the ratio. It also omits the frame factors (sec. 4), which are 1.01-1.1 and matter at x1.137. |
| "K = 0.18 fails / is too thin", hence raise to 0.25 | Not needed. At K = 0.18: K2 x1.69, K3 x1.93. |
| K3 per member, (L_c/L_m)^3 C3 t, worst 0.1115, x2.23 | Scaling correct. It omits the member data on each piece, the k^2 frame term, and the (50/L_m)^4 gain in C3 t. My value at K = 0.25 is x1.88-1.94. |
| Far family must be re-certified at small (M, K) | Correct, and this is the essential change. |
| Far data meet M = 0.1 at x2.34 (slope); far H4 x1.93 | Correct. Far H4 x1.93 is thin; (0.15, 0.3) gives x4.37. |
| Far inner strip from cone, x1.21 | Correct in size for their constants (my accounting: far K2 x1.21, K3 x1.20 at cone 0.25 / far 0.2). This is the thinnest point of their proposal and is caused by raising the cone K. |
| Inner strip K2 passes via M1core <= 0.1844 (x1.29) | Right conclusion, wrong bound. M1core is a bound on rescaled y in [28, 60]; the physical strip sits at y = r/lambda(tau) >= 50 and leaves that range. The windows give y\|f'\| <= 0.057 there, so K2 strip <= 0.050. |
| Inner strip K3 fails; needs E2 <= 1.0e-3 | Arithmetic reproduced (2.106; 1.05e-3) and the units are right **at tau_1**. Two errors: (i) 0.248 is not the bound step 3 provides for these windows; (ii) a tau-uniform bound on the omega(tau)-normalised E2 is not sufficient anyway. In lambda_1 units int v_rr^2 dr = (lambda_1/lambda) (E2+E1)/(14 r_w^2), which grows like lambda_1/lambda. What is needed is the bound in the weight frozen at entry, which is exactly the physical-window statement. |
| H8 dropped | Agreed (step3_proof l.97, l.105, l.111). |

## 3. Candidates

- **(a) Sharper window bound on [50, 150).** No new estimate is needed: it is step 3's own D5 input ("D4 on the strip [50, 150]"). What is needed is (i) to say so in the H6 proof, with the per-window sups, and (ii) the frozen-profile wording of F-U. **Adopted.**
- **(b) Fixed L = 50, centres >= 300.** Rejected.
  - Time ratio 1.6e-7 (x2.30), self-covering ratio 1, H4 error part about 1.1e-3 M^2: all fine.
  - But the uncovered strip becomes [150, 250). Under the agent's premise (E2 <= 0.248) it still fails: 0.50 in L2 against 0.18. Under the correct premise it is unnecessary.
  - It would also move R_1 and the outer ramp to >= 250, remove the cone cover of the far inner part [150, 263), and force a rerun of every cone check in step 4.
- **(c) Recommended:** below.

## 4. The H6 argument

Units lambda_1; r is the distance along the S1 ray. t ranges over [t_1, T), T - t_1 <= lambda_1^4/(1 - 4 abar).

**Families.**
- Cone: centres r_m in [200, 380], unit L_m = r_m/4, X = 3, cylinder [r_m/4, 7r_m/4], conclusion region [3r_m/4, 5r_m/4]. M = 0.12, K2 = K3 = 0.18 (unchanged).
- Far: centres >= 376, unit 113, X = 2, conclusion region +-113, so the conclusion regions cover r >= 263. **M = 0.15, K2 = K3 = 0.3** (was 1, 1.5).
- Strip: r in [50, 150), the (B) windows present at tau_1.

**Covering of a cone cylinder c.**
- r_p in [r_c/4, 150): windows.
- r_p in [150, 285]: the cone member r_m = 4r_p/3 (L_m = r_p/3 >= 50).
- r_p in (285, 475]: the member r_m = 380.
- r_p in (475, 7r_c/4]: the far member centred at p.

**Covering of a far cylinder.** r_p in [150, 263): cone members as above (ratio <= 113/50). r_p >= 263: far members (ratio 1).

**Scaling.** C_j do not depend on t. For a member of unit L_m the horizon ratio is 50^-4 (50/L_m)^4, so its deviations are those at 50^-4 times (50/L_m)^2. Into the units of c: |h_xx| times L_c/L_m, the L2 norm of h_xxx times (L_c/L_m)^(3/2).

**Frames.** With k the curvature, phi the tangent angle in a frame and s a bound for |tan phi|:
- h_xx = sec^3(phi) k, h_xxx = sec^4(phi) (k_sigma + 3 tan(phi) k^2), dx = cos(phi) d sigma.
- So |h^c_xx| <= (1+s_c^2)^(3/2) |k| with |k| <= |h^m_xx|, and |h^c_xxx| <= (1+s_c^2)^2 (|h^m_xxx| + 3(s_c+s_m) k^2), dx_c <= (1+s_m^2)^(1/2) dx_m.
- From the ray chart: |k| <= |f'|, |k_sigma| <= |f''| + 3|f| f'^2.

**Order of the derivation** (under BS, i.e. (A) in every cylinder up to the current time):
1. Lemma O in every member gives |g_x| <= C1 sqrt t, |g_xx| <= C2 sqrt t, int_{-1}^{1} g_xxx^2 <= C3 t on its conclusion region.
2. Slope. The tangent angle at a covered point changes by <= C1 sqrt t, so the slope in any frame is <= tan(atan s_0 + C1 sqrt t). On the strip: <= tan(atan sup|f| + atan |f(r_c, tau_1)|).
3. K2. L_c|k(p,t)| <= (1+s_0^2)^(3/2) m_2 + (L_c/L_m)(50/L_m)^2 C2 sqrt t; worst at r_p = 150, r_c = 380.
4. K3. Sum over the pieces (strip; members 200, 333.3, 380; one far piece) of the squares of (1+s_c^2)^2 [ (1+s_m^2)^(1/4) (L_c/L_m)^(3/2) (n_3 + sqrt(C3 t)_m) + 3(s_c+s_m) (L_c k)^2 |P|^(1/2) ].
5. All four are strictly below 1/4, K2, K3; the families are compact in the centre, so the continuity argument closes.

**New data thresholds ("H6-data", own cylinder units, whole cylinder), to be added to (H5):**

| | slope s_0 | sup\|h_0''\| m_2 | \|\|h_0'''\|\|_L2 n_3 |
|---|---|---|---|
| cone | 0.08 (data 0.0573, x1.40) | 0.03 (0.0189, x1.59) | 0.004 (1.20e-3, x3.3) |
| far | 0.08 (0.0428, x1.87) | 0.04 (0.0208, x1.93) | 0.004 (7.2e-4, x5.6) |

The closure genuinely needs these: with the data only bounded by M = 0.12 it fails (0.12 + 1.9*0.0395 > 0.18). The old H5 thresholds (0.14 cone; 1.0 far) are replaced by n_3.

## 5. Certified margins (`h6f_margins_0.18_0.15_0.3_380.json`)

| item | bound | against | margin |
|---|---|---|---|
| cone t* (unchanged; grid 14) | 3.683e-7 | 1.6001e-7 | x2.30 |
| far t* at (0.15, 0.3) | 9.511e-8 | 6.134e-9 | x15.5 |
| cone slope (strip is worst; covered points 0.0875) | 0.1354 | 0.25 | x1.85 |
| far slope | 0.0875 | 0.25 | x2.86 |
| **cone K2** (covered, r_c = 380, r_p = 150; strip 0.0497) | 0.1067 | 0.18 | **x1.69** |
| cone K3 (pieces 0.018, 0.080, 0.023, 0.015, 0.034) | 0.0935 | 0.18 | x1.93 |
| far K2 | 0.1313 | 0.3 | x2.29 |
| far K3 (pieces 0.105, 0.023, 0.046, 0.046) | 0.1255 | 0.3 | x2.39 |
| cone H4 (unchanged) | 0.304 M^2 | 1 | x3.28 |
| far H4 at M = 0.15 | 0.229 M^2 | 1 | x4.37 |
| far jets vs M = 0.15: slope, j = 3, j = 4 | 0.0428, 3.6e-4, 0.0121 | 0.15 | x3.5, x418, x12.4 |

Deviations used: cone C1, C2 sqrt t = 0.00748, 0.03955, sqrt(C3 t) = 0.02392 (C0 sqrt t = 0.00227); far 0.00438, 0.03810, 0.02878 (C0 sqrt t = 9.3e-4).

- **Time-equivalent margin of the closure.** With all deviations taken at mu times the horizon, every item still closes at mu = 3.5 (cone K2 x1.04, K3 x1.04). So the closure is not the binding time constraint; t* x2.30 is.
- **Comparison runs** (same script):
  - current far (1, 1.5): cone K2 x0.54, K3 x0.44. Fails.
  - cone 0.18, far (0.1, 0.2): min x1.52 (far K2).
  - assembly proposal, cone 0.25, far (0.1, 0.2): min x1.20 (far K3).
  - cone 0.25, far (0.15, 0.4): min x1.85, cone t* x5.69. This is an optional upgrade if a larger t* margin is wanted; it costs the extra step-3 edits to B4, W4 and the table.

Far physical thresholds for (H5), lambda_1 units: |h_xx| <= 3.53e-4, |h_xxx| <= 1.17e-5, |h_xxxx| <= 1.03e-7, ||h_xxx||_L2(|x| <= 2L) <= 3.32e-6, slope <= 0.08, int (h_0^(5))^2 <= 0.0225 in cylinder units.

## 6. Risks and findings that need the step-3 owner

- **F-U (apparent gap in Lemma S(B), strip).** `lemmaS_Bbox.py` propagates the deviation v = f - U(alpha(tau)) with forcing r0 = 1.4(beta + 1.02 alpha abar) and the forcing integral IB evaluated in entry coordinates.
  - At fixed physical r, d/dtau [alpha lambda^-p] = (beta - p a alpha) lambda^-p, so the forcing grows like e^(0.256 tau). Step 3 certifies only sup|beta| <= 3.34e-9 and a floor for alpha, not decay at that rate.
  - Estimate: the integrated forcing on the window [50, 75] is about 1e-6 e^(0.256 (tau - tau_1)); it reaches the strip sup 0.012 after about 37 units of tau. So "deviation energies <= 2x entry" (B3) is not derivable for all tau as written.
  - The full slope f at fixed r is nearly frozen, so the statement H6 needs is true. Repair: define the strip deviation against the profile frozen at tau_1, U_1(r) = 1.3056 alpha(tau_1) (r/lambda_1)^p. Its forcing is only the b-translation term (about 4 bbar p U/y_E, 2e-10) and the y^-4 residual of U_1. The K formulas in `lemmaS_Bbox.py` already use U at entry coordinates, so the numbers should not move; `lemmaS_Bbox.py` must be rerun with the changed forcing.
  - I have not re-derived Lemma S; this is from reading the script and sec. 3.7. It affects D4/D5 independently of H6.
- **Outer ramp [150, 175].** The windows use B4 there "only through R_1^-4" (lemmaS l.154); `lemmaS_Bbox.py` has no ramp term. Unchanged by this amendment since the cone K stays 0.18, but the strip bound rests on it.
- **Window values are step-3 floats.** The per-window sups are read from the Bbox JSON (inflated by 1e-6). They match the binding statement (K2 <= 0.024, K3 <= 0.025 at L = 50). The strip terms have large room (K2 x3.6; K3 part 0.018 of 0.18).
- **Stale deviation numbers in step 3.** D5 and the Lemma O table quote C2 sqrt t = 0.034 and C0 sqrt t = 0.0021; those are at the withdrawn ratio 1.2e-7. At 1.6e-7 they are 0.0395 and 0.0023. `t9_lemmaO.py` uses C0 = 0.0021 for the cone H1.
- **Thresholds.** Cone slope x1.40 and cone h'' x1.59 against the new thresholds are my choices; they trade against the closure margins.

## 7. Knock-on changes

**Step 3 (needs approval; nothing edited).**
- `writeup/lemmaO/lemmaO_v3.tex`: H5 and H6 texts (sec. 4 above); Revision paragraph, far case (M, K) = (0.15, 0.3), t* >= 9.51e-8, deviations; H7 far margin; needed-time deviations at 1.6e-7.
- `writeup/lemmaO/lemmaO_v3_consts.json` (and `lemmaO_v3_interval.py` main): key `far_margin2` regenerated at (0.15, 0.3).
- `writeup/step3_proof/step3_proof.tex`:
  - (H5) l.97: far constants and physical thresholds; add the six H6-data thresholds; drop 0.14 and 1.0 L^-3/2.
  - (W4) l.105: "resp. 1.5" -> 0.3.
  - Sec. 3.8, l.240-246: far column of the table; H2/H3/H5 far numbers; H6 sentence replaced by the covering statement.
  - B4 l.263: far K = 0.3, M = 0.15.
  - D5 l.311: far t* and factor; the recovery sentence replaced by steps 1-5 of sec. 4; 0.034 -> 0.0395.
  - Sec. 3.7 and B3: frozen-profile wording (F-U).
- `writeup/lemmaS/lemmaS.tex` (H-H6, Part (B)) and `lemmaS_Bbox.py`: F-U wording and rerun.
- Unchanged: the cone family, the verifier, Lemma S(A), all (B) boxes, the time ratios, sec. 5 (smoothness uses only bounded K2, K3).

**Step 4 reruns.**
- `p2/t9_far.py`: new thresholds (h2, h3, h4, slope, L2).
- `p2/t8se/h4loc.py`: far M = 0.15 (ratio 0.229).
- `p2/t9.py`: add the three cone H6-data checks (values already computed: slope, K2, N3).
- `p2/t9_h6.py`: replace by the closure of `h6f_margins.py`.
- `p2/t9_lemmaO.py`: cone C0 = 0.0023 (H1); `t9_lemmaO_far.py` unchanged (C0 = 0.05 covers 9.3e-4).
- `STEP4_RESULT.md` sec. 2 (H5 rows) and sec. 3 (H6 item), `H6_PLAN.md`.
- Not affected: the chain, T7, T8, H1-H4, the cone jets, cone H4.
