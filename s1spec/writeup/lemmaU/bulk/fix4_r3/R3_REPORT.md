# R3 correction (tag v36pre): honest forcing constant of Lemma S(A) and Lemma U's c_dU, trial run

Nov 19, 2026. `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `N` = `B/fix4_r3`. Labels: [cert] = python-flint arb with outward rounding, all constants exact rationals; [est] = not a bound; [float] = a float quotient of certified values, not outward-rounded (as in v35).

## Result in one paragraph
The R3 correction is in place and certified in new files only, and the whole t2 -> t4 -> zone III chain was run with it in `N`. **No margin fails.** The largest relative loss on any margin is x0.99944 (flux FLUX0 margin at beta = 2, probe hand-off setting). The smallest margins are unchanged: ramp M1 x1.18473, ramp z1 x1.20143, Ghat22 x1.91741 and H6 far K3 x1.86087. The forcing norms F_j go up by x1.0696, x1.0746 and x1.0892 for j = 0, 1, 2 [cert], which matches the U2_TRACE estimate. Ehat_j goes up by at most x1.000394 (level 0) [cert]. c_dU goes up by x1.14336 to 1.96498e-8 [cert]. The thinnest new item is that **c_dU is now only x1.0178 below the typed bound 2e-8** [float of cert values; `c_dU < 2e-8` itself is asserted in arb]. u1_ghat's M1z check and zone3's assertion both still use that typed bound; before R3 it had x1.164 of room. Two runs of the R3 chain give byte-identical JSON outputs (all 32). Every file in the proof tree is unchanged (sha256, section 6). **This chain still consumes the fin12f4 verifier,** whose E_I inequality failed audit (`B/fix3/audit/AUDIT_RAMP.md`). The numbers below are therefore a preparation run; they are not a closure of step 3.

## 1. What was done (files)
- `N/t2_A/`: a copy of `B/fix2/n12/t2_A_v35/` (base_ibp, base_v31H2, prov, ref_phase2a, rigor, run_t2.sh). Its `out*` and probe directories were not copied; `out/` is fresh.
- `N/t4_cert/`: a copy of the `B/t4_cert/` scripts. The `run_all*.sh` files were removed from the copy because they `cd` into the original directory.
- `N/t5_zone3/`: a copy of the `B/t5_zone3/` scripts. `out/ramp_U.json` is a byte copy of the v35 one: an input, not rerun (see section 7).
- **Path-only edits**, needed so that nothing writes outside `N`. The diffs are in `N/patches/paths_*.diff`:
  - `t2_A/rigor/arlib.py` ROOT;
  - `t2_A/run_t2.sh` cd;
  - `t4_cert/t4lib.py` T2/T5;
  - `t5_zone3/zone3.py`: the `--ucert_full` default is now None, which forces an explicit path. The old default pointed at the v34 U_cert.
- **R3 edits.** The pre-R3 backups are `*_pre_r3.py`, each byte-identical to the v35 script (checked by `cmp`). The diffs are in `N/patches/r3_*.diff`; reapply them on top of the re-derived verifier chain.
  - `t2_A/rigor/assembly_rep.py`:
    - r0 is now r0';
    - F_j^2 = (r0' p^j sqrt(I_int) + sqrt(Feps_j^2))^2;
    - c_F in `f4b_honest` uses the honest F_1;
    - each setting records `r3 = {r0, FE2, F2}`;
    - the R3 module is added to the args inputs.
  - `t2_A/base_ibp/rigor/assembly_arb.py` (the copy; imported by assembly_rep): the honest forcing goes into three places, each with its own eps-part on its own y-range and weight:
    - `window_loop` (the v3.1 (B) windows, j = 0..3);
    - `f4b` c_F (I_3 budget, [14, 350]);
    - the `f4b` (B)-exit (y_E = 48, 50, up to 350).
  - `t4_cert/u34_sweeps.py`: c_dU now uses betabar' = betabar + DL alpha_1, and the R3 module is added to the args inputs.
- **New:**
  - `N/r3/r3_forcing.py`: the certified constants, importable, with the derivation in its docstring. It writes `r3_forcing.json` with an args block.
  - `N/r3/r3_symcheck.py`: the exact sympy checks of the derivation (all residuals 0). It writes `r3_symcheck.json`.
  - `N/run_chain.sh`: the driver.
  - `N/cmpjson.py`, `N/r3/compare_r3.py`, `N/r3/key_table.py`: read-only comparison tools.
- **Runs**, all in `N/runs/`, with `stages.txt` giving wall time and max RSS per stage:
  - `base`: the paths-only copy before the R3 edit;
  - `r3run1`, `r3run2`: the R3 chain;
  - `failed_r3run1_symcheck_sign`: the first attempt. It stopped at stage 0 because my symbolic check had the wrong sign on the b-term target (+ c_U p b cos(theta) alpha y^(p-1), not -). The sign is irrelevant to the bound, which uses absolute values. I fixed the check and the docstring; the record is kept.

## 2. Derivation of the corrected constants
These are exact; the same text is in the `N/r3/r3_forcing.py` docstring. Notation:
- (A)'s variable is u = f - U, with U = c_U alpha(tau) y^p, c_U = 1.3056 and p = 1.0246, both exact.
- s = ln y and eps = y^-4.
- The slope equation (`lemmaS.tex` section 1; `lemmaS/lemmaS_symbolic.py`) is V_tau = -(1/4 - a) y V' + b cos(theta) V' - d_y^2 (A(V) V'' + B(V) V'^2), with A(x) = (1+x^2)^-2 and B(x) = -3x(1+x^2)^-3.

Since d_y^2 (y^-2 G) = y^-4 (d_s^2 - 5 d_s + 6) G, in s this reads

  V_tau + c V_s = -eps (d^2 - 5d + 6) G(V),  G(V) = A(V)(V_ss - V_s) + B(V) V_s^2,  c = 1/4 - a - b cos(theta) e^-s.

With V = U + u, G(V) = G(U) + A_2 u_ss + A_1 u_s + A_0 u exactly, so the forcing of (A) is **R = -U_tau - c U_s - eps (d^2 - 5d + 6) G(U)**.

**(i) Non-eps part.** U_s = pU, and alpha_tau = lambda~ alpha + beta with |beta| <= betabar (D1). Then

  -U_tau - c U_s = -c_U y^p [beta + (lambda~ + p/4) alpha - p a alpha] + c_U p b cos(theta) alpha y^(p-1)   (r3_symcheck item 3: residual 0).

- d_s^j y^p = p^j y^p.
- |d_s^j y^(p-1)| = |p-1|^j y^(p-1) <= p^j y^p / 14 on y >= 14. Every range used is >= 14.

So |d_s^j (non-eps part)| <= r0' p^j y^p, with

  **r0' = c_U (betabar + DL alpha_1) + c_U p alpha_1 (abar + bbar/14)**,  DL = max(|lambda~ + p/4|, |LAMT + p/4|) = 3.2224577463e-6 (exact).

The constants in DL:
- lambda~ = -0.2561532224577463 (`step3nl/h2md.py:13`, `psi3beta.py:5`);
- LAMT = -0.2561532 (`step3nl_v33/verifier/constants.py:14`);
- taking the maximum covers either reading of beta.

The inputs are the ones the assembly reads (`base_ibp/lemmaS_F10p_M0.52.json` args, one ulp up): alpha_1 = 1.5e-4, abar = bbar = 2.19103e-7 and betabar = 3.33811e-9. The parts of r0' [cert]:

| part | value |
|---|---|
| c_U betabar | 4.3582e-9 |
| c_U DL alpha_1 | 6.3109e-10 |
| c_U p alpha_1 abar | 4.3965e-11 |
| c_U p alpha_1 bbar/14 | 3.1403e-12 |
| **r0'** | **5.03643e-9** (old 1.4(betabar + 1.02 alpha_1 abar) = 4.72029e-9; ratio x1.066975) |

This also covers U2_TRACE R3.2 item (iii), the a/b coefficient for y < 14.8.

**(ii) Eps part.** On U = c_U alpha y^p, d_s acts on functions of U as pD, where D = U d/dU. Also d_s^j (eps h(U)) = eps (pD - 4)^j h(U). Hence

  d_s^j R_eps = -y^-4 Phi_j(U),  Phi_j = (pD - 4)^j (p^2 D^2 - 5pD + 6) g,  g(U) = (p^2 - p) U (1+U^2)^-2 - 3 p^2 U^3 (1+U^2)^-3   (r3_symcheck items 1, 2: residual 0 for j = 0..3).

Phi_j is odd. h_j = Phi_j / U is an exact rational function of w = U^2: sympy computes its numerator and denominator with rational coefficients, and arb evaluates them on the ball [0, Um^2], with the denominator certified > 0.

- **Linear coefficient** (exact): c_eps,j = c_U |p-4|^j (p^2-5p+6)(p^2-p) = 0.063407, 0.18866, 0.56134, 1.6702 for j = 0..3. This is the U2_TRACE value 0.0634 at j = 0.
- **Cubic coefficient** of Phi_j (exact, record only): g_1 (3p-4)^j (9p^2-15p+6), with g_1 = -2(p^2-p) - 3p^2. Its value over c_U at j = 0 is -0.25357.
- **No series is used for the bound.** On each geometric y-cell [ya, yb], |U| <= Um = c_U alpha_1 yb^p, so |d_s^j R_eps| <= c_U alpha_1 y^(p-4) k_cell with k_cell = sup_{0 <= w <= Um^2} |h_j(w)| [cert].
- At y = 700, where Um = 0.162, the cubic terms are not small. The sup of |h_j| there is 0.0906, 0.173, 0.568, 1.738 (j = 0..3), against linear parts 0.0486, 0.1445, 0.4300, 1.2793 (= c_eps,j / c_U). They are covered because k_cell is a sup over the whole cell.

**(iii) L^2 forcing.** F_j^2 = int omega |d_s^j R|^2 ds, with omega = (y_ref/y)^gamma, gamma = 1 (asserted), and ds = dy/y. By Minkowski,

  **F_j <= r0' p^j sqrt(I_int) + sqrt(Feps_j^2)**,  I_int = y_ref int y^(2p-2) dy,  Feps_j^2 = sum_cells (c_U alpha_1 k_cell)^2 y_ref int_cell y^(2p-10) dy

The cell integrals are taken over [ya.lower, yb.upper]. The ranges and cell counts are:
- (A): [14, 700], y_ref = 28, 400 cells;
- windows: [y_E, 1.5 y_E], y_ref = y_E;
- (B) exit: [y_E, 350];
- f4b c_F: [14, 350].
The last three use 100 cells.

For (A) on [14, 700] [cert]:

| j | Feps_j^2 | F_j old | F_j new | new/old |
|---|---|---|---|---|
| 0 | 3.936e-18 | 7.51699e-7 | 8.04028e-7 | x1.06961 |
| 1 | 3.485e-17 | 7.70191e-7 | 8.27678e-7 | x1.07464 |
| 2 | 3.085e-16 | 7.89138e-7 | 8.59556e-7 | x1.08923 |

The old F_j^2 reproduce `numbers_A_v33.tex:40` (5.6506e-13, 5.9320e-13).

**(iv) Lemma U, delta-U lemma (LU:55-60).** With the hypothesis replaced by alpha_tau = lambda~ alpha + beta,

  d/dtau[alpha lambda^-p] = (beta + (lambda~ + p/4) alpha - p a alpha) lambda^-p,

so **c_dU = c_U (betabar + DL alpha_1 + p abar alpha_1) / (p (1/4 - abar))**, with betabar' = betabar + DL alpha_1 = 3.82148e-9.
- Inputs: the verifier strings sup_ab = 2.19103e-7 and sup_beta = 3.33811e-9 (t4lib.verifier_sups).
- Result: c_dU = 1.96498e-8, against the old 1.71860e-8 (x1.14336) [cert]. `c_dU < 2e-8` is asserted in arb in both `u34_sweeps.py` and `r3_forcing.py`.
- c_dU enters Delta_in (x1.30727, which is (x1.14336)^2), dU_j, Lemma N's Nbar, and zone III's n_III.

## 3. Where the old constant was used, and what is now used
| site | old | new |
|---|---|---|
| `assembly_rep.assemble` (E_j forcing; export C1 = 2F^2/gp + ...) | (1.4(betabar + 1.02 alpha_1 abar) p^j)^2 I_int | (r0' p^j sqrt(I_int) + sqrt(Feps_j^2))^2 |
| `assembly_rep.f4b_honest` c_F (I_3 budget) | 1.4(...) p sqrt(I_int)/sqrt(Lambda) | F_1(honest)/sqrt(Lambda) |
| `assembly_arb.f4b` c_F, (B)-exit EB_j | the same 1.4 form ([14, 350], [y_E, 350]) | honest, with Feps on the same ranges |
| `assembly_arb.window_loop` EB_j (v3.1 (B) windows) | 1.4 form | honest, with Feps per window |
| `u34_sweeps` c_dU | betabar | betabar + DL alpha_1 |
| U2's "P32:235-236" | the line `d/dtau[alpha lambda^-p] = (beta - p a alpha) lambda^-p (lambda~ = -p/4)` in `writeup/step3_proof_v32/step3_proof_v32.tex` (Prop. B' discussion; a closed text) | not edited; the same identity, corrected, is item (iv) |

No other use of lambda~ = -p/4 was found in the chain scripts (searched `p/4`, `1.4'`, `1.02`, `r0`). Lemma U's own forcing R^ (frozen U_1) has no alpha_tau, so it is unaffected; it was not re-audited.

## 4. Runs and reproducibility
| check | result |
|---|---|
| `base` (paths-only copy) against the pinned v35 outputs (`cmpjson.py`, args/timing ignored) | `A_interface.json`, hjw, ic, a3, weights, ramp, outer, u2_K, u34 **identical**. The asm_* files, u1_ghat, U_cert passes and zone3 differ **only in stored path strings** (outer file, ramp_U src, A_interface path, Bhat_source, note). The copy reproduces v35 exactly. |
| `r3run1` vs `r3run2` | **all 32 JSON outputs byte-identical** (sha256), including args. The only differences are 5 logs, in timing lines. |
| R3 vs base, byte level | unchanged: a3_signs, hjw, ic_*, outer_L2, ramp_rep*, weights_t2*, u1_ghat, u2_K, ramp_U. Changed: asm_* (7), A_interface, u34_primary, U_cert passes, zone3 interfaces. So the R3 change is isolated to the intended stages. |
| Lemma U fixed point | `U_cert_r3_pass3` == `U_cert_r3_pass2` apart from args (as in v35). The zone3 pass-1 and final interfaces differ at ~1e-10 relative, the same pattern as v35. |
| N35 assertions (not inside the pinned scripts; values reported) | ramp M1 x1.18473, ramp z1 x1.20143, ramp M0 x7.02 (> 1; unchanged). Nbar^2/(4 N2ref) <= 0.41264 (<= 1). |
| failing legacy claims in `asm_*.log` (v3.1 page-1 comparisons) | the same set of names as in base/v35; they are legacy and not consumed |
| resources | 13 stages, single process, nice 19, PYTHONHASHSEED=0, no .pyc. Max RSS 92.8 MB. Wall time about 2.5 min per chain: t2 48 s, u34 46 s, u2_K 49 s. |

## 5. Old (v35) vs new (R3)
Old = the `base` run, whose numbers equal v35. New = `r3run1`. Values are arb midpoints (radii < 1e-20 relative). Margins are certified lower bounds where the source exports them as such. The table was generated by `N/r3/key_table.py runs/base runs/r3run1` (`N/r3/key_table_base_vs_r3.md`).

The **complete** leaf-by-leaf comparison of all 16 chain outputs is in `N/r3/compare_base_vs_r3.md`, generated by `r3/compare_r3.py`. It has two parts: every margin leaf (954 rows, changed or not) and every other changed number (2534 rows). Of all margin leaves, the largest relative decrease is x0.999442.

Labels:
- zone III `at_fine_est` rows are [est] (as in v35);
- zone III `box_380` margins are float quotients 0.15/b0, 0.52/b1, 0.16/bz of certified bounds [float];
- H6 margins are arb lower bounds (`LO`) of certified quantities;
- the U_cert `margins` are `flo` of arb [cert].

One extra row, not in the generated table: **c_dU against the typed bound 2e-8** (used by u1_ghat's M1z check and asserted by zone3) goes from x1.16373 to x1.01782 [float of cert; the inequality is asserted in arb].

| quantity | old (v35) | new (R3) | new/old | label |
|---|---|---|---|---|
| F_j^2 (A, [14,700]) j=0 | 5.65052e-13 | 6.46462e-13 | 1.144075 | [cert] |
| F_j^2 (A, [14,700]) j=1 | 5.93194e-13 | 6.85051e-13 | 1.154851 | [cert] |
| F_j^2 (A, [14,700]) j=2 | 6.22738e-13 | 7.38836e-13 | 1.186431 | [cert] |
| forcing share of E_j, j=0 | 6.78362e-09 | 7.25729e-09 | 1.069825 | [cert] |
| forcing share of E_j, j=1 | 6.52847e-08 | 7.01592e-08 | 1.074665 | [cert] |
| forcing share of E_j, j=2 | 3.14157e-06 | 3.42192e-06 | 1.089237 | [cert] |
| E_j (A, tight, N2ref) j=0 | 1.2015e-06 | 1.20198e-06 | 1.000394 | [cert] |
| E_j (A, tight, N2ref) j=1 | 0.000104384 | 0.000104388 | 1.000047 | [cert] |
| E_j (A, tight, N2ref) j=2 | 0.0566892 | 0.0566895 | 1.000005 | [cert] |
| Ehat = EA_ref j=0 | 1.2015e-06 | 1.20198e-06 | 1.000394 | [cert] |
| Ehat = EA_ref j=1 | 0.000104384 | 0.000104388 | 1.000047 | [cert] |
| Ehat = EA_ref j=2 | 0.0566892 | 0.0566895 | 1.000005 | [cert] |
| EA_ref4 (4 N2ref) j=0 | 1.21296e-06 | 1.21343e-06 | 1.000392 | [cert] |
| EA_ref4 (4 N2ref) j=1 | 0.000104511 | 0.000104516 | 1.000047 | [cert] |
| EA_ref4 (4 N2ref) j=2 | 0.0582464 | 0.0582466 | 1.000005 | [cert] |
| EA0 j=0 | 1.19768e-06 | 1.19816e-06 | 1.000395 | [cert] |
| EA0 j=1 | 0.000104341 | 0.000104346 | 1.000047 | [cert] |
| EA0 j=2 | 0.0561702 | 0.0561705 | 1.000005 | [cert] |
| C0_ref j=0 | 1.13376e-06 | 1.13415e-06 | 1.000343 | [cert] |
| C0_ref j=1 | 0.000166208 | 0.000166212 | 1.000024 | [cert] |
| C1_0 j=0 | 7.20596e-10 | 7.21146e-10 | 1.000763 | [cert] |
| C1_0 j=1 | 4.31427e-08 | 4.31433e-08 | 1.000014 | [cert] |
| C1_ref j=0 | 1.47907e-09 | 1.47962e-09 | 1.000372 | [cert] |
| C1_ref j=1 | 5.15665e-08 | 5.15671e-08 | 1.000012 | [cert] |
| S0v (core sup|w'|) | 0.00716755 | 0.00716838 | 1.000116 | [cert] |
| S1v | 0.103578 | 0.103579 | 1.000013 | [cert] |
| K0 at E (flux) | 2.43409e-21 | 2.43502e-21 | 1.000383 | [cert] |
| K2 at E (flux) | 6.70519e-21 | 6.70756e-21 | 1.000354 | [cert] |
| margin S0v (tight, N2ref) | 2.77401 | 2.77369 | 0.999885 | [cert] |
| margin S0v_vs_0.0153 (tight, N2ref) | 2.12212 | 2.12187 | 0.999885 | [cert] |
| margin S1v_vs_0.6 (tight, N2ref) | 5.78795 | 5.78788 | 0.999987 | [cert] |
| margin M0core (tight, N2ref) | 5.45551 | 5.45528 | 0.999959 | [cert] |
| margin M1core (tight, N2ref) | 4.43598 | 4.43592 | 0.999988 | [cert] |
| margin ramp M1 (tight, N2ref) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin ramp z1 (tight, N2ref) | 1.20143 | 1.20143 | 1.000000 | [cert] |
| margin ramp M0 (tight, N2ref) | 7.02045 | 7.02045 | 1.000000 | [cert] |
| margin Ghat22 (tight, N2ref) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| margin win M1 (tight, N2ref) | 2.05158 | 2.05158 | 1.000000 | [cert] |
| margin win M0 (tight, N2ref) | 2.54378 | 2.54378 | 0.999999 | [cert] |
| margin win z2 (tight, N2ref) | 4.26836 | 4.26836 | 1.000000 | [cert] |
| margin |f|<=1/2 (tight, N2ref) | 5.38808 | 5.38808 | 0.999999 | [cert] |
| margin K2 (tight, N2ref) | 5.08042 | 5.08042 | 1.000000 | [cert] |
| margin K3 (tight, N2ref) | 4.82451 | 4.82451 | 1.000000 | [cert] |
| margin S0v (tight, 4 N2ref) | 2.76626 | 2.76594 | 0.999885 | [cert] |
| margin S0v_vs_0.0153 (tight, 4 N2ref) | 2.11619 | 2.11595 | 0.999885 | [cert] |
| margin S1v_vs_0.6 (tight, 4 N2ref) | 5.74823 | 5.74815 | 0.999987 | [cert] |
| margin M0core (tight, 4 N2ref) | 5.45005 | 5.44982 | 0.999959 | [cert] |
| margin M1core (tight, 4 N2ref) | 4.40903 | 4.40898 | 0.999988 | [cert] |
| margin ramp M1 (tight, 4 N2ref) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin ramp z1 (tight, 4 N2ref) | 1.20143 | 1.20143 | 1.000000 | [cert] |
| margin ramp M0 (tight, 4 N2ref) | 7.02045 | 7.02045 | 1.000000 | [cert] |
| margin Ghat22 (tight, 4 N2ref) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| margin win M1 (tight, 4 N2ref) | 2.05158 | 2.05158 | 1.000000 | [cert] |
| margin win M0 (tight, 4 N2ref) | 2.54378 | 2.54378 | 0.999999 | [cert] |
| margin win z2 (tight, 4 N2ref) | 4.26836 | 4.26836 | 1.000000 | [cert] |
| margin |f|<=1/2 (tight, 4 N2ref) | 5.38808 | 5.38808 | 0.999999 | [cert] |
| margin K2 (tight, 4 N2ref) | 5.08042 | 5.08042 | 1.000000 | [cert] |
| margin K3 (tight, 4 N2ref) | 4.82451 | 4.82451 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=5/4 (tight, N2ref) | 5.4061 | 5.40403 | 0.999616 | [cert] |
| flux margin_FLUX0 beta=5/4 (tight, 4 N2ref) | 5.35602 | 5.35398 | 0.999618 | [cert] |
| flux margin_FLUX2 beta=5/4 (tight, N2ref) | 9.57203 | 9.56863 | 0.999645 | [cert] |
| flux margin_FLUX2 beta=5/4 (tight, 4 N2ref) | 9.48113 | 9.47778 | 0.999647 | [cert] |
| flux margin_FLUX0 beta=2 (tight, N2ref) | 3.39666 | 3.39535 | 0.999615 | [cert] |
| flux margin_FLUX0 beta=2 (tight, 4 N2ref) | 3.36511 | 3.36382 | 0.999617 | [cert] |
| flux margin_FLUX2 beta=2 (tight, N2ref) | 6.01283 | 6.0107 | 0.999644 | [cert] |
| flux margin_FLUX2 beta=2 (tight, 4 N2ref) | 5.95559 | 5.95348 | 0.999646 | [cert] |
| margin S0v (v31, N2ref) | 2.28265 | 2.28242 | 0.999896 | [cert] |
| margin S0v_vs_0.0153 (v31, N2ref) | 1.74623 | 1.74605 | 0.999896 | [cert] |
| margin S1v_vs_0.6 (v31, N2ref) | 3.37942 | 3.37939 | 0.999991 | [cert] |
| margin M0core (v31, N2ref) | 5.06561 | 5.06539 | 0.999958 | [cert] |
| margin M1core (v31, N2ref) | 2.72101 | 2.72099 | 0.999992 | [cert] |
| margin ramp M1 (v31, N2ref) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin ramp z1 (v31, N2ref) | 1.20143 | 1.20143 | 1.000000 | [cert] |
| margin ramp M0 (v31, N2ref) | 7.02045 | 7.02045 | 1.000000 | [cert] |
| margin Ghat22 (v31, N2ref) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| margin win M1 (v31, N2ref) | 2.05158 | 2.05158 | 1.000000 | [cert] |
| margin win M0 (v31, N2ref) | 2.54378 | 2.54378 | 0.999999 | [cert] |
| margin win z2 (v31, N2ref) | 4.26836 | 4.26836 | 1.000000 | [cert] |
| margin |f|<=1/2 (v31, N2ref) | 5.38808 | 5.38808 | 0.999999 | [cert] |
| margin K2 (v31, N2ref) | 5.08042 | 5.08042 | 1.000000 | [cert] |
| margin K3 (v31, N2ref) | 4.82451 | 4.82451 | 1.000000 | [cert] |
| margin S0v (v31, 4 N2ref) | 2.27743 | 2.27719 | 0.999897 | [cert] |
| margin S0v_vs_0.0153 (v31, 4 N2ref) | 1.74223 | 1.74205 | 0.999897 | [cert] |
| margin S1v_vs_0.6 (v31, 4 N2ref) | 3.37385 | 3.37382 | 0.999991 | [cert] |
| margin M0core (v31, 4 N2ref) | 5.06092 | 5.06071 | 0.999958 | [cert] |
| margin M1core (v31, 4 N2ref) | 2.71685 | 2.71683 | 0.999992 | [cert] |
| margin ramp M1 (v31, 4 N2ref) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin ramp z1 (v31, 4 N2ref) | 1.20143 | 1.20143 | 1.000000 | [cert] |
| margin ramp M0 (v31, 4 N2ref) | 7.02045 | 7.02045 | 1.000000 | [cert] |
| margin Ghat22 (v31, 4 N2ref) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| margin win M1 (v31, 4 N2ref) | 2.05158 | 2.05158 | 1.000000 | [cert] |
| margin win M0 (v31, 4 N2ref) | 2.54378 | 2.54378 | 0.999999 | [cert] |
| margin win z2 (v31, 4 N2ref) | 4.26836 | 4.26836 | 1.000000 | [cert] |
| margin |f|<=1/2 (v31, 4 N2ref) | 5.38808 | 5.38808 | 0.999999 | [cert] |
| margin K2 (v31, 4 N2ref) | 5.08042 | 5.08042 | 1.000000 | [cert] |
| margin K3 (v31, 4 N2ref) | 4.82451 | 4.82451 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=5/4 (v31, N2ref) | 4.60062 | 4.59901 | 0.999649 | [cert] |
| flux margin_FLUX0 beta=5/4 (v31, 4 N2ref) | 4.56435 | 4.56275 | 0.999650 | [cert] |
| flux margin_FLUX2 beta=5/4 (v31, N2ref) | 7.22166 | 7.21958 | 0.999712 | [cert] |
| flux margin_FLUX2 beta=5/4 (v31, 4 N2ref) | 7.16989 | 7.16783 | 0.999713 | [cert] |
| flux margin_FLUX0 beta=2 (v31, N2ref) | 2.88939 | 2.88838 | 0.999648 | [cert] |
| flux margin_FLUX0 beta=2 (v31, 4 N2ref) | 2.86655 | 2.86555 | 0.999649 | [cert] |
| flux margin_FLUX2 beta=2 (v31, N2ref) | 4.5334 | 4.53209 | 0.999711 | [cert] |
| flux margin_FLUX2 beta=2 (v31, 4 N2ref) | 4.50083 | 4.49953 | 0.999712 | [cert] |
| margin S0v_vs_0.0153 (tight, lo) | 2.22703 | 2.22675 | 0.999875 | [cert] |
| margin M1core (tight, lo) | 4.43598 | 4.43592 | 0.999988 | [cert] |
| margin ramp M1 (tight, lo) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin Ghat22 (tight, lo) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=2 (tight, lo) | 4.05707 | 4.05537 | 0.999580 | [cert] |
| flux margin_FLUX2 beta=2 (tight, lo) | 7.06956 | 7.06686 | 0.999618 | [cert] |
| margin S0v_vs_0.0153 (tight, hi) | 2.0369 | 2.03668 | 0.999892 | [cert] |
| margin M1core (tight, hi) | 4.43598 | 4.43592 | 0.999988 | [cert] |
| margin ramp M1 (tight, hi) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin Ghat22 (tight, hi) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=2 (tight, hi) | 2.92158 | 2.92054 | 0.999643 | [cert] |
| flux margin_FLUX2 beta=2 (tight, hi) | 5.23159 | 5.22984 | 0.999666 | [cert] |
| margin S0v_vs_0.0153 (tight, probe) | 2.62879 | 2.62836 | 0.999837 | [cert] |
| margin M1core (tight, probe) | 4.47465 | 4.47459 | 0.999988 | [cert] |
| margin ramp M1 (tight, probe) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin Ghat22 (tight, probe) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=2 (tight, probe) | 7.23471 | 7.23068 | 0.999442 | [cert] |
| flux margin_FLUX2 beta=2 (tight, probe) | 11.7621 | 11.7565 | 0.999525 | [cert] |
| margin S0v_vs_0.0153 (tight, N2=0) | 2.12411 | 2.12386 | 0.999885 | [cert] |
| margin M1core (tight, N2=0) | 4.44515 | 4.44509 | 0.999988 | [cert] |
| margin ramp M1 (tight, N2=0) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin Ghat22 (tight, N2=0) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=2 (tight, N2=0) | 3.4073 | 3.40599 | 0.999615 | [cert] |
| flux margin_FLUX2 beta=2 (tight, N2=0) | 6.03216 | 6.03001 | 0.999644 | [cert] |
| margin S0v_vs_0.0153 (tight, pm3ibp) | 2.13972 | 2.13947 | 0.999884 | [cert] |
| margin M1core (tight, pm3ibp) | 4.46561 | 4.46555 | 0.999988 | [cert] |
| margin ramp M1 (tight, pm3ibp) | 1.18473 | 1.18473 | 1.000000 | [cert] |
| margin Ghat22 (tight, pm3ibp) | 1.91741 | 1.91741 | 1.000000 | [cert] |
| flux margin_FLUX0 beta=2 (tight, pm3ibp) | 3.41622 | 3.4149 | 0.999614 | [cert] |
| flux margin_FLUX2 beta=2 (tight, pm3ibp) | 6.05286 | 6.0507 | 0.999643 | [cert] |
| c_dU | 1.7186e-08 | 1.96498e-08 | 1.143357 | [cert] |
| Delta_in j=0 | 3.13092e-12 | 4.09294e-12 | 1.307265 | [cert] |
| Delta_in j=1 | 3.28686e-12 | 4.29679e-12 | 1.307265 | [cert] |
| Delta_in j=2 | 3.45056e-12 | 4.5108e-12 | 1.307265 | [cert] |
| Nbar^2 (Lemma N, beta=2) k=0 | 0.0714713 | 0.0714723 | 1.000013 | [cert] |
| Nbar^2 (Lemma N, beta=2) k=1 | 0.330111 | 0.330114 | 1.000007 | [cert] |
| Nbar^2 (Lemma N, beta=2) k=2 | 5.04054 | 5.04058 | 1.000008 | [cert] |
| Nbar^2 / (4 N2ref) k=0 | 0.178678 | 0.178681 | 1.000013 | [cert] |
| Nbar^2 / (4 N2ref) k=1 | 0.412639 | 0.412642 | 1.000007 | [cert] |
| Nbar^2 / (4 N2ref) k=2 | 0.252027 | 0.252029 | 1.000008 | [cert] |
| Bhat (supersolution) k=0 | 4.87347e-06 | 4.87349e-06 | 1.000004 | [cert] |
| Bhat (supersolution) k=1 | 0.000559323 | 0.000559356 | 1.000060 | [cert] |
| Bhat (supersolution) k=2 | 0.366742 | 0.366747 | 1.000013 | [cert] |
| U-fine Bhat k=0 | 4.87221e-06 | 4.87224e-06 | 1.000006 | [cert] |
| U-fine Bhat k=1 | 0.00026266 | 0.00026266 | 1.000002 | [cert] |
| U-fine Bhat k=2 | 0.0258413 | 0.0258416 | 1.000011 | [cert] |
| EA (Lemma U fixed point) k=0 | 1.20173e-06 | 1.2022e-06 | 1.000394 | [cert] |
| EA (Lemma U fixed point) k=1 | 0.000104402 | 0.000104407 | 1.000047 | [cert] |
| EA (Lemma U fixed point) k=2 | 0.0566935 | 0.0566938 | 1.000005 | [cert] |
| C0 (U) k=0 | 1.13395e-06 | 1.13434e-06 | 1.000343 | [cert] |
| C1 (U) k=0 | 1.52408e-09 | 1.52463e-09 | 1.000366 | [cert] |
| C0 (U) k=1 | 0.000166223 | 0.000166227 | 1.000024 | [cert] |
| C1 (U) k=1 | 5.52548e-08 | 5.52555e-08 | 1.000013 | [cert] |
| margin late_M1 (Lemma U) | 2.35381 | 2.35377 | 0.999983 | [cert] |
| margin late_M0 (Lemma U) | 4.89246 | 4.89209 | 0.999924 | [cert] |
| margin late_M0z (Lemma U) | 12.9929 | 12.9904 | 0.999812 | [cert] |
| margin present_M1_min (Lemma U) | 2.7708 | 2.77077 | 0.999989 | [cert] |
| margin present_M0_min (Lemma U) | 2.22514 | 2.22514 | 1.000000 | [cert] |
| margin present_M0z_min (Lemma U) | 14.3774 | 14.3763 | 0.999925 | [cert] |
| R<2 margin j=0 (Lemma U) | 4.43231 | 4.43231 | 1.000000 | [cert] |
| R<2 margin j=1 (Lemma U) | 6.69783 | 6.69783 | 1.000000 | [cert] |
| rate_margin (Lemma U) | 28.2313 | 28.2313 | 1.000000 | [cert] |
| loop gain total j=0 | 1.22518e-07 | 1.22571e-07 | 1.000428 | [cert] |
| loop gain total j=1 | 0.0371686 | 0.0371686 | 1.000000 | [cert] |
| loop gain total j=2 | 0.0390342 | 0.0390342 | 1.000000 | [cert] |
| M1z check: zone III b1 + pU1 (vs M1z 0.534) | 0.50786 | 0.50786 | 1.000000 | [cert] |
| M1z check: S1(2 Bhat) zones I/II | 0.285844 | 0.285849 | 1.000019 | [cert] |
| zone III n_III k=0 | 0.0714447 | 0.0714456 | 1.000013 | [cert] |
| zone III n_III k=1 | 0.327066 | 0.327068 | 1.000006 | [cert] |
| zone III n_III k=2 | 3.04441 | 3.04442 | 1.000002 | [cert] |
| zone III DeltaR k=0 | 4.49852e-08 | 4.49852e-08 | 1.000000 | [cert] |
| zone III DeltaR k=1 | 6.79715e-06 | 6.79715e-06 | 1.000000 | [cert] |
| zone III DeltaR k=2 | 0.01148 | 0.01148 | 1.000000 | [cert] |
| zone III DeltaR_b k=0 | 2.65313e-12 | 2.65313e-12 | 1.000000 | [cert] |
| zone III DeltaR_b k=1 | 3.25182e-09 | 3.25182e-09 | 1.000000 | [cert] |
| zone III DeltaR_b k=2 | 4.34718e-08 | 4.34718e-08 | 1.000000 | [cert] |
| zone III box_380 b0 | 0.100314 | 0.100314 | 1.000000 | [cert] |
| zone III box_380 margin b0 | 1.4953 | 1.4953 | 1.000000 | [float] |
| zone III box_380 b1 | 0.419731 | 0.419731 | 1.000000 | [cert] |
| zone III box_380 margin b1 | 1.23889 | 1.23889 | 1.000000 | [float] |
| zone III box_380 bz | 0.0142902 | 0.0142902 | 1.000000 | [cert] |
| zone III box_380 margin bz | 11.1965 | 11.1965 | 1.000000 | [float] |
| zone III box_700 b0 | 0.227569 | 0.227569 | 1.000000 | [cert] |
| zone III box_700 b1 | 0.542338 | 0.542338 | 1.000000 | [cert] |
| H6 cone K2 margin (at Bhat_ref) | 2.54818 | 2.54818 | 1.000000 | [cert] |
| H6 cone K2 margin (at fine est) | 2.54818 | 2.54818 | 1.000000 | [est] |
| H6 cone K3 margin (at Bhat_ref) | 2.51066 | 2.51065 | 0.999997 | [cert] |
| H6 cone K3 margin (at fine est) | 2.27459 | 2.27459 | 1.000000 | [est] |
| H6 cone slope margin (at Bhat_ref) | 2.46299 | 2.46299 | 1.000000 | [cert] |
| H6 cone slope margin (at fine est) | 2.21043 | 2.21043 | 1.000000 | [est] |
| H6 far K2 margin (at Bhat_ref) | 2.28566 | 2.28566 | 1.000000 | [cert] |
| H6 far K2 margin (at fine est) | 2.28566 | 2.28566 | 1.000000 | [est] |
| H6 far K3 margin (at Bhat_ref) | 1.86087 | 1.86087 | 1.000000 | [cert] |
| H6 far K3 margin (at fine est) | 1.86087 | 1.86087 | 1.000000 | [est] |
| H6 far slope margin (at Bhat_ref) | 2.85617 | 2.85617 | 1.000000 | [cert] |
| H6 far slope margin (at fine est) | 2.85617 | 2.85617 | 1.000000 | [est] |
| zone III tstar margin | 2.30194 | 2.30194 | 1.000000 | [cert] |
| zone III box margin min | 1.23889 | 1.23889 | 1.000000 | [float] |
| zone III DeltaR2 (acceptance) | 0.01148 | 0.01148 | 1.000000 | [cert] |
| c_dU vs typed bound 2e-8 (u1_ghat M1z, zone3 assert) | 1.16373 | 1.01782 | 0.874624 | [float of cert] |

Summary (key_table.py): margins listed: 132; largest relative decrease: flux margin_FLUX0 beta=2 (tight, probe) x0.999442; smallest new margin: margin ramp M1 (tight, N2ref) = 1.18473

## 6. Proof tree untouched
- **The full tree is unchanged.** I took a sha256 snapshot of 4292 files before the work and again after it. The files cover `writeup/lemmaU` (excluding `fix4_r3`), `step4/p2/out_v33`, `step4/p2/v33`, `step3nl_v33/verifier`, `writeup/lemmaO_v32` and `writeup/lemmaS`. The two snapshots are identical.
- **Nothing was written outside the new directory.** `find /home/claude -newer <start marker>` outside `N` lists no file.
- **Existing SHA256SUMS** (`sha256sum -c`):

| file | OK | FAILED |
|---|---|---|
| `B/fix2/n12/t2_A_v35/out/SHA256SUMS` | 20 | 0 |
| `B/t4_cert/out/SHA256SUMS` | 34 | 0 |
| `B/t5_zone3/out/SHA256SUMS` | 10 | 0 |
| `B/t5_zone3/SHA256SUMS_scripts` | 4 | 0 |
| `B/t4_cert/out/SHA256SUMS_scripts_m14` (checked from `t4_cert/`) | 6 | 0 |

- **`B/fix2/n12/t2_A_v35/out/SHA256SUMS_scripts_m14` fails 2 of 2 (rigor/assembly_rep.py, rigor/export_t2.py).** This is pre-existing. The scripts were edited after M14, in N1, and the snapshot above shows I changed nothing. The file is stale, not damaged.
- **Pinned hashes are unchanged:**
  - A_interface 78fe2f8e...;
  - U_cert_v35 960809fd...;
  - zone3 45e8d4dc...;
  - t2_A_v35 SHA256SUMS 7791c6af...;
  - closed `FINAL_manifest.json` e01cb4da...;
  - `FINAL_manifest_v33.json` 8bc47e83....

## 7. What I did not do, and caveats
1. **Inputs from the failed verifier.** The chain consumes the fin12f4 verifier:
   - its slopes, J's and sqrt(Lambda), through `base_ibp/lemmaS_F10p_M0.52.json` args;
   - its sups, D1 (betabar, abar);
   - its rate table and ramp sups.

   All of these change once the verifier is re-derived. R3 is kept separable:
   - the `r3_forcing.py` module, with inputs passed in;
   - three small patches.

   Reapply the patches on the new chain and rerun.
2. **Texts not edited.** These are documents (new versions to be made later):
   - `lemmaS_A_v33.tex` section "Step 6: forcing" (r0^F, A:266);
   - `numbers_A_v33.tex`;
   - `lemmaU_v34.tex` Lemma dU (hypothesis alpha_tau = -(p/4) alpha + beta, and the value 1.7187e-8);
   - `step3_proof_v32.tex:235` (closed);
   - `main_v3.tex`, the erratum, the V33_INDEX.
3. **Not rerun:** step 4, the kits and the manifests (as instructed). The step-4 stages that pin `zone3_interface.json` (stage 9 and later) will need re-pinning when the final chain is made.
4. **`ramp_U.py` was not rerun.** It is unaffected: it reads u1_ghat, which is byte-identical with and without R3, and it uses no c_dU. Its v35 output was copied as an input (byte-identical).
5. **Unused stale constant.** The unused function `assembly_arb.assemble` (line 106, run only by that module's `__main__`) still contains the old 1.4 constant. It is not called in the chain. I left it unchanged, to keep the patch minimal.
6. **The N35 assertions are not inside the scripts.** These are the ramp-box check and N2 <= 4 N2ref. I used the pinned scripts, as v35 does. The asserted quantities are reported in section 4, and both hold.
7. **Remaining dependence on the typed bound 2e-8.** u1_ghat still uses the typed c_dU <= 2e-8 in its M1z check, and its room is now x1.018 instead of x1.164. This holds [cert], but any further growth of betabar or abar from the re-derived verifier could break it. The bound would then have to be raised; the check itself has room: M1z check 0.50786 against 0.534.
8. **Not in scope and unchanged:**
   - the flux-lemma gap (U2 R3.3, fix3/F1);
   - G0's verifier findings;
   - Lemma U's c'_j of Prop. B' ([paper input]);
   - the eps-part of Lemma U's own forcing R^.

## 8. SHA256SUMS of the new outputs
The outputs in `N/{t2_A,t4_cert,t5_zone3}/out` are the r3run2 outputs. They are byte-identical to `N/runs/r3run1`; for each JSON this was checked by `cmp`. The list is `N/SHA256SUMS` (sha256 of the list: 288043a8d43a765dfef0b31459f6a60b3a98f3e037114b48cd8ce50fb95ec66e). The scripts are in `N/SHA256SUMS_scripts` (eb701df700a69220b8d742df562d0a812708331cc3e16a3b5213ca8dd25ae28b).
```
1c463361805799a73f22330428e8ef74619dbbe0e4c8118821df7feb141832c6  r3/r3_forcing.json
4820937c2ffcae116931725367812f0aba6c5e32b0d83411e291c6cddbb9e1c5  r3/r3_symcheck.json
d8d41fbebbadea175a0b375cce978301f8209ea65311192c798c8d47e1503f94  t2_A/out/A_interface.json
1e9debc352ff9e0489e4d5fe4fdcd63cf0e9cd38def06601685e13d6c3d3b9d9  t2_A/out/a3_signs.json
0c93031c13d78c1b898fb067ab9a92be54737b9024a45786c49b6b90a8cc2f19  t2_A/out/asm_design_N20.json
37f03127a7ccd21b2f886ca3dbb4f82d3fadf9b026b83c3f8aeba35138a42b0f  t2_A/out/asm_design_N2ref.json
14242e6dfc00e95e0661ced5beac4e558fdd78a70dd46e0861a4154f686743ea  t2_A/out/asm_design_N2ref4.json
a3d4aa0a6d8814b1c055b92eb0b5d1c105abbfc70b882650bd9e301760258884  t2_A/out/asm_design_pm3ibp.json
385bce154bafe1421888027e63ab3a250078879280210338abaf87648f8803d8  t2_A/out/asm_hi_N2ref.json
ee29344b38d432f52575e7d1b7c5fad99cd493818bf275d526c507580ef9899a  t2_A/out/asm_lo_N2ref.json
936755df2fe50cadc322b696459ae110a24963da78371f62d621f138a16eb54b  t2_A/out/asm_probe_N2ref.json
3d37330abe6fbf9700f9451564d8530bed559ff7c88b2e0f1cfbefaf0fd4e345  t2_A/out/hjw_rep.json
b89253a8ad439e7becbae992f9c927fe8c8858f4e727e1061324046751d2d140  t2_A/out/ic_design.json
cecb268c6e4bb29af5af9eff929c0bc4ae372252913f0d04c75662a0d88d694e  t2_A/out/ic_hi.json
bb62d7472d8a96419f0a1c171796f57564b59ed9d31b5bbaa80a93b0f4d73f00  t2_A/out/ic_lo.json
eb74f968ca62a6e837f26af7875417c51ee9b6d924c71c5b9ed26968997cfddc  t2_A/out/ic_probe.json
6d846bf0636b3d0af410c8a4917b9a134f58e26970144cc25c6ce2692cafe21d  t2_A/out/outer_L2.json
214fcfaa25111a3755fc4c62503722e8dd9147e9f09fd133bcd0c94e339d0333  t2_A/out/ramp_rep.json
619e2764150f076af2cffd6129a1196a2c149b551c1bb982d6ad13dc4bd5afcf  t2_A/out/ramp_rep_pm3ibp.json
68cac85443ef6c49190d16f02e0a87b9b1cebb395678c08f121fb76a3ec8de61  t2_A/out/weights_t2.json
c6998893f992ec83c7683a4bf076bf2c6fc20ce1563b3f6289291125e74826f9  t2_A/out/weights_t2_fb.json
6fced8aca9cc07f3405cb96e74aabbe25490b6db281f3313870b71070b8135f1  t4_cert/out/U_cert_for_zone3_r3.json
71362358e1c8541c522d15dfdf789d729fbaa30b54e4ff12073597b54836843e  t4_cert/out/U_cert_for_zone3_r3_pass1.json
fb877fdae4fca9bc6fd12700c54f96811a007df499a3dbc4ec39366522a9a117  t4_cert/out/U_cert_r3_pass1.json
466ac45c7ff4c22afacb95ac9dd0e5704a32aae4ae977bbe2bab148b2861f289  t4_cert/out/U_cert_r3_pass2.json
df4499068c8fff6741103a9753a590d9ff2a8dc1e0bf82d93d15f2bb74e52935  t4_cert/out/U_cert_r3_pass3.json
e66446a17629c2f3c3072c45b270ecfb3bd7dafc004992cfcafc08f5410fd0bc  t4_cert/out/u1_ghat.json
8aa0a99178e71dd9e81d30e1353a2b6fd93de7783e4b6007d240cc64acb714c5  t4_cert/out/u2_K.json
0449711acda75fdae3cac593650a5257b587a01c93830408019a354f633c42b9  t4_cert/out/u34_primary.json
55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de  t5_zone3/out/ramp_U.json
58bc892540ce59d051c270e97a941844f21880759b408be855d55295bfaca1f7  t5_zone3/out/zone3_interface.json
82713f5c10b027db3ed098fa9499313735c4eaec6ca994537c69b92f4b6ea742  t5_zone3/out/zone3_interface_pass1.json
```
