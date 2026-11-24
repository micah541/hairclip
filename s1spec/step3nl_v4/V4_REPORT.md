# V4_REPORT: verifier v4 (tag fin13v4), the corrected E_I inequality of EI_SPEC.md

Nov 19, 2026. Work order: `writeup/lemmaU/bulk/fix5/ei_spec/EI_SPEC.md` §10. Base: `step3nl_v33` (fin12f4), unchanged. Everything new is in `/home/claude/s1spec/step3nl_v4/`.
Running log: `V4_LOG.md`. Hashes: `SHA256SUMS` (130 files).

## Verdict
**fin13v4 closes at both points, alpha1 = 1.5e-4 and alpha1 = 1.25e-4, with every term of EI_SPEC §5 in it and every constant read from raw certified JSON and rounded up.** No margin fails. The only lever used is L1: the far coefficients are recertified at the alpha-level 2.5e-4, by zone. No stop condition was hit:
- the U-3 restatement goes through;
- J moves by at most 0.14%;
- sup|beta| goes down by 29%, so c_dU is not threatened.

These results are certified for the E_I block. Step 3 is **not** closed by this. See "What remains unjustified" below, especially item 1, a side finding in the E0 block.

Checks:
- **Regression T(i): PASS.** With NEW=0 the v4 code reproduces fin12f4 bit for bit, apart from timing and the new keys. This covers both rate tables, final_alpha_points and the tail_cert log. It also reproduces the two cells quoted in the spec: first cell U = 0.03034405538749633 and tail anchor U = -0.0026930418625762047 (`runs/reg_fin12f4/t_regression.json`).
- **Double run: PASS.** The full chain (certification stages, then production) was run twice: `v4cert/out_A`+`runs/P_A` and `v4cert/out`+`runs/P`. The outputs are byte-identical apart from timing fields (`runtime_s`, `secs`, `cell ...s`) and the sha256 pins of upstream outputs that carry those fields (`runs/t_double_run.json`: 8 identical, 17 equal apart from timing).
- **tail_cert: PASS** on the production code path (`runs/P/tail_cert_log.txt`): radial H (0,2), L (0,1); s-mode H (2,5); joint with floors off H (2,5), L (2,5); f1 = 0.
- **T(iii) float point test: PASS.** An independent plain-float evaluation of §5.1, §5.2, the E0 strip term and L, written from the spec text, lies inside the interval rhs at 400 random points (`verifier/t_float_check.json`).
- **Unchanged files (sha256): step3nl_v33** (555 files), A_interface v35 78fe2f8e…, U_cert_v35 960809fd…, zone3 45e8d4dc…, manifest_v33 8bc47e83…, and the closed FINAL_manifest e01cb4da… (`audit/*_before.sha` = `audit/*_after.sha`).
- **Resources:** one process at a time, nice 19, PYTHONHASHSEED=0. Peak RSS: 571 MB (S1), 214 MB (C2), 79 MB (production). Production took 11 min per pass.
- **One RSS overrun:** an early version of C2 peaked at 1.1 GB for about 10 min. It has been fixed (see the log).

## Old vs new [cert; doubles as stored; ratio = fin13v4/fin12f4]
alpha1 = 1.5e-4 (sqrt(Lambda_1) = 5e-8):

| quantity | fin12f4 | fin13v4 | ratio |
|---|---|---|---|
| steps to the tail | 10 | 10 | |
| peak sqrt(Lambda) | 5.307458e-08 | 5.307698e-08 | 1.000045 |
| tail rho | 2.693042e-03 | 2.684366e-03 | 0.99678 |
| rho / rho_min (0.002) | 1.3465 | 1.3422 | |
| sqrt(L_F) | 9.871704e-10 | 9.887644e-10 | 1.0016 |
| alpha_F | 6.866519e-07 | 6.876746e-07 | 1.0015 |
| alpha-condition a*/alpha_F | 53.66 | 53.58 | |
| V1 delta_t | 0.2 | 0.2 | |
| J0 (V0^2) const / slope | 6.027237e-12 / 3.454743e-18 | 6.020438e-12 / 3.452361e-18 | 0.99887 / 0.99931 |
| J1 (V1^2) | 8.467211e-12 / 4.853308e-18 | 8.461143e-12 / 4.85196e-18 | 0.99928 / 0.99972 |
| J2 (V2^2) | 1.870781e-12 / 1.07231e-18 | 1.872555e-12 / 1.073798e-18 | 1.00095 / 1.00139 |
| J4 (V4^2) | 6.236181e-08 / 3.574508e-14 | 6.242096e-08 / 3.579469e-14 | 1.00095 / 1.00139 |
| sup\|(a,b)\| (6 digits up) | 2.19103e-07 | 2.14024e-07 | 0.9768 |
| sup\|beta\| | 3.33811e-09 | 2.35741e-09 | 0.7062 |
| sI1, sI2 | 1.160338e-3, 1.545618e-2 | 1.160390e-3, 1.545688e-2 | 1.000045 |
| M1_ramp, S_ramp | 0.432773, 7.473592e-4 | 0.4327926, 7.473931e-4 | 1.000045 |
| first cell rate U | +3.034406e-2 | +3.035649e-2 | 1.00041 |
| floor F0 | 1.312192e-21 | 1.312192e-21 | 1 |

alpha1 = 1.25e-4 (remark point):
- steps 9 / 9; peak 5.213518e-8 → 5.213696e-8 (×1.000034);
- rho 2.860963e-3 → 2.855715e-3 (×0.99817; rho/rho_min 1.428);
- sqrt(L_F) ×1.0009; a*/alpha_F 54.34 → 54.29; delta_t 0.2;
- J0..J4 const ×0.99849 / 0.99891 / 1.00057 / 1.00057;
- sup|(a,b)| 2.1537e-7 → 2.10289e-7; sup|beta| 3.20505e-9 → 2.22438e-9.

The full tables are in `runs/v4_compare.md` and `runs/v4_compare.json`.

Box hypotheses are asserted at every step and at the tail anchor. The worst value/bound ratios along the run (1.5e-4) are:
- |w| on [34.25,60] (≤ 0.5): 0.789;
- S2 on the strip (≤ 0.02): 0.773;
- S0V (≤ 0.02): 0.765;
- IR jet (≤ 1e-3): 0.747;
- modulation-lemma alpha (≤ 5e-4): 0.32;
- everything else ≤ 0.06.

Reading the numbers:
- The new terms cost 0.3% of rho, as the probe predicted (spec 0.2%).
- J0 and J1 fall slightly, because MU1 now comes from the certified 0.1555226 instead of the typed 0.1552.
- sup|beta| falls because of the U-3 strip constant (below).
- Downstream, (A) uses kappa_J = 0.04235 against a limit of 0.13, so J changes of 0.1% are harmless there. They do exceed the 6e-7 rounding reserve of ROUTE C6, so the (A) rerun with the fin13v4 J's (already planned, §10.4) is required.
- c_dU scales with betabar, which fell by 29%. It is not threatened, and it gains room.

## Spec items: status
| item | status |
|---|---|
| S0 setup | done. `cp -a` of v33; old outputs removed; external inputs pinned (`v4cert/v4common.py: EXT`) |
| S1 identities (U-2) | done [sym]: chain rule D X = X^a + X^b g'; c = 1 - J^-4; D(d000) = E1c g' + E0c with E0c(0,0) = 0 (generic G); N_low' table; full v-equation T1..T12 = zeta g_tau. All residuals are exactly 0. The FD cross-check of the table against N_low built from F itself agrees to 3e-10 relative, limited by the mpmath ODE solver (`v4cert/out/s1_identities.json`) |
| C1 cut_sups | done [cert]: PSI1 1.304871, PSI2 1.695438, RWR 0.655792, KWIR 70258.13, KAP IR/Z1/OR 0.53397/1.28274/1.58832, k4, k3, Z2i 0.21132, T0 37.18449, T2 39.18420 (now typed honestly; c62 deleted), TC0 0.35992, TC2 1.27089, Q2B 2.10e-4. sup\|zeta'\| = 35/96 and 35/48 exactly (sympy) |
| C2 t3d (incl. lever L1) | done [cert]: zones IR, [16.25,34.25], [34.25,60] at levels 2.5e-4 and 1e-3, for c, CA, CB, d, d^a, d^b, A0, A1, B0, B1, kU, DQ0, DQ1, EA0, EA1, using C^4 jets and a5. At 2.5e-4: CBAR 3.365e-3 (typed 2.27e-2), CB 0.165, d110 0.411, d300 3.05, CBAR_IR 6.92e-5. q-dependence is split off by an exact q-Taylor split. Same code path as r60_t3/pb_ramp: 30 of 37 table entries are bitwise equal; the other 7 are lower because v4 evaluates at 200 bits instead of 53 |
| C3 mod_x (U-6) | done [cert], C^4: NP3 0.2595063, NH2 1.795356, NAB 0.273467, NP12 1.7168e-7, C0P 9.875e-12 |
| C4 hardy | done [cert]: hA0/HA0, hB0/HB0, hEA1/HEA1 (Phi upper sums), C0Hb 1.9e-42 |
| C5 constants | done: `v4cert/out/constants_v4.json`, 269 entries with (raw, value, dir, source, key) and an assert that every value is outward. MU1 and MU2 are rounded down, from v2_T4_P2_0.25.json. Corrected down-roundings: AF, CF0/2/3, DS1/2, SB, SA, IA, FN, C3, the far TQ, and (new) KREG[12,16.25] = 3.0208e4 < raw 30208.006, now the arb 1/rho(16.25). BETA0 and CAB0 now include the lemma's dlead |
| R rhs | done: `verifier/rhs_iv.py` (`L_list`, `_ei_v4`, NEW branch), every term of §5.1/§5.2 with zone-wise sups; c62 deleted; FnW→Fn in the (a,b) bound; FNW2 for level 2 |
| unc/rhs share one L list (SF-3) | done: `unc.lk_v4` linearises `rhs_iv.L_list`, the same list rhs_pn sums |
| global sup\|v'\| includes OR | done: sG = max(sI2r, S1V/32, S1V/54.75 + 35/48·S0V, S2r + 35/96·S1r) |
| box hypotheses asserted (SF-4, U-4) | done: B-b, B-c, B-d (incl. \|w(34.25)\| + 25.75 S0V ≤ 0.5), strip S0..S2 ≤ 0.02, and the modulation box are asserted in `closure.box_assert` at every step and at the tail anchor. B-e (S0V, S1V from (A)) remains an input |
| E0 beta-strip term (SF-1) | done: + 7.52·BTs·A in E0, with BTs in state form |
| U-3 strip as a linear form | done [cert/read]: `v4cert/u3_strip_form.py`; it goes through, so the fallback (F0 ×7.4) is not needed. Two corrections found: (i) `linear_alpha_S0` multiplies \|alpha\|·S0 (ibp_run.py), so the old STRIP_BETA counted it with alpha = 1; (ii) the lemma omitted the Neumann feedback on the strip sources (ibp_final.py:75-81), which v4 now adds. New STRIP_BETA 8.04e-11 (was 1.06e-9) and STRIP_AB 4.37e-10 (was 5.52e-9) are used in the alpha-equation and mod_sups |
| strip hypothesis sup\|w'''\| ≤ 0.02 (SF-2) | **open.** In E_I/E0, S3 is replaced by supv2 (a state sup), so no hypothesis is needed there. The constants STRIP_BETA and STRIP_AB (alpha-equation, mod_sups → abar, betabar) still assume S3 ≤ 0.02, and this is not implied by Lambda. Its weight is small: the S3 coefficient is 1.98e-10 of the 4.0e-9 per-unit-S beta form |
| T tests | (i) PASS; (ii) PASS; (iii) PASS; (iv) sanity.py and sanity_unc.py were **not run** (they assume the old configuration) |
| P production + double run | done, PASS |
| D record | this report, plus SHA256SUMS. `V4_VERIFIER.md` was not written separately; this report covers it |
| W write-up (`ei_lemma_v4.tex`, U-1 regularity lemma) | **not done** |
| V Fable review | not done (not in my scope) |

## What remains unjustified
1. **E0 beta-loop constants (side finding, outside EI_SPEC, not fixed).** E0 uses PI3 = 0.01607 = n3·eps1/G33 (the old lemma, diagonal only) and 5.72 = 1 + n3·0.6261 (the pre-feedback cPhi). The IBP lemma used everywhere else has BETA0 = 0.0024398 and BETAF = 0.6565. If E0's beta coupling is n3·|beta|·A, the right values are 0.01835 and 5.937.
   - With those values the tail anchor rate goes from -2.684e-3 to **-2.85e-4**, above the -0.002 the tail lemma needs at this anchor. The first cell goes to +3.340e-2. This is [probe] (`probe_E0/e0probe.json`).
   - So either E0 has a structure that justifies the old constants (the "exact modulation alpha-form" R1, MODA = 25), or the tail must move to smaller alpha (lever L5).
   - **This needs an E0 audit and a decision by Micah before step 3 can be called closed.** The E and E2 blocks were not audited in general. The E2 modulation term 43.33 A² also has no strip or forcing part.
2. SF-2 (above).
3. U-1: the regularity lemma for the energy identities is not written.
4. The outer ramp: K0^OR and K2^OR (EI_SPEC §4.10) are still the inputs FLUX0 and FLUX2 (hypothesis I3). Flux lemma v2 (G2) and the time-integrated g''' parts (G3) are open.
5. Legacy typed literals that EI_SPEC does not cover: LAMT, PI3, MODA, MOD_E2, ALIN, nd*, nf*, r143/r163/r107, t22, t12, f572, 7.52. They are listed in constants_v4.json under `legacy_typed` and were not re-derived.
6. Partly unchecked inputs:
   - the B0/B1 mean-value bound uses the (U0,U1) box; its derivation is checked symbolically (I14), and the sups are [cert];
   - the T12 terms use |cos|, |sin| ≤ 1;
   - the third-party raw JSONs (coupling, sob, v3, c4switch, mod_ibp) are taken as certified and get one extra ulp.
7. The (A) chain, step 4, kits and manifests are untouched. Re-pinning them with fin13v4 needs Micah's go. That is where J, peak, abar, betabar, rho, L_F and alpha_F change; the R3 patches (fix4_r3) must be reapplied there.
8. Bytecode caches: the first runs wrote `__pycache__` files into `step3nl/`, `step3/` and `s1_tight/`. They were deleted (only the directory mtimes changed). The run scripts now set PYTHONDONTWRITEBYTECODE=1; this was added after the double run, and its only effect is on cache writing.

## Files
- certification stages: `v4cert/{s1_identities,c1_cut_sups,u3_strip_form,c3_modx,c2_t3d,c4_hardy,c5_constants}.py`, driver `v4cert/run_cert_chain.sh`, outputs `v4cert/out/` (pass A: `v4cert/out_A/`);
- verifier: `verifier/{constants,rhs_iv,unc,closure,mod_sups}.py` (backups `*_pre_v4.py`, byte copies of v33), `run_v4.sh`, `t_regression.py`, `t_float_check.py`, `t_double_run.py`, `v4_compare.py`;
- runs: `runs/P` (production, pass B), `runs/P_A`, `runs/reg_fin12f4` (regression), `runs/v4_compare.{md,json}`, `runs/t_double_run.json`;
- probe: `probe_E0/`; audit hashes: `audit/`.
