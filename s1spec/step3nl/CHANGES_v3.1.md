# CHANGES v3.1 (step-3 final configuration: outward-rounding pass)

Cause: the kit reproduction check found decimals typed into the certificate on the wrong side of the certified value. All fixed in the originals; the old outputs are kept (`*_v3.0*`, `fin10` tables, `writeup/lemmaS/v3.0_backup/`). No conclusion changes: every page-1 inequality of lemmaS.tex v3 still holds (59/59 in arb), peak sqrt(Lambda), rho, floor, J unchanged to 6 digits, minimum margin still 2.05.

## 1. Verifier (`step3nl/verifier/`)
| file : line | quantity | old | new |
|---|---|---|---|
| constants.py:72 | BETA1 | 66.5 | 66.53 (certified 66.52383) |
| constants.py:73 | CAB0 | 3.708 | 3.7083 (certified 3.708257) |
| constants.py:73 | QAB | 51.3 | 51.34 (certified 51.33974) |
| constants.py:73 | PHIAB | 3.84 | 3.8404 (certified 3.840306) |
| constants.py:77 | STRIP_AB per unit S | 2.5e-7 (max of a,b) | 2.76e-7 (Euclidean hypot(1.2923e-7, 2.4303e-7) = 2.7525e-7) |
Checked and fine (>= certified): BETA0 0.00244 (0.00243977), BETAF 0.657 (0.656528), QB 5.53 (5.526374), MAB 622.4 (622.3777), STRIP_BETA 5.3e-8 (5.2563e-8), FLUX0/FLUX2 4.1e-21/2.0e-20 (6.83e-22/3.26e-21), S0V 0.0153, S1V 0.6. The remaining constants (FN, TQ, CF, K tables, planB ramp table) come from other certified tables that were not recomputed here.
New files: `mod_sups.py` (certified sup|(a,b)|, sup|beta| over all march steps and the tail, outward-rounded to 6 digits, writes `final_alpha_points_5e-8.json`), `run_final.sh` (tag fin11), tables `rate_table_fin11_a0.00015.json`, `rate_table_fin11_a0.000125.json`, logs `verify_log_fin11_*`. Old summary: `final_alpha_points_5e-8.v3.0.json`. tail_cert: PASS (unchanged log).

Regenerated sups (the old values had no generator script; the old sup_ab did NOT dominate the lemma formula):
| alpha1 | quantity | old (reference) | new (mod_sups.py) |
|---|---|---|---|
| 1.5e-4 | sup\|(a,b)\| | 2.18321e-7 (typed 2.183e-7) | **2.19103e-7** (+0.36 %) |
| 1.5e-4 | sup\|beta\| | 3.33786e-9 (typed 3.34e-9) | 3.33811e-9 |
| 1.25e-4 | sup\|(a,b)\| | 2.14589e-7 (typed 2.146e-7) | **2.15369e-7** |
| 1.25e-4 | sup\|beta\| | 3.20485e-9 (typed 3.20e-9) | **3.20505e-9** (old typed 3.20e-9 was below) |
Table values moved by at most 1e-6 relative (peak 5.307447863e-8 -> 5.307447915e-8, rho 2.69303671e-3 -> 2.69303560e-3, floor 4.9358569e-10 -> 4.9358579e-10, J0 6.027230549e-12 -> 6.027231224e-12, J1 8.467201596e-12 -> 8.467202545e-12, J2 1.870778785e-12 -> 1.870778994e-12, J4 6.23617454e-8 -> 6.23617524e-8, slopes 8.6368780e-19 -> 8.6368794e-19, 1.21332985e-18 -> 1.21333004e-18, 2.6807815e-19 -> 2.6807819e-19, 8.93628979e-15 -> 8.93629121e-15; ramp M1 0.43277215437 -> 0.43277215856; sI1 1.160335217e-3 -> 1.160335229e-3; sI2 1.5456148370e-2 -> 1.5456148520e-2).

## 2. Lemma S (`writeup/lemmaS/`), final chain rerun (`run_final.sh`, `run_final_125.sh`, `rigor/capture_F2.py`, `rigor/assembly_arb.py` now reads fin11)
Typed inputs of lemmaS_F2.py (all rounded up, rho down): J 6.0272e-12 -> 6.0273e-12, 8.4672e-12 -> 8.4673e-12 (1.8708e-12 kept); sqLam 5.3074e-8 -> 5.3075e-8; slopes 8.637e-19 -> 8.6369e-19, 1.213e-18 -> 1.2134e-18, 2.681e-19 -> 2.6808e-19; sJ4 8.936e-15 -> 8.9363e-15; abar (--abar/ABAR) 2.183e-7 -> 2.19103e-7; beta 3.34e-9 -> 3.33811e-9; rho 2.693e-3 kept (certified 2.6930356e-3, used as a lower bound in F4B: intA = 2 sq/rho). Bbox uses sq = 8.19e-8 (default SQ, conservative, unchanged).
1.25e-4 point (run_final_125.sh): J 5.6765e-12, 7.9745e-12, 1.7620e-12 (was 1.7619e-12); sqLam 5.2135e-8 -> 5.2136e-8; slopes 8.365e-19 -> 8.3648e-19, 1.175e-18 -> 1.1752e-18, 2.596e-19 -> 2.5964e-19; sJ4 8.655e-15 -> 8.6548e-15; rho 2.861e-3 -> 2.8609e-3; abar 2.146e-7 -> 2.15369e-7; beta 3.2e-9 -> 3.20505e-9.
Page-1 inequalities: all 59 hold with the new inputs; arb values moved by <= 6e-6 relative except I3 slope (1.38239e-27 -> 1.38170e-27, down 5e-4). Changed outputs (arb upper bounds): E0 3.1109872e-7 -> 3.1110055e-7; E1 1.807408e-4 -> 1.807418e-4; E2 0.24704506 -> 0.24704507; S0v 0.0057455489 -> 0.0057455649; S1v 0.17075255 -> 0.17075278; core M0 0.018741114 -> 0.018741130; core M1 0.18439789 -> 0.18439813 (page 0.1844, slack 1e-5); E_B0 3.1441547e-7 -> 3.1441729e-7; K0 6.8285837e-22 -> 6.8286219e-22; K2 flux 3.2515668e-21 -> 3.2515794e-21; I3 const 1.8399834e-21 -> 1.8399376e-21.

### lemmaS.tex page 1 (lines of `lemmaS.tex`)
| line | old | new |
|---|---|---|
| 15 (H-B1) | modulation sups \|(a,b)\| <= 2.183e-7; table file fin10 | **2.192e-7**; fin11 |
| 16 (H-J) | J_{0,1,2} = (6.03, 8.47, **1.87**)e-12 | (6.03, 8.47, **1.88**)e-12 (certified 1.870779e-12) |
| 16 (H-J) | slopes (8.64, **12.1**, **2.68**)e-19 | (8.64, **12.2**, **2.69**)e-19 (certified 1.21333e-18, 2.68078e-19) |
| 17 (H-R) | cert. 1.16e-3 + U(28) | 1.161e-3 + U(28) (sI1 = 1.1603e-3) |
| 25 | (6) Verifier: fin10 | fin11 |
| 27 (Margins) | M0core vs 0.11: **5.87** | **5.86** (certified lower bound 5.8694) |
| 29 (remark 1.25e-4) | peak 5.21e-8; rho 2.861e-3; J (5.68, 7.97, 1.76)e-12; \|(a,b)\| 2.146e-7; \|beta\| 3.20e-9 | 5.22e-8; 2.860e-3; (5.68, 7.98, 1.77)e-12; 2.154e-7; 3.21e-9 |
| 29 | E = (3.05e-7, 1.76e-4, 0.244); S1v 0.169; M0 0.0165; ramp 0.430; boxes >= 2.21; \|f\| <= 0.078 | (3.05e-7, 1.761e-4, 0.245); 0.170; 0.0166; 0.431; >= 2.2; 0.079 |
| 29 | K2 3.19e-21; I3 <= 1.78e-21 + 1.29e-27 | 3.20e-21; 1.79e-21 + 1.30e-27 |
| 31 (Rigor notes) | - | v3.1 sentence added |
`lemmaS_consts.json` key `v3_nov10`: point.peak 5.3074e-8 -> 5.3075e-8; inputs J/J4/slopes/sJ4/abar/beta as above, sI1 0.00116 -> 0.001161, sI2 0.015456 -> 0.015457, verifier_files fin10 -> fin11; core.S1_v 0.1708 kept; alt_alpha1_1.25e-4 block as the remark; `rigor_pass.v3_1` note. (The values E, M0, M1, F4 etc. of page 1 did not change at the printed digits.) lemmaS.pdf rebuilt, 0 LaTeX errors.

## 3. Lemma O (`writeup/lemmaO/lemmaO_v3.tex`; computed numbers unchanged, 9 more quotes than the 4 reported were not outward)
| line | quantity | old | new (certified) |
|---|---|---|---|
| 22 | C2, C3 (X=3) | 98.62, 4451 | 98.63, 4452 (98.6224, 4451.46) |
| 22 | C2 (X=2.4), C3 | 171.5, 1.26e4 | 171.6, 1.261e4 (171.502, 12605) |
| 24, 26, 53? | t2*, certified t* (X=2.4) | 1.224e-7 | 1.223e-7 (1.22396e-7) |
| 28 | margin X=2.4 | 2.45 | 2.44 (2.4479) |
| 29 | C1 sqrt t (X=2.4); C2 sqrt t | 0.0058; 0.034, 0.038 | 0.0059; 0.035, 0.039 |
| 30 | global L2 term (X=3); slope margin (X=2.4) | 0.0400; 0.1242 | 0.0401; 0.1241 |
| 34 | sweep t* | X=3: 3.83, 3.04; X=2.4: 1.08, 0.89 (e-7) | 3.82, 3.03; 1.07, 0.88 |
| 53, 58 | far t* | 1.232e-8 | 1.231e-8 (1.23185e-8) |
| 58 | far C0 sqrt t, C1 sqrt t | 0.0088, 0.036 | 0.0089, 0.037 (0.00880819, 0.0363257) |
| 58 | far physical limit \|h_xx\| (step-4 hand-off) | 8.9e-3 | 8.8e-3 (1/113 = 8.85e-3; it is an allowed bound, so rounds down) |
The same far numbers (1.232e-8, 0.0088, 0.036, 8.9e-3) occur in `writeup/step3_proof/step3_proof.tex` (not edited here). lemmaO_v3.pdf rebuilt.

## 4. To update in step3_proof.tex afterwards
sup|(a,b)| 2.183e-7 -> 2.192e-7 and sup|beta| 3.34e-9 (unchanged; certified 3.33811e-9) at 1.5e-4; 2.154e-7 and 3.21e-9 at 1.25e-4; J_2 1.87e-12 -> 1.88e-12 and slopes 12.1e-19, 2.68e-19 -> 12.2e-19, 2.69e-19; M0core margin 5.87 -> 5.86; the Lemma O items of section 3; table file names fin10 -> fin11; the 1.25e-4 remark numbers above; modulation constants BETA1 66.5 -> 66.53, CAB0 3.708 -> 3.7083, QAB 51.3 -> 51.34, PHIAB 3.84 -> 3.8404 (step3_proof.tex 3.4/5.1 if quoted).

## 5. Far-curve jet hand-off thresholds of lemmaO_v3.tex (Revision section): derived from quoted or certified values?
Checked: none of them was derived from the under-quoted decimals.
* |h''| <= 8.9e-3: this is M/L = 1/113 = 8.8496e-3 rounded to nearest, i.e. 0.5 percent ABOVE the allowed value (an allowed bound must round down). **Changed to 8.8e-3** (lemmaO_v3.tex:58; the same 8.9e-3 occurs in step3_proof.tex and in step 4's hand-off requirement list: use 8.8e-3).
* |h'''| <= 7.0e-5 (exact M/L^2 = 7.83e-5) and |h''''| <= 6.2e-7 (M/L^3 = 6.93e-7): already below the exact limits, safe, unchanged.
* slope <= 0.21: from the H2 condition 1/4 - C1 sqrt(t) with the certified C1 sqrt(t) = 0.0363257, i.e. 0.21367; the quote 0.036 would give 0.214, so 0.21 is safe either way, unchanged.
* K3 slack 1.008 (certified 1.00833), margin 2.008 (certified 2.0085), C2 sqrt(t) 0.353 (0.352804): quoted on the safe side, unchanged. The margins 2.31 and 2.44 are computed from the certified t* (2.4 case corrected from 2.45).
So the only threshold that changes for step 4 is |h''| <= 8.8e-3 (lambda_1 units).
