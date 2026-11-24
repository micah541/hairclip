# kit_step3 v4 test run (WSL, Nov 20 2026): 80 PASS, 0 FAIL, 1 WARN, 0 SKIP
- zip `kit_step3_v4.zip` sha256 0f899eb2de74cf51e25aafb9ac0c401343279a9a34f1dbf0dc7a2286809ce853 (16296685 bytes), extracted into this new empty directory; `PYTHONHASHSEED=0 nice -n 19 python3 -B run_all.py` then `python3 -B check.py`. Python 3.14.4, flint 0.9.0, numpy 2.5.3, sympy 1.14.0. Wall finished in 1675s (all parts, one process at a time).
- WARN: ibp final_5e-4.json (unchanged kit-v2 part): numbers differ at 5.7e-16 relative (float noise), within rule C's 1e-9.
- Every v4 output (verifier rate tables, mod_sups, cert strip_cubic / A2 budget, n1, base chain, R3, 19 t2_A files, Lemma U/zone III passes, 8 flux JSON + pkl + 6 texts, G1/U1/GS_boot) equals the proof tree's reference under rule C (provenance normalised).
- Write guard: no step tried to write / chdir / start a child outside the kit; 0 files read outside the kit.
- Outside-write check: `find /home/claude -newer START_MARKER -type f` lists nothing outside this run directory apart from my own build/report files (kit builders, zips, p2/v38 log); the 29-stage v38 manifest pins and the closed manifests e01cb4da / 8bc47e83 verify.
- Earlier debugging runs (kept, not the result): kit3v4_run_wsl1..wsl4 (zips 7ad217fb, f730ad6a, 17e4850e, ad294522); fixes recorded in step4/p2/v38/V38_LOG.md (missing flux inputs, copy-onto-itself, checker rules, unnecessary headers on verifier files).

## check.py SUMMARY
```
SUMMARY: 80 PASS, 0 FAIL, 1 WARN, 0 SKIP
```

## check.py output
```
REPRODUCTION CHECK kit_step3 v4  (python 3.14.4, python-flint 0.9.0, Linux-6.18.58.1-microsoft-standard-WSL2-x86_64-with-glibc2.43)

PASS  verifier: tail_cert PASS, log identical to the reference
PASS  verifier: rate_table_fin13v4c_a0.00015.json (closure alpha1 = 1.5e-4)
        equal to the reference (rule C)
PASS  verifier closes at alpha1 = 1.5e-4: tail certified, rho > 0, every box-hypothesis ratio <= 1
        peak sqrtL 5.308499e-08, rho 2.679225e-03, max box ratio 0.7891314071814564
PASS  verifier: rate_table_fin13v4c_a0.000125.json (closure alpha1 = 1.25e-4)
        equal to the reference (rule C)
PASS  verifier closes at alpha1 = 1.25e-4: tail certified, rho > 0, every box-hypothesis ratio <= 1
        peak sqrtL 5.214218e-08, rho 2.850927e-03, max box ratio 0.7891104250359056
PASS  verifier: mod_sups final_alpha_points_5e-8.json (sup|(a,b)|, sup|beta|, J, rho)
        equal to the reference (rule C)
PASS  verifier: sup|(a,b)| = 2.14039e-7, sup|beta| = 2.35488e-9 at alpha1 = 1.5e-4 (v38 values)
        2.14039e-07, 2.35488e-09
PASS  cert: strip_cubic_v4c.json
        equal to the reference (rule C)
PASS  cert: e0_A2_budget_v4c.json
        equal to the reference (rule C)
PASS  n1: n1_chart_fin13v4c.json
        equal to the reference (rule C)
PASS  n1: ramp e < 2^-19 and sup chi_add_0 <= 2e-7
        e 4.388e-10, sup_add0 2.086e-08
PASS  SA base: lemmaS_F10p_M0.52.json
        equal to the reference (rule C)
PASS  SA base: cap_F10p.json
        equal to the reference (rule C)
PASS  SA base: lemmaS_rigcap_M0.52.json
        equal to the reference (rule C)
PASS  SA base: weights_exact.json
        equal to the reference (rule C)
PASS  SA base: assembly_arb.json
        equal to the reference (rule C)
PASS  SA R3: r3_forcing.json
        equal to the reference (rule C)
PASS  SA R3: r3_symcheck.json
        equal to the reference (rule C)
PASS  SA t2_A: hjw_rep.json
        equal to the reference (rule C)
PASS  SA t2_A: ic_design.json
        equal to the reference (rule C)
PASS  SA t2_A: ic_lo.json
        equal to the reference (rule C)
PASS  SA t2_A: ic_hi.json
        equal to the reference (rule C)
PASS  SA t2_A: ic_probe.json
        equal to the reference (rule C)
PASS  SA t2_A: a3_signs.json
        equal to the reference (rule C)
PASS  SA t2_A: weights_t2.json
        equal to the reference (rule C)
PASS  SA t2_A: weights_t2_fb.json
        equal to the reference (rule C)
PASS  SA t2_A: ramp_rep.json
        equal to the reference (rule C)
PASS  SA t2_A: ramp_rep_pm3ibp.json
        equal to the reference (rule C)
PASS  SA t2_A: outer_L2.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_design_N2ref.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_design_N2ref4.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_design_N20.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_lo_N2ref.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_hi_N2ref.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_probe_N2ref.json
        equal to the reference (rule C)
PASS  SA t2_A: asm_design_pm3ibp.json
        equal to the reference (rule C)
PASS  SA t2_A: A_interface.json
        equal to the reference (rule C)
PASS  SA: A_interface rule_reproduces_EA_ref (tight, v31)
PASS  SA: asm_design_N2ref tight: failing claims are exactly those of the reference (v3.1 page-1 quotes and the v3.1 S0_B quote; not v3.3+ claims)
        11: ['E_0 (v3.1 page 1)', 'E_B0 (y_E=50) (v3.1 page 1)', 'I3 const (v3.1 page 1.84e-21)']
PASS  SA: asm_design_N2ref4 tight: failing claims are exactly those of the reference (v3.1 page-1 quotes and the v3.1 S0_B quote; not v3.3+ claims)
        11: ['E_0 (v3.1 page 1)', 'E_B0 (y_E=50) (v3.1 page 1)', 'I3 const (v3.1 page 1.84e-21)']
PASS  U: u1_ghat.json
        equal to the reference (rule C)
PASS  U: u2_K.json
        equal to the reference (rule C)
PASS  U: u34_primary.json
        equal to the reference (rule C)
PASS  U: U_cert_v36_pass1.json
        equal to the reference (rule C)
PASS  U: U_cert_for_zone3_v36_pass1.json
        equal to the reference (rule C)
PASS  U: U_cert_v36_pass2.json
        equal to the reference (rule C)
PASS  U: U_cert_for_zone3_v36.json
        equal to the reference (rule C)
PASS  U: U_cert_v36_pass3.json
        equal to the reference (rule C)
PASS  zone III: zone3_interface_pass1.json
        equal to the reference (rule C)
PASS  zone III: zone3_interface.json
        equal to the reference (rule C)
PASS  zone III: ramp_U.json (an input of chain v37, the v35 byte copy; not regenerated)
        byte-identical to the reference
PASS  U: fixed point, pass 3 == pass 2 (rule C)
PASS  U: supersolution K1, K2, K3, N2 <= 4 N2_ref, U-fine requirement
PASS  zone III: H6 min margin > 1
        1.8608671374173327
PASS  flux: chart4_shell.json
        equal to the reference (rule C)
PASS  flux: flux2_cert.json
        equal to the reference (rule C)
PASS  flux: flux2_coef.json
        equal to the reference (rule C)
PASS  flux: flux2_data.json
        equal to the reference (rule C)
PASS  flux: flux2_forms.json
        equal to the reference (rule C)
PASS  flux: flux2_sanity.json
        equal to the reference (rule C)
PASS  flux: gen_texts.json
        equal to the reference (rule C)
PASS  flux: inputs.json
        equal to the reference (rule C)
PASS  flux: flux2_forms.pkl
        byte-identical to the reference
PASS  flux texts: BPhi_continuation_v36.tex
        byte-identical to the reference
PASS  flux texts: absorption_lemma_v2.tex
        byte-identical to the reference
PASS  flux texts: erratum_flux2.tex
        byte-identical to the reference
PASS  flux texts: flux2_lemma_v2.tex
        byte-identical to the reference
PASS  flux texts: lemmaD_k2_sketch.tex
        byte-identical to the reference
PASS  flux texts: numbers_flux2.json
        byte-identical to the reference
PASS  flux: flux2_cert verdict PASS, no failing condition (C8 = U-1 by hand)
        combined margin None
PASS  proofs: G1_check.json
        equal to the reference (rule C)
PASS  proofs: G1_check all checks pass
PASS  proofs: U1_sym.json
        equal to the reference (rule C)
PASS  proofs: U1_sym all checks pass
PASS  proofs: GS_boot_v4c.json
        equal to the reference (rule C)
PASS  proofs: GS_boot_v4c all checks pass
PASS  ibp: verify_symbolic.json
        equal to the reference (rule C)
WARN  ibp: final_5e-4.json
        only numbers differ, max relative difference 5.68e-16 (/final/core/b/gamma_ref: 0.048848138358923356 vs 0.04884813835892333)
PASS  lemmaO: lemmaO_v3_consts.json + far_margin2.json (merged, as the reference)
        equal to the reference (rule C)
PASS  guard: files read outside the kit and outside the Python installation / OS directories
        0
PASS  guard: no step tried to write / chdir / start a child outside the kit

SUMMARY: 80 PASS, 0 FAIL, 1 WARN, 0 SKIP
```
