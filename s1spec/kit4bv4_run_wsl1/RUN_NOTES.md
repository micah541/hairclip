# kit_step4b v4 test run (WSL, Nov 20 2026): TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]
- zip `kit_step4b_v4.zip` sha256 ab47b5416e8a6d58d5267e2a93c154bf979809005de0abb4e992dd989df7f0fa (43912237 bytes), extracted into this new empty directory; `PYTHONHASHSEED=0 nice -n 19 python3 -B run_all.py --jobs 2` (check.py run by run_all at the end), 15:08-21:29, 6 h 20 min, no intervention. Python 3.14.4, flint 0.9.0, numpy 2.5.3.
- All 40 run_all stages ok (v3.2 chain, 21 v3.3 stages, stage v38 = the 29-stage v38 driver). v38: 29 stages, 16 carried by hash from this run's v3.3 manifest, all_items_ok; every v38 result item ok; f12_sliver 339/339; flux combined x4.8045; H7, far H4 x15.59, N1, C1_ke (x287) hold.
- WARN (benign, by design): `iface_compare_v38.json` has 2 fewer checks than the proof tree's, because in the kit C1_ke/GS_tau1 read this run's regenerated step-4 data, so their equality with the proof-tree references is recorded, not asserted (kit patch). Their numbers differ from the proof tree's by <= 1.0e-4 relative (C1_ke) and ~1e-6 (GS_tau1).
- Write guard (kit_guard.py, every process): 0 attempts to write / chdir / start a child outside the kit; 0 files read outside the kit.
- Outside-write check: `find /home/claude -newer START_MARKER -type f` outside this directory lists only files I wrote meanwhile (kit3v4_run_wsl5/RUN_NOTES.md, p2/out_v38/REPIN_REPORT.md, p2/v38/V38_LOG.md, p2/v38/repin_report_v38*.py). Proof-tree manifests after the run: v38 be909bbf (992 pins, 0 mismatches); closed v3.2 e01cb4da verifies with the B8 exception only; v3.3 8bc47e83 verifies apart from the known F34 constants_needed_v34.json pin (documented in REPIN_REPORT).

## check.py SUMMARY
```
======================================== SUMMARY (paste this) ========================================
TOTAL PASS   [159 PASS, 0 FAIL, 1 WARN]
env: Python 3.14.4 | Linux-6.18.58.1-microsoft-standard-WSL2-x86_64-with-glibc2.43 | flint 0.9.0 numpy 2.5.3 | jobs 2
outputs vs proof tree: byte-identical / timing-fields-only / provenance-only / different / missing: 30 / 0 / 0 / 468 / 0  of 498 reference files
  T7: closure constant rho_c (min over certificates) >= 0.0471: yours 0.04711 (ours 0.04713) margin x1
  switch (A): Y(tau_a) <= typed Y0 = 6.21e-13 (Ybar = 1e-12): yours 6.203e-13 (ours 6.203e-13) margin x1
  (Y): comparison level on the direct residuals <= 3.22e-14: yours 3.214e-14 (ours 3.214e-14) margin x1
  (F) item 3 E3 (window, binding): Ybar limit >= 1e-12: yours 1.576e-11 (ours 1.576e-11) margin x15.8
  (F) item 4 (K4 at the good time): Ybar limit >= 1e-12: yours 2.891e-11 (ours 2.891e-11) margin x28.9
  (F) Ymax (min over items) >= the 1.7e-12 box level of N3: yours 1.576e-11 (ours 1.576e-11) margin x9.27
  boxes: far core x^2|e'''| <= b2 = 1e-2: yours 0.008335 (ours 0.008335) margin x1.2
  boxes: sqrt(F0) ratio <= fac = 1.7: yours 1 (ours 1) margin x1.7
  (H1) sqrt(Lambda) over slab 180 incl. C_Lambda Ybar <= 5e-8: yours 2.284e-08 (ours 2.284e-08) margin x2.19
  (H1) |alpha| over slab 180 + C_alpha Ybar <= 1.5e-4: yours 0.0001174 (ours 0.0001174) margin x1.28
  (H3) window energy E3 (base, slab 180) <= threshold: yours 9.706e-07 (ours 9.706e-07) margin x10.3
  KCH (chart factor, computed) <= 1.188: yours 1.188 (ours 1.188) margin x1
  (H5) cone K4 <= 0.12 (r_c in [200, 380], base + error at Ybar): yours 0.0123 (ours 0.0123) margin x9.76
  (H5) Lemma O H4: worst cylinder of the final window, int (h0^(5))^2 / M^2 <= 1: yours 0.2392 (ours 0.2392) margin x4.18
  v3.2: closure margin cone_K2 >= 1: yours 1.688 (ours 1.688) margin x1.69
  v3.3 (H2f) j=0: with-error value <= 1e-06: yours 3.136e-07 (ours 3.136e-07) margin x3.19
  v3.3 (H2f) j=0 [margin vs addendum x3.188]: yours x3.189 (ours x3.188) margin 
  v3.3 (H2f) j=1: with-error value <= 5e-06: yours 1.003e-06 (ours 1.003e-06) margin x4.98
  v3.3 (H2f) j=1 [margin vs addendum x4.983]: yours x4.983 (ours x4.983) margin 
  v3.3 (H2f) j=2: with-error value <= 1e-05: yours 2.435e-06 (ours 2.435e-06) margin x4.11
  v3.3 (H2f) j=2 [margin vs addendum x4.106]: yours x4.107 (ours x4.106) margin 
  v3.3 tightened (H2), w' j=0: with-error value <= 1e-08: yours 8.013e-11 (ours 8.013e-11) margin x125
  v3.3 tightened (H2) w' j=0 [margin vs addendum x124.78]: yours x124.8 (ours x124.78) margin 
  v3.3 tightened (H2), w' j=1: with-error value <= 1e-06: yours 4.706e-10 (ours 4.706e-10) margin x2.12e+03
  v3.3 tightened (H2) w' j=1 [margin vs addendum x2124.9]: yours x2125 (ours x2124.9) margin 
  v3.3 tightened (H2), w' j=2: with-error value <= 0.001: yours 1.559e-06 (ours 1.559e-06) margin x642
  v3.3 tightened (H2) w' j=2 [margin vs addendum x641.5]: yours x641.5 (ours x641.5) margin 
  v3.3 tightened (H2), f - U_1 (B0) j=0: bound <= 1e-08: yours 3.657e-10 (ours 3.657e-10) margin x27.3
  v3.3 tightened (H2) f-U_1 j=0 [margin vs addendum x27.34]: yours x27.34 (ours x27.34) margin 
  v3.3 tightened (H2), f - U_1 (B0) j=1: bound <= 1e-06: yours 1.777e-09 (ours 1.777e-09) margin x563
  v3.3 tightened (H2) f-U_1 j=1 [margin vs addendum x562.8]: yours x562.8 (ours x562.8) margin 
  v3.3 tightened (H2), f - U_1 (B0) j=2: bound <= 0.001: yours 1.664e-06 (ours 1.664e-06) margin x601
  v3.3 tightened (H2) f-U_1 j=2 [margin vs addendum x600.9]: yours x601 (ours x600.9) margin 
  v3.3 (H3) window j=0: with-error value <= 2e-07: yours 3.509e-08 (ours 3.509e-08) margin x5.7
  v3.3 (H3) window j=0 [margin vs addendum x5.7]: yours x5.7 (ours x5.7) margin 
  v3.3 (H3) window j=1: with-error value <= 1e-05: yours 1.262e-07 (ours 1.262e-07) margin x79.3
  v3.3 (H3) window j=1 [margin vs addendum x79.26]: yours x79.27 (ours x79.26) margin 
  v3.3 (H3) window j=2: with-error value <= 1e-05: yours 3.772e-07 (ours 3.772e-07) margin x26.5
  v3.3 (H3) window j=2 [margin vs addendum x26.5]: yours x26.51 (ours x26.5) margin 
  v3.3 (H3) window j=3: with-error value <= 1e-05: yours 1.282e-06 (ours 1.282e-06) margin x7.8
  v3.3 (H3) window j=3 [margin vs addendum x7.799]: yours x7.8 (ours x7.799) margin 
  v3.3 (H3) Ybar limit, binding (E_3) >= 1e-12: yours 1.576e-11 (ours 1.576e-11) margin x15.8
  v3.3 (H3) Ybar limit, binding (E_3): margin >= 15.75 (addendum x15.75): yours 15.76 (ours 15.76) margin x1
  v3.3 far-family Lemma O H4 at every centre >= 490: min margin >= 15.58 (addendum x15.58): yours 15.59 (ours 15.59) margin x1
  v3.3 fixed-L cone slope datum (H5) <= 0.08: yours 0.03029 (ours 0.03029) margin x2.64
  v3.3 fixed-L cone slope datum (H5) [margin vs addendum x2.64]: yours x2.641 (ours x2.64) margin 
  v3.3 fixed-L cone sup|h_0''| = K2 (H5) <= 0.03: yours 0.009876 (ours 0.009876) margin x3.04
  v3.3 fixed-L cone sup|h_0''| = K2 (H5) [margin vs addendum x3.03]: yours x3.038 (ours x3.03) margin 
  v3.3 fixed-L cone ||h_0'''||_L2 = N3 (H5) <= 0.004: yours 0.0001461 (ours 0.0001461) margin x27.4
  v3.3 fixed-L cone ||h_0'''||_L2 = N3 (H5) [margin vs addendum x27.3]: yours x27.38 (ours x27.3) margin 
  v3.3 fixed-L cone H3 jet K3 <= 0.12: yours 5.964e-05 (ours 5.964e-05) margin x2.01e+03
  v3.3 fixed-L cone H3 jet K3 [margin vs addendum x2012]: yours x2012 (ours x2012) margin 
  v3.3 fixed-L cone H3 jet K4 <= 0.12: yours 0.0008363 (ours 0.0008363) margin x143
  v3.3 fixed-L cone H3 jet K4 [margin vs addendum x143]: yours x143.5 (ours x143) margin 
  v3.3 fixed-L cone Lemma O H4: int (h_0^(5))^2 <= 0.0144: yours 1.286e-05 (ours 1.286e-05) margin x1.12e+03
  v3.3 fixed-L cone H4 [margin vs addendum x1119]: yours x1119 (ours x1119) margin 
  v3.3 fixed-L cone Lemma O H1 (graph): sup|h_0| <= H - C_0 sqrt(t_1): yours 0 (ours ) margin xnan
  v3.3 fixed-L cone H7: ratio <= needed: yours 1.6e-07 (ours 1.6e-07) margin x1
  v3.3 fixed-L cone H7 t* [margin vs addendum x2.3]: yours x2.302 (ours x2.3) margin 
  v3.3 N1 chart: |e| (ramp) < 2^-19: yours 4.388e-10 (ours 4.388e-10) margin x4.35e+03
  v3.3 N1 chart: sup|chi^add_0| on the ramp <= 2e-7: yours 2.086e-08 (ours 2.086e-08) margin x9.59
  v3.3 f12_sliver: all checks pass and no quoted number fails: yours 326 checks; failing quoted numbers: [] (ours 326 checks) margin 
  v3.3 f12 ramp M1 (A) >= 1.1847 (quoted x1.1847): yours 1.185 (ours 1.185) margin x1
  v3.3 f12 K6 >= 1.189: yours 1.19 (ours 1.19) margin x1
  v3.3 f12 H6 min margin >= 1.860: yours 1.861 (ours 1.861) margin x1
  v3.3 f12 flux FLUX0/K0 (beta 2, 4 N2_ref) >= 3.365: yours 3.365 (ours 3.365) margin x1
  v3.3 f12 flux FLUX2/K2 (beta 2, 4 N2_ref) >= 5.955: yours 5.956 (ours 5.956) margin x1
  v3.3 f12 M1' = 28 sI2' <= 0.4328 (what (A) reads): yours 0.4328 (ours 0.4328) margin x1
  v38 result item H7_fixedL: ok: yours True (ours True) margin 
  v38 result item notes: ok: yours True (ours True) margin 
  v38 result item M9_H4_far: ok: yours True (ours True) margin 
  v38 result item N1_chart: ok: yours True (ours True) margin 
  v38 result item F12_sliver: ok: yours True (ours True) margin 
  v38 result item flux_v37: ok: yours True (ours True) margin 
  v38 result item v4c_runs_pins: ok: yours True (ours True) margin 
  v38 result item interface_proofs: ok: yours True (ours True) margin 
  v38 H7: (T-t1)/L^4 <= 1.6001e-7 (fin13v4c abar): yours 1.6000013698507732e-07 (ours 1.6000013698507732e-07) margin x1
  v38 far-family H4 min margin >= 1: yours 15.587791462957798 (ours 15.58779119601459) margin x15.59
  v38 N1 ramp e < 2^-19: yours 4.3881233015929833e-10 (ours 4.3881233015929833e-10) margin x4347
  v38 f12_sliver: all checks pass: yours 339/339 (ours 339/339) margin 
  v38 f12 ramp M1 (A) > 1: yours 1.1844642357292894 (ours 1.1844642357292894) margin x1.184
  v38 f12 K6 > 1: yours 1.1884203179916526 (ours 1.1884203179916526) margin x1.188
  v38 flux v2 combined pointwise margin > 1: yours 4.804534359260718 (ours 4.804534359260718) margin x4.805
  v38 C1_ke sup|k~ e| < 1/2: yours 0.0017428461585192827 (ours 0.0017428478070338894) margin x286.9
  v38 outputs vs proof tree: byte-identical / equal with args+timing dropped / numbers differ / structure differs / missing: yours 7 / 42 / 3 / 1 / 0 (ours of 53 reference files) margin 
  v38 outputs with differing numbers: largest relative difference: yours 1.01e-04 (ours ) margin 
  v38 outputs whose structure/text differs from the proof tree: yours 1: p2/out_v38/iface/iface_compare_v38.json (ours ) margin 
  guard: files read outside the kit and outside the Python installation / OS directories: yours 0 (ours ) margin 
  guard: no process tried to write / chdir / start a child outside the kit: yours 0 (ours ) margin 
=======================================================================================================
```
