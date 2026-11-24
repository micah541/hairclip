# kit_step4b v3 first end-to-end run (WSL, Nov 15 2026)
zip 646f2552..., fresh unzip, `nice -n 19 python3 run_all.py --jobs 2`, 13:07-19:31 (6 h 24 min), Python 3.14.4, flint 0.9.0, numpy 2.5.3. No stage failed, no manual intervention.
- check.py: **137 PASS, 0 FAIL, 1 WARN**. All v3.2 claims and all v3.3 addendum claims (sections 1, 5, 7) hold with the proof tree's values (e.g. (H2) f-U_1 x563/x601, (H3) x5.7, cone slope x2.64, H7 x2.30, N1 x9.59, f12 326 checks, M1 x1.1847, K6 x1.19, H6 x1.861, FLUX x3.365/x5.956).
- WARN: f12_sliver.json "structure differs": only its `outputs` dict is keyed by absolute paths (kit vs proof tree root); provenance, no number involved. Checker fix: normalise keys.
- v3.3 outputs: 0 byte-identical / 63 equal with args+timing dropped / 109 numbers differ / 1 structure (above) / 0 missing of 173. Largest "relative difference 1.0" (h3_arb.json) is a fidelity residual 0.0 vs -2.2e-16; real values differ at ~1e-8 relative (regenerated v3.2 chain).
- v3.2 part: 30 byte-identical / 468 numerically different / 0 missing of 498 (B8 now runs inside the kit).
- Proof tree check (find -newermt run start): only my own CLAUDE.md and the two exception files changed; the kit wrote nothing outside its directory.
