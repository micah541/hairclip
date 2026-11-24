# Reproduction kit step 3, v5 (verifier fin13v4c, (A) chain v37, flux lemma v2, interface / proof checks)

**v5 = v4 + the Windows fixes** of the first Windows run (Nov 24 2026; `CHANGES_v5.md`): the write guard splits the Windows command line, `chart4_shell.py` runs without the Unix-only `resource` module, and check.py compares the generated flux texts modulo CRLF line endings. The certified computations and the reference outputs are those of v4.

Recomputes the step-3 certificates that the step-4 v38 manifest (`reference/step4/FINAL_manifest_v38.json`) pins, and compares every output with the
proof tree's (`reference/`). Plain Python + python-flint (arb), one process at a time, no compiler, no network.

```
python -m venv kitenv
kitenv/bin/python -m pip install -r requirements.txt        # Windows: kitenv\Scripts\python.exe
kitenv/bin/python run_all.py                                 # all parts, about 30-35 min
kitenv/bin/python check.py                                   # paste the SUMMARY (and any FAIL/WARN lines) back
```

**Windows 11 (PowerShell).** Tested with CPython 3.13.7, python-flint 0.9.0, numpy 2.5.3, sympy 1.14.0, mpmath 1.3.0 (kit v4 + these fixes). Unzip to a short path outside OneDrive (e.g. `C:\kitruns\kit3v5`), then:
```powershell
cd C:\kitruns\kit3v5\kit_step3_v5
py -3.13 -m venv kitenv
.\kitenv\Scripts\python.exe -m pip install --upgrade pip
.\kitenv\Scripts\python.exe -m pip install -r requirements.txt
.\kitenv\Scripts\python.exe run_all.py                  # about 35-40 min
.\kitenv\Scripts\python.exe check.py
```
Expected on Windows: every certified result PASS; a few WARN for last-bit float differences (relative <= 1e-15) in diagnostic fields (tau grid, gamma_ref) are platform noise, not failures. Flux texts written with CRLF line endings PASS as "identical modulo CRLF".

Parts (`python run_all.py --list`): smoke, verifier (tail_cert, closure 1.5e-4 and 1.25e-4, mod_sups: about 12 min), cert (strip cubic D-3, E0 A^2 budget A-1),
n1, lemmaSA (base chain, R3, t2_A), lemmaU (Lemma U and zone III, then a snapshot `runs/prod2` of this run's chain), flux, proofs (G1, U-1 sympy, GS_boot),
ibp, lemmaO. `python run_all.py lemmaSA lemmaU` runs only those parts (in that order; later parts need the earlier ones).

**Write guard.** Every process the kit starts runs through `kit_guard.py`: writing, deleting or renaming a file, changing directory, or starting a
child process OUTSIDE this folder raises an error (a 2026-11-15 kit run once wrote into the proof tree). `check.py` lists any file read outside the kit.

**What is regenerated, what is pinned.** Regenerated: everything above. Pinned inputs (shipped, SHA-256 in `inputs/MANIFEST.sha256`): the verifier's
certified constant files (`step3nl_v4c/v4cert/out`, `v4bcert/out`, `cubic_E0fine...json`, `planB_ramp.json`), the production verifier run
`step3nl_v4c/runs/P` that the chain reads (as in the proof tree; part `verifier` recomputes it and check.py compares), step 4's `b0_dprof.json`,
`ramp_U.json` (kept by the chain, review B3) and the v35 zone III interface that seeds assemble pass 1.
Comparison rule: see `python check.py --help` (provenance fields such as `args`, paths and input hashes are normalised; numbers are compared exactly).
