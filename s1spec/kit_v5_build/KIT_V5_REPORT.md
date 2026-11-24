# KIT_V5_REPORT (Nov 24 2026): reproduction kits v5 = v4 + Windows fixes

| kit | zip | bytes | sha256 | files |
|---|---|---|---|---|
| kit_step3 v5 | `s1spec/kit_step3_v5.zip` (entries under `kit_step3_v5/`) | 16264692 | `f0ea28743d340af8a440b4280e7ccc7ed71b5306635c577da0d568b619ee413e` | 549 |
| kit_step4b v5 | `s1spec/kit_step4b_v5.zip` (flat, as v4) | 43917654 | `49c4b7227069090d0b3573e6c6a10d42d77e94a437b3306101c3e303dbdf8ea2` | 1612 |

Build: `kit_step3_v5_build/build_v5.py`, `kit_step4b_v5_build/build_v5.py`, shared `kit_v5_build/v5common.py`. Each copies the v4 TREE (asserted == v4 zip, sha 0f899eb2 / ab47b541), applies string patches (each must match exactly once; `[kit patch v5, winfix N]` comments), regenerates the changed pins in `inputs/MANIFEST.sha256` (then all lines verified), asserts every other file byte-identical to v4, parses every .py, writes a deterministic zip (rebuilt twice: same sha, same tree). Full change list with before/after sha256 and the Windows path audit: `CHANGES_v5.md` in each kit and build dir. v4 trees, zips and build dirs unchanged (hashes checked before/after). Proof tree not touched.

## Fixes
1. `kit_guard.py` (both): CommandLineToArgvW split of the joined Windows Popen command line (code verbatim from the tested run copies).
2. `chart4_shell.py` (both): `resource` import guarded; Windows PeakWorkingSetSize via psapi, else -1 (RSS report only, timing sidecar). Pins regenerated: kit3 MANIFEST `kit_guard.py` 970e4e50->c977ce51, `chart4_shell.py` 81418ebb->bd545115; kit4b MANIFEST `chart4_shell.py` 81418ebb->bd545115 (kit_guard is not in kit4b's manifest). The env shim is no longer needed. `c2_t3d.py` also imports `resource` but is never executed (hashed input only): left unchanged.
3. CRLF, compare side (proof-tree scripts and outputs unchanged): kit3 `check.py` cmp_text (flux texts: bytes or modulo CRLF); kit4b `flux_v38.py` (tested run-copy patch: texts + gen_texts.json output hashes modulo CRLF), `verify_v33.py` (non-JSON outputs equal modulo CRLF -> EXACT; removes the Windows WARN on the 5 .tex), `check.py` labels.
4. `verify_v4c_runs_with_exceptions.py`: relpath -> '/' at :43 (tested) and the output keys :44, :51. Audit of all other relpath/sep/split('/') uses in both kits: none needs a change (listed in CHANGES_v5.md).
5. READMEs: v5 note, Windows 11 PowerShell instructions (venv, CPython 3.13.7 tested, python-flint 0.9.0), expected WARNs; kit4b: the "outputs vs proof tree ... different" line is INFO.

## Tests (WSL)
- kit_step3 v5 full run `s1spec/kit3v5_run_wsl1/` (1696 s): **80 PASS, 0 FAIL, 1 WARN, 0 SKIP**, status lines identical to kit3v4_run_wsl5 (WARN = ibp float noise 5.7e-16); input manifest ok with the new pins; guard clean.
- kit_step4b v5 partial `s1spec/kit4bv5_test_wsl1/` (copy of kit4bv4_run_wsl1 + v5 overlay; TEST_NOTES_v5.md): parse ok; inputs on pristine v5 extraction 1351/0 problems; v38 stage rerun 5m21s, flux_v38 8/8, v4c pins 6/6; **check.py TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]** (= v4; WARN iface_compare_v38.json by design).
- Unit tests `kit_v5_build/tests/unit_v5.py`: 9 PASS (CRLF branches, rss() without `resource`).

## Not tested
- No Windows run of either v5 kit. The Windows-only branches (guard split, psapi RSS, CRLF compares, '\\' relpath keys) were exercised on Windows only in their v4 run-copy form (guard, flux_v38, verify_v4c:43 verbatim; psapi code from the env shim, now inside chart4_shell.py and reporting MB); the CRLF compare branches are unit-tested on WSL with CRLF files; verify_v4c :44/:51 key changes are untested on Windows (no effect on POSIX).
- kit4b v5 was not run end to end (stages before v38 not rerun; they are unchanged from v4 except kit_guard.py and chart4_shell.py, both exercised by the v38 rerun).
