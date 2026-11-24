# CHANGES_v5: kit_step4b v4 -> v5 (Windows fixes)

Base: `kit_step4b_v4.zip` sha256 `ab47b5416e8a6d58d5267e2a93c154bf979809005de0abb4e992dd989df7f0fa` (tree verified equal to the zip). Source of the fixes: the kit4bv4 Windows run (`s1spec/kit4bv4_run_win1/RUN_NOTES.md`, `*.diff`; run copy `C:\Users\micahw\kitruns\kit4bv4_win1`). Every changed code line carries a `[kit patch v5 ...]` comment (README / version strings excepted). Fix 2 is made in the script (not by an environment shim), so its pin in `inputs/MANIFEST.sha256` is regenerated (table below); no other pin of the kit refers to a changed file (flux_v38.py, verify_v4c_runs_with_exceptions.py, verify_v33.py, check.py, run_all.py, kit_guard.py, README.md are not in the manifest). Choice for CRLF (fix 3): compare side (flux_v38.py, verify_v33.py), so the proof-tree scripts (gen_texts.py) and their outputs are unchanged. All other files are byte-identical to v4 (checked by the builder).

## Changed kit files (sha256 before = v4, after = v5)

| file | v4 sha256 | v5 sha256 | patches | why |
|---|---|---|---|---|
| `kit_guard.py` | `970e4e50d7de034944c49143e5f4fdf95532b0b08778c305a38e7c90dc6cf390` | `c977ce5143896a417db18538ff2efee93cd8dd4b4261c8a812948a3720c5814c` | docstring, winfix1-def, winfix1-split | winfix 1: on Windows the subprocess.Popen audit event carries ONE joined command line, so the guard denied the kit's own guarded Python children; split it with CommandLineToArgvW (code of the tested Windows run copies). No effect on POSIX. |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | `81418ebb2e200731e3b219f6065bd3edd002ac64a5518c08c6d8a7fe7933b282` | `bd5451152c00086d7deec2242f752879eda7853554746dd6700be4acc8a40820` | winfix2-import, winfix2-rss | winfix 2: `resource` is Unix-only (RSS report only, written to the timing sidecar / log, never compared): import guarded, Windows peak working set via psapi, else -1. Replaces the env shim of the kit4bv4 Windows run. Certified computation unchanged. |
| `step4/p2/v38/flux_v38.py` | `05af4627d35e88f81799a22b1579ae2f528c5d87a7ddcd77c91fc27f7cfd23b5` | `6f096f1192dfaf3e833824a3981c9a1ff615d6c927933ede80909845c38e15a0` | winfix3-lf, winfix3-gen_texts, winfix3-texts | winfix 3: the 6 rerun flux texts and the output hashes recorded in gen_texts.json are compared modulo CRLF (Python text mode writes CRLF on Windows); the tested run-copy patch (kit4bv4_run_win1/flux_v38_winfix.diff). Compare side: gen_texts.py (proof-tree script) unchanged. Identical result on POSIX. |
| `step4/p2/v38/verify_v4c_runs_with_exceptions.py` | `59a4064f4bc4e28ea701f4cf0919862d25e49ff9bf3f32a38db8ec3b022262cc` | `33b96308aff25741146fc7a031b9fc41772e8b309bd96f0f6be9686950be1315` | winfix4-affected, winfix4-rowkey, winfix4-comment, winfix4-eqkey | winfix 4: the exception record lists affected records as runs/<run>/<file> with '/'; os.path.relpath gives backslashes on Windows, so the documented exception was not recognised (stage 28 failed). relpath joined with '/' (tested line) and the two output dict keys made separator-independent. Identical on POSIX. |
| `verify_v33.py` | `326546c14fe92ade5a0c925e31ae6dd2edb572422cb924acd19274c835213acd` | `7d22961294b703603820228af7fb2ade4b8bc2709c1c73396d89c027fdb3e8c1` | winfix3-exact, winfix3-helper, doc | winfix 3 (checker): a non-JSON output (the 5 flux .tex) equal to the reference modulo CRLF is EXACT (was DIFFERENT -> WARN in the Windows run). |
| `check.py` | `e0e4d2c09e04bbd3063aa657f1bfb48428489cc7e7c0e74ec6a7665a69c121d2` | `6811e6d8151c3be05061fce0cdafebee0c7aebfd613d6067cd1d0d1b4f57228d` | version, label, label, version | version strings v4 -> v5; INFO labels say that text outputs count as byte-identical modulo CRLF (verify_v33.py change). |
| `run_all.py` | `50cfd5c6f48fa6a492bca3efe3669f94508e949cd0d8574abe0e723fce835f65` | `a2824e0b8cea14471963c4772e32c24af64c2cfb46fc0187e64b74255bab3dbb` | version | version string (docstring only). |
| `README.md` | `9483475f21edae12249466d46971deddf802c0487a147a239bf8fb0abf4c1d4e` | `deecc9a67b54013a9df91e97c6d1622c36c0a9375fe6b6c1ef442ec5b8bb5d72` | v5-head | winfix 5: README: v5 note, Windows instructions (venv, Python 3.13.7 tested, flint 0.9.0), expected WARN, the "outputs vs proof tree ... different" line is INFO. |
| `inputs/MANIFEST.sha256` | `3e33e549f9b2457d7c09fc9c30f1e924602fddcda2f35864e2b247b80843517f` | `776bce86656d0b52b1f97a8284aaa6b07bc778aa47c84447573b49c9f57d4b4a` | pins | pins of the changed kit files regenerated (kit v5) |

## Pins regenerated in `inputs/MANIFEST.sha256`

| pinned file | old pin (v4) | new pin (v5) | why |
|---|---|---|---|
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | `81418ebb2e200731e3b219f6065bd3edd002ac64a5518c08c6d8a7fe7933b282` | `bd5451152c00086d7deec2242f752879eda7853554746dd6700be4acc8a40820` | winfix 2: `resource` is Unix-only (RSS report only, written to the timing sidecar / log, never compared): import guarded, Windows peak working set via psapi, else -1. Replaces the env shim of the kit4bv4 Windows run. Certified computation unchanged. |

## Added files

- `CHANGES_v5.md` (this file)

## Windows path audit (fix 4), kit_step4b

Searched every `.py` of the kit for `os.path.relpath`, `os.sep`, `.split('/')`, `startswith(... + '/')` and backslash literals.
- FIXED `step4/p2/v38/verify_v4c_runs_with_exceptions.py:43` (the run-copy patch: membership test against the exception record), `:44` and `:51` (dict keys of its output: now '/' on every OS, so `v4c_runs_pins.json` is the same JSON on Windows).
- `step4/p2/v38/verify_v4c_runs_with_exceptions.py:50` `os.path.join(b, os.path.relpath(...))`: a file path, correct on both.
- `step4/p2/v33/final_tail_v33.py:49-50`, `step4/p2/v38/final_tail_v38.py:46-47`: relpath already followed by `replace('\\', '/')`; `startswith(S4 + '/')`, `startswith(OV + '/')`: S4/OV are built from `_S1` (forward slashes) and the stage paths the drivers compare are built the same way (the Windows v4 run passed every stage).
- `verify_v33.py:65/67`, `verify_outputs.py:29/33`, `run_all.py:191/380/382`: relpath with `replace(os.sep, '/')` / manifest paths split on '/' then `os.path.join` / `os.sep`: correct on both.
- `step4/p2/v38/iface_compare_v38.py:67`: relpath in a check NAME only (text); not compared as a key (iface_compare_v38.json is the known benign WARN in both runs).
- step-4 v3.2 provenance (`step4/t7/se7.py:68-116`, `merge.py:125/139`, `nlterms7.py:596`, `provenance.py:22`, `t6c/runc.py:34`, `t7/run_all.py:63`): relpath written into args/provenance fields of outputs; check.py / verify_outputs.py drop args and normalise separators; proof-tree reproductions, left unchanged (the Windows v4 run passed with them).
- `t1/hotab.py:11`, `t1/cmp_hls.py:12` `f.split('/')`: diagnostics, not run by the kit.
- `writeup/.../flux/code/inputs_manifest.py:23`: relpath is one of three alternatives; the third (`'runs/' + basename`) matches on Windows. Proof-tree script, unchanged.
- `step3nl_v4c/v4cert/c2_t3d.py` imports `resource` but is NOT executed by the kit (pinned input hashed by inputs_manifest.py); left byte-identical.
- `kit_guard.py`: realpath + `os.sep` comparisons, correct; command-line split is fix 1.

