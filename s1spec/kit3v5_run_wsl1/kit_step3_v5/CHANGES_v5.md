# CHANGES_v5: kit_step3 v4 -> v5 (Windows fixes)

Base: `kit_step3_v4.zip` sha256 `0f899eb2de74cf51e25aafb9ac0c401343279a9a34f1dbf0dc7a2286809ce853` (tree verified equal to the zip). Source of the fixes: the kit3v4 Windows run (`s1spec/kit3v4_run_win1/RUN_NOTES.md`; run copy `C:\Users\micahw\kitruns\kit3v4_win2`). Every changed line carries a `[kit patch v5 ...]` comment, except README/version strings. Choice for CRLF (fix 3): compare side (check.py), so the proof-tree scripts (gen_texts.py) and their outputs stay as in the proof tree. All other files are byte-identical to v4 (checked by the builder).

## Changed kit files (sha256 before = v4, after = v5)

| file | v4 sha256 | v5 sha256 | patches | why |
|---|---|---|---|---|
| `kit_guard.py` | `970e4e50d7de034944c49143e5f4fdf95532b0b08778c305a38e7c90dc6cf390` | `c977ce5143896a417db18538ff2efee93cd8dd4b4261c8a812948a3720c5814c` | docstring, winfix1-def, winfix1-split | winfix 1: on Windows the subprocess.Popen audit event carries ONE joined command line, so the guard denied the kit's own guarded Python children; split it with CommandLineToArgvW (code of the tested Windows run copies). No effect on POSIX. |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | `81418ebb2e200731e3b219f6065bd3edd002ac64a5518c08c6d8a7fe7933b282` | `bd5451152c00086d7deec2242f752879eda7853554746dd6700be4acc8a40820` | winfix2-import, winfix2-rss | winfix 2: `resource` is Unix-only (RSS report only, written to the timing sidecar / log, never compared): import guarded, Windows peak working set via psapi, else -1. Replaces the env shim of the kit4bv4 Windows run. Certified computation unchanged. |
| `check.py` | `9405af602e2519d84acf209fd476196ab910e3e755c11d58df19fe635efdcc0f` | `831648b00d641620cf53c1436546b45f5d1c6f25b9668caf133ee1e16d5f5b03` | version, winfix3-cmp_text, winfix3-use, version | winfix 3: the 6 generated flux texts (5 .tex + numbers_flux2.json) are compared byte for byte OR modulo CRLF (Windows text-mode writes; the 6 FAIL of the kit3v4 Windows run); compare-side fix, the proof-tree script gen_texts.py and its outputs are unchanged. Version strings v4 -> v5. |
| `run_all.py` | `f093b0a6eb107cc514a996a59069d15d1843fc383d31981e25fea097830051ea` | `d293291c708171600a5aed1950e25588e1fb8d9beadca168cb6e44dc5fdaf073` | version | version string v4 -> v5 (docstring only). |
| `README.md` | `df62dead025604985eed3b4cd56c5213dad01bb2da3a9473262b46ffdec7e45e` | `fa7ffebf31e1eb6ea093b9d0d0483c39353c4bd343fbd07a9f988705318c78c5` | v5-title, v5-windows | winfix 5: README: v5 note and Windows instructions (venv, Python 3.13.7 tested, python-flint 0.9.0), expected WARNs. |
| `inputs/MANIFEST.sha256` | `ecbd977e7022155d6d7846858aaa40b184fc9eb412df5261ffa90d3c124e2e82` | `08750160f86a4c16b1b4137c82bf449d804ec3cb4ab54fd13ce13566f3073cba` | pins | pins of the changed kit files regenerated (kit v5) |

## Pins regenerated in `inputs/MANIFEST.sha256`

| pinned file | old pin (v4) | new pin (v5) | why |
|---|---|---|---|
| `kit_guard.py` | `970e4e50d7de034944c49143e5f4fdf95532b0b08778c305a38e7c90dc6cf390` | `c977ce5143896a417db18538ff2efee93cd8dd4b4261c8a812948a3720c5814c` | winfix 1: on Windows the subprocess.Popen audit event carries ONE joined command line, so the guard denied the kit's own guarded Python children; split it with CommandLineToArgvW (code of the tested Windows run copies). No effect on POSIX. |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | `81418ebb2e200731e3b219f6065bd3edd002ac64a5518c08c6d8a7fe7933b282` | `bd5451152c00086d7deec2242f752879eda7853554746dd6700be4acc8a40820` | winfix 2: `resource` is Unix-only (RSS report only, written to the timing sidecar / log, never compared): import guarded, Windows peak working set via psapi, else -1. Replaces the env shim of the kit4bv4 Windows run. Certified computation unchanged. |

## Added files

- `CHANGES_v5.md` (this file)

## Windows path audit (fix 4), kit_step3

Searched every `.py` of the kit for `os.path.relpath`, `os.sep`, `.split('/')`, `startswith(... + '/')` and backslash literals.  Nothing to change in kit 3:
- `check.py:64/67/75` relpath: in FAIL messages only.  `check.py` already reads path strings/keys by file name with `[/\\]` (rule C) and normalises absolute paths with `replace('\\', '/')`.
- `run_all.py:39` cwd assertion uses `os.sep` (correct on both); `:72` manifest paths split on '/' then `os.path.join` (correct); `:78` `os.sep + 'kit_out'` (correct).
- `kit_guard.py`: paths compared after `realpath`, with `os.sep` (correct); the command-line issue is fix 1.
- `flux/code/inputs_manifest.py:23` `rel = os.path.relpath(run, FC.N36)` is only one of three alternatives (`run in rep`, `rel in rep`, `'runs/' + basename in rep`); the last matches on Windows (the run passed).  Proof-tree script, left unchanged.
- `r3/compare_r3.py:46` `p.split('/')`: not run by the kit.
- `step3nl_v4c/v4cert/c2_t3d.py` imports `resource` but is NOT executed by the kit (a pinned input hashed by inputs_manifest.py); left byte-identical.

