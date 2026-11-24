#!/usr/bin/env python3
"""Build kit_step4b v5 = kit_step4b v4 (tree == zip ab47b541...) + the Windows fixes of the kit4bv4 Windows run (s1spec/kit4bv4_run_win1/RUN_NOTES.md).
Writes /home/claude/s1spec/kit_step4b_v5/ (must not exist), kit_step4b_v5.zip (deterministic, flat like v4), CHANGES_v5.md (kit root and here), BUILD_REPORT_v5.md here.
Reads kit_step4b_v4/ and kit_step4b_v4.zip only.   usage: nice -n 19 python3 -B build_v5.py"""
import os, sys, json, shutil
sys.path.insert(0, '/home/claude/s1spec/kit_v5_build'); import v5common as C   # noqa: E402
S = C.S; HERE = os.path.dirname(os.path.abspath(__file__))
B = C.Build(S + '/kit_step4b_v4', S + '/kit_step4b_v4.zip', 'ab47b5416e8a6d58d5267e2a93c154bf979809005de0abb4e992dd989df7f0fa', '',
            S + '/kit_step4b_v5', S + '/kit_step4b_v5.zip', '', (2026, 11, 24, 0, 0, 0))
B.check_v4_tree_equals_zip(); B.copy()

# 1 + 2: common
C.patch_guard(B); C.patch_chart4(B)

# 3: flux_v38.py text / gen_texts.json output-hash compares modulo CRLF (the tested run-copy patch, flux_v38_winfix.diff)
FV = 'step4/p2/v38/flux_v38.py'
B.patch(FV, [
    ("CMP = {}\n",
     "_lf = lambda q: hashlib.sha256(open(q, 'rb').read().replace(b'\\r\\n', b'\\n')).hexdigest()   # [kit patch v5, winfix 3] CRLF-insensitive (texts written in text mode on Windows)\n"
     "CMP = {}\n", 'winfix3-lf'),
    ("    CMP[f] = 'equal apart from args' if a == b else 'DIFFERS'\n",
     "    if f == 'gen_texts.json' and isinstance(b.get('outputs'), dict):   # [kit patch v5, winfix 3] output hashes compared modulo CRLF\n"
     "        b = dict(b, outputs={k: (a.get('outputs', {}).get(k) if os.path.exists(f'{OT}/{k}') and _lf(f'{OT}/{k}') == a.get('outputs', {}).get(k) else v) for k, v in b['outputs'].items()})\n"
     "    CMP[f] = 'equal apart from args' if a == b else 'DIFFERS'\n", 'winfix3-gen_texts'),
    ("CMP_T = {f: ('byte-identical' if os.path.exists(f'{OT}/{f}') and sha(f'{FX}/texts/{f}') == sha(f'{OT}/{f}') else 'DIFFERS') for f in TEXTS}\n",
     "CMP_T = {f: ('byte-identical' if os.path.exists(f'{OT}/{f}') and _lf(f'{FX}/texts/{f}') == _lf(f'{OT}/{f}') else 'DIFFERS') for f in TEXTS}   # [kit patch v5, winfix 3] 'byte-identical' = modulo CRLF (label kept: flux_v38.json is compared with the proof tree's)\n",
     'winfix3-texts'),
], 'winfix 3: the 6 rerun flux texts and the output hashes recorded in gen_texts.json are compared modulo CRLF (Python text mode writes CRLF on Windows); '
   'the tested run-copy patch (kit4bv4_run_win1/flux_v38_winfix.diff). Compare side: gen_texts.py (proof-tree script) unchanged. Identical result on POSIX.')

# 4: verify_v4c_runs_with_exceptions.py: relpath -> '/' (the tested line 43 patch, plus the two record keys of its output)
VV = 'step4/p2/v38/verify_v4c_runs_with_exceptions.py'
B.patch(VV, [
    ("and 'runs/' + os.path.relpath(f, R) in E['affected_records']}\n",
     "and 'runs/' + os.path.relpath(f, R).replace(os.sep, '/') in E['affected_records']}   # [kit patch v5, winfix 4] separator: relpath joined with '/'\n", 'winfix4-affected'),
    ("        rows[os.path.relpath(f, R)] = dict(",
     "        rows[os.path.relpath(f, R).replace(os.sep, '/')] = dict(", 'winfix4-rowkey'),
    ("pairs and all(eq.values()), eq)\n", "pairs and all(eq.values()), eq)   # [kit patch v5, winfix 4] keys of rows / eq: relpath joined with '/' (same JSON on Windows as on POSIX)\n", 'winfix4-comment'),
    ("    eq = {os.path.relpath(f, R): os.path.exists(g)", "    eq = {os.path.relpath(f, R).replace(os.sep, '/'): os.path.exists(g)", 'winfix4-eqkey'),
], 'winfix 4: the exception record lists affected records as runs/<run>/<file> with \'/\'; os.path.relpath gives backslashes on Windows, so the documented '
   'exception was not recognised (stage 28 failed). relpath joined with \'/\' (tested line) and the two output dict keys made separator-independent. Identical on POSIX.')

# 3 (checker side): verify_v33.compare: a non-JSON output equal modulo CRLF counts as EXACT (the 5 flux .tex WARN of the Windows run)
B.patch('verify_v33.py', [
    ("            if not rel.endswith('.json'): out.append((rel, 'DIFFERENT', 'not JSON')); continue\n",
     "            if not rel.endswith('.json') and _crlf_eq(p, rp): out.append((rel, 'EXACT', 'identical modulo CRLF line endings')); continue   # [kit patch v5, winfix 3]\n"
     "            if not rel.endswith('.json'): out.append((rel, 'DIFFERENT', 'not JSON')); continue\n", 'winfix3-exact'),
    ("def compare(REFD=REFD, skip_manifest='FINAL_manifest_v33.json'):\n",
     "def _crlf_eq(p, q):\n"
     "    \"\"\"[kit patch v5, winfix 3] text files equal modulo CRLF line endings (Python text mode writes CRLF on Windows)\"\"\"\n"
     "    with open(p, 'rb') as f1, open(q, 'rb') as f2: return f1.read().replace(b'\\r\\n', b'\\n') == f2.read().replace(b'\\r\\n', b'\\n')\n\n\n"
     "def compare(REFD=REFD, skip_manifest='FINAL_manifest_v33.json'):\n", 'winfix3-helper'),
    ("  EXACT      same bytes\n", "  EXACT      same bytes  ([v5] or, for a non-JSON text file, same bytes modulo CRLF line endings)\n", 'doc'),
], 'winfix 3 (checker): a non-JSON output (the 5 flux .tex) equal to the reference modulo CRLF is EXACT (was DIFFERENT -> WARN in the Windows run).')
B.patch('check.py', [
    ('"""Kit B v4 check (', '"""Kit B v5 check ([v5] Windows fixes, CHANGES_v5.md; text outputs compared modulo CRLF).  Kit B v4 check (', 'version'),
    ("    add('INFO', 'v3.3 outputs vs proof tree: byte-identical / equal",
     "    add('INFO', 'v3.3 outputs vs proof tree: byte-identical (text: modulo CRLF) / equal", 'label'),
    ("    add('INFO', 'v38 outputs vs proof tree: byte-identical / equal",
     "    add('INFO', 'v38 outputs vs proof tree: byte-identical (text: modulo CRLF) / equal", 'label'),
    ("'KIT 4B check: certified", "'KIT 4B check (kit v5): certified", 'version'),
], 'version strings v4 -> v5; INFO labels say that text outputs count as byte-identical modulo CRLF (verify_v33.py change).')
B.patch('run_all.py', [('"""Reproduction kit B v3: ', '"""Reproduction kit B v5 (= v4 + Windows fixes, CHANGES_v5.md; v4 = v3 + stage v38 + write guard).  Kit B v3: ', 'version')],
        'version string (docstring only).')

# 5: README
B.patch('README.md', [
    ('# Reproduction kit 4B v4: step 4 (v3.2 chain), the 21 v3.3 stages and the 29-stage v38 re-pin chain\n',
     '# Reproduction kit 4B v5: step 4 (v3.2 chain), the 21 v3.3 stages and the 29-stage v38 re-pin chain\n\n'
     '**New in v5 (against v4): the Windows fixes** of the first Windows run of v4 (Nov 24 2026, 159 PASS / 0 FAIL / 1 WARN after run-copy patches; `CHANGES_v5.md`): '
     '(1) the write guard splits the joined Windows command line (CommandLineToArgvW); (2) `chart4_shell.py` runs without the Unix-only `resource` module '
     '(its pin in `inputs/MANIFEST.sha256` is regenerated); (3) the flux texts and their recorded hashes are compared modulo CRLF line endings '
     '(`flux_v38.py`, `verify_v33.py`); (4) `verify_v4c_runs_with_exceptions.py` joins relative paths with `/`. No certified computation changed.\n\n'
     '**Windows 11 (PowerShell).** Tested: CPython 3.13.7 from python.org, a venv, python-flint 0.9.0, numpy 2.5.4, sympy 1.14.0, mpmath 1.3.0, scipy 1.18.1 '
     '(`requirements.txt`). Unzip to a short path outside OneDrive, e.g. `C:\\kitruns\\kit4bv5`:\n'
     '```powershell\n'
     'Expand-Archive -Path .\\kit_step4b_v5.zip -DestinationPath C:\\kitruns\\kit4bv5\n'
     'cd C:\\kitruns\\kit4bv5\n'
     'py -3.13 -m venv kitenv\n'
     '.\\kitenv\\Scripts\\python.exe -m pip install --upgrade pip\n'
     '.\\kitenv\\Scripts\\python.exe -m pip install -r requirements.txt\n'
     '.\\kitenv\\Scripts\\python.exe run_all.py --jobs 4          # 5-7 h; resumable: run the same command again\n'
     '.\\kitenv\\Scripts\\python.exe check.py\n'
     '```\n'
     'Expected: `TOTAL PASS`; 1 WARN (`iface_compare_v38.json` has 2 fewer checks than the proof tree\'s, by design: C1_ke / GS_tau1 read this run\'s own step-4 data). '
     '**The INFO line "outputs vs proof tree: byte-identical / timing-fields-only / provenance-only / different / missing" is INFO, not a failure**: with other floats '
     '(another platform / libm) almost every regenerated v3.2 file differs from the proof tree\'s in the last digits (Windows v4 run: 0/0/0/498 different; WSL 30/0/0/468); '
     'every claim is recomputed on this run\'s own outputs and checked as an inequality.\n', 'v5-head'),
], 'winfix 5: README: v5 note, Windows instructions (venv, Python 3.13.7 tested, flint 0.9.0), expected WARN, the "outputs vs proof tree ... different" line is INFO.')

B.regen_manifest_pins('inputs/MANIFEST.sha256')

AUDIT = '''## Windows path audit (fix 4), kit_step4b

Searched every `.py` of the kit for `os.path.relpath`, `os.sep`, `.split('/')`, `startswith(... + '/')` and backslash literals.
- FIXED `step4/p2/v38/verify_v4c_runs_with_exceptions.py:43` (the run-copy patch: membership test against the exception record), `:44` and `:51` (dict keys of its output: now '/' on every OS, so `v4c_runs_pins.json` is the same JSON on Windows).
- `step4/p2/v38/verify_v4c_runs_with_exceptions.py:50` `os.path.join(b, os.path.relpath(...))`: a file path, correct on both.
- `step4/p2/v33/final_tail_v33.py:49-50`, `step4/p2/v38/final_tail_v38.py:46-47`: relpath already followed by `replace('\\\\', '/')`; `startswith(S4 + '/')`, `startswith(OV + '/')`: S4/OV are built from `_S1` (forward slashes) and the stage paths the drivers compare are built the same way (the Windows v4 run passed every stage).
- `verify_v33.py:65/67`, `verify_outputs.py:29/33`, `run_all.py:191/380/382`: relpath with `replace(os.sep, '/')` / manifest paths split on '/' then `os.path.join` / `os.sep`: correct on both.
- `step4/p2/v38/iface_compare_v38.py:67`: relpath in a check NAME only (text); not compared as a key (iface_compare_v38.json is the known benign WARN in both runs).
- step-4 v3.2 provenance (`step4/t7/se7.py:68-116`, `merge.py:125/139`, `nlterms7.py:596`, `provenance.py:22`, `t6c/runc.py:34`, `t7/run_all.py:63`): relpath written into args/provenance fields of outputs; check.py / verify_outputs.py drop args and normalise separators; proof-tree reproductions, left unchanged (the Windows v4 run passed with them).
- `t1/hotab.py:11`, `t1/cmp_hls.py:12` `f.split('/')`: diagnostics, not run by the kit.
- `writeup/.../flux/code/inputs_manifest.py:23`: relpath is one of three alternatives; the third (`'runs/' + basename`) matches on Windows. Proof-tree script, unchanged.
- `step3nl_v4c/v4cert/c2_t3d.py` imports `resource` but is NOT executed by the kit (pinned input hashed by inputs_manifest.py); left byte-identical.
- `kit_guard.py`: realpath + `os.sep` comparisons, correct; command-line split is fix 1.
'''
txt = B.changes_md('CHANGES_v5: kit_step4b v4 -> v5 (Windows fixes)',
                   'Base: `kit_step4b_v4.zip` sha256 `ab47b5416e8a6d58d5267e2a93c154bf979809005de0abb4e992dd989df7f0fa` (tree verified equal to the zip). '
                   'Source of the fixes: the kit4bv4 Windows run (`s1spec/kit4bv4_run_win1/RUN_NOTES.md`, `*.diff`; run copy `C:\\Users\\micahw\\kitruns\\kit4bv4_win1`). '
                   'Every changed code line carries a `[kit patch v5 ...]` comment (README / version strings excepted). '
                   'Fix 2 is made in the script (not by an environment shim), so its pin in `inputs/MANIFEST.sha256` is regenerated (table below); no other pin of the kit '
                   'refers to a changed file (flux_v38.py, verify_v4c_runs_with_exceptions.py, verify_v33.py, check.py, run_all.py, kit_guard.py, README.md are not in the manifest). '
                   'Choice for CRLF (fix 3): compare side (flux_v38.py, verify_v33.py), so the proof-tree scripts (gen_texts.py) and their outputs are unchanged. '
                   'All other files are byte-identical to v4 (checked by the builder).', AUDIT)
B.add_file('CHANGES_v5.md', txt, 'this list')
B.verify_unchanged_rest(); B.all_py_parse(); B.zip()
shutil.copy(B.v5 + '/CHANGES_v5.md', HERE + '/CHANGES_v5.md')
with open(HERE + '/BUILD_REPORT_v5.md', 'w', encoding='utf-8', newline='\n') as f:
    f.write(f'# BUILD_REPORT_v5: kit_step4b v5\n\nBuilt by `kit_step4b_v5_build/build_v5.py` (shared code `kit_v5_build/v5common.py`) from `kit_step4b_v4/` (== `kit_step4b_v4.zip`).\n\n'
            f'**Zip `kit_step4b_v5.zip`: {B.zsize} bytes, sha256 `{B.zsha}`.** Tree `kit_step4b_v5/`: {B.nfiles} files (zip is flat, as v4; date 2026-11-24, mode 644).\n\n'
            f'Changed files: {len(B.changes)}; manifest pins regenerated: {len(B.pins)}; added: {len(B.added)}. See CHANGES_v5.md (copy next to this file).\n')
json.dump(dict(zip=B.zip5, sha256=B.zsha, size=B.zsize, nfiles=B.nfiles, changes=B.changes, pins=B.pins, added=B.added), open(HERE + '/build_v5.json', 'w'), indent=1)
print('done')
