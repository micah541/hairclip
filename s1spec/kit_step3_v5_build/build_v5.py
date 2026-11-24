#!/usr/bin/env python3
"""Build kit_step3 v5 = kit_step3 v4 (tree == zip 0f899eb2...) + the Windows fixes of the kit3v4 Windows run (s1spec/kit3v4_run_win1/RUN_NOTES.md).
Writes /home/claude/s1spec/kit_step3_v5/ (must not exist), kit_step3_v5.zip (deterministic), CHANGES_v5.md (kit root and here), BUILD_REPORT_v5.md here.
Reads kit_step3_v4/ and kit_step3_v4.zip only.   usage: nice -n 19 python3 -B build_v5.py"""
import os, sys, json, shutil
sys.path.insert(0, '/home/claude/s1spec/kit_v5_build'); import v5common as C   # noqa: E402
S = C.S; HERE = os.path.dirname(os.path.abspath(__file__))
B = C.Build(S + '/kit_step3_v4', S + '/kit_step3_v4.zip', '0f899eb2de74cf51e25aafb9ac0c401343279a9a34f1dbf0dc7a2286809ce853', 'kit_step3_v4/',
            S + '/kit_step3_v5', S + '/kit_step3_v5.zip', 'kit_step3_v5/', (2026, 11, 24, 12, 0, 0))
B.check_v4_tree_equals_zip(); B.copy()

# 1 + 2: common
C.patch_guard(B); C.patch_chart4(B)

# 3: check.py text compares modulo CRLF (compare side; the proof-tree scripts and their outputs are unchanged)
CMP_TEXT = '''def cmp_text(label, newp, refp):
    """[kit patch v5, winfix 3] a generated TEXT file: bytes first; else equal modulo CRLF line endings (Python text mode writes CRLF on Windows) -> PASS"""
    if not os.path.exists(newp): add('FAIL', label, 'missing ' + os.path.relpath(newp, KIT)); return
    a, b = open(newp, 'rb').read(), open(refp, 'rb').read()
    if a == b: add('PASS', label, 'byte-identical to the reference'); return
    same = a.replace(b'\\r\\n', b'\\n') == b.replace(b'\\r\\n', b'\\n')
    add('PASS' if same else 'FAIL', label, 'identical to the reference modulo CRLF line endings' if same else 'bytes differ (also modulo CRLF line endings)')
'''
B.patch('check.py', [
    ('SUMMARY.   (kit_step3 v4)\n', 'SUMMARY.   (kit_step3 v5 = v4 + Windows fixes, CHANGES_v5.md)\n', 'version'),
    ("def ran(*logs):", CMP_TEXT + "def ran(*logs):", 'winfix3-cmp_text'),
    ("for f in sorted(os.listdir(os.path.join(REF, 'flux', 'texts'))): cmp_bytes('flux texts: ' + f,",
     "for f in sorted(os.listdir(os.path.join(REF, 'flux', 'texts'))): cmp_text('flux texts: ' + f,", 'winfix3-use'),
    ("'REPRODUCTION CHECK kit_step3 v4  (python", "'REPRODUCTION CHECK kit_step3 v5  (python", 'version'),
], 'winfix 3: the 6 generated flux texts (5 .tex + numbers_flux2.json) are compared byte for byte OR modulo CRLF (Windows text-mode writes; the 6 FAIL of the '
   'kit3v4 Windows run); compare-side fix, the proof-tree script gen_texts.py and its outputs are unchanged. Version strings v4 -> v5.')
B.patch('run_all.py', [('(kit_step3 v4).\n', '(kit_step3 v5 = v4 + Windows fixes, CHANGES_v5.md).\n', 'version')], 'version string v4 -> v5 (docstring only).')

# 5: README
B.patch('README.md', [
    ('# Reproduction kit step 3, v4 (verifier fin13v4c, (A) chain v37, flux lemma v2, interface / proof checks)\n',
     '# Reproduction kit step 3, v5 (verifier fin13v4c, (A) chain v37, flux lemma v2, interface / proof checks)\n\n'
     '**v5 = v4 + the Windows fixes** of the first Windows run (Nov 24 2026; `CHANGES_v5.md`): the write guard splits the Windows command line, '
     '`chart4_shell.py` runs without the Unix-only `resource` module, and check.py compares the generated flux texts modulo CRLF line endings. '
     'The certified computations and the reference outputs are those of v4.\n', 'v5-title'),
    ('kitenv/bin/python check.py                                   # paste the SUMMARY (and any FAIL/WARN lines) back\n```\n',
     'kitenv/bin/python check.py                                   # paste the SUMMARY (and any FAIL/WARN lines) back\n```\n\n'
     '**Windows 11 (PowerShell).** Tested with CPython 3.13.7, python-flint 0.9.0, numpy 2.5.3, sympy 1.14.0, mpmath 1.3.0 (kit v4 + these fixes). '
     'Unzip to a short path outside OneDrive (e.g. `C:\\kitruns\\kit3v5`), then:\n'
     '```powershell\n'
     'cd C:\\kitruns\\kit3v5\\kit_step3_v5\n'
     'py -3.13 -m venv kitenv\n'
     '.\\kitenv\\Scripts\\python.exe -m pip install --upgrade pip\n'
     '.\\kitenv\\Scripts\\python.exe -m pip install -r requirements.txt\n'
     '.\\kitenv\\Scripts\\python.exe run_all.py                  # about 35-40 min\n'
     '.\\kitenv\\Scripts\\python.exe check.py\n'
     '```\n'
     'Expected on Windows: every certified result PASS; a few WARN for last-bit float differences (relative <= 1e-15) in diagnostic fields '
     '(tau grid, gamma_ref) are platform noise, not failures. Flux texts written with CRLF line endings PASS as "identical modulo CRLF".\n\n', 'v5-windows'),
], 'winfix 5: README: v5 note and Windows instructions (venv, Python 3.13.7 tested, python-flint 0.9.0), expected WARNs.')

# manifest pins of the changed kit files (kit_guard.py, chart4_shell.py are pinned; check.py / run_all.py / README.md are not)
B.regen_manifest_pins('inputs/MANIFEST.sha256')

AUDIT = '''## Windows path audit (fix 4), kit_step3

Searched every `.py` of the kit for `os.path.relpath`, `os.sep`, `.split('/')`, `startswith(... + '/')` and backslash literals.  Nothing to change in kit 3:
- `check.py:64/67/75` relpath: in FAIL messages only.  `check.py` already reads path strings/keys by file name with `[/\\\\]` (rule C) and normalises absolute paths with `replace('\\\\', '/')`.
- `run_all.py:39` cwd assertion uses `os.sep` (correct on both); `:72` manifest paths split on '/' then `os.path.join` (correct); `:78` `os.sep + 'kit_out'` (correct).
- `kit_guard.py`: paths compared after `realpath`, with `os.sep` (correct); the command-line issue is fix 1.
- `flux/code/inputs_manifest.py:23` `rel = os.path.relpath(run, FC.N36)` is only one of three alternatives (`run in rep`, `rel in rep`, `'runs/' + basename in rep`); the last matches on Windows (the run passed).  Proof-tree script, left unchanged.
- `r3/compare_r3.py:46` `p.split('/')`: not run by the kit.
- `step3nl_v4c/v4cert/c2_t3d.py` imports `resource` but is NOT executed by the kit (a pinned input hashed by inputs_manifest.py); left byte-identical.
'''
txt = B.changes_md('CHANGES_v5: kit_step3 v4 -> v5 (Windows fixes)',
                   'Base: `kit_step3_v4.zip` sha256 `0f899eb2de74cf51e25aafb9ac0c401343279a9a34f1dbf0dc7a2286809ce853` (tree verified equal to the zip). '
                   'Source of the fixes: the kit3v4 Windows run (`s1spec/kit3v4_run_win1/RUN_NOTES.md`; run copy `C:\\Users\\micahw\\kitruns\\kit3v4_win2`). '
                   'Every changed line carries a `[kit patch v5 ...]` comment, except README/version strings. '
                   'Choice for CRLF (fix 3): compare side (check.py), so the proof-tree scripts (gen_texts.py) and their outputs stay as in the proof tree. '
                   'All other files are byte-identical to v4 (checked by the builder).', AUDIT)
B.add_file('CHANGES_v5.md', txt, 'this list')
B.verify_unchanged_rest(); B.all_py_parse(); B.zip()
shutil.copy(B.v5 + '/CHANGES_v5.md', HERE + '/CHANGES_v5.md')
with open(HERE + '/BUILD_REPORT_v5.md', 'w', encoding='utf-8', newline='\n') as f:
    f.write(f'# BUILD_REPORT_v5: kit_step3 v5\n\nBuilt by `kit_step3_v5_build/build_v5.py` (shared code `kit_v5_build/v5common.py`) from `kit_step3_v4/` (== `kit_step3_v4.zip`).\n\n'
            f'**Zip `kit_step3_v5.zip`: {B.zsize} bytes, sha256 `{B.zsha}`.** Tree `kit_step3_v5/`: {B.nfiles} files. Zip entries under `kit_step3_v5/`, date 2026-11-24 12:00, mode 644.\n\n'
            f'Changed files: {len(B.changes)}; manifest pins regenerated: {len(B.pins)}; added: {len(B.added)}. See CHANGES_v5.md (copy next to this file).\n')
json.dump(dict(zip=B.zip5, sha256=B.zsha, size=B.zsize, nfiles=B.nfiles, changes=B.changes, pins=B.pins, added=B.added), open(HERE + '/build_v5.json', 'w'), indent=1)
print('done')
