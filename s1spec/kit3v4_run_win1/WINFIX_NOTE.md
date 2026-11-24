Run copy fix (Nov 24, main session): kit_guard.py on Windows denied the kit's own guarded Python children, because CPython's
subprocess.Popen audit event on Windows carries the joined command line (one string), so _is_py() saw the whole string.
Patch: split it with CommandLineToArgvW before the checks (see diff vs kit_guard_pre_winfix.py). No other file changed.
First attempt (unpatched): C:\Users\micahw\kitruns\kit3v4_win1 (2 PASS, 4 FAIL, 8 SKIP; all FAILs from this guard error).

Also: writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py imported Unix-only `resource` (RSS report only); guarded (diff vs chart4_shell_pre_winfix.py). Parts flux+proofs rerun after.
