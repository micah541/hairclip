# kit_step4b_v4 — Windows run 1 (Micah's i7-10700, Nov 24 2026)
Run dir: C:\Users\micahw\kitruns\kit4bv4_win1 (Python 3.13.7, flint 0.9.0, numpy 2.5.4, Windows 11 26100, --jobs 4, below-normal priority).
Started 00:03:30 Windows time; stages 1-39 ran unattended (~5 h); stage 40 (v38) retried three times after run-copy fixes.

**Result: check.py TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]** (identical counts to the WSL run kit4bv4_run_wsl1).
- WARN = v38 outputs whose text differs: the 5 flux .tex files (CRLF only; LF-normalized sha256 equal the reference, checked)
  and iface_compare_v38.json (also differs in the WSL run).
- guard: no reads outside kit/Python/OS dirs, no write/chdir/child outside the kit.
- INFO "outputs vs proof tree 0/0/0/498 different": expected (platform floats; WSL run 30/0/0/468); every claim is recomputed on this run's own outputs.

## Run-copy fixes (all Windows-only; kit v5 must include them)
1. kit_guard.py: Popen audit event gets one joined command line on Windows -> split with CommandLineToArgvW (backup kit_guard_pre_winfix.py; WINFIX_NOTE.md).
2. `resource` module missing on Windows: env-only shim kitenv/Lib/site-packages/resource.py (peak working set); chart4_shell.py left byte-identical because the input manifest pins it.
3. step4/p2/v38/flux_v38.py: rerun texts compared byte for byte, but Python writes CRLF on Windows -> compare with CRLF normalized, and gen_texts.json output hashes compared modulo CRLF (diff: flux_v38_winfix.diff; backup flux_v38_pre_winfix.py).
4. step4/p2/v38/verify_v4c_runs_with_exceptions.py: affected_records matched against os.path.relpath (backslashes on Windows) -> joined with '/' (diff: verify_v4c_winfix.diff; backup *_pre_winfix.py). (First attempt of this patch had a syntax error of mine; corrected before try 3.)

Logs: run_all_try1/2.log, check_try1/2.log (1 FAIL each: v38 did not finish), run_all.log, check.log (final).
