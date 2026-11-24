# kit_step3 v5 test run (WSL, Nov 24 2026): 80 PASS, 0 FAIL, 1 WARN, 0 SKIP  (= kit3v4_run_wsl5, line for line)
- zip `kit_step3_v5.zip` sha256 f0ea28743d340af8a440b4280e7ccc7ed71b5306635c577da0d568b619ee413e (16264692 bytes), extracted into this new empty directory; `PYTHONHASHSEED=0 nice -n 19 python3 -B run_all.py` then `python3 -B check.py`. Python 3.14.4, flint 0.9.0, numpy 2.5.3, sympy 1.14.0. run_all finished in 1696 s, failed steps none.
- The PASS/FAIL/WARN/SKIP status lines are identical to kit3v4_run_wsl5's (diff empty). WARN = ibp final_5e-4.json float noise 5.68e-16 (unchanged kit-v2 part), as in v4.
- smoke: input manifest ok (the 2 regenerated pins kit_guard.py, chart4_shell.py verify). Flux texts PASS "byte-identical" (LF on POSIX; the CRLF branch is unit-tested in kit_v5_build/tests/unit_v5.py).
- Guard: 0 reads outside the kit; no write/chdir/child outside. Outside-write check (find -newer START_MARKER): only /home/claude/CLAUDE.md and writeup/lemmaU/bulk/fix27_texts_v44/* (other sessions' work), and my v5 build/test dirs.
Logs: run_all.log, check.log.
