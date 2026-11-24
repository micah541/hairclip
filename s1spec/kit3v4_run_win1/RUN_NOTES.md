# kit_step3_v4 (zip 0f899eb2...) on Micah's Windows 11 machine, Nov 24 2026 (run by main session from WSL via cmd.exe)
Python 3.13.7 (Windows), venv per README, python-flint 0.9.0, numpy 2.5.3, sympy 1.14.0, mpmath 1.3.0. Run dir C:\Users\micahw\kitruns\kit3v4_win2.
Two kit bugs on Windows, fixed in the run copy only (WINFIX_NOTE.md): (1) kit_guard.py: Popen audit gives joined command line -> split with CommandLineToArgvW; (2) chart4_shell.py imports Unix-only `resource` (RSS report) -> guarded. Run 1 without fix (kitruns\kit3v4_win1) failed at the guard.
Result (check2 after rerunning parts flux, proofs): 72 PASS, 6 FAIL, 3 WARN.
- 6 FAIL "flux texts ... bytes differ": all identical to the reference modulo CRLF line endings (checked with tr -d '\r'; Windows text-mode writes).
- 3 WARN: last-bit float differences (<= 5.7e-16 rel) in diagnostic fields (tau grid, gamma_ref); certified outputs equal.
- smoke flags kit_guard.py sha differs (the fix).
- Certified results equal to reference: both closures (rho 2.679225e-3 at 1.5e-4; 2.850927e-3 at 1.25e-4; max box ratio 0.78913), mod_sups (rule C), flux2_cert PASS, G1, U-1 sym, GS_boot, lemmaSA/lemmaU chain, lemmaO, ibp.
- Guard: nothing written outside the kit.
Kit v5 should include both fixes and a CRLF-insensitive text comparison in check.py.
