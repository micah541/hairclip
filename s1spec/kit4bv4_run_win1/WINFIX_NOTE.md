Run copy fixes (Nov 24, main session), same as kit3v4 Windows run: kit_guard.py (Windows Popen audit joined cmdline -> CommandLineToArgvW split) and
writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py (Unix-only `resource` guarded). Originals: *_pre_winfix.py.
Update: chart4_shell.py patch REVERTED (kit4b input manifest pins its hash; inputs stage failed). Only kit_guard.py changed.
Instead: shim kitenv/Lib/site-packages/resource.py (environment, not a kit file) provides getrusage().ru_maxrss on Windows.
