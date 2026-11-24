Refresh kit B after code changes in step4/p2, step4/t7, ...:  python3 build.py  (rebuilds /home/claude/s1spec/kit_step4b from the proof tree:
copies the .py sources, re-applies the path patches, regenerates inputs/MANIFEST.sha256, copies run_all.py / check.py / README.md from this folder).
Per-stage commands are in run_all.py (functions st_*; one place); thresholds and "ours" numbers in check.py.  Then re-test with  python run_all.py --quick  and zip with files at the root.
One command: sh refresh.sh [zip]  (build.py; "zip" also writes kit_step4b.zip). Spec of the current reference JSONs: chain.json. Quick test afterwards: python run_all.py --quick
