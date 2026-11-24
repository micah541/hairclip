#!/bin/bash
# T6 rerun on the rebuilt base (polished_v2, slabs 0-180); restartable (run.py skips finished slabs)
cd /home/claude/s1spec/step4/t6
nice -n 10 python3 -u run.py --slabs 0-180 --out out_v2 --jobs 2 --poldir ../t1/polished_v2 >> run_v2.log 2>&1 || exit 1
python3 -u merge.py --in out_v2 --first 0 --last 180 >> run_v2.log 2>&1 || exit 1
date > T6V2_DONE
