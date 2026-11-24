#!/bin/sh
# E2 audit: nonlinear coefficient recertification, sequential, one process, nice 19.
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/e2_audit
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
cd $D || exit 1
R="R1:0:8:1e-4:1e-4 R1w:0:8:1e-3:1e-3 R2a:8:12:1e-3:1e-3 R2b:12:16.25:1e-3:1e-3 R3:16.25:34.25:0.02:0.02"
nice -n 19 python3 -B $D/e2_nlin_regions.py 2.5e-4 1 nl_regions_a2.5e-4.json $R > nl_regions_a2.5e-4.log 2>&1; echo "a2.5e-4 rc=$?" > nl_status.txt
nice -n 19 python3 -B $D/e2_nlin_regions.py 1e-3 1 nl_regions_a1e-3.json $R > nl_regions_a1e-3.log 2>&1; echo "a1e-3 rc=$?" >> nl_status.txt
echo done >> nl_status.txt
