#!/bin/sh
# E2 audit probe: one process, nice 19, no bytecode.  Outputs under e2_audit/runs only.
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/e2_audit
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=1 THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20 LY1SQ=5e-8
mkdir -p $D/runs; cd $D/runs || exit 1
nice -n 19 python3 -B $D/e2_probe.py 5e-4 > probe_5e-4.out 2>&1; echo "probe rc=$?" > status.txt
