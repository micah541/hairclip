#!/bin/sh
cd /tmp/claude-1000/-home-claude/1625b396-7223-44a0-a7db-c079a426b542/scratchpad/probe/verifier || exit 1
export THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01 FLUX0=1.64e-20 FLUX2=8.0e-20
PROBE=1 PBOX=1e-3 PFAC=1 TAG=probe_b1e-3_f1 LY1SQ=5e-8 /usr/bin/time -v nice -n 19 python3 closure.py 1.5e-4 > probe_b1e-3_f1.out 2> probe_b1e-3_f1.err
PROBE=1 PBOX=1e-3 PFAC=3 TAG=probe_b1e-3_f3 LY1SQ=5e-8 /usr/bin/time -v nice -n 19 python3 closure.py 1.5e-4 > probe_b1e-3_f3.out 2> probe_b1e-3_f3.err
echo done > probe_done.txt
