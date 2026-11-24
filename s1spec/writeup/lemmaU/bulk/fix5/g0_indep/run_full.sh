#!/bin/sh
# G0_INDEP: full march + tail + unconditional pass with ALL reviewer terms switched on (G0FIX=bcfmpqsx, RAMPMULT=4), same environment as run_v33.sh (fin12f4).
# Runs in the COPY vp/ (outputs land in vp/); the proof tree is only read.
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix5/g0_indep/vp || exit 1
export THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01 FLUX0=1.64e-20 FLUX2=8.0e-20
export G0FIX=bcfmpqsx RAMPMULT=4
nice -n 19 python3 tail_cert.py > ../g0all_tail.out 2>&1; echo "tail rc=$?" > ../g0all_status.txt
TAG=g0all LY1SQ=5e-8 nice -n 19 python3 closure.py 1.5e-4 > ../g0all_1.5e-4.out 2>&1; echo "closure 1.5e-4 rc=$?" >> ../g0all_status.txt
