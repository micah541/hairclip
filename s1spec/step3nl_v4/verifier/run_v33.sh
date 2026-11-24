#!/bin/sh
# v3.3 verifier rerun (tag fin12f4): FLUX0=1.64e-20 FLUX2=8.0e-20 (x4). Same environment as run_final.sh. Sequential, nice 19.
cd /home/claude/s1spec/step3nl_v33/verifier || exit 1
export THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20
nice -n 19 python3 tail_cert.py > fin12f4_tail.out 2>&1; echo "tail rc=$?" > fin12f4_status.txt
TAG=fin12f4 LY1SQ=5e-8 nice -n 19 python3 closure.py 1.5e-4 > fin12f4_1.5e-4.out 2>&1; echo "closure 1.5e-4 rc=$?" >> fin12f4_status.txt
# optional remark point, run after the first succeeded
TAG=fin12f4 LY1SQ=5e-8 nice -n 19 python3 closure.py 1.25e-4 > fin12f4_1.25e-4.out 2>&1; echo "closure 1.25e-4 rc=$?" >> fin12f4_status.txt
