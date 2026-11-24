#!/bin/sh
# second batch: sharp constants with a later tail anchor.  Waits for run_e0.sh to finish (one process at a time).
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/e0_audit
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=1 THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20 LY1SQ=5e-8
cd $D/runs || exit 1
until grep -q "^done" status.txt; do sleep 5; done
echo $$ > run_e0b.pid
mkdir -p c_sharp_rho75; cd c_sharp_rho75
nice -n 19 python3 -B $D/e0_closure_rho.py sharp_rho75 0.0104 40 0.0075 1.5e-4 > closure_1.5e-4.out 2>&1; echo "closure sharp_rho75 rc=$?" >> ../status.txt
cd ..; echo done2 >> status.txt
