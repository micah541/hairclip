#!/bin/sh
# E0 audit runs: sequential, one process, nice 19, no bytecode.  Outputs under e0_audit/ only.
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/e0_audit
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=1 THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20 LY1SQ=5e-8
echo $$ > $D/runs/run_e0.pid
cd $D/runs || exit 1
nice -n 19 python3 -B $D/e0_probe.py 5e-4 > probe_5e-4.out 2>&1; echo "probe 5e-4 rc=$?" > status.txt
nice -n 19 python3 -B $D/e0_probe.py 1e-4 > probe_1e-4.out 2>&1; echo "probe 1e-4 rc=$?" >> status.txt
for spec in "conservative 0.01607 25" "sharp 0.0104 40" "nosym 0.01574642 38.4"; do
  set -- $spec
  mkdir -p c_$1; cd c_$1
  nice -n 19 python3 -B $D/e0_closure.py $1 $2 $3 1.5e-4 > closure_1.5e-4.out 2>&1; echo "closure $1 rc=$?" >> ../status.txt
  cd ..
done
echo done >> status.txt
