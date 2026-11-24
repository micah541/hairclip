#!/bin/sh
# G0_INDEP sensitivity: certified B&B rate bound at (1) the fin12f4 tail anchor cell (floors off), (2) the first march cell, per fix flag.
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix5/g0_indep || exit 1
export THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01 FLUX0=1.64e-20 FLUX2=8.0e-20
: > cells.jsonl
for f in "" x b p m q f s c xbpmqfsc; do
  G0FIX=$f nice -n 19 python3 cell.py 3.6849e-05 2.8065e-15 2.834565e-15 0 >> cells.jsonl 2>&1
  G0FIX=$f nice -n 19 python3 cell.py 1.5e-04 2.4975e-15 2.5150e-15 1 >> cells.jsonl 2>&1
done
G0FIX="" RAMPMULT=4 nice -n 19 python3 cell.py 3.6849e-05 2.8065e-15 2.834565e-15 0 >> cells.jsonl 2>&1
G0FIX=xbpmqfsc RAMPMULT=4 nice -n 19 python3 cell.py 3.6849e-05 2.8065e-15 2.834565e-15 0 >> cells.jsonl 2>&1
G0FIX=xbpmqfsc RAMPMULT=4 nice -n 19 python3 cell.py 1.5e-04 2.4975e-15 2.5150e-15 1 >> cells.jsonl 2>&1
echo done >> cells.jsonl
