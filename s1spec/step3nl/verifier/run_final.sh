#!/bin/sh
# v3.1 final verifier runs (tag fin11): tail_cert, then both alpha points, sequentially.  ~13 min.
export THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
python3 tail_cert.py > /tmp/fin11_tail.out 2>&1
for a in 1.5e-4 1.25e-4; do TAG=fin11 LY1SQ=5e-8 python3 closure.py $a > /tmp/fin11_$a.out 2>&1; done
