#!/bin/sh
# Reviewer rerun of fin13v4 production from a COPY of the verifier (v4_review/copy). One process, nice 19.
V=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/v4_review/copy/step3nl_v4/verifier
OUT=/home/claude/s1spec/writeup/lemmaU/bulk/fix5/v4_review/rerun/P_R; T=fin13v4
mkdir -p "$OUT"; cd "$OUT" || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=1 THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20
nice -n 19 python3 -B $V/tail_cert.py > ${T}_tail.out 2>&1; echo "tail rc=$?" > ${T}_status.txt
TAG=$T LY1SQ=5e-8 nice -n 19 python3 -B $V/closure.py 1.5e-4 > ${T}_1.5e-4.out 2>&1; echo "closure 1.5e-4 rc=$?" >> ${T}_status.txt
TAG=$T LY1SQ=5e-8 nice -n 19 python3 -B $V/closure.py 1.25e-4 > ${T}_1.25e-4.out 2>&1; echo "closure 1.25e-4 rc=$?" >> ${T}_status.txt
mkdir -p modsups_$T; nice -n 19 python3 -B $V/mod_sups.py $T modsups_$T > modsups_$T/mod_sups.out 2>&1; echo "mod_sups rc=$?" >> ${T}_status.txt
