#!/bin/sh
# fin13v4b production (verifier v4 + E0 fix F-B, E0_LEMMA_v4b.tex). Same environment as run_v4.sh (fin13v4) + V4B=1.
# Usage: run_v4b.sh OUTDIR [V4BFLAG] [TAG]  -- V4BFLAG=0 with TAG=fin13v4 is the regression against fin13v4.  Sequential, nice 19, no bytecode.
V=/home/claude/s1spec/step3nl_v4c/verifier
OUT=${1:?outdir}; V4BF=${2:-1}; T=${3:-fin13v4b}
mkdir -p "$OUT"; cd "$OUT" || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=1 V4B=$V4BF THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20
TM="/usr/bin/time -a -o ${T}_status.txt -f"
: > ${T}_status.txt
$TM "tail maxrss_kb=%M wall_s=%e" nice -n 19 python3 -B $V/tail_cert.py > ${T}_tail.out 2>&1; echo "tail rc=$?" >> ${T}_status.txt
TAG=$T LY1SQ=5e-8 $TM "c1.5 maxrss_kb=%M wall_s=%e" nice -n 19 python3 -B $V/closure.py 1.5e-4 > ${T}_1.5e-4.out 2>&1; echo "closure 1.5e-4 rc=$?" >> ${T}_status.txt
TAG=$T LY1SQ=5e-8 $TM "c1.25 maxrss_kb=%M wall_s=%e" nice -n 19 python3 -B $V/closure.py 1.25e-4 > ${T}_1.25e-4.out 2>&1; echo "closure 1.25e-4 rc=$?" >> ${T}_status.txt
mkdir -p modsups_$T; $TM "modsups maxrss_kb=%M wall_s=%e" nice -n 19 python3 -B $V/mod_sups.py $T modsups_$T > modsups_$T/mod_sups.out 2>&1; echo "mod_sups rc=$?" >> ${T}_status.txt
