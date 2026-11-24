#!/bin/sh
# fin13v4 production (verifier v4, EI_SPEC.md). Same environment as run_v33.sh (fin12f4) + NEW=1 (v4 E_I inequality, constants_v4.json).
# Usage: run_v4.sh OUTDIR [NEWFLAG] [TAG]  -- scripts are run from verifier/ by absolute path, all outputs go to OUTDIR (cwd). Sequential, nice 19.
V=/home/claude/s1spec/step3nl_v4/verifier
OUT=${1:?outdir}; NEWF=${2:-1}; T=${3:-fin13v4}
mkdir -p "$OUT"; cd "$OUT" || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 NEW=$NEWF THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 DELTA=0.2 MOD=ibp MODA=25 ETA=0.01
export FLUX0=1.64e-20 FLUX2=8.0e-20
nice -n 19 python3 $V/tail_cert.py > ${T}_tail.out 2>&1; echo "tail rc=$?" > ${T}_status.txt
TAG=$T LY1SQ=5e-8 nice -n 19 python3 $V/closure.py 1.5e-4 > ${T}_1.5e-4.out 2>&1; echo "closure 1.5e-4 rc=$?" >> ${T}_status.txt
TAG=$T LY1SQ=5e-8 nice -n 19 python3 $V/closure.py 1.25e-4 > ${T}_1.25e-4.out 2>&1; echo "closure 1.25e-4 rc=$?" >> ${T}_status.txt
mkdir -p modsups_$T; nice -n 19 python3 $V/mod_sups.py $T modsups_$T > modsups_$T/mod_sups.out 2>&1; echo "mod_sups rc=$?" >> ${T}_status.txt
