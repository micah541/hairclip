#!/bin/bash
# usage: bash code/run_stage.sh <script.py> [args...]   -- one process, nice 19, PYTHONHASHSEED=0, python3 -B, maxrss logged
set -e
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix8_flux2_impl
S=$1; shift
name=$(basename $S .py)
cd $D/code
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
/usr/bin/time -v -o $D/out/$name.time nice -n 19 python3 -B $S "$@" 2>&1 | tee $D/out/$name.log
grep -E "Maximum resident|Elapsed" $D/out/$name.time
