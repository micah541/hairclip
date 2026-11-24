#!/bin/sh
# Task 4 chain: one process at a time, nice 19. Usage: sh run_all.sh [A_interface path]
set -e
cd /home/claude/s1spec/writeup/lemmaU/bulk/t4_cert
AI=${1:-/home/claude/s1spec/writeup/lemmaU/bulk/t2_A/out/A_interface.json}
R="env PYTHONHASHSEED=0 nice -n 19 python3 -B"
$R u1_ghat.py > out/u1_ghat.log 2>&1
$R u2_K.py --fresh > out/u2_K.log 2>&1
$R u2_K.py --fresh --nophi --tag _nophi > out/u2_K_nophi.log 2>&1
$R u2_K.py --fresh --case fallback --glo 620000 --ghi 770000 --tag _fallback > out/u2_K_fallback.log 2>&1
$R u34_sweeps.py --case primary > out/u34_primary.log 2>&1
$R u34_sweeps.py --case fallback > out/u34_fallback.log 2>&1
$R u5_enum.py > out/u5_enum.log 2>&1
$R assemble_U.py --A "$AI" > out/assemble_U.log 2>&1
$R assemble_U.py --A "$AI" --case fallback > out/assemble_U_fallback.log 2>&1
rm -f out/*.partial.json
cd out && sha256sum *.json > SHA256SUMS
