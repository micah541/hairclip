#!/bin/bash
# v3.3 rerun of the base_ibp reference chain with the fin12f4 verifier inputs (slopes x4). Sequential, nice 19. Writes only inside this copy.
set -e; cd /home/claude/s1spec/writeup/lemmaU/bulk/t2_A_v33/base_ibp
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
nice -n 19 python3 run_F10p_ibp.py > F2_v33.log 2>&1;              echo "F10p done $(date)"
nice -n 19 python3 rigor/capture_F2.py > rigor_capture_F2_v33.log 2>&1; echo "capture done $(date)"
nice -n 19 python3 rigor/weights_exact.py > rigor_weights_exact_v33.log 2>&1; echo "weights done $(date)"
nice -n 19 python3 rigor/assembly_arb.py > rigor_assembly_arb_v33.log 2>&1;  echo "assembly done $(date)"
