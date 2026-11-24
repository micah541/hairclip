#!/bin/bash
# Lemma S(A), tightened (H2), core {2,3,3} by IBP. Sequential, low priority (Micah: low RAM, no hurry).
set -e; cd /home/claude/s1spec/writeup/step3_gap_v32/gonogo/lemmaS_A_ibp_v31H2
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
nice -n 19 python3 run_F10p_ibp.py > F2_ibp.log 2>&1;              echo "F10p done $(date)"
nice -n 19 python3 rigor/capture_F2.py > rigor_capture_F2.log 2>&1; echo "capture done $(date)"
nice -n 19 python3 rigor/weights_exact.py > rigor_weights_exact.log 2>&1; echo "weights done $(date)"
nice -n 19 python3 rigor/assembly_arb.py > rigor_assembly_arb.log 2>&1;  echo "assembly done $(date)"
