#!/bin/bash
# Phase-2a repair of Lemma S(A): one driver for the whole repaired chain, both (H2) settings (one process at a time, low priority; ~20 s total).
# Stages: weights_rep (R3, R4, R5 adds) -> ramp_rep (R1, R2) -> outer_rep (R5, hypotheses H-O1/H-O2) -> assembly_rep (both settings + export R4).
set -e; cd /home/claude/s1spec/writeup/lemmaU/phase2a/A_repair
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
nice -n 19 python3 -B rigor/weights_rep.py  > out/weights_rep.log 2>&1;  echo "weights  done $(date)"
nice -n 19 python3 -B rigor/ramp_rep.py --N 1200 > out/ramp_rep.log 2>&1; echo "ramp     done $(date)"
nice -n 19 python3 -B rigor/outer_rep.py    > out/outer_rep.log 2>&1;    echo "outer    done $(date)"
nice -n 19 python3 -B rigor/assembly_rep.py > out/assembly_rep.log 2>&1; echo "assembly done $(date)"

nice -n 19 python3 -B rigor/export_for_U.py > out/export_for_U.log 2>&1; echo "export   done $(date)"
sha256sum out/*.json export_for_U.json > out/SHA256SUMS
