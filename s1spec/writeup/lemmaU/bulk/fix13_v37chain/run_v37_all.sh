#!/bin/bash
# fix13_v37chain (tag v37): regression (fin13v4b inputs -> must reproduce fix6_v36/runs/prod2), then two production passes on fin13v4c.
# Resumable: a finished snapshot is skipped (run_chain_v37.sh refuses existing names).  One process at a time, nice 19 (inside the chain).
N=/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain
echo $$ > $N/run.pid
for spec in "reg_v36 fin13v4b" "prod1 fin13v4c" "prod2 fin13v4c"; do set -- $spec
  if [ -e $N/runs/$1/env.txt ] && grep -q "^end " $N/runs/$1/env.txt; then echo "skip $1"; continue; fi
  rm -rf $N/runs/$1
  V36_VER=$2 V36_R3=1 V36_R2=1 bash $N/run_chain_v37.sh $1 > $N/runs_$1.log 2>&1 || { echo "FAIL $1"; exit 1; }
  echo "done $1 $(date -Is)"
done
python3 -I -B $N/cmp_runs.py $N/runs/reg_v36 /home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36/runs/prod2 --json $N/runs/cmp_reg_vs_v36prod2.json > $N/runs/cmp_reg_vs_v36prod2.txt 2>&1
python3 -I -B $N/cmp_runs.py $N/runs/prod1 $N/runs/prod2 --json $N/runs/cmp_prod1_prod2.json > $N/runs/cmp_prod1_prod2.txt 2>&1
tail -1 $N/runs/cmp_reg_vs_v36prod2.txt; tail -1 $N/runs/cmp_prod1_prod2.txt
echo ALLDONE
