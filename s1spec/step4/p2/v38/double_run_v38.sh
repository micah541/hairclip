#!/bin/sh
# two clean complete passes of final_tail_v38.py (nice 19, one process); copies of the manifest after each
cd /home/claude/s1spec/step4/p2/v38 || exit 1
O=/home/claude/s1spec/step4/p2/out_v38
for P in A B; do
  rm -f $O/FINAL_manifest_v38.partial.json
  nice -n 19 python3 -I -B final_tail_v38.py > $O/final_tail_v38_pass$P.log 2>&1 || { echo "pass $P FAILED"; exit 1; }
  cp $O/FINAL_manifest_v38.json $O/FINAL_manifest_v38_pass$P.json
done
echo DOUBLE_RUN_DONE
