#!/bin/sh
# fin13v4 certification chain (stages S1, C1, U3, C3, C2, C4, C5), sequential, nice 19, fresh outputs in v4cert/out (checkpoints removed).
cd /home/claude/s1spec/step3nl_v4c/v4cert || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
rm -rf out/c2_ckpt out/c2_ckpt_test
for s in s1_identities c1_cut_sups u3_strip_form c3_modx c2_t3d c4_hardy c5_constants; do
  /usr/bin/time -v nice -n 19 python3 -u $s.py > out/$s.log 2> out/$s.time; rc=$?
  echo "$s rc=$rc maxrss_kB=$(grep 'Maximum resident' out/$s.time | awk '{print $NF}')"
  [ $rc -eq 0 ] || exit $rc
done
