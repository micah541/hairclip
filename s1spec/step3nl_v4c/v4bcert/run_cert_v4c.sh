#!/bin/sh
# fin13v4c certification chain = the fin13v4b chain (run_cert_v4b.sh, repointed to step3nl_v4c) + two new stages:
#   strip_cubic_v4c   (D-3: certified strip cubic coefficients C3S_abc, read by the verifier when V4C=1)
#   e0_A2_budget_v4c  (A-1: E0 A^2 budget with the (a,b) lead mismatch and the third-order loop piece; asserts the absorption)
# Sequential, one process, nice 19, no bytecode.  Outputs in v4bcert/out.
B=/home/claude/s1spec/step3nl_v4c
sh $B/v4bcert/run_cert_v4b.sh
cd $B/v4bcert || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
for st in strip_cubic_v4c e0_A2_budget_v4c; do
  /usr/bin/time -f "$st maxrss_kb=%M wall_s=%e" -a -o out/chain_status.txt nice -n 19 python3 -B $st.py > out/$st.log 2>&1; echo "$st rc=$?" >> out/chain_status.txt
done
