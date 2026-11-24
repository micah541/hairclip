#!/bin/sh
# fin13v4b: two complete passes (cert chain + production) for byte-identity.  (Regression V4B=0 against fin13v4 was run before: runs/reg_fin13v4.)
B=/home/claude/s1spec/step3nl_v4c
echo $$ > $B/runs_v4b.pid; : > $B/runs/all_status.txt
rm -rf $B/v4bcert/out $B/v4bcert/out_A
sh $B/v4bcert/run_cert_v4b.sh; echo "cert pass A done" >> $B/runs/all_status.txt
sh $B/verifier/run_v4b.sh $B/runs/P_A 1 fin13v4b; echo "production pass A done" >> $B/runs/all_status.txt
mv $B/v4bcert/out $B/v4bcert/out_A
sh $B/v4bcert/run_cert_v4b.sh; echo "cert pass B done" >> $B/runs/all_status.txt
sh $B/verifier/run_v4b.sh $B/runs/P 1 fin13v4b; echo "production pass B done" >> $B/runs/all_status.txt
echo ALLDONE >> $B/runs/all_status.txt
