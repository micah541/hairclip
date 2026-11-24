#!/bin/sh
# fin13v4c: regression (V4C=0, TAG fin13v4b, on cert pass A) + two complete passes (cert chain + production) for byte-identity.
# Cert pass A is run first by hand (runs/certA.out); this driver starts after it.  Resumable: each finished stage is listed in runs/all_status.txt and skipped.
B=/home/claude/s1spec/step3nl_v4c
echo $$ > $B/runs_v4c.pid; touch $B/runs/all_status.txt
done_() { grep -qx "$1" $B/runs/all_status.txt; }
grep -q "^done$" $B/runs/certA.out || { echo "cert pass A not finished" >&2; exit 1; }
done_ "regression done" || { rm -rf $B/runs/reg_fin13v4b; sh $B/verifier/run_v4c.sh $B/runs/reg_fin13v4b 0 fin13v4b; echo "regression done" >> $B/runs/all_status.txt; }
done_ "production pass A done" || { rm -rf $B/runs/P_A; sh $B/verifier/run_v4c.sh $B/runs/P_A 1 fin13v4c; echo "production pass A done" >> $B/runs/all_status.txt; }
done_ "cert pass B done" || { [ -d $B/v4bcert/out_A ] || mv $B/v4bcert/out $B/v4bcert/out_A; rm -rf $B/v4bcert/out; sh $B/v4bcert/run_cert_v4c.sh; echo "cert pass B done" >> $B/runs/all_status.txt; }
done_ "production pass B done" || { rm -rf $B/runs/P; sh $B/verifier/run_v4c.sh $B/runs/P 1 fin13v4c; echo "production pass B done" >> $B/runs/all_status.txt; }
echo ALLDONE >> $B/runs/all_status.txt
