#!/bin/bash
N=/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36; cd $N
V36_VER=fin12f4 V36_R3=0 V36_R2=0 bash run_chain_v36.sh reg_v35_final > runs_reg_v35_final.log 2>&1 || echo FAIL1
V36_VER=fin12f4 V36_R3=1 V36_R2=0 bash run_chain_v36.sh reg_r3_final > runs_reg_r3_final.log 2>&1 || echo FAIL2
echo DONE
