#!/bin/bash
# v37 flux reassembly: regression (F2_VTAG=fin13v4b on fix6_v36/runs/prod2 -> must reproduce fix8_flux2_impl/out), then two passes on fin13v4c + v37 prod2.
set -e
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/flux
B=/home/claude/s1spec/writeup/lemmaU/bulk
F2_VTAG=fin13v4b F2_V36RUN=$B/fix6_v36/runs/prod2 bash $D/code/run_all.sh reg_v36 > $D/snap_reg_v36.log 2>&1
F2_VTAG=fin13v4c F2_V36RUN=$B/fix13_v37chain/runs/prod2 bash $D/code/run_all.sh v37_pass1 > $D/snap_v37_pass1.log 2>&1
F2_VTAG=fin13v4c F2_V36RUN=$B/fix13_v37chain/runs/prod2 bash $D/code/run_all.sh v37_pass2 > $D/snap_v37_pass2.log 2>&1
echo ALLDONE
