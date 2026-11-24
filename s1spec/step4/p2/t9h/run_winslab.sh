#!/bin/bash
# windows and shells (win.py) on the 8 t_c balls of radius 1/8 covering slab 180; usage: run_winslab.sh M0 M1
cd /home/claude/s1spec/step4/p2/t9h
for m in $(seq $1 $2); do
  line=$(sed -n "$((m + 1))p" ../out/h1slab/win_cmds.txt)
  [ -f ../out/h1slab/win_${m}.json ] || env $line nice -n 10 python3 -I win.py > ../out/h1slab/win_${m}.log 2>&1
done
