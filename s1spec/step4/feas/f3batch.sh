#!/bin/bash
cd /home/claude/s1spec/step4/feas
export OMP_NUM_THREADS=1
python3 f3run.py 100 400 1600 30 10 > logs/f3_main_tip.log 2>&1
python3 f3run.py --frame fit 100 400 1600 > logs/f3_main_fit.log 2>&1
python3 f3run.py --notau 100 400 1600 30 > logs/f3_notau.log 2>&1
python3 f3run.py --sigma 28 --tag _S28 100 1600 > logs/f3_S28.log 2>&1
python3 f3run.py --sigma 45 --tag _S45 100 1600 > logs/f3_S45.log 2>&1
python3 f3run.py --eps2 1e-5 100 400 1600 30 > logs/f3_h2.log 2>&1
