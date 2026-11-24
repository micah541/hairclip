#!/bin/bash
# usage: test_cfg.sh TAG P JT h0 cT  : short float run to k=37 then polished chain over slabs 10..15
cd /home/claude/s1spec/step4/t1
export OMP_NUM_THREADS=1
python3 run2.py $1 37 $2 $3 $4 $5 80 0.12 0.02 1.5 8 0.15
python3 polish_chain.py $1 10 6 4 > chain_$1.out 2>&1
