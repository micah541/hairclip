#!/bin/bash
cd /home/claude/s1spec/step4/feas2; export PYTHONPYCACHEPREFIX=$PWD/pycache
python3 p1run.py k1600,k100,k10,t4.59 0.1 3 100 > logs/p1_S100_nc3.log 2>&1
python3 p1run.py t4.59,t4.5,t2.0,t0.5 0.1 4 100 > logs/p1_S100_nc4.log 2>&1
