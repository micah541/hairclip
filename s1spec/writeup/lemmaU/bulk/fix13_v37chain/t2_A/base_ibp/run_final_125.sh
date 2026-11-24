#!/bin/sh
# v3.1 alternate point alpha1 = 1.25e-4 (remark of lemmaS.tex), inputs rounded outward from the fin11 tables; ~7 min
export PYTHONHASHSEED=0
RT=../../step3nl/verifier/rate_table_fin11_a0.000125.json; ABAR=2.15369e-7; BETA=3.20505e-9
python3 lemmaS_F2.py --J 5.6765e-12 7.9745e-12 1.7620e-12 --handoff 2e-7 1e-4 0.2 --sqLam 5.2136e-8 --alpha1 1.25e-4 --core 0.11 0.52 --ramp 0.05 0.9 --mhr 0 --slopes 8.3648e-19 1.1752e-18 2.5964e-19 --rho 2.8609e-3 --M1s 0.52 --J4 5.8733e-8 --sJ4 8.6548e-15 --rampl4 --rs33 ibp2 --l4 --abar $ABAR --beta $BETA --tag F10 > F10_run.out 2>&1
RT=$RT ABAR=$ABAR BETA=$BETA python3 lemmaS_F4B.py F10_M0.52 > F4B_F10.out 2>&1
AL=1.25e-4 ABAR=$ABAR BETA=$BETA python3 lemmaS_Bbox.py 1e-5 0.15 > B10_a1.25e-4_e1e-5.out 2>&1
