#!/bin/sh
# v3.1 final Lemma S chain at alpha1 = 1.5e-4 (inputs rounded in the safe direction from the fin11 verifier tables; ~16 min)
export PYTHONHASHSEED=0
RT=../../step3nl/verifier/rate_table_fin11_a0.00015.json; ABAR=2.19103e-7; BETA=3.33811e-9
python3 lemmaS_F2.py --J 6.0273e-12 8.4673e-12 1.8708e-12 --handoff 2e-7 1e-4 0.2 --sqLam 5.3075e-8 --alpha1 1.5e-4 --core 0.11 0.52 --ramp 0.05 0.9 --mhr 0 --slopes 8.6369e-19 1.2134e-18 2.6808e-19 --rho 2.693e-3 --M1s 0.52 --J4 6.2362e-8 --sJ4 8.9363e-15 --rampl4 --rs33 ibp2 --l4 --abar $ABAR --beta $BETA --tag F10p > F10p_run.out 2>&1
RT=$RT ABAR=$ABAR BETA=$BETA python3 lemmaS_F4B.py F10p_M0.52 > F4B_F10p.out 2>&1
AL=1.5e-4 ABAR=$ABAR BETA=$BETA python3 lemmaS_Bbox.py 1e-5 0.15 > B10_a1.5e-4_e1e-5.out 2>&1
python3 rigor/capture_F2.py > rigor/cap_run.out 2>&1
python3 rigor/weights_exact.py > rigor/weights_exact.out 2>&1
python3 rigor/ramp_arb.py > rigor/ramp_arb.out 2>&1
python3 rigor/assembly_arb.py > rigor/assembly_arb.out 2>&1
