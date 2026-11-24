#!/bin/bash
# H1 uniformly on slab 180 (t_c in [-1, 1]): 16 exact points t_k = -1 + (2k+1)/16 (Krawczyk IFT + Lambda), and for each the two
# derivative balls [t_k - 1/16, t_k] and [t_k, t_k + 1/16] (h1dt.py, reference t_k).  usage: run_h1slab.sh K0 K1  (k range of one worker)
cd /home/claude/s1spec/step4/p2/t9h
export COHERENT=1 S1PROF=/home/claude/s1_tight/ext_tip/prof_tip.txt S1TT=10.249077651120778,10.249077651121361
for k in $(seq $1 $2); do
  tk=$(python3 -c "print(-1 + (2 * $k + 1) / 16)")
  ref=$(python3 -c "
a = (1.04655427e-4, 0.99999491111, -1.17386340e-4); b = (1.00712248e-4, 0.99999510274, -1.12963590e-4); f = ($tk + 1) / 1.9990234375
print(','.join(repr(x + f * (y - x)) for x, y in zip(a, b)))")
  out=../out/t9_h1_p${tk}_coh_s180.json
  [ -f $out ] || H1TAG=_s180 H1REF=$ref nice -n 10 python3 -I h1.py p$tk > ../out/h1slab/h1_p${tk}.log 2>&1
  for sgn in -1 1; do
    tc=$(python3 -c "print($tk + $sgn * 0.03125)")
    [ -f ../out/h1slab/dt_${tc}_0.03125.json ] || [ -n "$NODT" ] || TREF=$tk nice -n 10 python3 -I h1dt.py $tc 0.03125 $out > ../out/h1slab/dt_${tc}.log 2>&1
  done
done
