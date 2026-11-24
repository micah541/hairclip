#!/usr/bin/env bash
# Reproduces the tightened step-1 enclosure of S1 (needs CAPD 6.x, python-flint >= 0.9; ../s1 holds tail.py, geometry.py, s1proof, s1curve).
set -euo pipefail
CAPDCFG=${CAPDCFG:-capd-config}
S1DIR=${S1DIR:-../s1}
g++ -O2 $S1DIR/s1proof.cpp $($CAPDCFG --cflags --libs) -o s1proof_bin
g++ -O2 $S1DIR/s1curve.cpp $($CAPDCFG --cflags --libs) -o s1curve_bin
g++ -O2 prof/s1prof.cpp $($CAPDCFG --cflags --libs) -o s1prof_bin
echo "== 1. exact linear tail, cubic remainder, Jhat"
python3 tail2.py
python3 jhat_rig.py
echo "== 2. tight bisection (each sign certified by CAPD)"
sed -i "s|^BIN = .*|BIN = \"./s1proof_bin\"|" bracket2.py
python3 bracket2.py
echo "== 3. trajectory enclosure + geometry on the tight bracket set"
python3 -c "
import json
from gen_input2 import build, write
b = json.load(open('bracket2.json'))
write('in_br2.txt', *build('8', repr(b['eps_a']), repr(b['eps_b']), 8.5, emax=1.62e-4*(1+1e-12))[:6])"
./s1curve_bin in_br2.txt 10.30 8000 > curve_br2.txt
(cd $S1DIR && python3 geometry.py $OLDPWD/curve_br2.txt 100)
echo "== 4. tight profile boxes for the downstream certificates"
./s1prof_bin in_br2.txt 10.249077651519443 51246 > prof/prof_tight.txt
