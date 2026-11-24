#!/usr/bin/env bash
# Reproduces every rigorous step of the computer-assisted proof of the S1 shrinker.
# Requirements: CAPD 6.x (capd-config on PATH or CAPDCFG=/path/to/capd-config), python-flint >= 0.9, numpy.
set -euo pipefail
CAPDCFG=${CAPDCFG:-capd-config}
g++ -O2 s1proof.cpp $($CAPDCFG --cflags --libs) -o s1proof
g++ -O2 s1curve.cpp $($CAPDCFG --cflags --libs) -o s1curve
pass=1
check() { if grep -q "$2" <<< "$1"; then echo "   ok: $3"; else echo "   FAILED: $3"; pass=0; fi; }

echo "== 1. Tail theorem constants (S=8, |eps| <= 1.62e-4)"
python3 -c "
from flint import arb
from tail import radii, gamma1
T=radii(arb(8), arb(1.62e-4*(1+1e-12)))
print('   verified:', T['verified'], ' kappa=%.6g r0=%.4g r1=%.4g delta1=%.4g L=%.4f Jbar=%.4g' % tuple(float(T[k].mid()) for k in ['kappa','r0','r1','delta1','L','Jbar']))
print('   gamma1 =', gamma1(8))
assert T['verified']"

echo "== 2. M(eps_a) > 0 at eps_a = 1.45e-4"
python3 gen_input.py 8 1.45e-4 1.45e-4 8.5 in_a.txt > /dev/null
out=$(./s1proof in_a.txt); echo "$out" | grep -E "^  m =|return time"; check "$out" "m sign: POSITIVE" "M(eps_a) > 0"; check "$out" "RESULT OK" "CAPD validation"

echo "== 3. M(eps_b) < 0 at eps_b = 1.62e-4"
python3 gen_input.py 8 1.62e-4 1.62e-4 8.5 in_b.txt > /dev/null
out=$(./s1proof in_b.txt); echo "$out" | grep -E "^  m =|return time"; check "$out" "m sign: NEGATIVE" "M(eps_b) < 0"; check "$out" "RESULT OK" "CAPD validation"

echo "== 4. Whole interval [eps_a, eps_b]: q(tau1)<0 and transversal first return (continuity of M)"
python3 gen_input.py 8 1.45e-4 1.62e-4 8.5 in_full.txt > /dev/null
out=$(./s1proof in_full.txt); echo "$out" | grep -E "y1\[1\]|return time"; check "$out" "q(tau1) < 0 : YES" "q<0 at tau1 on the whole set"; check "$out" "RESULT OK" "Poincare map validated on the whole set"

echo "== 5. Tight bracket for eps* (signs re-verified at both ends)"
EA=0.00015345849609375; EB=0.00015346679687499998
python3 gen_input.py 8 $EA $EA 8.5 in_ea.txt > /dev/null; out=$(./s1proof in_ea.txt); check "$out" "m sign: POSITIVE" "M > 0 at $EA"
python3 gen_input.py 8 $EB $EB 8.5 in_eb.txt > /dev/null; out=$(./s1proof in_eb.txt); check "$out" "m sign: NEGATIVE" "M < 0 at $EB"

echo "== 6. Geometry of S1 on the bracket set"
python3 gen_input.py 8 $EA $EB 8.5 in_br.txt > /dev/null
./s1curve in_br.txt 10.30 8000 > curve_br.txt
out=$(python3 geometry.py curve_br.txt 100); echo "$out" | sed 's/^/   /'
check "$out" "ALL GEOMETRIC CHECKS PASSED" "double point / embedding / mirror checks"
check "$out" "TAIL/RADIUS CHECKS PASSED" "tail and radius checks"

[ $pass = 1 ] && echo "ALL STEPS PASSED" || { echo "SOME STEP FAILED"; exit 1; }
