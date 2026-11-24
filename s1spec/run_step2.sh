#!/usr/bin/env bash
# Reproduces every rigorous step of the computer-assisted proof of linear stability of S1 (step 2).
# Requirements: CAPD 6.x (CAPDCFG=/path/to/capd-config), python-flint >= 0.9, numpy; cvxpy+clarabel only to
# regenerate the certificates (not part of the trust base).  Expects the step-1 bundle in ../s1_cap.
set -euo pipefail
CAPDCFG=${CAPDCFG:-capd-config}
ROOT=$(cd "$(dirname "$0")" && pwd)
S1=${S1:-$ROOT/../s1_cap}
export PYTHONPATH="$S1:$ROOT/lin:$ROOT/cap:${PYTHONPATH:-}"

echo "== 1. Profile enclosures of S1 on the shooting bracket (30000 boxes)"
cd "$ROOT/apr"
g++ -O2 "$S1/s1curve.cpp" $($CAPDCFG --cflags --libs) -o s1curve
python3 "$S1/gen_input.py" 8 0.00015345849609375 0.00015346679687499998 8.5 in_br.txt > /dev/null
./s1curve in_br.txt 10.30 30000 > curve_fine.txt

echo "== 2. Numerical-range certificates (a priori trapezoid)"
for c in cert_b0.25_th0.0_tau0.03.json cert_b0.25_th0.48_tau0.03.json; do
  python3 verify.py "$c" | tee "ver_$c.log" | grep -E "Lam|FAILURES"
  grep -q "FAILURES: 0" "ver_$c.log" || { echo "certificate $c FAILED"; exit 1; }
done

echo "== 3. Evans-function winding numbers along the contour"
cd "$ROOT/cap"
python3 gen_vf.py
g++ -O2 evans_seg.cpp $($CAPDCFG --cflags --libs) -o evans_seg
rm -f contour_ckpt.json
python3 contour.py | tee contour_rerun.log | tail -20
python3 winding_cert.py contour_result.json

echo "== 4. Lyapunov-Perron bounds on the interior (analyticity)"
python3 interior.py | tail -3
echo "DONE"
