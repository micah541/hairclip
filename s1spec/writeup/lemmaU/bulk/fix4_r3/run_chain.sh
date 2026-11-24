#!/bin/bash
# fix4_r3 (tag v36pre) driver: t2 -> t4 -> zone3 chain of v35 (V33_INDEX.md section 4.6), run entirely inside fix4_r3/.
# Usage: bash run_chain.sh <snapshot-name>      (outputs in fix4_r3/{t2_A,t4_cert,t5_zone3}/out; afterwards copied to fix4_r3/runs/<name>/)
# One process at a time, nice 19, PYTHONHASHSEED=0, no .pyc; /usr/bin/time records wall time and max RSS of every stage in runs/<name>/stages.txt.
# Reads (never writes) outside fix4_r3: step-4 v33 outputs b0_dprof.json, n1_chart.json; verifier fin12f4 tables; phase1 / phase2a Lemma-U inputs;
# lemmaO_v32 constants; the v35 zone3_interface.json (initial zone-III input of assemble pass 1, as in the v35 chain).
set -e
NAME=${1:?snapshot name}
N=/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3
Z3V35=/home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
SNAP=$N/runs/$NAME; test ! -e $SNAP; mkdir -p $SNAP
ST=$SNAP/stages.txt
T="/usr/bin/time -a -o $ST -f %C|%e_s|%M_KB"
R="$T nice -n 19 python3 -B"
rm -f $N/t2_A/out/* ; find $N/t4_cert/out $N/t5_zone3/out -type f ! -name ramp_U.json -delete
# 0. R3 constants (only present after the R3 edit)
if [ -f $N/r3/r3_forcing.py ]; then (cd $N/r3 && rm -f r3_forcing.json r3_symcheck.json && $R r3_forcing.py > r3_forcing.log 2>&1 && $R r3_symcheck.py > r3_symcheck.log 2>&1); echo "r3       done $(date +%T)"; fi
# 1. (A): t2 chain (run_t2.sh, every python call through nice 19; time per call recorded by wrapping the whole script)
$T bash $N/t2_A/run_t2.sh
# 2. Lemma U first pieces
cd $N/t4_cert
$R u1_ghat.py > out/u1_ghat.log 2>&1
$R u34_sweeps.py --case primary > out/u34_primary.log 2>&1
$R u2_K.py --fresh > out/u2_K.log 2>&1
rm -f out/*.partial.json
AI=$N/t2_A/out/A_interface.json
# 3. assemble pass 1 (v35 zone3 interface), zone3, pass 2, zone3, pass 3
$R assemble_U.py --A $AI --zone3 $Z3V35 --out U_cert_r3_pass1.json > out/assemble_pass1.log 2>&1
$R make_ucert_for_zone3.py out/U_cert_r3_pass1.json out/U_cert_for_zone3_r3_pass1.json
(cd $N/t5_zone3 && $R zone3.py --ucert ../t4_cert/out/U_cert_for_zone3_r3_pass1.json --ucert_full $N/t4_cert/out/U_cert_r3_pass1.json > out/zone3_pass1.log 2>&1 && cp out/zone3_interface.json out/zone3_interface_pass1.json)
$R assemble_U.py --A $AI --zone3 $N/t5_zone3/out/zone3_interface_pass1.json --out U_cert_r3_pass2.json > out/assemble_pass2.log 2>&1
$R make_ucert_for_zone3.py out/U_cert_r3_pass2.json out/U_cert_for_zone3_r3.json
(cd $N/t5_zone3 && $R zone3.py --ucert ../t4_cert/out/U_cert_for_zone3_r3.json --ucert_full $N/t4_cert/out/U_cert_r3_pass2.json > out/zone3.log 2>&1)
$R assemble_U.py --A $AI --zone3 $N/t5_zone3/out/zone3_interface.json --out U_cert_r3_pass3.json > out/assemble_pass3.log 2>&1
echo "chain    done $(date +%T)"
for d in t2_A t4_cert t5_zone3; do (cd $N/$d/out && sha256sum *.json > SHA256SUMS); mkdir -p $SNAP/$d; cp -a $N/$d/out/. $SNAP/$d/; done
if [ -f $N/r3/r3_forcing.json ]; then mkdir -p $SNAP/r3; cp -a $N/r3/r3_forcing.json $N/r3/r3_forcing.log $N/r3/r3_symcheck.json $N/r3/r3_symcheck.log $SNAP/r3/; fi
echo "snapshot $SNAP"
