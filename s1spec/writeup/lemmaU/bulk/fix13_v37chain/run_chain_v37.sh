#!/bin/bash
# fix13_v37chain (tag v37; copy of fix6_v36/run_chain_v36.sh, default verifier fin13v4c, n1 chart per V36_VER) driver: [n1_chart on fin13v4b] -> base_ibp reference chain (verifier argv) -> R3 constants -> t2 (A) -> t4 Lemma U -> zone III,
# run entirely inside fix6_v36/.  Modelled on fix4_r3/run_chain.sh (V33_INDEX.md section 4.6 order).
# Usage: V36_VER=fin13v4b|fin12f4 V36_R3=1|0 V36_R2=1|0 bash run_chain_v36.sh <snapshot-name>
#   outputs in fix6_v36/{t2_A/out,t2_A/base_ibp,t4_cert/out,t5_zone3/out,r3,n1/out}; snapshot copied to fix6_v36/runs/<name>/ (refuses an existing name).
# One process at a time, nice 19, PYTHONHASHSEED=0, python3 -B (no .pyc); /usr/bin/time records wall time and max RSS per stage in runs/<name>/stages.txt.
# Reads (never writes) outside fix6_v36: the verifier outputs (v36cfg.py), step-4 v33 b0_dprof.json (and n1_chart.json in fin12f4 mode), phase1/phase2a
# Lemma-U inputs, lemmaO_v32 constants, step3nl/step3/step4-t9h modules (n1_chart), the v35 zone3_interface.json (initial input of assemble pass 1).
set -e
NAME=${1:?snapshot name}
N=/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain
Z3V35=/home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export V36_VER=${V36_VER:-fin13v4c} V36_R3=${V36_R3:-1} V36_R2=${V36_R2:-1}
SNAP=$N/runs/$NAME; test ! -e $SNAP; mkdir -p $SNAP
ST=$SNAP/stages.txt
echo "V36_VER=$V36_VER V36_R3=$V36_R3 V36_R2=$V36_R2 start $(date -Is)" > $SNAP/env.txt
T="/usr/bin/time -a -o $ST -f %C|%e_s|%M_KB"
R="$T nice -n 19 python3 -B"
python3 -I -B $N/v36cfg.py > $SNAP/v36cfg.json
rm -f $N/t2_A/out/* ; find $N/t4_cert/out $N/t5_zone3/out -type f ! -name ramp_U.json -delete
# 0. N1 chart terms on the fin13v4b sups (fin12f4 mode reads the pinned step-4 file instead)
if [ "$V36_VER" != fin12f4 ]; then mkdir -p $N/n1/out; (cd $N/n1 && rm -f out/n1_chart_${V36_VER}.json && $R n1_chart.py > out/n1_chart_${V36_VER}.log 2>&1); echo "n1       done $(date +%T)"; fi
# 1. base_ibp reference chain (lemmaS_F10p_M0.52.json args = the verifier inputs (A) consumes; cap, weights_exact, assembly_arb reference)
(cd $N/t2_A/base_ibp && rm -f lemmaS_F10p_M0.52.json rigor/cap_F10p.json rigor/lemmaS_rigcap_M0.52.json rigor/weights_exact.json rigor/assembly_arb.json \
  && $R run_F10p_ibp.py > F2_v36.log 2>&1 && $R rigor/capture_F2.py > rigor_capture_F2_v36.log 2>&1 \
  && $R rigor/weights_exact.py > rigor_weights_exact_v36.log 2>&1 && $R rigor/assembly_arb.py > rigor_assembly_arb_v36.log 2>&1)
echo "base_ibp done $(date +%T)"
# 2. R3 constants (record; the chain imports the module)
if [ "$V36_R3" = 1 ]; then (cd $N/r3 && rm -f r3_forcing.json r3_symcheck.json && $R r3_forcing.py > r3_forcing.log 2>&1 && $R r3_symcheck.py > r3_symcheck.log 2>&1); echo "r3       done $(date +%T)"; fi
# 3. (A): t2 chain
$T bash $N/t2_A/run_t2.sh
# 4. Lemma U first pieces
cd $N/t4_cert
$R u1_ghat.py > out/u1_ghat.log 2>&1
$R u34_sweeps.py --case primary > out/u34_primary.log 2>&1
$R u2_K.py --fresh > out/u2_K.log 2>&1
rm -f out/*.partial.json
AI=$N/t2_A/out/A_interface.json
# 5. assemble pass 1 (v35 zone3 interface), zone3, pass 2, zone3, pass 3
$R assemble_U.py --A $AI --zone3 $Z3V35 --out U_cert_v36_pass1.json > out/assemble_pass1.log 2>&1
$R make_ucert_for_zone3.py out/U_cert_v36_pass1.json out/U_cert_for_zone3_v36_pass1.json
(cd $N/t5_zone3 && $R zone3.py --ucert ../t4_cert/out/U_cert_for_zone3_v36_pass1.json --ucert_full $N/t4_cert/out/U_cert_v36_pass1.json > out/zone3_pass1.log 2>&1 && cp out/zone3_interface.json out/zone3_interface_pass1.json)
$R assemble_U.py --A $AI --zone3 $N/t5_zone3/out/zone3_interface_pass1.json --out U_cert_v36_pass2.json > out/assemble_pass2.log 2>&1
$R make_ucert_for_zone3.py out/U_cert_v36_pass2.json out/U_cert_for_zone3_v36.json
(cd $N/t5_zone3 && $R zone3.py --ucert ../t4_cert/out/U_cert_for_zone3_v36.json --ucert_full $N/t4_cert/out/U_cert_v36_pass2.json > out/zone3.log 2>&1)
$R assemble_U.py --A $AI --zone3 $N/t5_zone3/out/zone3_interface.json --out U_cert_v36_pass3.json > out/assemble_pass3.log 2>&1
echo "chain    done $(date +%T)"
for d in t2_A t4_cert t5_zone3; do (cd $N/$d/out && sha256sum *.json > SHA256SUMS); mkdir -p $SNAP/$d; cp -a $N/$d/out/. $SNAP/$d/; done
mkdir -p $SNAP/base_ibp/rigor; (cd $N/t2_A/base_ibp && cp -a lemmaS_F10p_M0.52.json F2_v36.log rigor_capture_F2_v36.log rigor_weights_exact_v36.log rigor_assembly_arb_v36.log $SNAP/base_ibp/ \
  && cp -a rigor/cap_F10p.json rigor/lemmaS_rigcap_M0.52.json rigor/weights_exact.json rigor/assembly_arb.json $SNAP/base_ibp/rigor/)
if [ "$V36_R3" = 1 ]; then mkdir -p $SNAP/r3; cp -a $N/r3/r3_forcing.json $N/r3/r3_forcing.log $N/r3/r3_symcheck.json $N/r3/r3_symcheck.log $SNAP/r3/; fi
if [ "$V36_VER" != fin12f4 ]; then mkdir -p $SNAP/n1; cp -a $N/n1/out/n1_chart_${V36_VER}.json $N/n1/out/n1_chart_${V36_VER}.log $SNAP/n1/; fi
echo "end $(date -Is)" >> $SNAP/env.txt
echo "snapshot $SNAP"
