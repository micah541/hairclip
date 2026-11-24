#!/bin/sh
# fin13v4b certification chain (sequential, one process, nice 19, no bytecode). Outputs in v4bcert/out.
cd /home/claude/s1spec/step3nl_v4b/v4bcert || exit 1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
O=/home/claude/s1spec/step3nl_v4b/v4bcert/out
mkdir -p out; : > out/chain_status.txt
for st in tc_weighted e0_strip_drift e0_lemma_consts e0_identities r2_strip_ibp; do
  /usr/bin/time -f "$st maxrss_kb=%M wall_s=%e" -a -o out/chain_status.txt nice -n 19 python3 -B $st.py > out/$st.log 2>&1; echo "$st rc=$?" >> out/chain_status.txt
done
# E2 audit fix (promoted from fix5/e2_audit): nonlinear constants at the asserted boxes, modulation/strip constants, verifier JSON, identity check
R="R1:0:8:1e-4:1e-4 R1w:0:8:1e-3:1e-3 R2a:8:12:1e-3:1e-3 R2b:12:16.25:1e-3:1e-3 R3:16.25:34.25:0.02:0.02"
/usr/bin/time -f "e2_nlin_regions maxrss_kb=%M wall_s=%e" -a -o out/chain_status.txt nice -n 19 python3 -B e2_nlin_regions.py 2.5e-4 1 $O/nl_regions_a2.5e-4.json $R > out/e2_nlin_regions.log 2>&1; echo "e2_nlin_regions rc=$?" >> out/chain_status.txt
for st in e2_consts e2_v4b_constants e2_identity_check; do
  /usr/bin/time -f "$st maxrss_kb=%M wall_s=%e" -a -o out/chain_status.txt nice -n 19 python3 -B $st.py > out/$st.log 2>&1; echo "$st rc=$?" >> out/chain_status.txt
done
