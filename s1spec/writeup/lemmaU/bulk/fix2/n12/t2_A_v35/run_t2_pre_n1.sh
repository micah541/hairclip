#!/bin/bash
# Fix B0 copy (v34) of the task 2 (r2) driver; assembly_rep adds d = (f - U) - w' (step-4 v33 b0_dprof.json).
# Task 2 (r2) driver: (A) with R5' (L^2 form), H-J-W, the (H2f) datum, point terms, localised identities; writes out/A_interface.json.
# One process at a time, nice 19, < 100 MB.  Stages:
#   hjw_rep -> ic_rep (4 H2f settings) -> a3_signs -> weights_t2 (box (0.15,0.52); fallback (0.20,0.60)) -> ramp_rep (A3 route; and the phase-2a route for dE)
#   -> outer_L2 (R5') -> assembly_rep (design at N2ref, 4 N2ref, 0; lo, hi, probe at N2ref; design with the phase-2a ramp route) -> export_t2
set -e; cd /home/claude/s1spec/writeup/lemmaU/bulk/fix1/b0/t2_A_v34
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
R="nice -n 19 python3 -B"
$R rigor/hjw_rep.py > out/hjw_rep.log 2>&1;                                   echo "hjw      done $(date +%T)"
for s in "design 1e-6,5e-6,1e-5" "lo 8e-7,5e-6,1e-5" "hi 1.2e-6,5e-6,1e-5" "probe 3.5e-7,1.1e-6,2.7e-6"; do set -- $s
  $R rigor/ic_rep.py --h2f $2 --name $1 > out/ic_$1.log 2>&1; done;             echo "ic       done $(date +%T)"
$R rigor/a3_signs.py > out/a3_signs.log 2>&1;                                 echo "signs    done $(date +%T)"
$R rigor/weights_t2.py > out/weights_t2.log 2>&1
$R rigor/weights_t2.py --seg-box 20/100 60/100 --tag fb > out/weights_t2_fb.log 2>&1; echo "weights  done $(date +%T)"
$R rigor/ramp_rep.py --N 1200 > out/ramp_rep.log 2>&1
$R rigor/ramp_rep.py --N 1200 --pm3-ibp --out out/ramp_rep_pm3ibp.json > out/ramp_rep_pm3ibp.log 2>&1; echo "ramp     done $(date +%T)"
$R rigor/outer_L2.py > out/outer_L2.log 2>&1;                                 echo "outer    done $(date +%T)"
A="rigor/assembly_rep.py --hjw out/hjw_rep.json --dprof /home/claude/s1spec/step4/p2/out_v33/b0_dprof.json"   # fix B0
$R $A --ic out/ic_design.json --h3key design --N2 0.1 0.2 5.0 --out out/asm_design_N2ref.json  > out/asm_design_N2ref.log 2>&1
$R $A --ic out/ic_design.json --h3key design --N2 0.4 0.8 20.0 --out out/asm_design_N2ref4.json > out/asm_design_N2ref4.log 2>&1
$R $A --ic out/ic_design.json --h3key design --N2 0 0 0 --out out/asm_design_N20.json          > out/asm_design_N20.log 2>&1
for s in lo hi probe; do $R $A --ic out/ic_$s.json --h3key $s --N2 0.1 0.2 5.0 --out out/asm_${s}_N2ref.json > out/asm_${s}_N2ref.log 2>&1; done
$R $A --ic out/ic_design.json --h3key design --N2 0.1 0.2 5.0 --ramp out/ramp_rep_pm3ibp.json --out out/asm_design_pm3ibp.json > out/asm_design_pm3ibp.log 2>&1
echo "assembly done $(date +%T)"
$R rigor/export_t2.py > out/export_t2.log 2>&1;                               echo "export   done $(date +%T)"
(cd out && sha256sum *.json > SHA256SUMS)
