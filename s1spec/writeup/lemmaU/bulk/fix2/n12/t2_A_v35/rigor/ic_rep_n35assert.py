#!/usr/bin/env python3
"""Task 2 / A4 (revision r2, CONTRACTS_r2.md): the initial datum of (A) from the new hand-off threshold (H2f).
Audit (r1, unchanged): step 4 certifies (H2) on the shells [14, 60] only (STEP4_RESULT.md l.33, camp item 2); the (H3) windows reach 379.7;
nothing on [379.5, 700].  r2 (STATEMENTS_r2.tex Remark 'initial datum', change R2):
   (H2f)  G_j = int_{ln 60}^{ln 700} (28/y) (d_s^j (f - U_1))^2 ds <= H2f_j  at t_1,     E^A_j(tau_1) <= H2_j + H2f_j   (eta_A <= 1, z^A_j(tau_1) = d_s^j(f - U_1)).
No H3 term, no zone-III term.  (H2f) is a HYPOTHESIS here (label 'hyp'): step 4 has not certified it (task 6).
Writes E_tau1 under '<H2 setting>/<name>' for H2 settings tight (1e-8, 1e-6, 1e-3) and v31 (2e-7, 1e-4, 0.2).
The r1 version (H3 + zone III) is kept as rigor/ic_rep_pre_r2.py; its outputs are in probe/.
Output: --out (default out/ic_<name>.json).   Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/ic_rep.py --h2f 5e-7,5e-6,1e-5 --name design"""
import os, sys, json, argparse
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, Q_, arb
H2S = {'tight': ('1e-8', '1e-6', '1e-3'), 'v31': ('2e-7', '1e-4', '0.2')}
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--h2f', required=True, help='a,b,c: (H2f) thresholds, levels 0,1,2')
    ap.add_argument('--name', required=True); ap.add_argument('--out', default=None); a = ap.parse_args()
    h2f = tuple(x.strip() for x in a.h2f.split(',')); assert len(h2f) == 3
    out = dict(step4_H2_range='[14, 60] for f - U_1 (fix B0: step-4 v33 stage h2u_v33_arb.json (arb row; supersedes the float h2u_v33.json) = v3.2 shells of w\' + certified d, b0_dprof.json); (H3) windows on [50, 379.7]; (H2f) on [60, 700] (step-4 v33 h2f_result_arb.json, arb row; supersedes the float h2f_result.json)',
               name=a.name, H2f=list(h2f), H2f_label='hyp', datum_rule='E^A_j(tau_1) <= H2_j + H2f_j', E_tau1={})
    for s, h2 in H2S.items():
        E = [A(Q_(h2[j])) + A(Q_(h2f[j])) for j in range(3)]
        out['E_tau1']['%s/%s' % (s, a.name)] = dict(E=E, H2=list(h2), H2f=list(h2f))
        print('E^A(tau_1) %-5s %-7s = %s  (H2 %s + H2f %s)' % (s, a.name, ['%.4g' % L.fup(e) for e in E], h2, h2f))
    L.dump(out, a.out or (L.OUT + '/ic_%s.json' % a.name), [os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py'])
