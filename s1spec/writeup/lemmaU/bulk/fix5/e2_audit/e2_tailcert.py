"""E2 audit: the structural (degree) conditions of tail_cert.py on vp/rhs_iv.py with the E2 additions switched on (both variants).
Mechanical: the fin13v4 tail_cert.run is executed unmodified on the patched rhs.  usage: python3 -B e2_tailcert.py   (env as run_v4.sh)"""
import sys, os, json
D = os.path.dirname(os.path.abspath(__file__)); V = '/home/claude/s1spec/step3nl_v4/verifier'
sys.path.insert(0, V); sys.path.insert(0, D + '/vp')
import rhs_iv as R, tail_cert as TC
assert R.__file__.startswith(D + '/vp') and TC.__file__.startswith(V) and TC.R is R
_c = R._c; out = {}
base = dict(dab=0.013114, n3=7.69699, nh2=13.1085, ngab=8.0286, srw=4.4292e-3, cw0=1.00114, cw1=33.0322, r2=32.4022, xwc=5.71e-6)
for name, extra in (('fin13v4', None), ('E2fix_CS', dict(nf12=2.06197, fe2=13.99)), ('E2fix_IBP', dict(ndf=12.34096, fe2=6.07, fe2A=47.4))):
    R.E2X = None if extra is None else {k: _c(v) for k, v in {**base, **extra}.items()}
    res = {}
    for mode in ('radial', 's', 'joint'):
        H, L, N = TC.run(mode); res[mode] = dict(H=H, L=L, NEG=N)
    H2, L2, N2 = TC.run('joint', floors=False); res['joint_nofloor'] = dict(H=H2, L=L2, NEG=N2)
    ok = (res['radial']['H'][1] <= 2 and res['radial']['L'][1] <= 2 and res['s']['H'][0] >= 2 and res['s']['L'][1] <= 2 and res['s']['NEG'] == (2, 2)
          and H2[0] >= 2 and L2[0] >= 2 and N2 == (2, 2))
    res['PASS'] = bool(ok); out[name] = res; print(name, res)
json.dump(out, open('e2_tailcert.json', 'w'), indent=1)
