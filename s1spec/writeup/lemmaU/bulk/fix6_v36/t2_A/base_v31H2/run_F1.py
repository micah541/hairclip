#!/usr/bin/env python3
"""CLI assembly of the Lemma S shell estimates (fix plan F1/F4).
  python3 run_F1.py --J J0 J1 J2 --J4 4.3e-7 --handoff 5e-7 3e-4 20 [--sqLam 2.49e-7] [--alpha1 5e-4]
                    [--mode lin|quasi] [--box M2 M3 M4] [--gam 1] [--tag name]
J3 is computed as sqrt(J2*J4)+0.8477*J2 unless --J3 is given.  Output: table + lemmaS_F1_<tag>.json.
Entries with gamma' <= 0 are INVALID (the energy inequality does not close) and are flagged."""
import argparse, json, math, sys
sys.path.insert(0, '.')
import lemmaS_F1 as F

def box(M2, M3, M4, Umax=0.05, M0=0.1, M1=1.1):
    return dict(M0=M0, M1=M1, M2=M2, M3=M3, M4=M4, M5=M4*1.5, M6=M4*2, M7=M4*2.5, Umax=Umax, p=1.0246)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--J', nargs=3, type=float, default=[4.1e-11, 5.8e-11, 1.3e-11])
    ap.add_argument('--J3', type=float, default=None)
    ap.add_argument('--J4', type=float, default=4.3e-7)
    ap.add_argument('--handoff', nargs=3, type=float, default=[5e-7, 3e-4, 20.0])
    ap.add_argument('--sqLam', type=float, default=2.49e-7)
    ap.add_argument('--alpha1', type=float, default=5e-4)
    ap.add_argument('--mode', choices=['lin', 'quasi'], default='lin')
    ap.add_argument('--box', nargs=3, type=float, default=[1.1, 2, 4])
    ap.add_argument('--M0', type=float, default=0.1)
    ap.add_argument('--M1', type=float, default=1.1)
    ap.add_argument('--gam', type=float, nargs='+', default=[1.0])
    ap.add_argument('--tag', default='run')
    a = ap.parse_args()
    J3 = a.J3 if a.J3 is not None else math.sqrt(a.J[2]*a.J4) + 0.8477*a.J[2]
    Jl = tuple(a.J) + (J3,)
    tables = F.build_tables(); out = {}
    for gam in a.gam:
        c, G = F.constants(gam, box(*a.box, M0=a.M0, M1=a.M1), a.mode, tables=tables)
        r = F.assemble(G, Jl, a.sqLam, gam, a.alpha1, handoff=tuple(a.handoff))
        valid = all(float(p['gamma_prime']) > 0 for p in r['parts'])
        print(f"mode={a.mode} gam={gam} box={a.box} J={a.J} J3={J3:.3g} J4={a.J4}")
        print('  Ghat_jj =', [round(float(c[f'Ghat_{j}'][str(j)]), 1) for j in range(3)], " gamma' =", [round(float(p['gamma_prime']), 4) for p in r['parts']])
        print('  E =', [f'{e:.3g}' for e in r['E']], ' S0_v = %.4f  S1_v = %.3f' % (r['S0'], r['S1']), '' if valid else '  ** INVALID: gamma_prime<=0 at some level **')
        out[f'gam{gam}'] = dict(consts=c, result=r, valid=valid, inputs=vars(a), J3=J3)
    json.dump(out, open(f'lemmaS_F1_{a.tag}.json', 'w'), indent=1, default=str)
if __name__ == '__main__': main()
