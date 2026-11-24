#!/usr/bin/env python3
"""(B) level 3 on the physically fixed strip: Ghat_33 (zero-jet + L4/interpolation extras) and gamma'_3 at eps_B = yE^-4; requirement on hand-off E3(tau1) for K3."""
import sys, pickle, math, json
sys.path.insert(0, '.')
import sympy as sp
import lemmaS_F1 as F, lemmaS_prop1_sym as P, lemmaS_prop1_bounds as PB, lemmaS_L4terms as L4, run_F1 as R
def div_form3(Q, n):
    out = {}
    def add(k, v): k = tuple(sorted(k)); out[k] = out.get(k, 0) + v
    for (a, b), c in Q.items():
        if b < n or (a, b) == (n, n): add((a, b), c)
        else:
            if a == n - 1: add((a, a), -sp.Rational(1, 2)*P.D(c))
            else: add((a + 1, n - 1), -c); add((a, n - 1), -P.D(c))
    return {k: sp.expand(v) for k, v in out.items() if v != 0}
if __name__ == '__main__':
    M1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6; M0 = 0.11; gam = 1.0
    Q = pickle.load(open('prop1_forms3.pkl', 'rb'))[3][0]
    Qd = div_form3(Q, 5); tab3 = PB.table(Qd)
    box = R.box(0, 0, 0, M0=M0, M1=M1); Ab = F.A_bound_set('quasi', box); amin = (1 + M0**2)**-2
    YC = {k: 0.0 for k in range(6)}
    for (a, b), mons in tab3.items():
        if (a, b) == (5, 5): continue
        gC = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
    T = {3: tab3}
    add, log, dis, spare = L4.l4_extra_YC(T, 3, M0, M1, gam, fallback={2: float(sys.argv[2]) if len(sys.argv) > 2 else 3.0, 3: float(sys.argv[3]) if len(sys.argv) > 3 else 3.0})
    YC2 = {k: YC[k] + add.get(k, 0) for k in YC}
    G0 = L4.ghat_from_YC(YC, 3, amin); G1 = L4.ghat_from_YC(YC2, 3, amin)
    print('level 3: Ghat_33 zero-jet %.0f ; with L4 extras %.0f ; dissipation %.3f of %.3f' % (G0[3], G1[3], dis, spare))
    sq = 1.45e-7; abar = 7.57*sq; cmin = 0.25 - abar - abar/50
    for yE in (50., 48.):
        gp = gam*cmin - 4*abar/yE - yE**-4*G1[3]
        print(f'  yE={yE:g}: gamma_3\' = {gp:.4f}', '' if gp > 0 else 'INVALID')
    json.dump(dict(G33_zero=float(G0[3]), G33_L4=float(G1[3]), YC={k: float(v) for k, v in YC2.items()}, M1=M1), open('lemmaS_B3.json', 'w'), indent=1)
