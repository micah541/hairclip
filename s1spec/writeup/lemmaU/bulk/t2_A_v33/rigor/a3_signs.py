#!/usr/bin/env python3
"""Task 2 / A3 (part 1): signs of the table entries used by Lemma P on the STORED symbolic forms (lemmaS_F1.div_form(j), base_ibp copy).
weights_exact.table_exact stores |coef|; here the same expansion WITHOUT abs, for (2,2,0) on (j+1,j+1) and (3,2,0) on (j,j+1), j = 0,1,2.
Lemma P needs: +4 int h z_{j+1}^2 and +2 int g z_j z_{j+1}, h = (Xi a)'', g = (Xi a)''' (Xi = eps*omega*eta^2, a = A2), i.e. both signed coefficients
positive in the convention of the stored form, which IS the contribution to dE_j/dtau (lemmaS_prop1_sym.py l.3: T_j = -2 om z_j d^j(eps L G));
check of the convention: the dissipation entry (0,2,0) on (n,n) is -2.  Then case (i) of Lemma P gives the net point term -g_c z_j(y*)^2 <= 0 in dE_j/dtau,
case (ii) leaves +g_c z_j(y*)^2 (unfavourable), as stated.
Output: out/a3_signs.json.  Exact sympy Rationals (no floating point)."""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX
import sympy as sp
F, P = WX.F, WX.P
def signed(Q, pair, want):
    c = Q[pair]; s = P.s
    e = sp.expand(sp.simplify(c.subs(P.ep, sp.exp(-4*s)).doit()*sp.exp(4*s)))
    gens = list(e.atoms(sp.Derivative)) + [f for f in e.atoms(sp.Function) if not isinstance(f, sp.Derivative) and f.func.__name__ in ('om', 'A0', 'A1', 'A2')]
    P_ = sp.Poly(e, *gens); tot = sp.Integer(0); allm = {}
    for mon, co in P_.terms():
        m = 0; i = r = None
        for g, pw in zip(P_.gens, mon):
            if pw == 0: continue
            if isinstance(g, sp.Derivative): nm, k = g.expr.func.__name__, g.derivative_count
            else: nm, k = g.func.__name__, 0
            if nm == 'om': m = k
            else: i, r = int(nm[1:]), k
        allm[(m, i, r)] = allm.get((m, i, r), 0) + co
    return allm.get(want), {str(k): str(v) for k, v in allm.items()}
if __name__ == '__main__':
    out = {}
    for j in (0, 1, 2):
        Qd, ok = F.div_form(j); assert ok
        h, hall = signed(Qd, (j + 1, j + 1), (2, 2, 0)); g, gall = signed(Qd, (j, j + 1), (3, 2, 0)); dn, _ = signed(Qd, (j + 2, j + 2), (0, 2, 0))
        out[str(j)] = dict(h_220_on_jp1jp1=str(h), g_320_on_jjp1=str(g), dissipation_020_on_nn=str(dn), signed_pair_jp1=hall, signed_pair_jjp1=gall,
                           lemmaP_signs_ok=bool(h == 4 and g == 2 and dn == -2))
        print('level', j, ': (2,2,0) on (%d,%d) = %s ; (3,2,0) on (%d,%d) = %s' % (j + 1, j + 1, h, j, j + 1, g), flush=True)
    L.dump(dict(res=out), L.OUT + '/a3_signs.json', [L.BASE + '/rigor/weights_exact.py', L.BASE + '/lemmaS_F1.py', L.BASE + '/lemmaS_prop1_sym.py', os.path.abspath(__file__)])
