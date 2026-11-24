#!/usr/bin/env python3
"""FIXPLAN2: one shell system levels 0-2, core [28,60] (box M0c,M1c), ramp [14,28] (box M0r,M1r), no jets M2..M4 (divergence form).
Ramp commutators at order j (L1): pairs with z_{j+1} are majorised by m = L0*eps*omega*t^3, then
   int m z_{j+1}^2 = -int m z_j z_{j+2} + 1/2 int m'' z_j^2 <= 0.5*amin*int Theta z_{j+2}^2 + int (m^2/(2 amin Theta) + |m''|/2) z_j^2.
Only J0,J1,J2 are used.  python3 lemmaS_F2.py --J J0 J1 J2 --handoff E0 E1 E2 --sqLam x --core M0 M1 --ramp M0 M1"""
import sys, math, json, argparse
sys.path.insert(0, '.')
import lemmaS_F1 as F, run_F1 as R
from flint import arb
STIR, BELL = dict(F.STIR), dict(F.BELL)
STIR[4] = {1: 1, 2: 7, 3: 6, 4: 1}; BELL[4] = 15
import os
DIAG = bool(os.environ.get('DIAG'))
import sympy as sp
_t = F.tt
def _remq_fun(gam, L):
    e2 = F.Sx**2; out = {}
    for m in (1, 2, 3):
        P = sum(sp.binomial(m, l)*sp.diff(e2, _t, l)/sp.Float(L)**l*sp.Float(-gam)**(m - l) for l in range(1, m + 1))
        Q = sp.cancel(sp.expand(P)/_t**3)
        out[m] = sp.lambdify(_t, Q, modules=[{}])
    return out
from math import factorial
def _expand(a, b, d):
    """exact IBP reduction of int c^{(d)} z_a z_b to squares z_k^2 (a<=b): returns {(k,d'):weight}"""
    if a == b: return {(a, d): 1.0}
    if b == a + 1: return {(a, d + 1): 0.5}
    out = {}
    for key, v in list(_expand(a + 1, b - 1, d).items()) + list(_expand(a, b - 1, d + 1).items()): out[key] = out.get(key, 0) + v
    return out
def ramp_level(j, T, Abr, Ghr, YCr, Rm, gam, tmin=1e-3):
    """L1 in table form: rem-parts of pairs touching z_{j+1} are reduced by repeated exact IBP (z_a z_b -> z_{a+1}z_{b-1}, c'), NO Young, only J_k, k<=j.
    Derivatives of the rem coefficient by Leibniz: d^d[eps*rem_m*A^(r)] bounded termwise (|(d/ds - gam)^q R_m| * omega, (-4)^p, A^(r+u))."""
    n = j + 2; L = Rm.L; e2 = F.Sx**2
    pairs = [(a, b, m) for (a, b), m in T[j].items() if (a, b) != (n, n)]
    Dq = {}
    for m in (1, 2, 3, 4):
        R = sum(sp.binomial(m, l)*sp.diff(e2, _t, l)/sp.Float(L)**l*sp.Float(-gam)**(m - l) for l in range(1, m + 1))
        f = sp.expand(R); Dq[m] = []
        for q_ in range(4):
            Dq[m].append(sp.lambdify(_t, f, modules=[{}])); f = sp.expand(sp.diff(f, _t)/sp.Float(L) - sp.Float(gam)*f)
    def cd(mons, d, t):
        tot = arb(0)
        for (m, i, r), co in mons.items():
            if m == 0: continue
            for p_ in range(d + 1):
                for q_ in range(d + 1 - p_):
                    u = d - p_ - q_; mult = factorial(d)//(factorial(p_)*factorial(q_)*factorial(u))
                    tot += arb(abs(float(co))*mult*4**p_*Abr[(i, r + u)])*abs(arb(Dq[m][q_](t)))*Rm.om(t)
        return tot*Rm.eps(t)
    wts = []
    for (a, b, mons) in pairs:
        if b <= j:
            if a == b: wts.append((mons, {(a, 0): 1.0}))
            else: wts.append((mons, {(a, 0): 0.5, (b, 0): 0.5}))
        else: wts.append((mons, _expand(a, b, 0)))
    def q(k, t, part=None):
        tot = arb(0)
        if part in (None, 'pairs'):
            for mons, w in wts:
                for (kk, d), v in w.items():
                    if kk == k: tot += arb(v)*cd(mons, d, t)
        if part in (None, 'G'): tot += arb(Ghr[k])*Rm.eps(t)*Rm.e2(t, 0)*Rm.om(t)
        if k == j and part in (None, 'rho'): tot += arb(YCr[j + 1])*Rm.rho(t)
        return tot
    def sup_on(f, lo=tmin, n_=200, tol=0.02, depth=12):
        def ev(a_, b_):
            u = float(abs(f(arb((a_ + b_)/2, (b_ - a_)/2))).upper()); return u if u == u else 1e300
        def pt(a_, b_): return float(abs(f(arb((a_ + b_)/2))).mid())
        cells = [(lo + (1 - lo)*i/n_, lo + (1 - lo)*(i + 1)/n_) for i in range(n_)]
        lb = max(pt(a_, b_) for a_, b_ in cells); best = 0.0; stack = [(a_, b_, 0) for a_, b_ in cells]
        while stack:
            a_, b_, dd = stack.pop(); u = ev(a_, b_)
            if u <= lb*(1 + tol) or dd >= depth: best = max(best, u); continue
            m_ = (a_ + b_)/2; lb = max(lb, pt(a_, b_)); stack += [(a_, m_, dd + 1), (m_, b_, dd + 1)]
        return best
    out = {}
    for k in range(j + 1):
        out[k] = {m: sup_on(lambda t, m=m, k=k: q(k, t)*Rm.y(t)**(2*m - 1)/Rm.W(t)) for m in STIR[k]}
        if DIAG: print('  level', j, 'k', k, {p: {m: '%.2g' % sup_on(lambda t, m=m, k=k, p=p: q(k, t, p)*Rm.y(t)**(2*m - 1)/Rm.W(t), n_=100) for m in STIR[k]} for p in ('pairs', 'G', 'rho')}, flush=True)
    out[j + 1] = {m: 0.0 for m in STIR[j + 1]}
    return out, 0.0

def ramp_extras(j, T, M0r, M1r, gam, Jm, Rm, M2star=3.0, M0z=0.02, kap=None):
    """Space-time L^4 treatment of the extra jet factors on the ramp (Theta = eps*omega*eta^2, eta = S(t)).
    Jm = {m: J_m^eff, m=0..4}.  JTh_k = int int Theta z_k^2 <= B_k sum_m S(k,m) sup(Theta y^{2m-1}/W) J_m  (free: Theta -> eps*omega).
    L^4 identity (exact): Y = int int Theta z_2^4 = -3 int int Theta z_1 z_2^2 z_3 - int int Theta' z_1 z_2^3, Theta' = -kap*Theta + 2 eps*omega*eta*eta'.
      int int eps*omega*eta*eta'|z_2|^3 <= (int int eps*omega z_2^2)^{1/2} (eta'_max) sqrt(Y)   [homogeneous: eta'^2 eta^2 <= eta'_max^2 eta^2]
    => sqrt(Y) <= M1 [3 sqrt(JTh_3) + kap sqrt(JTh_2) + 2 eta'_max sqrt(Jfree_2)]."""
    kap = kap or gam + 4
    th = lambda t: Rm.eps(t)*Rm.om(t)*Rm.e2(t, 0); fr = lambda t: Rm.eps(t)*Rm.om(t)
    def Jint(wfun, k):
        return BELL[k]*sum(STIR[k][m]*F.sup_grid(lambda t, m=m: wfun(t)*Rm.y(t)**(2*m - 1)/Rm.W(t), 400)*Jm[m] for m in STIR[k])
    JTh = {k: Jint(th, k) for k in range(5)}; Jfree2 = Jint(fr, 2)
    etap = 1.875/Rm.L                                   # max_t S'(t) / L  (d eta/ds)
    sqY = M1r*(3*math.sqrt(JTh[3]) + kap*math.sqrt(JTh[2]) + 2*etap*math.sqrt(Jfree2))
    items = classify_ramp(T, j, M0r, M1r, gam, M0z)
    tot = 0.0; det = {}
    for it in items:
        if it[0] == 'P': v = it[3]*math.sqrt(JTh[it[1]]*JTh[it[2]])
        elif it[0] == 'L4': v = it[2]*sqY*math.sqrt(JTh[it[1]])
        else: v = it[1]*M2star*JTh[3]
        tot += v; det[it[0]] = det.get(it[0], 0) + v
    return tot, det, JTh, sqY
def classify_ramp(T, j, M0r, M1r, gam, M0z):
    import lemmaS_L4terms as L4
    return L4.classify(T, j, M0r, M1r, gam, M0z)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--J', nargs=3, type=float, default=[1.17e-11, 1.66e-11, 3.7e-12])
    ap.add_argument('--handoff', nargs=3, type=float, default=[2e-7, 1e-4, 0.2])
    ap.add_argument('--sqLam', type=float, default=1.45e-7)
    ap.add_argument('--alpha1', type=float, default=5e-4)
    ap.add_argument('--core', nargs=2, type=float, default=[0.11, 0.75])
    ap.add_argument('--ramp', nargs=2, type=float, default=[0.05, 2.5])
    ap.add_argument('--gam', type=float, default=1.0)
    ap.add_argument('--mh', type=float, default=0.0)
    ap.add_argument('--mhr', type=float, default=None)
    ap.add_argument('--slopes', nargs=3, type=float, default=[0, 0, 0])
    ap.add_argument('--rho', type=float, default=2.6e-3)
    ap.add_argument('--M1s', nargs='+', type=float, default=None)
    ap.add_argument('--J4', type=float, default=None)
    ap.add_argument('--sJ4', type=float, default=0.0)
    ap.add_argument('--c1W', type=float, default=0.8477)
    ap.add_argument('--rampl4', action='store_true')
    ap.add_argument('--s33', default='assumed'); ap.add_argument('--rs33', default='ibp'); ap.add_argument('--m2r', type=float, default=3.0)
    ap.add_argument('--l4', action='store_true')
    ap.add_argument('--tmin', type=float, default=0.0)
    ap.add_argument('--tag', default='F2')
    ap.add_argument('--abar', type=float, default=None); ap.add_argument('--beta', type=float, default=None)
    a = ap.parse_args(); gam = a.gam
    T = F.build_tables()
    M1s = a.M1s or [a.core[1]]
    br = R.box(*(3*[a.mh if a.mhr is None else a.mhr]), M0=a.ramp[0], M1=a.ramp[1])
    cr, Gr = F.constants(gam, br, 'quasi', tables=T)
    import lemmaS_prop1_bounds as PB
    Abr = PB.A_bounds_s(br, rmax=6)[0]
    Rm = F.Ramp(14.0, 28.0, gam, True); rampd = {}; L0s = {}
    footw = {}; W33 = {}
    if a.rampl4:
        import lemmaS_L4terms as L4
        amin_r = (1 + a.ramp[0]**2)**-2; DEL = 0.1; etap = 1.875/Rm.L
        _L = arb(Rm.L)
        def _ce(t):
            S = t**3*(10 - 15*t + 6*t**2); S1 = 30*t**2 - 60*t**3 + 30*t**4; S2 = 60*t - 180*t**2 + 120*t**3
            return -5*S*S1/_L + S1**2/(_L*_L) + S*S2/(_L*_L)
        CETA = float(F.sup_grid(_ce, 4000))   # arb enclosure over 4000 cells (rigorous upper bound)
        print('C_eta =', CETA)
        for j in (0, 1, 2):
            add, log, dis, spare = L4.l4_extra_YC(T, j, a.ramp[0], a.ramp[1], gam, M2star=a.m2r, s33=a.rs33, dfrac=0.08)
            items = L4.classify(T, j, a.ramp[0], a.ramp[1], gam)
            fw = 0.0
            w33 = sum(it[1] for it in items if it[0] == 'S33')
            if a.rs33 in ('ibp', 'ibp2') and w33 > 0:
                c1 = w33*a.ramp[1]*etap; kap = 5.0
                add[3] = add.get(3, 0) + c1*kap*DEL; fw += c1*kap/DEL
                dd = 0.1*amin_r; add[2] = add.get(2, 0); dis += dd; fw += c1**2/(4*dd)
                if a.rs33 == 'ibp': fw += w33*CETA/3*a.m2r
            W33[j] = w33 if a.rs33 == 'ibp2' else 0.0
            assert dis <= spare*(1+1e-9), ('ramp dissipation budget exceeded', j, dis, spare)
            for it in items:
                if it[0] != 'L4': continue
                w = it[2]
                if it[1] == 3: add[3] = add.get(3, 0) + w*a.ramp[1]*etap*DEL; fw += w*a.ramp[1]*etap/DEL
                elif it[1] == 2: add[2] = add.get(2, 0) + w*a.ramp[1]*etap; fw += w*a.ramp[1]*etap
                else: raise ValueError(('foot r', it))
            YC = {int(k): float(v) for k, v in Gr[j]['YC'].items()}; YC = {k: YC[k] + add.get(k, 0) for k in YC}
            Gr[j]['YC'] = YC; Gr[j]['Ghat'] = L4.ghat_from_YC(YC, j, amin_r); footw[j] = fw
            cr[f'Ghat_{j}'] = {str(k): v for k, v in Gr[j]['Ghat'].items()}
    for j in (0, 1, 2): rampd[j], L0s[j] = ramp_level(j, T, Abr, Gr[j]['Ghat'], Gr[j]['YC'], Rm, gam, tmin=a.tmin)
    for M1c in M1s:
        a.core = [a.core[0], M1c]
        bc = R.box(a.mh, a.mh, a.mh, M0=a.core[0], M1=M1c)
        cc, G = F.constants(gam, bc, 'quasi', tables=T)
        if a.l4:
            import lemmaS_L4terms as L4
            amin_c = (1 + a.core[0]**2)**-2
            for j in (0, 1, 2):
                add, log, dis, spare = L4.l4_extra_YC(T, j, a.core[0], M1c, gam, s33=a.s33)
                assert dis <= spare*(1+1e-9), ('dissipation budget exceeded', j, dis, spare)
                YC = {int(k): float(v) for k, v in G[j]['YC'].items()}; YC = {k: YC[k] + add.get(k, 0) for k in YC}
                Gn = L4.ghat_from_YC(YC, j, amin_c); G[j]['YC'] = YC; G[j]['Ghat'] = Gn
                cc[f'Ghat_{j}'] = {str(k): v for k, v in Gn.items()}
        for j in (0, 1, 2): G[j]['ramp'] = rampd[j]
        r = F.assemble(G, tuple(a.J) + (0.0,), a.sqLam, gam, a.alpha1, handoff=tuple(a.handoff), slopes=tuple(a.slopes) + (0.0,), abar_c=a.abar, beta_c=a.beta)
        extra = (0.0, 0.0, 0.0); exinfo = {}
        if a.rampl4:
            ex = []
            for j in (0, 1, 2):
                gpj = float(r['parts'][j]['gamma_prime'])
                Jm = {m: a.J[m] + a.slopes[m]/gpj for m in range(3)}
                fr = lambda t: Rm.eps(t)*Rm.om(t)
                Jf2 = BELL[2]*sum(STIR[2][m]*F.sup_grid(lambda t, m=m: fr(t)*Rm.y(t)**(2*m - 1)/Rm.W(t), 400)*Jm[m] for m in STIR[2])
                v = footw[j]*Jf2; exinfo[j] = dict(foot=v, Jfree2=Jf2, footw=footw[j])
                if W33.get(j, 0) > 0:
                    import math as _m
                    fr = lambda t: Rm.eps(t)*Rm.om(t)
                    cv = lambda k, Jd: BELL[k]*sum(STIR[k][m]*float(F.sup_grid(lambda t, m=m: fr(t)*Rm.y(t)**(2*m - 1)/Rm.W(t), 400))*Jd[m] for m in STIR[k])
                    J4e = a.J4 + a.sJ4/gpj; J3e = _m.sqrt(Jm[2]*J4e) + a.c1W*Jm[2]
                    Jd3 = dict(Jm); Jd3[3] = J3e
                    STIR.setdefault(3, {1: 1, 2: 3, 3: 1}); BELL.setdefault(3, 5)
                    JG = 2*(cv(3, Jd3) + 6.25*cv(2, Jd3))
                    L2 = a.sqLam**2; V0s = L2/0.08; V2s = L2/(0.08*3e-5); V1s = _m.sqrt(V0s*V2s) + a.c1W*V0s
                    Fmax = cv(2, {0: V0s, 1: V1s, 2: V2s})
                    frmin = float(fr(arb(1)))
                    Kf = W33[j]*CETA/3/_m.sqrt(frmin); ell = float(Rm.L)
                    D = Kf*_m.sqrt(Fmax)*(Jf2/_m.sqrt(ell) + _m.sqrt(2)*Jf2**0.75*JG**0.25)
                    exinfo[j].update(S33data=D, Kf=Kf, Fmax=Fmax, JG=JG, J3e=J3e); v += D
                    print('    ramp S33 sup-free data: Kf %.3g Fmax %.3g JG %.3g J3 %.3g -> %.3g' % (Kf, Fmax, JG, J3e, D), flush=True)
                ex.append(v)
                print('    ramp-L4 level', j, 'foot data term = %.3g (Jfree2 %.3g, weight %.3g)' % (v, Jf2, footw[j]), flush=True)
            extra = tuple(ex)
            r = F.assemble(G, tuple(a.J) + (0.0,), a.sqLam, gam, a.alpha1, handoff=tuple(a.handoff), slopes=tuple(a.slopes) + (0.0,), abar_c=a.abar, beta_c=a.beta, extra_outer=extra)
            r['extra_outer'] = exinfo
        gp = [float(p['gamma_prime']) for p in r['parts']]; valid = all(g > 0 for g in gp)
        print('M1core', M1c, 'Ghat core', [round(float(cc[f'Ghat_{j}'][str(j)])) for j in range(3)], ' Ghat ramp', [round(float(cr[f'Ghat_{j}'][str(j)])) for j in range(3)])
        print("  gamma'", [round(g, 4) for g in gp], 'E', [f'{e:.3g}' for e in r['E']], 'S0_v %.4f S1_v %.4f' % (float(r['S0']), float(r['S1'])), '' if valid else '** INVALID **', flush=True)
        for j in range(3): print('    level', j, {k: '%.3g' % float(v) for k, v in r['parts'][j].items()})
        json.dump(dict(args=dict(vars(a)), Ghat_core=cc, Ghat_ramp=cr, L0=L0s, result=r, valid=valid), open(f'lemmaS_{a.tag}_M{M1c}.json', 'w'), indent=1, default=str)
if __name__ == '__main__': main()
