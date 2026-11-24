#!/usr/bin/env python3
"""Lemma R (level-4 companion of Lemma O), certified constants.  python-flint arb, 200 bits, outward rounding at 1e-6.
Formulas: lemmaR.tex.  Everything of Lemma O is taken from lemmaO_v32_interval.py (same cutoffs zeta_E, zeta_q, same polynomial f of
v1 Prop. 3 / v2 Theorem), with two changes that are proved in lemmaR.tex:
  (K)   K := kappa_R = K3 + min(lambda M, K3, n3), n3 >= ||h_0'''||_L2(-X,X)   (any K >= sup_t ||g_xxx||_L2(supp zeta) is admissible in Prop. 3);
  (th)  f_theta(u) := f(u) + theta a_min u^8,  C_theta := 2 max_{u>=0} f_theta(u)   (branch and bound, same routine as C_*).
Constants:
  C_*      = 2 max f                                   E_3(t) <= C_* t
  D5:      int_0^t ||zeta^2 g_5||^2 <= C_theta t / (2 theta a_min)
  R4(t)    = [ sqrt(C_* C_theta t / (4 theta a_min)) + Z_q^2 C_* t / (2 delta) ] / (1 - delta)
             >= (1/t) int_0^t int_{|x|<=1+w} g_xxxx^2 dx dt'        (0 < theta, delta < 1 free)
usage: python3 -B lemmaR_interval.py        (writes lemmaR_consts.json)"""
import json, sys, os
from fractions import Fraction as F
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lemmaO_v32_interval as O
from flint import arb
R, up, down = O.R, O.up, O.down

def check_identity():
    """h_xxxx = s^(5/2) [k_ss + 10 p k k_s + (3 + 15 p^2) k^3],  s = 1 + p^2, p = h_x, k curvature, subscript s = arclength (sympy, exact)."""
    import sympy as sp
    x = sp.Symbol('x'); h = sp.Function('h')(x); p = sp.diff(h, x); s = 1 + p**2
    k = sp.diff(h, x, 2)/s**sp.Rational(3, 2); ds = lambda u: sp.diff(u, x)/sp.sqrt(s)
    rhs = s**sp.Rational(5, 2)*(ds(ds(k)) + 10*p*k*ds(k) + (3 + 15*p**2)*k**3)
    assert sp.simplify(sp.diff(h, x, 4) - rhs) == 0
    rhs3 = s**2*(ds(k) + 3*p*k**2); assert sp.simplify(sp.diff(h, x, 3) - rhs3) == 0
    return True

def build_f(X, M, K2, K3, w, n3):
    """the polynomial f of Prop. 3 (coefficients as in lemmaO_v32_interval.consts) with K = kappa_R; returns (f, kappa_R, Zq)"""
    X, M, K2, K3, w, n3 = map(R, (X, M, K2, K3, w, n3))
    lam = (2*X).sqrt(); wE = X - 1 - 2*w; assert wE.lower() > 0 and w.lower() > 0
    Z = O.S1/wE; Z2 = O.S2/(wE*wE); q = 12*Z*Z + 4*Z2; lm = lam*M
    kap = K3 + O.amin(O.amin(lm, K3), n3)
    K = kap; sK = K.sqrt()
    Yb = [2*Z*K, arb(0), sK]
    om2 = [K*sK*(3*Z*K).sqrt(), K*sK*sK.sqrt()]
    nu0 = M + O.CC1*K2*lm + 2*O.CC2*lam*M*M + O.CC3*K2*K2*(lm + K) + lam*O.CC4*K2**4
    nu = O.padd(O.padd([nu0], O.pscale(Yb, O.CC1*K2)), O.pscale(om2, 2*O.CC2))
    lin = O.padd(O.pscale(Yb, 8*Z), [q*K])
    f = O.padd([arb(0)]*8 + [-O.A_MIN], O.pmul(lin, [arb(0)]*4 + [arb(1)]))
    f = O.padd(f, O.pmul(O.padd([arb(0)]*4 + [arb(1)], lin), nu))
    f = f + [arb(0)]*(9 - len(f))
    return f, kap, O.S1/w

def lemmaR(X, M, K2, K3, w, n3, theta, delta, t):
    f, kap, Zq = build_f(X, M, K2, K3, w, n3)
    Cs = 2*O.fmax_upper(f)
    th, de, tt = R(theta), R(delta), R(t)
    fth = list(f); fth[8] = -(1 - th)*O.A_MIN
    Cth = 2*O.fmax_upper(fth)
    D5 = Cth/(2*th*O.A_MIN)                                   # int_0^t ||zeta^2 g_5||^2 <= D5 * t
    R4 = ((Cs*Cth*tt/(4*th*O.A_MIN)).sqrt() + Zq*Zq*Cs*tt/(2*de))/(1 - de)
    return dict(kappa_R=kap, Cstar=Cs, Ctheta=Cth, D5rate=D5, R4=R4, Zq=Zq)

def best(X, M, K2, K3, n3, t, wgrid, thgrid, degrid):
    b = None
    for w in wgrid:
        for th in thgrid:
            for de in degrid:
                c = lemmaR(X, M, K2, K3, w, n3, th, de, t)
                if b is None or c['R4'].upper() < b[0]['R4'].upper(): b = (c, w, th, de)
    return b

if __name__ == '__main__':
    assert check_identity(); print('frame identities verified (sympy)', flush=True)
    out = dict(note='Lemma R constants; *_upper rigorous (arb, 200 bits, rounded outward 1e-6). R4_upper >= time average over [0,t] of int_{|x|<=1+w} g_xxxx^2 (member units).', cases=[])
    X, M, K = F(3), F(12, 100), F(18, 100); n3 = F(4, 1000)
    wg = [F(i, 100) for i in (6, 8, 10, 12, 16, 20)]; tg = [F(i, 100) for i in (25, 30, 35)]; dg = [F(1, 10), F(15, 100), F(2, 10), F(3, 10)]
    # members used by the outer ramp [190, 380] (lemmaS_v32/ramp_cert.py): r_m = 4*190/3, 4*sqrt(190*285)/3, 380; and the innermost member
    for Lm, tag in ((F(50), 'cone member L=50 (r_m=200)'), (F(190, 3), 'cone member L=63.33 (r_m=253.33)'), (F(7757, 100), 'cone member L=77.57 (r_m=310.28)'), (F(95), 'cone member L=95 (r_m=380)')):
        t = float(F(16001, 10**11)*(F(50)/Lm)**4) * (1 + 1e-12)
        c, w, th, de = best(X, M, K, K, n3, t, wg, tg, dg)
        c0 = lemmaR(X, M, K, K, w, K, th, de, t)            # n3 := K3: the kappa of Lemma O v3 (no data threshold)
        row = dict(name=tag, X='3', M='0.12', K2='0.18', K3='0.18', n3=str(n3), L_m=float(Lm), t_needed=t, w=str(w), theta=str(th), delta=str(de),
                   kappa_R_upper=up(c['kappa_R']), Cstar_upper=up(c['Cstar']), Ctheta_upper=up(c['Ctheta']), D5rate_upper=up(c['D5rate']),
                   R4_upper=up(c['R4']), R4_upper_without_n3=up(c0['R4']), Cstar_without_n3=up(c0['Cstar']))
        out['cases'].append(row); print(row, flush=True)
    json.dump(out, open(os.path.join(HERE, 'lemmaR_consts.json'), 'w'), indent=1)
