#!/usr/bin/env python3
"""Explicit constants for Lemma O (short-time stability of curve diffusion flow, graph form).
Every formula here is the one proved in lemmaO.tex (section/eq numbers in comments).
Usage: python3 lemmaO_consts.py K2 K3 M [t1]
"""
import sys, math
import numpy as np

# ---- universal constants (|h_x|<=1/4) -------------------------------------------------
r = (16/17)**3                    # sup_{|p|<=1/4} (1+p^2)^-3 * ... : p(1+p^2)^-3 increasing on [0,1/sqrt5]
A_MAX = 1.0                       # sup a
A_MIN = (17/16)**-2               # inf a = 0.8858
A1  = r                           # sup|a'|  = 4|p|(1+p^2)^-3 <= 4*(1/4)*r
B0  = 0.75*r                      # sup|b|   = 3|p|(1+p^2)^-3
B1  = 3.0                         # sup|b'|  = |-3+21p^2... | <= 3  (b' = -3(1-5p^2)(1+p^2)^-4... see tex)
CC1 = 3.5*r                       # sup|c1| = 14|p|(1+p^2)^-3
CC2 = 2.5*r                       # sup|c2| = 10|p|(1+p^2)^-3   (also = |a'+2b|)
CC3 = 19.0                        # sup|c3| = 19(1-5p^2)(1+p^2)^-4
CC4 = 13.5                        # sup|c4| = 18|p|(3-5p^2)(1+p^2)^-5 <= 54/4
# quintic smoothstep S(s)=10s^3-15s^4+6s^5 : max S'=15/8, max |S''|=10/sqrt3, int_0^1 S'^2 = 10/7
S1, S2, S1SQ = 15/8, 10/math.sqrt(3), 10/7

def b_prime_sup():
    p = np.linspace(0, .25, 100001)
    return np.max(np.abs(-3*(1-5*p**2)/(1+p**2)**4))   # sympy-checked closed form of b'

def E1_const(K2, kappa, M, Z, Z2, mE):
    """Lambda1 with E1(t) = int zeta^2 g_x^2 <= Lambda1 * t   (tex Prop. 3).
    kappa >= sup_t ||g_xxx||_{L2(supp zeta)};  mE = |{zeta' != 0}|."""
    p1 = 0.25 + M; k2 = K2 + M; sm = math.sqrt(mE)
    rho0 = M*(1 + CC2*K2) + B1*K2**3
    beta = (4*Z*k2*sm + 2*A1*K2*k2 + 2*Z*A1*K2*p1*sm + 2*k2*CC2*K2 + 2*Z*p1*CC2*K2*sm)
    gamma = 2*(Z**2 + Z2)*p1*sm*kappa + 4*k2*rho0 + 2*Z*p1*rho0*mE
    return 2*gamma + beta**2/(2*A_MIN)

def E3_const(K2, kappa, M, Z, q):
    """Cstar with E3(t) = int zeta^4 g_xxx^2 <= Cstar * t   (tex Prop. 4).
    kappa >= sup_t ||g_xxx||_{L2(supp zeta)}; Z>=|zeta'|, q = 12 Z^2 + 4 sup|zeta''|."""
    K = kappa
    # polynomial in u = s^(1/4):   s = u^4 ,  s = ||zeta^2 g_5||
    P = np.polynomial.polynomial
    Yb = np.array([2*Z*K, 0, K**0.5])                      # 2ZK + K^(1/2) u^2
    om2 = np.array([K**1.5*math.sqrt(3*Z*K), K**1.5*K**0.25])   # <= K^(3/2)[(3ZK)^(1/2)+K^(1/4)u]
    nu0 = (M + 2*CC1*K2*M + CC2*4*M**2 + CC3*K2**2*(2*M + K) + 2*CC4*K2**4)
    nu = P.polyadd(P.polyadd([nu0], CC1*K2*Yb), 2*CC2*om2)
    c8 = P.polymul([0,0,0,0,1], [1])                        # u^4 = s
    lin = P.polyadd(8*Z*Yb, [q*K])                          # 8 Z Ybar + q K
    f = P.polyadd(-A_MIN*np.array([0]*8+[1.0]), P.polymul(lin, [0,0,0,0,1.0]))
    f = P.polyadd(f, P.polymul(P.polyadd([0,0,0,0,1.0], lin), nu))
    df = P.polyder(f)
    rts = P.polyroots(df)
    cand = [0.0] + [z.real for z in rts if abs(z.imag) < 1e-9*max(1,abs(z)) and z.real > 0]
    fmax = max(P.polyval(c, f) for c in cand)
    uu = np.geomspace(1e-6, 1e6, 4001)
    assert fmax >= P.polyval(uu, f).max()*(1-1e-9) - 1e-9   # root-finding sanity
    return 2*fmax*(1+1e-9)

def level(K2, K3, M, w):
    """Single energy cutoff zeta_E: =1 on |x|<=1+2w, =0 for |x|>=2 (quintic smoothstep, width 1-2w)."""
    wE = 1 - 2*w
    Z = S1/wE; Z2 = S2/wE**2; q = 12*Z**2 + 4*Z2
    kap = K3 + 2*M                        # ||g_xxx||_{L2(-2,2)} <= K3 + ||h0_xxx||_{L2(-2,2)} <= K3 + 2M
    return E1_const(K2, kap, M, Z, Z2, 2*wE), E3_const(K2, kap, M, Z, q), kap

def constants(K2, K3, M, t1=1.0, w=None):
    """Return dict with C0..C3 (valid for all 0<=t<=t1) and intermediates."""
    def one(w):
        L1, Cs, kap = level(K2, K3, M, w)
        Zq = Zp = S1/w
        Yc = 2*Zq*math.sqrt(L1) + (L1*Cs)**0.25                    # ||zeta_q g_xx|| <= Yc sqrt(t)
        C1 = math.sqrt(math.sqrt(L1)*(Yc + Zq*math.sqrt(L1)))
        C2 = math.sqrt(Yc*(math.sqrt(Cs) + Zp*Yc))
        Gam = K3 + 2*B0*K2**2
        thp = math.sqrt(2*S1SQ/w)                                    # ||zeta_p'||_2
        C0 = math.sqrt(2 + w)*math.sqrt(L1) + thp*Gam*math.sqrt(t1)/2
        return dict(C0=C0, C1=C1, C2=C2, C3=Cs, Lambda1=L1, w=w, kappa=kap, Yc=Yc, Gamma=Gam)
    if w is not None:
        return one(w)
    best = None
    for ww in np.linspace(0.04, 0.46, 22):
        c = one(ww)
        sc = max(c['C0'], c['C1'], c['C2'], math.sqrt(c['C3']))
        if best is None or sc < best[0]:
            best = (sc, c)
    return best[1]

def usable_time(K2, K3, M, eps=0.1):
    """Largest t such that all four deviation bounds are <= eps (C0,C1,C2 sqrt t <= eps, sqrt(C3 t) <= eps)."""
    lo, hi = 1e-30, 1.0
    for _ in range(70):
        t = math.sqrt(lo*hi); c = constants(K2, K3, M, t)
        ok = max(c['C0'], c['C1'], c['C2'], math.sqrt(c['C3']))*math.sqrt(t) <= eps
        lo, hi = (t, hi) if ok else (lo, t)
    return lo

def delta0(K2, K3, M):
    """The proof needs NO smallness of t1 (the constants are valid for every t1; only C0 depends on t1
    through the term thp*Gamma*sqrt(t1)/2).  So delta0 = infinity formally.  See usable_time() for the
    time below which the bounds say something."""
    return float('inf')

if __name__ == '__main__':
    K2, K3, M = map(float, sys.argv[1:4]); t1 = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    assert abs(b_prime_sup() - B1) < 1e-12 or b_prime_sup() <= B1
    c = constants(K2, K3, M, t1)
    print({k: float(f'{v:.6g}') for k, v in c.items()})
    for k in ('C0', 'C1', 'C2'):
        print(f'{k}*sqrt(t1) = {c[k]*math.sqrt(t1):.4g}')
    print(f"sqrt(C3*t1) = {math.sqrt(c['C3']*t1):.4g}")
    print('usable_time(eps=0.1) =', usable_time(K2, K3, M, 0.1), '  trivial a-priori sizes: |g_x|<=%.3g |g_xx|<=%.3g' % (0.25+M, K2+M))
