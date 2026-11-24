#!/usr/bin/env python3
"""Lemma O, version 2: cylinder [-X,X] (eta=1 on [-1-2w,1+2w], transition to X), small-data regime,
additive slack thresholds t*.  Same proof as v1 with |I0|=2X; formulas are those of lemmaO_v2.tex.
Usage: python3 lemmaO_consts_v2.py K2 K3 M X [slack_factor=0.5]
"""
import sys, math
import numpy as np
from lemmaO_consts import (A_MAX, A_MIN, A1, B0, B1, CC1, CC2, CC3, CC4, S1, S2, S1SQ)

def E1_const(K2, kappa, M, X, Z, Z2, mE, drop=False):
    lam = math.sqrt(2*X)
    p1 = 0.25 + M; k2 = K2 + M; sm = math.sqrt(mE)
    rho0 = M*(1 + CC2*K2) + B1*K2**3
    if drop: rho0 = 0.0
    beta = (4*Z*k2*sm + lam*A1*K2*k2 + 2*Z*A1*K2*p1*sm + lam*k2*CC2*K2 + 2*Z*p1*CC2*K2*sm)
    gamma = 2*(Z**2 + Z2)*p1*sm*kappa + 2*X*k2*rho0 + 2*Z*p1*rho0*mE
    return 2*gamma + beta**2/(2*A_MIN)

def E3_const(K2, kappa, M, X, Z, q, drop=False):
    lam = math.sqrt(2*X); K = kappa
    P = np.polynomial.polynomial
    Yb = np.array([2*Z*K, 0, K**0.5])
    om2 = np.array([K**1.5*math.sqrt(3*Z*K), K**1.5*K**0.25])
    nu0 = (M + 2*CC1*K2*lam*M/2*1.0 + 2*CC2*lam*M**2 + CC3*K2**2*(lam*M + K) + lam*CC4*K2**4)
    # note: ||zeta^2 N2|| <= C1 K2 (lam M + Y): the constant part is C1*K2*lam*M
    nu0 = (M + CC1*K2*lam*M + 2*CC2*lam*M**2 + CC3*K2**2*(lam*M + K) + lam*CC4*K2**4)
    if drop: nu0 = 0.0
    nu = P.polyadd(P.polyadd([nu0], CC1*K2*Yb), 2*CC2*om2)
    lin = P.polyadd(8*Z*Yb, [q*K])
    f = P.polyadd(-A_MIN*np.array([0]*8+[1.0]), P.polymul(lin, [0,0,0,0,1.0]))
    f = P.polyadd(f, P.polymul(P.polyadd([0,0,0,0,1.0], lin), nu))
    rts = P.polyroots(P.polyder(f))
    cand = [0.0] + [z.real for z in rts if abs(z.imag) < 1e-9*max(1,abs(z)) and z.real > 0]
    fmax = max(P.polyval(c, f) for c in cand)
    uu = np.geomspace(1e-6, 1e6, 4001)
    assert fmax >= P.polyval(uu, f).max()*(1-1e-9) - 1e-12
    return 2*fmax*(1+1e-9)

def constants(K2, K3, M, X, t1=1.0, w=None, drop=False):
    lam = math.sqrt(2*X)
    def one(w):
        wE = X - 1 - 2*w
        if wE <= 0.05: return None
        Z = S1/wE; Z2 = S2/wE**2; q = 12*Z**2 + 4*Z2
        kap = K3 + lam*M
        L1 = E1_const(K2, kap, M, X, Z, Z2, 2*wE, drop); Cs = E3_const(K2, kap, M, X, Z, q, drop)
        Zq = Zp = S1/w
        Yc = 2*Zq*math.sqrt(L1) + (L1*Cs)**0.25
        ell = 2 + 2*w          # J = [-1-w, 1+w]; sup_J u^2 <= ||u||_J^2/ell + 2||u||_J ||u'||_J   (Lemma 3 (S'))
        C1 = math.sqrt(min(math.sqrt(L1)*(Yc + Zq*math.sqrt(L1)), L1/ell + 2*math.sqrt(L1)*Yc))
        C2 = math.sqrt(min(Yc*(math.sqrt(Cs) + Zp*Yc), Yc**2/ell + 2*Yc*math.sqrt(Cs)))
        Gam = K3 + lam*B0*K2**2
        C0 = math.sqrt(2 + w)*math.sqrt(L1) + math.sqrt(2*S1SQ/w)*Gam*math.sqrt(t1)/2
        return dict(C0=C0, C1=C1, C2=C2, C3=Cs, Lambda1=L1, w=w, wE=wE, X=X, kappa=kap, Yc=Yc, Gamma=Gam)
    if w is not None: return one(w)
    best = None
    for ww in np.linspace(0.03, (X-1.05)/2, 40):
        c = one(ww)
        if c is None: continue
        sc = max(c['C1']/1.0, c['C2'], math.sqrt(c['C3'])/math.sqrt(2))
        if best is None or sc < best[0]: best = (sc, c)
    return best[1]

def t_star(K2, K3, M, X, f=0.5, drop=False):
    """Thresholds where C_j sqrt(t) = slack (C1,C2,C3 do not depend on t1):
    slack_1 = slack_2 = f*M  (pointwise g_x, g_xx on [-1,1]);  L2 of g_xxx over [-1,1]: sqrt(C3 t) = f*M*sqrt(2).
    w is chosen to maximise min(t1,t2,t3)."""
    best = None
    for ww in np.linspace(0.03, (X-1.05)/2, 40):
        c = constants(K2, K3, M, X, w=ww, drop=drop)
        if c is None: continue
        out = dict(t1=(f*M/c['C1'])**2, t2=(f*M/c['C2'])**2, t3=(f*M*math.sqrt(2))**2/c['C3'])
        m = min(out.values())
        if best is None or m > best[0]: best = (m, out, c)
    return best[1], best[2]

if __name__ == '__main__':
    K2, K3, M, X = map(float, sys.argv[1:5]); f = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
    ts, c = t_star(K2, K3, M, X, f)
    print({k: float(f'{v:.4g}') for k, v in c.items()}); print({k: v for k, v in ts.items()})
