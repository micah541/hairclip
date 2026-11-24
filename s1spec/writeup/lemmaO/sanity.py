#!/usr/bin/env python3
"""Numerical sanity check for Lemma O.
 (1) sympy: verify graph equation and the expansions of G_x, G_xx used in the proof;
     verify numerically the sup bounds of the coefficient functions on |p|<=1/4.
 (2) solve h_t = -(a(h_x)h_xxx + b(h_x)h_xx^2)_x on the periodic interval [-4,4) (Fourier collocation,
     stiff BDF), for smooth data; measure K2,K3,M on the box [-2,2] (hypotheses (A),(I)), compute the
     constants of lemmaO_consts.py and compare with the observed deviations on [-1,1], and the
     energies E1, E3 with the cutoff zeta_E of the proof.
"""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp
from lemmaO_consts import *

def check_algebra():
    import sympy as sp
    x = sp.symbols('x'); h = sp.Function('h')(x); p = h.diff(x); s = 1+p**2
    kap = h.diff(x,2)/s**sp.Rational(3,2); ds = lambda f: f.diff(x)/sp.sqrt(s)
    ht = sp.sqrt(s)*(-ds(ds(kap)))
    a = lambda q: (1+q**2)**-2; b = lambda q: -3*q*(1+q**2)**-3
    G = a(p)*h.diff(x,3)+b(p)*h.diff(x,2)**2
    assert sp.simplify(ht+G.diff(x)) == 0
    q = sp.symbols('q'); A = a(q); Bq = b(q)
    a1, a2, b1, b2 = A.diff(q), A.diff(q,2), Bq.diff(q), Bq.diff(q,2)
    c1, c2, c3, c4 = 2*a1+2*Bq, a1+2*Bq, a2+5*b1, b2
    P1,P2,P3,P4 = [h.diff(x,k) for k in (2,3,4,5)]
    sub = lambda e: e.subs(q,p)
    assert sp.simplify(sub(A)*P3+sub(c2)*P1*P2+sub(b1)*P1**3 - G.diff(x)) == 0
    assert sp.simplify(sub(A)*P4+sub(c1)*P1*P3+sub(c2)*P2**2+sub(c3)*P1**2*P2+sub(c4)*P1**4 - G.diff(x,2)) == 0
    pp = np.linspace(-.25, .25, 200001)
    f = lambda e: sp.lambdify(q, e, 'numpy')(pp)
    sup = lambda e: np.max(np.abs(f(e)))
    got = dict(a=(sup(A), A_MAX), amin=(np.min(f(A)), A_MIN), a1=(sup(a1),A1), b=(sup(Bq),B0), b1=(sup(b1),B1),
               c1=(sup(c1),CC1), c2=(sup(c2),CC2), c3=(sup(c3),CC3), c4=(sup(c4),CC4))
    for k,(v,bd) in got.items():
        ok = (v >= bd-1e-9) if k=='amin' else (v <= bd+1e-9)
        print(f'  coeff {k:5s}: true sup = {v:.5f}  used bound = {bd:.5f}  {"OK" if ok else "FAIL"}'); assert ok
    print('algebra check passed')

L = 8.0; N = 128
xg = -L/2 + L*np.arange(N)/N
kk = 2*np.pi*np.fft.fftfreq(N, d=L/N)
def D(u, m=1): return np.real(np.fft.ifft((1j*kk)**m*np.fft.fft(u)))
a = lambda p: (1+p**2)**-2; b = lambda p: -3*p*(1+p**2)**-3
def rhs(t, h):
    p = D(h); return -D(a(p)*D(h,3) + b(p)*D(h,2)**2)

def zeta_E(w):
    wE = 1-2*w; r = (2-np.abs(xg))/wE; r = np.clip(r, 0, 1)
    return 10*r**3-15*r**4+6*r**5

def run(name, modes, Tend, slope=0.2):
    rng = np.random.default_rng(1)
    h0 = sum(c*np.cos(2*np.pi*k*xg/L + ph) for k,c,ph in modes)
    h0 *= slope/np.max(np.abs(D(h0)))
    ins = np.abs(xg) <= 2
    dx = L/N
    M = max([np.max(np.abs(D(h0,j)[ins])) for j in range(1,5)] + [math.sqrt(np.sum(D(h0,5)[ins]**2)*dx)])
    ts = np.concatenate([[0], np.geomspace(Tend*1e-4, Tend, 60)])
    sol = solve_ivp(rhs, (0, Tend), h0, method='BDF', t_eval=ts, rtol=1e-11, atol=1e-13)
    H = sol.y.T
    K2 = max(np.max(np.abs(D(h,2)[ins])) for h in H)
    K3 = max(math.sqrt(np.sum(D(h,3)[ins]**2)*dx) for h in H)
    sl = max(np.max(np.abs(D(h)[ins])) for h in H)
    assert sl <= 0.25, sl
    c = constants(K2, K3, M, Tend); w = c['w']; ze = zeta_E(w)
    print(f'\n[{name}] M={M:.3f} K2={K2:.3f} K3={K3:.3f} max|h_x|={sl:.3f}  Tend={Tend:g}  (w={w:.2f})')
    print('  constants:', {k: float(f"{v:.4g}") for k,v in c.items() if k in ('C0','C1','C2','C3','Lambda1')})
    inn = np.abs(xg) <= 1
    worst = np.zeros(6)
    for t, h in zip(ts[1:], H[1:]):
        g = h-h0
        obs = [np.max(np.abs(g[inn])), np.max(np.abs(D(g)[inn])), np.max(np.abs(D(g,2)[inn])),
               np.sum(D(g,3)[inn]**2)*dx, np.sum(ze**2*D(g)**2*ins)*dx, np.sum(ze**4*D(g,3)**2*ins)*dx]
        bnd = [c['C0']*math.sqrt(t), c['C1']*math.sqrt(t), c['C2']*math.sqrt(t), c['C3']*t, c['Lambda1']*t, c['C3']*t]
        worst = np.maximum(worst, np.array(obs)/np.array(bnd))
    lab = ['|g|', '|g_x|', '|g_xx|', 'int g_xxx^2 [-1,1]', 'E1', 'E3']
    for l, wv in zip(lab, worst): print(f'  max over t of observed/bound  {l:20s} = {wv:.3e}')
    assert np.all(worst < 1)
    # also the unweighted observed rate: g ~ t ?  g_xx ~ t^(1/2)? print fitted exponents
    g = [np.max(np.abs(D(h-h0,2)[inn])) for h in H[1:]]
    sl_ = np.polyfit(np.log(ts[1:20]), np.log(g[:19]), 1)[0]; print(f'  observed small-t exponent of sup|g_xx|: {sl_:.2f} (bound exponent 0.5)')
    return worst

if __name__ == '__main__':
    check_algebra()
    run('smooth', [(1,1.0,0.3),(2,0.6,1.1),(3,0.4,2.0)], 1e-2)
    run('rough-ish (k<=8)', [(k, 1.0/k, 0.7*k) for k in range(1,9)], 1e-3)
    run('large amplitude curvature', [(3,1.0,0.0),(5,0.7,1.0),(7,0.5,2.0)], 1e-4)
    rng = np.random.default_rng(5)
    run('H5-marginal (spectrum k^-6, k<=40)', [(k, k**-6.0, 6.28*rng.random()) for k in range(1,41)], 1e-3)
    print('\nALL SANITY CHECKS PASSED')
