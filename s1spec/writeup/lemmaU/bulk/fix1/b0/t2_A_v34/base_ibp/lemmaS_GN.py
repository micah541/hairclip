"""Certified constants + numerical sanity test of the GN inequality for E_I (sup|h'| in terms of sup|h| and ||h''||_L2).
 h = w' (slope perturbation), h' = w'', h'' = w'''.
 (GN1) int eta^2 h'^2 <= K S0^2 + S0 Leta N2,   K = ||(eta^2)''||_L1 / 2,  Leta = ||eta||_L2,  S0=sup_{I'}|h|,  N2=||h''||_L2(I')
 (GN2) sup_I |h'|^2 <= 4 S0^2/|I|^2 + 2 sqrt(int_I h'^2) N2
 =>   sup_I|h'|^2 <= 4 S0^2/|I|^2 + 2 N2 sqrt(K S0^2 + S0 Leta N2)."""
import numpy as np, sympy as sp, json
t = sp.Symbol('t'); S = t**3*(10 - 15*t + 6*t**2)
def eta_consts(delta, ends=2, power=1):
    # eta^2 = S(t)^power on a ramp of width delta (t=(x-a)/delta), ends = number of ramps
    e2 = S**power; d2 = sp.lambdify(t, sp.diff(e2, t, 2)); xs = np.linspace(0, 1, 200001)
    K = 0.5*ends*np.trapezoid(np.abs(d2(xs)), xs)/delta           # int |(eta^2)''| dx = int |d2|/delta^2 * delta dt
    return K
def consts(Ilen, delta, ends=2):
    K = eta_consts(delta, ends); Leta = np.sqrt(Ilen + ends*delta*0.5)   # int eta^2 <= |I| + ends*delta (<=; 0.5 not rigorous -> use 1)
    Leta = np.sqrt(Ilen + ends*delta)
    return K, Leta
def bound(S0, N2, Ilen, delta, ends=2):
    K, Leta = consts(Ilen, delta, ends)
    return np.sqrt(4*S0**2/Ilen**2 + 2*N2*np.sqrt(K*S0**2 + S0*Leta*N2))
if __name__ == '__main__':
    rng = np.random.default_rng(0); worst = 0
    a, b, delta = 34.25, 60.0, 5.0; Ilen = b - a
    x = np.linspace(a - delta, b + delta, 40001)
    for trial in range(300):
        k = rng.integers(1, 40, size=4); ph = rng.uniform(0, 6.3, 4); am = rng.normal(size=4)/k**rng.uniform(0, 2)
        h = sum(am[i]*np.sin(k[i]*x*0.3 + ph[i]) for i in range(4))*10**rng.uniform(-3, 0)
        hp = np.gradient(h, x); hpp = np.gradient(hp, x)
        S0 = np.abs(h).max(); N2 = np.sqrt(np.trapezoid(hpp**2, x)); m = (x >= a) & (x <= b)
        lhs = np.abs(hp[m]).max(); rhs = bound(S0, N2, Ilen, delta)
        worst = max(worst, lhs/rhs)
    K, Leta = consts(Ilen, delta)
    print('K=%.4f Leta=%.3f  worst lhs/rhs over random tests = %.3f (must be <= 1)' % (K, Leta, worst))
    # example: S0 = 0.0153, N2^2 = ||w''''||... illustrative N2 = 1e-3
    for N2 in (1e-4, 1e-3, 1e-2): print('S0=0.0153 N2=%g -> sup|w\'\'| <= %.4g' % (N2, bound(0.0153, N2, Ilen, delta)))
