"""Weights of Lemma F as 4-jets in x (= jets in y):  T_j = 28 eta(y)^2 Yl(y)^(2j-2) Om(Yl / Yh(tau)).
eta  : S(t), t = ln(y/ya)/ln(yb/ya) (Lemma S ramp: ya, yb = 14, 28), S(u) = u^4(35-84u+70u^2-20u^3).
Yl   : similarity flattening at the far end, Yl = ly * Phi(y/ly), ly = ell - X0; Phi(s) = s (s <= s0), Phi' = 1 - psi(u),
       u = (s-s0)/(1-s0), psi = u^5 (21/4 - 49/8 u + 15/8 u^2): Phi' odd about s = 1  ->  T_j even about the far end (C^5).
Om   : optional outer factor (hand-off energies): Om(s) = exp(-p L(ln s)), L(u) = 0 (u<=0), smooth convex, L = u - m (u >= ln2);
       Yh(tau) moves with the characteristics (Yh' = c_h Yh, c_h <= c).  Om = 1 for the global energies.
Jets: arrays (5, n) of f, f', f'', f''', f''''.  tau-derivative: dlogT/dtau."""
import numpy as np
from math import comb

def jmul(a, b):
    out = np.zeros_like(a)
    for n in range(5):
        for m in range(n + 1): out[n] += comb(n, m) * a[m] * b[n - m]
    return out

def jcomp(F, g):
    """F: (5,n) derivatives of the outer function at g0; g: (5,n) jet. Faa di Bruno to order 4"""
    g1, g2, g3, g4 = g[1], g[2], g[3], g[4]
    out = np.empty_like(g)
    out[0] = F[0]; out[1] = F[1] * g1
    out[2] = F[2] * g1**2 + F[1] * g2
    out[3] = F[3] * g1**3 + 3 * F[2] * g1 * g2 + F[1] * g3
    out[4] = F[4] * g1**4 + 6 * F[3] * g1**2 * g2 + F[2] * (3 * g2**2 + 4 * g1 * g3) + F[1] * g4
    return out

def S_der(u):
    """S and derivatives 1..4 at u (clipped: 0 below, 1 above)"""
    u = np.asarray(u, float); uc = np.clip(u, 0, 1)
    P = np.array([0, 0, 0, 0, 35, -84, 70, -20], float)
    out = np.zeros((5,) + u.shape)
    cp = P.copy()
    for m in range(5):
        out[m] = np.polyval(cp[::-1], uc)
        cp = np.array([i * cp[i] for i in range(1, len(cp))]) if len(cp) > 1 else np.zeros(1)
    out[:, u <= 0] = 0.0
    out[1:, u >= 1] = 0.0; out[0, u >= 1] = 1.0
    return out

def eta2_jet(y, ya=14.0, yb=28.0, qr=0.0):
    """eta = S(t); qr = 0: t = ln(y/ya)/ln(yb/ya) (Lemma S ramp); qr > 0: t = (1-(ya/y)^qr)/(1-(ya/yb)^qr) (inflow pushed to the foot)"""
    t = np.zeros((5, len(y)))
    if qr == 0:
        L = np.log(yb / ya)
        t[0] = np.log(y / ya) / L; t[1] = 1 / (y * L); t[2] = -1 / (y**2 * L); t[3] = 2 / (y**3 * L); t[4] = -6 / (y**4 * L)
    else:
        Bq = ya**qr / (1 - (ya / yb)**qr)
        t[0] = (1 - (ya / y)**qr) / (1 - (ya / yb)**qr)
        f = qr
        for m in range(1, 5):
            t[m] = Bq * f * (-1) ** (m + 1) * y ** (-qr - m); f *= (qr + m)
    e = jcomp(S_der(t[0]), t)
    return jmul(e, e), e

PSI = (21 / 4, -49 / 8, 15 / 8)      # psi(u) = u^5 (a + b u + c u^2)
def psi_der(u):
    a, b, c = PSI; u = np.clip(u, 0, 1)
    p = [a * u**5 + b * u**6 + c * u**7, 5 * a * u**4 + 6 * b * u**5 + 7 * c * u**6,
         20 * a * u**3 + 30 * b * u**4 + 42 * c * u**5, 60 * a * u**2 + 120 * b * u**3 + 210 * c * u**4,
         120 * a * u + 360 * b * u**2 + 840 * c * u**3]
    return np.array(p)
def Psi_int(u):
    a, b, c = PSI; u = np.clip(u, 0, 1); return a * u**6 / 6 + b * u**7 / 7 + c * u**8 / 8

def Yl_jet(y, ly, s0=0.5):
    """Yl and x-derivatives; also (Phi - s Phi') for the tau-derivative"""
    s = y / ly; u = (s - s0) / (1 - s0); m = s > s0
    Ph = np.where(m, s0 + (1 - s0) * (np.clip(u, 0, 1) - Psi_int(u)), s)
    ps = psi_der(u)
    d1 = np.where(m, 1 - ps[0], 1.0); d2 = np.where(m, -ps[1] / (1 - s0), 0.0)
    d3 = np.where(m, -ps[2] / (1 - s0)**2, 0.0); d4 = np.where(m, -ps[3] / (1 - s0)**3, 0.0)
    J = np.array([ly * Ph, d1, d2 / ly, d3 / ly**2, d4 / ly**3])
    return J, Ph - s * d1

def power_jet(J, p):
    F = np.array([J[0]**p, p * J[0]**(p - 1), p * (p - 1) * J[0]**(p - 2), p * (p - 1) * (p - 2) * J[0]**(p - 3),
                  p * (p - 1) * (p - 2) * (p - 3) * J[0]**(p - 4)])
    return jcomp(F, J)

def Lsm(u, m=0.25):
    """smooth convex L: 0 for u <= 0, u - m for u >= 2m (C^4), via integral of S"""
    t = u / (2 * m); Sd = S_der(t)
    # L' = S(t)  -> L = 2m * int_0^t S
    from numpy.polynomial import polynomial as Pn
    Pint = Pn.polyint([0, 0, 0, 0, 35, -84, 70, -20])
    tc = np.clip(t, 0, 1)
    L = np.where(t <= 0, 0.0, np.where(t >= 1, u - m, 2 * m * Pn.polyval(tc, Pint)))
    out = np.array([L, Sd[0], Sd[1] / (2 * m), Sd[2] / (2 * m)**2, Sd[3] / (2 * m)**3])
    return out

def weight(j, y, ly, ellt, ramp=(14.0, 28.0), s0=0.5, om=None):
    """ramp = (ya, yb) or (ya, yb, qr)"""
    """jet of T_j (5, n) and dT/dtau (n,).  om = (p, Yh, Yh_tau) for the outer factor, or None."""
    e2, _ = eta2_jet(y, *ramp)
    Jl, ph = Yl_jet(y, ly, s0)
    T = 28.0 * jmul(e2, power_jet(Jl, 2 * j - 2))
    dlnYl = ellt * ph / Jl[0]                       # d ln Yl / dtau  (ly' = ell_tau)
    dlogT = (2 * j - 2) * dlnYl
    if om is not None:
        p, Yh, Yht = om
        lnS = np.zeros_like(Jl); lnS[0] = np.log(Jl[0] / Yh)
        # jet of ln(Yl/Yh) in x
        inv = np.array([1 / Jl[0], -1 / Jl[0]**2, 2 / Jl[0]**3, -6 / Jl[0]**4, 24 / Jl[0]**5])
        lnS = jcomp(np.array([np.log(Jl[0] / Yh), inv[0], inv[1], inv[2], inv[3]]), Jl)
        Lj = jcomp(Lsm(lnS[0]), lnS)                    # L(ln s) jet
        Om = jcomp(np.array([np.exp(-p * Lj[0])] * 1 + [None] * 0).repeat(1, 0) if False else
                   np.array([np.exp(-p * Lj[0]), -p * np.exp(-p * Lj[0]), p**2 * np.exp(-p * Lj[0]),
                             -p**3 * np.exp(-p * Lj[0]), p**4 * np.exp(-p * Lj[0])]), Lj)
        T = jmul(T, Om)
        Lp = Lsm(lnS[0])[1]                             # L'(ln s)
        dlogT += -p * Lp * (dlnYl - Yht / Yh)
    return T, T[0] * dlogT

if __name__ == '__main__':
    # self-test: jets against finite differences
    y = np.linspace(14.5, 900.0, 4001); ly = 1000.0
    for j in (0, 2, 4):
        for om in (None, (4.0, 150.0, 37.0)):
            T, dT = weight(j, y, ly, 250.0, om=om)
            h = y[1] - y[0]
            fd = [np.gradient(T[m], h) for m in range(4)]
            err = max(np.max(np.abs(fd[m][5:-5] - T[m + 1][5:-5]) / (np.abs(T[m + 1][5:-5]).max() + 1e-300)) for m in range(4))
            # tau-derivative by finite difference in ly (ell) and Yh
            dt = 1e-6
            om2 = None if om is None else (om[0], om[1] + om[2] * dt, om[2])
            T2, _ = weight(j, y, ly + 250.0 * dt, 250.0, om=om2)
            errt = np.max(np.abs((T2[0] - T[0]) / dt - dT)) / np.abs(dT).max() if np.abs(dT).max() > 0 else 0
            print(f'j={j} om={om is not None}: jet rel err {err:.1e}, tau-derivative rel err {errt:.1e}')
