"""T5 residual norms of a (polished) slab.
r := V~ + k~_xx - (1/4 - a~) p~ - b~ n~ (normal-velocity defect of the reconstructed base). With rho_kappa = kappa_tau|_x - Phi_x
and rho_theta = theta_tau - Phi = int_0^x rho_kappa (both exact polynomials per element, computed in arb), the reconstruction gives
   r(x) = Re[ e^{i theta(x)} int_0^x e^{-i theta} rho_theta ],   A = Im[...],   r' = rho_theta - k A,  r'' = rho_kappa - k' A - k^2 r,
   r''' = rho_kappa' - k'' A - 3 k k' r - k^2 r'      (r(0) = 0 by the choice of b~, r'(0) = 0 by symmetry),
   S1 = int_0^x r,  S2 = int_0^x S1.  Norms as in T5 (R): R2, Rs on [0,34.25] with rho; R1, RIs with zeta_inf, W_inf on the
whole half curve; R^F_j = (28 int eta^2 x^(2j-2) (d^(j-1) r)^2)^(1/2), d^-1 r = S1, eta = S(ln(y/14)/ln 2), y = x + STIP."""
import numpy as np
from numpy.polynomial import legendre as LG
from flint import arb, arb_poly
from rep import K0, DS1, STIP

NQ = 40
g, w = LG.leggauss(NQ); UQ = (g + 1) / 2; WQ = w / 2
_V = LG.legvander(g, NQ - 1); _Vi = np.linalg.inv(_V)
_I = np.array([LG.legval(g, LG.legint(np.eye(NQ)[j], lbnd=-1) / 2.0) for j in range(NQ)]).T @ _Vi    # cumulative on [0,u]
UQA = [arb(float(v)) for v in UQ]

def smooth(u):
    u = np.clip(u, 0, 1); return u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)

def samples(SP, theta):
    """arb -> float samples at Gauss points of every element at relative time theta"""
    Y, Yd = SP.Y_at(theta)
    els, gl = SP.op.elements(Y, Yd)
    N = len(els); X = np.zeros((N, NQ)); H = np.zeros(N)
    rt, rk, rk1, k0, k1, k2, tt = [np.zeros((N, NQ)) for _ in range(7)]
    acc = arb(0); acct = arb(0)
    for el in els:
        e = el['e']; he = el['h']; ih = 1 / he
        rho_k = el['kt'] - el['Phi'].derivative() * ih                  # kappa_tau|_x - Phi_x
        rho_t = arb_poly([acc]) + rho_k.integral() * he                  # theta_tau - Phi
        th_t = arb_poly([acct]) + el['kt'].integral() * he               # theta_tau (diagnostic)
        d1 = rho_k.derivative() * ih
        kp = el['kp']; ku = el['ku']; kuu = ku.derivative()
        for q in range(NQ):
            uq = UQA[q]
            rt[e, q] = float(rho_t(uq).mid()); rk[e, q] = float(rho_k(uq).mid()); rk1[e, q] = float(d1(uq).mid())
            k0[e, q] = float(kp(uq).mid()); k1[e, q] = float((ku(uq) * ih).mid()); k2[e, q] = float((kuu(uq) * ih * ih).mid())
            tt[e, q] = float(th_t(uq).mid())
        acc = rho_t(arb(1)); acct = th_t(arb(1))
        X[e] = float(el['x'].mid()) + float(he.mid()) * UQ; H[e] = float(he.mid())
    return dict(X=X, H=H, rt=rt, rk=rk, rk1=rk1, k0=k0, k1=k1, k2=k2, tt=tt, b=float((gl['k2t'] + (arb(1) / 4 - gl['a']) * arb(DS1)).mid()), ell=float(gl['ell'].mid()), a=float(gl['a'].mid()))

def cumint(f, H):
    """cumulative integral from x=0 at all Gauss points (element-wise Gauss + exact-polynomial cumulative matrix)"""
    loc = H[:, None] * (f @ _I.T)
    tot = np.r_[0.0, np.cumsum(H * (f @ WQ))] if np.isrealobj(f) else np.r_[0.0, np.cumsum(H * (f @ WQ))]
    return tot[:-1, None] + loc

def norms(sm):
    X, H = sm['X'], sm['H']; k0 = sm['k0']
    th = -np.pi / 2 + cumint(k0, H)
    z = np.exp(1j * th) * cumint(np.exp(-1j * th) * sm['rt'], H)
    r, Aim = z.real, z.imag
    r1 = sm['rt'] - k0 * Aim
    r2 = sm['rk'] - sm['k1'] * Aim - k0**2 * r
    r3 = sm['rk1'] - sm['k2'] * Aim - 3 * k0 * sm['k1'] * r - k0**2 * r1
    S1 = cumint(r, H); S2 = cumint(S1, H)
    Wq = H[:, None] * WQ[None, :]
    rho = np.exp(-0.25 * (1 + X**2) ** (2 / 3)); mr = X <= 34.25
    Winf = (1 + X**2) ** -2.0; zeta = smooth((X - 10.25) / 6)
    y = X + STIP; eta = np.where(y > 14, smooth(np.log(np.maximum(y, 14) / 14) / np.log(2)), 0.0)
    I = lambda f, m=None: float(np.sqrt(np.sum((Wq * f)[m] if m is not None else Wq * f)))
    out = dict(R2=I(rho * S2**2, mr), Rs=I(rho * r**2, mr), R1=I(Winf * (zeta * S1)**2), RIs=I(Winf * (zeta * r1)**2))
    ders = [S1, r, r1, r2, r3]
    for j in range(5):
        out[f'RF{j}'] = float(np.sqrt(28 * np.sum(Wq * eta**2 * X**(2 * j - 2) * ders[j]**2)))
    # where the weak residual sits, and pointwise sizes in the tip zone
    out['supS2_tip'] = float(np.max(np.abs(S2[X <= 20]))); out['supr_tip'] = float(np.max(np.abs(r[X <= 20])))
    out['supS1_W'] = float(np.max(np.abs(S1[(X >= 10) & (X <= 100)])))
    return out, dict(r=r, r1=r1, r2=r2, r3=r3, S1=S1, S2=S2, X=X, Wq=Wq, eta=eta, rt=sm['rt'], rk=sm['rk'], rk1=sm['rk1'])
