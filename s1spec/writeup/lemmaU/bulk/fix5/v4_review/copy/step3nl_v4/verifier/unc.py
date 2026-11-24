"""V1: unconditional time-integral machinery.
 Phi(alpha; delta, theta) := sum_k l_k(alpha)^2 / (4 delta n_k)  >=  sup_{s>0, shapes} ( L(alpha, s x_hat) - delta N(s x_hat) ),
 where L <= sum_k l_k x_k is the linear (x-degree 1) part of tot (forcing by alpha), majorised variable by variable (AM-GM for W3, V1, V3 with free
 parameters kappa), N = NEG - theta*NEGV (theta-reduced explicit quadratic, n_k its coefficients).  Pure arb arithmetic (200 bit), no search needed
 for rigor: any kappa>0 gives a valid bound; kappa is chosen by a float coordinate search.
 qbl_sup(): rigorous sup over shapes of (A+W1+B)^2/Lyap  (for the placeholder q_beta term in the alpha-equation)."""
import numpy as np
from flint import arb
import constants as K
A_ = lambda v: arb(str(v))
def Fn_(al): return A_(K.FN[0]) * al + A_(K.FN[1]) * al ** 2 + A_(K.FN[2]) * al ** 3
def FnW_(al): return A_(K.FNW[0]) * al + A_(K.FNW[1]) * al ** 2 + A_(K.FNW[2]) * al ** 3
def nk(theta):
    th = A_(theta); mu = A_(K.MU); eps = A_(K.EPS); eI = A_(K.EPSI)
    return dict(A=arb(1) / 32 - A_(K.PI3), B=arb(1) / 9, D4=eps / 2, V0=mu * A_(K.MU1) * (1 - th), V2=(1 - th) * (mu / 2 + mu * eI * A_(K.MU2)), V4=(1 - th) * mu * eI / 2)
def lk_v4(al, kap):
    """fin13v4 (SF-3): linear majorant of L generated from rhs_iv.L_list, the SAME list that rhs_pn sums (no hand-written copy).
    W3 = sqrt(B (D4 + c1 B)) <= (k3 + c1/k3)/2 B + D4/(2 k3);  V1 <= (k1 + C1W/k1)/2 V0 + V2/(2 k1);  V3 <= (k2 + C1W/k2)/2 V2 + V4/(2 k2)."""
    import rhs_iv as R
    k3, k1, k2 = [A_(x) for x in kap]
    kd = {n: arb(repr(v)) for n, v in K.V4C.items()}
    kd['two'] = arb(2); kd['f572'] = A_('5.72'); kd['t22'] = A_('2.2'); kd['t12'] = A_('1.2')
    mu = A_(K.MU); eI = A_(K.EPSI); eps = A_(K.EPS); c1 = arb(repr(K.c1)); C1W = arb(repr(K.C1W))
    l = {n: arb(0) for n in ('A', 'B', 'D4', 'V0', 'V2', 'V4')}
    for cf, mono in R.L_list(al, kd, mu, mu * eI, eps):
        if mono in l: l[mono] += cf
        elif mono == 'W3': l['B'] += cf * (k3 + c1 / k3) / 2; l['D4'] += cf / (2 * k3)
        elif mono == 'V1': l['V0'] += cf * (k1 + C1W / k1) / 2; l['V2'] += cf / (2 * k1)
        elif mono == 'V3': l['V2'] += cf * (k2 + C1W / k2) / 2; l['V4'] += cf / (2 * k2)
        else: raise KeyError(mono)
    return l
def lk(al, kap):
    """linear majorant coefficients l_k(alpha) of L, alpha an arb (use the upper end of the alpha cell)"""
    if K.NEW: return lk_v4(al, kap)
    k3, k1, k2 = [A_(x) for x in kap]; eps = A_(K.EPS); mu = A_(K.MU); eI = A_(K.EPSI); c1 = A_(K.c1); C1W = A_(K.C1W)
    fn = Fn_(al); fw = FnW_(al); a3 = al ** 3
    l = dict(A=A_('5.72') * fn)
    l['B'] = eps * fn * (A_('1.2') + A_('1.1') * (k3 + c1 / k3))
    l['D4'] = eps * fn * (1 + A_('1.1') / k3)
    l['V0'] = mu * fw * (A_('1.67') + (k1 + C1W / k1) / 2)
    l['V2'] = mu * fw / (2 * k1) + mu * eI * a3 * (A_('1.7') + A_('1.305') * (k2 + C1W / k2))
    l['V4'] = mu * eI * a3 * (1 + A_('1.305') / k2)
    if K.FIX and K.CLIN > 0: l['A'] = l['A'] + A_(K.CLIN)            # absolute linear forcing c_j A (F3, pending)
    return l
def phi_val(al, delta, theta, kap):
    n = nk(theta); l = lk(al, kap); d = A_(delta)
    return sum((l[k] * l[k] / (4 * d * n[k]) for k in n), arb(0))
_cache = {}
def Phi(a_hi, delta, theta):
    """returns float upper bound of Phi(a_hi) and the kappas used"""
    key = (a_hi, delta, theta)
    if key in _cache: return _cache[key]
    al = A_(a_hi)
    best = None
    grid = [2.0 ** e for e in range(-6, 12)]
    kap = [1.0, 1.0, 1.0]
    f = lambda kp: float(phi_val(al, delta, theta, kp).upper())
    cur = f(kap)
    for _ in range(3):
        for i in range(3):
            for g in grid:
                kp = list(kap); kp[i] = g; v = f(kp)
                if v < cur: cur, kap = v, kp
    out = (float(np.nextafter(float(phi_val(al, delta, theta, kap).upper()), np.inf)), kap)
    _cache[key] = out
    return out
def qbl_sup(n=4000):
    """rigorous float upper bound of sup_shapes (A+W1+B)^2 / Lyap  (maximal at cw=1; sweep of the angle t, interval arithmetic)"""
    import rhs_iv as R
    from ival import cc, sqrt
    c = R.consts(); t0 = np.arange(n) / n; t1 = (np.arange(n) + 1) / n
    ct, st = R.circ(t0, t1); one = cc(1.0)
    Aa = ct; Bb = st * c.inv_sqrt_eps; W1 = sqrt(Aa * (Bb + c.c1 * Aa)); q = Aa + W1 + Bb
    return float(np.nextafter(float(np.max((q * q).hi)), np.inf))
