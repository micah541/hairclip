"""T6 weights in arb: rigorous upper bounds of every weight over an x-interval [xa, xb] (and ell in [la, lb]).
T5 (R):      rho = exp(-(1+x^2)^(2/3)/4) on [0, 34.25];  W = (1+x^2)^-2 zeta^2, zeta = S((x-10.25)/6) (clipped);
             chi = chi_32 = S((32 - y)/3), y = x - X0 (T7's constraint cutoff; used for ||chi r||_rho).
Lemma F:     T_j = 28 eta_j(y)^2 Yl^(2j-2), y = x - X0, eta_j = S(ln(y/ya)/ln(yb/ya)) clipped, Yl = yl Phi(y/yl), yl = ell - X0,
             Phi(s) = s (s <= 1/2), Phi' = 1 - psi((s - 1/2)/(1/2)), psi(u) = u^5 (21/4 - 49 u/8 + 15 u^2/8), clipped at u = 1.
             hats: T_j'^2/T_j and |T_j''| (residual parts of Lemma F's R^_j).  Designs: 'T5' ramps (14,28) all j;
             'opt4n' ramps (14,22) j <= 3, (22,30) j = 4 (T8 recommended).
T5 literal:  28 eta^2 x^(2j-2), eta = S(ln(y'/14)/ln 2), y' = x + STIP (T1 resid.py definition).
Monotonicity used: rho, W decreasing; zeta, eta, Yl (in y and in yl) increasing; S increasing on [0, 1]."""
from flint import arb

X0 = arb('2.2490776515')
STIP = arb(-2.249077651121112)
S0 = arb(1) / 2
PA, PB, PC = arb(21) / 4, arb(-49) / 8, arb(15) / 8
DESIGNS = {'T5': [(14, 28)] * 5, 'opt4n': [(14, 22)] * 4 + [(22, 30)]}


def S(u):
    return u ** 4 * (35 - 84 * u + 70 * u ** 2 - 20 * u ** 3)


def Sp(u):
    return 140 * u ** 3 * (1 - u) ** 3


def Spp(u):
    return 420 * u ** 2 * (1 - u) ** 2 * (1 - 2 * u)


def clip01(lo, hi):
    """clip an interval [lo, hi] (exact arbs) to [0,1]"""
    return (max(lo, arb(0)), min(hi, arb(1)))


def iv(lo, hi):
    """ball covering [lo, hi]"""
    return arb((lo + hi) / 2, ((hi - lo) / 2).abs_upper())


def up(x):
    return x.upper()


def lo_(x):
    return x.lower()


# ------------------------------------------------------------------------------------------------ T5 weights
def rho_sup(xa, xb):
    xa = max(xa, arb(0))
    return up((-((1 + xa ** 2) ** (arb(2) / 3)) / 4).exp())


def rho_chi2_sup(xa, xb):
    """rho chi_32^2: chi = S((32 - y)/3) is decreasing in x"""
    y = xa - X0
    u = (32 - y) / 3
    c = arb(1) if u >= 1 else (arb(0) if u <= 0 else S(u))
    return up(rho_sup(xa, xb) * c ** 2)


def W_sup(xa, xb):
    return up((1 + xa ** 2) ** -2)


def zeta_sup(xb):
    u = (xb - arb('10.25')) / 6
    if u <= 0: return arb(0)
    if u >= 1: return arb(1)
    return up(S(u))


def Wzeta2_sup(xa, xb):
    return up(W_sup(xa, xb) * zeta_sup(xb) ** 2)


# ------------------------------------------------------------------------------------------------ Lemma F weights
def eta_t(y, ya, yb):
    return (y / ya).log() / (arb(yb) / ya).log()


def eta_sup(ylo, yhi, ya, yb):
    """sup of eta over [ylo, yhi] (eta increasing)"""
    if yhi <= ya: return arb(0)
    t = eta_t(yhi, ya, yb)
    return arb(1) if t >= 1 else up(S(t))


def Phi(s):
    """Phi for an exact or ball s >= 0 (monotone increasing); evaluated by branch"""
    if s <= S0: return s
    u = (s - S0) / (1 - S0)
    if u >= 1: u = arb(1)
    Psi = PA * u ** 6 / 6 + PB * u ** 7 / 7 + PC * u ** 8 / 8
    return S0 + (1 - S0) * (u - Psi)


def Yl(y, yl):
    return yl * Phi(y / yl)


def Ybounds(ylo, yhi, lla, llb):
    """[inf, sup] of Yl over y in [ylo, yhi], yl in [lla, llb] (Yl increasing in both)"""
    return lo_(Yl(ylo, lla)), up(Yl(yhi, llb))


def T_sup(j, ylo, yhi, lla, llb, ramp):
    ya, yb = ramp
    e = eta_sup(ylo, yhi, ya, yb)
    if e == 0: return arb(0)
    Ylo, Yhi = Ybounds(max(ylo, arb(ya)), yhi, lla, llb)
    p = 2 * j - 2
    Yp = Yhi ** p if p >= 0 else Ylo ** p
    return up(28 * e ** 2 * Yp)


def _iv_sup_abs(x):
    return x.abs_upper()


def hats_sup(j, ylo, yhi, lla, llb, ramp):
    """(sup T'^2/T, sup |T''|) over the cell; 0 below the ramp foot"""
    ya, yb = ramp
    if yhi <= ya: return arb(0), arb(0)
    ylo = max(ylo, arb(ya))
    L = (arb(yb) / ya).log()
    tlo, thi = eta_t(ylo, ya, yb), eta_t(yhi, ya, yb)
    tl, th = clip01(lo_(tlo), up(thi))
    if tl >= 1:
        e, e1, e2 = arb(1), arb(0), arb(0)
    else:
        tb = iv(tl, th)
        yb_ = iv(ylo, yhi)
        e = arb(1) if th >= 1 and tl >= 1 else up(S(th)) if th < 1 else arb(1)
        t1 = 1 / (yb_ * L); t2 = -1 / (yb_ ** 2 * L)
        e1 = (Sp(tb) * t1).abs_upper()
        e2 = (Spp(tb) * t1 ** 2 + Sp(tb) * t2).abs_upper()
    Ylo, Yhi = Ybounds(ylo, yhi, lla, llb)
    # Y' = Phi'(s) in [0, 1], Y'' = Phi''(s)/yl
    slo, shi = lo_(ylo / llb), up(yhi / lla)
    if shi <= S0:
        Y1, Y2 = arb(1), arb(0)
    else:
        ulo = max(lo_((slo - S0) / (1 - S0)), arb(0)); uhi = min(up((shi - S0) / (1 - S0)), arb(1))
        ub_ = iv(ulo, uhi)
        psip = 5 * PA * ub_ ** 4 + 6 * PB * ub_ ** 5 + 7 * PC * ub_ ** 6
        Y1 = arb(1)
        Y2 = (psip / (1 - S0) / lla).abs_upper()
    p = 2 * j - 2
    Yp = Yhi ** p if p >= 0 else Ylo ** p                 # sup Y^p
    Ypm1 = (Yhi ** (p - 1) if p - 1 >= 0 else Ylo ** (p - 1))
    Ypm2 = (Yhi ** (p - 2) if p - 2 >= 0 else Ylo ** (p - 2))
    ap = abs(p)
    # T'^2/T = 28 Y^p (2 eta' + p eta Y'/Y)^2
    br = 2 * e1 + ap * e * Y1 / Ylo
    hA = up(28 * Yp * br ** 2)
    # T'' = 28 [ (eta^2)'' Y^p + 2 (eta^2)' (Y^p)' + eta^2 (Y^p)'' ]
    d1 = 2 * e * e1; d2 = 2 * e1 ** 2 + 2 * e * e2
    Yp1 = ap * Ypm1 * Y1
    Yp2 = abs(p * (p - 1)) * Ypm2 * Y1 ** 2 + ap * Ypm1 * Y2
    hB = up(28 * (d2 * Yp + 2 * d1 * Yp1 + e ** 2 * Yp2))
    return hA, hB


# ------------------------------------------------------------------------------------------------ T5 literal (resid.py)
def Tlit_sup(j, xa, xb):
    yhi = xb + STIP
    if yhi <= 14: return arb(0)
    t = (yhi / 14).log() / arb(2).log()
    e = arb(1) if t >= 1 else up(S(t))
    p = 2 * j - 2
    xp = xb ** p if p >= 0 else max(xa, arb('0.5')) ** p
    return up(28 * e ** 2 * xp)
