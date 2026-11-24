#!/usr/bin/env python3
"""Rigorous (arb ball arithmetic) constants for weighted interpolation inequalities
and cutoff-commutator bounds.  Writes sob_constants.json next to this file."""
import json, os, sys
from fractions import Fraction
from flint import arb, acb, ctx

ctx.prec = 100
HERE = os.path.dirname(os.path.abspath(__file__))
X = 34.25
A_W = 6.0
H = 2.0 ** -12           # cell width (dyadic => all breakpoints exact)
OUT = {}


# ---------------------------------------------------------------- helpers
def chk(v):
    assert not v.is_nan(), "NaN ball"
    return v


def to_frac(b):
    """exact Fraction of the midpoint of a zero-radius ball"""
    m, e = b.mid().man_exp()
    m = int(m); e = int(e)
    return Fraction(m) * (Fraction(2) ** e)


def _round(fr, sig, up):
    import math
    if fr == 0:
        return Fraction(0)
    neg = fr < 0
    a = -fr if neg else fr
    e = 0
    # find e with 10^e <= a < 10^(e+1)
    e = int(math.floor(math.log10(float(a)))) if a > 0 and 1e-300 < float(a) < 1e300 else None
    if e is None:
        e = len(str(a.numerator)) - len(str(a.denominator))
    while Fraction(10) ** e > a:
        e -= 1
    while Fraction(10) ** (e + 1) <= a:
        e += 1
    unit = Fraction(10) ** (e - sig + 1)
    q = a / unit
    n = q.numerator // q.denominator
    if up and n * unit < a:
        n += 1
    if (not up) and neg:      # rounding the magnitude up to round the value down
        if n * unit < a:
            n += 1
    if neg:
        return -(n * unit)
    return n * unit


def up4(b, sig=4):
    """round an upper bound of ball b UP to sig significant digits (as float)"""
    chk(b)
    return float(_round(to_frac(b.upper()), sig, True))


def dn4(b, sig=4):
    """round a lower bound of ball b DOWN to sig significant digits (as float)"""
    chk(b)
    lo = to_frac(b.lower())
    return float(-_round(-lo, sig, True)) if lo != 0 else 0.0


def ball(lo, hi):
    return arb((lo + hi) / 2, (hi - lo) / 2)    # lo,hi exact dyadic floats


def cells(lo, hi, h=H):
    n = int(round((hi - lo) / h))
    assert abs(lo + n * h - hi) < 1e-12, "grid does not hit endpoint"
    for i in range(n):
        yield lo + i * h, lo + (i + 1) * h


# ---------------------------------------------------------------- rho side
def phi(x):
    return (arb(1) + x * x) ** (arb(2) / 3) / 4


def dphi(x):
    q = arb(1) + x * x
    return x * (q ** (-arb(1) / 3)) / 3


def ddphi(x):
    q = arb(1) + x * x
    q13 = q ** (-arb(1) / 3)
    return q13 / 3 - arb(2) / 9 * (x * x) * q13 / q


def sup_cells(lo, hi, f):
    best = None
    for l, h in cells(lo, hi):
        v = chk(f(ball(l, h))).upper()
        if best is None or v > best:
            best = v
    return best


def min_lower_cells(lo, hi, f):
    best = None
    for l, h in cells(lo, hi):
        v = chk(f(ball(l, h))).lower()
        if best is None or v < best:
            best = v
    return best


def rho_part():
    R = {}
    # phi' > 0 on (0,X] and phi'' > 0 (so phi' increasing, phi increasing, rho decreasing)
    R["phi_dd_min_lower_bound"] = float(min_lower_cells(0.0, X, ddphi).mid())
    assert R["phi_dd_min_lower_bound"] > 0
    # (phi'(x) > 0 for x > 0 is evident from the closed form; check on cells away from 0)
    assert min_lower_cells(H, X, dphi) > 0
    sup_expr = sup_cells(0.0, X, lambda x: dphi(x) * dphi(x) - ddphi(x))
    if sup_expr < 0:
        sup_expr = arb(0)
    c1 = sup_expr / 2
    R["sup_phi1sq_minus_phi2_plus"] = up4(sup_expr)
    R["c1"] = up4(c1)
    Phi1X = sup_cells(0.0, X, dphi)      # sup over [0,X] >= sup over [x,X]; equals phi'(X) by monotonicity
    PhiX_pt = dphi(arb(X))
    R["Phi1_X"] = up4(PhiX_pt)
    R["Phi1_X_cellsup"] = up4(Phi1X)
    regions = [(0, 1), (1, 2), (2, 4), (4, 8), (8, 12), (12, 20), (20, X)]
    K = {}
    for (l, r) in regions:
        # rho decreasing => min rho on [l,r] = rho(r); K = exp(phi(r))
        Kv = chk(phi(arb(r)).exp())
        K["[%g,%g]" % (l, r)] = up4(Kv)
    R["K_region"] = K
    # I4
    for b in (2, 12):
        R["K_b%d" % b] = up4(chk(phi(arb(b)).exp()))
    R["rho_at_X"] = up4(chk((-phi(arb(X))).exp()))
    return R


# ---------------------------------------------------------------- W side
def Spoly(t):
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    return t4 * (35 - 84 * t + 70 * t2 - 20 * t3)


def Sprime(t):          # S'(t) = 140 t^3 (1-t)^3
    u = arb(1) - t
    return 140 * (t * t * t) * (u * u * u)


def psi_pp(x, lo, hi):
    """returns (psi', psi'') on the ball x (cell aligned to 30.25 / 38.25 breakpoints)"""
    q = arb(1) + x * x
    q13 = q ** (-arb(1) / 3)
    A = 6 * x / q
    Ap = 6 * (arb(1) - x * x) / (q * q)
    B = x * q13 / 3
    Bp = q13 / 3 - arb(2) / 9 * (x * x) * q13 / q
    if hi <= 30.25:
        return A, Ap
    if lo >= 38.25:
        return B, Bp
    assert lo >= 30.25 and hi <= 38.25
    t = (x - arb(30.25)) / 8
    s = Spoly(t)
    sp = Sprime(t) / 8
    d1 = (arb(1) - s) * A + s * B
    d2 = (arb(1) - s) * Ap + s * Bp + sp * (B - A)
    return d1, d2


def psi_point_deriv(xv):
    """psi'(x) for a single exact x"""
    x = arb(xv)
    return psi_pp(x, xv, xv)[0]


def psi_value(xv):
    """psi(x) normalised so that exp(-psi) = (1+x^2)^-3 for x<=30.25"""
    x = arb(xv)
    q = arb(1) + x * x
    if xv <= 30.25:
        return 3 * q.log()
    base = 3 * (arb(1) + arb(30.25) * arb(30.25)).log()

    def integrand(z, _an):
        z2 = z * z
        qq = acb(1) + z2
        q13 = (qq.log() * (-acb(1) / 3)).exp()
        A = 6 * z / qq
        B = z * q13 / 3
        t = (z - acb("30.25")) / 8
        t2 = t * t; t3 = t2 * t; t4 = t2 * t2
        s = t4 * (35 - 84 * t + 70 * t2 - 20 * t3)
        return (1 - s) * A + s * B

    if xv <= 38.25:
        I = acb.integral(integrand, acb(30.25), acb(xv))
        return base + I.real
    I = acb.integral(integrand, acb(30.25), acb(38.25))
    p38 = base + I.real
    return p38 + phi(x) - phi(arb(38.25))


def W_part():
    R = {}
    RS = [40, 60, 80]
    # cell data on [6, 80]
    cell_lo = []
    up_d1 = []
    up_c1 = []
    min_d1 = None
    n = int(round((80 - A_W) / H))
    for i in range(n):
        l = A_W + i * H; h = l + H
        xb = ball(l, h)
        d1, d2 = psi_pp(xb, l, h)
        chk(d1); chk(d2)
        cell_lo.append(l)
        up_d1.append(d1.upper())
        up_c1.append((d1 * d1 - d2).upper())
        lw = d1.lower()
        if min_d1 is None or lw < min_d1:
            min_d1 = lw
    R["psi1_min_lower_bound_on_[6,80]"] = float(min_d1.mid())
    assert min_d1 > 0
    idx = lambda x: int(round((x - A_W) / H))
    R["c1W"] = {}
    R["Psi1"] = {}
    for Rr in RS:
        m = max(up_c1[: idx(Rr)])
        if m < 0:
            m = arb(0)
        R["c1W"]["R=%d" % Rr] = up4(m / 2)
        R["Psi1"]["R=%d" % Rr] = {}
        for left in (6, 12, 20, 30, 40, 60):
            if left >= Rr:
                continue
            mm = max(up_d1[idx(left): idx(Rr)])
            R["Psi1"]["R=%d" % Rr]["from_%d" % left] = up4(mm)
    # K_W per region: psi increasing (psi'>0 certified) => min W on [l,r] = W(r), K_W = exp(psi(r))
    regs = [(12, 20), (20, 30), (30, 40), (40, 60), (60, 80)]
    KW = {}
    psiv = {}
    for (l, r) in regs:
        pv = chk(psi_value(float(r)))
        psiv[r] = pv
        KW["[%d,%d]" % (l, r)] = up4(pv.exp())
        KW["[%d,%d]_psi(right)" % (l, r)] = up4(pv)
    R["K_W_region"] = KW
    R["W_normalisation"] = "W = (1+x^2)^(-3) for x <= 30.25 (psi(x) = 3 log(1+x^2) there)"
    # cross-check: Riemann-sum upper bound for psi(38.25)-psi(30.25)
    s = arb(0)
    for i in range(idx(30.25), idx(38.25)):
        s += up_d1[i] * H
    R["check_psi38.25_minus_psi30.25_riemann_upper"] = float(s.mid())
    R["check_psi38.25_minus_psi30.25_integral"] = float((psi_value(38.25) - psi_value(30.25)).mid())
    return R


# ---------------------------------------------------------------- part (b)
def part_b():
    import sympy as sp
    t = sp.symbols("t")
    S = t ** 4 * (35 - 84 * t + 70 * t ** 2 - 20 * t ** 3)
    Sd = [sp.Poly(sp.diff(S, t, j), t) for j in range(5)]
    supS = {}
    # rigorous: cells + interval Horner of exact rational coefficients
    for j in range(1, 5):
        coeffs = [arb(int(c.p)) / int(c.q) for c in Sd[j].all_coeffs()]
        N = 1 << 19
        best = None
        for i in range(N):
            tb = arb((2 * i + 1) / (2.0 * N), 1.0 / (2 * N))
            v = arb(0)
            for c in coeffs:
                v = v * tb + c
            chk(v)
            u = max(v.upper().abs_upper(), (-v).upper().abs_upper()) if False else None
            a_up = v.abs_upper()
            if best is None or a_up > best:
                best = a_up
        supS[j] = best
    # exact check via critical points
    exact = {}
    for j in range(1, 5):
        crit = [sp.Rational(0), sp.Rational(1)] + [r for r in sp.real_roots(Sd[j + 1]) if 0 <= r <= 1] if j < 4 else [sp.Rational(0), sp.Rational(1)]
        if j == 4:
            crit += [r for r in sp.real_roots(sp.Poly(sp.diff(S, t, 5), t)) if 0 <= r <= 1]
        exact[j] = max(abs(sp.N(Sd[j].eval(r), 30)) for r in crit)
    Rb = {"supS_cells_upper": {j: float(supS[j].mid()) for j in supS},
          "supS_exact_check": {j: float(exact[j]) for j in exact}}
    # chi^(j) = (-1/3)^j S^(j)(u)
    supchi = {j: supS[j] / (arb(3) ** j) if False else supS[j] / arb(3 ** j) for j in range(1, 5)}
    Rb["sup_chi_j"] = {j: up4(supchi[j]) for j in supchi}
    eps = arb("1e-20")
    k2 = eps * eps
    qb = arb(32) + eps                    # |q| <= sigma + 1e-20 <= 32 + 1e-20
    coef_a = 5 * eps * eps + qb / 4       # |5km + q/4|
    C3 = 4 * supchi[1]
    C2 = 6 * supchi[2]
    C1 = 4 * supchi[3] + 2 * k2 * supchi[1]
    C0 = supchi[4] + k2 * supchi[2] + coef_a * supchi[1]
    Rb["C"] = {"C0": up4(C0), "C1": up4(C1), "C2": up4(C2), "C3": up4(C3)}
    Rb["C_unrounded_upper"] = {"C0": float(C0.upper().mid()), "C1": float(C1.upper().mid()),
                               "C2": float(C2.upper().mid()), "C3": float(C3.upper().mid())}
    Rb["bound_on_|5km+q/4|"] = up4(coef_a)
    # strip, rho_max
    tip_lo, tip_hi = arb("-2.2490776516"), arb("-2.2490776505")
    # sigma in [29,32] => x = sigma - tip in [29 - tip_hi, 32 - tip_lo]
    x_lo_enc = (arb(29) - tip_hi).lower()      # rigorous lower bound of x on supp(chi')
    x_hi_enc = (arb(32) - tip_lo).upper()
    Rb["x_strip_enclosure"] = [float(x_lo_enc.mid()), float(x_hi_enc.mid())]
    assert float(x_lo_enc.mid()) >= 31.249 and float(x_hi_enc.mid()) <= 34.250
    rho_max_tight = (-phi(x_lo_enc)).exp()
    rho_max_stated = (-phi(arb(31.249))).exp()
    Rb["rho_max"] = up4(rho_max_stated)
    Rb["rho_max_tight_from_enclosure"] = up4(rho_max_tight)
    Rb["len_exact_support_of_chi_prime"] = 3.0
    Rb["len_stated_strip"] = 3.001
    Rb["rho_max_times_len"] = up4(rho_max_stated * arb("3.001"))
    return Rb


def main():
    OUT["setting"] = {"X": X, "a": A_W, "cell_width": H}
    OUT["rho"] = rho_part()
    OUT["W"] = W_part()
    OUT["commutator"] = part_b()
    with open(os.path.join(HERE, "sob_constants.json"), "w") as f:
        json.dump(OUT, f, indent=2)
    print(json.dumps(OUT, indent=2))


if __name__ == "__main__":
    main()
