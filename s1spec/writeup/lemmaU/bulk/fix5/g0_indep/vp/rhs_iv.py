"""Interval version of rhs() of ../inner_closure_float2.py, written as  tot = POS - NEG  (every summand of POS and NEG is a
nonnegative interval).  MONOTONICITY (needed because W1,W3,V1,V3 are replaced by their interpolation upper bounds, and the
sup-bounds S[..], supv1, supv2 by their upper bounds): every one of these quantities, and |alpha|, enters POS only through
sums/products/maxima with NONNEGATIVE coefficients (asserted structurally: POS is built only from P-arithmetic, which has no
subtraction, and all constants are asserted >= 0 by cc()).  The only negative terms are the explicit quadratics in NEG, which
involve only the free variables A,B,D4,V0,V2,V4 (never W1,W3,V1,V3,S,supv).  See VERIFIER.md.
Variables: parametrisation  A = s cw ct, B = s cw st/sqrt(eps), V0 = s sw cI/sqrt(mu), V2 = s sw sI/sqrt(mu epsI),
r = s q,  D4 = r cu/sqrt(eps), V4 = r su/sqrt(mu epsI),   s = sqrt(Ly);  cx^2+sx^2 = 1  (so Lyap = Ly exactly)."""
import numpy as np
from ival import *
import constants as K

def _c(x):
    """constant -> interval.  arb (derived constants, computed at 200 bit) or decimal literal (float(), widened by 1 ulp)"""
    return arb2P(x) if isinstance(x, arb) else cc(x)
def _A(v): return arb(str(v))      # decimal string -> exact arb (never the binary double)
class Consts:
    def __init__(self):
        k = K
        self.c1 = _c(k.c1); self.PHX = _c(k.PHX); self.KW = _c(k.KW); self.PSI1W = _c(k.PSI1W); self.C1W = _c(k.C1W)
        self.KR16 = arb2P(k.KR16_arb); self.K = {r: (self.KR16 if v is None else _c(v)) for r, v in k.KREG.items()}
        self.K2 = _c(k.K2); self.K12 = _c(k.K12)
        self.sqrtPHI34 = arb2P(_A(k.PHI34).sqrt())
        self.sw3125 = arb2P((6 * k.W3125_arb).sqrt())
        self.inv_sqrt_eps = arb2P(1 / _A(k.EPS).sqrt()); self.inv_sqrt_mu = arb2P(1 / _A(k.MU).sqrt())
        self.inv_sqrt_muepsI = arb2P(1 / (_A(k.MU) * _A(k.EPSI)).sqrt())
        self.EPS = _c(k.EPS); self.MU = _c(k.MU); self.EPSI = _c(k.EPSI)
        self.MUEPSI = _c(_A(k.MU) * _A(k.EPSI))
        # cubic groups: per region, per d: list of (k0,k1,coef)
        self.cub = {}
        for reg, Kr in k.KREG.items():
            groups = {}
            for key, C in k.cubf[reg].items():
                assert C >= 0
                ab = eval(key); mono = [0] * ab[0] + [1] * ab[1] + [2] * ab[2]; mono.sort(reverse=True)
                kk = (mono[0], mono[1], tuple(mono[2:])); groups[kk] = groups.get(kk, 0.0) + 0.5 * C
            self.cub[reg] = [(a, b, d, P(np.float64(c * (1 - 1e-13)), np.float64(c * (1 + 1e-13)))) for (a, b, d), c in groups.items()]   # d: tuple of 'rest' indices
        self.cl = {n: _c(v) for n, v in dict(
            lin1=arb(1) / 32 - _A(k.PI3), one9=1 / 9, alin_sq=k.ALIN["sq"], alin_AB=k.ALIN["AB"], alin_W1B=k.ALIN["W1B"], alin_BB=k.ALIN["BB"],
            mu1t=_A(k.MU1) * (1 - _A(k.THETA)), halft=(1 - _A(k.THETA)) / 2, mu2t=_A(k.MU2) * (1 - _A(k.THETA)), mod=k.MOD_ALPHA, f572=5.72, SB_SA=_A(k.SB) + _A(k.SA), twoIA=2 * _A(k.IA), two=2.0, half=0.5,
            nd1=6.9e-3, nd2=1.2e-2, nd3=4.2e-2, nd4=0.54, nd5=0.62, nd6=4.75, nd7=0.136, nd8=3.05,
            nf1=4.3e-3, nf2=6.9e-4, nf3=2.9e-5, nf4=0.014, nf5=0.46, nf6=3.06, nf7=8.6e-6,
            DS0=k.DS[0], DS1=k.DS[1], DS2=k.DS[2], DS3=k.DS[3], t22=2.2, t12=1.2, r143=14.34, r163=16.3, r107=10.7, modE2=k.MOD_E2,
            MU1=k.MU1, MU2=k.MU2, AF=k.AF, BF=k.BF, d001=k.TQ['d001'], d010=k.TQ['d010'], d100=k.TQ['d100'], A1=k.TQ['A1'], A0=k.TQ['A0'],
            d200=k.TQ['d200'], d110=k.TQ['d110'], d300=k.TQ['d300'], t167=1.67, CABQ=_A(k.CAB) * (_A('0.273') + _A('0.0159')), CAB3985=_A(k.CAB) * (_A('38.5') + _A('0.66')),
            CAB25=_A(k.CAB) * _A('2.5'), CF0=k.CF[0], CF1=k.CF[1], CF2=k.CF[2], CF3=k.CF[3], t261=2.61, t17=1.7, one=1.0, c62=62.0, CAB03=_A(k.CAB) * _A('0.3'),
            F0=k.FN[0], F1=k.FN[1], F2=k.FN[2], G0=k.FNW[0], G1=k.FNW[1], G2=k.FNW[2],
            ar0=k.ADD_I0_RAMP[0], ar1=k.ADD_I0_RAMP[1], ar2=k.ADD_I0_RAMP[2], ar3=k.ADD_I0_RAMP[3], ar4=k.ADD_I0_RAMP[4], ar5=k.ADD_I0_RAMP[5],
            am_a=_A(k.CAB0) / _A(k.CAB), am_al=_A(k.MAB) / _A(k.CAB), am_q=(_A(k.QAB) if k.MODIBP else _A(k.QM)) / _A(k.CAB), am_phi=_A(k.PHIAB) / _A(k.CAB), am_strip=_A(k.STRIP_AB) / _A(k.CAB),
            negscale=_A(1) - _A(k.DELTA_NEG), s0v=k.S0V, s1v32=_A(k.S1V) / 32, qb752=_A('7.52') * _A(k.QB), qm_cab=_A(k.QM) / _A(k.CAB), aflux0a=k.ABS_FLUX_I0, aflux2a=k.ABS_FLUX_I2,
            ac0w=k.ADD_I0_C0W, c0pn=k.C0P_NORM, cflux=k.CFLUX, aflux0=k.ADD_FLUX_I0, aflux2=k.ADD_FLUX_I2, astrip=k.ADD_STRIP_CUBIC).items()}
        self.g0 = {n: _c(v) for n, v in k.G0.items()}; self.CABP = _c(_A(k.CAB) * _A(k.G0['cpsi']))
        self.ramp = {lev: [(_c(c0), _c(dw)) for (c0, dw) in tab] for lev, tab in k.RAMP_TAB.items()}
def set_mode(delta=0.0, dropL=False):
    global _CC
    K.DELTA_NEG = float(delta); K.NEGSCALE = 1.0 - float(delta); K.DROPL = bool(dropL); _CC = None
def radial_factor():
    """NEG >= radial_factor * (eps D4^2 + mu epsI V4^2)/2 : lower-rounded float, used in find_qstar/qcut"""
    return (1.0 - K.DELTA_NEG) * (1.0 - K.THETA) * (1 - 1e-12)
def set_theta(th):
    global _CC
    K.THETA = float(th); _CC = None
def gm(a, b):
    """sqrt(a*b) for nonneg intervals.  In gradient (PG) mode an AM-GM majorant 1/2 (th a + b/th) with th^2 = mid(b)/mid(a) is
    used instead (valid upper bound for ANY th>0; Lipschitz, no 1/sqrt singularity)."""
    if isinstance(a, PG) and isinstance(b, PG) or (isinstance(a, PG) or isinstance(b, PG)):
        am = 0.5 * (a.lo + a.hi); bm = 0.5 * (b.lo + b.hi)
        ok = (am > 0) & (bm > 0)
        th = np.sqrt(np.where(ok, bm, 1.0) / np.where(ok, am, 1.0))
        thP = P(fd(th), fu(th)); ithP = P(fd(1.0 / th), fu(1.0 / th))      # th is a float; thP contains th, ithP contains 1/th (valid majorant for any th>0)
        return cc(0.5) * (thP * a + ithP * b)
    return sqrt(a * b)

_CC = None
def consts():
    global _CC
    if _CC is None: _CC = Consts()
    return _CC

def Fn(al, c):
    return c.cl['F0'] * al + c.cl['F1'] * (al * al) + c.cl['F2'] * (al * al * al)
def FnW(al, c):
    return c.cl['G0'] * al + c.cl['G1'] * (al * al) + c.cl['G2'] * (al * al * al)

def rhs_pn(al, A, B, D4, V0, V2, V4, parts=False, floors=True):
    """al: P (scalar, >=0). A..V4: P arrays (physical units, i.e. s * normalised).  Returns (H, L, NEG):
       H   = POS terms of x-degree >= 2  (every term scales like s^d with d>=2, so H/s^2 is nondecreasing in s=sqrt(Ly))
       L   = POS terms of x-degree 1 and 0 (alpha-forcing Fn*A etc., the c0' w coupling, the flux constants; L/s^2 nonincreasing in s)
       NEG = explicit negative quadratics (NEG/s^2 independent of s).   tot = H + L - NEG.
       floors=False drops the terms of joint (alpha,x)-degree < 2 (the added b,c terms): used only for the scale-invariant tail lemma."""
    c = consts(); cl = c.cl; ADD = K.ADD
    FL = floors  # floors=False drops the degree-<2 terms (ac0w V0 and the two flux constants): used for the scale-invariant tail lemma
    AA = A * A; BB = B * B; V0V0 = V0 * V0; V2V2 = V2 * V2
    W1 = gm(A, B + c.c1 * A); W3 = gm(B, D4 + c.c1 * B)
    V1 = gm(V0, V2 + c.C1W * V0); V3 = gm(V2, V4 + c.C1W * V2)
    W1W1 = W1 * W1; AW1 = A * W1; W1B = W1 * B; BW3 = B * W3
    def sset(Kc):
        return (gm(A, Kc * (c.PHX * A + cl['two'] * W1)), gm(W1, Kc * (c.PHX * W1 + cl['two'] * B)), gm(B, Kc * (c.PHX * B + cl['two'] * W3)))
    mid = sset(c.K12); tip = sset(c.K2)
    sr = sset(c.KR16)
    sI1 = gm(V0, c.KW * (c.PSI1W * V0 + cl['two'] * V1)); sI2 = gm(V1, c.KW * (c.PSI1W * V1 + cl['two'] * V2))
    if K.FIX:       # F2: on [34.25,60] the pointwise data are the bootstrap constants S0_v = sup|w'|, S1_v/32 >= sup|w''| (Lemma S); sup over [16.25,60] <= max
        sI1 = mx(sI1, cl['s0v']); sI2 = mx(sI2, cl['s1v32'])
    sI0 = sr[0] + c.sqrtPHI34 * V0
    far = (mx(sr[0], sI0), mx(sr[1], sI1), mx(sr[2], sI2))
    L2 = (A, W1, B)
    # ---- E0 ----
    nlin = cl['lin1'] * AA + cl['one9'] * BB
    alin = al * (cl['alin_sq'] * (AA + W1W1) + cl['alin_AB'] * (A * B) + cl['alin_W1B'] * W1B + cl['alin_BB'] * BB)
    mod = cl['mod'] * al * AA
    fn = Fn(al, c)
    forc = cl['f572'] * fn * A       # [L group]
    cub = None; cubfar = None; regS = {}
    for reg in K.KREG:
        if K.KREG[reg] is None: s = (sI0, sI1, sI2)
        else: s = sset(c.K[reg])
        regS[reg] = s
        reg_t = None
        for (k0, k1, d, coef) in c.cub[reg]:
            t = coef * L2[k0] * L2[k1]
            for dd in d: t = t * s[dd]
            reg_t = t if reg_t is None else reg_t + t
        cub = reg_t if cub is None else cub + reg_t
        if K.KREG[reg] is None: cubfar = reg_t
    if 'd' in ADD: cub = cub + cl['astrip'] * cubfar
    strip = cl['SB_SA'] * V0V0 + cl['twoIA'] * c.KR16 * (AA + cl['two'] * AW1)
    R0p = alin + mod + cub + strip
    _q = A + W1 + B
    if K.FIX: R0p = R0p + cl['qb752'] * (_q * _q) * A        # F3 placeholder: quadratic modulation form q_beta (A+W1+B)^2 * A in E0
    L = forc
    # ---- E2 ----
    supw2 = mx(tip[2], mid[2]); supw2f = far[2]
    ND = (cl['nd1'] * D4 + cl['nd2'] * W3 + cl['nd3'] * B + cl['nd4'] * A + cl['nd5'] * W1 + cl['nd6'] * supw2 * B + cl['nd7'] * supw2 * W3
          + cl['nd8'] * (supw2 * supw2) * B + cl['nf1'] * D4 + cl['nf2'] * W3 + cl['nf3'] * B + cl['nf4'] * supw2f * B + cl['nf5'] * supw2f * W3
          + cl['nf6'] * (supw2f * supw2f) * B + cl['nf7'] * A + cl['nf7'] * W1)
    wstrip = c.sw3125 * (sr[0] + c.sqrtPHI34 * V0)
    Dq = D4 + cl['t22'] * W3 + cl['t12'] * B
    R2p = ((cl['DS0'] * V2 + cl['DS1'] * V1 + cl['DS2'] * V0 + cl['DS3'] * wstrip) * Dq + cl['r143'] * BB + cl['r163'] * W1B + cl['r107'] * (A * B)
           + ND * Dq + cl['modE2'] * AA)
    L = L + c.EPS * (fn * Dq)
    # ---- I0 ----
    supv1 = sI2
    Vq = V1 + cl['t167'] * V0
    Amx = (cl['am_a'] * A + cl['am_al'] * al * A + cl['am_q'] * (_q * _q)) if K.MODIBP else (A + cl['am_q'] * (_q * _q) if K.FIX else A)
    Amc = (cl['am_phi'] * FnW(al, c) + cl['am_strip']) if K.MODIBP else None      # x-degree 0 part of the (a,b) bound: its products with al*V0, al*D2v are x-degree 1 -> L group
    Am = Amx + Amc if Amc is not None else Amx
    #      # F3 placeholder: CAB*A -> CAB*A + q_(a,b) (A+W1+B)^2 in the a,b-modulation terms of I0, I2
    I0p = (c.cl['AF'] * W1W1 + cl['BF'] * BB + cl['d001'] * V2V2 + (cl['d010'] * V2 + cl['d100'] * V1 + cl['A1'] * V0) * Vq
           + (cl['d200'] * supv1 * V1 + cl['d110'] * supv1 * V2 + cl['d300'] * (supv1 * supv1) * V1) * Vq
           + cl['CABQ'] * Amx * al * V0 + cl['CAB3985'] * Am * V0V0 + cl['CAB25'] * Am * W1W1)
    L = L + c.MU * (FnW(al, c) * Vq)
    if Amc is not None: L = L + c.MU * (cl['CABQ'] * Amc * al * V0)
    if 'a' in ADD:
        if K.RAMP == "table":
            lev = K.ramp_level(float(np.max(al.hi)))
            S_ramp = mx(mx(mid[0], mid[1]), mx(regS["12<x<=16.25"][0], regS["12<x<=16.25"][1]))     # sup(|w|,|w'|) on [10.25,16.25] (regions 8<x<=12, 12<x<=16.25)
            tab = c.ramp[lev]; xs5 = (D4, W3, B, W1, A)
            ramp = None
            for (c0, dw), xj in zip(tab, xs5):
                t = (c0 + dw * S_ramp) * xj
                ramp = t if ramp is None else ramp + t
            I0p = I0p + cl['ar0'] * ramp * Vq
        else:
            I0p = I0p + cl['ar0'] * (cl['ar1'] * D4 + cl['ar2'] * W3 + cl['ar3'] * B + cl['ar4'] * W1 + cl['ar5'] * A) * Vq
    if 'b' in ADD:
        if K.QUAD:
            S_all = None
            for reg in K.KREG: S_all = regS[reg][0] if S_all is None else mx(S_all, regS[reg][0])      # sup|w| over every region (0..60)
            I0p = I0p + cl['c0pn'] * S_all * V0
        elif FL: L = L + c.MU * (cl['ac0w'] * V0)
    if 'c' in ADD:
        if K.CMODE == "abs":
            if FL: L = L + c.MU * cl['aflux0a']                       # F4: absolute forcing (degree 0)
        elif K.QUAD: I0p = I0p + cl['cflux'] * (V0V0 + V1 * V1 + V2V2)
        elif FL: L = L + c.MU * cl['aflux0']
    # ---- I2 ----
    D2v = V4 + cl['t261'] * V3 + cl['t17'] * V2
    supv2 = gm(V2, c.KW * (c.PSI1W * V2 + cl['two'] * V3))
    I2p = ((cl['CF0'] * D4 + cl['CF1'] * W3 + cl['CF2'] * B + cl['CF3'] * W1) * D2v + cl['d001'] * (V4 * V4)
           + (cl['d010'] * V3 + cl['d110'] * supv1 * V3 + cl['two'] * cl['d200'] * supv1 * V3 + (cl['d110'] * supv1 * (cc(3.0) * V3 + cc(1.31) * V2) if K.GN else cl['d110'] * supv2 * V2) + cl['d300'] * (supv1 * supv1) * V2) * D2v
           + cl['CAB3985'] * Am * V2V2 + cl['c62'] * Am * V2 * V3 + cl['CAB03'] * Amx * al * D2v)
    L = L + c.MUEPSI * ((al * al * al) * D2v)
    if Amc is not None: L = L + c.MUEPSI * (cl['CAB03'] * Amc * al * D2v)
    if 'c' in ADD:
        if K.CMODE == "abs":
            if FL: L = L + c.MUEPSI * cl['aflux2a']
        elif K.QUAD: I2p = I2p + cl['cflux'] * (V2V2 + V3 * V3 + V4 * V4)
        elif FL: L = L + c.MUEPSI * cl['aflux2']
    # ======================= G0_INDEP review patch =======================
    G0F = K.G0FIX
    if G0F:
        g = c.g0; three = cc(3.0)
        Sr2 = mx(mid[2], regS["12<x<=16.25"][2])                 # sup|w_xx| on the inner ramp [10.25,16.25] (rho-Sobolev, regions 8..12, 12..16.25)
        betaH = (g['b0'] + g['b1'] * al) * A + g['bq'] * (_q * _q)   # x-degree>=1 part of the |beta| bound
        betaL = g['bF'] * fn                                     # x-degree 0 part (strip constant NOT included: see report)
        cb = g['cb0'] + g['cb1'] * al
        if 'x' in G0F:    # level 0: cross terms of the twice-integrated top-order term (1-J^-4) v_xxx
            I0p = I0p + (g['xa1'] + g['xb1'] * supv1) * (V1 * V2) + (g['xa0'] + g['xb0'] * supv1) * (V0 * V2)
        if 'b' in G0F:    # beta coupling, levels 0 and 2
            I0p = I0p + cb * betaH * V0; L = L + c.MU * (cb * betaL * V0)
            I2p = I2p + cb * betaH * D2v; L = L + c.MUEPSI * (cb * betaL * D2v)
        if 'p' in G0F:    # -a psi1_x - b psi2_x
            I0p = I0p + c.CABP * Amx * V0; I2p = I2p + c.CABP * Amx * D2v
            if Amc is not None:
                AmcF = cl['am_phi'] * FnW(al, c)        # joint (alpha,x)-degree 2 with V0: homogeneous part
                L = L + c.MU * (c.CABP * AmcF * V0) + c.MUEPSI * (c.CABP * AmcF * D2v)
                if FL: L = L + c.MU * (c.CABP * cl['am_strip'] * V0) + c.MUEPSI * (c.CABP * cl['am_strip'] * D2v)   # strip constant of |(a,b)|: degree-1 FLOOR term (f1-type), NOT in closure.tail's floor formula (size 1e-31, see report)
        if 'm' in G0F:    # level 0: quadratic/cubic monomials of D on the inner ramp
            I0p = I0p + cl['ar0'] * (g['d200r'] * Sr2 * B + g['d110r'] * Sr2 * W3 + g['d300r'] * (Sr2 * Sr2) * B) * Vq
        if 'q' in G0F:    # level 2, far region: missing cross terms / factors
            I2p = I2p + (cl['d001'] * V4 * (cl['t261'] * V3 + cl['t17'] * V2) + (g['cp1'] + g['cp2'] * supv1) * V3 * D2v
                         + (cl['two'] * cl['d200'] * supv1 * V2 + cl['two'] * cl['d300'] * (supv1 * supv1) * V2) * D2v
                         + (cl['two'] * cl['d100'] * V2 + cl['two'] * cl['A1'] * V1 + cl['two'] * cl['A1'] * V0) * D2v + g['e1'] * (V1 + V0) * V2)
        if 'f' in G0F:    # level 2 forcing: alpha- and alpha^2-residual derivative norms (v3 F1b)
            L = L + c.MUEPSI * ((g['f21'] * al + g['f22'] * (al * al)) * D2v)
        if 's' in G0F and K.RAMP == "table" and 'a' in ADD:    # level 2 on the inner ramp (linear deviation + monomials), [est] model
            I2p = I2p + cl['ar0'] * g['r2lin'] * ramp * D2v
            I2p = I2p + (g['d110r'] * Sr2 * (V3 + cl['ar0'] * (g['r3a'] * W3 + g['r3b'] * B + g['r3c'] * W1)) + cl['ar0'] * g['d110r'] * supv2 * W3
                         + cl['ar0'] * (cl['two'] * g['d200r'] * Sr2 * W3 + three * g['d300r'] * (Sr2 * Sr2) * W3)) * D2v
    H = R0p + c.EPS * R2p + c.MU * (I0p + c.EPSI * I2p)
    NEG = nlin + c.EPS * (cl['half'] * (D4 * D4)) + c.MU * (cl['mu1t'] * V0V0 + cl['halft'] * V2V2) + c.MUEPSI * (cl['halft'] * (V4 * V4) + cl['mu2t'] * V2V2)
    if K.DELTA_NEG > 0: NEG = NEG * cl['negscale']
    if K.DROPL: L = L * cc(0.0)
    if parts:
        return H, L, NEG, dict(nlin=nlin, alin=alin, mod=mod, forc=forc, cub=cub, strip=strip, R0p=R0p, E2p=c.EPS * R2p, E2n=c.EPS * (cl['half'] * (D4 * D4)),
                              I0p=c.MU * I0p, I0n=c.MU * (cl['MU1'] * V0V0 + cl['half'] * V2V2), I2p=c.MUEPSI * I2p, I2n=c.MUEPSI * (cl['half'] * (V4 * V4) + cl['MU2'] * V2V2))
    return H, L, NEG

def circ(t0, t1):
    """enclosure of (c,s) = ((1-t^2)/(1+t^2), 2t/(1+t^2)) for t in [t0,t1] subset [0,1]; c decreasing, s increasing. abs pad 1e-15."""
    c = lambda t: (1 - t * t) / (1 + t * t); s = lambda t: 2 * t / (1 + t * t)
    PADA = 1e-15
    clo = np.maximum(c(t1) - PADA, 0.0); chi = np.minimum(c(t0) + PADA, 1.0)
    slo = np.maximum(s(t0) - PADA, 0.0); shi = np.minimum(s(t1) + PADA, 1.0)
    return P(clo, chi), P(slo, shi)

def xhat_ranges(lo, hi):
    """NORMALISED physical-variable ranges (A,B,D4,V0,V2,V4)/s (P intervals) for the shape boxes lo,hi (N,5)=(w,t,tI,u,q)"""
    c = consts()
    cw, sw = circ(lo[:, 0], hi[:, 0]); ct, st = circ(lo[:, 1], hi[:, 1]); cI, sI = circ(lo[:, 2], hi[:, 2]); cu, su = circ(lo[:, 3], hi[:, 3])
    q = P(lo[:, 4], hi[:, 4])
    A = cw * ct; B = cw * st * c.inv_sqrt_eps
    V0 = sw * cI * c.inv_sqrt_mu; V2 = sw * sI * c.inv_sqrt_muepsI
    D4 = q * cu * c.inv_sqrt_eps; V4 = q * su * c.inv_sqrt_muepsI
    return A, B, D4, V0, V2, V4

def rate_parts(al, s, xh, parts=False, floors=True):
    """(H/s^2, L/s^2, NEG/s^2) at the POINT scale s (P scalar) for normalised ranges xh (list of 6 P/PG)."""
    xs = [s * x for x in xh]
    H, L, N = rhs_pn(al, *xs, parts=False, floors=floors)
    s2 = s * s
    return H / s2, L / s2, N / s2

def two_scales(al, sl, sh, xh, floors=True):
    """((H,L,N)/s^2 at s_hi, (H,L,N)/s^2 at s_lo); evaluates rhs only once if the cell is a point"""
    ph = rate_parts(al, sh, xh, floors=floors)
    same = (float(np.max(sl.lo)) == float(np.max(sh.lo)))
    pl = ph if same else rate_parts(al, sl, xh, floors=floors)
    return ph, pl

def rate_ub_lo(al, sl, sh, lo, hi, floors=True):
    """rigorous enclosure [rlo, rhi] of rate = tot/Ly over the shape boxes, for all Ly in [s_lo^2, s_hi^2], alpha in al.
    rate = (H-NEG)/s^2 + L/s^2 ; (H-NEG)/s^2 is nondecreasing in s, L/s^2 nonincreasing."""
    xh = xhat_ranges(lo, hi)
    (Hh, Lh, Nh), (Hl, Ll, Nl) = two_scales(al, sl, sh, xh, floors)
    _, ub = sub_signed(Hh + Ll, Nh)          # upper end uses H@s_hi, L@s_lo
    lbv, _ = sub_signed(Hl + Lh, Nl)         # lower end uses H@s_lo, L@s_hi
    return lbv, ub

def nvec():
    """lower-rounded diagonal coefficients n_i of NEG = sum n_i x_i^2, x = (A,B,D4,V0,V2,V4)"""
    k = K
    n = np.array([1 / 32 - k.PI3, 1 / 9, k.EPS / 2, k.MU * k.MU1 * (1 - k.THETA), (1 - k.THETA) * (k.MU / 2 + k.MU * k.EPSI * k.MU2), k.MU * k.EPSI / 2 * (1 - k.THETA)])
    return n * k.NEGSCALE * (1 - 1e-13)
