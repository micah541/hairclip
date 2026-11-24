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
        self.KR16 = arb2P(k.KR16_arb); self.K = {r: (self.KR16 if (v is None or v == 'KR16') else _c(v)) for r, v in k.KREG.items()}
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
        self.ramp = {lev: [(_c(c0), _c(dw)) for (c0, dw) in tab] for lev, tab in k.RAMP_TAB.items()}
        if getattr(k, 'V4B', False):      # E2 fix: recertified nonlinear constants at the asserted boxes (H-E2a..d) and the certified Gaarding r143 (E2-3, E2-4)
            for n_ in ('nd1', 'nd2', 'nd3', 'nd4', 'nd5', 'nd6', 'nd7', 'nd8', 'nf1', 'nf2', 'nf3', 'nf4', 'nf5', 'nf6', 'nf7'): self.cl[n_] = _c(k.E0B[n_.upper()])
            self.cl['r143'] = _c(k.E0B['R143'])
        self.v4 = None; self.e0b = None
        if getattr(k, 'V4B', False):      # fin13v4b: E0 constants of E0_LEMMA_v4b (values already rounded up by e0_lemma_consts.py; _c widens 1 ulp more)
            self.e0b = {n: _c(v) for n, v in k.E0B.items()}; self.e0b['one'] = _c(1.0)
        if k.NEW:      # fin13v4: every v4 constant as an outward interval (values already rounded up by c5_constants.py; _c widens 1 ulp more)
            self.v4 = {n: _c(v) for n, v in k.V4C.items()}
            for n, v in dict(two=2.0, three=3.0, b752=7.52, s1v54=_A(k.S1V) / _A('54.75'), z35_48=arb(35) / 48, z35_96=arb(35) / 96, z35_48s0=arb(35) / 48 * _A(k.S0V),
                             s1v32=_A(k.S1V) / 32, s0v=k.S0V, twoKWIR=2 * _A(k.V4C['KWIR']), boxS=k.BOX_STRIP, boxS2=2 * _A(k.BOX_STRIP)).items():
                self.v4[n] = _c(v)
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

# =====================================================================================================================================
# fin13v4: the E_I inequality of EI_SPEC.md §5 (level 0: §5.1, level 2: §5.2) and the shared L list (SF-3).
# Every coefficient is an outward interval from constants_v4.json; every state quantity enters with a nonnegative coefficient.
def Fn_gen(al, k):
    return k['FN0'] * al + k['FN1'] * (al * al) + k['C3'] * (al * al * al)
def L_list(al, k, mu, mueps, eps):
    """The L group (x-degree <= 1, every term carries a factor alpha) as a list of (coefficient, monomial); monomials in
    {A, B, D4, W3, V0, V1, V2, V3, V4}.  Used verbatim by rhs_pn (P intervals) and by unc.lk (arb, linear majorant): one list, no hand copy.
    k: dict of constants (P or arb), with keys of constants_v4.json and 'f572', 't22', 't12'."""
    fn = Fn_gen(al, k)
    fnw = k['FNW0'] * al + k['FNW1'] * (al * al) + k['FNW2'] * (al * al * al)          # ||Phi_F||_W (level 0, IBP form)
    fnwd = k['FNWd0'] * al + k['FNWd1'] * (al * al) + k['FNWd2'] * (al * al * al)      # ||zeta Phi_F'||_W (level 2)
    ABc = k['PHIAB'] * fn; BTc = k['BETAF'] * fn
    np3e = k['NP3'] + k['two'] * al * k['NH2']
    c0 = fnw * k['KAP'] + ABc * al * k['NAB'] + BTc * np3e + ABc * k['NP12']          # coefficient of V0 at level 0 (with mu)
    c2 = fnwd + ABc * al * k['NAB'] + BTc * np3e + ABc * k['NP12']                    # coefficient of DV at level 2 (with mu eps_I)
    fe = eps * fn
    if 'N3E2' in k:     # fin13v4b E2 fix (F3): Fn parts of beta and (a,b) at the E2 level: FE2 = 1 + BETAF (n3 + 2 al ||h2||) + PHIAB (NF12 + al (GABN + al (RA2 + RB2)))
        fe = eps * fn * (k['one'] + k['BETAF'] * (k['N3E2'] + k['two'] * al * k['NH2R']) + k['PHIAB'] * (k['NF12'] + al * (k['GABN'] + al * (k['RA2'] + k['RB2']))))
    out = [(k['f572'] * fn, 'A'), (fe, 'D4'), (fe * k['t22'], 'W3'), (fe * k['t12'], 'B'),      # E0 forcing, E2 forcing (fin12f4 forms; v4b: times FE2)
           (mu * fnw, 'V1'), (mu * c0, 'V0'),
           (mueps * c2, 'V4'), (mueps * c2 * k['two_PSI1'], 'V3'), (mueps * c2 * k['PSI2'], 'V2')]
    return out
def _ei_v4(al, A, B, D4, V0, V2, V4, W1, W3, V1, V3, mid, regS, sr, sI1r, sI2r, _q, c, FL, fn, e2ctx=None):
    v = c.v4; cl = c.cl
    lev = K.v4_level(float(np.max(al.hi)))
    lk_ = K.lev_key(lev)
    T = lambda z, key: v[f'TQ[{lk_}][{z}][{key}]']
    TF = lambda key: mx(T('ZA', key), T('ZB', key))           # far: sup over [16.25,60] (contains Z1 and OR)
    TG = lambda key: mx(TF(key), T('IR', key))                # global: sup over [10.25,60]
    Ti = lambda key: T('IR', key)
    HD = lambda key: v[f'HD[{lk_}][{key}]']
    QUF = v[f'QUF[{lk_}]']
    AA = A * A; BB = B * B; V0V0 = V0 * V0; V2V2 = V2 * V2; W1W1 = W1 * W1
    # ---- derived sups (EI_SPEC 1.2)
    S1r = mx(mid[1], regS["12<x<=16.25"][1]); S2r = mx(mid[2], regS["12<x<=16.25"][2])      # sup over IR of |w'|, |w''|
    sZ = mx(sI2r, v['s1v32'])                                                                 # sup over Z1 of |g'|
    sOR = v['s1v54'] + v['z35_48s0']                                                          # sup over OR of |v'|
    sg1 = mx(sZ, S2r)                                                                         # sup over [10.25,60] of |g'|
    sG = mx(mx(sZ, sOR), S2r + v['z35_96'] * S1r)                                             # sup over [10.25,60] of |v'|
    S1z = gm(V1, v['twoKWIR'] * V2) + v['z35_96'] * S1r                                       # sup over IR of |zeta g'|
    S2z = gm(V2, v['twoKWIR'] * V3) + v['two'] * v['z35_96'] * S2r + v['Z2i'] * S1r           # sup over IR of |zeta g''|
    sw16 = sr[0]                                                                              # |w(16.25)|
    S0s = sr[0] + c.sqrtPHI34 * V0                                                            # sup|w| on [16.25,34.25]
    supv2 = gm(V2, c.KW * (c.PSI1W * V2 + v['two'] * V3))                                     # sup|v''| = sup|w'''| on [16.25,34.25]
    S_all34 = None
    for reg in K.KREG: S_all34 = regS[reg][0] if S_all34 is None else mx(S_all34, regS[reg][0])   # sup|w| on [0,34.25]
    DV = V4 + v['two_PSI1'] * V3 + v['PSI2'] * V2
    # ---- coefficients
    CBAR = TG('c'); C1 = TG('c_a') + TG('c_b') * sg1
    RQ = Ti('c') * (v['k4_1'] * D4 + v['k4_2'] * W3 + v['k4_3'] * B + v['k4_4'] * W1) + (Ti('c_a') + Ti('c_b') * S2r) * (v['k3_1'] * W3 + v['k3_2'] * B + v['k3_3'] * W1)
    C0W = v['C0P'] * S_all34 + v['C0Hb'] * V0
    # ---- modulation (EI_SPEC 4.9; U-3 state form of the strip parts, sups on sigma in [29,32] subset x in [16.25,34.25])
    Q = _q * _q
    ABs = (v['SF_ab_0'] + v['SF_ab_al'] * al) * S0s + v['SF_ab_1'] * sI1r + v['SF_ab_2'] * sI2r + v['SF_ab_3'] * supv2
    BTs = (v['SF_beta_0'] + v['SF_beta_al'] * al) * S0s + v['SF_beta_1'] * sI1r + v['SF_beta_2'] * sI2r + v['SF_beta_3'] * supv2
    ABx = (v['CAB0'] + v['CAB0_dlead'] + v['MAB'] * al) * A + v['QAB'] * Q + ABs
    BTx = (v['BETA0'] + v['BETA0_dlead'] + v['BETA1'] * al) * A + v['QB'] * Q + BTs
    ABc = v['PHIAB'] * fn
    AB = ABx + ABc
    np3e = v['NP3'] + v['two'] * al * v['NH2']
    EAt = TG('EA0') * V0 + v['RWR'] * Ti('EA1') * A + HD('hEA1') * sw16 + HD('HEA1') * V0
    # ================= level 0 (EI_SPEC 5.1) =================
    I0 = (v['AF'] * W1W1 + v['BF'] * BB                                                                         # 4.2
          + CBAR * V2V2 + (C1 + v['two_PSI1'] * CBAR) * (V1 * V2) + (v['PSI1'] * C1 + v['PSI2'] * CBAR) * (V0 * V2)  # 4.3
          + V0 * RQ                                                                                             # 4.4
          + (V1 + v['KAP_Z1'] * V0) * (TF('d010') * V2 + TF('d100') * V1 + TF('A1') * V0 + TF('d200') * sZ * V1 + TF('d110') * sZ * V2
                                       + TF('d300') * (sZ * sZ) * V1 + HD('hA0') * sw16 + HD('HA0') * V0)        # 4.5
          + (V1 + v['KAP_IR'] * V0) * v['RWR'] * (Ti('d010') * W3 + Ti('d100') * B + Ti('A1') * W1 + Ti('A0') * A + Ti('d200') * S2r * B
                                                  + Ti('d110') * S2r * W3 + Ti('d300') * (S2r * S2r) * B)       # 4.6
          + V0 * C0W                                                                                            # 4.7
          + AB * v['T0'] * V0V0 + AB * v['TC0'] * QUF * W1W1 + ABx * al * v['NAB'] * V0                         # 4.9 T8, T9, T10
          + BTx * np3e * V0 + ABx * v['NP12'] * V0                                                              # 4.9 T11
          + AB * (TG('DQ0') * (V0 * V1) + TG('DQ1') * V0V0) + AB * V0 * EAt)                                    # 4.9 T12
    # ================= level 2 (EI_SPEC 5.2) =================
    NZ = (TF('d010') * V3 + (TF('d010_a') + TF('d100')) * V2 + TF('d010_b') * sZ * V2 + TF('d100_a') * V1 + TF('d100_b') * sZ * V1
          + v['two'] * TF('d200') * sZ * V2 + TF('d200_a') * sZ * V1 + TF('d200_b') * (sZ * sZ) * V1
          + TF('d110') * sG * (v['three'] * V3 + v['PSI1'] * V2) + TF('d110') * sZ * V3 + TF('d110_a') * sZ * V2 + TF('d110_b') * (sZ * sZ) * V2
          + v['three'] * TF('d300') * (sZ * sZ) * V2 + TF('d300_a') * (sZ * sZ) * V1 + TF('d300_b') * (sZ * sZ * sZ) * V1
          + TF('A1') * V1 + TF('B1') * V0 + HD('hB0') * sw16 + HD('HB0') * V0)
    R5 = (Ti('d010') * D4 + (Ti('d010_a') + Ti('d100')) * W3 + Ti('d010_b') * S2r * W3 + Ti('d100_a') * B + Ti('d100_b') * S2r * B
          + v['two'] * Ti('d200') * S2r * W3 + Ti('d200_a') * S2r * B + Ti('d200_b') * (S2r * S2r) * B
          + Ti('d110') * (S1z * D4 + S2z * W3) + Ti('d110_a') * S2r * W3 + Ti('d110_b') * (S2r * S2r) * W3
          + v['three'] * Ti('d300') * (S2r * S2r) * W3 + Ti('d300_a') * (S2r * S2r) * B + Ti('d300_b') * (S2r * S2r * S2r) * B
          + Ti('A1') * B + Ti('B1') * W1 + Ti('B0') * A)
    I2 = (v['E1'] * (V1 * V2) + v['E0'] * (V0 * V2)                                                            # 4.1
          + (v['CF0'] * D4 + v['CF1'] * W3 + v['CF2'] * B + v['CF3'] * W1) * DV                                 # 4.2
          + CBAR * V4 * DV + C1 * V3 * DV                                                                       # 4.3
          + RQ * DV                                                                                             # 4.4
          + NZ * DV                                                                                             # 4.5
          + v['RWR'] * R5 * DV                                                                                  # 4.6
          + C0W * DV                                                                                            # 4.7
          + AB * (v['T2'] * V2V2 + v['two'] * v['Q2B'] * (V1 * V2) + v['Q2B'] * V0 * (V3 + v['PSI1'] * V2))      # 4.9 T8
          + AB * v['TC2'] * QUF * W1 * DV + ABx * al * v['NAB'] * DV                                            # T9, T10
          + BTx * np3e * DV + ABx * v['NP12'] * DV                                                              # T11
          + AB * (TG('DQ0') * V1 + TG('DQ1') * V0) * DV + AB * EAt * DV)                                        # T12
    # ================= E0: strip part of beta in state form (SF-1): n3 |beta_strip| A, n3 = ||chi psi3^|| = 7.52 =================
    if c.e0b is None:
        R0add = v['b752'] * BTs * A
    else:
        # fin13v4b (E0_LEMMA_v4b.tex): -beta_S <g0, w~> = -beta_S <Pi g0, w~> (w~ orthogonal to f1, f2, f3): coupling n3_123 = ||Pi g0|| (not 7.52)
        e = c.e0b
        R0add = e['N3P'] * BTs * A
        XW = e['XW_COEF'] * S0s                     # ||chi' w||_rho <= sup|chi'| (int_strip rho)^1/2 sup_strip |w|
        W1X = W1 + XW                               # ||chi w'||_rho <= ||(chi w)'|| + ||chi' w||
        V1V1 = V1 * V1
        R0add = R0add + al * (e['M_STRIP_00'] * e['IR'] * (S0s * S0s) + (e['M_STRIP_01'] * V0 + e['M_STRIP_02'] * V1) * e['IRW'] * S0s
                              + e['RW'] * (e['M_STRIP_11'] * V0V0 + e['M_STRIP_12'] * (V0 * V1) + e['M_STRIP_22'] * V1V1))       # (E0-8) strip alpha-linear form
        R0add = R0add + S_all34 * (e['CD1'] * A * W1X + e['CD0'] * AA)                                                              # cubic drift
        R0add = R0add + AB * ((e['one'] + e['CA0'] + e['CB0']) * AA + (e['CA1'] + e['CB1']) * A * W1X)                              # (a,b) w-parts
        R0add = R0add + al * (e['GABN'] + e['LOOPC']) * (v['QAB'] * Q + ABs) * A                                                     # (a,b) beyond lead: Q1, loop
        R0add = R0add + ABx * (al * al) * (e['RA2'] + e['RB2']) * A                                                                  # alpha^2 remainders of A1, B1 in Q1
        R0add = R0add + al * e['HOB'] * (v['BETA1'] * al * A + v['QB'] * Q + BTs) * A                                               # beta beyond lead: h2, E33
        R0add = R0add + al * (e['H12AB'] * ABx + e['H12B'] * BTx) * A                                                                # H31, H32 part of beta_E
    # ================= fin13v4b E2 fix (E2_AUDIT.md section 5, F2): modulation beyond the Q2 lead + strip part of chi N(w), times Dq (inside eps*[...])
    E2add = None
    if c.e0b is not None and e2ctx is not None:
        e = c.e0b; Dq, supw2f, wstrip = e2ctx
        RABx = (e['DAB'] + v['CAB0_dlead'] + v['MAB'] * al) * A + v['QAB'] * Q + ABs                    # |(a,b) - Minv l|, x-part
        n3e = e['N3E2'] + v['two'] * al * e['NH2R']
        wprt = (e['one'] + e['CA0'] + e['CB0']) * A + (e['CA1'] + e['CB1']) * (W1 + e['XW_COEF'] * S0s)  # ||chi (A1,B1)(u) - chi (A1,B1)(u0)||
        u0pr = al * (e['GABN'] + al * (e['RA2'] + e['RB2']))                                          # ||chi (A1,B1)(u0)||
        RM = RABx * e['NF12'] + BTx * n3e + ABx * (u0pr + wprt) + ABc * wprt
        NS = e['SRW'] * (cl['nf1'] * V3 + cl['nf2'] * V2 + cl['nf3'] * V1 + cl['nf7'] * (V0 + wstrip) + supw2f * (cl['nf4'] * V1 + cl['nf5'] * V2) + cl['nf6'] * (supw2f * supw2f) * V1)
        E2add = (RM + NS) * Dq
    # ================= L (shared list) + floors =================
    kd = dict(v); kd['f572'] = cl['f572']; kd['t22'] = cl['t22']; kd['t12'] = cl['t12']
    if c.e0b is not None:
        for n_ in ('N3E2', 'NH2R', 'NF12', 'GABN', 'RA2', 'RB2', 'one'): kd[n_] = c.e0b[n_]
    X = dict(A=A, B=B, D4=D4, W3=W3, V0=V0, V1=V1, V2=V2, V3=V3, V4=V4)
    L = None
    for cf, mono in L_list(al, kd, c.MU, c.MUEPSI, c.EPS):
        t = cf * X[mono]; L = t if L is None else L + t
    if FL and 'c' in K.ADD and K.CMODE == "abs":
        L = L + c.MU * cl['aflux0a'] + c.MUEPSI * cl['aflux2a']         # FLUX0, FLUX2 (inputs; K^OR of EI_SPEC 4.10)
    return I0, I2, R0add, L, E2add

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
    sI1r, sI2r = sI1, sI2          # raw E_I sups on [16.25,34.25] (x-valuation 1), before the max with S0V, S1V/32
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
    if getattr(K, 'V4B', False): strip = cl['SB_SA'] * V0V0 + cl['twoIA'] * c.KR16 * (c.PHX * AA + cl['two'] * AW1)      # fin13v4b: (I2) has Phi1 (E0-6)
    else: strip = cl['SB_SA'] * V0V0 + cl['twoIA'] * c.KR16 * (AA + cl['two'] * AW1)
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
    if K.NEW:
        I0p, I2p, R0add, Lv4, E2add = _ei_v4(al, A, B, D4, V0, V2, V4, W1, W3, V1, V3, mid, regS, sr, sI1r, sI2r, _q, c, FL, fn, e2ctx=(Dq, supw2f, wstrip))
        R0p = R0p + R0add
        if E2add is not None: R2p = R2p + E2add
        H = R0p + c.EPS * R2p + c.MU * (I0p + c.EPSI * I2p)
        L = Lv4
        NEG = nlin + c.EPS * (cl['half'] * (D4 * D4)) + c.MU * (cl['mu1t'] * V0V0 + cl['halft'] * V2V2) + c.MUEPSI * (cl['halft'] * (V4 * V4) + cl['mu2t'] * V2V2)
        if K.DELTA_NEG > 0: NEG = NEG * cl['negscale']
        if K.DROPL: L = L * cc(0.0)
        if parts:
            return H, L, NEG, dict(nlin=nlin, alin=alin, mod=mod, forc=forc, cub=cub, strip=strip, R0p=R0p, E2p=c.EPS * R2p, I0p=c.MU * I0p, I2p=c.MUEPSI * I2p, R0add=R0add)
        return H, L, NEG
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
