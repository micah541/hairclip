"""SINGLE SOURCE of all constants of the rigorous verifier.  Values copied verbatim from ../inner_closure_float2.py
(CUBFINE=1, ALIN=spec, BLOOP=0, CAB=7.57 configuration); provenance tags as in ../ASSEMBLY.md.
Change a constant HERE and rerun (closure.py / run_all.sh).  check_vs_float.py asserts that the interval rhs
reproduces ../inner_closure_float2.py (plus the explicitly listed ADDED terms, switchable by ADDED=0)."""
import json, os
from flint import arb, ctx
ctx.prec = 200
SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../step3nl  (read only)
PARAMS = dict(EPS=1e-5, MU=0.08, EPSI=3e-5)
for k in PARAMS: PARAMS[k] = float(os.environ.get(k, PARAMS[k]))
EPS, MU, EPSI = PARAMS["EPS"], PARAMS["MU"], PARAMS["EPSI"]
ADD = set(os.environ.get("ADD", "abcd"))        # which of the terms missing in the float rhs are included: a=ramp coupling (I0), b=c0' w coupling (I0), c=flux constants (I0,I2), d=strip cubic. ADD="" reproduces the float rhs
c1 = 0.5804; PHX = 1.083; K2 = 2.078; K12 = 992.4
LAMT = -0.2561532; PI3 = 0.01607
MOD_ALPHA = 6.19; MOD_E2 = 43.33; C3 = 12.104; CAB = float(os.environ.get("CAB", 7.57))
ALIN = dict(sq=59.3, AB=35.9, W1B=7.9, BB=6.2)     # alin = al*(sq (A^2+W1^2) + AB A B + W1B W1 B + BB B^2)
KREG = {"0<x<=1": 1.488, "1<x<=2": 2.078, "2<x<=4": 5.222, "4<x<=8": 56.92, "8<x<=12": 992.4, "12<x<=16.25": 3.0208e4, "16.25<x<=31.25": None}
cubf = json.load(open(os.path.join(SRC, "cubic_E0fine_0.001_0.01_0.01.json")))["sup"]
MU1, MU2 = 0.1552, 0.6555
KW = 2.323e6; PSI1W = 1.305; C1W = 0.8477
PHI34 = 9.8547e6
SB, SA, IA = 3.5955e-5, 1.3372e-4, 9.5997e-12
AF, BF = 6.9555e-2, 3.9359e-2
CF = (0.4421, 0.5259, 0.4235, 0.4250)
DS = (6.841e-3, 1.5418e-2, 2.3648e-2, 4.5627e-2)
# recertified far-field T3 maxima for slope box |w'|<=0.02 (agent_int/planB/r60_T3_boxw1_0.02.json, max over regions 16.25..60); valid only if sup|w'|(x>=16.25)<=0.02 (checked in closure.py: far_max)
# previous (box 1e-2): d001=1.85e-2, d010=1.5e-3, d200=1.5e-2, d110=0.98, d300=3.26.   d100 (1.1e-4, includes the ramp), A0, A1 unchanged.
TQ = dict(d001=2.27e-2, d010=1.67e-3, d100=1.1e-4, d200=1.53e-2, d110=1.086, d300=3.314, A0=1.7e-5, A1=7.1e-5)
QPSI = 0.273 + 0.0159
FN = (7.93e-7, 2.769e-6, C3)         # Fn(al)  = FN[0] al + FN[1] al^2 + FN[2] al^3
FNW = (1.09e-7, 6.5e-9, 0.1823)      # FnW(al) = ...
# slow mode / modulation beta bound (agent_mod/lemma): beta = (BETA0 + BETA1|al|) sqrt(Ly) + BETAF Fn(al)
BETA0, BETA1, BETAF = 0.00244, 63.5, 0.626
# ---- ADDED TERMS (not in the float rhs; ASSEMBLY.md sec. 5, 6 and task statement) ----
ADD_I0_RAMP = (float(os.environ.get('AR0', 0.66)), float(os.environ.get('AR1', 1.85e-2)), 1.5e-3, 1.1e-4, 7.1e-5, 1.7e-5)     # 0.66*(d001 D4 + d010 W3 + d100 B + A1 W1 + A0 A)*(V1+1.67 V0)  in I0
ADD_I0_C0W = 9.9e-12 * 1e-2                                     # c0' w coupling: 9.9e-12*1e-2*V0 in I0
ADD_FLUX_I0 = 1e-17; ADD_FLUX_I2 = 1e-17                        # (old placeholders; FIX=1 uses ABS_FLUX below)                        # outer ramp flux constants (I0 and I2)
ADD_STRIP_CUBIC = 2.6e-3                                        # cubic on strip [31.25,34.25] = 2.6e-3 * (16.25<x<=31.25 region term)
# ---- certified inner-ramp deviation coefficients (agent_int/planB/planB_ramp.json 'table'), replace the fixed ADD_I0_RAMP[1:] when RAMP=table (default)
#      coefficient(alpha<=lev, S) = c0(lev) + dw(lev) * S ,  c0 = value at w-box eta=0, dw = (eta=1e-2 value - eta=0 value)/1e-2 , S = sup(|w|,|w'|) on [10.25,16.25]
RAMP = os.environ.get("RAMP", "table")
RAMP_NAMES = ("d001_w4", "d010_w3", "d100_w2", "A1_w1", "A0_w0")       # multiply D4, W3, B, W1, A  (derivative order 4,3,2,1,0 of w)
RAMPMULT = float(os.environ.get("RAMPMULT", 1.0))      # sensitivity only (default 1)
RAMP_TAB = {}
if RAMP == "table":
    _rt = json.load(open(os.path.join(SRC, "agent_int", "planB", "planB_ramp.json")))["table"]
    for _am in (2.5e-4, 5e-4, 1e-3):
        RAMP_TAB[_am] = [(RAMPMULT * _rt[f"alpha={_am:g},eta=0"][n], RAMPMULT * max(_rt[f"alpha={_am:g},eta=0.01"][n] - _rt[f"alpha={_am:g},eta=0"][n], 0.0) / 1e-2) for n in RAMP_NAMES]
def ramp_level(a_hi):
    for lev in sorted(RAMP_TAB):
        if a_hi <= lev * (1 + 1e-12): return lev
    raise ValueError(f"alpha {a_hi:g} above the certified ramp table (1e-3)")
THETA = 0.0       # fraction of the V-dissipation moved to the left side (time-integral bounds); set by rhs_iv.set_theta
# ---- state-dependent recasting of (b),(c) (QUAD=1, default): (b) 9.9e-12 * sup|w|(all regions, from shape variables) * V0 in I0 (the 1e-2 box placeholder removed);
#      (c) CFLUX*(V0^2+V1^2+V2^2) in I0 and CFLUX*(V2^2+V3^2+V4^2) in I2 (replaces the constant fluxes 1e-17).  CFLUX is an INPUT (not derived here); default = the old placeholder magnitude.
QUAD = os.environ.get("QUAD", "1") == "1"
C0P_NORM = 9.9e-12
CFLUX = float(os.environ.get("CFLUX", 1e-17))
# ---- FIX PLAN (STEP3_FIXPLAN.md): FIX=1 (default) switches on F2, F3, F4; FIX=0 reproduces the 'final pass' configuration
FIX = os.environ.get("FIX", "1") == "1"
S0V = float(os.environ.get("S0V", 0.011)); S1V = float(os.environ.get("S1V", 0.87))     # F2 design inputs (Lemma S will certify): sup|w'|, sup y|w''| on [34.25,60]
GN = FIX and os.environ.get("GN", "1") == "1"                                           # F2 identity ||v''^2||_W <= sup|v'| (3 V3 + 1.31 V2)
QB = float(os.environ.get("QB", 50)); QM = float(os.environ.get("QM", 30))              # F3 placeholders q_beta, q_(a,b)  (swap in agent_mod/ibp/mod_ibp.json when it exists)
CMODE = os.environ.get("CMODE", "abs" if FIX else "quad")                                # F4: term (c) as absolute forcing (abs) or quadratic placeholder (quad)
ABS_FLUX_I0 = float(os.environ.get('FLUX0', 4.1e-21)); ABS_FLUX_I2 = float(os.environ.get('FLUX2', 2.0e-20))   # certified lemmaS_consts.json K0<=4.1e-21, K2<=2.0e-20 per unit tau                                            # lemmaS_cflux.json abs_box K0_abs=6.093e-18, K2_abs=1.544e-17 (placeholders until Lemma S)
CLIN = float(os.environ.get('CLIN', 0.0))      # F3: absolute linear forcing c_j*A coefficient (pending mod_ibp.json)
DELTA_NEG = 0.0; NEGSCALE = 1.0; DROPL = False                                                            # V1 pass: NEG -> (1-delta) NEG, L dropped (set by rhs_iv.set_mode)
# F3 certified (agent_mod/ibp/mod_ibp.json, alpha1=5e-4, eta=0.02, valid for A,W1<=1e-3, |alpha|<=5e-4): MOD=ibp (default with FIX)
MODIBP = FIX and os.environ.get("MOD", "ibp") == "ibp"
CAB0, MAB, QAB, PHIAB = CAB, 0.0, 0.0, 0.0
if MODIBP:
    BETA1 = 66.53; BETAF = 0.657; QB = 5.53; MOD_ALPHA = float(os.environ.get('MODA', 6.19 * 66.5 / 63.5));     # |beta| <= (0.00244+66.53|a|)A + 5.53 Q + 0.657||Phi_F||  (v3.1: BETA1 was 66.5 < certified 66.5238); E0 alpha-form scaled with BETA1 (assumption)
    CAB0, MAB, QAB, PHIAB = 3.7083, 622.4, 51.34, 3.8404                        # |(a,b)| <= (3.7083+622.4|a|)A + 51.34 Q + 3.8404||Phi_F||  (v3.1: outward-rounded final_5e-4.json core.ab: 3.708257, 622.378, 51.3397, 3.84031; were 3.708/51.3/3.84)  (+strip terms <=1e-8/unit S: negligible, see VERIFIER.md)
# strip terms (mod_ibp.json 'strip', sums of all |coefficients| per unit S, rounded up) times the strip hypothesis S <= 0.02 (|w|,|w'|,.. on sigma in [29,32])
STRIP_S = 0.02; STRIP_ON = MODIBP and os.environ.get("STRIP", "1") == "1"
STRIP_BETA = STRIP_S * 5.3e-8 if STRIP_ON else 0.0          # beta: 5.256e-8 per unit S
STRIP_AB = STRIP_S * 2.76e-7 if STRIP_ON else 0.0           # (a,b) Euclidean: hypot(1.2923e-7, 2.4303e-7) = 2.7525e-7 per unit S (v3.1; was 2.5e-7 = max, too small for the Euclidean norm)
ETA_ANCHOR = float(os.environ.get('ETA', 1e-9))     # relative width of the tail anchor cell Ly in [L*, (1+eta)L*]
# ---- derived rigorous constants (arb balls) ----
KR16_arb = (((1 + arb(16.25) ** 2) ** (arb(2) / 3)) / 4).exp()  # 1/rho(16.25)
W3125_arb = (1 + arb(31.25) ** 2) ** -2
# =====================================================================================================================================
# ---- fin13v4 (verifier v4, EI_SPEC.md): NEW=1 (default) takes every certified constant below from ../v4cert/out/constants_v4.json
#      (raw certified values rounded outward by v4cert/c5_constants.py); NEW=0 keeps the fin12f4 typed values above (regression only).
NEW = os.environ.get("NEW", "1") == "1"
V4C = None
if NEW:
    _cj = json.load(open(os.path.join(SRC, "v4cert", "out", "constants_v4.json")))["constants"]
    V4C = {k: v["value"] for k, v in _cj.items()}
    def _v(n): return V4C[n]
    SB, SA, IA = _v('SB'), _v('SA'), _v('IA')
    AF, BF = _v('AF'), _v('BF')
    CF = (_v('CF0'), _v('CF1'), _v('CF2'), _v('CF3'))
    DS = (_v('DS0'), _v('DS1'), _v('DS2'), _v('DS3'))
    C3 = _v('C3'); FN = (_v('FN0'), _v('FN1'), C3)
    FNW = (_v('FNW0'), _v('FNW1'), _v('FNW2'))           # ||Phi_F||_W, level 0 (C^4 functions)
    FNWD = (_v('FNWd0'), _v('FNWd1'), _v('FNWd2'))       # ||zeta Phi_F'||_W, level 2 (C^4 functions)
    MU1, MU2 = _v('MU1'), _v('MU2')                      # lower bounds (rounded down)
    PHI34 = min(_v('PHI34'), _v('PHI34_c4'))          # both certified upper bounds of Phi(34.25) (coupling: tight; c4_hardy: coarse upper sum)
    c1 = _v('c1'); PHX = _v('PHX'); K2 = _v('K2'); K12 = _v('K12'); KW = _v('KW'); PSI1W = _v('PSI1W'); C1W = _v('C1W')
    KREG = {"0<x<=1": _v('KREG[0,1]'), "1<x<=2": _v('KREG[1,2]'), "2<x<=4": _v('KREG[2,4]'), "4<x<=8": _v('KREG[4,8]'), "8<x<=12": _v('KREG[8,12]'),
            "12<x<=16.25": "KR16", "16.25<x<=31.25": None}       # 'KR16': the arb 1/rho(16.25) (typed 3.0208e4 was below the raw 30208.006)
    if MODIBP:
        BETA0 = _v('BETA0') + _v('BETA0_dlead'); BETA1 = _v('BETA1'); BETAF = _v('BETAF'); QB = _v('QB')
        CAB0 = _v('CAB0') + _v('CAB0_dlead'); MAB = _v('MAB'); QAB = _v('QAB'); PHIAB = _v('PHIAB')
        STRIP_BETA = _v('STRIP_BETA') if STRIP_ON else 0.0       # U-3 linear form at S_r <= 0.02, |alpha| <= 5e-4, with Neumann feedback (v4cert/u3_strip_form.py)
        STRIP_AB = _v('STRIP_AB') if STRIP_ON else 0.0
    C0P_NORM = _v('C0P')
    AMAX_LEVELS = (2.5e-4, 1e-3)
    def lev_key(lev): return {2.5e-4: '0.00025', 1e-3: '0.001'}[lev]
    def TQv(lev, zone, key): return V4C[f'TQ[{lev_key(lev)}][{zone}][{key}]']
    def HDv(lev, key): return V4C[f'HD[{lev_key(lev)}][{key}]']
    def v4_level(a_hi):
        for lev in AMAX_LEVELS:
            if a_hi <= lev * (1 + 1e-12): return lev
        raise ValueError(f"alpha {a_hi:g} above the certified v4 coefficient tables (1e-3)")
    # standing box hypotheses (EI_SPEC 1.3), asserted in closure.py at every step and at the tail anchor
    ETA_IR = 1e-3          # B-b: S0r, S1r <= ETA_IR (IR jet box of the C2 tables)
    BOX_W1 = 0.02          # B-c: |w'| <= 0.02 on [16.25,60] (far tables)
    BOX_W0A = 0.02; BOX_W0B = 0.5      # B-d: |w| <= 0.02 on [16.25,26], <= 0.5 on [26,60]
    BOX_STRIP = 0.02       # U-3: S0, S1, S2 <= 0.02 on the strip sigma in [29,32]
