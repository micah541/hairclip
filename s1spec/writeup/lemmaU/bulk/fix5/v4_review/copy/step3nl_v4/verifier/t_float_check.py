#!/usr/bin/env python3
"""T(iii) (fin13v4): point test of the v4 E_I terms.  An independent plain-float evaluation of EI_SPEC.md §5.1 (level 0), §5.2 (level 2),
the E0 strip term (SF-1) and the L group, written from the spec text (not from rhs_iv._ei_v4), is compared with the interval rhs_pn(parts=True)
at random points (degenerate boxes): the float value must lie in the interval (up to 1e-12 relative for the float's own rounding).
Run in the verifier environment: NEW=1 PYTHONHASHSEED=0 python3 t_float_check.py  -> t_float_check.json"""
import os, sys, json, math, time, hashlib
import numpy as np
os.environ.setdefault("NEW", "1")
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import constants as K, ival, rhs_iv as R
assert K.NEW
V = {k: float(v) for k, v in K.V4C.items()}
sq = math.sqrt
KR16 = float(K.KR16_arb.mid()); KREGf = {r: (KR16 if v in (None, 'KR16') else float(v)) for r, v in K.KREG.items()}
def spec_eval(al, A, B, D4, V0, V2, V4, lev='0.00025'):
    T = lambda z, k: V[f'TQ[{lev}][{z}][{k}]']
    F = lambda k: max(T('ZA', k), T('ZB', k)); G = lambda k: max(F(k), T('IR', k)); I = lambda k: T('IR', k)
    HD = lambda k: V[f'HD[{lev}][{k}]']
    c1, PHX, C1W, KW, PSI1W = K.c1, K.PHX, K.C1W, K.KW, K.PSI1W
    W1 = sq(A * (B + c1 * A)); W3 = sq(B * (D4 + c1 * B)); V1 = sq(V0 * (V2 + C1W * V0)); V3 = sq(V2 * (V4 + C1W * V2))
    ss = lambda Kc: (sq(A * Kc * (PHX * A + 2 * W1)), sq(W1 * Kc * (PHX * W1 + 2 * B)), sq(B * Kc * (PHX * B + 2 * W3)))
    mid = ss(K.K12); r12 = ss(KR16); sr = ss(KR16)
    S1r = max(mid[1], r12[1]); S2r = max(mid[2], r12[2])
    sI1r = sq(V0 * KW * (PSI1W * V0 + 2 * V1)); sI2r = sq(V1 * KW * (PSI1W * V1 + 2 * V2)); supv2 = sq(V2 * KW * (PSI1W * V2 + 2 * V3))
    sw16 = sr[0]; S0s = sr[0] + sq(K.PHI34) * V0
    swall = [ss(KREGf[r])[0] for r in KREGf if KREGf[r] is not None and r != "16.25<x<=31.25"] + [S0s]
    S_all34 = max(swall)
    sZ = max(sI2r, K.S1V / 32); sOR = K.S1V / 54.75 + 35 / 48 * K.S0V
    sg1 = max(sZ, S2r); sG = max(sI2r, K.S1V / 32, sOR, S2r + 35 / 96 * S1r)
    S1z = sq(2 * V['KWIR'] * V1 * V2) + 35 / 96 * S1r
    S2z = sq(2 * V['KWIR'] * V2 * V3) + 2 * 35 / 96 * S2r + V['Z2i'] * S1r
    PSI1, PSI2 = V['PSI1'], V['PSI2']
    DV = V4 + V['two_PSI1'] * V3 + PSI2 * V2
    CBAR = G('c'); CA = G('c_a'); CB = G('c_b'); C1 = CA + CB * sg1
    k4 = [V[f'k4_{j}'] for j in (1, 2, 3, 4)]; k3 = [V[f'k3_{j}'] for j in (1, 2, 3)]
    RQ = I('c') * (k4[0] * D4 + k4[1] * W3 + k4[2] * B + k4[3] * W1) + (I('c_a') + I('c_b') * S2r) * (k3[0] * W3 + k3[1] * B + k3[2] * W1)
    Fn = V['FN0'] * al + V['FN1'] * al ** 2 + V['C3'] * al ** 3
    Q = (A + W1 + B) ** 2
    st_ab = (V['SF_ab_0'] + V['SF_ab_al'] * al) * S0s + V['SF_ab_1'] * sI1r + V['SF_ab_2'] * sI2r + V['SF_ab_3'] * supv2
    st_b = (V['SF_beta_0'] + V['SF_beta_al'] * al) * S0s + V['SF_beta_1'] * sI1r + V['SF_beta_2'] * sI2r + V['SF_beta_3'] * supv2
    ABx = (V['CAB0'] + V['CAB0_dlead'] + V['MAB'] * al) * A + V['QAB'] * Q + st_ab; ABc = V['PHIAB'] * Fn; AB = ABx + ABc
    BTx = (V['BETA0'] + V['BETA0_dlead'] + V['BETA1'] * al) * A + V['QB'] * Q + st_b; BTc = V['BETAF'] * Fn
    C0W = V['C0P'] * S_all34 + V['C0Hb'] * V0
    np3 = V['NP3'] + 2 * al * V['NH2']; QUF = V[f'QUF[{lev}]']
    T12 = G('EA0') * V0 + V['RWR'] * I('EA1') * A + HD('hEA1') * sw16 + HD('HEA1') * V0
    I0 = (V['AF'] * W1 ** 2 + V['BF'] * B ** 2
          + CBAR * V2 ** 2 + (C1 + 2 * PSI1 * CBAR) * V1 * V2 + (PSI1 * C1 + PSI2 * CBAR) * V0 * V2
          + V0 * RQ
          + (V1 + V['KAP_Z1'] * V0) * (F('d010') * V2 + F('d100') * V1 + F('A1') * V0 + F('d200') * sZ * V1 + F('d110') * sZ * V2 + F('d300') * sZ ** 2 * V1 + HD('hA0') * sw16 + HD('HA0') * V0)
          + (V1 + V['KAP_IR'] * V0) * V['RWR'] * (I('d010') * W3 + I('d100') * B + I('A1') * W1 + I('A0') * A + I('d200') * S2r * B + I('d110') * S2r * W3 + I('d300') * S2r ** 2 * B)
          + V0 * C0W
          + AB * V['T0'] * V0 ** 2 + AB * V['TC0'] * QUF * W1 ** 2 + ABx * al * V['NAB'] * V0
          + BTx * np3 * V0 + ABx * V['NP12'] * V0
          + AB * (G('DQ0') * V0 * V1 + G('DQ1') * V0 ** 2) + AB * V0 * T12)
    tab = (F('d010') * V3 + (F('d010_a') + F('d100')) * V2 + F('d010_b') * sZ * V2 + F('d100_a') * V1 + F('d100_b') * sZ * V1 + 2 * F('d200') * sZ * V2
           + F('d200_a') * sZ * V1 + F('d200_b') * sZ ** 2 * V1 + F('d110') * sG * (3 * V3 + PSI1 * V2) + F('d110') * sZ * V3 + F('d110_a') * sZ * V2
           + F('d110_b') * sZ ** 2 * V2 + 3 * F('d300') * sZ ** 2 * V2 + F('d300_a') * sZ ** 2 * V1 + F('d300_b') * sZ ** 3 * V1 + F('A1') * V1 + F('B1') * V0
           + HD('hB0') * sw16 + HD('HB0') * V0)
    R5 = (I('d010') * D4 + (I('d010_a') + I('d100')) * W3 + I('d010_b') * S2r * W3 + I('d100_a') * B + I('d100_b') * S2r * B + 2 * I('d200') * S2r * W3
          + I('d200_a') * S2r * B + I('d200_b') * S2r ** 2 * B + I('d110') * (S1z * D4 + S2z * W3) + I('d110_a') * S2r * W3 + I('d110_b') * S2r ** 2 * W3
          + 3 * I('d300') * S2r ** 2 * W3 + I('d300_a') * S2r ** 2 * B + I('d300_b') * S2r ** 3 * B + I('A1') * B + I('B1') * W1 + I('B0') * A)
    I2 = (V['E1'] * V1 * V2 + V['E0'] * V0 * V2 + (V['CF0'] * D4 + V['CF1'] * W3 + V['CF2'] * B + V['CF3'] * W1) * DV
          + CBAR * V4 * DV + C1 * V3 * DV + RQ * DV + tab * DV + V['RWR'] * R5 * DV + C0W * DV
          + AB * (V['T2'] * V2 ** 2 + 2 * V['Q2B'] * V1 * V2 + V['Q2B'] * V0 * (V3 + PSI1 * V2)) + AB * V['TC2'] * QUF * W1 * DV + ABx * al * V['NAB'] * DV
          + BTx * np3 * DV + ABx * V['NP12'] * DV + AB * (G('DQ0') * V1 + G('DQ1') * V0) * DV + AB * T12 * DV)
    R0add = 7.52 * st_b * A
    fnw = V['FNW0'] * al + V['FNW1'] * al ** 2 + V['FNW2'] * al ** 3; fnwd = V['FNWd0'] * al + V['FNWd1'] * al ** 2 + V['FNWd2'] * al ** 3
    L0 = fnw * (V1 + V['KAP'] * V0) + ABc * al * V['NAB'] * V0 + BTc * np3 * V0 + ABc * V['NP12'] * V0
    L2 = fnwd * DV + ABc * al * V['NAB'] * DV + BTc * np3 * DV + ABc * V['NP12'] * DV
    LE = 5.72 * Fn * A + K.EPS * Fn * (D4 + 2.2 * W3 + 1.2 * B)
    L = LE + K.MU * L0 + K.MU * K.EPSI * L2 + K.MU * K.ABS_FLUX_I0 + K.MU * K.EPSI * K.ABS_FLUX_I2
    return dict(I0=I0, I2=I2, R0add=R0add, L=L)
rng = np.random.default_rng(7); worst = {}; n = 0
for trial in range(400):
    al = 10 ** rng.uniform(-7, math.log10(2.4e-4)); s = 10 ** rng.uniform(-9, -6.5)
    A, B, D4, V0, V2, V4 = (s * rng.uniform(0, 1), s * rng.uniform(0, 1) / math.sqrt(K.EPS), s * 10 ** rng.uniform(-2, 1.5) / math.sqrt(K.EPS),
                            s * rng.uniform(0, 1) / math.sqrt(K.MU), s * rng.uniform(0, 1) / math.sqrt(K.MU * K.EPSI), s * 10 ** rng.uniform(-2, 1.5) / math.sqrt(K.MU * K.EPSI))
    P = lambda x: ival.P(np.array([x]), np.array([x]))
    H, Lp, N, parts = R.rhs_pn(ival.cc(al), P(A), P(B), P(D4), P(V0), P(V2), P(V4), parts=True)
    f = spec_eval(al, A, B, D4, V0, V2, V4)
    for nm, (iv, fv) in dict(I0=(parts['I0p'], K.MU * f['I0']), I2=(parts['I2p'], K.MU * K.EPSI * f['I2']), R0add=(parts['R0add'], f['R0add']), L=(Lp, f['L'])).items():
        lo, hi = float(iv.lo[0]), float(iv.hi[0])
        ok = lo * (1 - 1e-12) - 1e-300 <= fv <= hi * (1 + 1e-12) + 1e-300
        assert ok, (nm, trial, lo, fv, hi)
        worst[nm] = max(worst.get(nm, 0.0), (hi - lo) / max(abs(fv), 1e-300))
    n += 1
out = dict(n_points=n, all_enclosed=True, max_rel_width=worst,
           args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('NEW', 'PYTHONHASHSEED', 'S0V', 'S1V', 'FLUX0', 'FLUX2', 'MOD', 'FIX', 'ADD')},
                     inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in [os.path.abspath(__file__), os.path.join(here, 'rhs_iv.py'), os.path.join(here, 'constants.py'), os.path.join(K.SRC, 'v4cert', 'out', 'constants_v4.json')]}))
json.dump(out, open(os.path.join(here, 't_float_check.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps(out['max_rel_width'])); print(f"T(iii) PASS: float spec evaluation enclosed by the interval rhs at {n} random points")
