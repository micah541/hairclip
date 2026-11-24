"""Step 3 amendment v3.2 (draft), script for Theorem H6: closure of the Lemma O a-priori bounds (A) for the cone and far families, arc by arc.
python-flint arb; every printed margin is a lower bound (rounded down), every bound an upper bound (rounded up).
Nothing is sampled: the cone centre r_c runs over [200, 380] in cells [a, a+1] with monotone enclosures.

GEOMETRY (lambda_1 units, r = distance along the S1 ray).
  cone member m: centre r_m in [200, r_max], unit L_m = r_m/4, cylinder [r_m/4, 7 r_m/4], conclusion region |x| <= 1: [3r_m/4, 5r_m/4].
  far member:    centre >= 376, unit 113, cylinder +-226, conclusion region +-113  (so conclusion regions cover r >= 263).
  strip:         r in [50, 150): the physically fixed (B) windows present at tau_1, deviation against the FROZEN profile U_1 (Prop. B'),
                 levels 0-2 only: sup|v|, sup y|v'| by Sobolev (S[0], S[1]); the third derivative enters in L2 through E_1, E_2.
                 env STRIPX = x multiplies the strip energies by x (budget for the outer-ramp lemma R).
SCALING.  C_j do not depend on t, so for a member of unit L_m the deviations at the horizon are
  C1 sqrt t, C2 sqrt t at ratio 50^-4, times (50/L_m)^2;  sqrt(C3 t) likewise.  In the units of a cylinder c:
  |h_xx| scales by L_c/L_m,  ||h_xxx||_L2 by (L_c/L_m)^(3/2).
FRAMES.  k = curvature, phi = tangent angle in a frame, s = bound for |tan phi| in that frame:
  h_xx = sec^3(phi) k,   h_xxx = sec^4(phi) (k_sigma + 3 tan(phi) k^2),   dx = cos(phi) dsigma.   Hence for frames c (target), m (source)
  |h^c_xx|  <= (1+s_c^2)^(3/2) |k|,  |k| <= |h^m_xx|,
  |h^c_xxx| <= (1+s_c^2)^2 (|h^m_xxx| + 3 (s_c + s_m) k^2),   dx_c <= (1+s_m^2)^(1/2) dx_m.
  From the ray chart (slope f over the ray, ' = d/dy): |k| <= |f'|, |k_sigma| <= |f''| + 3 |f| f'^2, dx_c <= (1+f^2)^(1/2) dy.
DATA THRESHOLDS (new hand-off hypotheses "H6-data", own cylinder units, on the whole cylinder):
  cone: slope <= S0C, |h0''| <= M2C, ||h0'''||_L2(-3,3) <= N3C;   far: slope <= S0F, |h0''| <= M2F, ||h0'''||_L2(-2,2) <= N3F.
usage: python3 -I v32_h6.py [far_tag] [380] [cone_tag]     (env MU = mu: all deviations taken at mu times the horizon)"""
import sys, os, json, math
from flint import arb, ctx
ctx.prec = 128
HERE = '/home/claude/s1spec/writeup/step3_amend_v32/scripts'
O = '/home/claude/s1spec/step4/p2/out'
FAR_TAG = sys.argv[1] if len(sys.argv) > 1 else 'far X=2 M=0.15 K=0.3'
RMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 380
assert RMAX == 380, 'the piece lists below are written for cone centres [200, 380]'
CONE_TAG = sys.argv[3] if len(sys.argv) > 3 else 'cone X=3 M=0.12 K=0.18'
A = lambda v: arb(repr(v)) if isinstance(v, float) else arb(v)
U = lambda v: math.nextafter(float(arb(v).upper()), math.inf)
Lw = lambda v: math.nextafter(float(arb(v).lower()), -math.inf)
def exmax(a, b): return a if a.lower() >= b.upper() else (b if b.lower() >= a.upper() else a.union(b))   # enclosure of max(a, b)
def amax(*xs):
    m = xs[0]
    for x in xs[1:]: m = m if m.upper() >= x.upper() else x
    return arb(U(m))                                    # an upper bound of the maximum
C = json.load(open(f'{HERE}/v32_cases.json'))
cone, far = C[CONE_TAG], C[FAR_TAG]
KC = arb(CONE_TAG.split('K=')[1]); KF = arb(FAR_TAG.split('K=')[1]); MF = arb(FAR_TAG.split('M=')[1].split()[0])
dc, df = cone['deviation_upper_at_needed_time'], far['deviation_upper_at_needed_time']
Q = arb(os.environ.get('MU', '1')).sqrt()                               # MU > 1: all deviations at mu times the horizon (C_j sqrt(mu t))
C1c, C2c, G3c = Q * A(dc['C1']), Q * A(dc['C2']), Q * A(dc['sqrtC3t'])  # cone, at ratio 50^-4 (1.6001e-7)
C1f, C2f, G3f = Q * A(df['C1']), Q * A(df['C2']), Q * A(df['sqrtC3t'])  # far, at ratio 113^-4 (6.134e-9)
# ---- hand-off data (step 4, certified) and the proposed thresholds
t9 = json.load(open(f'{O}/t9_cone_kch.json'))['cone_max']; tf = json.load(open(f'{O}/t9_far.json'))
h4 = json.load(open(f'{O}/h4loc_A1_w180.json'))['windows']['constants/amp_A1_w180.json']['cylinders']
TH = dict(S0C='0.08', M2C='0.03', N3C='0.004', S0F='0.08', M2F='0.04', N3F='0.004')
S0C, M2C, N3C, S0F, M2F, N3F = [arb(TH[k]) for k in ('S0C', 'M2C', 'N3C', 'S0F', 'M2F', 'N3F')]
data = dict(cone_slope=A(t9['slope']), cone_h2=A(t9['K2']), cone_N3=A(t9['N3']),
            far_slope=A(tf['slope']), far_h2=113 * A(tf['h2']), far_h3=arb(113) ** 2 * A(tf['h3']), far_h4=arb(113) ** 3 * A(tf['h4']),
            far_N3=arb(113) ** 3 * A(tf['L2']) ** 2)
data['far_N3'] = data['far_N3'].sqrt()
far_h4_ratio_M1 = max(c['ratio'] for c in h4 if c['name'].startswith('far'))            # int (h0^(5))^2 in units M^2 = 1
cone_h4_ratio = max(c['ratio'] for c in h4 if c['name'].startswith('cone') and int(c['name'].split('-')[1]) <= RMAX)
# ---- windows: scratch rerun of lemmaS_Bbox with the frozen profile (alpha_1 = 1.5e-4, hand-off E = 1e-5), floats inflated by 1e-6
import os
WBF = os.environ.get('WINDOWS', f'{HERE}/../scratch_lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15_frozen.json')
WB = json.load(open(WBF)); SX = arb(os.environ.get('STRIPX', '1'))
infl = arb('1.000001'); p = arb('1.0246'); al = arb('1.5e-4')
Ufun = lambda y: arb('1.3056') * al * arb(y) ** p
wins = []
for w in WB:
    yE = arb(repr(w['yE'])); yb = yE * 3 / 2; S = [A(x) * infl * SX.sqrt() for x in w['S']]
    Eu = [yb / 28 * A(x) * infl * SX for x in w['EB']]  # unweighted window energies int z_j^2 ds (omega = 28/y >= 28/yb)
    F0 = Ufun(yb) + S[0]                               # |f|
    D1 = arb('1.05') * Ufun(yb) + S[1]                 # y |f'|
    T2 = S[2] + S[1] + (p * p - p) * Ufun(yb)          # y^2 |f''|
    wins.append(dict(yE=yE, yb=yb, F0=F0, D1=D1, T2=T2, Eu=Eu, Ub=Ufun(yb)))
strip_w = [w for w in wins if w['yE'].lower() < 150]   # windows meeting [50, 150)
F0s = amax(*[w['F0'] for w in strip_w])
sq = lambda s: 1 + s * s
tan_add = lambda s0, d: (arb(U(s0)).atan() + arb(U(d))).tan()
# ---- slopes during [t_1, T)
s_c_cov = tan_add(S0C, amax(C1c, C1f))                 # cone frame, at points covered by a member (angle change <= C1 sqrt t)
FAX = amax(*[w['F0'] for w in wins])                   # |f(r_c, tau_1)| <= derived |f| on the hand-off windows up to y = 380: axis angle of cylinder c
s_c_str_a = tan_add(S0C, 2 * F0s.atan())               # strip points: data slope + change of the angle to the ray (<= 2 atan sup|f|)
s_c_str_b = (F0s.atan() + FAX.atan()).tan()            # strip points: angle to the ray at (p, t) + angle of the axis of c to the ray
s_c_str = s_c_str_a if s_c_str_a.upper() <= s_c_str_b.upper() else s_c_str_b
s_f_cov = tan_add(S0F, amax(C1c, C1f))
s_c = amax(s_c_cov, s_c_str)
R = {}
R['slope'] = dict(cone_covered=U(s_c_cov), cone_strip=U(s_c_str), cone_strip_a=U(s_c_str_a), cone_strip_b=U(s_c_str_b), far=U(s_f_cov),
                  margin_cone=Lw(arb('0.25') / s_c), margin_far=Lw(arb('0.25') / s_f_cov))
# ---- cone cylinders, r_c in [a, a+1]
def cone_cell(a, b):
    Lc = arb(b) / 4; Lc_lo = arb(a) / 4; lo = arb(a) / 4; hi = arb(7 * b) / 4
    # K2: data part (own frame, own units) + deviation part
    dev_cov = (Lc / 50) * C2c                          # worst covered point r_p = 150: member r_m = 200, L_m = 50, ratio L_c/50
    dev_far = (Lc / 113) * C2f if 7 * b / 4 > 475 else arb(0)
    k_cov = sq(S0C) ** arb('1.5') * M2C + amax(dev_cov, dev_far)      # bound for L_c |k| at covered points
    K2_cov = sq(s_c_cov) ** arb('1.5') * k_cov
    k_str = arb(0)
    for w in strip_w:
        if w['yb'].lower() <= a / 4: continue
        y0 = exmax(w['yE'], lo); k_str = amax(k_str, Lc * w['D1'] / y0)
    K2_str = sq(s_c_str) ** arb('1.5') * k_str
    K2 = amax(K2_cov, K2_str); kap = amax(k_cov, k_str)
    # K3: pieces
    Fc = sq(s_c_cov) ** 2
    def piece(Lm, n3m, g3m, sm, length):               # member piece, in c-units (covered points: c-slope <= s_c_cov, L_c|k| <= k_cov)
        return Fc * (sq(sm) ** arb('0.25') * (Lc / Lm) ** arb('1.5') * (n3m + g3m) + 3 * (s_c_cov + sm) * k_cov ** 2 * (arb(length) / Lc_lo).sqrt())
    tot = arb(0); parts = {}
    # strip [r_c/4, 150)
    acc = arb(0)
    for w in strip_w:
        y0 = exmax(w['yE'], lo); y1 = w['yb'] if w['yb'].upper() < 150 else arb(150)
        if y1.lower() <= y0.upper(): continue
        # ||f''||_L2(dy) <= ||v''|| + ||U_1''|| ; int v''^2 dy = int (z_2 - z_1)^2 y^-3 ds <= y0^-3 (sqrt Eu_2 + sqrt Eu_1)^2 ; U_1'' = p(p-1) U_1 / y^2
        pw = ((y0 ** -3 - y1 ** -3) / 3).sqrt()
        n_w = y0 ** arb('-1.5') * (w['Eu'][2].sqrt() + w['Eu'][1].sqrt()) + ((p * p - p) * w['Ub'] + 3 * (w['F0'] + s_c_str) * w['D1'] ** 2) * pw
        acc += n_w ** 2
    parts['strip'] = sq(s_c_str) ** 2 * sq(F0s) ** arb('0.25') * Lc ** arb('1.5') * acc.sqrt()
    sm_c = tan_add(S0C, C1c); sm_f = tan_add(S0F, C1f)
    parts['m200'] = piece(arb(50), N3C, G3c, sm_c, 100)                                    # [150, 250]
    L333 = arb(250) / 3
    parts['m333'] = piece(L333, N3C, G3c * (50 / L333) ** 2, sm_c, arb(500) / 3)           # [250, 416.7]
    if 7 * b / 4 > 1250 / 3: parts['m380'] = piece(arb(95), N3C, G3c * (arb(50) / 95) ** 2, sm_c, arb(175) / 3)   # [416.7, 475]
    if 7 * b / 4 > 475: parts['far'] = piece(arb(113), N3F, G3f, sm_f, 190)               # (475, 665] inside one far conclusion region
    for v in parts.values(): tot += v * v
    return K2, tot.sqrt(), {k: U(v) for k, v in parts.items()}, U(K2_cov), U(K2_str)
worst2 = worst3 = None
for a in range(200, RMAX):
    K2, K3, parts, kc, ks = cone_cell(a, a + 1)
    if worst2 is None or K2.upper() > worst2[0].upper(): worst2 = (K2, a, kc, ks)
    if worst3 is None or K3.upper() > worst3[0].upper(): worst3 = (K3, a, parts)
dev2 = (arb(RMAX) / 200) * C2c
R['cone'] = dict(K=float(KC.mid()), r_max=RMAX,
    K2_recovered=U(worst2[0]), K2_cell=worst2[1], K2_covered=worst2[2], K2_strip=worst2[3], K2_margin=Lw(KC / worst2[0]),
    K2_time_margin=Lw(((KC / sq(s_c_cov) ** arb('1.5') - sq(S0C) ** arb('1.5') * M2C) / dev2) ** 2),
    K3_recovered=U(worst3[0]), K3_cell=worst3[1], K3_parts=worst3[2], K3_margin=Lw(KC / worst3[0]),
    tstar_margin=cone['margin_factor_lower'], H4_ratio=cone_h4_ratio, H4_margin=1 / cone_h4_ratio)
# ---- far cylinders (centre >= 376): inner points [150, 263) from cone members, the rest from far members
Lf = arb(113)
kf = sq(S0F) ** arb('1.5') * M2F + amax(C2f, (Lf / 50) * C2c)
K2f = sq(s_f_cov) ** arb('1.5') * kf
Ff = sq(s_f_cov) ** 2; sm_c = tan_add(S0C, C1c); sm_f = tan_add(S0F, C1f)
def fpiece(Lm, n3m, g3m, sm, length): return Ff * (sq(sm) ** arb('0.25') * (Lf / Lm) ** arb('1.5') * (n3m + g3m) + 3 * (s_f_cov + sm) * kf ** 2 * (arb(length) / Lf).sqrt())
L333 = arb(250) / 3
fp = dict(m200=fpiece(arb(50), N3C, G3c, sm_c, 100), m333=fpiece(L333, N3C, G3c * (50 / L333) ** 2, sm_c, 13),
          far1=fpiece(Lf, N3F, G3f, sm_f, 226), far2=fpiece(Lf, N3F, G3f, sm_f, 226))
K3f = sum((v * v for v in fp.values()), arb(0)).sqrt()
R['far'] = dict(tag=FAR_TAG, K2_recovered=U(K2f), K2_margin=Lw(KF / K2f), K3_recovered=U(K3f), K3_parts={k: U(v) for k, v in fp.items()},
                K3_margin=Lw(KF / K3f), tstar_margin=far['margin_factor_lower'],
                H4_ratio_at_M=U(arb(repr(far_h4_ratio_M1)) / MF ** 2), H4_margin=Lw(MF ** 2 / arb(repr(far_h4_ratio_M1))))
# ---- data against thresholds
R['data_vs_thresholds'] = dict(
    cone_slope=Lw(S0C / data['cone_slope']), cone_h2=Lw(M2C / data['cone_h2']), cone_N3=Lw(N3C / data['cone_N3']),
    far_slope=Lw(S0F / data['far_slope']), far_h2=Lw(M2F / data['far_h2']), far_N3=Lw(N3F / data['far_N3']),
    far_h3_vs_M=Lw(MF / data['far_h3']), far_h4_vs_M=Lw(MF / data['far_h4']), far_slope_vs_M=Lw(MF / data['far_slope']))
R['data'] = {k: U(v) for k, v in data.items()}
R['thresholds'] = dict(TH, far_M=float(MF.mid()), far_K=float(KF.mid()), cone_M=0.12, cone_K=float(KC.mid()))
R['windows'] = [dict(yE=U(w['yE']), F0=U(w['F0']), D1=U(w['D1']), Eu1=U(w['Eu'][1]), Eu2=U(w['Eu'][2])) for w in wins]
R['far_physical_thresholds_lambda1_units'] = dict(h2=Lw(M2F / 113), h3=Lw(MF / arb(113) ** 2), h4=Lw(MF / arb(113) ** 3),
                                                   L2_h3=Lw(N3F / arb(113) ** arb('1.5')), slope=float(S0F.mid()))
R['args'] = dict(cone=CONE_TAG, far=FAR_TAG, r_max=RMAX, cases='v32_cases.json', t9=f'{O}/t9_cone_kch.json', t9_far=f'{O}/t9_far.json',
                 h4=f'{O}/h4loc_A1_w180.json', windows=WBF, STRIPX=os.environ.get('STRIPX', '1'))
tag = f"{CONE_TAG.split('K=')[1]}_{FAR_TAG.split('M=')[1].replace(' K=', '_')}_{RMAX}"
R['MU'] = os.environ.get('MU', '1')
mins = dict(slope_cone=R['slope']['margin_cone'], slope_far=R['slope']['margin_far'], cone_K2=R['cone']['K2_margin'], cone_K3=R['cone']['K3_margin'],
            far_K2=R['far']['K2_margin'], far_K3=R['far']['K3_margin'])
R['closure_margins'] = mins; R['min_closure_margin'] = min(mins.values())
if os.environ.get('MU') is None and os.environ.get('STRIPX') is None and os.environ.get('WINDOWS') is None and len(sys.argv) == 1:
    json.dump(R, open(f'{HERE}/v32_h6.json', 'w'), indent=1); print(json.dumps(R, indent=1))
else: print('MU', R['MU'], 'STRIPX', os.environ.get('STRIPX', '1'), tag, {k: round(v, 3) for k, v in mins.items()}, 'K3 strip part', R['cone']['K3_parts']['strip'], 'K2 strip', R['cone']['K2_strip'])
