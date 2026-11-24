#!/usr/bin/env python3
"""Task 5, Z1/Z2 (zone-III box from D5a), PRELIMINARY: python-flint arb, outward rounding; cell covers of width 1/4 lambda_1, monotone enclosures.
Written to quantify the zone-III box under the three coverings that are possible (see T5_REPORT.md: the contract's coverage of [190, 250) is a design question).

Definitions (lambda_1 units; r = coordinate along the S1 ray; f = slope over the ray = tan(psi), psi = tangent angle relative to the (straight) ray).
  D5a (Lemma O, under B4): on the conclusion region |x| <= 1 of a member of unit L (frame: tangent at the member centre at t_1),
      |g| <= d0, |g_x| <= d1, |g_xx| <= d2 (member units), g = h - h_0, with d_j = C_j sqrt(t) at the needed ratio times (50/L)^2 (cone), or the far values (L = 113).
  Data at t_1 (hand-off THRESHOLDS only): (H3) windows [y_w, 1.5 y_w], y_w = 50, 75, 112.5, 168.75, 253.125, weight 28/y:
      E_0 <= 2e-7, E_1, E_2 <= 1e-5;  |alpha(tau_1)| <= 1.5e-4 (B1 / hand-off), U_1(y) = 1.3056 alpha_1 y^p, p = 1.0246.
      Unweighted Eu_j <= (1.5 y_w/28) E_j;  sup|z_j| <= S_j = sqrt(Eu_j/ln1.5 + 2 sqrt(Eu_j Eu_{j+1}))  (z_{j+1} = d_s z_j);
      |v(t_1)| <= S_0,  |f(t_1)| <= F0 = U_1(y_b) + S_0,  |y f'(t_1)| <= D1 = p U_1(y_b) + S_1.
      (H5) far (own units, L = 113, whole cylinder): sup|h_0''| <= 0.04  (used only beyond the last window, Z2).
  Tilt slack: a curve point at ray coordinate r lies in the conclusion region of a member with |r - r_c| <= L (1 - e),
      e = (1 - cos th) + Hm sin th, th = atan(F0max) (axis vs ray), Hm = X tan(atan s0 + d1) + d0 (height in the member frame).
Bounds at a point r, time t, member m (unit L), t_1-point r' with |r - r'| <= d0 L (same x in the member frame):
  psi(r,t) <= atan F0(r') + d1                                    -> b0 = tan(.)
  |f - U_1| <= (1 + b0^2) d1 + S_0(r') + p U_1(r+) d0 L / r-          -> bz
  r|f'|    <= (1 + b0^2)^{3/2} [ D1(r') r/r' + (r/L)(d2 + m2c 1.5 d1 (2 s0 + d1)) ]     (k = h_xx sigma^{-3/2}, |k| <= |f'|)
      with m2c = 0.03 the (H5) cone curvature threshold (0.04 far) in the factor sigma^{-3/2} change only.
Z2 beyond the last window (r > 379.6875): psi(r,t_1) <= atan F0_last + (0.04/113) sqrt(1 + F^2)(r - 379.5) (curvature threshold, chained), then as above.
usage: PYTHONHASHSEED=0 nice -n 19 python3 -B z1_box.py     -> out/z1_box_prelim.json"""
import json, sys, os, math, hashlib, platform
from flint import arb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
CONSTS = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json'
A = lambda v: arb(repr(v)) if isinstance(v, float) else arb(v)
UP = lambda v: math.nextafter(float(arb(v).upper()), math.inf)
LO = lambda v: math.nextafter(float(arb(v).lower()), -math.inf)
def amax(xs):
    m = xs[0]
    for x in xs[1:]: m = x if x.upper() > m.upper() else m
    return arb(UP(m))
J = json.load(open(CONSTS))
dc = J['cases'][0]['deviation_upper_at_needed_time']; df = J['far_v32']['case']['deviation_upper_at_needed_time']
assert J['cases'][0]['name'] == 'S1 cone, X=3' and J['far_v32']['case']['X'] in ('2', 2)
# deviation values are upper bounds stored as doubles (already rounded up by lemmaO_v32_interval.py)
CONE = dict(d0=A(dc['C0']), d1=A(dc['C1']), d2=A(dc['C2']), X=arb(3), m2=arb('0.03'), s0=arb('0.08'))
FAR = dict(d0=A(df['C0']), d1=A(df['C1']), d2=A(df['C2']), X=arb(2), m2=arb('0.04'), s0=arb('0.08'), L=arb(113))
p = arb('1.0246'); al = arb('1.5e-4'); c13 = arb('1.3056'); ell = arb(1.5).log()
U1 = lambda y: c13 * al * arb(y) ** p
E = [arb('2e-7'), arb('1e-5'), arb('1e-5')]
WIN = []
for yw in ('50', '75', '112.5', '168.75', '253.125'):
    a = arb(yw); b = a * 3 / 2; Eu = [b / 28 * e for e in E]
    S0 = (Eu[0] / ell + 2 * (Eu[0] * Eu[1]).sqrt()).sqrt(); S1 = (Eu[1] / ell + 2 * (Eu[1] * Eu[2]).sqrt()).sqrt()
    WIN.append(dict(a=a, b=b, S0=S0, S1=S1, F0=U1(b) + S0, D1=p * U1(b) + S1))
RWEND = WIN[-1]['b']                                         # 379.6875
F0MAX = amax([w['F0'] for w in WIN])
def wins_meeting(lo, hi):
    ws = [w for w in WIN if w['a'].lower() <= hi.upper() and w['b'].upper() >= lo.lower()]
    return ws
def tilt(fam, scale):
    th = F0MAX.atan(); d0, d1 = fam['d0'] * scale, fam['d1'] * scale
    Hm = fam['X'] * (fam['s0'].atan() + d1).tan() + d0
    return (1 - th.cos()) + Hm * th.sin()
def point_bounds(lo, hi, L, d0, d1, d2, m2, s0, chained=False):
    """bounds on the cell r in [lo, hi] for a member of unit L with deviations d_j (member units)"""
    rl, rh = lo - d0 * L, hi + d0 * L                          # range of t_1-points r'
    if not chained:
        ws = wins_meeting(rl, rh); assert ws and rh.upper() <= RWEND.lower(), (float(lo.mid()), float(hi.mid()))
        F0 = amax([w['F0'] for w in ws]); S0 = amax([w['S0'] for w in ws]); D1 = amax([w['D1'] for w in ws])
        psi0 = F0.atan()
    else:                                                     # beyond the windows: chain the curvature threshold from 379.5
        wl = WIN[-1]; k0 = arb('0.04') / 113
        F = arb('0.5')                                        # a-priori |f| <= 1/2 used only inside sqrt(1+F^2); checked a posteriori
        psi0 = wl['F0'].atan() + k0 * (1 + F * F).sqrt() * (rh - arb('379.5'))
        F0 = psi0.tan(); S0 = None; D1 = k0 * (1 + F0 * F0) ** arb('1.5') * rh
        assert UP(F0) < 0.5
    b0 = (psi0 + d1).tan()
    sig = 1 + b0 * b0
    bz = sig * d1 + S0 + p * U1(rh) * d0 * L / lo if S0 is not None else None
    kdev = d2 + m2 * arb('1.5') * d1 * (2 * s0 + d1)
    b1 = sig ** arb('1.5') * (D1 * (hi / rl if not chained else arb(1)) + (hi / L) * kdev)
    return b0, b1, bz
def scan(lo, hi, members, step=arb(1) / 4, chained_from=None):
    """members(r_lo, r_hi) -> list of (L, scale, fam) admissible on the cell; take the best (min) bound per quantity, then the max over cells"""
    worst = dict(b0=None, b1=None, bz=None); n = 0
    r = arb(lo)
    while r.upper() < float(arb(hi).lower()):
        r2 = r + step
        if r2.upper() > float(arb(hi).upper()): r2 = arb(hi)
        ch = chained_from is not None and r2.upper() > float(chained_from)
        cands = []
        for (L, sc, fam) in members(r, r2):
            cands.append(point_bounds(r, r2, L, fam['d0'] * sc, fam['d1'] * sc, fam['d2'] * sc, fam['m2'], fam['s0'], chained=ch))
        assert cands, ('uncovered cell', float(r.mid()))
        best = [amax([c[i] for c in cands]) if False else min((c[i] for c in cands if c[i] is not None), key=lambda x: x.upper(), default=None) for i in range(3)]
        for k, v in zip(('b0', 'b1', 'bz'), best):
            if v is None: continue
            if worst[k] is None or v.upper() > worst[k][0].upper(): worst[k] = (v, float(r.mid()))
        r = r2; n += 1
    return {k: (dict(value=UP(v[0]), at_r=v[1]) if v else None) for k, v in worst.items()}, n
# ---- coverings
eC = tilt(CONE, arb(1)); eF = tilt(FAR, arb(1))
def far_members(lo, hi):                                       # far centres >= 376, unit 113: conclusion |r - r_c| <= 113 (1 - eF)
    return [(FAR['L'], arb(1), FAR)] if lo.lower() >= (376 - 113 * (1 - eF)).upper() else []
def fixedL(cmin):
    def m(lo, hi):
        L = arb(50); half = L * (1 - eC)
        ok = lo.lower() >= (arb(cmin) - half).upper() and hi.upper() <= (arb(380) + half).lower()
        return ([(L, arb(1), CONE)] if ok else []) + far_members(lo, hi)
    return m
def v32(lo, hi):                                               # v3.2 cone members r_m in [200, 380], L_m = r_m/4, conclusion |r - r_m| <= L_m (1 - e)
    out = []
    rm = arb(min(LO(4 * lo / (3 + eC)) * (1 - 1e-12), 380.0))   # (just below) the largest centre with r_m - L_m(1-e) <= lo, capped at 380
    Lm = rm / 4
    if lo.lower() >= (rm - Lm * (1 - eC)).upper() and hi.upper() <= (rm + Lm * (1 - eC)).lower() and rm.lower() >= 200:
        sc = (arb(50) / Lm) ** 2; out.append((Lm, sc, CONE))
    return out + far_members(lo, hi)
res = dict(); cells = {}
OPTS = {'A_fixedL_c300_zone[251.6,379.5]': (fixedL(300), None), 'B_fixedL_c238_zone[190,379.5]': (fixedL(238), None), 'C_v32members_zone[190,379.5]': (v32, None)}
lowA = UP(arb(300) - 50 * (1 - eC))
ZONES = {'A_fixedL_c300_zone[251.6,379.5]': (lowA, 379.5), 'B_fixedL_c238_zone[190,379.5]': (190, 379.5), 'C_v32members_zone[190,379.5]': (190, 379.5)}
for k, (mem, _) in OPTS.items():
    lo, hi = ZONES[k]
    res[k], cells[k] = scan(arb(repr(lo)) if isinstance(lo, float) else arb(lo), arb(repr(hi)), mem)
    b = res[k]
    res[k]['decision'] = 'primary' if (b['b0']['value'] <= 0.15 and b['b1']['value'] <= 0.52 and b['bz']['value'] <= 0.16) else \
        ('fallback' if (b['b0']['value'] <= 0.20 and b['b1']['value'] <= 0.60 and b['bz']['value'] <= 0.20) else 'fail')
    res[k]['margins_primary'] = dict(b0=0.15 / b['b0']['value'], b1=0.52 / b['b1']['value'], bz=0.16 / b['bz']['value'])
    res[k]['zone'] = [lo, hi]; res[k]['cells'] = cells[k]
    print(k, json.dumps(res[k]), flush=True)
# ---- Z2: [379.5, 700] with far members, curvature chained beyond the windows (option-independent part)
z2, n2 = scan(arb('379.5'), arb(700), far_members, step=arb(1), chained_from=RWEND.lower())
z2['HO1_ok_beyond_379.5'] = z2['b0']['value'] <= 0.5 and z2['b1']['value'] <= 1.0
print('Z2 [379.5,700]', json.dumps(z2), flush=True)
def sha(fn): return hashlib.sha256(open(fn, 'rb').read()).hexdigest()
out = dict(label='cert (arb) given B4 for the members used, D1 horizon, and the hand-off THRESHOLDS (H3) windows, alpha_1 <= 1.5e-4, (H5) far m2 = 0.04, cone m2 = 0.03, s0 = 0.08',
           provisional=True, note='PRELIMINARY: options A/B/C are the three possible coverings; the contract (zone III = [190, 379.5] with the fixed-L family, centres >= 300) is not realisable as written (T5_REPORT.md).',
           tilt=dict(e_cone=UP(eC), e_far=UP(eF), F0max=UP(F0MAX)),
           windows=[dict(yE=UP(w['a']), S0=UP(w['S0']), S1=UP(w['S1']), F0=UP(w['F0']), D1=UP(w['D1'])) for w in WIN],
           Z1=res, Z2_beyond=z2,
           args=dict(argv=sys.argv, env=dict(PYTHONHASHSEED=os.environ.get('PYTHONHASHSEED')), python=platform.python_version(), prec=ctx.prec,
                     inputs={CONSTS: sha(CONSTS)}, script_sha256=sha(os.path.abspath(__file__))))
json.dump(out, open(os.path.join(HERE, 'out', 'z1_box_prelim.json'), 'w'), indent=1, sort_keys=True)
