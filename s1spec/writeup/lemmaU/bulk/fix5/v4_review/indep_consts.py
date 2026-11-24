#!/usr/bin/env python3
"""Reviewer's independent recomputation (python-flint arb, 200 bits, outward) of the most important fin13v4 constants.
Written from EI_SPEC.md definitions, not from v4cert/*.py.  Read-only on the proof tree.  Output: indep_consts.json (this directory).
Each entry: mine (certified upper bound, or lower bound for 'down'), v4 (value in constants_v4.json), ok = (v4 covers the true value as far as
my enclosure can tell: for 'up' constants ok means v4 >= my LOWER enclosure of the sup; 'tight' = v4 / my upper bound).
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B indep_consts.py"""
import json, math, os, sys
from flint import arb, acb, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
V4 = '/home/claude/s1spec/step3nl_v4'
C = json.load(open(V4 + '/v4cert/out/constants_v4.json'))['constants']
C2 = json.load(open(V4 + '/v4cert/out/c2_t3d.json'))
C1j = json.load(open(V4 + '/v4cert/out/c1_cut_sups.json'))
res = {}
def up(z): return math.nextafter(float(arb(z).upper()), math.inf)
def dn(z): return math.nextafter(float(arb(z).lower()), -math.inf)
def rec(name, lo, hi, key=None, note=''):
    """lo/hi: my rigorous lower/upper enclosure of the true sup; v4 value must be >= lo (necessary) ; v4 >= hi is sufficient"""
    v = C[key or name]['value']
    res[name] = dict(mine_lower=lo, mine_upper=hi, v4=v, v4_ge_mine_lower=bool(v >= lo), v4_ge_mine_upper=bool(v >= hi), v4_over_mine_upper=v / hi if hi else None, note=note)
# ------------------------------------------------------------------ cutoff and weights
def S7(t): return t ** 4 * (35 - 84 * t + 70 * t ** 2 - 20 * t ** 3)
def S7_1(t): return 140 * t ** 3 * (1 - t) ** 3
def S7_2(t): return 420 * t ** 2 * (1 - t) ** 2 * (1 - 2 * t)
def S7_3(t): return 840 * t * (1 - t) * (1 - 5 * t + 5 * t * t)
def S7_4(t): return 840 * (1 - 2 * t) * (1 - 10 * t + 10 * t * t)
import sympy as sp
_t = sp.symbols('t'); _S = _t ** 4 * (35 - 84 * _t + 70 * _t ** 2 - 20 * _t ** 3)
assert sp.expand(sp.diff(_S, _t) - 140 * _t ** 3 * (1 - _t) ** 3) == 0
assert sp.expand(sp.diff(_S, _t, 2) - 420 * _t ** 2 * (1 - _t) ** 2 * (1 - 2 * _t)) == 0
assert sp.expand(sp.diff(_S, _t, 3) - 840 * _t * (1 - _t) * (1 - 5 * _t + 5 * _t ** 2)) == 0
assert sp.expand(sp.diff(_S, _t, 4) - 840 * (1 - 2 * _t) * (1 - 10 * _t + 10 * _t ** 2)) == 0
def zin(X, j):          # zeta = 1 - S7((16.25 - x)/6) on IR
    u = (arb('16.25') - X) / 6
    if j == 0: return 1 - S7(u)
    f = (S7_1, S7_2, S7_3, S7_4)[j - 1]
    return -f(u) * (arb(-1) / 6) ** j
def zout(X, j):
    u = (arb(60) - X) / 3
    if j == 0: return S7(u)
    f = (S7_1, S7_2, S7_3, S7_4)[j - 1]
    return f(u) * (arb(-1) / 3) ** j
def psiA(X): return 4 * X / (1 + X * X)
def psiB(X): return X * (1 + X * X) ** (-arb(1) / 3) / 3          # phi' , phi = (1+x^2)^(2/3)/4
def psip(X, a, b):
    if b <= 30.25: return psiA(X)
    if a >= 38.25: return psiB(X)
    s = S7((X - arb('30.25')) / 8); return (1 - s) * psiA(X) + s * psiB(X)
def psipp(X, a, b):
    A2 = 4 * (1 - X * X) / (1 + X * X) ** 2
    B2 = (1 + X * X) ** (-arb(1) / 3) / 3 - 2 * X * X * (1 + X * X) ** (-arb(4) / 3) / 9
    if b <= 30.25: return A2
    if a >= 38.25: return B2
    t = (X - arb('30.25')) / 8; s = S7(t); s1 = S7_1(t) / 8
    return (1 - s) * A2 + s * B2 + s1 * (psiB(X) - psiA(X))
def WR(X): return ((1 + X * X) ** (arb(2) / 3) / 4 - 2 * (1 + X * X).log()).exp()        # W/rho, x <= 30.25
def supinf(f, a, b, n):
    """(max of lower endpoints at midpoints = rigorous LOWER bound of the sup, max of upper endpoints over the cover = UPPER bound)"""
    hi = None; lo = None; h = (b - a) / n
    for i in range(n):
        xa = a + i * h; xb = b if i == n - 1 else a + (i + 1) * h
        X = arb(xa).union(arb(xb)); v = f(X, xa, xb); assert v.is_finite()
        hi = v.upper() if hi is None or v.upper() > hi else hi
        xm = arb(0.5 * (xa + xb)); vm = f(xm, xa, xb)
        lo = vm.lower() if lo is None or vm.lower() > lo else lo
    return dn(lo), up(hi)
def over(f, a, b, per=300):
    pts = [a] + [c for c in (16.25, 30.25, 38.25, 57.0) if a < c < b] + [b]
    L = H = None
    for u, v in zip(pts[:-1], pts[1:]):
        l, h = supinf(f, u, v, max(16, int((v - u) * per)))
        L = l if L is None else max(L, l); H = h if H is None else max(H, h)
    return L, H
NI = 4000
rec('PSI1', *over(lambda X, a, b: psip(X, a, b), 10.25, 60.0))
rec('PSI2', *over(lambda X, a, b: abs(psip(X, a, b) ** 2 - psipp(X, a, b)), 10.25, 60.0))
l, h = supinf(lambda X, a, b: WR(X).sqrt(), 10.25, 16.25, NI); rec('RWR', l, h)
kw = (1 + arb('16.25') ** 2) ** 2; rec('KWIR', dn(kw), up(kw))
rec('KAP_IR', *supinf(lambda X, a, b: abs(zin(X, 1)) + psip(X, a, b) * zin(X, 0), 10.25, 16.25, NI))
rec('KAP_Z1', *over(lambda X, a, b: psip(X, a, b), 16.25, 57.0))
rec('KAP_OR', *supinf(lambda X, a, b: abs(zout(X, 1)) + psip(X, a, b) * zout(X, 0), 57.0, 60.0, NI))
for j, bn in ((1, 4), (2, 6), (3, 4), (4, 1)):
    rec(f'k4_{j}', *supinf(lambda X, a, b, j=j, bn=bn: bn * abs(zin(X, j)) * WR(X).sqrt(), 10.25, 16.25, NI))
for j, bn in ((1, 3), (2, 3), (3, 1)):
    rec(f'k3_{j}', *supinf(lambda X, a, b, j=j, bn=bn: bn * abs(zin(X, j)) * WR(X).sqrt(), 10.25, 16.25, NI))
rec('Z2i', *supinf(lambda X, a, b: abs(zin(X, 2)), 10.25, 16.25, NI))
# psi >= 0 increasing? (needed for the left-end Sobolev S1z, S2z and for 1/W increasing): inf psi' on [10.25,60] > 0
l_, h_ = over(lambda X, a, b: -psip(X, a, b), 10.25, 60.0)
res['inf_psi_prime_[10.25,60]'] = dict(lower_bound=-h_, positive=bool(-h_ > 0))
# ------------------------------------------------------------------ W(34.25), Phi(34.25) (W = (1+x^2)^-2 on x <= 30.25; blend beyond)
def psip_c(z, _):
    s = S7((z - acb('30.25')) / 8)
    return (1 - s) * 4 * z / (1 + z * z) + s * z * ((1 + z * z).log() * (-acb(1) / 3)).exp() / 3
def psi_at(x):
    if x <= 30.25: return 2 * (1 + arb(x) ** 2).log()
    return 2 * (1 + arb('30.25') ** 2).log() + acb.integral(psip_c, acb('30.25'), acb(x)).real
KW = psi_at(34.25).exp()
rec('KW', dn(KW), up(KW), note='1/W(34.25)')
# Phi(34.25) = int_16.25^30.25 (1+x^2)^2 dx + int_30.25^34.25 exp(psi)
P1 = lambda x: x + 2 * x ** 3 / 3 + x ** 5 / 5
I1 = P1(arb('30.25')) - P1(arb('16.25'))
n = 800; I2lo = arb(0); I2hi = arb(0); hh = arb(4) / n; prev = psi_at(30.25)
for i in range(n):
    xb = 30.25 + 4.0 * (i + 1) / n
    cur = psi_at(xb)
    I2lo += hh * prev.exp(); I2hi += hh * cur.exp(); prev = cur       # 1/W = e^psi increasing
PHI_lo, PHI_hi = I1 + I2lo, I1 + I2hi
rec('PHI34', dn(PHI_lo), up(PHI_hi), note='coupling raw; constants.py uses min(PHI34, PHI34_c4)')
# ------------------------------------------------------------------ K_rho(b) = 1/rho(b) = exp((1+b^2)^(2/3)/4)
for b_, key in ((1, 'KREG[0,1]'), (2, 'KREG[1,2]'), (4, 'KREG[2,4]'), (8, 'KREG[4,8]'), (12, 'KREG[8,12]'), (16.25, 'K_rho_16.25')):
    kr = ((1 + arb(b_) ** 2) ** (arb(2) / 3) / 4).exp(); rec(key, dn(kr), up(kr), note='1/rho(b)')
# ------------------------------------------------------------------ transport T0, T2 (profile bounds at sigma >= 8 taken from c1 json [their cert])
pr = C1j['profile_sigma8']; KB, MB, PB, QS = (arb(pr[k_]) for k_ in ('kb', 'mb', 'pb', 'qs'))
TIP = arb('-2.2490776515195').union(arb('-2.2490776505754'))
def Tn(lev):
    c5 = 1 if lev == 0 else 5
    def f(X, a, b):
        qq = X + TIP + arb(0, 1) * QS; qp = 1 + arb(0, 1) * KB * PB
        ca = psip(X, a, b) * qq / 2 + c5 * qp / 2 - 1; cb = psip(X, a, b) / 2 + c5 * KB / 2
        return (ca * ca + cb * cb).sqrt()
    return over(f, 10.25, 60.0)
rec('T0', *Tn(0)); rec('T2', *Tn(2))
# ------------------------------------------------------------------ C2 top-order coefficients from closed forms (I13 [sym]): c = 1 - J^-4,
# c_U1 = 4 U1 J^-6, d110 = 10 U1 J^-6, d300 = 3 J^-6 (1 - 6 U1^2 J^-2);  J^2 = (1 - k U0)^2 + U1^2,  |U1| <= sup|u0'| + e1, |k U0| <= kU (their table)
for lev in ('0.00025', '0.001'):
    for zn, e1 in (('IR', 1e-3), ('ZA', 0.02), ('ZB', 0.02)):
        T = C2['tables'][lev][zn]['sup']
        u1 = arb(T["sup|u0'|"]) + arb(e1); kU = arb(T['kU'])
        J2lo = (1 - kU) ** 2; J2hi = (1 + kU) ** 2 + u1 * u1
        c_hi = max(abs(1 - 1 / J2lo ** 2), abs(1 - 1 / J2hi ** 2), key=lambda z: float(z.upper()))
        cb_hi = 4 * u1 / J2lo ** 3; d110_hi = 10 * u1 / J2lo ** 3; d300_hi = 3 / J2lo ** 3
        for nm, mine in (('c', c_hi), ('c_b', cb_hi), ('d110', d110_hi), ('d300', d300_hi)):
            key = f'TQ[{lev}][{zn}][{nm}]'; v = C[key]['value']
            # 'mine' is a valid upper bound GIVEN their sup|u0'| and kU; it is attained up to the box structure (zone-wide sup of |u0'| vs per-box), so v4 may be slightly smaller
            res[key] = dict(mine_upper_from_closed_form=up(mine), v4=v, v4_over_mine=v / up(mine), sup_u0p=T["sup|u0'|"], kU=T['kU'])
# ------------------------------------------------------------------ U-3 strip form, recomputed from mod_ibp.json (my own code)
Mj = json.load(open(V4 + '/agent_mod/ibp/mod_ibp.json'))['alpha_max_5e-4']
A_ = lambda v: arb(repr(v)); SBx = arb('0.02'); ETA = A_(Mj['eta'])
comps = ['a', 'b', 'beta']; d = [A_(x) for x in Mj['neumann']['d']]; gam = [A_(Mj['final']['core'][c]['gamma_ref']) for c in comps]
form = {}
for c in comps:
    s = Mj['strip'][c]; g = {k_: A_(v) for k_, v in s['quad_gamma_L1'].items()}; assert set(g) == {'00', '01', '11', '02', '12', '22'}
    kap = [A_(v) for v in s['commutator_S0_S3']]; M1, M2, M3 = (A_(s[f'cubic_M{i}_L1']) for i in (1, 2, 3))
    # quadratic: g_ab S_a S_b <= SB g_ab S_min(a,b);  cubic: eta[M3/3 (S0+S1)^2 + M2 (S0+S1) S2 + 2 M1 S2^2] <= eta[(2 SB M3/3 + SB M2)(S0+S1) + 2 SB M1 S2]
    form[c] = [kap[0] + SBx * (g['00'] + g['01'] + g['02']) + ETA * (2 * SBx * M3 / 3 + SBx * M2),
               kap[1] + SBx * (g['11'] + g['12']) + ETA * (2 * SBx * M3 / 3 + SBx * M2),
               kap[2] + SBx * g['22'] + 2 * ETA * SBx * M1, kap[3], A_(s['linear_alpha_S0'])]
def amax(l): return arb(max(float(z.upper()) for z in l)) * (1 + arb(2) ** -50)
fb = {c: [form[c][r] + gam[i] * amax([form[k_][r] / d[j] for j, k_ in enumerate(comps)]) for r in range(5)] for i, c in enumerate(comps)}
for r, key in enumerate(('SF_beta_0', 'SF_beta_1', 'SF_beta_2', 'SF_beta_3', 'SF_beta_al')):
    rec(key, dn(fb['beta'][r]), up(fb['beta'][r]))
for r, key in enumerate(('SF_ab_0', 'SF_ab_1', 'SF_ab_2', 'SF_ab_3', 'SF_ab_al')):
    z = (fb['a'][r] ** 2 + fb['b'][r] ** 2).sqrt(); rec(key, dn(z), up(z))
sbeta = sum((fb['beta'][r] for r in range(4)), arb(0)) * SBx + arb('5e-4') * fb['beta'][4] * SBx
rec('STRIP_BETA', dn(sbeta), up(sbeta), note='S0..S3 <= 0.02 (S3: hypothesis), |alpha| <= 5e-4')
res['STRIP_BETA_share_of_S3_term'] = float((fb['beta'][3] * SBx / sbeta).mid())
# ------------------------------------------------------------------ raw read-backs (source JSON vs constants_v4.json, direction)
G = json.load(open(V4 + '/agent_int/v2/v2_T4_P2_0.25.json'))['sup']
mu1 = min(-v['-mu1'] for v in G.values()); mu2 = min(-v['-mu2'] for v in G.values())
res['MU1'] = dict(raw=mu1, v4=C['MU1']['value'], ok_down=bool(C['MU1']['value'] <= mu1)); res['MU2'] = dict(raw=mu2, v4=C['MU2']['value'], ok_down=bool(C['MU2']['value'] <= mu2))
cp = json.load(open(V4 + '/agent_coup/coupling_constants.json'))['P']['2']['rigorous']
for nm, k_ in (('AF', 'a_f'), ('BF', 'b_f'), ('CF0', 'c3'), ('CF1', 'c2'), ('CF2', 'c1'), ('CF3', 'c0')):
    res[nm] = dict(raw=cp[k_], v4=C[nm]['value'], ok_up=bool(C[nm]['value'] >= cp[k_]))
bad = [k_ for k_, e in res.items() if e is not None and isinstance(e, dict) and (e.get('v4_ge_mine_lower') is False or e.get('ok_up') is False or e.get('ok_down') is False)]
res['_summary'] = dict(n=len(res), v4_below_my_lower_bound=bad, not_covering_my_upper=[k_ for k_, e in res.items() if isinstance(e, dict) and e.get('v4_ge_mine_upper') is False])
json.dump(res, open(os.path.join(HERE, 'indep_consts.json'), 'w'), indent=1)
for k_, e in res.items(): print(k_, e)
