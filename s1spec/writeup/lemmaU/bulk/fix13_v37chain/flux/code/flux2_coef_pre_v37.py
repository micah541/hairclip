#!/usr/bin/env python3
"""Stage X2 (FLUX2_DESIGN.md sec. 8, constants #1-#7): certified coefficient constants of flux lemma v2 on OR = [57, 60].
python-flint arb, 200 bits; every sup rounded UP (one ulp, f2common.fup), every lower bound DOWN.  Labels: [cert] unless stated.

#1  z_i = sup_OR |zeta^(i)|, i = 0..4 (zeta = S7((60 - x)/3)), NCELL cells, mean-value form on the clipped cell of u = (60 - x)/3
    (F1's method, fix3/flux/flux_cert.py); lower checks z_1 >= 35/48, z_4 >= 280/27 (exact values at u = 1/2 and at the ends).
#2  frak a_i = sup_OR |coefficient of g^(i) in D2v| (X1 form (c)), KAP2 = sup_OR |2 zeta zeta' - psi' zeta^2|, sup psi', inf psi',
    sup |psi2|, psi2 = psi'^2 - psi''.  On x >= 38.25, psi' = phi' = x/(3(1+x^2)^(1/3)) exactly (verifier weight, v4cert/c1_cut_sups.py:4).
#3  W(57) = exp(-psi(57)) with the verifier's psi: psi = 2 log(1+x^2) on x <= 30.25 (W = (1+x^2)^-2), psi' = (1-s) 4x/(1+x^2) + s phi',
    s = S7((x - 30.25)/8) on [30.25, 38.25], psi' = phi' beyond (so psi(57) = psi(38.25) + phi(57) - phi(38.25)).  The blend integral is
    enclosed by lower/upper Riemann sums of psi' over NB cells (ball evaluation).  omega_W^2 = int_OR W <= W(57)/inf_OR psi'.
    Regression: the same code with F1-era quintic blend and W = (1+x^2)^-2 reproduces lemmaS_cflux.json W57_upper (lower or equal, close).
#4  profile balls on OR (step-1 tail theorem, step3/hstail.profile_bounds at sigma_min = 57 + STIP lower): k, m, p, q - sigma;
    C = 7km + q/4 per cell; k^2; sup|c0'|, c0' = -(5/2) m p + k^2 q; |Q_u| <= |q|/(1 - kU), |C_u| <= 1/(1 - kU).
#5  coefficient sups on OR from the certified C2 stage of the verifier (step3nl_v4b/v4cert/out/c2_t3d.json, tables/0.00025/ZB/sub/[57,60]:
    the code path of c2_t3d.py restricted to the boxes meeting (57, 60), boxes B-a (alpha <= 2.5e-4), B-c (|w'| <= 0.02), B-d (|w| <= 0.5));
    asserted <= the ZB values typed in constants_v4.json.  CBAR_OR = c, c_a, c_b; d010, d100, d200, d110, d300 with ^a, ^b; A0, A1, B0, B1;
    EA1 = sup |(E_A', E_B')|_2; kU.  Each raw value is rounded up by one more ulp.
#6  Gamma^(0)_ab and kappa^(0)_l: the reduced level-0 forms of X1 (out/flux2_forms.pkl) compiled to arb and evaluated on the NCELL cells with
    j in [1 - CBAR_OR, 1 + CBAR_OR], j' in [-C1_OR, C1_OR], C1_OR = c_a + c_b s1, C_eff = C + [-1, 1] abar (|Q_u| + |C_u|) as balls;
    off-diagonal Gamma_ab = sup |c_ab|/W, diagonal one-sided Gamma_aa = sup max(c_aa/W, 0); weighted AM-GM with the fixed rationals of F1
    (lambda_10 = 3/4, lambda_20 = 3/4, lambda_21 = 1/2): kappa_a += Gamma_ab lambda/2, kappa_b += Gamma_ab/(2 lambda).
    Regression: j = 1, j' = 0, abar = 0 reproduces F1's level-0 Gamma and kappa^(0) = (8.4244, 4.1592, 2.2596) (fix3/flux/flux_cert.json).
#7  s0 = S0V, s1 = S1V / 54.75 (S0V, S1V from the env of the fin13v4b rate tables; 54.75 <= inf of y, sigma, x on OR).
Level-2 table (sec. 5): pointwise row coefficients of frak X_pt (multipliers of N_0, N_1, N_2, wbar omega_W), x_3, and frak D coefficients.
abar: typed hypothesis ABAR_HYP = 1e-6 >= every bootstrap level in use (v36cfg abar 2.14041e-7; zone3 assert <= 2.19103e-7) [hyp].
Output out/flux2_coef.json."""
import sys, os, json, pickle, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail
T = FC.Timer()
NCELL = int(os.environ.get('F2_NCELL', '3000'))
NB = 8000
ABAR_HYP = arb('1e-6')
RT = {a: os.path.join(FC.V4B, 'runs', 'P', 'rate_table_fin13v4b_a%s.json' % a) for a in ('0.00015', '0.000125')}
C2 = os.path.join(FC.V4B, 'v4cert', 'out', 'c2_t3d.json')
CV4 = os.path.join(FC.V4B, 'v4cert', 'out', 'constants_v4.json')
OLDCF = '/home/claude/s1spec/writeup/lemmaS/lemmaS_cflux.json'
F1CERT = os.path.join(FC.F1, 'flux_cert.json')
FORMS = os.path.join(FC.OUT, 'flux2_forms.pkl')
R = {}
CHK = {}


def ok(name, val):
    CHK[name] = bool(val); print('%-78s %s' % (name, bool(val)), flush=True); assert val, name


# ------------------------------------------------------------------ #7 boxes from the rate tables
envs = [json.load(open(p))['args']['env'] for p in RT.values()]
ok('rate tables: env (FLUX0, FLUX2, S0V, S1V) = (1.64e-20, 8.0e-20, 0.0153, 0.6) at both alpha points',
   all((e['FLUX0'], e['FLUX2'], e['S0V'], e['S1V']) == ('1.64e-20', '8.0e-20', '0.0153', '0.6') for e in envs))
S0V = arb(envs[0]['S0V']); S1V = arb(envs[0]['S1V'])          # decimal strings -> balls containing the exact decimals
s0 = FC.aup(S0V); s1 = FC.aup(S1V / arb('54.75'))
R['s0'] = FC.rec(s0, 'cert', 'sup_OR |g| <= S0V (box B-e)'); R['s1'] = FC.rec(s1, 'cert', 'sup_OR |g\'| <= S1V/54.75 (box B-e; 54.75 <= inf y on OR)')
WBAR = arb('0.5')

# ------------------------------------------------------------------ #5 coefficient sups on [57, 60]
c2 = json.load(open(C2)); sub = c2['tables']['0.00025']['ZB']['sub']['[57,60]']; zb = c2['tables']['0.00025']['ZB']['sup']
cv4 = json.load(open(CV4))['constants']
KEYS = ['c', 'c_a', 'c_b', 'd010', 'd010_a', 'd010_b', 'd100', 'd100_a', 'd100_b', 'd200', 'd200_a', 'd200_b', 'd110', 'd110_a', 'd110_b',
        'd300', 'd300_a', 'd300_b', 'A0', 'A1', 'B0', 'B1', 'EA1', 'EA0', 'kU', 'DQ0', 'DQ1']
CO = {}
for k in KEYS:
    CO[k] = FC.from_raw_up(sub[k])
    typed = cv4['TQ[0.00025][ZB][%s]' % k]['value']
    ok('OR sup %-7s = %.6e <= ZB value %.6e (constants_v4.json)' % (k, sub[k], typed), sub[k] <= typed and sub[k] <= zb[k])
R['coef_OR'] = {k: FC.rec(v, 'cert', 'c2_t3d.json tables/0.00025/ZB/sub/[57,60], raw +1 ulp') for k, v in CO.items()}
CBAR = CO['c']; C1OR = FC.aup(CO['c_a'] + CO['c_b'] * s1)
R['CBAR_OR'] = FC.rec(CBAR, 'cert'); R['C1_OR'] = FC.rec(C1OR, 'cert', 'sup_OR |j\'| = sup |c\'| <= c_a + c_b s1')
JB = arb(1) + arb(0, CBAR.upper())             # j ball
JPB = arb(0, C1OR.upper())                      # j' ball
jbar = FC.aup(1 + CBAR)

# ------------------------------------------------------------------ #4 profile balls
STIP = hstail.STIP
sig_lo = FC.adn(arb(57) + STIP); sig_hi = FC.aup(arb(60) + STIP)
kb, mb, pb, qs = hstail.profile_bounds(sig_lo)
kb, mb, pb, qs = [arb(0, abs(a).upper()) for a in (kb, mb, pb, qs)]
kU = CO['kU']
ok('profile: kU < 1e-30 on OR (so 1/(1 - kU) <= 1 + 2e-30)', kU.upper() < 1e-30)
R['profile'] = dict(sigma_min=FC.rec_dn(sig_lo, 'cert'), sigma_max=FC.rec(sig_hi, 'cert'), kb=FC.rec(kb, 'cert'), mb=FC.rec(mb, 'cert'), pb=FC.rec(pb, 'cert'),
                    qs=FC.rec(qs, 'cert', '|q - sigma|'), STIP=FC.astr(STIP))
c0p = abs(-arb(5) / 2 * mb * pb + kb * kb * (arb(sig_hi) + qs))
R['sup_c0p'] = FC.rec(c0p, 'cert', 'sup_OR |c0\'|, c0\' = -(5/2) m p + k^2 q')

# ------------------------------------------------------------------ forms of X1 compiled
forms = pickle.load(open(FORMS, 'rb'))
x = sp.Symbol('x', real=True)
Wf, Zf, Cf, Kf, Jf = [sp.Function(n)(x) for n in ('W', 'Z', 'C', 'K', 'J')]
nW, nZ = 6, 8
rsym = sp.symbols('r0:%d' % nW); zsym = sp.symbols('z0:%d' % nZ)
Js, Jps, Cs, Ks = sp.symbols('Js Jps Cs Ks')


class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return "R(%d,%d)" % (e.p, e.q)
    def _print_Integer(self, e): return "arb(%d)" % e.p

    def _print_Pow(self, e, rational=False):
        if e.exp.is_Integer: return "IP(%s, %d)" % (self._print(e.base), int(e.exp))
        raise ValueError(str(e))


def R_(p, q): return arb(p) / arb(q)


def IP(b, n):
    r = arb(1)
    for _ in range(abs(n)): r = r * b
    return r if n >= 0 else 1 / r


def compile_expr(e):
    code = compile(_Pr().doprint(e), '<expr>', 'eval')
    return lambda env: eval(code, {'arb': arb, 'R': R_, 'IP': IP}, env)


def to_eval(e):
    rep = {}
    for d in e.atoms(sp.Derivative):
        n = d.derivative_count
        if d.expr == Wf: rep[d] = rsym[n]
        elif d.expr == Zf: rep[d] = zsym[n]
        elif d.expr == Jf: assert n == 1; rep[d] = Jps
        else: raise AssertionError('unexpected derivative %s' % d)
    e = e.xreplace(rep).xreplace({Wf: rsym[0], Zf: zsym[0], Cf: Cs, Kf: Ks, Jf: Js})
    return sp.expand(e)


ns = sp.__dict__.copy(); ns['x'] = x
L0 = {k: to_eval(eval(s, ns)) for k, s in forms['level0'].items()}
L0c = {k: compile_expr(e) for k, e in L0.items()}
# F1-regression versions: j = 1, j' = 0 (exact substitution, then expand: same canonical form as F1's to_eval)
L0r = {k: compile_expr(sp.expand(e.subs({Js: 1, Jps: 0}))) for k, e in L0.items()}
# W ratios r_n = W^(n)/W on x >= 38.25 (psi = phi + const), exact in x and S = (1+x^2)^(1/3)
psi = sp.Function('psi')(x); Ssym = sp.Symbol('S', positive=True); phi = (1 + x**2)**sp.Rational(2, 3) / 4


def wratio(n):
    e = sp.simplify(sp.diff(sp.exp(-psi), x, n) / sp.exp(-psi))
    for j in range(n, 0, -1): e = e.subs(sp.Derivative(psi, (x, j)), sp.diff(phi, x, j))
    e = sp.expand(e)
    e = e.replace(lambda t: t.is_Pow and t.base == 1 + x**2 and t.exp.is_Rational and (t.exp * 3).is_Integer, lambda t: Ssym**int(t.exp * 3))
    assert not e.has(sp.Pow(1 + x**2, sp.Rational(1, 3))), e
    return e


WR = [compile_expr(wratio(n)) for n in range(nW)]
uS = sp.Symbol('u'); S7 = uS**4 * (35 - 84 * uS + 70 * uS**2 - 20 * uS**3)
S7d = [compile_expr(sp.horner(sp.expand(sp.diff(S7, uS, i)), uS)) for i in range(nZ)]
ok('S7^(8) == 0 (mean-value form exact with nZ = 8)', sp.diff(S7, uS, 8) == 0)
# D2v coefficients (X1 (c)) with psi' = -r1, psi'' = -(r2 - r1^2)  (W'/W = -psi', W''/W = psi'^2 - psi'')
p1s, p2s = sp.symbols('p1s p2s')   # psi', psi2 = psi'^2 - psi''
D2C = {}
for i, srep in forms['D2v'].items():
    e = eval(srep, ns)
    rep = {}
    for d in e.atoms(sp.Derivative):
        if d.expr == Zf: rep[d] = zsym[d.derivative_count]
    e = e.xreplace(rep)
    e = e.xreplace({sp.Derivative(psi, (x, 2)): p1s**2 - p2s}).xreplace({sp.Derivative(psi, x): p1s}).xreplace({Zf: zsym[0]})
    D2C[int(i)] = compile_expr(sp.expand(e))


# ------------------------------------------------------------------ cells
def ball(lo, hi): return arb((lo + hi) / 2, (hi - lo) / 2)


sup = {}
def upd(name, val, onesided=False):
    cand = (arb(val.upper()) if val.upper() > 0 else arb(0)) if onesided else arb(abs(val).upper())
    if name not in sup or cand.upper() > sup[name].upper(): sup[name] = cand


inf_p1 = None; inf_Ceff = None
a0, a1 = arb(57), arb(60)
for i in range(NCELL):
    lo = a0 + (a1 - a0) * i / NCELL; hi = a0 + (a1 - a0) * (i + 1) / NCELL
    X = ball(lo, hi)
    Sx = (1 + X * X) ** (arb(1) / 3)
    env = {'x': X, 'S': Sx}
    rv = [f(env) for f in WR]
    u = (arb(60) - X) / 3
    ulo = max(0.0, float(u.lower())); uhi = min(1.0, float(u.upper()))
    uL = arb(math.nextafter(ulo, -math.inf)) if ulo > 0 else arb(0)
    uH = arb(math.nextafter(uhi, math.inf)) if uhi < 1 else arb(1)
    ub = uL.union(uH)
    assert ub.lower() <= max(u.lower(), arb(0)).lower() and ub.upper() >= min(u.upper(), arb(1)).upper()
    uc = arb(ub.mid())
    zv = []
    for j in range(nZ):
        sj = S7d[j]({'u': uc}) + (S7d[j + 1]({'u': ub}) * (ub - uc) if j + 1 < len(S7d) else arb(0))
        zv.append((arb(-1) / 3) ** j * sj)
    sig = X + STIP
    qb = sig + qs
    Cb = 7 * kb * mb + qb / 4                                 # C = 7 k m + q/4
    Kb = kb * kb
    mT = arb(0, (ABAR_HYP * ((abs(qb) + 1) / (1 - kU))).upper())   # |a Q_u + b C_u| <= abar (|Q_u| + |C_u|)
    Ceff = Cb + mT
    p1 = -rv[1]                                                  # psi' = -W'/W
    p2 = rv[2]                                                   # W''/W = psi'^2 - psi'' = psi2 exactly
    zabs = [abs(z) for z in zv]
    for j in range(5): upd('z%d' % j, zv[j])
    upd('psi1', p1); upd('psi2', p2)
    inf_p1 = p1 if inf_p1 is None or p1.lower() < inf_p1.lower() else inf_p1
    inf_Ceff = Ceff if inf_Ceff is None or Ceff.lower() < inf_Ceff.lower() else inf_Ceff
    envD = {'p1s': p1, 'p2s': p2, **{'z%d' % j: zv[j] for j in range(nZ)}}
    for k_, f in D2C.items(): upd('a%d' % k_, f(envD))
    upd('KAP2', 2 * zv[0] * zv[1] - p1 * zv[0] * zv[0])
    upd('X0lin', jbar * zabs[4] + abs(Ceff * zv[1]))           # sup (jbar |zeta''''| + |C_eff zeta'|)
    upd('k2', Kb)
    env2 = {**{'r%d' % j: rv[j] for j in range(nW)}, **{'z%d' % j: zv[j] for j in range(nZ)}, 'Js': JB, 'Jps': JPB, 'Cs': Ceff, 'Ks': Kb}
    env1 = {**{'r%d' % j: rv[j] for j in range(nW)}, **{'z%d' % j: zv[j] for j in range(nZ)}, 'Cs': Cb, 'Ks': Kb}
    for key in L0c:
        a, b = map(int, key.split(','))
        upd('G:' + key, L0c[key](env2), onesided=(a == b))
        upd('Greg:' + key, L0r[key](env1), onesided=(a == b))
    if i % 1000 == 0: print('  cell', i, '/', NCELL, '%.0fs' % T.s(), flush=True)

# the compiled level-0 coefficients already carry W through r0 = W/W = 1: check that r0 == 1
ok('W ratio r0 == 1 (coefficients c_ab/W are evaluated directly)', str(sp.simplify(wratio(0))) == '1')
Z = {j: FC.aup(sup['z%d' % j]) for j in range(5)}
ok('z_1 >= 35/48 and z_4 >= 280/27 (exact maxima; mean-value form encloses them)', Z[1].upper() >= arb(35) / 48 and Z[4].upper() >= arb(280) / 27)
ok('z_0 <= 1.001 (cell enclosure; exactly 0 <= zeta <= 1 since S7 is monotone on [0, 1])', Z[0].upper() <= 1.001)
R['z'] = {str(j): FC.rec(Z[j], 'cert') for j in range(5)}
AF = {j: FC.aup(sup['a%d' % j]) for j in range(5)}
ok('frak a_4 = sup |zeta| (coefficient of g\'\'\'\') <= 1.001 (exactly <= 1)', AF[4].upper() <= 1.001)
R['frak_a'] = {str(j): FC.rec(AF[j], 'cert', 'sup_OR |coefficient of g^(%d) in D2v|' % j) for j in range(4)}
KAP2 = FC.aup(sup['KAP2']); R['KAP2'] = FC.rec(KAP2, 'cert', 'sup_OR |2 zeta zeta\' - psi\' zeta^2|')
R['sup_psi1'] = FC.rec(sup['psi1'], 'cert'); R['inf_psi1'] = FC.rec_dn(inf_p1, 'cert'); R['sup_abs_psi2'] = FC.rec(sup['psi2'], 'cert')
R['inf_Ceff'] = FC.rec_dn(inf_Ceff, 'cert', 'C_eff >= this on OR (diagnostic: drift sign)')
ok('C_eff > 13 on OR (drift sign used once, inside the one-sided diagonal)', inf_Ceff.lower() > 13)
R['sup_k2'] = FC.rec(sup['k2'], 'cert')
R['X0lin'] = FC.rec(sup['X0lin'], 'cert', 'sup_OR (jbar |zeta\'\'\'\'| + |C_eff zeta\'|)')
GAM = {k[2:]: FC.aup(v) for k, v in sup.items() if k.startswith('G:')}
GREG = {k[5:]: FC.aup(v) for k, v in sup.items() if k.startswith('Greg:')}
R['Gamma0'] = {k: FC.rec(v, 'cert', 'one-sided' if k.split(',')[0] == k.split(',')[1] else 'sup |c_ab|/W') for k, v in GAM.items()}
LAM = {'1,0': FC.frac_arb('3/4'), '2,0': FC.frac_arb('3/4'), '2,1': FC.frac_arb('1/2')}


def kappa(Gm):
    k = [arb(0)] * 3
    for key, s_ in Gm.items():
        a, b = map(int, key.split(','))
        assert max(a, b) <= 2
        if a == b: k[a] += s_
        else: k[a] += s_ * LAM[key] / 2; k[b] += s_ / (2 * LAM[key])
    return [FC.aup(z) for z in k]


K0K = kappa(GAM); K0R = kappa(GREG)
R['kappa0'] = [FC.rec(z, 'cert', 'kappa^(0)_%d' % l) for l, z in enumerate(K0K)]
R['amgm_weights_level0'] = {k: str(v) for k, v in (('1,0', '3/4'), ('2,0', '3/4'), ('2,1', '1/2'))}
# regression vs F1
f1 = json.load(open(F1CERT))
reg = {}
for key, v in GREG.items():
    old = f1['Gamma_up']['0:' + key]
    reg[key] = dict(new=FC.fup(v), F1=old, rel=abs(FC.fup(v) / old - 1) if old else abs(FC.fup(v)))
regk = [dict(new=FC.fup(K0R[l]), F1=f1['kappa0_new'][l]) for l in range(3)]
ok('regression: j = 1, j\' = 0, abar = 0 reproduces F1 level-0 Gamma (rel <= 1e-12)', all(r['rel'] <= 1e-12 for r in reg.values()))
ok('regression: kappa^(0) reproduces F1 (8.4244, 4.1592, 2.2596) (rel <= 1e-12)', all(abs(r['new'] / r['F1'] - 1) <= 1e-12 for r in regk))
R['regression_F1'] = dict(Gamma=reg, kappa0=regk, bitwise_equal=sum(r['new'] == r['F1'] for r in reg.values()) + sum(r['new'] == r['F1'] for r in regk))
R['kappa0_vs_F1'] = [FC.fup(K0K[l] / K0R[l]) for l in range(3)]

# ------------------------------------------------------------------ #3 W(57)
def S7a(t): return t**4 * (35 - 84 * t + 70 * t * t - 20 * t**3)
def S5a(t): return t**3 * (10 - 15 * t + 6 * t * t)
def phi_a(xx): return (1 + xx * xx) ** (arb(2) / 3) / 4


def psi57(blend, nb):
    """(lower, upper) of psi(57) with psi = 2 log(1 + x^2) on x <= 30.25 and the given blend on [30.25, 38.25]"""
    a = arb('30.25'); b = arb('38.25'); base = 2 * (1 + a * a).log()
    lo = arb(0); hi = arb(0)
    for i in range(nb):
        xl = a + (b - a) * i / nb; xh = a + (b - a) * (i + 1) / nb; X = ball(xl, xh)
        t = (X - a) / 8
        s = blend(t)
        d = (1 - s) * 4 * X / (1 + X * X) + s * X / (3 * (1 + X * X) ** (arb(1) / 3))
        lo += arb(d.lower()) * (xh - xl); hi += arb(d.upper()) * (xh - xl)
    tail = phi_a(arb(57)) - phi_a(b)
    return base + lo + tail, base + hi + tail


def psi_lower_F1(xend, n=3000):
    """verbatim mirror of writeup/lemmaS/lemmaS_cflux_num.py psi_lower (quintic blend, one Riemann sum over [30.25, xend]); regression only"""
    a = arb(30.25); psi0 = 2*(1 + a*a).log()
    S = lambda t: t**3*(10 - 15*t + 6*t*t)
    tot = arb(0); L = arb(xend) - a
    def dpsi(xb):
        t = (xb - a)/8
        s = arb(1) if float(t.lower()) >= 1 else S(t)
        if float(t.lower()) < 1 and float(t.upper()) > 1: s = arb(0, 1) + arb(0.5)
        phip = xb/(3*(1 + xb*xb)**(arb(1)/3))
        return (1 - s)*4*xb/(1 + xb*xb) + s*phip
    for i in range(n):
        lo = a + L*i/n; hi = a + L*(i + 1)/n
        b = arb((lo + hi)/2, (hi - lo)/2); d = dpsi(b)
        tot += arb(float(d.lower()))*(hi - lo)
    return psi0 + tot


from flint import acb
def _blend_integrand(z, _an):
    qq = acb(1) + z * z
    q13 = (qq.log() * (-acb(1) / 3)).exp()          # (1 + z^2)^(-1/3), principal branch (Re(1 + z^2) > 0 near [30.25, 38.25])
    t = (z - acb('30.25')) / 8
    s = t**4 * (35 - 84 * t + 70 * t * t - 20 * t**3)
    return (1 - s) * 4 * z / qq + s * z * q13 / 3
_I = acb.integral(_blend_integrand, acb('30.25'), acb('38.25')).real
psi57_ball = 2 * (1 + arb('30.25')**2).log() + _I + phi_a(arb(57)) - phi_a(arb('38.25'))
pl, ph = psi57(S7a, NB)                               # Riemann lower/upper sums: containment cross-check only
ok('W(57): acb.integral enclosure of psi(57) lies inside the Riemann lower/upper sums', psi57_ball.lower() >= pl.lower() and psi57_ball.upper() <= ph.upper())
ok('W(57): acb.integral enclosure radius < 1e-30', float(psi57_ball.rad()) < 1e-30)
W57 = FC.aup((-psi57_ball).exp()); W57lo = FC.adn((-psi57_ball).exp())
W57old_typed = json.load(open(OLDCF))['W57_upper']
ctx.prec = 100                     # lemmaS_cflux_num.py:5
_pF1 = psi_lower_F1(57)
_W57F1 = float(arb(-_pF1).exp().upper())
ctx.prec = 200
ok('W(57) regression: mirror of the F1-era psi_lower reproduces lemmaS_cflux W57_upper bitwise', _W57F1 == W57old_typed)
ol, oh = psi57(S5a, NB)
W57old = FC.aup((-ol).exp())
ok('W(57) method check: same quintic definition, finer sums: W57 <= the F1 value', W57old.upper() <= W57old_typed)
R['W57'] = FC.rec(W57, 'cert', 'W(57) = exp(-psi(57)), verifier psi (S7 blend, W = (1+x^2)^-2 on x <= 30.25)')
R['W57_lower'] = FC.rec_dn(W57lo, 'cert'); R['W57_quintic_def'] = FC.rec(W57old, 'cert', 'F1-era quintic-blend definition with this code (finer sums); regression only')
R['W57_ratio_S7_over_F1typed'] = FC.fup(W57 / arb(W57old_typed)); R['W57_upper_typed_old'] = W57old_typed
omW2 = FC.aup(W57 / FC.adn(inf_p1)); omW = FC.aup(omW2.sqrt())
R['omega_W2'] = FC.rec(omW2, 'cert', 'int_OR W <= W(57)/inf psi\''); R['omega_W'] = FC.rec(omW, 'cert')
ok('W(57) <= 1.6e-18', W57.upper() <= 1.6e-18)

# ------------------------------------------------------------------ level-2 table coefficients (sec. 5)
d = CO
s1a = s1
X2 = {}                       # multipliers of N_2, N_1, N_0, wbar*omega_W (pointwise part frak X_pt), each [cert]
X2['N2'] = (jbar * 6 * Z[2] + C1OR * 3 * Z[1] + (d['d010_a'] + d['d010_b'] * s1a + d['d100']) + 2 * d['d200'] * s1a
            + (d['d110_a'] + d['d110_b'] * s1a) * s1a + 3 * d['d300'] * s1a**2)
X2['N1'] = (jbar * 4 * Z[3] + C1OR * 3 * Z[2] + sup['k2'] * 2 * Z[1] + (d['d100_a'] + d['d100_b'] * s1a) + (d['d200_a'] + d['d200_b'] * s1a) * s1a
            + (d['d300_a'] + d['d300_b'] * s1a) * s1a**2 + d['A1'])
X2['N0'] = sup['X0lin'] + C1OR * Z[3] + sup['k2'] * Z[2] + d['B1']
X2['w'] = d['B0'] + c0p + ABAR_HYP * d['EA1']
X2 = {k: FC.aup(v) for k, v in X2.items()}
x3 = FC.aup(4 * Z[1] * jbar + d['d010'] + d['d110'] * s1a)
R['X2_pt_coef'] = {k: FC.rec(v, 'cert', 'frak X_pt = sum coef * (N_2, N_1, N_0, wbar omega_W)') for k, v in X2.items()}
R['x3'] = FC.rec(x3, 'cert', 'x_3 = 4 z_1 jbar + d010 + d110 s1')
R['d110'] = FC.rec(d['d110'], 'cert')
# level-0 nonlinear and w coefficients
R['K0nl_coef'] = dict(N2=FC.rec(d['d010'] + d['d110'] * s1a, 'cert'), N1=FC.rec(d['d100'] + d['d200'] * s1a + d['d300'] * s1a**2, 'cert'),
                      N0=FC.rec(d['A1'], 'cert'), w=FC.rec(d['A0'], 'cert'), note='K0nl = (N_1 + KAP2 N_0) (coef.N2 N_2 + coef.N1 N_1 + coef.N0 N_0 + coef.w wbar omega_W)')
R['K0w_coef'] = FC.rec(c0p + ABAR_HYP * d['EA1'], 'cert', 'K0w = N_0 wbar omega_W * this')
R['wbar'] = '0.5'; R['ABAR_HYP'] = dict(value='1e-6', label='hyp', note='typed modulation level; >= 2.19103e-7 (v36 zone3 assert) and 2.14041e-7 (v36cfg)')
R['jbar'] = FC.rec(jbar, 'cert')
out = dict(stage='X2 flux2_coef', ncell=NCELL, nb_psi=NB, prec=ctx.prec, res=R, checks=CHK, all_pass=all(CHK.values()),
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), FORMS, os.path.join(FC.OUT, 'flux2_forms.json'), C2, CV4, OLDCF, F1CERT,
                               hstail.__file__, '/home/claude/s1spec/lin/lintail2.py'] + list(RT.values())))
out['timing'] = dict(runtime_s=T.s())
h = FC.dump(os.path.join(FC.OUT, 'flux2_coef.json'), out)
print('kappa0', [FC.fup(z) for z in K0K], 'F1-regression', [FC.fup(z) for z in K0R])
print('z', [FC.fup(Z[j]) for j in range(5)], 'frak a', [FC.fup(AF[j]) for j in range(4)], 'KAP2', FC.fup(KAP2))
print('W57', FC.fup(W57), 'omega_W', FC.fup(omW), 'x3', FC.fup(x3), 'X2 coef', {k: FC.fup(v) for k, v in X2.items()})
print('X2 done sha', h, '(%.1fs)' % T.s())
