#!/usr/bin/env python3
"""Stage X3 (FLUX2_DESIGN.md sec. 3.4, constants #8): the chart identity at level 4 on the shell S, and the chart/profile data on S for levels 0..4.
Notation and construction exactly as fix6_v36/n1/n1_chart.py (= step4/p2/v33/n1_chart.py with v36 wiring) and the chart lemma
(fix1/b12_text/chart_lemma_body.tex): Y = q cos d - (p + u) sin d, e = Y - sigma, eta = 2 sin^2(d/2) + b sin d + u k cos d, Y' = 1 - eta,
f(Y) = b + Phi, T = (1 + eps) D, 1 + eps = (1 + e/sigma)/(1 - eta), chi_k = (T^k - D^k) b + T^k Phi + p^k (U(sigma) - U(Y)); b = u', b_i = d^i b.

[sym] Built in sympy for k = 0..4 (n1_chart builds k <= 3).  ASSERTED for k = 4 (the stop rule of the task):
   (i)  chi_4 contains no b_5 (= u^(6)) and no unresolved derivative;
   (ii) chi_4 is a polynomial of degree <= 2 in (b_2, b_3, b_4) whose ONLY quadratic monomial is b_2^2 (i.e. (u''')^2):
        d^2/db2 db3 = d^2/db2 db4 = d^2/db3^2 = d^2/db3 db4 = d^2/db4^2 = 0 and d^3/db2^3 = 0, all exactly;
   (iii) chi_4 = 0 for a ray (d_i = 0, e = 0);  (iv) the k <= 3 structure of n1_chart (affine in b2, b3; B_02 = B_03 = B_12 = B_13 = B_23 = 0)
        is re-asserted here.  So chi_4 = A_4 + B_42 b_2 + B_43 b_3 + B_44 b_4 + B_422 b_2^2 with A, B functions of (sigma, e, d_0..d_5, u, b, b_1).
[cert] arb on cells of the shell sigma in [54.75, 57.76] (contains S = [57 + STIP, 60 + STIP]) with the boxes
   |alpha| <= 2.5e-4 (B-a; P = alpha psi3^' + alpha^2 h2^', c4fun jets to order 6), |w| <= 1/2 (B-d), |w'| <= S0V (B-e), |w''| <= S1V/54.75 (B-e),
   profile envelopes of hstail.profile_bounds and |d| <= kb/mu (n1_chart), |k''''| <= (mb (s1 + qs) + kb (1 + kb pb))/4:
     Abar_k := sup |A_k| + p^k sup|U(sigma) - U(Y)|, Bbar_ki := sup |B_ki|, Bbar_422 := sup |B_422|  (k = 0..4);
     sup |P^(i)| and int_S (P^(i))^2 ds (i = 0..4), sup|P''^2| ; sup |D^k d|, int_S (D^k d)^2 ds (k = 0..4), d = P - U(sigma), U = 1.3056 alpha sigma^1.0246;
     eps_c, e on S, rho_meas = max_{m <= 8} Y' (sigma/Y)^m - 1 (weight transfer ds <-> dy/y with any power y^{-m}, m <= 8).
   Regression: the k <= 3 sups on S are <= the n1_chart core values (fix6_v36/n1/out/n1_chart_fin13v4b.json if present, else step4 out_v33).
[float, not part of the proof] spot check of the level-4 identity z^A_4 o Y = D^4(w' + d) + chi_4 on a synthetic curve (finite differences).
Output out/chart4_shell.json."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_int/v3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
sys.path.insert(0, '/home/claude/s1spec/step3')
import sympy as sp
from flint import arb, ctx
import s1chart
import hstail as HT
import lintail2 as LT
import c4fun as C4
import resource
ctx.prec = 128
T = FC.Timer()
CHK = {}


def ok(name, val):
    CHK[name] = bool(val); print('%-80s %s' % (name, bool(val)), flush=True); assert val, name


def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


STIP = s1chart.SIG_TIP
P = arb('1.0246'); C13 = arb('1.3056'); AL = arb('2.5e-4')
RTF = os.path.join(FC.V4B, 'runs', 'P', 'rate_table_%s_a0.00015.json' % FC.VTAG)
env = json.load(open(RTF))['args']['env']
S0V = FC.aup(arb(env['S0V'])); S1V = arb(env['S1V'])
WSUP = arb('0.5'); W1SUP = S0V; W2SUP = FC.aup(S1V / arb('54.75'))
NMAX = 4

# ---------------------------------------------------------------- symbolic chi_k, k = 0..4 (n1_chart construction, extended)
s = sp.Symbol('s', positive=True)
dlf = sp.Function('dl')(s); uf = sp.Function('u')(s); ef = sp.Function('e')(s); Ef = sp.Function('E')(s)
kf = dlf.diff(s); bf = uf.diff(s)
ETA = 2 * sp.sin(dlf / 2) ** 2 + bf * sp.sin(dlf) + uf * kf * sp.cos(dlf)
EPS = (ef / s + ETA) / (1 - ETA)
NUM = (1 - kf * uf + bf ** 2) * sp.sin(dlf) + kf * uf * bf * sp.cos(dlf); DEN = (1 - kf * uf) * sp.cos(dlf) - bf * sp.sin(dlf); PHI = NUM / DEN
NJ = 7
DS = sp.symbols('d0:%d' % NJ); US = sp.symbols('u0 b0 b1 b2 b3 b4 b5'); E0 = sp.Symbol('e0')


def to_syms(ex, K):
    for i in range(K, -1, -1): ex = ex.subs(Ef.diff(s, i) if i else Ef, EPS.diff(s, i) if i else EPS)
    for n in range(NJ, 0, -1): ex = ex.subs(ef.diff(s, n), -ETA.diff(s, n - 1))
    for n in range(NJ - 1, -1, -1): ex = ex.subs(uf.diff(s, n) if n else uf, US[n])
    for n in range(NJ - 1, -1, -1): ex = ex.subs(dlf.diff(s, n) if n else dlf, DS[n])
    return ex.subs(ef, E0)


Top = lambda F: (1 + Ef) * s * F.diff(s); Dop = lambda F: s * F.diff(s)
CHI = {}; Fb = bf; Fp = PHI; Db = bf
for kk in range(NMAX + 1):
    if kk: Fb, Fp, Db = Top(Fb), Top(Fp), Dop(Db)
    CHI[kk] = to_syms(sp.expand(Fb - Db), kk) + to_syms(Fp, kk)
    print('chi_%d built (%.1fs, %d MB)' % (kk, T.s(), rss()), flush=True)
b2, b3, b4, b5 = US[3], US[4], US[5], US[6]
ZERO = {**{DS[i]: 0 for i in range(NJ)}, E0: 0}
for kk in range(NMAX + 1):
    c = CHI[kk]
    ok('chi_%d: no unresolved function/derivative' % kk, not c.atoms(sp.Derivative) and not any(isinstance(f, sp.core.function.AppliedUndef) for f in c.atoms(sp.Function)))
    ok('chi_%d: no b_5 (u^(6))%s' % (kk, '' if kk == 4 else ''), b5 not in c.free_symbols)
    if kk <= 3: ok('chi_%d: no b_4' % kk, b4 not in c.free_symbols)
    if kk <= 3: ok('chi_%d: no k\'\'\'\' (d_5)' % kk, DS[5] not in c.free_symbols)


def zero(e):
    if e == 0: return True
    e2 = sp.expand(e)
    if e2 == 0: return True
    return sp.simplify(e2) == 0


for kk in range(NMAX + 1):
    c = CHI[kk]
    pairs = [(b2, b3), (b2, b4), (b3, b3), (b3, b4), (b4, b4)] if kk == 4 else [(b2, b2), (b2, b3), (b3, b3)]
    for (p_, q_) in pairs:
        ok('chi_%d: d^2/d%s d%s == 0' % (kk, p_, q_), zero(sp.diff(c, p_, q_)))
    if kk == 4:
        ok('chi_4: d^3/db2^3 == 0 (degree <= 2 in b2)', zero(sp.diff(c, b2, 3)))
    ok('chi_%d == 0 for a ray (d_i = 0, e = 0)' % kk, zero(c.subs(ZERO)))
    print('  structure k=%d done (%.1fs, %d MB)' % (kk, T.s(), rss()), flush=True)
for kk, bb in ((0, b2), (0, b3), (1, b2), (1, b3), (2, b3)):
    ok('chi_%d: B_%d%s == 0 (n1_chart lemma (e))' % (kk, kk, str(bb)[1]), zero(sp.diff(CHI[kk], bb)))
SYM = {}
for kk in range(NMAX + 1):
    c = CHI[kk]
    A = c.subs({b2: 0, b3: 0, b4: 0})
    Bs = {2: sp.diff(c, b2).subs({b2: 0}), 3: sp.diff(c, b3), 4: sp.diff(c, b4)}
    B22q = sp.diff(c, b2, 2) / 2
    SYM[kk] = dict(A=A, B=Bs, Q=B22q)
# quadratic coefficient: only k = 4 may have it
for kk in range(4):
    ok('chi_%d: no b2^2 term' % kk, zero(SYM[kk]['Q']))
# the B coefficients and A depend only on (s, e0, d_0..d_5, u0, b0, b1)
ALLOW = {s, E0, *DS[:6], US[0], US[1], US[2]}
for kk in range(NMAX + 1):
    fs = set(SYM[kk]['A'].free_symbols)
    for e in list(SYM[kk]['B'].values()) + [SYM[kk]['Q']]: fs |= e.free_symbols
    ok('chi_%d: A, B, Q depend on (sigma, e, d_0..d_5, u, b, b_1) only' % kk, fs <= ALLOW)
print('symbolic part done (%.1fs, %d MB)' % (T.s(), rss()), flush=True)

# ---------------------------------------------------------------- arb evaluation
IP = sp.Function('ipow')


def ipow(a, n):
    n = int(n)
    if n < 0: return 1 / ipow(a, -n)
    r = arb(1)
    for _ in range(n): r = r * a
    return r


fixpow = lambda ex: ex.replace(lambda q: q.is_Pow and q.exp.is_Integer, lambda q: IP(q.base, q.exp))
MODL = [{'sin': lambda a: a.sin(), 'cos': lambda a: a.cos(), 'ipow': ipow}]
ARGS = (s, E0, *DS[:6], US[0], US[1], US[2])
LAM = lambda ex: sp.lambdify(ARGS, fixpow(sp.sympify(ex)), modules=MODL)
FA = {kk: LAM(SYM[kk]['A']) for kk in range(NMAX + 1)}
BKEYS = [(2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4)]          # all other B_ki vanish (asserted above: B_02, B_03, B_12, B_13, B_23; b4 absent for k <= 3)
FB = {(kk, i): LAM(SYM[kk]['B'][i]) for (kk, i) in BKEYS}
FQ = {4: LAM(SYM[4]['Q'])}
ETA_S = ETA.subs({bf: US[1], uf: US[0], kf: DS[1]}).subs(dlf, DS[0])
FETA = LAM(ETA_S); FEPS = LAM((E0 / s + ETA_S) / (1 - ETA_S))
print('lambdified: B keys', sorted(FB), '(%.1fs, %d MB)' % (T.s(), rss()), flush=True)
UPf = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
ev = lambda f, a: arb(UPf(arb(0) + f(*a)))
STIR = {0: {}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}, 4: {1: 1, 2: 7, 3: 6, 4: 1}}


def boxes(s0, s1):
    S = s0.union(s1)
    kb, mb, pb, qs = HT.profile_bounds(s0); mu = LT.mu(s0)
    dlb = kb / mu
    X = S - STIP; Jp = C4.PSI.jets_any(X, float(X.mid()), 7); Jh = C4.H2.jets_any(X, float(X.mid()), 7)
    pj = lambda j: AL * abs(Jp[j]) + AL * AL * abs(Jh[j])
    Ub = WSUP + pj(0); B0b = W1SUP + pj(1); B1b = W2SUP + pj(2)
    sn = (dlb / 2).sin()
    eb = qs + 2 * s1 * sn * sn + (pb + Ub) * dlb
    k4b = (mb * (s1 + qs) + kb * (1 + kb * pb)) / 4
    bx = lambda v: arb(0, UPf(v))
    a = (S, bx(eb), bx(dlb), bx(kb), bx(mb), bx(pb / 4), bx(kb * (s1 + qs) / 4), bx(k4b), bx(Ub), bx(B0b), bx(B1b))
    Pj = [pj(i + 1) for i in range(6)]                                    # |P^(i)|, P = alpha psi3' + alpha^2 h2'
    # D^k d, d = P - U(sigma), U = C13 alpha sigma^P
    pw = lambda m: math.prod([P - i for i in range(m)], start=arb(1)) * S ** (P - m)
    Dk = lambda J, k: (J[1] if k == 0 else sum((arb(c) * S ** m * J[m + 1] for m, c in STIR[k].items()), arb(0)))
    dD = [AL * abs(Dk(Jp, k) - C13 * P ** k * S ** P) + AL * AL * abs(Dk(Jh, k)) for k in range(5)]
    return dict(S=S, a=a, e=arb(UPf(eb)), Pj=Pj, dD=dD)


def cell(s0, s1):
    bx = boxes(s0, s1); a = bx['a']; S = bx['S']
    A = [ev(FA[k], a) for k in range(NMAX + 1)]
    Bv = {key: ev(f, a) for key, f in FB.items()}
    Q4 = ev(FQ[4], a)
    eta = ev(FETA, a); eps = ev(FEPS, a); e = bx['e']
    dU = e * C13 * AL * P * (s1 + e) ** (P - 1)                           # |U(sigma) - U(sigma + e)|
    Abar = [arb(UPf(A[k] + P ** k * dU)) for k in range(NMAX + 1)]
    Yp = 1 + arb(0, eta.upper()); E = arb(0, e.upper()); r = 1 / (1 + E / S)
    meas = max((Yp * r ** m for m in range(9)), key=lambda v: float(v.upper()))
    meas2 = max((r ** (-m) / (1 - arb(0, eta.upper())) for m in range(9)), key=lambda v: float(v.upper()))   # reverse direction
    return dict(s0=s0, s1=s1, Abar=Abar, B=Bv, Q4=Q4, eta=eta, eps=eps, e=e, meas=meas, meas2=meas2, Pj=bx['Pj'], dD=bx['dD'])


NC = 49
SLO, SHI = arb('54.75'), arb('57.76')
ok('shell [54.75, 57.76] contains S = [57 + STIP, 60 + STIP]', (arb(57) + STIP).lower() >= SLO.lower() and (arb(60) + STIP).upper() <= SHI.upper())
edges = [SLO + (SHI - SLO) * i / NC for i in range(NC + 1)]
cells = [cell(edges[i], edges[i + 1]) for i in range(NC)]
mx = lambda L: arb(UPf(max(L, key=lambda v: float(v.upper()))))
Abar = [mx([c['Abar'][k] for c in cells]) for k in range(NMAX + 1)]
Bbar = {'%d%d' % key: mx([c['B'][key] for c in cells]) for key in FB}
Q4 = mx([c['Q4'] for c in cells])
Pbar = [mx([c['Pj'][i] for c in cells]) for i in range(5)]
dbar = [mx([c['dD'][k] for c in cells]) for k in range(5)]
ds_int = lambda f: arb(UPf(sum((f(c) * ((c['s1'] / c['s0']).log()) for c in cells), arb(0))))      # int_S f ds, ds = dsigma/sigma
Pint = [ds_int(lambda c, i=i: c['Pj'][i] ** 2) for i in range(5)]
P2sq_int = ds_int(lambda c: c['Pj'][2] ** 4)
dint = [ds_int(lambda c, k=k: c['dD'][k] ** 2) for k in range(5)]
Sds = arb(UPf((SHI / SLO).log()))
rho_meas = arb(UPf(max(mx([c['meas'] for c in cells]), mx([c['meas2'] for c in cells])) - 1))
eps_c = mx([c['eps'] for c in cells]); e_sup = mx([c['e'] for c in cells])
print('Abar', ['%.3e' % float(v.upper()) for v in Abar], 'Bbar', {k: '%.3e' % float(v.upper()) for k, v in Bbar.items()}, 'Q4 %.3e' % float(Q4.upper()), flush=True)
print('Pbar', ['%.3e' % float(v.upper()) for v in Pbar], 'dbar', ['%.3e' % float(v.upper()) for v in dbar], 'rho_meas %.3e' % float(rho_meas.upper()), flush=True)
# regression: k <= 3 sups on S <= n1_chart core values
n1p = os.path.join(FC.N36, 'n1', 'out', 'n1_chart_%s.json' % FC.VTAG)
if not os.path.exists(n1p): n1p = '/home/claude/s1spec/step4/p2/out_v33/n1_chart.json'
n1 = json.load(open(n1p))['core']
ok('regression: Abar_k (k <= 3, chart part incl. U-shift) on S <= n1_chart core sup_add_k', all(float(Abar[k].upper()) <= n1['sup_add'][k] for k in range(4)))
ok('regression: chart e on S <= n1_chart core e', float(e_sup.upper()) <= n1['e'])
ok('chart4 coefficients small: max(Abar_4, Bbar_4i, Q4) <= 1e-30', max(float(Abar[4].upper()), max(float(v.upper()) for k, v in Bbar.items() if k[0] == '4'), float(Q4.upper())) <= 1e-30)

# ---------------------------------------------------------------- float spot check of the level-4 identity (not part of the proof)
def spot():
    """geometric route: Gamma = Gamma1 + u i T1, T1 = e^{i theta1}, theta1 = Theta + delta  =>  e^{-i Theta} Gamma' = e^{i delta}((1 - k u) + i b);
    f = Im/Re of that (the slope dh/dy), y' = its real part; T = (y/y') d/dsigma with y a symbol (y' known in closed form).
    Synthetic delta(s) = 0.02 exp(-s/7), u(s) = 0.003 sin(s/3) + 1e-4 s; U = 0.  Checks T^4 f - D^4 b == chi_4 at s0 = 31, e0 = 1e-3."""
    import mpmath as mp
    mp.mp.dps = 40
    Ysym = sp.Symbol('Ysym')
    dls = sp.Rational(2, 100) * sp.exp(-s / 7); us = sp.Rational(3, 1000) * sp.sin(s / 3) + sp.Rational(1, 10000) * s
    ks = sp.diff(dls, s); bs = sp.diff(us, s)
    re_ = sp.cos(dls) * (1 - ks * us) - sp.sin(dls) * bs; im_ = sp.sin(dls) * (1 - ks * us) + sp.cos(dls) * bs
    yp = re_; f = im_ / re_
    Tg = lambda G: Ysym / yp * (sp.diff(G, s) + sp.diff(G, Ysym) * yp)
    G = f
    for _ in range(4): G = Tg(G)
    Db = bs
    for _ in range(4): Db = s * sp.diff(Db, s)
    s0_ = sp.Integer(31); e0v = sp.Rational(1, 1000)
    lhs = sp.N((G - Db).subs({Ysym: s0_ + e0v}).subs(s, s0_), 35)
    vals = {s: s0_, E0: e0v}
    for i in range(6): vals[DS[i]] = sp.N(sp.diff(dls, s, i).subs(s, s0_), 40)
    for i in range(6): vals[US[i]] = sp.N(sp.diff(us, s, i).subs(s, s0_), 40)
    fchi = sp.lambdify(list(vals.keys()), CHI[4], modules='mpmath')
    rhs = fchi(*[mp.mpf(str(v)) for v in vals.values()])
    return float(lhs), float(rhs)


try:
    l_, r_ = spot()
    spot_rel = abs(l_ - r_) / max(abs(l_), 1e-300)
except Exception as ex:                       # the spot check is diagnostic only
    l_, r_, spot_rel = None, None, 'failed: %s' % ex
print('float spot check level 4 (diagnostic):', l_, r_, spot_rel, flush=True)

R = dict(shell=['54.75', '57.76'], cells=NC, boxes=dict(alpha='2.5e-4', w='0.5', w1=FC.fup(W1SUP), w2=FC.fup(W2SUP), note='standing boxes B-a, B-d, B-e of %s' % FC.VTAG),
         Abar=[FC.rec(v, 'cert', 'sup_S |A_%d| + p^%d sup |U(sigma) - U(Y)|' % (k, k)) for k, v in enumerate(Abar)],
         Bbar={k: FC.rec(v, 'cert', 'sup_S |B_%s| (coefficient of b_%s in chi_%s)' % (k, k[1], k[0])) for k, v in Bbar.items()},
         B422=FC.rec(Q4, 'cert', 'sup_S |B_422| (coefficient of b_2^2 in chi_4)'),
         Pbar=[FC.rec(v, 'cert', 'sup_S |P^(%d)|' % i) for i, v in enumerate(Pbar)],
         P_int_ds=[FC.rec(v, 'cert', 'int_S (P^(%d))^2 ds' % i) for i, v in enumerate(Pint)], P2sq_int_ds=FC.rec(P2sq_int, 'cert', 'int_S (P\'\')^4 ds'),
         Dd_sup=[FC.rec(v, 'cert', 'sup_S |D^%d d|' % k) for k, v in enumerate(dbar)], Dd_int_ds=[FC.rec(v, 'cert', 'int_S (D^%d d)^2 ds' % k) for k, v in enumerate(dint)],
         S_ds=FC.rec(Sds, 'cert', '|S|_s = log(57.76/54.75)'), rho_meas=FC.rec(rho_meas, 'cert', 'weight transfer ds <-> dy/y incl. y^-m, m <= 8, both directions'),
         eps_c=FC.rec(eps_c, 'cert'), e=FC.rec(e_sup, 'cert'),
         spot_check=dict(lhs=l_, rhs=r_, rel=spot_rel, label='float diagonstic, not part of the proof'),
         n1_core_reference=dict(file=n1p, sup_add=n1['sup_add'], e=n1['e']))
out = dict(stage='X3 chart4_shell', res=R, checks=CHK, all_pass=all(CHK.values()), prec=ctx.prec,
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), RTF, n1p, s1chart.__file__, HT.__file__, LT.__file__, C4.__file__,
                               '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py']))
out['timing'] = dict(runtime_s=T.s(), maxrss_MB=rss())
h = FC.dump(os.path.join(FC.OUT, 'chart4_shell.json'), out)
print('X3 done sha', h, '(%.1fs, %d MB)' % (T.s(), rss()))
