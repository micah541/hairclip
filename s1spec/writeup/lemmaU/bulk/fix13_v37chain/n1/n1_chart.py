"""v3.3 re-review N1 (REREVIEW_v33.md): certified bounds on the chart terms between the S1-chart variable of the verifier and the
ray-graph variable of Lemma S(A), on (A)'s ramp sigma in [14 - WS, 28 + WS] and on the core sigma in [28, 60].  Lemma: writeup/lemmaU/bulk/
fix2/n12/chart_lemma.tex.  Labels: every number written is [cert] (python-flint arb, outward rounding) unless its key says 'hyp'.
Setting (rescaled frame, fixed tau).  S1: Gamma1(sigma), arclength sigma (CAPD convention, x = sigma - SIG_TIP), T1 = e^{i theta1},
theta1' = k, k' = m, m' = p/4, p' = -k q, q' = 1 + k p, Gamma1 = (q + i p) e^{i theta1}; Theta = lim theta1 (sigma -> oo);
delta := theta1 - Theta = -int_sigma^oo k.  Curve: Gamma = Gamma1 + u N1, N1 = i T1 (step3nl graph_eq.py), u = alpha psi3^ + alpha^2 h2^ + w,
b := u' = w' + P, P := alpha psi3^' + alpha^2 h2^'.  Ray chart (A): y = Re(e^{-i Theta} Gamma), h = Im(e^{-i Theta} Gamma), f = dh/dy.
Exact identities (chart_lemma.tex, Lemma C1):  Y(sigma) := y at the curve point sigma = q cos d - p sin d - u sin d  (d = delta),
  Y' = 1 - eta,  eta = 2 sin^2(d/2) + b sin d + u k cos d;   e := Y - sigma = (q - sigma) cos d - 2 sigma sin^2(d/2) - (p + u) sin d;
  f(Y(sigma)) = b + Phi,  Phi = [(1 - k u + b^2) sin d + k u b cos d] / [(1 - k u) cos d - b sin d];
  with T := y d/dy = (1 + eps) D,  D = sigma d/dsigma,  1 + eps = Y/(sigma Y') = (1 + e/sigma)/(1 - eta):
  z^A_k o Y = D^k (w' + d_B0) + chi_k,   chi_k = (T^k - D^k) b + T^k Phi + p^k (U(sigma) - U(Y))     (d_B0 = P - U(sigma), stage b0_dprof).
chi_k (k <= 3) is computed symbolically (sympy) as A_k + B_k2 b'' ' + B_k3 b''''  (x-derivatives: b2 = u''', b3 = u''''; affine, asserted),
all coefficients carrying explicit small factors (sin d, k, k', k'', k''', eta, e), evaluated in arb on cells with boxes:
  |d| <= kb(s0)/mu(s0) (int_s^oo kb <= kb(s)/mu(s): E convex, E' = mu), |k| <= kb, |k'| <= mb, |k''| <= pb/4, |k'''| = |k q|/4 <= kb (s1 + qs)/4,
  (kb, mb, pb, qs) = step3/hstail.profile_bounds(s0) (step-1 tail theorem, decreasing envelopes, valid for sigma >= s0);
  |e| <= qs + 2 s1 sin^2(d/2) + (pb + |u|)|d|;  |u|, |b|, |b'| from the verifier fin12f4 sups (alpha_1 = 1.5e-4 row: max(far_w, S_ramp),
  max(sI1, S_ramp), sI2) plus |alpha| <= alpha_1 [hyp] times the certified C^4 profile jets (step3nl/agent_int/v3/c4fun.py).
Split (D-form): b2 = (D^2 w' - D w')/sigma^2 + P'', b3 = (D^3 w' - 3 D^2 w' + 2 D w')/sigma^3 + P''', |D w'| = sigma |w''| <= sigma sI2:
  chi_k = chi_k^add + beta_k2 D^2 w' + beta_k3 D^3 w',  |chi_k^add| <= add_k := |A_k| + |B_k2||P''| + |B_k3||P'''| + (|B_k2|/s^2 + 2|B_k3|/s^3) s sI2
  + p^k |U(s) - U(Y)|,   |beta_k2| <= |B_k2|/s^2 + 3|B_k3|/s^3,  |beta_k3| <= |B_k3|/s^3.  For Phi of (A): beta22 = B_22/s^2 (signed) and D beta22.
Measure: (A)'s ramp weight Theta(y) dy/y at y = Y(sigma) against the verifier's W(x) dx: factor <= 1 + rho_meas (max over m <= 3 of
  Y' (sigma/Y)^{2m} W(Y + X0)/W(sigma + X0)); core (28/y) dy/y: Y' max_m (sigma/Y)^m, m <= 6.
Also: d_B0 on the ramp end slivers [14 - WS, 14], [28, 28 + WS] (the ramp sigma-domain Y^{-1}([14, 28]) lies in [14 - WS, 28 + WS], asserted
e <= WS), which b0_dprof's [14, 28] tables do not cover.
Writes p2/out_v33/n1_chart.json.   usage: python3 -I -B n1_chart.py   [--out PATH: test runs only]"""
import sys, os, json, math, time, hashlib
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_int/v3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
sys.path.insert(0, '/home/claude/s1spec/step3')
import sympy as sp
from flint import arb, ctx
import s1chart
import hstail as HT
import lintail2 as LT
import c4fun as C4
ctx.prec = 128
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33'
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain'); import v36cfg as CFG   # fix6_v36 (v36): verifier source + switches V36_VER/V36_R3/V36_R2
VER = CFG.RT                                         # v36: verifier per V36_VER (was fin12f4)
MODS = ['/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
        f'{P2}/t9h/s1chart.py', '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/lin/lintail2.py']
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
UPf = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
UPd = lambda x: math.nextafter(float(x), math.inf)
STIP = s1chart.SIG_TIP
P = arb('1.0246'); C13 = arb('1.3056'); AL1 = arb('1.5e-4')
X0 = arb('2.2490776515').union(-STIP)
WS = arb(2) ** -19                                  # ramp end slivers (exact dyadic, 1.9073e-6)
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}

# ---------------------------------------------------------------- verifier sups (alpha_1 = 1.5e-4), outward (one ulp up)
V = json.load(open(VER)); assert V['alpha1'] == 0.00015
FM = V['far_max']; SR = V['S_ramp_max_along_trajectory']
WSUP = arb(UPd(max(FM['far_w'], SR)))               # sup |w|   on x >= 10.25 (far_w: x >= 16.25; S_ramp: [10.25, 16.25])
W1SUP = arb(UPd(max(FM['sI1'], SR)))                # sup |w'|  on x >= 10.25
W2SUP = arb(UPd(FM['sI2']))                         # sup |w''| on x >= 16.25 (on x in [16.2491, 16.25] as (A)'s ramp sups: inherited convention)

# ---------------------------------------------------------------- symbolic chi_k
s = sp.Symbol('s', positive=True)
dlf = sp.Function('dl')(s); uf = sp.Function('u')(s); ef = sp.Function('e')(s); Ef = sp.Function('E')(s)
kf = dlf.diff(s); bf = uf.diff(s)
ETA = 2 * sp.sin(dlf / 2) ** 2 + bf * sp.sin(dlf) + uf * kf * sp.cos(dlf)
EPS = (ef / s + ETA) / (1 - ETA)
NUM = (1 - kf * uf + bf ** 2) * sp.sin(dlf) + kf * uf * bf * sp.cos(dlf); DEN = (1 - kf * uf) * sp.cos(dlf) - bf * sp.sin(dlf); PHI = NUM / DEN
DS = sp.symbols('d0:6'); US = sp.symbols('u0 b0 b1 b2 b3 b4'); E0 = sp.Symbol('e0')
def to_syms(ex, K):
    for i in range(K, -1, -1): ex = ex.subs(Ef.diff(s, i) if i else Ef, EPS.diff(s, i) if i else EPS)
    for n in range(6, 0, -1): ex = ex.subs(ef.diff(s, n), -ETA.diff(s, n - 1))     # e' = Y' - 1 = -eta
    for n in range(5, -1, -1): ex = ex.subs(uf.diff(s, n) if n else uf, US[n])
    for n in range(5, -1, -1): ex = ex.subs(dlf.diff(s, n) if n else dlf, DS[n])
    return ex.subs(ef, E0)
Top = lambda F: (1 + Ef) * s * F.diff(s); Dop = lambda F: s * F.diff(s)
CHI = {}; Fb = bf; Fp = PHI; Db = bf
for kk in range(4):
    if kk: Fb, Fp, Db = Top(Fb), Top(Fp), Dop(Db)
    CHI[kk] = to_syms(sp.expand(Fb - Db), kk) + to_syms(Fp, kk)
b2, b3, b4 = US[3], US[4], US[5]
ZERO = {**{DS[i]: 0 for i in range(6)}, E0: 0}
for kk in range(4):
    c = CHI[kk]; assert b4 not in c.free_symbols
    assert sp.simplify(sp.diff(c, b2, 2)) == 0 and sp.simplify(sp.diff(c, b2, b3)) == 0 and sp.simplify(sp.diff(c, b3, 2)) == 0, kk   # affine in b2, b3
    assert sp.simplify(c.subs(ZERO)) == 0, kk                                                                                       # vanishes for a ray
SYM_A = {kk: CHI[kk].subs({b2: 0, b3: 0}) for kk in range(4)}; SYM_B2 = {kk: sp.diff(CHI[kk], b2) for kk in range(4)}; SYM_B3 = {kk: sp.diff(CHI[kk], b3) for kk in range(4)}
B22 = SYM_B2[2]; assert US[2] not in B22.free_symbols            # beta22 independent of b1 = u'', so D beta22 contains no u'''
ETA_S = ETA.subs({bf: US[1], uf: US[0], kf: DS[1]}).subs(dlf, DS[0])
BETA22 = B22 / s ** 2
def dsig(ex):                                                     # total d/dsigma of an expression in (s, e0, d_i, u0, b0, b1)
    return (sp.diff(ex, s) + sp.diff(ex, E0) * (-ETA_S) + sum(sp.diff(ex, DS[i]) * DS[i + 1] for i in range(5)) + sp.diff(ex, US[0]) * US[1]
            + sp.diff(ex, US[1]) * US[2] + sp.diff(ex, US[2]) * US[3])
DBETA22 = s * dsig(BETA22); assert b2 not in DBETA22.free_symbols
EPS_S = (E0 / s + ETA_S) / (1 - ETA_S)
IP = sp.Function('ipow')
def ipow(a, n):
    n = int(n)
    if n < 0: return 1 / ipow(a, -n)
    r = arb(1)
    for _ in range(n): r = r * a
    return r
fixpow = lambda ex: ex.replace(lambda q: q.is_Pow and q.exp.is_Integer, lambda q: IP(q.base, q.exp))
MODL = [{'sin': lambda a: a.sin(), 'cos': lambda a: a.cos(), 'ipow': ipow}]
ARGS = (s, E0, *DS[:5], US[0], US[1], US[2])
LAM = lambda ex: sp.lambdify(ARGS, fixpow(sp.sympify(ex)), modules=MODL)
FA = {kk: LAM(SYM_A[kk]) for kk in range(4)}; FB2 = {kk: LAM(SYM_B2[kk]) for kk in range(4)}; FB3 = {kk: LAM(SYM_B3[kk]) for kk in range(4)}
FDB22 = LAM(DBETA22); FETA = LAM(ETA_S); FEPS = LAM(EPS_S)
ev = lambda f, a: arb(UPf(arb(0) + f(*a)))

def boxes(s0, s1):
    """profile boxes for sigma in [s0, s1] (arb), at |alpha| <= alpha_1"""
    S = s0.union(s1)
    kb, mb, pb, qs = HT.profile_bounds(s0); mu = LT.mu(s0)
    dlb = kb / mu
    X = S - STIP; Jp = C4.PSI.jets_any(X, float(X.mid()), 5); Jh = C4.H2.jets_any(X, float(X.mid()), 5)
    pj = lambda j: AL1 * abs(Jp[j]) + AL1 * AL1 * abs(Jh[j])
    Ub = WSUP + pj(0); B0b = W1SUP + pj(1); B1b = W2SUP + pj(2); Pdd = pj(3); Pddd = pj(4)
    sn = (dlb / 2).sin()
    eb = qs + 2 * s1 * sn * sn + (pb + Ub) * dlb
    bx = lambda v: arb(0, UPf(v))
    a = (S, bx(eb), bx(dlb), bx(kb), bx(mb), bx(pb / 4), bx(kb * (s1 + qs) / 4), bx(Ub), bx(B0b), bx(B1b))
    return dict(S=S, a=a, e=arb(UPf(eb)), dl=arb(UPf(dlb)), kb=kb, mb=mb, pb=pb, qs=qs, Ub=Ub, B0=B0b, B1=B1b, Pdd=Pdd, Pddd=Pddd)

def cell(s0, s1, mode):
    bx = boxes(s0, s1); a = bx['a']; S = bx['S']
    A = [ev(FA[k], a) for k in range(4)]; B2 = [ev(FB2[k], a) for k in range(4)]; B3 = [ev(FB3[k], a) for k in range(4)]
    eta = ev(FETA, a); eps = ev(FEPS, a)
    e = bx['e']
    dU = e * C13 * AL1 * P * (s1 + e) ** (P - 1)                                   # |U(s) - U(s + e)|, U' = p U/y increasing
    add = [arb(UPf(A[k] + B2[k] * bx['Pdd'] + B3[k] * bx['Pddd'] + (B2[k] / s0 ** 2 + 2 * B3[k] / s0 ** 3) * s1 * W2SUP + P ** k * dU)) for k in range(4)]
    beta = {'22': arb(UPf(B2[2] / s0 ** 2)), '32': arb(UPf(B2[3] / s0 ** 2 + 3 * B3[3] / s0 ** 3)), '33': arb(UPf(B3[3] / s0 ** 3)), 'D22': ev(FDB22, a)}
    # measure factors
    Yp = 1 + arb(0, eta.upper()); E = arb(0, e.upper()); r = 1 / (1 + E / S)       # sigma / Y
    if mode == 'ramp':
        x = S + X0; q = (2 * x * E + E * E) / (1 + x * x); wr = (1 / (1 + q)) ** 2   # W(Y + X0)/W(sigma + X0)
        meas = max((Yp * r ** (2 * m) * wr for m in range(4)), key=lambda v: float(v.upper()))
    else:
        meas = max((Yp * r ** m for m in range(7)), key=lambda v: float(v.upper()))
    W = (1 + (S + X0) ** 2) ** -2
    return dict(s0=s0, s1=s1, add=add, beta=beta, eta=eta, eps=eps, e=e, meas=arb(UPf(meas)), W=arb(UPf(W)), bx=bx, A=A)

def dsliver(s0, s1):
    """d_B0 = alpha (psi3^' - 1.3056 y^p) + alpha^2 h2^' on [s0, s1] at |alpha| <= alpha_1: sups of |d_x^m d| (m <= 3) and |D^k d| (k <= 3), naive"""
    S = s0.union(s1); X = S - STIP
    Jp = C4.PSI.jets_any(X, float(X.mid()), 5); Jh = C4.H2.jets_any(X, float(X.mid()), 5)
    pw = lambda m: math.prod([P - i for i in range(m)], start=arb(1)) * S ** (P - m)  # d^m y^p
    g = [Jp[m + 1] - C13 * pw(m) for m in range(4)]; h = [Jh[m + 1] for m in range(4)]
    dx = [arb(UPf(AL1 * abs(g[m]) + AL1 * AL1 * abs(h[m]))) for m in range(4)]
    Dk = lambda J, k: sum((arb(c) * S ** m * J[m + 1] for m, c in STIR[k].items()), arb(0))
    gD = [Dk(Jp, k) - C13 * P ** k * S ** P for k in range(4)]; hD = [Dk(Jh, k) for k in range(4)]
    dD = [arb(UPf(AL1 * abs(gD[k]) + AL1 * AL1 * abs(hD[k]))) for k in range(4)]
    W = (1 + (S + X0) ** 2) ** -2
    return dict(dx=dx, dD=dD, W=arb(UPf(W)), w=s1 - s0)

if __name__ == '__main__':
    t0 = time.time()
    # ---- ramp: [14 - WS, 14], [14, 28] in steps 1/64, [28, 28 + WS]
    edges = [arb(14) - WS] + [arb(14) + arb(i) / 64 for i in range(14 * 64 + 1)] + [arb(28) + WS]
    ramp = [cell(edges[i], edges[i + 1], 'ramp') for i in range(len(edges) - 1)]
    assert ramp[0]['s0'] == arb(14) - WS and ramp[-1]['s1'] == arb(28) + WS and all(ramp[i]['s1'] == ramp[i + 1]['s0'] for i in range(len(ramp) - 1))
    e_ramp = max((c['e'] for c in ramp), key=lambda v: float(v.upper())); assert (e_ramp - WS).upper() < 0, 'shift exceeds the end sliver'
    core_edges = [arb(28) + arb(i) / 4 for i in range(32 * 4 + 1)]
    core = [cell(core_edges[i], core_edges[i + 1], 'core') for i in range(len(core_edges) - 1)]
    assert core[0]['s0'] == 28 and core[-1]['s1'] == 60
    mx = lambda L: arb(UPf(max(L, key=lambda v: float(v.upper()))))
    X = [arb(UPf(sum((c['W'] * c['add'][k] ** 2 * (c['s1'] - c['s0']) for c in ramp), arb(0)))) for k in range(4)]    # int W (chi^add_k)^2 dsigma
    Ec = [arb(UPf(sum((28 / c['s0'] ** 2 * c['add'][k] ** 2 * (c['s1'] - c['s0']) for c in core), arb(0)))) for k in range(4)]   # int (28/s) chi^2 ds
    R = dict(sup_add=[UPf(mx([c['add'][k] for c in ramp])) for k in range(4)], X_W=[UPf(v) for v in X],
             beta22=UPf(mx([c['beta']['22'] for c in ramp])), beta32=UPf(mx([c['beta']['32'] for c in ramp])), beta33=UPf(mx([c['beta']['33'] for c in ramp])),
             Dbeta22=UPf(mx([c['beta']['D22'] for c in ramp])), eps=UPf(mx([c['eps'] for c in ramp])), eta=UPf(mx([c['eta'] for c in ramp])),
             e=UPf(e_ramp), rho_meas=UPf(mx([c['meas'] for c in ramp]) - 1), domain=['14 - 2^-19', '28 + 2^-19'], cells=len(ramp))
    e_core = mx([c['e'] for c in core])
    C = dict(sup_add=[UPf(mx([c['add'][k] for c in core])) for k in range(4)], E=[UPf(v) for v in Ec],
             beta22=UPf(mx([c['beta']['22'] for c in core])), beta3=UPf(mx([c['beta']['32'] + c['beta']['33'] for c in core])), e=UPf(e_core),
             rho_meas=UPf(mx([c['meas'] for c in core]) - 1), domain=['28', '60'], cells=len(core))
    # ---- d_B0 on the ramp end slivers
    sl = [dsliver(arb(14) - WS, arb(14)), dsliver(arb(28), arb(28) + WS)]
    DS_ = dict(W=[UPf(sum((q['W'] * q['dx'][m] ** 2 * q['w'] for q in sl), arb(0))) for m in range(4)],
               sup=[UPf(max((q['dD'][k] for q in sl), key=lambda v: float(v.upper()))) for k in range(4)], intervals=['[14 - 2^-19, 14]', '[28, 28 + 2^-19]'])
    # ---- comparison with the reviewer's float probe (scratch2/chartprobe.py): delta-jet bounds and add_k at sigma = 14, 16, 20, 28
    probe = {}
    for s_ in (14, 16, 20, 28):
        c = next(c for c in ramp + core if c['s0'] == arb(s_))
        b = c['bx']; S0 = arb(s_)
        probe[str(s_)] = dict(delta=UPf(b['dl']), sigma_k=UPf(S0 * b['kb']), sigma2_kp=UPf(S0 ** 2 * b['mb']), sigma3_kpp=UPf(S0 ** 3 * b['pb'] / 4),
                              add=[UPf(v) for v in c['add']], A=[UPf(v) for v in c['A']], cell=[str(c['s0']), str(c['s1'])])
    out = dict(item='N1: chart terms chi_k = z^A_k o Y - D^k (w\' + d_B0), k <= 3 (lemma chart_lemma.tex), on (A)\'s ramp and core',
               label='certified (arb, outward); hypothesis |alpha| <= alpha_1 = 1.5e-4 (v3.3 bootstrap), verifier %s sups (alpha_1 row)' % CFG.VER,
               alpha_1='1.5e-4', verifier=dict(file=VER, w_sup=float(WSUP), w1_sup=float(W1SUP), w2_sup=float(W2SUP)),
               ramp=R, core=C, dB0_slivers=DS_, probe=probe,
               definitions=dict(X_W='int_{14-2^-19}^{28+2^-19} W(sigma + X0) add_k^2 dsigma (rate per unit tau)', E_core='int_28^60 (28/sigma) add_k^2 dsigma/sigma',
                                add='|chi_k^add| (D-form level k, includes p^k|U(sigma) - U(Y)| and the D w\' parts via sigma sI2)',
                                beta='|chi_k - chi_k^add| <= beta_k2 |D^2 w\'| + beta_k3 |D^3 w\'|; Dbeta22 = sup |D (B_22/sigma^2)|; eps: T = (1 + eps) D',
                                rho_meas='weight transfer y -> sigma, see docstring'),
               args=dict(argv=sys.argv[1:], env={k: os.environ.get(k) for k in ('S1TT', 'S1PROF', 'FORM_MULT', 'PYTHONHASHSEED')},
                         inputs={p: sha(p) for p in [os.path.abspath(__file__), VER] + MODS}, v36_env=CFG.env_record()))
    dst = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else f'{CFG.N}/n1/out/n1_chart_{CFG.VER}.json'   # v36: never writes into step4
    tmp = dst + '.tmp'; json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, dst)
    print('ramp sup add_k', ['%.4e' % v for v in R['sup_add']], 'X_W', ['%.4e' % v for v in R['X_W']])
    print('ramp beta22 %.3e beta32 %.3e beta33 %.3e Dbeta22 %.3e eps %.3e eta %.3e e %.3e rho_meas %.3e' % tuple(R[k] for k in ('beta22', 'beta32', 'beta33', 'Dbeta22', 'eps', 'eta', 'e', 'rho_meas')))
    print('core sup add', ['%.3e' % v for v in C['sup_add']], 'E', ['%.3e' % v for v in C['E']], 'beta22 %.3e beta3 %.3e e %.3e rho %.3e' % (C['beta22'], C['beta3'], C['e'], C['rho_meas']))
    print('d slivers W', ['%.3e' % v for v in DS_['W']], 'sup', ['%.3e' % v for v in DS_['sup']])
    for k_, v in probe.items(): print('probe', k_, {q: (('%.3e' % w) if isinstance(w, float) else w) for q, w in v.items() if q != 'cell'})
    print(f'({time.time() - t0:.1f}s)')
