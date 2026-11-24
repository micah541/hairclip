"""fix11 review, part D (1): theta-certificate post-processing, Gaarding identity + constants, Q2 Gram form.
Read-only on the proof tree (cwd = this directory; nothing is written elsewhere).  [cert] arb outward, [num] float/mpmath check.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B d_checks.py"""
import sys, os, json, hashlib, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
from flint import arb, ctx
ctx.prec = 200
out = {}; checks = []
def chk(name, ok, **kw): checks.append(dict(name=name, ok=bool(ok), **kw)); print(('PASS ' if ok else 'FAIL ') + name, kw, flush=True)
up = lambda z: float(arb(z).upper()); lo = lambda z: float(arb(z).lower())
S3 = '/home/claude/s1spec/step3'; V = '/home/claude/s1spec/step3nl_v4b'
# ============================================================ (1) theta-certificate
C = json.load(open(S3 + '/certificate_step3_theta_T1.json'))
G = [[arb(C['G3_mid'][i][j], C['G3_rad'][i][j] * (1 + 1e-12)) for j in range(3)] for i in range(3)]
g11, g12, g22, g13, g23, g33 = G[0][0], G[0][1], G[1][1], G[0][2], G[1][2], G[2][2]
det2 = g11 * g22 - g12 * g12
Sch = g33 - (g22 * g13 * g13 - 2 * g12 * g13 * g23 + g11 * g23 * g23) / det2
out['theta'] = dict(S_from_stored_mid_rad=[lo(Sch), up(Sch)], S_stored=C['S'], g11=up(g11), det2_lo=lo(det2), form_mult=C['form_mult'], omega=C['omega'])
chk('theta: G3 negative definite from the stored mid/rad (plain interval, no shared symbols)', up(g11) < 0 and lo(det2) > 0 and up(Sch) < 0, S=[lo(Sch), up(Sch)])
# re-run the certificate's own post-processing on the T1 CAPD outputs (read only; absolute prefix)
try:
    sys.path.insert(0, S3)
    os.environ['FORM_MULT'] = '9/8'; os.environ['OMEGA_Q'] = '1/32'
    import gram_post8 as GP
    r = GP.certify(S3 + '/outT1', verbose=False)
    Sr = [float(r['S'].lower()), float(r['S'].upper())]
    out['theta']['S_recomputed_gram_post8'] = Sr
    out['theta']['G3_recomputed_mid'] = [[float(v.mid()) for v in row] for row in r['G3']]
    same = all(abs(float(r['G3'][i][j].mid()) - C['G3_mid'][i][j]) <= 1e-12 * (1 + abs(C['G3_mid'][i][j])) for i in range(3) for j in range(3))
    chk('theta: gram_post8.certify(outT1) reproduces the stored certificate', r['ok'] and same and abs(Sr[0] - C['S'][0]) < 1e-12 and abs(Sr[1] - C['S'][1]) < 1e-12, S=Sr)
except Exception as e:
    out['theta']['recompute_error'] = repr(e); chk('theta: gram_post8.certify(outT1) rerun', False, err=repr(e))
n_ok = 0
import glob
fl = sorted(glob.glob(S3 + '/outT1*.txt')) + [S3 + '/out2P_T1.txt']
for f in fl: n_ok += ('RESULT OK' in open(f).read())
P = json.load(open(S3 + '/positivity_cert_T1.json'))
chk('theta: 34 CAPD outputs + Morse run say RESULT OK; positivity cert ok with the same (form_mult, omega)', n_ok == 35 and P['ok'] and P['form_mult'] == '9/8' and P['omega'] == '1/32' and C['form_mult'] == '9/8' and C['omega'] == '1/32', n_ok=n_ok)
out['theta']['pins'] = {os.path.basename(f): hashlib.sha256(open(f, 'rb').read()).hexdigest() for f in fl + [S3 + '/certificate_step3_theta_T1.json', S3 + '/positivity_cert_T1.json']}
# ============================================================ (2) Gaarding: identity [num] and constants [cert]
import mpmath as mp
mp.mp.dps = 30
x = mp.mpf
phi = lambda t: (1 + t * t) ** (mp.mpf(2) / 3) / 4
rho = lambda t: mp.exp(-phi(t))
d = lambda f, n=1: (lambda t: mp.diff(f, t, n))
# generic smooth data: v even, c2 even, c1 odd, c0 even
v = lambda t: (1 + mp.mpf('0.7') * t * t - mp.mpf('0.2') * t ** 4) * mp.exp(-t * t / 5) + mp.mpf('0.3') * mp.cos(t) * mp.exp(-t * t / 9)
c2 = lambda t: -mp.mpf('1.3') / (1 + t * t) + mp.mpf('0.2') * mp.cos(t)
c1 = lambda t: -mp.mpf('0.25') * t - mp.mpf('2') * t / (1 + t ** 4)
c0 = lambda t: mp.mpf('0.25') - 3 / (1 + t * t) ** 2
rr2 = lambda f: (lambda t: mp.diff(f, t, 2) - 2 * mp.diff(phi, t) * mp.diff(f, t) + (mp.diff(phi, t) ** 2 - mp.diff(phi, t, 2)) * f(t))
Af = lambda t: 2 * (mp.diff(phi, t) ** 2 - mp.diff(phi, t, 2)) - c2(t)
g2 = lambda t: mp.diff(phi, t) ** 2 - mp.diff(phi, t, 2)
Vf = lambda t: c0(t) - (mp.diff(c1, t) - mp.diff(phi, t) * c1(t)) / 2 + rr2(c2)(t) / 2 - rr2(g2)(t) / 2
Mv = lambda t: -mp.diff(v, t, 4) + c2(t) * mp.diff(v, t, 2) + c1(t) * mp.diff(v, t) + c0(t) * v(t)
pts = [0, 1, 2, 4, 7, 11, 16, 22]
lhs = mp.quad(lambda t: rho(t) * Mv(t) * v(t), pts)
rhs = mp.quad(lambda t: rho(t) * (-mp.diff(v, t, 2) ** 2 + Af(t) * mp.diff(v, t) ** 2 + Vf(t) * v(t) ** 2), pts)
muf = lambda t: -Af(t) ** 2 / 2 - rr2(Af)(t) / 2 - Vf(t)
ub = mp.quad(lambda t: rho(t) * (-mp.diff(v, t, 2) ** 2 / 2 - muf(t) * v(t) ** 2), pts)
out['garding_identity'] = dict(lhs=float(lhs), rhs=float(rhs), abs_diff=float(abs(lhs - rhs)), upper_form=float(ub))
chk("Gaarding: <Mv,v>_rho = int rho(-v''^2 + A v'^2 + V v^2) on generic even data [num]", abs(lhs - rhs) < 1e-12 * (1 + abs(lhs)), diff=float(abs(lhs - rhs)))
chk("Gaarding: <Mv,v>_rho <= -1/2||v''||^2 - int rho mu v^2 [num]", lhs <= ub + 1e-15, lhs=float(lhs), ub=float(ub))
# level-2 reduction (L w)'' = Mt w'' + e1 w' + e0 w and the constants, with my own sympy expressions, on the same hull boxes (finer grouping)
import sympy as sp
X, k, m, p, q = sp.symbols('X k m p q')
ph = (1 + X ** 2) ** sp.Rational(2, 3) / 4
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4}
D = lambda e: sp.diff(e, X) + sum(sp.diff(e, s) * r for s, r in rules.items())
a2, a1, a0 = -k ** 2, -(5 * k * m + q / 4), sp.Rational(1, 4) - 3 * m ** 2 - k * p
wf = sp.Function('w')(X)
# check (L w)'' expansion symbolically with k,m,p,q as functions through D on coefficients only
w0, w1, w2, w3, w4, w5, w6 = sp.symbols('w0:7'); wj = [w0, w1, w2, w3, w4, w5, w6]
def Dw(e): return D(e) + sum(sp.diff(e, wj[i]) * wj[i + 1] for i in range(6))
Lw = -w4 + a2 * w2 + a1 * w1 + a0 * w0
Lw2 = Dw(Dw(Lw))
c1t = 2 * D(a2) + a1; c0t = D(D(a2)) + 2 * D(a1) + a0; e1 = D(D(a1)) + 2 * D(a0); e0 = D(D(a0))
res = sp.expand(Lw2 - (-w6 + a2 * w4 + c1t * w3 + c0t * w2 + e1 * w1 + e0 * w0))
chk("Gaarding: (Lw)'' = Mt w'' + e1 w' + e0 w, Mt = -d^4 + a2 d^2 + (2a2'+a1) d + (a2''+2a1'+a0) [sym]", res == 0 and sp.expand(c1t + 9 * k * m + q / 4) == 0)
chk('Gaarding: boundary data: a1, c1t odd-type (vanish at the tip m = q = 0)', sp.expand(a1.subs({m: 0, q: 0})) == 0 and sp.expand(c1t.subs({m: 0, q: 0})) == 0)
P1 = sp.diff(ph, X)
R2 = lambda f: D(D(f)) - 2 * P1 * D(f) + (P1 ** 2 - D(P1)) * f
A2 = 2 * (P1 ** 2 - D(P1)) - a2
V2 = c0t - (D(c1t) - P1 * c1t) / 2 + R2(a2) / 2 - R2(P1 ** 2 - D(P1)) / 2
mu2 = -A2 ** 2 / 2 - R2(A2) / 2 - V2
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n)); return f"({prod})" if n > 0 else f"(1/({prod}))"
        return f"({bb}**(arb({ex.p})/arb({ex.q})))"
mk = lambda e: eval("lambda X, k, m, p, q: " + _Pr().doprint(e), {"arb": arb})
ctx.prec = 100
f2, f1e, f0e = mk(mu2), mk(sp.expand(e1)), mk(sp.expand(e0))
TT = arb('10.249073'); mx = dict(mu2=-1e300, e1=0.0, e0=0.0); nb = 0
# prof_tight.txt rows: t_lo t_hi ... hull columns 13.. = (p, q, k, m) lo/hi  (same reading as garding_cert.py; x = c (T0 - t) with the reader's own TT from h2res)
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, S3)
import h2res as HR
lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
GROUP = 10
for g0 in range(0, len(lines), GROUP):
    grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
    ta, tb = min(v_[0] for v_ in grp), max(v_[1] for v_ in grp)
    pp, qq, kk, mm = [HR.I(min(v_[13 + 2 * j] for v_ in grp), max(v_[14 + 2 * j] for v_ in grp)) for j in range(4)]
    Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
    if Xb.upper() < 0: continue
    if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
    mx['mu2'] = max(mx['mu2'], up(-f2(Xb, kk, mm, pp, qq))); mx['e1'] = max(mx['e1'], up(abs(f1e(Xb, kk, mm, pp, qq)))); mx['e0'] = max(mx['e0'], up(abs(f0e(Xb, kk, mm, pp, qq)))); nb += 1
for ci in range(8, 32):
    for i in range(50):
        sa = arb(ci) + arb(i) / 50; sb = arb(ci) + arb(i + 1) / 50
        kb, mb, pb, qs = HR.HS.profile_bounds(sa); Xb = (sa - HR.ST).union(sb - HR.ST)
        a_ = (Xb, arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper()), sa.union(sb) + arb(0, qs.upper()))
        mx['mu2'] = max(mx['mu2'], up(-f2(*a_))); mx['e1'] = max(mx['e1'], up(abs(f1e(*a_)))); mx['e0'] = max(mx['e0'], up(abs(f0e(*a_)))); nb += 1
ctx.prec = 200
E2c = {k_: v_['value'] for k_, v_ in json.load(open(V + '/v4bcert/out/e2v4b_constants.json'))['constants'].items()}
out['garding_consts'] = dict(mine=mx, boxes=nb, coded=dict(r143=E2c['R143'], r163=16.3, r107=10.7))
chk('Gaarding constants: coded r143, r163, r107 >= my sups (own expressions, finer hull grouping) [cert]', E2c['R143'] >= mx['mu2'] and 16.3 >= mx['e1'] and 10.7 >= mx['e0'], mine=mx)
# D_rho constants t22 = 2.2, t12 = 1.2 on [0, 34.25]: sup 2|phi'|, sup |phi'^2 - phi''|
s1 = 0.0; s2 = 0.0; N = 3425
for i in range(N):
    Xb = arb(i) / 100; Xb = Xb.union(arb(i + 1) / 100); E = (1 + Xb * Xb) ** (arb(-1) / 3); p1 = Xb * E / 3; p2 = E / 3 - 2 * Xb * Xb * E ** 4 / 9
    s1 = max(s1, up(abs(2 * p1))); s2 = max(s2, up(abs(p1 * p1 - p2)))
out['Drho'] = dict(sup_2phi1=s1, sup_phi1sq_minus_phi2=s2)
chk('D_rho: 2.2 >= sup 2|phi\'| and 1.2 >= sup|phi\'^2 - phi\'\'| on [0,34.25] [cert]', s1 <= 2.2 and s2 <= 1.2, s1=s1, s2=s2)
# ============================================================ (3) Q2 Gram form: own eigenvalue bound
Gj = json.load(open('/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json'))
assert Gj['names'] == ['Lf1', 'Lf2', 'gA', 'gB', 'Df1', 'Df2', 'f1', 'f2']
B8 = [[arb(repr(Gj['G'][i][j][0])).union(arb(repr(Gj['G'][i][j][1]))) for j in range(8)] for i in range(8)]
sym_ok = all(B8[i][j].overlaps(B8[j][i]) for i in range(8) for j in range(8))
M = [[B8[6][6], B8[6][7]], [B8[7][6], B8[7][7]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; assert det.lower() > 0
Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
idx = [0, 1, 4, 5]; F = [6, 7]
Gs = [[B8[a][b] - sum((B8[a][F[i]] * Mi[i][j] * B8[F[j]][b] for i in range(2) for j in range(2)), arb(0)) for b in idx] for a in idx]
# Q2(w) = -a <w, Df1> - b <w, Df2>, (a, b) = Minv (l1, l2):  Q = l^T K l with K (4x4) = -(1/2)(N + N^T), N[r][2+c]... build directly
K4 = [[arb(0)] * 4 for _ in range(4)]
for i in range(2):          # a_i = sum_j Mi[i][j] l_j ; Q = - sum_i a_i l_{2+i}
    for j in range(2):
        K4[j][2 + i] += -Mi[i][j] / 2; K4[2 + i][j] += -Mi[i][j] / 2
# sup over the ellipsoid {l = G^{1/2} b}: lambda_max(K4 Gs).  Certified bound: c is an upper bound iff c G^-1 - K4 >= 0 on the range, i.e. iff  c Gs - Gs K4 Gs is PSD (congruence with Gs).
def psd_cert(c):
    Mx = [[c * Gs[i][j] - sum((Gs[i][a] * K4[a][b] * Gs[b][j] for a in range(4) for b in range(4)), arb(0)) for j in range(4)] for i in range(4)]
    # symmetric interval matrix -> mid - (row radius sums) I must be PD (arb Cholesky on exact mids)
    mid = [[(arb(Mx[i][j].mid()) + arb(Mx[j][i].mid())) / 2 for j in range(4)] for i in range(4)]
    rad = [sum((arb(Mx[i][j].rad()) + abs(arb(Mx[i][j].mid()) - arb(Mx[j][i].mid())) / 2 for j in range(4)), arb(0)) for i in range(4)]
    A_ = [[mid[i][j] - (arb(rad[i].upper()) if i == j else 0) for j in range(4)] for i in range(4)]
    R = [[arb(0)] * 4 for _ in range(4)]
    for i in range(4):
        s = A_[i][i] - sum((R[k_][i] * R[k_][i] for k_ in range(i)), arb(0))
        if not (s.lower() > 0): return False
        R[i][i] = s.sqrt()
        for j in range(i + 1, 4): R[i][j] = (A_[i][j] - sum((R[k_][i] * R[k_][j] for k_ in range(i)), arb(0))) / R[i][i]
    return True
import numpy as np
Gf = np.array([[float(Gs[i][j].mid()) for j in range(4)] for i in range(4)]); Kf = np.array([[float(K4[i][j].mid()) for j in range(4)] for i in range(4)])
ev = np.sort(np.linalg.eigvals(Kf @ Gf).real)
cands = [43.33, 43.3271, 43.3275, 43.33, 43.4, 43.5, 44.0]
cert = {repr(c): psd_cert(arb(repr(c))) for c in cands}
out['Q2'] = dict(float_eigs=ev.tolist(), certified_upper={k_: v_ for k_, v_ in cert.items()}, reported=43.327097178333574, coded=43.33, gram_symmetric_overlap=sym_ok,
                 Gs_rel_rad=max(float(Gs[i][j].rad()) / abs(float(Gs[i][j].mid())) for i in range(4) for j in range(4)))
chk('Q2: coded 43.33 >= float lambda_max of the Gram form (own construction)', 43.33 >= ev[-1], lam=float(ev[-1]))
out['Q2']['direct_congruence_test_c_G_minus_GKG'] = cert      # informational: inconclusive (Gs is nearly singular, the interval radii swamp the small pivots)
# Loewner-majorant route (own code): every true Gram G in [Gs] satisfies G <= G+ := mid + diag(row radius sums) (Gershgorin on G - mid), E(G) subset E(G+),
# sup_{l in E(G+)} l^T K l = lambda_max(R K R^T), R^T R = G+ (R upper).  K is an interval matrix (Minv): lambda_max <= c iff c I - mid(B) - diag(row rad sums) PD.
def chol(A_, n):
    R = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        s_ = A_[i][i] - sum((R[k_][i] * R[k_][i] for k_ in range(i)), arb(0))
        if not (s_.lower() > 0): return None
        R[i][i] = s_.sqrt()
        for j in range(i + 1, n): R[i][j] = (A_[i][j] - sum((R[k_][i] * R[k_][j] for k_ in range(i)), arb(0))) / R[i][i]
    return R
Gm = [[(arb(Gs[i][j].mid()) + arb(Gs[j][i].mid())) / 2 for j in range(4)] for i in range(4)]
Gr = [sum((arb((arb(Gs[i][j].rad()) + abs(arb(Gs[i][j].mid()) - arb(Gs[j][i].mid())) / 2).upper()) for j in range(4)), arb(0)) for i in range(4)]
Gp = [[Gm[i][j] + (Gr[i] if i == j else 0) for j in range(4)] for i in range(4)]
Rc = chol(Gp, 4); assert Rc is not None
Ks = [[K4[i][j].union(K4[j][i]) for j in range(4)] for i in range(4)]
Bm = [[sum((Rc[i][a] * Ks[a][b] * Rc[j][b] for a in range(4) for b in range(4)), arb(0)) for j in range(4)] for i in range(4)]
def lam_le(c):
    rad = [sum((arb((arb(Bm[i][j].rad()) + abs(arb(Bm[i][j].mid()) - arb(Bm[j][i].mid())) / 2).upper()) for j in range(4)), arb(0)) for i in range(4)]
    A_ = [[(c if i == j else 0) - (arb(Bm[i][j].mid()) + arb(Bm[j][i].mid())) / 2 - (rad[i] if i == j else 0) for j in range(4)] for i in range(4)]
    return chol(A_, 4) is not None
cert2 = {repr(c): lam_le(arb(repr(c))) for c in (43.3271, 43.3272, 43.328, 43.33)}
out['Q2']['certified_upper_majorant_route'] = cert2
chk('Q2: lambda_max <= 43.33 certified by my own Loewner-majorant + Cholesky code [cert]', cert2[repr(43.33)], cert=cert2)
out['checks'] = checks; out['n_fail'] = sum(1 for c in checks if not c['ok'])
out['args'] = dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                   inputs={p_: hashlib.sha256(open(p_, 'rb').read()).hexdigest() for p_ in (S3 + '/certificate_step3_theta_T1.json', S3 + '/gram_post8.py', S3 + '/gram_post7.py', V + '/garding_cert.py', V + '/garding_cert.json',
                                                                                           '/home/claude/s1_tight/prof/prof_tight.txt', '/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json', V + '/agent_mod/mod_forms.py', V + '/v4bcert/out/e2v4b_constants.json')},
                   script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
json.dump(out, open(os.path.join(HERE, 'd_checks.json'), 'w'), indent=1, sort_keys=True)
print('n_fail', out['n_fail'])
