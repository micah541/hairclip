"""T8 certified (arb): per-slab enclosures of the Lemma F constants (LEMMA_F.md sections 2-3 and Lemma B) for levels
j = 0..3, design opt4n (ramp [14, 22] for every level, delta_P = (1/4, 1/4, 1/4, 3/2)), on the slabs of the polished P11 base
that T6 certified (v2 base t1/polished_v2: slabs 0..180; v1: 5..187), from T6's certified per-element hulls (t6/jets) and
T6's certified scalar defects (t6/residuals.json).

  python cconst.py --slabs 0-180 --out constants/slabs --jobs 2        (restartable; one JSON per slab)
  python merge.py --in constants/slabs --out constants/constants.json  (slab range: T6's residuals.json)

Method (no sampling):
  * Every base quantity is the T6 hull [lo, hi] over (element sub-cell) x (whole slab) -> an arb ball.  On the fixed
    elements (x <= 40) the sub-cell x-interval is exact, so it is split further (ramp, tip norms); on the moving outer
    elements the x-hull already contains the motion of the knots over the slab.
  * Coefficients: beta = Q - c X0 - 6 k k' + A (true transport, T6 sec. 5), b2 = -k^2, c0 = -3 k'^2 - 4 k k'', all
    derivatives by Leibniz on hulls; A' = k r and higher from the hulls of r, r', r'', r'''.
  * The identity forms come from cforms.py (generic beta, b2, c0; Euler-checked).  Core (eta = 1): every term is divided
    by T analytically (log-derivative jets tau_m = T^(m)/T of the flattened Y_l), and the transport/time-derivative
    cancellation  d_tau ln T + c y tau_1 + (2-2j) c = (2j-2)(D - beta(l))(Phi - s Phi')/Y,  D := b1(l) + l_tau,
    is done exactly before any interval evaluation.  Ramp: the foot ratios P^2/T, T'^2/T are evaluated as bounded
    rational functions of t (S = t^4 sigma0, S^(l) = t^(4-l) sigma_l, sigma0 >= 1).
  * Sups are upper bounds of arb balls (rounded outward); integrals are sums of (cell length) x (cell sup).
Inputs tagged [T7] (K_e, mu_e, rho_c) come from t7/constants.json when it has them (certified closure block and K_e),
otherwise they are T7 prototype values."""
import os, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import json, math, time, argparse, subprocess
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
T6DIR = (os.path.normpath(os.path.abspath(os.environ['T6DIR'])) if os.environ.get('T6DIR')    # env T6DIR: another T6 data dir
         else os.path.normpath(os.path.join(HERE, '..', 't6')))                               # (e.g. p2/t6c, review B1)
T7DIR = os.path.normpath(os.path.join(HERE, '..', 't7'))
ENVKEYS = ('TOPLOSS', 'NLFAC', 'T6DIR', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')   # env read by this script
PVKEYS = ('t6dir', 'jets_sha256', 'residuals_sha256', 't7_inputs_sha256', 'env')   # a cached slab file is reused only if these match


def provenance(i, outdir, t7):
    """args block of a slab file (kit B): argv, env and digests of every input the slab depends on"""
    import hashlib
    dig = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    jf = os.path.join(T6DIR, 'jets', f'jets_{i:03d}.npz')
    t7d = hashlib.sha256(json.dumps({k: v for k, v in t7.items() if k != '_cert_status'}, sort_keys=True).encode()).hexdigest()
    return dict(argv=sys.argv[1:], slab=i, out=os.path.normpath(os.path.abspath(outdir)), t6dir=T6DIR, jets=jf, jets_sha256=dig(jf),
                residuals_sha256=dig(os.path.join(T6DIR, 'residuals.json')), t7const=os.path.join(T7DIR, 'constants.json'),
                t7_inputs_sha256=t7d, env={k: os.environ[k] for k in ENVKEYS if k in os.environ})

from flint import arb, ctx
ctx.prec = 64

JMAX = 3
X0 = arb('2.2490776515')
YA, YB = 14.0, 22.0                       # opt4n ramp, levels 0..3 (LEMMA_F 1.3; T6 RF*_opt4n uses the same weights)
DP = (0.25, 0.25, 0.25, 1.5)
DS = (0.5, 0.5, 0.5, 0.25)            # (S6): -2 int T w_{j+2} rho <= DS_j D_j + ||rho||_T^2 / DS_j  (theta_3 > 0 needs DS_3 < 1/2)
A_S3, NU = arb('0.01'), arb('0.001')
ZETA_FAC = 1.0 + 1e-12                    # zeta_inf >= 1 - 2e-14 on [x_a, 16.25]: V_m(W only) <= (1 + 1e-13) V_m(T5)
NSUB_RAMP, NSUB_TIP = 64, 4
XTIP = 34.25                              # supp chi_32 ends at y = 32 (x = 34.249)


def t7_inputs():
    """T7 constants for Lemma J / Corollary C.  From t7/constants.json (certified, some entries pending/provisional) where it
    provides them: mu1, mu2, a_f, b_f, c_I, Psi1, c1W, C_gamma0, C_gamma2 (worst case over all its slabs), C_x (max of the
    available Cx_full), nV_min (per slab), and, when present, the certified closure (mu_e, min rho_c) and K_e.  Otherwise
    mu_e, rho_c are the t7/t8req.json P=2 best row and K_e = 1.4e4, flagged PLACEHOLDER.  theta_0, eps, eps_I, c1, kappa_D:
    design values."""
    q = json.load(open(os.path.join(T7DIR, 't8req.json')))['P=2']
    b = q['best']; w = q['worst']
    d = dict(muE=b['muE'], rho=b['rho'], th0=0.025, eps=1e-5, epsI=3e-5, c1=0.5804, Ke=1.4e4, kD=0.5)
    src = dict(muE='PLACEHOLDER t8req', rho='PLACEHOLDER t8req (float closure)', th0='T7 theta_I', eps='T5 default',
               epsI='T5 default', c1='S3 interpolation', Ke='PLACEHOLDER (T5 2e3/sqrt(mu_e) scale)', kD='Lemma J (T7 closure reserves 1/2)')
    p = os.path.join(T7DIR, 'constants.json')
    nVmin = {}
    if os.path.exists(p):
        C = json.load(open(p)); S = C['slabs']
        vals = lambda f: [f(v) for v in S.values()]
        d.update(mu1=min(vals(lambda v: v['E_I']['mu1'])), mu2=min(vals(lambda v: v['E_I']['mu2'])),
                 af=max(vals(lambda v: v['ramp']['a_f'])), bf=max(vals(lambda v: v['ramp']['b_f'])),
                 cI=[max(vals(lambda v, k=k: v['ramp'][k])) for k in ('c0', 'c1', 'c2', 'c3')],
                 c1W=C['meta']['weight_consts']['c1W'], Psi1=C['meta']['weight_consts']['Psi1'],
                 Cg0=max(vals(lambda v: v['modulation']['Cgamma0'])), Cg2=max(vals(lambda v: v['modulation']['Cgamma2'])))
        cx = [v['C3']['Cx_full'] for v in S.values() if 'C3' in v]
        d['Cx'] = max(cx + [7.5]) if cx else 7.5
        nVmin = {int(k): v['C3']['nV_min'] for k, v in S.items() if 'C3' in v}
        npend = sum(1 for v in S.values() if 'C3' not in v)
        src.update(mu1='t7/constants.json', mu2='t7/constants.json', af='t7/constants.json', bf='t7/constants.json',
                   cI='t7/constants.json', c1W='t7/constants.json', Psi1='t7/constants.json',
                   Cg0='t7/constants.json (PROVISIONAL: C_x at tau_c until c2r3)', Cg2='t7/constants.json (PROVISIONAL)',
                   Cx=f'max(Cx_full of the {len(cx)} slabs with C3, 7.5) ({npend} slabs pending c2r3)')
        d['_cert_status'] = C['meta']['summary'].get('rank3_status', {})
        cl = C['meta'].get('closure')
        if cl and cl.get('all_certified'):
            d.update(muE=cl['mu_e'], rho=cl['rho_c_min'])
            src.update(muE='t7/constants.json closure (certified mu_e)', rho='t7/constants.json closure (min certified rho_c)')
        if 'K_e' in C['meta']:
            d['Ke'] = C['meta']['K_e']['value']
            src['Ke'] = 't7/constants.json K_e (certified: |e| <= K_e Ybar for x <= 14 + X0)'
    else:
        d.update(mu1=w['mu1'], mu2=w['mu2'], af=w['af'], bf=w['bf'], cI=list(w['cI']), c1W=w['c1W'], Psi1=w['Psi1'],
                 Cg0=w['Cg0'], Cg2=w['Cg2'], Cx=7.5)
        for k in ('mu1', 'mu2', 'af', 'bf', 'cI', 'c1W', 'Psi1', 'Cg0', 'Cg2', 'Cx'): src[k] = 'PLACEHOLDER T7 prototype'
    d['gam0'] = d['Cg0'] / d['Cx']; d['gam2'] = d['Cg2'] / d['Cx']
    d['nVmin'] = nVmin; d['_sources'] = src
    d['_source'] = 't7/constants.json + placeholders' if os.path.exists(p) else 'PLACEHOLDER: T7 prototype'
    return d


# ---------------------------------------------------------------------------------------------------------- arb helpers
def H(v): return arb(float(v[0])).union(arb(float(v[1])))
def iv(lo, hi): return arb(lo).union(arb(hi))
def U(a):
    f = float(a.upper())
    if not math.isfinite(f): raise ValueError(f'non-finite upper bound {a}')
    return math.nextafter(f, math.inf)
def L(a):
    f = float(a.lower())
    if not math.isfinite(f): raise ValueError(f'non-finite lower bound {a}')
    return math.nextafter(f, -math.inf)
def pos(a):                                # enclosure of max(a, 0)
    lo, hi = float(a.lower()), float(a.upper())
    if hi <= 0: return arb(0)
    return iv(max(lo, 0.0), math.nextafter(hi, math.inf))
def absb(a):                               # enclosure of |a| with a correct lower end
    lo, hi = float(a.lower()), float(a.upper())
    if lo >= 0: return a
    if hi <= 0: return -a
    return iv(0.0, math.nextafter(max(-lo, hi), math.inf))
def sqr(a): b = absb(a); return b * b
def sqrtp(a):                              # sqrt of a nonnegative quantity (lower end clipped at 0)
    lo, hi = max(0.0, float(a.lower())), float(a.upper())
    if not math.isfinite(hi): raise ValueError(f'non-finite {a}')
    return arb(lo).sqrt().lower().union(arb(math.nextafter(hi, math.inf)).sqrt().upper())
def jmul(a, b, n):
    return [sum((comb(m, l) * a[l] * b[m - l] for l in range(m + 1)), arb(0)) for m in range(n + 1)]
def U2(a, b): return max(a, b)


# ------------------------------------------------------------------------------------------------------- forms (cforms)
_FORMS = None
def forms():
    global _FORMS
    if _FORMS is None:
        F = json.load(open(os.path.join(HERE, 'forms.json')))
        _FORMS = {}
        for j, d in F.items():
            j = int(j)
            qn = [m for m in d['q'] if 'c' not in m[2]]
            qc = sorted((m[0], m[1], tuple(sorted(m[2].items()))) for m in d['q'] if 'c' in m[2])
            assert qc == sorted([(1, 1, (('T1', 1), ('c', 1), ('y', 1))), (1 - 2 * j, 1, (('T0', 1), ('c', 1)))]), qc
            _FORMS[j] = dict(q=_compile(d['q']), qn=_compile(qn), P=_compile(d['P']),
                             pairs={tuple(int(t) for t in k.split(',')): _compile(v) for k, v in d['pairs'].items()})
    return _FORMS

def _compile(mons):
    """polynomial -> python function of a dict of arb values (exact rational coefficients)"""
    terms = []
    for num, den, pw in mons:
        cf = 'arb(%d)' % num if den == 1 else '(arb(%d)/%d)' % (num, den)
        fac = '*'.join(('V["%s"]' % s) if p == 1 else ('V["%s"]**%d' % (s, p)) for s, p in sorted(pw.items()))
        terms.append(cf + ('*' + fac if fac else ''))
    src = 'def f(V):\n    return ' + (' + '.join(terms) if terms else 'arb(0)') + '\n'
    ns = dict(arb=arb); exec(src, ns); return ns['f']


# --------------------------------------------------------------------------------------------------------- weights
SIG_CHECKED = False
def _check_sig():
    global SIG_CHECKED
    if SIG_CHECKED: return
    import sympy as sp
    u = sp.Symbol('u'); S = 35 * u**4 - 84 * u**5 + 70 * u**6 - 20 * u**7
    sig = [35 - 84 * u + 70 * u**2 - 20 * u**3, 140 * (1 - u)**3, 420 * (1 - u)**2 * (1 - 2 * u),
           840 * (1 - u) * (1 - 5 * u + 5 * u**2), 840 * (1 - 2 * u) * (1 - 10 * u + 10 * u**2)]
    for l in range(5): assert sp.expand(sp.diff(S, u, l) - u**(4 - l) * sig[l]) == 0
    SIG_CHECKED = True

def sigmas(u):
    """sigma_l(u) with S^(l) = u^(4-l) sigma_l; sigma0 by monotonicity (decreasing on [0,1], sigma0 >= 1)"""
    lo, hi = float(u.lower()), float(u.upper())
    s0f = lambda z: 35 - 84 * z + 70 * z * z - 20 * z ** 3
    s0 = (arb(35) - 84 * arb(hi) + 70 * arb(hi) ** 2 - 20 * arb(hi) ** 3).union(arb(35) - 84 * arb(lo) + 70 * arb(lo) ** 2 - 20 * arb(lo) ** 3)
    v = 1 - u
    return [s0, 140 * v ** 3, 420 * v ** 2 * (1 - 2 * u), 840 * v * (1 - 5 * u + 5 * u ** 2), 840 * (1 - 2 * u) * (1 - 10 * u + 10 * u ** 2)]

def ramp_weight(y, p, LN):
    """level weight T = 28 eta^2 y^p on the ramp, eta = S(t), t = ln(y/ya)/ln(yb/ya).  Returns dict with
    T[m] = T^(m) (m <= 4), Tt[m] = T^(m)/eta (regular), Z = y^p, eta, eta1 = eta'."""
    t = (y / YA).log() / LN
    lo, hi = max(0.0, float(t.lower())), min(1.0, float(t.upper()))
    t = iv(lo, max(lo, hi))
    sg = sigmas(t)
    Sl = [t ** (4 - l) * sg[l] if l < 4 else sg[4] for l in range(5)]
    R = {}
    def Rl(l1, l2):
        key = (min(l1, l2), max(l1, l2))
        if key not in R: R[key] = t ** (4 - l1 - l2) * sg[l1] * sg[l2] / sg[0]
        return R[key]
    t1 = 1 / (y * LN); t2 = -1 / (y * y * LN); t3 = 2 / (y ** 3 * LN); t4 = -6 / (y ** 4 * LN)
    cf = [{0: arb(1)}, {1: t1}, {2: t1 * t1, 1: t2}, {3: t1 ** 3, 2: 3 * t1 * t2, 1: t3},
          {4: t1 ** 4, 3: 6 * t1 * t1 * t2, 2: 3 * t2 * t2 + 4 * t1 * t3, 1: t4}]
    eta_d = [sum((c * Sl[l] for l, c in cf[k].items()), arb(0)) for k in range(5)]
    def E(a, b):               # eta^(a) eta^(b) / eta
        if a == 0: return eta_d[b]
        if b == 0: return eta_d[a]
        return sum((ca * cb * Rl(la, lb) for la, ca in cf[a].items() for lb, cb in cf[b].items()), arb(0))
    e2 = [sum((comb(m, k) * E(k, m - k) for k in range(m + 1)), arb(0)) for m in range(5)]     # (eta^2)^(m)/eta
    fall = lambda k: math.prod(p - i for i in range(k))
    Zd = [fall(k) * y ** (p - k) for k in range(5)]
    Tt = [28 * sum((comb(m, a) * e2[a] * Zd[m - a] for a in range(m + 1)), arb(0)) for m in range(5)]
    eta = eta_d[0]
    return dict(T=[eta * v for v in Tt], Tt=Tt, Z=Zd[0], eta=eta, eta1=eta_d[1])

PSIC = (arb(21) / 4, arb(-49) / 8, arb(15) / 8)
def _psi(u):
    a, b, c = PSIC
    return [u ** 5 * (a + b * u + c * u * u), 5 * a * u ** 4 + 6 * b * u ** 5 + 7 * c * u ** 6,
            20 * a * u ** 3 + 30 * b * u ** 4 + 42 * c * u ** 5, 60 * a * u * u + 120 * b * u ** 3 + 210 * c * u ** 4,
            a * u ** 6 / 6 + b * u ** 7 / 7 + c * u ** 8 / 8]

def core_weight(y, yl, p):
    """T = 28 Y^p, Y = yl Phi(y/yl): returns tau[1..4] = T^(m)/T, T0 = T, Y, and xi := (Phi - s Phi')/Y"""
    slo = max(0.0, float((y / yl).lower())); shi = min(1.0, float((y / yl).upper()))
    branches = []
    if slo < 0.5:                                   # Phi = s: Y = y
        branches.append(dict(Y=y, L=[1 / y, arb(0), arb(0), arb(0)], xi=arb(0)))
    if shi > 0.5:
        s = iv(max(slo, 0.5), shi); u = 2 * s - 1
        lo_u = max(0.0, float(u.lower())); u = iv(lo_u, max(lo_u, float(u.upper())))
        ps = _psi(u)
        Ph = arb('0.5') + (u - ps[4]) / 2; d1 = 1 - ps[0]; d2 = -2 * ps[1]; d3 = -4 * ps[2]; d4 = -8 * ps[3]
        Y = yl * Ph
        branches.append(dict(Y=Y, L=[d1 / Y, d2 / (yl * Y), d3 / (yl * yl * Y), d4 / (yl ** 3 * Y)], xi=(Ph - s * d1) / Y))
    out = None
    for b in branches:
        L1, L2, L3, L4 = b['L']
        tau = [None, p * L1, p * (p - 1) * L1 * L1 + p * L2,
               p * (p - 1) * (p - 2) * L1 ** 3 + 3 * p * (p - 1) * L1 * L2 + p * L3,
               p * (p - 1) * (p - 2) * (p - 3) * L1 ** 4 + 6 * p * (p - 1) * (p - 2) * L1 * L1 * L2 + p * (p - 1) * (3 * L2 * L2 + 4 * L1 * L3) + p * L4]
        r = dict(tau=tau, T0=28 * b['Y'] ** p, Y=b['Y'], xi=b['xi'], L1=L1)
        if out is None: out = r
        else:
            out = dict(tau=[None] + [out['tau'][m].union(r['tau'][m]) for m in range(1, 5)], T0=out['T0'].union(r['T0']),
                       Y=out['Y'].union(r['Y']), xi=out['xi'].union(r['xi']), L1=out['L1'].union(r['L1']))
    return out


# ------------------------------------------------------------------------------------------------------------ data
def load_slab(i):
    import numpy as np
    z = np.load(os.path.join(T6DIR, 'jets', f'jets_{i:03d}.npz'))
    J = {k: z[k] for k in z.files}
    R = json.load(open(os.path.join(T6DIR, 'residuals.json')))
    rec = [r for r in R['per_slab'] if r['slab'] == i][0]
    return J, rec


def cell_data(J, e, c, sc):
    """arb hulls of the base coefficients on sub-cell (e, c) over the slab"""
    k = [H(J['k%d' % m][e, c]) for m in range(6)]
    r = [H(J[nm][e, c]) for nm in ('r', 'r1', 'r2', 'r3')]
    Q, A, A1 = H(J['Q'][e, c]), H(J['A'][e, c]), H(J['A1'][e, c])
    k1 = k[1:]; k2 = k[2:]
    b2 = [-v for v in jmul(k, k, 4)]
    c0 = [-3 * a - 4 * b for a, b in zip(jmul(k1, k1, 3), jmul(k, k2, 3))]
    kk1 = jmul(k, k1, 4); k1k1 = jmul(k1, k1, 3); kr = jmul(k, r, 3)
    beta = [Q - sc['c'] * X0 - 6 * kk1[0] + A, k1k1[0] - 6 * kk1[1] + A1]
    for m in range(2, 5): beta.append(k1k1[m - 1] - 6 * kk1[m] + kr[m - 1])
    V = {}
    for m in range(5): V['B%d' % m] = beta[m]; V['D%d' % m] = b2[m]
    for m in range(4): V['C%d' % m] = c0[m]
    # Frenet jets of the modulation directions (p, n) up to order 4
    P = [H(J['p'][e, c])]; Qq = [H(J['q'][e, c])]; N = [H(J['n'][e, c])]; T1 = [H(J['t1'][e, c])]
    for m in range(4):
        P.append(-sum((comb(m, l) * k[l] * Qq[m - l] for l in range(m + 1)), arb(0)))
        Qq.append((arb(1) if m == 0 else arb(0)) + sum((comb(m, l) * k[l] * P[m - l] for l in range(m + 1)), arb(0)))
        N.append(-sum((comb(m, l) * k[l] * T1[m - l] for l in range(m + 1)), arb(0)))
        T1.append(sum((comb(m, l) * k[l] * N[m - l] for l in range(m + 1)), arb(0)))
    return dict(V=V, k=k, P=P, N=N, beta=beta)


def pieces(J, xmin, xmax, splits, nsub_in, nsub_region):
    """yield (e, c, xint(arb), length(float, upper), fixed) for sub-cells meeting [xmin, xmax]; fixed sub-cells are split
    at the points in `splits` and, inside nsub_region = (a, b), into nsub_in equal parts"""
    import numpy as np
    X = J['x']; Ne = X.shape[0]
    for e in range(Ne):
        ew = X[e, 3, 1] - X[e, 0, 0]; w = X[e, :, 1] - X[e, :, 0]
        fixed = abs(w.sum() - ew) < 1e-6 * ew + 1e-8
        for c in range(4):
            lo, hi = float(X[e, c, 0]), float(X[e, c, 1])
            if hi <= xmin or lo >= xmax: continue
            if not fixed:
                # moving element: knots XT (1 + C l^3)^(1/3) increase with tau, so x(s, tau) is increasing in tau and the hull
                # ends are x(s_c, tau0), x(s_{c+1}, tau1); the true sub-cell length at any tau is h_e(tau)/4 <= h_e(tau1)/4
                hs = [X[e, cc + 1, 1] - X[e, cc, 1] for cc in range(3)] + [X[e, cc + 1, 0] - X[e, cc, 0] for cc in range(3)]
                ltrue = max(hs) * (1 + 1e-8)
                assert ltrue > 0 and ltrue <= (hi - lo) * (1 + 1e-8), (e, c, ltrue, hi - lo)
                u0, u1 = max(lo, xmin), min(hi, xmax)
                yield e, c, iv(u0, u1), math.nextafter(min(ltrue, u1 - u0), math.inf), False, u0, u1
                continue
            pts = [max(lo, xmin), min(hi, xmax)]
            pts += [s for s in splits if pts[0] < s < pts[1]]
            a, b = nsub_region
            if hi > a and lo < b:
                pts += [lo + (hi - lo) * q / nsub_in for q in range(1, nsub_in)]
            pts = sorted(set(p for p in pts if pts[0] <= p <= pts[1]))
            for u0, u1 in zip(pts[:-1], pts[1:]):
                if u1 > u0: yield e, c, iv(u0, u1), math.nextafter(u1 - u0, math.inf), True, u0, u1


# ------------------------------------------------------------------------------------------------------- per slab
def run_slab(i, t7=None, n3=True):
    _check_sig()
    t7 = t7 or t7_inputs()
    J, rec = load_slab(i)
    S = rec['scalars']
    sc = dict(c=iv(*S['range_c']), b=iv(*S['range_b']), ell=iv(*S['range_ell']), Qend=iv(*S['range_Q_end']))
    Dabs = S['sup_A_end'] + S['sup_eps_ell']                          # |b1(l) + l_tau| = |A(l) + eps_l|
    sc['D'] = iv(-Dabs, Dabs)
    sc['beta_end'] = sc['Qend'] - sc['c'] * X0 + iv(-S['sup_A_end'], S['sup_A_end'])
    yl = sc['ell'] - X0
    yl_lo = L(yl)
    LN = arb(YB / YA).log()
    x_a = arb(YA) + X0; x_b = arb(YB) + X0
    xa_f, xb_f = U(x_a), U(x_b)
    F = forms()
    muE = arb(t7['muE']); Ke = arb(t7['Ke'])
    Fx = lambda x: x + 2 * x ** 3 / 3 + x ** 5 / 5
    Fxa = Fx(x_a)

    # ---- rho-norm lower bounds of chi g, g = (p, n, V*) on [0, 34.25] (tip; fixed elements)
    nrm2 = [arb(0), arb(0), arb(0)]
    cache = {}
    def cd(e, c):
        if (e, c) not in cache: cache[(e, c)] = cell_data(J, e, c, sc)
        return cache[(e, c)]
    for e, c, xi, ln, fixed, u0, u1 in pieces(J, 0.0, XTIP, [], NSUB_TIP, (0.0, XTIP)):
        xh = arb(u1)
        rho = (-(1 + xh * xh) ** (arb(2) / 3) / 4).exp()                 # rho decreasing: inf at the right end
        yv = u1 - 2.2490776515 + 1e-9
        u = max(0.0, min(1.0, (32.0 - yv) / 3.0))
        chi = arb(u) ** 4 * (35 - 84 * arb(u) + 70 * arb(u) ** 2 - 20 * arb(u) ** 3) if u < 1 else arb(1)
        lnl = pos(arb(u1) - arb(u0) - arb('1e-9'))
        for q, nm in enumerate(('p', 'n', 'Vstar')):
            g = H(J[nm][e, c]); glo = absb(g).lower()
            nrm2[q] += lnl * rho.lower() * chi.lower() ** 2 * glo * glo
    Rs_chi = rec['norms']['Rs_chi']
    Np, Nn, Nv = (L(sqrtp(nrm2[0])), L(sqrtp(nrm2[1])), L(sqrtp(nrm2[2])))
    nV7 = t7.get('nVmin', {}).get(i)                                    # T7's certified lower bound of ||chi V_s||, if present
    nV = max(Nv, nV7 or 0.0) - math.nextafter(Rs_chi, math.inf)        # n = ||chi V~|| >= ||chi V_s|| - ||chi r||
    assert nV > 0
    Np_a, Nn_a, nV_a = arb(Np), arb(Nn), arb(nV)

    # ---- far field
    out = dict(slab=i, tau0=rec['tau0'], k=rec['k'], norms=dict(Np=Np, Nn=Nn, nV=nV, nV_own=Nv, nV_T7=nV7, Rs_chi=Rs_chi), lev={},
               scalars=dict(Dabs=Dabs, ell=S['range_ell'], c=S['range_c'], b=S['range_b'],
                            closure_q=S['sup_closure_q'], turning_t1=S['sup_turning_t1'], D_end_T6=S['sup_D_end']))
    xmax = float(sc['ell'].upper())
    xb_s = math.nextafter(xb_f, math.inf)                              # split point >= x_b: left pieces are ramp
    xa_s = L(x_a)                                                      # far field starts at x_a (T = 0 below)
    far = [(e, c, xi, ln, fx, ('core' if u0 >= xb_s else 'ramp'), u1)
           for e, c, xi, ln, fx, u0, u1 in pieces(J, xa_s, xmax + 1.0, [xa_s, xb_s], NSUB_RAMP, (xa_s, xb_s))]
    assert all(u1 <= xb_s for e, c, xi, ln, fx, u0, u1 in pieces(J, xa_s, xmax + 1.0, [xa_s, xb_s], NSUB_RAMP, (xa_s, xb_s)) if u0 < xb_s)
    Kbar = 0.0
    for e, c, xi, ln, fixed, kind, u1 in far: Kbar = max(Kbar, U(absb(H(J['k0'][e, c]))))
    # e-bound pieces
    eps0 = Ke + sqrtp(Fx(x_b) - Fxa) / muE.sqrt()
    phi1 = arb('0.8828125')                                            # Phi(1) = 1/2 + (1 - 15/64)/2
    for j in range(JMAX + 1):
        p = 2 * j - 2
        Fj = F[j]
        npairs = sum(1 for (a, b) in Fj['pairs'] if a == j + 1)
        nj = npairs + 1
        mj = (arb('0.5') / phi1) ** p if p > 0 else arb(1)
        epsj = arb(Dabs) / mj                                          # Lemma B
        dFj = A_S3 * nj / 8 + epsj * (2 / yl + arb(5) / 4)
        kN = arb(0); kF = arb(0); n3d = {}
        lev = dict(nj=nj, k={}, c={}, r={}, d={}, eq=0.0, A1=0.0, A2=0.0, re=0.0, gam=0.0, gr=0.0, gha=0.0, ghb=0.0)
        # ------------------------------------------------------------------ N3 (level 3 nonlinear classes, placeholders)
        if j == 3 and n3:
            etaF = arb('0.01'); mu = arb('1e-5')
            NLF = arb(os.environ.get('NLFAC', '1'))      # sensitivity factor on the N3 majorants (default 1)
            a0 = 20 * etaF * NLF; a1 = 40 * (etaF / x_a + arb(Kbar)) * NLF
            s0 = a0 * etaF / x_a ** 2; s1 = a1 * etaF / x_a ** 2
            ell0 = ell1 = r0 = r1 = 0.0
        # ------------------------------------------------------------------ core pass
        Ssup = -math.inf; Ssup_y = None
        I_pe = arb(0); I_peH = arb(0); I_ne = arb(0); I_neH = arb(0)
        gam2 = arb(0); gr2 = arb(0); gha2 = arb(0); ghb2 = arb(0); cmax = {}
        tau2p = 0.0
        for e, c, xi, ln, fixed, kind, u1 in far:
            if kind != 'core': continue
            y = xi - X0
            D = cd(e, c); V = dict(D['V'])
            w = core_weight(y, yl, p)
            tau = w['tau']
            V['T0'] = arb(1)
            for m in range(1, 5): V['T%d' % m] = tau[m]
            V['c'] = sc['c']; V['y'] = y
            qn = Fj['qn'](V)
            xi0 = p * (sc['D'] - sc['beta_end']) * w['xi']
            PT = 4 * tau[2] - 2 * V['D0']
            P2T = 4 * tau[4] - 2 * (V['D2'] + 2 * V['D1'] * tau[1] + V['D0'] * tau[2])
            t2p = pos(tau[2])
            rate = xi0 + qn + P2T / 2 + sqr(PT) / (4 * arb(DP[j])) + (A_S3 * nj / 4 + epsj / 2) * t2p + NU * absb(tau[2])
            ru = U(rate)
            if ru > Ssup: Ssup, Ssup_y = ru, (L(y), U(y))
            tau2p = max(tau2p, U(t2p))
            Y = w['Y']; T0 = w['T0']; lnA = arb(ln)
            for (a, b), f in Fj['pairs'].items():
                pT = f(V)                            # p / T_j  (T-symbols = tau)
                if a == j and b >= 0:
                    cmax[b] = max(cmax.get(b, 0.0), U(absb(pT) * Y ** (j - b)))
                elif a == j and b == -1:
                    yy = y
                    H2 = (yy ** 3 - arb(YB) ** 3) / 84
                    I_pe += lnA * sqr(pT) * T0; I_peH += lnA * sqr(pT) * T0 * pos(H2)
                elif a == j + 1 and b == -1:
                    H2 = (y ** 3 - arb(YB) ** 3) / 84
                    I_ne += lnA * sqr(pT) * T0; I_neH += lnA * sqr(pT) * T0 * pos(H2)
                else:
                    raise RuntimeError(f'unexpected pair {(a, b)} at level {j}')
            # modulation
            g1 = D['P'][j + 1] / Np_a; g2 = D['N'][j + 1] / Nn_a; g3 = (sc['c'] * D['P'][j + 1] + sc['b'] * D['N'][j + 1]) / nV_a
            gam2 += lnA * T0 * (sqr(g1) + sqr(g2) + sqr(g3))
            grf = (D['k'][1] if j == 0 else D['k'][j + 1]) / nV_a
            gr2 += lnA * T0 * sqr(grf); gha2 += lnA * sqr(tau[1]) * T0 * sqr(grf); ghb2 += lnA * absb(tau[2]) * T0 * sqr(grf)
            if j == 3 and n3:
                ell0 = max(ell0, U(4 * absb(w['L1']) * a0 * etaF / xi ** 2)); ell1 = max(ell1, U(4 * absb(w['L1']) * a1 * etaF / xi ** 2))
        # ------------------------------------------------------------------ ramp pass
        Kramp = -math.inf; Kat = None; rmax = {}; supPe = 0.0
        ramp_cells = [(e, c, xi, ln, u1) for e, c, xi, ln, fixed, kind, u1 in far if kind == 'ramp']
        for e, c, xi, ln, u1 in ramp_cells:
            y = xi - X0
            D = cd(e, c); V = dict(D['V'])
            w = ramp_weight(y, p, LN)
            T = w['T']; Tt = w['Tt']
            for m in range(5): V['T%d' % m] = T[m]
            V['c'] = sc['c']; V['y'] = y
            Winv = (1 + xi * xi) ** 2
            lnA = arb(ln)
            eh = Ke + sqrtp(Fx(arb(u1)) - Fxa) / muE.sqrt()
            # pairs on the ramp
            for (a, b), f in Fj['pairs'].items():
                pv = f(V)
                if a == j and b >= 0:
                    rmax[b] = max(rmax.get(b, 0.0), U(absb(pv) * Winv))
                elif a == j and b == -1:
                    supPe = max(supPe, U(absb(pv) * Winv))
                    lev.setdefault('_Iramp_e', arb(0))
                    lev['_Iramp_e'] = lev['_Iramp_e'] + lnA * absb(pv) * eh * eh
                elif a == j + 1 and b == -1:
                    Vt = dict(V)
                    for m in range(5): Vt['T%d' % m] = Tt[m]
                    p2T = sqr(f(Vt)) / (28 * w['Z'])
                    lev.setdefault('_Iramp_ne', arb(0))
                    lev['_Iramp_ne'] = lev['_Iramp_ne'] + lnA * p2T * eh * eh
            # modulation (ramp part)
            g1 = D['P'][j + 1] / Np_a; g2 = D['N'][j + 1] / Nn_a; g3 = (sc['c'] * D['P'][j + 1] + sc['b'] * D['N'][j + 1]) / nV_a
            T0 = T[0]
            gam2 += lnA * T0 * (sqr(g1) + sqr(g2) + sqr(g3))
            grf = (D['k'][1] if j == 0 else D['k'][j + 1]) / nV_a
            gr2 += lnA * T0 * sqr(grf); gha2 += lnA * sqr(Tt[1]) / (28 * w['Z']) * sqr(grf); ghb2 += lnA * absb(T[2]) * sqr(grf)
            if j == 3 and n3:
                ell0 = max(ell0, U(4 / y * a0 * etaF / xi ** 2)); ell1 = max(ell1, U(4 / y * a1 * etaF / xi ** 2))
                r0 = max(r0, U(y * y * a0 * etaF / xi ** 2 * absb(w['eta1']) * (1 + xi * xi)))
                r1 = max(r1, U(y * y * a1 * etaF / xi ** 2 * absb(w['eta1']) * (1 + xi * xi)))
        if j == 3 and n3:
            ell0 = arb(2).sqrt() * arb(ell0); ell1 = arb(2).sqrt() * arb(ell1)
            r0 = arb(224).sqrt() * arb(r0); r1 = arb(224).sqrt() * arb(r1)
            kN = ell0 / 2 + r0 * mu / 2 + 3 * s0 + 3 * s1 / 2
            kF = ell0 / 2 + ell1 + r1 * mu / 2 + 3 * s1 / 2
            dFj = dFj + kF + kN / 4
            n3d = dict(s0=U(s0), s1=U(s1), ell0=U(ell0), ell1=U(ell1), r0=U(r0), r1=U(r1), kN=U(kN), kF=U(kF), a0=U(a0), a1=U(a1),
                       Kbar=Kbar, k33_add=U((r0 + r1) / (2 * mu)), note='placeholders a0 = 20 eta_F, a1 = 40(eta_F/x_a + sup|k|), eta_F = 1e-2, mu = 1e-5 [T7: n2-n4]')
            Ssup = Ssup + U(kN / 2 * arb(tau2p))           # core: kN/2 (T'')^+ / T  <= kN/2 sup tau2^+
        lam = L(sc['c']) - Ssup - U(dFj)
        lamA = arb(lam) + dFj
        for e, c, xi, ln, u1 in ramp_cells:
            y = xi - X0
            D = cd(e, c); V = dict(D['V'])
            w = ramp_weight(y, p, LN)
            T = w['T']; Tt = w['Tt']
            for m in range(5): V['T%d' % m] = T[m]
            V['c'] = sc['c']; V['y'] = y
            Winv = (1 + xi * xi) ** 2
            q = Fj['q'](V)
            P2 = 4 * T[4] - 2 * (V['D2'] * T[0] + 2 * V['D1'] * T[1] + V['D0'] * T[2])
            P2T = 16 * sqr(Tt[2]) / (28 * w['Z']) - 16 * V['D0'] * T[2] + 4 * sqr(V['D0']) * T[0]
            qh = q + P2 / 2 + P2T / (4 * arb(DP[j])) + (A_S3 * nj / 4 + epsj / 2 + kN / 2) * pos(T[2]) + NU * absb(T[2])
            kv = U(pos(qh + lamA * T[0]) * Winv)
            if kv > Kramp: Kramp, Kat = kv, (L(y), U(y))
        # ------------------------------------------------------------------ assemble
        Iq = lambda key: lev.pop(key, arb(0))
        I_re = Iq('_Iramp_e'); I_rne = Iq('_Iramp_ne')
        A1 = eps0 * sqrtp(I_pe); A2 = sqrtp(I_peH)
        lev['A1'] = U(A1); lev['A2'] = U(A2)
        lev['re'] = U(sqrtp(arb(supPe) * I_re)) * ZETA_FAC if supPe > 0 else 0.0
        if npairs:
            eq = (2 * eps0 ** 2 * I_ne + I_rne) / (2 * A_S3); dn = 2 * I_neH / (2 * A_S3)
            lev['eq'] = U(eq); lev['d'] = {0: U(dn)}
        else:
            lev['d'] = {}
        lev['c'] = {b: v for b, v in cmax.items()}
        if j >= 1: lev['c'][0] = lev['c'].get(0, 0.0) + 2 * lev['A2']
        lev['r'] = {b: v * ZETA_FAC for b, v in rmax.items()}
        lev['k'] = {j: Kramp * ZETA_FAC}
        if j == 3 and n3: lev['k'][j] += n3d['k33_add'] * ZETA_FAC
        lev['k_at_y'] = Kat
        lev['sigma'] = Ssup + U(dFj); lev['lam'] = lam; lev['rate_at_y'] = Ssup_y
        # TOPLOSS: top-coefficient loss 2 sup|A4 - A4(0)| from the certified bend box (p2/boxes.py; default 0)
        lev['Dcoef'] = L(arb(2) - arb(DP[j]) - A_S3 * nj / 2 - arb(DS[j]) - epsj - kN - arb(os.environ.get('TOPLOSS', '0')))
        lev['DS'] = DS[j]; lev['DP'] = DP[j]
        assert lev['Dcoef'] > 0, (j, lev['Dcoef'])
        lev['gam'] = U(sqrtp(gam2)); lev['gr'] = U(sqrtp(gr2)); lev['gha'] = U(sqrtp(gha2)); lev['ghb'] = U(sqrtp(ghb2))
        lev['bdry'] = dict(eps=U(epsj), m=L(mj), tau2p_core=tau2p)
        if n3d: lev['n3'] = n3d
        nm = rec['norms']
        lev['RF'] = nm[f'RF{j}_opt4n']; lev['RFa'] = nm[f'RFa{j}_opt4n']; lev['RFb'] = nm[f'RFb{j}_opt4n']
        lev['RF_int'] = rec['norms_int'].get(f'RF{j}_opt4n')
        out['lev'][j] = lev
    out['t7_source'] = t7['_source']; out['t7_sources'] = t7['_sources']
    return out


# ------------------------------------------------------------------------------------------------------------ CLI
def parse_slabs(spec):
    out = []
    for part in spec.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        elif part: out.append(int(part))
    return out


def work(slabs, outdir):
    os.makedirs(outdir, exist_ok=True)
    t7 = t7_inputs()
    for i in slabs:
        fj = os.path.join(outdir, f'slab_{i:03d}.json')
        pv = provenance(i, outdir, t7)
        if os.path.exists(fj) and all((json.load(open(fj)).get('args') or {}).get(k) == pv[k] for k in PVKEYS): continue
        t0 = time.time(); c0 = time.process_time()
        o = run_slab(i, t7)
        o['args'] = pv
        o['cpu'] = time.process_time() - c0; o['wall'] = time.time() - t0
        tmp = fj + '.tmp'
        with open(tmp, 'w') as f: json.dump(o, f, indent=1)
        os.replace(tmp, fj)
        lv = o['lev']
        print(f"slab {i}: k={o['k']:.4g} cpu={o['cpu']:.1f}s lam=" + ','.join(f"{lv[j]['lam']:.4f}" for j in lv) +
              ' k_jj=' + ','.join(f"{lv[j]['k'][j]:.3e}" for j in lv) + f" n={o['norms']['nV']:.3e}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slabs', required=True)
    ap.add_argument('--out', default=os.path.join(HERE, 'constants', 'slabs'))
    ap.add_argument('--jobs', type=int, default=1)
    a = ap.parse_args()
    slabs = parse_slabs(a.slabs)
    if a.jobs <= 1:
        work(slabs, a.out); return
    os.makedirs(a.out, exist_ok=True)
    procs = []
    for w in range(a.jobs):
        sub = slabs[w::a.jobs]
        if not sub: continue
        cmd = [sys.executable, os.path.abspath(__file__), '--slabs', ','.join(map(str, sub)), '--out', a.out]
        log = open(os.path.join(a.out, f'worker{w}.log'), 'a')
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT))
    for pr in procs: pr.wait()


if __name__ == '__main__':
    main()
