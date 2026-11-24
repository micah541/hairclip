"""T7: certified closure rate rho_c of T5 (E0) + (E2) + mu_e [(I0) + (I2)] in linear mode (default) or, with --nl FILE, in
nonlinear mode: the terms of nlterms7.py (T7 Lemma 4c: N, Q and the nonlinear part of (M), reduced to quadratic terms on the
bootstrap box) are added to the 'bad' side, so that Q(z) - Ds - NL(z) >= rho_c L(z) with the same reserve Ds; per certificate interval, with the
certified T7 constants (worst case over the interval's slabs).  T5 2.4 form (closure.py), (E2) modulation in AB form (G02 + G_x) A B,
dissipation reserve Ds (kappa_D = 1/2 for Lemma J):
   Q(z) := -(E0 + E2 + mu_e (I0 + I2)) - Ds  >=  rho_c L(z),   L = A^2 + eps B^2 + mu_e (V0^2 + eps_I V2^2),
for all shapes z = (A, B, D4, V0, V2, V4) >= 0 with W1^2 <= A B + c1 A^2, W3^2 <= B D4 + c1 B^2, V1^2 <= V0 V2 + c1W V0^2,
V3^2 <= V2 V4 + c1W V2^2.
Proof: every product x y of nonnegative atoms with a sqrt atom is bounded by (t x^2 + y^2/t)/2 (t > 0 fixed, chosen by float
optimisation), squares of sqrt atoms by their interpolation bounds; products of basic atoms are kept.  Then
E0 + E2 + mu_e(I0 + I2) + Ds + rho L <= z^T M z, and -M - rho diag(L) is certified positive definite (arb LDL^T) on R^6.
One mu_e for the whole run (maximises min_I rho_c).  usage: python3 closure_cert.py  -> closure_cert.json (and t8req.json update)"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import minimize
from flint import arb, ctx
ctx.prec = 128
EPS, EPSI, C1, PHI1, PHI2 = 1e-5, 3e-5, 0.5804, 1.083, 1.3
NAMES = ['A', 'B', 'D4', 'V0', 'V2', 'V4']; IX = {n: i for i, n in enumerate(NAMES)}
SQRT = {'W1': (('A', 'B', 0.5), ('A', 'A', 'C1')), 'W3': (('B', 'D4', 0.5), ('B', 'B', 'C1')),
        'V1': (('V0', 'V2', 0.5), ('V0', 'V0', 'c1W')), 'V3': (('V2', 'V4', 0.5), ('V2', 'V2', 'c1W'))}


def terms(c, mu):
    """list of (coef, atom1, atom2) with coef >= 0 for the 'bad' part and the exact negative squares; c = constants"""
    T = []
    add = lambda k, a, b: T.append((k, a, b))
    # E0
    add(-c['om'], 'A', 'A'); add(-c['th0'], 'B', 'B'); add(c['strip'], 'V0', 'V0'); add(c['sA0'], 'A', 'A'); add(2 * c['sA0'], 'A', 'W1')
    # E2 (times eps)
    e = EPS
    add(-0.5 * e, 'D4', 'D4'); add(e * c['G22'], 'B', 'B'); add(e * c['G12'], 'W1', 'B'); add(e * (c['G02'] + c['Gx']), 'A', 'B')
    for dk, v in ((c['d3'], 'V2'), (c['d2'], 'V1'), (c['d1'], 'V0')):
        for fk, w in ((1.0, 'D4'), (2 * PHI1, 'W3'), (PHI2, 'B')): add(e * dk * fk, v, w)
    # I0 (times mu)
    add(-mu * c['mu1'], 'V0', 'V0'); add(-0.5 * mu, 'V2', 'V2'); add(mu * c['a_f'], 'W1', 'W1'); add(mu * c['b_f'], 'B', 'B'); add(mu * c['Cg0'], 'A', 'V0')
    # I2 (times mu eps_I)
    m2 = mu * EPSI
    add(-0.5 * m2, 'V4', 'V4'); add(-m2 * c['mu2'], 'V2', 'V2')
    for lk, a in ((c['c3'], 'D4'), (c['c2'], 'W3'), (c['c1'], 'B'), (c['c0'], 'W1'), (c['Cg2'], 'A')):
        for rk, b in ((1.0, 'V4'), (2 * c['Psi1'], 'V3'), (2 * c['c1W'], 'V2')): add(m2 * lk * rk, a, b)
    # dissipation reserve
    add(c['th0'] / 4, 'B', 'B'); add(EPS / 8, 'D4', 'D4'); add(mu / 8, 'V2', 'V2'); add(mu * EPSI / 8, 'V4', 'V4')
    return T


def amgm_list(T):
    return [i for i, (k, a, b) in enumerate(T) if (a in SQRT or b in SQRT) and not (a == b)]


def build(T, c, tpar, R=float):
    """matrix M (list of lists, entries of type R) with z^T M z >= sum of terms, AM-GM parameters tpar (dict term index -> t)"""
    M = [[R(0) for _ in range(6)] for _ in range(6)]
    def addsq(k, a):                          # k * a^2 (a basic or sqrt atom: square -> interpolation bound)
        if a in SQRT:
            for (u, v, w) in SQRT[a]:
                kk = k * (R(c[w]) if isinstance(w, str) else R(w) * 2 if (u != v) else R(w))
                if u == v: M[IX[u]][IX[u]] += kk
                else:
                    M[IX[u]][IX[v]] += kk / 2; M[IX[v]][IX[u]] += kk / 2
        else: M[IX[a]][IX[a]] += k
    for i, (k, a, b) in enumerate(T):
        k = R(k)
        if a == b: addsq(k, a); continue
        if a in SQRT or b in SQRT:
            t = R(tpar[i]); addsq(k * t / 2, a); addsq(k / (2 * t), b)
        else:
            M[IX[a]][IX[b]] += k / 2; M[IX[b]][IX[a]] += k / 2
    return M


def Lvec(mu, R=float): return [R(1), R(EPS), R(0), R(mu), R(mu * EPSI), R(0)]


def rho_max(M, mu):
    """largest rho with -M - rho L >= 0 (float, bisection)"""
    Mm = -np.array(M); L = np.diag(Lvec(mu))
    if np.linalg.eigvalsh(Mm).min() <= 0: return -1.0
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if np.linalg.eigvalsh(Mm - mid * L).min() > 0: lo = mid
        else: hi = mid
    return lo


def optimise(c, mu, x0=None, extra=()):
    T = terms(c, mu) + [tuple(t) for t in extra]; idx = amgm_list(T)
    def f(lt):
        tp = {i: math.exp(v) for i, v in zip(idx, lt)}
        return -rho_max(build(T, c, tp), mu)
    x0 = np.zeros(len(idx)) if x0 is None else x0
    best = None
    for start in (x0, np.full(len(idx), -2.0), np.full(len(idx), 2.0)):
        r = minimize(f, start, method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-3, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    return -best.fun, best.x, T, idx


def ldl_pd(S):
    """arb LDL^T: True if the symmetric ball matrix S is certified positive definite"""
    n = len(S); Lm = [[arb(0)] * n for _ in range(n)]; D = [arb(0)] * n
    for i in range(n):
        s = S[i][i] - sum((Lm[i][k] ** 2 * D[k] for k in range(i)), arb(0))
        if not float(s.lower()) > 0: return False
        D[i] = s
        for j in range(i + 1, n):
            Lm[j][i] = (S[j][i] - sum((Lm[j][k] * Lm[i][k] * D[k] for k in range(i)), arb(0))) / s
    return True


def certify(c, mu, lt, T, idx, rho):
    tp = {i: float(math.exp(v)) for i, v in zip(idx, lt)}
    ca = {k: (arb(v) if isinstance(v, (int, float)) else v) for k, v in c.items()}
    M = build(T, ca, tp, R=arb)
    L = Lvec(mu, R=arb)
    S = [[-M[i][j] - (arb(rho) * L[i] if i == j else arb(0)) for j in range(6)] for i in range(6)]
    return ldl_pd(S)


def interval_constants(C, cert):
    lo, hi = cert['slabs']; S = [C['slabs'][str(s)] for s in range(lo, hi + 1)]
    W = C['meta']['weight_consts']
    mx = lambda f: max(f(e) for e in S); mn = lambda f: min(f(e) for e in S)
    return dict(om=cert['c_I'], th0=cert['theta'], strip=mx(lambda e: e['strip']['sA'] + e['strip']['sB']), sA0=mx(lambda e: e['strip']['sA0']),
                G22=mx(lambda e: e['E2']['G22']), G12=mx(lambda e: e['E2']['G12']), G02=mx(lambda e: e['E2']['G02']),
                d3=mx(lambda e: e['strip']['d3']), d2=mx(lambda e: e['strip']['d2']), d1=mx(lambda e: e['strip']['d1']),
                mu1=mn(lambda e: e['E_I']['mu1']), mu2=mn(lambda e: e['E_I']['mu2']), a_f=mx(lambda e: e['ramp']['a_f']),
                b_f=mx(lambda e: e['ramp']['b_f']), c3=mx(lambda e: e['ramp']['c3']), c2=mx(lambda e: e['ramp']['c2']),
                c1=mx(lambda e: e['ramp']['c1']), c0=mx(lambda e: e['ramp']['c0']), Psi1=W['Psi1'], c1W=W['c1W'], C1=C1,
                Cg0=mx(lambda e: e['modulation']['Cgamma0']), Cg2=mx(lambda e: e['modulation']['Cgamma2']), Gx=mx(lambda e: e['modulation']['Gx']))


def main(mu=None, out_name='closure_cert.json', work=HERE, nl=None):
    C = json.load(open(os.path.join(work, 'constants.json')))
    certs = [cc for cc in C['certificates']]
    NL = None
    if nl:                                    # nonlinear mode (T7 Lemma 4c): extra terms per interval from nlterms7.py
        D = json.load(open(nl)); NL = {r['id']: r for r in D['intervals']}
        certs = [cc for cc in certs if cc['id'] in NL]
        assert mu is not None, 'nonlinear mode needs the --mu the terms were generated with'
        assert all(len(t) == 3 and t[1] in IX or t[1] in SQRT for r in NL.values() for t in r['terms'])
    consts = [interval_constants(C, cc) for cc in certs]
    if mu is None:
        # common mu_e: maximise the minimum over intervals of the (float) relaxed rate
        def minrate(mu):
            return min(optimise(c, mu)[0] for c in consts[::6])
        grid = [0.01, 0.015, 0.02, 0.025, 0.03, 0.036, 0.045, 0.055]
        vals = [(minrate(m), m) for m in grid]
        print('mu_e scan (subset of intervals):', [(m, round(v, 5)) for v, m in vals], flush=True)
        mu = max(vals)[1]
    out = []; x0 = None
    for cc, c in zip(certs, consts):
        extra = [(float(k), a, b) for k, a, b in NL[cc['id']]['terms']] if NL else ()
        r, lt, T, idx = optimise(c, mu, x0 if not NL else None, extra); x0 = lt
        rc = math.floor(r * 0.995 * 1e5) / 1e5
        ok = certify(c, mu, lt, T, idx, rc)
        out.append(dict(id=cc['id'], slabs=cc['slabs'], k_c=cc['k_c'], mu_e=mu, rho_c=rc if ok else None, rho_float_relaxed=r, certified=ok,
                        amgm_t=[float(math.exp(v)) for v in lt], consts=c, n_nl_terms=len(extra)))
        print(f"cert {cc['id']:2d} (k_c {cc['k_c']:8.1f}): rho_c = {rc:.5f} certified={ok} (relaxed float {r:.5f}); om {c['om']:.4f} Cg0 {c['Cg0']:.3g}", flush=True)
    import provenance as PV
    args = dict(argv=sys.argv[1:], mu_e=mu, work=PV.rel(work), inputs=dict(constants=PV.entry(os.path.join(work, 'constants.json'))),
                code=PV.code([os.path.abspath(__file__), os.path.join(HERE, 'provenance.py')]))
    if nl: args['inputs']['nl_terms'] = PV.entry(nl)
    res = dict(mu_e=mu, rho_c_min=min(o['rho_c'] for o in out if o['rho_c'] is not None), all_certified=all(o['certified'] for o in out), per_interval=out,
               nonlinear=(os.path.basename(nl) if nl else None), args=args)
    json.dump(res, open(os.path.join(work, out_name), 'w'), indent=1)
    print('mu_e', mu, 'min rho_c', res['rho_c_min'], 'all certified', res['all_certified'])
    return res


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--mu', type=float, default=None); ap.add_argument('--out', default='closure_cert.json')
    ap.add_argument('--work', default='.'); ap.add_argument('--nl', default=None, help='nl_terms.json of nlterms7.py (nonlinear mode)')
    A = ap.parse_args(); W = A.work if os.path.isabs(A.work) else os.path.normpath(os.path.join(HERE, A.work))
    main(A.mu, A.out, W, A.nl if (A.nl is None or os.path.isabs(A.nl)) else os.path.join(HERE, A.nl))
