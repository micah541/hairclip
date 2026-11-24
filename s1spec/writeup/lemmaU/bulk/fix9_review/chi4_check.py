#!/usr/bin/env python3
"""fix9_review/chi4_check.py -- independent exact check of the level-4 chart structure (FLUX2_DESIGN 3.4; fix8 chart4_shell.py X3).
Independent construction: jet symbols and a total-derivative operator (no sympy Function objects), trigonometric functions replaced by
algebraic symbols Sn = sin(delta), Cn = cos(delta) (2 sin^2(delta/2) = 1 - Cn) with d/dsigma Sn = k Cn, d/dsigma Cn = -k Sn.  All zero tests are
exact identities of rational functions in independent symbols (sufficient: if an expression vanishes with Sn, Cn independent it vanishes on Sn^2+Cn^2=1).
chi_k = (T^k - D^k) b + T^k Phi, T = (1 + eps) D, D = sigma d/dsigma (the p^k (U(sigma) - U(Y)) part has no jets of u beyond u itself).
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B chi4_check.py"""
import sympy as sp, json, os, sys, time, resource
t0 = time.time()
s, e, Sn, Cn, u0 = sp.symbols('s e Sn Cn u0')
kk = sp.symbols('k0:6')            # k, k', ..., k^(5)
b = sp.symbols('b0:7')             # b = u', b_i = d^i b
eta = (1 - Cn) + b[0] * Sn + u0 * kk[0] * Cn
eps = (e / s + eta) / (1 - eta)
Phi = ((1 - kk[0] * u0 + b[0] ** 2) * Sn + kk[0] * u0 * b[0] * Cn) / ((1 - kk[0] * u0) * Cn - b[0] * Sn)
def Dt(F):
    r = sp.diff(F, s) + sp.diff(F, e) * (-eta) + sp.diff(F, Sn) * kk[0] * Cn + sp.diff(F, Cn) * (-kk[0] * Sn) + sp.diff(F, u0) * b[0]
    for i in range(5): r += sp.diff(F, kk[i]) * kk[i + 1]
    for i in range(6): r += sp.diff(F, b[i]) * b[i + 1]
    return r
T = lambda F: (1 + eps) * s * Dt(F); D = lambda F: s * Dt(F)
res = {}; chi = {}
Fb, Fp, Db = b[0], Phi, b[0]
for n in range(5):
    if n: Fb, Fp, Db = T(Fb), T(Fp), D(Db)
    chi[n] = Fb - Db + Fp
def iszero(x): return sp.cancel(sp.together(x)) == 0
ok = True
def chk(name, val):
    global ok
    res[name] = bool(val); ok = ok and bool(val); print('%-70s %s  (%.0fs, %d MB)' % (name, bool(val), time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024), flush=True)
c4 = chi[4]
chk('chi_4 free of b5, b6 (u^(6), u^(7)) and of k^(5)', not ({b[5], b[6], kk[5]} & c4.free_symbols))
for i, j in ((2, 3), (2, 4), (3, 3), (3, 4), (4, 4)):
    chk('d^2 chi_4 / db%d db%d == 0' % (i, j), iszero(sp.diff(c4, b[i], b[j])))
Q = sp.diff(c4, b[2], 2)
chk('d^3 chi_4 / db2^3 == 0', iszero(sp.diff(Q, b[2])))
chk('B_422 = (1/2) d^2 chi_4/db2^2 is NOT identically zero (the quadratic term exists)', not iszero(Q))
chk('B_422, B_42, B_43, B_44, A_4 depend on (sigma, e, delta, k..k\'\'\'\', u, b, b_1) only',
    not ({b[2], b[3], b[4]} & (sp.cancel(sp.together(Q)).free_symbols | sp.diff(c4, b[3]).free_symbols | sp.diff(c4, b[4]).free_symbols | sp.diff(c4, b[2]).subs(b[2], 0).free_symbols)))
ray = {Sn: 0, Cn: 1, e: 0, **{kk[i]: 0 for i in range(6)}}
for n in range(5): chk('chi_%d == 0 on a ray' % n, iszero(chi[n].subs(ray)))
for n in range(4):
    chk('chi_%d (k<=3): affine in (b2, b3), no b4' % n, b[4] not in chi[n].free_symbols and all(iszero(sp.diff(chi[n], b[i], b[j])) for i, j in ((2, 2), (2, 3), (3, 3))))
# explicit leading structure of the quadratic coefficient at a ray-like expansion: B_422 to first order in Sn (diagnostic, exact)
Qs = sp.cancel(sp.together(Q / 2))
lead = sp.series(Qs.subs({Cn: 1, e: 0, kk[0]: 0}), Sn, 0, 2).removeO()
res['B422_first_order_in_sin_delta_at_k=0_e=0'] = str(sp.factor(lead))
print('B_422 ~', res['B422_first_order_in_sin_delta_at_k=0_e=0'])
res['all_pass'] = ok; res['args'] = dict(argv=sys.argv, sympy=sp.__version__, env={'PYTHONHASHSEED': os.environ.get('PYTHONHASHSEED')})
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chi4_check.json'), 'w'), indent=1, sort_keys=True)
print('ALL PASS' if ok else 'FAIL')
