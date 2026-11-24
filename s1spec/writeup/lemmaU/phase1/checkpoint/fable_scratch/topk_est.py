# Reviewer scratch (read-only on project): port of lemmaS/rigor/ramp_arb.py ramp_level, plus the k = j+1 table it does not compute.
import sys, json, math, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/phase1/ghat')
import sympy as sp
from flint import arb, ctx
import ghat_arb as GA, arbio_local as AL
ctx.prec = 64
A = GA.A; Q_ = sp.Rational
LS = '/home/claude/s1spec/writeup/lemmaS'
T = pickle.load(open(LS + '/rigor/tables_exact.pkl', 'rb'))
WE = AL.load(json.load(open(LS + '/rigor/weights_exact.json'))['res'])
REF = AL.load(json.load(open(LS + '/rigor/ramp_tables.json')))
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}; BELL = {0: 1, 1: 1, 2: 2, 3: 5}
LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); gam = 1
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2s = Sx**2
Dq = {}
for m in (1, 2, 3, 4):
    R = sum(sp.binomial(m, l)*sp.diff(e2s, t_, l)/Ls**l*(-gam)**(m - l) for l in range(1, m + 1)); f = sp.expand(R); Dq[m] = []
    for q_ in range(4):
        Dq[m].append(sp.lambdify((t_, Ls), f, modules=[{}])); f = sp.expand(sp.diff(f, t_)/Ls - gam*f)
def yy(t): return 14*(t*LN2).exp()
X0 = arb(22490776515)/10**10
om = lambda t: (-(arb(gam))*((yy(t)/28).log())).exp()
eps = lambda t: yy(t)**-4
Wt = lambda t: (1 + (yy(t) + X0)**2)**-2
def expand_ab(a, b, d):
    if a == b: return {(a, d): 1}
    if b == a + 1: return {(a, d + 1): arb(1)/2}
    out = {}
    for key, v in list(expand_ab(a + 1, b - 1, d).items()) + list(expand_ab(a, b - 1, d + 1).items()): out[key] = out.get(key, 0) + v
    return out
Abr = GA.ab_exact(Q_(5, 100), Q_(9, 10), 6)
def level(j, ks):
    n = j + 2; pairs = [(a, b, m) for (a, b), m in T[j].items() if (a, b) != (n, n)]
    def cd(mons, d, t):
        tot = arb(0)
        for (m, i, r), co in mons.items():
            if m == 0: continue
            for p_ in range(d + 1):
                for q_ in range(d + 1 - p_):
                    u = d - p_ - q_; mult = math.factorial(d)//(math.factorial(p_)*math.factorial(q_)*math.factorial(u))
                    tot += A(co)*mult*4**p_*Abr[(i, r + u)]*abs(Dq[m][q_](t, LN2))*om(t)
        return tot*eps(t)
    wts = []
    for (a, b, mons) in pairs:
        if b <= j: wts.append((mons, {(a, 0): 1} if a == b else {(a, 0): arb(1)/2, (b, 0): arb(1)/2}))
        else: wts.append((mons, expand_ab(a, b, 0)))
    keys = sorted({kd for _, w in wts for kd in w})
    def q(k, t, pairs_only=False):
        tot = arb(0)
        for mons, w in wts:
            for (kk, d), v in w.items():
                if kk == k: tot += arb(v)*cd(mons, d, t)
        if not pairs_only and k <= j:
            tot += AL.up(WE['ramp'][str(j)]['Ghat'][str(k)])*eps(t)*(Sx.subs(t_, 0) + arb(1)*0 + (t**3*(10 - 15*t + 6*t*t))**2)*om(t)
            if k == j: tot += AL.up(WE['ramp'][str(j)]['YC'][str(j + 1)])*abs(((sp.lambdify(t_, sp.diff(e2s, t_, 2), modules=[{}])(t))/LN2**2 - 2*(gam + 4)*(sp.lambdify(t_, sp.diff(e2s, t_, 1), modules=[{}])(t))/LN2)*om(t)*eps(t))/2
        return tot
    def sup_on(f, n_=1500):
        best = 0.0
        for i in range(n_):
            v = float(abs(f(arb((i + 0.5)/n_))).mid())
            if v > best: best = v
        return arb(best), best
    out = {}
    for k in ks:
        out[k] = {m: sup_on(lambda t, m=m, k=k: q(k, t)*yy(t)**(2*m - 1)/Wt(t)) for m in STIR[k]}
    return out, keys
# final-run data (lemmaS_F10p_M0.52.json args; gamma' from RIGOR_PASS)
J = [6.0273e-12, 8.4673e-12, 1.8708e-12]; sl = [8.6369e-19, 1.2134e-18, 2.6808e-19]; gp = [0.24628897, 0.24617403, 0.14306229]
J4, sJ4, c1W = 6.2362e-08, 8.9363e-15, 0.8477
res = {}
for j in (0, 1, 2):
    Je = [J[m] + sl[m]/gp[j] for m in range(3)]; J4e = J4 + sJ4/gp[j]; Je.append(math.sqrt(Je[2]*J4e) + c1W*Je[2])
    ks = [j, j + 1] if len(sys.argv) < 2 else [j + 1]
    out, keys = level(j, ks)
    print('level', j, 'IBP-reduced keys (k,d):', keys, flush=True)
    for k in ks:
        tab = {m: float(v[0]) for m, v in out[k].items()}; lbs = {m: v[1] for m, v in out[k].items()}
        val = BELL[k]*sum(STIR[k][m]*tab[m]*Je[m] for m in STIR[k])
        ref = {m: float(v) for m, v in REF[str(j)][str(k)].items()} if k <= j and str(k) in REF[str(j)] else None
        print('  k =', k, 'table', {m: '%.6g' % v for m, v in tab.items()}, 'pt-lb', {m: '%.4g' % v for m, v in lbs.items()}, '| ramp_tables.json:', ref and {m: '%.6g' % v for m, v in ref.items()},
              '| contribution to E_%d: %.4g' % (j, val), 'per m:', {m: '%.3g' % (BELL[k]*STIR[k][m]*tab[m]*Je[m]) for m in STIR[k]}, flush=True)
    res[j] = dict(Je=Je)
print('J3e-type data', res)
