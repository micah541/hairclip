#!/usr/bin/env python3
"""Item 2 of the rigor pass: the L^4 class weights and the constants hat-G, recomputed WITHOUT sympy Floats / double arithmetic.
Exact sympy Rationals for the boxes (M0,M1,Umax,p), the table coefficients and the Bell/Faa di Bruno polynomials; the only irrational inputs are the suprema
Fa_k = sup_{|p|<=M0}|d^k (1+p^2)^-2|, Fb_k (arb ball upper endpoints on 600 cells, as in lemmaS_prop1_bounds.sup_rat); all later arithmetic is arb (prec 200), outward.
Ports (line by line) of: lemmaS_prop1_bounds.{table,A_bounds_s}, lemmaS_F1.constants (YC part), lemmaS_L4terms.{Aterms,classify,l4_extra_YC,_old_fallback_impl,ghat_from_YC},
the ramp-L4 block of lemmaS_F2.main (add[], footw, W33, C_eta), the core block of lemmaS_F2.main, and the per-window block of lemmaS_Bbox.py.
Results (arb) -> rigor/weights_exact.json, comparison with the double values of the final run -> printed + weights_exact.json['compare'].

Rerun (single core, a few minutes; run alone):
    cd /home/claude/s1spec/writeup/lemmaS && PYTHONHASHSEED=0 python3 rigor/weights_exact.py
(needs rigor/cap_F10p.json from rigor/capture_F2.py for the ramp comparison; the core/window comparisons use lemmaS_F10p_M0.52.json and lemmaS_Bbox_1e-05_a0.00015_dv0.15.json)"""
import os, sys, json, pickle, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); PAR = os.path.dirname(HERE)
os.chdir(PAR); sys.path.insert(0, PAR); sys.path.insert(0, HERE)
import sympy as sp
from sympy.functions.combinatorial.numbers import bell as bellpoly
from flint import arb, ctx
import lemmaS_F1 as F, lemmaS_prop1_sym as P, lemmaS_prop1_bounds as PB, lemmaS_B3 as B
import arbio
ctx.prec = 200
Q_ = sp.Rational
s = P.s
GAM = 1
PQ = Q_(10246, 10000); UQ = Q_(1, 20); M0Z = Q_(1, 50); KAP = 5; M2STAR = 3

def A(q):
    q = Q_(q); return arb(int(q.p)) / arb(int(q.q))
def ub(x): return x.upper()

# ------------------------------------------------------------------ exact tables (port of PB.table with Rational |coef|)
def table_exact(Q):
    out = {}
    for (a, b), c in Q.items():
        e = sp.expand(sp.simplify(c.subs(P.ep, sp.exp(-4*s)).doit()*sp.exp(4*s)))
        gens = list(e.atoms(sp.Derivative)) + [f for f in e.atoms(sp.Function) if not isinstance(f, sp.Derivative) and f.func.__name__ in ('om', 'A0', 'A1', 'A2')]
        P_ = sp.Poly(e, *gens); tab = {}
        for mon, co in P_.terms():
            m = i = r = None; m = 0
            for gsym, pw in zip(P_.gens, mon):
                if pw == 0: continue
                assert pw == 1
                if isinstance(gsym, sp.Derivative): nm, k = gsym.expr.func.__name__, gsym.derivative_count
                else: nm, k = gsym.func.__name__, 0
                if nm == 'om': m = k
                else: i, r = int(nm[1:]), k
            assert co.is_Rational, co
            tab[(m, i, r)] = tab.get((m, i, r), 0) + abs(co)
        out[(a, b)] = tab
    return out

def build_tables(cache=os.path.join(HERE, 'tables_exact.pkl')):
    if os.path.exists(cache): return pickle.load(open(cache, 'rb'))
    T = {}
    for j in (0, 1, 2):
        Qd, ok = F.div_form(j); assert ok; T[j] = table_exact(Qd)
    Q = pickle.load(open('prop1_forms3.pkl', 'rb'))[3][0]; T[3] = table_exact(B.div_form3(Q, 5))
    pickle.dump(T, open(cache, 'wb')); return T

# ------------------------------------------------------------------ Ab = A_bounds_s (exact)
_pp = sp.Symbol('p'); _a = (1 + _pp**2)**-2; _b = -3*_pp*(1 + _pp**2)**-3
_lam = {}
def _fun(kind, k):
    if (kind, k) not in _lam:
        _lam[(kind, k)] = sp.lambdify(_pp, sp.diff(_a if kind == 'a' else _b, _pp, k), modules=[{}])
    return _lam[(kind, k)]
_supc = {}
def sup_rat_arb(kind, k, M0q, N=600):
    key = (kind, k, M0q)
    if key in _supc: return _supc[key]
    f = _fun(kind, k); M0 = A(M0q); best = arb(0)
    for i in range(N):
        lo = M0*(2*i - N)/N; hi = M0*(2*i + 2 - N)/N
        v = abs(f(lo.union(hi))); u = v.upper()
        if u > best: best = u
    _supc[key] = best; return best

def faa_num(F_, r, V):
    if r == 0: return F_[0]
    return sum((F_[k]*A(bellpoly(r, k, list(V[1:r - k + 2]))) for k in range(1, r + 1)), arb(0))

def ab_exact(M0q, M1q, rmax):
    """A_bounds_s for the box R.box(0,0,0,M0,M1): V = [M0, M1, 0, 0, ...] (M2..M7 = 0)"""
    V = [Q_(M0q), Q_(M1q)] + [Q_(0)]*6
    Fa = [sup_rat_arb('a', k, M0q) for k in range(rmax + 3)]; Fb = [sup_rat_arb('b', k, M0q) for k in range(rmax + 3)]
    pp, Um = PQ, UQ; Ab = {}
    for r in range(rmax + 1):
        Ab[(2, r)] = faa_num(Fa, r, V)
        Ab[(1, r)] = Ab[(2, r)] + sum((A(math.comb(r, q))*faa_num(Fb, q, V)*A(V[r - q + 1] + pp**(r - q + 1)*Um) for q in range(r + 1)), arb(0))
        Ab[(0, r)] = sum((A(math.comb(r, q))*(faa_num(Fa[1:], q, V)*A((pp**2 - pp)*pp**(r - q)*Um) + faa_num(Fb[1:], q, V)*A((2*pp)**(r - q)*pp**2*Um**2)) for q in range(r + 1)), arb(0))
    return Ab

# ------------------------------------------------------------------ YC base (port of F.constants)
def YC_base(tab, j, Ab):
    n = j + 2; YC = {k: arb(0) for k in range(n + 1)}
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        gC = sum((A(co)*Ab[(i, r)] for (m, i, r), co in mons.items()), arb(0))      # gam = 1
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
    return YC

# ------------------------------------------------------------------ Aterms (symbolic once, evaluated per box)
def build_At_sym(rmax=6):
    FA = sp.symbols('fa0:%d' % (rmax + 3)); FB = sp.symbols('fb0:%d' % (rmax + 3)); m1 = sp.Symbol('m1'); Z = sp.symbols('Z0:10')
    def Vj(j): return m1 if j == 1 else PQ**j*UQ + Z[j]
    def faa(F_, r):
        if r == 0: return F_[0]
        return sum(F_[k]*bellpoly(r, k, [Vj(j) for j in range(1, r - k + 2)]) for k in range(1, r + 1))
    def collect(expr, scale=1):
        expr = sp.expand(expr*scale); out = {}
        for term in sp.Add.make_args(expr):
            c, mon = term.as_independent(*Z); key = []
            for f in sp.Mul.make_args(mon):
                if f == 1: continue
                base, ex = f.as_base_exp(); key += [int(str(base)[1:])]*int(ex)
            out.setdefault(tuple(sorted(key)), []).append(c)
        return out
    def merge(*ds):
        o = {}
        for d in ds:
            for k, v in d.items(): o.setdefault(k, []).extend(v)
        return o
    At = {}
    for r in range(rmax + 1):
        At[(2, r)] = collect(faa(FA, r))
        t1 = merge(*[collect(faa(FB, q)*(Vj(r - q + 1) + PQ**(r - q + 1)*UQ), sp.binomial(r, q)) for q in range(r + 1)])
        At[(1, r)] = merge(At[(2, r)], t1)
        At[(0, r)] = merge(*[merge(collect(faa(FA[1:], q), sp.binomial(r, q)*(PQ**2 - PQ)*PQ**(r - q)*UQ), collect(faa(FB[1:], q), sp.binomial(r, q)*(2*PQ)**(r - q)*PQ**2*UQ**2)) for q in range(r + 1)])
    return At, FA, FB, m1

def eval_At(Asym, M0q, M1q, rmax=6):
    At, FA, FB, m1 = Asym
    vals = {m1: A(M1q)}
    for k in range(rmax + 3): vals[FA[k]] = sup_rat_arb('a', k, M0q); vals[FB[k]] = sup_rat_arb('b', k, M0q)
    out = {}
    for key, d in At.items():
        o = {}
        for zk, lst in d.items():
            tot = arb(0)
            for c in lst:
                coef, rest = c.as_coeff_Mul(); v = A(abs(coef))
                for f in sp.Mul.make_args(rest):
                    if f == 1: continue
                    base, ex = f.as_base_exp(); v = v*vals[base]**int(ex)
                tot += v
            o[zk] = tot
        out[key] = o
    return out

# ------------------------------------------------------------------ classify / l4_extra_YC / fallback / ghat (arb ports)
def classify_x(T, j, M1q, At, M0zq=M0Z):
    n = j + 2; sup = {0: A(M0zq), 1: A(M1q)}; supq = {0: M0zq, 1: Q_(M1q)}; items = []
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, arb(0)) + A(co)*v          # gam = 1
        for X, w in W.items():
            lows = [x for x in (a, b) if x < 2]; R = sorted([x for x in (a, b) if x >= 2] + list(X)); Ml = arb(1)
            for l in lows: Ml = Ml*sup[l]
            if len(R) == 1:
                if len(lows) != 2: raise ValueError(('linear term', (a, b), X))
                lo = sorted(lows, key=lambda l: supq[l]); items.append(('P', lo[1], R[0], w*sup[lo[0]]))
            elif len(R) == 2: items.append(('P', R[0], R[1], w*Ml))
            elif len(R) == 3 and R.count(2) >= 2:
                r = list(R); r.remove(2); r.remove(2); items.append(('L4', r[0], w*Ml))
            elif len(R) == 3 and R.count(3) >= 2 and 2 in R: items.append(('S33', w*Ml))
            else: raise ValueError(('unhandled', (a, b), X, R))
    return items

def l4_extra_YC_x(T, j, M0q, M1q, At, dfrac=Q_(1, 10), fallback=None, s33='assumed', kap=KAP):
    n = j + 2; amin = (1 + A(M0q)**2)**-2
    add = {k: arb(0) for k in range(n + 1)}; dissip = [arb(0)]; ndp = [0]; M1 = A(M1q)
    def put(x, y, wt, d_fr):
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y))
        if y == n and x == n: dissip[0] += wt
        elif y == n:
            d = A(d_fr)*amin; add[x] += wt**2/(4*d); ndp[0] += 1
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
    try:
        items = classify_x(T, j, M1q, At)
    except ValueError:
        if fallback is None: raise
        return fallback_x(T, j, M0q, M1q, At, kap, fallback)
    for it in items:
        if it[0] == 'P': put(it[1], it[2], it[3], dfrac)
        elif it[0] == 'L4': put(3, it[1], 3*M1*it[2], dfrac); put(2, it[1], kap*M1*it[2], dfrac)
        elif s33 == 'ibp' and n >= 4:
            w = it[1]; put(3, 3, A(Q_(3, 2)*kap)*M1*w, dfrac); put(3, 4, A(Q_(3, 2))*M1*w, dfrac); put(2, 3, A(Q_(kap**2, 2))*M1*w, dfrac); put(2, 4, A(Q_(kap, 2))*M1*w, dfrac)
        else: put(3, 3, it[1]*M2STAR, dfrac)
    return add, dissip[0], amin/2, ndp[0]*dfrac, items          # dissipation used = dissip[0] + ndp*dfrac*amin

def fallback_x(T, j, M0q, M1q, At, kap, fallback):
    """port of _old_fallback_impl (all index-beyond puts use d = 0.1*amin)"""
    n = j + 2; amin = (1 + A(M0q)**2)**-2; dfr = Q_(1, 10)
    add = {k: arb(0) for k in range(n + 1)}; dissip = [arb(0)]; ndp = [0]
    sup = {0: A(M0Z), 1: A(M1q)}; supq = {0: M0Z, 1: Q_(M1q)}
    def put(x, y, wt):
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y))
        if y == n and x == n: dissipation_add(wt)
        elif y == n:
            d = A(dfr)*amin; add[x] += wt**2/(4*d); ndp[0] += 1
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
    def dissipation_add(wt): dissip[0] += wt
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, arb(0)) + A(co)*v
        for X, w in W.items():
            lows = [x for x in (a, b) if x < 2]; R = sorted([x for x in (a, b) if x >= 2] + list(X)); Ml = arb(1)
            for l in lows: Ml = Ml*sup[l]
            if len(R) == 1:
                if len(lows) != 2: raise ValueError('linear')
                lo = sorted(lows, key=lambda l: supq[l]); put(lo[1], R[0], w*sup[lo[0]]); continue
            if len(R) == 2: put(R[0], R[1], w*Ml)
            elif len(R) == 3 and R.count(2) >= 2:
                r = list(R); r.remove(2); r.remove(2); r = r[0]
                put(3, r, 3*A(M1q)*w*Ml); put(2, r, kap*A(M1q)*w*Ml)
            elif len(R) == 3 and R.count(3) >= 2 and 2 in R: put(3, 3, w*Ml*M2STAR)
            elif len(R) >= 3:
                rest = R[-2:]; ww = w*Ml
                for k in R[:-2]: ww = ww*A(Q_(fallback[k]))
                put(rest[0], rest[1], ww)
            else: raise ValueError(('unhandled', (a, b), X, R))
    return add, dissip[0], amin/2, ndp[0]*dfr, None

def budget_ok(nn, ndfr, M0q, extra=Q_(0)):
    """dissipation actually used = nn + (ndfr+extra)*amin  must be <= amin/2.  nn = sum of (n,n)-pair weights (arb).  Returns (ok, slack/amin as float)"""
    amin = (1 + A(M0q)**2)**-2; fr = Q_(ndfr) + extra
    if nn == 0: return fr <= Q_(1, 2), float(Q_(1, 2) - fr)
    sl = (Q_(1, 2) - A(fr)) - nn/amin
    return sl.lower() >= 0, float(sl.lower())

def ghat_x(YC, j, M0q, kap=KAP):
    amin = (1 + A(M0q)**2)**-2; Bc = A(Q_(3, 2))*amin; G = {k: YC[k] for k in range(j + 1)}
    G[j] = G[j] + YC[j + 1]*(YC[j + 1]/(4*Bc) + A(Q_(kap**2, 2))); return G

# ------------------------------------------------------------------ cases
def core(T, Asym, M0q, M1q, s33='assumed'):
    Ab = ab_exact(M0q, M1q, 4); At = eval_At(Asym, M0q, M1q); out = {}
    for j in (0, 1, 2):
        YCb = YC_base(T[j], j, Ab); add, dis, spare, ndfr, _ = l4_extra_YC_x(T, j, M0q, M1q, At, s33=s33)
        YC = {k: YCb[k] + add.get(k, 0) for k in YCb}; ok, sl = budget_ok(dis, ndfr, M0q); assert ok, ('core budget', j)
        out[j] = dict(YC=YC, Ghat=ghat_x(YC, j, M0q), slack=sl, YCb=YCb)
    return out

def ramp(T, Asym, M0q, M1q, rs33='ibp2', dfrac=Q_(2, 25)):
    """port of the --rampl4 block of F2.main (rs33 'ibp2'), box (M0r,M1r) = (0.05, 0.9)"""
    Ab = ab_exact(M0q, M1q, 4); At = eval_At(Asym, M0q, M1q); L = arb(2).log()
    etap = A(Q_(15, 8))/L; DEL = A(Q_(1, 10)); amin_r = (1 + A(M0q)**2)**-2; M1 = A(M1q); out = {}
    for j in (0, 1, 2):
        YCb = YC_base(T[j], j, Ab)
        add, dis, spare, ndfr, items = l4_extra_YC_x(T, j, M0q, M1q, At, dfrac=dfrac, s33=rs33)
        fw = arb(0); extra = Q_(0); w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0)); has33 = any(it[0] == 'S33' for it in items)
        if rs33 in ('ibp', 'ibp2') and has33:
            c1 = w33*M1*etap; kap = 5
            add[3] = add.get(3, arb(0)) + c1*kap*DEL; fw += c1*kap/DEL
            dd = A(Q_(1, 10))*amin_r; extra = Q_(1, 10); fw += c1**2/(4*dd)
        for it in items:
            if it[0] != 'L4': continue
            w = it[2]
            if it[1] == 3: add[3] = add.get(3, arb(0)) + w*M1*etap*DEL; fw += w*M1*etap/DEL
            elif it[1] == 2: add[2] = add.get(2, arb(0)) + w*M1*etap; fw += w*M1*etap
            else: raise ValueError(('foot r', it))
        ok, sl = budget_ok(dis, ndfr, M0q, extra); assert ok, ('ramp budget', j)
        YC = {k: YCb[k] + add.get(k, arb(0)) for k in YCb}
        out[j] = dict(YC=YC, Ghat=ghat_x(YC, j, M0q), footw=fw, W33=(w33 if rs33 == 'ibp2' else arb(0)), slack=sl, has33=has33, YCb=YCb)
    return out

def window(T, Asym, M0q, M1q, f2=Q_(1, 2)):
    """per-window block of lemmaS_Bbox.py, levels 0..3; returns Ghat_jj / Ghat rows and YC3"""
    Ab = ab_exact(M0q, M1q, 4); At = eval_At(Asym, M0q, M1q); Gk = {}
    for j in range(3):
        YCb = YC_base(T[j], j, Ab); add, dis, spare, nd, _ = l4_extra_YC_x(T, j, M0q, M1q, At, s33='ibp'); assert budget_ok(dis, nd, M0q)[0], ('window budget', j)
        YC = {k: YCb[k] + add.get(k, 0) for k in YCb}; Gk[j] = ghat_x(YC, j, M0q)
    tab3 = T[3]; YCb = YC_base(tab3, 3, Ab)       # same loop as in Bbox (range(6) keys 0..5)
    add, dis, spare, nd, _ = l4_extra_YC_x({3: tab3}, 3, M0q, M1q, At, fallback={2: f2, 3: Q_(0)}, s33='ibp')
    assert budget_ok(dis, nd, M0q)[0], 'level 3 dissipation'
    YC3 = {k: YCb[k] + add.get(k, 0) for k in YCb}; Gk[3] = ghat_x(YC3, 3, M0q)
    return Gk, YC3

def up_rational(x, den=10**15):
    return Q_(int((x.upper()*den).ceil().unique_fmpz()), den)

def window_boxes(al=Q_(15, 100000), dv=Q_(15, 100)):
    out = []
    for yE in (Q_(50), Q_(75), Q_(225, 2), Q_(675, 4), Q_(2025, 8)):
        yb = A(Q_(3, 2)*yE); Uyb = A(Q_(13056, 10000))*A(al)*(yb.log()*A(PQ)).exp()
        out.append((yE, up_rational(Uyb + A(dv)), up_rational(A(Q_(105, 100))*Uyb + A(dv))))
    return out

if __name__ == '__main__':
    T = build_tables(); Asym = build_At_sym(); res = {}; cmp = {}
    def rel(new, old):          # (arb upper - old)/old
        return float(((new.upper() - arb(old))/arb(old)).upper())
    # ---- table coefficients: exact vs the double tables of PB.table
    Tf = {j: PB.table(F.div_form(j)[0]) for j in (0, 1, 2)}; md = 0.0
    for j in (0, 1, 2):
        for ab, mons in Tf[j].items():
            for k, v in mons.items(): md = max(md, abs(float(T[j][ab][k]) - v)/v)
    cmp['table_coef_max_rel_diff'] = md; print('table coefficients: max rel diff exact vs float', md)
    # ---- core
    c = core(T, Asym, Q_(11, 100), Q_(52, 100)); d = json.load(open('lemmaS_F10p_M0.52.json'))['Ghat_core']; res['core'] = {}; cmp['core'] = {}
    for j in (0, 1, 2):
        res['core'][j] = {k: v for k, v in c[j]['Ghat'].items()}; cmp['core'][j] = {}
        for k, v in c[j]['Ghat'].items():
            old = float(d[f'Ghat_{j}'][str(k)]); cmp['core'][j][k] = (old, float(v.upper()), rel(v, old))
        print('core level', j, {k: (round(o, 6), '%.3e' % r) for k, (o, u, r) in cmp['core'][j].items()}, 'budget slack/amin', c[j]['slack'])
    # ---- ramp
    r_ = ramp(T, Asym, Q_(5, 100), Q_(9, 10)); cap = json.load(open(os.path.join(HERE, 'cap_F10p.json')))['ramp_level']
    res['ramp'] = {}; cmp['ramp'] = {}; ramp_excess = 0.0
    for j in (0, 1, 2):
        res['ramp'][j] = dict(Ghat=r_[j]['Ghat'], YC=r_[j]['YC'], footw=r_[j]['footw'], W33=r_[j]['W33']); cmp['ramp'][j] = {}
        Ghr = cap[j]['Ghr']; YCr = cap[j]['YCr']
        for k, v in r_[j]['Ghat'].items():
            old = float(Ghr[str(k)]); rr = rel(v, old); cmp['ramp'][j][f'Ghat{k}'] = (old, float(v.upper()), rr); ramp_excess = max(ramp_excess, rr)
        for k, v in r_[j]['YC'].items():
            if str(k) in YCr: old = float(YCr[str(k)]); rr = rel(v, old) if old else 0.0; cmp['ramp'][j][f'YC{k}'] = (old, float(v.upper()), rr); ramp_excess = max(ramp_excess, rr)
        print('ramp level', j, {k: (round(o, 4), '%.3e' % r) for k, (o, u, r) in cmp['ramp'][j].items() if k.startswith('Ghat')}, 'footw', float(r_[j]['footw'].upper()), 'W33', float(r_[j]['W33'].upper()), 'budget slack/amin', r_[j]['slack'])
    # Abr used in ramp_level (doubles) vs exact
    Abx = ab_exact(Q_(5, 100), Q_(9, 10), 6); Abr = cap[0]['Abr']; mx = 0.0
    for key, v in Abx.items():
        old = float(Abr[f'({key[0]}, {key[1]})']);
        if old > 0: mx = max(mx, rel(v, old))
    cmp['ramp_Abr_max_rel_excess'] = mx; cmp['ramp_Ghat_YC_max_rel_excess'] = ramp_excess
    print('ramp: max relative excess (exact arb upper / double - 1): Abr %.3e, Ghat/YC %.3e' % (mx, ramp_excess))
    # C_eta (arb, exact ln 2)
    L = arb(2).log(); N = 4000; best = arb(0)
    for i in range(N):
        t = (arb(i)/N).union(arb(i + 1)/N); S = t**3*(10 - 15*t + 6*t**2); S1 = 30*t**2 - 60*t**3 + 30*t**4; S2 = 60*t - 180*t**2 + 120*t**3
        v = abs(-5*S*S1/L + S1**2/(L*L) + S*S2/(L*L)).upper()
        if v > best: best = v
    res['C_eta'] = best; print('C_eta arb upper', float(best), '(double run 15.83531249835125)')
    # ---- windows
    bb = json.load(open('lemmaS_Bbox_1e-05_a0.00015_dv0.15.json')); res['windows'] = {}; cmp['windows'] = {}
    for (yE, M0q, M1q), row in zip(window_boxes(), bb):
        Gk, YC3 = window(T, Asym, M0q, M1q); key = str(yE)
        res['windows'][key] = dict(M0=float(M0q), M1=float(M1q), M0q=str(M0q), M1q=str(M1q), Ghat={j: Gk[j] for j in range(4)}, YC3=YC3)
        cmp['windows'][key] = [(row['Ghat'][j], float(Gk[j][j].upper()), rel(Gk[j][j], row['Ghat'][j])) for j in range(4)]
        print('window', key, 'box', float(M0q), float(M1q), 'Ghat_jj exact-arb-upper vs run:', ['%.6g/%.6g (%.1e)' % (o, u, r) for o, u, r in cmp['windows'][key]])
    json.dump(dict(res=arbio.dump(res), compare=cmp), open(os.path.join(HERE, 'weights_exact.json'), 'w'), default=str)
    worst = max([abs(v[2]) for c_ in cmp['core'].values() for v in c_.values()] + [ramp_excess, mx] + [abs(v[2]) for w in cmp['windows'].values() for v in w])
    print('WORST relative deviation of the exact arb upper value above the double used: %.3e' % worst)
