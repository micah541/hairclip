#!/usr/bin/env python3
"""P1 (Lemma U, phase 1): CERTIFIED (arb, outward) rate constants hat-G_jj (j = 0,1,2) for the Lemma-U weight with the general-weight L^4 classes,
the certified damping fixed point, and the sanity reproduction of the power-weight values of lemmaS/rigor/weights_exact.json.

Port (copied, NOT imported: weights_exact.py chdirs at import) of lemmaS/rigor/weights_exact.py {ab_exact, sup_rat_arb, faa_num, YC_base, build_At_sym,
eval_At, classify_x, l4_extra_YC_x, budget_ok, ghat_x, core}, with exactly these changes (L4_general_weight.md sec. 2-3):
  * weight-derivative factor gamma^m (there gam = 1, i.e. omitted)  ->  Gamma_m (certified arb upper bounds, weight_consts_arb.certify); Gamma_0 = 1.
    Tables store |coef| per monomial (m,i,r) (exact Rationals, checked >= 0 below), so |c om^(m) A| <= |c| Gamma_m Theta |A| term by term:
    no sign cancellation between different m (or different monomials) is used.  The coefficient of one monomial is a weight-independent
    rational (om is a generic sympy Function in lemmaS_prop1_sym), so its exact value is valid for every weight.
  * kappa in the L^4 classes ('L4': 3 M1 N3 + kappa M1 N2; 'S33' ibp: 3/2 kappa, 3/2, kappa^2/2, kappa/2)  ->  kappa~ = sup|(ln eps Theta)'| (arb).
    (kappa^2 there is a product of two sup|L| factors, so kappa~^2, not kappa2~.)
  * interpolation in ghat: kappa^2/2  ->  kappa2~/2, kappa2~ = sup((eps Theta)''/(eps Theta))_+  (L4_general_weight.md (d)).
Everything else (M0z = 1/50, p = 1.0246, Umax = 1/20, dfrac = 1/10, Young splitting, budget check) is unchanged; inherited hypotheses, see P1_REPORT.md.
Run:  cd <this dir> && PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 nice -n 19 python3 ghat_arb.py      (single process, < 200 MB, ~1-3 min)"""
import os, sys, json, pickle, math, hashlib, time, platform
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from sympy.functions.combinatorial.numbers import bell as bellpoly
from flint import arb, ctx
import weight_consts_arb as WC
ctx.prec = 200
Q_ = sp.Rational
LS = '/home/claude/s1spec/writeup/lemmaS'
TABLES = os.path.join(LS, 'rigor/tables_exact.pkl'); WEXACT = os.path.join(LS, 'rigor/weights_exact.json')
AARB = '/home/claude/s1spec/writeup/step3_gap_v32/gonogo/lemmaS_A_rerun/rigor/assembly_arb.json'
T8 = '/home/claude/s1spec/writeup/step3_gap_v32/gonogo/t8_L4_cost.json'
PQ = Q_(10246, 10000); UQ = Q_(1, 20); M0Z = Q_(1, 50); M2STAR = 3

def A(q):
    if isinstance(q, arb): return q
    q = Q_(q); return arb(int(q.p))/arb(int(q.q))

# ------------------------------------------------------------------ Ab (copy of weights_exact.ab_exact & helpers)
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
        u = abs(f(lo.union(hi))).upper()
        if u > best: best = u
    _supc[key] = best; return best

def faa_num(F_, r, V):
    if r == 0: return F_[0]
    return sum((F_[k]*A(bellpoly(r, k, list(V[1:r - k + 2]))) for k in range(1, r + 1)), arb(0))

def ab_exact(M0q, M1q, rmax, Umq=UQ):
    V = [Q_(M0q), Q_(M1q)] + [Q_(0)]*6
    Fa = [sup_rat_arb('a', k, M0q) for k in range(rmax + 3)]; Fb = [sup_rat_arb('b', k, M0q) for k in range(rmax + 3)]
    pp, Um = PQ, Umq; Ab = {}
    for r in range(rmax + 1):
        Ab[(2, r)] = faa_num(Fa, r, V)
        Ab[(1, r)] = Ab[(2, r)] + sum((A(math.comb(r, q))*faa_num(Fb, q, V)*A(V[r - q + 1] + pp**(r - q + 1)*Um) for q in range(r + 1)), arb(0))
        Ab[(0, r)] = sum((A(math.comb(r, q))*(faa_num(Fa[1:], q, V)*A((pp**2 - pp)*pp**(r - q)*Um) + faa_num(Fb[1:], q, V)*A((2*pp)**(r - q)*pp**2*Um**2)) for q in range(r + 1)), arb(0))
    return Ab

def YC_base(tab, j, Ab, Gam):
    n = j + 2; YC = {k: arb(0) for k in range(n + 1)}
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        gC = sum((A(co)*Gam[m]*Ab[(i, r)] for (m, i, r), co in mons.items()), arb(0))      # gamma^m -> Gamma_m
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
    return YC

_ASYM = {}
def build_At_sym(Umq=UQ, rmax=6):
    if Umq in _ASYM: return _ASYM[Umq]
    FA = sp.symbols('fa0:%d' % (rmax + 3)); FB = sp.symbols('fb0:%d' % (rmax + 3)); m1 = sp.Symbol('m1'); Z = sp.symbols('Z0:10')
    def Vj(j): return m1 if j == 1 else PQ**j*Umq + Z[j]
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
        t1 = merge(*[collect(faa(FB, q)*(Vj(r - q + 1) + PQ**(r - q + 1)*Umq), sp.binomial(r, q)) for q in range(r + 1)])
        At[(1, r)] = merge(At[(2, r)], t1)
        At[(0, r)] = merge(*[merge(collect(faa(FA[1:], q), sp.binomial(r, q)*(PQ**2 - PQ)*PQ**(r - q)*Umq), collect(faa(FB[1:], q), sp.binomial(r, q)*(2*PQ)**(r - q)*PQ**2*Umq**2)) for q in range(r + 1)])
    _ASYM[Umq] = (At, FA, FB, m1); return _ASYM[Umq]

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

def classify_x(T, j, M1q, At, Gam, M0zq=M0Z):
    n = j + 2; sup = {0: A(M0zq), 1: A(M1q)}; supq = {0: M0zq, 1: Q_(M1q)}; items = []
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, arb(0)) + A(co)*Gam[m]*v          # gamma^m -> Gamma_m
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

def l4_extra_YC_x(T, j, M0q, M1q, At, Gam, kap, dfrac=Q_(1, 10), s33='ibp', M0zq=M0Z):
    n = j + 2; amin = (1 + A(M0q)**2)**-2; kap = A(kap)
    add = {k: arb(0) for k in range(n + 1)}; dissip = [arb(0)]; ndp = [0]; M1 = A(M1q)
    def put(x, y, wt, d_fr):
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y))
        if y == n and x == n: dissip[0] += wt
        elif y == n:
            d = A(d_fr)*amin; add[x] += wt**2/(4*d); ndp[0] += 1
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
    items = classify_x(T, j, M1q, At, Gam, M0zq)
    for it in items:
        if it[0] == 'P': put(it[1], it[2], it[3], dfrac)
        elif it[0] == 'L4': put(3, it[1], 3*M1*it[2], dfrac); put(2, it[1], kap*M1*it[2], dfrac)        # kappa -> kappa~
        elif s33 == 'ibp' and n >= 4:
            w = it[1]; put(3, 3, A(Q_(3, 2))*kap*M1*w, dfrac); put(3, 4, A(Q_(3, 2))*M1*w, dfrac); put(2, 3, kap*kap/2*M1*w, dfrac); put(2, 4, kap/2*M1*w, dfrac)
        else: put(3, 3, it[1]*M2STAR, dfrac)                                                           # 'assumed' sup|z2| <= 3 (reference only)
    return add, dissip[0], amin/2, ndp[0]*dfrac, items

def budget_ok(nn, ndfr, M0q, extra=Q_(0)):
    amin = (1 + A(M0q)**2)**-2; fr = Q_(ndfr) + extra
    if nn == 0: return fr <= Q_(1, 2), float(Q_(1, 2) - fr)
    sl = (Q_(1, 2) - A(fr)) - nn/amin
    return sl.lower() >= 0, float(sl.lower())

def ghat_x(YC, j, M0q, kap2):
    amin = (1 + A(M0q)**2)**-2; Bc = A(Q_(3, 2))*amin; G = {k: YC[k] for k in range(j + 1)}
    G[j] = G[j] + YC[j + 1]*(YC[j + 1]/(4*Bc) + A(kap2)/2); return G                                   # kappa^2 -> kappa2~

def core(T, M0q, M1q, Gam, kap, kap2, s33='ibp', Umq=UQ, M0zq=M0Z):
    Asym = build_At_sym(Umq); Ab = ab_exact(M0q, M1q, 4, Umq); At = eval_At(Asym, M0q, M1q); out = {}
    for j in (0, 1, 2):
        YCb = YC_base(T[j], j, Ab, Gam); add, dis, spare, ndfr, items = l4_extra_YC_x(T, j, M0q, M1q, At, Gam, kap, s33=s33, M0zq=M0zq)
        YC = {k: YCb[k] + add.get(k, 0) for k in YCb}; ok, sl = budget_ok(dis, ndfr, M0q)
        GL = ghat_x(YC, j, M0q, kap2); G0 = ghat_x(YCb, j, M0q, kap2)
        out[j] = dict(Ghat=GL, Ghat_noL4=G0, budget_ok=bool(ok), slack=sl, n_items=len(items))
    return out

# ------------------------------------------------------------------ helpers
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def fup(x):
    x = x.upper() if isinstance(x, arb) else arb(x); v = float(x)
    return v if arb(v) >= x else math.nextafter(v, math.inf)
def flo(x):
    x = x.lower() if isinstance(x, arb) else arb(x); v = float(x)
    return v if arb(v) <= x else math.nextafter(v, -math.inf)
def up_rat(x, den=10**12):
    return Q_(int((x.upper()*den).ceil().unique_fmpz()), den)
def rows(G):
    return [[flo(G[j]['Ghat'][j]), fup(G[j]['Ghat'][j])] for j in (0, 1, 2)]

if __name__ == '__main__':
    t0 = time.time(); out = dict(sanity={}, consts={}, cases={}, notes=[])
    T = pickle.load(open(TABLES, 'rb'))
    for j in (0, 1, 2):
        for mons in T[j].values():
            for (m, i, r), co in mons.items(): assert Q_(co) >= 0 and m <= 3
    print('tables: |coef| exact Rationals, all >= 0, weight-derivative orders m <= 3', flush=True)
    # ---------------- 4. sanity: power weight Gamma_m = 1, kappa = 5, kappa2 = 25
    one = {m: arb(1) for m in range(4)}
    we = json.load(open(WEXACT))['res']['core']
    import arbio_local as AL
    for s33 in ('assumed', 'ibp'):
        c = core(T, Q_(11, 100), Q_(52, 100), one, 5, 25, s33=s33)
        out['sanity'][s33] = dict(Ghat_jj=rows(c), Ghat_noL4_jj=[fup(c[j]['Ghat_noL4'][j]) for j in range(3)], budget_ok=[c[j]['budget_ok'] for j in range(3)], slack=[c[j]['slack'] for j in range(3)])
        print(f'power weight (gamma=1, kappa=5), crude box, s33={s33}: Ghat_jj', rows(c), 'budget', [c[j]['budget_ok'] for j in range(3)], flush=True)
        if s33 == 'assumed':
            ref = {j: AL.load(we[str(j)]) for j in range(3)}; dev = []
            for j in range(3):
                for k in range(j + 1):
                    r_ = ref[j][str(k)]; dev.append(float(abs(c[j]['Ghat'][k] - r_).upper()/r_.lower()))
            out['sanity']['max_rel_dev_vs_weights_exact'] = max(dev)
            print('   vs lemmaS/rigor/weights_exact.json core Ghat (all k <= j): max rel deviation %.3e' % max(dev), flush=True)
            assert max(dev) < 1e-12, 'sanity reproduction FAILED'
    t8 = json.load(open(T8)); out['sanity']['t8_float_power_ref_ibp'] = t8['refs']['power weight, Lemma S(A) core (gamma=1, kappa=5)']['Ghat_L4']
    # ---------------- 1. weight constants: certify once at theta_a = 3/5 (bounds are monotone in theta_a -> valid for all theta_a <= 3/5)
    WCONST = {}
    for (mu, nu, alo, ahi) in ((4, 4, Q_(0), Q_(3, 10)), (4, 2, Q_(1, 2), Q_(1, 2))):
        for yL in (56, 150):
            r = WC.certify(mu, nu, alo, ahi, Q_(3, 5), yL)
            assert not r['unresolved_cells']
            WCONST[(mu, nu, yL)] = r
            out['consts'][f'mu={mu},nu={nu},a=[{alo},{ahi}],yL={yL}'] = dict(theta_a_max=str(Q_(3, 5)), Gamma_m_up=[r['bound_float_up'][k] for k in ('G1', 'G2', 'G3')],
                kappa_up=r['bound_float_up']['kap'], kappa2_up=r['bound_float_up']['kap2'], limits=r['limits'], leaves=r['leaves'], geometry=r['geometry'])
            print(f'consts mu={mu} nu={nu} a=[{alo},{ahi}] yL={yL} (all theta_a <= 3/5):', {k: '%.15g' % v for k, v in r['bound_float_up'].items()}, flush=True)
    # ---------------- 2./3. Ghat and fixed point
    def ghat_case(mu, nu, yL, M0q, M1q, Umq=UQ, M0zq=M0Z):
        b = WCONST[(mu, nu, yL)]['bound']; Gam = {0: arb(1), 1: b['G1'], 2: b['G2'], 3: b['G3']}
        return core(T, M0q, M1q, Gam, b['kap'], b['kap2'], s33='ibp', Umq=Umq, M0zq=M0zq)
    boxes = dict(crude=(Q_(11, 100), Q_(52, 100), 56), fine=(Q_(21, 100), Q_(21, 100), 150))
    # (A)'s certified core box on [28,60] (tightened (H2) rerun): read the arb values, round up to 1e-7
    AA = json.load(open(AARB)); Abox = {}
    for row in AA['claims']:
        nm = str(row[0]) if isinstance(row, list) else str(row.get('name', ''))
        if nm.startswith('core M0 (page') or nm.startswith('core M1 (page'): Abox[nm[5:7]] = row
    out['A_core_box_rows'] = Abox
    boxes['A_core'] = (Q_(167137, 10**7), Q_(1070809, 10**7), 56)
    assert Q_(Abox['M0'][3]) <= boxes['A_core'][0] and Q_(Abox['M1'][3]) <= boxes['A_core'][1], 'A core box rounding'
    for (mu, nu) in ((4, 4), (4, 2)):
        for name, (M0q, M1q, yL) in boxes.items():
            c = ghat_case(mu, nu, yL, M0q, M1q); G22 = c[2]['Ghat'][2]
            assert all(c[j]['budget_ok'] for j in range(3)), ('budget', mu, nu, name)
            # certified fixed point: theta_a := up(G22/28^4) (constants do not depend on theta_a <= 3/5); check G22 <= 28^4 theta_a
            th = up_rat(G22/A(28**4), 10**9); assert th <= Q_(3, 5) and A(28**4)*A(th) >= G22.upper()
            g = th*28**4/(1 - th); Gw = 1/(1 - A(th)); rm = A(28**4)/G22
            out['cases'][f'{name}|mu={mu},nu={nu}'] = dict(box=[str(M0q), str(M1q)], yL=yL, Ghat_jj=rows(c), Ghat_noL4_jj=[fup(c[j]['Ghat_noL4'][j]) for j in range(3)],
                Ghat_full={j: {k: [flo(v), fup(v)] for k, v in c[j]['Ghat'].items()} for j in range(3)},
                budget_slack=[c[j]['slack'] for j in range(3)], L4_factor_22=[flo(G22/c[2]['Ghat_noL4'][2]), fup(G22/c[2]['Ghat_noL4'][2])],
                theta_a=str(th), theta_a_float=float(th), g=float(g), Gw=[flo(Gw), fup(Gw)], rate_margin_28=[flo(rm), fup(rm)], rate_margin_28_lower=flo(rm))
            print(f'{name:7s} mu={mu} nu={nu} box ({float(M0q)}, {float(M1q)}): Ghat_jj', ['%.7g' % x[1] for x in rows(c)], 'noL4', ['%.6g' % fup(c[j]['Ghat_noL4'][j]) for j in range(3)],
                  f'| theta_a = {float(th):.9f}, Gw <= {fup(Gw):.6f}, 28^4/Ghat22 >= {flo(rm):.5f}', 'slack', [round(c[j]['slack'], 3) for j in range(3)], flush=True)
    # ---------------- radius-dependent (hybrid) option, CONDITIONAL (see report): small box on [28,60], crude box beyond 60
    for (mu, nu) in ((4, 4), (4, 2)):
        Gs = A(fup(A(out['cases'][f'A_core|mu={mu},nu={nu}']['Ghat_jj'][2][1]))); Gc = A(out['cases'][f'crude|mu={mu},nu={nu}']['Ghat_jj'][2][1])
        # need g/(1+g e^{-4s}) >= Ghat(box at s) for all s >= ln28:  g >= Gs/(1 - Gs/28^4)  and  g >= Gc/(1 - Gc/60^4)
        g1 = Gs/(1 - Gs/A(28**4)); g2 = Gc/(1 - Gc/A(60**4)); gq = g1 if g1 > g2 else g2
        th = gq/(A(28**4) + gq); Gw = 1/(1 - th)
        out['cases'][f'hybrid_CONDITIONAL|mu={mu},nu={nu}'] = dict(g_from_small=fup(g1), g_from_crude_beyond_60=fup(g2), binding='crude beyond 60' if g2 > g1 else 'small box at 28',
            theta_a_up=fup(th), Gw_up=fup(Gw), equivalent_rate_margin=flo(1/th), status='CONDITIONAL: needs a localised (position-dependent box) Ghat derivation, not available; see P1_REPORT.md')
        print(f'hybrid (CONDITIONAL) mu={mu} nu={nu}: g >= max({fup(g1):.5g} [small @28], {fup(g2):.5g} [crude @60]) -> theta_a <= {fup(th):.6f}, Gw <= {fup(Gw):.5f}', flush=True)
    # ---------------- sensitivity: Umax (inherited 1/20) -> 9/100 (covers U(380 lambda1) = 0.086 at alpha1 = 1.5e-4)
    c = ghat_case(4, 4, 56, Q_(11, 100), Q_(52, 100), Umq=Q_(9, 100)); G22 = c[2]['Ghat'][2]
    out['sensitivity_Umax_0.09'] = dict(Ghat_jj=rows(c), budget_ok=[c[j]['budget_ok'] for j in range(3)], rate_margin_28_lower=flo(A(28**4)/G22), Gw_up=fup(1/(1 - A(up_rat(G22/A(28**4), 10**9)))))
    print('sensitivity Umax = 0.09 (crude, 4,4):', rows(c), 'rate margin >= %.4f' % flo(A(28**4)/G22), flush=True)
    # sensitivity: M0z (inherited sup|z_0| = sup|v| <= 1/50) -> 11/100 and 13/100 (|v| <= |f| + |U|; crude box), crude (4,4)
    for mz in (Q_(11, 100), Q_(13, 100)):
        c = ghat_case(4, 4, 56, Q_(11, 100), Q_(52, 100), M0zq=mz); G22 = c[2]['Ghat'][2]
        out[f'sensitivity_M0z_{float(mz)}'] = dict(Ghat_jj=rows(c), budget_ok=[c[j]['budget_ok'] for j in range(3)], rate_margin_28_lower=flo(A(28**4)/G22))
        print(f'sensitivity M0z = {float(mz)} (crude, 4,4):', rows(c), 'rate margin >= %.4f' % flo(A(28**4)/G22), flush=True)
    out['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), flint=__import__('flint').__version__, sympy=sp.__version__,
                       env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, prec=ctx.prec,
                       inputs={p: sha(p) for p in (TABLES, WEXACT, AARB, T8, os.path.join(HERE, 'weight_consts_arb.py'), os.path.join(HERE, 'arbio_local.py'))},
                       script_sha256=sha(os.path.abspath(__file__)), wall_s=round(time.time() - t0, 1))
    p = os.path.join(HERE, 'ghat_cert.json')
    json.dump(out, open(p + '.tmp', 'w'), indent=1, default=str); os.replace(p + '.tmp', p)
    print('wrote', p, 'in %.1f s' % (time.time() - t0))
