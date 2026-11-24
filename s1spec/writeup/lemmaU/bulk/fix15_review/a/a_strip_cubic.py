"""fix15 review, part A: independent recomputation of the strip cubic coefficients C3S_abc of fin13v4c (D-3).
Differences from v4bcert/strip_cubic_v4c.py: NS = 300 cells per unit sigma (theirs 200); x-cells built in arb from the sigma cell (no float midpoint/radius);
my own arb evaluation through sympy.lambdify-free code generation of the second derivatives (same printer idea, re-derived d^2C from the pickle);
plus (b) a pointwise float sanity test of the monomial bound and (c) the H8 Hoelder assembly constants.
[cert] = arb, outward; inputs (tc_weighted.pkl, profile hulls, jets) taken as given.  Writes only a_strip_cubic.json next to this file."""
import sys, os, json, pickle, math, hashlib, random, time
sys.dont_write_bytecode = True
from flint import arb, ctx
import sympy as sp
ctx.prec = 140
t0 = time.time()
V = '/home/claude/s1spec/step3nl_v4c'
PKL = V + '/v4bcert/out/tc_weighted.pkl'; SCJ = V + '/v4bcert/out/strip_cubic_v4c.json'; E0J = V + '/v4bcert/out/e0v4b_constants.json'
C = pickle.load(open(PKL, 'rb'))['C']
U = sp.symbols('U0 U1 U2'); names = "phi1, phi2, k, m, p, X, X1, X2, U0, U1, U2"
def gen(e):
    # arb code: rationals exact, sqrt via .sqrt(), integer powers by repeated multiplication
    def pr(x):
        if x.is_Symbol: return str(x)
        if x.is_Integer: return f"arb({int(x)})"
        if x.is_Rational: return f"(arb({x.p})/arb({x.q}))"
        if x.is_Add: return "(" + "+".join(pr(a) for a in x.args) + ")"
        if x.is_Mul: return "(" + "*".join(pr(a) for a in x.args) + ")"
        if x.is_Pow:
            b, ex = x.as_base_exp(); bb = pr(b)
            if ex.is_Integer:
                n = int(ex); s = "*".join([bb] * abs(n)); return f"({s})" if n > 0 else f"(arb(1)/({s}))"
            if ex.is_Rational and ex.q == 2:
                n = int(ex.p); s = "*".join([f"({bb}).sqrt()"] * abs(n)); return f"({s})" if n > 0 else f"(arb(1)/({s}))"
        raise ValueError(x)
    return eval(f"lambda {names}: " + pr(e), {"arb": arb})
D2 = {}; D2f = {}
for l in range(3):
    for i in range(3):
        for j in range(i, 3):
            e = sp.diff(C[l], U[i], U[j])
            if e != 0:
                D2[(l, i, j)] = gen(e); D2f[(l, i, j)] = sp.lambdify(sp.symbols(names.replace(' ', '')), e, 'math')
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
AM = arb('2.5e-4'); Ab = arb(0, AM); SB = arb(0, arb('0.02'))
def Sm(u): return u**4 * (35 - 84*u + 70*u*u - 20*u**3)
def Sm1(u): return 140 * u**3 * (1 - u)**3
def Sm2(u): return 420 * u*u * (1 - u)**2 * (1 - 2*u)
NS = 300; mono = {}
for ci in range(29, 32):
    for i in range(NS):
        sa = arb(ci) + arb(i) / NS; sb = arb(ci) + arb(i + 1) / NS
        kb, mb, pb, qs = HR.HS.profile_bounds(sa)
        Xv = (sa - HR.ST).union(sb - HR.ST)                      # x-cell, arb hull
        xc = float(Xv.mid())
        E = (1 + Xv * Xv) ** (arb(-1) / 3); ph1 = Xv * E / 3; ph2 = E / 3 - 2 * Xv * Xv * E**4 / 9
        y = HR.jets_any(HR.F, Xv, xc)[:3]; h = HR.jets_any(HR.Hf, Xv, xc)[:3]
        T = ((32 - sb) / 3).union((32 - sa) / 3)
        ch, ch1, ch2 = Sm(T), -Sm1(T) / 3, Sm2(T) / 9
        Xw, X1w, X2w = ch * ch, 2 * ch * ch1, 2 * ch1 * ch1 + 2 * ch * ch2
        Ub = [Ab * y[j] + Ab * Ab * h[j] + SB for j in range(3)]
        kk, mm, pp = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())
        for key, f in D2.items():
            v = abs(f(ph1, ph2, kk, mm, pp, Xw, X1w, X2w, *Ub)); assert v.is_finite()
            cnt = [0, 0, 0]
            for o in key: cnt[o] += 1
            mult = 1 if key[1] == key[2] else 2
            val = arb(v.upper()) * mult / 2
            kk2 = (tuple(cnt), key)
            if kk2 not in mono or val.upper() > mono[kk2].upper(): mono[kk2] = val
tot = {}
for (cnt, key), v in mono.items(): tot[cnt] = tot.get(cnt, arb(0)) + v
sc = json.load(open(SCJ)); coded = {k: v['value'] for k, v in sc['constants'].items()}
res = {}; fails = []
for cnt in sorted(tot):
    n = 'C3S_%d%d%d' % cnt; mine = float(tot[cnt].upper()); ok = bool(arb(coded[n]) >= tot[cnt].upper())
    res[n] = dict(mine_NS300=mine, coded=coded[n], coded_over_mine=coded[n] / mine if mine else None, covered=ok)
    if not ok: fails.append(n)
res['keys_equal'] = sorted('C3S_%d%d%d' % c for c in tot) == sorted(k for k, v in coded.items() if v > 0) or [sorted('C3S_%d%d%d' % c for c in tot), sorted(coded)]
# (b) pointwise sanity [num]: random sigma, random jets in the box, exact t-integral by Gauss-Legendre (float)
random.seed(1); worst = 0.0
gl = [(0.5 - 0.5 * 0.7745966692414834, 5 / 18), (0.5, 8 / 18), (0.5 + 0.5 * 0.7745966692414834, 5 / 18)]
for trial in range(400):
    s = 29 + 3 * random.random(); sa = arb(s)
    kb, mb, pb, qs = HR.HS.profile_bounds(sa)
    x = float((sa - HR.ST).mid()); Xv = arb(x)
    Ef = (1 + x * x) ** (-1 / 3); ph1 = x * Ef / 3; ph2 = Ef / 3 - 2 * x * x * Ef**4 / 9
    y = [float(z.mid()) for z in HR.jets_any(HR.F, Xv, x)[:3]]; h = [float(z.mid()) for z in HR.jets_any(HR.Hf, Xv, x)[:3]]
    kf, mf, pf = [random.uniform(-1, 1) * float(b.upper()) for b in (kb, mb, pb)]
    t = (32 - s) / 3; ch, ch1, ch2 = Sm(t), -Sm1(t) / 3, Sm2(t) / 9
    Xw, X1w, X2w = ch * ch, 2 * ch * ch1, 2 * ch1 * ch1 + 2 * ch * ch2
    al = random.uniform(0, 2.5e-4); w = [random.uniform(-0.02, 0.02) for _ in range(3)]
    # 3-point GL is not exact; use 3 points on (1-t) weight as an approximation and a 200-point midpoint rule as reference
    R3 = 0.0; NT = 60
    for q in range(NT):
        tt = (q + 0.5) / NT; Uv = [al * y[j] + al * al * h[j] + tt * w[j] for j in range(3)]
        sm = 0.0
        for (l, i_, j_), f in D2f.items():
            sm += (1 if i_ == j_ else 2) * f(ph1, ph2, kf, mf, pf, Xw, X1w, X2w, *Uv) * w[i_] * w[j_] * w[l]
        R3 += (1 - tt) * sm / NT
    bnd = sum(coded['C3S_%d%d%d' % (a, b, c)] * abs(w[0])**a * abs(w[1])**b * abs(w[2])**c for (a, b, c) in ((0,1,2),(0,2,1),(0,3,0),(1,0,2),(1,1,1),(1,2,0),(2,0,1),(2,1,0),(3,0,0)))
    worst = max(worst, abs(R3) / bnd)
res['pointwise_sanity_max_ratio_400_samples'] = worst
if not worst <= 1.0: fails.append('pointwise sanity')
# (c) Hoelder constants read by the verifier
E0 = json.load(open(E0J))['constants']
res['strip_integrals'] = {n: E0[n]['value'] for n in ('IR', 'IRW', 'RW')}
res['xcheck_in_stage'] = sc['strip_integrals_crosscheck']
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
json.dump(dict(res=res, fails=fails, args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, prec=ctx.prec, NS=NS,
          inputs={p: sha(p) for p in (PKL, SCJ, E0J, os.path.abspath(__file__))}, runtime_s=round(time.time() - t0, 1))),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a_strip_cubic.json'), 'w'), indent=1)
for k, v in res.items(): print(k, v)
print('FAILS', fails)
