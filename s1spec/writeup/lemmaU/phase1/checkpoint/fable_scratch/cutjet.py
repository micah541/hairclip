# Reviewer scratch [est]: cutoff x jet cubic classes on (A)'s ramp [14,28] (box 0.05, 0.9; Umax 1/20), which lemmaS_F2.ramp_level / rigor/ramp_arb.py do not bound
# (they use the jet-free Abr).  Float midpoint grid; sizes only.
import sys, math, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/phase1/ghat')
import sympy as sp
from flint import arb, ctx
import ghat_arb as GA
ctx.prec = 64
A = GA.A; Q_ = sp.Rational
T = pickle.load(open('/home/claude/s1spec/writeup/lemmaS/rigor/tables_exact.pkl', 'rb'))
LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); gam = 1
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2s = Sx**2; E2 = sp.lambdify(t_, e2s, modules=[{}])
REM = {m: sp.lambdify((t_, Ls), sp.expand(sum(sp.binomial(m, l)*sp.diff(e2s, t_, l)/Ls**l*(-gam)**(m - l) for l in range(1, m + 1))), modules=[{}]) for m in (1, 2, 3, 4)}
def yy(t): return 14*(t*LN2).exp()
X0 = arb(22490776515)/10**10
om = lambda t: (-(arb(gam))*((yy(t)/28).log())).exp(); eps = lambda t: yy(t)**-4; Wt = lambda t: (1 + (yy(t) + X0)**2)**-2
M0r, M1r = Q_(5, 100), Q_(9, 10)
At = GA.eval_At(GA.build_At_sym(Q_(1, 20)), M0r, M1r); Abr = GA.ab_exact(M0r, M1r, 6); amin = float((1 + A(M0r)**2)**-2)
J = [6.0273e-12, 8.4673e-12, 1.8708e-12, 3.4315e-10]; STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}; BELL = {0: 1, 1: 1, 2: 2, 3: 5}
N = 1500; ts = [arb((i + 0.5)/N) for i in range(N)]
remv = {m: [abs(float(REM[m](t, LN2).mid())) for t in ts] for m in REM}
eo = [float((eps(t)*om(t)).mid()) for t in ts]; th = [float((eps(t)*om(t)*E2(t)).mid()) for t in ts]
cw = {m: [float((yy(t)**(2*m - 1)/Wt(t)).mid()) for t in ts] for m in (1, 2, 3)}; cw[0] = [float((1/(yy(t)*Wt(t))).mid()) for t in ts]
SUP = {0: 0.05, 1: 0.9}
def conv(k, coef):   # data int int coef z_k^2 via the verifier J's
    return BELL[k]*sum(STIR[k][m]*max(coef[i]*cw[m][i] for i in range(N))*J[m] for m in STIR[k])
for j in (1, 2):
    n = j + 2; cls = {}; jf = {}
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        for (m, i, r), co in mons.items():
            if m == 0: continue
            for X, v in At[(i, r)].items():
                c = float((A(co)*v).mid())
                if not X: continue
                key = tuple(sorted([a, b] + list(X))); cls.setdefault(key, {}); cls[key][m] = cls[key].get(m, 0.0) + c
    print('level', j, '(A) ramp: cutoff x jet classes (index multiset: sup over ramp of weight per eps*omega)')
    tot_data = 0.0
    for key, cm in sorted(cls.items()):
        f = [sum(c*remv[m][i] for m, c in cm.items()) for i in range(N)]
        idx = list(key); w = 1.0
        while len(idx) > 2 and any(x < 2 for x in idx):
            x = min(x for x in idx if x < 2); idx.remove(x); w *= SUP[x]
        note = ''
        if len(idx) == 2 and max(idx) <= j:
            a_, b_ = idx; val = sum(conv(k, [0.5*w*f[i]*eo[i] for i in range(N)]) for k in (a_, b_)); tot_data += val; note = 'folded to z_%d z_%d: data + %.3g' % (a_, b_, val)
        elif len(idx) == 2 and sorted(idx) == [j, j + 1]:
            lam = 50.0; coef = [(w*f[i]*eo[i])**2/(4*lam*th[i]) for i in range(N)]; val = conv(j, coef); tot_data += val
            note = 'z_%d z_%d: Young vs Theta z_%d^2 (lambda = 50 into ramp Y_%d): data + %.3g' % (j, j + 1, j + 1, j + 1, val)
        else: note = 'GENUINELY CUBIC or top-top, residual indices %s, weight x lows %.3g' % (idx, w*max(f))
        print('   %-12s sup %.4g  %s' % (key, max(f), note))
    print('   -> sum of the data additions listed: + %.3g to E_%d  [est]' % (tot_data, j))
