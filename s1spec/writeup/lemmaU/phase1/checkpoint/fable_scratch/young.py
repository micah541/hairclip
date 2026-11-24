# Reviewer scratch [est]: size of the (j+1,j+1) cutoff-derivative pair on (A)'s ramp [14,28] at level 2 when treated by
#   int c z3^2 = -int c z2 z4 + 1/2 int c'' z2^2,  |c z2 z4| <= delta*amin*Theta z4^2 + c^2/(4 delta amin Theta) z2^2,  Theta = eps*om*eta^2,
# i.e. data only in z2 (J_1, J_2 of the verifier).  Float grid (midpoints), not certified.
import sys, json, math, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/phase1/ghat')
import sympy as sp
from flint import arb, ctx
import ghat_arb as GA
ctx.prec = 64
A = GA.A; Q_ = sp.Rational
T = pickle.load(open('/home/claude/s1spec/writeup/lemmaS/rigor/tables_exact.pkl', 'rb'))
LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); gam = 1
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2s = Sx**2
E2 = sp.lambdify(t_, e2s, modules=[{}])
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
Abr = GA.ab_exact(Q_(5, 100), Q_(9, 10), 6); amin = float((1 + A(Q_(5, 100))**2)**-2)
def cd(mons, d, t):
    tot = arb(0)
    for (m, i, r), co in mons.items():
        if m == 0: continue
        for p_ in range(d + 1):
            for q_ in range(d + 1 - p_):
                u = d - p_ - q_; mult = math.factorial(d)//(math.factorial(p_)*math.factorial(q_)*math.factorial(u))
                tot += A(co)*mult*4**p_*Abr[(i, r + u)]*abs(Dq[m][q_](t, LN2))*om(t)
    return tot*eps(t)
J = [6.0273e-12, 8.4673e-12, 1.8708e-12]; J3e = 3.4315e-10
STIR = {1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}; BELL = {1: 1, 2: 2, 3: 5}
N = 3000; ts = [arb((i + 0.5)/N) for i in range(N)]
for j in (0, 1, 2):
    mons = T[j][(j + 1, j + 1)]
    c0 = [float(cd(mons, 0, t).mid()) for t in ts]; c2 = [float(cd(mons, 2, t).mid()) for t in ts]
    th = [float((eps(t)*om(t)*E2(t)).mid()) for t in ts]; cw = {m: [float((yy(t)**(2*m - 1)/Wt(t)).mid()) for t in ts] for m in (1, 2, 3)}
    print('level', j, 'pair', (j + 1, j + 1), 'rem monomials', {k: str(v) for k, v in mons.items() if k[0] >= 1})
    # direct route (data z_{j+1}): what F1 did
    k = j + 1; direct = BELL[k]*sum(STIR[k][m]*max(c0[i]*cw[m][i] for i in range(N))*(J + [J3e])[m] for m in STIR[k])
    print('   direct J_%d route: +%.4g' % (k, direct))
    if j >= 1:
        for delta in (0.2, 0.1, 0.05, 0.02):
            coef = [c0[i]**2/(4*delta*amin*th[i]) + c2[i]/2 for i in range(N)]
            val = BELL[j]*sum(STIR[j][m]*max(coef[i]*cw[m][i] for i in range(N))*J[m] for m in STIR[j])
            v2 = BELL[j]*sum(STIR[j][m]*max(c2[i]/2*cw[m][i] for i in range(N))*J[m] for m in STIR[j])
            print('   Young route delta = %.2f amin: + %.4g to E_%d (of which c\'\'/2 part %.3g)' % (delta, val, j, v2))
    # split route for level 2: t <= t0 by J3, t >= t0 by retained level-(j-1) dissipation  D = int int Theta z_{j+1}^2 <= Bud/(kappa*amin)
    if j == 2:
        for t0 in (0.05, 0.1, 0.15, 0.2, 0.3):
            i0 = int(t0*N)
            foot = BELL[3]*sum(STIR[3][m]*max(c0[i]*cw[m][i] for i in range(i0))*([0]+J[1:]+[J3e])[m] for m in STIR[3])
            ratio = max(c0[i]/th[i] for i in range(i0, N))
            print('   split t0 = %.2f: foot (J3) %.4g ; bulk ratio sup c/Theta = %.4g  -> x D1 (retained level-1 dissipation, <= Bud1/(kappa amin))' % (t0, foot, ratio))
