# Reviewer scratch: reconcile P1/P2 Ghat conventions with P1's certified pipeline (ghat_arb functions, P1's Gamma/kappa enclosures).
import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/phase1/ghat')
import sympy as sp
from flint import arb
import ghat_arb as GA, weight_consts_arb as WC
Q_ = sp.Rational; A = GA.A
T = pickle.load(open(GA.TABLES, 'rb'))
r = WC.certify(4, 4, Q_(0), Q_(3, 10), Q_(3, 5), 56); b = r['bound']
Gam = {0: arb(1), 1: b['G1'], 2: b['G2'], 3: b['G3']}; kap = b['kap']; kap2 = b['kap2']
LN2 = arb(2).log(); etap = A(Q_(15, 8))/LN2; DEL = A(Q_(1, 10))
def run(M0q, M1q, Umq, M0zq, dfrac, foot, honest=False, Gam=Gam, kap=kap, kap2=kap2, s33='ibp'):
    Asym = GA.build_At_sym(Umq); Ab = GA.ab_exact(M0q, M1q, 4, Umq); At = GA.eval_At(Asym, M0q, M1q); out = []
    if honest: Ab = {k: d.get((), arb(0)) for k, d in At.items()}
    for j in (0, 1, 2):
        YCb = GA.YC_base(T[j], j, Ab, Gam); add, dis, spare, ndfr, items = GA.l4_extra_YC_x(T, j, M0q, M1q, At, Gam, kap, dfrac=dfrac, s33=s33, M0zq=M0zq)
        fadd = {k: arb(0) for k in add}; extra = Q_(0); M1 = A(M1q)
        if foot:
            w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0))
            if any(it[0] == 'S33' for it in items):
                c1 = w33*M1*etap; fadd[3] += c1*A(kap)*DEL; extra = Q_(1, 10)
            for it in items:
                if it[0] != 'L4': continue
                w = it[2]
                if it[1] == 3: fadd[3] += w*M1*etap*DEL
                elif it[1] == 2: fadd[2] += w*M1*etap
        ok, sl = GA.budget_ok(dis, ndfr, M0q, extra)
        YC = {k: YCb[k] + add.get(k, 0) + fadd.get(k, 0) for k in YCb}
        out.append((float(GA.ghat_x(YC, j, M0q, kap2)[j].upper()), bool(ok), round(sl, 3)))
    return out
cases = [('P1 as is', (Q_(11,100), Q_(52,100), Q_(1,20), Q_(1,50), Q_(1,10), False)),
         ('P1 Umax .09', (Q_(11,100), Q_(52,100), Q_(9,100), Q_(1,50), Q_(1,10), False)),
         ('P2 conv dfrac 1/10 no foot', (Q_(11,100), Q_(52,100), Q_(11,100), Q_(11,100), Q_(1,10), False)),
         ('P2 conv dfrac 2/25 + foot', (Q_(11,100), Q_(52,100), Q_(11,100), Q_(11,100), Q_(2,25), True)),
         ('support: U .09, M0z .16, 2/25 + foot', (Q_(11,100), Q_(52,100), Q_(9,100), Q_(16,100), Q_(2,25), True)),
         ('zones I/II: U .05, M0z .16, 2/25 + foot', (Q_(11,100), Q_(52,100), Q_(1,20), Q_(16,100), Q_(2,25), True)),
         ('M0 box .15: U .09, M0z .16, 2/25 + foot', (Q_(15,100), Q_(52,100), Q_(9,100), Q_(16,100), Q_(2,25), True)),
         ('M1 box .45: U .09, M0z .16, 2/25 + foot', (Q_(11,100), Q_(45,100), Q_(9,100), Q_(16,100), Q_(2,25), True))]
one = {m: arb(1) for m in range(4)}
for name, c in cases + [(n + ' [Ab incl. U_j]', c + (True,)) for n, c in cases] + [('(A) core power wt ibp', (Q_(11,100), Q_(52,100), Q_(1,20), Q_(1,50), Q_(1,10), False, False, one, 5, 25)), ('(A) core power wt ibp [Ab incl. U_j]', (Q_(11,100), Q_(52,100), Q_(1,20), Q_(1,50), Q_(1,10), False, True, one, 5, 25))]:
    o = run(*c); G22 = o[2][0]
    print('%-42s Ghat_jj = %.1f %.1f %.1f  budget %s slack %s | rate margin %.3f Gw %.3f' % (name, o[0][0], o[1][0], G22, [x[1] for x in o], [x[2] for x in o], 28**4/G22, 1/(1 - G22/28**4)), flush=True)
