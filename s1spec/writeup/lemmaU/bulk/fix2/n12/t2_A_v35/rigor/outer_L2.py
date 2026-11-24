#!/usr/bin/env python3
"""Task 2 / A2: R5' (STATEMENTS_r2.tex Prop. 'R5''), the cutoff terms of (A) on the outer annulus O = [350, 700] in L^2 form, certified (arb).
Replaces rigor/outer_rep.py.  No pointwise bound on z_2, z_3 is used; data are the L^2 norms Nbar_k^2 >= sup_tau int_O (z^A_k)^2 ds, k = 0,1,2 (H-O2'),
the box H-O1 (T1: |z_0| <= Z0 = 331/500, |z_1| <= Z1 = 1166/1000; coefficient box (1/2, 1), U <= 0.162) and the reserved outer dissipation
delta_o a_o int Theta z_n^2 (weights_t2.py lowers B_c on O by delta_o).
Cutoff: eta_o = 1 - S(t), t = ln(y/350)/ln 2, Theta = eps*omega*eta_o^2, a_o = (1 + (1/2)^2)^-2.
Class enumerator: that of ramp_rep.py, mirrored to eta_o (cutoff part of a level-j pair = eps*omega*co*rem_m*A_i^(r), m >= 1), with the routes
  levels 0, 1: no integration by parts at all (Lemma P case (iii));  level 2: top pair (3,3) and pair (2,3) integrated by parts (Lemma P case (i):
  net point term -g_c z_2(350)^2 <= 0, dropped), pairs (a,3), a < 2, one integration by parts (T4; m <= 1, so c is C^1).
Tools per class (asserted; anything else raises):
  pair (x,y), x,y <= 2         : T2   sup|w| (AM-GM split fixed at Nbar_ref)
  pair (x,n), x <= 2           : T3   Young against the reserve; needs sup w^2/Theta (asserted finite: vanishing order >= 3 at t = 1)
  pair (2,3)  [level 2]        : envelope w <= K e, e = eps*om*t(1-t)^4;  int e|z2||z3| <= mu/2 int theta_* z3^2 + 1/(2mu) sup(e^2/theta_*) Nbar_2^2;  T5 on theta_*
  (2,2,2)     [level 2]        : T6 with theta = K e (continuous, piecewise C^1, zero at both ends): <= Z1 (sup|theta'| Nbar_2^2 + 2 int theta|z2||z3|) -> (2,3) chain
  (2,2,3)     [level 2]        : T7 (top-jet coefficient a'(V) or b(V); kappa continuous): = -1/3 int kappa' z2^3, |kappa'| <= K' e -> (2,2,2)
  T5: int theta_* z3^2 <= delta' a_o int Theta z4^2 + sup[theta_*^2/(4 delta' a_o Theta) + |theta_*''|/2] Nbar_2^2,  theta_* = eps*om*t^2(1-t)^3.
Other cutoff sources on O: the interpolation remainder (1/2) Ytilde_o eps*om |(eta_o^2)'' - 2 kappa (eta_o^2)'| z_j^2 (T2), and the L4/S33 outer-foot
terms (computed in weights_t2.py: 'outer_foot_N2', 'outer_foot_reserve').
Output: out/outer_L2.json: kappa_o[j][k] (o_j = sum_k kappa_o[j][k] Nbar_k^2 per unit tau), the class list with tools, reserve usage.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/outer_L2.py [--weights out/weights_t2.json] [--N2ref 0.1 0.2 5.0]"""
import os, sys, json, math, argparse, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX, A, Q_, arb
import ramp_rep as RR
import sympy as sp

LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); KAP = 5
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2d = sp.expand((1 - Sx)**2)
M0o, M1o, Uo = Q_(1, 2), Q_(1), Q_(162, 1000)
Z0, Z1 = A(Q_(331, 500)), A(Q_(1166, 1000)); SUP = {0: Z0, 1: Z1}
ENV = t_*(1 - t_)**4                                  # envelope e/(eps*om)
THS = t_**2*(1 - t_)**3                               # theta_*/(eps*om)

def tri(g): return tuple(sp.lambdify((t_, Ls), sp.diff(g, t_, k), modules=[{}]) for k in range(3))
def tenc(F, lo, hi):
    m = (lo + hi)/2; h = (hi - lo)/2; cell = lo.union(hi)
    return (abs(F[0](m, LN2)) + abs(F[1](m, LN2))*h + abs(F[2](cell, LN2))*h*h/2).upper()

def build():
    D, DE, D3 = {}, {}, {}
    for m in (1, 2, 3, 4):
        f = sp.expand(sum(sp.binomial(m, l)*sp.diff(e2d, t_, l)/Ls**l*(-1)**(m - l) for l in range(1, m + 1)))
        for q in range(5):
            D[(m, q)] = tri(f)
            g = sp.cancel(f/ENV); DE[(m, q)] = tri(sp.expand(g)) if t_ not in sp.fraction(sp.together(g))[1].free_symbols else None
            g3 = sp.cancel(f/(1 - t_)**3); D3[(m, q)] = tri(sp.expand(g3)) if t_ not in sp.fraction(sp.together(g3))[1].free_symbols else None
            f = sp.expand(sp.diff(f, t_)/Ls - f)
    return D, DE, D3

class Cells:
    def __init__(self, N, D, DE, D3):
        self.N = N; self.B, self.BE, self.B3, self.eo, self.aux = [], [], [], [], []
        dds = lambda g: sp.diff(g, t_)/Ls
        Fe1 = tri(sp.expand(dds(ENV) - KAP*ENV))                          # (eps om e)'/(eps om)
        Fth2 = tri(sp.expand(dds(dds(THS)) - 2*KAP*dds(THS) + KAP**2*THS))  # (eps om theta_*)''/(eps om)
        Fth2Q = tri(sp.expand(t_**4))                                     # theta_*^2/(eps om)/Theta*(eps om) = t^4/(1+3t+6t^2)^2 (divide below)
        Frho = tri(sp.expand(sp.diff(e2d, t_, 2)/Ls**2 - 2*KAP*sp.diff(e2d, t_)/Ls))
        Fe2 = tri(sp.expand(ENV**2/THS))                                  # e^2/theta_* /(eps om) = (1-t)^5 after cancel
        Fe2 = tri(sp.cancel(ENV**2/THS))
        for i in range(N):
            lo, hi = arb(i)/N, arb(i + 1)/N; t = lo.union(hi)
            y = 350*(t*LN2).exp(); eo = (y**-4*28/y)
            self.eo.append(eo.upper())
            self.B.append({k: tenc(F, lo, hi) for k, F in D.items()})
            self.BE.append({k: (tenc(F, lo, hi) if F is not None else None) for k, F in DE.items()})
            self.B3.append({k: (tenc(F, lo, hi) if F is not None else None) for k, F in D3.items()})
            p3 = (1 + 3*t + 6*t*t).lower()                                 # (1 - S) = (1-t)^3 (1 + 3t + 6t^2)
            self.aux.append(dict(e1=tenc(Fe1, lo, hi), th2=tenc(Fth2, lo, hi), thQ=tenc(Fth2Q, lo, hi)/arb(p3)**2, rho=tenc(Frho, lo, hi)/2, e2=tenc(Fe2, lo, hi), p3=p3))
    def sup(self, f):
        best = arb(0)
        for i in range(self.N):
            v = f(i)
            if v is None: return None
            v = v.upper() if isinstance(v, arb) else v
            if v > best: best = v
        return best

def leib(d): yield from RR.leib(d)

def fold(I):
    I = sorted(I); fac = arb(1)
    while len(I) > 2 and I[0] < 2: fac = fac*SUP[I[0]]; I = I[1:]
    return tuple(I), fac

def level_classes(j, T, At):
    n = j + 2; cls = {}
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        if j < 2 or b <= j: wts = [((a, b), 0, Q_(1))]; route = 'direct'
        elif a == b: wts = [((j, n), 0, Q_(1)), ((j, j), 2, Q_(1, 2))]; route = 'P(i): IBP top pair'
        elif a == j: wts = [((j, j), 1, Q_(1, 2))]; route = 'P(i): IBP (j,j+1)'
        else: wts = [((a, j), 1, Q_(1)), ((a + 1, j), 0, Q_(1))]; route = 'T4: one IBP'
        for (m, i, r), co in mons.items():
            if m == 0: continue
            if route == 'T4: one IBP': assert m <= 1, ('T4 needs c in C^1', j, (a, b), m)
            for pr, d, v in wts:
                for p, q, u, mult in leib(d):
                    R = r + u; c0 = A(v*mult*co)
                    for X, val in At[(i, R)].items():
                        I, fac = fold(list(pr) + list(X))
                        ft = RR.top_jet(i, R, X) if sorted(I) == [2, 2, 3] and X == (3,) else None
                        cls.setdefault((I, ft), []).append(dict(c=c0*val*fac, c_nof=c0*fac, k=(m, q), i=i, R=R, X=X, pr=pr, src=(a, b), route=route))
    return cls

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--weights', default=L.OUT + '/weights_t2.json'); ap.add_argument('--N', type=int, default=1200)
    ap.add_argument('--N2ref', nargs=3, default=['0.1', '0.2', '5.0']); ap.add_argument('--out', default=L.OUT + '/outer_L2.json')
    a = ap.parse_args(); t0 = time.time()
    WR = L.load(a.weights); T = WX.build_tables(); At = L.at_eval(M0o, M1o, Uo)
    a_o = (1 + A(M0o)**2)**-2
    Fa1, Fa2 = WX.sup_rat_arb('a', 1, M0o), WX.sup_rat_arb('a', 2, M0o); Fb0, Fb1 = WX.sup_rat_arb('b', 0, M0o), WX.sup_rat_arb('b', 1, M0o)
    FSUP = {'a1': (Fa1, Fa2), 'b0': (Fb0, Fb1)}
    D, DE, D3 = build(); C = Cells(a.N, D, DE, D3); print('cells %.0fs' % (time.time() - t0), flush=True)
    Nref = [A(Q_(x)) for x in a.N2ref]
    def mu_xy(x, y): return A(Q_(math.sqrt(float(Q_(a.N2ref[y])/Q_(a.N2ref[x]))))).mid() if x != y else arb(1)
    supE = lambda comps: C.sup(lambda i: sum((c*C.B[i][k] for c, k in comps), arb(0))*C.eo[i])                 # sup |w|
    def supEnv(comps):                                                                                       # sup |w|/e  (K)
        def f(i):
            tot = arb(0)
            for c, k in comps:
                b = C.BE[i][k]
                if b is None: return None
                tot += c*b
            return tot
        return C.sup(f)
    def supW2Th(comps):                                                                                      # sup w^2/Theta
        def f(i):
            tot = arb(0)
            for c, k in comps:
                b = C.B3[i][k]
                if b is None: return None
                tot += c*b
            return (tot/C.aux[i]['p3'])**2*C.eo[i]
        return C.sup(f)
    sup_e1 = C.sup(lambda i: C.aux[i]['e1']*C.eo[i]); sup_e2 = C.sup(lambda i: C.aux[i]['e2']*C.eo[i])
    sup_thQ = C.sup(lambda i: C.aux[i]['thQ']*C.eo[i]); sup_th2 = C.sup(lambda i: C.aux[i]['th2']*C.eo[i]); sup_rho = C.sup(lambda i: C.aux[i]['rho']*C.eo[i])
    res = {}
    for j in (0, 1, 2):
        n = j + 2; LV = WR['levels'][str(j)]; dlo = A(Q_(LV['O']['delta_o'])); Yo = LV['Ytilde']['plateau_O']
        kap = {k: arb(0) for k in range(3)}; classes = []; T3 = []; chainA = arb(0)
        cls = level_classes(j, T, At)
        for (I, ft), lst in sorted(cls.items(), key=lambda kv: str(kv[0])):
            comps = [(e['c'], e['k']) for e in lst]; ent = dict(I=list(I), ft=ft, n_terms=len(lst), routes=sorted(set(e['route'] for e in lst)), mq=sorted(set(e['k'] for e in lst)))
            if len(I) == 2 and max(I) <= 2:
                s = supE(comps); x, y = I; ent.update(tool='T2', sup_w=s)
                if x == y: kap[x] += s; ent['cost'] = {str(x): s}
                else:
                    mu = mu_xy(x, y); kap[x] += s*mu/2; kap[y] += s/(2*mu); ent['cost'] = {str(x): s*mu/2, str(y): s/(2*mu)}; ent['mu'] = mu
            elif len(I) == 2 and I[1] == n and I[0] <= 2:
                s = supW2Th(comps); assert s is not None, ('T3: w^2/Theta unbounded', j, I)
                ent.update(tool='T3', sup_w2_over_Theta=s); T3.append((I[0], s, ent))
            elif j == 2 and I == (2, 3):
                K = supEnv(comps); assert K is not None, ('(2,3): no envelope t(1-t)^4', j)
                ent.update(tool='envelope + AM-GM -> T5 (theta_*) + T2', K=K); chainA += K
            elif j == 2 and I == (2, 2, 2):
                K = supEnv(comps); assert K is not None, ('(2,2,2): no continuous envelope', j)
                ent.update(tool='T6 (theta = K e) -> T2 + (2,3)-chain', K=K); kap[2] += Z1*K*sup_e1; chainA += 2*Z1*K
                ent['cost'] = {'2': Z1*K*sup_e1}
            elif j == 2 and I == (2, 2, 3):
                assert ft in ('a1', 'b0') and all(RR.top_jet(e['i'], e['R'], e['X']) == ft for e in lst), ('(2,2,3) non-top', lst)
                dcomps = []
                for e in lst:
                    f, fp = FSUP[ft]; m, q = e['k']; assert m + q <= 2, ('T7: kappa must be continuous', e)
                    dcomps += [(e['c_nof']*f, (m, q + 1)), (e['c_nof']*f*4, (m, q)), (e['c_nof']*fp*A(M1o), (m, q))]
                Kp = supEnv(dcomps); assert Kp is not None, ('(2,2,3): kappa\' has no envelope', j)
                ent.update(tool='T7 -> (2,2,2) [T6] -> T2 + (2,3)-chain', Kprime=Kp); kap[2] += Z1*Kp/3*sup_e1; chainA += 2*Z1*Kp/3
                ent['cost'] = {'2': Z1*Kp/3*sup_e1}
            else:
                raise SystemExit('REPORT: class %s at level %d is not closed by T1-T7, P (CONTRACTS r1: report, do not bound pointwise)' % (I, j))
            classes.append(ent)
        # interpolation remainder on O: (1/2) Ytilde_o eps om |(eta_o^2)'' - 2 kappa (eta_o^2)'| z_j^2   (T2)
        rho = Yo*sup_rho; kap[j] += rho
        # outer foot (weights_t2)
        kap[2] += LV['O']['outer_foot_N2']; res_foot = LV['O']['outer_foot_reserve']
        # reserve allocation: T3 classes and the chain share  dlo - res_foot  (units a_o)
        budget = dlo - res_foot; assert budget.lower() > 0
        parts = []
        for x, s, ent in T3: parts.append(('T3', x, s, ent))
        # float optimisation of the shares (any positive shares summing to <= budget are valid; values then evaluated in arb)
        Bf = float(budget.lower())*0.999
        def chain_cost(Dc, A_):            # A_ int e|z2||z3|: min over mu of A_[(mu/2)(c1/delta' + c2) + s_e/(2mu)], dissipation A_ mu delta'/2 = Dc
            c1 = float(sup_thQ.upper())/(4*float(a_o.lower())); c2 = float(sup_th2.upper())/2; se = float(sup_e2.upper())
            best = None
            for k in range(-60, 61):
                mu = 10**(k/10); dl = 2*Dc/(A_*mu)
                v = A_*(mu/2*(c1/dl + c2) + se/(2*mu))
                if best is None or v < best[0]: best = (v, mu, dl)
            return best
        Af = float(chainA.upper()) if chainA > 0 else 0.0
        ao = float(a_o.lower()); N_ = [float(x.upper()) for x in Nref]
        sq = [(float(s.upper())*N_[x])**0.5 for (_, x, s, _) in parts]
        best = None
        for fi in ([0.0] if Af == 0 else [k/1000 for k in range(1, 1000)]):    # chain fraction; T3 shares optimal (prop. to sqrt(s N)) in the rest
            rest = Bf*(1 - fi); tot3 = sum(sq) or 1.0
            c3 = sum((float(s.upper())*N_[x])/(4*ao*rest*q/tot3) for (_, x, s, _), q in zip(parts, sq)) if parts else 0.0
            cc = chain_cost(Bf*fi, Af)[0]*N_[2] if Af > 0 else 0.0
            if best is None or c3 + cc < best[0]: best = (c3 + cc, fi)
        fi = best[1]; tot3 = sum(sq) or 1.0
        shares = [Bf*(1 - fi)*q/tot3 for q in sq] + ([Bf*fi] if Af > 0 else [])
        used = arb(0); alloc_info = []
        for (tag, x, s, ent), sh in zip(parts, shares[:len(parts)]):
            dl = A(Q_(math.floor(sh*1e9), 10**9)); assert dl > 0
            c = s/(4*dl*a_o); kap[x] += c; used += dl; ent['delta'] = dl; ent['cost'] = {str(x): c}; alloc_info.append((tag, x, float(dl)))
        if Af > 0:
            Dc = shares[-1]; v, mu, dl = chain_cost(Dc, Af)
            muq = A(Q_(mu).limit_denominator(10**6)); dlq = A(Q_(math.floor(float(2*Dc/(Af*mu))*0.999*1e12), 10**12))
            dcost = chainA*(muq/2*(sup_thQ/(4*dlq*a_o) + sup_th2/2) + sup_e2/(2*muq)); kap[2] += dcost
            used += chainA*muq*dlq/2; alloc_info.append(('chain', 2, float((chainA*muq*dlq/2).upper())))
            chain = dict(A=chainA, mu=muq, delta_T5=dlq, cost_N2=dcost, sup_e2_over_theta=sup_e2, sup_thetaQ=sup_thQ, sup_theta_pp=sup_th2)
        else: chain = None
        assert (used + res_foot - dlo).upper() <= 0, ('reserve exceeded', j)
        res[str(j)] = dict(kappa=[kap[k] for k in range(3)], classes=classes, chain=chain, interp_remainder=dict(Ytilde_o=Yo, sup=sup_rho, cost=rho),
                           outer_foot=dict(N2=LV['O']['outer_foot_N2'], reserve=res_foot), reserve=dict(delta_o=dlo, used=used + res_foot, alloc=alloc_info),
                           pointterm='level %d: %s' % (j, 'Lemma P case (i) at y = 350 (net -g_c z_2(350)^2 <= 0, dropped)' if j == 2 else 'Lemma P case (iii) (no integration by parts on O)'))
        o_ref = sum((kap[k]*Nref[k] for k in range(3)), arb(0))
        print('level %d: kappa_o = %s ; o_j at N2ref = %.4g, at 4 N2ref = %.4g ; reserve used %.4f of %s' % (j, ['%.4g' % L.fup(kap[k]) for k in range(3)], L.fup(o_ref), L.fup(4*o_ref), L.fup(used + res_foot), LV['O']['delta_o']), flush=True)
        for e in classes: print('    class %-10s %-6s %-55s %s' % (e['I'], e['ft'] or '', e['tool'], {k: '%.3g' % L.fup(v) for k, v in e.get('cost', {}).items()}), flush=True)
    inputs = [a.weights, L.BASE + '/rigor/weights_exact.py', L.BASE + '/rigor/tables_exact.pkl', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py', L.ROOT + '/rigor/ramp_rep.py']
    L.dump(dict(res=res, N2ref=a.N2ref, box=dict(M0o=str(M0o), M1o=str(M1o), U=str(Uo), Z0='331/500', Z1='1166/1000'), N=a.N,
                sups=dict(e1=sup_e1, e2=sup_e2, thetaQ=sup_thQ, theta_pp=sup_th2, rho=sup_rho)), a.out, inputs)
    print('wrote', a.out, '%.0fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
