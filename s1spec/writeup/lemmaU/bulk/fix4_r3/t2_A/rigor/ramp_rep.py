#!/usr/bin/env python3
"""Phase-2a repair of Lemma S(A), stage 2: the inner-ramp [14,28] commutators, certified (exact Rational table coefficients, arb on a cover of [0,1]).
Ramp: y = 14*2^t, eta^2 = S(t)^2 (quintic smoothstep), Theta = eps*omega*eta^2, eps = y^-4, omega = 28/y, box (M0,M1) = (1/20, 9/10), U <= 1/20.
The cutoff part of a level-j pair coefficient is  c = eps*omega*co*rem_m(t)*A_i^(r)  (m >= 1),  rem_m*omega = (omega*eta^2)^(m) - (-1)^m omega*eta^2.

v3.1 (lemmaS_F2.ramp_level, rigor/ramp_arb.py) bounded only the JET-FREE part of A_i^(r+u) (Ab with V_j = 0 for j >= 2), reduced the pairs
(a, j+1), a < j, by repeated IBP to squares, and SET THE TOP PAIR (j+1, j+1) TO ZERO (lemmaS_F2.py:81, ramp_arb.py l.105, assembly_arb.py:111).
Here (all with the honest Ab of R3):
  R1 top pair (j+1,j+1): its cutoff coefficient c is jet-free (monomials (1,1,0),(1,2,0),(1,2,1),(2,2,0): A_1, A_2, A_2'), so either
       route J : int int c z_{j+1}^2 <= conv_{j+1}(c0)                          (verifier data J_{j+1}; J_3 = sqrt(J_2 J_4) + c_1W J_2)
       route IBP: int c z_{j+1}^2 = -int c z_j z_{j+2} + 1/2 int c'' z_j^2    (exact; c, c' vanish at t = 0, 1)
                  |c z_j z_{j+2}| <= delta a_min Theta z_{j+2}^2 + c0^2/(4 delta a_min Theta) z_j^2   (spare ramp dissipation)
                  c'' = Leibniz with A^(r+u), u <= 2: jet-free part -> data z_j^2; jet parts z_2, z_3 -> classes (j,j,2), (j,j,3).
  R2 cutoff x jet classes: every jet monomial X of A_i^(r+u) (incl. those created by the IBP derivatives) gives a class (indices of the pair) + X with an
       explicit weight function; low factors z_0, z_1 are folded by sup|z_0| <= 1/50, sup|z_1| <= 9/10 (as classify_x) down to two or three factors;
       handlers: 'data' (x,y <= 3: Cauchy-Schwarz in space-time against verifier data), 'young' (x, j+2: spare ramp dissipation, needs w/S bounded),
       'phi' (2,2,2): <= sup(w/(eps*omega)) * Phi,  Phi = bound of int int eps*omega|z_2|^3 (sup-free; the S33 data formula of v3.1),
       'ibp223' (2,2,3) whose z_3 coefficient is a TOP jet coefficient f(V) in {a'(V), b(V)}:  int kappa z_2^2 z_3 = -(1/3) int kappa' z_2^3,
       |kappa'| <= sum |coef|*eps*omega*[ sup|f| (|D_{q+1}| + 4|D_q|) + sup|f'| M1 |D_q| ]  -> 'phi'.
  Pairs (a, j+1), a < j: ONE integration by parts, c z_a z_{j+1} -> -c' z_a z_j - c z_{a+1} z_j (the full reduction of v3.1 creates jets z_5 of
       A_0^(5), A_1^(4) that no datum controls); (j, j+1): -1/2 c' z_j^2 as in v3.1.
Output: out/ramp_rep.json (sup tables against the conversion weights y^(2m-1)/W; the assembly turns them into energies).
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/ramp_rep.py [--N 1200]"""
import os, sys, json, math, argparse, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX, A, Q_, arb
import sympy as sp

STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}, 4: {1: 1, 2: 7, 3: 6, 4: 1}}
LN2 = arb(2).log(); X0 = arb(22490776515)/10**10
t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); GAM = 1; KAP = 5
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2s = Sx**2
M0R, M1R, UR = Q_(5, 100), Q_(9, 10), Q_(1, 20)
PM3_DIRECT = True        # task 2 / A3; --pm3-ibp restores the phase-2a convention (case (ii), point term missing)
SUPLOW = {0: A(Q_(1, 50)), 1: A(M1R)}

def build_D(qmax=5):
    """D[m][q] = (d/ds - gamma)^q rem_m (polynomials in t, Ls = ln 2), and D/t^3 when polynomial in t (for the Young weight w/S).
    Each entry is a triple (f, f_t, f_tt) of lambdified functions: cell enclosure |f| <= |f(mid)| + |f_t(mid)| h + |f_tt(cell)| h^2/2 (Taylor, rigorous)."""
    D, D3 = {}, {}
    def tri(g): return tuple(sp.lambdify((t_, Ls), sp.diff(g, t_, k), modules=[{}]) for k in range(3))
    for m in (1, 2, 3, 4):
        R = sum(sp.binomial(m, l)*sp.diff(e2s, t_, l)/Ls**l*(-GAM)**(m - l) for l in range(1, m + 1)); f = sp.expand(R)
        for q in range(qmax + 1):
            D[(m, q)] = tri(f)
            g = sp.cancel(f/t_**3); num, den = sp.fraction(sp.together(g))
            D3[(m, q)] = tri(sp.expand(g)) if t_ not in den.free_symbols else None
            f = sp.expand(sp.diff(f, t_)/Ls - GAM*f)
    return D, D3

def tenc(F, lo, hi, args=()):
    """rigorous upper bound of |f| on [lo,hi] by second-order Taylor about the midpoint (F = (f, f_t, f_tt))"""
    m = (lo + hi)/2; h = (hi - lo)/2; cell = lo.union(hi)
    return (abs(F[0](m, *args)) + abs(F[1](m, *args))*h + abs(F[2](cell, *args))*h*h/2).upper()

class Cells:
    """cover of [0,1] by N closed cells; per cell arb enclosures of all basis functions"""
    def __init__(self, N, D, D3):
        self.N = N; self.B = []; self.B3 = []; self.cw = []; self.cy = []; self.eo = []; self.S3 = []; self.G = []; self.rho = []
        tri1 = lambda g: tuple(sp.lambdify((t_, Ls), sp.diff(g, t_, k), modules=[{}]) for k in range(3))
        FG = tri1(e2s)                                                 # eta^2
        FR = tri1(sp.expand(sp.diff(e2s, t_, 2)/Ls**2 - 2*KAP*sp.diff(e2s, t_, 1)/Ls))   # (eta^2)'' - 2 kappa (eta^2)'  (s-derivatives; Ls = ln 2 ball)
        for i in range(N):
            lo, hi = arb(i)/N, arb(i + 1)/N
            t = lo.union(hi)
            y = 14*(t*LN2).exp(); om = 28/y; eps = y**-4; W = (1 + (y + X0)**2)**-2
            self.B.append({k: tenc(F, lo, hi, (LN2,)) for k, F in D.items()})
            self.B3.append({k: (tenc(F, lo, hi, (LN2,)) if F is not None else None) for k, F in D3.items()})
            self.eo.append((eps*om).upper())
            self.cw.append({m: (eps*om*y**(2*m - 1)/W).upper() for m in (1, 2, 3, 4)})
            self.cw[-1][0] = (eps*om/(y*W)).upper()
            self.cy.append({m: (y**(2*m - 1)/W).upper() for m in (0, 1, 2, 3, 4)})
            self.S3.append((10 - 15*t + 6*t*t).lower())          # S/t^3 >= 1 on [0,1]
            self.G.append((eps*om).upper()*tenc(FG, lo, hi, (LN2,)))
            self.rho.append((om*eps).upper()*tenc(FR, lo, hi, (LN2,))/2)
    def prof(self, i, comps):
        """sum coef*|D[m,q]| on cell i (per eps*omega)"""
        return sum((c*self.B[i][k] for c, k in comps), arb(0))
    def prof3(self, i, comps):
        tot = arb(0)
        for c, k in comps:
            b = self.B3[i][k]
            if b is None: return None
            tot += c*b
        return tot/self.S3[i]

def leib(d):
    for p in range(d + 1):
        for q in range(d + 1 - p):
            u = d - p - q; yield p, q, u, math.factorial(d)//(math.factorial(p)*math.factorial(q)*math.factorial(u))

def top_jet(i, R, X):
    """is the single jet X = (k,) the TOP jet of A_i^(R)?  -> f-type ('a1' = a'(V), 'b0' = b(V)) or None"""
    if len(X) != 1: return None
    if i == 2 and X[0] == R and R >= 2: return 'a1'
    if i == 1 and X[0] == R + 1: return 'b0'
    return None

def fold(I):
    """fold low factors z_0 (then z_1) by sup while more than two factors remain; returns (folded indices, sup factor)"""
    I = sorted(I); fac = arb(1)
    while len(I) > 2 and I[0] < 2:
        fac = fac*SUPLOW[I[0]]; I = I[1:]
    return tuple(I), fac

def level_terms(j, T, At, Abr, top_mode):
    """returns (jet-free data components per index pair, jet classes).  Components (coef, (m,q))."""
    n = j + 2; jf = {}; cls = {}
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        if b <= j: wts = [((a, b), 0, Q_(1))]
        elif a == b:                                 # top pair (j+1, j+1)
            if top_mode == 'J': wts = [((b, b), 0, Q_(1))]
            else: wts = [((j, n), 0, Q_(1)), ((j, j), 2, Q_(1, 2))]
        elif a == j: wts = [((j, j), 1, Q_(1, 2))]
        else: wts = [((a, j), 1, Q_(1)), ((a + 1, j), 0, Q_(1))]
        tag = 'top' if a == b == j + 1 else 'pairs'
        for (m, i, r), co in mons.items():
            if m == 0: continue
            if PM3_DIRECT and a == j and b == j + 1 and m == 3 and j in (0, 1):
                # task 2 / A3 (Lemma P case (iii)): the m = 3 monomial of pair (j, j+1) is NOT integrated by parts at levels 0, 1
                # (its weight jumps at y = 28, the IBP would leave +g_c z_j(28)^2); bounded directly by space-time Cauchy-Schwarz (conv_j, conv_{j+1})
                c0 = A(co)
                for X, val in At[(i, r)].items():
                    I, fac = fold([j, j + 1] + list(X))
                    cls.setdefault(('pairs', I, None), []).append(dict(c=c0*val*fac if X else c0*Abr[(i, r)], c_nof=c0*fac, k=(m, 0), i=i, R=r, X=X, pr=(j, j + 1), pm3=True))
                continue
            for pr, d, v in wts:
                for p, q, u, mult in leib(d):
                    R = r + u; c0 = A(v*mult*4**p*co)
                    if (i, R) not in At: raise ValueError(('A order beyond rmax', i, R))
                    for X, val in At[(i, R)].items():
                        if not X:
                            jf.setdefault((tag, pr), []).append((c0*Abr[(i, R)], (m, q)))
                            continue
                        if tag == 'top' and pr == (j, n): raise ValueError('jets in the top-pair coefficient itself')
                        I, fac = fold(list(pr) + list(X))
                        ft = top_jet(i, R, X) if sorted(I) == [2, 2, 3] and X == (3,) else None
                        key = (tag, I, ft)
                        cls.setdefault(key, []).append(dict(c=c0*val*fac, c_nof=c0*fac, k=(m, q), i=i, R=R, X=X, pr=pr))
    return jf, cls

def sup_tab(C, comps, idx, mode='data'):
    """sup over cells of profile * eps*omega*y^(2m-1)/W for m in STIR[idx]; mode 'young': (profile/S)^2 * eps*omega*y^(2m-1)/W"""
    out = {}
    for m in STIR[idx]:
        best = arb(0)
        for i in range(C.N):
            if mode == 'data': v = C.prof(i, comps)
            else:
                v = C.prof3(i, comps)
                if v is None: return None
                v = v*v
            v = (v*C.cw[i][m]).upper()
            if v > best: best = v
        out[m] = best
    return out

def sup_plain(C, comps):
    best = arb(0)
    for i in range(C.N):
        v = C.prof(i, comps).upper()
        if v > best: best = v
    return best

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--N', type=int, default=1200); ap.add_argument('--pm3-ibp', action='store_true'); ap.add_argument('--weights', default=L.OUT + '/weights_t2.json'); ap.add_argument('--out', default=None); a = ap.parse_args()
    PM3_DIRECT = not a.pm3_ibp
    t0 = time.time()
    T = WX.build_tables(); At = L.at_eval(M0R, M1R, UR); Abr = L.ab_honest(At)
    WR = L.load(a.weights)
    Fa2 = WX.sup_rat_arb('a', 2, M0R); Fb1 = WX.sup_rat_arb('b', 1, M0R); Fa1 = WX.sup_rat_arb('a', 1, M0R); Fb0 = WX.sup_rat_arb('b', 0, M0R)
    FSUP = {'a1': (Fa1, Fa2), 'b0': (Fb0, Fb1)}
    D, D3 = build_D(); C = Cells(a.N, D, D3); print('cells built %.0fs' % (time.time() - t0), flush=True)
    res = {}
    for j in (0, 1, 2):
        n = j + 2; Ghr = {int(k): v.upper() for k, v in WR['ramp'][str(j)]['Ghat'].items()}; YCr = {int(k): v.upper() for k, v in WR['ramp'][str(j)]['YC'].items()}
        R = dict()
        for mode in ('J', 'IBP'):
            jf, cls = level_terms(j, T, At, Abr, mode)
            # ---- data tables q_k (k <= j): jet-free 'pairs' parts, 1/2-1/2 split of mixed pairs, + ramp share of the rate + interpolation remainder
            if mode == 'J':
                qk = {}
                for k in range(j + 1):
                    comps = []
                    for (tag, pr), cl in jf.items():
                        if tag != 'pairs': continue
                        x, y = pr
                        if x == y == k: comps += cl
                        elif k in pr: comps += [(c/2, kk) for c, kk in cl]
                    tab = {}
                    for m in STIR[k]:
                        best = arb(0)
                        for i in range(C.N):
                            v = (C.prof(i, comps)*C.cw[i][m] + (Ghr[k]*C.G[i] + (YCr[j + 1]*C.rho[i] if k == j else 0))*C.cy[i][m]).upper()
                            if v > best: best = v
                        tab[m] = best
                    qk[k] = tab
                R['q'] = qk
                for (tag, pr), cl in jf.items():
                    if tag == 'pairs' and max(pr) > j: raise ValueError(('jet-free pair beyond j', pr))
                topc = jf.get(('top', (j + 1, j + 1)), [])
                R['topJ'] = sup_tab(C, topc, j + 1)
                R['topJ_sup_per_eo'] = sup_plain(C, topc)
            else:
                c0 = jf.get(('top', (j, n)), []); c2 = jf.get(('top', (j, j)), [])
                R['topY'] = sup_tab(C, c0, j, 'young')            # (c0/(eps om S))^2 eps om y^(2m-1)/W ; cost conv_j(.)/(4 delta a_min)
                R['top2'] = sup_tab(C, c2, j)                     # 1/2 |c''| jet-free (weight 1/2 included)
            # ---- jet classes
            if mode == 'J': assert not any(k[0] == 'top' for k in cls), ('jets in the top-pair coefficient (route J would drop them)', j)
            out = []
            for (tag, I, ft), lst in sorted(cls.items(), key=lambda kv: (kv[0][0], kv[0][1], str(kv[0][2]))):
                if mode == 'J' and tag == 'top': continue
                if mode == 'IBP' and tag != 'top': continue
                comps = [(e['c'], e['k']) for e in lst]; ent = dict(tag=tag, I=list(I), ft=ft, n_terms=len(lst), sup_w=sup_plain(C, comps))
                srcs = {}
                for e in lst: srcs.setdefault(str(e['pr']) + ' X=' + str(e['X']), []).append((e['c'], e['k']))
                ent['sup_w_by_source'] = {k: sup_plain(C, v) for k, v in srcs.items()}
                if len(I) == 2:
                    x, y = I; ent['handler'] = 'pair'
                    ent['data'] = {str(k): sup_tab(C, comps, k) for k in set(I) if k <= 3}
                    if y == n: ent['young'] = sup_tab(C, comps, x, 'young')
                    if y > 3 and y != n: raise ValueError(('class beyond data and dissipation', j, I))
                elif sorted(I) == [2, 2, 2]: ent['handler'] = 'phi'
                elif sorted(I) == [2, 2, 3]:
                    if any(top_jet(e['i'], e['R'], e['X']) is None for e in lst): raise ValueError(('(2,2,3) with a non-top z_3 coefficient', j, lst))
                    ent['handler'] = 'ibp223'
                    dcomps = []
                    for e in lst:
                        f, fp = FSUP[top_jet(e['i'], e['R'], e['X'])]; m, q = e['k']
                        dcomps += [(e['c_nof']*f, (m, q + 1)), (e['c_nof']*f*4, (m, q)), (e['c_nof']*fp*A(M1R), (m, q))]
                    ent['sup_dw'] = sup_plain(C, dcomps)          # sup |kappa'|/(eps omega)
                    ent['check_sup_w_vs_top'] = sup_plain(C, [(e['c_nof']*FSUP[top_jet(e['i'], e['R'], e['X'])][0], e['k']) for e in lst])
                else: raise ValueError(('unhandled class', j, I))
                out.append(ent)
            R['classes_' + mode] = out
        res[str(j)] = R
        print('level', j, 'q', {k: {m: '%.6g' % float(v) for m, v in d.items()} for k, d in R['q'].items()}, flush=True)
        print('   topJ', {m: '%.4g' % float(v) for m, v in R['topJ'].items()}, 'topY', R['topY'] and {m: '%.4g' % float(v) for m, v in R['topY'].items()},
              'top2', {m: '%.4g' % float(v) for m, v in R['top2'].items()}, flush=True)
        for md in ('J', 'IBP'):
            for e in R['classes_' + md]: print('   class[%s] %-4s %-10s %-8s sup w/(eps om) %.4g %s' % (md, e['tag'], e['I'], e['handler'], float(e['sup_w']), ('sup dw %.4g' % float(e['sup_dw'])) if 'sup_dw' in e else ''), flush=True)
    inputs = [L.BASE + '/rigor/weights_exact.py', L.BASE + '/rigor/tables_exact.pkl', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py', a.weights]
    gc28 = A(28)**-4*A(120)/LN2**3        # Lemma P: g_c at y = 28 (eps*omega*a*J_in with omega(28) = 1, a <= 1), per unit tau
    pt = dict(gc_28=gc28, routes={'0': 'case (iii): pair (0,1) m=3 by data (no IBP); top pair route J', '1': 'case (iii): pair (1,2) m=3 by data (no IBP); top pair route J',
              '2': 'case (i): top pair (3,3) and pair (2,3) both by IBP; net point term -g_c z_2(28)^2 <= 0 dropped'} if PM3_DIRECT else 'phase-2a convention (case (ii) at levels 0,1)')
    L.dump(dict(res=res, N=a.N, box=dict(M0=str(M0R), M1=str(M1R), U=str(UR)), pointterm=pt), a.out or (L.OUT + '/ramp_rep.json'), inputs)
    print('wrote out/ramp_rep.json  %.0fs' % (time.time() - t0))
