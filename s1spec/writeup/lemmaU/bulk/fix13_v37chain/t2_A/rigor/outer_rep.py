#!/usr/bin/env python3
"""Phase-2a repair of Lemma S(A), stage 3 (R5): the outer cutoff on y in [350, 700], certified given EXPLICIT HYPOTHESES (symbolic thresholds).
Energy weight: omega*eta^2, eta^2 = eta_in^2 * eta_out^2, eta_out^2 = (1 - S(t))^2, t = ln(y/350)/ln 2 (decreasing quintic cutoff), Theta = eps*omega*eta^2.
Split (pointwise integrands; the quadratic-form identity of Prop. 1 is global, its total derivatives vanish at both ends of the support):
  * the core identities (z_{j+1}-interpolation, L4 identity, S33 integration by parts) are taken over [14, 350] (eta_out = 1 there); their boundary
    terms at y = 350 are bounded pointwise by the hypotheses below (R_b, R_e) or, inside the L4 square root, by lam*N_r^2 + w^2 O_b/(4 lam)
    (YC adds lam, dissipation eps4: weights_rep.py) (R_d), O_b = Theta(350)*M1bar*M2bar^3;
  * EVERY quadratic, cubic and quartic term on [350, 700] is bounded term by term: cutoff (rem) parts of all pairs incl. the top pair (R_a); base parts
    (weight Theta) of all pairs except (n,n) (R_c); a factor z_n (n = j+2) is absorbed by the outer dissipation 2 a_min^o Theta z_n^2 (Young, delta_o),
    every other factor z_k by its pointwise bound.
HYPOTHESES (for all tau >= tau_1, on 350 <= y <= 700; NOT certified here, to be supplied by Lemma U zone III / Lemma O / (B)):
  H-O1  coefficient box |f| <= M0o, y|f'| <= M1o (and |U| <= 1.3056 alpha_1 700^p < 0.162, certified);
  H-O2  deviation jets |z_k| <= Mbar_k, k = 0..3 (z_k = (y d/dy)^k v).
The rates are per unit tau and enter E_j as rate/gamma'_j.  Output: out/outer_rep.json with the rate as an explicit function of (Mbar_0..3) (list of
monomials coef * prod Mbar^e) plus its value at the recorded placeholder values.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/outer_rep.py [--box 1/2 1] [--Mbar 1 1 1 1]"""
import os, sys, json, math, argparse, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX, A, Q_, arb
import sympy as sp

LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); GAM = 1; KAP = 5
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2d = (1 - Sx)**2
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}

def tri(g): return tuple(sp.lambdify((t_, Ls), sp.diff(g, t_, k), modules=[{}]) for k in range(3))
def tenc(F, lo, hi):
    m = (lo + hi)/2; h = (hi - lo)/2; cell = lo.union(hi)
    return (abs(F[0](m, LN2)) + abs(F[1](m, LN2))*h + abs(F[2](cell, LN2))*h*h/2).upper()

def sups(N=600):
    """sup over [350,700] of eps*omega*|rem_m| (m=1..4), eps*omega*eta^2, eps*omega*|(eta^2)'' - 2 kappa (eta^2)'|/2"""
    Fr = {m: tri(sp.expand(sum(sp.binomial(m, l)*sp.diff(e2d, t_, l)/Ls**l*(-GAM)**(m - l) for l in range(1, m + 1)))) for m in (1, 2, 3, 4)}
    Fe = tri(sp.expand(e2d)); Fh = tri(sp.expand(sp.diff(e2d, t_, 2)/Ls**2 - 2*KAP*sp.diff(e2d, t_, 1)/Ls))
    out = {('rem', m): arb(0) for m in Fr}; out['Theta'] = arb(0); out['rho'] = arb(0)
    for i in range(N):
        lo, hi = arb(i)/N, arb(i + 1)/N; t = lo.union(hi); y = 350*(t*LN2).exp(); eo = (y**-4*28/y).upper()
        for m, F in Fr.items(): out[('rem', m)] = max(out[('rem', m)], eo*tenc(F, lo, hi))
        out['Theta'] = max(out['Theta'], eo*tenc(Fe, lo, hi)); out['rho'] = max(out['rho'], eo*tenc(Fh, lo, hi)/2)
    return out

def mono(idx):
    e = [0, 0, 0, 0]
    for k in idx: e[k] += 1
    return tuple(e)

def addm(P, e, c): P[e] = P.get(e, arb(0)) + c

def evalP(P, Mb):
    tot = arb(0)
    for e, c in P.items():
        v = c
        for k in range(4): v = v*A(Mb[k])**e[k]
        tot += v
    return tot

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--box', nargs=2, default=['1/2', '1'], help='H-O1 placeholder values M0o, M1o')
    ap.add_argument('--Mbar', nargs=4, default=['1', '1', '1', '1'], help='H-O2 placeholder values Mbar_0..Mbar_3')
    ap.add_argument('--delta-o', default='1/4', help='outer dissipation fraction (units a_min^o) for factors z_n')
    a = ap.parse_args(); t0 = time.time()
    M0o, M1o = Q_(a.box[0]), Q_(a.box[1]); Mb = [Q_(x) for x in a.Mbar]; DO = Q_(a.delta_o)
    T = WX.build_tables(); Ato = L.at_eval(M0o, M1o, Q_(162, 1000)); amin_o = (1 + A(M0o)**2)**-2
    assert (L.U_of_y(700) - A(Q_(162, 1000))).upper() < 0
    S = sups(); ln2 = LN2
    WR = L.load(L.OUT + '/weights_rep.json'); LAM = A(Q_(WR['lam'])); EPS4 = A(Q_(WR['eps4'])); amin_c = (1 + A(Q_(11, 100))**2)**-2
    Th350 = A(350)**-4*A(28)/350
    res = {}
    for j in (0, 1, 2):
        n = j + 2; P = {'a': {}, 'b': {}, 'c': {}, 'cY': {}, 'd': {}, 'e': {}}
        for (a_, b_), mons in T[j].items():
            if (a_, b_) == (n, n): continue
            for (m, i, r), co in mons.items():
                for X, val in Ato[(i, r)].items():
                    idx = [a_, b_] + list(X)
                    if max(idx) > 4 or idx.count(4) > 1 or (4 in idx and n != 4): raise ValueError(('outer index', j, idx))
                    # R_a: cutoff part (m >= 1)
                    if m >= 1:
                        if 4 in idx: raise ValueError(('z_4 in a cutoff term', j, idx))
                        addm(P['a'], mono(idx), ln2*A(co)*val*S[('rem', m)])
                    # R_c: base part, weight gamma^m * Theta
                    w = A(co)*val
                    if 4 in idx:
                        rest = [k for k in idx if k != 4]       # Young: w|z_rest||z_4| Theta <= delta_o a_min^o Theta z_4^2 + w^2 |z_rest|^2 Theta/(4 delta_o a_min^o)
                        addm(P['cY'], mono(rest + rest), ln2*S['Theta']*w*w/(4*A(DO)*amin_o))
                    else: addm(P['c'], mono(idx), ln2*S['Theta']*w)
        # rho on [350,700]: interpolation remainder for the (global) z_{j+1} interpolation with the outer box YC: use the term-by-term bound instead:
        # the outer base parts are already in R_c, so the interpolation is NOT used on [350,700]; its boundary terms at 350 (core YC) are R_b.
        YCc = {int(k): v.upper() for k, v in WR['core'][str(j)]['YC'].items()}
        if j + 1 <= 3:
            addm(P['b'], mono([j, j + 1]), YCc[j + 1]*Th350); addm(P['b'], mono([j, j]), YCc[j + 1]*A(KAP)*Th350/2)
        # R_d, R_e: L4 / S33 boundary terms at 350 (core coefficient weights w_c)
        l4w = {k: v.upper() for k, v in WR['core'][str(j)]['l4w'].items()}
        Ob = {(0, 1, 3, 0): Th350}                       # O_b = Theta(350) * Mbar_1 * Mbar_2^3   (monomial exponents in Mbar)
        for r, w in l4w.items():
            if r == 'S33':
                c1 = (w/2*KAP)**2/(4*LAM) + (w/2)**2/(4*EPS4*amin_c)
                for e, c in Ob.items(): addm(P['d'], e, c*c1)
                addm(P['e'], mono([3, 2, 2]), w*Th350/2)
            else:
                for e, c in Ob.items(): addm(P['d'], e, c*w*w/(4*LAM))
        tot = {}
        for part, Pp in P.items():
            for e, c in Pp.items(): addm(tot, e, c)
        val = evalP(tot, Mb); parts = {k: evalP(v, Mb) for k, v in P.items()}
        res[str(j)] = dict(poly=[dict(exp=list(e), coef=c) for e, c in sorted(tot.items())], rate_at_placeholder=val, parts=parts,
                           rate_at_x10=evalP(tot, [10*x for x in Mb]))
        print('level %d: outer rate per unit tau at the placeholder Mbar = %s: %.4g  (parts %s); at 10x Mbar: %.4g' % (j, a.Mbar, L.fup(val), {k: '%.2g' % L.fup(v) for k, v in parts.items()}, L.fup(res[str(j)]['rate_at_x10'])), flush=True)
    inputs = [L.BASE + '/rigor/weights_exact.py', L.BASE + '/rigor/tables_exact.pkl', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py', L.OUT + '/weights_rep.json']
    L.dump(dict(res=res, hyp=dict(M0o=str(M0o), M1o=str(M1o), Mbar=[str(x) for x in Mb], delta_o=str(DO), status='HYPOTHESES H-O1, H-O2: placeholder values, not certified'),
                sups={str(k): v for k, v in S.items()}), L.OUT + '/outer_rep.json', inputs)
    print('wrote out/outer_rep.json  %.0fs' % (time.time() - t0))
