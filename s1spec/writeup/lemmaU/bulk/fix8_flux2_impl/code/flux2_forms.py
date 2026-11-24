#!/usr/bin/env python3
"""Stage X1 (FLUX2_DESIGN.md sec. 10): exact sympy forms for flux lemma v2.  [sym] throughout; no numbers.

(a) The outer-ramp identity of sec. 2.1 from the v-equation (EI_SPEC sec. 2.3), with generic functions g, zeta (Z), j (J), k^2 (K),
    C = 7km + q/4, e0, and the modulation coefficient m_T9 := a Q_u + b C_u:
      zeta (L1 g + (c g''')') - [L1 v + (c v''')'] - zeta' m_T9 g  ==  J K4g + J' K3g + K K2g + C_eff zeta' g,
    with c = 1 - J, C_eff = C - m_T9, L1 f = -f'''' - K f'' - C f' + e0 f, v = zeta g;  plus SPEC I2 and the j-commutator form.
(b) Level-0 reduction of  W zeta g (J K4g + J' K3g + K K2g + C_eff zeta' g)  by F1's algorithm (fix3/flux/flux_reduce.py, reduce_level(0, 2)):
    a monomial c g^(a) g^(b), a > 2, is integrated by parts only if c vanishes identically at 57+ and 60- (exact one-sided zeta^(i) values;
    W, C_eff, K, J and their derivatives generic).  Asserted: Euler-operator check; NO J'' (and no derivative of C_eff, K) in the result;
    the boundary coefficient of every integration by parts vanishes at both ends; only the monomial 4 W zeta zeta' J g g''' is moved.
    Regression: J = 1, J' = 0 reproduces fix3/flux/flux_forms.pkl level 0 exactly (pair by pair, sympy expand of the difference == 0).
    Hand form of sec. 4 item 1 checked as well (Euler).
(c) D2v = (W v'')''/W coefficients (sec. 5): zeta g'''' + (4Z' - 2 psi' Z) g''' + (6Z'' - 6psi' Z' + psi2 Z) g'' + ... , psi2 = psi'^2 - psi''.
    Also (W zeta v)' = W (zeta^2 g' + (2 zeta zeta' - psi' zeta^2) g)  (level-0 nonlinear term, sec. 4 item 2).
(d) N_low in g-jets from step3nl_v4b/graph_struct.pkl: F is a polynomial in (u'', u''', u'''') with exactly the 7 monomials of EI_SPEC 2.2;
    substituting u'' = a2 + g', u''' = a3 + g'' gives N_low = d010 g'' + d100 g' + d200 g'^2 + d110 g'g'' + d300 g'^3 + d000 (exactly these
    monomials; d110 = G110 = 10 Ji^6 u1, d300 = G300, d200 = G200 + 3 a2 G300, d010 = G010 + a2 G110); and N_low' = the 12-monomial list of
    EI_SPEC 4.5 (generic d's as functions of x).
(e) The conversions sigma^3 d^3 = D^3 - 3D^2 + 2D, sigma^4 d^4 = D^4 - 6D^3 + 11D^2 - 6D (D = sigma d/dsigma), and Cauchy-Schwarz weights.
(f) Total-derivative residue of the trace route (sec. 6.2): (g'g'')' g'''' + (5/2) g'' g'''^2 is a total derivative (and 0, 2 are not).
Outputs: out/flux2_forms.pkl (srepr of the reduced level-0 forms, D2v coefficients), out/flux2_forms.json (summary, assertions, args)."""
import sys, os, json, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
import sympy as sp
from sympy.calculus.euler import euler_equations

T = FC.Timer()
x = sp.Symbol('x', real=True)
Wf = sp.Function('W')(x); Zf = sp.Function('Z')(x); Cf = sp.Function('C')(x); Kf = sp.Function('K')(x); Jf = sp.Function('J')(x)
NG = 8
G = sp.symbols('G0:%d' % NG)
uu = sp.Symbol('u')
S7 = uu**4*(35 - 84*uu + 70*uu**2 - 20*uu**3)
CHECKS = {}


def ok(name, val):
    CHECKS[name] = bool(val)
    print('%-70s %s' % (name, bool(val)), flush=True)
    assert val, name


def D(e, n=1): return sp.diff(e, x, n)


def zeta_end(i, end):
    u0 = 1 if end == 57 else 0
    return sp.Rational(-1, 3)**i*sp.diff(S7, uu, i).subs(uu, u0)


def td(expr, gfun):
    """exact: expr (polynomial in the jets of gfun) is a total x-derivative  <=>  Euler operator vanishes"""
    expr = sp.expand(expr)
    if expr == 0:
        return True
    eq = euler_equations(expr, [gfun], x)
    return (not eq) or eq[0] == True or sp.simplify(sp.expand(eq[0].lhs)) == 0


# ======================================================================= (a) identities from the v-equation
gF = sp.Function('g')(x); e0F = sp.Function('e0')(x); mT = sp.Function('mT9')(x)
Zd = lambda i: D(Zf, i) if i else Zf
gd = lambda i: D(gF, i) if i else gF
K2g = 2*Zd(1)*gd(1) + Zd(2)*gd(0)
K3g = 3*Zd(1)*gd(2) + 3*Zd(2)*gd(1) + Zd(3)*gd(0)
K4g = 4*Zd(1)*gd(3) + 6*Zd(2)*gd(2) + 4*Zd(3)*gd(1) + Zd(4)*gd(0)
L1 = lambda f: -D(f, 4) - Kf*D(f, 2) - Cf*D(f) + e0F*f
cF = 1 - Jf
v = Zf*gF
lhs = Zf*(L1(gF) + D(cF*D(gF, 3))) - (L1(v) + D(cF*D(v, 3))) - Zd(1)*mT*gF
rhs = Jf*K4g + D(Jf)*K3g + Kf*K2g + (Cf - mT)*Zd(1)*gF
ok('(a1) zeta(L1 g + (c g\'\'\')\') - L1 v - (c v\'\'\')\' - zeta\' mT9 g == j K4 + j\' K3 + k^2 K2 + C_eff zeta\' g', sp.expand(lhs - rhs) == 0)
ok('(a2) SPEC I2: zeta (c g\'\'\')\' == (c v\'\'\')\' - c K4 - c\' K3', sp.expand(Zf*D(cF*D(gF, 3)) - (D(cF*D(v, 3)) - cF*K4g - D(cF)*K3g)) == 0)
ok('(a3) zeta (-(j g\'\'\')\') == -(j v\'\'\')\' + j K4 + j\' K3', sp.expand(-Zf*D(Jf*D(gF, 3)) - (-D(Jf*D(v, 3)) + Jf*K4g + D(Jf)*K3g)) == 0)
ok('(a4) -[L1, zeta] g == K4 + k^2 K2 + C zeta\' g', sp.expand(-(L1(v) - Zf*L1(gF)) - (K4g + Kf*K2g + Cf*Zd(1)*gF)) == 0)
ok('(a5) -d^4 + (c d^3)\' == -(j d^3)\' on g', sp.expand(-D(gF, 4) + D(cF*D(gF, 3)) + D(Jf*D(gF, 3))) == 0)


# ======================================================================= (b) level-0 reduction with generic j
def Hc_of(gs):
    z = lambda i: D(Zf, i) if i else Zf
    k4 = 4*z(1)*gs[3] + 6*z(2)*gs[2] + 4*z(3)*gs[1] + z(4)*gs[0]
    k3 = 3*z(1)*gs[2] + 3*z(2)*gs[1] + z(3)*gs[0]
    k2 = 2*z(1)*gs[1] + z(2)*gs[0]
    return Jf*k4 + D(Jf)*k3 + Kf*k2 + Cf*z(1)*gs[0]          # Cf stands for C_eff here


def integrand0():
    return sp.expand(Wf*Zf*G[0]*Hc_of(G))


def monomials(expr):
    out = []
    for term in sp.Add.make_args(sp.expand(expr)):
        idx = []; coef = term
        for a in range(NG):
            p = sp.degree(term, G[a]) if term.has(G[a]) else 0
            idx += [a]*p
            if p:
                coef = coef/G[a]**p
        assert len(idx) == 2, (term, idx)
        b, a = sorted(idx)
        out.append((a, b, sp.expand(coef)))
    return out


def vanishes_at_ends(coef):
    for end in (57, 60):
        e = coef
        rep = {}
        for d in sorted(e.atoms(sp.Derivative), key=lambda d: -d.derivative_count):
            if d.expr == Zf:
                rep[d] = zeta_end(d.derivative_count, end)
        e = e.xreplace(rep).xreplace({Zf: zeta_end(0, end)})
        rep2 = {}; cnt = 0
        for d in sorted(e.atoms(sp.Derivative), key=str):
            rep2[d] = sp.Symbol('t%d' % cnt); cnt += 1
        for f in (Wf, Cf, Kf, Jf):
            rep2[f] = sp.Symbol('t%d' % cnt); cnt += 1
        if sp.expand(e.xreplace(rep2)) != 0:
            return False
    return True


MOVED = []


def reduce_level0(target=2):
    terms = {}

    def add(a, b, c):
        if b > a:
            a, b = b, a
        terms.setdefault((a, b), []).append(c)
    for a, b, c in monomials(integrand0()):
        add(a, b, c)
    stuck = {}
    changed = True
    while changed:
        changed = False
        for key in sorted(terms, reverse=True):
            a, b = key
            if a <= target or a == b:
                continue
            lst = terms.pop(key); keep = []
            for c in lst:
                if not vanishes_at_ends(c):
                    keep.append(c); continue
                changed = True
                MOVED.append(((a, b), sp.srepr(c)))
                if a == b + 1:
                    for t in sp.Add.make_args(sp.expand(-D(c)/2)):
                        add(b, b, t)
                else:
                    for t in sp.Add.make_args(sp.expand(-D(c))):
                        add(a - 1, b, t)
                    add(a - 1, b + 1, -c)
            if keep:
                terms[key] = keep
                stuck[key] = len(keep)
            if changed:
                break
    red = {k: sp.expand(sp.Add(*v)) for k, v in terms.items()}
    red = {k: v for k, v in red.items() if v != 0}
    return red, {str(k): n for k, n in stuck.items() if k in red}


def to_fun(e, gfun):
    for a in range(NG - 1, -1, -1):
        e = e.xreplace({G[a]: (D(gfun, a) if a else gfun)})
    return e


red0, stuck0 = reduce_level0(2)
print('level-0 reduced pairs', sorted(red0), 'stuck', stuck0, flush=True)
gfun = sp.Function('g')(x)
ok('(b1) Euler: integrand - reduced level-0 form is a total derivative', td(to_fun(integrand0(), gfun) - to_fun(sum(c*G[a]*G[b] for (a, b), c in red0.items()), gfun), gfun))
ok('(b2) no stuck monomial (every pair has a <= 2)', not stuck0 and all(a <= 2 for (a, b) in red0))
derivs = set()
for c in red0.values():
    for d in c.atoms(sp.Derivative):
        derivs.add((str(d.expr), d.derivative_count))
ok('(b3) no j\'\' in the level-0 form (J derivative order <= 1)', all(n <= 1 for f, n in derivs if f == 'J(x)'))
ok('(b4) no derivative of C_eff or of k^2 in the level-0 form', not any(f in ('C(x)', 'K(x)') for f, n in derivs))
ok('(b5) exactly one monomial moved, and it is 4 W zeta zeta\' J g g\'\'\'', len(MOVED) == 1 and MOVED[0][0] == (3, 0)
   and sp.expand(sp.sympify(eval(MOVED[0][1], sp.__dict__.copy() | {'x': x})) - 4*Wf*Zf*D(Zf)*Jf) == 0)
ok('(b6) the moved monomial\'s coefficient vanishes at 57+ and 60-', vanishes_at_ends(4*Wf*Zf*D(Zf)*Jf))
# affine in (J, J'): each coefficient is linear in {J, J'} jointly (no J*J', no J^2)
Js, Jps = sp.symbols('Jsym Jpsym')
aff = True
for c in red0.values():
    cs = c.xreplace({D(Jf): Jps}).xreplace({Jf: Js})
    P = sp.Poly(cs, Js, Jps)
    aff &= P.total_degree() <= 1
ok('(b7) reduced coefficients are affine (degree <= 1) in (j, j\')', aff)
# regression against F1 level 0
f1forms = pickle.load(open(os.path.join(FC.F1, 'flux_forms.pkl'), 'rb'))[0]
ns = sp.__dict__.copy(); ns['x'] = x
regr = {}
for key, srep in f1forms.items():
    a, b = map(int, key.split(','))
    mine = red0.get((a, b), sp.Integer(0))
    mine1 = sp.expand(mine.xreplace({D(Jf): 0}).xreplace({Jf: 1}))
    regr[key] = sp.expand(mine1 - eval(srep, ns)) == 0
extra = [k for k in red0 if '%d,%d' % k not in f1forms and sp.expand(red0[k].xreplace({D(Jf): 0}).xreplace({Jf: 1})) != 0]
ok('(b8) regression: j = 1, j\' = 0 reproduces F1 flux_forms.pkl level 0 pair by pair', all(regr.values()) and not extra)
# hand form of sec. 4 item 1: kappa = 4 W Z Z' J;  int kappa g g''' = -int kappa' g g'' + 1/2 int kappa' g'^2
kap = 4*Wf*Zf*D(Zf)*Jf
ok('(b9) hand reduction int kappa g g\'\'\' = -int kappa\' g g\'\' + 1/2 int kappa\' g\'^2 (Euler)',
   td(kap*D(gfun, 3)*gfun - (-D(kap)*gfun*D(gfun, 2) + sp.Rational(1, 2)*D(kap)*D(gfun)**2), gfun))
# the hand pair list (sec. 4 item 1) vs the reduced form: equal modulo total derivatives
hand = ((-D(kap) + 6*Wf*Zf*D(Zf, 2)*Jf + 3*Wf*Zf*D(Zf)*D(Jf))*G[2]*G[0] + sp.Rational(1, 2)*D(kap)*G[1]**2
        + (4*Wf*Zf*D(Zf, 3)*Jf + 3*Wf*Zf*D(Zf, 2)*D(Jf) + 2*Kf*Wf*Zf*D(Zf))*G[1]*G[0]
        + (Wf*Zf*D(Zf, 4)*Jf + Wf*Zf*D(Zf, 3)*D(Jf) + Kf*Wf*Zf*D(Zf, 2) + Wf*Cf*Zf*D(Zf))*G[0]**2)
ok('(b10) hand pair list of sec. 4 item 1 == reduced form modulo total derivatives',
   td(to_fun(hand, gfun) - to_fun(sum(c*G[a]*G[b] for (a, b), c in red0.items()), gfun), gfun))
# boundary-vanishing diagnostic of the reduced pairs (as F1)
ends0 = {'%d,%d' % k: vanishes_at_ends(v) for k, v in red0.items()}

# ======================================================================= (c) D2v and (W zeta v)'
psi = sp.Function('psi')(x)
Wexp = sp.exp(-psi); vv = Zf*gfun
D2 = sp.expand(sp.simplify((Wexp*D(vv, 2)).diff(x, 2)/Wexp))
p1 = D(psi); p2 = p1**2 - D(psi, 2)
D2C = {4: Zf, 3: 4*Zd(1) - 2*p1*Zf, 2: 6*Zd(2) - 6*p1*Zd(1) + p2*Zf, 1: 4*Zd(3) - 6*p1*Zd(2) + 2*p2*Zd(1), 0: Zd(4) - 2*p1*Zd(3) + p2*Zd(2)}
ok('(c1) D2v = (W v\'\')\'\'/W coefficients of g^(i), i = 0..4 (sec. 5)', sp.simplify(sp.expand(D2 - sum(c*D(gfun, i) if i else c*gfun for i, c in D2C.items()))) == 0)
ok('(c2) D2v = v\'\'\'\' - 2 psi\' v\'\'\' + psi2 v\'\'', sp.simplify(sp.expand(D2 - (D(vv, 4) - 2*p1*D(vv, 3) + p2*D(vv, 2)))) == 0)
WZv = sp.simplify(D(Wexp*Zf*vv)/Wexp)
ok('(c3) (W zeta v)\'/W = zeta^2 g\' + (2 zeta zeta\' - psi\' zeta^2) g', sp.simplify(sp.expand(WZv - (Zf**2*D(gfun) + (2*Zf*Zd(1) - p1*Zf**2)*gfun))) == 0)
ok('(c4) level-0 IBP: W v zeta f\' + (W zeta v)\' f is a total derivative', sp.simplify(sp.expand(Wexp*vv*Zf*D(sp.Function('f')(x)) + D(Wexp*Zf*vv)*sp.Function('f')(x) - D(Wexp*Zf*vv*sp.Function('f')(x)))) == 0)

# ======================================================================= (d) N_low from graph_struct.pkl
GS = pickle.load(open(os.path.join(FC.V4B, 'graph_struct.pkl'), 'rb'))
keys = sorted(GS.keys())
ok('(d1) graph_struct monomials = {(0,0,0),(0,0,1),(0,1,0),(1,0,0),(1,1,0),(2,0,0),(3,0,0)}', set(keys) == {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 0), (2, 0, 0), (3, 0, 0)})
Ji, u1 = sp.symbols('Ji u1')
ok('(d2) u\'\'\'\' coefficient = -Ji^4 (so c = 1 - J^-4)', sp.expand(GS[(0, 0, 1)] + Ji**4) == 0)
ok('(d3) d110 = G110 = 10 Ji^6 u1', sp.expand(GS[(1, 1, 0)] - 10*Ji**6*u1) == 0)
ok('(d4) d300 = G300 = 3 Ji^6 (1 - 6 Ji^2 u1^2)', sp.expand(GS[(3, 0, 0)] - 3*Ji**6*(1 - 6*Ji**2*u1**2)) == 0)
# the coefficients depend on (U0, U1) and profile/alpha symbols only, never on u'', u''', u'''' (polynomial structure)
a2s, a3s, a4s, g1s, g2s, g3s = sp.symbols('a2 a3 a4 g1 g2 g3')
Gc = {k: sp.Symbol('G%d%d%d' % k) for k in keys}
Fpoly = sum(Gc[k]*(a2s + g1s)**k[0]*(a3s + g2s)**k[1]*(a4s + g3s)**k[2] for k in keys)
Npoly = sp.expand(Fpoly - Fpoly.subs({g1s: 0, g2s: 0, g3s: 0}))
Pn = sp.Poly(Npoly, g1s, g2s, g3s)
mons = set(Pn.monoms())
ok('(d5) N(w) - lin. part in g-jets: monomials of (g\', g\'\', g\'\'\') = {g\'\'\', g\'\', g\', g\'^2, g\'g\'\', g\'^3}', mons == {(0, 0, 1), (0, 1, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (3, 0, 0)})
dl = {m: Pn.coeff_monomial(m) for m in mons}
ok('(d6) d001 (c) = G001; d110 = G110; d300 = G300; d200 = G200 + 3 a2 G300; d010 = G010 + a2 G110',
   sp.expand(dl[(0, 0, 1)] - Gc[(0, 0, 1)]) == 0 and sp.expand(dl[(1, 1, 0)] - Gc[(1, 1, 0)]) == 0 and sp.expand(dl[(3, 0, 0)] - Gc[(3, 0, 0)]) == 0
   and sp.expand(dl[(2, 0, 0)] - Gc[(2, 0, 0)] - 3*a2s*Gc[(3, 0, 0)]) == 0 and sp.expand(dl[(0, 1, 0)] - Gc[(0, 1, 0)] - a2s*Gc[(1, 1, 0)]) == 0)
# N_low' = 12-monomial list (generic d's as functions of x)
dn = {n: sp.Function(n)(x) for n in ('d010', 'd100', 'd200', 'd110', 'd300', 'd000')}
g1, g2, g3 = D(gfun), D(gfun, 2), D(gfun, 3)
Nlow = dn['d010']*g2 + dn['d100']*g1 + dn['d200']*g1**2 + dn['d110']*g1*g2 + dn['d300']*g1**3 + dn['d000']
Dd = {n: D(f) for n, f in dn.items()}
LIST = (dn['d010']*g3 + (Dd['d010'] + dn['d100'])*g2 + Dd['d100']*g1 + 2*dn['d200']*g1*g2 + Dd['d200']*g1**2
        + dn['d110']*(g2**2 + g1*g3) + Dd['d110']*g1*g2 + 3*dn['d300']*g1**2*g2 + Dd['d300']*g1**3 + Dd['d000'])
ok('(d7) N_low\' == the twelve-monomial list of EI_SPEC 4.5 / FLUX2_DESIGN sec. 5 table', sp.expand(D(Nlow) - LIST) == 0)
# chain rule behind Dd = d^a + d^b g' (d depends on x, U0, U1; U0' = U1, U1' = a2 + g'): linear in g' with coefficient d_U1
xs, U0s, U1s = sp.symbols('xs U0s U1s'); dfun = sp.Function('dd')(xs, U0s, U1s)
U0f = sp.Function('U0')(x); U1f = sp.Function('U1')(x); a2f = sp.Function('a2')(x)
dcomp = sp.Function('dd')(x, U0f, U1f)
tot = D(dcomp).subs({D(U0f): U1f, D(U1f): a2f + g1})
ok('(d8) D d = d_x + d_U0 U1 + d_U1 (a2 + g\'): affine in g\' with slope d_U1', sp.Poly(sp.expand(tot), g1).degree() == 1)

# ======================================================================= (e) conversions
s = sp.Symbol('s', positive=True); f = sp.Function('f')(s); Dop = lambda e: s*sp.diff(e, s)
D1 = Dop(f); D2_ = Dop(D1); D3 = Dop(D2_); D4 = Dop(D3)
ok('(e1) sigma^3 d^3 = D^3 - 3D^2 + 2D', sp.simplify(s**3*f.diff(s, 3) - (D3 - 3*D2_ + 2*D1)) == 0)
ok('(e2) sigma^4 d^4 = D^4 - 6D^3 + 11D^2 - 6D', sp.simplify(s**4*f.diff(s, 4) - (D4 - 6*D3 + 11*D2_ - 6*D1)) == 0)
ok('(e3) sigma^2 d^2 = D^2 - D; sigma d = D', sp.simplify(s**2*f.diff(s, 2) - (D2_ - D1)) == 0)
ok('(e4) D^4 = sigma d + 7 sigma^2 d^2 + 6 sigma^3 d^3 + sigma^4 d^4 (Stirling, used for D^4 d)',
   sp.simplify(D4 - (s*f.diff(s) + 7*s**2*f.diff(s, 2) + 6*s**3*f.diff(s, 3) + s**4*f.diff(s, 4))) == 0)

# ======================================================================= (f) trace-route residue
G_ = [gfun] + [D(gfun, i) for i in range(1, 6)]
ex = D(G_[1]*G_[2])*G_[4]
ok('(f1) (g\'g\'\')\' g\'\'\'\' + (5/2) g\'\' g\'\'\'^2 is a total derivative', td(ex + sp.Rational(5, 2)*G_[2]*G_[3]**2, gfun))
ok('(f2) ... with 0 in place of 5/2: not a total derivative', not td(ex, gfun))
ok('(f3) ... with 2 in place of 5/2: not a total derivative', not td(ex + 2*G_[2]*G_[3]**2, gfun))

# ======================================================================= outputs
forms = dict(level0={'%d,%d' % k: sp.srepr(v) for k, v in sorted(red0.items())},
             D2v={str(i): sp.srepr(c) for i, c in D2C.items()},
             symbols=dict(W='W(x)', zeta='Z(x)', C_eff='C(x)', k2='K(x)', j='J(x)', G='G0..G7 = g^(a)', psi='psi(x)'))
pk = os.path.join(FC.OUT, 'flux2_forms.pkl')
tmp = pk + '.tmp'
with open(tmp, 'wb') as fh:
    pickle.dump(forms, fh, protocol=4)
os.replace(tmp, pk)
res = dict(stage='X1 flux2_forms', label='sym (exact sympy identities; no numbers)', checks=CHECKS, all_pass=all(CHECKS.values()),
           level0_pairs=sorted('%d,%d' % k for k in red0), level0_stuck=stuck0, level0_ends_vanish=ends0,
           level0_moved=[dict(pair='%d,%d' % k, coef=c) for k, c in MOVED],
           regression_F1_level0=regr, outputs={'flux2_forms.pkl': FC.sha(pk)},
           zeta_end_values={str(e): [str(zeta_end(i, e)) for i in range(8)] for e in (57, 60)},
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), os.path.join(FC.F1, 'flux_forms.pkl'),
                               os.path.join(FC.V4B, 'graph_struct.pkl')]))
res['timing'] = dict(runtime_s=T.s())
h = FC.dump(os.path.join(FC.OUT, 'flux2_forms.json'), res)
print('X1 done: all checks pass; pairs', res['level0_pairs'], 'sha', h, '(%.1fs)' % T.s())
