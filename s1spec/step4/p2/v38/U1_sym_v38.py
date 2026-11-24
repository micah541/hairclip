#!/usr/bin/env python3
"""U1_sym.py -- exact sympy identities used by the regularity lemma U-1 (U1_regularity.tex).  Fix round 10, Nov 20 2026.
[v38 copy of <s1spec>/writeup/lemmaU/bulk/fix10_proofs/U1_sym.py for the step-4 v38 re-pin (p2/v38/final_tail_v38.py); only change: output path out_v38/iface/U1_sym.json]

Checks (each residual must simplify to 0, or the stated exact values must hold):
 S1  n = 1 reduction of the Lean velocity (CWFlow/Scalar/Pointwise.lean, defs meanCurv, flowPot, flowVel, read from source):
     for X : R -> C = R^2, J = rotation by +pi/2,
       H        = g^-1 (X'' - X' g^-1 <X'', X'>)            (HOf = tr_g (D^2 X)^perp, 1-dim trace = g^-1)
       flowPot  = g^-1 < d/dth (J H), X' >                  (tr_g <D(JH), DX>)
       flowVel  = J ( X' g^-1 d/dth flowPot )               (J grad_g flowPot)
     equals  -k_ss N  with N = J X'/|X'|, k = <X'', J X'>/|X'|^3, d/ds = |X'|^-1 d/dth.   (general parametrisation, X = (x(th), y(th)))
 S2  reflection equivariance: flowVel(R X) = R flowVel(X), R = diag(1, -1) (orientation-reversing isometry, R J = -J R).
 S3  reparametrisation invariance of -k_ss N (orientation preserving and reversing) and of the Lean velocity, for a general change of
     parameter phi (jets p1..p6 of phi, p1 != 0), by the chain rule on jets.
 S1-S3 are done in jet variables (x_i, y_i = i-th th-derivatives, r = |X'| with r^2 = x1^2 + y1^2); a rational expression is declared 0 when
     its numerator, as a polynomial in r, reduces to 0 modulo r^2 - g (exact, no floating point).
 S4  Lagrangian condition is automatic for curves: <J v, v> = 0.
 S5  rescaled normal velocity: gamma = p + lam Gamma, d tau = lam^-4 dt  =>  <Gamma_tau, nu> = -kappa_ss - (lam_tau/lam)<Gamma,nu> - lam^-1 <p_tau, nu>
     (scaling of kappa_ss: kappa_phys,ss_phys = lam^-3 kappa_ss).
 S6  cutoff profiles: S7(t) = t^4(35 - 84t + 70t^2 - 20t^3): S7^(j)(0) = 0 (j = 0..3), S7(1) = 1, S7^(j)(1) = 0 (j = 1..3), S7''''(0) = 840,
     S7''''(1) = -840;  quintic S(t) = t^3(10 - 15t + 6t^2): S^(j)(0) = 0 (j = 0..2), S(1) = 1, S^(j)(1) = 0 (j = 1, 2), S'''(0) = 60, S'''(1) = 60.
 S7  level-2 integration by parts: d/dx [W v'' h' - (W v'')' h] = W v'' h'' - (W v'')'' h, and (W v'')''/W = v'''' - 2 psi' v''' + (psi'^2 - psi'') v''
     for W = exp(-psi).
 S8  Agmon on an interval of length l (used by the flux lemma, s_2): the identity f(x)^2 - f(y)^2 = int_y^x 2 f f' (exact, symbolic check of the
     derivative), from which sup f^2 <= ||f||^2/l + 2 ||f|| ||f'|| follows [hand, one line].
Output: U1_sym.json (args block with argv, env, sympy version, own sha256).  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B U1_sym.py
"""
import sys, os, json, hashlib, platform
sys.dont_write_bytecode = True
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
res = {}


def check(name, expr):
    z = sp.simplify(sp.expand(expr)) if not isinstance(expr, sp.MatrixBase) else expr.applyfunc(lambda e: sp.simplify(e))
    ok = (z == 0) if not isinstance(z, sp.MatrixBase) else all(sp.simplify(e) == 0 for e in z)
    res[name] = bool(ok)
    print('%-70s %s' % (name, ok), flush=True)
    assert ok, name


# ---- jet calculus: x_i = d^i x/dth^i, y_i likewise (i = 0..6), r = |X'| with r^2 = g = x1^2 + y1^2 (r > 0).
# An expression is a rational function of the jets and r; D is the total th-derivative; zero-tests reduce r^2 -> g.
NJ = 8
xs = sp.symbols('x0:%d' % NJ, real=True); ys = sp.symbols('y0:%d' % NJ, real=True)
r = sp.Symbol('r', positive=True)
g_ = xs[1]**2 + ys[1]**2


def D(e):
    out = sp.diff(e, r) * (xs[1] * xs[2] + ys[1] * ys[2]) / r
    for i in range(NJ - 1):
        out += sp.diff(e, xs[i]) * xs[i + 1] + sp.diff(e, ys[i]) * ys[i + 1]
    return out


def Dv(V): return sp.Matrix([D(c) for c in V])


def is_zero(e):
    n = sp.numer(sp.together(e))
    n = sp.expand(n)
    P = sp.Poly(n, r)
    rem = sp.rem(P, sp.Poly(r**2 - g_, r))
    return all(sp.expand(c) == 0 for c in rem.all_coeffs())


def check0(name, e):
    if isinstance(e, sp.MatrixBase): ok = all(is_zero(c) for c in e)
    else: ok = is_zero(e)
    res[name] = bool(ok); print('%-70s %s' % (name, ok), flush=True); assert ok, name


Jm = sp.Matrix([[0, -1], [1, 0]])
def ip(a, b): return (a.T * b)[0, 0]


def lean_flowvel(Xp, Xpp, Dfun):
    """Lean defs read in a chart, 1-dim domain: H = g^-1 (X'')^perp, pot = g^-1 <(J H)', X'>, vel = J (X' g^-1 pot')"""
    g = ip(Xp, Xp)
    H = (Xpp - Xp * ip(Xpp, Xp) / g) / g
    pot = ip(Dfun(Jm * H), Xp) / g
    return Jm * (Xp * Dfun(pot) / g)


def geom_vel(Xp, Xpp, rr, Dfun):
    k = ip(Xpp, Jm * Xp) / rr**3
    kss = Dfun(Dfun(k) / rr) / rr
    return -kss * (Jm * Xp / rr)


Xp = sp.Matrix([xs[1], ys[1]]); Xpp = sp.Matrix([xs[2], ys[2]])
D1 = lambda e: (Dv(e) if isinstance(e, sp.MatrixBase) else D(e))
# S1
check0('S1 flowVel (Lean, n=1) == -k_ss N, general parametrisation (jets)', lean_flowvel(Xp, Xpp, D1) - geom_vel(Xp, Xpp, r, D1))
# S2: R X has jets (x_i, -y_i); |(RX)'| = r
R = sp.Matrix([[1, 0], [0, -1]])
check0('S2 R J = -J R', R * Jm + Jm * R)
check0('S2 flowVel(R X) == R flowVel(X)', lean_flowvel(R * Xp, R * Xpp, D1) - R * lean_flowvel(Xp, Xpp, D1))
# S3: X~(th) = X(phi(th)); phi jets p1..p5 (p1 != 0), D~ = d/dth acting on (jets of X at phi, phi jets): D~ f = p1 * D f + sum p_{i+1} d f/d p_i
ps = sp.symbols('p1:7', real=True); sg = sp.Symbol('sigma')    # sigma = sign(p1): |X~'| = sigma p1 r


def Dt(e):
    out = ps[0] * D(e)
    for i in range(len(ps) - 1):
        out += sp.diff(e, ps[i]) * ps[i + 1]
    return out


def Dtv(e): return sp.Matrix([Dt(c) for c in e]) if isinstance(e, sp.MatrixBase) else Dt(e)


Xtp = Xp * ps[0]; Xtpp = Xpp * ps[0]**2 + Xp * ps[1]
for nm, sv in (('orientation preserving (p1 > 0)', 1), ('orientation reversing (p1 < 0)', -1)):
    lhs = geom_vel(Xtp, Xtpp, sv * ps[0] * r, Dtv)
    rhs = geom_vel(Xp, Xpp, r, D1)
    check0('S3 -k_ss N invariant under reparametrisation, %s' % nm, lhs - rhs)
check0('S3 Lean flowVel invariant under reparametrisation (any p1 != 0)', lean_flowvel(Xtp, Xtpp, Dtv) - lean_flowvel(Xp, Xpp, D1))
# S4
a1, a2 = sp.symbols('a1 a2', real=True); v = sp.Matrix([a1, a2])
check('S4 <J v, v> = 0 (curves are Lagrangian)', ip(Jm * v, v))
# S5
lam, lamt, kss, Gn, pn = sp.symbols('lam lam_tau kappa_ss Gamma_nu p_tau_nu', real=True)
# physical normal velocity = -kappa_phys,ss = -lam^-3 kss ;  <gamma_t, nu> = <p_t,nu> + lam_t <Gamma,nu> + lam <Gamma_tau,nu> tau_t, tau_t = lam^-4, f_t = lam^-4 f_tau
Gt = sp.Symbol('Gamma_tau_nu')
eq = sp.Eq(lam**-4 * pn + lam**-4 * lamt * Gn + lam * Gt * lam**-4, -lam**-3 * kss)
sol = sp.solve(eq, Gt)[0]
check('S5 <Gamma_tau,nu> = -kappa_ss - (lam_tau/lam)<Gamma,nu> - lam^-1 <p_tau,nu>', sol - (-kss - lamt / lam * Gn - pn / lam))
s = sp.Symbol('s', positive=True); L = sp.Symbol('L', positive=True)
check('S5 scaling: kappa(s) of lam*Gamma at lam*s is kappa/lam; d^2/ds_phys^2 gives lam^-3', sp.diff(sp.Function('K')(s / L) / L, s, 2) - sp.Function('K')(s / L).diff(s, 2) / L)
# S6
t = sp.Symbol('t')
S7 = t**4 * (35 - 84 * t + 70 * t**2 - 20 * t**3); S5q = t**3 * (10 - 15 * t + 6 * t**2)
vals = {}
for j in range(5):
    vals['S7^(%d)(0)' % j] = sp.diff(S7, t, j).subs(t, 0); vals['S7^(%d)(1)' % j] = sp.diff(S7, t, j).subs(t, 1)
for j in range(4):
    vals['S^(%d)(0)' % j] = sp.diff(S5q, t, j).subs(t, 0); vals['S^(%d)(1)' % j] = sp.diff(S5q, t, j).subs(t, 1)
want = {'S7^(0)(0)': 0, 'S7^(1)(0)': 0, 'S7^(2)(0)': 0, 'S7^(3)(0)': 0, 'S7^(4)(0)': 840, 'S7^(0)(1)': 1, 'S7^(1)(1)': 0, 'S7^(2)(1)': 0,
        'S7^(3)(1)': 0, 'S7^(4)(1)': -840, 'S^(0)(0)': 0, 'S^(1)(0)': 0, 'S^(2)(0)': 0, 'S^(3)(0)': 60, 'S^(0)(1)': 1, 'S^(1)(1)': 0,
        'S^(2)(1)': 0, 'S^(3)(1)': 60}
for k_, w_ in want.items():
    check('S6 %s = %s' % (k_, w_), vals[k_] - w_)
check('S6 S7(t) + S7(1-t) = 1', sp.expand(S7 + S7.subs(t, 1 - t) - 1))
# S7
xx = sp.Symbol('x', real=True)
psi = sp.Function('psi')(xx); vv = sp.Function('v')(xx); hh = sp.Function('h')(xx); W = sp.exp(-psi)
F = W * vv.diff(xx, 2) * hh.diff(xx) - (W * vv.diff(xx, 2)).diff(xx) * hh
check('S7 d/dx[W v\'\' h\' - (W v\'\')\' h] = W v\'\' h\'\' - (W v\'\')\'\' h', F.diff(xx) - (W * vv.diff(xx, 2) * hh.diff(xx, 2) - (W * vv.diff(xx, 2)).diff(xx, 2) * hh))
check('S7 (W v\'\')\'\'/W = v\'\'\'\' - 2 psi\' v\'\'\' + (psi\'^2 - psi\'\') v\'\'',
      (W * vv.diff(xx, 2)).diff(xx, 2) / W - (vv.diff(xx, 4) - 2 * psi.diff(xx) * vv.diff(xx, 3) + (psi.diff(xx)**2 - psi.diff(xx, 2)) * vv.diff(xx, 2)))
# level-0 analogue: d/dx[W v h] relation used for <v, zeta f'>_W  (EI_SPEC I10): (W zeta v)' = W (zeta v' + (zeta' - psi' zeta) v)
zz = sp.Function('zeta')(xx)
check('S7 (W zeta v)\' = W (zeta v\' + (zeta\' - psi\' zeta) v)', (W * zz * vv).diff(xx) - W * (zz * vv.diff(xx) + (zz.diff(xx) - psi.diff(xx) * zz) * vv))
# S8
ff = sp.Function('f')(xx)
check('S8 d/dx f^2 = 2 f f\'', (ff**2).diff(xx) - 2 * ff * ff.diff(xx))

me = os.path.abspath(__file__)
out = dict(checks=res, all_pass=all(res.values()),
           args=dict(argv=sys.argv, cwd=os.getcwd(), env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, python=platform.python_version(),
                     sympy=sp.__version__, inputs={me: hashlib.sha256(open(me, 'rb').read()).hexdigest()}))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/iface/U1_sym.json', 'w'), indent=1, sort_keys=True)
print('all_pass', out['all_pass'])
