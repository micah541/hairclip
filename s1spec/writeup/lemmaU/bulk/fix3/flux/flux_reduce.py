#!/usr/bin/env python3
"""F1 flux lemma, step 1 (exact, sympy): the outer-cutoff commutator forms K0, K2 of the verifier's E_I energies,
reduced by integration by parts on [57, 60] ONLY where the boundary coefficient vanishes identically at x = 57+ and x = 60-.

Objects (verifier step3nl_v33, ASSEMBLY.md sec. 0, 5, 6; CERT_PLAN.md:192; agent_sob/SOB_NOTES.md:43-45):
  g = w' (uncut slope deviation), v = zeta g, zeta = S7((60-x)/3) on [57,60] (=1 left of 57, =0 right of 60),
  S7(u) = u^4 (35 - 84u + 70u^2 - 20u^3)   (C^3 smoothstep; int_f23.py:14-18).
  L1 = -d^4 - k^2 d^2 - C d + e0,  C = 7km + q/4  (linear operator on v = w').
  Hc := -[L1, zeta] g = [d^4, zeta] g + k^2 [d^2, zeta] g + C zeta' g        (supported on [57,60]).
  K0 := int W v Hc dx                (contribution to (1/2) d/dtau V0^2)
  K2 := int (W v'')'' Hc dx          (contribution to (1/2) d/dtau V2^2, H^4-continuous form; see FLUX_REPORT.md sec. 2)
Each integrand is expanded into monomials  coef(x) * g^(a) g^(b)  (a >= b); coef is a product of W^(n), zeta^(i), C^(j), (k^2)^(j).
A monomial with a > target is integrated by parts only if coef vanishes at BOTH endpoints (exact check with the one-sided
values of zeta^(i) at 57+ and 60-; W, C, k^2 factors are treated as nonzero).  Otherwise it is kept and reported.
Output: flux_forms.pkl (sympy srepr of the reduced coefficient sums per pair, per level), flux_reduce.json (summary + args)."""
import sys, os, json, time, pickle, hashlib, platform
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol('x', real=True)
Wf = sp.Function('W')(x); Zf = sp.Function('Z')(x); Cf = sp.Function('C')(x); Kf = sp.Function('K')(x)   # K = k^2
NG = 8
G = sp.symbols('G0:%d' % NG)            # G_a stands for g^(a)
uu = sp.Symbol('u')
S7 = uu**4*(35 - 84*uu + 70*uu**2 - 20*uu**3)

def D(e, n=1): return sp.diff(e, x, n)

def zeta_end(i, end):
    """exact one-sided value of zeta^(i) at x = 57+ (u = 1) or 60- (u = 0), zeta(x) = S7((60-x)/3)."""
    u0 = 1 if end == 57 else 0
    return sp.Rational(-1, 3)**i*sp.diff(S7, uu, i).subs(uu, u0)

def Hc_of(gs):
    """Hc as a linear form in g^(j) (gs[j] = symbol for g^(j))"""
    z = lambda i: D(Zf, i) if i else Zf
    d4 = 4*z(1)*gs[3] + 6*z(2)*gs[2] + 4*z(3)*gs[1] + z(4)*gs[0]
    d2 = 2*z(1)*gs[1] + z(2)*gs[0]
    return d4 + Kf*d2 + Cf*z(1)*gs[0]

def dx_lin(expr):
    """x-derivative of an expression linear/bilinear in G (G_a -> G_{a+1})"""
    out = D(expr)
    for a in range(NG - 1):
        out += sp.diff(expr, G[a])*G[a + 1]
    return out

def integrand(level):
    v = Zf*G[0]
    if level == 0:
        return sp.expand(Wf*v*Hc_of(G))
    v2 = dx_lin(dx_lin(v))
    Wv2 = Wf*v2
    T = dx_lin(dx_lin(Wv2))
    return sp.expand(T*Hc_of(G))

def monomials(expr):
    """list of (a, b, coef) with a >= b, coef a product (no G)"""
    out = []
    for term in sp.Add.make_args(sp.expand(expr)):
        idx = []; coef = term
        for a in range(NG):
            p = sp.degree(term, G[a]) if term.has(G[a]) else 0
            idx += [a]*p
            if p: coef = coef/G[a]**p
        assert len(idx) == 2, (term, idx)
        b, a = sorted(idx)
        out.append((a, b, sp.expand(coef)))
    return out

_generic = None
def vanishes_at_ends(coef):
    """exact: coef (a sum of products) is identically zero at 57+ and at 60- (W, C, K and their derivatives arbitrary)."""
    for end in (57, 60):
        e = coef
        derivs = sorted(e.atoms(sp.Derivative), key=lambda d: -d.derivative_count)
        rep = {}
        for d in derivs:
            if d.expr == Zf: rep[d] = zeta_end(d.derivative_count, end)
        e = e.xreplace(rep).xreplace({Zf: zeta_end(0, end)})
        # remaining: W, C, K and derivatives -> independent symbols
        rep2 = {}; cnt = 0
        for d in sorted(e.atoms(sp.Derivative), key=str):
            rep2[d] = sp.Symbol('t%d' % cnt); cnt += 1
        for f in (Wf, Cf, Kf):
            rep2[f] = sp.Symbol('t%d' % cnt); cnt += 1
        if sp.expand(e.xreplace(rep2)) != 0: return False
    return True

def reduce_level(level, target):
    terms = {}
    def add(a, b, c):
        if b > a: a, b = b, a
        terms.setdefault((a, b), []).append(c)
    for a, b, c in monomials(integrand(level)): add(a, b, c)
    stuck = {}; log = []
    changed = True
    while changed:
        changed = False
        for key in sorted(terms, reverse=True):
            a, b = key
            if a <= target or a == b: continue
            lst = terms.pop(key); keep = []
            for c in lst:
                if not vanishes_at_ends(c):
                    keep.append(c); continue
                changed = True
                if a == b + 1:          # c g_{b+1} g_b = (1/2)(c g_b^2)' - (1/2) c' g_b^2
                    for t in sp.Add.make_args(sp.expand(-D(c)/2)): add(b, b, t)
                else:                   # c g_a g_b = (c g_{a-1} g_b)' - c' g_{a-1} g_b - c g_{a-1} g_{b+1}
                    for t in sp.Add.make_args(sp.expand(-D(c))): add(a - 1, b, t)
                    add(a - 1, b + 1, -c)
            if keep:
                terms[key] = terms.get(key, []) + keep if key in terms else keep
                stuck[key] = len(keep)
            if changed: break
    red = {k: sp.expand(sp.Add(*v)) for k, v in terms.items()}
    red = {k: v for k, v in red.items() if v != 0}
    # per-pair: does the summed coefficient vanish at the ends (diagnostic only)
    return red, {str(k): n for k, n in stuck.items() if k in red}

def is_td_check(level, red):
    """machine check: integrand - sum_red c_ab g_a g_b is a total derivative of a form whose coefficients vanish at 57+, 60-.
    Done constructively: we verify  integrand - reduced  ==  d/dx( boundary form ) by recomputing the boundary form along the reduction
    is not stored; instead use the Euler-operator test (exact): E[integrand - reduced] == 0 with g a generic function."""
    gfun = sp.Function('g')(x)
    def to_fun(e):
        for a in range(NG - 1, -1, -1): e = e.xreplace({G[a]: (D(gfun, a) if a else gfun)})
        return e
    diff = to_fun(integrand(level)) - to_fun(sum(c*G[a]*G[b] for (a, b), c in red.items()))
    from sympy.calculus.euler import euler_equations
    diff = sp.expand(diff)
    if diff == 0: return True
    eq = euler_equations(diff, [gfun], x)
    return (not eq) or sp.simplify(sp.expand(eq[0].lhs)) == 0

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

if __name__ == '__main__':
    t0 = time.time()
    res = {}; forms = {}
    for level, target in ((0, 2), (2, 2)):
        red, stuck = reduce_level(level, target)
        chk = is_td_check(level, red)
        print('level', level, 'pairs', sorted(red), 'stuck', stuck, 'Euler check', chk, flush=True)
        assert chk
        forms[level] = {'%d,%d' % k: sp.srepr(v) for k, v in red.items()}
        res[level] = dict(pairs=sorted('%d,%d' % k for k in red), stuck_monomials=stuck, euler_check=bool(chk),
                          ends_vanish={'%d,%d' % k: vanishes_at_ends(v) for k, v in red.items()})
    ze = {str(e): [str(zeta_end(i, e)) for i in range(8)] for e in (57, 60)}
    pickle.dump(forms, open(os.path.join(HERE, 'flux_forms.pkl'), 'wb'))
    args = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), sympy=sp.__version__,
                env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, inputs={os.path.abspath(__file__): sha(os.path.abspath(__file__))},
                runtime_s=round(time.time() - t0, 1))
    json.dump(dict(res=res, zeta_end_values=ze, outputs={'flux_forms.pkl': sha(os.path.join(HERE, 'flux_forms.pkl'))}, args=args),
              open(os.path.join(HERE, 'flux_reduce.json'), 'w'), indent=1)
    print('zeta^(i) at 57+:', ze['57']); print('zeta^(i) at 60-:', ze['60']); print('done %.1fs' % (time.time() - t0))
