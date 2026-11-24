"""Expression families from ibp_sym, compiled to python-flint arb functions (integer powers -> products; cse)."""
import sympy as sp, time, pickle, sys
from flint import arb
import ibp_sym as S
t0 = time.time()
SOURCES = []
class Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp()
        assert ex.is_Integer, e
        n = int(ex); bb = f"({self._print(b)})"
        if n == 0: return "arb(1)"
        prod = "*".join([bb] * abs(n))
        return f"({prod})" if n > 0 else f"(1/({prod}))"
def make_func(exprs, args, name='f'):
    """exprs: dict name->sympy expr. returns python function(*args)->dict of arb"""
    keys = list(exprs)
    reps, red = sp.cse([exprs[k_] for k_ in keys], optimizations='basic')
    pr = Pr()
    lines = [f"def {name}({', '.join(map(str, args))}):"]
    for s, e in reps: lines.append(f"    {s} = {pr.doprint(e)}")
    lines.append("    return {" + ", ".join(f"{k_!r}: {pr.doprint(e)}" for k_, e in zip(keys, red)) + "}")
    src = "\n".join(lines); SOURCES.append(src)
    ns = {"arb": arb}; exec(src, ns)
    return ns[name], len(reps)
# ---------------------------------------------------------------- linear part:  Lam(alpha)/..., d/dalpha
psi = sp.symbols('ps0:5'); hh = sp.symbols('hh0:5'); al_ = sp.Symbol('alpha')
def d_alpha(e):
    r = 0
    for i in range(5):
        if e.has(S.a[i]) or True:
            pass
    # chain rule through a_i (explicit), alb, Jb
    da_ = [psi[i] + 2 * al_ * hh[i] for i in range(5)]
    r = sum(sp.diff(e, S.a[i]) * da_[i] for i in range(5))
    a0p, a1p = da_[0], da_[1]
    r += sp.diff(e, S.alb) * (-S.k * a0p) + sp.diff(e, S.Jb) * (S.Jb**3 * S.alb * S.k * a0p - S.Jb**3 * S.a[1] * a1p)
    return sp.expand(r)
LamA = d_alpha(S.Lam)
LIN_ARGS = [*S.Y[:5], *S.a[:5], *psi, *hh, al_, S.k, S.m, S.p, S.q, S.alb, S.Jb]
LinF, nl = make_func({'LamA': LamA, 'Lam': S.Lam}, LIN_ARGS, 'linf')
print(f"[{time.time()-t0:.0f}s] lin codegen, cse {nl}", flush=True)
# ---------------------------------------------------------------- quadratic part families
w0, w1, w2 = S.w[0], S.w[1], S.w[2]
fam = {}
for r in range(3):
    Qr = sp.expand(sp.diff(S.QR, S.Y[r]))
    P_ = sp.Poly(Qr, w2); assert P_.degree() <= 2
    Qj = [sp.expand(P_.coeff_monomial(w2**j)) for j in range(3)]
    def D(e, n0, n1):
        for _ in range(n0): e = S.dw(e, 0)
        for _ in range(n1): e = S.dw(e, 1)
        return e
    z = {S.alf: S.alb, S.Jf: S.Jb, w0: 0, w1: 0}
    at0 = lambda e: sp.expand(e.subs(z, simultaneous=True))
    # quadratic coefficients at w = 0
    fam[f'c00_{r}'] = at0(D(Qj[0], 2, 0)) / 2; fam[f'c01_{r}'] = at0(D(Qj[0], 1, 1)); fam[f'c11_{r}'] = at0(D(Qj[0], 0, 2)) / 2
    fam[f'c02_{r}'] = at0(D(Qj[1], 1, 0)); fam[f'c12_{r}'] = at0(D(Qj[1], 0, 1)); fam[f'c22_{r}'] = at0(Qj[2])
    # vanishing checks (symbolic): Q0(0), grad Q0(0), Q1(0)
    assert at0(Qj[0]) == 0 and at0(D(Qj[0], 1, 0)) == 0 and at0(D(Qj[0], 0, 1)) == 0 and at0(Qj[1]) == 0
    # remainder families (evaluated over the jet box)
    for s in range(4): fam[f't3_{s}_{r}'] = D(Qj[0], 3 - s, s)
    for s in range(3): fam[f't2_{s}_{r}'] = D(Qj[1], 2 - s, s)
    for s in range(2): fam[f't1_{s}_{r}'] = D(Qj[2], 1 - s, s)
    print(f"[{time.time()-t0:.0f}s] family r={r}", flush=True)
QA = [S.k, S.m, S.p, S.q, *S.a[:5], w0, w1, S.alf, S.Jf, S.alb, S.Jb]
famc = {k_: v for k_, v in fam.items() if k_[0] == 'c'}; famt = {k_: v for k_, v in fam.items() if k_[0] == 't'}
for v in famc.values(): assert not (v.has(w0) or v.has(w1) or v.has(S.alf) or v.has(S.Jf))
CoefF, nq = make_func(famc, [S.k, S.m, S.p, S.q, *S.a[:5], S.alb, S.Jb], 'coeff')
RemF, nq2 = make_func(famt, QA, 'remf')
print(f"[{time.time()-t0:.0f}s] quad codegen, cse {nq},{nq2}, {len(fam)} expressions", flush=True)
# Full QR (to test)
QRF, _ = make_func({'QR': S.QR}, [*S.Y[:3], S.k, S.m, S.p, S.q, *S.a[:5], w0, w1, S.w[2], S.alf, S.Jf, S.alb, S.Jb], 'qrf')
if __name__ == "__main__":
    open('ibp_gen.py', 'w').write("from flint import arb\n" + "\n\n".join(SOURCES) + "\n")
    print("wrote ibp_gen.py (functions linf, coeff, remf, qrf)")
