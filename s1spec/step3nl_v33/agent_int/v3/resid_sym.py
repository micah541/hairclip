"""Symbolic residuals and their total x-derivatives (profile ODE p' = -kq, q' = 1+kp, k' = m, m' = p/4; jets Y_i = psi3~^(i), H_i = h2~^(i)):
  R_psi = (L - lamt) psi3~ = -Y4 - k^2 Y2 - (5km + q/4) Y1 + (1/4 - LAM - 3m^2 - kp) Y0                      (LAM = lamt)
  R_h2  = (L - 2 lamt) h2~ + Q2(psi3~) = same operator on H with LAM = 2 lamt, + Q2(Y; p, k, m)              (Q2 from q2_derivs.pkl)
  DR^j = j-th total derivative (exact, sympy).  Lambdified to arb ball arithmetic (explicit products, no ball ** n)."""
import pickle, sympy as sp
from flint import arb
NY = 9
Y = sp.symbols('Y0:%d' % NY); H = sp.symbols('H0:%d' % NY); p, q, k, m, LAM = sp.symbols('p q k m LAM')
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4}
def D(e):
    out = sum(sp.diff(e, v) * r for v, r in rules.items())
    for i in range(NY - 1):
        out += sp.diff(e, Y[i]) * Y[i + 1] + sp.diff(e, H[i]) * H[i + 1]
    return sp.expand(out)
Q = pickle.load(open('/home/claude/s1spec/step3nl/q2_derivs.pkl', 'rb'))
_u = sp.symbols('u0:7'); _pp, _qq, _kk, _mm = sp.symbols('pp qq kk mm')
Q2s = Q['Q2'].subs({s: Y[i] for i, s in enumerate(_u)}, simultaneous=True)
Q2s = Q2s.subs({sp.Symbol('pp'): p, sp.Symbol('qq'): q, sp.Symbol('kk'): k, sp.Symbol('mm'): m}, simultaneous=True)
def op(U): return -U[4] - k**2 * U[2] - (5 * k * m + q / 4) * U[1] + (sp.Rational(1, 4) - LAM - 3 * m**2 - k * p) * U[0]
Rpsi = sp.expand(op(Y)); Rh = sp.expand(op(H) + Q2s)
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
ARGS = list(Y[:8]) + list(H[:8]) + [p, q, k, m, LAM]
def _mk(e):
    assert not (e.free_symbols - set(ARGS)), e.free_symbols - set(ARGS)
    return eval("lambda " + ",".join(map(str, ARGS)) + ": " + _Pr().doprint(e), {"arb": arb})
def derivs(R, n=3):
    out = [R]
    for _ in range(n): out.append(D(out[-1]))
    return out
EXPR = {'psi': derivs(Rpsi), 'h2': derivs(Rh)}
FN = {w: [_mk(e) for e in EXPR[w]] for w in EXPR}
def call(which, order, Yj, Hj, pp, qq, kk, mm, lam):
    Yj = list(Yj) + [arb(0)] * (8 - len(Yj)); Hj = list(Hj) + [arb(0)] * (8 - len(Hj))
    return FN[which][order](*Yj[:8], *Hj[:8], pp, qq, kk, mm, lam)
if __name__ == "__main__":
    print({w: [sp.count_ops(e) for e in EXPR[w]] for w in EXPR})
