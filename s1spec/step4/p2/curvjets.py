"""Review B5: the true curvature jets of curve = base + e N (normal graph over the base, base arclength s), exactly, for T9.
  X' = (1 - k e) T + e' N,  X'' = (-k' e - 2 k e') T + (k (1 - k e) + e'') N,  J = |X'| = sqrt((1 - k e)^2 + e'^2),
  kappa = X' x X'' / J^3 = [(1 - k e)(k (1 - k e) + e'') + e'(k' e + 2 k e')] / J^3,   K_j = (J^-1 d/ds)^j kappa  (true arclength).
K_j(k_0..k_{j+1}, e_0..e_{j+2}) by sympy, evaluated on arb boxes (no typed allowances).  kmajor(...) returns sup|K_j|, j = 0..3."""
import math
import sympy as sp
from flint import arb
from sympy.printing.pycode import PythonCodePrinter

class _P(PythonCodePrinter):
    def _print_Pow(self, expr, rational=False):
        b, e = expr.base, expr.exp
        if e.is_Rational and e.q in (1, 2):
            n = abs(int(e.p)); unit = '(' + self._print(b) + ')' if e.q == 1 else 'sqrt(' + self._print(b) + ')'
            prod = '*'.join([unit] * n) if n > 0 else '1'
            return '(' + prod + ')' if e > 0 else '(1/(' + prod + '))'
        return super()._print_Pow(expr)

_s = sp.symbols('s'); _k = sp.Function('k')(_s); _e = sp.Function('e')(_s)
_J = sp.sqrt((1 - _k * _e) ** 2 + sp.diff(_e, _s) ** 2)
_kap = ((1 - _k * _e) * (_k * (1 - _k * _e) + sp.diff(_e, _s, 2)) + sp.diff(_e, _s) * (sp.diff(_k, _s) * _e + 2 * _k * sp.diff(_e, _s))) / _J ** 3
KS = sp.symbols('k0:5'); ES = sp.symbols('e0:6')
def _sub(ex):
    for j in range(5, 0, -1): ex = ex.subs(sp.Derivative(_e, (_s, j)), ES[j])
    for j in range(4, 0, -1): ex = ex.subs(sp.Derivative(_k, (_s, j)), KS[j])
    return ex.subs(_e, ES[0]).subs(_k, KS[0])
_Ks = [_kap]
for _ in range(3): _Ks.append(sp.diff(_Ks[-1], _s) / _J)
_mods = [{'sqrt': lambda a: a.sqrt()}, 'math']
FK = [sp.lambdify(list(KS) + list(ES), _sub(ex), modules=_mods, printer=_P) for ex in _Ks]
FA = sp.lambdify(list(KS) + list(ES), sp.diff(_sub(_Ks[3]), ES[5]), modules=_mods, printer=_P)   # coefficient of e^(5) in K_3
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)

def kmajor(kb, eb, n=4):
    """kb: arb hulls of base k, k', k'', k''' (k'''' optional); eb: bounds E_0..E_5 >= 0 of |e^(j)| (E_5 only for the e5-part, see A3)
    returns sup|K_j| for j < n (K_3 without its e^(5) term) and sup|dK_3/de5|"""
    k = list(kb) + [arb(0)] * (5 - len(kb)); e = [arb(0, float(v)) for v in eb] + [arb(0)] * (6 - len(eb))
    out = [U(FK[j](*k, *e)) for j in range(min(n, 3))]
    if n > 3:
        e5 = list(e); e5[5] = arb(0)
        out.append(U(FK[3](*k, *e5))); out.append(U(FA(*k, *e)))
    return out
