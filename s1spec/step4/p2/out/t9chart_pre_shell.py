"""Review B5 (T9): exact variation of the graph slope over the S1 ray, in arb (python-flint), no sympy.

Setting (as p2/t9.py).  Base Gamma~(x), x = arclength from the tip, theta~ = tangent angle, k = theta~' (k0..k3 = k, k', k'', k'''),
p = <Gamma~, N>, q = <Gamma~, T> (N = iT).  Ray R = e^{i Theta}, phi = theta~ - Theta, so R = cos(phi) T - sin(phi) N and
  Y_b = <Gamma~, R> = q cos(phi) - p sin(phi),  Y_b' = cos(phi),  base slope over the ray F_b = tan(phi).
True curve gamma = Gamma~ + e N (e = T8's error, x the base arclength): gamma' = a T + b N with a = 1 - k e, b = e', hence at the same x
  Y_t = Y_b - e sin(phi),  M_t := Y_t' = a cos(phi) - b sin(phi),  F_t = tan(phi + beta) = (a sin(phi) + b cos(phi)) / M_t.
D = y d/dy is D_t = (Y_t/M_t) d/dx along gamma and D_b = (Y_b/cos phi) d/dx along Gamma~.  With Q^t_j = D_t^j F_t, Q^b_j = D_b^j F_b and
X_t, X_b the inverses of Y_t, Y_b, exactly:
  D^j f_t(y) - D^j f_b(y) = [Q^t_j - Q^b_j](X_t(y)) + [Q^b_j(X_t(y)) - Q^b_j(X_b(y))].                                     (*)
Second bracket: Y_b(X_t) - Y_b(X_b) = e(X_t) sin(phi(X_t)), so |X_t - X_b| <= |e| sup|sin phi| / inf cos phi.
First bracket: Taylor mode in u = (x - x0)/y0, y0 = x0 - X0 (so every series coefficient is O(1)), in the dimensionless jets
  eta_i = y0^(i-1) e^(i)(x0),  i = 0..4   (e/y0, e', y0 e'', y0^2 e''', y0^3 e'''').
Each Taylor coefficient is an 'E':  c + sum_i l_i eta_i + NL  with |NL| <= sum_i r_i |eta_i|;  c and l are exact (balls), r >= 0.
Products and the reciprocal charge a pair eta_i eta_k (k <= i) to eta_i with the pointwise box |eta_k| <= B_k (k <= 3); truncation
at order 3 never pairs eta_4 (or eta_3) with itself, so eta_4 needs no box (asserted: B_4 = None raises).
At first order (*) is D^j(sec^2 phi e' + sec^2 phi tan phi k e), the exact first variation (review B5); this module does not use that
form, but p2/t9chart_check.py compares the two (and tests the remainder bounds)."""
from flint import arb

NV = 5
BOX = [None] * NV                 # |eta_i| <= BOX[i] at the current point (set by the caller); BOX[4] stays None
ONE = arb(1)
FACT = [1, 1, 2, 6, 24, 120]


def absu(v):
    """an arb >= |v|"""
    return abs(v).upper()


class E:
    __slots__ = ('c', 'l', 'r')

    def __init__(s, c, l=None, r=None):
        s.c = c if isinstance(c, arb) else arb(c)
        s.l = l if l is not None else {}
        s.r = r if r is not None else {}


def _acc(d, i, v):
    d[i] = d[i] + v if i in d else v


def eadd(X, Y, sgn=1):
    l = dict(X.l); r = dict(X.r)
    for i, v in Y.l.items(): _acc(l, i, v if sgn > 0 else -v)
    for i, v in Y.r.items(): _acc(r, i, v)
    return E(X.c + Y.c if sgn > 0 else X.c - Y.c, l, r)


def escal(X, k):
    ka = absu(arb(k))
    return E(X.c * k, {i: v * k for i, v in X.l.items()}, {i: v * ka for i, v in X.r.items()})


def epart(X):
    """a_i with |X - X.c| <= sum_i a_i |eta_i|"""
    a = {i: absu(v) for i, v in X.l.items()}
    for i, v in X.r.items(): _acc(a, i, v)
    return a


def _cross(ax, ay, r, f=None):
    """r += (sum_i ax_i |eta_i|)(sum_k ay_k |eta_k|), the pair (i, k) charged to max(i, k) with the box of min(i, k)"""
    for i, u in ax.items():
        for k, v in ay.items():
            hi, lo = (i, k) if i >= k else (k, i)
            if BOX[lo] is None:
                raise ValueError(f'pair (eta_{i}, eta_{k}) needs the box of eta_{lo}')
            t = u * v * BOX[lo]
            _acc(r, hi, t * f if f is not None else t)


def emul(X, Y):
    l, r = {}, {}
    for i, v in Y.l.items(): _acc(l, i, X.c * v)
    for i, v in X.l.items(): _acc(l, i, Y.c * v)
    xa, ya = absu(X.c), absu(Y.c)
    for i, v in Y.r.items(): _acc(r, i, xa * v)
    for i, v in X.r.items(): _acc(r, i, ya * v)
    ax, ay = epart(X), epart(Y)
    if ax and ay: _cross(ax, ay, r)
    return E(X.c * Y.c, l, r)


def enclose(X):
    """ball enclosing X for every eta in the boxes"""
    a = epart(X)
    if not a: return X.c
    D = arb(0)
    for i, v in a.items():
        if BOX[i] is None: raise ValueError(f'enclose needs the box of eta_{i}')
        D += v * BOX[i]
    return X.c + D * arb(0, 1)


def erecip(X):
    """1/X = 1/c - d/c^2 + d^2/(c^2 X), d = X - c"""
    c = 1 / X.c; c2 = c * c; c2a = absu(c2)
    a = epart(X)
    l = {i: -v * c2 for i, v in X.l.items()}
    r = {i: v * c2a for i, v in X.r.items()}
    if a:
        D = arb(0)
        for i, v in a.items():
            if BOX[i] is None: raise ValueError(f'reciprocal needs the box of eta_{i}')
            D += v * BOX[i]
        xmin = abs(X.c).lower() - D.upper()
        if not xmin > 0: raise ValueError('reciprocal: X may vanish on the boxes')
        _cross(a, a, r, c2a / xmin)
    return E(c, l, r)


# ---- truncated Taylor series (lists of E, index = order)
def s_mul(X, Y, n):
    out = []
    for k in range(n + 1):
        acc = emul(X[0], Y[k])
        for i in range(1, k + 1): acc = eadd(acc, emul(X[i], Y[k - i]))
        out.append(acc)
    return out


def s_add(X, Y, n, sgn=1): return [eadd(X[k], Y[k], sgn) for k in range(n + 1)]


def s_der(X): return [escal(X[k + 1], k + 1) for k in range(len(X) - 1)]


def s_recip(X, n):
    R = [erecip(X[0])]
    for k in range(1, n + 1):
        acc = emul(X[1], R[k - 1])
        for i in range(2, k + 1): acc = eadd(acc, emul(X[i], R[k - i]))
        R.append(escal(emul(R[0], acc), -1))
    return R


# ---- plain arb series (base only)
def a_mul(X, Y, n): return [sum((X[i] * Y[k - i] for i in range(k + 1)), arb(0)) for k in range(n + 1)]


def a_recip(X, n):
    R = [1 / X[0]]
    for k in range(1, n + 1): R.append(-R[0] * sum((X[i] * R[k - i] for i in range(1, k + 1)), arb(0)))
    return R


def a_der(X): return [X[k + 1] * (k + 1) for k in range(len(X) - 1)]


def a_sincos(ph, n):
    S, C = [ph[0].sin()], [ph[0].cos()]
    for k in range(1, n + 1):
        S.append(sum((i * ph[i] * C[k - i] for i in range(1, k + 1)), arb(0)) / k)
        C.append(-sum((i * ph[i] * S[k - i] for i in range(1, k + 1)), arb(0)) / k)
    return S, C


def base_series(y0, phi0, k, ratio, n=4):
    """phi(u), kappa-hat(u) = d phi/du = y0 k(x0 + y0 u), sin, cos, Yhat_b(u) = Y_b/y0 (dYhat_b/du = cos phi), orders n, n-1, n, n, n"""
    ph = [phi0] + [y0 ** (i + 1) * k[i] / FACT[i + 1] for i in range(n)]
    kh = [y0 ** (i + 1) * k[i] / FACT[i] for i in range(n)]
    S, C = a_sincos(ph, n)
    Yb = [ratio] + [C[i] / (i + 1) for i in range(n)]
    return ph, kh, S, C, Yb


def base_Q(S, C, Yb, n=4):
    """Q^b_j (j = 0..3) as u-series (orders n, n-1, n-2, n-3) from F_b = S/C and D_b = (Yhat_b / C) d/du"""
    iC = a_recip(C, n)
    Q = [a_mul(S, iC, n)]
    Rb = a_mul(Yb, iC, n - 1)
    for j in range(3): Q.append(a_mul(Rb, a_der(Q[j]), n - 1 - j))
    return Q


def true_Q(kh, S, C, Yb, eta=None):
    """Q^t_j (j = 0..3) at u = 0 as E's, with M_t(0) and Yhat_t(0).  eta: None -> symbolic jets (E variables);
    a list of 5 balls -> plain evaluation at those jets (for tests)."""
    if eta is None:
        eh = [E(0, {i: ONE / FACT[i]}) for i in range(4)]              # e/y0 = sum eta_i u^i / i!
        ed = [E(0, {i + 1: ONE / FACT[i]}) for i in range(4)]          # d(e/y0)/du = e' = sum eta_{i+1} u^i / i!
    else:
        eh = [E(eta[i] / FACT[i]) for i in range(4)]
        ed = [E(eta[i + 1] / FACT[i]) for i in range(4)]
    Se, Ce, Ke, Ye = ([E(v) for v in X] for X in (S, C, kh, Yb))
    a = s_add([E(1), E(0), E(0), E(0)], s_mul(Ke, eh, 3), 3, -1)       # 1 - k e
    Mt = s_add(s_mul(a, Ce, 3), s_mul(ed, Se, 3), 3, -1)                # a cos - b sin
    Nt = s_add(s_mul(a, Se, 3), s_mul(ed, Ce, 3), 3)                     # a sin + b cos
    iM = s_recip(Mt, 3)
    Q = [s_mul(Nt, iM, 3)]                                               # F_t
    Yt = s_add(Ye, s_mul(eh, Se, 2), 2, -1)                              # Yhat_t = Yhat_b - (e/y0) sin
    R = s_mul(Yt, iM, 2)
    for j in range(3): Q.append(s_mul(R, s_der(Q[j]), 2 - j))
    return [Qj[0] for Qj in Q], Mt[0], Yt[0]
