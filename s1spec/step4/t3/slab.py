"""slab.py -- T3 prototype: validated slabs for HLS curve-diffusion flow in the sine class.

u(s) = sum_n u_n sin(n a), theta = -pi/2 + a + u;  clock tau' (s), lam_n = n^4 - n^2:
    u_n' = -lam_n u_n + N_n(u),   N = -d(v^2 + v^3/3) + (1+v) dinv(w^2 - <w^2>),  v = u_a, w = u_aa.
Method (T3_DESIGN.md): Chebyshev-in-time base (float collocation + arb polish), error e = P(s) rho + z(s),
P = float Chebyshev propagator with arb defect, z by Duhamel self-consistent bounds (componentwise for
n <= Nc, tail T n^-q beyond), Lohner QR block update.  All rigorous quantities in arb (PREC bits).
"""
import numpy as np, time, math
from numpy.fft import rfft, irfft
from scipy.linalg import lu_factor, lu_solve, qr as sqr
from scipy.sparse.linalg import gmres, LinearOperator
from flint import arb, arb_poly, arb_mat, ctx

PREC = 128
ctx.prec = PREC
A0 = arb(0)
QT = 6                    # tail decay exponent q
DEBUG = False
MU = 0.03                 # relative margin of the a priori bounds
ETACAP = 1e-2             # a priori cap (relative to w) on the intra-slab absorbed radii


# ============================================================ float model (sine-class Galerkin, M modes)
class FM:
    def __init__(self, M):
        self.M = M
        Ng = 4 * M + 8
        self.Ng = Ng + (Ng % 2)
        self.k = np.arange(self.Ng // 2 + 1, dtype=float)
        self.n = np.arange(1, M + 1, dtype=float)
        self.lam = self.n ** 4 - self.n ** 2
        kk = self.k.copy(); kk[0] = 1.0
        self.inv_ik = 1.0 / (1j * kk); self.inv_ik[0] = 0.0

    def spec(self, U):
        X = np.zeros(U.shape[:-1] + (self.Ng // 2 + 1,), complex)
        X[..., 1:self.M + 1] = -0.5j * self.Ng * U
        return X

    def sin_of(self, f):
        return -2.0 * rfft(f)[..., 1:self.M + 1].imag / self.Ng

    def fields(self, U):
        X = self.spec(U)
        v = irfft(1j * self.k * X, self.Ng); w = irfft(-(self.k ** 2) * X, self.Ng)
        Xw = rfft(w * w); Xw[..., 0] = 0
        W = irfft(Xw * self.inv_ik, self.Ng)
        return v, w, W

    def N(self, U):
        v, w, W = self.fields(U)
        dA = irfft(1j * self.k * rfft(v * v + v ** 3 / 3), self.Ng)
        return self.sin_of(-dA + (1 + v) * W)

    def Japply(self, flds, E):
        """flds=(v,w,W) with shape (nodes,Ng); E shape (nodes,cols,M) -> (nodes,cols,M)"""
        v, w, W = (f[:, None, :] for f in flds)
        X = self.spec(E)
        eps = irfft(1j * self.k * X, self.Ng); eta = irfft(-(self.k ** 2) * X, self.Ng)
        dA = irfft(1j * self.k * rfft((2 * v + v * v) * eps), self.Ng)
        Xy = rfft(2 * w * eta); Xy[..., 0] = 0
        Y = irfft(Xy * self.inv_ik, self.Ng)
        return self.sin_of(-dA + eps * W + (1 + v) * Y)

    def thmax(self, U):
        X = self.spec(U)
        return np.abs(1 + irfft(1j * self.k * X, self.Ng)).max(axis=-1)


def cheb(K):
    """Chebyshev-Lobatto nodes x_j = cos(pi j/K) (x_0 = 1) and differentiation matrix (Trefethen)."""
    x = np.cos(np.pi * np.arange(K + 1) / K)
    c = np.hstack([2.0, np.ones(K - 1), 2.0]) * (-1.0) ** np.arange(K + 1)
    X = np.tile(x, (K + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(K + 1))
    return x, D - np.diag(D.sum(axis=1))


def vinv_f(K):
    j = np.arange(K + 1); cb = np.ones(K + 1); cb[0] = cb[-1] = 2
    return (2.0 / K) * np.cos(np.pi * np.outer(j, j) / K) / np.outer(cb, cb)


def teval_f(K, x):
    return np.cos(np.outer(np.arccos(np.clip(x, -1, 1)), np.arange(K + 1)))


class Colloc:
    """Float Lobatto collocation on [0,h] for L d = (2/h) D d + Lam d - J(s_i) d (rows i<K), d_K (row K: s=0).
    GMRES, preconditioner: dense block for modes <= Mf, (2/h D + lam_n) per mode above."""

    def __init__(self, fm, h, K, Unodes, Mf, Ms=None):
        """initial condition imposed for modes <= Ms; stiffer modes: ODE also at s=0 (slaved start)"""
        Ms = fm.M if Ms is None else Ms
        self.fm, self.h, self.K, self.Mf, self.Ms = fm, h, K, Mf, Ms
        self.x, self.D = cheb(K)
        self.flds = fm.fields(Unodes)
        M = fm.M; D2 = (2.0 / h) * self.D
        E = np.broadcast_to(np.eye(M)[:Mf], (K + 1, Mf, M))
        JE = fm.Japply(self.flds, E)
        Jll = JE[:, :, :Mf].transpose(0, 2, 1)
        nn = K + 1
        A = np.kron(D2, np.eye(Mf))
        for i in range(K):
            A[i * Mf:(i + 1) * Mf, i * Mf:(i + 1) * Mf] += np.diag(fm.lam[:Mf]) - Jll[i]
        A[K * Mf:, :] += 0
        for n in range(Mf):
            r = K * Mf + n
            if n < Ms: A[r, :] = 0; A[r, r] = 1
            else: A[r, r] += fm.lam[n] - 0; A[r, K * Mf:(K + 1) * Mf] -= Jll[K][n]
        self.lu = lu_factor(A)
        Ah = np.broadcast_to(D2, (M - Mf, nn, nn)).copy()
        Ah[:, np.arange(nn), np.arange(nn)] += fm.lam[Mf:, None]
        ic = np.arange(Mf, M) < Ms
        Ah[ic, K, :] = 0; Ah[ic, K, K] = 1
        self.Ahinv = np.linalg.inv(Ah)

    def apply(self, d):
        out = (2.0 / self.h) * (self.D @ d) + self.fm.lam * d - self.fm.Japply(self.flds, d[:, None, :])[:, 0, :]
        out[self.K, :self.Ms] = d[self.K, :self.Ms]
        return out

    def prec(self, r):
        K, Mf = self.K, self.Mf
        out = np.empty_like(r)
        out[:, :Mf] = lu_solve(self.lu, r[:, :Mf].reshape(-1)).reshape(K + 1, Mf)
        out[:, Mf:] = np.einsum('nij,jn->in', self.Ahinv, r[:, Mf:])
        return out

    def apply_blk(self, d):
        """d: (K+1, ncol, M)"""
        out = (2.0 / self.h) * np.einsum('ij,jcn->icn', self.D, d) + self.fm.lam * d - self.fm.Japply(self.flds, d)
        out[self.K, :, :self.Ms] = d[self.K, :, :self.Ms]
        return out

    def prec_blk(self, r):
        K, Mf = self.K, self.Mf; nc = r.shape[1]
        out = np.empty_like(r)
        lo = r[:, :, :Mf].transpose(0, 2, 1).reshape((K + 1) * Mf, nc)
        out[:, :, :Mf] = lu_solve(self.lu, lo).reshape(K + 1, Mf, nc).transpose(0, 2, 1)
        out[:, :, Mf:] = np.einsum('nij,jcn->icn', self.Ahinv, r[:, :, Mf:])
        return out

    def solve_blk(self, rhs, tol=1e-12, maxit=40):
        """preconditioned Richardson for many columns at once (float floor ~1e-10 relative in the true residual)"""
        X = np.zeros_like(rhs); sc = np.abs(rhs).max(); hist = []
        for it in range(maxit):
            R = rhs - self.apply_blk(X); dX = self.prec_blk(R); X = X + dX
            hist.append(np.abs(dX).max() / max(np.abs(X).max(), 1e-300))
            if hist[-1] < tol: break
        return X, hist

    def solve(self, rhs, tol=1e-12):
        sh = rhs.shape; N = rhs.size
        Lop = LinearOperator((N, N), matvec=lambda v: self.apply(v.reshape(sh)).ravel())
        Pop = LinearOperator((N, N), matvec=lambda v: self.prec(v.reshape(sh)).ravel())
        sc = np.abs(rhs).max()
        if sc == 0: return np.zeros(sh), 0
        x, info = gmres(Lop, rhs.ravel() / sc, M=Pop, rtol=tol, atol=0, restart=40, maxiter=2)
        return x.reshape(sh) * sc, info


def base_float(fm, u0, h, K, Uguess, Mf, Ms=None, maxit=12, verbose=False):
    Ms = fm.M if Ms is None else Ms
    x, D = cheb(K)
    U = Uguess.copy(); U[K] = u0
    for it in range(maxit):
        F = (2.0 / h) * (D @ U) + fm.lam * U - fm.N(U); F[K, :Ms] = U[K, :Ms] - u0[:Ms]
        col = Colloc(fm, h, K, U, Mf, Ms)
        dU, info = col.solve(-F)
        U = U + dU
        if verbose: print(f'   float Newton {it}: |F|={np.abs(F[:K]).max():.2e} |dU|={np.abs(dU).max():.2e} gmres={info}')
        if np.abs(dU).max() < 1e-14 * np.abs(U).max(): break
    return U



def solve_base(fm, x0arb, xf, h, K, Mf, Ms, Uguess=None, npol=3):
    """float Lobatto-collocation Newton + arb polish (residual at the collocation nodes in arb).
    Returns nodal arb values Ua, Chebyshev coefficient matrix C (arb), its rows, polish history, float nodal U."""
    M = fm.M; harb = arb(h)
    lam_arb = [arb(n) ** 4 - arb(n) ** 2 for n in range(1, M + 1)]
    if Uguess is None: Uguess = np.tile(xf, (K + 1, 1))
    U = base_float(fm, xf, h, K, Uguess, Mf, Ms, verbose=False)
    Varb = cheb_arb(K); Darb = teval_arb(K, K, True) * Varb
    Ua = [arb_vec(U[i]) for i in range(K + 1)]
    for n in range(Ms): Ua[K][n] = x0arb[n]
    col = Colloc(fm, h, K, U, Mf, Ms)
    pol = []
    for it in range(npol):
        DU = rows_of(Darb * mat_from_rows(Ua))
        F = []
        for i in range(K + 1):
            Ni = N_arb(Ua[i])
            F.append([DU[i][n] * 2 / harb + Ua[i][n] * lam_arb[n] - Ni[n] for n in range(M)])
        for n in range(Ms): F[K][n] = Ua[K][n] - x0arb[n]
        Ff = np.array([[float(x.mid()) for x in r] for r in F])
        pol.append(np.abs(Ff[:K]).max())
        dU, info = col.solve(-Ff)
        Ua = [[Ua[i][n] + float(dU[i, n]) for n in range(M)] for i in range(K + 1)]
    C = Varb * mat_from_rows(Ua)
    return Ua, C, rows_of(C), pol, U


def prelim_slab(x0arb, fm, h, K=16, Mf=60):
    """non-validated arb base solve; returns exact (128-bit midpoint) end values u(h) -- a consistent datum"""
    Ms = int(np.sum(fm.lam * h <= 8.0 * K * K))
    xf = np.array([float(x.mid()) for x in x0arb])
    Ua, C, Crows, pol, U = solve_base(fm, x0arb, xf, h, K, Mf, Ms)
    u_h = [sum((Crows[k][n] for k in range(K + 1)), A0) for n in range(fm.M)]
    return [arb(x.mid()) for x in u_h], pol


# ============================================================ arb Chebyshev grids
_cache = {}


def arb_cos_pi(p, q):
    return (arb(p) / q).cos_pi()


def cheb_arb(K):
    """arb matrices: Vinv (values at Lobatto nodes -> Chebyshev coefficients)."""
    key = ('vinv', K)
    if key not in _cache:
        rows = []
        for k in range(K + 1):
            ck = 2 if k in (0, K) else 1
            for j in range(K + 1):
                cj = 2 if j in (0, K) else 1
                rows.append(arb_cos_pi(j * k, K) * 2 / (K * ck * cj))
        _cache[key] = arb_mat(K + 1, K + 1, rows)
    return _cache[key]


def teval_arb(K, K2, deriv=False):
    """(K2+1) x (K+1): T_k(x_i) (or T_k'(x_i)) at Lobatto nodes x_i = cos(pi i/K2)."""
    key = ('te', K, K2, deriv)
    if key not in _cache:
        ent = []
        for i in range(K2 + 1):
            for k in range(K + 1):
                if not deriv:
                    ent.append(arb_cos_pi(i * k, K2))
                else:
                    if i == 0: ent.append(arb(k * k))
                    elif i == K2: ent.append(arb((-1) ** (k + 1) * k * k))
                    else: ent.append(k * (arb(i * k) / K2).sin_pi() / (arb(i) / K2).sin_pi())
        _cache[key] = arb_mat(K2 + 1, K + 1, ent)
    return _cache[key]


def sup_sub(CM, Kd, a, b):
    """rigorous sup over s/h in [a,b] of each column of a Chebyshev series (coefficient matrix CM, (Kd+1) x ncol,
    variable x = 2s/h-1): re-expand on [a,b] (exact: same degree) and sum |coefficients|."""
    key = ('sub', Kd, a, b)
    if key not in _cache:
        xa, xb = 2 * arb(a) - 1, 2 * arb(b) - 1
        ent = []
        for i in range(Kd + 1):
            xi = (xa + xb) / 2 + (xb - xa) / 2 * arb_cos_pi(i, Kd)
            T0, T1 = arb(1), xi
            row = [T0, T1]
            for k in range(2, Kd + 1): row.append(2 * xi * row[-1] - row[-2])
            ent.extend(row[:Kd + 1])
        _cache[key] = cheb_arb(Kd) * arb_mat(Kd + 1, Kd + 1, ent)
    loc = rows_of(_cache[key] * CM)
    return [_ub(sum((abs(loc[k][n]) for k in range(Kd + 1)), A0)) for n in range(len(loc[0]))]


def cheb_two(a):
    return [x / 2 for x in reversed(a[1:])] + [a[0]] + [x / 2 for x in a[1:]]


def cheb_mul(a, b):
    """product of Chebyshev series (coefficient lists, arb): T_i T_j = (T_{i+j} + T_{|i-j|})/2, exact"""
    C = conv(cheb_two(a), cheb_two(b)); o = len(a) + len(b) - 2
    return [C[o]] + [2 * C[o + k] for k in range(1, o + 1)]


def cheb_add(a, b):
    n = max(len(a), len(b)); a = a + [A0] * (n - len(a)); b = b + [A0] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def cheb_intdef(a):
    """int_{-1}^{1} sum a_k T_k dx"""
    return sum((a[k] * 2 / (1 - k * k) for k in range(0, len(a), 2)), A0)


def cheb_antider(a):
    """coefficients of int_{-1}^{x} sum a_k T_k"""
    out = [A0] * (len(a) + 1)
    out[1] = out[1] + a[0]
    if len(a) > 1: out[2] = out[2] + a[1] / 4
    for k in range(2, len(a)):
        out[k + 1] = out[k + 1] + a[k] / (2 * (k + 1)); out[k - 1] = out[k - 1] - a[k] / (2 * (k - 1))
    out[0] = out[0] - sum((out[k] * (-1) ** k for k in range(len(out))), A0)
    return out


EB = arb(0, 1)          # the unit ball [-1, 1]


def kmax_enclosure(xarb, e1, L, Ng=2048, chunk=256):
    """rigorous enclosure of k_max = (2pi/L) max_alpha |theta_alpha| for theta = -pi/2 + alpha + sum x_n sin(n alpha) + e,
    |e|_1 <= e1.  Grid alpha_j = (j+1/2) pi/Ng on [0,pi] (theta_alpha is even and 2pi-periodic); on each cell
    |f| <= |f_j| + |f'_j| d + |f''_j| d^2/2 + (sum n^4|x_n|) d^3/6, d = pi/(2Ng).  cos(n alpha_j) = cos(pi k/(2Ng)),
    k = n(2j+1) mod 4Ng, from an exact table (no recurrence: interval recurrences inflate radii)."""
    M = len(xarb); P4 = 4 * Ng
    tc = [(arb(k) / (2 * Ng)).cos_pi() for k in range(P4)]; ts = [(arb(k) / (2 * Ng)).sin_pi() for k in range(P4)]
    co = arb_mat(M, 2, [v for n in range(1, M + 1) for v in (xarb[n - 1] * n, -xarb[n - 1] * (n * n * n))])
    si = arb_mat(M, 1, [-xarb[n - 1] * (n * n) for n in range(1, M + 1)])
    S4 = sum((abs(xarb[n - 1]) * arb(n) ** 4 for n in range(1, M + 1)), A0)
    dlt = arb.pi() / (2 * Ng)
    up = A0; lo = A0; jmax = 0; fmaxf = -1.0; cells = []
    for j0 in range(0, Ng, chunk):
        js = range(j0, min(j0 + chunk, Ng))
        Cm = arb_mat(len(js), M, [tc[(n * (2 * j + 1)) % P4] for j in js for n in range(1, M + 1)])
        Sm = arb_mat(len(js), M, [ts[(n * (2 * j + 1)) % P4] for j in js for n in range(1, M + 1)])
        F = rows_of(Cm * co); F1 = rows_of(Sm * si)
        for r_, j in enumerate(js):
            f = 1 + F[r_][0]; f1 = F1[r_][0]; f2 = F[r_][1]
            u = _ub(abs(f) + abs(f1) * dlt + abs(f2) * dlt ** 2 / 2 + S4 * dlt ** 3 / 6)
            cells.append(((arb(2 * j + 1) / (2 * Ng)) * arb.pi(), dlt, u))
            l_ = abs(f).lower()
            if l_ > lo: lo = l_
            if float(abs(f).mid()) > fmaxf: fmaxf = float(abs(f).mid()); jmax = j
    # branch and bound: keep cells whose upper bound reaches the best lower bound, bisect, evaluate directly
    def ev(a):
        cs = [(a * n).cos() for n in range(1, M + 1)]; ss = [(a * n).sin() for n in range(1, M + 1)]
        f = 1 + sum((xarb[n - 1] * n * cs[n - 1] for n in range(1, M + 1)), A0)
        f1 = -sum((xarb[n - 1] * (n * n) * ss[n - 1] for n in range(1, M + 1)), A0)
        f2 = -sum((xarb[n - 1] * (n * n * n) * cs[n - 1] for n in range(1, M + 1)), A0)
        return f, f1, f2
    for rnd in range(40):
        cells = [c for c in cells if not (c[2] < lo)]
        up = max((c[2] for c in cells), key=lambda z: float(z.mid()))
        if float((up - lo).mid()) < 1e-13 * float(up.mid()) or len(cells) > 64: break
        new = []
        for (cm, dl, u) in cells:
            for sgn in (-1, 1):
                a = arb((cm + sgn * dl / 2).mid()); d2 = dl / 2 + abs(a - cm - sgn * dl / 2)
                f, f1, f2 = ev(a)
                new.append((a, d2, _ub(abs(f) + abs(f1) * d2 + abs(f2) * d2 ** 2 / 2 + S4 * d2 ** 3 / 6)))
                if abs(f).lower() > lo: lo = abs(f).lower()
        cells = new
    up = max((c[2] for c in cells), key=lambda z: float(z.mid()))
    Kup = _ub(up + e1); Klo = lo - e1
    kup = (2 * arb.pi() * Kup / L.lower()).upper(); klo = (2 * arb.pi() * Klo / L.upper()).lower()
    return klo, kup, float((arb(2 * jmax + 1) / (2 * Ng) * arb.pi()).mid())


def resize_M(st, fm_old, M_new):
    """change the number of base modes: new modes enter with centre 0 (their error stays in the box);
    dropped modes move to the box (|xbar_n| + |G_n| w + box_n); the componentwise range Nc = 4M follows."""
    M = fm_old.M; Nc, Nc2 = 4 * M, 4 * M_new; d = st.G.shape[1]
    wA = [arb(float(x)) for x in st.w]
    if M_new >= M:
        xbar = np.concatenate([st.xbar, np.zeros(M_new - M)]); xa = list(st.xbar_arb) + [A0] * (M_new - M)
        G = np.vstack([st.G, np.zeros((M_new - M, d))])
        n = np.arange(Nc + 1, Nc2 + 1, dtype=float)
        box = np.concatenate([st.box, st.T * n ** (-QT) * (1 + 1e-12)]); T = st.T
    else:
        box = st.box.copy()
        for i in range(M_new, M):
            ex = abs(st.xbar_arb[i]) + sum((abs(arb(float(st.G[i, j]))) * wA[j] for j in range(d)), A0)
            box[i] = float(_ub(ex + arb(float(box[i]))).mid()) * (1 + 1e-12)
        tl = [float(_ub(arb(float(box[i])) * arb(i + 1) ** QT).mid()) for i in range(Nc2, Nc)]
        T = max(st.T, max(tl)) * (1 + 1e-12)
        box = box[:Nc2]; xbar = st.xbar[:M_new].copy(); xa = list(st.xbar_arb[:M_new]); G = st.G[:M_new].copy()
    new = SlabState(xbar, G, st.w.copy(), box, T, st.tau, st.Lint, xa, st.L, st.t, st.X0)
    return new


def mat_from_rows(rows):
    return arb_mat(len(rows), len(rows[0]), [x for r in rows for x in r])


def rows_of(m):
    return m.tolist()


def arb_vec(x):
    return [arb(float(v)) for v in x]


# ============================================================ arb series products (sigma/gamma two-sided forms)
def conv(a, b):
    c = (arb_poly(a) * arb_poly(b)).coeffs()
    L = len(a) + len(b) - 1
    if len(c) < L: c = c + [A0] * (L - len(c))
    return c


def sig2(u):
    return [-x for x in reversed(u)] + [A0] + list(u)


def node_fields(u):
    """base u (M sine coeffs) -> gv (gamma of v, len 2M+1), sw (sigma of w), gv2, sW (sigma of W), gc (gamma of 2v+v^2), len 4M+1"""
    M = len(u); s = sig2(u)
    gv = [s[i] * ((i - M) * 0.5) for i in range(2 * M + 1)]
    sw = [s[i] * (-(i - M) * (i - M)) for i in range(2 * M + 1)]
    gv2 = conv(gv, gv)
    gw2 = conv(sw, sw)
    sW = [A0 if i == 2 * M else gw2[i] / (-2 * (i - 2 * M)) for i in range(4 * M + 1)]
    gc = list(gv2)
    for i in range(2 * M + 1): gc[M + i] = gc[M + i] + 2 * gv[i]
    return gv, sw, gv2, sW, gc


def N_arb(u, fld=None):
    M = len(u)
    gv, sw, gv2, sW, gc = fld if fld is not None else node_fields(u)
    gv3 = conv(gv2, gv); gvsW = conv(gv, sW)
    out = []
    for k in range(1, 3 * M + 1):
        t = gv3[3 * M + k] / 3
        if k <= 2 * M: t = t + gv2[2 * M + k]
        t = t * (2 * k) + gvsW[3 * M + k]
        if k <= 2 * M: t = t + sW[2 * M + k]
        out.append(t)
    return out


def J_arb(fld, p):
    """J(u)p for base fields fld=(gv,sw,_,sW,gc) (two-sided lists) and p (M sine coeffs) -> 3M sine coeffs"""
    gv, sw, _, sW, gc = fld
    M = len(p); s = sig2(p)
    ge = [s[i] * ((i - M) * 0.5) for i in range(2 * M + 1)]
    se = [s[i] * (-(i - M) * (i - M)) for i in range(2 * M + 1)]
    Mc = (len(gc) - 1) // 2; Mw = (len(sW) - 1) // 2; Mv = (len(gv) - 1) // 2; Mws = (len(sw) - 1) // 2
    A = conv(gc, ge); oA = Mc + M
    B = conv(ge, sW); oB = M + Mw
    C = conv(sw, se); oC = Mws + M
    Y = [A0 if i == oC else -C[i] / (i - oC) for i in range(len(C))]
    Dd = conv(gv, Y); oD = Mv + oC
    out = []
    for k in range(1, 3 * M + 1):
        t = A[oA + k] * (2 * k) if k <= oA else A0
        if k <= oB: t = t + B[oB + k]
        if k <= oC: t = t + Y[oC + k]
        if k <= oD: t = t + Dd[oD + k]
        out.append(t)
    return out


# ============================================================ majorant sequences (nonneg, even, one-sided m=0..Nx)
class MS:
    """a(m) <= v[m] (m<=Nx);  a(m) <= F m^-p (m>Nx).  backend 'f' (numpy) or 'a' (exact arb upper bounds)."""

    def __init__(self, v, F=0.0, p=math.inf, be='f'):
        self.v, self.F, self.p, self.be = v, F, p, be

    @property
    def Nx(self): return len(self.v) - 1


def _ub(x): return x.abs_upper()


def ms_from(vals, be, F=0.0, p=math.inf):
    if be == 'f': return MS(np.asarray(vals, float), float(F), p, 'f')
    return MS([_ub(arb(x)) if not isinstance(x, arb) else _ub(x) for x in vals], _ub(arb(F)) if not isinstance(F, arb) else _ub(F), p, 'a')


def ms_l1(a):
    if a.be == 'f':
        s = a.v[0] + 2 * a.v[1:].sum()
        if a.F > 0: s += 2 * a.F * a.Nx ** (1 - a.p) / (a.p - 1)
        return s
    s = a.v[0] + 2 * sum(a.v[1:], A0)
    if a.F != 0: s = s + 2 * a.F * arb(a.Nx) ** (1 - a.p) / (a.p - 1)
    return _ub(s)


def ms_farsup(a):
    if (a.F == 0) if a.be == 'f' else a.F.is_zero(): return 0.0 if a.be == 'f' else A0
    if a.be == 'f': return a.F * (a.Nx + 1.0) ** (-a.p)
    return _ub(a.F * arb(a.Nx + 1) ** (-a.p))


def ms_env(a, R, p):
    """sup_{m>R} a(m) m^p (p <= a.p)"""
    R = int(R)
    if a.be == 'f':
        m = np.arange(R + 1, a.Nx + 1, dtype=float)
        e = (a.v[R + 1:] * m ** p).max() if len(m) else 0.0
        if a.F > 0: e = max(e, a.F * (a.Nx + 1.0) ** (p - a.p))
        return e
    e = A0
    for m in range(R + 1, a.Nx + 1):
        if not a.v[m].is_zero(): e = e.max(_ub(a.v[m] * arb(m) ** p))
    if not a.F.is_zero(): e = e.max(_ub(a.F * arb(a.Nx + 1) ** (p - a.p)))
    return e


def ms_conv(a, b, Nout):
    be = a.be
    if be == 'f':
        ta = np.concatenate([a.v[:0:-1], a.v]); tb = np.concatenate([b.v[:0:-1], b.v])
        c = np.convolve(ta, tb)[a.Nx + b.Nx:]
    else:
        ta = a.v[:0:-1] + a.v; tb = b.v[:0:-1] + b.v
        cc = conv(ta, tb); c = cc[a.Nx + b.Nx:]
    p = min(a.p, b.p)
    if p == math.inf:                       # finite support: keep everything (trim trailing zeros)
        L = len(c)
        while L > 1 and ((c[L - 1] == 0) if be == 'f' else c[L - 1].is_zero()): L -= 1
        c = c[:max(L, Nout + 1)]
    if len(c) < Nout + 1:
        c = (np.concatenate([c, np.zeros(Nout + 1 - len(c))]) if be == 'f' else c + [A0] * (Nout + 1 - len(c)))
    if p < math.inf: c = c[:Nout + 1]
    U = ms_farsup(a) * ms_l1(b) + ms_farsup(b) * ms_l1(a)
    if be == 'f':
        c = c + U
    else:
        U = _ub(U); c = [_ub(x + U) for x in c]
    if p == math.inf:
        F = 0.0 if be == 'f' else A0
    else:
        F = 2 ** p * (ms_env(a, Nout / 2, p) * ms_l1(b) + ms_env(b, Nout / 2, p) * ms_l1(a))
        if be == 'a': F = _ub(F)
    return MS(c, F, p, be)


def ms_pow(a, j, fac=1):
    """m^j * fac * a  (j in {-1,1,2}); m=0 entry set to 0 for j<0"""
    if a.be == 'f':
        m = np.arange(a.Nx + 1, dtype=float); mm = m.copy(); mm[0] = 1
        v = a.v * mm ** j * fac
        if j < 0: v[0] = 0
        return MS(v, a.F * fac, a.p - j, 'f')
    v = []
    for m, x in enumerate(a.v):
        if j < 0: v.append(A0 if m == 0 else _ub(x * fac / m))
        else: v.append(_ub(x * (m ** j) * fac))
    return MS(v, _ub(a.F * fac), a.p - j, 'a')


def ms_add(*ss):
    be = ss[0].be; N = max(s.Nx for s in ss); p = min(s.p for s in ss)
    if be == 'f':
        v = np.zeros(N + 1); F = 0.0
        for s in ss:
            v[:s.Nx + 1] += s.v
            if s.Nx < N and s.F > 0: v[s.Nx + 1:] += s.F * np.arange(s.Nx + 1, N + 1, dtype=float) ** (-s.p)
            if s.F > 0: F += s.F * (N + 1.0) ** (p - s.p)
        return MS(v, F, p, 'f')
    v = [A0] * (N + 1); F = A0
    for s in ss:
        for m in range(s.Nx + 1): v[m] = v[m] + s.v[m]
        if not s.F.is_zero():
            for m in range(s.Nx + 1, N + 1): v[m] = v[m] + s.F * arb(m) ** (-s.p)
            F = F + s.F * arb(N + 1) ** (p - s.p)
    return MS([_ub(x) for x in v], _ub(F), p, 'a')


def ms_trunc(a, N):
    """restrict explicit part to m<=N (moving the rest into the far bound)"""
    if a.Nx <= N: return a
    if a.be == 'f':
        m = np.arange(N + 1, a.Nx + 1, dtype=float)
        F = max(a.F * (a.Nx + 1.0) ** 0, (a.v[N + 1:] * m ** a.p).max()) if a.p < math.inf else 0.0
        return MS(a.v[:N + 1].copy(), max(F, a.F), a.p, 'f')
    F = ms_env(a, N, a.p) if a.p < math.inf else A0
    return MS(a.v[:N + 1], _ub(F.max(a.F)), a.p, 'a')


def ms_extend(a, N):
    """make explicit up to N using the far bound"""
    if a.Nx >= N: return a
    if a.be == 'f':
        m = np.arange(a.Nx + 1, N + 1, dtype=float)
        return MS(np.concatenate([a.v, a.F * m ** (-a.p) if a.F > 0 else 0 * m]), a.F, a.p, 'f')
    ext = [(_ub(a.F * arb(m) ** (-a.p)) if not a.F.is_zero() else A0) for m in range(a.Nx + 1, N + 1)]
    return MS(a.v + ext, a.F, a.p, 'a')


class BaseSup:
    """sup-over-slab coefficient sequences of the base (one-sided):
    GV=|gamma v|, SW=|sigma w|, GC=|gamma c|, SWW=|sigma W|; *_r: remainders beyond time degree KJ"""
    pass


def J_major(bs, y, Nout, M, rem=False):
    """majorant of |J e| given |e|<=y (MS), for J (rem=False) or J - J_trunc (rem=True)"""
    E = ms_pow(y, 1, 0.5); H = ms_pow(y, 2)
    Ni = Nout + 2 * M
    if not rem:
        t1 = ms_pow(ms_conv(bs.GC, E, Nout), 1, 2)
        t2 = ms_conv(E, bs.SWW, Nout)
        Yv = ms_pow(ms_conv(bs.SW, H, Ni), -1)
        t3 = ms_conv(bs.G1V, Yv, Nout)
        return ms_add(t1, t2, t3)
    t1 = ms_pow(ms_conv(bs.GC_r, E, Nout), 1, 2)
    t2 = ms_conv(E, bs.SWW_r, Nout)
    Yr = ms_pow(ms_conv(bs.SW_r, H, Ni), -1)
    Yf = ms_pow(ms_conv(bs.SW, H, Ni), -1)
    t3 = ms_conv(bs.G1V, Yr, Nout); t4 = ms_conv(bs.GV_r, Yf, Nout)
    return ms_add(t1, t2, t3, t4)


def Q_major(bs, y, Nout, M):
    E = ms_pow(y, 1, 0.5); H = ms_pow(y, 2)
    Ni = Nout + 2 * M
    EE = ms_conv(E, E, Ni)
    ta = ms_pow(ms_add(ms_conv(bs.G1V, EE, Nout), ms_pow(ms_conv(EE, E, Nout), 0, 1.0 / 3 if y.be == 'f' else arb(1) / 3)), 1, 2)
    Yv = ms_pow(ms_conv(bs.SW, H, Ni), -1)
    tb = ms_conv(E, Yv, Nout)
    Z2 = ms_pow(ms_conv(H, H, Ni), -1, 0.5)
    tc = ms_conv(bs.G1V, Z2, Nout); td = ms_conv(E, Z2, Nout)
    return ms_add(ta, tb, tc, td)


def bs_convert(bs, be):
    out = BaseSup()
    for k, v in bs.__dict__.items():
        if isinstance(v, MS):
            if be == 'f': setattr(out, k, MS(np.array([float(x.mid()) + float(x.rad()) for x in v.v]), 0.0, math.inf, 'f'))
            else: setattr(out, k, v)
        else: setattr(out, k, v)
    return out


# ============================================================ the slab
class SlabState:
    """chainable set: u = xbar + G rho + beta,  |rho_i|<=w_i,  |beta_n|<=box_n (n<=Nc), <= T n^-q (n>Nc)."""

    def __init__(self, xbar, G, w, box, T, tau=0.0, Lint=None, xbar_arb=None, L=None, t=None, X0=None):
        self.xbar, self.G, self.w, self.box, self.T, self.tau = xbar, G, w, box, T, tau
        self.Lint = Lint; self.xbar_arb = xbar_arb
        self.L, self.t, self.X0 = L, t, X0        # arb enclosures of length, time and base point x0 (X0 = (x0, 0))


def slaving_quasistatic(xbar, Ma, fm):
    M = fm.M
    flds = fm.fields(xbar[None, :])
    JE = fm.Japply(flds, np.eye(M)[None])[0].T          # J[n,m]
    return np.linalg.solve(np.diag(fm.lam[Ma:]) - JE[Ma:, Ma:], JE[Ma:, :Ma])


def slaving_propagated(xprev, fm, h, Ma, K=16, KP=16, Mf=60):
    """slaved high part of the block generators at the END of a (preliminary) slab [0,h] started at xprev:
    float propagator P with P(0) = [I; S_qs], S = P_r(h) P_a(h)^-1 (consistent with the flow)."""
    Ms = int(np.sum(fm.lam * h <= 8.0 * K * K))
    U = base_float(fm, xprev, h, K, np.tile(xprev, (K + 1, 1)), Mf, Ms)
    Cf = vinv_f(K) @ U
    xP, _ = cheb(KP); UP = teval_f(K, xP) @ Cf
    MsP = int(np.sum(fm.lam * h <= 8.0 * KP * KP))
    colP = Colloc(fm, h, KP, UP, Mf, MsP)
    G = np.vstack([np.eye(Ma), slaving_quasistatic(xprev, Ma, fm)])
    rhs = np.zeros((KP + 1, Ma, fm.M)); rhs[KP, :, :MsP] = G[:MsP, :].T
    sol, _ = colP.solve_blk(rhs)
    Ph = sol[0].T                                        # node 0 is s = h; (M, Ma)
    return Ph[Ma:] @ np.linalg.inv(Ph[:Ma])


def initial_state(xbar, Ma, Nc, eps0, fm, S0=None):
    """exact data xbar + ball |e_n| <= eps0 n^-q; generators [I; S0] (S0: slaved high part, default quasi-static)."""
    M = fm.M
    if S0 is None: S0 = slaving_quasistatic(xbar, Ma, fm)
    G = np.vstack([np.eye(Ma), S0])
    n = np.arange(1, Nc + 1, dtype=float)
    w = eps0 * n[:Ma] ** (-QT)
    box = eps0 * n ** (-QT)
    box[:Ma] = 0.0
    box[Ma:M] += np.abs(S0) @ w
    return SlabState(xbar.copy(), G, w, box * (1 + 1e-12), eps0)


def validate_slab(st, fm, h, K=16, KP=16, KJ=6, Ma=None, Mf=None, msub=8, Uguess=None, verbose=True, timing=None, lamcut=None,
                   grade=(1 / 512, 1 / 64), etacap_in=None, plan_next=None, kmax_ng=0):
    """validate one slab [0,h] in tau'.  Returns (new_state, info).
    plan_next(u_h_float) -> (h_next, Ma_next): block size of the NEXT slab (block resizing at the Lohner step)."""
    T0 = time.time(); tm = {} if timing is None else timing
    if etacap_in is None: validate_slab._depth = 0
    M = fm.M; Nc = 4 * M; d = st.G.shape[1]
    Ma = d if Ma is None else Ma
    Mf = Mf or min(max(Ma, 60), 140)
    harb = arb(h)
    lam_arb = [arb(n) ** 4 - arb(n) ** 2 for n in range(1, Nc + 2)]

    # ---------------- 1. base: float collocation + arb polish
    t = time.time()
    lamcut = 8.0 * K * K if lamcut is None else lamcut   # slaved start only where well conditioned: lam h > 8K^2
    Ms = int(np.sum(fm.lam * h <= lamcut))           # modes with imposed initial value
    x0arb = st.xbar_arb if getattr(st, 'xbar_arb', None) is not None else [arb(float(v)) for v in st.xbar]
    Ua, C, Crows, pol, U = solve_base(fm, x0arb, st.xbar, h, K, Mf, Ms, Uguess)
    tm['base'] = time.time() - t
    Varb = cheb_arb(K)
    u_h = [sum((Crows[k][n] for k in range(K + 1)), A0) for n in range(M)]
    u_0 = [sum((Crows[k][n] * (-1) ** k for k in range(K + 1)), A0) for n in range(M)]
    Cf = np.array([[float(x.mid()) for x in r] for r in Crows])

    # ---------------- 2. rigorous residual (degree 3K)
    t = time.time()
    Kf = 3 * K
    Ef = teval_arb(K, Kf); Edf = teval_arb(K, Kf, True)
    Uf = rows_of(Ef * C); dUf = rows_of(Edf * C)
    R = []
    for i in range(Kf + 1):
        Ni = N_arb(Uf[i])
        R.append([(dUf[i][n] * 2 / harb + Uf[i][n] * lam_arb[n] if n < M else A0) - Ni[n] for n in range(3 * M)])
    if DEBUG:
        import pickle
        pickle.dump(dict(R=[[float(x.mid()) for x in r] for r in R], pol=pol, x=[math.cos(math.pi*i/Kf) for i in range(Kf+1)]), open('/tmp/resid_dump.pkl','wb'))
    RcM = cheb_arb(Kf) * mat_from_rows(R)
    Rc = rows_of(RcM)
    rhat = [_ub(sum((abs(Rc[k][n]) for k in range(Kf + 1)), A0)) for n in range(3 * M)]
    # sub-intervals (graded towards s=0: initial layer) and local sup bounds by re-expansion
    edges = sorted(set([0.0] + [grade[i] for i in range(len(grade))] + [j / msub for j in range(1, msub + 1)]))
    subs = list(zip(edges[:-1], edges[1:]))
    rhat_sub = [sup_sub(RcM, Kf, a_, b_) for a_, b_ in subs]
    tm['resid'] = time.time() - t

    # ---------------- 3. coefficient functions at 2K+1 nodes, sup sequences, truncation at KJ
    t = time.time()
    Kc = 2 * K
    Uc = rows_of(teval_arb(K, Kc) * C)
    flds = [node_fields(Uc[i]) for i in range(Kc + 1)]
    Vc = cheb_arb(Kc)
    bs = BaseSup(); trunc = {}
    Lq = []
    for name, idx, half in (('GV', 0, M), ('SW', 1, M), ('GC', 4, 2 * M), ('SWW', 3, 2 * M)):
        L2 = 2 * half + 1
        coef = rows_of(Vc * mat_from_rows([flds[i][idx] for i in range(Kc + 1)]))
        full = [_ub(sum((abs(coef[k][half + m]) for k in range(Kc + 1)), A0)) for m in range(half + 1)]
        remv = [_ub(sum((abs(coef[k][half + m]) for k in range(KJ + 1, Kc + 1)), A0)) for m in range(half + 1)]
        setattr(bs, name, MS(full, A0, math.inf, 'a')); setattr(bs, name + '_r', MS(remv, A0, math.inf, 'a'))
        trunc[name] = coef[:KJ + 1]
    g1 = list(bs.GV.v); g1[0] = _ub(g1[0] + 1)
    bs.G1V = MS(g1, A0, math.inf, 'a')
    # quadrature: int_0^h 1/2 sum n^4 u_n^2 ds  (degree 2K, exact on 2K+1 nodes)
    gq = [sum((Uc[i][n] * Uc[i][n] * (n + 1) ** 4 for n in range(M)), A0) / 2 for i in range(Kc + 1)]
    gqc = rows_of(Vc * arb_mat(Kc + 1, 1, gq))
    Lbase = sum((gqc[k][0] * arb(2) / (1 - k * k) for k in range(0, Kc + 1, 2)), A0) * harb / 2
    tm['coef'] = time.time() - t

    # ---------------- 4. propagator P (float collocation, degree KP) and arb defect
    t = time.time()
    xP, _ = cheb(KP)
    UP = (teval_f(K, xP) @ Cf)
    MsP = int(np.sum(fm.lam * h <= 8.0 * KP * KP))
    colP = Colloc(fm, h, KP, UP, Mf, MsP)
    h_next, a1 = (None, Ma)
    if plan_next is not None:
        h_next, a1 = plan_next(Cf.sum(axis=0))
    pas = max(a1 - d, 0)                                            # passenger columns (rho = 0): new slow coordinates
    Gall = st.G
    if pas:
        # passenger seed for new slow mode m: e_m + slaved rest, slaving SHIFTED by the column's own decay rate:
        # (Lam_r - lam_m - J_rr) S = J_rm   (a column dominated by e^{-lam_m s}; unshifted quasi-static slaving
        # leaves an O(1) relative inconsistency in the modes just above the block -> defect feedback)
        JE = fm.Japply(fm.fields(st.xbar[None, :]), np.eye(M)[None])[0].T
        Gp = np.zeros((M, pas))
        for i, m in enumerate(range(d, a1)):
            Ar = np.diag(fm.lam[a1:] - fm.lam[m]) - JE[a1:, a1:]
            Gp[m, i] = 1.0; Gp[a1:, i] = np.linalg.solve(Ar, JE[a1:, m])
        Gall = np.hstack([st.G, Gp])
    dd = Gall.shape[1]
    rhs = np.zeros((KP + 1, dd, M)); rhs[KP, :, :MsP] = Gall[:MsP, :].T
    sol, ghist = colP.solve_blk(rhs); gm = [len(ghist), ghist[-1]]
    Pn = sol.transpose(0, 2, 1)                                    # (KP+1, M, dd)
    Pk = np.einsum('kj,jnd->knd', vinv_f(KP), Pn)                 # Chebyshev coefficients, taken as exact
    tm['Pfloat'] = time.time() - t
    t = time.time()
    KD = 2 * KJ + KP          # (1+v_t) dinv(w_t eta) is bilinear in the base: degree 2KJ+KP
    trunc_at = {}
    EJ = teval_arb(KJ, KD)
    for name in ('GV', 'SW', 'GC', 'SWW'):
        trunc_at[name] = rows_of(EJ * mat_from_rows(trunc[name]))
    fldD = [(trunc_at['GV'][i], trunc_at['SW'][i], None, trunc_at['SWW'][i], trunc_at['GC'][i]) for i in range(KD + 1)]
    EP = teval_arb(KP, KD); EdP = teval_arb(KP, KD, True)
    VD = cheb_arb(KD)
    Dsup_a = []
    Pabs = [[A0] * d for _ in range(M)]       # sum_k |P_k| (sup over slab of |P_nj(s)|)
    EPmis = []                                # P(0) - G0
    for j in range(d):
        Pkj = mat_from_rows([arb_vec(Pk[k, :, j]) for k in range(KP + 1)])
        Pj = rows_of(EP * Pkj); dPj = rows_of(EdP * Pkj)
        Dv = []
        for i in range(KD + 1):
            Jp = J_arb(fldD[i], Pj[i])
            Dv.append([(dPj[i][n] * 2 / harb + Pj[i][n] * lam_arb[n] if n < M else A0) - Jp[n] for n in range(3 * M)])
        DcM = VD * mat_from_rows(Dv); Dc = rows_of(DcM)
        Dsup_a.append([_ub(sum((abs(Dc[k][n]) for k in range(KD + 1)), A0)) for n in range(3 * M)])

        pr = rows_of(Pkj)
        for n in range(M):
            Pabs[n][j] = _ub(sum((abs(pr[k][n]) for k in range(KP + 1)), A0))
        EPmis.append([_ub(abs(sum((pr[k][n] * (-1) ** k for k in range(KP + 1)), A0) - float(st.G[n, j]))) for n in range(M)])
    tm['defect'] = time.time() - t

    # ---------------- 5. forcing, a priori bounds on sub-intervals
    t = time.time()
    wA0 = [_ub(arb(float(x))) for x in st.w]
    # radii used inside the slab: w + cap for the intra-slab absorptions (verified below)
    wmax = max((float(x.mid()) for x in wA0), default=0.0)
    etacap = [_ub(arb(wmax * ETACAP + 1e-40)) for x in wA0] if etacap_in is None else etacap_in
    wA = [_ub(wA0[j] + etacap[j]) for j in range(d)]
    eP = [_ub(sum((Pabs[n][j] * wA[j] for j in range(d)), A0)) for n in range(M)]
    eP_ms = MS(eP + [A0] * (Nc - M), A0, math.inf, 'a')
    Jrem = J_major(bs, MS(eP + [A0], A0, math.inf, 'a'), 3 * M, M, rem=True)
    dhat = [_ub(sum((Dsup_a[j][n] * wA[j] for j in range(d)), A0) + Jrem.v[n + 1]) for n in range(3 * M)]
    f0_sub = [[_ub(rhat_sub[i][n] + dhat[n]) if n < 3 * M else A0 for n in range(Nc)] for i in range(len(subs))]
    a = []
    for n in range(Nc):
        x = arb(float(st.box[n]))
        if n < M:
            x = x + abs(x0arb[n] - u_0[n]) + sum((EPmis[j][n] * wA0[j] for j in range(d)), A0)
        a.append(_ub(x))
    aT = _ub(arb(float(st.T)))
    bsf = bs_convert(bs, 'f')
    tailfac = arb(Nc + 1) ** QT / lam_arb[Nc]     # times F n^-(q-2) far forcing:  n^q phi_n F n^{2-q} <= F (Nc+1)^2/lam
    Zsup = [A0] * Nc; TZsup = A0; subinfo = []; eta_tot = [A0] * d; absorbed = []
    pass
    for sub, (sa_, sb_) in enumerate(subs):
        Delta = harb * (arb(sb_) - arb(sa_))
        ex = [(-lam_arb[n] * Delta).exp() for n in range(Nc + 1)]
        phi = [(-(-lam_arb[n] * Delta).expm1() / lam_arb[n]) if n > 0 else Delta for n in range(Nc + 1)]
        exf = np.array([float(x.mid()) for x in ex]); phif = np.array([float(x.mid()) for x in phi])
        f0 = f0_sub[sub]
        # float iteration (margin mu in every mode); Q (quadratic) enters as fixed forcing after a failed check
        af = np.array([float(x.mid()) for x in a]); aTf = float(aT.mid())
        f0f = np.array([float(x.mid()) for x in f0]); ePf = np.array([float(x.mid()) for x in eP_ms.v])
        Qx = np.zeros(Nc); QxT = 0.0
        Z = af + phif[:Nc] * f0f; TZ = aTf
        for vround in range(5):
          for it in range(40):
            y = ms_extend(MS(np.concatenate([[0.0], Z]), TZ, QT, 'f'), Nc + 2 * M)
            Fm = J_major(bsf, y, Nc, M)
            Fv = Fm.v[1:] + f0f + Qx
            Zn = np.maximum(af, exf[:Nc] * af + phif[:Nc] * Fv) / (1 - MU)
            TZn = (aTf + (Fm.F + QxT) * (Nc + 1.0) ** (QT - Fm.p) / float(lam_arb[Nc].mid())) / (1 - MU)
            ch = np.max(np.abs(Zn - Z) / np.maximum(Zn, 1e-300))
            Z, TZ = np.maximum(Zn, Z) if vround else Zn, max(TZn, TZ) if vround else TZn
            if ch < 1e-5: break
          Zs = Z * (1 + 1e-6) + 1e-250; TZs = TZ * (1 + 1e-6) + 1e-250
          # arb verification (strict)
          if True:
            Zsa = [arb(float(x)) for x in Zs]; TZsa = arb(float(TZs))
            y = ms_extend(MS([A0] + Zsa, TZsa, QT, 'a'), Nc + 2 * M)
            ye = ms_extend(MS([A0] + [Zsa[n] + eP_ms.v[n] for n in range(Nc)], TZsa, QT, 'a'), Nc + 2 * M)
            Qa = Q_major(bs, ye, Nc, M)
            Fa = ms_add(J_major(bs, y, Nc, M), Qa)
            Fv = [Fa.v[n + 1] + f0[n] for n in range(Nc)]
            rhs = [a[n].max(ex[n] * a[n] + phi[n] * Fv[n]) for n in range(Nc)]
            if Fa.p < QT - 2: raise ValueError('tail power')
            tail_rhs = aT + Fa.F * arb(Nc + 1) ** (QT - Fa.p) / lam_arb[Nc]
            bad = [n for n in range(Nc) if not (rhs[n] < Zsa[n])]
            tbad = not (tail_rhs < TZsa)
            if not bad and not tbad: break
            Qx = np.maximum(Qx, 1.5 * np.array([float(x.mid()) for x in Qa.v[1:Nc + 1]])); QxT = max(QxT, 1.5 * float(Qa.F.mid()))
            if DEBUG or vround >= 2:
                print(f'   a priori sub {sub} round {vround}: {len(bad)} bad modes {bad[:8]} max rhs/Z '
                      f'{max(float(rhs[n].upper().mid()) / Zs[n] for n in bad) if bad else 0:.4f}; phi*Q/Z at bad '
                      f'{max(float(phi[n].mid()) * float(Qa.v[n + 1].mid()) / Zs[n] for n in bad) if bad else 0:.2e}; tail {tbad}')
            if DEBUG:
                rr = [float(rhs[n].upper().mid()) / Zs[n] for n in bad]
                print(f'   sub {sub} round {vround}: float it {it} ch {ch:.1e}; bad {len(bad)} modes {bad[:12]}... rhs/Z max {max(rr):.3e}; tail {tbad}')
            for n in bad: Zs[n] = max(Zs[n], float(rhs[n].upper().mid()) * 1.05 + 1e-250)
            if tbad: TZs = max(TZs, float(tail_rhs.upper().mid()) * 1.05)
            Z = np.maximum(Z, Zs); TZ = max(TZ, TZs)
        else:
            raise RuntimeError(f'a priori bound not verified on sub-interval {sub}: {len(bad)} modes, tail {tbad}')
        fl = lambda L, i0, i1: max(float(x.mid()) for x in L[i0:i1])
        subinfo.append((it, vround, len(bad)))
        if DEBUG:
            print(f'   sub {sub} [{sa_:.4f},{sb_:.4f}]: a_blk {fl(a,0,Ma):.2e} a_rest {fl(a,Ma,Nc):.2e} | f0_blk {fl(f0,0,Ma):.2e} f0_rest {fl(f0,Ma,3*M):.2e} '
                  f'| r_blk {fl(rhat_sub[sub],0,Ma):.2e} d_blk {fl(dhat,0,Ma):.2e} Jrem_blk {max(float(x.mid()) for x in Jrem.v[1:Ma+1]):.2e} | Z_blk {max(Zs[:Ma]):.2e} Z_rest {max(Zs[Ma:]):.2e}'
                  f' | J_blk {max(Fm.v[1:Ma + 1]):.2e} (n={int(np.argmax(Fm.v[1:Ma + 1])) + 1}) Q_blk {max(float(x.mid()) for x in Qa.v[1:Ma + 1]):.2e} argmax a_rest n={Ma + int(np.argmax([float(x.mid()) for x in a[Ma:Nc]])) + 1}')
        for n in range(Nc): Zsup[n] = Zsup[n].max(Zsa[n])
        TZsup = TZsup.max(TZsa)
        # end-of-sub-interval bounds
        a = [_ub(ex[n] * a[n] + phi[n] * Fv[n]) for n in range(Nc)]
        aT = _ub(ex[Nc] * aT + Fa.F * arb(Nc + 1) ** (QT - Fa.p) / lam_arb[Nc])
        # (ex[Nc] = exp(-lam_{Nc+1} Delta) is the max over n > Nc)
        if sub < len(subs) - 1 and d > 0:
            # intra-slab absorption at s_i: z_a = P_a(s_i) eta, rho <- rho + eta, z_r <- z_r - P_r(s_i) eta
            xi = 2 * arb(sb_) - 1
            Tk = [arb(1), xi]
            for k in range(2, KP + 1): Tk.append(2 * xi * Tk[-1] - Tk[-2])
            Ps = [[sum((Tk[k] * float(Pk[k, n, j]) for k in range(KP + 1)), A0) for j in range(d)] for n in range(M)]
            Ai = mat_from_rows(Ps[:Ma]).inv()                         # eta_i = A_i zeta_i, |zeta_i| <= Z_a(s_i)
            Si = rows_of(mat_from_rows(Ps[Ma:]) * Ai)                  # slaving at s_i (product first: no cond loss)
            Za = [a[k] for k in range(Ma)]
            for r_ in range(M - Ma):
                a[Ma + r_] = _ub(a[Ma + r_] + sum((abs(Si[r_][k]) * Za[k] for k in range(Ma)), A0))
            for n in range(Ma): a[n] = A0
            Al = rows_of(Ai)
            eta_tot = [_ub(eta_tot[j] + sum((abs(Al[j][k]) * Za[k] for k in range(Ma)), A0)) for j in range(d)]
            absorbed.append((Ai, Za))
    tm['apriori'] = time.time() - t
    Zend = a; Tend = aT
    capok = all(eta_tot[j] < etacap[j] for j in range(d))
    if not capok:
        if (etacap_in is not None) and getattr(validate_slab, '_depth', 0) >= 3: raise RuntimeError('eta cap not met after 3 redos')
        validate_slab._depth = getattr(validate_slab, '_depth', 0) + 1
        newcap = [_ub(2 * eta_tot[j] + etacap[j]) for j in range(d)]
        print(f'   eta cap exceeded (max eta/cap {max(float(eta_tot[j].mid()) / float(etacap[j].mid()) for j in range(d)):.2e}): redo')
        if verbose: print('   eta cap exceeded: redo with larger cap')
        return validate_slab(st, fm, h, K=K, KP=KP, KJ=KJ, Ma=Ma, Mf=Mf, msub=msub, Uguess=Uguess, verbose=verbose,
                             timing=timing, lamcut=lamcut, grade=grade, etacap_in=newcap, plan_next=plan_next, kmax_ng=kmax_ng)
    wA = wA0                                                        # rho' = rho + sum_i A_i zeta_i (kept separately)

    # ---------------- 6. Lohner step
    t = time.time()
    # block resizing: new block = modes 1..a1 (a1 > d: passengers supply the new columns; a1 < d: wide QR, the
    # direction annihilated by rows 1..a1 goes to the box through the slaving mismatch Rmis)
    Ph = Pk.sum(axis=0)                                            # P(h), float approx of the exact sum (M, dd)
    Ph_rows = [[sum((arb(float(Pk[k, n, j])) for k in range(KP + 1)), A0) for j in range(dd)] for n in range(M)]
    w0f = np.asarray(st.w, float)
    wext = np.concatenate([w0f, np.full(pas, max(w0f.max(), 1e-300) * 1e-3)])
    Asc = Ph[:a1] * wext[None, :]
    order = np.argsort(-np.linalg.norm(Asc, axis=0))
    Qm, Rm = np.linalg.qr(Asc[:, order])
    Qm = Qm[:, :a1]
    Gr = Ph[a1:] @ (np.linalg.solve(Ph[:a1], Qm) if dd == a1 else np.linalg.pinv(Ph[:a1]) @ Qm)
    Gnew = np.vstack([Qm, Gr])
    Qa = arb_mat(a1, a1, [arb(float(x)) for x in Qm.ravel()])
    Qinv = Qa.inv()
    PhA = mat_from_rows(Ph_rows[:a1])
    Bm = Qinv * PhA                                                # a1 x dd; columns >= d (passengers) carry rho = 0
    Bl = rows_of(Bm); Ql = rows_of(Qinv)
    eta = [_ub(sum((abs(Ql[i][k]) * Zend[k] for k in range(a1)), A0)) for i in range(a1)]
    wnew = [_ub(sum((abs(Bl[i][j]) * wA[j] for j in range(d)), A0) + eta[i]) for i in range(a1)]
    Bd = mat_from_rows([r_[:d] for r_ in Bl])
    for (Ai, Za) in absorbed:
        BA = rows_of(Bd * Ai)
        wnew = [_ub(wnew[i] + sum((abs(BA[i][k]) * Za[k] for k in range(Ma)), A0)) for i in range(a1)]
    GrA = arb_mat(M - a1, a1, [arb(float(x)) for x in Gr.ravel()])
    RmisM = mat_from_rows([r_[:d] for r_ in Ph_rows[a1:]]) - GrA * Bd
    Rmis = rows_of(RmisM)
    Rabs = [A0] * (M - a1)
    for (Ai, Za) in absorbed:
        RA = rows_of(RmisM * Ai)
        Rabs = [_ub(Rabs[r_] + sum((abs(RA[r_][k]) * Za[k] for k in range(Ma)), A0)) for r_ in range(M - a1)]
    newbox = [A0] * Nc
    for r in range(M - a1):
        n = a1 + r
        s = sum((abs(Rmis[r][j]) * wA[j] for j in range(d)), A0) + Zend[n] + Rabs[r]
        s = s + sum((abs(float(Gr[r, i])) * eta[i] for i in range(a1)), A0)
        newbox[n] = s
    for n in range(M, Nc): newbox[n] = Zend[n]
    xnew_arb = [arb(x.mid()) for x in u_h]                         # exact 128-bit centre (consistent datum)
    xnew = np.array([float(x.mid()) for x in u_h])
    for n in range(M):
        newbox[n] = _ub(newbox[n] + u_h[n].rad())
    tm['lohner'] = time.time() - t

    # ---------------- 7. slab-wide outputs
    eh = [_ub(eP_ms.v[n] + Zsup[n]) for n in range(Nc)]
    tailsum = lambda j: TZsup * arb(Nc) ** (j + 1 - QT) / (QT - j - 1)
    supnorm = [float((sum((eh[n] * arb(n + 1) ** j for n in range(Nc)), A0) + tailsum(j)).upper().mid()) for j in range(4)]
    Ubar = [_ub(sum((abs(Crows[k][n]) for k in range(K + 1)), A0)) for n in range(M)]
    Lerr = harb * (sum((Ubar[n] * eh[n] * arb(n + 1) ** 4 for n in range(M)), A0)
                   + sum((eh[n] * eh[n] * arb(n + 1) ** 4 for n in range(Nc)), A0) / 2 + TZsup ** 2 * arb(Nc) ** (5 - 2 * QT))
    Lint = Lbase + _ub(Lerr) * EB
    # ---------------- 8. quadratures: L, t, X0 (slaved quantities), exact Chebyshev algebra + error terms
    t = time.time(); quadinfo = None; Lnew = tnew = X0new = None
    if st.L is not None:
        gco = [r_[0] for r_ in gqc]                                  # gbar(x) = 1/2 sum n^4 ubar_n^2 (degree 2K)
        Gb = [c * harb / 2 for c in cheb_antider(gco)]               # G_base(s) = int_0^s gbar >= 0
        Gmax = _ub(sum((abs(c) for c in Gb), A0))
        e3 = _ub(sum((eh[n] * arb(n + 1) ** 3 for n in range(Nc)), A0) + tailsum(3))
        Jx = 14
        def expser(fac):
            out = [arb(1)]; term = [arb(1)]
            for jj in range(1, Jx + 1):
                term = [c * fac / jj for c in cheb_mul(term, Gb)]
                out = cheb_add(out, term)
            x_ = abs(fac) * Gmax
            return out, _ub(x_ ** (Jx + 1) / arb.fac_ui(Jx + 1) * x_.exp())
        E4, R4 = expser(arb(-4))
        I4 = (cheb_intdef(E4) * harb / 2 + _ub(R4 * harb) * EB) * (4 * _ub(Lerr) * EB).exp()
        Lr0 = st.L / (2 * arb.pi())
        dt_ = Lr0 ** 4 * I4
        E1, R1 = expser(arb(-1))
        pk = [sum((Crows[k][n] * (n + 1) ** 3 for n in range(M)), A0) for k in range(K + 1)]
        supp = _ub(sum((abs(c) for c in pk), A0))
        mainX = cheb_intdef(cheb_mul(E1, pk)) * harb / 2
        errX = _ub(harb * (supp * (R1 + _ub(Lerr).exp() - 1) + e3))
        dX0 = Lr0 * (mainX + errX * EB)
        Lnew = st.L * (-Lint).exp(); tnew = st.t + dt_; X0new = st.X0 + dX0
        quadinfo = dict(L=Lnew.str(20), t=tnew.str(20), X0=X0new.str(20), dt=dt_.str(12), dX0=dX0.str(12),
                        Lrad=float(Lnew.rad()), trad=float(tnew.rad()), X0rad=float(X0new.rad()))
    tm['quad'] = time.time() - t
    newst = SlabState(xnew, Gnew, np.array([float(x.mid()) + float(x.rad()) for x in wnew]) * (1 + 1e-12),
                      np.array([float(x.upper().mid()) for x in newbox]) * (1 + 1e-12) + 1e-300,
                      float(Tend.upper().mid()) * (1 + 1e-12), st.tau + h, Lint, xnew_arb, Lnew, tnew, X0new)
    if kmax_ng and Lnew is not None:
        t = time.time()
        nn_ = np.arange(1, Nc + 1, dtype=float)
        e1end = arb(float(((np.abs(Gnew) * nn_[:M, None]).sum(axis=0) @ newst.w + (newst.box * nn_).sum()
                           + newst.T * Nc ** (2.0 - QT) / (QT - 2))) * (1 + 1e-9))
        klo, kup, amax = kmax_enclosure(xnew_arb, e1end, Lnew, kmax_ng)
        quadinfo.update(kmax=[float(klo.mid()), float(kup.mid())], alpha_kmax=amax, e1_end=float(e1end.mid()))
        tm['kmax'] = time.time() - t
    # ---------------- diagnostics
    n = np.arange(1, Nc + 1, dtype=float)
    f = lambda L: np.array([float(x.mid()) + float(x.rad()) for x in L])
    rh = f(rhat); dh = f(dhat); Ze = f(Zend); Zs_ = f(Zsup)
    info = dict(Ms=Ms, h=h, K=K, KP=KP, KJ=KJ, M=M, Ma=Ma, d=d, Nc=Nc, msub=msub, polish=pol, Prich=gm,
                resid_max=rh.max(), resid_w3=float((rh * np.arange(1, 3 * M + 1) ** 3).sum()),
                resid_ind=float((rh / np.maximum(f(lam_arb[:3 * M]), 1.0 / h) * np.arange(1, 3 * M + 1) ** 3).sum()),
                defect_w3=float((dh * np.arange(1, 3 * M + 1) ** 3).sum()),
                w_in=float(np.max(st.w)), w_out=float(max(f(wnew))),
                wnorm_in=float(np.abs(st.G).sum(axis=0) @ st.w), wnorm_out=float(np.abs(Gnew).sum(axis=0) @ f(wnew)),
                box_in_w3=float((st.box * n ** 3).sum()), box_out_w3=float((newst.box * n ** 3).sum()),
                Zend_a_w3=float((Ze[:Ma] * n[:Ma] ** 3).sum()), Zend_r_w3=float((Ze[Ma:] * n[Ma:] ** 3).sum()),
                eta_w3=float((f(eta)).sum()), eta_intra=float(max(f(eta_tot))) if d else 0.0,
                hull3_in=float((np.abs(st.G) * (n[:M, None] ** 3)).sum(axis=0) @ st.w),
                hull3_out=float((np.abs(Gnew) * (n[:M, None] ** 3)).sum(axis=0) @ newst.w),
                hull1_in=float((np.abs(st.G) * n[:M, None]).sum(axis=0) @ st.w),
                hull1_out=float((np.abs(Gnew) * n[:M, None]).sum(axis=0) @ newst.w),
                box1_out=float((newst.box * n).sum()),
                wB_diag=float(max(abs(float(Bl[i][i].mid())) * st.w[i] for i in range(min(a1, d)))), a1=a1, h_next=h_next, passengers=pas,
                wB_full=float(max(sum(abs(float(Bl[i][j].mid())) * st.w[j] for j in range(d)) for i in range(a1))),
                w_eta=float(max(f(eta))), sv_Pa=[float(x) for x in np.linalg.svd(Ph[:a1, :] , compute_uv=False)[[0, 1, 2, -1]]],
                Rdiag_ratio=float(np.max(np.abs(np.diag(Rm))) / max(np.max(w0f * 0 + 1e-300), 1e-300)),
                T_in=float(st.T), T_out=newst.T, supnorm=supnorm, Lint=(float(Lint.mid()), float(Lint.rad())),
                theta_end=float(np.abs(1 + np.sum(np.arange(1, M + 1) * xnew)) if False else fm.thmax(xnew[None])[0]),
                umis=float(max(abs(x0arb[n] - u_0[n]).upper().mid() for n in range(M))), Mslaved=M - Ms,
                cond_Pa=float(np.linalg.cond(Ph[:a1, :])), sub=subinfo, Cf=Cf, time=time.time() - T0, tm=tm, quad=quadinfo)
    if verbose:
        print(f"  slab h={h:.3e}: resid max {info['resid_max']:.2e} | w {info['w_in']:.2e}->{info['w_out']:.2e} "
              f"| box|3 {info['box_in_w3']:.2e}->{info['box_out_w3']:.2e} | T {info['T_in']:.1e}->{info['T_out']:.1e} | {info['time']:.0f}s {tm}")
    return newst, info
