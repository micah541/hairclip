"""
Rigorous verification of a numerical-range certificate for the linearisation L at S1 (step 2, a priori bounds).

Claim certified: for all x in R,  P(x) := H_th(x) + N_W(x) - Lam E33  is negative semidefinite, where for the jet
J = (u'', u', u) and rho = exp(-phi), phi = beta (1+x^2)^{2/3},
    Re(e^{i th} <L u, u>_rho) = int J^* H_th J rho,     int J^* N_W J rho = 0   (N_W a null form built from W).
Consequently Re(e^{i th} z) <= Lam for every z in the numerical range of L in L^2(rho), in particular for every
eigenvalue.  By the reflection symmetry (parities of the coefficients and of W) it suffices to take x >= 0.

H_11 = -cos th, H_21 = e^{i th} phi', H_31 = -e^{i th}(g + k^2)/2, H_32 = e^{i th} b/2, H_33 = cos th * d,
   g = phi'^2 - phi'', b = -(5 k m + q/4), d = 1/4 - 3 m^2 - k p   (profile (p,q,k,m) of S1 at arclength x from the tip)
N = [[0, w11, w12], [w11, w11' + 2 wc - phi' w11, w12' + w22 - phi' w12], [w21, w21' + w22 - phi' w21, w22' - phi' w22]],
   w12 = wc + i wd, w21 = conj(w12);  w11, w22 odd, wc, wd even (explicit basis functions, certificate coefficients).

Regions:  (1) core x in [0, x_S] (x_S = arclength from the tip to sigma = 8): CAPD enclosures of S1 (curve_br.txt);
          (2) tail sigma in [8, sigma_v]: step-1 tail theorem bounds;   (3) x >= X_v: u = x^{-1/3} in [0, X_v^{-1/3}],
              scaled matrix  Phat = u^4 D P D,  D = diag(u^-2, u^-1, 1), analytic in u up to u = 0.
NSD test on each box: midpoint matrix Mc of -P, Frobenius radius rho of the ball matrix; Sylvester's criterion for
Mc - rho I in ball arithmetic.
"""
import sys, json, math
from flint import arb, acb, ctx
sys.path.insert(0, '/home/claude/s1_cap')
from tail import radii

ctx.prec = 96
THIRD = arb(1) / 3


def ahull(lo, hi):
    lo, hi = arb(lo), arb(hi)
    l = lo.lower(); h = hi.upper()
    return arb((l + h) / 2) + arb(0, 1) * ((h - l) / 2)


class Cert:
    def __init__(self, path):
        d = json.load(open(path))
        self.d = d
        self.beta = arb(repr(d["beta"])); self.th = arb(repr(d["th"])); self.Lam = arb(repr(d["Lam"]))
        self.s = arb(repr(d["s"]))
        self.nodes = [arb(repr(c)) for c in d["nodes"]]
        self.t11 = [arb(repr(v)) for v in d["t11"]]; self.t22 = [arb(repr(v)) for v in d["t22"]]
        self.tc = [arb(repr(v)) for v in d["tc"]]; self.td = [arb(repr(v)) for v in d["td"]]
        self.e = acb(self.th.cos(), self.th.sin())
        self.cth = self.th.cos()

    # ---- basis (same ordering as sdp.py) ----
    def gauss(self, x, c):
        s2 = 2 * self.s * self.s
        v = (-(x - c) ** 2 / s2).exp()
        return v, -(x - c) / (self.s * self.s) * v

    def W(self, x):
        """(w11, w11', w22, w22', wc, wc', wd, wd') at x (arb ball)"""
        u = 1 + x * x
        f_odd1 = (x * u ** (-THIRD), u ** (-THIRD) - (arb(2) / 3) * x * x * u ** (-4 * THIRD))
        f_x = (x, arb(1))
        f_e2 = (u ** THIRD, (arb(2) / 3) * x * u ** (-2 * THIRD))
        f_1 = (arb(1), arb(0))
        odd = []; even = []
        for c in self.nodes:
            gp, dgp = self.gauss(x, c); gm, dgm = self.gauss(x, -c)
            if not (c == 0):
                odd.append((gp - gm, dgp - dgm))
            even.append((gp + gm, dgp + dgm))
        F11 = [f_odd1] + odd
        F22 = [f_x, f_odd1] + odd
        Fc = [f_e2, f_1] + even
        def comb(th, F):
            v = arb(0); dv = arb(0)
            for t, (f, df) in zip(th, F):
                v += t * f; dv += t * df
            return v, dv
        w11 = comb(self.t11, F11); w22 = comb(self.t22, F22); wc = comb(self.tc, Fc); wd = comb(self.td, Fc)
        return w11 + w22 + wc + wd

    def phi(self, x):
        u = 1 + x * x
        b = self.beta
        p1 = b * (arb(4) / 3) * x * u ** (-THIRD)
        p2 = b * (arb(4) / 3) * (u ** (-THIRD) - (arb(2) / 3) * x * x * u ** (-4 * THIRD))
        return p1, p2

    def P(self, x, p, q, k, m):
        """entries of P = H + N - Lam E33 at arclength x with profile balls (p, q, k, m)"""
        w11, dw11, w22, dw22, wc, dwc, wd, dwd = self.W(x)
        f1, f2 = self.phi(x)
        g = f1 * f1 - f2
        b = -(5 * k * m + q / 4)
        d = arb(1) / 4 - 3 * m * m - k * p
        e = self.e
        w12 = acb(wc, wd); w21 = acb(wc, -wd); dw21 = acb(dwc, -dwd)
        P11 = -self.cth
        P21 = e * f1 + w11
        P31 = -e * (g + k * k) / 2 + w21
        P22 = dw11 + 2 * wc - f1 * w11
        P32 = e * b / 2 + dw21 + w22 - f1 * w21
        P33 = self.cth * d + dw22 - f1 * w22 - self.Lam
        return P11, P21, P31, P22, P32, P33

    def Phat_u(self, u, sig_tip, tiny):
        """scaled matrix u^4 D P D (D = diag(u^-2, u^-1, 1)) for x = u^-3 >= X_v, cone coefficients;
        'tiny' = bound for all profile corrections and Gaussian terms (added as balls)"""
        b_ = self.beta
        u2 = u * u; u3 = u2 * u; u4 = u2 * u2; u6 = u3 * u3
        v = 1 + u6
        vm13 = v ** (-THIRD); vm23 = v ** (-2 * THIRD); vm43 = v ** (-4 * THIRD); v13 = v ** THIRD
        A1 = self.t11[0]; A2 = self.t22[0]; A3 = self.t22[1]
        C1 = self.tc[0]; C2 = self.tc[1]; D1 = self.td[0]; D2 = self.td[1]
        e = self.e
        T = acb(arb(0, 1) * tiny, arb(0, 1) * tiny)
        # scaled pieces
        uphi1 = (arb(4) / 3) * b_ * vm13                                  # u * phi'
        u2phi2 = (arb(4) / 3) * b_ * u4 * (vm13 - (arb(2) / 3) * vm43)    # u^2 * phi''
        u2g = uphi1 * uphi1 - u2phi2                                       # u^2 * g
        uw11 = A1 * vm13                                                   # u * w11
        u2dw11 = A1 * u4 * (vm13 - (arb(2) / 3) * vm43)                    # u^2 * w11'
        u3w22 = A2 + A3 * u2 * vm13                                        # u^3 * w22
        u4dw22 = A2 * u4 + A3 * u6 * (vm13 - (arb(2) / 3) * vm43)          # u^4 * w22'
        u2wc = C1 * v13 + C2 * u2; u2wd = D1 * v13 + D2 * u2               # u^2 * wc, u^2 * wd
        u3dwc = (arb(2) / 3) * C1 * u4 * vm23; u3dwd = (arb(2) / 3) * D1 * u4 * vm23
        u3b = -(1 + sig_tip * u3) / 4                                      # u^3 * b (cone: q = x + sig_tip)
        u4d = u4 / 4
        u2w21 = acb(u2wc, -u2wd)
        P11 = -self.cth
        P21 = e * uphi1 + uw11 + T
        P31 = -e * u2g / 2 + u2w21 + T
        P22 = u2dw11 + 2 * u2wc - uphi1 * uw11 + T.real
        P32 = e * u3b / 2 + acb(u3dwc, -u3dwd) + u3w22 - uphi1 * u2w21 + T
        P33 = self.cth * u4d + u4dw22 - uphi1 * u3w22 - self.Lam * u4 + T.real
        return P11, P21, P31, P22, P32, P33


def nsd_check(P11, P21, P31, P22, P32, P33):
    """rigorous: P (Hermitian, given by lower triangle; diagonal real) <= 0 ?  returns (ok, margin)"""
    M = [[-acb(P11), -acb(P21).conjugate(), -acb(P31).conjugate()],
         [-acb(P21), -acb(P22), -acb(P32).conjugate()],
         [-acb(P31), -acb(P32), -acb(P33)]]
    rad2 = arb(0)
    Mc = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            z = M[i][j]
            mr, mi = float(z.real.mid()), float(z.imag.mid())
            if i == j:
                mi = 0.0
            Mc[i][j] = acb(mr, mi)
            dr = abs(z.real - arb(mr)).upper(); di = abs(z.imag - arb(mi)).upper()
            rad2 += arb(dr) ** 2 + arb(di) ** 2
    rho = rad2.sqrt()
    A = [[Mc[i][j] - (rho if i == j else 0) for j in range(3)] for i in range(3)]
    m1 = A[0][0].real
    m2 = (A[0][0] * A[1][1] - A[0][1] * A[1][0]).real
    det = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
           + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])).real
    ok = bool(m1 > 0 and m2 > 0 and det > 0)
    return ok, (float(m1.lower()), float(m2.lower()), float(det.lower()), float(rho.upper()))


def read_curve(path):
    L = open(path).read().split("\n")
    tip = tuple(map(float, L[0].split()[1:3]))
    rows = []
    for l in L[6:]:
        if l.startswith("END") or not l.strip():
            break
        v = list(map(float, l.split()))
        rows.append(v)
    return tip, rows


def verify(cert_path, curve_path="/home/claude/s1spec/apr/curve_fine.txt", sig_v=40.0, dsig=0.02, nu=400, verbose=True, maxdepth=6):
    C = Cert(cert_path)
    fails = []
    # ---------------- (1) core
    tip, rows = read_curve(curve_path)
    tlo, thi = tip
    ncore = 0
    xS_lo = None
    for v in rows:
        t0, t1 = v[0], v[1]
        xlo, xhi = tlo - t1, thi - t0          # x = tau_tip - tau
        if xhi < 0:
            continue
        p = ahull(v[2], v[3]); q = ahull(v[4], v[5]); k = ahull(v[6], v[7]); m = ahull(v[8], v[9])
        stack = [(max(xlo, -1e-3), xhi, 0)]
        while stack:
            a, b_, dep = stack.pop()
            x = ahull(arb(a), arb(b_))
            ok, info = nsd_check(*C.P(x, p, q, k, m))
            ncore += 1
            if not ok:
                if dep < maxdepth:
                    mid = 0.5 * (a + b_); stack += [(a, mid, dep + 1), (mid, b_, dep + 1)]
                else:
                    fails.append(("core", (a, b_), info))
    if verbose:
        print(f"core: {ncore} boxes, {len([f for f in fails if f[0]=='core'])} failures", flush=True)
    # ---------------- (2) tail sigma in [8, sig_v]: x = sigma - sig_tip
    T1 = radii(arb(8), arb(1.62e-4 * (1 + 1e-12)))
    R = T1["R"]
    S = arb(8)
    A_E = arb(3) / 4 * arb(4) ** (-THIRD)
    sig_tip = ahull(arb(8) - arb(thi), arb(8) - arb(tlo))
    n = int(math.ceil((sig_v - 8.0) / dsig))
    ntail = 0
    tail_stack = [(8.0 + i * dsig, 8.0 + (i + 1) * dsig, 0) for i in range(n)]
    while tail_stack:
        a_, b_, dep = tail_stack.pop()
        s0 = arb(repr(a_)); s1 = arb(repr(b_))
        sig = ahull(s0, s1)
        mu0 = (s0 / 4) ** THIRD
        K1 = R * (-A_E * (s0 ** (4 * THIRD) - S ** (4 * THIRD))).exp() * (S / s0) ** THIRD
        kb, mb, pb = K1, mu0 * K1, 4 * mu0 * mu0 * K1      # decreasing in sigma (see lintail2.profile_da)
        J = R * R * (mu0 / 2 + 1 / (12 * s0)) * (-2 * A_E * (s0 ** (4 * THIRD) - S ** (4 * THIRD))).exp() * (S / s0) ** (2 * THIRD)
        k = arb(0, 1) * kb; m = arb(0, 1) * mb; p = arb(0, 1) * pb
        q = sig + arb(0, 1) * (4 * kb * mb) + ahull(arb(0), 4 * J)
        x = sig - sig_tip
        ok, info = nsd_check(*C.P(x, p, q, k, m))
        ntail += 1
        if not ok:
            if dep < maxdepth:
                mm = 0.5 * (a_ + b_); tail_stack += [(a_, mm, dep + 1), (mm, b_, dep + 1)]
            else:
                fails.append(("tail", (a_, b_), info))
    if verbose:
        print(f"tail: {ntail} boxes, {len([f for f in fails if f[0]=='tail'])} failures", flush=True)
    # ---------------- (3) asymptotic: x >= X_v = sig_v - sig_tip_hi
    Xv = arb(sig_v) - sig_tip.upper()
    Xv = arb(Xv.lower())
    uv = Xv ** (-THIRD)
    uv = arb(uv.upper())
    # tiny: Gaussian terms and profile corrections for x >= X_v (crude, rigorous upper bounds)
    smax = max(float(c.mid()) for c in C.nodes)
    dist = Xv - arb(smax)
    gb = 2 * (-(dist * dist) / (2 * C.s * C.s)).exp() * (1 + dist / (C.s * C.s) + 1)
    tsum = sum((abs(t) for t in C.t11[1:] + C.t22[2:] + C.tc[2:] + C.td[2:]), arb(0))
    s_v = arb(sig_v)
    K1 = R * (-A_E * (s_v ** (4 * THIRD) - S ** (4 * THIRD))).exp()
    # profile corrections in scaled entries are bounded by (poly in x) * K1 with x <= anything: use K1 * 1e6 slack
    tiny = arb((tsum * gb * Xv ** 4 + K1 * arb(10) ** 6).upper())
    nasy = 0
    for i in range(nu):
        u = ahull(uv * i / nu, uv * (i + 1) / nu)
        ok, info = nsd_check(*C.Phat_u(u, sig_tip, tiny))
        nasy += 1
        if not ok:
            fails.append(("asym", (float((uv * i / nu).mid()), float((uv * (i + 1) / nu).mid())), info))
    if verbose:
        print(f"asymptotic: {nasy} boxes in u in [0, {float(uv.mid()):.4f}], {len([f for f in fails if f[0]=='asym'])} failures; tiny = {float(tiny.mid()):.2e}", flush=True)
    return fails, C


if __name__ == "__main__":
    fails, C = verify(sys.argv[1])
    print("Lam =", C.d["Lam"], " th =", C.d["th"], " beta =", C.d["beta"])
    print("FAILURES:", len(fails))
    for f in fails[:15]:
        print(f)
