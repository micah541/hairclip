"""
Rigorous admissible plane at sigma = S = 8 for the eigenvalue ODE of the linearisation at S1
(conical end), uniformly for lam in a complex box.  Ball arithmetic: python-flint / Arb.

Eigen-ODE (tail normalisation, sigma >= 8):  u'''' = a2 u'' + a1 u' + a0 u,
    a0 = c - 3 m^2 - k p,  c = 1/4 - lam,   a1 = -sigma/4 - 5 k m - (q - sigma)/4,   a2 = -k^2,
with the profile bounds of the step-1 tail theorem (|k|, |k'|/mu, |k''|/mu^2 <= Rbar/w, w = e^E (sigma/S)^{1/3},
J <= Rbar^2 (mu/2 + 1/(12 sigma)) w^-2, q - sigma = 4 k k' + 4 J, p = 4 k'').

Frozen block basis.  mu = (sigma/4)^{1/3}, eps = c/mu^4, zeta^4 + zeta = eps.  G-roots z2, z3 (-> e^{+-i pi/3}),
S-pair (slow, decaying) through s = z0+z1 = -(z2+z3), p = z0 z1 = (z2+z3)^2 - z2 z3.
    U = [e+, e-, v(z2), v(z3)],  e+ = (1, s/2, (s^2-2p)/2, (s^3-3sp)/2),  e- = (0, 1, s, s^2-p),  v(z) = (1,z,z^2,z^3),
    Y = (u,u',u'',u''') = T Z,  T = D_mu U Wsc,  D_mu = diag(1,mu,mu^2,mu^3),  Wsc = diag(1,w,1,1).
    Z' = B Z,  B = Wsc^{-1}[ mu blockdiag(Bs, z2, z3) - Mt/(3 sigma) + P ] Wsc,
    Bs = [[s/2, 1], [(s^2-4p)/4, s/2]],  Mt = U^{-1} diag(0,1,2,3) U - 4 eps U^{-1} dU/deps,
    P = U^{-1} e4 (da0/mu^3, da1/mu^2, da2/mu, 0) U   (profile perturbation).
(The block form and Mt were checked against direct computation in check_alg.py; the proof only uses
 the exact identities  U^{-1} A_zeta U = blockdiag(Bs,z2,z3)  and  T^{-1}T' = Mt/(3 sigma).)

Lyapunov-Perron: the admissible plane is {Z_G = Gamma Z_S} with Gamma the unique bounded solution of
    Gamma' = B_GS + B_GG Gamma - Gamma B_SS - Gamma B_SG Gamma,
and sup ||Gamma|| <= eta whenever  (a + eta^2 b)/g <= eta  on every sigma-piece, where a >= ||B_GS||,
b >= ||B_SG||, g <= lambda_min(Herm B_GG) - lambda_max(Herm B_SS); uniqueness/contraction: 2 eta b/g < 1.
"""
import sys, math
import numpy as np
from flint import arb, acb, acb_mat, ctx

sys.path.insert(0, '/home/claude/s1_cap')
from tail import radii   # step-1 tail theorem (verified there)

ctx.prec = 128
S = arb(8)            # base point of the step-1 profile tail theorem
S0 = arb(32)          # start of the Lyapunov-Perron interval for the eigen-ODE
THIRD = arb(1) / 3
A_E = arb(3) / 4 * arb(4) ** (-THIRD)
WSC = arb('0.5')
_T1 = radii(arb(8), arb(1.62e-4 * (1 + 1e-12)))
assert _T1["verified"]
RBAR = _T1["R"]
G0 = np.array([np.exp(1j * np.pi / 3), np.exp(-1j * np.pi / 3)])


def mu(s):
    return (arb(s) / 4) ** THIRD


def Efun(s):
    s = arb(s)
    return A_E * (s ** (arb(4) / 3) - S ** (arb(4) / 3))


def hull(a, b):
    """arb ball containing a and b (arb)"""
    lo = min(a.lower(), b.lower()); hi = max(a.upper(), b.upper())
    return arb((lo + hi) / 2, (hi - lo) / 2) if True else None


def ahull(lo, hi):
    lo, hi = arb(lo), arb(hi)
    l = lo.lower(); h = hi.upper()
    return arb((l + h) / 2) + arb(0, 1) * ((h - l) / 2)


def sq(x):
    """rigorous square of a real ball, with nonnegative lower bound"""
    a = abs(x)
    lo, hi = a.lower(), a.upper()
    if lo < 0:
        lo = arb(0)
    l2, h2 = lo * lo, hi * hi
    l = l2.lower()
    if l < 0:
        l = arb(0)
    return ahull(l, h2.upper())


def sqrtpos(x):
    lo, hi = x.lower(), x.upper()
    if not (lo > 0):
        lo = arb(0)
    return ahull(lo.sqrt().lower() if lo > 0 else arb(0), hi.sqrt().upper())


def cbox(re_lo, re_hi, im_lo, im_hi):
    return acb(ahull(re_lo, re_hi), ahull(im_lo, im_hi))


def numroots(eps, n=40):
    z = G0.copy()
    for t in np.linspace(0, 1, n + 1)[1:]:
        for _ in range(8):
            z = z - (z ** 4 + z - t * eps) / (4 * z ** 3 + 1)
    for _ in range(30):
        z = z - (z ** 4 + z - eps) / (4 * z ** 3 + 1)
    return z


def krawczyk(zt, eps_ball):
    """enclosure of the unique root of z^4 + z - eps in a disc around zt, for every eps in eps_ball"""
    zt = acb(complex(zt).real, complex(zt).imag)
    fp = 4 * zt ** 3 + 1
    Y = 1 / acb(fp.real.mid(), fp.imag.mid())
    f0 = zt ** 4 + zt - eps_ball
    r = 4 * float(abs(f0).upper()) / float(abs(fp).lower()) + 1e-30
    for _ in range(40):
        Z = acb(arb(zt.real, r), arb(zt.imag, r))
        K = zt - Y * f0 + (1 - Y * (4 * Z ** 3 + 1)) * (Z - zt)
        if Z.contains_interior(K):
            return K
        r *= 1.6
    raise RuntimeError("Krawczyk failed")


def frozen(eps_ball, eps_num):
    """rigorous frozen-basis objects for all eps in eps_ball; eps_num = numerical centre for root tracking"""
    zn = numroots(eps_num)
    z2 = krawczyk(zn[0], eps_ball); z3 = krawczyk(zn[1], eps_ball)
    sp_ = z2 + z3; pp = z2 * z3
    s = -sp_; p = -eps_ball / pp          # z0 z1 z2 z3 = -eps  (no cancellation)
    dz2 = 1 / (4 * z2 ** 3 + 1); dz3 = 1 / (4 * z3 ** 3 + 1)
    dsp = dz2 + dz3; dpp = dz2 * z3 + z2 * dz3
    ds = -dsp; dp = -1 / pp + eps_ball * dpp / (pp * pp)
    ep = [acb(1), s / 2, (s * s - 2 * p) / 2, (s ** 3 - 3 * s * p) / 2]
    em = [acb(0), acb(1), s, s * s - p]
    dep = [acb(0), ds / 2, (2 * s * ds - 2 * dp) / 2, (3 * s * s * ds - 3 * ds * p - 3 * s * dp) / 2]
    dem = [acb(0), acb(0), ds, 2 * s * ds - dp]
    v = lambda z: [acb(1), z, z * z, z ** 3]
    dv = lambda z, d: [acb(0), d, 2 * z * d, 3 * z * z * d]
    cols = [ep, em, v(z2), v(z3)]
    dcols = [dep, dem, dv(z2, dz2), dv(z3, dz3)]
    U = acb_mat([[cols[j][i] for j in range(4)] for i in range(4)])
    dU = acb_mat([[dcols[j][i] for j in range(4)] for i in range(4)])
    Ui = U.inv()
    Dg = acb_mat([[acb(i if i == j else 0) for j in range(4)] for i in range(4)])
    Mt = Ui * Dg * U - (4 * eps_ball) * (Ui * dU)
    Bs = [[s / 2, acb(1)], [(s * s - 4 * p) / 4, s / 2]]
    return dict(z2=z2, z3=z3, s=s, p=p, U=U, Ui=Ui, Mt=Mt, Bs=Bs, ds=ds, dp=dp, dz2=dz2, dz3=dz3, dU=dU)


def profile_da(s1):
    """bounds (arb, >=0) for |da0|, |da1|, |da2| valid on [s1, oo) (all are decreasing in sigma)"""
    s1 = arb(s1)
    m1 = mu(s1)
    K1 = RBAR * (-Efun(s1)).exp() * (S / s1) ** THIRD
    kb, mb, pb = K1, m1 * K1, 4 * m1 * m1 * K1
    J = RBAR * RBAR * (m1 / 2 + 1 / (12 * s1)) * (-2 * Efun(s1)).exp() * (S / s1) ** (2 * THIRD)
    qs = 4 * kb * mb + 4 * J
    return 3 * mb * mb + kb * pb, 5 * kb * mb + qs / 4, kb * kb


def herm2_max(a, b, c, d):
    """upper bound of lambda_max of the Hermitian part of [[a,b],[c,d]] (acb entries)"""
    off = (b + c.conjugate()) / 2
    h = (a.real + d.real) / 2 + sqrtpos(sq((a.real - d.real) / 2) + sq(abs(off)))
    return h


def herm2_min(a, b, c, d):
    off = (b + c.conjugate()) / 2
    return (a.real + d.real) / 2 - sqrtpos(sq((a.real - d.real) / 2) + sq(abs(off)))


def fro(M, rows, cols):
    """upper bound of the spectral norm of the 2x2 block M[rows, cols]:
       ||M||^2 = (F^2 + sqrt(F^4 - 4|det|^2))/2, increasing in F^2 and decreasing in |det|"""
    r0, r1 = rows; c0, c1 = cols
    F2 = sq(abs(M[r0, c0])) + sq(abs(M[r0, c1])) + sq(abs(M[r1, c0])) + sq(abs(M[r1, c1]))
    F2u = arb(F2.upper())
    D = abs(M[r0, c0] * M[r1, c1] - M[r0, c1] * M[r1, c0])
    Dl = D.lower()
    if not (Dl > 0):
        Dl = arb(0)
    inner = sq(F2u) - 4 * Dl * Dl
    n2 = (F2u + sqrtpos(inner)) / 2
    return arb(sqrtpos(n2).upper())


def scaled(M):
    """Wsc^{-1} M Wsc,  Wsc = diag(1, w, 1, 1)"""
    w = [arb(1), WSC, arb(1), arb(1)]
    return acb_mat([[M[i, j] * w[j] / w[i] for j in range(4)] for i in range(4)])


def pert_matrix(F, s1):
    """P = U^{-1} e4 row U as a ball matrix, row = (da0/mu^3, da1/mu^2, da2/mu, 0) with da as [-b,b] balls"""
    d0, d1, d2 = profile_da(s1)
    m1 = mu(s1)
    row = [acb(arb(0, 1) * d0 / m1 ** 3), acb(arb(0, 1) * d1 / m1 ** 2), acb(arb(0, 1) * d2 / m1), acb(0)]
    U, Ui = F["U"], F["Ui"]
    rU = [sum((row[l] * U[l, j] for l in range(4)), acb(0)) for j in range(4)]
    return acb_mat([[Ui[i, 3] * rU[j] for j in range(4)] for i in range(4)])


def piece(c, c_num, s1, s2):
    """bounds (a, b, g) on the piece sigma in [s1, s2] (s2 = None: final piece [s1, oo))"""
    s1 = arb(s1)
    t1 = (4 / s1) ** (4 * THIRD)
    if s2 is None:
        t = ahull(arb(0), t1)
    else:
        s2 = arb(s2); t = ahull((4 / s2) ** (4 * THIRD), t1)
    eps = c * t
    tn = float(t.mid())
    F = frozen(eps, c_num * tn)
    Bs, z2, z3, Mt = F["Bs"], F["z2"], F["z3"], F["Mt"]
    P = pert_matrix(F, s1)
    w = WSC
    Bss = [[Bs[0][0], Bs[0][1] * w], [Bs[1][0] / w, Bs[1][1]]]       # scaled S-block (zeta units)
    Mts, Ps = scaled(Mt), scaled(P)
    mx = max(float(abs(Mts[i, j]).upper()) for i in range(4) for j in range(4))
    rq = max(float(Mts[i, j].real.rad()) + float(Mts[i, j].imag.rad()) for i in range(4) for j in range(4))
    rq = max(rq / mx, (float(z2.real.rad()) + float(z2.imag.rad())) / 0.5)
    piece.quality = rq
    nMSS, nMGG = fro(Mts, (0, 1), (0, 1)), fro(Mts, (2, 3), (2, 3))
    nMGS, nMSG = fro(Mts, (2, 3), (0, 1)), fro(Mts, (0, 1), (2, 3))
    nPSS, nPGG = fro(Ps, (0, 1), (0, 1)), fro(Ps, (2, 3), (2, 3))
    nPGS, nPSG = fro(Ps, (2, 3), (0, 1)), fro(Ps, (0, 1), (2, 3))
    lS = herm2_max(Bss[0][0], Bss[0][1], Bss[1][0], Bss[1][1])      # zeta units
    mG = herm2_min(z2, acb(0), acb(0), z3)
    if s2 is None:
        # final piece: g(sigma) >= mu (mG - lS) - (nMSS+nMGG)/(3 sigma) - (nPSS+nPGG) is increasing in sigma
        # provided mG - lS > 0 on the eps-ball (checked); a, b decreasing.  Evaluate at sigma = s1.
        gz = mG - lS
        assert gz > 0, "final piece: zeta-gap not positive"
        m1 = mu(s1)
        g = m1 * gz - (nMSS + nMGG) / (3 * s1) - (nPSS + nPGG)
        a = nMGS / (3 * s1) + nPGS
        b = nMSG / (3 * s1) + nPSG
        return a, b, g
    mu_iv = ahull(mu(s1), mu(s2))
    inv3s = ahull(1 / (3 * s2), 1 / (3 * s1))
    g = mu_iv * mG - mu_iv * lS
    # mu*mG - mu*lS as interval products (mG, lS may have any sign): take worst case
    g = (mu_iv * mG).lower() - (mu_iv * lS).upper()
    g = arb(g) - (nMSS + nMGG) * inv3s.upper() - (nPSS + nPGG)
    a = nMGS * inv3s.upper() + nPGS
    b = nMSG * inv3s.upper() + nPSG
    return a, b, g


def sigma_grid(N=48, smax=400.0, s0=None):
    s0 = float(S0.mid()) if s0 is None else s0
    xs = np.exp(np.linspace(np.log(s0), np.log(smax), N + 1))
    return [arb(repr(float(x))) for x in xs[:-1]] + [arb(repr(float(smax)))]


def eta_bound(c, c_num, N=120, smax=600.0):
    """constant-radius Lyapunov-Perron on [S0, oo): eta with (a + eta^2 b)/g <= eta on all pieces, 2 eta b/g < 1"""
    grid = sigma_grid(N=N, smax=smax)
    def refined(s_lo, s_hi, depth=0):
        try:
            a, b, g = piece(c, c_num, s_lo, s_hi)
            good = (g > 0) and float((4 * a * b / (g * g)).upper()) < 0.5
        except (ZeroDivisionError, RuntimeError):
            if depth >= 8:
                raise
            good = False
        if good or depth >= 8:
            return [(a, b, g)]
        sm = arb(repr(float((s_lo * s_hi).sqrt().mid())))
        return refined(s_lo, sm, depth + 1) + refined(sm, s_hi, depth + 1)
    rows = []
    for i in range(len(grid) - 1):
        rows += refined(grid[i], grid[i + 1])
    rows.append(piece(c, c_num, grid[-1], None))
    eta_lo = arb(0); eta_hi = arb(1e30)
    for a, b, g in rows:
        if not (g > 0):
            raise RuntimeError("gap not positive")
        disc = sq(g) - 4 * a * b
        if not (disc > 0):
            raise RuntimeError("no invariant ball")
        r_lo = 2 * a / (g + sqrtpos(disc))
        r_hi = (g + sqrtpos(disc)) / (2 * b) if b > 0 else arb(1e30)
        eta_lo = arb(max(eta_lo, arb(r_lo.upper())).upper())
        eta_hi = arb(min(eta_hi, arb(r_hi.lower())).lower())
    eta = arb(eta_lo.upper() * (1 + 2.0 ** -30))
    if not (eta < eta_hi):
        raise RuntimeError("eta interval empty")
    Lip = arb(0)
    for a, b, g in rows:
        if not ((a + eta * eta * b) / g <= eta):
            raise RuntimeError("invariance check failed")
        Lip = arb(max(Lip, arb((2 * eta * b / g).upper())).upper())
    if not (Lip < 1):
        raise RuntimeError("no contraction")
    return eta, Lip, min(float(g.lower()) for a, b, g in rows)


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def plucker(ya, yb):
    return [ya[i] * yb[j] - ya[j] * yb[i] for (i, j) in PAIRS]


def plane_at_S(c, c_num):
    """T(S0) = D_mu U Wsc for all c in the ball; returns columns as lists of acb"""
    m = mu(S0)
    eps = c / m ** 4
    F = frozen(eps, c_num / float(m.mid()) ** 4)
    U = F["U"]
    w = [arb(1), WSC, arb(1), arb(1)]
    Dm = [arb(1), m, m * m, m ** 3]
    T = [[U[i, j] * Dm[i] * w[j] for j in range(4)] for i in range(4)]
    return T


def admissible_W(lam_box):
    """lam_box = (re_lo, re_hi, im_lo, im_hi) (python floats).  Returns a dict with
       W0 (6 acb, valid over the box), cols (list of 8 complex 6-vectors at the centre, as python complex),
       dcols (6 acb per column: enclosure of col(lam)-col_centre over the box), eta, W2 bound (arb)."""
    rl, rh, il, ih = lam_box
    lam = cbox(rl, rh, il, ih)
    c = 0.25 - lam
    lc = complex((rl + rh) / 2, (il + ih) / 2)
    c_num = 0.25 - lc
    eta, Lip, gmin = eta_bound(c, c_num)
    T = plane_at_S(c, c_num)
    lamc = acb(lc.real, lc.imag)
    Tc = plane_at_S(0.25 - lamc, c_num)
    col = lambda TT, j: [TT[i][j] for i in range(4)]
    W0 = plucker(col(T, 0), col(T, 1))
    # W1 = d/dGamma: y_a = T_S0 + T_G Gamma[:,0], y_b = T_S1 + T_G Gamma[:,1];  Gamma[l][j], l in G (2,3), j in S (0,1)
    def w1cols(TT):
        out = []
        for j in (0, 1):            # which S column is perturbed
            for l in (2, 3):        # by which G column
                if j == 0:
                    v = plucker(col(TT, l), col(TT, 1))
                else:
                    v = plucker(col(TT, 0), col(TT, l))
                out.append(v)
        return out
    C_box = w1cols(T); C_c = w1cols(Tc)
    cols, dcols = [], []
    for v_box, v_c in zip(C_box, C_c):
        vc = [complex(float(x.real.mid()), float(x.imag.mid())) for x in v_c]
        cols.append(vc)
        dcols.append([vb - acb(z.real, z.imag) for vb, z in zip(v_box, vc)])
    W2dir = plucker(col(T, 2), col(T, 3))           # times det(Gamma), |det Gamma| <= eta^2
    return dict(W0=W0, cols=cols, dcols=dcols, eta=eta, Lip=Lip, gmin=gmin, W2dir=W2dir)


if False:
    import time
    for box in [(1.0, 1.0, 0.0, 0.0), (-0.125, -0.125, 0.0, 0.0), (-0.125, -0.12, 3.9, 3.95), (1.75, 1.8, 0.3, 0.35),
                (0.3, 0.31, 0.2, 0.21)]:
        t0 = time.time()
        R = admissible_W(box)
        print(f"box {box}: pieces={eta_bound.npieces} eta = {float(R['eta'].upper()):.4f}  Lip = {float(R['Lip'].upper()):.3f}  gmin = {R['gmin']:.3f}"
              f"  |W0_13| rad = {float(R['W0'][4].rad()):.2e}  ({time.time()-t0:.1f}s)", flush=True)


def profile_box_at(s0, eps_lo, eps_hi):
    """enclosure of the profile state (p, q, k, m) at sigma = s0 >= S from the step-1 tail theorem, for every
       eps in [eps_lo, eps_hi]:  Z0(s0) = (eps + [-r0,r0]) e^{-E}(S/s0)^{1/3},  |Z1(s0)| <= r1 e^{-E}(S/s0)^{5/3},
       (k, k', k'') = (Z0 + 2Re Z1, -mu (Z0 + 2Re(w1 Z1)), mu^2 (Z0 + 2Re(w2 Z1))),  q = s0 + 4 k k' + 4 J."""
    from tail import W1, W2
    s0 = arb(s0)
    m = mu(s0)
    fac = (-Efun(s0)).exp() * (S / s0) ** THIRD
    r0, r1 = _T1["r0"], _T1["r1"]
    e = ahull(arb(eps_lo), arb(eps_hi))
    Z0 = (e + arb(0, 1) * r0) * fac
    rz1 = r1 * fac * (S / s0) ** (4 * THIRD)
    Z1 = acb(arb(0, 1) * rz1, arb(0, 1) * rz1)
    k = Z0 + 2 * Z1.real
    kp = -m * (Z0 + 2 * (W1 * Z1).real)
    kpp = m * m * (Z0 + 2 * (W2 * Z1).real)
    Jb = RBAR * RBAR * (m / 2 + 1 / (12 * s0)) * (-2 * Efun(s0)).exp() * (S / s0) ** (2 * THIRD)
    J = ahull(arb(0), Jb)
    q = s0 + 4 * k * kp + 4 * J
    return [4 * kpp, q, k, kp]


if __name__ == "__main__":
    import time
    for box in [(1.0, 1.0, 0.0, 0.0), (-0.125, -0.125, 0.0, 0.0), (-0.125, -0.125, 1.0, 1.0), (-0.125, -0.12, 3.9, 3.95),
                (1.75, 1.8, 0.3, 0.35), (0.3, 0.31, 0.2, 0.21), (1.8, 1.85, -0.05, 0.05), (-0.2, 0.05, 0.75, 1.0)]:
        t0 = time.time()
        try:
            R = admissible_W(box)
            print(f"box {box}: eta={float(R['eta'].upper()):.5f} Lip={float(R['Lip'].upper()):.4f} gmin={R['gmin']:.3f} ({time.time()-t0:.1f}s)", flush=True)
        except Exception as ex:
            print(f"box {box}: FAILED {ex!r} ({time.time()-t0:.1f}s)", flush=True)
    print(profile_box_at(32, 1.5345849609375e-4, 1.5346679687499998e-4))


def plane_and_deriv(c, c_num):
    """T(S0) and dT/dlam (dlam = -dc; d eps/d lam = -1/mu^4) for all c in the ball"""
    m = mu(S0)
    m4 = m ** 4
    F = frozen(c / m4, c_num / float(m4.mid()))
    U, dU = F["U"], F["dU"]
    w = [arb(1), WSC, arb(1), arb(1)]
    Dm = [arb(1), m, m * m, m ** 3]
    T = [[U[i, j] * Dm[i] * w[j] for j in range(4)] for i in range(4)]
    dT = [[-dU[i, j] * Dm[i] * w[j] / m4 for j in range(4)] for i in range(4)]
    return T, dT


def segment_W(lam_a, lam_b, eta=None):
    """Data for the CAPD run on the segment lam(t) = lam_a + t (lam_b - lam_a), t in [0,1]:
       W0c = W0(lam_c) (6 acb, lam_c = midpoint), dWt = dW0/dt at lam_c (6 acb),
       rem = enclosure of W0(lam(t)) - W0c - dWt (t - 1/2) over the segment (6 acb),
       cols = Gamma-columns at lam_c (4 complex 6-vectors), dcols = their variation over the segment box,
       eta = LP radius valid on the segment box, W2dir (6 acb over the box)."""
    la, lb = complex(lam_a), complex(lam_b)
    D = lb - la
    lc = (la + lb) / 2
    box = (min(la.real, lb.real), max(la.real, lb.real), min(la.imag, lb.imag), max(la.imag, lb.imag))
    lam_box = cbox(*box)
    c_box = 0.25 - lam_box
    c_num = 0.25 - lc
    if eta is None:
        eta, Lip, gmin = eta_bound(c_box, c_num)
    else:
        Lip, gmin = None, None
    lamc = acb(lc.real, lc.imag)
    Dl = acb(D.real, D.imag)
    T_c, dT_c = plane_and_deriv(0.25 - lamc, c_num)
    T_b, dT_b = plane_and_deriv(c_box, c_num)
    col = lambda TT, j: [TT[i][j] for i in range(4)]
    W0c = plucker(col(T_c, 0), col(T_c, 1))
    def dW(TT, dTT):
        a = plucker(col(dTT, 0), col(TT, 1)); b = plucker(col(TT, 0), col(dTT, 1))
        return [x + y for x, y in zip(a, b)]
    dW_c = dW(T_c, dT_c); dW_b = dW(T_b, dT_b)
    dWt = [x * Dl for x in dW_c]
    half = acb(arb(0, 0.5))
    rem = [(xb - xc) * Dl * half for xb, xc in zip(dW_b, dW_c)]
    def w1cols(TT):
        out = []
        for j in (0, 1):
            for l in (2, 3):
                out.append(plucker(col(TT, l), col(TT, 1)) if j == 0 else plucker(col(TT, 0), col(TT, l)))
        return out
    C_b, C_c = w1cols(T_b), w1cols(T_c)
    cols, dcols = [], []
    for vb, vc in zip(C_b, C_c):
        z = [complex(float(x.real.mid()), float(x.imag.mid())) for x in vc]
        cols.append(z)
        dcols.append([x - acb(y.real, y.imag) for x, y in zip(vb, z)])
    W2dir = plucker(col(T_b, 2), col(T_b, 3))
    W2c = plucker(col(T_c, 2), col(T_c, 3))
    W2dir_c = [complex(float(x.real.mid()), float(x.imag.mid())) for x in W2c]
    return dict(W0c=W0c, dWt=dWt, rem=rem, cols=cols, dcols=dcols, eta=eta, Lip=Lip, gmin=gmin, W2dir=W2dir,
                W2dir_c=W2dir_c, box=box)


def box_W(box):
    """Data for a lambda-box run (lam an interval parameter): W0 over the box as centre + mean-value remainder,
       Gamma-columns at the centre and their variation, eta valid on the box."""
    rl, rh, il, ih = box
    lc = complex((rl + rh) / 2, (il + ih) / 2)
    lam_box = cbox(rl, rh, il, ih)
    c_box = 0.25 - lam_box
    c_num = 0.25 - lc
    eta, Lip, gmin = eta_bound(c_box, c_num)
    lamc = acb(lc.real, lc.imag)
    T_c, dT_c = plane_and_deriv(0.25 - lamc, c_num)
    T_b, dT_b = plane_and_deriv(c_box, c_num)
    col = lambda TT, j: [TT[i][j] for i in range(4)]
    W0c = plucker(col(T_c, 0), col(T_c, 1))
    a = plucker(col(dT_b, 0), col(T_b, 1)); b = plucker(col(T_b, 0), col(dT_b, 1))
    dW_b = [x + y for x, y in zip(a, b)]
    dlam = lam_box - lamc
    rem = [x * dlam for x in dW_b]
    def w1cols(TT):
        out = []
        for j in (0, 1):
            for l in (2, 3):
                out.append(plucker(col(TT, l), col(TT, 1)) if j == 0 else plucker(col(TT, 0), col(TT, l)))
        return out
    C_b, C_c = w1cols(T_b), w1cols(T_c)
    cols, dcols = [], []
    for vb, vc in zip(C_b, C_c):
        z = [complex(float(x.real.mid()), float(x.imag.mid())) for x in vc]
        cols.append(z)
        dcols.append([x - acb(y.real, y.imag) for x, y in zip(vb, z)])
    W2dir = plucker(col(T_b, 2), col(T_b, 3))
    zero = [acb(0)] * 6
    return dict(W0c=W0c, dWt=zero, rem=rem, cols=cols, dcols=dcols, eta=eta, Lip=Lip, gmin=gmin, W2dir=W2dir, box=box)
