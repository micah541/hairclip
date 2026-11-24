"""Rigorous conical-end analysis for the Jacobi system of A = -H_sym - omega (lam = -omega), sigma >= 32.

Jacobi system (forward in x = sigma - stip):  y' = M y,  y = (u, u', u'', y4),
   M = [[0,1,0,0],[0,0,1,0],[0,-A,P1,-1],[-(V-lam),0,0,P1]],  A = 2(P1^2-P2)+k^2,  V = V0 + dV.
Scaled frozen basis T = D U (D = diag(1,P1,P1^2,P1^3)), U = [Re v(z_u), Im v(z_u), Re v(z_s), Im v(z_s)],
v(z) = (1, z, z^2, -z(z^2 - z + a)),  z(z-1)(z^2-z+a) = c.  With y = T Z:
   Z' = (P1 blockdiag(R(z_u), R(z_s)) + E) Z,  E = -U^{-1} diag(0,r,2r,3r) U - U^{-1} U' + U^{-1} D^{-1} M1 D U.
Lyapunov-Perron (as in step 2): the L^2(rho)-admissible plane is {Z_u = Gamma Z_s} with sup ||Gamma||_2 <= eta,
provided (a + eta^2 b)/g <= eta and 2 eta b/g < 1 on every piece, a >= ||E_us||, b >= ||E_su||,
g <= P1 (Re z_u - Re z_s) - ||E_uu|| - ||E_ss||.  Pieces are taken in u = x^{-1/3} in [0, u0] (u = 0 is x = oo).
Also certified: the plane stays transverse to the Dirichlet plane {u = u' = 0} for all x >= x0
(no conjugate points beyond sigma = 32).
"""
import sys, pickle
import sympy as sp
from flint import arb, acb, arb_mat, acb_mat, ctx
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap')
import lintail2 as LT

ctx.prec = 200
STIP = LT.ahull(arb('-2.249127'), arb('-2.249019'))
SIG0 = arb(32)

class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e):
        return f"R({e.p},{e.q})"
    def _print_Integer(self, e):
        return f"arb({e.p})"
_raw = pickle.load(open('/home/claude/s1spec/step3/hs_uforms.pkl', 'rb'))
_src = {}
for k, v in _raw.items():
    e = eval(v, {**sp.__dict__})
    if k == 'P1u':
        e = sp.Rational(1, 3)/sp.Symbol('w6', positive=True)**sp.Rational(1, 3)   # P1 * x^{-1/3} (hand-simplified)
    import re as _re
    _src[k] = _re.sub(r"\*\*arb\((-?\d+)\)", r"**(\1)", _Pr().doprint(e))

def R(p, q):
    return arb(p)/arb(q)

import os as _os
_FM = _os.environ.get('FORM_MULT', '')
FMULT = (arb(int(_FM.split('/')[0])) / int(_FM.split('/')[1])) if _FM else arb(1)
def ev(name, u, stip, lam):
    w6 = 1 + u**6
    v = eval(_src[name], {'arb': arb, 'R': R, 'u': u, 'w6': w6, 'stip': stip, 'lam': lam})
    return v * FMULT if name in ('a', 'c', 'dadx', 'dcdx') else v

def profile_bounds(s1):
    """|k|,|m|,|p|,|q - sigma| <= (kb, mb, pb, qs) for sigma >= s1 (step-1 tail theorem; decreasing)."""
    s1 = arb(s1)
    m1 = LT.mu(s1)
    K1 = LT.RBAR * (-LT.Efun(s1)).exp() * (LT.S / s1) ** LT.THIRD
    kb, mb, pb = K1, m1 * K1, 4 * m1 * m1 * K1
    J = LT.RBAR * LT.RBAR * (m1 / 2 + 1 / (12 * s1)) * (-2 * LT.Efun(s1)).exp() * (LT.S / s1) ** (2 * LT.THIRD)
    qs = 4 * kb * mb + 4 * J
    return kb, mb, pb, qs

def ball(r):
    return arb(0, arb(r).upper())

def frozen(a, c):
    d = -a * a - 4 * c
    assert d > 0, "frozen discriminant"
    sd = d.sqrt()
    w = acb(-a / 2, sd / 2)
    st = (1 + 4 * w).sqrt()
    assert st.real > 0
    zu, zs = (1 + st) / 2, (1 - st) / 2
    return zu, zs

def vec(z, a):
    return [acb(1), z, z * z, -z * (z * z - z + a)]

def dvec(z, a, dz, da):
    return [acb(0), dz, 2 * z * dz, -(3 * z * z - 2 * z + a) * dz - z * da]

def Umat(zu, zs, a):
    vu, vs = vec(zu, a), vec(zs, a)
    return arb_mat([[vu[i].real, vu[i].imag, vs[i].real, vs[i].imag] for i in range(4)])

def fro(M, rows, cols):
    s = arb(0)
    for i in rows:
        for j in cols:
            t = arb(abs(M[i, j]).upper())
            s += t * t
    return arb(s.sqrt().upper())

def minor(M, r, c):
    return M[r[0], c[0]] * M[r[1], c[1]] - M[r[0], c[1]] * M[r[1], c[0]]

def piece(u_lo, u_hi, lam, stip=STIP):
    """(a_LP, b_LP, g, transversality data) on x in [u_hi^{-3}, u_lo^{-3}] (u_lo = 0 allowed)."""
    u = LT.ahull(arb(u_lo), arb(u_hi))
    a, c = ev('a', u, stip, lam), ev('c', u, stip, lam)
    da, dc = ev('dadx', u, stip, lam), ev('dcdx', u, stip, lam)
    r = ev('r', u, stip, lam)
    iP1, iP13 = ev('invP1', u, stip, lam), ev('invP13', u, stip, lam)
    zu, zs = frozen(a, c)
    U = Umat(zu, zs, a)
    Ui = U.inv()
    def dz(z):
        return (dc - z * (z - 1) * da) / ((2 * z - 1) * (2 * z * z - 2 * z + a))
    vu1, vs1 = dvec(zu, a, dz(zu), da), dvec(zs, a, dz(zs), da)
    dU = arb_mat([[vu1[i].real, vu1[i].imag, vs1[i].real, vs1[i].imag] for i in range(4)])
    DD = arb_mat([[0, 0, 0, 0], [0, r, 0, 0], [0, 0, 2 * r, 0], [0, 0, 0, 3 * r]])
    # profile perturbation (scaled): entries (2,1) = -k^2/P1, (3,0) = -dV/P1^3
    x_min_sig = arb(u_hi) ** -3 + stip.lower()          # smallest sigma on the piece
    kb, mb, pb, qs = profile_bounds(arb(x_min_sig.lower()))
    # |dV|/P1^3 <= k^2/(2P1) + |r| k^2/(2P1^2) + |k m|/(2P1^2) + (|k p| + 3 m^2)/(2 P1^3) + |q - sigma|/(8 P1^2)
    dVs = (kb * kb / 2 * iP1 + kb * kb / 2 * abs(r) * iP1 * iP1 + kb * mb / 2 * iP1 * iP1
           + (kb * pb / 2 + 3 * mb * mb / 2) * iP13 + qs / 8 * iP1 * iP1)
    k2s = kb * kb * iP1
    dVs, k2s = dVs * FMULT, k2s * FMULT
    P1lo = 1 / (3 * arb(u_hi) * (1 + arb(u_hi) ** 6) ** (arb(1) / 3))
    M1 = arb_mat([[0, 0, 0, 0], [0, 0, 0, 0], [0, ball(abs(k2s)), 0, 0], [ball(abs(dVs)), 0, 0, 0]])
    E = -(Ui * DD * U) - Ui * dU + Ui * M1 * U
    nuu, nss = fro(E, (0, 1), (0, 1)), fro(E, (2, 3), (2, 3))
    nus, nsu = fro(E, (0, 1), (2, 3)), fro(E, (2, 3), (0, 1))
    gz = (zu.real - zs.real)
    g = arb((P1lo * gz).lower()) - nuu - nss
    # transversality to the Dirichlet plane: minors of U rows (0,1)
    rows = (0, 1)
    main = minor(U, rows, (2, 3))
    others = [minor(U, rows, cc) for cc in [(0, 2), (0, 3), (1, 2), (1, 3)]]
    m01 = minor(U, rows, (0, 1))
    return dict(a=nus, b=nsu, g=g, main=main, others=others, m01=m01, zu=zu, zs=zs, U=U)

def u_grid(u0, N):
    # geometric-ish grid in u from u0 down to small, then final piece [0, u_min]
    import math
    us = [float(u0) * (0.97 ** i) for i in range(N)]
    return us

def tail_theorem(lam, N=140):
    lam = arb(lam)
    x0 = SIG0 - STIP                         # x0 interval
    u0 = arb((x0.lower()) ** (-arb(1) / 3)).upper()   # u0 >= x0^{-1/3} for every stip
    u0 = arb(u0)
    grid = u_grid(float(u0.upper()) * (1 + 1e-12), N)
    pieces = []
    for i in range(len(grid) - 1):
        pieces.append(piece(arb(repr(grid[i + 1])), arb(repr(grid[i])), lam))
    pieces.append(piece(arb(0), arb(repr(grid[-1])), lam))
    eta_lo = arb(0); eta_hi = arb(1e30)
    for P in pieces:
        a, b, g = P['a'], P['b'], P['g']
        assert g > 0, "gap"
        disc = g * g - 4 * a * b
        assert disc > 0, "no invariant ball"
        r_lo = 2 * a / (g + disc.sqrt())
        r_hi = (g + disc.sqrt()) / (2 * b) if b > 0 else arb(1e30)
        eta_lo = arb(max(eta_lo, arb(r_lo.upper())).upper())
        eta_hi = arb(min(eta_hi, arb(r_hi.lower())).lower())
    eta = arb(eta_lo.upper() * (1 + 2.0 ** -30))
    assert eta < eta_hi
    Lip = arb(0)
    trans_margin = arb(1e30)
    for P in pieces:
        a, b, g = P['a'], P['b'], P['g']
        assert (a + eta * eta * b) / g <= eta
        Lip = arb(max(Lip, arb((2 * eta * b / g).upper())).upper())
        lhs = abs(P['main'])
        rhs = eta * sum(abs(o) for o in P['others']) + eta * eta * abs(P['m01'])
        assert lhs > rhs, "transversality"
        trans_margin = arb(min(trans_margin, arb((lhs - rhs).lower())).lower())
    assert Lip < 1
    return dict(eta=eta, Lip=Lip, gmin=min(float(P['g'].lower()) for P in pieces), trans=trans_margin,
                npieces=len(pieces), u0=u0)

def plane_at_x0(lam, eta):
    """frozen basis T(x0) = D U at sigma = 32 (all stip), for building CAPD inputs."""
    lam = arb(lam)
    x0 = SIG0 - STIP
    u = x0 ** (-arb(1) / 3)
    a, c = ev('a', u, STIP, lam), ev('c', u, STIP, lam)
    zu, zs = frozen(a, c)
    U = Umat(zu, zs, a)
    P1 = 1 / (3 * u * (1 + u ** 6) ** (arb(1) / 3))
    D = [arb(1), P1, P1 * P1, P1 ** 3]
    T = arb_mat([[D[i] * U[i, j] for j in range(4)] for i in range(4)])
    return T, zu, zs, P1

if __name__ == "__main__":
    import time
    lam = sys.argv[1] if len(sys.argv) > 1 else '-0.05'
    t0 = time.time()
    res = tail_theorem(lam)
    print(f"lam={lam}: eta={float(res['eta'].upper()):.4e}  Lip={float(res['Lip'].upper()):.3e}  gmin={res['gmin']:.4f}"
          f"  transversality margin={float(res['trans'].lower()):.4e}  pieces={res['npieces']}  ({time.time()-t0:.1f}s)")
    T, zu, zs, P1 = plane_at_x0(lam, res['eta'])
    print("z_u =", zu, "\nz_s =", zs, "\nP1(x0) =", P1)
