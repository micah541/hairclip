"""float (mpmath) check of T4: mu1, mu2 (and the exact relation mu2 = mu1 + 1/2 + O(profile)) at points, with P_j = psi^(j) by numerical differentiation of psi'."""
import sys; sys.argv = ['x', '0.5']
import mpmath as mp
from int_garding import *
mp.mp.dps = 40
def dpsi(x):
    t = (x - mp.mpf('30.25')) / 8
    s = 0 if t <= 0 else (1 if t >= 1 else t**4 * (35 - 84*t + 70*t**2 - 20*t**3))
    return (1 - s) * 6 * x / (1 + x * x) + s * x * (1 + x * x) ** (-mp.mpf(1) / 3) / 3
f2 = sp.lambdify((P1, P2, P3, P4, k, m, p, q), mu2, 'mpmath'); f1 = sp.lambdify((P1, P2, P3, P4, k, m, p, q), mu1, 'mpmath')
for x in (12.0, 16.0, 25.0, 33.0, 36.0, 45.0, 59.0):
    xm = mp.mpf(x); P = [dpsi(xm)] + [mp.diff(dpsi, xm, j) for j in (1, 2, 3)]
    sg = xm + mp.mpf(float(ST.mid()))
    v2, v1 = f2(*P, 0, 0, 0, sg), f1(*P, 0, 0, 0, sg)
    for sa_, sb_ in make_segments([('x', 12.0), ('x', 60.0)], [0.01]):
        b = make_box(sa_, sb_)
        if b.xa <= x <= b.xb: break
    Pb = Pj_box(b.X)
    r2, r1 = fmu2(*Pb, b.k, b.m, b.p, b.q), fmu1(*Pb, b.k, b.m, b.p, b.q)
    print(f"x={x}: float mu1={mp.nstr(v1,8)} mu2={mp.nstr(v2,8)} (mu2-mu1={mp.nstr(v2-v1,6)});  rigorous box mu1 in [{float(r1.lower().mid()):.6f},{float(r1.upper().mid()):.6f}] mu2 in [{float(r2.lower().mid()):.6f},{float(r2.upper().mid()):.6f}]")
