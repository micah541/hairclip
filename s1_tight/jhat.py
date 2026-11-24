"""Jhat = int_S^oo K1(sigma)^2 dsigma via the double contour formula
   int_S^oo K1^2 = -(1/(2 pi i))^2 int_C int_C t s exp(t^4 + s^4 + S(t+s)) / (t+s) dt ds   (Re(t+s) < 0 on C x C).
C = upper ray (t0 + r e3, outward) minus lower ray (t0 + r conj(e3), outward)."""
import sys
from flint import arb, acb, ctx
from laplace import K
ctx.prec = 80
E3 = acb(-1, 1) / arb(2).sqrt(); E3c = E3.conjugate()
def jhat(S, R=6):
    S = arb(S); t0 = acb(-(S / 4) ** (arb(1) / 3))
    rays = [(E3, 1), (E3c, -1)]
    tot = acb(0)
    for ea, sa in rays:
        for eb, sb in rays:
            def inner(r, analytic):
                t = t0 + r * ea
                def f(rho, an2):
                    s = t0 + rho * eb
                    return t * s * (t ** 4 + s ** 4 + S * (t + s)).exp() / (t + s) * eb
                return acb.integral(f, acb(0), acb(R), rel_tol=arb(2) ** -40) * ea
            tot += sa * sb * acb.integral(inner, acb(0), acb(R), rel_tol=arb(2) ** -36)
    val = -tot / (2 * arb.pi() * acb(0, 1)) ** 2
    return val
if __name__ == "__main__":
    import time
    t = time.time(); J = jhat(8); print("int_8^oo K1^2 =", J, f"[{time.time()-t:.0f}s]")
    # float cross-check: trapezoid in sigma with K1 from laplace.K
    import numpy as np
    sig = np.linspace(8, 20, 241); v = np.array([float(K(s, 1).mid())**2 for s in sig])
    print("float check:", np.trapezoid(v, sig))
