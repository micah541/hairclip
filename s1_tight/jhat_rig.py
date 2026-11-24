"""Rigorous Jhat = int_S^oo K1^2 : double contour integral on [0,R]^2 (jhat.py) + explicit remainder for r > R or rho > R.
On both rays t = t0 + r e^{+-3 i pi/4}, t0 = -a, a = (S/4)^(1/3):
   Re(t^4 + S t) = a^4 + 2 sqrt2 a^3 r - 2 sqrt2 a r^3 - r^4 - S a - S r/sqrt2  <=  c0 - r^4,   c0 = a^4 - S a
(the r-linear terms cancel because a^3 = S/4), |t| <= a + r, |t + s| >= |Re(t + s)| >= 2a.
So |integrand| <= g(r) g(rho)/(2a), g(r) = (a + r) e^{c0 - r^4}, and the remainder over the 4 ray pairs is
   <= (1/(4 pi^2)) * 4 * 2 G_R G_0/(2a),   G_0 = int_0^oo g = e^{c0}(a Gamma(5/4) + sqrt(pi)/4),
   G_R = int_R^oo g <= e^{c0 - phi(R)}/phi'(R),  phi(r) = r^4 - log(a + r),  phi' >= 4R^3 - 1/(a + R)."""
import json, sys
from flint import arb, acb, ctx
from jhat import jhat
ctx.prec = 80
def remainder(S, R):
    S = arb(S); R = arb(R); a = (S / 4) ** (arb(1) / 3); c0 = a ** 4 - S * a
    G0 = c0.exp() * (a * arb(5 / 4).gamma() + arb.pi().sqrt() / 4)
    GR = (c0 - R ** 4 + (a + R).log()).exp() / (4 * R ** 3 - 1 / (a + R))
    return 4 * 2 * G0 * GR / (2 * a) / (4 * arb.pi() ** 2)
if __name__ == "__main__":
    S, R = 8, 6
    J = jhat(S, R); rem = remainder(S, R)
    val = J.real + arb(0, (abs(J.imag).upper() + rem.upper()))      # Jhat real; imaginary part is enclosure noise
    print("Jhat =", val, "  remainder bound", rem)
    json.dump({"int_K1sq_str": val.str(30, radius=True), "int_K1sq": [float(val.mid()), float(val.rad()) * 1.01 + 1e-30],
               "S": S, "R": R, "remainder": float(rem.upper())}, open("jhat8.json", "w"), indent=1)
