"""Exact decaying solution of the linear conical-end equation  k''' + (sigma/4) k = 0:
   K_j(sigma) = (1/(2 pi i)) int_C t^j exp(sigma t + t^4) dt  = d^j/dsigma^j K_0,
C: from oo e^{-3 i pi/4} through t0 = -(sigma/4)^{1/3} to oo e^{3 i pi/4} (straight rays from t0).
(Check: K''' + (sigma/4) K = (1/2pi i) int (t^3 + sigma/4) e^{...} = (1/2pi i) int (1/4) d/dt e^{...} = 0.)
Rigorous: arb acb.integral on r in [0, R] + explicit tail bound for r > R."""
import sys
sys.path.insert(0, __import__('os').environ.get('S1DIR', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 's1_cap')))
from flint import arb, acb, ctx
ctx.prec = 200
E3 = acb(-1, 1) / arb(2).sqrt()            # e^{3 i pi/4}
def K(sig, j, R=8):
    sig = arb(sig)
    t0 = -(sig / 4) ** (arb(1) / 3)
    def f(r, analytic):
        t = acb(t0) + r * E3
        return t ** j * (sig * t + t ** 4).exp() * E3
    I = acb.integral(f, acb(0), acb(R), rel_tol=arb(2) ** -150)
    # tail r > R (Lemma 1 of s1_tight.tex): Re(sig t + t^4) = a^4 + 2 sqrt2 a^3 r - 2 sqrt2 a r^3 - r^4 - sig a - sig r/sqrt2
    #   <= -(3/4) sig a - r^4  (a^3 = sig/4), |t| <= a + r; the exponent phi(r) = r^4 - j log(a + r) has
    #   phi' >= 4R^3 - j/(a + R) > 0 on [R, oo), so int_R^oo (a+r)^j e^{-(3/4) sig a - r^4} dr <= e^{-(3/4)sig a - phi(R)}/phi'(R).
    rR = arb(R); a0 = abs(t0)
    tail = (-(arb(3) / 4) * sig * a0 - rR ** 4).exp() * (a0 + rR) ** j / (4 * rR ** 3 - j / (a0 + rR))
    I = I + acb(arb(0, tail.upper()), arb(0, tail.upper()))
    # (1/(2 pi i)) (I - conj(I)) = Im(I)/pi
    return I.imag / arb.pi()
if __name__ == "__main__":
    import tail as T
    S = 8
    Ks = [K(S, j) for j in range(3)]
    print("K0, K1, K2 at S=8:", Ks)
    print("ratios K1/K0, K2/K0:", Ks[1] / Ks[0], Ks[2] / Ks[0])
    g = T.gamma1(S); m = T.mu(S)
    k = 1 + 2 * g.real; kp = -m * (1 + 2 * (T.W1 * g).real); kpp = m * m * (1 + 2 * (T.W2 * g).real)
    print("step-1 first-order ratios k'/k, k''/k:", kp / k, kpp / k, "  (should agree to O(kappa^2) ~ 1e-4)")
    # check the ODE: K0''' + (S/4) K0 = K3 + (S/4)K0
    print("ODE check K3 + S/4 K0 =", K(S, 3) + arb(S) / 4 * Ks[0])
