"""
Rigorous tail (conical end) enclosure for the CDF shrinker system
    p' = -k q,  q' = 1 + k p,  k' = m,  m' = p/4        (theta' = k)
on [S, oo), normalised so that q(sigma) - sigma -> 0.

Reduction (see writeup):  k''' + (sigma/4) k = N,  N = -k (k k' + J),  J = int_sigma^oo k'^2,
with m = k', p = 4 k'', q = sigma + 4 k k' + 4 J.

Diagonalised coordinates: Y=(k,k',k'') = P(sigma) Z,  P = diag(1,mu,mu^2) V,  V_ij = (-w_j)^i,
mu = (sigma/4)^(1/3),  w_j = exp(2 pi i j/3).  Then
    Z_j' = (lambda_j - 1/(3 sigma)) Z_j - (1/(3 sigma)) sum_{l!=j} C_jl Z_l + (w_j/(3 mu^2)) N,
lambda_j = -w_j mu,  C_jl = (1/3)(zeta + 2 zeta^2), zeta = w_l / w_j, |C_jl| = 1/sqrt(3).

Fixed point in the weighted norms
    ||Z0||_0 = sup w |Z0|,  ||Z1||_1 = sup w g |Z1|,  w = e^{E}(sigma/S)^{1/3},  g = (sigma/S)^{4/3},
    E(sigma) = int_S^sigma mu = a (sigma^{4/3} - S^{4/3}),  a = (3/4) 4^{-1/3}.
All constants below are evaluated in ball arithmetic (python-flint / arb).
"""
from flint import arb, acb, ctx

ctx.prec = 256

THIRD = arb(1) / 3
A_E = arb(3) / 4 * arb(4) ** (-THIRD)          # a = (3/4) 4^{-1/3}
C_ABS = 1 / arb(3).sqrt()                       # |C_jl|
W1 = acb(-arb(1) / 2, arb(3).sqrt() / 2)        # omega_1
W2 = acb(-arb(1) / 2, -arb(3).sqrt() / 2)       # omega_2 = omega_1^2
C10 = acb(-arb(1) / 2, arb(3).sqrt() / 6)       # C_10 (l-j = 2 mod 3)


def mu(s):
    return (arb(s) / 4) ** THIRD


def constants(S):
    S = arb(S)
    m = mu(S)
    kappa = 2 * C_ABS / (9 * S * m)
    n = 1 / (2 * m) + 1 / (36 * S * m * m)
    nu0 = n / (2 * m)
    nu1 = arb(2) / 7 * n / m
    return dict(S=S, mu=m, kappa=kappa, n=n, nu0=nu0, nu1=nu1)


def radii(S, eps_max, inflate=arb("1.0001")):
    """Find r0, r1 with  r0 >= (c/2) r1 + nu0 R^3,  r1 >= kappa (eps+r0+r1) + nu1 R^3,
    R = eps + r0 + 2 r1, and contraction constant L < 1.  Verified in ball arithmetic."""
    K = constants(S)
    e = arb(eps_max)
    r1 = K["kappa"] * e / (1 - K["kappa"])
    r0 = C_ABS / 2 * r1
    for _ in range(60):
        R = e + r0 + 2 * r1
        r1n = (K["kappa"] * (e + r0 + r1) + K["nu1"] * R ** 3) * inflate
        r0n = (C_ABS / 2 * r1n + K["nu0"] * R ** 3) * inflate
        r0, r1 = r0n.upper(), r1n.upper()   # use exact upper bounds as the radii
        r0, r1 = arb(r0), arb(r1)
    R = e + r0 + 2 * r1
    ok0 = (C_ABS / 2 * r1 + K["nu0"] * R ** 3) <= r0
    ok1 = (K["kappa"] * (e + r0 + r1) + K["nu1"] * R ** 3) <= r1
    L = C_ABS / 2 + K["kappa"] + 6 * (K["nu0"] + K["nu1"]) * R ** 2
    okL = L < 1
    delta1 = K["kappa"] * (r0 + r1) + K["nu1"] * R ** 3
    Jbar = R ** 2 * (K["mu"] / 2 + 1 / (12 * K["S"]))
    return dict(r0=r0, r1=r1, R=R, L=L, delta1=delta1, Jbar=Jbar,
                verified=bool(ok0 and ok1 and okL), **K)


def gamma1(S, cut=None):
    """gamma^(1) = (C_10/3) int_S^oo tau^{-1} exp(-(1-omega_1) E(tau)) dtau, rigorous.
    The integral over [S, cut] uses arb's validated integration; the tail beyond cut is bounded
    by (1/cut) int_cut^oo e^{-(3/2)E} <= (2/3) e^{-(3/2)E(cut)} / (cut mu(cut))."""
    S = arb(S)
    if cut is None:
        cut = S + 30
    cut = arb(cut)
    s43 = S ** (arb(4) / 3)
    beta = acb(arb(3) / 2, -arb(3).sqrt() / 2)      # 1 - omega_1

    def f(x, analytic):
        # x^{4/3} via exp((4/3) log x); holomorphic for Re x > 0
        if analytic and not (x.real > 0):
            return acb("nan")
        E = A_E * ((acb(x).log() * (arb(4) / 3)).exp() - s43)
        return (-beta * E).exp() / x

    I = acb.integral(f, acb(S), acb(cut), rel_tol=arb(2) ** (-200))
    Ecut = A_E * (cut ** (arb(4) / 3) - s43)
    tail = arb(2) / 3 * (-(arb(3) / 2) * Ecut).exp() / (cut * mu(cut))
    I = I + acb(arb(0, tail), arb(0, tail))
    return C10 / 3 * I


def initial_state(S, eps, g1, delta1, Jbar):
    """Return (center_state, directions, radii) describing the enclosure of the conical-end
    state at sigma=S for parameter eps (an arb point or ball):
        state = center + d1*v_re + d2*v_im + dq*e_q,  |d1|,|d2| <= delta1',  dq in [0, 4 Jbar] + q-slack.
    Ordering of state: (p, q, k, m, theta)."""
    S = arb(S)
    m_ = mu(S)
    Z1 = eps * g1                                 # acb
    k = eps + 2 * Z1.real
    kp = -m_ * (eps + 2 * (W1 * Z1).real)
    kpp = m_ ** 2 * (eps + 2 * (W2 * Z1).real)
    p = 4 * kpp
    # direction vectors for d1 = Re dZ1, d2 = Im dZ1 (exact linear dependence of (k,k',k''))
    #  2 Re(w (d1 + i d2)) = 2 Re(w) d1 - 2 Im(w) d2
    vk = (2, 0)
    vkp = (-m_ * 2 * W1.real, m_ * 2 * W1.imag)
    vkpp = (m_ ** 2 * 2 * W2.real, -m_ ** 2 * 2 * W2.imag)
    v_re = [4 * vkpp[0], arb(0), arb(vk[0]), vkp[0], arb(0)]
    v_im = [4 * vkpp[1], arb(0), arb(vk[1]), vkp[1], arb(0)]
    return dict(p=p, k=k, m=kp, v_re=v_re, v_im=v_im)
