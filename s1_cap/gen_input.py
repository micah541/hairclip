"""
Generate a CAPD input file whose initial set  x + C r0  contains the rigorous conical-end
states at sigma = S for every eps in [eps_lo, eps_hi]   (tail theorem in tail.py).

Usage: python3 gen_input.py S eps_lo eps_hi tau1 outfile
"""
import sys, math
from flint import arb, acb, ctx
from tail import radii, gamma1, mu, W1, W2  # noqa

ctx.prec = 256


def dn(a):
    """largest double <= every point of ball a"""
    lo = a.lower()
    v = float(lo)
    while arb(v) > lo:
        v = math.nextafter(v, -math.inf)
    return v


def up(a):
    hi = a.upper()
    v = float(hi)
    while arb(v) < hi:
        v = math.nextafter(v, math.inf)
    return v


def dbl(a):
    """an exact double near the midpoint of a"""
    return float(a.mid())


def build(S, eps_lo, eps_hi, tau1):
    S = arb(S)
    elo, ehi = arb(eps_lo), arb(eps_hi)
    emax = max(abs(float(eps_lo)), abs(float(eps_hi))) * (1 + 1e-12)
    T = radii(S, arb(emax))
    assert T["verified"], "tail fixed point not verified"
    g = gamma1(S)
    gc = acb(arb(dbl(g.real)), arb(dbl(g.imag)))
    rg = max(up(abs(g.real - gc.real)), up(abs(g.imag - gc.imag)))
    deltap = T["delta1"] + arb(emax) * arb(rg) * 2          # >= delta1 + emax*|gamma-gc|
    deltap = arb(up(deltap))
    m_ = mu(S)
    # derivative of (p,k,m) w.r.t. eps at Z1 = eps*gc
    u_p = 4 * m_ ** 2 * (1 + 2 * (W2 * gc).real)
    u_k = 1 + 2 * gc.real
    u_m = -m_ * (1 + 2 * (W1 * gc).real)
    # directions for d1 = Re dZ1, d2 = Im dZ1
    vre = [4 * m_ ** 2 * 2 * W2.real, arb(0), arb(2), -m_ * 2 * W1.real, arb(0)]
    vim = [-4 * m_ ** 2 * 2 * W2.imag, arb(0), arb(0), m_ * 2 * W1.imag, arb(0)]
    ec = arb(dbl((elo + ehi) / 2))
    Delta = arb(up(max(abs(ehi - ec), abs(elo - ec))))
    # box enclosures of k, m over the whole parameter set (for q)
    E = arb(0, 1)  # [-1,1]
    eps_box = ec + Delta * E
    d_box = deltap * E
    k_box = eps_box * u_k + d_box * vre[2] + d_box * vim[2]
    m_box = eps_box * u_m + d_box * vre[3] + d_box * vim[3]
    q_box = S + 4 * k_box * m_box + 4 * T["Jbar"] * (arb(1) + E) / 2
    qc = arb(dbl(q_box))
    dq_lo, dq_hi = q_box.lower() - qc, q_box.upper() - qc
    # exact centre (balls) and exact columns (balls)
    xc = [ec * u_p, qc, ec * u_k, ec * u_m, arb(0)]
    cols = [[u_p, arb(0), u_k, u_m, arb(0)], vre, vim,
            [arb(0), arb(1), arb(0), arb(0), arb(0)],
            [arb(0), arb(0), arb(0), arb(0), arb(1)]]
    r0 = [arb(0, 1) * Delta, d_box, d_box, None, arb(0)]
    r0abs = [Delta, deltap, deltap, max(abs(dq_lo), abs(dq_hi)), arb(0)]
    # double-rounded matrix and the induced error added to the centre interval
    Ct = [[dbl(cols[j][i]) for j in range(5)] for i in range(5)]
    xlo, xhi = [], []
    for i in range(5):
        err = arb(0)
        for j in range(5):
            err += abs(cols[j][i] - arb(Ct[i][j])) * r0abs[j]
        xi = xc[i] + arb(0, 1) * err
        xlo.append(dn(xi)); xhi.append(up(xi))
    r0lo = [-up(Delta), -up(deltap), -up(deltap), dn(dq_lo), 0.0]
    r0hi = [up(Delta), up(deltap), up(deltap), up(dq_hi), 0.0]
    info = dict(S=S, eps_c=ec, Delta=Delta, deltap=deltap, kappa=T["kappa"], r0=T["r0"], r1=T["r1"],
                L=T["L"], Jbar=T["Jbar"], gamma=g, q_box=q_box)
    return tau1, xlo, xhi, Ct, r0lo, r0hi, info


def write(path, tau1, xlo, xhi, Ct, r0lo, r0hi):
    with open(path, "w") as f:
        f.write(repr(float(tau1)) + "\n")
        for a, b in zip(xlo, xhi):
            f.write(f"{a!r} {b!r}\n")
        for i in range(5):
            f.write(" ".join(repr(Ct[i][j]) for j in range(5)) + "\n")
        for a, b in zip(r0lo, r0hi):
            f.write(f"{a!r} {b!r}\n")


if __name__ == "__main__":
    S, elo, ehi, tau1, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
    tau1, xlo, xhi, Ct, r0lo, r0hi, info = build(S, elo, ehi, tau1)
    write(out, tau1, xlo, xhi, Ct, r0lo, r0hi)
    print("wrote", out, " eps_c=", float(info["eps_c"].mid()), " Delta=", float(info["Delta"].mid()),
          " delta'=", float(info["deltap"].mid()), " q range=", float(info["q_box"].lower()), float(info["q_box"].upper()))
