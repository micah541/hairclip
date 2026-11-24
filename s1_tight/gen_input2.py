"""
Tight CAPD input: parallelepiped x + C r0 containing the conical-end states at sigma = S for every
eps in [eps_lo, eps_hi], using the EXACT linear tail solution (tail2.py / laplace.py) and the cubic
Lyapunov-Perron remainder only.

State (p, q, k, m, th) at sigma = S:
   (k, k', k'') = eps * Yhat(S) + P(S) (0, d, conj d),   |d| <= rn1   (d = Z_nl,1(S), Z_nl,0(S) = 0)
   p = 4 k'',  m = k',  th = 0 (convention),
   q = S + 4 k k' + 4 J,   J = int_S^oo k'^2 = eps^2 Jhat + J',
   Jhat = cn^2 int_S^oo K1^2  (jhat.py),  |J'| <= 2|eps| sqrt(Jhat) B + B^2,
   B = (rn0 + 2 rn1) sqrt(mu(S)/2 + 1/(12 S))  >=  || k'_nl ||_{L^2(S,oo)}.
Columns of C: eps-direction u (including dq/deps at the centre), Re d, Im d (with their q-components at
the centre), e_q (quadratic remainder + J'), e_th.

Usage: python3 gen_input2.py S eps_lo eps_hi tau1 outfile
"""
import sys, math, json, os
sys.path.insert(0, __import__('os').environ.get('S1DIR', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 's1_cap')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import tail as T
from tail2 import nl_radii, yhat, Pmat
from laplace import K

ctx.prec = 256


def dn(a):
    lo = a.lower(); v = float(lo)
    while arb(v) > lo:
        v = math.nextafter(v, -math.inf)
    return v


def up(a):
    hi = a.upper(); v = float(hi)
    while arb(v) < hi:
        v = math.nextafter(v, math.inf)
    return v


def dbl(a):
    return float(a.mid())


_cache = {}


def tail_data(S, emax):
    key = (str(S), repr(emax))
    if key in _cache:
        return _cache[key]
    Sa = arb(S)
    R = nl_radii(Sa, arb(emax))
    assert R["verified"], "nonlinear LP remainder not verified"
    Y, Z = yhat(Sa)
    # Yhat is real: imaginary parts are pure enclosure noise; real parts with the imaginary radius absorbed
    Yr = [y.real + arb(0, abs(y.imag).upper()) for y in Y]
    # normalisation constant cn = 1/Zhat_0 used inside yhat: recompute from K0
    K0 = acb(K(Sa, 0))
    cn = (Y[0] / K0)
    cnr = cn.real + arb(0, abs(cn.imag).upper())
    jdat = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "jhat%d.json" % int(S))))
    I1 = arb(jdat["int_K1sq"][0], max(jdat["int_K1sq"][1], 1e-18))   # radius widened: covers decimal rounding of the stored midpoint
    Jhat = cnr ** 2 * I1
    m_ = T.mu(Sa)
    B = (R["rn0"] + 2 * R["rn1"]) * (m_ / 2 + 1 / (12 * Sa)).sqrt()
    out = dict(R=R, Y=Yr, cn=cnr, Jhat=Jhat, B=B, mu=m_)
    _cache[key] = out
    return out


def build(S, eps_lo, eps_hi, tau1, emax=None):
    Sa = arb(S)
    elo, ehi = arb(eps_lo), arb(eps_hi)
    em = max(abs(float(eps_lo)), abs(float(eps_hi))) * (1 + 1e-12)
    if emax is None:
        emax = em
    assert emax >= em
    D = tail_data(S, emax)
    Y0, Y1, Y2 = D["Y"]
    Jhat, B, m_ = D["Jhat"], D["B"], D["mu"]
    rn1 = D["R"]["rn1"]
    dmax = arb(up(rn1))
    ec = arb(dbl((elo + ehi) / 2))
    Delta = arb(up(max(abs(ehi - ec), abs(elo - ec))))
    E = arb(0, 1)
    # d-directions (exact linear dependence of (k, k', k'') on d1 = Re d, d2 = Im d)
    vk = (arb(2), arb(0))
    vm = (-m_ * 2 * T.W1.real, m_ * 2 * T.W1.imag)
    vpp = (m_ ** 2 * 2 * T.W2.real, -m_ ** 2 * 2 * T.W2.imag)
    vp = (4 * vpp[0], 4 * vpp[1])
    G = Y0 * Y1 + Jhat                      # q = S + 4 eps^2 G + (d-terms) + 4 J'
    u_q = 8 * ec * G
    vre_q = 4 * ec * (Y0 * vm[0] + Y1 * vk[0])
    vim_q = 4 * ec * (Y0 * vm[1] + Y1 * vk[1])
    u = [4 * Y2, u_q, Y0, Y1, arb(0)]
    vre = [vp[0], vre_q, vk[0], vm[0], arb(0)]
    vim = [vp[1], vim_q, vk[1], vm[1], arb(0)]
    # remainder in q over the parameter box: e in [-Delta, Delta], d1, d2 in [-dmax, dmax], eps = ec + e
    e = Delta * E; d1 = dmax * E; d2 = dmax * E
    eps_box = ec + e
    Jerr = 2 * abs(eps_box).upper() * Jhat.sqrt() * B + B * B
    dk = vk[0] * d1 + vk[1] * d2
    dm = vm[0] * d1 + vm[1] * d2
    # 4 (eps Y0 + dk)(eps Y1 + dm) = 4 eps^2 Y0 Y1 + 4 eps (Y0 dm + Y1 dk) + 4 dk dm
    rem = (4 * e * e * G
           + 4 * e * (Y0 * dm + Y1 * dk)          # (eps - ec) part of the linear d-terms
           + 4 * dk * dm
           + 4 * arb(0, 1) * Jerr)
    qc_ball = Sa + 4 * ec * ec * G
    qc = arb(dbl(qc_ball))
    dq = rem + (qc_ball - qc)
    xc = [ec * u[0], qc, ec * u[2], ec * u[3], arb(0)]
    cols = [u, vre, vim, [arb(0), arb(1), arb(0), arb(0), arb(0)], [arb(0)] * 4 + [arb(1)]]
    r0abs = [Delta, dmax, dmax, max(abs(dq.lower()), abs(dq.upper())), arb(0)]
    Ct = [[dbl(cols[j][i]) for j in range(5)] for i in range(5)]
    xlo, xhi = [], []
    for i in range(5):
        err = arb(0)
        for j in range(5):
            err += abs(cols[j][i] - arb(Ct[i][j])) * r0abs[j]
        xi = xc[i] + arb(0, 1) * err
        xlo.append(dn(xi)); xhi.append(up(xi))
    r0lo = [-up(Delta), -up(dmax), -up(dmax), dn(dq), 0.0]
    r0hi = [up(Delta), up(dmax), up(dmax), up(dq), 0.0]
    info = dict(eps_c=ec, Delta=Delta, dmax=dmax, dq=dq, Jhat=Jhat, Jerr=Jerr, Y=D["Y"], cn=D["cn"],
                rn0=D["R"]["rn0"], rn1=rn1, L=D["R"]["L"])
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
    print("wrote", out, " eps_c=", float(info["eps_c"].mid()), " Delta=%.3e" % float(info["Delta"].mid()),
          " dmax=%.3e" % float(info["dmax"].mid()), " dq=[%.3e,%.3e]" % (float(info["dq"].lower()), float(info["dq"].upper())),
          " Jhat=", info["Jhat"], " Jerr=%.3e" % float(info["Jerr"].upper()))
