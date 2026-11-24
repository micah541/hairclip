"""Certified slow eigenfunction psi3^ (CERT_PLAN.md, 13:55): rows u0..u3 solve (L - lamt)u = 0 inward from
sigma = 32, with data = exact jet of g(sigma) = sigma^A sum_{j<4} a_j sigma^(-4j) + beta e^{S(sigma)},
S' = -(sigma/4)^(1/3). Integral rows give <psi3^, f_i>-type quantities on [8,32] (B*) and [tip,8] (C*), as in build_A2.
psi3^ := U on [tip, 32] and g beyond (C^3 at 32), minus tip bumps (from the certified tip values)."""
import json, numpy as np
from flint import arb
import gen_hs3 as G3
from gen_hs3 import coeffs, write, ltr, wstar_init

D = json.load(open('/home/claude/s1spec/step3nl/psi3beta.json'))
BETA = D["beta"]; LT = G3.LAMT


def g_jet():
    """arb enclosure of the jet (g, g', g'', g''') at sigma = 32."""
    s = arb(32); A = 1 - 4 * arb(LT)
    a = [arb(1)]
    for j in range(1, 4):
        a.append(a[-1] * (A - 4 * j + 4) * (A - 4 * j + 3) * (A - 4 * j + 2) * (A - 4 * j + 1) / j)
    jet = [arb(0)] * 4
    for j in range(4):
        e = A - 4 * j; f = arb(1)
        for n in range(4):
            jet[n] += a[j] * f * s ** (e - n); f = f * (e - n)
    z = -(s / 4) ** (arb(1) / 3); z1 = -(arb(1) / 12) * (s / 4) ** (-arb(2) / 3)
    z2 = (arb(1) / 72) * (s / 4) ** (-arb(5) / 3) * arb(2) / 3 * 3 / 2      # d/ds of z1 = (1/18)(1/4)(s/4)^(-5/3)
    z2 = (arb(1) / 72) * (s / 4) ** (-arb(5) / 3)
    ex = [arb(1), z, z1 + z * z, z2 + 3 * z * z1 + z ** 3]
    b = arb(BETA)
    return [jet[n] + b * ex[n] for n in range(4)]


def lrow(c, u):
    u0, u1, u2, u3 = u
    k, m, q, p = c['k'], c['m'], c['q'], c['p']
    return [f"-{u1}", f"-{u2}", f"-{u3}",
            f"-(-{k}^2*{u2}-(5*{k}*{m}+{q}/4)*{u1}+(1/4-3*{m}^2-{k}*{p}-lamt)*{u0})"]


P3_INT = [f"B{i}" for i in range(7)] + [f"C{i}" for i in range(7)]


def build_P3(path):
    aux = ["C", "S", "Ct2", "St2", "ws0", "ws1", "ws2", "ws3", "u0", "u1", "u2", "u3"]

    def fun(c, stage):
        k = c['k']; chi = c['chi']; R = c['R']; P = c['p']
        if stage == "B":
            out = [f"-S*{k}", f"C*{k}", "0", "0"]; cs = ("C", "S")
        else:
            out = ["0", "0", f"-St2*{k}", f"Ct2*{k}"]; cs = ("Ct2", "St2")
        out += ltr(c, ["ws0", "ws1", "ws2", "ws3"]) + lrow(c, ["u0", "u1", "u2", "u3"])
        ch = f"({chi})*" if chi != "1" else ""
        Cv, Sv = cs
        e1 = c.get('e1', "0"); e3 = c.get('e3', "0")
        ints = [f"{ch}ws0*u0*{R}", f"{ch}{P}*u0*{R}", f"{ch}{Cv}*u0*{R}", f"{ch}{Sv}*u0*{R}", f"u0^2*{R}",
                "0" if e1 == "0" else f"{e1}*u0*{R}", "0" if e3 == "0" else f"{e3}*u0*{R}"]
        z = ["0"] * 7
        return out + (ints + z if stage == "B" else z + ints)

    gj = g_jet()
    x = ([(1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (0.0, 0.0)] + [(cc, cc) for cc in wstar_init()]
         + [(float(v.lower()), float(v.upper())) for v in gj] + [(0.0, 0.0)] * len(P3_INT))
    return write(path, aux + P3_INT, fun(coeffs("B"), "B"), fun(coeffs("C"), "C"), fun(coeffs("C2"), "C"), x, [], [],
                 [(0.0, 0.0)] * len(x))


if __name__ == "__main__":
    import time
    print("g jet at 32:", [str(v) for v in g_jet()])
    build_P3("in11P3.txt")
    t0 = time.time(); out = G3.run("in11P3.txt", env={"REORG": "1e300"})
    open("out11P3.txt", "w").write(out)
    print("RESULT OK" in out, f"{time.time() - t0:.0f}s")
    for line in out.splitlines():
        if line.startswith(("TIP", "TX", "RESULT")): print(line[:300])
