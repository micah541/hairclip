"""Global Chebyshev refit of the multi-domain psi3^ on x in [0, XC] (stage C range) and its float residual."""
import sys, json, numpy as np
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
J = json.load(open('psi3md_14.json')); lam = J["lam"]
def val(x):
    out = np.empty_like(x)
    for pc in J["pieces"]:
        sel = (x >= pc["a"]) & (x <= pc["b"])
        out[sel] = Ch.chebval(1 - 2 * (x[sel] - pc["a"]) / (pc["b"] - pc["a"]), pc["coef"])
    return out
XC = 8 - E.STIP + 0.01
for deg in [40, 50, 60, 70]:
    t = np.cos(np.pi * (np.arange(4 * deg) + 0.5) / (4 * deg)); x = (1 - t) * XC / 2
    cf = Ch.chebfit(t, val(x), deg)
    xf = np.linspace(0, XC, 8001); tf = 1 - 2 * xf / XC
    dd = [Ch.chebval(tf, Ch.chebder(cf, n)) * (-2 / XC) ** n for n in range(5)]
    p, q, k, m = np.array([E.PROF(s + E.STIP) for s in xf]).T
    r = -dd[4] - k * k * dd[2] - (5 * k * m + q / 4) * dd[1] + (0.25 - 3 * m * m - k * p - lam) * dd[0]
    rho = np.exp(-0.25 * (1 + xf**2) ** (2 / 3))
    print(f"deg {deg}: ||r||_rho on [0,{XC:.2f}] = {np.sqrt(2*np.trapezoid(rho*r*r, xf)):.2e}, max weighted {np.max(np.abs(r)*np.sqrt(rho)):.2e}, |c_last| = {abs(cf[-1]):.1e}")
    json.dump({"XC": XC, "lam": lam, "coef": cf.tolist()}, open(f'psi3glob_{deg}.json', 'w'))
