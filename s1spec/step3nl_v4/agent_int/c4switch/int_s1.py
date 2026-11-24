"""S1: perturbation bounds d = psi3~ - psi3^ and d = h2~ - h2^ (C^4 minus C^3 version):  both share the raw Chebyshev pieces and the extension, so
d = (degree-9 Hermite correction) - (degree-7 Hermite correction) on every piece, and d == 0 for x >= 28.
Outputs: ||d^(j)||_rho (half line), ||d^(j)||_{W,[10.25,60]} (P = 2 and 1.5), j = 0..4;  sup |d^(j)| over [0, 34.25] for j = 0..4.
Boxes: every piece cut into 400 sub-boxes; norm^2 = sum len * weight(x_lo) * sup(d^2) (weights decreasing)."""
import sys, json
from int_common import *
from c4fun import PSI, H2
import psi3fun as F3, h2fun as H3
NB = 400
def diff_jets(Cf, e, X):
    P = Cf.M.PIECES[e]; hv = (X - P.a) / P.Lx
    old = [P.Hd[d](hv) / P.Lx ** d for d in range(5)]          # C^3 correction (HR.Hd, degree 7)
    new = Cf.corr(e, X, 5)
    return [chk(a - b) for a, b in zip(new, old)]
def run(name, Cf):
    out = {"rho": [arb(0)] * 5, "W2": [arb(0)] * 5, "W1.5": [arb(0)] * 5}; sup = [0.0] * 5; rec = []
    for e, P in enumerate(Cf.pieces):
        a, b = P.af, P.bf; xmax = 34.25
        for i in range(NB):
            xa = a + (b - a) * i / NB; xb = a + (b - a) * (i + 1) / NB
            X = I(xa, xb); dj = diff_jets(Cf, e, X); ln = arb(xb) - arb(xa)
            rho = HR.rho_of(xa)
            for j in range(5):
                s2 = arb(up(abs(dj[j]))) ** 2
                out["rho"][j] += ln * rho * s2
                if xb > 10.25:
                    out["W2"][j] += ln * W_up('2', xa) * s2; out["W1.5"][j] += ln * W_up('1.5', xa) * s2
                sup[j] = max(sup[j], up(abs(dj[j])))
    res = {k: [up(v.sqrt()) for v in vs] for k, vs in out.items()}; res["sup_abs"] = sup
    return res
if __name__ == "__main__":
    R = {"psi3": run("psi3", PSI), "h2": run("h2", H2)}
    # note: the first W box starting below 10.25 is handled by the test xa >= 10.25-1e-12 (boxes straddling 10.25 are counted fully only if xa >= 10.25): add the straddling box conservatively
    json.dump(R, open('c4switch_S1.json', 'w'), indent=1)
    for k, v in R.items():
        for kk, vv in v.items(): print(k, kk, ["%.3e" % t for t in vv])
