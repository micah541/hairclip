"""Extension of gen_hs3 for the spectral-complement certification (CERT_PLAN.md): two more adjoint generators,
wa (lambda = 1) and wb (lambda = 1/4), alongside w (lambda = lamt). The symmetry eigenvalues are exact rationals,
so they are inlined in the vector field (no new CAPD parameters; parameter indices of gen_hs3 are unchanged).
Initial data at sigma = 32: alpha*Re(v) - Im(v), v = (1, z, z^2, z^3) the decaying WKB pair, alpha from
wstar_all.json (float shooting with even tip conditions)."""
import json, numpy as np
import gen_hs3 as G3

LAMS = {"wa": ("1", 1.0), "wb": ("1/4", 0.25)}
ALPHA = {k: json.load(open('/home/claude/s1spec/step3/wstar_all.json'))[str(v[1])]["alpha"] for k, v in LAMS.items()}


def frame(lam):
    roots = np.roots([1, 0, 0, -32.0 / 4, -(0.5 - lam)])
    z = [r for r in roots if r.real < -0.5][0]
    v = np.array([1, z, z * z, z ** 3])
    return v.real, v.imag


def init(g):
    wa, wb = frame(LAMS[g][1])
    return [float(ALPHA[g] * a - b) for a, b in zip(wa, wb)]


def ltr_l(c, w, lamtxt):
    w0, w1, w2, w3 = w
    k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
    return [f"-{w1}-{P1}*{w0}", f"-{w2}-{P1}*{w1}", f"-{w3}-{P1}*{w2}",
            f"-(-{k}^2*{w2}+({k}*{m}+{q}/4)*{w1}+(1/2-({lamtxt}))*{w0})-{P1}*{w3}"]


_aux_for, _aux_fun, _aux_init, _gexpr, _fexpr = G3.aux_for, G3.aux_fun, G3.aux_init, G3.gexpr, G3.fexpr


def aux_for(gs):
    rows = _aux_for(gs)
    for g in ("wa", "wb"):
        if g in gs: rows += [f"{g}{i}" for i in range(4)]
    return rows


def aux_fun(c, rows):
    out = _aux_fun(c, rows)
    for g in ("wa", "wb"):
        if f"{g}0" in rows: out += ltr_l(c, [f"{g}{i}" for i in range(4)], LAMS[g][0])
    return out


def aux_init(rows):
    x = _aux_init(rows)
    for g in ("wa", "wb"):
        if f"{g}0" in rows: x += [(cc, cc) for cc in init(g)]
    return x


def gexpr(c, g):
    if g in ("wa", "wb"):
        chi = f"{c['chi']}*" if c['chi'] != "1" else ""
        return chi + f"{g}0*{c['R']}"
    return _gexpr(c, g)


def fexpr(c, g):
    if g in ("wa", "wb"):
        chi = f"{c['chi']}*" if c['chi'] != "1" else ""
        return chi + f"{g}0"
    return _fexpr(c, g)


G3.aux_for, G3.aux_fun, G3.aux_init, G3.gexpr, G3.fexpr = aux_for, aux_fun, aux_init, gexpr, fexpr

if __name__ == "__main__":
    import sys, time
    g = sys.argv[1] if len(sys.argv) > 1 else "wa"
    path = f"in9Z{g}.txt"
    G3.build_Z(path, g)
    t0 = time.time(); out = G3.run(path, env={"REORG": "1e300"})
    open(f"out9Z{g}.txt", "w").write(out)
    print(g, "RESULT OK" in out, f"{time.time() - t0:.0f}s")
