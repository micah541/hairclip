"""F2: sup|1 - q psi'|, sup|(W cos theta)'/W| = sup|-psi' cos(theta) - k sin(theta)| <= sup sqrt(psi'^2 + k^2) (frame angle arbitrary), and over the inner ramp
    [10.25, 16.25] sup(q zeta' W/rho)  (zeta' = S'((x-10.25)/6)/6, S'(t) = 140 t^3 (1-t)^3, W/rho = exp(phi - psi) increasing there: certified per box).
F3: W-norms of zeta*A1(u0)'''/|alpha|, zeta*B1(u0)'''/|alpha|, zeta*psi3~''', zeta*h2~''' on [10.25, R], R = 40, 60 (P = 2).
  A1(u) = -(u - u' q/al_g), B1(u) = u' cos(theta)/al_g, al_g = 1 - k u, theta' = k.  u0^(n) = A yh_n, yh_n = psi3~^(n) + A h2~^(n), A = [-amax, amax];
  Taylor series in the x-increment eps (order 3): profile jets from the ODE (exact recursion, balls), cos(theta0 + Theta) = c0 cos Theta - s0 sin Theta with
  Theta = int k, |B1'''| <= 6 sqrt(X3^2 + Y3^2), X = (u'/al_g) cos Theta, Y = (u'/al_g) sin Theta  (|c0|^2 + |s0|^2 = 1).  (A1', B1' coefficient-1 values are also reported
  as a check against v2 V5.)   Usage: python3 int_f23.py [width_scale]"""
import sys
from int_common import *
from flint import arb_series
from c4fun import PSI, H2
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
PK = '2'; AMAX = 1e-3; A = arb(0, AMAX)
def Spoly_f(t):
    t = min(max(t, 0.0), 1.0); a = arb(t); return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb, R):
    zi = up(Spoly_f((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = up(Spoly_f((R - xa) / 3)) if xa > R - 3 else 1.0
    return min(1.0, zi) * min(1.0, zo)
def profile_coeffs(k0, m0, p0, q0):
    k1 = m0; m1 = p0 / 4; p1 = -k0 * q0; q1 = 1 + k0 * p0
    k2 = m1 / 2; m2 = p1 / 8; p2 = -(k0 * q1 + k1 * q0) / 2; q2 = (k0 * p1 + k1 * p0) / 2
    k3 = m2 / 3; m3 = p2 / 12; p3 = -(k0 * q2 + k1 * q1 + k2 * q0) / 3; q3 = (k0 * p2 + k1 * p1 + k2 * p0) / 3
    return [k0, k1, k2, k3], [m0, m1, m2, m3], [p0, p1, p2, p3], [q0, q1, q2, q3]
def f3_quantities(b, y, h):
    fact = [1, 1, 2, 6]
    yh = [y[n] + A * h[n] for n in range(5)]
    u = arb_series([yh[n] / fact[n] for n in range(4)], prec=4); up_ = arb_series([yh[n + 1] / fact[n] for n in range(4)], prec=4)
    kc, mc, pc, qc = profile_coeffs(b.k, b.m, b.p, b.q)
    K = arb_series(kc, prec=4); Q = arb_series(qc, prec=4)
    ag = 1 - K * A * u; assert ag.coeffs()[0].lower() > 0
    Z = up_ / ag
    A1 = -(u - Z * Q)
    Th = arb_series([arb(0), kc[0], kc[1] / 2, kc[2] / 3], prec=4)
    Xs = Z * Th.cos(); Ys = Z * Th.sin()
    cA = A1.coeffs() + [arb(0)] * 4; cX = Xs.coeffs() + [arb(0)] * 4; cY = Ys.coeffs() + [arb(0)] * 4
    def mag(a_, c_):
        ua, uc = arb(up(abs(a_))), arb(up(abs(c_))); return (ua * ua + uc * uc).sqrt()
    d = {"A1'''_over_alpha": 6 * abs(cA[3]), "B1'''_over_alpha": 6 * mag(cX[3], cY[3]), "psi3'''": abs(y[3]), "h2'''": abs(h[3]),
         "A1'_over_alpha": abs(cA[1]), "B1'_over_alpha": mag(cX[1], cY[1])}
    return {kk: up(v) for kk, v in d.items()}
def f2_quantities(b):
    P = Pj_box(b.X, PK)
    d = {"|1-q psi'|": up(abs(1 - b.q * P[0]))}
    u1 = arb(up(abs(P[0]))); kb = arb(b.kb.upper())
    d["sqrt(psi'^2+k^2)"] = up((u1 * u1 + kb * kb).sqrt()); d["psi'"] = up(P[0])
    if b.xa < 16.25 + 1e-9:
        xb = b.xb; X = b.X
        phi1 = X * ((arb(1) + X * X) ** (-arb(1) / 3)) / 3
        assert (phi1 - P[0]).lower() > 0, 'W/rho not monotone on box'
        t = (X - arb(10.25)) / 6; t = t.union(t) if False else t
        tt = arb(max(0.0, min(1.0, float((arb(b.xa) - arb(10.25)).lower()) / 6))).union(arb(max(0.0, min(1.0, float((arb(b.xb) - arb(10.25)).upper()) / 6))))
        om = 1 - tt
        zp = 140 * (tt * tt * tt) * (om * om * om) / 6
        phi_xb = (arb(1) + arb(xb) * arb(xb)) ** (arb(2) / 3) / 4
        wr = (phi_xb - psi_value(PK, math.nextafter(xb, UP))).exp()      # >= W/rho at xb' >= xb, increasing => sup over the box
        d["zeta'"] = up(zp); d["W/rho"] = up(wr); d["q zeta' W/rho"] = up(abs(b.q) * zp * wr)
    return d
def run():
    t0 = time.time(); res = []
    bl = make_segments([('x', 10.25), ('x', 16.25), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)],
                       [0.01 * WS, 0.01 * WS, 0.0125 * WS, 0.0125 * WS, 0.025 * WS])
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
        for v in y + h: chk(v)
        res.append(dict(xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)), W=float(W_up(PK, b.xa).mid()), f3=f3_quantities(b, y, h), f2=f2_quantities(b),
                        z={str(R): zeta_sup(b.xa, b.xb, R) for R in (40.0, 60.0)}))
    print(f"[{time.time()-t0:.0f}s] {len(res)} boxes"); return res
def summarize(res):
    out = {}
    for Rr in (40.0, 60.0):
        o = {}
        for nm in res[0]['f3']:
            tz = arb(0); tn = arb(0)
            for r in res:
                if r['xb'] > Rr + 1e-6: continue
                v = arb(r['f3'][nm]); base = arb(r['len']) * arb(r['W']); z = arb(r['z'][str(Rr)])
                tz += base * z * z * v * v; tn += base * v * v
            o[nm] = dict(norm_zeta=up(tz.sqrt()), norm_nozeta=up(tn.sqrt()))
        out[f"R={Rr:g}"] = o
    sup = {}
    for r in res:
        for ri in region_of(r['xa'], r['xb']):
            s = sup.setdefault(f"[{REG[ri][0]:g},{REG[ri][1]:g}]", {})
            for dd in (r['f3'], r['f2']):
                for nm, v in dd.items(): s[nm] = max(s.get(nm, -1e300), v)
    out["sup_per_region"] = sup
    # F2 cumulative sups
    for Rr in (40.0, 60.0):
        c = {}
        for r in res:
            if r['xb'] > Rr + 1e-6: continue
            for nm, v in r['f2'].items(): c[nm] = max(c.get(nm, -1e300), v)
        out[f"F2_sup_R={Rr:g}"] = c
    return out
if __name__ == "__main__":
    res = run(); S = summarize(res)
    json.dump({"width_scale": WS, **S}, open(f'v3_f23_{WS:g}.json', 'w'), indent=1)
    for k_ in ("R=40", "R=60"): print(k_, json.dumps({a: (f"{b['norm_zeta']:.4g}", f"{b['norm_nozeta']:.4g}") for a, b in S[k_].items()}))
    print(json.dumps(S["F2_sup_R=40"]), json.dumps(S["F2_sup_R=60"]))
    for r_, v in S["sup_per_region"].items(): print(r_, {a: float(f"{b:.4g}") for a, b in v.items()})
