"""Plan B (c)/(e): modulation sources, F2/F3-type constants on [10.25, 260] (|alpha| <= 1e-3), W-norms with zeta (box sup; inner S((x-10.25)/6), outer S((260-x)/10)).
 psi1' = -k q ; psi2' <= |k| ; psi3~' ; h2~' ; (A1(u0))'/alpha, (B1(u0))'/alpha, (A1(u0))'''/alpha, (B1(u0))'''/alpha ; psi3~''', h2~''' ; c0' = -(5/2) m p + k^2 q ;
 q psi3~'' ;  |1 - q psi'| ; sqrt(psi'^2 + k^2) (>= |(W cos th)'/W|) ; q zeta' W/rho on the inner ramp.   A1 = -(u - u' q/al_g), B1 = u' cos th/al_g, theta' = k.
 Taylor series of order 3 in the x-increment (profile jets from the ODE), cos/sin via (c0,s0) frame: |B1'''| <= 6 sqrt(X3^2 + Y3^2)."""
import sys
from pb_common import *
from flint import arb_series
from c4fun import PSI, H2
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
AMAX = 1e-3; A = arb(0, AMAX)
def profile_coeffs(k0, m0, p0, q0):
    k1 = m0; m1 = p0 / 4; p1 = -k0 * q0; q1 = 1 + k0 * p0
    k2 = m1 / 2; m2 = p1 / 8; p2 = -(k0 * q1 + k1 * q0) / 2; q2 = (k0 * p1 + k1 * p0) / 2
    k3 = m2 / 3; m3 = p2 / 12; p3 = -(k0 * q2 + k1 * q1 + k2 * q0) / 3; q3 = (k0 * p2 + k1 * p1 + k2 * p0) / 3
    return [k0, k1, k2, k3], [m0, m1, m2, m3], [p0, p1, p2, p3], [q0, q1, q2, q3]
def mag(a_, c_):
    ua, uc = arb(up(abs(a_))), arb(up(abs(c_))); return (ua * ua + uc * uc).sqrt()
def quantities(b, y, h):
    fact = [1, 1, 2, 6]; k, m, p, q = b.k, b.m, b.p, b.q
    yh = [y[n] + A * h[n] for n in range(5)]
    u = arb_series([yh[n] / fact[n] for n in range(4)], prec=4); up_ = arb_series([yh[n + 1] / fact[n] for n in range(4)], prec=4)
    kc, mc, pc, qc = profile_coeffs(k, m, p, q); K = arb_series(kc, prec=4); Q = arb_series(qc, prec=4)
    ag = 1 - K * A * u; assert ag.coeffs()[0].lower() > 0
    Z = up_ / ag; A1 = -(u - Z * Q); Th = arb_series([arb(0), kc[0], kc[1] / 2, kc[2] / 3], prec=4)
    Xs = Z * Th.cos(); Ys = Z * Th.sin()
    cA = A1.coeffs() + [arb(0)] * 4; cX = Xs.coeffs() + [arb(0)] * 4; cY = Ys.coeffs() + [arb(0)] * 4
    P = Pj_box(b.X, 0)
    d = {"psi1p": abs(k * q), "psi2p": abs(k), "psi3p": abs(y[1]), "h2p": abs(h[1]), "A1p_over_alpha": abs(cA[1]), "B1p_over_alpha": mag(cX[1], cY[1]),
         "A1ppp_over_alpha": 6 * abs(cA[3]), "B1ppp_over_alpha": 6 * mag(cX[3], cY[3]), "psi3ppp": abs(y[3]), "h2ppp": abs(h[3]),
         "c0p": abs(-arb(5) / 2 * m * p + k * k * q), "q_psi3pp": abs(q * y[2]), "q_psi3pp_alpha": abs(q * (y[2] + A * h[2])),
         "1-q_psi_p": abs(1 - q * P[0]), "sqrt(psi'^2+k^2)": (arb(up(abs(P[0]))) ** 2 + arb(kb_up(b)) ** 2).sqrt()}
    out = {kk: up(v) for kk, v in d.items()}
    if b.xa < 16.25 + 1e-9:
        X = b.X; phi1 = X * ((arb(1) + X * X) ** (-arb(1) / 3)) / 3
        assert (phi1 - P[0]).lower() > 0, 'W/rho not monotone on box'
        tt = arb(max(0.0, min(1.0, float((arb(b.xa) - arb(10.25)).lower()) / 6))).union(arb(max(0.0, min(1.0, float((arb(b.xb) - arb(10.25)).upper()) / 6))))
        om = 1 - tt; zp = 140 * (tt * tt * tt) * (om * om * om) / 6
        phi_xb = (arb(1) + arb(b.xb) * arb(b.xb)) ** (arb(2) / 3) / 4
        wr = (phi_xb - psi_value(math.nextafter(b.xb, UP))).exp()
        out["inner_ramp_q_zeta_p_W_over_rho"] = up(abs(b.q) * zp * wr); out["inner_ramp_W_over_rho"] = up(wr); out["inner_ramp_zeta_p"] = up(zp)
    return out
def kb_up(b): return up(b.k.upper() if False else arb(b.kb.upper()))
def run():
    t0 = time.time(); res = []
    bl = make_segments(BPS, [w * WS for w in (0.01, 0.01, 0.0125, 0.025, 0.05, 0.1, 0.1, 0.1)])
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb); y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
        for v in y + h: chk(v)
        res.append(dict(xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)), W=float(W_up(b.xa).mid()), q=quantities(b, y, h), z=zeta_sup(b.xa, b.xb)))
    print(f"[{time.time()-t0:.0f}s] {len(res)} boxes"); return res
def summarize(res):
    out = {"norms": {}, "sup_per_region": {}}
    names = [n for n in res[0]['q'] if not n.startswith('inner_ramp') and n not in ('1-q_psi_p', "sqrt(psi'^2+k^2)")]
    for nm in names:
        tz = arb(0); tn = arb(0); reg = {}
        for r in res:
            v = arb(r['q'][nm]); base = arb(r['len']) * arb(r['W']); z = arb(r['z']); tz += base * z * z * v * v; tn += base * v * v
            for i in region_of(r['xa'], r['xb']): reg[i] = reg.get(i, arb(0)) + base * z * z * v * v
        out["norms"][nm] = dict(norm_zeta=up(tz.sqrt()), norm_nozeta=up(tn.sqrt()), per_region_norm_zeta={reg_name(i): up(v.sqrt()) for i, v in sorted(reg.items())})
    for r in res:
        for ri in region_of(r['xa'], r['xb']):
            s = out["sup_per_region"].setdefault(reg_name(ri), {})
            for nm, v in r['q'].items(): s[nm] = max(s.get(nm, -1e300), v)
    out["sup_all"] = {nm: max(r['q'][nm] for r in res if nm in r['q']) for nm in res[0]['q']}
    return out
if __name__ == "__main__":
    res = run(); S = summarize(res); json.dump({"width_scale": WS, **S}, open(f'planB_mod_{WS:g}.json', 'w'), indent=1)
    for nm, v in S["norms"].items(): print(f"{nm:22s} zeta {v['norm_zeta']:.4g}  nozeta {v['norm_nozeta']:.4g}")
    print(json.dumps(S["sup_all"]))
