"""V5: W-norms of the modulation / beta source terms of the v-equation (v = zeta w'), V6: outer-ramp bookkeeping.
Pointwise upper bounds per box (profile = HS.profile_bounds(sigma_lo); |cos theta|, |sin theta| <= 1 used for psi2' and B1):
  psi1' = -k q ;  psi2' = -k cos(theta) (bounded by |k|) ;  psi3^' ;  h2^' ;
  (A1(u0))'/alpha,  A1(u) = -(u - u' q/al_g), al_g = 1 - k u :  A1' = -(u' - u'' q/al_g - u' (q/al_g)'),  (q/al_g)' = q'/al_g - q al_g'/al_g^2,
       q' = 1 + k p, al_g' = -(m u + k u'); u = A yh0, yh_j = psi3^(j) + A h2^(j), A = [-amax, amax] (u0^(j) = A yh_j, so the ratio is taken exactly);
  (B1(u0))'/alpha,  B1(u) = u' cos(theta)/al_g, theta' = k:  |B1'| <= (|u''| + |u'||k|)/al_g + |u'||al_g'|/al_g^2 ;
  c0' = -(5/2) m p + k^2 q.
zeta = S((x-10.25)/6) (inner ramp, 0 below 10.25), outer ramp S((R-x)/3) on [R-3, R]; box sup of zeta used (S increasing): norms
  ||zeta f||_W^2 = sum len * W(x_lo) * sup(zeta^2) * sup f^2.   Also the zeta-free norm over [10.25, R] is reported.
Usage: python3 int_v56.py [width_scale]"""
import sys
from int_common import *
AMAX = 1e-3
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
def Spoly_f(t):
    t = min(max(t, 0.0), 1.0); a = arb(t)
    return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb, R):
    zi = up(Spoly_f((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = up(Spoly_f((R - xa) / 3)) if xa > R - 3 else 1.0
    return min(1.0, zi) * min(1.0, zo)
A = arb(0, AMAX)
def quantities(b, y, h):
    k, m, p, q = b.k, b.m, b.p, b.q
    yh = [y[j] + A * h[j] for j in range(3)]
    ag = 1 - k * A * yh[0]; assert ag.lower() > 0
    agp = -(m * A * yh[0] + k * A * yh[1])
    qa = q / ag; qap = (1 + k * p) / ag - q * agp / (ag * ag)
    A1 = -(yh[1] - yh[2] * qa - yh[1] * qap)
    B1 = (abs(yh[2]) + abs(yh[1]) * abs(k)) / ag + abs(yh[1]) * abs(agp) / (ag * ag)
    c0p = -arb(5) / 2 * m * p + k * k * q
    d = {"psi1p": abs(k * q), "psi2p": abs(k), "psi3p": abs(y[1]), "h2p": abs(h[1]), "A1p_over_alpha": abs(A1),
         "B1p_over_alpha": B1, "c0p": abs(c0p)}
    return {kk: up(v) for kk, v in d.items()}
def run():
    t0 = time.time()
    bl = make_segments([('x', 10.25), ('x', 16.25), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)],
                       [0.01 * WS, 0.01 * WS, 0.0125 * WS, 0.0125 * WS, 0.025 * WS])
    res = []
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = jets_any(F, b.X, b.xc); h = jets_any(Hf, b.X, b.xc)
        for v in y + h: chk(v)
        r = dict(xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)), q=quantities(b, y, h),
                 W={pk: float(W_up(pk, b.xa).mid()) for pk in PXS}, z={str(R): zeta_sup(b.xa, b.xb, R) for R in (40.0, 60.0)})
        res.append(r)
    print(f"[{time.time()-t0:.0f}s] {len(res)} boxes"); return res
def summarize(res):
    out = {}
    names = list(res[0]['q'].keys())
    for pk in PXS:
        for Rr in (40.0, 60.0):
            o = {}
            for nm in names:
                tz = arb(0); tn = arb(0); supz = 0.0
                for r in res:
                    if r['xb'] > Rr + 1e-6: continue
                    v = arb(r['q'][nm]); Wv = arb(r['W'][pk]); z = arb(r['z'][str(Rr)])
                    tz += arb(r['len']) * Wv * z * z * v * v; tn += arb(r['len']) * Wv * v * v
                o[nm] = dict(norm_zeta=up(tz.sqrt()), norm_nozeta=up(tn.sqrt()))
            out[f"P={pk}_R={Rr:g}"] = o
    sup = {}
    regs = REG
    for r in res:
        for ri in region_of(r['xa'], r['xb']):
            s = sup.setdefault(f"[{regs[ri][0]:g},{regs[ri][1]:g}]", {})
            for nm, v in r['q'].items(): s[nm] = max(s.get(nm, 0.0), v)
    out["sup_per_region"] = sup
    return out
def v6():
    out = {}
    for pk in PXS:
        for Rr in (40.0, 60.0):
            o = {"W(R-3)": up(arb(0) + W_up(pk, Rr - 3)), "W(R)": up(W_up(pk, Rr)), "W(R)_lower": dn(W_lo(pk, Rr)), "W(R-3)_lower": dn(W_lo(pk, Rr - 3))}
            n = 3000; up_i = arb(0); lo_i = arb(0); h_ = 3.0 / n
            for i in range(n):
                a, c = Rr - 3 + i * h_, Rr - 3 + (i + 1) * h_
                up_i += arb(h_) * W_up(pk, a); lo_i += arb(h_) * W_lo(pk, c)
            o["int_W_outer_upper"] = up(up_i); o["int_W_outer_lower"] = dn(lo_i)
            o["psi(R)"] = [float(psi_value(pk, Rr).mid()), float(psi_value(pk, Rr).rad())]
            out[f"P={pk}_R={Rr:g}"] = o
    return out
if __name__ == "__main__":
    res = run(); S = summarize(res); V6 = v6()
    json.dump({"width_scale": WS, "V5": S, "V6": V6}, open(f'v2_V56_{WS:g}.json', 'w'), indent=1)
    for k_ in ("P=2_R=60", "P=1.5_R=60"): print(k_, json.dumps({a: b['norm_zeta'] for a, b in S[k_].items()}))
    print(json.dumps(S["sup_per_region"], indent=0)); print(json.dumps(V6, indent=0))
