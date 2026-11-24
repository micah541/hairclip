"""T1: third-order forcing in the intermediate weight W on [12, R] (adapted from n3_forcing.py).
 g(alpha) = N0(alpha psi3^ + alpha^2 h2^), g(0)=g'(0)=0, g''(0)/2 = Q2(psi3^);
 |g(alpha) - alpha^2 Q2| <= |alpha|^3 (|c3(0)| + amax |c4(0)| + amax^2 sup_A |c5|)  pointwise (c_j = g^(j)/j!, series at the EXACT point 0
 for c3, c4; interval point A = [-amax, amax] only for c5).  Norm: C3W(R)^2 = sum_boxes len * W(x_lo) * bound^2  (W decreasing).
Also sup |c3(0)| = sup |g'''/6| (pointwise coefficient).  Usage: python3 int_n3.py [AMAX] [width_scale]"""
import sys, pickle
from int_common import *
from flint import arb_series
AMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
WS = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
PREC = 6
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
class _Pr2(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"({e.p}/{e.q})"
    def _print_Integer(self, e): return f"({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
CF = {mon: eval("lambda Ji, al, k, m, p, q, u0, u1: " + _Pr2().doprint(cf)) for mon, cf in GS.items()}
def N0(u, p, q, k, m):
    al = 1 - k * u[0]
    J2 = al * al + u[1] * u[1]
    Ji = J2.rsqrt() if hasattr(J2, 'rsqrt') else 1 / J2.sqrt()
    Fc = 0
    for (i, j, l), f in CF.items():
        term = f(Ji, al, k, m, p, q, u[0], u[1])
        for _ in range(i): term = term * u[2]
        for _ in range(j): term = term * u[3]
        for _ in range(l): term = term * u[4]
        Fc = Fc + term
    drift = (p + u[0] - u[1] * q / al) / 4
    lin = -u[4] - k * k * u[2] - (5 * k * m + q / 4) * u[1] + (arb(1) / 4 - 3 * m * m - k * p) * u[0]
    return Fc + drift - lin
A = arb(0, AMAX)
def g3_bound(y, h, p, q, k, m):
    z = arb(0)
    u0s = [arb_series([z, y[j], h[j]], prec=PREC) for j in range(5)]
    c0 = N0(u0s, p, q, k, m).coeffs() + [arb(0)] * PREC
    uAs = [arb_series([A * y[j] + A * A * h[j], y[j] + 2 * A * h[j], h[j]], prec=PREC) for j in range(5)]
    cA = N0(uAs, p, q, k, m).coeffs() + [arb(0)] * PREC
    c3a = abs(c0[3]); tot = c3a + AMAX * abs(c0[4]) + AMAX * AMAX * abs(cA[5])
    for v in (c3a, tot): chk(v)
    return tot, c3a, c0
def check_Q2():
    y = [arb(0.9), arb(0.3), arb(1.2), arb(-0.4), arb(0.7)]; h = [arb(0.1)] * 5
    p, q, k, m = arb(-2.4), arb(0.3), arb(1.5), arb(0.05)
    us = [arb_series([arb(0), y[j], h[j]], prec=PREC) for j in range(5)]
    c = N0(us, p, q, k, m).coeffs(); q2 = HR.q2('Q2', y, p, q, k, m)
    print("sanity: c2 =", float(c[2].mid()), " Q2(y) =", float(q2.mid()), " c0,c1 =", float(c[0].mid()), float(c[1].mid()))
    assert abs(float((c[2] - q2).mid())) < 1e-9 and abs(float(c[0].mid())) < 1e-12 and abs(float(c[1].mid())) < 1e-12
def boxlist():
    bps = [('x', 12.0), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)]
    wd = [0.02 * WS, 0.025 * WS, 0.025 * WS, 0.05 * WS]
    return make_segments(bps, wd), [sig_of_x(v) for v in (12.0, 20.0, 30.25, 40.0, 60.0)]
def run():
    t0 = time.time(); bl, edges = boxlist(); res = []
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = jets_any(F, b.X, b.xc); h = jets_any(Hf, b.X, b.xc)
        for v in y + h: chk(v)
        tot, c3a, c0 = g3_bound(y, h, b.p, b.q, b.k, b.m)
        W = W_up(b.xa)
        res.append(dict(sa=sa, sb=sb, xa=b.xa, xb=b.xb, best=up(tot), c3=up(c3a), W=float(W.mid()), len=up(arb(b.xb) - arb(b.xa))))
        if n % 100 == 0: print(n, len(bl), b.xa, up(tot), f"{time.time()-t0:.0f}s", flush=True)
    return res
def summarize(res, tag=""):
    out = {}
    for Rr in (40.0, 60.0):
        tot = arb(0); reg = {}; sup3 = {}; wmax = (0.0, 0.0); sc3 = (0.0, 0.0)
        for r in res:
            if r['xb'] > Rr + 1e-6: continue
            c = arb(r['len']) * arb(r['W']) * arb(r['best']) ** 2
            tot += c
            key = [i for i, (a, b) in enumerate(REG) if r['xa'] >= a - 1e-6 and r['xa'] < b - 1e-6][0]
            reg[key] = reg.get(key, arb(0)) + c
            sup3[key] = max(sup3.get(key, 0.0), r['c3'])
            w = r['best'] * r['W'] ** .5
            if w > wmax[0]: wmax = (w, r['xa'])
            if r['c3'] > sc3[0]: sc3 = (r['c3'], r['xa'])
        out[f"R={Rr:g}"] = dict(C3W=up(tot.sqrt()), sup_c3=sc3[0], sup_c3_at=sc3[1], max_weighted_bound=wmax[0], at=wmax[1],
                               per_region_sq={f"[{REG[i][0]},{REG[i][1]}]": up(v) for i, v in reg.items()},
                               per_region_sup_c3={f"[{REG[i][0]},{REG[i][1]}]": v for i, v in sup3.items()})
    return out
if __name__ == "__main__":
    check_Q2()
    res = run(); S = summarize(res)
    json.dump({"alpha_max": AMAX, "width_scale": WS, "boxes": len(res), **S}, open(f'int_T1_{AMAX:g}_{WS:g}.json', 'w'), indent=1)
    json.dump(res, open(f'int_T1_boxes_{AMAX:g}_{WS:g}.json', 'w'))
    for k_, v in S.items(): print(k_, json.dumps(v, indent=1))
