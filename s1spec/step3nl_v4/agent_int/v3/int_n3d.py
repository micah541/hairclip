"""F1a: DERIVATIVE-level third-order forcing for the C^4 functions psi3~, h2~:
   g(alpha; x) = N0(alpha psi3~ + alpha^2 h2~) - alpha^2 Q2(psi3~),  c_j(x) = alpha^j-coefficient;  c0 = c1 = 0, c2 = Q2 (value) and (c2)' = DQ2 (checked).
   |g'(alpha; x)| <= |alpha|^3 ( |c3'| + amax |c4'| + amax^2 sup_{xi in A} |c5'| ),   A = [-amax, amax],
 with c_j' = d/dx c_j(x), computed by forward-mode differentiation along x: dual numbers over python-flint arb_series in alpha, with
   u_j' = u_{j+1} (jets of the C^4 functions to order 5), profile ' from the ODE (k' = m, m' = p/4, p' = -kq, q' = 1+kp).
 Series at the EXACT point alpha = 0 for c3', c4' (interval point only for c5'); the value parts give the level-0 constant for the C^4 functions as a by-product.
 Norms: ||zeta g'||_W (zeta = box sup, inner ramp S((x-10.25)/6), outer ramp S((R-x)/3)), W = (1+x^2)^-2 for x <= 30.25 (P=2).
Usage: python3 int_n3d.py [AMAX] [width_scale]"""
import sys, pickle
from int_common import *
from flint import arb_series
from c4fun import PSI, H2
AMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
WS = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
PREC = 6; PK = '2'
class Dual:
    __slots__ = ('a', 'b')
    def __init__(s, a, b): s.a, s.b = a, b
    @staticmethod
    def lift(o): return o if isinstance(o, Dual) else Dual(o, 0)
    def __add__(s, o): o = Dual.lift(o); return Dual(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = Dual.lift(o); return Dual(s.a - o.a, s.b - o.b)
    def __rsub__(s, o): o = Dual.lift(o); return Dual(o.a - s.a, o.b - s.b)
    def __neg__(s): return Dual(-s.a, -s.b)
    def __mul__(s, o):
        o = Dual.lift(o)
        if isinstance(o.b, int) and o.b == 0: return Dual(s.a * o.a, s.b * o.a)
        return Dual(s.a * o.a, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = Dual.lift(o)
        if isinstance(o.b, int) and o.b == 0: return Dual(s.a / o.a, s.b / o.a)
        return Dual(s.a / o.a, (s.b * o.a - s.a * o.b) / (o.a * o.a))
    def __rtruediv__(s, o): return Dual.lift(o) / s
    def rsqrt(s):
        r = s.a.rsqrt(); return Dual(r, -s.b * r * r * r / 2)
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
class _Pr2(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
CF = {mon: eval("lambda Ji, al, k, m, p, q, u0, u1: " + _Pr2().doprint(cf), {"arb": arb}) for mon, cf in GS.items()}
def N0(u, p, q, k, m):
    al = 1 - k * u[0]; J2 = al * al + u[1] * u[1]; Ji = J2.rsqrt(); Fc = 0
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
def prof(p, q, k, m):
    return (Dual(p, -k * q), Dual(q, 1 + k * p), Dual(k, m), Dual(m, p / 4))
def dseries(y, h, p, q, k, m, exact):
    if exact:
        z = arb(0); S = lambda j: arb_series([z, y[j], h[j]], prec=PREC)
    else:
        S = lambda j: arb_series([A * y[j] + A * A * h[j], y[j] + 2 * A * h[j], h[j]], prec=PREC)
    us = [Dual(S(j), S(j + 1)) for j in range(5)]
    P, Q, K, M = prof(p, q, k, m)
    return N0(us, P, Q, K, M)
def gd_bound(y, h, p, q, k, m):
    r0 = dseries(y, h, p, q, k, m, True); rA = dseries(y, h, p, q, k, m, False)
    cv = r0.a.coeffs() + [arb(0)] * PREC; cd = r0.b.coeffs() + [arb(0)] * PREC; cA = rA.b.coeffs() + [arb(0)] * PREC
    d3 = abs(cd[3]); tot = d3 + AMAX * abs(cd[4]) + AMAX * AMAX * abs(cA[5])
    v3 = abs(cv[3]); totv = v3 + AMAX * abs(cv[4]) + AMAX * AMAX * abs(rA.a.coeffs()[5] if len(rA.a.coeffs()) > 5 else arb(0))
    return chk(tot), chk(d3), chk(totv), cv, cd
def check():
    y = [arb(0.9), arb(0.3), arb(1.2), arb(-0.4), arb(0.7), arb(0.2)]; h = [arb(0.1)] * 6
    p, q, k, m = arb(-2.4), arb(0.3), arb(1.5), arb(0.05)
    r = dseries(y, h, p, q, k, m, True); cv = r.a.coeffs(); cd = r.b.coeffs()
    import h2res as HR2
    q2 = HR2.q2('Q2', y[:5], p, q, k, m); dq2 = HR2.q2('DQ2', y[:6], p, q, k, m)
    print("sanity: c2 =", float(cv[2].mid()), "Q2 =", float(q2.mid()), "; (c2)' =", float(cd[2].mid()), "DQ2 =", float(dq2.mid()), "; c0',c1' =", float(cd[0].mid()), float(cd[1].mid()))
    assert abs(float((cv[2] - q2).mid())) < 1e-9 and abs(float((cd[2] - dq2).mid())) < 1e-9 and abs(float(cd[0].mid())) < 1e-12 and abs(float(cd[1].mid())) < 1e-12
def Spoly_f(t):
    t = min(max(t, 0.0), 1.0); a = arb(t); return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb, R):
    zi = up(Spoly_f((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = up(Spoly_f((R - xa) / 3)) if xa > R - 3 else 1.0
    return min(1.0, zi) * min(1.0, zo)
def run():
    t0 = time.time()
    bl = make_segments([('x', 10.25), ('x', 16.25), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)], [0.02 * WS, 0.02 * WS, 0.025 * WS, 0.025 * WS, 0.05 * WS])
    res = []
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = PSI.jets_any(b.X, b.xc, 6); h = H2.jets_any(b.X, b.xc, 6)
        for v in y + h: chk(v)
        tot, d3, totv, cv, cd = gd_bound(y, h, b.p, b.q, b.k, b.m)
        res.append(dict(xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)), W=float(W_up(PK, b.xa).mid()), best=up(tot), d3=up(d3), bestv=up(totv),
                        z={str(R): zeta_sup(b.xa, b.xb, R) for R in (40.0, 60.0)}))
        if n % 200 == 0: print(n, len(bl), b.xa, up(tot), f"{time.time()-t0:.0f}s", flush=True)
    return res
def summarize(res):
    out = {}
    for Rr in (40.0, 60.0):
        tz = arb(0); tn = arb(0); tv = arb(0); reg = {}; sup = {}
        for r in res:
            if r['xb'] > Rr + 1e-6: continue
            z = arb(r['z'][str(Rr)]); base = arb(r['len']) * arb(r['W'])
            c = base * z * z * arb(r['best']) ** 2; tz += c; tn += base * arb(r['best']) ** 2; tv += base * arb(r['bestv']) ** 2
            key = [i for i, (a, b_) in enumerate(REG) if r['xa'] >= a - 1e-6 and r['xa'] < b_ - 1e-6][0]
            reg[key] = reg.get(key, arb(0)) + c; sup[key] = max(sup.get(key, 0.0), r['d3'])
        out[f"R={Rr:g}"] = dict(C3W_deriv_zeta=up(tz.sqrt()), C3W_deriv_nozeta=up(tn.sqrt()), C3W_level0_C4=up(tv.sqrt()),
                              sup_c3prime=max(sup.values()),
                              per_region_norm={f"[{REG[i][0]},{REG[i][1]}]": up(v.sqrt()) for i, v in sorted(reg.items())},
                              per_region_sup_c3prime={f"[{REG[i][0]},{REG[i][1]}]": v for i, v in sorted(sup.items())})
    return out
if __name__ == "__main__":
    check(); res = run(); S = summarize(res)
    json.dump({"alpha_max": AMAX, "width_scale": WS, **S}, open(f'v3_n3d_{AMAX:g}_{WS:g}.json', 'w'), indent=1)
    for k_, v in S.items(): print(k_, json.dumps(v, indent=1))
