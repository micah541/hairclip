"""F1b: level-0 and DERIVATIVE-level residual bounds for the C^4 functions psi3~, h2~ (c4fun.py) on [10.25, 60], weight W (P=2).
  r_psi = (L - lamt) psi3~,  r_h2 = (L - 2lamt) h2~ + Q2(psi3~);   level 1 = d/dx of these.
sigma <= 20 (Stage Bt, BOTH functions now): tight-tail profile at the box centre sigma_c (HR.tail_state), exact total derivatives R, R', R'', R'''
  (resid_sym.py):   level 0:  R(xc) + R'(xc) h + R''(box) h^2/2 ;   level 1:  R'(xc) + R''(xc) h + R'''(box) h^2/2   (h = x - xc, hull of box profile for the box terms).
sigma >= 20 (Stage B): ray residual (Taylor form on the polynomial pieces; closed form (exact telescoping + WKB term) on the extension) + tip term
  dST/4 |u'| (level 1: |u''|) + profile terms (chunk/tail bounds of k, m, p, q-sigma; level 1 uses the differentiated coefficient bounds) + |Q2| (level 1: |DQ2|) for h2.
Norms: level 0 plain W-norm, level 1 with zeta (box sup).  Usage: python3 int_res3.py [width_scale]"""
import sys
from int_common import *
import c4fun, resid_sym as RS
from c4fun import PSI, H2, ext_jets, _zjets
import psi3fun as F3, h2fun as H3
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
PK = '2'
LAMS = {'psi': F3.LAM, 'h2': H3.LAM}
SB_END = HR.SB_END
q2 = HR.q2
assert all(abs(float((a.a - b.a).mid())) < 1e-15 for a, b in zip(PSI.pieces, H2.pieces))
def Spoly_f(t):
    t = min(max(t, 0.0), 1.0); a = arb(t); return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb, R):
    zi = up(Spoly_f((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = up(Spoly_f((R - xa) / 3)) if xa > R - 3 else 1.0
    return min(1.0, zi) * min(1.0, zo)
def bt(b):
    sa, sb = arb(b.sa), arb(b.sb); sc = (sa + sb) / 2
    p, q, k, m = HR.tail_state(sc)
    Xc = sc - ST; Hh = (sa - sc).union(sb - sc); Xb = Xc + Hh
    P_, K_, M_ = b.p, b.k, b.m; Q_ = b.sig + arb(0, b.qs.upper())
    out = {w: [None, None] for w in ('psi', 'h2')}
    for e, P in enumerate(PSI.pieces):
        if Xb.upper() < P.a.lower() or Xb.lower() > P.b.upper(): continue
        yc = PSI.piece_jets(e, Xc, 8); hc = H2.piece_jets(e, Xc, 8); Yb = PSI.piece_jets(e, Xb, 8); Hb = H2.piece_jets(e, Xb, 8)
        for w in ('psi', 'h2'):
            lam = LAMS[w]
            R0 = RS.call(w, 0, yc, hc, p, q, k, m, lam); R1 = RS.call(w, 1, yc, hc, p, q, k, m, lam); R2c = RS.call(w, 2, yc, hc, p, q, k, m, lam)
            R2b = RS.call(w, 2, Yb, Hb, P_, Q_, K_, M_, lam); R3b = RS.call(w, 3, Yb, Hb, P_, Q_, K_, M_, lam)
            l0 = chk(abs(R0 + R1 * Hh + R2b * Hh * Hh / 2)); l1 = chk(abs(R1 + R2c * Hh + R3b * Hh * Hh / 2))
            for i, v in enumerate((l0, l1)):
                out[w][i] = arb(v.upper()) if out[w][i] is None else arb(max(out[w][i].upper(), v.upper()))
    assert out['psi'][0] is not None
    return {w: (up(v[0]), up(v[1])) for w, v in out.items()}
def rpoly(C, lam, X, xc):
    best0 = best1 = None; c = arb(xc); hh = X - c
    for e, P in enumerate(C.pieces):
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = C.raw(e, c, 8); jX = C.raw(e, X, 8); Hc = C.corr(e, X, 6)
        sc = c + C.M.STIP; sX = X + C.M.STIP; kap = arb(1) / 4 - lam
        r0 = -jc[4] - sc / 4 * jc[1] + kap * jc[0]
        r1 = -jc[5] - jc[1] / 4 - sc / 4 * jc[2] + kap * jc[1]
        r2 = -jX[6] - jX[2] / 2 - sX / 4 * jX[3] + kap * jX[2]
        rH0 = -Hc[4] - sX / 4 * Hc[1] + kap * Hc[0]
        a1c = -jc[6] - jc[2] / 2 - sc / 4 * jc[3] + kap * jc[2]
        a2X = -jX[7] - arb(3) / 4 * jX[3] - sX / 4 * jX[4] + kap * jX[3]
        rH1 = -Hc[5] - Hc[1] / 4 - sX / 4 * Hc[2] + kap * Hc[1]
        v0 = abs(r0 + r1 * hh + r2 * hh * hh / 2) + abs(rH0); v1 = abs(r1 + a1c * hh + a2X * hh * hh / 2) + abs(rH1)
        best0 = v0 if best0 is None else best0.union(v0); best1 = v1 if best1 is None else best1.union(v1)
    return best0, best1
def rext(M, lam, X):
    s = X + M.STIP; z = _zjets(s)
    E = (-3 * (s / 4) ** (arb(4) / 3) + 3 * ((M.XC + M.STIP) / 4) ** (arb(4) / 3)).exp()
    kap = arb(1) / 4 - lam
    rho = -(z[3] + 4 * z[0] * z[2] + 3 * z[1] * z[1] + 6 * z[0] * z[0] * z[1]) + kap
    drho = -(z[4] + 10 * z[1] * z[2] + 4 * z[0] * z[3] + 12 * z[0] * z[1] * z[1] + 6 * z[0] * z[0] * z[2])
    g0 = abs(M.C1 * 4 * M.AJ[4] * s ** (M.A - 16)); g1 = abs(M.C1 * 4 * M.AJ[4] * (M.A - 16) * s ** (M.A - 17))
    return g0 + abs(M.C2 * E * rho), g1 + abs(M.C2 * E * (z[0] * rho + drho))
def ray(C, lam, b):
    M = C.M
    if b.xb <= float(M.XC.lower()): return rpoly(C, lam, b.X, b.xc)
    if b.xa >= float(M.XC.upper()): return rext(M, lam, b.X)
    a = rpoly(C, lam, b.X, b.xc); c = rext(M, lam, b.X); return a[0].union(c[0]), a[1].union(c[1])
def stageB(b):
    X, xc = b.X, b.xc; y = PSI.jets_any(X, xc, 7); h = H2.jets_any(X, xc, 7)
    for v in y + h: chk(v)
    Bk, Bm, Bp, Bq = b.kb.upper(), b.mb.upper(), b.pb.upper(), b.qs.upper(); qa = arb(up(abs(b.q)))
    out = {}
    for w, C, u in (('psi', PSI, y), ('h2', H2, h)):
        lam = LAMS[w]; rr0, rr1 = ray(C, lam, b)
        rs0 = dST / 4 * abs(u[1]); rs1 = dST / 4 * abs(u[2])
        dl0 = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
        dl1 = (2 * Bk * Bm * abs(u[2]) + Bk * Bk * abs(u[3]) + (5 * Bm * Bm + arb(3) / 2 * Bk * Bp) * abs(u[1]) + (5 * Bk * Bm + Bq / 4) * abs(u[2])
               + (arb(5) / 2 * Bm * Bp + Bk * Bk * qa) * abs(u[0]) + (3 * Bm * Bm + Bk * Bp) * abs(u[1]))
        t0 = rr0 + rs0 + dl0; t1 = rr1 + rs1 + dl1
        if w == 'h2':
            t0 = t0 + abs(q2('Q2', y, arb(0, Bp), X + STIP, arb(0, Bk), arb(0, Bm)))
            t1 = t1 + abs(q2('DQ2', y[:6], arb(0, Bp), X + STIP, arb(0, Bk), arb(0, Bm)))
        out[w] = (up(chk(t0)), up(chk(t1)))
    return out
def boxlist():
    bps = [('x', 10.25), ('x', 16.25), ('x', 20.0), ('s', float(SB_END)), ('x', 28.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)]
    return make_segments(bps, [0.005 * WS, 0.005 * WS, 0.005 * WS, 0.01 * WS, 0.01 * WS, 0.01 * WS, 0.02 * WS])
SUBREG = [(10.25, 16.25), (16.25, 20.0), (20.0, 28.0), (28.0, 30.25), (30.25, 40.0), (40.0, 60.0)]
def run():
    t0 = time.time(); bl = boxlist(); res = []
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        r = bt(b) if sb <= SB_END else stageB(b)
        xlo = b.xa
        res.append(dict(xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)), W=float(W_up(PK, xlo).mid()), stage='Bt' if sb <= SB_END else 'B',
                        z={str(R): zeta_sup(b.xa, b.xb, R) for R in (40.0, 60.0)}, r=r))
        if n % 500 == 0: print(n, len(bl), f"x={b.xa:.3f}", r, f"{time.time()-t0:.0f}s", flush=True)
    return res
def summarize(res):
    out = {}
    for w in ('psi', 'h2'):
        for lev in (0, 1):
            for Rr in (40.0, 60.0):
                tot = arb(0); reg = {}; wmax = (0.0, 0.0)
                for r in res:
                    if r['xb'] > Rr + 1e-6: continue
                    z = arb(r['z'][str(Rr)]) if lev == 1 else arb(1)
                    c = arb(r['len']) * arb(r['W']) * z * z * arb(r['r'][w][lev]) ** 2; tot += c
                    key = [i for i, (a, b_) in enumerate(SUBREG) if r['xa'] >= a - 1e-6 and r['xa'] < b_ - 1e-6][0]
                    reg[key] = reg.get(key, arb(0)) + c
                    wv = r['r'][w][lev] * (r['W'] * float(z.mid()) ** 2) ** .5
                    if wv > wmax[0]: wmax = (wv, r['xa'])
                out[f"{w}_level{lev}_R={Rr:g}"] = dict(norm=up(tot.sqrt()), max_weighted_pointwise=wmax[0], at=wmax[1],
                    per_region_norm={f"[{SUBREG[i][0]:g},{SUBREG[i][1]:g}]": up(v.sqrt()) for i, v in sorted(reg.items())})
    return out
if __name__ == "__main__":
    res = run(); S = summarize(res)
    json.dump({"width_scale": WS, **S}, open(f'v3_res_{WS:g}.json', 'w'), indent=1)
    json.dump(res, open(f'v3_res_boxes_{WS:g}.json', 'w'))
    for k_, v in S.items(): print(k_, f"{v['norm']:.4e}", {a: f"{b:.3e}" for a, b in v['per_region_norm'].items()}, v['max_weighted_pointwise'], v['at'])
