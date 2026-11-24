"""T2: residual norms in the intermediate weight W on [12, R] (adapted from psi3res3.py stage B and h2res.py stages Bt, B, T).
  r_psi = (L - lamt) psi3^ ,   r_h2 = (L - 2 lamt) h2^ + Q2(psi3^),   norms ||.||_{W,[12,R]}^2 = int_12^R W r^2.
sigma boxes [sa, sb] (x = sigma - ST, ST = tip enclosure; box lower x >= 12 rigorously; last box ends at R + 1e-9):
  psi3: ray + chunk bounds (psi3res3 stage B) with per-box step-1 tail bounds HS.profile_bounds(sa) on ALL boxes (also sigma >= 32);
        ray residual from the poly pieces (x <= 28), the explicit extension formula r_ray_ext (x >= 28), union across x = 28.
  h2:   sigma <= 20 : tight-tail second-order Taylor form (h2res stage_Bt) ;  sigma >= 20 : ray + chunk bounds + |Q2| (h2res stage_B).
 W-weight per box: W(x_lo) (W decreasing).   Usage: python3 int_res.py [width_scale] [weight: W|rho]   (rho = parity check with the rho-certificates)"""
import sys
from int_common import *
WS = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
WEIGHT = 'W'
lamF, lamH = F.LAM, Hf.LAM
q2, coef0, coef1, coef2, tail_state = HR.q2, HR.coef0, HR.coef1, HR.coef2, HR.tail_state
SB_END = HR.SB_END
def psi_res(b, X, xc):
    u = jets_any(F, X, xc)
    if b.xb <= float(F.XC.lower()): rr = HR.r_ray_poly(F, lamF, X, xc)
    elif b.xa >= float(F.XC.upper()): rr = HR.r_ray_ext(F, lamF, X)
    else: rr = HR.r_ray_poly(F, lamF, X, xc).union(HR.r_ray_ext(F, lamF, X))
    Bk, Bm, Bp, Bq = b.kb.upper(), b.mb.upper(), b.pb.upper(), b.qs.upper()
    rs = dST / 4 * abs(u[1])
    dl = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
    parts = dict(ray=rr, st=rs, prof=dl)
    return chk(rr + rs + dl), parts
def h2_res_B(b, X, xc):
    u = jets_any(Hf, X, xc)
    if b.xb <= float(Hf.XC.lower()): rr = HR.r_ray_poly(Hf, lamH, X, xc)
    elif b.xa >= float(Hf.XC.upper()): rr = HR.r_ray_ext(Hf, lamH, X)
    else: rr = HR.r_ray_poly(Hf, lamH, X, xc).union(HR.r_ray_ext(Hf, lamH, X))
    y = jets_any(F, X, xc)
    Bk, Bm, Bp, Bq = b.kb.upper(), b.mb.upper(), b.pb.upper(), b.qs.upper()
    q2v = abs(q2('Q2', y, arb(0, Bp), X + STIP, arb(0, Bk), arb(0, Bm)))
    rs = dST / 4 * abs(u[1])
    dl = Bk * Bk * abs(u[2]) + (5 * Bk * Bm + Bq / 4) * abs(u[1]) + (3 * Bm * Bm + Bk * Bp) * abs(u[0])
    parts = dict(ray=rr, st=rs, prof=dl, Q2=q2v)
    return chk(rr + rs + dl + q2v), parts
def h2_res_Bt(b):
    sa, sb = arb(b.sa), arb(b.sb); sc = (sa + sb) / 2
    p, q, k, m = tail_state(sc)
    Xc = sc - ST; Hh = (sa - sc).union(sb - sc); Xb = Xc + Hh
    P_, K_, M_ = b.p, b.k, b.m
    Q_ = b.sig + arb(0, b.qs.upper())
    a2, a1, a0 = coef0(p, q, k, m); a2p, a1p, a0p = coef1(p, q, k, m)
    b2, b1, b0 = coef0(P_, Q_, K_, M_); b2p, b1p, b0p = coef1(P_, Q_, K_, M_); b2pp, b1pp, b0pp = coef2(P_, Q_, K_, M_)
    best = None; parts = dict(R0=arb(0), R1h=arb(0), R2h2=arb(0))
    for Pf, Ph in zip(F.PIECES, Hf.PIECES):
        if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
        h = HR.jets(Hf, Ph, Xc); Hb = HR.jets(Hf, Ph, Xb); y = HR.jets(F, Pf, Xc); Yb = HR.jets(F, Pf, Xb)
        R0 = -h[4] - a2 * h[2] - a1 * h[1] + a0 * h[0] + q2('Q2', y[:5], p, q, k, m)
        R1 = (-h[5] - a2p * h[2] - a2 * h[3] - a1p * h[1] - a1 * h[2] + a0p * h[0] + a0 * h[1] + q2('DQ2', y[:6], p, q, k, m))
        R2 = (-Hb[6] - b2pp * Hb[2] - 2 * b2p * Hb[3] - b2 * Hb[4] - b1pp * Hb[1] - 2 * b1p * Hb[2] - b1 * Hb[3]
              + b0pp * Hb[0] + 2 * b0p * Hb[1] + b0 * Hb[2] + q2('D2Q2', Yb[:7], P_, Q_, K_, M_))
        r = abs(R0 + R1 * Hh + R2 * Hh * Hh / 2); chk(r)
        best = r if best is None else arb(max(best.upper(), r.upper()))
        hm = abs(Hh).upper()
        for key, v in (('R0', abs(R0)), ('R1h', abs(R1) * hm), ('R2h2', abs(R2) * hm * hm / 2)):
            if v.upper() > parts[key]: parts[key] = arb(v.upper())
    assert best is not None
    return arb(best.upper()), parts, float(Xb.lower())
def boxlist():
    bps = [('x', 10.25), ('x', 16.25), ('x', 20.0), ('s', float(SB_END)), ('x', 28.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)]
    wd = [0.005 * WS, 0.005 * WS, 0.005 * WS, 0.01 * WS, 0.01 * WS, 0.01 * WS, 0.02 * WS]
    return make_segments(bps, wd)
def run(stride=1, which=('psi', 'h2')):
    t0 = time.time(); bl = boxlist(); res = []
    for n, (sa, sb) in enumerate(bl):
        if n % stride: continue
        b = make_box(sa, sb); xc = b.xc
        rec = dict(sa=sa, sb=sb, xa=b.xa, xb=b.xb, len=up(arb(b.xb) - arb(b.xa)))
        xlo = b.xa
        if 'psi' in which:
            rp, pp_ = psi_res(b, b.X, xc); rec['psi'] = up(rp); rec['psi_parts'] = {k_: up(v) for k_, v in pp_.items()}
        if 'h2' in which:
            if sb <= SB_END:
                rh, ph, xl2 = h2_res_Bt(b); xlo = min(xlo, xl2); rec['h2'] = up(rh); rec['h2_parts'] = {k_: up(v) for k_, v in ph.items()}; rec['h2_stage'] = 'Bt'
            else:
                rh, ph = h2_res_B(b, b.X, xc); rec['h2'] = up(rh); rec['h2_parts'] = {k_: up(v) for k_, v in ph.items()}; rec['h2_stage'] = 'B'
        rec['W'] = {pk: float(W_up(pk, xlo).mid()) for pk in PXS}; rec['xlo_used'] = xlo
        res.append(rec)
        if n % 400 == 0: print(n, len(bl), f"x={b.xa:.3f}", {k_: rec.get(k_) for k_ in ('psi', 'h2')}, f"{time.time()-t0:.0f}s", flush=True)
    return res
SUBREG = [(10.25, 16.25), (16.25, 20.0), (20.0, 28.0), (28.0, 30.25), (30.25, 40.0), (40.0, 60.0)]
def summarize(res, stride=1):
    out = {}
    for pk in PXS:
      for which in ('psi', 'h2'):
        for Rr in (40.0, 60.0):
            tot = arb(0); reg = {}; wmax = (0.0, 0.0); pmax = {}
            for r in res:
                if which not in r or r['xb'] > Rr + 1e-6: continue
                Wv = r['W'][pk]
                c = arb(r['len']) * arb(Wv) * arb(r[which]) ** 2 * stride
                tot += c
                key = [i for i, (a, b) in enumerate(SUBREG) if r['xa'] >= a - 1e-6 and r['xa'] < b - 1e-6][0]
                reg[key] = reg.get(key, arb(0)) + c
                w = r[which] * Wv ** .5
                if w > wmax[0]: wmax = (w, r['xa'])
                for kk, v in r[which + '_parts'].items():
                    pm = pmax.setdefault(f"[{SUBREG[key][0]:g},{SUBREG[key][1]:g}]", {})
                    pm[kk] = max(pm.get(kk, 0.0), v * Wv ** .5)
            out[f"P={pk}_{which}_R={Rr:g}"] = dict(norm=up(tot.sqrt()), norm_sq=up(tot), max_weighted_pointwise=wmax[0], at=wmax[1],
                per_region_norm={f"[{SUBREG[i][0]:g},{SUBREG[i][1]:g}]": up(v.sqrt()) for i, v in sorted(reg.items())},
                max_weighted_parts_by_region=pmax)
    return out
if __name__ == "__main__":
    stride = 1
    res = run(stride); S = summarize(res, stride)
    tag = f"{WEIGHT}_{WS:g}" + ("" if stride == 1 else f"_stride{stride}")
    json.dump({"width_scale": WS, "weight": WEIGHT, "stride": stride, "boxes": len(res), **S}, open(f'v2_T2_{tag}.json', 'w'), indent=1)
    json.dump(res, open(f'v2_T2_boxes_{tag}.json', 'w'))
    for k_, v in S.items():
        print(k_, json.dumps({a: b for a, b in v.items() if a != 'per_region_norm_sq'}, indent=1))
