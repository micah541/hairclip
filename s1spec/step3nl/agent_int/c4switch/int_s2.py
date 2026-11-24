"""S2: rho-residual certificates for the C^4 functions psi3~, h2~ (half line, weight rho = exp(-phi), phi = (1+x^2)^(2/3)/4):
   ||(L - lamt) psi3~||_rho   (cf. psi3res3.py: 1.2772e-6 for psi3^)      ||(L - 2lamt) h2~ + Q2(psi3~)||_rho  (cf. h2res.py: 1.5719e-6 for h2^)
 Stages as in the originals, with the generic symbolic Taylor forms (resid_sym.py: exact R, R', R''):
   C  : x in [0, tau_tip] (51,246 CAPD boxes of prof_tight.txt): R(xc) + R'(xc) h + R''(box) h^2/2   (both functions; Q2 evaluated on psi3~ jets)
   Bt : sigma in [8, 20], tight-tail profile at the box centre, same Taylor form (psi3 too; the original psi3res3 used ray+chunk bounds there: ALSO reported as 'psi_rayB')
   B  : sigma in [20, 32] ray + chunk bounds (+ |Q2| for h2)
   T  : sigma >= 32: the extension is untouched by the C^4 corrections, so the ORIGINAL analytic tail bounds (psi3res3.stage_T, h2res.stage_T) apply verbatim.
 Usage: python3 int_s2.py [width_scale] [stageC_stride]    (stride > 1: diagnostic only)"""
import sys, json, time
WS_ARGV = sys.argv[1:]
sys.argv = ['x', WS_ARGV[0] if WS_ARGV else '1']
from int_common import *
import int_res3_base as B3
from c4fun import PSI, H2
import resid_sym as RS
WS = float(WS_ARGV[0]) if WS_ARGV else 1.0
STRIDE = int(WS_ARGV[1]) if len(WS_ARGV) > 1 else 1
LAMS = B3.LAMS
rho_of = HR.rho_of
def stage_C(path='/home/claude/s1_tight/prof/prof_tight.txt'):
    t0 = time.time(); tot = {'psi': arb(0), 'h2': arb(0)}; worst = {'psi': (0.0, 0.0), 'h2': (0.0, 0.0)}; nb = 0
    for li, line in enumerate(open(path)):
        if line.startswith('END'): break
        if li % STRIDE: continue
        v = list(map(float, line.split())); ta, tb, tc = v[0], v[1], v[2]
        c = [HR.I(v[3 + 2 * j], v[4 + 2 * j]) for j in range(4)]; hl = [HR.I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(4)]
        Xc = HR.TT - arb(tc); Hh = arb(tc - tb).union(arb(tc - ta)); Xb = Xc + Hh
        if Xb.upper() < 0: continue
        xlo = max(float(Xb.lower()), 0.0); xhi = float(Xb.upper())
        p, q, k, m = c; P_, Q_, K_, M_ = hl
        best = {'psi': None, 'h2': None}
        for e, P in enumerate(PSI.pieces):
            if Xb.upper() < P.a.lower() or Xb.lower() > P.b.upper(): continue
            yc = PSI.piece_jets(e, Xc, 7); hc = H2.piece_jets(e, Xc, 7); Yb = PSI.piece_jets(e, Xb, 7); Hb = H2.piece_jets(e, Xb, 7)
            for w in ('psi', 'h2'):
                lam = LAMS[w]
                R0 = RS.call(w, 0, yc, hc, p, q, k, m, lam); R1 = RS.call(w, 1, yc, hc, p, q, k, m, lam); R2 = RS.call(w, 2, Yb, Hb, P_, Q_, K_, M_, lam)
                r = chk(abs(R0 + R1 * Hh + R2 * Hh * Hh / 2))
                best[w] = arb(r.upper()) if best[w] is None else arb(max(best[w].upper(), r.upper()))
        rho = rho_of(xlo)
        for w in best:
            tot[w] += arb(xhi - xlo) * rho * best[w] * best[w]
            wv = float(best[w].upper()) * float(rho.sqrt().upper())
            if wv > worst[w][0]: worst[w] = (wv, xlo)
        nb += 1
    print(f"C: {nb} boxes: psi ||r||^2 <= {float(tot['psi'].upper()):.4e}, h2 ||r||^2 <= {float(tot['h2'].upper()):.4e}; worst weighted {worst} [{time.time()-t0:.0f}s]", flush=True)
    return tot, worst
def boxes_BtB():
    return make_segments([('s', 8.0), ('s', 20.0), ('s', 32.0)], [0.005 * WS, 0.01 * WS])
def stage_BtB():
    t0 = time.time(); tot = {k: {'Bt': arb(0), 'B': arb(0)} for k in ('psi', 'h2', 'psi_rayB')}; worst = {}
    for n, (sa, sb) in enumerate(boxes_BtB()):
        b = make_box(sa, sb); rho = rho_of(b.xa - 1e-12); ln = arb(b.xb) - arb(b.xa) + arb(2e-12)
        if sb <= 20.0:
            r = B3.bt(b); rb = B3.stageB(b); st = 'Bt'
            vals = {'psi': r['psi'][0], 'h2': r['h2'][0], 'psi_rayB': rb['psi'][0]}
        else:
            rb = B3.stageB(b); st = 'B'; vals = {'psi': rb['psi'][0], 'h2': rb['h2'][0], 'psi_rayB': rb['psi'][0]}
        for kk, v in vals.items():
            tot[kk][st] += ln * rho * arb(v) ** 2
            wv = v * float(rho.sqrt().upper())
            if wv > worst.get(kk, (0.0, 0.0))[0]: worst[kk] = (wv, b.xa)
    print(f"Bt/B done [{time.time()-t0:.0f}s]", {k: {s: float(v.upper()) for s, v in d.items()} for k, d in tot.items()}, worst, flush=True)
    return tot, worst
def stage_T():
    sys.argv = ['x', '1']
    import psi3res3 as P3
    return P3.stage_T(), HR.stage_T()
if __name__ == "__main__":
    print("tip interval ST =", ST)
    totC, worstC = stage_C()
    if STRIDE > 1: sys.exit(0)
    totB, worstB = stage_BtB(); Tpsi, Th2 = stage_T()
    res = {}
    for w, T in (('psi', Tpsi), ('h2', Th2)):
        tot = totC[w] + totB[w]['Bt'] + totB[w]['B'] + T
        res[w] = dict(C2=float(totC[w].upper()), Bt2=float(totB[w]['Bt'].upper()), B2=float(totB[w]['B'].upper()), T2=float(T.upper()), r_half=float(tot.upper()) ** .5,
                      worstC=worstC[w], worstBtB=worstB[w])
    tot = totC['psi'] + totB['psi_rayB']['Bt'] + totB['psi_rayB']['B'] + Tpsi
    res['psi_rayB_for_sigma_8_20'] = dict(r_half=float(tot.upper()) ** .5, Bt_as_rayB2=float(totB['psi_rayB']['Bt'].upper()))
    json.dump({"width_scale": WS, **res}, open(f'c4switch_S2_{WS:g}.json', 'w'), indent=1)
    print(json.dumps(res, indent=1))
