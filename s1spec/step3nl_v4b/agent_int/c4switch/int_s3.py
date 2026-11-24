"""S3: ||chi (N0(alpha psi3~ + alpha^2 h2~) - alpha^2 Q2(psi3~))||_rho <= C3 |alpha|^3, |alpha| <= AMAX (n3_forcing.py with the C^4 jets).
Series at the exact point alpha = 0 (c3, c4) + interval c5 remainder; stage C (CAPD boxes) and stage B (sigma 8..32, chunk bounds).
Usage: python3 int_s3.py AMAX [stride]"""
import sys, json, time
AM = sys.argv[1] if len(sys.argv) > 1 else '1e-3'
sys.argv = ['x', AM, '1']
sys.path.insert(0, '/home/claude/s1spec/step3nl')
from int_common import *
import n3_forcing as N3
from c4fun import PSI, H2
def stage_C(path='/home/claude/s1_tight/prof/prof_tight.txt'):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0); nb = 0
    for line in open(path):
        if line.startswith('END'): break
        v = list(map(float, line.split())); ta, tb = v[0], v[1]
        hl = [HR.I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        xlo = max(float(Xb.lower()), 0.0); xhi = float(Xb.upper()); p, q, k, m = hl
        best = arb(0)
        for e, P in enumerate(PSI.pieces):
            if Xb.upper() < P.a.lower() or Xb.lower() > P.b.upper(): continue
            y = PSI.piece_jets(e, Xb, 5); h = H2.piece_jets(e, Xb, 5)
            g = N3.g3_bound(y, h, p, q, k, m)
            chk(g); best = arb(max(best.upper(), g.upper()))
        rho = HR.rho_of(xlo); tot += arb(xhi - xlo) * rho * best * best; nb += 1
        w = float(best.upper()) * float(rho.sqrt().upper())
        if w > worst[0]: worst = (w, xlo)
    print(f"C: {nb} boxes {float(tot.upper()):.4e} worst {worst} [{time.time()-t0:.0f}s]", flush=True); return tot, worst
def stage_B(Nn=50):
    t0 = time.time(); tot = arb(0); worst = (0.0, 0.0)
    for ci in range(8, 32):
        for i in range(Nn):
            sa = arb(ci) + arb(i) / Nn; sb = arb(ci) + arb(i + 1) / Nn
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            p, k, m = arb(0, pb.upper()), arb(0, kb.upper()), arb(0, mb.upper()); q = sa.union(sb) + arb(0, qs.upper())
            xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2; X = arb(xc, (xb - xa) / 2)
            y = PSI.jets_any(X, xc, 5); h = H2.jets_any(X, xc, 5)
            g = chk(N3.g3_bound(y, h, p, q, k, m)); rho = HR.rho_of(xa)
            tot += arb(xb - xa) * rho * g * g
            w = float(g.upper()) * float(rho.sqrt().upper())
            if w > worst[0]: worst = (w, xa)
    print(f"B: {float(tot.upper()):.4e} worst {worst} [{time.time()-t0:.0f}s]", flush=True); return tot, worst
if __name__ == "__main__":
    N3.check_Q2(); tC, wC = stage_C(); tB, wB = stage_B(); C3 = (tC + tB).sqrt()
    out = {"alpha_max": float(AM), "C3": float(C3.upper()), "C2part": float(tC.upper()), "B2part": float(tB.upper()), "worstC": wC, "worstB": wB}
    json.dump(out, open(f'c4switch_S3_{AM}.json', 'w'), indent=1); print(out)
