"""Box generator for the IBP modulation lemma.  Boxes = the hull boxes of outA5.txt (2400 'MB' rows, sigma in [8,32], and 1025 'MC' rows near the tip),
exactly those of lemma/lem_integrals.py, but with: C^4 jets (0..4) of psi3^, h2^ (agent_int/v3/c4fun.py), the profile (k,m,p,q), the frame (C,S,Ct,St),
the jets 0..4 of chi, of the base functions b1 = p, b2 = psi2 = Ct C + St S, b3 = F = e^phi w* (ws0..ws3 from the hull, F'''' from L*F = lamt F),
tip bumps of f3, and psi_hat_{j,r} = rho^-1 (rho chi f_j)^(r), r = 0..4  (f1 = chi p, f2 = chi psi2, f3 = chi F - a eta1 - b eta3)."""
import sys, json, time
from flint import arb
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod'); sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_int/v3')
import mod_gram as MG
sys.argv = sys.argv_saved
import sympy as sp
from math import comb
import gen_hs3 as G3
import c4fun as C4
ab = lambda z: arb(abs(z).upper())
A_, hull, centered, dSm_pts, _pw, cS, dp, dC, dS, weight, chunks = MG.A, MG.hull, MG.centered, MG.dSm_pts, MG._pw, MG.cS, MG.dp, MG.dC, MG.dS, MG.weight, MG.chunks
LEMX = json.load(open('/home/claude/s1spec/step3nl/agent_mod/lemma/lyap_extra.json'))
A_MAG, B_MAG = LEMX['a_'] * (1 + 1e-12), LEMX['b_'] * (1 + 1e-12)
LAMT = arb(G3.LAMT)
# numeric evaluators of Em, Bp (polynomials in P1..P4)
P1s, P2s, P3s, P4s = MG.P1, MG.P2, MG.P3, MG.P4
Emf = [MG.fn(MG.Em[n], (P1s, P2s, P3s, P4s)) for n in range(5)]
Bpf = [MG.fn(MG.Bp[n], (P1s, P2s, P3s, P4s)) for n in range(5)]
# bumps
sv = sp.Symbol('s'); zeta = (1 - sv**2) ** 4; LBs = sp.Rational(G3.LB).limit_denominator(10)
e1s = LBs * sv * zeta; e3s = LBs**3 * sv**3 * zeta / 6
de1_pts = [MG.fn(sp.diff(e1s, sv, j), (sv,)) for j in range(10)]; de3_pts = [MG.fn(sp.diff(e3s, sv, j), (sv,)) for j in range(12)]
def bump_ders(Xv, cc, which):
    s_ = Xv / (cc * arb(G3.LB)); f = de1_pts if which == 1 else de3_pts
    return [centered(f, s_, j) / _pw(cc * arb(G3.LB), j) for j in range(5)]
def boxes():
    CtB, StB, Ct2h, St2h = MG.frame_data(); T0 = arb(G3.T0); CH = chunks()
    for line in open('/home/claude/s1spec/step3/outA5.txt'):
        tag = line[:2]
        if tag not in ("MB", "MC"): continue
        v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
        h = [A_(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        bump = False; cc = None
        if tag == "MB":
            st, C_, S_ = h[0], h[1], h[2]; Ct, St = CtB, StB; ws = h[5:9]
            c = min(int(ta + 1e-9), len(CH) - 1)
            kv, mv, pv_, dq = [A_(lo, hi) for lo, hi in CH[c][2]]
            sig = A_(32 - tb, 32 - ta); Xv = sig - st; qv = sig + dq; dX = arb(tb) - arb(ta)
            if tb <= 3.0 + 1e-12:
                uu = A_(ta / 3, tb / 3); chi = [centered(dSm_pts, uu, j) * _pw(arb(-1) / 3, j) for j in range(5)]; strip = True
            else: chi = [arb(1)] + [arb(0)] * 4; strip = False
        else:
            pv_, qv, kv, mv, st, C_, S_ = h[0], h[1], h[2], h[3], h[4], h[7], h[8]; Ct, St = Ct2h, St2h; ws = h[9:13]
            cc = (8 - st) / T0
            Xv = hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta))); dX = cc * (arb(tb) - arb(ta)); chi = [arb(1)] + [arb(0)] * 4; strip = False
            bump = ta >= float(T0) - G3.LB - 1e-12
        lo = max(float(Xv.lower()), 0.0); hi = float(Xv.upper())
        Xv = arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)
        P, rho = weight(Xv)
        pvs = (kv, mv, qv, pv_, C_, S_)
        b1 = [dp[n](*pvs) for n in range(5)]
        b2 = [Ct * dC[n](*pvs) + St * dS[n](*pvs) for n in range(5)]
        # F^(n) = sum_j C(n,j) ws_j B_{n-j}, n <= 3 ;  F'''' from  L* F = lamt F  (aS[4] = -1)
        Bv = [Bpf[n](*P) for n in range(5)]
        F = [sum((comb(n, j) * ws[j] * Bv[n - j] for j in range(n + 1)), arb(0)) for n in range(4)]
        aSv = {n: cS[n](P[0], P[1], P[2], P[3], kv, mv, qv, pv_) for n in range(5)}
        F4 = sum((aSv[n] * F[n] for n in range(4)), arb(0)) - LAMT * F[0]
        F.append(F4)
        # g_j = chi f_j = chi^2 b_j (j = 1,2,3) minus chi*bumps for j = 3
        chi2 = [sum((comb(n, s) * chi[s] * chi[n - s] for s in range(n + 1)), arb(0)) for n in range(5)]
        g = {}
        for j, bj in ((1, b1), (2, b2), (3, F)): g[j] = [sum((comb(n, s) * chi2[s] * bj[n - s] for s in range(n + 1))) if False else sum((comb(n, s) * chi2[s] * bj[n - s] for s in range(n + 1)), arb(0)) for n in range(5)]
        if bump:
            ed1 = bump_ders(Xv, cc, 1); ed3 = bump_ders(Xv, cc, 3)
            for n in range(5): g[3][n] = g[3][n] - arb(0, A_MAG) * ab(ed1[n]) - arb(0, B_MAG) * ab(ed3[n])      # chi = 1 on the bump support (sign of a, b irrelevant)
            # (ball: a*eta1 in [-|a||eta1|, |a||eta1|]; the arb(0, r) * |eta| form encloses it)
        psih = {j: [sum((comb(r, s) * Emf[r - s](*P) * g[j][s] for s in range(r + 1)), arb(0)) for r in range(5)] for j in (1, 2, 3)}
        xc = float(Xv.mid())
        psi_j = C4.PSI.jets_any(Xv, xc, 5); h2_j = C4.H2.jets_any(Xv, xc, 5)
        yield dict(tag=tag, ta=ta, tb=tb, strip=strip, bump=bump, X=Xv, dX=dX, rho=rho, w=dX * rho, k=kv, m=mv, p=pv_, q=qv, chi=chi, ps=psi_j, hh=h2_j, psih=psih, g=g, F=F, b1=b1, b2=b2)
if __name__ == "__main__":
    t0 = time.time(); n = 0
    for b in boxes():
        n += 1
        if n in (1, 500, 1500, 2300, 2420, 3000, 3400):
            print(n, b['tag'], float(b['X'].mid()), "psih", [float(b['psih'][j][0].mid()) for j in (1, 2, 3)], "rad", [float(b['psih'][j][0].rad()) for j in (1, 2, 3)], "strip", b['strip'], "bump", b['bump'])
    print(n, "boxes", time.time() - t0)
