"""Rigorous 8x8 Gram matrix of V = (L*f1, L*f2, gA, gB, D f1, D f2) and F = (f1, f2) in half-line L^2(rho), using
second-order Taylor models per box (centre-point enclosure + hull enclosure of the second Taylor coefficient):
   int_box H = len*H(c) + H'(c) int(x-c) + (1/2 H''(hull)) int (x-c)^2,    H = rho u_i u_j .
Regions:  x in [0, tau_tip]  (sigma <= 8):  tight profile boxes /home/claude/s1_tight/prof/prof_tight.txt  (centre + hull of (p,q,k,m,X5),
          X5 = theta - theta(sigma=8); theta(x) = -pi/2 + X5(tau_tip - x) - X5(tau_tip));
          x >= tau_tip (sigma in [8,32]):  boxes of outA5.txt ('MB' rows: chunk profile bounds, C,S hulls, as in lyap_consts2.py).
Everything is expressed through truncated Taylor series in s = x - c (python-flint arb_series), with the ODE
p' = -kq, q' = 1 + kp, k' = m, m' = p/4, c' = -k s, s' = k c (c = cos theta, s = sin theta; psi2 = -s, cos(theta) = c)."""
import sys, json, time
from math import comb
from flint import arb, arb_series, ctx
sys.argv_saved = sys.argv; ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod')
import mod_gram as MG
sys.argv = sys.argv_saved
HR, F = MG.HR, MG.F
ctx.prec = 100
N = 7                                   # series length for state / psi derivatives (need coefficient 2 of 4th derivatives of f)
A_, hull, centered, dSm_pts, _pw, cS = MG.A, MG.hull, MG.centered, MG.dSm_pts, MG._pw, MG.cS
NAMES = ["Lf1", "Lf2", "gA", "gB", "g0", "f1", "f2"]
def S_(x, n=N): return arb_series([x], prec=n)
def tr(sr, n):
    if isinstance(sr, arb): sr = arb_series([sr], prec=n)
    c = sr.coeffs(); c = (c + [arb(0)] * n)[:n]; return arb_series(c, prec=n)
def state_series(p0, q0, k0, m0, c0, s0, n=N):
    P, Q, K, M, C, S = [S_(v, n) for v in (p0, q0, k0, m0, c0, s0)]
    one = S_(arb(1), n)
    for _ in range(n + 1):
        Pn = S_(p0, n) + tr((-(K * Q)).integral(), n); Qn = S_(q0, n) + tr((one + K * P).integral(), n)
        Kn = S_(k0, n) + tr(M.integral(), n); Mn = S_(m0, n) + tr((P * arb(1) / 4).integral(), n)
        Cn = S_(c0, n) + tr((-(K * S)).integral(), n); Sn = S_(s0, n) + tr((K * C).integral(), n)
        P, Q, K, M, C, S = Pn, Qn, Kn, Mn, Cn, Sn
    return P, Q, K, M, C, S
def weight_series(Xs, n=N):
    one = S_(arb(1), n); X2 = Xs * Xs; L = (one + X2).log()
    E = (L * (arb(-1) / 3)).exp(); E2 = E * E; E4 = E2 * E2; E7 = E4 * E2 * E; E10 = E7 * E2 * E
    r = lambda a, b: arb(a) / arb(b)
    P1 = Xs * E * r(1, 3)
    P2 = E * r(1, 3) - X2 * E4 * r(2, 9)
    P3 = Xs * E4 * r(-2, 3) + X2 * Xs * E7 * r(16, 27)
    P4 = E4 * r(-2, 3) + X2 * E7 * r(32, 9) - X2 * X2 * E10 * r(224, 81)
    rho = ((-((L * (arb(2) / 3)).exp())) * (arb(1) / 4)).exp()
    return [P1, P2, P3, P4], rho
def dser(sr, j, n):
    for _ in range(j): sr = sr.derivative()
    return tr(sr, n)
def mul3(*a):
    r = a[0]
    for x in a[1:]: r = r * x
    return r
def functions(st, chi, Xs, yser, n=N, nout=3):
    """series (length nout) of the 8 functions.  st = state series (p,q,k,m,c,s); chi: series of chi (length n); Xs: series of x;
    yser: (y, y') series of psi3^ (length >= nout)."""
    p_, q_, k_, m_, c_, s_ = st
    Pw, rho = weight_series(Xs, n)
    a = {j: cS[j](Pw[0], Pw[1], Pw[2], Pw[3], k_, m_, q_, p_) for j in range(5)}
    out = {}
    for key, psi in (("1", p_), ("2", -s_)):
        g = tr(chi * psi, n)
        dg = [dser(g, j, nout) for j in range(5)]
        out["Lf" + key] = tr(sum((tr(a[j], nout) * dg[j] for j in range(1, 5)), tr(a[0], nout) * dg[0]), nout)
        out["f" + key] = tr(g, nout)
    y, y1 = yser
    c3 = tr(chi, nout)
    out["gA"] = tr(c3 * (tr(y, nout) - tr(q_, nout) * tr(y1, nout)), nout)
    out["gB"] = tr(c3 * tr(y1, nout) * tr(c_, nout), nout)
    out["g0"] = tr(c3 * tr(y, nout), nout)
    return out, tr(rho, nout)
def psi_series(Xb, n=3):
    """series of psi3^ and psi3^' about the (hull or point) ball Xb: coefficients y^(j)/j! (enclosures valid at every point of Xb)."""
    xc = float(Xb.mid()); J = HR.jets_any(F, Xb, xc)
    fac = [1, 1, 2, 6, 24]
    y = arb_series([J[j] / fac[j] for j in range(n)], prec=n)
    y1 = arb_series([J[j + 1] / fac[j] for j in range(n)], prec=n)
    return y, y1
def clampX(Xv):
    lo = max(float(Xv.lower()), 0.0); hi = float(Xv.upper())
    return arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def chi_series(uu, n=N):
    """chi(sigma) = S(u), u = (32 - sigma)/3 as series in s = x - c;  uu: ball for u over the region (point or hull)"""
    cs = []
    fact = 1
    for j in range(n):
        if j > 0: fact *= j
        cs.append(centered(dSm_pts, uu, j) * _pw(arb(-1) / 3, j) / fact)
    return arb_series(cs, prec=n)
NF = len(NAMES)
def accumulate(G, uc, uh, Rc, Rh, L, I1, I2, extra=None):
    """G[i][j] += L*H_c0 + I1*H_c1 + I2*h2  with  H = rho u_i u_j"""
    for i in range(NF):
        ci = uc[NAMES[i]].coeffs() + [arb(0)] * 3; hi_ = uh[NAMES[i]]
        for j in range(i, NF):
            cj = uc[NAMES[j]].coeffs() + [arb(0)] * 3
            Hc = Rc * uc[NAMES[i]] * uc[NAMES[j]]; hc = Hc.coeffs() + [arb(0)] * 3
            Hh = (Rh * uh[NAMES[i]] * uh[NAMES[j]]).coeffs() + [arb(0)] * 3
            t = L * hc[0] + I1 * hc[1] + I2 * Hh[2]
            G[i][j] += t
def tip_region(G, GROUP=3, stride_limit=None, verbose=True):
    t0 = time.time()
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    rows = [list(map(float, l.split())) for l in lines]
    TT = HR.TT
    def hullv(grp, j): return HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp))
    last = rows[-1]; kl = HR.I(last[13 + 4], last[14 + 4])
    X5c_last = HR.I(last[3 + 8], last[4 + 8])
    X5_tip = X5c_last - kl * (TT - arb(last[2]))
    nb = 0; extra_tot = arb(0)
    assert len(rows) % GROUP == 0 and GROUP % 2 == 1
    for g0 in range(0, len(rows), GROUP):
        grp = rows[g0:g0 + GROUP]; mid = grp[GROUP // 2]
        ta, tb, tc = arb(grp[0][0]), arb(grp[-1][1]), arb(mid[2])
        Xb = (TT - tb).union(TT - ta); Xc = TT - tc
        if Xb.upper() < 0: continue
        Xb = clampX(Xb)
        Xcl = Xc if Xc.lower() >= 0 else clampX(Xc)
        cs = [HR.I(mid[3 + 2 * j], mid[4 + 2 * j]) for j in range(5)]
        hs = [hullv(grp, j) for j in range(5)]
        def cs_(v5):
            th = -arb.pi() / 2 + (v5 - X5_tip); return th.cos(), th.sin()
        cc, sc = cs_(cs[4]); ch_, sh = cs_(hs[4])
        stc = state_series(cs[0], cs[1], cs[2], cs[3], cc, sc, 6)
        sth = state_series(hs[0], hs[1], hs[2], hs[3], ch_, sh, N)
        one6 = S_(arb(1), 6); oneN = S_(arb(1), N)
        Xs_c = arb_series([Xcl, arb(1)], prec=6); Xs_h = arb_series([Xb, arb(1)], prec=N)
        yc = psi_series(Xcl); yh = psi_series(Xb)
        uc, Rc = functions(stc, one6, Xs_c, yc, 6, 3); uh, Rh = functions(sth, oneN, Xs_h, yh, N, 3)
        # integration in tau: x - c = -(tau - tau_c); s in [tc - tb, tc - ta]
        lo_, hi_ = tc - tb, tc - ta
        L = tb - ta; I1 = (hi_ * hi_ - lo_ * lo_) / 2; I2 = (hi_ * hi_ * hi_ - lo_ * lo_ * lo_) / 3
        accumulate(G, uc, uh, Rc, Rh, L, I1, I2)
        # last box sticks out below x = 0 by at most (tb - TT.lower()): add +-(that)*sup|H|
        if TT.lower() < tb:
            over = tb - TT.lower()
            for i in range(NF):
                for j in range(i, NF):
                    Hb = (Rh * uh[NAMES[i]] * uh[NAMES[j]]).coeffs()
                    mag = sum((abs(c_) * arb(max(abs(float(lo_.lower())), abs(float(hi_.upper())))) ** 0 for c_ in Hb[:1]), arb(0)) + abs(Hb[1] if len(Hb) > 1 else arb(0)) * arb(abs(float(lo_.lower())) + 1e-3) + abs(Hb[2] if len(Hb) > 2 else arb(0)) * arb((abs(float(lo_.lower())) + 1e-3) ** 2)
                    G[i][j] += arb(0, float((over.upper()) * mag.upper()))
        nb += 1
    if verbose: print(f"tip region: {nb} boxes [{time.time()-t0:.0f}s]", flush=True)
    return nb
def far_region(G, NSUB=1, verbose=True):
    t0 = time.time()
    CtB, StB, Ct2h, St2h = MG.frame_data(); CH = MG.chunks()
    nb = 0
    for line in open('/home/claude/s1spec/step3/outA5.txt'):
        if line[:2] != "MB": continue
        v = list(map(float, line.split()[1:])); ta0, tb0 = v[0], v[1]
        h = [A_(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        st, C_, S_v = h[0], h[1], h[2]
        c = min(int(ta0 + 1e-9), len(CH) - 1)
        kv, mv, pv_, dq = [A_(lo, hi) for lo, hi in CH[c][2]]
        s0h = -(CtB * C_ + StB * S_v); c0h = StB * C_ - CtB * S_v           # sin(theta), cos(theta) hulls over the row
        for isub in range(NSUB):
            ta = arb(ta0) + (arb(tb0) - arb(ta0)) * isub / NSUB; tb = arb(ta0) + (arb(tb0) - arb(ta0)) * (isub + 1) / NSUB
            tc = (ta + tb) / 2
            sig_b = 32 - (ta.union(tb)); sig_c = 32 - tc
            Xb = clampX(sig_b - st); Xc = clampX(sig_c - st)
            q_c = sig_c + dq; q_b = sig_b + dq
            if float(tb.upper()) <= 3.0 + 1e-12:
                chic = chi_series(arb((32 - sig_c.mid()) / 3 * 0 + (tc / 3).mid(), 0) if False else (tc / 3), N)
                chih = chi_series((ta / 3).union(tb / 3), N)
            else:
                chic = chih = arb_series([arb(1)], prec=N)
            stc = state_series(pv_, q_c, kv, mv, c0h, s0h, N); sth = state_series(pv_, q_b, kv, mv, c0h, s0h, N)
            Xs_c = arb_series([Xc, arb(1)], prec=N); Xs_h = arb_series([Xb, arb(1)], prec=N)
            uc, Rc = functions(stc, chic, Xs_c, psi_series(Xc), N, 3); uh, Rh = functions(sth, chih, Xs_h, psi_series(Xb), N, 3)
            # x = sigma - sigma_tip: s = x - x_c = sigma - sigma_c = -(tau - tau_c); tau in [ta, tb]
            lo_, hi_ = tc - tb, tc - ta
            L = tb - ta; I1 = (hi_ * hi_ - lo_ * lo_) / 2; I2 = (hi_ * hi_ * hi_ - lo_ * lo_ * lo_) / 3
            accumulate(G, uc, uh, Rc, Rh, L, I1, I2)
            nb += 1
    if verbose: print(f"far region: {nb} boxes [{time.time()-t0:.0f}s]", flush=True)
    return nb
if __name__ == "__main__":
    GROUP = int(ARGS[0]) if ARGS else 3; NSUB = int(ARGS[1]) if len(ARGS) > 1 else 1; tag = ARGS[2] if len(ARGS) > 2 else ""
    G = [[arb(0) for _ in range(NF)] for _ in range(NF)]
    t0 = time.time()
    far_region(G, NSUB); tip_region(G, GROUP)
    for i in range(NF):
        for j in range(i):  G[i][j] = G[j][i]
        for j in range(NF): assert G[i][j].is_finite(), (i, j)
    json.dump({"names": NAMES, "G": [[[float(G[i][j].lower()), float(G[i][j].upper())] for j in range(NF)] for i in range(NF)],
               "group": GROUP, "nsub": NSUB, "runtime_s": time.time() - t0}, open(f'gram7{tag}.json', 'w'), indent=1)
    for i in range(NF): print(f"{NAMES[i]:4s}", "  ".join(f"{float(G[i][j].mid()):+.6g}(r{float(G[i][j].rad()):.0e})" for j in range(NF)))
    print(f"[{time.time()-t0:.0f}s]")
