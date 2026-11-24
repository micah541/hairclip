"""Hull-box (outA5.txt, 0.01 boxes) upper bounds of the weighted integrals needed by the modulation lemma (crude is fine: they only enter
perturbation terms).  Reuses the box generator logic of ../mod_gram.py (copied with extra fields).  All integrals: sum len*rho(box)*|hull|,
absolute values taken as arb(abs(z).upper()).  Regions: core = sigma <= 29 (chi = 1), strip = sigma in [29,32]; x<=12 / x>12 split for norms."""
import sys, json, time
from flint import arb
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import mod_gram as MG
sys.argv = sys.argv_saved
HR = MG.HR; F = HR.F; Hf = HR.Hf
A_, hull, centered, dSm_pts, _pw, cS, dp, dC, dS, weight, leib, chunks = MG.A, MG.hull, MG.centered, MG.dSm_pts, MG._pw, MG.cS, MG.dp, MG.dC, MG.dS, MG.weight, MG.leib, MG.chunks
import gen_hs3 as G3
AMAX = 1e-3; ETA0 = 1e-2
ab = lambda z: arb(abs(z).upper())
def boxes():
    CtB, StB, Ct2h, St2h = MG.frame_data(); T0 = arb(G3.T0); CH = chunks()
    for line in open('outA5.txt'):
        tag = line[:2]
        if tag not in ("MB", "MC"): continue
        v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
        h = [A_(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        bump = False
        if tag == "MB":
            st, C_, S_ = h[0], h[1], h[2]; Ct, St = CtB, StB; ws0 = h[5]
            c = min(int(ta + 1e-9), len(CH) - 1)
            kv, mv, pv_, dq = [A_(lo, hi) for lo, hi in CH[c][2]]
            sig = A_(32 - tb, 32 - ta); Xv = sig - st; qv = sig + dq; dX = arb(tb) - arb(ta)
            if tb <= 3.0 + 1e-12:
                uu = A_(ta / 3, tb / 3); chi = [centered(dSm_pts, uu, j) * _pw(arb(-1) / 3, j) for j in range(5)]; strip = True
            else: chi = [arb(1)] + [arb(0)] * 4; strip = False
        else:
            pv_, qv, kv, mv, st, C_, S_ = h[0], h[1], h[2], h[3], h[4], h[7], h[8]; Ct, St = Ct2h, St2h; ws0 = h[9]
            cc = (8 - st) / T0
            Xv = MG.hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta))); dX = cc * (arb(tb) - arb(ta)); chi = [arb(1)] + [arb(0)] * 4; strip = False
            bump = ta >= float(T0) - G3.LB - 1e-12
        lo = max(float(Xv.lower()), 0.0); hi = float(Xv.upper())
        Xv = arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)
        P, rho = weight(Xv)
        pvs = (kv, mv, qv, pv_, C_, S_)
        ps1 = [dp[n](*pvs) for n in range(5)]; ps2 = [Ct * dC[n](*pvs) + St * dS[n](*pvs) for n in range(5)]
        xc = float(Xv.mid()); y = HR.jets_any(F, Xv, xc)[:2]; hh = HR.jets_any(Hf, Xv, xc)[:2]
        cth = Ct and (St * C_ - Ct * S_)
        yield dict(strip=strip, bump=bump, X=Xv, w=dX * rho, chi=chi, k=kv, m=mv, q=qv, p=pv_, f1=chi[0] * ps1[0], f2=chi[0] * ps2[0], F=chi[0] * ws0,
                   cth=cth, y=y, hh=hh, tag=tag, psi1=ps1[0], psi2=ps2[0], ws0=ws0)
def run():
    t0 = time.time()
    BUMP = 0.8324 * 9.0923e-8 + 0.4245 * 2.5707e-7          # |a_c eta1| + |b_c eta3| <= 1.9e-7 on the tip bump support (psi3ip.py constants, rounded up)
    BUMP = BUMP * 1.01
    # first pass: global alpha_g lower bound
    kmax = 0.0; kpsi = 0.0
    data = list(boxes())
    for b in data:
        kk = float(ab(b['k']).upper()); kmax = max(kmax, kk)
        kpsi = max(kpsi, kk * (AMAX * float(ab(b['y'][0]).upper()) + AMAX ** 2 * float(ab(b['hh'][0]).upper())))
    ag_min = 1 - kpsi - kmax * ETA0
    print(f"boxes {len(data)}; sup|k| = {kmax:.5f}; sup|k||u0| = {kpsi:.3e}; alpha_g >= {ag_min:.5f}  [{time.time()-t0:.0f}s]", flush=True)
    AGM = arb(ag_min)
    Z = lambda: arb(0)
    regs = ("all", "x<=12", "x>12", "core")
    nf2 = {j: {r: Z() for r in regs} for j in (1, 2, 3)}
    # E^alpha (abs integrals over all boxes): coefficient of |alpha| (h2 parts with alpha<=AMAX)
    e1 = {j: Z() for j in (1, 2, 3)}; e2 = {j: Z() for j in (1, 2, 3)}
    # E^w (core): squared L2 norms of q f_j, cos f_j, and with the k u0' factor
    q2n = {j: Z() for j in (1, 2, 3)}; c2n = {j: Z() for j in (1, 2, 3)}; q2k = {j: Z() for j in (1, 2, 3)}; c2k = {j: Z() for j in (1, 2, 3)}
    # strip L1 integrals
    IS = {j: Z() for j in (1, 2, 3)}; ISq = {j: Z() for j in (1, 2, 3)}; ISc = {j: Z() for j in (1, 2, 3)}
    ISqk = {j: Z() for j in (1, 2, 3)}; ISck = {j: Z() for j in (1, 2, 3)}
    kap = {j: [Z() for _ in range(4)] for j in (1, 2, 3)}
    SP = {'y_1': Z(), 'y_2': Z(), 'y_F': Z(), '1_F': Z(), '2_F': Z()}
    for b in data:
        W = arb(b['w'].upper() if True else 0)
        W = ab(b['w']); X = b['X']
        fj = {1: ab(b['f1']), 2: ab(b['f2']), 3: ab(b['F']) + (arb(BUMP) if b['bump'] else arb(0))}
        xhi = float(X.upper())
        k_, q_, c_ = ab(b['k']), ab(b['q']), ab(b['cth'])
        y0, y1, h0, h1 = ab(b['y'][0]), ab(b['y'][1]), ab(b['hh'][0]), ab(b['hh'][1])
        u0a = y0 + arb(AMAX) * h0; u1a = y1 + arb(AMAX) * h1                 # |u0| / |alpha|, |u0'| / |alpha|
        if b['strip']:
            Wy = W * ab(b['y'][0]); SP['y_1'] += Wy * ab(b['psi1']); SP['y_2'] += Wy * ab(b['psi2']); SP['y_F'] += Wy * ab(b['ws0'])
            SP['1_F'] += W * ab(b['psi1']) * ab(b['ws0']); SP['2_F'] += W * ab(b['psi2']) * ab(b['ws0'])
        for j in (1, 2, 3):
            f = fj[j]
            for r in regs:
                if r == "all" or (r == "x<=12" and xhi <= 12) or (r == "x>12" and xhi > 12) or (r == "core" and not b['strip']):
                    nf2[j][r] += W * f * f
            e1[j] += W * f * (u0a + q_ * u1a / AGM); e2[j] += W * f * c_ * u1a / AGM
            if not b['strip']:
                q2n[j] += W * (q_ * f / AGM) ** 2 if False else W * (q_ * f / AGM) * (q_ * f / AGM)
                c2n[j] += W * (c_ * f / AGM) * (c_ * f / AGM)
                tq = q_ * f * k_ * u1a / (AGM * AGM); tc = c_ * f * k_ * u1a / (AGM * AGM)
                q2k[j] += W * tq * tq; c2k[j] += W * tc * tc
            else:
                IS[j] += W * f; ISq[j] += W * f * q_ / AGM; ISc[j] += W * f * c_ / AGM
                ISqk[j] += W * f * q_ * k_ * u1a / (AGM * AGM); ISck[j] += W * f * c_ * k_ * u1a / (AGM * AGM)
                chi = b['chi']
                ch = [ab(chi[i]) for i in range(5)]
                k2 = k_ * k_
                cc = [ch[4] + k2 * ch[2] + (arb(5) * k_ * ab(b['m']) + q_ / 4) * ch[1], arb(4) * ch[3] + 2 * k2 * ch[1], arb(6) * ch[2], arb(4) * ch[1]]
                for i in range(4): kap[j][i] += W * f * cc[i]
    out = {"ag_min": ag_min, "sup_k": kmax, "AMAX": AMAX, "ETA0": ETA0, "BUMP": BUMP}
    sq = lambda z: float(z.sqrt().upper()) if z.upper() > 0 else 0.0
    up = lambda z: float(z.upper())
    for j in (1, 2, 3):
        out[f"nf_{j}"] = {r: sq(nf2[j][r]) for r in regs}
        out[f"e1_{j}"] = up(e1[j]); out[f"e2_{j}"] = up(e2[j])
        out[f"Pq1_{j}"] = sq(q2n[j]); out[f"Pc1_{j}"] = sq(c2n[j]); out[f"Pq0_{j}"] = sq(q2k[j]); out[f"Pc0_{j}"] = sq(c2k[j])
        out[f"IS_{j}"] = up(IS[j]); out[f"ISq_{j}"] = up(ISq[j]); out[f"ISc_{j}"] = up(ISc[j]); out[f"ISqk_{j}"] = up(ISqk[j]); out[f"ISck_{j}"] = up(ISck[j])
        out[f"kappa_{j}"] = [up(z) for z in kap[j]]
    out['SP'] = {k_: up(v_) for k_, v_ in SP.items()}
    for kk_, v in out.items():
        if isinstance(v, float): assert v == v
    json.dump(out, open('lem_integrals.json', 'w'), indent=1)
    for kk_, v in out.items(): print(kk_, v)
    print(f"[{time.time()-t0:.0f}s]")
if __name__ == "__main__": run()
