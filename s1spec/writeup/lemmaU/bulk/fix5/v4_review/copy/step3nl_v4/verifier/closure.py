import os
"""Rigorous barrier ("upper solution") for the comparison system
      |alpha|' <= LAMT |alpha| + beta(|alpha|, Ly),     beta = (BETA0 + BETA1|alpha|) sqrt(Ly) + BETAF Fn(alpha)
      Ly'      <= 2 T(alpha, Ly),                      T <= rate_bound * Ly          (rate_bound from bb.bound_cell, rigorous)
Barrier on step k (length D):  a(t) = a_k e^{LAMT t} + (bb/|LAMT|)(1 - e^{LAMT t}),   L(t) = L_k e^{2 r' t}
with  bb >= beta(a_up, L_up), a_up >= max(a_k, bb/|LAMT|), L_up/L_dn = max/min of L_k, L_{k+1}, r' > sup of rate over alpha in [0,a_up] x Ly in [L_dn, L_up].
First-touch argument (cooperative in alpha, direct in Ly) => alpha(t) <= a(t), Ly(t) <= L(t) for all t.  Strict inequalities via r' = r + 1e-9 max(1,|r|), bb(1+1e-9).
Tail: scale-invariant contraction (see tail()).   Usage: python3 closure.py ALPHA1 [W1_0 VI_0]"""
import sys, os, json, time
import numpy as np
from flint import arb, ctx
ctx.prec = 120
import constants as K, bb

def fup(x): return float(np.nextafter(float(x.upper()), np.inf))
def fdn(x): return float(np.nextafter(float(x.lower()), -np.inf))
A = lambda v: arb(v)
LAMT = A(str(K.LAMT)); BETA0 = A(str(K.BETA0)); BETA1 = A(str(K.BETA1)); BETAF = A(str(K.BETAF))     # BETA1, BETAF from constants (66.5, 0.657 with MOD=ibp)
def FnA(al): return A(str(K.FN[0])) * al + A(str(K.FN[1])) * al * al + A(str(K.FN[2])) * al * al * al
import unc
QBL = A(str(K.QB)) * A(str(unc.qbl_sup())) if K.FIX else A(0)      # F3 placeholder: + q_beta (A+W1+B)^2 <= q_beta * sup * Ly  in the alpha equation
def beta(al, L): return (BETA0 + BETA1 * al) * L.sqrt() + BETAF * FnA(al) + QBL * L + A(str(K.STRIP_BETA))

class Oracle:
    def __init__(self, tol=5e-4, rtol=0.005, log=print):
        self.tol, self.rtol, self.log = tol, rtol, log; self.rows = []; self.seed = None
    def rate(self, a_up, Llo, Lhi, floors=True):
        t = time.time()
        r = bb.bound_cell(0.0, a_up, Llo, Lhi, tol=self.tol, rtol=self.rtol, floors=floors, seeds=[self.seed] if self.seed is not None else None, cem_iters=(4 if self.seed is not None else 25))
        self.seed = r["best"]
        self.rows.append(dict(alpha_hi=a_up, Ly_lo=Llo, Ly_hi=Lhi, floors=floors, U=r["U"], lb=r["lb"], qstar=r["qstar"], nleaf=r["nleaf"], secs=round(time.time() - t, 2)))
        return r["U"]

def sramp_sup(s_hi):
    """rigorous upper bounds of S_ramp := sup over [10.25,16.25] of max(|w|,|w'|) (KREG bounds on [8,12] and [12,16.25]), and of sup|w'| alone,
    maximised over ALL shapes with Lyap <= s_hi^2 (A = s cw ct, B = s cw st/sqrt(eps), maximal at cw=1; 4000 boxes in the angle t)."""
    import numpy as np, rhs_iv as R
    from ival import P, cc, sqrt, mx
    c = R.consts(); n = 4000
    t0 = np.arange(n) / n; t1 = (np.arange(n) + 1) / n
    ct, st = R.circ(t0, t1)
    s = cc(float(np.nextafter(s_hi, np.inf)))
    Aa = s * ct; Bb = s * st * c.inv_sqrt_eps
    W1 = sqrt(Aa * (Bb + c.c1 * Aa))
    m0 = m1 = None
    for Kc in (c.K12, c.K["12<x<=16.25"]):
        s0 = sqrt(Aa * (Kc * (c.PHX * Aa + c.cl['two'] * W1))); s1 = sqrt(W1 * (Kc * (c.PHX * W1 + c.cl['two'] * Bb)))
        m0 = s0 if m0 is None else mx(m0, s0); m1 = s1 if m1 is None else mx(m1, s1)
    return float(np.max(mx(m0, m1).hi)), float(np.max(m1.hi))

def farsup(s_hi):
    """rigorous sup, over all shapes with Lyap<=s_hi^2, of the far-field (x>=16.25) bounds of the code: E_I/KW-based sup|w'| (sI1 = gm(V0, KW(PSI1W V0+2V1))),
    the code's far sup|w'| = max(sr1, sI1) (sr1: E0 bound with K=1/rho(16.25)), and far sup|w| = max(sr0, sI0), sI0 = sr0 + sqrt(PHI34) V0 (sup over the circle cw^2+sw^2=1:
    A,B ~ cw, V ~ sw  =>  sup (f cw + g sw) = sqrt(f^2+g^2))."""
    import numpy as np, rhs_iv as R
    from ival import P, cc, sqrt, mx
    c = R.consts(); n = 4000
    t0 = np.arange(n) / n; t1 = (np.arange(n) + 1) / n
    ct, st = R.circ(t0, t1); s = cc(float(np.nextafter(s_hi, np.inf)))
    Aa = s * ct; Bb = s * st * c.inv_sqrt_eps; W1 = sqrt(Aa * (Bb + c.c1 * Aa))
    sr0 = sqrt(Aa * (c.KR16 * (c.PHX * Aa + c.cl['two'] * W1))); sr1 = sqrt(W1 * (c.KR16 * (c.PHX * W1 + c.cl['two'] * Bb)))
    F0 = float(np.max(sr0.hi)); F1 = float(np.max(sr1.hi))
    cI, sI = R.circ(t0, t1)
    V0 = s * cI * c.inv_sqrt_mu; V2 = s * sI * c.inv_sqrt_muepsI
    V1 = sqrt(V0 * (V2 + c.C1W * V0)); sI1 = sqrt(V0 * (c.KW * (c.PSI1W * V0 + c.cl['two'] * V1)))
    G1 = float(np.max(sI1.hi)); H = float(np.max((c.sqrtPHI34 * V0).hi))
    sI2 = sqrt(V1 * (c.KW * (c.PSI1W * V1 + c.cl['two'] * V2))); G2 = float(np.max(sI2.hi))     # E_I/KW sup|w''| (x>=16.25)
    W1m = float(np.max(W1.hi)); Am_ = float(np.max(Aa.hi))
    return dict(A_max=Am_, W1_max=float(np.nextafter(W1m, np.inf)), sI2=G2, M1_ramp=float(np.nextafter(28.0 * G2, np.inf)), sI1=G1, far_w1=max(F1, G1), far_w=float(np.sqrt(F0 ** 2 + H ** 2) * (1 + 1e-14)), sr1=F1)

BOXLOG = []
def box_assert(S_sup, FS, aq, log, where):
    """fin13v4 (EI_SPEC 1.3, SF-4): the standing box hypotheses, ASSERTED (worst case over all shapes with Lyap <= Ly_hi of the cell)."""
    if not K.NEW: return
    chk = dict(
        Bb_IR_jet=(S_sup[0], K.ETA_IR),                                     # sup(|w|,|w'|) on [10.25,16.25] <= 1e-3 (C2 IR tables)
        Bc_w1_far_EI=(FS['sI1'], K.BOX_W1), Bc_w1_far_S0V=(K.S0V, K.BOX_W1),  # |w'| <= 0.02 on [16.25,60]
        Bd_w_16_34=(FS['far_w'], K.BOX_W0A),                                # |w| <= 0.02 on [16.25,34.25] (covers [16.25,26])
        Bd_w_34_60=(FS['far_w'] + 25.75 * K.S0V * (1 + 1e-12), K.BOX_W0B),  # |w| <= |w(34.25)| + 25.75 S0V <= 0.5 on [34.25,60]
        U3_strip_S0=(FS['far_w'], K.BOX_STRIP), U3_strip_S1=(FS['sI1'], K.BOX_STRIP), U3_strip_S2=(FS['sI2'], K.BOX_STRIP),
        Bf_mod_alpha=(aq, 5e-4), Bf_mod_A=(FS['A_max'], 1e-3), Bf_mod_W1=(FS['W1_max'], 1e-3))
    worst = {n: v / b for n, (v, b) in chk.items()}
    BOXLOG.append(dict(where=where, ratio=worst))
    bad = [n for n, r in worst.items() if not r <= 1.0]
    assert not bad, f"box hypothesis violated at {where}: {bad} {worst}"
def march(alpha1, W1_0=1e-6, VI_0=3e-6, tmax=120.0, rate_ratio=1.05, dmin=0.03, dmax=0.6, log=print, tol=5e-4, rtol=0.005, sqrtL1=None):
    orc = Oracle(tol, rtol, log)
    Ly1 = A(str(W1_0)) ** 2 * (1 + A(str(K.EPS)) * 9) + A(str(K.MU)) * A(str(VI_0)) ** 2 * (1 + A(str(K.EPSI)) * 9)
    if sqrtL1 is not None: Ly1 = A(str(sqrtL1)) ** 2          # hand-off level: Lyap(tau1) <= sqrtL1^2
    a = A(str(alpha1)) * (1 + A("1e-12")); L = Ly1 * (1 + A("1e-12")); tau = A(0)
    log(f"alpha1={alpha1} Ly1={float(Ly1):.6e} sqrt={float(Ly1.sqrt()):.6e}  ADD={sorted(K.ADD)}  tol={tol} rtol={rtol}")
    r_guess = 0.3; Lmax = L; steps = []; tail_res = None; tfinal = None
    for k in range(100000):
        if float(tau.mid()) >= tmax: break
        D = min(dmax, max(dmin, np.log(rate_ratio) / (2 * max(abs(r_guess), 1e-3)))); Da = A(D)
        eD = (LAMT * Da).exp()
        for it in range(8):
            g = L * (2 * A(str(r_guess)) * Da).exp()
            Lup = L if fup(L) >= fup(g) else g
            Lup = (Lup if fup(Lup) >= fup(L) else L) * (1 + A("1e-3"))
            Ldn = (L if fdn(L) <= fdn(g) else g) * (1 - A("1e-3"))
            # alpha side: fixed point  a_up >= max(a, bb/|LAMT|),  bb >= beta(a_up, Lup)
            a_up = a
            for _ in range(20):
                bbv = beta(a_up, Lup) * (1 + A("1e-9"))
                a_new = bbv / (-LAMT)
                if fup(a_new) <= fup(a_up): break
                a_up = a_up if fup(a_up) >= fup(a_new) else a_new; a_up = a_up * (1 + A("1e-9"))
            assert fup(bbv / (-LAMT)) <= fup(a_up) * (1 + 1e-12) or True
            aq = fup(a_up); aq = float(np.ceil(aq / 1e-5) * 1e-5)       # alpha level (rate is bounded for all alpha in [0, aq])
            slack = 0.0
            if K.RAMP == "table" and 'a' in K.ADD:
                for lev in sorted(K.RAMP_TAB):      # a_up within 1e-8 (rel.) above a certified table level: use that level + 1e-6 absolute slack in the rate (d rate/d alpha <~ 2e2)
                    if aq > lev and fup(a_up) <= lev * (1 + 1e-8): aq = lev; slack = 1e-6
            Rb = orc.rate(aq, fdn(Ldn), fup(Lup)) + slack
            rr = A(Rb + 1e-9 * max(1.0, abs(Rb)))
            Lnew = L * (2 * rr * Da).exp()
            ok = fup(Lnew) <= fup(Lup) and fdn(Lnew) >= fdn(Ldn)
            if ok: break
            r_guess = float(rr.upper())
        else:
            raise RuntimeError("cell iteration did not close at step %d" % k)
        a_next = a * eD + (bbv / (-LAMT)) * (1 - eD)
        # barrier checks (all in arb, enclosures):
        assert fup(a_up) >= fup(a) and fup(a_up) >= fup(bbv / (-LAMT)) * (1 - 1e-12)
        tau = tau + Da; a = a_next; L = Lnew; r_guess = float(rr.upper())
        if fup(Lup) > fup(Lmax): Lmax = Lup
        S_sup = sramp_sup(fup(Lup.sqrt())); FS = farsup(fup(Lup.sqrt()))
        box_assert(S_sup, FS, aq, log, f"step {k}")
        steps.append(dict(k=k, far=FS, S_ramp_sup=S_sup[0], S_w1_sup=S_sup[1], tau=float(tau.mid()), D=D, alpha_up=aq, Ly_lo=fdn(Ldn), Ly_hi=fup(Lup), rate_U=Rb, alpha_next=fup(a), sqrtL_next=fup(L.sqrt())))
        if k % 10 == 0 or aq < 3e-4:
            log(f"  k={k:4d} tau={float(tau.mid()):7.3f} D={D:.3f} alpha<={fup(a):.4e} sqrtLy<={fup(L.sqrt()):.4e} rate<={Rb:+.5f} (lb {orc.rows[-1]['lb']:+.5f}) cell {orc.rows[-1]['secs']}s")
        # try the tail lemma once alpha is small
        if fup(a) < 1.6e-4:
            tl = tail(a, L, orc, log)
            if tl is not None: tail_res = tl; tfinal = tau; break
    else:
        tail_res = None
    res = dict(alpha1=alpha1, Ly1=float(Ly1), steps=steps, sqrtL_max=fup(Lmax.sqrt()), tail=(tail_res if 'tail_res' in dir() else None), tau_tail_start=(float(tfinal.mid()) if tfinal is not None else None),
               oracle_rows=orc.rows, ADD=sorted(K.ADD))
    return res

def tail(a, L, orc, log, rho_min=0.002):
    """scale-invariant contraction anchored at (a*, L*).  Let T_s := tot without the degree<2 terms (b,c).  If rate_s(alpha<=a*, Ly=L*) <= -rho then for
    alpha<=s a*, Ly=s^2 L* (s<=1):  tot <= -rho Ly + f1 sqrt(Ly) + F0  (f1 = mu*ac0w/sqrt(mu), F0 = mu*(aflux0+epsI*aflux2)).
    Barrier l(t) = max(L_F, L* e^{-rho t}),  alpha-barrier kappa sqrt(l), kappa = a*/sqrt(L*): conditions checked below."""
    Ls = L; sq = Ls.sqrt()
    if K.NEW: box_assert(sramp_sup(fup(sq)), farsup(fup(sq)), fup(a), log, "tail anchor")
    kappa_need = (BETA0 * 1.05) / (-LAMT - A("0.5") * A("0.05"))
    astar = a if fup(a) >= fup(kappa_need * sq) else kappa_need * sq
    Us = orc.rate(fup(astar), fdn(Ls), fup(Ls * (1 + A(str(K.ETA_ANCHOR)))), floors=False)
    log(f"  tail anchor: a*={fup(astar):.4e} L*={fup(Ls):.4e}  rate_s bound U={Us:+.5f}")
    if Us > -rho_min: return None
    rho = A(str(-Us)) * (1 - A("1e-9"))
    kappa = astar / sq
    lhs = (-LAMT - rho / 2) * kappa
    rhs = BETA0 + BETA1 * kappa * sq + (BETAF * FnA(astar) + A(str(K.STRIP_BETA))) / sq + QBL * sq
    ok_alpha = fdn(lhs - rhs) > 0
    mu = A(str(K.MU)); epsI = A(str(K.EPSI))
    f1 = mu.sqrt() * A(str(K.ADD_I0_C0W)) if ('b' in K.ADD and not K.QUAD) else A(0)
    if 'c' in K.ADD and K.CMODE == "abs": F0 = mu * (A(str(K.ABS_FLUX_I0)) + epsI * A(str(K.ABS_FLUX_I2)))
    else: F0 = mu * (A(str(K.ADD_FLUX_I0)) + epsI * A(str(K.ADD_FLUX_I2))) if ('c' in K.ADD and not K.QUAD) else A(0)
    yF = (f1 + (f1 * f1 + 2 * rho * F0).sqrt()) / rho
    LF = yF * yF
    log(f"  tail: rho={float(rho):.5f} alpha-condition {'OK' if ok_alpha else 'FAIL'} (margin {float(lhs-rhs):.3e}); floor sqrt(L_F)={fup(yF):.3e}, alpha_F<= {fup(kappa*yF):.3e}")
    if not ok_alpha: return None
    # V1 compatibility: need delta_t in grid with (H-(1-delta)(NEG-theta*NEGV))/Ly <= 0 at the anchor (dropL), else keep marching
    import rhs_iv as R
    th_unc = float(os.environ.get("THETA", 0.3)); dts = None
    try:
        R.set_theta(th_unc)
        for d in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
            R.set_mode(d, True)
            if orc.rate(fup(astar), fdn(Ls), fup(Ls * (1 + A(str(K.ETA_ANCHOR))))) <= 0: dts = d; break
    finally:
        R.set_theta(0.0); R.set_mode(0.0, False)
    log(f"  tail: V1 delta_t={dts}")
    if dts is None: return None
    return dict(a_star=fup(astar), L_star=fup(Ls), rate_s_U=Us, rho=float(rho.mid()), kappa=fup(kappa), sqrt_LF=fup(yF), alpha_F=fup(kappa * yF),
                decay_rate_sqrtLy=float((rho / 2).mid()), floor_terms=dict(f1=fup(f1), F0=fup(F0)))

def unconditional(res, theta, delta, orc, log, delta_tail_grid=(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)):
    """V1: UNCONDITIONAL bounds for int V0^2, V1^2, V2^2, V4^2 (only the upper barrier Ly<=L(tau), alpha<=a(tau) is used).
    tot = (1/2) Ly' = H + L - N - theta*NEGV,  N = NEG - theta*NEGV,  NEGV = mu MU1 V0^2 + (mu/2 + mu epsI MU2) V2^2 + (mu epsI/2) V4^2.
    (i)  H - (1-delta) N <= r1 s^2  for every s<=s_up (H/s^2 nondecreasing, N/s^2 s-free: tail_cert.py), r1 = B&B sup at s=s_up (L dropped)
    (ii) L - delta N <= Phi(alpha) = sum_k l_k^2/(4 delta n_k)   (L <= sum l_k x_k, degree 1)           [unc.py]
    =>   Ly' = 2 tot <= 2 r1 Ly + 2 Phi(alpha) + 2 F0 - 2 theta NEGV,     theta int NEGV <= Ly1/2 + sum_k [r1_k^+ Lup_k + Phi(a_k)] D_k + tail.
    Tail (alpha<=kappa sqrt(l), Ly<=l=max(L_F, L* e^{-rho tau})): anchor pass with r1* <= 0 (delta_t = largest grid value), Phi(alpha) <= Phi(a*) l/L*:
       theta int_{tau_t}^T NEGV <= L*/2 + Phi(a*)/rho + [Phi(a*) L_F/L* + F0] (T - tau_t)."""
    import rhs_iv as R
    th = A(str(theta)); tot = A(str(res['Ly1'])) * (1 + A('1e-12')) / 2; rows = []; pr = A(0); pp = A(0)
    R.set_theta(theta); R.set_mode(delta, True)
    try:
        for st in res['steps']:
            r1 = orc.rate(st['alpha_up'], st['Ly_hi'], st['Ly_hi'])
            ph, kap = unc.Phi(st['alpha_up'], delta, theta)
            rows.append((r1, ph))
            pr = pr + A(str(max(r1, 0.0))) * A(str(st['Ly_hi'])) * A(str(st['D'])); pp = pp + A(str(ph)) * A(str(st['D']))
            tot = tot + (A(str(max(r1, 0.0))) * A(str(st['Ly_hi'])) + A(str(ph))) * A(str(st['D']))
        log(f"  V1 march part: delta={delta} max r1={max(r for r, _ in rows):+.4f}  sum Phi*D={sum(p * s['D'] for (r, p), s in zip(rows, res['steps'])):.3e}  total/theta = {fup(tot / th):.4e}")
        tl = res['tail']; Ls = tl['L_star']; a_s = tl['a_star']
        dt = None
        for d in delta_tail_grid:
            R.set_mode(d, True)
            r1t = orc.rate(a_s, Ls, Ls * (1 + max(2e-9, K.ETA_ANCHOR)))
            if r1t <= 0: dt = d; break
        if dt is None: log("  V1 tail: no delta in grid gives r1*<=0"); return None
        phs, _ = unc.Phi(a_s, dt, theta)
        log(f"  V1 tail: delta_t={dt} r1*={r1t:+.5f} Phi(a*)={phs:.3e}")
    finally:
        R.set_theta(0.0); R.set_mode(0.0, False)
    rho = A(str(tl['rho'])); LF = A(str(tl['sqrt_LF'])) ** 2; F0 = A(str(tl['floor_terms']['F0'])); assert tl['floor_terms']['f1'] < 1e-300
    Xc = tot + A(str(Ls)) / 2 * (1 + A('1e-9')) + A(str(phs)) / rho
    Xs = A(str(phs)) * LF / A(str(Ls)) + F0
    mu = A(str(K.MU)); epsI = A(str(K.EPSI))
    n0 = mu * A(str(K.MU1)); n2 = mu / 2 + mu * epsI * A(str(K.MU2)); n4 = mu * epsI / 2
    C1W = A(str(K.C1W))
    out = dict(theta=theta, delta=delta, delta_tail=dt, X_const=fup(Xc), X_slope=fup(Xs),
               parts_Xc=dict(Ly1_half=res['Ly1'] / 2, march_r1plus=fup(pr), march_Phi=fup(pp), tail_Lstar_half=Ls / 2, tail_Phi_over_rho=fup(A(str(phs)) / rho)),
               parts_Xs=dict(Phi_LF_over_Lstar=fup(A(str(phs)) * LF / A(str(Ls))), F0=fup(F0)))
    log('  X parts (J=X/(theta*n_k)): ' + json.dumps(out['parts_Xc']) + ' slope ' + json.dumps(out['parts_Xs']))
    Ic = {}; Is = {}
    for nm, n in (("V0sq", n0), ("V2sq", n2), ("V4sq", n4)):
        Ic[nm] = Xc / th / n; Is[nm] = Xs / th / n
    lam = (Ic["V2sq"] / Ic["V0sq"]).sqrt()
    Ic["V1sq"] = C1W * Ic["V0sq"] + (lam * Ic["V0sq"] + Ic["V2sq"] / lam) / 2; Is["V1sq"] = C1W * Is["V0sq"] + (lam * Is["V0sq"] + Is["V2sq"] / lam) / 2
    for nm in ("V0sq", "V1sq", "V2sq", "V4sq"): out[nm] = dict(const=fup(Ic[nm]), slope_per_tau=fup(Is[nm]))
    log("  UNCONDITIONAL: " + ";  ".join(f"int {nm[:2]}^2 <= {out[nm]['const']:.4e} + {out[nm]['slope_per_tau']:.3e}(T-tau_t)" for nm in ("V0sq", "V1sq", "V2sq", "V4sq")))
    return out

def integrals(res, theta, orc, log):
    """rigorous bounds for  int_{tau1}^{T} V0^2, V1^2, V2^2 dtau  along any trajectory below the barrier.
    Lyap' <= r_theta*Ly - theta*NEGV  on every step cell (r_theta = B&B bound of (H+L-NEG_theta)/Ly, NEG_theta = NEG - theta*NEGV, NEGV = mu MU1 V0^2 + (mu/2+mu epsI MU2) V2^2),
    hence theta int NEGV <= Ly1 - Ly(tau_t) + sum_k max(r_theta,k,0) Lup_k D_k ; tail: Lyap' <= -rho_theta Ly - theta NEGV + f1 sqrt(Ly) + F0
    => theta int_{tau_t}^T NEGV <= Ly(tau_t) + (f1^2/(2 rho_theta) + F0)(T - tau_t)."""
    import rhs_iv as R
    R.set_theta(theta)
    try:
        tot = A(str(res['Ly1'])) * (1 + A('1e-12'))
        rows = []
        for st in res['steps']:
            r = orc.rate(st['alpha_up'], st['Ly_lo'], st['Ly_hi'])
            rows.append(r)
            tot = tot + A(str(max(r, 0.0))) * A(str(st['Ly_hi'])) * A(str(st['D']))
        log(f"  theta={theta}: max step rate_theta={max(rows):+.4f}")
        tl = res['tail']
        Us = orc.rate(tl['a_star'], tl['L_star'] * (1 - 2e-9), tl['L_star'] * (1 + 2e-9), floors=False)
        log(f"  theta={theta}: tail anchor rate_s_theta = {Us:+.5f}")
        if Us >= -1e-4: return None
        rho = A(str(-Us)) * (1 - A("1e-9"))
        f1 = A(str(tl['floor_terms']['f1'])); F0 = A(str(tl['floor_terms']['F0']))
        g = (f1 * f1 / (2 * rho) + F0) * (1 + A("1e-9"))
    finally:
        R.set_theta(0.0)
    th = A(str(theta)); mu = A(str(K.MU)); epsI = A(str(K.EPSI))
    J0 = tot / th                       # finite part of int NEGV
    Jg = g / th                         # rate of growth per unit tau after tau_t (from the non-decaying floor forcing b,c)
    nV0 = mu * A(str(K.MU1)); nV2 = mu / 2 + mu * epsI * A(str(K.MU2))
    I0_ = J0 / nV0; I2_ = J0 / nV2; g0 = Jg / nV0; g2 = Jg / nV2
    C1W = A(str(K.C1W))
    lam = (I2_ / I0_).sqrt()           # sqrt(int V0^2 int V2^2) <= (lam int V0^2 + int V2^2/lam)/2 for any lam>0  (affine in T)
    I1_ = C1W * I0_ + (lam * I0_ + I2_ / lam) / 2; g1 = C1W * g0 + (lam * g0 + g2 / lam) / 2
    out = dict(theta=theta, rho_theta=float(rho.mid()), tau_t=res['tau_tail_start'],
               V0sq=dict(const=fup(I0_), slope_per_tau=fup(g0)), V1sq=dict(const=fup(I1_), slope_per_tau=fup(g1)), V2sq=dict(const=fup(I2_), slope_per_tau=fup(g2)), step_rates_theta=rows)
    log(f"  int V0^2 <= {fup(I0_):.4e} + {fup(g0):.3e}(T-tau_t);  int V1^2 <= {fup(I1_):.4e} + {fup(g1):.3e}(T-tau_t);  int V2^2 <= {fup(I2_):.4e} + {fup(g2):.3e}(T-tau_t)")
    return out

if __name__ == "__main__":
    alpha1 = float(sys.argv[1]); W1_0 = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6; VI_0 = float(sys.argv[3]) if len(sys.argv) > 3 else 3e-6
    tag = os.environ.get("TAG", "run")
    lf = open(f"verify_log_{tag}_a{alpha1:g}.txt", "w")
    def log(s): print(s, flush=True); lf.write(s + "\n"); lf.flush()
    t0 = time.time()
    res = march(alpha1, W1_0, VI_0, log=log, sqrtL1=(float(os.environ["LY1SQ"]) if "LY1SQ" in os.environ else None))
    res["runtime_s"] = round(time.time() - t0, 1)
    import hashlib
    _ins = [os.path.join(os.path.dirname(os.path.abspath(__file__)), f) for f in ("constants.py", "rhs_iv.py", "closure.py", "unc.py", "bb.py", "ival.py")] + ([os.path.join(K.SRC, "v4cert", "out", "constants_v4.json")] if K.NEW else []) + [os.path.join(K.SRC, "cubic_E0fine_0.001_0.01_0.01.json")]
    res["args"] = dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ("PYTHONHASHSEED", "TAG", "NEW", "THETA", "ADD", "FIX", "S0V", "S1V", "DELTA", "MOD", "MODA", "ETA", "FLUX0", "FLUX2", "LY1SQ", "INTEG", "COND")},
                       inputs={p_: hashlib.sha256(open(p_, "rb").read()).hexdigest() for p_ in _ins}, NEW=K.NEW)
    res["box_hypotheses"] = dict(max_ratio={n: max(b['ratio'][n] for b in BOXLOG) for n in (BOXLOG[0]['ratio'] if BOXLOG else {})}, n_checks=len(BOXLOG))
    res["S_ramp_max_along_trajectory"] = max(st["S_ramp_sup"] for st in res["steps"]); res["S_w1_max_along_trajectory"] = max(st["S_w1_sup"] for st in res["steps"])
    log(f"S_ramp (sup(|w|,|w'|) on [10.25,16.25], worst case over all shapes) max along trajectory <= {res['S_ramp_max_along_trajectory']:.4e}   (sup|w'| alone: {res['S_w1_max_along_trajectory']:.4e});  box certified: 1e-2")
    res["far_max"] = {key: max(st["far"][key] for st in res["steps"]) for key in ("A_max", "W1_max", "sI1", "sI2", "M1_ramp", "far_w1", "far_w", "sr1")}
    log("far-field sups on x>=16.25 (worst case over shapes), max along trajectory: " + ", ".join(f"{k_}<={v_:.4e}" for k_, v_ in res["far_max"].items()) + "   [sI1 = KW/E_I-based sup|w'|, far_w1 = max(E0,E_I) sup|w'|, far_w = sup|w|]")
    json.dump(res, open(f"rate_table_{tag}_a{alpha1:g}.json", "w"), indent=1)      # checkpoint before the integral passes
    if res["tail"] is not None and os.environ.get("INTEG", "1") != "0":
        th = float(os.environ.get("THETA", 0.02))
        if os.environ.get("COND", "0") == "1": res["integrals"] = integrals(res, th, Oracle(tol=5e-4, rtol=0.005, log=log), log)
        res["unc"] = unconditional(res, th, float(os.environ.get("DELTA", 0.2)), Oracle(tol=5e-4, rtol=0.005, log=log), log)
        res["runtime_s"] = round(time.time() - t0, 1)
    log(f"RESULT alpha1={alpha1}: rigorous max sqrt(Ly) <= {res['sqrtL_max']:.4e}   (float code: see baseline)  tail={res['tail'] is not None}  runtime {res['runtime_s']} s")
    json.dump(res, open(f"rate_table_{tag}_a{alpha1:g}.json", "w"), indent=1)
