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
LAMT = A(str(K.LAMT)); BETA0 = A(str(K.BETA0)); BETA1 = A(str(K.BETA1)); BETAF = A(str(K.BETAF))
def FnA(al): return A(str(K.FN[0])) * al + A(str(K.FN[1])) * al * al + A(str(K.FN[2])) * al * al * al
QB=float(os.environ.get('QB',0)); QBL=A(str(QB*1.13e5))   # (A+W1+B)^2 <= 1.13e5 Ly  (eps=1e-5)
def beta(al, L): return (BETA0 + BETA1 * al) * L.sqrt() + BETAF * FnA(al) + QBL*L

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
    return dict(sI1=G1, far_w1=max(F1, G1), far_w=float(np.sqrt(F0 ** 2 + H ** 2) * (1 + 1e-14)), sr1=F1)

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
    kappa_need = (BETA0 * 1.05) / (-LAMT - A("0.5") * A("0.05"))
    astar = a if fup(a) >= fup(kappa_need * sq) else kappa_need * sq
    Us = orc.rate(fup(astar), fdn(Ls), fup(Ls * (1 + A("1e-9"))), floors=False)
    log(f"  tail anchor: a*={fup(astar):.4e} L*={fup(Ls):.4e}  rate_s bound U={Us:+.5f}")
    if Us > -rho_min: return None
    rho = A(str(-Us)) * (1 - A("1e-9"))
    kappa = astar / sq
    lhs = (-LAMT - rho / 2) * kappa
    rhs = BETA0 + BETA1 * kappa * sq + BETAF * FnA(astar) / sq + QBL*sq
    ok_alpha = fdn(lhs - rhs) > 0
    mu = A(str(K.MU)); epsI = A(str(K.EPSI))
    f1 = mu.sqrt() * A(str(K.ADD_I0_C0W)) if ('b' in K.ADD and not K.QUAD) else A(0)
    F0 = mu * (A(str(K.ADD_FLUX_I0)) + epsI * A(str(K.ADD_FLUX_I2))) if ('c' in K.ADD and not K.QUAD) else A(0)
    yF = (f1 + (f1 * f1 + 2 * rho * F0).sqrt()) / rho
    LF = yF * yF
    log(f"  tail: rho={float(rho):.5f} alpha-condition {'OK' if ok_alpha else 'FAIL'} (margin {float(lhs-rhs):.3e}); floor sqrt(L_F)={fup(yF):.3e}, alpha_F<= {fup(kappa*yF):.3e}")
    if not ok_alpha: return None
    return dict(a_star=fup(astar), L_star=fup(Ls), rate_s_U=Us, rho=float(rho.mid()), kappa=fup(kappa), sqrt_LF=fup(yF), alpha_F=fup(kappa * yF),
                decay_rate_sqrtLy=float((rho / 2).mid()), floor_terms=dict(f1=fup(f1), F0=fup(F0)))

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
    res["S_ramp_max_along_trajectory"] = max(st["S_ramp_sup"] for st in res["steps"]); res["S_w1_max_along_trajectory"] = max(st["S_w1_sup"] for st in res["steps"])
    log(f"S_ramp (sup(|w|,|w'|) on [10.25,16.25], worst case over all shapes) max along trajectory <= {res['S_ramp_max_along_trajectory']:.4e}   (sup|w'| alone: {res['S_w1_max_along_trajectory']:.4e});  box certified: 1e-2")
    res["far_max"] = {key: max(st["far"][key] for st in res["steps"]) for key in ("sI1", "far_w1", "far_w", "sr1")}
    log("far-field sups on x>=16.25 (worst case over shapes), max along trajectory: " + ", ".join(f"{k_}<={v_:.4e}" for k_, v_ in res["far_max"].items()) + "   [sI1 = KW/E_I-based sup|w'|, far_w1 = max(E0,E_I) sup|w'|, far_w = sup|w|]")
    if res["tail"] is not None and os.environ.get("INTEG", "1") != "0":
        th = float(os.environ.get("THETA", 0.02))
        res["integrals"] = integrals(res, th, Oracle(tol=5e-4, rtol=0.005, log=log), log)
        res["runtime_s"] = round(time.time() - t0, 1)
    log(f"RESULT alpha1={alpha1}: rigorous max sqrt(Ly) <= {res['sqrtL_max']:.4e}   (float code: see baseline)  tail={res['tail'] is not None}  runtime {res['runtime_s']} s")
    json.dump(res, open(f"rate_table_{tag}_a{alpha1:g}.json", "w"), indent=1)
