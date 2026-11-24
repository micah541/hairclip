#!/usr/bin/env python3
"""Certified (python-flint arb) evaluation of the constants of Lemma O (formulas of lemmaO_v2.tex / lemmaO_v3.tex)
and certified lower bounds for t*(X,M,K2,K3).  Inputs are exact rationals.  Every returned 'upper'/'lower' is a
rigorous bound given the formulas (the only non-closed-form step, max_u f(u), is done by branch and bound with
ball arithmetic plus an analytic tail argument)."""
import json, sys
from fractions import Fraction
from flint import arb, fmpq, ctx
ctx.prec = 200

def R(x):                                   # exact rational -> arb
    f = Fraction(x).limit_denominator(10**12) if not isinstance(x, Fraction) else x
    return arb(fmpq(f.numerator, f.denominator))

r      = (arb(16)/17)**3
A_MIN  = (arb(16)/17)**2                    # (17/16)^-2
A1, B0, B1 = r, arb(3)/4*r, arb(3)
CC1, CC2, CC3, CC4 = arb(7)/2*r, arb(5)/2*r, arb(19), arb(27)/2
S1, S2, S1SQ = arb(15)/8, 10/arb(3).sqrt(), arb(10)/7

def padd(a, b):
    n = max(len(a), len(b)); z = arb(0)
    return [(a[i] if i < len(a) else z) + (b[i] if i < len(b) else z) for i in range(n)]
def pmul(a, b):
    out = [arb(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out
def pscale(a, c): return [c*x for x in a]
def peval(p, u):
    s = arb(0)
    for c in reversed(p): s = s*u + c
    return s
def amin(a, b): return a if a.upper() <= b.upper() else b     # min of two enclosures of upper bounds

def fmax_upper(f, tol=1e-6):
    """Rigorous upper bound for max_{u>=0} f(u), f = -a u^8 + sum_{k<8} r_k u^k, r_k >= 0 (r_k = coefficient balls)."""
    assert len(f) == 9 and all(f[k].lower() >= 0 for k in range(8))
    a = -f[8]; assert a.lower() > 0
    U = arb(1)
    while True:                              # tail: for u>=U, f(u)/u^8 <= -a + sum r_k U^(k-8) < 0
        s = sum((f[k]/(U**(8-k)) for k in range(8)), arb(0))
        if (s - a).upper() < 0: break
        U = U*2
    best_lo = max(peval(f, U*arb(i)/64).lower() for i in range(65))
    stack = [(arb(0), U)]
    maxup = best_lo
    while stack:
        lo, hi = stack.pop()
        mid = (lo+hi)/2; rad = (hi-lo)/2
        up = peval(f, arb(mid.mid(), (rad+mid.rad()).upper())).upper()
        m_lo = peval(f, mid).lower()
        if m_lo > best_lo: best_lo = m_lo
        if up <= best_lo*(1+tol) + 1e-300:
            if up > maxup: maxup = up
            continue
        if rad.upper() < 1e-12:
            if up > maxup: maxup = up
            continue
        stack.append((lo, mid)); stack.append((mid, hi))
    return maxup if maxup > best_lo else best_lo

def consts(X, M, K2, K3, w, kappa_mode='min', t1=None):
    """All inputs Fractions. Returns dict of arb enclosures (upper bounds of the true constants)."""
    X, M, K2, K3, w = map(R, (X, M, K2, K3, w))
    lam = (2*X).sqrt(); wE = X-1-2*w
    assert wE.lower() > 0 and w.lower() > 0
    Z = S1/wE; Z2 = S2/(wE*wE); m = 2*wE; q = 12*Z*Z + 4*Z2
    lm = lam*M
    kappa = K3 + (lm if (kappa_mode == 'old' or lm.upper() < K3.lower()) else K3)   # ||h0'''||<=K3 by (A) at t=0
    # --- E1
    p1 = R(Fraction(1, 4)) + M; k2 = K2+M; sm = m.sqrt()
    rho0 = M*(1+CC2*K2) + B1*K2**3
    beta = 4*Z*k2*sm + lam*A1*K2*k2 + 2*Z*A1*K2*p1*sm + lam*k2*CC2*K2 + 2*Z*p1*CC2*K2*sm
    gamma = 2*(Z*Z+Z2)*p1*sm*kappa + 2*X*k2*rho0 + 2*Z*p1*rho0*m
    L1 = beta*beta/(2*A_MIN) + 2*gamma
    # --- E3
    K = kappa; sK = K.sqrt()
    Yb = [2*Z*K, arb(0), sK]
    om2 = [K*sK*(3*Z*K).sqrt(), K*sK*sK.sqrt()]
    nu0 = M + CC1*K2*lm + 2*CC2*lam*M*M + CC3*K2*K2*(lm+K) + lam*CC4*K2**4
    nu = padd(padd([nu0], pscale(Yb, CC1*K2)), pscale(om2, 2*CC2))
    lin = padd(pscale(Yb, 8*Z), [q*K])
    f = padd([arb(0)]*8+[-A_MIN], pmul(lin, [arb(0)]*4+[arb(1)]))
    f = padd(f, pmul(padd([arb(0)]*4+[arb(1)], lin), nu))
    f = f + [arb(0)]*(9-len(f))
    Cs = 2*fmax_upper(f)
    # --- interpolation
    Zq = S1/w; sL = L1.sqrt(); sC = Cs.sqrt(); ell = 2+2*w
    Yc = 2*Zq*sL + (L1*Cs).sqrt().sqrt()
    C1sq = amin(sL*(Yc+Zq*sL), L1/ell + 2*sL*Yc)
    C2sq = amin(Yc*(sC+Zq*Yc), Yc*Yc/ell + 2*Yc*sC)
    out = dict(kappa=kappa, Lambda1=L1, Cstar=Cs, Yc=Yc, C1=C1sq.sqrt(), C2=C2sq.sqrt(), C3=Cs, w=w, q=q, Z=Z)
    if t1 is not None:
        Gam = K3 + lam*B0*K2*K2
        out['C0'] = (2+w).sqrt()*sL + (2*S1SQ/w).sqrt()*Gam*R(t1).sqrt()/2
    return {k: v for k, v in out.items()}

def tstar(X, M, K2, K3, w, f=Fraction(1, 2), **kw):
    c = consts(X, M, K2, K3, w, **kw); fM = R(f)*R(M)
    t1 = fM*fM/(c['C1']*c['C1']); t2 = fM*fM/(c['C2']*c['C2']); t3 = 2*fM*fM/c['C3']
    lows = [t.lower() for t in (t1, t2, t3)]
    return dict(t1=lows[0], t2=lows[1], t3=lows[2], tmin=min(lows)), c

def best_w(X, M, K2, K3, grid=24, **kw):
    X = Fraction(X); best = None
    for i in range(1, grid+1):
        w = Fraction(i, grid) * (X-1)/2 * Fraction(95, 100)
        w = Fraction(round(w*1000), 1000)
        if w <= 0 or X-1-2*w <= Fraction(1, 20): continue
        ts, c = tstar(X, M, K2, K3, w, **kw)
        if best is None or ts['tmin'] > best[0]['tmin']: best = (ts, c, w)
    return best

def fl(x, n=5): return float(f'{float(x.mid()):.{n}g}')   # display only (midpoint of an enclosure)

def down(x): return float(f"{float(x.lower())*(1-1e-6):.6g}")        # safe decimal lower bound
def up(x):   return float(f"{float(x.upper())*(1+1e-6):.6g}")          # safe decimal upper bound

def case(name, X, M, K2, K3, f, need, grid=24):
    ts, c, w = best_w(X, M, K2, K3, f=f, grid=grid)
    t_need = R(need); fM = R(f)*R(M)
    cc = consts(X, M, K2, K3, w, t1=need)
    dev = {j: up(cc[j]*t_need.sqrt()) for j in ('C0', 'C1', 'C2')}
    dev['sqrtC3t'] = up((cc['C3']*t_need).sqrt())
    import math
    ncov = math.ceil(X)
    g3glob = up((ncov*cc['C3']*t_need).sqrt())      # ||g_xxx||_{L2(-X,X)} <= sqrt(ceil(X) C3 t) by covering with intervals [c-1,c+1]
    return dict(name=name, X=str(X), M=str(M), K2=str(K2), K3=str(K3), slack_factor=str(f), w=str(w),
        kappa_upper=up(cc['kappa']), C1_upper=up(cc['C1']), C2_upper=up(cc['C2']), C3_upper=up(cc['C3']),
        Lambda1_upper=up(cc['Lambda1']), t1_star_lower=down(ts['t1']), t2_star_lower=down(ts['t2']),
        t3_star_lower=down(ts['t3']), t_star_lower=down(ts['tmin']), needed=need,
        margin_factor_lower=down(ts['tmin'])/need,
        g3_global_term_at_needed_time=g3glob, K3_minus_g3_global_term=float(R(K3))-g3glob, covering_number=ncov,
        slope_actual_margin_lower=float(Fraction(1, 4))-float(M)-dev['C1'],
        slack_pointwise=float(R(f)*R(M)), slack_L2_g3=float(fM*arb(2).sqrt()),
        deviation_upper_at_needed_time=dev,
        slope_margin_to_quarter_lower=float(Fraction(1, 4)-(1+f)*M))

def main():
    F = Fraction; out = dict(note="all *_lower / *_upper are rigorous enclosure bounds (python-flint arb, 200 bits), rounded outward at 1e-6; "
        "slack = f*M pointwise (f=1/2 => K=1.5M); t* = min(t1*,t2*,t3*) of lemmaO_v3.tex Sec.3", cases=[], sweeps=[], far=[])
    M = F(12, 100)
    out['cases'].append(case('S1 cone, X=3', F(3), M, F(3, 2)*M, F(3, 2)*M, F(1, 2), 1.2e-7))
    out['cases'].append(case('S1 cone, X=2.4', F(12, 5), M, F(3, 2)*M, F(3, 2)*M, F(1, 2), 5e-8))
    for X in (F(3), F(12, 5)):
        for Mi in (F(8, 100), F(1, 10), F(12, 100), F(15, 100), F(2, 10), F(3, 10)):
            ts, c, w = best_w(X, Mi, F(3, 2)*Mi, F(3, 2)*Mi, grid=14)
            out['sweeps'].append(dict(X=str(X), M=str(Mi), K="1.5M", w=str(w), t_star_lower=down(ts['tmin'])))
    tau = 1e-9
    for X in (F(2), F(12, 5)):
        for f in (F(1, 2), F(1)):
            lo, hi = F(1, 100), F(10)
            for _ in range(10):
                mid = (lo+hi)/2
                if best_w(X, mid, (1+f)*mid, (1+f)*mid, f=f, grid=10)[0]['tmin'] >= tau: lo = mid
                else: hi = mid
            lo = F(round(lo*100), 100)       # round down-ish to a 2-decimal rational and re-certify below
            while best_w(X, lo, (1+f)*lo, (1+f)*lo, f=f, grid=10)[0]['tmin'] < tau: lo -= F(1, 100)
            out['far'].append(case(f'far, X={X}, f={f}, largest certified M', X, lo, (1+f)*lo, (1+f)*lo, f, tau, grid=10))
    for X in (F(2), F(12, 5)):
        for Mi, Ki in ((F(20, 3), F(10)), (F(5), F(10)), (F(10), F(15))):
            f = (Ki-Mi)/Mi
            ts, c, w = best_w(X, Mi, Ki, Ki, f=f, grid=10)
            out['far'].append(dict(name=f'far corner X={X} M={float(Mi):.3g} K={float(Ki)}', t_star_lower=down(ts['tmin']),
                                   factor_below_1em9=1e-9/down(ts['tmin'])))
    json.dump(out, open('lemmaO_v3_consts.json', 'w'), indent=1); print(json.dumps(out, indent=1))

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'run':
    main()
elif __name__ == '__main__':
    F = Fraction
    # cross-check against the floating-point v2 script (old kappa)
    import lemmaO_consts_v2 as v2
    M, K, X, w = F(1, 10), F(15, 100), F(12, 5), F(1, 5)
    c = consts(X, M, K, K, w, kappa_mode='old'); d = v2.constants(float(K), float(K), float(M), float(X), w=float(w))
    for k, kk in (('Lambda1', 'Lambda1'), ('C3', 'C3'), ('C1', 'C1'), ('C2', 'C2')):
        print(f'check {k}: interval upper {float(c[k].upper()):.6g}  float v2 {d[kk]:.6g}  ratio {float(c[k].upper())/d[kk]:.6f}')
