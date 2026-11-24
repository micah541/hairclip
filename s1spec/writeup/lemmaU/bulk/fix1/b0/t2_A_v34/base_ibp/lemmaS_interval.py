#!/usr/bin/env python3
"""Lemma S constants (python-flint arb).  `python3 lemmaS_interval.py` writes lemmaS_consts.json.
Coordinates: y = sigma = x - X0 (X0 = 2.2490776515 tip offset; W(x)=(1+x^2)^-2 is the closure weight in x).
Shell energies (log variable s=ln y, v~(s)=v(e^s)):  E0=int v^2 ds, E1=int v~_s^2 ds (=int y v'^2 dy), E2=int v~_ss^2 ds.
All sups of explicit functions are rigorous upper bounds (arb branch and bound)."""
import json, math
from flint import arb, ctx
ctx.prec = 100
X0 = arb('2.2490776515')

def sup_arb(f, a, b, nmax=4000):
    """rigorous upper bound of f on [a,b] (f takes arb, returns arb)"""
    stack = [(arb(a), arb(b))]; best_low = arb(-1e300); ub = arb(-1e300); n = 0
    work = [(a, b)]
    lo_best = -1e300; res = -1e300
    import heapq
    h = []
    def ev(x0, x1):
        m = (arb(x0) + arb(x1)) / 2; r = (arb(x1) - arb(x0)) / 2
        return f(arb(m, r)), f(arb(m))
    UB = lambda z: (lambda u: 1e300 if u != u else u)(float(z.upper()))
    box = ev(a, b); lo_best = float(box[1].lower())
    heapq.heappush(h, (-UB(box[0]), a, b))
    while h and n < nmax:
        negu, x0, x1 = heapq.heappop(h)
        if -negu <= lo_best * (1 + 1e-9) + 1e-300 and n > 0:
            res = -negu; break
        mid = (x0 + x1) / 2; n += 1
        for (u0, u1) in ((x0, mid), (mid, x1)):
            e, p = ev(u0, u1); lo_best = max(lo_best, float(p.lower()) if float(p.lower()) == float(p.lower()) else -1e300)
            heapq.heappush(h, (-UB(e), u0, u1))
    top = -h[0][0] if h else lo_best
    return arb(top)   # all remaining boxes have upper <= top

# quintic smoothstep
S  = lambda t: t**3 * (10 - 15*t + 6*t*t)
Sp = lambda t: 30 * t*t * (1 - t)**2
Spp= lambda t: 60 * t * (1 - t) * (1 - 2*t)

def eta2p_ramp(y, y0, y1):          # (eta^2)' on inner ramp [y0,y1]
    t = (y - y0) / (y1 - y0)
    return 2 * S(t) * Sp(t) / (y1 - y0)

def W(y):  # closure weight at x = y + X0
    return (1 + (y + X0)**2) ** -2

def inflow_K(y0, y1):
    out = {}
    for j, p in ((0, 0), (1, 2), (2, 4)):
        out['K%d' % j] = sup_arb(lambda y: eta2p_ramp(y, y0, y1) * y**p / W(y), y0, y1)
        out['Ku%d' % j] = sup_arb(lambda y: eta2p_ramp(y, y0, y1) * y**p, y0, y1)  # unweighted (W=1)
    return out

def coeff_sups(pmax):
    pm = arb(pmax)
    a = lambda p: (1 + p*p) ** -2
    b = lambda p: -3*p*(1 + p*p) ** -3
    amin = a(pm)
    bmax = sup_arb(lambda p: 3*p*(1 + p*p)**-3, 0, pmax)
    apmax = sup_arb(lambda p: 4*p*(1 + p*p)**-3, 0, pmax)
    return dict(a_min=amin, b_max=bmax, ap_max=apmax)

def cutoff_norms():
    # sup S', |S''|, S'^2/S   (t in (0,1)); S'^2/S -> 90 t at 0
    d = {}
    d['Sp_max'] = sup_arb(lambda t: Sp(t), 0, 1)
    d['Spp_max'] = sup_arb(lambda t: abs(Spp(t)), 0, 1)
    d['SpSq_over_S_max'] = sup_arb(lambda t: 900*t*(1-t)**4/(10 - 15*t + 6*t*t), 0, 1)   # S'^2/S, cancelled t^3
    d['int_Sp2'] = arb(10)/7
    d['eta2p_max_unit'] = sup_arb(lambda t: 2*S(t)*Sp(t), 0, 1)   # per unit ramp width
    return d

def f(x): return float(x.upper()) if isinstance(x, arb) else x

def main():
    out = {}
    out['ramp'] = [19, 24]
    K = inflow_K(19, 24)
    out['inflow_K'] = {k: f(v) for k, v in K.items()}
    out['cutoff'] = {k: f(v) for k, v in cutoff_norms().items()}
    out['coeff_p<=1/2'] = {k: f(v) for k, v in coeff_sups(0.5).items()}
    # alternative ramps (inside x in [20,30] <=> y in [17.75,27.75])
    alt = {}
    for (y0, y1) in ((19, 24), (20, 24), (21, 24), (18, 24)):
        kk = inflow_K(y0, y1)
        alt['%d-%d' % (y0, y1)] = {k: f(v) for k, v in kk.items() if k in ('K0', 'K1', 'K2')}
    out['ramp_options'] = alt
    # ---------- reference inputs ----------
    ref = dict(J0=3e-10, J1=1e-8, J2=2e-6, sqrtLam=2e-6, alpha1=1e-3)
    ell = math.log(60/24)
    def level_inflow(Kd, J0, J1, J2, abar):
        c = 0.25 + abar + 0 # c/y at the ramp (b/y <= bbar/19 added below)
        E0 = c * Kd['K0'] * J0
        E1 = c * Kd['K1'] * J1
        E2 = c * (2*Kd['K1'] * J1 + 2*Kd['K2'] * J2)
        return E0, E1, E2
    Kd = out['inflow_K']
    tab = []
    for sl in (2e-6, 1e-6, 5e-7, 2.5e-7):
        sc = (sl / 2e-6) ** 2
        J0, J1, J2 = ref['J0']*sc, ref['J1']*sc, ref['J2']*sc
        abar = 7.57 * sl + 7.57*sl/19
        E0, E1, E2 = level_inflow(Kd, J0, J1, J2, abar)
        # forcing (beta, a, b) -- see tex Sec. 4
        beta = (0.00244 + 63.5*ref['alpha1']) * sl
        psi3p60 = 93.0   # 1.4*60^1.0246 (coordinator quotes 87)
        R0 = beta * psi3p60 + 7.57*sl * (60 + 1) * 1.43*ref['alpha1'] + 0.2*ref['alpha1']*24**-2.97
        F0 = R0 * math.sqrt(math.log(60/19))
        gam = 1.0; gp = gam * (0.25 - abar) - 0.0
        Ef = 4 * F0**2 / gp**2
        conv = (60/24) ** gam * (24/19) ** gam * math.exp(1e5 * 24.0**-4 * 0.0)  # eps-factor handled separately
        S0 = math.sqrt(conv * (E0 + Ef) / ell + 2*math.sqrt(conv**2 * (E0 + Ef) * E1))
        S1_lower = math.sqrt(conv * E1 / ell)      # even if E2 = 0
        S1 = math.sqrt(conv * E1 / ell + 2*math.sqrt(conv**2 * E1 * max(E2, 1e-30)))
        tab.append(dict(sqrtLam1=sl, J0=J0, J1=J1, J2=J2, E0_in=E0, E1_in=E1, E2_in=E2, forcing_R0=R0, E0_forcing=Ef,
                        S0_v=S0, S1_v_floor_E2_zero=S1_lower, S1_v=S1))
    out['reference_table'] = tab
    # thresholds: what J1 (resp. J2) would give S1<=0.25 if E2 were zero, and the unweighted alternative
    c = 0.25
    out['needed'] = dict(
        J1_for_S1_le_0p25_floor=0.25**2 * ell / (c * Kd['K1'] * 1.5625),
        J1_unweighted_equiv_E1=c*out['inflow_K']['Ku1'],
    )
    json.dump(out, open('lemmaS_consts.json', 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))

if __name__ == '__main__':
    main()
