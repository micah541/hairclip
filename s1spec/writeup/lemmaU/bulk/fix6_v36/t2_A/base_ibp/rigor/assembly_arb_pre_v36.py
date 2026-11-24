#!/usr/bin/env python3
"""Rigor pass, items 1 and 3: the final Lemma S assembly in python-flint arb (outward rounding) + the zeta<1 sliver.

Rerun (about 1 minute, after the two steps below have produced rigor/cap_F10p.json and rigor/weights_exact.json):
    cd /home/claude/s1spec/writeup/lemmaS
    PYTHONHASHSEED=0 python3 rigor/capture_F2.py          # ~2 min, reruns the final F2 command, records the ramp tables that F2 does not store
    PYTHONHASHSEED=0 python3 rigor/weights_exact.py       # exact-Rational / arb Ghat (core, ramp, 5 windows)
    PYTHONHASHSEED=0 python3 rigor/ramp_arb.py            # ~3 min, ramp_level in arb
    PYTHONHASHSEED=0 python3 rigor/assembly_arb.py        # writes rigor/assembly_arb.json, prints the verification table

What is ported (line by line, then checked against the doubles of lemmaS_F10p_M0.52.json / lemmaS_F4B_F10p_M0.52.json / lemmaS_Bbox_1e-05_a0.00015_dv0.15.json):
  F1.assemble (gamma', inflow, ramp commutators, forcing, coupling, S0/S1), the extra_outer (foot + S33) block of F2.main,
  the whole of lemmaS_Bbox.py (window loop, K2/K3/K4), the F4 and (B)-exit/H6 parts of lemmaS_F4B.py.
Inputs: every double read from the run JSONs / cflux / rate table is inflated by one ulp in the safe direction (upper bounds up, lower bounds down);
these doubles are themselves rounded certified values (hypotheses typed as decimals; verifier sups printed to <=17 digits; arb upper endpoints converted with float()),
hence the ulp.  Decimal constants of the proof (1.4, 1.02, 1.3056, 1.0246, 0.25, 0.08, 3e-5, 0.8477, ...) are taken as exact rationals.
Ghat = arb upper endpoints from rigor/weights_exact.py (exact Rationals).  Conversion-weight suprema: recomputed here in arb with exact ln 2 (400/1500 cells, upper endpoints);
the ramp-commutator tables come from rigor/ramp_arb.py (ramp_level ported to exact Rationals / arb, L = ln 2 ball): no doubles, no blanket factor.
"""
import os, sys, json, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); PAR = os.path.dirname(HERE)
os.chdir(PAR); sys.path.insert(0, HERE)
import sympy as sp
from flint import arb, ctx
import arbio
ctx.prec = 200

def D(s):
    f = Fraction(str(s)); return arb(f.numerator)/arb(f.denominator)
def UP(x): return D(repr(math.nextafter(float(x), math.inf)))      # double read from json -> upper bound (exact double, one ulp up)
def DN(x): return D(repr(math.nextafter(float(x), -math.inf)))
def fl(x): return float(x.upper())
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}; BELL = {0: 1, 1: 1, 2: 2, 3: 5}
STIRN = lambda k, m: STIR[k][m]
GAM = arb(1); P = D('1.0246'); X0 = D('2.2490776515'); LN2 = arb(2).log()
def U(y): return D('1.3056')*AL*(y.log()*P).exp()         # |U(y)| bound, AL set below
def sqrt(x): return x.sqrt()

# ----------------------------------------------------------------- arb ramp objects (exact ln 2)
t_ = sp.Symbol('t'); Sx_ = t_**3*(10 - 15*t_ + 6*t_**2)
E2F = [sp.lambdify(t_, sp.diff(Sx_**2, t_, l), modules=[{}]) for l in range(3)]      # unexpanded (factored) derivatives: far less cancellation on the cells than Horner
def e2(t, l): return E2F[l](t)/LN2**l
def yy(t): return 14*(t*LN2).exp()
def om(t): return (-(GAM)*((yy(t)/28).log())).exp()
def eps(t): return yy(t)**-4
def Wt(t): return (1 + (yy(t) + X0)**2)**-2
def sup_grid(f, n):
    best = arb(0)
    for i in range(n):
        u = abs(f((arb(i)/n).union(arb(i + 1)/n))).upper()
        if u > best: best = u
    return best
def frm_tab(n=400): return {m: sup_grid(lambda t, m=m: eps(t)*om(t)*yy(t)**(2*m - 1)/Wt(t), n) for m in (1, 2, 3)}
def infl_tab(j, n=1500): return {m: sup_grid(lambda t, m=m: abs(e2(t, 1)*om(t))*yy(t)**(2*m - 1)/Wt(t), n) for m in STIR[j]}

# ----------------------------------------------------------------- zeta<1 sliver (item 3)
def sliver_G(nsub=64):
    """G_i = sup over the sliver y in [14, 16.25-X0] of |D^i g - delta_{i0}|, g = 1/zeta, D = y d/dx, zeta(x) = S7((x-10.25)/6), S7(u)=u^4(35-84u+70u^2-20u^3) (C^3 smoothstep, verifier).
    With v=(16.25-x)/6 in [0, vmax]: zeta = 1 - S7(v) (S7(u)+S7(1-u)=1), d/dx = -(1/6) d/dv, D^i g = sum_m Stirling2(i,m) y^m g_x^(m)."""
    v = sp.Symbol('v'); h = v**4*(35 - 84*v + 70*v**2 - 20*v**3); g = 1/(1 - h)
    der = [sp.lambdify(v, sp.diff(g, v, m), modules=[{}]) for m in range(4)]
    vmax = (D('16.25') - (X0 + 14))/6; ymax = D('16.25') - X0
    # cells in v to limit the dependency overestimate
    sup = {m: arb(0) for m in range(4)}
    for c in range(nsub):
        ball = (vmax*c/nsub).union(vmax*(c + 1)/nsub)
        for m in range(4):
            val = abs(der[m](ball)).upper() if m else abs(der[0](ball) - 1).upper()
            val = val/6**m if m else val
            if val > sup[m]: sup[m] = val
    S2 = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
    G = {0: sup[0]}
    for i in (1, 2, 3): G[i] = sum((arb(S2[i][m])*ymax**m*sup[m] for m in S2[i]), arb(0))
    return G, vmax, ymax, sup

def convk(k, a, Jd, G=None):
    """BELL_k sum_m S(k,m) a_m Jd_m  (the code's conversion).  With G (sliver sups): (sqrt(base_k) + sum_i C(k,i) G_i sqrt(base_{k-i}))^2, see RIGOR_PASS.md item 3."""
    def base(kk):
        if kk == 0: a0 = a[0] if 0 in a else a[1]/196; return a0*Jd[0]       # a_0 <= a_1 / 14^2 (y >= 14)
        return BELL[kk]*sum((arb(STIR[kk][m])*a[m]*Jd[m] for m in STIR[kk]), arb(0))
    b = base(k)
    if G is None: return b
    d = sum((arb(math.comb(k, i))*G[i]*sqrt(base(k - i)) for i in range(k + 1)), arb(0))
    return (sqrt(b) + d)**2

# ----------------------------------------------------------------- inputs
FIN = json.load(open('lemmaS_F10p_M0.52.json')); CAP = json.load(open(os.path.join(HERE, 'cap_F10p.json'))); a_ = FIN['args']
WE = arbio.load(json.load(open(os.path.join(HERE, 'weights_exact.json')))['res'])
AL = UP(a_['alpha1']); SQ = UP(a_['sqLam']); ABAR = UP(a_['abar']); BETA = UP(a_['beta'])
J = [UP(x) for x in a_['J']] + [arb(0)]; SL = [UP(x) for x in a_['slopes']] + [arb(0)]; HO = [UP(x) for x in a_['handoff']]
J4 = UP(a_['J4']); SJ4 = UP(a_['sJ4']); C1W = D('0.8477')
assert [float(x) for x in a_['J']] == CAP['assemble'][1]['J'][:3] and a_['sqLam'] == CAP['assemble'][1]['sqLam']
CE = arbio.up(WE['C_eta'])

RT_ = arbio.load(json.load(open(os.path.join(HERE, 'ramp_tables.json'))))      # rigor/ramp_arb.py: ramp_level re-done in arb, no doubles, no blanket factor
RAMPT = {j: {int(k): {int(m): arbio.up(v) for m, v in d.items()} for k, d in RT_[str(j)].items()} for j in range(3)}
GHAT = {j: {int(k): arbio.up(v) for k, v in WE['core'][str(j)].items()} for j in range(3)}

def assemble(sliver=None, tag=''):
    """port of lemmaS_F1.assemble + F2.main extra_outer block (core box 0.11/0.52, abar/beta certified).  sliver = dict G_i or None."""
    y0, y1, YO, Rout = arb(14), arb(28), arb(60), arb(350)
    cbar = D('0.25') + ABAR + ABAR/y0; cmin = D('0.25') - ABAR - ABAR/y0
    eps_core = y1**-4; ee = 2*P - GAM; Iint = y1**GAM*((Rout**ee - y0**ee)/ee)
    r0 = D('1.4')*(BETA + D('1.02')*AL*ABAR)
    INFL = {j: infl_tab(j) for j in range(3)}; FR = frm_tab(400); E = []; parts = []; gps = []
    for j in range(3):
        gp = GAM*cmin - (abs(1 - 2*j) + 1)*ABAR/y0 - eps_core*GHAT[j][j]; assert gp.lower() > 0; gps.append(gp)
        Je = [J[m] + SL[m]/gp for m in range(4)]
        infl = cbar*convk(j, INFL[j], Je, sliver)
        rampc = sum((convk(k, RAMPT[j][k], Je, sliver) for k in range(j + 1)), arb(0))     # k=j+1 table is 0
        coup = sum((GHAT[j].get(k, arb(0))*eps_core*E[k] for k in range(j)), arb(0))
        F2 = (r0*P**j)**2*Iint
        # ---- extra_outer (F2.main, --rampl4 --rs33 ibp2)
        Jm = {m: J[m] + SL[m]/gp for m in range(3)}
        Jf2 = convk(2, FR, Jm, sliver); foot = arbio.up(WE['ramp'][str(j)]['footw'])*Jf2; extra = foot; info = dict(foot=foot, Jf2=Jf2)
        W33 = arbio.up(WE['ramp'][str(j)]['W33'])
        if W33 > 0:
            J4e = J4 + SJ4/gp; J3e = sqrt(Jm[2]*J4e) + C1W*Jm[2]; Jd3 = dict(Jm); Jd3[3] = J3e
            JG = 2*(convk(3, FR, Jd3, sliver) + D('6.25')*convk(2, FR, Jd3, sliver))
            L2 = SQ**2; V0s = L2/D('0.08'); V2s = L2/(D('0.08')*D('3e-5')); V1s = sqrt(V0s*V2s) + C1W*V0s
            Fmax = convk(2, FR, {0: V0s, 1: V1s, 2: V2s}, sliver)
            frmin = eps(arb(1))*om(arb(1)); Kf = W33*CE/3/sqrt(frmin)
            Dd = Kf*sqrt(Fmax)*(Jf2/sqrt(LN2) + sqrt(arb(2))*Jf2**D('0.75')*JG**D('0.25'))
            extra = extra + Dd; info.update(S33=Dd, JG=JG, Fmax=Fmax, Kf=Kf)
        Ej = HO[j] + infl + rampc + 4*F2/gp**2 + 2*coup/gp + extra
        E.append(Ej); parts.append(dict(gamma_prime=gp, inflow=infl, ramp_comm=rampc, forcing=4*F2/gp**2, coupling=2*coup/gp, handoff=HO[j], extra=extra, **info))
    conv = (YO/y1)**GAM; ell = (YO/y1).log(); Eu = [conv*e for e in E]
    S0 = sqrt(Eu[0]/ell + 2*sqrt(Eu[0]*Eu[1])); S1 = sqrt(Eu[1]/ell + 2*sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts, gp=gps, cmin=cmin, cbar=cbar, r0=r0, Iint=Iint)

def sup1(a0, a1, a2, ell):
    return sqrt(a0/ell + 2*sqrt(a0*a1)), sqrt(a1/ell + 2*sqrt(a1*(a1 + a2)))

def window_loop():
    """port of lemmaS_Bbox.py (e0=1e-5, dv0=0.15, al=1.5e-4, ABAR, BETA); Ghat rows from weights_exact.window()"""
    e0 = D('1e-5'); dv = D('0.15'); out = []
    r0 = D('1.4')*(BETA + D('1.02')*AL*ABAR); cmin = D('0.25') - ABAR - ABAR/50
    for key, w in WE['windows'].items():
        yE = arb(Fraction(key).numerator)/arb(Fraction(key).denominator); yb = D('1.5')*yE
        M0b = arb(Fraction(w['M0q']).numerator)/arb(Fraction(w['M0q']).denominator); M1b = arb(Fraction(w['M1q']).numerator)/arb(Fraction(w['M1q']).denominator)
        assert (U(yb) + dv - M0b).upper() <= 0 and (D('1.05')*U(yb) + dv - M1b).upper() <= 0     # rational box >= exact box
        Gk = {j: {int(k): arbio.up(v) for k, v in w['Ghat'][str(j)].items()} for j in range(4)}; YC3 = {int(k): arbio.up(v) for k, v in w['YC3'].items()}
        ep = yE**-4; gp = [GAM*cmin - (abs(1 - 2*j) + 1)*ABAR/yE - ep*Gk[j][j] for j in range(4)]
        ee = 2*P - GAM; IB = yE**GAM*((yb**ee - yE**ee)/ee)
        Es = [D('2e-7') if e0 > D('2e-7') else e0, e0, e0, e0]; Cexp = max([Gk[j][j] for j in range(4)], key=lambda x: float(x.upper()))*ep; EB = []
        for j in range(4):
            coup = sum((Gk[j].get(k, arb(0))*ep*EB[k] for k in range(j)), arb(0)) if j < 3 else sum((YC3[k]*ep*EB[k] for k in range(3)), arb(0))
            EB.append(Es[j]*Cexp.exp()*(8*ABAR/yE).exp() + 4*(r0*P**j)**2*IB/gp[j]**2 + 2*coup/gp[j])
        cv = yb/28; ell = D('1.5').log(); Eu = [cv*x for x in EB]
        def sup(a_, b_): return sqrt(a_/ell + 2*sqrt(a_*(a_ + b_)))          # Bbox.sup (note: a*(a+b))
        S = [sup(Eu[0], Eu[1]), sup(Eu[1], Eu[2]), sup(Eu[2], Eu[3]), sup(Eu[3], Eu[3])]
        der = (U(yb) + S[0], D('1.05')*U(yb) + S[1], P*P*U(yb) + S[2]); y = yE; L = arb(50)
        K2 = (L/y)*(D('1.05')*U(y) + S[1]); K3 = (L/y)**2*(S[2] + S[1] + U(y)*abs(P*P - P)); K4 = (L/y)**3*(S[3] + 3*S[2] + 2*S[1] + U(y)*abs(P*(P - 1)*(P - 2)))
        K2b = (L/yb)*(D('1.05')*U(yb) + S[1]); K2 = K2 if K2.upper() >= K2b.upper() else K2b
        out.append(dict(yE=key, box=(M0b, M1b), Ghat=[Gk[j][j] for j in range(4)], gp=gp, EB=EB, S=S, der=der, K=[K2, K3, K4], Ubox=(U(yb)+dv, D('1.05')*U(yb)+dv)))
    return out

def f4b(R):
    """port of lemmaS_F4B.py (F4 affine constants, (B)-exit y=48,50, core S_f, H6 strip)"""
    c = json.load(open('lemmaS_cflux.json')); RT = json.load(open('/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json')); tl = RT['tail']
    E = R['E']; parts = R['parts']; gp = R['gp']; Rout = arb(350); y1 = arb(28)
    W = UP(c['W57_upper']); k0 = [UP(x) for x in c['kappa_I0']]; k2 = [UP(x) for x in c['kappa_I2']]
    conv = (arb(60)/y1)**GAM; Eu = [conv*x for x in E]
    ee = 2*P - GAM; Iint = y1**GAM*((Rout**ee - arb(14)**ee)/ee)
    cF = D('1.4')*(BETA + D('1.02')*AL*ABAR)/SQ*P*sqrt(Iint)
    PHI_C = sum((max_(sqrt(UP(x['Ly_hi'])), UP(x['sqrtL_next']))*UP(x['D']) for x in RT['steps']), arb(0)) + 2*sqrt(UP(tl['L_star']))/DN(tl['rho']); PHI_S = UP(tl['sqrt_LF'])
    pr = [(p['inflow'] + p['ramp_comm'], p['gamma_prime'], p['forcing'] + p['coupling'], p['handoff']) for p in parts]
    Es, Ec, ssl = [], [], []
    for j in range(3):
        ir, g, fc, ho_ = pr[j]; r = max_(*[SL[m]/(g*J[m]) for m in range(3)]); fr = r/(1 + r)
        Es.append(fr*ir + fc); Ec.append(E[j] - Es[-1]); ssl.append(g*fr*ir)
    IE_c = [Ec[j]/pr[j][1] for j in range(3)]; IE_s = Es
    pa = parts[1]; bud_c = pa['handoff'] + pa['inflow'] + pa['ramp_comm'] + 2*cF*sqrt(E[1])*PHI_C; bud_s = ssl[1] + 2*cF*sqrt(E[1])*PHI_S
    amin = (1 + D('0.11')**2)**-2; yw = D('57.75'); fac = 1/((amin*yw**-4*(28/yw))*2)*(yw/28)**GAM
    Eu3c, Eu3s = bud_c*fac, bud_s*fac
    N0 = W*yw*Eu[0]; N1 = W*Eu[1]/D('54.75'); N2 = W*2*(Eu[2] + Eu[1])/D('54.75')**3
    N3c = W*3*(Eu3c + conv*(9*IE_c[2] + 4*IE_c[1]))/D('54.75')**5; N3s = W*3*(Eu3s + conv*(9*IE_s[2] + 4*IE_s[1]))/D('54.75')**5
    K0 = k0[0]*N0 + k0[1]*N1 + k0[2]*N2; K2 = k2[0]*N0 + k2[1]*N1 + k2[2]*N2
    out = dict(K0=K0, K2=K2, I3c=k2[3]*N3c, I3s=k2[3]*N3s, PHI_C=PHI_C, PHI_S=PHI_S, IE_c=IE_c, IE_s=IE_s)
    # (B) exit
    cmin = D('0.25') - ABAR - ABAR/50; r0 = R['r0']; Gc = [GHAT[j][j] for j in range(3)]; Bx = {}
    for yEi in (48, 50):
        yE = arb(yEi); ep = yE**-4; IB = yE**GAM*((Rout**ee - yE**ee)/ee); Cexp = max_(*Gc)*ep; EB = []; gpB = []
        for j in range(3):
            g = GAM*cmin - (abs(1 - 2*j) + 1)*ABAR/yE - ep*Gc[j]
            coup = sum((GHAT[j].get(k, arb(0))*ep*EB[k] for k in range(j)), arb(0))
            EB.append(E[j]*Cexp.exp()*(4*2*ABAR/yE).exp() + 4*(r0*P**j)**2*IB/g**2 + 2*coup/g); gpB.append(g)
        cv = (Rout/y1)**GAM; S0, S1 = sup1(cv*EB[0], cv*EB[1], cv*EB[2], (Rout/yE).log())
        Bx[yEi] = dict(gp=gpB, EB=EB, S0=S0, S1=S1)
    out['B'] = Bx
    # core S_f and H6 strip
    Um60 = U(arb(60)); out['S0f'] = Um60 + R['S0']; out['S1f'] = D('1.05')*Um60 + R['S1']
    eC = (max_(*Gc)*arb(50)**-4).exp(); cv = (arb(150)/y1)**GAM; ell = (arb(150)/50).log()
    S0h, S1h = sup1(cv*HO[0]*eC, cv*HO[1]*eC, cv*HO[2]*eC, ell); out['H6'] = dict(S0h=S0h, S1h=S1h, curvL={yv: arb(50)/yv*(D('1.05')*U(arb(yv))) + arb(50)/yv*S1h for yv in (50, 100, 150)})
    return out

def max_(*xs):
    b = xs[0]
    for x in xs[1:]:
        if x.upper() > b.upper(): b = x
    return b

if __name__ == '__main__':
    base = assemble(None); G, vmax, ymax, supg = sliver_G(); sl = assemble(G)
    print('sliver: v in [0, %.6e], y in [14, %.7f]; G_0..G_3 (arb upper) = %s' % (fl(vmax), fl(ymax), ['%.3e' % fl(G[i]) for i in range(4)]))
    # ----- F1.assemble numbers vs doubles
    stored = FIN['result']; sE = [float(x) for x in stored['E']]; rows = []
    print('\n== core assembly (arb) vs stored double ==')
    print('%-22s %-14s %-14s %-14s %s' % ('quantity', 'stored double', 'arb upper', 'arb+sliver up', 'rel.diff arb/double-1'))
    def line(name, dbl, val, vsl=None, lower=False):
        up_ = float(val.lower()) if lower else fl(val); rel = up_/dbl - 1 if dbl else 0.0
        print('%-22s %-14.8g %-14.8g %-14s %+.2e %s' % (name, dbl, up_, ('%.8g' % (float(vsl.lower()) if lower else fl(vsl))) if vsl is not None else '', rel, '<-- FLAG >1e-6' if abs(rel) > 1e-6 else ''))
        rows.append((name, dbl, up_, rel))
    for j in range(3):
        line('E_%d' % j, sE[j], base['E'][j], sl['E'][j])
        line("gamma'_%d (lower)" % j, float(stored['parts'][j]['gamma_prime']), base['gp'][j], sl['gp'][j], lower=True)
        for nm in ('inflow', 'ramp_comm', 'forcing', 'coupling'):
            dv_ = float(stored['parts'][j][nm])
            if dv_ > 0: line('  %s_%d' % (nm, j), dv_, base['parts'][j][nm], sl['parts'][j][nm])
        eo = stored['extra_outer'][str(j)]; line('  extra_outer_%d' % j, float(eo['foot']) + float(eo.get('S33data', 0.0)), base['parts'][j]['extra'], sl['parts'][j]['extra'])
    line('S0_v', stored['S0'], base['S0'], sl['S0']); line('S1_v', stored['S1'], base['S1'], sl['S1'])
    R = sl     # all downstream claims use the sliver-inclusive E (the larger one)
    print('\nsliver contribution to E_j (arb upper, absolute / relative):', ['%.3e / %.3e' % (fl(sl['E'][j] - base['E'][j]), fl((sl['E'][j] - base['E'][j])/base['E'][j])) for j in range(3)])
    # ----- (B) windows
    W_ = window_loop(); Bb = json.load(open('lemmaS_Bbox_1e-05_a0.00015_dv0.15.json'))
    print('\n== (B) windows (arb) vs Bbox json ==')
    for w, d in zip(W_, Bb):
        print('yE=%s' % w['yE'], "gamma' lower", ['%.6f/%.6f' % (x, float(g.lower())) for x, g in zip(d['gp'], w['gp'])], 'box', ['%.4f/%.4f' % (b, fl(c)) for b, c in zip(d['box'], w['Ubox'])])
        print('   EB', ['%.5g/%.5g' % (x, fl(y)) for x, y in zip(d['EB'], w['EB'])], 'derived', ['%.5f/%.5f' % (x, fl(y)) for x, y in zip(d['derived'], w['der'])], 'K', ['%.5f/%.5f' % (x, fl(y)) for x, y in zip(d['K'], w['K'])])
    F4 = f4b(R); Fj = json.load(open('lemmaS_F4B_F10p_M0.52.json'))
    print('\n== F4 / (B)-exit / H6 (arb, sliver-inclusive E) vs F4B json ==')
    f4d = Fj['F4_affine']
    line('K0_rate', f4d['K0_rate'], F4['K0']); line('K2_rate', f4d['K2_rate'], F4['K2']); line('I3_const', f4d['I3_const'], F4['I3c']); line('I3_slope', f4d['I3_slope'], F4['I3s'])
    line('PHI_C', f4d['PHI_C'], F4['PHI_C']); line('PHI_S', f4d['PHI_S'], F4['PHI_S'])
    for yv in (48, 50):
        b_ = Fj['B_y%d' % yv]; x = F4['B'][yv]
        for j in range(3): line('B y%d gp_%d (lower)' % (yv, j), b_['gamma_prime'][j], x['gp'][j], lower=True); line('B y%d E_B%d' % (yv, j), b_['E_B'][j], x['EB'][j])
        line('B y%d S0_B' % yv, b_['S0_B'], x['S0']); line('B y%d S1_B' % yv, b_['S1_B'], x['S1'])
    line('H6 S0h', Fj['H6']['S0h'], F4['H6']['S0h']); line('H6 S1h', Fj['H6']['S1h'], F4['H6']['S1h'])
    for yv in (50, 100, 150): line('H6 curvL(%d)' % yv, Fj['H6']['curvL'][str(float(yv))], F4['H6']['curvL'][yv])
    # ----- page-1 checks
    print('\n== page 1 (Summary v3) inequalities: quoted value vs arb ==')
    chk = []
    def claim(name, quoted, rel, val, note=''):
        v = val.upper() if rel == '<=' else val.lower(); ok = (D(quoted) - val).lower() >= 0 if rel == '<=' else (val - D(quoted)).lower() >= 0
        chk.append((name, quoted, rel, float(v), ok, note)); print('%-44s %s %-10s arb %-14.8g %s %s' % (name, rel, quoted, float(v), 'OK' if ok else '**QUOTE TOO TIGHT**', note))
    for j, q in enumerate(('3.12e-7', '1.81e-4', '0.248')): claim('E_%d' % j, q, '<=', R['E'][j])
    for j, q in enumerate(('0.246', '0.246', '0.143')): claim("gamma'_%d" % j, q, '>=', R['gp'][j])
    claim('core sup|v| = S0_v (page: 0.00575)', '0.00575', '<=', R['S0']); claim('core sup y|v\'| = S1_v (page: 0.171)', '0.171', '<=', R['S1'])
    claim('core sup|v| <= 0.02 (H-B2)', '0.02', '<=', R['S0']); claim('S0_v <= S0V = 0.0153 (verifier input)', '0.0153', '<=', R['S0']); claim('S1_v <= S1V = 0.6 (verifier input)', '0.6', '<=', R['S1'])
    M0c = F4['S0f']; M1c = F4['S1f']; claim('core M0 (page: 0.0188)', '0.0188', '<=', M0c); claim('core M1 (page: 0.1844)', '0.1844', '<=', M1c)
    claim('core M0 <= 0.11', '0.11', '<=', M0c); claim('core M1 <= 0.52', '0.52', '<=', M1c)
    for j, q in enumerate(('3.15e-7', '1.83e-4', '0.250')): claim('E_B%d (y_E=50)' % j, q, '<=', F4['B'][50]['EB'][j])
    claim("gamma'_B (page: >=0.2394, y_E=50, min over j)", '0.2394', '>=', min(F4['B'][50]['gp'], key=lambda x: float(x.lower())))
    claim('S0_B (y_E=50) <= S0V=0.0153', '0.0153', '<=', F4['B'][50]['S0']); claim('S1_B (y_E=50) <= S1V=0.6', '0.6', '<=', F4['B'][50]['S1'])
    wmax = lambda f: max((f(w) for w in W_), key=lambda x: float(x.upper()))
    claim('windows |f| (page 0.093)', '0.093', '<=', wmax(lambda w: w['der'][0])); claim("windows y|f'| (page 0.118)", '0.118', '<=', wmax(lambda w: w['der'][1])); claim('windows |z2| (page 0.118)', '0.118', '<=', wmax(lambda w: w['der'][2]))
    claim("windows gamma' (page >=0.181)", '0.181', '>=', min((g for w in W_ for g in w['gp']), key=lambda x: float(x.lower())))
    claim('windows |f| <= 1/2', '0.5', '<=', wmax(lambda w: w['der'][0])); claim('windows |z2| <= 0.5', '0.5', '<=', wmax(lambda w: w['der'][2]))
    for w in W_:
        claim('  window %s: |f| <= U(1.5y)+0.15 (box M0)' % w['yE'], '%.12f' % fl(w['Ubox'][0]), '<=', w['der'][0]); claim('  window %s: y|f\'| <= 1.05U+0.15 (box M1)' % w['yE'], '%.12f' % fl(w['Ubox'][1]), '<=', w['der'][1])
    claim('box M0 on [253,380] (page <=0.237)', '0.237', '<=', W_[-1]['Ubox'][0]); claim('box M1 on [253,380] (page <=0.241)', '0.241', '<=', W_[-1]['Ubox'][1])
    claim('K2 (page 0.024)', '0.024', '<=', wmax(lambda w: w['K'][0])); claim('K3 (page 0.025)', '0.025', '<=', wmax(lambda w: w['K'][1])); claim('K2 <= 0.12', '0.12', '<=', wmax(lambda w: w['K'][0])); claim('K3 <= 0.12', '0.12', '<=', wmax(lambda w: w['K'][1]))
    print('   (window K4 = %.4f > 0.06 is NOT used: page 1 certifies K4 <= 0.06 directly at the hand-off)' % fl(wmax(lambda w: w['K'][2])))
    claim('K0 (page 6.83e-22)', '6.83e-22', '<=', F4['K0']); claim('K2 flux (page 3.26e-21)', '3.26e-21', '<=', F4['K2']); claim('I3 const (page 1.84e-21)', '1.84e-21', '<=', F4['I3c']); claim('I3 slope (page 1.39e-27)', '1.39e-27', '<=', F4['I3s'])
    claim('Phi const, Lambda-integral (page 3.97e-5)', '3.97e-5', '<=', F4['PHI_C']); claim('Phi slope (page 4.94e-10)', '4.94e-10', '<=', F4['PHI_S'])
    claim('K0 <= 1e-3 * placeholder 6.093e-18', '6.093e-21', '<=', F4['K0']); claim('K2 <= 1e-3 * placeholder 1.5445e-17', '1.5445e-20', '<=', F4['K2'])
    # ramp boxes
    U28 = U(arb(28)); a2h = AL**2*D('785.7'); rM1 = D('0.4328') + P*U28 + D('2e-5'); rM0 = D('0.00116') + U28
    claim('ramp M1 cert 0.4328+pU(28)+2e-5 (page 0.4390)', '0.4390', '<=', rM1); claim('ramp M1 <= 0.9', '0.9', '<=', rM1); claim('  (alpha1^2 * 785.7 <= 2e-5 used above)', '0.00002', '<=', a2h)
    claim('ramp M0 cert 1.16e-3+U(28) (page 7.12e-3)', '0.00712', '<=', rM0); claim('ramp M0 <= 0.05', '0.05', '<=', rM0)
    claim('Ghat_22 core (page 6.58e4)', '65800', '<=', GHAT[2][2]); claim("Ghat_22 < cmin/eps_core = 1.54e5 (page)", '153664', '<=', GHAT[2][2])
    claim('ramp Ghat_22 (page 2.36e5)', '235593', '<=', arbio.up(WE['ramp']['2']['Ghat']['2']))
    print('sliver increments per part (arb upper, absolute):')
    for j in range(3): print('   level %d:' % j, {nm: '%.2e' % fl(sl['parts'][j][nm] - base['parts'][j][nm]) for nm in ('inflow', 'ramp_comm', 'extra')})
    print('\nmargins (arb, certified lower bounds): S0v: %.4f, M0core: %.4f, M1core: %.4f, ramp M1: %.4f, ramp M0: %.4f, Ghat22: %.4f, windows M1/M0/z2: %.4f/%.4f/%.4f, |f|<=1/2: %.4f, K2/K3: %.4f/%.4f' % (
        float((D('0.02')/R['S0']).lower()), float((D('0.11')/M0c).lower()), float((D('0.52')/M1c).lower()), float((D('0.9')/rM1).lower()), float((D('0.05')/rM0).lower()), float((D('153664')/GHAT[2][2]).lower()),
        min(float((w['Ubox'][1]/w['der'][1]).lower()) for w in W_), min(float((w['Ubox'][0]/w['der'][0]).lower()) for w in W_), float((D('0.5')/wmax(lambda w: w['der'][2])).lower()),
        float((D('0.5')/wmax(lambda w: w['der'][0])).lower()), float((D('0.12')/wmax(lambda w: w['K'][0])).lower()), float((D('0.12')/wmax(lambda w: w['K'][1])).lower())))
    json.dump(dict(rows=rows, claims=chk, sliver=dict(G=[fl(G[i]) for i in range(4)], vmax=fl(vmax), ymax=fl(ymax), dE=[fl(sl['E'][j] - base['E'][j]) for j in range(3)], dE_rel=[fl((sl['E'][j] - base['E'][j])/base['E'][j]) for j in range(3)]),
                   E=[fl(x) for x in sl['E']], gp=[float(x.lower()) for x in sl['gp']], E_nosliver=[fl(x) for x in base['E']]), open(os.path.join(HERE, 'assembly_arb.json'), 'w'), indent=1)
    bad = [c for c in chk if not c[4]]; print('\nclaims failing as quoted (rounded-to-nearest quotes below the arb bound):', [(c[0], c[1], c[3]) for c in bad])
