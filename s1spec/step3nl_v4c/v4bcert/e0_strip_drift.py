"""v4b E0 certificate, box part (python-flint arb, outward).  [cert]: every output is a rigorous upper bound over the boxes, given the
certified inputs (prof_tight hull boxes, step-1 tail profile bounds hstail.profile_bounds, the psi3^ / h2^ jet modules) -- the same inputs and
the same box covers as the closed alin_E0_spec.py (stage C: prof_tight hull groups on x <= tau_tip; stage B: sigma chunks of 1/N).
Parameters: AMAX = 2.5e-4 (|alpha| level; asserted in closure.py), WB = 0.02 (|w| <= WB on [0, 34.25]; asserted in closure.py).

Outputs (out/e0_strip_drift.json), all per unit |alpha| where stated:
 core (sigma <= 29, chi = 1):
   C_DR_CORE  = sup |c_dr| / |alpha|, c_dr = (1/8)[q' G - phi' q G + q (m u0 - k u0')/(1-k u0)^2], G = k u0/(1-k u0)
                (drift alpha-form: int_core rho w Dlin(w) = int_core rho c_dr w^2 <= C_DR_CORE |alpha| A^2;  identity: tc_weighted.py (3))
 strip (sigma in [29, 32], weight X = chi^2):
   M_STRIP[ij] = sup over the strip of the alpha-linear coefficients of the weighted divergence form T_c^X (tc_weighted.py (1)) plus, for ij = 00,
                the drift coefficient c_dr (with X, X1):  m00 = |b00| + |c_dr|, m01 = |b01 + b10|, m02 = |b02 + b20|, m11 = |b11|,
                m12 = |b12 + b21|, m22 = |b22|  (b_ij = a^X_ij(u0) - a^X_ij(0), |b_ij|/|alpha| <= |b'_ij(0)| + AMAX |b''_ij(A)|/2 as in alin_E0_spec.py)
   strip integrals: IR = int_strip rho, IRW = (int_strip rho^2/W)^(1/2), RW = sup_strip rho/W  (W >= 1/KW on the strip, W decreasing)
 global (x in [0, 34.25], |w| <= WB, |alpha| <= AMAX):
   AG_MIN  = inf (1 - |k| (|u0| + WB))   (lower bound of alpha_g(u0 + w))
   cubic drift  N_dr2 = -(q/4)[w'(G(u0+w) - G(u0)) + u0'(G(u0+w) - G(u0) - G_u(u0) w)]:  |w N_dr2| <= |w| (CD1 |w||w'| + CD0 w^2),
       CD1 = sup |q||k| / (4 ag^2),  CD0 = sup |q||u0'| k^2 / (4 ag^3)
   (a,b) w-parts: |A1(u) - A1(u0)| <= |w| + CA1 |w'| + CA0 |w|,  |B1(u) - B1(u0)| <= CB1 |w'| + CB0 |w|,
       CA1 = sup |q|/ag, CA0 = sup |q||u0'||k|/ag^2, CB1 = 1/AG_MIN (|cos theta| <= 1), CB0 = sup |u0'||k|/ag^2
   alpha^2 remainders of the (a,b) modulation directions (Q1 uses their first order alpha gA, alpha gB):
       RA2 = ||chi (A1(u0) + alpha gA)||_rho / alpha^2 (first order of A1(u0) is -alpha gA),  A1(u0) + alpha gA = alpha^2 (q h2' - h2) + q u0' k u0/alpha_g(u0),  gA = chi (psi - q psi')
       RB2 = sup_theta ||chi (B1(u0) - alpha gB)||_rho / alpha^2 <= || h2' + u0' k u0/(alpha^2 alpha_g(u0)) ||  (|cos theta| <= 1), gB = chi psi' cos theta
       (L2(rho) norms over x in [0, 34.25] as box sums  sum len * sup rho * sup |.|^2, chi <= 1)
Run: cd v4bcert; PYTHONHASHSEED=0 nice -n 19 python3 -B e0_strip_drift.py"""
import sys, os, json, time, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump, fup, fdn, EXT, V4B
from flint import arb, arb_series, ctx
import sympy as sp
ctx.prec = 120
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
t0 = time.time()
AMAX = arb('2.5e-4'); WB = arb('0.02'); NB = 50; NS = 200
Aball = arb(0, AMAX)
ab = lambda z: arb(abs(z).upper())
C4 = json.load(open(os.path.join(V4B, 'v4cert', 'out', 'constants_v4.json')))['constants']
KW = arb(repr(C4['KW']['value'])); assert C4['KW']['dir'] == 'up'
# ---------------------------------------------------------------- symbolic coefficient functions (from tc_weighted.py)
PK = pickle.load(open(os.path.join(OUT, 'tc_weighted.pkl'), 'rb'))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        if ex.is_Rational and ex.q == 2:
            n = int(ex.p); base = f"({bb}.sqrt())"; prod = "*".join([base] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        raise ValueError(f"power {ex}")
ARGS = "phi1, phi2, k, m, p, X, X1, X2, U0, U1, U2"
allowed = set(sp.symbols(ARGS.replace(' ', '')))
TF = {}
for ij, e in PK['TF'].items():
    if e == 0: continue
    assert e.free_symbols <= allowed, (ij, e.free_symbols - allowed)
    TF[ij] = eval(f"lambda {ARGS}: " + _Pr().doprint(e), {"arb": arb})
# ---------------------------------------------------------------- per-box quantities
def phis(X):
    E = (1 + X * X) ** (arb(-1) / 3)
    return X * E / 3, E / 3 - 2 * X * X * E * E * E * E / 9
def drift_coef(Xw, X1w, ph1, k, m, p, q, y0, y1, h0, h1):
    """c_dr / alpha for all alpha in Aball (direct enclosure; c_dr is O(alpha) with c_dr/alpha a rational function of alpha)"""
    u0 = Aball * y0 + Aball * Aball * h0
    den = 1 - k * u0
    assert den.lower() > 0
    Ga = k * (y0 + Aball * h0) / den                          # G(u0)/alpha
    q1 = 1 + k * p                                            # q'
    bra = q1 * Ga - ph1 * q * Ga + q * (m * (y0 + Aball * h0) - k * (y1 + Aball * h1)) / (den * den)
    return Xw * bra / 8 + X1w * q * Ga / 8
def glob_coefs(k, q, y0, y1, h0, h1):
    u0a = ab(AMAX * ab(y0) + AMAX * AMAX * ab(h0)); u1a = ab(AMAX * ab(y1) + AMAX * AMAX * ab(h1))
    kk, qq = ab(k), ab(q)
    ag = 1 - kk * (u0a + WB); assert ag.lower() > 0
    ag = arb(ag.lower())
    return dict(AG=ag, CD1=qq * kk / (4 * ag * ag), CD0=qq * u1a * kk * kk / (4 * ag * ag * ag), CA1=qq / ag, CA0=qq * u1a * kk / (ag * ag), CB0=u1a * kk / (ag * ag))
def rem_sq(k, q, y0, y1, h0, h1):
    """sup over the box of |A1(u0) - alpha gA|^2/alpha^4 and |u0'/alpha_g - alpha psi'|^2/alpha^4 (both for all |alpha| <= AMAX)"""
    u0 = Aball * y0 + Aball * Aball * h0; den = 1 - k * u0; assert den.lower() > 0
    t = q * (y1 + Aball * h1) * k * (y0 + Aball * h0) / den                 # q u0' k u0 / (alpha^2 alpha_g)
    ra = ab(q * h1 - h0 + t); rb = ab(h1 + (y1 + Aball * h1) * k * (y0 + Aball * h0) / den)
    return ra * ra, rb * rb
def curv_strip(y, h, ph1, ph2, kk, mm, pp, Xw, X1w, X2w):
    z = arb(0); c = {}
    u0s = [arb_series([z, y[j], h[j]], prec=3) for j in range(3)]
    uAs = [arb_series([AMAX * 0 + Aball * y[j] + Aball * Aball * h[j], y[j] + 2 * Aball * h[j], h[j]], prec=3) for j in range(3)]
    for ij, f in TF.items():
        c0 = f(ph1, ph2, kk, mm, pp, Xw, X1w, X2w, *u0s).coeffs() + [arb(0)] * 3
        cA = f(ph1, ph2, kk, mm, pp, Xw, X1w, X2w, *uAs).coeffs() + [arb(0)] * 3
        c[ij] = c0[1] + arb(0, AMAX) * cA[2]
    g = lambda i, j: c.get((i, j), arb(0))
    return {"00": ab(g(0, 0)), "01": ab(g(0, 1) + g(1, 0)), "02": ab(g(0, 2) + g(2, 0)), "11": ab(g(1, 1)), "12": ab(g(1, 2) + g(2, 1)), "22": ab(g(2, 2))}
Su = sp.symbols('u'); Sm = Su ** 4 * (35 - 84 * Su + 70 * Su ** 2 - 20 * Su ** 3)            # chi_32(sigma) = Sm((32 - sigma)/3) (agent_mod/mod_gram.py:60)
dSm = [eval("lambda u: " + _Pr().doprint(sp.expand(sp.diff(Sm, Su, j))), {"arb": arb}) for j in range(3)]
def chi_jets(sa, sb):
    """enclosures of chi, chi', chi'' (d/dsigma = d/dx) for sigma in [sa, sb] subset [29, 32]"""
    tl = (32 - sb) / 3; th = (32 - sa) / 3
    T = arb((tl + th) / 2, 0).union(tl).union(th)
    return dSm[0](T), -dSm[1](T) / 3, dSm[2](T) / 9
# exact sup |chi'| = 35/48 at t = 1/2 (Sm' = 140 t^3 (1-t)^3)
assert sp.simplify(sp.diff(Sm, Su) - 140 * Su ** 3 * (1 - Su) ** 3) == 0
CHI1 = arb(35) / 48
# ---------------------------------------------------------------- run
sup = {}
def upd(key, v):
    vu = v.upper(); assert vu.is_finite()
    if key not in sup or vu > sup[key]: sup[key] = arb(vu)
def updmin(key, v):
    vl = v.lower()
    if key not in sup or vl < sup[key]: sup[key] = arb(vl)
argmax = {}
# stage C (prof_tight hull groups, x <= tau_tip ~ 10.25), as alin_E0_spec.py
lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
GROUP = 10; nC = 0; RA = arb(0); RB = arb(0)
for g0 in range(0, len(lines), GROUP):
    grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
    ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
    hl = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
    Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
    if Xb.upper() < 0: continue
    Xb = arb(0).union(Xb) if Xb.lower() < 0 else Xb
    pp_, qq_, kk_, mm_ = hl
    ph1, ph2 = phis(Xb)
    for Pf, Ph in zip(F.PIECES, Hf.PIECES):
        if Xb.upper() < Ph.a.lower() or Xb.lower() > Ph.b.upper(): continue
        y = HR.jets(F, Pf, Xb, 3); h = HR.jets(Hf, Ph, Xb, 3)
        cd = drift_coef(arb(1), arb(0), ph1, kk_, mm_, pp_, qq_, y[0], y[1], h[0], h[1])
        upd('C_DR_CORE', ab(cd)); nC += 1
        gc = glob_coefs(kk_, qq_, y[0], y[1], h[0], h[1])
        for kq, v in gc.items(): (updmin if kq == 'AG' else upd)(kq if kq == 'AG' else kq, v)
        ra, rb = rem_sq(kk_, qq_, y[0], y[1], h[0], h[1]); Lx = 2 * arb(Xb.rad()); rs = arb(HR.rho_of(max(float(Xb.lower()), 0.0)).upper())
        RA += Lx * rs * ra; RB += Lx * rs * rb
print(f"stage C: {nC} boxes, C_DR_CORE <= {float(sup['C_DR_CORE'].upper()):.6g}  [{time.time()-t0:.0f}s]", flush=True)
# stage B (sigma chunks [8, 29], 1/NB), as alin_E0_spec.py; then the strip [29, 32] (1/NS) with chi
IR = arb(0); IR2 = arb(0); RWs = arb(0); nB = 0; nS = 0
for ci in range(8, 32):
    strip = ci >= 29; N = NS if strip else NB
    for i in range(N):
        sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
        kb, mb, pb, qs = HR.HS.profile_bounds(sa)
        xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
        X = arb(xc, (xb - xa) / 2); ph1, ph2 = phis(X)
        y = HR.jets_any(F, X, xc)[:3]; h = HR.jets_any(Hf, X, xc)[:3]
        kk_, mm_, pp_ = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())
        sig = arb((sa + sb) / 2, 0).union(sa).union(sb); qq_ = sig + arb(0, qs.upper())
        gc = glob_coefs(kk_, qq_, y[0], y[1], h[0], h[1])
        for kq, v in gc.items(): (updmin if kq == 'AG' else upd)(kq, v)
        ra, rb = rem_sq(kk_, qq_, y[0], y[1], h[0], h[1]); Lx = arb(xb) - arb(xa); rs = arb(HR.rho_of(xa).upper())
        RA += Lx * rs * ra; RB += Lx * rs * rb
        if not strip:
            cd = drift_coef(arb(1), arb(0), ph1, kk_, mm_, pp_, qq_, y[0], y[1], h[0], h[1]); upd('C_DR_CORE', ab(cd)); nB += 1
            continue
        ch, ch1, ch2 = chi_jets(sa, sb)
        Xw = ch * ch; X1w = 2 * ch * ch1; X2w = 2 * ch1 * ch1 + 2 * ch * ch2
        cs = curv_strip(y, h, ph1, ph2, kk_, mm_, pp_, Xw, X1w, X2w)
        cd = ab(drift_coef(Xw, X1w, ph1, kk_, mm_, pp_, qq_, y[0], y[1], h[0], h[1]))
        upd('DR_STRIP', cd)
        for kq, v in cs.items():
            upd('CURV_STRIP_' + kq, v)
            vv = v + cd if kq == '00' else v
            if 'M_STRIP_' + kq not in sup or vv.upper() > sup['M_STRIP_' + kq]: argmax['M_STRIP_' + kq] = xa
            upd('M_STRIP_' + kq, vv)
        rho_lo = HR.rho_of(xa)                                 # rho decreasing: sup over the box
        L = arb(xb) - arb(xa)
        IR += L * arb(rho_lo.upper()); IR2 += L * arb((rho_lo * rho_lo).upper())
        if rho_lo.upper() > RWs: RWs = arb(rho_lo.upper())
        nS += 1
print(f"stage B: {nB} core boxes, strip {nS} boxes  [{time.time()-t0:.0f}s]", flush=True)
assert xb <= 34.25
RW = RWs * KW; IRW = (IR2 * KW).sqrt()
res = {kq: fup(v) for kq, v in sup.items() if kq != 'AG'}
res['AG_MIN'] = fdn(sup['AG']); res['CB1'] = fup(1 / sup['AG'])
res.update(RA2=fup(RA.sqrt()), RB2=fup(RB.sqrt()), IR=fup(IR), IRW=fup(IRW), RW=fup(RW), CHI1=fup(CHI1), XW_COEF=fup(CHI1 * IR.sqrt()), KW=float(KW.upper()))
res['strip_argmax_x'] = argmax
res['params'] = dict(AMAX=float(AMAX.upper()), WB=float(WB.upper()), NB=NB, NS=NS, GROUP=GROUP, prec=ctx.prec, strip_sigma=[29, 32], strip_x_hi=xb)
res['labels'] = 'all values [cert] (arb, outward) given the certified inputs'
res['runtime_s'] = round(time.time() - t0, 1)
res['args'] = args_block(EXT + [os.path.join(OUT, 'tc_weighted.pkl'), os.path.join(V4B, 'v4cert', 'out', 'constants_v4.json')])
print(json.dumps({k: v for k, v in res.items() if k != 'args'}, indent=1))
print('sha256', dump(os.path.join(OUT, 'e0_strip_drift.json'), res))
