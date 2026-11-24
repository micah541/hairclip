#!/usr/bin/env python3
"""Phase-2a repair of Lemma S(A), stage 4: the assembly in arb, at BOTH (H2) settings, with the dissipation export (R4).
Reuses (by import, unmodified) the certified assembly of the reference chain  base_ibp/rigor/assembly_arb.py  (copy of gonogo/lemmaS_A_ibp):
  sliver_G, convk, frm_tab, infl_tab, window_loop, f4b, U, sup1, max_  and its inputs (J, slopes, alpha_1, sqrt Lambda, abar, beta, J4, sJ4 from
  base_ibp/lemmaS_F10p_M0.52.json args; both reference chains have identical args except the hand-off (H2)).
Replaced inputs:  core/ramp/window hat G (out/weights_rep.json: R3 honest U_j, R4 re-split, R5 lam adds),  ramp tables (out/ramp_rep.json: R1, R2),
  outer rates (out/outer_rep.json: R5, HYPOTHESES H-O1/H-O2 at placeholder values),  forcing integral to Rout = 700 (outer cutoff end).
New parts of E_j (all time-integrated data, added directly; the outer rate is per unit tau and enters as rate/gamma'_j):
  top   : R1, route J (levels 0,1) / IBP + Young against the spare ramp dissipation (level 2)
  cls   : R2, cutoff x jet classes ('pair': Cauchy-Schwarz data or Young;  'phi': sup_w * Phi;  'ibp223': sup_dw/3 * Phi)
  outer : R5
Phi_j = (eps*omega)_min^{-1/2} sqrt(Fmax) (Jf2_j/sqrt(ln2) + sqrt2 Jf2_j^{3/4} JG_j^{1/4}) bounds int int_{ramp} eps*omega|z_2|^3 (the sup-free S33 data
  formula of v3.1, lemmaS_F2.py:204-208 / assembly_arb.py:119-124).
Export (R4), for j = 0, 1 and any tau_a >= tau_1, T >= 0:
  X_j(T) := int_{tau_a}^{tau_a+T} int_{ln 28}^{ln 350} eps*omega z_{j+2}^2 ds dtau <= C0_j + C1_j T,
  C0_j = (E_j^A + D_j)/(r_j a_min^c),  C1_j = (2F_j^2/gamma'_j + coup_j + outer_j + rho_J D_j)/(r_j a_min^c),
  D_j = time-integrated data sources of level j, r_j = exported retained fraction (weights_rep), rho_J = max slope/J ratio.
Output: out/assembly_rep.json, printed tables.  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/assembly_rep.py"""
import os, sys, json, math, re
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, Q_, arb
_cwd = os.getcwd()
import assembly_arb as AS                       # base_ibp/rigor (copy); chdirs to base_ibp, loads the reference inputs
os.chdir(_cwd)
from fractions import Fraction
D = AS.D; sqrt = AS.sqrt; fl = AS.fl
STIR = AS.STIR
SETTINGS = {'tight': ('1e-8', '1e-6', '1e-3'), 'v31': ('2e-7', '1e-4', '0.2')}

import argparse
_ap = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
_ap.add_argument('--weights', default=L.OUT + '/weights_t2.json'); _ap.add_argument('--ramp', default=L.OUT + '/ramp_rep.json')
_ap.add_argument('--outer', default=L.OUT + '/outer_L2.json'); _ap.add_argument('--N2', nargs=3, default=['0.1', '0.2', '5.0'], help="outer data Nbar_k^2 (task 2 r2: R5', o_j = sum_k kappa_jk Nbar_k^2)")
_ap.add_argument('--dprof', default=None, help="fix B0: step-4 v33 p2/out_v33/b0_dprof.json (d = (f - U) - w'): adds d at |alpha| <= alpha_1 to every ramp data "
                 "conversion (verifier J's are for w'), to the verifier fluxes K0, K2, I3 and to the verifier inputs S0V, S1V (both for w')")
PRE, _ = _ap.parse_known_args()
WR = L.load(PRE.weights); RR = L.load(PRE.ramp)['res']; OR = L.load(PRE.outer)
N2 = [A(Q_(x)) for x in PRE.N2]
GH = {j: {int(k): v.upper() for k, v in WR['core'][str(j)]['Ghat'].items()} for j in range(3)}
SLACK = {j: WR['ramp'][str(j)]['slack'].lower() for j in range(3)}
FOOTW = {j: WR['ramp'][str(j)]['footw'].upper() for j in range(3)}; W33 = {j: WR['ramp'][str(j)]['W33'].upper() for j in range(3)}
REXP = [A(Q_(x)) for x in WR['r_export']]
AMIN_R = (1 + A(Q_(5, 100))**2)**-2; AMIN_C = (1 + A(Q_(11, 100))**2)**-2
if 'kappa' in OR['res']['0']:                       # outer_L2.json (R5'): o_j = sum_k kappa_o[j][k] Nbar_k^2
    OUTER = {j: sum((OR['res'][str(j)]['kappa'][k]*N2[k] for k in range(3)), arb(0)).upper() for j in range(3)}
else:                                               # phase-2a outer_rep.json (pointwise placeholders)
    OUTER = {j: OR['res'][str(j)]['rate_at_placeholder'].upper() for j in range(3)}
TOP_MODE = {0: 'J', 1: 'J', 2: 'IBP'}

# ---- fix B0 (REVIEW_v33.md B0).  (A)'s variable is z^A_k = d_s^k (f - U(alpha(tau), y)) on [14, 700] (forcing, core box and outputs are for it);
# the verifier's J_m, sqrt(Lambda) and its inputs FLUX0/FLUX2/S0V/S1V are for w' = f - alpha psi3^' - alpha^2 h2^'.  z^A = (w'-part) + D^k d, with
# d = alpha (psi3^' - 1.3056 y^p) + alpha^2 h2^' (certified tables b0_dprof.json at |alpha| <= alpha_1, plus the step-4 chart allowance 2e-7 at
# level 0).  Every quadratic ramp data term t splits as t <= (1 + eta) t_w + (1 + 1/eta) t_d (per term, eta = sqrt(t_d/t_w) as a point value);
# cubic terms (Phi, int int eps*omega |z_2|^3) as (1 + eta)^2 t_w + (1 + 1/eta)^2 t_d; Cauchy-Schwarz products sqrt(X Y) by AM-GM in the window form.
# t_d is a rate per unit tau (|alpha| <= alpha_1 at every tau): damped (pointwise E^A) contribution t_d/gp_j, window (export) contribution t_d*T.
DP = None
if PRE.dprof:
    _dp = json.load(open(PRE.dprof)); _t = _dp['at_alpha_1']; assert _dp['alpha_1'] == '1.5e-4' and _dp['label'].startswith('certified')
    _ch = A(Q_(repr(_dp['chart_allowance_level0']['value'])))
    _x0 = A(Q_('2249', '1000'))                                              # <= X0 = 2.2490776515, -SIG_TIP: int_14^28 W <= int (y + 2.249)^-4
    _iW = ((14 + _x0)**-3 - (28 + _x0)**-3)/3
    _q = lambda v: A(Q_(repr(v)))
    DW = [_q(v) for v in _t['W_ramp_14_28']]; DW[0] = (DW[0].sqrt() + _ch*_iW.sqrt())**2
    DE = [_q(v) for v in _t['E_28_60']]; DE[0] = (DE[0].sqrt() + _ch*(28*(A(1)/28 - A(1)/60)).sqrt())**2
    DP = dict(path=PRE.dprof, W=DW, E=DE, sup_ramp=[_q(_t['sup_14_28'][0]) + _ch] + [_q(v) for v in _t['sup_14_28'][1:]],
              sup_core=[_q(_t['sup_28_60'][0]) + _ch] + [_q(v) for v in _t['sup_28_60'][1:]], chart=_ch)

class Tr:
    """data term: d = damped value (pointwise E^A), c = window constant, s = window slope (export X_j(T) <= C0 + C1 T)"""
    def __init__(self, d, c, s): self.d, self.c, self.s = d, c, s
    def __add__(self, o): return Tr(self.d + o.d, self.c + o.c, self.s + o.s)
    def sc(self, k): return Tr(k*self.d, k*self.c, k*self.s)
TR0 = Tr(arb(0), arb(0), arb(0))
def _pt(x): return A(Q_(repr(float(x.mid()))))                              # exact point value (any eta > 0 is admissible)
def trq(tw, td, gp, rho):
    """quadratic term: w-part tw (damped = window const, slope rho tw), d-part rate td"""
    if td is None or td.upper() == 0: return Tr(tw, tw, rho*tw)
    e = _pt((td/gp/tw).sqrt()); return Tr((1 + e)*tw + (1 + 1/e)*td/gp, (1 + e)*tw, (1 + e)*rho*tw + (1 + 1/e)*td)
def trc(tw, td, gp, rho):
    """cubic term ((a + b)^3 <= (1 + eta)^2 a^3 + (1 + 1/eta)^2 b^3)"""
    if td is None or td.upper() == 0: return Tr(tw, tw, rho*tw)
    e = _pt(((td/gp/tw).log()/3).exp()); return Tr((1 + e)**2*tw + (1 + 1/e)**2*td/gp, (1 + e)**2*tw, (1 + e)**2*rho*tw + (1 + 1/e)**2*td)
def trp(X, Y):
    """Cauchy-Schwarz product: damped sqrt(Xd Yd); window sqrt((Xc + Xs T)(Yc + Ys T)) <= (l (Xc + Xs T) + (Yc + Ys T)/l)/2, l = sqrt(Yc/Xc)"""
    l = _pt((Y.c/X.c).sqrt()); return Tr((X.d*Y.d).sqrt(), (l*X.c + Y.c/l)/2, (l*X.s + Y.s/l)/2)
def tabs(d): return {int(m): v.upper() for m, v in d.items()}

def rat_floor(x, den=10**4):
    return A(Q_(math.floor(float(x.lower())*den), den))

def assemble(HO, G, Rout=700):
    y0, y1, YO = arb(14), arb(28), arb(60); Ro = arb(Rout)
    ABAR, BETA, AL, P, GAM = AS.ABAR, AS.BETA, AS.AL, AS.P, AS.GAM
    cbar = D('0.25') + ABAR + ABAR/y0; cmin = D('0.25') - ABAR - ABAR/y0
    eps_core = y1**-4; ee = 2*P - GAM; Iint = y1**GAM*((Ro**ee - y0**ee)/ee)
    r0 = D('1.4')*(BETA + D('1.02')*AL*ABAR)
    INFL = {j: AS.infl_tab(j) for j in range(3)}; FR = AS.frm_tab(400); E = []; parts = []; gps = []; info = []
    frmin = AS.eps(arb(1))*AS.om(arb(1))
    rho = AS.max_(*([AS.SL[m]/AS.J[m] for m in range(3)] + [AS.SJ4/AS.J4]))      # window slope/constant ratio of the verifier data (as export)
    cvd = (lambda k, t: AS.convk(k, t, {m: DP['W'][m] for m in range(4)}, None)) if DP else (lambda k, t: None)
    for j in range(3):
        gp = GAM*cmin - (abs(1 - 2*j) + 1)*ABAR/y0 - eps_core*GH[j][j]; assert gp.lower() > 0; gps.append(gp)
        J4e = AS.J4 + AS.SJ4/gp; Je = [AS.J[m] + AS.SL[m]/gp for m in range(3)]; Je.append(sqrt(Je[2]*J4e) + AS.C1W*Je[2])
        cvw = lambda k, t: AS.convk(k, t, {m: Je[m] for m in range(4)}, G)
        cv = lambda k, t: trq(cvw(k, t), cvd(k, t), gp, rho)
        infl = cv(j, INFL[j]).sc(cbar)
        rampc = TR0
        for k in range(j + 1): rampc = rampc + cv(k, tabs(RR[str(j)]['q'][str(k)]))
        # ---- Phi (sup-free int int eps*omega |z_2|^3 on the ramp); w-part as before, d-part sup|D^2 d| int eps*omega (D^2 d)^2 per unit tau
        Jm = {m: Je[m] for m in range(4)}
        Jf2 = AS.convk(2, FR, Jm, G); JG = 2*(AS.convk(3, FR, Jm, G) + D('6.25')*AS.convk(2, FR, Jm, G))
        L2 = AS.SQ**2; V0s = L2/D('0.08'); V2s = L2/(D('0.08')*D('3e-5')); V1s = sqrt(V0s*V2s) + AS.C1W*V0s
        Fmax = AS.convk(2, FR, {0: V0s, 1: V1s, 2: V2s}, G)
        Phiw = sqrt(Fmax)*(Jf2/sqrt(AS.LN2) + sqrt(arb(2))*Jf2**D('0.75')*JG**D('0.25'))/sqrt(frmin)
        Phi = trc(Phiw, DP['sup_ramp'][2]*cvd(2, FR) if DP else None, gp, rho)
        Jf2T = trq(Jf2, cvd(2, FR), gp, rho)
        # ---- young users and data alternatives
        users = []                                    # (name, A_i = cost*delta (Tr), data alternative (Tr) or None)
        R = RR[str(j)]
        if TOP_MODE[j] == 'IBP':
            users.append(('top', cv(j, tabs(R['topY'])).sc(1/(4*AMIN_R)), None))
        clist = []
        for e in R['classes_J'] + (R['classes_IBP'] if TOP_MODE[j] == 'IBP' else []):
            I = [int(x) for x in e['I']]; h = e['handler']; c = dict(I=I, tag=e['tag'], handler=h)
            if h == 'phi': c['cost'] = Phi.sc(e['sup_w'].upper())
            elif h == 'ibp223': c['cost'] = Phi.sc(e['sup_dw'].upper()/3)
            else:
                x, y = I; dat = None
                if all(str(k) in e['data'] for k in (x, y)):
                    dx = cv(x, tabs(e['data'][str(x)])); dat = dx if x == y else trp(dx, cv(y, tabs(e['data'][str(y)])))
                if e.get('young') is not None: users.append((len(clist), cv(x, tabs(e['young'])).sc(1/(4*AMIN_R)), dat)); c['cost'] = None
                else:
                    assert dat is not None, ('class without a route', j, I); c['cost'] = dat; c['route'] = 'data'
            clist.append(c)
        # delta allocation over the young users, proportional to sqrt(A_i) (damped values), on the 1e-4 rational grid, total <= ramp slack
        dsum = arb(0); deltas = {}
        if users:
            sq = [sqrt(u[1].d) for u in users]; tot = sum(sq, arb(0))
            for u, s_ in zip(users, sq):
                d_ = rat_floor(SLACK[j]*D("0.999")*s_/tot); deltas[u[0]] = d_; dsum += d_
            assert (dsum - SLACK[j]).upper() <= 0
        top = TR0; topinfo = {}
        if TOP_MODE[j] == 'J':
            top = cv(j + 1, tabs(R['topJ'])); topinfo = dict(route='J', J=top.d)
        for u in users:
            yc = u[1].sc(1/deltas[u[0]]) if deltas[u[0]] > 0 else None
            if u[0] == 'top':
                t2 = cv(j, tabs(R['top2'])); top = yc + t2; topinfo = dict(route='IBP', young=yc.d, cpp=t2.d, delta=deltas['top'], J_route_for_comparison=cv(j + 1, tabs(R['topJ'])).d)
            else:
                c = clist[u[0]]; dat = u[2]
                if yc is not None and (dat is None or yc.d.upper() <= dat.d.upper()): c['cost'] = yc; c['route'] = 'young'; c['delta'] = deltas[u[0]]
                else: assert dat is not None; c['cost'] = dat; c['route'] = 'data'
        cls_top = TR0; cls_pairs = TR0
        for c in clist:
            if c['tag'] == 'top': cls_top = cls_top + c['cost']
            elif c['tag'] == 'pairs': cls_pairs = cls_pairs + c['cost']
        if TOP_MODE[j] == 'IBP': top = top + cls_top
        cls = cls_pairs
        # ---- extra_outer of v3.1 (foot + sup-free S33 on the ramp), with the R3 weights
        foot = Jf2T.sc(FOOTW[j]); S33d = TR0
        if W33[j] > 0: S33d = Phi.sc(W33[j]*AS.CE/3)
        extra = foot + S33d
        coup = sum((GH[j].get(k, arb(0))*eps_core*E[k] for k in range(j)), arb(0))
        F2 = (r0*P**j)**2*Iint
        outer = OUTER[j]/gp
        data = infl + rampc + top + cls + extra
        Ej = HO[j] + data.d + 4*F2/gp**2 + 2*coup/gp + outer
        E.append(Ej)
        parts.append(dict(gamma_prime=gp, inflow=infl.d, ramp_comm=rampc.d, forcing=4*F2/gp**2, coupling=2*coup/gp, handoff=HO[j], extra=extra.d,
                          top=top.d, cls=cls.d, outer=outer, foot=foot.d, S33=S33d.d))
        for c in clist: c['cost'] = c['cost'].d
        dw = dict(Phi_w=Phiw, Phi_d_rate=(DP['sup_ramp'][2]*cvd(2, FR)) if DP else arb(0), Jf2_d_rate=cvd(2, FR) if DP else arb(0))
        info.append(dict(Phi=Phi.d, Jf2=Jf2T.d, JG=JG, Fmax=Fmax, top=topinfo, classes=clist, deltas={str(k): v for k, v in deltas.items()}, F2=F2, coup_rate=coup,
                         D=data.d, D_window_const=data.c, D_window_slope=data.s, rho_w=rho, b0=dw))
    conv = (YO/y1)**AS.GAM; ell = (YO/y1).log(); Eu = [conv*e for e in E]
    S0 = sqrt(Eu[0]/ell + 2*sqrt(Eu[0]*Eu[1])); S1 = sqrt(Eu[1]/ell + 2*sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts, gp=gps, info=info, cmin=cmin, cbar=cbar, r0=r0, Iint=Iint)

def export(R):
    rhoJ = max([AS.SL[m]/AS.J[m] for m in range(3)] + [AS.SJ4/AS.J4], key=lambda x: float(x.upper()))
    out = {}
    for j in (0, 1):
        gp = R['gp'][j]; inf = R['info'][j]; Dj = inf['D']
        C0 = (R['E'][j] + Dj)/(REXP[j]*AMIN_C)
        C1 = (2*inf['F2']/gp + inf['coup_rate'] + OUTER[j] + inf['D_window_slope'])/(REXP[j]*AMIN_C)   # B0: = rhoJ*Dj without d; with d the d-rates
        assert (inf['D_window_const'] - Dj).upper() <= 0                                                   # C0 uses the damped Dj >= window constant
        out[j] = dict(C0=C0, C1=C1, r=REXP[j], amin_c=AMIN_C, D=Dj, E=R['E'][j], rhoJ=rhoJ, D_window_slope=inf['D_window_slope'],
                      statement='int_{tau_a}^{tau_a+T} int_{ln28}^{ln350} y^-4 (28/y) z_%d^2 ds dtau <= C0 + C1*T  (z_k = d_s^k v, v = f - U, s = ln y)' % (j + 2))
    return out

def S01w(R):
    """core sup|w'|, sup|y w''| (verifier variable) = (A)'s S0v, S1v plus sup|d|, sup|D d| on [28, 60] (B0); unchanged without --dprof"""
    return (R['S0'] + DP['sup_core'][0], R['S1'] + DP['sup_core'][1]) if DP else (R['S0'], R['S1'])

def margins(R, W_, F4):
    # M14 (M3): verifier ramp sups read from fin12f4 (alpha_1 = 1.5e-4 row), rounded UP to 4 significant digits, asserted
    import decimal as _dc
    _V = json.load(open('/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'))['1.5e-4']
    def _up4(x):
        d = _dc.Decimal(repr(float(x))); q = _dc.Decimal(1).scaleb(d.adjusted() - 3)
        r = d.quantize(q, rounding=_dc.ROUND_CEILING); assert r >= _dc.Decimal(repr(float(x))); return str(r)
    VM1 = _up4(_V['M1']); VI1 = _up4(_V['sI1']); print('M14 verifier ramp sups (rounded up):', 'M1', VM1, 'sI1', VI1, flush=True)
    if DP:   # B0: ramp |z^A_1| <= 0.4328 (verifier, w') + sup|D d|; |f| <= 0.00116 (w') + sup|d| + U(28); y f_y = D w' + D d + p U
        assert (DP['sup_ramp'][1] - D('2e-5')).upper() <= 0      # the typed 2e-5 of v3.1 covers sup_ramp |D d| (certified)
        rM1 = D(VM1) + AS.P*AS.U(arb(28)) + DP['sup_ramp'][1]; rM0 = D(VI1) + AS.U(arb(28)) + DP['sup_ramp'][0]; rz1 = D(VM1) + DP['sup_ramp'][1]
    else:
        rM1 = D(VM1) + AS.P*AS.U(arb(28)) + D('2e-5'); rM0 = D(VI1) + AS.U(arb(28)); rz1 = D(VM1) + D('2e-5')
    wmax = lambda f: AS.max_(*[f(w) for w in W_])
    S0w, S1w = S01w(R)
    return {'S0v': D('0.02')/S0w, 'S0v_vs_0.0153': D('0.0153')/S0w, 'S1v_vs_0.6': D('0.6')/S1w, 'M0core': D('0.11')/F4['S0f'], 'M1core': D('0.52')/F4['S1f'],
            'ramp M1': D('0.52')/rM1, 'ramp z1': D('0.52')/rz1, 'ramp M0': D('0.05')/rM0,
            'Ghat22': D('153664')/GH[2][2], 'win M1': min((w['Ubox'][1]/w['der'][1] for w in W_), key=lambda x: float(x.lower())),
            'win M0': min((w['Ubox'][0]/w['der'][0] for w in W_), key=lambda x: float(x.lower())), 'win z2': D('0.5')/wmax(lambda w: w['der'][2]),
            '|f|<=1/2': D('0.5')/wmax(lambda w: w['der'][0]), 'K2': D('0.12')/wmax(lambda w: w['K'][0]), 'K3': D('0.12')/wmax(lambda w: w['K'][1])}

def claims(R, W_, F4):
    chk = []
    def claim(name, quoted, rel, val):
        ok = (D(quoted) - val).lower() >= 0 if rel == '<=' else (val - D(quoted)).lower() >= 0
        chk.append((name, quoted, rel, float(val.upper() if rel == '<=' else val.lower()), ok))
    for j, q in enumerate(('3.12e-7', '1.81e-4', '0.248')): claim('E_%d (v3.1 page 1)' % j, q, '<=', R['E'][j])
    for j, q in enumerate(('0.246', '0.246', '0.143')): claim("gamma'_%d (v3.1 page 1)" % j, q, '>=', R['gp'][j])
    S0w, S1w = S01w(R)
    claim('core sup|v| <= 0.02 (H-B2)', '0.02', '<=', S0w); claim('S0_v <= S0V = 0.0153 (verifier input)', '0.0153', '<=', S0w); claim('S1_v <= S1V = 0.6 (verifier input)', '0.6', '<=', S1w)
    claim('core M0 <= 0.11', '0.11', '<=', F4['S0f']); claim('core M1 <= 0.52', '0.52', '<=', F4['S1f'])
    for j, q in enumerate(('3.15e-7', '1.83e-4', '0.250')): claim('E_B%d (y_E=50) (v3.1 page 1)' % j, q, '<=', F4['B'][50]['EB'][j])
    claim("gamma'_B >= 0.2394 (v3.1 page 1)", '0.2394', '>=', min(F4['B'][50]['gp'], key=lambda x: float(x.lower())))
    claim('S0_B (y_E=50) <= S0V=0.0153', '0.0153', '<=', F4['B'][50]['S0']); claim('S1_B (y_E=50) <= S1V=0.6', '0.6', '<=', F4['B'][50]['S1'])
    wmax = lambda f: AS.max_(*[f(w) for w in W_])
    claim('windows |f| <= 1/2', '0.5', '<=', wmax(lambda w: w['der'][0])); claim('windows |z2| <= 0.5', '0.5', '<=', wmax(lambda w: w['der'][2]))
    for w in W_:
        claim('  window %s: |f| <= box M0' % w['yE'], '%.12f' % fl(w['Ubox'][0]), '<=', w['der'][0]); claim('  window %s: y|f\'| <= box M1' % w['yE'], '%.12f' % fl(w['Ubox'][1]), '<=', w['der'][1])
    claim("windows gamma' >= 0.181 (v3.1 page 1)", '0.181', '>=', min((g for w in W_ for g in w['gp']), key=lambda x: float(x.lower())))
    claim('K2 <= 0.12', '0.12', '<=', wmax(lambda w: w['K'][0])); claim('K3 <= 0.12', '0.12', '<=', wmax(lambda w: w['K'][1]))
    claim('K0 <= 6.093e-21 (1e-3 x verifier placeholder)', '6.093e-21', '<=', F4['K0']); claim('K2 flux <= 1.5445e-20', '1.5445e-20', '<=', F4['K2'])
    claim('K0 <= FLUX0 = 1.64e-20 (verifier constants.py:65)', '1.64e-20', '<=', F4['K0']); claim('K2 flux <= FLUX2 = 8.0e-20 (verifier constants.py:65)', '8.0e-20', '<=', F4['K2'])
    claim('K0 (v3.1 page 6.83e-22)', '6.83e-22', '<=', F4['K0']); claim('K2 flux (v3.1 page 3.26e-21)', '3.26e-21', '<=', F4['K2'])
    claim('I3 const (v3.1 page 1.84e-21)', '1.84e-21', '<=', F4['I3c']); claim('I3 slope (v3.1 page 1.39e-27)', '1.39e-27', '<=', F4['I3s'])
    claim('Ghat_22 core < cmin/eps_core = 153664', '153664', '<=', GH[2][2])
    return chk

def f4b_honest(R, F4):
    """I_3 with the level-1 budget including every new data part and the RETAINED fraction r_1 (v3.1 used 2 a_min, assembly_arb.py:175)"""
    pa = R['parts'][1]; c = json.load(open(L.BASE + '/lemmaS_cflux.json')); k2 = [AS.UP(x) for x in c['kappa_I2']]; W = AS.UP(c['W57_upper'])
    conv = (arb(60)/arb(28))**AS.GAM; Eu = [conv*x for x in R['E']]
    cF = D('1.4')*(AS.BETA + D('1.02')*AS.AL*AS.ABAR)/AS.SQ*AS.P*sqrt(R['Iint'])
    bud_c = pa['handoff'] + R['info'][1]['D'] + 2*cF*sqrt(R['E'][1])*F4['PHI_C']
    bud_s = 2*cF*sqrt(R['E'][1])*F4['PHI_S'] + R['info'][1]['coup_rate'] + OUTER[1] + 2*R['info'][1]['F2']/R['gp'][1]
    yw = D('57.75'); fac = 1/((AMIN_C*REXP[1]*yw**-4*(28/yw)))*(yw/28)**AS.GAM
    IE = [R['E'][j]/R['gp'][j] for j in range(3)]
    N3c = W*3*(bud_c*fac + conv*(9*IE[2] + 4*IE[1]))/D('54.75')**5; N3s = W*3*(bud_s*fac)/D('54.75')**5
    return dict(I3c=k2[3]*N3c, I3s=k2[3]*N3s)

def ref_margins(path):
    try:
        txt = open(path).read(); m = re.search(r'margins \(arb, certified lower bounds\): (.*)', txt); return m.group(1) if m else None
    except Exception: return None

if __name__ == '__main__':
    ap = argparse.ArgumentParser(parents=[_ap], allow_abbrev=False)
    ap.add_argument('--betas', nargs='*', default=['5/4', '2'])
    ap.add_argument('--ic', default=None, help='task 2 / A4: out/ic_rep.json; E^A(tau_1) per setting from its E_tau1[<setting>/<h3key>]')
    ap.add_argument('--h3key', default='h3_v32')
    ap.add_argument('--hjw', default=None, help='task 2 / A1: out/hjw_rep.json; J_m -> J_m (1 + kappa_J), m = 0,1,2,4')
    ap.add_argument('--out', default=L.OUT + '/assembly_rep.json')
    CL = ap.parse_args(); EXTRA_IN = []
    ICE = None
    if CL.ic:
        ICE = L.load(CL.ic)['E_tau1']; EXTRA_IN.append(CL.ic)
    if CL.hjw:
        kJ = A(Q_(repr(json.load(open(CL.hjw))['kappa_J']))); EXTRA_IN.append(CL.hjw)
        AS.J = [x*(1 + kJ) for x in AS.J]; AS.J4 = AS.J4*(1 + kJ)
        print('H-J-W fallback: J_m scaled by 1 + kappa_J =', 1 + float(kJ.upper()))
    G, vmax, ymax, supg = AS.sliver_G()
    out = {}; refs = {'tight': (L.BASE + '/rigor/assembly_arb.json', L.BASE + '/rigor_assembly_arb.log'), 'v31': (L.BASE31 + '/rigor/assembly_arb.json', L.BASE31 + '/rigor_assembly_arb.log')}
    V31 = ('/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.json', '/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.out')
    AS.WE = dict(AS.WE); AS.WE['windows'] = WR['windows']; AS.GHAT = GH
    for name, ho in SETTINGS.items():
        HO = [D(x) for x in ho]; AS.HO = HO          # AS.HO: hand-off shells (H2) only (H6 diagnostic of assembly_arb)
        if ICE is not None: HO = [x.upper() for x in ICE['%s/%s' % (name, CL.h3key)]['E']]; print('E^A(tau_1) from', CL.ic, name, CL.h3key, [float(x) for x in HO])
        R = assemble(HO, G)
        os.chdir(L.BASE); W_ = AS.window_loop(); F4 = AS.f4b(R); os.chdir(_cwd)
        F4h = f4b_honest(R, F4); F4['I3c'] = AS.max_(F4['I3c'], F4h['I3c']); F4['I3s'] = AS.max_(F4['I3s'], F4h['I3s'])
        KD = None
        if DP:   # B0: the verifier's fluxes are for w' = z^A - d: K_i <= (sqrt(K_i(E)) + sqrt(K_i(E_d)))^2 (K_i = sum_l k_l N_l, k_l >= 0, N_l linear in E)
            cfl = json.load(open(L.BASE + '/lemmaS_cflux.json')); assert all(float(x) >= 0 for x in cfl['kappa_I0'] + cfl['kappa_I2'])
            os.chdir(L.BASE); F4d = AS.f4b(dict(R, E=DP['E'][:3])); os.chdir(_cwd)
            KD = dict(K0=F4d['K0'], K2=F4d['K2'], K0_w=F4['K0'], K2_w=F4['K2'])
            F4['K0'] = (sqrt(F4['K0']) + sqrt(F4d['K0']))**2; F4['K2'] = (sqrt(F4['K2']) + sqrt(F4d['K2']))**2
            k2 = [AS.UP(x) for x in cfl['kappa_I2']]; Wc = AS.UP(cfl['W57_upper']); cvc = (arb(60)/arb(28))**AS.GAM
            I3d = k2[3]*Wc*3*cvc*(DP['E'][3] + 9*DP['E'][2] + 4*DP['E'][1])/D('54.75')**5; et = A(Q_(1, 8))
            KD.update(I3_d_rate=I3d, I3c_w=F4['I3c'], I3s_w=F4['I3s'], I3_eta=et)
            F4['I3c'] = (1 + et)*F4['I3c']; F4['I3s'] = (1 + et)*F4['I3s'] + (1 + 1/et)*I3d
        M = margins(R, W_, F4); C = claims(R, W_, F4); X = export(R)
        # A4b: flux constants at the bootstrap level beta*E (K0, K2 are linear in E: recomputed by f4b on beta*E, asserted)
        FLX = dict(K0_at_E=F4['K0'], K2_at_E=F4['K2'], FLUX0=D('1.64e-20'), FLUX2=D('8.0e-20'), FLUX0_needed=4*F4['K0'], FLUX2_needed=4*F4['K2'], beta={})
        if KD: FLX['B0_d'] = KD
        for bq in CL.betas:
            b = A(Q_(bq)); Rb = dict(R); Rb['E'] = [b*e for e in R['E']]
            os.chdir(L.BASE); F4b = AS.f4b(Rb); os.chdir(_cwd)
            K0w, K2w = (KD['K0_w'], KD['K2_w']) if KD else (F4['K0'], F4['K2'])
            assert (F4b['K0'] - b*K0w).upper() <= 1e-30 and (F4b['K2'] - b*K2w).upper() <= 1e-30
            if KD: F4b['K0'] = (sqrt(F4b['K0']) + sqrt(KD['K0']))**2; F4b['K2'] = (sqrt(F4b['K2']) + sqrt(KD['K2']))**2
            FLX['beta'][bq] = dict(K0=F4b['K0'], K2=F4b['K2'], margin_FLUX0=D('1.64e-20')/F4b['K0'], margin_FLUX2=D('8.0e-20')/F4b['K2'])
            if bq == '2': FLX['FLUX0_needed'] = 2*F4b['K0']; FLX['FLUX2_needed'] = 2*F4b['K2']     # x2 at beta = 2 (= 4 K(E) without d)
        print('flux: K0(E) = %.4g, K2(E) = %.4g; ' % (fl(F4['K0']), fl(F4['K2'])) + '; '.join('beta %s: FLUX0/K0 x%.3f FLUX2/K2 x%.3f' % (bq, float(v['margin_FLUX0'].lower()), float(v['margin_FLUX2'].lower())) for bq, v in FLX['beta'].items()))
        ref = json.load(open(refs[name][0])); v31 = json.load(open(V31[0]))
        print('\n==================== setting %s: (H2) = %s ====================' % (name, ho))
        print('%-10s %-14s %-14s %-14s | %s' % ('', 'repaired', 'lemmaS_A_ibp' + ('' if name == 'tight' else '_v31H2'), 'v3.1', 'parts of the repaired E_j'))
        for j in range(3):
            p = R['parts'][j]
            print('E_%d       %-14.8g %-14.8g %-14.8g | handoff %.3g inflow %.3g ramp %.3g TOP %.3g CLS %.3g extra %.3g forcing %.2g coupling %.2g OUTER %.2g' % (
                j, fl(R['E'][j]), ref['E'][j], v31['E'][j], *(fl(p[k]) for k in ('handoff', 'inflow', 'ramp_comm', 'top', 'cls', 'extra', 'forcing', 'coupling', 'outer'))))
        for j in range(3): print("gamma'_%d  %-14.8g %-14.8g %-14.8g" % (j, float(R['gp'][j].lower()), ref['gp'][j], v31['gp'][j]))
        for j in range(3):
            inf = R['info'][j]; print('   level %d: Phi %.4g; top route %s %s; deltas %s' % (j, fl(inf['Phi']), inf['top'].get('route'), {k: '%.3g' % fl(v) for k, v in inf['top'].items() if isinstance(v, arb)}, {k: '%.4g' % fl(v) for k, v in inf['deltas'].items()}))
            for c in inf['classes']: print('       class %-4s %-10s %-7s %-6s %.3g' % (c['tag'], c['I'], c['handler'], c.get('route', ''), fl(c['cost'])))
        print('margins (arb lower):', {k: round(float(v.lower()), 4) for k, v in M.items()})
        print('   ref  :', ref_margins(refs[name][1])); print('   v3.1 :', ref_margins(V31[1]))
        print('export (R4):', {j: dict(C0='%.5g' % fl(x['C0']), C1='%.4g' % fl(x['C1']), r=str(WR['r_export'][j])) for j, x in X.items()})
        bad = [c for c in C if not c[4]]
        print('claims failing:', [(c[0], c[1], '%.6g' % c[3]) for c in bad])
        out[name] = dict(flux=FLX, OUTER={j: OUTER[j] for j in range(3)}, E=R['E'], gp=R['gp'], S0=R['S0'], S1=R['S1'], parts=R['parts'], info=R['info'], margins=M, claims=C, claims_failing=bad, export=X,
                         F4={k: v for k, v in F4.items() if k != 'B'}, B50=F4['B'][50], windows=[dict(yE=w['yE'], gp=w['gp'], der=w['der'], K=w['K'], Ubox=w['Ubox']) for w in W_],
                         ref=dict(E=ref['E'], gp=ref['gp'], margins=ref_margins(refs[name][1])), v31=dict(E=v31['E'], gp=v31['gp'], margins=ref_margins(V31[1])))
    inputs = [PRE.weights, PRE.ramp, PRE.outer, L.BASE + '/rigor/assembly_arb.py', L.BASE + '/lemmaS_F10p_M0.52.json',
              L.BASE + '/rigor/cap_F10p.json', L.BASE + '/lemmaS_cflux.json', '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json', os.path.abspath(__file__),
              refs['tight'][0], refs['v31'][0], V31[0]]
    b0 = None
    if DP:
        EXTRA_IN.append(PRE.dprof)
        b0 = dict(dprof=PRE.dprof, rule='data terms: (1+eta) t_w + (1+1/eta) t_d (quadratic), (1+eta)^2 t_w + (1+1/eta)^2 t_d (cubic); fluxes (sqrt K(E) + sqrt K(E_d))^2; S0V/S1V: + sup|d|, sup|D d| on [28, 60]',
                  W_ramp=DP['W'], E_core=DP['E'], sup_ramp=DP['sup_ramp'], sup_core=DP['sup_core'], chart_level0=DP['chart'])
    L.dump(dict(res=out, top_mode=TOP_MODE, outer=dict(file=PRE.outer, N2=PRE.N2, hyp=OR.get('hyp')), ic=dict(file=CL.ic, key=CL.h3key), hjw=CL.hjw, b0=b0), CL.out, inputs + EXTRA_IN)
    print('\nwrote', CL.out)
