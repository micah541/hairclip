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

WR = L.load(L.OUT + '/weights_rep.json'); RR = L.load(L.OUT + '/ramp_rep.json')['res']; OR = L.load(L.OUT + '/outer_rep.json')
GH = {j: {int(k): v.upper() for k, v in WR['core'][str(j)]['Ghat'].items()} for j in range(3)}
SLACK = {j: WR['ramp'][str(j)]['slack'].lower() for j in range(3)}
FOOTW = {j: WR['ramp'][str(j)]['footw'].upper() for j in range(3)}; W33 = {j: WR['ramp'][str(j)]['W33'].upper() for j in range(3)}
REXP = [A(Q_(x)) for x in WR['r_export']]
AMIN_R = (1 + A(Q_(5, 100))**2)**-2; AMIN_C = (1 + A(Q_(11, 100))**2)**-2
OUTER = {j: OR['res'][str(j)]['rate_at_placeholder'].upper() for j in range(3)}
TOP_MODE = {0: 'J', 1: 'J', 2: 'IBP'}
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
    for j in range(3):
        gp = GAM*cmin - (abs(1 - 2*j) + 1)*ABAR/y0 - eps_core*GH[j][j]; assert gp.lower() > 0; gps.append(gp)
        J4e = AS.J4 + AS.SJ4/gp; Je = [AS.J[m] + AS.SL[m]/gp for m in range(3)]; Je.append(sqrt(Je[2]*J4e) + AS.C1W*Je[2])
        cv = lambda k, t: AS.convk(k, t, {m: Je[m] for m in range(4)}, G)
        infl = cbar*cv(j, INFL[j])
        rampc = sum((cv(k, tabs(RR[str(j)]['q'][str(k)])) for k in range(j + 1)), arb(0))
        # ---- Phi (sup-free int int eps*omega |z_2|^3 on the ramp)
        Jm = {m: Je[m] for m in range(4)}
        Jf2 = AS.convk(2, FR, Jm, G); JG = 2*(AS.convk(3, FR, Jm, G) + D('6.25')*AS.convk(2, FR, Jm, G))
        L2 = AS.SQ**2; V0s = L2/D('0.08'); V2s = L2/(D('0.08')*D('3e-5')); V1s = sqrt(V0s*V2s) + AS.C1W*V0s
        Fmax = AS.convk(2, FR, {0: V0s, 1: V1s, 2: V2s}, G)
        Phi = sqrt(Fmax)*(Jf2/sqrt(AS.LN2) + sqrt(arb(2))*Jf2**D('0.75')*JG**D('0.25'))/sqrt(frmin)
        # ---- young users and data alternatives
        users = []                                    # (name, A_i = cost*delta, data alternative or None)
        R = RR[str(j)]
        if TOP_MODE[j] == 'IBP':
            users.append(('top', cv(j, tabs(R['topY']))/(4*AMIN_R), None))
        clist = []
        for e in R['classes_J'] + (R['classes_IBP'] if TOP_MODE[j] == 'IBP' else []):
            I = [int(x) for x in e['I']]; h = e['handler']; c = dict(I=I, tag=e['tag'], handler=h)
            if h == 'phi': c['cost'] = e['sup_w'].upper()*Phi
            elif h == 'ibp223': c['cost'] = e['sup_dw'].upper()*Phi/3
            else:
                x, y = I; dat = None
                if all(str(k) in e['data'] for k in (x, y)):
                    dx = cv(x, tabs(e['data'][str(x)])); dat = dx if x == y else sqrt(dx*cv(y, tabs(e['data'][str(y)])))
                if e.get('young') is not None: users.append((len(clist), cv(x, tabs(e['young']))/(4*AMIN_R), dat)); c['cost'] = None
                else:
                    assert dat is not None, ('class without a route', j, I); c['cost'] = dat; c['route'] = 'data'
            clist.append(c)
        # delta allocation over the young users, proportional to sqrt(A_i), on the 1e-4 rational grid, total <= ramp slack
        dsum = arb(0); deltas = {}
        if users:
            sq = [sqrt(u[1]) for u in users]; tot = sum(sq, arb(0))
            for u, s_ in zip(users, sq):
                d_ = rat_floor(SLACK[j]*D("0.999")*s_/tot); deltas[u[0]] = d_; dsum += d_
            assert (dsum - SLACK[j]).upper() <= 0
        top = arb(0); topinfo = {}
        if TOP_MODE[j] == 'J':
            top = cv(j + 1, tabs(R['topJ'])); topinfo = dict(route='J', J=top)
        for u in users:
            yc = u[1]/deltas[u[0]] if deltas[u[0]] > 0 else None
            if u[0] == 'top':
                t2 = cv(j, tabs(R['top2'])); top = yc + t2; topinfo = dict(route='IBP', young=yc, cpp=t2, delta=deltas['top'], J_route_for_comparison=cv(j + 1, tabs(R['topJ'])))
            else:
                c = clist[u[0]]; dat = u[2]
                if yc is not None and (dat is None or yc.upper() <= dat.upper()): c['cost'] = yc; c['route'] = 'young'; c['delta'] = deltas[u[0]]
                else: assert dat is not None; c['cost'] = dat; c['route'] = 'data'
        cls_top = sum((c['cost'] for c in clist if c['tag'] == 'top'), arb(0)); cls_pairs = sum((c['cost'] for c in clist if c['tag'] == 'pairs'), arb(0))
        if TOP_MODE[j] == 'IBP': top = top + cls_top
        cls = cls_pairs
        # ---- extra_outer of v3.1 (foot + sup-free S33 on the ramp), with the R3 weights
        foot = FOOTW[j]*Jf2; S33d = arb(0)
        if W33[j] > 0: S33d = W33[j]*AS.CE/3*Phi
        extra = foot + S33d
        coup = sum((GH[j].get(k, arb(0))*eps_core*E[k] for k in range(j)), arb(0))
        F2 = (r0*P**j)**2*Iint
        outer = OUTER[j]/gp
        Ej = HO[j] + infl + rampc + top + cls + extra + 4*F2/gp**2 + 2*coup/gp + outer
        E.append(Ej)
        parts.append(dict(gamma_prime=gp, inflow=infl, ramp_comm=rampc, forcing=4*F2/gp**2, coupling=2*coup/gp, handoff=HO[j], extra=extra,
                          top=top, cls=cls, outer=outer, foot=foot, S33=S33d))
        info.append(dict(Phi=Phi, Jf2=Jf2, JG=JG, Fmax=Fmax, top=topinfo, classes=clist, deltas={str(k): v for k, v in deltas.items()}, F2=F2, coup_rate=coup,
                         D=infl + rampc + top + cls + extra))
    conv = (YO/y1)**AS.GAM; ell = (YO/y1).log(); Eu = [conv*e for e in E]
    S0 = sqrt(Eu[0]/ell + 2*sqrt(Eu[0]*Eu[1])); S1 = sqrt(Eu[1]/ell + 2*sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts, gp=gps, info=info, cmin=cmin, cbar=cbar, r0=r0, Iint=Iint)

def export(R):
    rhoJ = max([AS.SL[m]/AS.J[m] for m in range(3)] + [AS.SJ4/AS.J4], key=lambda x: float(x.upper()))
    out = {}
    for j in (0, 1):
        gp = R['gp'][j]; inf = R['info'][j]; Dj = inf['D']
        C0 = (R['E'][j] + Dj)/(REXP[j]*AMIN_C)
        C1 = (2*inf['F2']/gp + inf['coup_rate'] + OUTER[j] + rhoJ*Dj)/(REXP[j]*AMIN_C)
        out[j] = dict(C0=C0, C1=C1, r=REXP[j], amin_c=AMIN_C, D=Dj, E=R['E'][j], rhoJ=rhoJ,
                      statement='int_{tau_a}^{tau_a+T} int_{ln28}^{ln350} y^-4 (28/y) z_%d^2 ds dtau <= C0 + C1*T  (z_k = d_s^k v, v = f - U, s = ln y)' % (j + 2))
    return out

def margins(R, W_, F4):
    rM1 = D('0.4328') + AS.P*AS.U(arb(28)) + D('2e-5'); rM0 = D('0.00116') + AS.U(arb(28))
    wmax = lambda f: AS.max_(*[f(w) for w in W_])
    return {'S0v': D('0.02')/R['S0'], 'M0core': D('0.11')/F4['S0f'], 'M1core': D('0.52')/F4['S1f'], 'ramp M1': D('0.9')/rM1, 'ramp M0': D('0.05')/rM0,
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
    claim('core sup|v| <= 0.02 (H-B2)', '0.02', '<=', R['S0']); claim('S0_v <= S0V = 0.0153 (verifier input)', '0.0153', '<=', R['S0']); claim('S1_v <= S1V = 0.6 (verifier input)', '0.6', '<=', R['S1'])
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
    claim('K0 <= FLUX0 = 4.1e-21 (verifier constants.py:65)', '4.1e-21', '<=', F4['K0']); claim('K2 flux <= FLUX2 = 2.0e-20 (verifier constants.py:65)', '2.0e-20', '<=', F4['K2'])
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
    G, vmax, ymax, supg = AS.sliver_G()
    out = {}; refs = {'tight': (L.BASE + '/rigor/assembly_arb.json', L.BASE + '/rigor_assembly_arb.log'), 'v31': (L.BASE31 + '/rigor/assembly_arb.json', L.BASE31 + '/rigor_assembly_arb.log')}
    V31 = ('/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.json', '/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.out')
    AS.WE = dict(AS.WE); AS.WE['windows'] = WR['windows']; AS.GHAT = GH
    for name, ho in SETTINGS.items():
        HO = [D(x) for x in ho]; AS.HO = HO
        R = assemble(HO, G)
        os.chdir(L.BASE); W_ = AS.window_loop(); F4 = AS.f4b(R); os.chdir(_cwd)
        F4h = f4b_honest(R, F4); F4['I3c'] = AS.max_(F4['I3c'], F4h['I3c']); F4['I3s'] = AS.max_(F4['I3s'], F4h['I3s'])
        M = margins(R, W_, F4); C = claims(R, W_, F4); X = export(R)
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
        out[name] = dict(E=R['E'], gp=R['gp'], S0=R['S0'], S1=R['S1'], parts=R['parts'], info=R['info'], margins=M, claims=C, claims_failing=bad, export=X,
                         F4={k: v for k, v in F4.items() if k != 'B'}, B50=F4['B'][50], windows=[dict(yE=w['yE'], gp=w['gp'], der=w['der'], K=w['K'], Ubox=w['Ubox']) for w in W_],
                         ref=dict(E=ref['E'], gp=ref['gp'], margins=ref_margins(refs[name][1])), v31=dict(E=v31['E'], gp=v31['gp'], margins=ref_margins(V31[1])))
    inputs = [L.OUT + '/weights_rep.json', L.OUT + '/ramp_rep.json', L.OUT + '/outer_rep.json', L.BASE + '/rigor/assembly_arb.py', L.BASE + '/lemmaS_F10p_M0.52.json',
              L.BASE + '/rigor/cap_F10p.json', L.BASE + '/lemmaS_cflux.json', '/home/claude/s1spec/step3nl/verifier/rate_table_fin11_a0.00015.json', os.path.abspath(__file__),
              refs['tight'][0], refs['v31'][0], V31[0]]
    L.dump(dict(res=out, top_mode=TOP_MODE, outer_hyp=OR['hyp']), L.OUT + '/assembly_rep.json', inputs)
    print('\nwrote out/assembly_rep.json')
