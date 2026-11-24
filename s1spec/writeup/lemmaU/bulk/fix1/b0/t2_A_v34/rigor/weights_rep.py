#!/usr/bin/env python3
"""Phase-2a repair of Lemma S(A), stage 1: quadratic coefficient bounds (hat G), certified (exact Rationals + arb, outward).
Repairs applied here (v3.1 files untouched; base functions imported from the copy base_ibp/rigor/weights_exact.py):
  R3  profile jets U_j (j>=2) in the jet-free coefficient bounds: Ab := Z-free part of the symbolic Aterms (weights_exact.ab_exact sets V_j = 0, j>=2,
      weights_exact.py:81-83); used for the core, the ramp [14,28] and the (B) windows.  Windows also use the honest profile bound U(1.5 y_E) instead of 1/20.
  R4  re-split of the core dissipation at levels 0,1: interpolation share Bc_j (v3.1: 3/2 a_min) lowered so that a fraction r_j*a_min of the core
      dissipation 2 a(V) eps*omega z_{j+2}^2 is RETAINED (not used); r_j is a certified lower bound (arb) and is exported by the assembly.
  core {2,3,3} by IBP (s33='ibp'), as in gonogo/lemmaS_A_ibp.
Segment checks (new): the core rate eps(28)*Ghat_core must dominate eps(y)*Ghat(box(y)) wherever the core gamma' accounting is used:
      [28, 220]: U <= 1/20 (U(220) < 0.0493) = core;  [220, 350]: same box (0.11, 0.52), U <= 0.08 (U(350) < 0.0792);
      [350, 700] (outer cutoff region of R5): HYPOTHESIS box H-O1 (M0o, M1o), U <= U(700) < 0.162, sup|v| <= M0zo.
Output: out/weights_rep.json.  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/weights_rep.py [--retain 1.25] [--outer-box 1/2 1 1/2]"""
import os, sys, json, argparse, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX, A, Q_, arb

def YC_base(tab, j, Ab):
    return WX.YC_base(tab, j, Ab)

def ghat_bc(YC, j, M0q, bc, kap=5):
    amin = (1 + A(M0q)**2)**-2; Bc = A(bc)*amin; G = {k: YC[k] for k in range(j + 1)}
    G[j] = G[j] + YC[j + 1]*(YC[j + 1]/(4*Bc) + A(Q_(kap**2, 2))); return G

def used_frac(dis, ndfr, M0q, extra=Q_(0)):
    """dissipation used by the L4/Young items, in units of a_min (arb upper): dis/amin + ndfr + extra"""
    amin = (1 + A(M0q)**2)**-2
    return dis/amin + A(Q_(ndfr) + extra)

def core(T, M0q, M1q, Uq, bcs, M0zq=WX.M0Z, honest=True, lam=None):
    """lam (R5): the L4 identity and the S33 IBP are taken over [14,350] (outer cutoff starts at 350); their boundary terms at y=350 enter as
    w*sqrt(O_b)*N_r, bounded by lam*N_r^2 + w^2 O_b/(4 lam): YC[r] += lam for every L4 class r and YC[3] += lam for S33 (kappa N_3 part);
    the S33 N_4 part uses the level-2 dissipation reserved by Bc_2 = 3/2 - eps4 (rates in outer_rep.py)."""
    At = L.at_eval(M0q, M1q, Uq); Ab = L.ab_honest(At) if honest else L.ab_v31(M0q, M1q, Uq, 4); out = {}
    for j in (0, 1, 2):
        YCb = YC_base(T[j], j, Ab)
        add, dis, spare, ndfr, items = _l4(T, j, M0q, M1q, At, M0zq)     # = weights_exact.l4_extra_YC_x(s33='ibp') with sup|z_0| <= M0zq
        YC = {k: YCb[k] + add.get(k, 0) for k in YCb}
        l4w = {}
        for it in items:
            if it[0] == 'L4': l4w[it[1]] = l4w.get(it[1], arb(0)) + it[2]
            elif it[0] == 'S33': l4w['S33'] = l4w.get('S33', arb(0)) + it[1]
        if lam is not None:
            for r in l4w: YC[3 if r == 'S33' else r] += A(lam)
        used = used_frac(dis, ndfr, M0q); assert (used - A(Q_(1, 2))).upper() <= 0, ('core L4 budget', j)
        ret = A(2) - A(bcs[j]) - used                 # retained dissipation / a_min (lower end certified)
        assert ret.lower() >= 0, ('negative retained dissipation', j)
        out[j] = dict(YC=YC, Ghat=ghat_bc(YC, j, M0q, bcs[j]), used=used, retained=ret, Bc=str(bcs[j]), YCb=YCb, l4w=l4w)
    return out

def _l4(T, j, M0q, M1q, At, M0zq):
    """classify with a custom sup|z_0| (M0zq) -- weights_exact.classify_x takes M0zq as keyword"""
    n = j + 2; amin = (1 + A(M0q)**2)**-2; dfrac = Q_(1, 10)
    add = {k: arb(0) for k in range(n + 1)}; dissip = [arb(0)]; ndp = [0]; M1 = A(M1q)
    def put(x, y, wt):
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y))
        if y == n and x == n: dissip[0] += wt
        elif y == n: d = A(dfrac)*amin; add[x] += wt**2/(4*d); ndp[0] += 1
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
    items = WX.classify_x(T, j, M1q, At, M0zq=M0zq); kap = 5
    for it in items:
        if it[0] == 'P': put(it[1], it[2], it[3])
        elif it[0] == 'L4': put(3, it[1], 3*M1*it[2]); put(2, it[1], kap*M1*it[2])
        elif n >= 4:
            w = it[1]; put(3, 3, A(Q_(3, 2)*kap)*M1*w); put(3, 4, A(Q_(3, 2))*M1*w); put(2, 3, A(Q_(kap**2, 2))*M1*w); put(2, 4, A(Q_(kap, 2))*M1*w)
        else: raise ValueError('S33 at n<4')
    return add, dissip[0], amin/2, ndp[0]*dfrac, items

def ramp(T, M0q, M1q, Uq, honest=True, dfrac=Q_(2, 25)):
    """port of weights_exact.ramp (rs33 'ibp2'), with honest Ab"""
    At = L.at_eval(M0q, M1q, Uq); Ab = L.ab_honest(At) if honest else L.ab_v31(M0q, M1q, Uq, 4); Lg = arb(2).log()
    etap = A(Q_(15, 8))/Lg; DEL = A(Q_(1, 10)); amin_r = (1 + A(M0q)**2)**-2; M1 = A(M1q); out = {}
    for j in (0, 1, 2):
        YCb = YC_base(T[j], j, Ab)
        add, dis, spare, ndfr, items = WX.l4_extra_YC_x(T, j, M0q, M1q, At, dfrac=dfrac, s33='ibp2')
        fw = arb(0); extra = Q_(0); w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0)); has33 = any(it[0] == 'S33' for it in items)
        if has33:
            c1 = w33*M1*etap; kap = 5
            add[3] = add.get(3, arb(0)) + c1*kap*DEL; fw += c1*kap/DEL
            dd = A(Q_(1, 10))*amin_r; extra = Q_(1, 10); fw += c1**2/(4*dd)
        for it in items:
            if it[0] != 'L4': continue
            w = it[2]
            if it[1] == 3: add[3] = add.get(3, arb(0)) + w*M1*etap*DEL; fw += w*M1*etap/DEL
            elif it[1] == 2: add[2] = add.get(2, arb(0)) + w*M1*etap; fw += w*M1*etap
            else: raise ValueError(('foot r', it))
        used = used_frac(dis, ndfr, M0q, extra); slack = A(Q_(1, 2)) - used; assert slack.lower() >= 0, ('ramp budget', j)
        YC = {k: YCb[k] + add.get(k, arb(0)) for k in YCb}
        out[j] = dict(YC=YC, Ghat=ghat_bc(YC, j, M0q, Q_(3, 2)), footw=fw, W33=w33, slack=slack, has33=has33, YCb=YCb)
    return out

def windows(T, al, dv, honest=True):
    res = {}
    for (yE, M0q, M1q) in WX.window_boxes(al, dv):
        yb = A(Q_(3, 2)*yE); Uq = WX.up_rational(L.U_of_y(yb, al)) if honest else Q_(1, 20)
        At = L.at_eval(M0q, M1q, Uq); Ab = L.ab_honest(At) if honest else L.ab_v31(M0q, M1q, Uq, 4); Gk = {}
        for j in range(3):
            YCb = YC_base(T[j], j, Ab); add, dis, spare, nd, _ = WX.l4_extra_YC_x(T, j, M0q, M1q, At, s33='ibp'); assert WX.budget_ok(dis, nd, M0q)[0]
            YC = {k: YCb[k] + add.get(k, 0) for k in YCb}; Gk[j] = WX.ghat_x(YC, j, M0q)
        YCb = YC_base(T[3], 3, Ab)
        add, dis, spare, nd, _ = WX.l4_extra_YC_x({3: T[3]}, 3, M0q, M1q, At, fallback={2: Q_(1, 2), 3: Q_(0)}, s33='ibp'); assert WX.budget_ok(dis, nd, M0q)[0]
        YC3 = {k: YCb[k] + add.get(k, 0) for k in YCb}; Gk[3] = WX.ghat_x(YC3, 3, M0q)
        res[str(yE)] = dict(M0q=str(M0q), M1q=str(M1q), Uq=str(Uq), Ghat={j: Gk[j] for j in range(4)}, YC3=YC3)
    return res

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--retain', default='5/4', help='target retained core dissipation (units a_min) at levels 0,1')
    ap.add_argument('--outer-box', nargs=3, default=['1/2', '1', '1'], help='HYPOTHESIS H-O1 on [350,700]: |f|<=M0o, y|f\'|<=M1o, |v|<=M0zo (info only)')
    ap.add_argument('--lam', default='1/10', help='R5: YC add for the boundary terms of the L4/S33 identities at y=350')
    ap.add_argument('--eps4', default='1/1000', help='R5: level-2 dissipation reserved for the S33 N_4 boundary part (units a_min)')
    a = ap.parse_args()
    T = WX.build_tables(); res = {}; M0c, M1c, Uc = Q_(11, 100), Q_(52, 100), Q_(1, 20)
    # ---- v3.1-convention core (for the comparison and to read off the L4 usage), Bc = 3/2 everywhere
    c31 = core(T, M0c, M1c, Uc, [Q_(3, 2)]*3, honest=False); ch = core(T, M0c, M1c, Uc, [Q_(3, 2)]*3, honest=True)
    print('core Ghat_jj  v3.1 Ab / honest Ab (R3):', [('%.2f' % L.fup(c31[j]['Ghat'][j]), '%.2f' % L.fup(ch[j]['Ghat'][j])) for j in range(3)], flush=True)
    print('core L4 dissipation used (a_min units, upper):', ['%.4f' % L.fup(ch[j]['used']) for j in range(3)])
    # ---- R4: choose Bc_j (j = 0,1) as a rational (1/100 grid, rounded down) so that the retained fraction >= target
    tgt = Q_(a.retain); bcs = []; LAM = Q_(a.lam); EPS4 = Q_(a.eps4)
    for j in range(3):
        if j == 2: bcs.append(Q_(3, 2) - EPS4); continue
        u = ch[j]['used'].upper(); import math; b = Q_(math.floor(float(((A(2) - A(tgt) - u)*100).lower())), 100); assert b > 0; bcs.append(min(b, Q_(3, 2)))
    cr = core(T, M0c, M1c, Uc, bcs, honest=True, lam=LAM)
    print('R4 re-split: Bc =', [str(b) for b in bcs], 'retained (a_min units, certified lower) =', ['%.4f' % L.flo(cr[j]['retained']) for j in range(3)])
    print('   Ghat_jj after re-split:', ['%.2f' % L.fup(cr[j]['Ghat'][j]) for j in range(3)], flush=True)
    res['core_v31Ab'] = {j: dict(Ghat=c31[j]['Ghat'], YC=c31[j]['YC']) for j in range(3)}
    res['core_honest_bc32'] = {j: dict(Ghat=ch[j]['Ghat'], YC=ch[j]['YC'], used=ch[j]['used'], retained=ch[j]['retained']) for j in range(3)}
    res['core'] = {j: dict(Ghat=cr[j]['Ghat'], YC=cr[j]['YC'], used=cr[j]['used'], retained=cr[j]['retained'], Bc=cr[j]['Bc'], l4w={str(k): v for k, v in cr[j]['l4w'].items()}) for j in range(3)}
    res['lam'] = str(LAM); res['eps4'] = str(EPS4)
    # ---- ramp [14,28], box (0.05, 0.9), U <= 1/20, honest Ab
    r31 = ramp(T, Q_(5, 100), Q_(9, 10), Uc, honest=False); rh = ramp(T, Q_(5, 100), Q_(9, 10), Uc, honest=True)
    print('ramp Ghat_jj v3.1 / honest:', [('%.2f' % L.fup(r31[j]['Ghat'][j]), '%.2f' % L.fup(rh[j]['Ghat'][j])) for j in range(3)],
          'slack (a_min):', ['%.4f' % L.flo(rh[j]['slack']) for j in range(3)], 'footw', ['%.6g' % L.fup(rh[j]['footw']) for j in range(3)], 'W33', ['%.6g' % L.fup(rh[j]['W33']) for j in range(3)], flush=True)
    res['ramp'] = {j: dict(Ghat=rh[j]['Ghat'], YC=rh[j]['YC'], footw=rh[j]['footw'], W33=rh[j]['W33'], slack=rh[j]['slack']) for j in range(3)}
    res['ramp_v31Ab'] = {j: dict(Ghat=r31[j]['Ghat'], YC=r31[j]['YC'], footw=r31[j]['footw'], W33=r31[j]['W33'], slack=r31[j]['slack']) for j in range(3)}
    # ---- windows (B), honest Ab and honest U(1.5 y_E)
    res['windows'] = windows(T, Q_(15, 100000), Q_(15, 100), honest=True)
    for k, w in res['windows'].items(): print('window', k, 'U<=%s' % w['Uq'], 'Ghat_jj', ['%.2f' % L.fup(w['Ghat'][j][j]) for j in range(4)], flush=True)
    # ---- segment checks for the core gamma' accounting
    seg = {}
    # exported retained fraction: rational, 1/1000 grid, below the certified lower ends (core); the segments must retain at least as much
    rexp = [Q_(math.floor(float(cr[j]['retained'].lower())*1000 - 1e-6), 1000) for j in range(2)] + [Q_(0)]      # level 2: eps4 reserved, not exported
    print('exported retained fractions r_j (rational, certified <= retained):', [str(r) for r in rexp])
    U220 = L.U_of_y(220); U350 = L.U_of_y(350); U700 = L.U_of_y(700)
    assert (U220 - A(Q_(1, 20))).upper() < 0 and (U350 - A(Q_(8, 100))).upper() < 0 and (U700 - A(Q_(162, 1000))).upper() < 0
    gs = core(T, M0c, M1c, Q_(8, 100), bcs, honest=True, lam=LAM)
    M0o, M1o, M0zo = (Q_(x) for x in a.outer_box)
    go = core(T, M0o, M1o, Q_(162, 1000), bcs, M0zq=M0zo, honest=True, lam=LAM)     # info only: [350,700] is bounded term by term in outer_rep.py
    e28 = A(28)**-4
    for name, g, y0 in (('[220,350] U<=0.08', gs, 220), ('[350,700] H-O1 box (info only)', go, 350)):
        ok = True; rows = {}
        for j in range(3):
            for k in range(j + 1):
                lhs = A(y0)**-4*g[j]['Ghat'][k]; rhs = e28*cr[j]['Ghat'][k]; rows[f'{j}{k}'] = (L.fup(lhs), L.flo(rhs)); ok = ok and (lhs - rhs).upper() <= 0
            ok = ok and (g[j]['retained'] - A(rexp[j])).lower() >= 0
        seg[name] = dict(ok=ok, rows=rows, retained=[L.flo(g[j]['retained']) for j in range(3)])
        print('segment', name, 'eps(y0)Ghat_seg <= eps(28)Ghat_core for all jk:', ok, 'retained', seg[name]['retained'], flush=True)
    res['segments'] = seg; res['r_export'] = [str(r) for r in rexp]; res['bcs'] = [str(b) for b in bcs]; res['outer_box_hyp'] = dict(M0o=str(M0o), M1o=str(M1o), M0zo=str(M0zo), U700=L.fup(U700))
    inputs = [L.BASE + '/rigor/weights_exact.py', L.BASE + '/rigor/tables_exact.pkl', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py']
    L.dump(res, L.OUT + '/weights_rep.json', inputs)
    print('wrote', L.OUT + '/weights_rep.json')
