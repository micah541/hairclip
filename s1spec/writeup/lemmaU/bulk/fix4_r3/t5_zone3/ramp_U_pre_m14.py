#!/usr/bin/env python3
"""Task 5 (r2), Z4 part 1: certified (arb) coefficients of the Lemma-U OUTER ramp chi on r in [r_chi, 1.5 r_chi] = [253, 379.5] lambda_1.

Why not a rerun of lemmaS_v32/ramp_out.py: that script bounds the (B) strip energies, weight omega = 28/y (weight-derivative factor gam = 1),
with Umax = 0.05 (run_F1.box default; U_1 reaches 0.0861 on this ramp) and jet-free A-bounds without the cutoff x jet (cubic) classes.
Lemma U's weight is Psi = chi^2 W(sigma) D phi, |W'/W| up to 4.  Here: the P2 machinery (phase1/ramp/ulib.py, exact |coef| tables, honest
Theta^(m)/Theta enclosed per cell) mirrored to the decreasing cutoff chi = 1 - S(t), t = (s - s_chi(tau))/ln1.5.

Normalisation (per unit eps*Theta at the ramp entry, eps = e^{-4s}): on the ramp at time tau, y = (r_chi lambda_1/lambda(tau)) e^{t ln1.5},
eps = eps_edge(tau) e^{-4 t ln1.5}, eps_edge = (lambda/(r_chi lambda_1))^4.  Centres rho <= 253 lambda_1, so sigma = ln(r/rho) = t ln1.5 + w, w >= 0.
  Theta = W(sigma) D(s) phi(s),  D = 1 + g e^{-4s} (theta = g e^{-4s}/D <= g/253^4),  phi = exp(K_b e^{-s}), |phi^(l)/phi| <= k_phi = x(1+3x+x^2), x = K_b/253.
  RT_i = Theta^(i)/Theta (W, D exact forms as P2; phi as a ball);  c2_l = (chi^2)^(l) (in s);
  rem_m = sum_{l=1..m} C(m,l) c2_l RT_{m-l}                  (cutoff part of (chi^2 Theta)^(m) / Theta)
  rho^  = (1/2)(c2_2 + 2 L c2_1),  L = (ln eps Theta)' = -4 + RT_1   (interpolation remainder per unit Y_{j+1})
  Cxi   = L chi chi' + chi'^2 + chi chi''                     ((eps Theta chi chi')' / (eps Theta): exact S33 foot, lemmaS.tex l.333)
  K_jk  = sup e^{-4 t ln1.5} [ sum_{pairs} w_k sum_{m>=1} |coef| Ab |rem_m| + [k = j] Y_{j+1} |rho^| ]       (pairs: k <= j+1)
  Kz_j[key] = sup e^{-4 t ln1.5} sum |coef| v_X |rem_m|      (cutoff x jet classes, key = sorted pair + jet multiset; folded in zone3.py)
  CXI   = sup e^{-4 t ln1.5} |Cxi|
Global Y_{j+1}: the widened design (0.15, 0.52), U .09, M0z .16, dfrac 2/25, inner-ramp foot adds (port of ghat_design_arb.run, returning YC);
  cross-check: Ghat_22 must reproduce 338258.77 (ghat_design_arb.json).
Outer L4 foot terms (port of the foot block with eta'_max -> chi'_max = 1.875/ln1.5, Young delta_o for the z_3 / global parts, dd = a_min/10 for z_4,
  which is admissible because the inner-ramp foot uses its share of the dissipation on [28, 56] only (disjoint supports)): footw_out_j (data
  weight on the space-time ramp norm of z_2) and Y_add_out (to be added to the global Y by task 4; tiny), W33_j.
usage: PYTHONHASHSEED=0 nice -n 19 python3 -B ramp_U.py [--maxeval N] [--tol x] [--fresh]   -> out/ramp_U.json (resumable: out/ramp_U.partial.json)"""
import sys, os, json, math, time, hashlib, argparse, platform
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
P2 = '/home/claude/s1spec/writeup/lemmaU/phase1/ramp'; GH = '/home/claude/s1spec/writeup/lemmaU/phase1/ghat'
GDA = '/home/claude/s1spec/writeup/lemmaU/phase2a/U_ramp/ghat_design_arb.json'
sys.path.insert(0, P2); sys.path.insert(0, GH)
import sympy as sp
from flint import arb, ctx
import ulib as U
import ghat_arb as GA, weight_consts_arb as WC
from ulib import A, Q_, hull
ctx.prec = 128
OUT = os.path.join(HERE, 'out', 'ramp_U.json'); PART = os.path.join(HERE, 'out', 'ramp_U.partial.json')
ap = argparse.ArgumentParser(); ap.add_argument('--maxeval', type=int, default=60000); ap.add_argument('--tol', type=float, default=0.003)
ap.add_argument('--fresh', action='store_true'); args = ap.parse_args()
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def fup(x):
    v = x.upper() if isinstance(x, arb) else arb(x); f = float(v)
    return f if arb(f) >= v else math.nextafter(f, math.inf)
R = {} if args.fresh or not os.path.exists(PART) else json.load(open(PART))
def save():
    json.dump(R, open(PART + '.tmp', 'w'), indent=1, sort_keys=True); os.replace(PART + '.tmp', PART)
def cached(key, fn):
    if key in R: return R[key]
    v = fn(); R[key] = v; save(); print(f'  {key}: ub {v["ub"]:.6g} lb {v["lb"]:.6g} ({v["nev"]} cells)', flush=True); return v
def bb(f, boxes):
    ub, lb, nev, left = U.bb_sup(f, boxes, tol=args.tol, maxeval=args.maxeval)
    return dict(ub=ub, lb=lb, nev=nev, left=left)

RCHI = 253; LN15 = (arb(3)/2).log(); M0q, M1q, UQ, M0zq, DFRAC = Q_(15, 100), Q_(52, 100), Q_(9, 100), Q_(4, 25), Q_(2, 25)
KB = A(Q_(1, 1000)); xph = KB/RCHI; KPHI = xph*(1 + 3*xph + xph*xph)
T = U.load_tables(); assert U.TABLES == GA.TABLES

# ---------------------------------------------------------------- global Y_{j+1} (widened design), Ghat cross-check
r = WC.certify(4, 4, Q_(0), Q_(3, 10), Q_(3, 5), 56); b = r['bound']
Gam = {0: arb(1), 1: b['G1'], 2: b['G2'], 3: b['G3']}; KAP = b['kap']; KAP2 = b['kap2']
etap_in = A(Q_(15, 8))/U.LN2; DEL = A(Q_(1, 10))
AsymG = GA.build_At_sym(UQ); AtG = GA.eval_At(AsymG, M0q, M1q); AbG = {k: d.get((), arb(0)) for k, d in AtG.items()}
YFULL, GHAT, ITEMS, DIS = {}, {}, {}, {}
for j in (0, 1, 2):
    YCb = GA.YC_base(T[j], j, AbG, Gam); add, dis, spare, ndfr, items = GA.l4_extra_YC_x(T, j, M0q, M1q, AtG, Gam, KAP, dfrac=DFRAC, s33='ibp', M0zq=M0zq)
    fadd = {k: arb(0) for k in add}; extra = Q_(0); M1 = A(M1q)
    w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0))
    if any(it[0] == 'S33' for it in items): c1 = w33*M1*etap_in; fadd[3] += c1*arb(KAP)*DEL; extra = Q_(1, 10)
    for it in items:
        if it[0] != 'L4': continue
        if it[1] == 3: fadd[3] += it[2]*M1*etap_in*DEL
        elif it[1] == 2: fadd[2] += it[2]*M1*etap_in
    ok, sl = GA.budget_ok(dis, ndfr, M0q, extra); assert ok, ('budget', j)
    YC = {k: YCb[k] + add.get(k, 0) + fadd.get(k, 0) for k in YCb}
    GHAT[j] = GA.ghat_x(YC, j, M0q, KAP2)[j]; YFULL[j] = YC[j + 1]; ITEMS[j] = items; DIS[j] = (dis, ndfr, extra)
g22ref = json.load(open(GDA))['cases']['widened (0.15,0.52), U .09, M0z .16']['Ghat_jj'][2]
assert abs(fup(GHAT[2])/g22ref[1] - 1) < 1e-9, (fup(GHAT[2]), g22ref)
print('global Y_{j+1}:', [fup(YFULL[j]) for j in range(3)], 'Ghat_jj', [fup(GHAT[j]) for j in range(3)], '(Ghat_22 reproduces ghat_design_arb.json)', flush=True)

# ---------------------------------------------------------------- outer L4 foot (port; chi'_max = 1.875/ln1.5)
CHIP = A(Q_(15, 8))/LN15; DELO = A(Q_(1, 1000)); FOOT = {}
for j in (0, 1, 2):
    M1 = A(M1q); amin = (1 + A(M0q)**2)**-2; fw = arb(0); yadd = {2: arb(0), 3: arb(0)}; items = ITEMS[j]
    w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0)); has33 = any(it[0] == 'S33' for it in items)
    if has33:
        c1 = w33*M1*CHIP; yadd[3] += c1*arb(KAP)*DELO; fw += c1*arb(KAP)/DELO
        dd = A(Q_(1, 10))*amin; fw += c1**2/(4*dd)
    for it in items:
        if it[0] != 'L4': continue
        w = it[2]
        if it[1] in (2, 3): yadd[it[1]] += w*M1*CHIP*DELO; fw += w*M1*CHIP/DELO
        else: raise ValueError(('foot r', it))
    dis, ndfr, extra = DIS[j]; ok, sl = GA.budget_ok(dis, ndfr, M0q, Q_(1, 10) if has33 else Q_(0))
    assert ok, ('outer foot budget', j)
    FOOT[j] = dict(footw_out=fup(fw), W33=fup(w33), Y_add_out={str(k): fup(v) for k, v in yadd.items()},
                   Y_add_rel={str(k): fup(v/YFULL[k - 1]) if k - 1 <= 2 and k - 1 >= 0 else None for k, v in yadd.items()}, budget_ok=bool(ok))
print('outer foot:', FOOT, flush=True)

# ---------------------------------------------------------------- ramp coefficient arrays (P2 port; box (0.15, 0.52), U = 9/100)
AsymR = U.build_At_sym(UQ); AtR = U.eval_At(AsymR, M0q, M1q); AbR = U.ab_from_At(AtR)
def coef_arrays(Ab):
    C = {}
    for j in (0, 1, 2):
        n = j + 2; C[j] = {k: {} for k in range(j + 2)}
        for (a, b_), mons in T[j].items():
            if (a, b_) == (n, n): continue
            for k in range(j + 2):
                w = 0
                if a == k: w += Q_(1) if a == b_ else Q_(1, 2)
                elif b_ == k: w += Q_(1, 2)
                if w == 0: continue
                for (m, i, rr), co in mons.items():
                    if m == 0: continue
                    C[j][k][m] = C[j][k].get(m, arb(0)) + A(w*co)*Ab[(i, rr)]
    return C
def coef_Z(At):
    Kz = {}
    for j in (0, 1, 2):
        n = j + 2; Kz[j] = {}
        for (a, b_), mons in T[j].items():
            if (a, b_) == (n, n): continue
            for (m, i, rr), co in mons.items():
                if m == 0: continue
                for X, v in At[(i, rr)].items():
                    if not X: continue
                    key = tuple(sorted([a, b_] + list(X))); Kz[j].setdefault(key, {}); Kz[j][key][m] = Kz[j][key].get(m, arb(0)) + A(co)*v
    return Kz
C = coef_arrays(AbR); KZ = coef_Z(AtR)
for j in (0, 1, 2):
    for (a, b_), mons in T[j].items():
        if (a, b_) != (j + 2, j + 2) and any(m >= 1 for (m, i, rr) in mons): assert max(a, b_) <= j + 1, (j, a, b_)

# ---------------------------------------------------------------- geometry on (t, w, a, g)
_t = sp.Symbol('t'); _S = _t**3*(10 - 15*_t + 6*_t**2); _chi = 1 - _S
C2P = [sp.Poly(sp.diff(_chi**2, _t, l), _t) for l in range(4)]; CP = [sp.Poly(sp.diff(_chi, _t, l), _t) for l in range(3)]
TC2 = [U.TP1(P) for P in C2P]; TC = [U.TP1(P) for P in CP]
def c2(t, l): return TC2[l](t)/U.ipow(LN15, l)
def ch(t, l): return TC[l](t)/U.ipow(LN15, l)
LNR = arb(RCHI).log(); WT = 16.0
def geom(box, point=False):
    tb, wb, ab, gb = box
    t = arb((tb[0] + tb[1])/2) if point else hull(tb[0], tb[1])
    a = arb((ab[0] + ab[1])/2) if point else hull(ab[0], ab[1])
    g = arb((gb[0] + gb[1])/2) if point else hull(gb[0], gb[1])
    if wb[1] == math.inf:
        sigT = (arb(tb[0]) if not point else t)*LN15 + wb[0]; p, q = U.pq_upper_tail(sigT, a)
    else:
        w = arb((wb[0] + wb[1])/2) if point else hull(wb[0], wb[1]); p, q = U.pq_from_sigma(t*LN15 + w, a)
    Rk = U.Rk(p, q)
    s = LNR + (arb(tb[0]) if not point else t)*LN15                         # smallest s on the cell at the earliest time (theta decreasing in s)
    th = hull(0, U.theta_s(s, g).upper())
    RT = {0: arb(1), 1: Rk[1] - 4*th, 2: Rk[2] + th*(16 - 8*Rk[1]), 3: Rk[3] + th*(-12*Rk[2] + 48*Rk[1] - 64)}
    RTp = {0: arb(1)}
    for i in (1, 2, 3):
        rad = sum((math.comb(i, l)*abs(RT[i - l]) for l in range(1, i + 1)), arb(0))*KPHI
        RTp[i] = RT[i] + hull(-rad.upper(), rad.upper())
    E = [c2(t, l) for l in range(4)]
    rem = {m: sum((math.comb(m, l)*E[l]*RTp[m - l] for l in range(1, m + 1)), arb(0)) for m in (1, 2, 3)}
    L = -4 + RTp[1]
    rho = (E[2] + 2*L*E[1])/2
    X = [ch(t, l) for l in range(3)]
    cxi = L*X[0]*X[1] + X[1]*X[1] + X[0]*X[2]
    dec = (-4*LN15*t).exp()
    return dict(rem=rem, rho=rho, cxi=cxi, dec=dec, L=L, RT1=RTp[1])
def mk(coef=None, Y=None, which='K'):
    def f(box, point=False):
        q = geom(box, point)
        if which == 'K':
            tot = arb(0)
            for m, c in (coef or {}).items(): tot += c*abs(q['rem'][m])
            if Y is not None: tot += Y*abs(q['rho'])
            return q['dec']*tot
        if which == 'cxi': return q['dec']*abs(q['cxi'])
        if which == 'L': return q['L']
        if which == 'RT1': return q['RT1']
        raise ValueError(which)
    f.scales = [0.4055, 1.0, 1.0, 1e-5]; return f
def boxes():
    ts = [(i/16, (i + 1)/16) for i in range(16)]
    ws = [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, WT), (WT, math.inf)]
    return [(tt, ww, (0.0, 0.3), (float(U.G_LO), float(U.G_HI))) for tt in ts for ww in ws]

print('ramp sups', flush=True)
GEO = {w: cached(f'geom|{w}', lambda w=w: bb(mk(which=w), boxes())) for w in ('L', 'RT1', 'cxi')}
K = {}
for j in (0, 1, 2):
    K[j] = {}
    for k in range(j + 2):
        Y = YFULL[j] if k == j else None
        K[j][k] = cached(f'K|{j}{k}', lambda j=j, k=k, Y=Y: bb(mk(C[j][k], Y), boxes()))['ub']
KPAIR = {}
for j in (0, 1, 2):
    KPAIR[j] = {k: cached(f'Kpair|{j}{k}', lambda j=j, k=k: bb(mk(C[j][k], None), boxes()))['ub'] for k in range(j + 2)}
KZV = {}
for j in (0, 1, 2):
    KZV[j] = {str(key): cached(f'Kz|{j}|{key}', lambda cm=cm: bb(mk(cm, None), boxes()))['ub'] for key, cm in sorted(KZ[j].items())}

res = dict(label='cert (arb, branch and bound over (t, w, a, g); upper endpoints rounded up)',
           ramp=dict(r=[RCHI, 379.5], cutoff='chi = 1 - S(t), t = ln(r/253)/ln1.5, quintic S', chi_prime_max=fup(CHIP)),
           weight=dict(W='(4,4,a), a in [0,0.3]', g=[str(U.G_LO), str(U.G_HI)], K_b='1/1000', k_phi=fup(KPHI), centres='rho <= 253 lambda_1 (sigma >= 0 on the ramp)'),
           box=dict(M0=str(M0q), M1=str(M1q), U=str(UQ), M0z=str(M0zq), dfrac=str(DFRAC)),
           normalisation='per unit eps_edge * Theta, eps_edge = (lambda/(253 lambda_1))^4; includes the factor e^{-4 t ln1.5}',
           Gamma={str(m): fup(Gam[m]) for m in (1, 2, 3)}, kappa=fup(arb(KAP)), kappa2=fup(arb(KAP2)),
           Y_next_global=[fup(YFULL[j]) for j in range(3)], Ghat_jj_global=[fup(GHAT[j]) for j in range(3)],
           K={str(j): {str(k): v for k, v in d.items()} for j, d in K.items()},
           K_pair_only={str(j): {str(k): v for k, v in d.items()} for j, d in KPAIR.items()},
           Kz={str(j): d for j, d in KZV.items()}, CXI=GEO['cxi']['ub'], L_sup=GEO['L']['ub'], RT1_sup=GEO['RT1']['ub'],
           foot_out={str(j): d for j, d in FOOT.items()}, delta_o=str(Q_(1, 1000)),
           bb={k: dict(lb=v['lb'], nev=v['nev'], left=v['left']) for k, v in sorted(R.items())},
           args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',) if os.environ.get(k)}, python=platform.python_version(),
                     flint=__import__('flint').__version__, prec=ctx.prec, maxeval=args.maxeval, tol=args.tol,
                     inputs={p: sha(p) for p in (U.TABLES, os.path.join(P2, 'ulib.py'), os.path.join(GH, 'ghat_arb.py'), os.path.join(GH, 'weight_consts_arb.py'), GDA, os.path.abspath(__file__))}))
json.dump(res, open(OUT + '.tmp', 'w'), indent=1, sort_keys=True); os.replace(OUT + '.tmp', OUT)
print('wrote', OUT, '\nK', res['K'], '\nKz', res['Kz'], '\nCXI', res['CXI'], flush=True)
