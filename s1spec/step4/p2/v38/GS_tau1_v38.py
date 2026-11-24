#!/usr/bin/env python3
"""GS_tau1.py -- [cert] (GS) at tau_1: the tip arc of the true curve at every good time of base slab 180 is a normal graph over S_1 on
[v38 copy of <s1spec>/writeup/lemmaU/bulk/fix15_review/c/GS_tau1.py for the step-4 v38 re-pin (p2/v38/final_tail_v38.py); changes: output path out_v38/iface/GS_tau1.json; cpu_s recorded as null (time printed to the log)]
|x| <= 60.25, with sup|k u| and sup|u'| bounded (fix round 14; GS_lemma.tex, Lemma GS-1 and the tube Lemma T).

Method (read-only; step-4 and step-1 code is imported with -B, nothing is written outside this directory):
 * base slab 180 at 16 t_c-balls of radius 1/16 covering t_c in [-1, 1] (p2/t9h/basegeo.Geo, the T6 engine on t1/polished_v2, the objects
   of step 4's H1 stage); each element is cut into NS s-balls (cells).  Per cell: hulls of the base position Gam, tangent angle th~,
   curvature k~ and base arclength x.
 * true curve c = (Gam + e nu~ - x_c)/lambda in the step-3 frame, (x_c, lambda) = P_hull of p2/out/h1slab/summary.json (the 32 IFT balls of
   slab 180) widened by 1e-9 (covers the true-curve shift <= 3.3e-11 of step4_proof.tex S9); |e| <= 1e-5 for y < 28 and on x <= 34.25,
   |e| <= 1.7 E(x) for y >= 28 with sqrt(F0) <= aF Ybar + bF of slab 180 (amp_A1_w180.json), Ybar = 1.7e-12 (box level, conservative);
   tangent of c: th~ + beta, |beta| <= |e'|/(1 - |k~ e|), |e'| <= 1e-5 (y < 28) or 1.5e-3 (box b0).
 * S_1 (p2/t9h/s1chart2: prof_tight + ext813 rows with the tip-position ball, tail envelopes beyond sigma = 13): hulls of position, angle and
   curvature on a fixed x-grid (width 1/64 up to x13, far chart beyond with width 1/8), mirrored for x < 0.
 * per cell: delta >= sup|c(xi) - S1(x(xi))|, eta >= sup|angle(c'(xi)) - theta_S1(x(xi))|, K >= sup|k_S1| on [x_lo - 2 delta, x_hi + 2 delta];
   asserted:  q1 := 3 K delta < 1/2,  eta + 2 K delta < pi/2,  q2 := (1 + 3 K delta) tan(eta + 2 K delta) < 1.
   Then Lemma T gives a unique foot X(xi) in [x - 2 delta, x + 2 delta], X' > 0, |u| <= 3 delta, |k u| <= q1, |u'| <= q2.
 * coverage: the last cell's x_lo - 2 delta >= 60.25 (asserted), and the tip cell contains xi = 0 (X(0) = 0 by symmetry).
Output GS_tau1.json.  Run: nice -n 19 python3 -B GS_tau1.py [NS]"""
import sys, os, json, hashlib, platform, math, time
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3')
from flint import arb, acb, ctx
ctx.prec = 128
import basegeo, s1chart2
import hstail as HT
HERE = os.path.dirname(os.path.abspath(__file__))
S4 = '/home/claude/s1spec/step4'
SUMMARY = f'{S4}/p2/out/h1slab/summary.json'
AMP = f'{S4}/p2/t8se/constants/amp_A1_w180.json'
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 16
NTC = 16
XS = 60.25
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
I = acb(0, 1)
t0 = time.time()


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def iv(lo, hi):
    return arb(lo).union(arb(hi))


summ = json.load(open(SUMMARY)); PH = summ['P_hull']
XC = iv(PH['xc'][0] - 1e-9, PH['xc'][1] + 1e-9); LAM = iv(PH['lam'][0] - 1e-9, PH['lam'][1] + 1e-9)
amp = json.load(open(AMP)); i180 = amp['slabs'].index(180)
YB = arb('1.7e-12'); SF0 = arb(repr(amp['per_slab']['0']['aF'][i180])) * YB + arb(repr(amp['per_slab']['0']['bF'][i180]))
CE = arb('1.6031e4'); X0 = arb('2.2490776515'); YCORE = 28.0 + 2.2490776515


def ebox(xl, xu):
    """(|e| bound, |e'| bound) on base arclength [xl, xu]"""
    e0, e1 = arb('1e-5'), arb('1e-5')
    if xu >= YCORE and xu > 34.25:
        yu = arb(xu) - X0
        Ef = arb('1.7') * (CE * YB + SF0 * (yu * yu * yu / 84).sqrt())
        e0 = arb(max(U(Ef), U(e0))); e1 = arb('1.5e-3')
    elif xu >= YCORE:                       # y >= 28 but still x <= 34.25: both boxes hold; take the larger (conservative)
        yu = arb(xu) - X0
        Ef = arb('1.7') * (CE * YB + SF0 * (yu * yu * yu / 84).sqrt())
        e0 = arb(max(U(Ef), U(e0))); e1 = arb('1.5e-3')
    return e0, e1


# ---------------- S1 hull cache
S = s1chart2.S1(); X13 = float((s1chart2.T0 + 5).lower()) - 1e-6
H_IN = 1.0 / 64; H_FAR = 1.0 / 8; XMAX = 80.0
cache = []                                   # (xlo, xhi, G (acb hull), th (arb hull), |k| upper)
x = 0.0
while x < X13 - 1e-12:
    xh = min(x + H_IN, X13)
    c = S.chart_ball(iv(x, xh)); cache.append((x, xh, c['G'], c['th'], U(c['k'])))
    x = xh
x = X13
while x < XMAX:
    xh = x + H_FAR
    c = S.far_chart(iv(x, xh), step=0.25)
    cache.append((x, xh, c['G'], c['th'], U(c['k'])))
    x = xh
c0 = cache[0]                                 # mirror of [0, H_IN] for x in [-H_IN, 0]: G -> conj G, theta -> -pi - theta
cache.insert(0, (-H_IN, 0.0, acb(c0[2].real, -c0[2].imag), -arb.pi() - c0[3], c0[4]))
KB_END = U(HT.profile_bounds(arb(XMAX) - s1chart2.T0 + 8 - arb(1e-6))[0])   # |k| beyond the cache (envelopes decreasing in sigma)
print(f'S1 cache: {len(cache)} intervals ({time.time() - t0:.0f}s); x13 = {X13:.6f}', flush=True)


def s1_hull(xl, xu):
    G = th = None; K = 0.0
    assert xl >= -H_IN and xu <= XMAX, (xl, xu)
    for (a, b, Gc, tc, kc) in cache:
        if b < xl or a > xu: continue
        G = Gc if G is None else acb(G.real.union(Gc.real), G.imag.union(Gc.imag)); th = tc if th is None else th.union(tc); K = max(K, kc)
    return G, th, K


def s1_kmax(xl, xu):
    lo, hi = xl, xu
    if lo < 0: hi = max(hi, -lo); lo = 0.0
    K = 0.0
    for (a, b, Gc, tc, kc) in cache:
        if b < lo or a > hi: continue
        K = max(K, kc)
    if hi > XMAX: K = max(K, KB_END)
    return K


# ---------------- cells
inputs = {SUMMARY: sha(SUMMARY), AMP: sha(AMP)}
for pth in ('/home/claude/s1_tight/prof/prof_tight.txt', '/home/claude/s1_tight/ext813/ext813_prof.txt'):
    inputs[pth] = sha(pth)
for m in ('basegeo', 's1chart2', 'hstail'):
    mod = sys.modules[m]; inputs[mod.__file__] = sha(mod.__file__)
worst = dict(q1=0.0, q2=0.0, delta_tip=0.0, eta=0.0); balls = []
last_cov = None
for jt in range(NTC):
    TC = arb(-1 + (2 * jt + 1) / NTC, 1.0 / NTC)
    g = basegeo.Geo(180, 64.0, TC)
    nel = len(g.els); cov = -1e9; bw = dict(q1=0.0, q2=0.0); cells = 0; prev_hi = None
    for i in range(nel):
        for j in range(NS):
            Sb = arb(-1 + (2 * j + 1) / NS, 1.0 / NS)
            q = g.at(i, Sb)
            xl, xu = float(q['x'].lower()), float(q['x'].upper())
            if xl > 75: break
            e0, e1 = ebox(xl, xu)
            kt = abs(q['k0'])
            cpt = (q['Gam'] - XC) / LAM
            cpt = acb(cpt.real + arb(0, U(e0 / LAM)), cpt.imag + arb(0, U(e0 / LAM)))
            beta = arb(0, U(e1 / (1 - kt * e0)))
            thc = q['th'] + beta
            G1, th1, _ = s1_hull(xl, xu)
            eta = U(abs(thc - th1))
            if xu <= X13:
                # inner cells: mean-value form about the cell centre (the S1 hull of a cell is too wide near the tip, |p| ~ 2.5):
                #   |D(xi)| <= |D(xi_c)| + |xi - xi_c| sup|D'|,  D = c - S1(x),  |D'| <= |th~ - th_S1| + |1/lam - 1| + (|k~ e| + |e'|)/lam
                qc = g.at(i, arb(Sb.mid())); xcb = qc['x']
                cc = S.chart_ball(xcb)
                Dc = (qc['Gam'] - XC) / LAM - cc['G']
                d0 = U(abs(Dc)) + U(e0 / LAM) * 1.5          # the e-disk (|re|, |im| <= e0/lam), modulus <= sqrt(2) e0/lam
                Dp = U(abs(q['th'] - th1)) + U(abs(1 / LAM - 1)) + U((kt * e0 + e1) / LAM)
                hw = U(abs(q['x'] - xcb))
                delta = math.nextafter(d0 + hw * Dp, math.inf) * (1 + 1e-12)
            else:
                delta = U(abs(cpt - G1))
            K = s1_kmax(xl - 2 * delta, xu + 2 * delta)
            q1 = 3 * K * delta; ang = eta + 2 * K * delta
            assert ang < math.pi / 2, ('angle', jt, i, j, ang)
            q2 = math.nextafter((1 + q1) * math.tan(ang), math.inf) * (1 + 1e-12)
            assert q1 < 0.5 and q2 < 1.0, ('cell fails', jt, i, j, xl, xu, delta, eta, K, q1, q2)
            bw['q1'] = max(bw['q1'], q1); bw['q2'] = max(bw['q2'], q2)
            if xu <= 1.0: worst['delta_tip'] = max(worst['delta_tip'], delta)
            worst['eta'] = max(worst['eta'], eta)
            if prev_hi is not None: assert xl <= prev_hi + 1e-9, ('cells not contiguous in x', xl, prev_hi)
            prev_hi = xu
            cov = max(cov, xl - 2 * delta); cells += 1
            if cov >= XS + 1.0: break
        if cov >= XS + 1.0: break
    assert cov >= XS, ('coverage short of 60.25', cov)
    # tip cell contains xi = 0: element 0, s = -1 is the base tip (x = 0, theta = -pi/2), asserted
    q0 = g.at(0, arb(-1)); assert float(abs(q0['x']).upper()) < 1e-6 and float(abs(q0['th'] + arb.pi() / 2).upper()) < 1e-6
    worst['q1'] = max(worst['q1'], bw['q1']); worst['q2'] = max(worst['q2'], bw['q2'])
    balls.append(dict(tc=[float(TC.lower()), float(TC.upper())], cells=cells, sup_ku=bw['q1'], sup_up=bw['q2'], covered_to=cov))
    print(f'  t_c ball {jt}: {cells} cells, sup|ku| <= {bw["q1"]:.3e}, sup|u\'| <= {bw["q2"]:.3e}, covered to x >= {cov:.3f} ({time.time() - t0:.0f}s)', flush=True)
ok = worst['q1'] < 0.5 and worst['q2'] < 1.0
print('GS at tau_1 on |x| <= 60.25:', worst, 'ok', ok)
me = os.path.abspath(__file__)
out = dict(label='[cert]', ok=bool(ok), XS=XS, sup_ku=worst['q1'], sup_uprime=worst['q2'], sup_delta_tip_x_le_1=worst['delta_tip'], sup_eta=worst['eta'],
           balls=balls, frame=dict(xc=[float(XC.lower()), float(XC.upper())], lam=[float(LAM.lower()), float(LAM.upper())]),
           sqrtF0_slab180=U(SF0),
           statement='for every t_c in [-1,1] of slab 180 and every S1 in the certified set, Lemma T holds cell by cell on the true tip arc up to x >= 60.25: (GS) with sup|k u| <= sup_ku < 1/2, sup|u\'| <= sup_uprime < 1',
           args=dict(argv=sys.argv, env={k_: os.environ[k_] for k_ in ('S1PROF', 'COHERENT') if k_ in os.environ}, python=platform.python_version(),
                     NS=NS, NTC=NTC, Ybar_box=str(YB), script=dict(path=me, sha256=sha(me)), inputs=inputs, cpu_s=None))   # v38: wall time printed to the stage log only (JSON reproducible byte for byte)
print(f'cpu_s {time.time() - t0:.1f}')
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/iface/GS_tau1.json', 'w'), indent=1)
