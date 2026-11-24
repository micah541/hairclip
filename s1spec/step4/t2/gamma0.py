"""T2: certified initial curve Gamma0 for the CDF blow-up proof.
Rerun:   cd /home/claude/s1spec/step4/t2 && python3 -I gamma0.py build      # coefficients + closure + embeddedness (writes gamma0_coeffs.json)
         cd /home/claude/s1spec/step4/t2 && python3 -I gamma0.py sanity     # float run of hls.py from Gamma0 for dt=0.009, compare p2_t0.01
         cd /home/claude/s1spec/step4/t2 && python3 -I gamma0.py final      # (T3) same curve, b_2 re-closed by arb Newton (exact 192-bit
                                                                            # centre) + Krawczyk -> tight interval; re-certified; gamma0_final.json
         G0SNAP=gamma0s_src.npz python3 -I gamma0.py final   # (T3, current Gamma0) the same for the curve with modes 186..1000 slaved
                                                            # (t3/run_proto.py slave); G0FINAL=... sets the output file
Convention (source: feas2/snaps/p2_t0.001.npz, t=0.001; alpha=0 already at an axis crossing, no shift):
  theta(alpha) = -pi/2 + phi(alpha),  phi = alpha + sum_{n=1..M} b_n sin(n alpha)     (orientation: theta(0) = -pi/2, i.e. the
  opposite of the pi/2 convention in the task; the half-curve alpha in [0,pi] lies in the LOWER half plane y<0).
  X(alpha) = X0 + (L/2pi) int_0^alpha e^{i theta},  X0 = (x0, 0) exact.  e^{i theta} = -i e^{i phi}:
  Re X = x0 + (L/2pi) int sin(phi),  Im X = -(L/2pi) int cos(phi).  Mirror axis = real axis.  Closure: C(b)=int_0^pi cos(phi)=0.
All certified quantities in arb (python-flint); Taylor models (degree K, even) with rigorous remainder per arc.
"""
import sys, os, json, time, numpy as np
from flint import arb, acb, arb_mat, ctx
ctx.prec = 192
K = 24                       # even Taylor degree
SNAP = os.environ.get('G0SNAP', '../feas2/snaps/p2_t0.001.npz')   # source state (npz: th_hat, N, L, X0, t)
FINAL = os.environ.get('G0FINAL', 'gamma0_final.json')
JF = 'gamma0_coeffs.json'
JJ = 2                       # closing coefficient b_2
TAU = arb(10) ** -22         # per-arc quadrature tolerance
PI = arb.pi()
FAC = [arb.fac_ui(k) for k in range(K + 2)]

def load():
    d = np.load(SNAP); h = d['th_hat']; N = int(d['N'])
    b = -2 * h.imag / N
    ev = 2 * h.real / N
    tail = lambda M: float((np.arange(M + 1, len(b)) ** 4 * np.abs(b[M + 1:])).sum())
    M = next(M for M in range(100, len(b)) if tail(M) <= 1e-12)
    return dict(b=b, M=M, tail=tail(M), even=float(np.abs(ev[1:]).max()), even_all=float(np.abs(ev[1:]).sum()),
                const=float(h[0].real / N + np.pi / 2), L=float(d['L']), x0=float(d['X0'][0]), y0=float(d['X0'][1]), t=float(d['t']))

# ---------------- series machinery
def Pmat(bs):
    M = len(bs)
    ent = []
    for k in range(K + 1):
        ent += [bs[n - 1] * arb(n ** k) / FAC[k] for n in range(1, M + 1)]
    return arb_mat(K + 1, M, ent)

def trig(c, M):
    ent = []
    for n in range(1, M + 1):
        s, co = (arb(n) * c).sin_cos(); ent += [s, co]
    return arb_mat(M, 2, ent)

def sin_series(P, c, M):
    """returns phi_k (k=0..K) for phi = alpha + sum b_n sin(n alpha) expanded at c (point or ball)"""
    A = P * trig(c, M); out = []
    for k in range(K + 1):
        r = k % 4; v = A[k, 0] if r == 0 else A[k, 1] if r == 1 else -A[k, 0] if r == 2 else -A[k, 1]
        out.append(v)
    out[0] = out[0] + c; out[1] = out[1] + 1
    return out

def exp_i(phi):
    s, co = phi[0].sin_cos(); e = [acb(co, s)]
    a = [acb(0, p) for p in phi]
    for k in range(1, K + 1):
        t = acb(0)
        for j in range(1, k + 1): t += j * a[j] * e[k - j]
        e.append(t / k)
    return e

def jsin_series(j, c):
    s, co = (arb(j) * c).sin_cos(); out = []
    for k in range(K + 1):
        r = k % 4; v = s if r == 0 else co if r == 1 else -s if r == 2 else -co
        out.append(v * arb(j ** k) / FAC[k])
    return out

def mul(p, q):
    return [sum((p[i] * q[k - i] for i in range(k + 1)), acb(0)) for k in range(K + 1)]

def Gint(co, boxK, xa, x):
    """int_{xa}^{x} of Taylor model (K even)"""
    s = acb(0)
    for k in range(K):
        s += co[k] * (x ** (k + 1) - xa ** (k + 1)) / (k + 1)
    return s + boxK * (x ** (K + 1) - xa ** (K + 1)) / (K + 1)

class Arc:
    """Taylor model of phi and E=e^{i phi} on alpha in [pi*ta, pi*tb] ; optionally times sin(j alpha)"""
    def __init__(s, ta, tb, P, M, jj=None):
        s.ta, s.tb = ta, tb
        a, b = PI * ta, PI * tb
        cm = ((a + b) / 2).mid()
        s.xa, s.xb = a - cm, b - cm
        s.xr = s.xa.union(s.xb)
        I = a.union(b)
        s.phi = sin_series(P, cm, M); phiI = sin_series(P, I, M)
        s.phiK = phiI[K]
        E = exp_i(s.phi); EI = exp_i(phiI)
        if jj is not None:
            E = mul(E, [acb(v) for v in jsin_series(jj, cm)]); EI = mul(EI, [acb(v) for v in jsin_series(jj, I)])
        s.E = E[:K]; s.EK = EI[K]
    def G(s, x=None): return Gint(s.E, s.EK, s.xa, s.xb if x is None else x)
    def phirange(s):
        x = s.xr; v = sum((s.phi[k] * x ** k for k in range(K)), arb(0)) + s.phiK * x ** K
        return v

def quad(P, M, leaves, jj=None, tau=TAU):
    """adaptive validated int over leaves (list of (ta,tb)); returns (G total acb, refined leaves)"""
    tot = acb(0); out = []
    stack = list(reversed(leaves))
    while stack:
        ta, tb = stack.pop(); A = Arc(ta, tb, P, M, jj); g = A.G()
        if g.real.is_finite() and g.imag.is_finite() and g.real.rad() < tau and g.imag.rad() < tau:
            tot += g; out.append((ta, tb))
        else:
            tm = (ta + tb) / 2; stack.append((tm, tb)); stack.append((ta, tm))
    return tot, out

def uniform(n): return [(i / n, (i + 1) / n) for i in range(n)]

def coeff_arb(D, M, bj):
    bs = [arb(float(x)) for x in D['b'][1:M + 1]]; bs[JJ - 1] = bj
    return bs

def Cfun(D, M, bj, leaves):
    P = Pmat(coeff_arb(D, M, bj)); G, lv = quad(P, M, leaves)
    return G.real, lv

def dCfun(D, M, bj, leaves, tau=TAU):
    P = Pmat(coeff_arb(D, M, bj)); G, lv = quad(P, M, leaves, jj=JJ, tau=tau)
    return -G.imag, lv

# ---------------- main certificates
def build(final=False):
    t00 = time.time(); D = load(); M = D['M']; L = D['L']; log = {}
    print('M', M, 'tail', D['tail'], 'even part max', D['even'], 'const dev', D['const'], flush=True)
    bj0 = float(D['b'][JJ])
    # closure: float Newton with validated C
    t0 = time.time(); leaves = uniform(64); m = bj0; cp = None
    for it in range(4):
        C, leaves = Cfun(D, M, arb(m), leaves)
        if cp is None:
            dC, _ = dCfun(D, M, arb(m), leaves); cp = float(dC.mid())
        print('  newton', it, 'm', repr(m), 'C(m)', C, "C'", cp, 'leaves', len(leaves), round(time.time() - t0), 's', flush=True)
        m = m - float(C.mid()) / cp
    ma = arb(m)
    if final:   # continue Newton in arb with exact (rounded-midpoint) iterates: centre no longer restricted to doubles
        for it in range(3):
            C, leaves = Cfun(D, M, ma, leaves)
            ma = arb((ma - C / arb(cp)).mid())
            print('  arb newton', it, 'C', C, flush=True)
    C, leaves = Cfun(D, M, ma, leaves)
    Cm = C; rho = 2 * abs(float(Cm.mid()) / cp) + 4 * float(Cm.rad()) / abs(cp) + 1e-300
    rho = float(np.nextafter(rho * 1.0000001, 1))
    X = ma + arb(0, rho)
    dCX, _ = dCfun(D, M, X, leaves, tau=arb(10) ** -9)
    cpa = arb(cp)
    Kr = ma - Cm / cpa + (1 - dCX / cpa) * (X - ma)
    ok = X.contains(Kr.mid()) and Kr.upper() < X.upper() and Kr.lower() > X.lower()
    contr = abs(float((1 - dCX / cpa).abs_upper()))
    print('Krawczyk: m', repr(m), 'rho', rho, "C(m)", Cm, "|1-C'(X)/c'|<=", contr, 'K(X) in int X:', ok, flush=True)
    assert ok
    lo, hi = X.lower(), X.upper(); bjX = X
    # embeddedness with b_j in [lo,hi]
    P = Pmat(coeff_arb(D, M, bjX)); Lf = arb(L); k2 = Lf / (2 * PI)
    x0 = arb(D['x0'])
    cache = {}
    nodes = sorted(set([t for lf in leaves for t in lf]))
    def refine_phi(nodes):  # enforce turning <= 0.05 later via build
        return nodes
    def build_arcs(nodes):
        arcs = []; X = acb(x0, 0)
        for ta, tb in zip(nodes[:-1], nodes[1:]):
            key = (ta, tb)
            if key not in cache: cache[key] = Arc(ta, tb, P, M)
            A = cache[key]; G = A.G(); Gx = A.G(A.xr)
            bx = acb(X.real + k2 * Gx.imag, X.imag - k2 * Gx.real)
            Xn = acb(X.real + k2 * G.imag, X.imag - k2 * G.real)
            pr = A.phirange()
            arcs.append(dict(ta=ta, tb=tb, box=bx, pr=pr, G=G, rad=max(float(G.real.rad()), float(G.imag.rad()))))
            X = Xn
        return arcs, X
    def need_split(nodes):
        arcs, Xend = build_arcs(nodes); new = set();
        for a in arcs:
            w = float(a['pr'].upper() - a['pr'].lower()) if a['pr'].is_finite() else 9.
            if w > 0.05 or a['rad'] > 1e-13:
                new.add((a['ta'], a['tb']))
        return arcs, Xend, new
    t1 = time.time()
    for rnd in range(12):
        arcs, Xend, bad = need_split(nodes)
        # box disjointness for non-adjacent
        n = len(arcs)
        def dbl(arbs, lowside):
            return np.array([float(v.lower() if lowside else v.upper()) for v in arbs])
        xlo = np.nextafter(dbl([a['box'].real for a in arcs], True), -np.inf); xhi = np.nextafter(dbl([a['box'].real for a in arcs], False), np.inf)
        ylo = np.nextafter(dbl([a['box'].imag for a in arcs], True), -np.inf); yhi = np.nextafter(dbl([a['box'].imag for a in arcs], False), np.inf)
        inter = set(); gapmin = np.inf
        for i in range(n):
            js = np.arange(i + 2, n)
            if len(js) == 0: continue
            gx = np.maximum(xlo[js] - xhi[i], xlo[i] - xhi[js]); gy = np.maximum(ylo[js] - yhi[i], ylo[i] - yhi[js])
            g = np.maximum(gx, gy)
            hit = js[g <= 0]
            for j in hit: inter.add(i); inter.add(int(j))
            if (g > 0).any(): gapmin = min(gapmin, g[g > 0].min())
        print(' round', rnd, 'arcs', n, 'turning/precision violations', len(bad), 'intersecting arcs', len(inter), 'min positive sep', gapmin, round(time.time() - t1), 's', flush=True)
        if not bad and not inter: break
        split = set(bad) | {(arcs[i]['ta'], arcs[i]['tb']) for i in inter}
        ns = set(nodes)
        for ta, tb in split: ns.add((ta + tb) / 2)
        nodes = sorted(ns)
    else:
        raise SystemExit('embedding refinement did not converge')
    # ---- adjacent pair turning, zones, y<0
    prl = [float(a['pr'].lower()) for a in arcs]; pru = [float(a['pr'].upper()) for a in arcs]
    maxturn = max(u - l for l, u in zip(prl, pru))
    maxpair = max(max(pru[i], pru[i + 1]) - min(prl[i], prl[i + 1]) for i in range(n - 1))
    assert maxpair < np.pi
    q0 = 0
    while q0 < n and prl[q0] > -1.5 and pru[q0] < 1.5: q0 += 1
    q1 = n
    while q1 > 0 and prl[q1 - 1] > 1.6 and pru[q1 - 1] < 4.7: q1 -= 1
    assert q0 >= 1 and q1 <= n - 1 and q0 <= q1
    mid_ok = all(arcs[i]['box'].imag.upper() < 0 for i in range(q0, q1))
    y0z = float(arcs[q0]['box'].imag.upper())
    # zone check for local argument: phi range of arcs in zones strictly inside cos>0 / cos<0 windows (done above by q0,q1 windows)
    alpha_mid = lambda a: float(PI * (a['ta'] + a['tb']) / 2)
    slot = [a for a in arcs[q0:q1] if 1.75 <= alpha_mid(a) <= 2.9]
    gap_cert = min(float(-a['box'].imag.upper()) for a in slot)       # lower bound for |y| on slot plateau
    # float counterpart
    b = D['b'][1:M + 1].copy(); b[JJ - 1] = m
    n_ = np.arange(1, M + 1)
    # y(al) = y(1.75) - integral; simpler: use dense cumulative over [0,pi]
    A2 = np.linspace(0, np.pi, 200001)
    ph2 = np.empty_like(A2)
    for s_ in range(0, len(A2), 4000): ph2[s_:s_ + 4000] = A2[s_:s_ + 4000] + np.sin(np.outer(A2[s_:s_ + 4000], n_)) @ b
    da = A2[1] - A2[0]; yf = -(L / 2 / np.pi) * np.concatenate([[0], np.cumsum((np.cos(ph2[1:]) + np.cos(ph2[:-1])) / 2 * da)])
    msk = (A2 >= 1.75) & (A2 <= 2.9); gap_float = float(-yf[msk].max())
    ends = (Xend.real, Xend.imag)
    print('arcs', n, 'zones q0', q0, 'q1', q1, 'middle y<0:', mid_ok, 'max arc turning', maxturn, 'max adjacent-pair turning', maxpair, '(<pi)')
    print('X(pi)', Xend, '| min |y| slot plateau (cert lower bd)', gap_cert, 'float', gap_float, flush=True)
    assert mid_ok
    # right zone cos<0 argument needs C=0 (closure) -> y(pi)=0 ; left zone: y(0)=0 exactly
    # JSON
    bhex = [float(x).hex() for x in D['b'][1:M + 1]]
    bhex[JJ - 1] = float(m).hex()
    if final:
        man, ex = ma.mid().man_exp()
        bhex[JJ - 1] = float(ma).hex()
        lo_s, hi_s = (X - X.mid()).lower(), (X - X.mid()).upper()
        json.dump(dict(source=SNAP, t=D['t'], M=M, j=JJ, b_hex=bhex,
                       b_j_center_exact=dict(man=str(man), exp=str(ex), note='b_j centre = man * 2^exp exactly (b_hex[j-1] is its nearest double)'),
                       b_j_halfwidth=repr(rho), b_j_closing_interval=[X.lower().str(45, radius=False), X.upper().str(45, radius=False)],
                       Cm=str(Cm), krawczyk_contraction=contr, L0_hex=float(L).hex(), L0=L, X0_hex=[float(D['x0']).hex(), (0.0).hex()],
                       convention='theta=-pi/2+phi, phi=alpha+sum b_n sin(n alpha), alpha=0 axis crossing, theta(0)=-pi/2, half-curve in y<0; closure int_0^pi cos(phi)=0; X0 on real axis; all b_n exact doubles except b_j (exact dyadic above)',
                       embedded=dict(arcs=n, gap_cert=gap_cert, gap_float=gap_float, maxturn=maxturn, maxpair=maxpair, minsep=float(gapmin), y_left_zone_end=y0z),
                       nodes_t=[float(t).hex() for t in nodes]), open(FINAL, 'w'), indent=1)
        print('wrote', FINAL, ': b_2 halfwidth', rho, flush=True)
        return
    json.dump(dict(source=SNAP, t=D['t'], M=M, j=JJ, b_hex=bhex, b_j_closing_interval=[lo.str(30, radius=False), hi.str(30, radius=False)],
                   b_j_center_hex=float(m).hex(), b_j_halfwidth_hex=float(rho).hex(), L0_hex=float(L).hex(), L0=L, X0_hex=[float(D['x0']).hex(), (0.0).hex()],
                   convention='theta=-pi/2+phi, phi=alpha+sum b_n sin(n alpha), alpha=0 axis crossing, theta(0)=-pi/2, half-curve in y<0; closure int_0^pi cos(phi)=0; X0 on real axis',
                   nodes_t=[float(t).hex() for t in nodes]), open(JF, 'w'), indent=1)
    res = dict(M=M, tail=D['tail'], even=D['even'], rho=rho, contr=contr, arcs=n, q0=q0, q1=q1, gap_cert=gap_cert, gap_float=gap_float,
               maxturn=maxturn, maxpair=maxpair, minsep=float(gapmin), Xend_imag=str(Xend.imag), Xend_real=str(Xend.real), time=time.time() - t00, closure_leaves=len(leaves),
               bj=repr(m), Cm=str(Cm), y_left_zone_end=y0z)
    json.dump(res, open('results.json', 'w'), indent=1); print(res)

def sanity():
    sys.path.insert(0, '/home/claude/s1spec/step4')
    from hls import HLS
    J = json.load(open(JF)); b = np.array([float.fromhex(h) for h in J['b_hex']]); M = len(b); L = float.fromhex(J['L0_hex'])
    N = 8192; al = 2 * np.pi * np.arange(N) / N
    S = np.zeros(N)
    for s_ in range(0, N, 1024): S[s_:s_ + 1024] = np.sin(np.outer(al[s_:s_ + 1024], np.arange(1, M + 1))) @ b
    th = -np.pi / 2 + S
    H = HLS(th, L, [float.fromhex(J['X0_hex'][0]), 0.0], omega=1, cdt=1e-2, tail_tol=1e-10, Nmax=1 << 21)
    t0 = time.time(); H.run_to(tend=0.009); print('ran to', H.t, 'in', round(time.time() - t0), 's N', H.N, flush=True)
    z = np.load('../feas2/snaps/p2_t0.01.npz'); R = HLS(np.zeros(int(z['N'])), float(z['L']), z['X0'], omega=1, cdt=1e-2, tail_tol=1e-10, Nmax=1 << 21)
    R.th_hat = z['th_hat']; R.t = float(z['t'])
    Nc = max(H.N, R.N); tg = lambda S_: np.fft.irfft(np.pad(S_.th_hat, (0, Nc // 2 + 1 - len(S_.th_hat))) * (Nc / S_.N), Nc)
    d = tg(H) - tg(R)
    print('t_run', H.t, 't_ref', R.t, 'max|dth|', np.abs(d).max(), 'dL', H.L - R.L, 'L', R.L, 'dX0', H.X0 - R.X0, 'kmax', H.kmax(), R.kmax())

if __name__ == '__main__':
    {'build': build, 'sanity': sanity, 'final': lambda: build(final=True)}[sys.argv[1]]()
