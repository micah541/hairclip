"""run_proto.py -- T3 prototype driver.
  python3 -I run_proto.py one   [c]          one slab at t=2 (p2_t2.0, M=360), h = c/max|theta_a|^4
  python3 -I run_proto.py chain [c] [n]      n consecutive slabs at t=2
  python3 -I run_proto.py gamma0 [c] [M]     one slab from exactly Gamma0 (t2/gamma0_coeffs.json, re-truncated to M modes)
  python3 -I run_proto.py gamma0c [c] [M]    same regime, Gamma0 made consistent first (relax + 2 preliminary arb slabs)
  python3 -I run_proto.py initial [n]        n graded slabs from the certified Gamma0 (t2/gamma0_final.json), dt = t/2, K=24
  python3 -I run_proto.py swing [n]          block-resizing test t=4 -> 4.59 (theta 16 -> 40), M=480, c=2
  python3 -I run_proto.py approach [n]       n slabs at k~12 (hls-regenerated resolved state), M=1000, k_max at the end
  python3 -I run_proto.py smoke              fast check (M=200): block add/drop + quadratures + k_max
Outputs: logs/*.json (per-slab info), printed summary.
"""
import sys, os, json, time, math
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import slab as S

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = '/home/claude/s1spec/step4/feas2/snaps'
os.makedirs(f'{HERE}/logs', exist_ok=True)


def load_sine(path, M):
    z = np.load(path)
    N = int(z['N']); th = z['th_hat']
    b = -2 * th.imag / N
    cos_drop = np.abs(2 * th.real[1:] / N).max()
    return b[1:M + 1].copy(), float(z['L']), float(z['t']), cos_drop, b[M + 1:]


def etd_coeffs(L, k):
    Mc = 32; r = np.exp(1j * np.pi * (np.arange(1, Mc + 1) - 0.5) / Mc); LR = k * L[:, None] + r[None, :]
    E = np.exp(k * L); E2 = np.exp(k * L / 2)
    Q = k * np.real(np.mean((np.exp(LR / 2) - 1) / LR, axis=1))
    f1 = k * np.real(np.mean((-4 - LR + np.exp(LR) * (4 - 3 * LR + LR ** 2)) / LR ** 3, axis=1))
    f2 = k * np.real(np.mean((2 + LR + np.exp(LR) * (-2 + LR)) / LR ** 3, axis=1))
    f3 = k * np.real(np.mean((-4 - 3 * LR - LR ** 2 + np.exp(LR) * (4 - LR)) / LR ** 3, axis=1))
    return E, E2, Q, f1, f2, f3


def relax(fm, b, tau, nsteps):
    """float ETDRK4 of the sine-class Galerkin system (own integrator): makes the datum consistent
    with the M-mode Galerkin slow manifold (stiff-mode initial layer of the projected snapshot)."""
    k = tau / nsteps; E, E2, Q, f1, f2, f3 = etd_coeffs(-fm.lam, k); u = b.copy(); F = lambda x: fm.N(x[None])[0]
    for _ in range(nsteps):
        Nu = F(u); a = E2 * u + Q * Nu; Na = F(a); bb = E2 * u + Q * Na; Nb = F(bb); c = E2 * a + Q * (2 * Nb - Nu); Nc = F(c)
        u = E * u + f1 * Nu + 2 * f2 * (Na + Nb) + f3 * Nc
    return u


def jsonable(info):
    out = {}
    for k, v in info.items():
        if isinstance(v, np.ndarray): continue
        if isinstance(v, (np.floating, np.integer)): v = v.item()
        out[k] = v
    return out


def summary(info, L):
    dtdtau = (L / (2 * np.pi)) ** 4
    s = (f"h={info['h']:.3e} (dt~{info['h'] * dtdtau:.3e})  M={info['M']} Ma={info['Ma']} Nc={info['Nc']} K={info['K']} "
         f"KP={info['KP']} KJ={info['KJ']} msub={info['msub']}\n"
         f"  residual: max_n sup|r_n| = {info['resid_max']:.2e}, sum n^3 sup|r_n| = {info['resid_w3']:.2e}, "
         f"induced sum n^3 r_n/max(lam_n,1/h) = {info['resid_ind']:.2e}; defect*w |.|_3 = {info['defect_w3']:.2e}\n"
         f"  polish |F| = {['%.1e' % x for x in info['polish']]}  u(0)-xbar = {info['umis']:.1e}  cond(P_a(h)) = {info['cond_Pa']:.1e}\n"
         f"  rho radii max: {info['w_in']:.3e} -> {info['w_out']:.3e} (|diag B|w {info['wB_diag']:.3e}, |B|w {info['wB_full']:.3e}, absorbed {info['w_eta']:.2e}, intra-slab {info['eta_intra']:.2e});  sv(P_a(h)) {['%.3g' % x for x in info['sv_Pa']]}\n"
         f"  parallelepiped hull sum_n n^j|G||w|: j=1 {info['hull1_in']:.3e} -> {info['hull1_out']:.3e}; j=3 {info['hull3_in']:.3e} -> {info['hull3_out']:.3e};  box |.|_1 out {info['box1_out']:.2e}\n"
         f"  box |.|_3: {info['box_in_w3']:.3e} -> {info['box_out_w3']:.3e}  (Zend block {info['Zend_a_w3']:.2e}, rest {info['Zend_r_w3']:.2e}); "
         f"tail T: {info['T_in']:.2e} -> {info['T_out']:.2e}\n"
         f"  slab-wide sup |e|_j (j=0..3): {['%.2e' % x for x in info['supnorm']]};  int <u_aa^2> = {info['Lint'][0]:.15e} +/- {info['Lint'][1]:.1e}\n"
         f"  time {info['time']:.1f}s  {{{', '.join('%s: %.1f' % kv for kv in info['tm'].items())}}}  sub-intervals {info['sub']}")
    return s


def run_one(c=2.0, nslab=1, M=360, Ma=None, eps0=1e-20, tag="t2", msub=8, npre=2, lam_slow=8.0, fixed_h=True):
    b, L, t0, cos_drop, rest = load_sine(f'{SNAP}/p2_t2.0.npz', M)
    print(f'p2_t2.0: t={t0} L={L:.6f}; dropped cosine part max {cos_drop:.1e}; dropped sine modes > {M}: max {np.abs(rest).max():.1e}')
    fm = S.FM(M)
    th4 = fm.thmax(b[None])[0] ** 4
    t0w = time.time(); taur = 4 * c / th4
    b = relax(fm, b, taur, 16000)
    b = relax(fm, b, taur / 16, 16000)          # fine phase: stiff-mode slaving error ~ (k theta^4)^p
    print(f'relaxed datum: own ETDRK4 over tau_r = {taur:.3e} (dt ~ {taur * (L / 2 / np.pi) ** 4:.2e}), 16000 steps, {time.time() - t0w:.0f}s')
    h0 = float(np.float64(c / th4)); xa = [S.arb(float(v)) for v in b]
    for k in range(npre):
        t0w = time.time(); xprev = np.array([float(x.mid()) for x in xa]); xa, pol = S.prelim_slab(xa, fm, h0)
        print(f'preliminary arb base slab {k + 1}: polish {["%.1e" % p for p in pol]}  ({time.time() - t0w:.1f}s)')
    b = np.array([float(x.mid()) for x in xa])
    if Ma is None: Ma = int(np.sum(fm.lam * float(np.float64(c / fm.thmax(b[None])[0] ** 4)) <= lam_slow))
    S0 = S.slaving_propagated(xprev, fm, h0, Ma)
    print(f'parallelepiped block: Ma = {Ma} slow modes (lam_n h <= {lam_slow}); initial generators slaved by a float propagator over the last preliminary slab')
    st = S.initial_state(b, Ma, 4 * M, eps0, fm, S0); st.xbar_arb = xa
    log = []
    Ug = None
    hfix = float(np.float64(c / fm.thmax(st.xbar[None])[0] ** 4))
    for k in range(nslab):
        th4 = fm.thmax(st.xbar[None])[0] ** 4
        h = hfix if fixed_h else float(np.float64(c / th4))       # fixed h keeps lam_{Ma} h (block conditioning) constant
        st_new, info = S.validate_slab(st, fm, h, Ma=Ma, Uguess=Ug, verbose=False, msub=msub)
        print(f'--- slab {k + 1}/{nslab}  (theta_max = {th4 ** 0.25:.3f}, tau = {st_new.tau:.4e})')
        print(summary(info, L))
        rec = jsonable(info); rec['slab'] = k + 1; rec['tau_end'] = st_new.tau; rec['thmax'] = th4 ** 0.25
        log.append(rec)
        # next guess: constant at the new centre
        Ug = None
        st = st_new
        with open(f'{HERE}/logs/{tag}_c{c}_n{nslab}.json', 'w') as f: json.dump(log, f, indent=1)
    return st, log


def run_gamma0(c=2.0, M=1000, Ma=None, eps0=1e-20, lam_slow=8.0, consistent=False):
    """one validated slab starting EXACTLY at Gamma0 (t2/gamma0_coeffs.json), re-truncated to M modes."""
    g = json.load(open('/home/claude/s1spec/step4/t2/gamma0_coeffs.json'))
    bfull = np.array([float.fromhex(x) for x in g['b_hex']])
    j = int(g['j']); hw = float.fromhex(g['b_j_halfwidth_hex'])
    L0 = float.fromhex(g['L0_hex'])
    nn = np.arange(1, len(bfull) + 1)
    print(f"Gamma0: {len(bfull)} coefficients, closing b_{j} = {bfull[j - 1]!r} +/- {hw:.3e}; L0={L0}")
    print(f"  dropped modes > {M}: max |b_n| {np.abs(bfull[M:]).max():.2e}, sum n^3|b_n| = {np.sum(nn[M:] ** 3 * np.abs(bfull[M:])):.2e}")
    b = bfull[:M].copy()
    fm = S.FM(M)
    th = fm.thmax(b[None])[0]
    h = float(np.float64(c / th ** 4))
    if Ma is None: Ma = int(np.sum(fm.lam * h <= lam_slow))
    Nc = 4 * M
    xa = None; S0p = None
    if consistent:
        # Gamma0 (hls t=0.001) is still inside the corner-clip initial layer in modes ~70-300 (O(1e-1) relative transient):
        # make it consistent with the M-mode sine Galerkin flow (own ETDRK4), then 2 preliminary arb base slabs.
        t1 = time.time(); b = relax(fm, b, 2 * h, 8000); b = relax(fm, b, h / 8, 8000)
        xa = [S.arb(float(v)) for v in b]
        for k in range(2):
            xprev = np.array([float(x.mid()) for x in xa]); xa, pol = S.prelim_slab(xa, fm, h)
        b = np.array([float(x.mid()) for x in xa])
        S0p = S.slaving_propagated(xprev, fm, h, Ma)
        print(f'  consistent datum: relax 2h+h/8 + 2 preliminary arb slabs ({time.time() - t1:.0f}s; datum at t ~ {0.001 + 4.125 * h * (L0 / 2 / np.pi) ** 4:.3f})')
    st = S.initial_state(b, Ma, Nc, eps0, fm, S0p)        # generators [I; S] (S quasi-static if no history)
    if xa is not None: st.xbar_arb = xa
    S0 = st.G[Ma:]
    w = st.w.copy(); w[j - 1] += hw * (1 + 1e-12)          # closing-coefficient interval -> parallelepiped
    box = st.box.copy(); box[Ma:M] = eps0 * np.arange(Ma + 1, M + 1, dtype=float) ** (-S.QT) + np.abs(S0) @ w
    m2 = min(len(bfull), Nc); box[M:m2] += np.abs(bfull[M:m2]); box *= (1 + 1e-12)
    st.w = w; st.box = box
    print(f'  M={M} Ma={Ma} Nc={Nc} theta_max={th:.3f} h={h:.3e} (dt~{h * (L0 / 2 / np.pi) ** 4:.2e})')
    t0 = time.time()
    st_new, info = S.validate_slab(st, fm, h, Ma=Ma, verbose=False)
    print(summary(info, L0))
    with open(f'{HERE}/logs/gamma0{"c" if consistent else ""}_c{c}_M{M}.json', 'w') as f: json.dump([jsonable(info)], f, indent=1)
    return st_new, info


T2DIR = '/home/claude/s1spec/step4/t2'


def qline(info, st):
    q = info.get('quad') or {}
    s = (f"  block {info['d']}->{info['a1']} (passengers {info['passengers']}), h_next={info['h_next']}, theta_end={info['theta_end']:.3f}\n"
         f"  quadratures: t={q.get('t')}  L={q.get('L')}  X0={q.get('X0')}")
    if 'kmax' in q: s += f"\n  k_max(s=h) in [{q['kmax'][0]:.12f}, {q['kmax'][1]:.12f}] at alpha={q['alpha_kmax']:.5f} (|e|_1 <= {q['e1_end']:.1e})"
    return s


def run_chain(st, fm, nslab, hrule, tag, K=16, kmax_last=0, stop=None, Mrule=None, Lfloat=None, logf=None, adapt=True, lamcut1=None, etacap1=None):
    """generic validated chain: h from hrule(x_float, t_float, L_float); next block size Ma = #(lam h_next <= 8)"""
    log = []
    h = hrule(st.xbar, float(st.t.mid()), float(st.L.mid()))
    rreq = 10 ** (24.0 / K)
    def hcheck(h):
        # adaptive slab length: float base Chebyshev decay ratio r (degrees 1..6) must satisfy r^K >= 1e24
        for it in range(4):
            Ms = int(np.sum(fm.lam * h <= 8.0 * K * K))
            U = S.base_float(fm, st.xbar, h, K, np.tile(st.xbar, (K + 1, 1)), min(max(int(np.sum(fm.lam * h <= 8)), 60), 140), Ms)
            cm = np.abs(S.vinv_f(K) @ U).max(axis=1)
            r = (cm[1] / cm[6]) ** 0.2
            if r >= 1.1 * rreq: return h, r
            h = h * r / (1.25 * rreq)
        return h, r
    for k in range(nslab):
        if adapt:
            h0_ = h; h, rdec = hcheck(h)
            print(f'  adaptive h: {h0_:.3e} -> {h:.3e} (Chebyshev decay ratio {rdec:.1f}, required {rreq:.1f})', flush=True)
        tcur, Lcur = float(st.t.mid()), float(st.L.mid())
        def plan(uh, h=h, tcur=tcur, Lcur=Lcur):
            hn = hrule(uh, tcur + h * (Lcur / 2 / np.pi) ** 4, Lcur)
            return hn, int(np.sum(fm.lam * hn <= 8.0))
        last = (k == nslab - 1)
        # intra-slab absorption cap (relative to max w): first slab as given, then 1.5x the last measured ratio
        # (a cap that is too small costs a full redo of the slab)
        S.ETACAP = (etacap1 or 1e-2) if k == 0 else min(1.0, max(1e-2, 1.5 * info['eta_intra'] / max(info['w_out'], 1e-300)))
        st_new, info = S.validate_slab(st, fm, h, K=K, Ma=st.G.shape[1], verbose=False, plan_next=plan,
                                       kmax_ng=kmax_last if last else 0, lamcut=lamcut1 if k == 0 else None)
        print(f'--- {tag} slab {k + 1}/{nslab}  (M={fm.M}, theta_max start {fm.thmax(st.xbar[None])[0]:.3f}, tau = {st_new.tau:.4e})', flush=True)
        print(summary(info, Lcur)); print(qline(info, st_new), flush=True)
        rec = jsonable(info); rec['slab'] = k + 1; log.append(rec)
        with open(f'{HERE}/logs/{tag}.json', 'w') as f: json.dump(log, f, indent=1)
        st = st_new; h = info['h_next']
        if Mrule is not None:
            Mn = Mrule(st.xbar)
            if Mn < fm.M - 40 or Mn > fm.M:
                st = S.resize_M(st, fm, Mn); fm = S.FM(Mn); print(f'  M resized -> {Mn}', flush=True)
        if stop is not None and stop(info): break
    return st, fm, log


def prep_consistent(fm, b, h0, npre=2):
    """short own-flow relaxation + preliminary arb slabs (all of length <= h0/4) -> exact consistent centre"""
    b = relax(fm, b, h0 / 4, 4000); b = relax(fm, b, h0 / 64, 4000)
    xa = [S.arb(float(v)) for v in b]
    for k in range(npre):
        xprev = np.array([float(x.mid()) for x in xa]); xa, pol = S.prelim_slab(xa, fm, h0 / 4)
    return xa, xprev


def run_slave(M=1000, K=24, hs_fac=1.8, nfree=186, src='gamma0_final_unslaved.json', two_pass=True, lamcut1=26.0):
    """option (c): put the fast modes nfree..M of the certified T2 curve on the slow manifold.  Free-start Lobatto
    collocation (ODE imposed at all K+1 nodes for n >= nfree, initial value only below) over [0, hs], hs = hs_fac*t0 in
    tau' (lam_nfree*hs >= 26: well conditioned), arb-polished.  Writes t2/gamma0s_src.npz (the snapshot with modes
    nfree..M replaced; doubles, exact) for  G0SNAP=gamma0s_src.npz python3 -I gamma0.py final  (re-closure + embedding)."""
    g = json.load(open(f'{T2DIR}/{src}'))
    bfull = np.array([float.fromhex(x) for x in g['b_hex']]); j = int(g['j'])
    bj = S.arb(int(g['b_j_center_exact']['man'])) * S.arb(2) ** int(g['b_j_center_exact']['exp'])
    L0 = float.fromhex(g['L0_hex']); tf0 = float(g['t'])
    T = tf0 / (L0 / 2 / np.pi) ** 4; hs = hs_fac * T
    fm = S.FM(M); b = bfull[:M].copy(); b[j - 1] = float(bj.mid())
    xa = [S.arb(float(v)) for v in b]; xa[j - 1] = bj
    t1 = time.time()
    Ua, C, Crows, pol, U = S.solve_base(fm, xa, b, hs, K, 140, nfree - 1)
    u0 = [sum((Crows[k][n] * (-1) ** k for k in range(K + 1)), S.A0) for n in range(M)]
    bn = np.array([float(x.mid()) for x in u0])
    print(f'slaving solve hs={hs:.3e} (= {hs_fac} t0 in tau), K={K}, free n>={nfree}: lam*hs at nfree = {fm.lam[nfree - 1] * hs:.1f}; '
          f'polish |F| {["%.1e" % p for p in pol]}  ({time.time() - t1:.0f}s)', flush=True)
    for a, e in [(nfree, 250), (250, 350), (350, 450), (450, 600), (600, 800), (800, M)]:
        sl = slice(a - 1, e); dl = np.abs(bn[sl] - b[sl])
        print(f'   n in [{a},{e}]: max|b_new-b| = {dl.max():.2e}  (rel {np.max(dl / (np.abs(b[sl]) + 1e-300)):.1e})  sum n^3|.| = {np.sum(fm.n[sl] ** 3 * dl):.2e}', flush=True)
    if two_pass:
        # second pass: the first validated slab (h1 = t0/2 in tau', free start for lam h1 > lamcut1) reproduces these
        # values exactly, so Gamma0 - u_base(0) is only the rounding to doubles (the one-pass mismatch, 5.5e-18 at
        # n = 250, was coupled into the block by the a priori bound: rho 8e-24 -> 1.3e-20 in slab 1)
        h1 = 0.5 * T; n2 = int(np.sum(fm.lam * h1 <= lamcut1)) + 1
        b2 = b.copy(); b2[nfree - 1:M] = bn[nfree - 1:M]
        xa2 = [S.arb(float(v)) for v in b2]; xa2[j - 1] = bj
        t1 = time.time(); Ua, C, Crows, pol, U = S.solve_base(fm, xa2, b2, h1, K, 140, n2 - 1)
        u0 = [sum((Crows[k][n] * (-1) ** k for k in range(K + 1)), S.A0) for n in range(M)]
        bn2 = np.array([float(x.mid()) for x in u0])
        dl = np.abs(bn2[n2 - 1:] - bn[n2 - 1:])
        print(f'second pass h1={h1:.3e}, free n>={n2}: max change {dl.max():.2e} (at n={n2 + int(np.argmax(dl))}); polish {["%.1e" % p for p in pol]} ({time.time() - t1:.0f}s)', flush=True)
        bn[n2 - 1:] = bn2[n2 - 1:]
    z = np.load(g['source'] if g['source'].startswith('/') else f'{T2DIR}/{g["source"]}'); N = int(z['N']); th = z['th_hat'].copy()
    for n in range(nfree, M + 1):
        th[n] = complex(th[n].real, -bn[n - 1] * N / 2)
    assert np.all(-2 * th.imag[nfree:M + 1] / N == bn[nfree - 1:M])
    np.savez(f'{T2DIR}/gamma0s_src.npz', th_hat=th, N=z['N'], L=z['L'], X0=z['X0'], t=z['t'], kmax=z['kmax'], hprev=z['hprev'])
    print(f'wrote {T2DIR}/gamma0s_src.npz (modes {nfree}..{M} slaved)', flush=True)


def cost_model(Ma, M):
    """1-core seconds per validated slab (fit: t~2 chain, Gamma0 slab M=1000/Ma=185 at 525 s CPU, swing M=480-720)"""
    return 0.645 * Ma * (M / 360) ** 1.35 + 10.0 * (M / 360) ** 1.63


def run_estimate(t0=0.001, rho_g=0.5, c=2.0, ceff_app=2.0, kend=12.0, Tb=4.598865):
    """P1 slab count and 1-core time: graded phase (dt = rho_g t) while t kmax^4 < c/rho_g, then h = c/theta^4,
    swing as measured (logs/swing.json), approach kmax^4 (T-t) = 5.55 -> int kmax^4 dt = 22.2 ln(k2/k1), dt = ceff/kmax^4."""
    import glob
    R = sorted((float(z['t']), float(z['kmax']), float(z['L'])) for z in (np.load(f) for f in glob.glob(f'{SNAP}/p2_t*.npz')))
    t, k, L = (np.array(a) for a in zip(*R))
    tk4 = t * k ** 4; t1 = float(np.exp(np.interp(c / rho_g, tk4, np.log(t))))
    nI = int(np.ceil(np.log(t1 / t0) / np.log(1 + rho_g))); tI = [t0 * (1 + rho_g) ** j for j in range(nI)]
    LI = np.interp(tI, t, L); MaI = [int((Lf / 2 / np.pi) * (8 / rho_g / tj) ** 0.25) for Lf, tj in zip(LI, tI)]
    cI = sum(cost_model(a, 1000) for a in MaI)
    # phase II: t1 -> 4.0 at h = c/theta^4; M: 930 (t<0.3), 700 (<1), 480 (<2), 360 (<4)
    tg = np.concatenate([[t1], t[(t > t1) & (t <= 4.0)]]); kg = np.interp(tg, t, k); Lg = np.interp(tg, t, L)
    nII, cII = 0.0, 0.0
    for a in range(len(tg) - 1):
        n = 0.5 * (kg[a] ** 4 + kg[a + 1] ** 4) * (tg[a + 1] - tg[a]) / c
        th = 0.5 * (kg[a] * Lg[a] + kg[a + 1] * Lg[a + 1]) / 2 / np.pi; tm = 0.5 * (tg[a] + tg[a + 1])
        M = 930 if tm < 0.3 else 700 if tm < 1 else 480 if tm < 2 else 360
        nII += n; cII += n * cost_model((8 / c) ** 0.25 * th, M)
    sw = json.load(open(f'{HERE}/logs/swing.json')); nIII = len(sw); cIII = sum(cost_model(r['d'], r['M']) for r in sw)
    k1 = 5.03; nIV = 22.2 * np.log(kend / k1) / ceff_app; LIV = 48.6
    # M along the approach: max(800, 100 k) (k ~ 13: 1300 for |x_M| ~ 1e-26; M = 1500 clean to k ~ 14.5)
    ks = np.exp(np.linspace(np.log(k1), np.log(kend), 200))
    cIV = nIV * np.mean([cost_model((8 / ceff_app) ** 0.25 * kk * LIV / 2 / np.pi, max(800, 100 * kk)) for kk in ks])
    print(f'phase I   graded dt={rho_g}t, t {t0} -> {t1:.3f}: {nI} slabs, Ma {MaI[0]} -> {MaI[-1]}, M=1000: {cI / 3600:.2f} h')
    print(f'phase II  c={c}, t {t1:.3f} -> 4.0: int k^4 dt = {nII * c:.1f}: {nII:.0f} slabs: {cII / 3600:.2f} h')
    print(f'phase III swing (measured, K=16): {nIII} slabs, theta 16 -> {sw[-1]["theta_end"]:.0f}: {cIII / 3600:.2f} h (model, 1 core)')
    print(f'phase IV  approach k {k1} -> {kend}, c_eff={ceff_app} (K=24): {nIV:.0f} slabs, M 800 -> {max(800, 100 * kend):.0f}: {cIV / 3600:.2f} h')
    print(f'total: {nI + nII + nIII + nIV:.0f} slabs, {(cI + cII + cIII + cIV) / 3600:.1f} h on one core')


def run_budget(ckdir, slab, kmarks=(12.0, 20.0), cdt=0.04, t_mid=0.03, srel=1e-6, Mmax=2400):
    """expected P1 error at the switch: float tangent propagation (ETDRK4 in tau') of the parallelepiped A = G diag(w)
    of a run_p1 checkpoint.  All columns to t_mid; then a truncated SVD A(t_mid) = U S V^T (S_i >= srel S_1) and only
    U S onward (the stable directions have mixed by then).  At k_max = kmarks: |e|_j = sum_n n^j sum_cols |(Phi A)_nj|
    (j = 0..4), physical C2 = (2pi/L)|e|_1 + kappa_max |e|_0, C5 = (2pi/L)^4 |e|_4.  The box is negligible (report)."""
    import run_p1 as P
    rec = json.load(open(f'{ckdir}/ckpt/slab_{slab:05d}.json')); z = np.load(f'{ckdir}/ckpt/slab_{slab:05d}.npz')
    u = np.array([float(S.arb(int(m)) * S.arb(2) ** int(e)) for m, e in rec['xbar']])
    V = (z['G'] * z['w'][None, :]).T.copy(); Vt = None                       # rows = columns G_j w_j
    L = float(P.j2a(rec['L']).mid()); t = 0.001 + float(P.j2a(rec['t']).mid()); M = len(u)
    print(f'budget from slab {slab}: t = {t:.6g}, M = {M}, {V.shape[0]} columns, max w {z["w"].max():.3e}', flush=True)
    def hull(V, Vt, M):
        A = V if Vt is None else (V.T @ Vt).T                                # (ncols, M)
        n = np.arange(1, M + 1.0); col = np.abs(A).sum(axis=0)
        return [float((n ** j * col).sum()) for j in range(5)]
    marks = list(kmarks); t1 = time.time(); nstep = 0; out = []; fm = S.FM(M)
    while marks:
        th = float(fm.thmax(u[None])[0]); k = 2 * np.pi / L * th
        if Vt is None and t >= t_mid:
            Uu, Ss, VVt = np.linalg.svd(V.T, full_matrices=False)              # V.T = (M, ncols)
            r = int(np.sum(Ss >= srel * Ss[0])); V = (Uu[:, :r] * Ss[:r]).T.copy(); Vt = VVt[:r]
            print(f'  t = {t:.4g}: SVD of the propagated set, kept {r} directions (S_1 = {Ss[0]:.2e}, S_2/S_1 = {Ss[1] / Ss[0]:.1e}); '
                  f'|e|_0..4 = ' + ' '.join('%.2e' % x for x in hull(V, Vt, M)), flush=True)
        if k >= marks[0]:
            e = hull(V, Vt, M); q = 2 * np.pi / L
            res = dict(k=k, t=t, L=L, M=M, e=e, C2=q * e[1] + q * th * e[0], C5=q ** 4 * e[4])
            print(f'  k = {k:.2f} (t = {t:.7f}, M = {M}, {nstep} steps, {time.time() - t1:.0f}s): |e|_0..4 = ' + ' '.join('%.2e' % x for x in e)
                  + f';  C2 = {res["C2"]:.2e}, C5 = {res["C5"]:.2e}', flush=True)
            out.append(res); marks.pop(0); continue
        if np.abs(u[M - 40:]).max() > 1e-30 and M < Mmax:
            u = np.concatenate([u, np.zeros(80)]); V = np.hstack([V, np.zeros((V.shape[0], 80))]); M += 80; fm = S.FM(M)
        hs = cdt / th ** 4
        E, E2, Qc, f1, f2, f3 = etd_coeffs(-fm.lam, hs)
        def F(x, Y):
            fl = fm.fields(x[None]); return fm.N(x[None])[0], fm.Japply(fl, Y[None])[0]
        for _ in range(10):
            Nu, NV = F(u, V)
            a = E2 * u + Qc * Nu; aV = E2 * V + Qc * NV; Na, NaV = F(a, aV)
            bq = E2 * u + Qc * Na; bV = E2 * V + Qc * NaV; Nb, NbV = F(bq, bV)
            c = E2 * a + Qc * (2 * Nb - Nu); cV = E2 * aV + Qc * (2 * NbV - NV); Nc, NcV = F(c, cV)
            t += (L / 2 / np.pi) ** 4 * hs; L *= np.exp(-0.5 * np.sum(fm.n ** 4 * u ** 2) * hs)
            u = E * u + f1 * Nu + 2 * f2 * (Na + Nb) + f3 * Nc
            V = E * V + f1 * NV + 2 * f2 * (NaV + NbV) + f3 * NcV
            nstep += 1
        if nstep % 500 == 0:
            print(f'    ... step {nstep}: t = {t:.6g}, k = {k:.3f}, M = {M} ({time.time() - t1:.0f}s); |e|_1 = {hull(V, Vt, M)[1]:.2e}', flush=True)
    json.dump(out, open(f'{ckdir}/budget_from_slab{slab}.json', 'w'), indent=1)
    return out


def run_initial(nslab=2, M=1000, K=24, rho_g=0.5, eps_reg=1e-30, c=2.0):
    """initial layer from the certified Gamma0 (gamma0_final.json): graded slabs dt_k = rho_g * t_float,k"""
    g = json.load(open(f'{T2DIR}/gamma0_final.json'))
    bfull = np.array([float.fromhex(x) for x in g['b_hex']]); j = int(g['j'])
    bj = S.arb(int(g['b_j_center_exact']['man'])) * S.arb(2) ** int(g['b_j_center_exact']['exp'])
    hw = float(g['b_j_halfwidth']); L0 = float.fromhex(g['L0_hex']); x0 = float.fromhex(g['X0_hex'][0]); tf0 = float(g['t'])
    xa = [S.arb(float(v)) for v in bfull[:M]]; xa[j - 1] = bj
    b = bfull[:M].copy(); b[j - 1] = float(bj.mid())
    fm = S.FM(M)
    hrule = lambda x, tp, Lf: min(rho_g * (tp + tf0) / (Lf / 2 / np.pi) ** 4, c / fm.thmax(x[None])[0] ** 4)
    h = hrule(b, 0.0, L0); Ma = int(np.sum(fm.lam * h <= 8.0))
    st = S.initial_state(b, Ma, 4 * M, eps_reg, fm)
    st.w[j - 1] += hw * (1 + 1e-12); st.box[Ma:M] += np.abs(st.G[Ma:, j - 1]) * hw * (1 + 1e-12)
    m2 = min(len(bfull), 4 * M); st.box[M:m2] += np.abs(bfull[M:m2]) * (1 + 1e-12)
    st.xbar_arb = xa; st.L = S.arb(L0); st.t = S.arb(0); st.X0 = S.arb(x0)
    print(f'Gamma0_final ({g["source"]}): b_{j} +/- {hw:.2e}; M={M} (dropped {len(bfull) - M} modes -> box), first h={h:.3e} (dt={h * (L0 / 2 / np.pi) ** 4:.2e}), Ma={Ma}, K={K}', flush=True)
    # first slab: slaved start (ODE at s=0) for lam h > 26 -- the slaved Gamma0 agrees with it to ~1e-15 relative; the
    # difference enters the initial box.  Later slabs start from the (consistent) base end values.
    return run_chain(st, fm, nslab, hrule, 'initial', K=K, lamcut1=26.0, etacap1=1e5)


_FMC = {}
def thmax_any(x):
    M = len(x)
    if M not in _FMC: _FMC[M] = S.FM(M)
    return _FMC[M].thmax(x[None])[0]


def mrule_decay(x, tol=1e-30, step=80, Mmax=1200):
    """grow M by `step` (never shrink during a run) while the top 40 modes of the centre exceed tol"""
    M = len(x)
    return int(min(Mmax, M + step)) if np.abs(x[M - 40:]).max() > tol else M


def run_swing(nslab=16, M=480, c=2.0, thstop=40.0):
    b, L, t0, cos_drop, rest = load_sine(f'{SNAP}/p2_t4.0.npz', M)
    x0 = float(np.load(f'{SNAP}/p2_t4.0.npz')['X0'][0])
    fm = S.FM(M); h0 = c / fm.thmax(b[None])[0] ** 4
    t1 = time.time(); xa, xprev = prep_consistent(fm, b, h0)
    b = np.array([float(x.mid()) for x in xa]); h = c / fm.thmax(b[None])[0] ** 4; Ma = int(np.sum(fm.lam * h <= 8.0))
    S0 = S.slaving_propagated(xprev, fm, h0 / 4, Ma)
    st = S.initial_state(b, Ma, 4 * M, 1e-20, fm, S0); st.xbar_arb = xa
    dtr = (h0 / 4 + h0 / 64 + 2 * h0 / 4) * (L / 2 / np.pi) ** 4
    st.L = S.arb(L); st.t = S.arb(t0 + dtr); st.X0 = S.arb(x0)       # test data (taken as exact)
    print(f'swing: start t~{t0 + dtr:.4f} theta={fm.thmax(b[None])[0]:.2f} Ma={Ma} M={M} (prep {time.time() - t1:.0f}s)', flush=True)
    hrule = lambda x, tp, Lf: c / thmax_any(x) ** 4
    return run_chain(st, fm, nslab, hrule, 'swing', stop=lambda info: info['theta_end'] >= thstop, Mrule=mrule_decay)


def run_approach(nslab=3, M=1000, c=2.0, kstop=12.0):
    tag = 'approach' if M == 1000 else f'approach_M{M}'
    sys.path.insert(0, '/home/claude/s1spec/step4/feas2'); sys.path.insert(0, '/home/claude/s1spec/step4')
    from run2 import load_state
    H = load_state(f'{SNAP}/p2_t4.59.npz'); H.tail_tol = 1e-15
    while H.N < 4096: H.refine()
    t1 = time.time(); H.run_to(kstop=kstop)
    N = H.N; th = H.th_hat; bb = -2 * th.imag / N
    print(f'hls regenerated: t={H.t:.9f} kmax={H.kmax():.4f} N={N} L={H.L:.6f} ({time.time() - t1:.0f}s); |b_n| at n=min(M,N/2): {abs(bb[min(M, len(bb) - 1)]):.1e}, max beyond: {np.abs(bb[M + 1:]).max() if len(bb) > M + 1 else 0.0:.1e}', flush=True)
    b = np.zeros(M); m = min(M, len(bb) - 1); b[:m] = bb[1:m + 1]; L = H.L; x0 = float(H.X0[0]); t0 = H.t
    fm = S.FM(M); h0 = c / fm.thmax(b[None])[0] ** 4
    t1 = time.time(); xa, xprev = prep_consistent(fm, b, h0)
    b = np.array([float(x.mid()) for x in xa]); h = c / fm.thmax(b[None])[0] ** 4; Ma = int(np.sum(fm.lam * h <= 8.0))
    S0 = S.slaving_propagated(xprev, fm, h0 / 4, Ma)
    st = S.initial_state(b, Ma, 4 * M, 1e-20, fm, S0); st.xbar_arb = xa
    dtr = (h0 / 4 + h0 / 64 + 2 * h0 / 4) * (L / 2 / np.pi) ** 4
    st.L = S.arb(L); st.t = S.arb(t0 + dtr); st.X0 = S.arb(x0)
    print(f'approach: start t~{t0 + dtr:.6f} theta={fm.thmax(b[None])[0]:.2f} Ma={Ma} M={M} (prep {time.time() - t1:.0f}s)', flush=True)
    hrule = lambda x, tp, Lf: c / fm.thmax(x[None])[0] ** 4
    # K = 24: near T the time scale is T - t; the decay test r^K >= 1e24 then allows dt ~ 0.33 (T-t) (K=16: 0.12 (T-t))
    return run_chain(st, fm, nslab, hrule, tag, K=24, kmax_last=2048)


def run_smoke():
    """fast check of block add/drop + quadratures (t=2, M=200)"""
    M = 300; b, L, t0, _, _ = load_sine(f'{SNAP}/p2_t2.0.npz', M); x0 = float(np.load(f'{SNAP}/p2_t2.0.npz')['X0'][0])
    fm = S.FM(M); h0 = 2.0 / fm.thmax(b[None])[0] ** 4
    xa, xprev = prep_consistent(fm, b, h0)
    b = np.array([float(x.mid()) for x in xa]); Ma = int(np.sum(fm.lam * h0 <= 8.0))
    st = S.initial_state(b, Ma, 4 * M, 1e-20, fm, S.slaving_propagated(xprev, fm, h0 / 4, Ma)); st.xbar_arb = xa
    st.L = S.arb(L); st.t = S.arb(t0); st.X0 = S.arb(x0)
    seq = [Ma + 4, Ma - 2, Ma - 2]
    for k, an in enumerate(seq):
        st, info = S.validate_slab(st, fm, h0, verbose=False, plan_next=lambda uh, an=an: (h0, an), kmax_ng=512 if k == 2 else 0)
        print(summary(info, L)); print(qline(info, st), flush=True)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'one'
    c = float(sys.argv[2]) if len(sys.argv) > 2 and mode not in ('initial', 'swing', 'approach', 'budget', 'slave', 'estimate', 'smoke') else 2.0
    if mode == 'one': run_one(c, 1)
    elif mode == 'chain': run_one(c, int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    elif mode == 'gamma0': run_gamma0(c, int(sys.argv[3]) if len(sys.argv) > 3 else 1000)
    elif mode == 'gamma0c': run_gamma0(c, int(sys.argv[3]) if len(sys.argv) > 3 else 1000, consistent=True)
    elif mode == 'initial': run_initial(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    elif mode == 'swing': run_swing(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    elif mode == 'approach': run_approach(int(sys.argv[2]) if len(sys.argv) > 2 else 3, int(sys.argv[3]) if len(sys.argv) > 3 else 1000)
    elif mode == 'smoke': run_smoke()
    elif mode == 'slave': run_slave()
    elif mode == 'estimate': run_estimate(kend=float(sys.argv[2]) if len(sys.argv) > 2 else 12.0)
    elif mode == 'budget': run_budget(sys.argv[2], int(sys.argv[3]), tuple(float(x) for x in sys.argv[4].split(',')) if len(sys.argv) > 4 else (12.0, 20.0))
