#!/usr/bin/env python3
"""P1 (Lemma U, phase 1): CERTIFIED (python-flint arb, outward) suprema of the weight constants of the Lemma-U weight
    Theta = W(sigma) D(s),  W = 1/(e^{-mu(sigma+a)} + 1 + e^{nu(sigma - ln1.5 - a)}),  D = 1 + g e^{-4s},  sigma = s - s_rho,  s = ln y,
over the geometry  s >= ln 28,  s_rho >= ln yL  (crude yL = 56, fine yL = 150),  all a in [a_lo, a_hi]:
    Gamma_m = sup|Theta^(m)/Theta| (m = 1,2,3),   kappa~ = sup|L|,  L = (ln eps Theta)' = -4 + Theta'/Theta,   kappa2~ = sup(L^2 + L')_+ .

Method (no sampling):
  * p = A/S, q = B/S (S = A + 1 + B) satisfy p' = -p(mu(1-p) + nu q) < 0, q' = q(nu(1-q) + mu p) > 0, and W^(k)/W = P_k(p,q) with
    P_1 = mu p - nu q, P_{k+1} = P_k' + P_k P_1 (polynomials; sympy, exact integer coefficients).
  * D^(j)/D = (-4)^j theta (j >= 1), theta = g e^{-4s}/D in [0, theta_max(sigma)],  theta_max(sigma) = g/(max(28, yL e^sigma)^4 + g)
    (every (sigma, theta) with 0 < theta <= theta_max(sigma) is realised by the geometry, and no other).
  * Theta^(m)/Theta = sum_k C(m,k) P_k (-4)^{m-k} theta [k<m] + P_m ;  L = -4 + P_1 - 4 theta ;
    L^2 + L' = (P_1-4)^2 + P_1' + theta (16 - 8(P_1-4))   (the theta^2 terms cancel exactly) : all AFFINE in theta, so theta is
    carried as an independent arb ball [0, theta_max(cell)] without dependency loss.
  * sigma in [SL, SR]: adaptive arb cells (sigma ball x a ball), p, q from the stable forms
        p = 1/(1 + e^{mu(sigma+a)} + e^{(mu+nu)sigma + (mu-nu)a - nu c}),  q = 1/(1 + e^{-(mu+nu)sigma - (mu-nu)a + nu c} + e^{-nu(sigma-c-a)}),  c = ln 1.5.
    Tails sigma <= SL and sigma >= SR (this covers sigma -> -inf, sigma -> +inf, s -> inf): by monotonicity p in [p(SL),1], q in [0,q(SL)],
    theta in [0, theta_a];  resp. p in [0,p(SR)], q in [q(SR),1], theta in [0, theta_max(SR)]  -> one arb evaluation each.
  The certified bound is the max of the arb upper endpoints over all leaves and both tails.  Cells whose bound exceeds the exact
  sigma -> +-inf limit (times 1 + 1e-12) are bisected (to depth DMAX); a leaf still above it would mean the float claim is wrong (reported).
Import-safe (no side effects); CLI writes weight_consts.json (args block)."""
import os, sys, json, math, hashlib, time, platform
sys.dont_write_bytecode = True
import sympy as sp
from flint import arb, ctx
ctx.prec = 120
Q_ = sp.Rational

def A(q):
    q = Q_(q); return arb(int(q.p))/arb(int(q.q))

_P = {}
def polys(mu, nu):
    """P_k (k = 1..3), P_1' as python callables of arb (p, q); exact integer coefficients"""
    key = (mu, nu)
    if key in _P: return _P[key]
    p, q = sp.symbols('p q'); mu_, nu_ = sp.Integer(mu), sp.Integer(nu)
    dp = -p*(mu_*(1 - p) + nu_*q); dq = q*(nu_*(1 - q) + mu_*p)
    d = lambda e: sp.expand(sp.diff(e, p)*dp + sp.diff(e, q)*dq)
    P = {0: sp.Integer(1), 1: sp.expand(mu_*p - nu_*q)}
    for k in (1, 2): P[k + 1] = sp.expand(d(P[k]) + P[k]*P[1])
    P1p = d(P[1])
    # sanity: exact limits P_k(1,0) = mu^k, P_k(0,1) = (-nu)^k ; and (ln W)'' = P1' <= 0 form
    for k in (1, 2, 3):
        assert P[k].subs({p: 1, q: 0}) == mu_**k and P[k].subs({p: 0, q: 1}) == (-nu_)**k
    f = {k: sp.lambdify((p, q), P[k], modules=[{}]) for k in (1, 2, 3)}
    f['P1p'] = sp.lambdify((p, q), P1p, modules=[{}])
    _P[key] = (f, {k: str(v) for k, v in P.items()}, str(P1p))
    return _P[key]

def quantities(f, p, q, th):
    """arb enclosures of (Gamma_1, Gamma_2, Gamma_3 integrands, L, L^2+L') on the box (p, q, theta)"""
    w = {0: arb(1), 1: f[1](p, q), 2: f[2](p, q), 3: f[3](p, q)}
    G = []
    for m in (1, 2, 3):
        cth = arb(0)
        for k in range(m): cth += math.comb(m, k)*w[k]*(arb(-4)**(m - k))
        G.append(w[m] + th*cth)            # theta multiplies ONE coefficient: affine, no dependency loss
    c0 = w[1] - 4
    L = c0 - 4*th
    K2 = c0*c0 + f['P1p'](p, q) + th*(16 - 8*c0)
    return G, L, K2

def pq(mu, nu, sig, a, c):
    e1 = (mu*(sig + a)).exp(); e2 = ((mu + nu)*sig + (mu - nu)*a - nu*c).exp()
    e3 = (-(mu + nu)*sig - (mu - nu)*a + nu*c).exp(); e4 = (-nu*(sig - c - a)).exp()
    return 1/(1 + e1 + e2), 1/(1 + e3 + e4)

def theta_max_up(g, yL, sig_lo):
    """certified upper bound of g/(max(28, yL e^sig)^4 + g) for sigma >= sig_lo (decreasing in sigma)"""
    y = arb(yL)*sig_lo.exp()
    ylow = max(arb(28), arb(y.lower()))
    return (g/(ylow**4 + g)).upper()

def certify(mu, nu, a_lo, a_hi, theta_aq, yL, SL=-12, SR=12, h0=Q_(1, 20), na=6, DMAX=14, tol=1e-12):
    """returns dict of certified upper bounds (arb exact upper endpoints) and diagnostics"""
    f, Pstr, P1pstr = polys(mu, nu)
    theta_aq = Q_(theta_aq); gq = theta_aq*28**4/(1 - theta_aq); g = A(gq); c = arb(Q_(3, 2).p)/arb(Q_(3, 2).q); c = c.log()
    lim = dict(G1=arb(max(mu, nu)), G2=arb(max(mu, nu)**2), G3=arb(max(mu, nu)**3), kap=arb(nu + 4), kap2=arb((nu + 4)**2))
    lim['kap'] = arb(max(nu + 4, abs(mu - 4)))
    names = ['G1', 'G2', 'G3', 'kap', 'kap2']
    best = {k: arb(0) for k in names}; where = {k: None for k in names}; nleaf = [0]; overs = []
    def vals(p, q, thmax):
        # all five quantities are AFFINE in theta (see doc): sup over theta in [0, thmax] of |.| (or (.)_+) is attained at an endpoint,
        # so evaluate at theta = 0 and theta = thmax (an exact arb upper endpoint) and take the max.
        out = {k: arb(0) for k in names}
        for th in (arb(0), arb(thmax)):
            G, L, K2 = quantities(f, p, q, th)
            for k, x in (('G1', abs(G[0]).upper()), ('G2', abs(G[1]).upper()), ('G3', abs(G[2]).upper()), ('kap', abs(L).upper()), ('kap2', max(arb(0), K2.upper()))):
                if x > out[k]: out[k] = x
        return out
    def upd(v, tag):
        for k in names:
            if v[k] > best[k]: best[k] = v[k]; where[k] = tag
    # abox = [(A(a_lo) + (A(a_hi) - A(a_lo))*i/na).union(A(a_lo) + (A(a_hi) - A(a_lo))*(i + 1)/na) for i in range(na)]
    aall = A(a_lo).union(A(a_hi))
    thA = (g/(28**4 + g)).upper()
    # inner tail sigma <= SL
    pL, qL = pq(mu, nu, arb(SL), aall, c)
    v = vals(arb(pL.lower()).union(arb(1)), arb(0).union(arb(qL.upper())), thA); upd(v, ('tail-', SL)); tails = {'inner': {k: float(x) for k, x in v.items()}}
    # outer tail sigma >= SR
    pR, qR = pq(mu, nu, arb(SR), aall, c)
    thR = theta_max_up(g, yL, arb(SR))
    v = vals(arb(0).union(arb(pR.upper())), arb(qR.lower()).union(arb(1)), thR); upd(v, ('tail+', SR)); tails['outer'] = {k: float(x) for k, x in v.items()}
    # middle: adaptive cells
    N = int((Q_(SR) - Q_(SL))/h0); a_lo, a_hi = Q_(a_lo), Q_(a_hi)
    da = (a_hi - a_lo)/na
    stack = [(Q_(SL) + h0*(N - 1 - i), Q_(SL) + h0*(N - i), a_lo + da*ia, a_lo + da*(ia + 1), 0) for i in range(N) for ia in range(na)]
    while stack:                               # depth-first: memory O(N*na + DMAX)
        lo, hi, al, ah, dep = stack.pop()
        sig = A(lo).union(A(hi)); ab = A(al).union(A(ah)); p, q = pq(mu, nu, sig, ab, c)
        v = vals(p, q, theta_max_up(g, yL, A(lo)))
        over = [k for k in names if v[k] > (lim[k]*(1 + arb(tol))).upper()]
        if over and dep < DMAX:
            if ah > al and (ah - al) > 4*(hi - lo):          # split a when its cell is relatively wide
                m_ = (al + ah)/2; stack.append((lo, hi, al, m_, dep + 1)); stack.append((lo, hi, m_, ah, dep + 1))
            else:
                m_ = (lo + hi)/2; stack.append((lo, m_, al, ah, dep + 1)); stack.append((m_, hi, al, ah, dep + 1))
            continue
        if over: overs.append((float(lo), float(hi), float(al), float(ah), over, {k: float(v[k]) for k in over}))
        nleaf[0] += 1; upd(v, ('cell', float(lo), float(hi), float(al), float(ah)))
    return dict(mu=mu, nu=nu, a=[str(Q_(a_lo)), str(Q_(a_hi))], yL=yL, theta_a=str(theta_aq), g=str(gq), g_float=float(gq),
                bound={k: best[k] for k in names}, bound_float_up={k: _fup(best[k]) for k in names}, limits={k: float(lim[k]) for k in names},
                argmax={k: where[k] for k in names}, leaves=nleaf[0], unresolved_cells=overs, tails=tails, P=Pstr, P1p=P1pstr,
                geometry=dict(SL=SL, SR=SR, h0=str(h0), na=na, DMAX=DMAX, tol=tol))

def _fup(x):
    v = float(x)
    return v if arb(v) >= x else math.nextafter(v, math.inf)

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

if __name__ == '__main__':
    t0 = time.time(); HERE = os.path.dirname(os.path.abspath(__file__)); out = dict(cases=[])
    # theta_a values at which the constants are evaluated: any theta_a in (0, 0.6] (bound reported per value); checked at 0.6 (worst) and the fixed points
    for (mu, nu, a_lo, a_hi) in ((4, 4, 0, Q_(3, 10)), (4, 2, Q_(1, 2), Q_(1, 2))):
        for yL in (56, 150):
            for th in (Q_(1, 10), Q_(4, 10), Q_(6, 10)):
                r = certify(mu, nu, a_lo, a_hi, th, yL)
                b = r['bound_float_up']
                print(f"mu={mu} nu={nu} a=[{a_lo},{a_hi}] yL={yL} theta_a={th}: Gamma_m <= {b['G1']:.12g}, {b['G2']:.12g}, {b['G3']:.12g}; "
                      f"kappa~ <= {b['kap']:.12g}; kappa2~ <= {b['kap2']:.12g}  (limits {r['limits']}); leaves {r['leaves']}, unresolved {len(r['unresolved_cells'])}", flush=True)
                r['bound'] = {k: str(v) for k, v in r['bound'].items()}
                out['cases'].append(r)
    out['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), flint=__import__('flint').__version__,
                       sympy=sp.__version__, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                       inputs={}, script_sha256=sha(os.path.abspath(__file__)), prec=ctx.prec, wall_s=time.time() - t0)
    p = os.path.join(HERE, 'weight_consts.json')
    json.dump(out, open(p + '.tmp', 'w'), indent=1, default=str); os.replace(p + '.tmp', p)
    print('wrote', p, 'in %.1f s' % (time.time() - t0))
