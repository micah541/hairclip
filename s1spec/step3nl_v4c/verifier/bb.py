"""Rigorous upper bound of  rate(cell) = sup_{shapes, alpha in [a_lo,a_hi], Ly in [L_lo,L_hi]} tot/Ly   (tot = H + L - NEG, see rhs_iv)
by interval branch-and-bound over the shape box (w,t,tI,u,q) in [0,1]^4 x [0,q*]:
  A = s cw ct, B = s cw st/sqrt(eps), V0 = s sw cI/sqrt(mu), V2 = s sw sI/sqrt(mu epsI), D4 = s q cu/sqrt(eps), V4 = s q su/sqrt(mu epsI)
  with (cx,sx) = ((1-t^2)/(1+t^2), 2t/(1+t^2))  (rational angle parametrisation, no libm).
SCALE (Ly) HANDLING: (H-NEG)/Ly is nondecreasing in s=sqrt(Ly) and L/Ly nonincreasing, so over a Ly-cell the sup is bounded by
  H/s_hi^2 - NEG/s_hi^2 + L/s_lo^2  -- no interval-dependency loss from Ly.
TAIL (unbounded D4,V4): see qcut / find_qstar (radial lemma: POS(r)/r^2 nonincreasing, POS(r*) <= r*^2/2  =>  tot(r) <= tot(r*) for r>=r*).
BOUNDING: (i) natural interval extension  P.hi - N.lo ;  (ii) lower-corner mean-value (monotone) form
  tot(x) <= tot(x_lo) + sum_k max_{0<=D<=w_k} [(G_k - 2 n_k x_lo,k) D - n_k D^2],   G_k >= dP/dx_k on the rectangle (forward-mode AD on nonneg
  intervals; sqrt(a b) replaced by the AM-GM majorant (th a + b/th)/2 in AD mode so that all gradients are finite).
Leaves are accepted when ub - lb <= max(tol, rtol |lb|); result U = max(leaf ub, lb) where lb is a certified point value."""
import numpy as np, time
from ival import *
import rhs_iv as R

def _s_pts(Ly_lo, Ly_hi):
    """outward point intervals for s_lo=sqrt(Ly_lo), s_hi=sqrt(Ly_hi)"""
    a = fd(np.sqrt(Ly_lo)); b = fu(np.sqrt(Ly_hi))
    return P(np.float64(a), np.float64(a)), P(np.float64(b), np.float64(b))

def plain_rate(al, sl, sh, lo, hi, floors=True):
    return R.rate_ub_lo(al, sl, sh, lo, hi, floors)

def pos_norm_hi(al, sl, sh, lo, hi, floors=True):
    xh = R.xhat_ranges(lo, hi)
    (Hh, _, _), (_, Ll, _) = R.two_scales(al, sl, sh, xh, floors)
    return (Hh + Ll).hi

def find_qstar(al, sl, sh, floors=True):
    q = 1.0
    for _ in range(60):
        lo = np.array([[0, 0, 0, 0, q]], float); hi = np.array([[1, 1, 1, 1, q]], float)
        if pos_norm_hi(al, sl, sh, lo, hi, floors)[0] <= 0.5 * q * q * 0.999 * R.radial_factor(): return q
        q *= 2
    raise RuntimeError("no q*")

def qcut(al, sl, sh, lo, hi, floors=True, maxhalv=40):
    """per-box tail truncation (radial lemma, valid for every s in the cell): shrink q1 while POS_norm(q)<=q^2/2 holds on the box."""
    lo = lo.copy(); hi = hi.copy()
    q0 = lo[:, 4]; q1 = hi[:, 4].copy(); act = np.ones(len(lo), bool)
    for _ in range(maxhalv):
        cand = np.maximum(0.5 * q1, q0)
        idx = np.nonzero(act & (cand < q1))[0]
        if len(idx) == 0: break
        l2 = lo[idx].copy(); h2 = hi[idx].copy(); l2[:, 4] = cand[idx]; h2[:, 4] = cand[idx]
        ok = pos_norm_hi(al, sl, sh, l2, h2, floors) <= 0.5 * cand[idx] ** 2 * 0.999 * R.radial_factor()
        q1[idx[ok]] = cand[idx[ok]]
        act[idx[~ok]] = False
        act[idx[ok & (cand[idx] <= q0[idx])]] = False
    hi[:, 4] = q1
    return lo, hi

def mvf_rate_ub(al, sl, sh, lo, hi, floors=True):
    """lower-corner mean-value upper bound of the rate over shape boxes; returns (ub, phi(6,N))."""
    n = len(lo)
    xh = R.xhat_ranges(lo, hi)
    xlo = np.array([x.lo for x in xh]); xhi = np.array([x.hi for x in xh])
    seeds = np.eye(6)
    xg = [PG(x.lo, x.hi, np.repeat(seeds[k][:, None], n, axis=1)) for k, x in enumerate(xh)]
    (Hg, _, _), (_, Lg, _) = R.two_scales(al, sl, sh, xg, floors)
    G = Hg.g + Lg.g
    pl = [P(xlo[k], xlo[k]) for k in range(6)]
    (Hc, _, Nc), (_, Lc, _) = R.two_scales(al, sl, sh, pl, floors)
    nv = R.nvec()[:, None]; w = xhi - xlo            # nv: LOWER bounds of the n_k (relative margin 1e-14 >> float error); w exact differences of floats -> rounded up below
    w = fu(xhi - xlo)
    with np.errstate(invalid="ignore", over="ignore"):
        t1 = fd(fd(2 * nv) * xlo)                    # lower bound of 2 n_k x_lo,k
        g = fu(G - t1)                               # upper bound of G_k - 2 n_k x_lo,k (G is an upper bound of the gradient)
        g = np.where(np.isfinite(G), g, np.inf)
        big = g >= fu(fu(2 * nv * w) * (1 + 1e-12))  # then the maximum of g D - n D^2 on [0,w] is at D = w (g/(2n) >= w); otherwise use g^2/(4n) >= max (valid always)
        quad = fu(fu(g * g) / fd(4 * nv))
        lin = fu(fu(g * w) - fd(fd(nv * w) * w))
        phi = np.where(g <= 0, 0.0, np.where(big, lin, quad))
        phi = np.where(w > 0, phi, 0.0)
    Nlo = fd(fd((nv * xlo * xlo).sum(axis=0)) * (1 - 1e-14))        # lower bound of sum n_k x_lo,k^2  (3 products + 5 additions: error << 1e-14 relative)
    corner = fu((Hc + Lc).hi - Nlo)
    return fu(corner + fu(phi.sum(axis=0) * (1 + 1e-14))), phi

def _rel(f):
    return (f.hi - f.lo) / np.maximum(f.hi, 1e-300)

def _split_phi(lo, hi, span, phi):
    """split along the beta-dim that most reduces the relative width of the physical variable with the largest mean-value penalty phi_k"""
    n = len(lo); k = phi.argmax(axis=0)
    cw, sw = R.circ(lo[:, 0], hi[:, 0]); ct, st = R.circ(lo[:, 1], hi[:, 1]); cI, sI = R.circ(lo[:, 2], hi[:, 2]); cu, su = R.circ(lo[:, 3], hi[:, 3])
    q = P(lo[:, 4], hi[:, 4])
    # candidates per physical variable: (dimA, relwidthA), (dimB, relwidthB)
    cand = {0: ((0, _rel(cw)), (1, _rel(ct))), 1: ((0, _rel(cw)), (1, _rel(st))), 2: ((3, _rel(cu)), (4, _rel(q))),
            3: ((0, _rel(sw)), (2, _rel(cI))), 4: ((0, _rel(sw)), (2, _rel(sI))), 5: ((3, _rel(su)), (4, _rel(q)))}
    d = np.zeros(n, int)
    for kk in range(6):
        m = k == kk
        if m.any():
            (da, ra), (db, rb) = cand[kk]
            d[m] = np.where(ra[m] >= rb[m], da, db)
    wd = (hi - lo) / span
    bad = (hi - lo)[np.arange(n), d] <= 1e-12 * span[d]
    d[bad] = wd[bad].argmax(axis=1)
    idx = np.arange(n); l = lo[idx, d]; h = hi[idx, d]
    m = np.where(l <= 0, h / 8.0, np.where(h > 4 * l, np.sqrt(np.maximum(l, 1e-300) * h), 0.5 * (l + h)))
    m = np.where((l <= 0) & (h > 0.30 * span[d]), 0.5 * h, m)
    lo2 = lo.copy(); hi1 = hi.copy(); hi1[idx, d] = m; lo2[idx, d] = m
    return np.vstack([lo, lo2]), np.vstack([hi1, hi])

def _cem(al, sl, sh, qs, floors, seeds, iters=25, n=4000):
    # lower bounds are certified point values: use the point alpha = upper end of the alpha interval
    rng = np.random.default_rng(12345)
    lg_lo = np.array([-6, -6, -6, -6, -4.0]); lg_hi = np.array([0, 0, 0, 0, np.log10(qs)])
    mu = 0.5 * (lg_lo + lg_hi); sd = 0.5 * (lg_hi - lg_lo); best = None; lb = -np.inf; nev = 0
    for p in (seeds or []):
        p = np.asarray(p, float)[None, :]
        r = max(plain_rate(al, sl, sl, p, p, floors)[0][0], plain_rate(al, sh, sh, p, p, floors)[0][0]); nev += 2
        if r > lb: lb = r; best = p[0].copy()
    for it in range(iters):
        y = np.clip(mu + sd * rng.standard_normal((n, 5)), lg_lo, lg_hi)
        if best is not None: y[:200] = np.clip(np.log10(np.maximum(best, 1e-12)) + 0.1 * sd * rng.standard_normal((200, 5)), lg_lo, lg_hi)
        p = 10 ** y; p[:, :4] = np.minimum(p[:, :4], 1.0); p[:, 4] = np.minimum(p[:, 4], qs)
        rl = np.maximum(plain_rate(al, sl, sl, p, p, floors)[0], plain_rate(al, sh, sh, p, p, floors)[0]); nev += 2 * n
        o = np.argsort(-rl)[:150]
        if rl[o[0]] > lb: lb = rl[o[0]]; best = p[o[0]].copy()
        mu = y[o].mean(axis=0); sd = np.maximum(y[o].std(axis=0), 0.02)
    return lb, best, nev

LEAVES = None
def bound_cell(a_lo, a_hi, L_lo, L_hi, tol=3e-3, rtol=0.03, floors=True, minw=1e-7, chunk=15000, maxeval=6_000_000, seeds=None, verbose=False, cem_iters=25):
    """rigorous bound for alpha in [a_lo,a_hi], Ly in [L_lo,L_hi].  returns dict(U=..., lb=..., ...)"""
    al = P(np.float64(a_lo), np.float64(a_hi)) if a_hi > 0 else P(np.float64(0), np.float64(0))
    if a_hi > 0: al = P(np.float64(fd(a_lo)), np.float64(fu(a_hi)))
    sl, sh = _s_pts(L_lo, L_hi); pointcell = (L_lo == L_hi)
    qs = find_qstar(al, sl, sh, floors)
    span = np.array([1, 1, 1, 1, qs])
    root = (np.zeros((1, 5)), np.concatenate([np.ones((1, 4)), [[qs]]], axis=1))
    alp = P(al.hi, al.hi)
    lb, best, nev = _cem(alp, sl, sh, qs, floors, seeds, iters=cem_iters)
    if verbose: print("phase1 lb", lb, best, flush=True)
    stack = [root]; U = -np.inf; nleaf = 0; nminw = 0
    while stack:
        lo, hi = stack.pop()
        if len(lo) > chunk:
            h = len(lo) // 2; stack.append((lo[h:], hi[h:])); lo, hi = lo[:h], hi[:h]
        lo, hi = qcut(al, sl, sh, lo, hi, floors)
        n = len(lo); mid = 0.5 * (lo + hi)
        (_, ub), = [plain_rate(al, sl, sh, lo, hi, floors)]
        mm = np.vstack([mid, mid])
        rl_s = plain_rate(alp, sh, sh, mid, mid, floors)[0]
        if not pointcell: rl_s = np.maximum(rl_s, plain_rate(alp, sl, sl, mid, mid, floors)[0])
        nev += 3 * n
        assert np.all(np.isfinite(ub)) and np.all(np.isfinite(rl_s))
        k = rl_s.argmax()
        if rl_s[k] > lb: lb = rl_s[k]; best = mid[k].copy()
        keep = ub > lb
        if not keep.any(): continue
        lo, hi, ub = lo[keep], hi[keep], ub[keep]
        m, phi = mvf_rate_ub(al, sl, sh, lo, hi, floors); nev += 2 * len(lo)
        ub2 = np.minimum(ub, m)
        keep = ub2 > lb
        if not keep.any(): continue
        lo, hi, ub2, phi = lo[keep], hi[keep], ub2[keep], phi[:, keep]
        done_tol = (ub2 - lb) <= max(tol, rtol * abs(lb))
        done_w = (((hi - lo) / span).max(axis=1) <= minw)
        done = done_tol | done_w
        nminw += int((done_w & ~done_tol).sum())
        if done.any():
            U = max(U, ub2[done].max()); nleaf += int(done.sum())
            if LEAVES is not None and len(LEAVES) < 200: LEAVES.append((lo[done].copy(), hi[done].copy(), ub2[done].copy()))
        lo, hi, phi = lo[~done], hi[~done], phi[:, ~done]
        if len(lo) == 0: continue
        if verbose and (nev // 1000) % 200 == 0: print('nev', nev, 'stack', sum(len(a) for a, b in stack), 'cur', len(lo), 'ubmax', ub2[~done].max(), 'leaf', nleaf, 'U', U, flush=True)
        stack.append(_split_phi(lo, hi, span, phi))
        if nev > maxeval: raise RuntimeError("maxeval exceeded; lb=%g U_so_far=%g" % (lb, U))
    U = max(U, lb)
    return dict(U=float(U), lb=float(lb), qstar=qs, nev=nev, nleaf=nleaf, nminw=nminw, best=[float(v) for v in best])

if __name__ == "__main__":
    import sys
    a = float(sys.argv[1]); L = float(sys.argv[2]); w = float(sys.argv[3]); tol = float(sys.argv[4])
    t0 = time.time()
    r = bound_cell(a, a, L, L * (1 + w), tol=tol, verbose=False, maxeval=int(float(sys.argv[5])) if len(sys.argv) > 5 else 6_000_000)
    print(r, round(time.time() - t0, 1), "s")
