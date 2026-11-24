"""(a),(b) Plan B: certified level-0 coercivity mu1 of L1 in W, and level-2 Gaarding quantities, on [10.25, 260].
 mu1 (lint_cert_p.py formula) ; mu2 for the PLAIN W and for the level-2 weight m2*W (m2 = 1 + (x/60)^4; psi -> psi - log m2 in the P_j) ; e1, e0.
 Adaptive bisection (dependency control only); every final piece is a rigorous enclosure; NaN asserts.  multiprocessing Pool(2).
 Output planB_lint.json: per-region minima, profile in x (min over 2.5-unit cells) and the global minima."""
import sys, json, time
from multiprocessing import Pool
sys.path.insert(0, '.')
from pb_common import *
from pb_common import _Pr, _lam
P1, P2, P3, P4, P5, k, m, p, q = sp.symbols('P1 P2 P3 P4 P5 k m p q')
rules = {P1: P2, P2: P3, P3: P4, P4: P5, p: -k*q, q: 1 + k*p, k: m, m: p/4}
def D(e): return sp.expand(sum(sp.diff(e, v)*r for v, r in rules.items()))
def rr2(f): return D(D(f)) - 2*P1*D(f) + (P1**2 - P2)*f
def forms(c2, c1, c0):
    A = 2*(P1**2 - P2) - c2
    V = c0 - (D(c1) - P1*c1)/2 + rr2(c2)/2 - rr2(P1**2 - P2)/2
    return sp.expand(-A**2/2 - rr2(A)/2 - V)
a2, a1, a0 = -k**2, -(7*k*m + q/4), -8*m**2 - sp.Rational(5, 2)*k*p
c1t = 2*D(a2) + a1; c0t = D(D(a2)) + 2*D(a1) + a0
e1 = sp.expand(D(D(a1)) + 2*D(a0)); e0 = sp.expand(D(D(a0)))
mu1 = forms(a2, a1, a0); mu2 = forms(a2, c1t, c0t)
assert not (mu1.has(P5) or mu2.has(P5))
mk = lambda e: eval("lambda P1, P2, P3, P4, k, m, p, q: " + _Pr().doprint(e), {"arb": arb})
MU1, MU2, E1, E0 = mk(mu1), mk(mu2), mk(e1), mk(e0)
WTOL = 0.01; DEPTH = 14
def piece(xa, xb):
    X = I(xa, xb); sig_lo = max(float((arb(xa) + STIP).lower()), 8.0)
    kb, mb, pb, qs = HR.HS.profile_bounds(arb(sig_lo))
    K_, M_, Pp = arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper())); Q_ = X + STIP + arb(0, float(qs.upper()))
    Pa = Pj_box(X, 0); Pb = Pj_box(X, 1)
    return dict(mu1=MU1(*Pa, K_, M_, Pp, Q_), mu2=MU2(*Pa, K_, M_, Pp, Q_), mu2w=MU2(*Pb, K_, M_, Pp, Q_), mu1w=MU1(*Pb, K_, M_, Pp, Q_),
                e1=E1(*Pa, K_, M_, Pp, Q_), e0=E0(*Pa, K_, M_, Pp, Q_), e1w=E1(*Pb, K_, M_, Pp, Q_), e0w=E0(*Pb, K_, M_, Pp, Q_))
def cert(a, b, out, d=0):
    r = piece(a, b)
    wid = max(float(r[n].upper()) - float(r[n].lower()) for n in ('mu1', 'mu2w'))
    if wid > WTOL and d < DEPTH:
        mid = 0.5 * (a + b); cert(a, mid, out, d + 1); cert(mid, b, out, d + 1); return
    for n, v in r.items(): chk(v)
    out.append((a, b, float(r['mu1'].lower()), float(r['mu2'].lower()), float(r['mu2w'].lower()), float(r['mu1w'].lower()),
                up(abs(r['e1'])), up(abs(r['e0'])), up(abs(r['e1w'])), up(abs(r['e0w']))))
def work(chunk):
    out = []
    for a, b in chunk: cert(a, b, out)
    return out
if __name__ == "__main__":
    t0 = time.time(); step_lo, step_hi = 0.01, 0.05
    grid = [10.25 + step_lo * i for i in range(int(round(49.75 / step_lo)) + 1)]
    grid += [60 + step_hi * i for i in range(1, int(round(200 / step_hi)) + 1)]
    pairs = list(zip(grid[:-1], grid[1:])); chunks = [pairs[i:i + 200] for i in range(0, len(pairs), 200)]
    with Pool(2) as pool: res = [r for part in pool.map(work, chunks) for r in part]
    res.sort(); print(f"{len(res)} final pieces [{time.time()-t0:.0f}s]")
    names = ['mu1', 'mu2', 'mu2w', 'mu1w']; out = {"n_pieces": len(res), "runtime_s": time.time() - t0, "regions": {}, "cell_profile": {}}
    for i, (ra, rb) in enumerate(REG[:-1] if False else REG):
        sel = [r for r in res if r[1] > ra + 1e-9 and r[0] < rb - 1e-9]
        out["regions"][reg_name(i)] = {n: min(r[2 + j] for r in sel) for j, n in enumerate(names)}
        out["regions"][reg_name(i)].update(e1=max(r[6] for r in sel), e0=max(r[7] for r in sel), e1w=max(r[8] for r in sel), e0w=max(r[9] for r in sel))
        out["regions"][reg_name(i)]["mu1_at"] = min(sel, key=lambda r: r[2])[0]
    for j, n in enumerate(names): g = min(res, key=lambda r: r[2 + j]); out[f"global_min_{n}"] = [g[2 + j], g[0]]
    out["global_max_e1"] = max(r[6] for r in res); out["global_max_e0"] = max(r[7] for r in res)
    out["global_max_e1w"] = max(r[8] for r in res); out["global_max_e0w"] = max(r[9] for r in res)
    cells = {}
    for r in res:
        c = int((r[0] - 10.25) // 2.5); cells.setdefault(c, []).append(r)
    out["cell_profile"] = {f"{10.25 + 2.5 * c:.2f}": {"mu1": min(r[2] for r in v), "mu2": min(r[3] for r in v), "mu2w": min(r[4] for r in v)} for c, v in sorted(cells.items())}
    json.dump(out, open('planB_lint.json', 'w'), indent=1)
    print(json.dumps({k_: v for k_, v in out.items() if k_ != 'cell_profile'}, indent=1))
