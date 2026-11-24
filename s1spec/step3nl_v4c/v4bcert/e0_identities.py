"""Machine checks of the algebraic identities used by the v4b E0 modulation lemma (E0_LEMMA_v4b.tex).  [sym] = exact sympy, residual 0;
[num] = random float test (supporting only).  Output out/e0_identities.json.
 I1 [sym] rank-one symmetrisation (S2): the symmetric form (u v^T + v u^T)/2 has eigenvalues (<u,v> +- |u||v|)/2 on span{u,v} (0 elsewhere),
          so sup_{|w|=1} |<w,u><w,v>| = (|u||v| + |<u,v>|)/2.        [num] random R^30.
 I2 [sym] oblique projection (S3): <g,h> = 1 => P = h g^T is idempotent and I - P has the singular value |h||g| (2-dim reduction), all others 1 <= |h||g|.
          [num] random R^30: ||I - h g^T||_2 = |h||g| when <g,h> = 1.
 I3 [sym] modulation algebra: x^T (G + E) = R  <=>  x = H R - H E^T x,  H = (G^-1)^T  (3x3 symbolic).
 I4 [sym] (a,b) w-parts: A1(u) - A1(u0) = -w + q [r(u0+w) - r(u0)],  B1(u) - B1(u0) = cos(th) [r(u0+w) - r(u0)],  r(u) = u'/(1 - k u),
          r(u0+w) - r(u0) = w'/ag(u0+w) + u0' k w/(ag(u0+w) ag(u0)).
 I5 [sym] cubic drift: Dr(u) = -(q/4) u' G(u), G = k u/(1-k u):  Dr(u0+w) - Dr(u0) - Dlin(w) = -(q/4)[w' dG + u0' (dG - G_u(u0) w)],
          dG = G(u0+w) - G(u0) = k w/(ag(u0) ag(u0+w)),  dG - G_u(u0) w = k^2 w^2/(ag(u0)^2 ag(u0+w)).
 I6 [sym] alpha^2 remainders (the first order of A1(u0) is -alpha (psi - q psi')): A1(u0) + alpha (psi - q psi') = alpha^2 (q h2' - h2) + q u0' k u0/ag(u0),  u0'/ag(u0) - alpha psi' = alpha^2 h2' + u0' k u0/ag(u0).
 I7 [sym] (reference) weighted divergence form and drift integration by parts: tc_weighted.py (exact Euler-operator residuals 0).
 I8 [sym] chi_32: Sm(t) = t^4 (35 - 84 t + 70 t^2 - 20 t^3), Sm' = 140 t^3 (1-t)^3, sup|chi'| = 35/48, Sm(0)=0, Sm(1)=1, Sm^(j)(0)=Sm^(j)(1)=0 (j=1..3)."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump
import sympy as sp, numpy as np
res = {}
# I1
a, b, c, s = sp.symbols('a b c s', real=True)
u = sp.Matrix([a, 0]); v = sp.Matrix([b * c, b * s])          # c^2 + s^2 = 1
Mf = (u * v.T + v * u.T) / 2
ev = list(Mf.eigenvals().keys())
target = [(a * b * c + a * b) / 2, (a * b * c - a * b) / 2]      # <u,v> = a b c, |u||v| = a b (a, b >= 0)
ok = all(any(sp.simplify((e - t).subs(s, sp.sqrt(1 - c ** 2))) == 0 for e in ev) for t in target)
res['I1_sym'] = bool(ok)
rng = np.random.default_rng(0); worst = 0.0
for _ in range(200):
    U, V = rng.normal(size=30), rng.normal(size=30)
    S = (np.outer(U, V) + np.outer(V, U)) / 2; lam = np.abs(np.linalg.eigvalsh(S)).max()
    worst = max(worst, abs(lam - (np.linalg.norm(U) * np.linalg.norm(V) + abs(U @ V)) / 2) / lam)
res['I1_num_rel'] = worst
# I2
g1, g2, h1, h2 = sp.symbols('g1 g2 h1 h2', real=True)
g = sp.Matrix([g1, g2]); h = sp.Matrix([h1, h2])
P = h * g.T
cons = {h2: (1 - g1 * h1) / g2}                                    # <g,h> = 1
idem = sp.simplify((P * P - P).subs(cons))
Q = (sp.eye(2) - P).subs(cons); QtQ = sp.simplify(Q.T * Q)
# singular values^2 of I - P: product = det^2 = 0 (I - P is a rank-1 projection in 2D), sum = trace = |h|^2|g|^2
tr = sp.simplify(QtQ.trace() - ((h1 ** 2 + h2 ** 2) * (g1 ** 2 + g2 ** 2)).subs(cons)); dt = sp.simplify(QtQ.det())
res['I2_sym'] = bool(idem == sp.zeros(2, 2) and tr == 0 and dt == 0)
worst = 0.0
for _ in range(200):
    G_, H_ = rng.normal(size=30), rng.normal(size=30); H_ = H_ / (G_ @ H_)
    nrm = np.linalg.norm(np.eye(30) - np.outer(H_, G_), 2); worst = max(worst, abs(nrm - np.linalg.norm(H_) * np.linalg.norm(G_)) / nrm)
res['I2_num_rel'] = worst
# I3
Gs = sp.Matrix(3, 3, sp.symbols('G0:9')); Es = sp.Matrix(3, 3, sp.symbols('E0:9')); Rs = sp.Matrix(3, 1, sp.symbols('R0:3'))
x = ((Gs + Es).T).inv() * Rs                                      # x^T (G+E) = R^T  <=>  (G+E)^T x = R
Hm = (Gs.inv()).T
r3 = sp.simplify(x - (Hm * Rs - Hm * Es.T * x))
res['I3_sym'] = bool(r3 == sp.zeros(3, 1))
# I4, I5, I6
k, q, th, al_ = sp.symbols('k q theta alpha')
U0, U1, W0, W1, y0, y1, H0, H1 = sp.symbols('U0 U1 W0 W1 y0 y1 H0 H1')
ag = lambda u_: 1 - k * u_
r = lambda u_, up_: up_ / ag(u_)
A1 = lambda u_, up_: -(u_ - up_ * q / ag(u_)); B1 = lambda u_, up_: up_ * sp.cos(th) / ag(u_)
dr = r(U0 + W0, U1 + W1) - r(U0, U1)
res['I4_sym'] = bool(sp.simplify(A1(U0 + W0, U1 + W1) - A1(U0, U1) - (-W0 + q * dr)) == 0 and sp.simplify(B1(U0 + W0, U1 + W1) - B1(U0, U1) - sp.cos(th) * dr) == 0
                     and sp.simplify(dr - (W1 / ag(U0 + W0) + U1 * k * W0 / (ag(U0 + W0) * ag(U0)))) == 0)
G = lambda u_: k * u_ / (1 - k * u_); Gu = sp.diff(G(U0), U0)
Dr = lambda u_, up_: -(q / 4) * up_ * G(u_)
Dlin = -(q / 4) * (G(U0) * W1 + U1 * Gu * W0)
dG = G(U0 + W0) - G(U0)
i5a = sp.simplify(Dr(U0 + W0, U1 + W1) - Dr(U0, U1) - Dlin - (-(q / 4) * (W1 * dG + U1 * (dG - Gu * W0))))
i5b = sp.simplify(dG - k * W0 / (ag(U0) * ag(U0 + W0))); i5c = sp.simplify(dG - Gu * W0 - k ** 2 * W0 ** 2 / (ag(U0) ** 2 * ag(U0 + W0)))
res['I5_sym'] = bool(i5a == 0 and i5b == 0 and i5c == 0)
u0 = al_ * y0 + al_ ** 2 * H0; u0p = al_ * y1 + al_ ** 2 * H1
i6a = sp.simplify(A1(u0, u0p) + al_ * (y0 - q * y1) - (al_ ** 2 * (q * H1 - H0) + q * u0p * k * u0 / ag(u0)))
i6b = sp.simplify(u0p / ag(u0) - al_ * y1 - (al_ ** 2 * H1 + u0p * k * u0 / ag(u0)))
res['I6_sym'] = bool(i6a == 0 and i6b == 0)
TW = json.load(open(os.path.join(OUT, 'tc_weighted.json')))
res['I7_ref_tc_weighted'] = dict(X1_reduction_ok=TW['X1_reduction_ok'], euler_curv_ok=TW['euler_curv_ok'], euler_drift_ok=TW['euler_drift_ok'],
                                 drift_split_ok=TW['drift_split_ok'], drift_lin_ok=TW['drift_lin_ok'])
t = sp.symbols('t'); Sm = t ** 4 * (35 - 84 * t + 70 * t ** 2 - 20 * t ** 3)
d1 = sp.diff(Sm, t)
ok8 = sp.expand(d1 - 140 * t ** 3 * (1 - t) ** 3) == 0 and Sm.subs(t, 0) == 0 and Sm.subs(t, 1) == 1
ok8 = ok8 and all(sp.diff(Sm, t, j).subs(t, 0) == 0 and sp.diff(Sm, t, j).subs(t, 1) == 0 for j in (1, 2, 3))
ok8 = ok8 and sp.Rational(140, 64) / 3 == sp.Rational(35, 48) and sp.solve(sp.diff(d1, t), t).count(sp.Rational(1, 2)) == 1
res['I8_sym'] = bool(ok8)
allok = all(v for k_, v in res.items() if k_.endswith('_sym')) and all(res['I7_ref_tc_weighted'].values()) and res['I1_num_rel'] < 1e-10 and res['I2_num_rel'] < 1e-8
res['all_ok'] = bool(allok)
res['args'] = args_block([os.path.join(OUT, 'tc_weighted.json')])
print(json.dumps({k_: v for k_, v in res.items() if k_ != 'args'}, indent=1)); assert allok
print('sha256', dump(os.path.join(OUT, 'e0_identities.json'), res))
