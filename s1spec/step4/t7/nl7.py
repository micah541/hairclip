"""T7 (C4): Hessian constants n4, n3, n2 of the normal-graph CDF operator (T5 2.5), per slab.
Graph over the base at the same label (a = a~, b = b~, clock fixed): Gamma = G~ + e nu~, alpha = 1 - k e, J = sqrt(alpha^2 + e'^2),
  kappa = (alpha e'' - e' alpha' + k (alpha^2 + e'^2)) / J^3,  kappa_ss = (kappa_x / J)_x / J,
  e_tau = F(e) := (J/alpha) [-kappa_ss + c (-q e' + (p + e) alpha)/J + b (-t1 e' + n alpha)/J] - V + e' (w - e th_t)/alpha,
(c = 1/4 - a~, w = varpi, th_t = theta~_tau, V = V~).  F(0) = -r, DF(0) = L (t1/signs.py), and exactly
  F = A4(e, e') e'''' + A3(e, e', e'') e''' + A2(e, e', e'')   (checked symbolically below),
so N(e) = F - F(0) - DF(0)e = [A4 - A4(0)] e'''' + [A3 - A3(0)] e''' + [A2 - A2(0) - grad A2(0) z],  z = (e, e', e''), and
  |N| <= n4 (|e| + |e'|)|e''''| + n3 (|e| + |e'| + |e''|)|e'''| + n2 (|e| + |e'| + |e''|)^2  with
  n4 = sup max(|d_e A4|, |d_e' A4|),  n3 = sup max_i |d_zi A3|,  n2 = (1/2) sup max_ij |d_zi d_zj A2|  over the box.
Tip x in [0, 34.25], NO condition on e'' (Gap 9): F is a polynomial in e'' (alpha, J do not contain e''), exactly
  A4 = A4(e, e'),  A3 = A30(e, e') + A31(e, e') e'',  A2 = A20 + A21 e'' + A22 e''^2 + A23 e''^3   (checked symbolically), so
  |N| <= n4 w|e''''| + n3 w|e'''| + m3|e''||e'''| + n2 w^2 + n21 w|e''| + m2 e''^2 + m23|e''|^3,   w = |e| + |e'|,
  n4 = sup|grad A4|, n3 = sup|grad A30|, m3 = sup|A31|, n2 = sup|Hess A20|/2, n21 = sup|grad A21|, m2 = sup|A22|, m23 = sup|A23|
  (grad, Hess in (e, e'); sups over the cells and |e|, |e'| <= 1e-5 only).  T5 form: n3* = max(n3, m3), n2* = max(n2, n21/2, m2),
  |N| <= n4 w|e''''| + n3*(w + |e''|)|e'''| + n2*(w + |e''|)^2 + m23 |e''|^3.
  Far x >= 10.25 (T5 box: x^(j-1)|e^(j)| <= eta_F = 1e-2, |e(10.25)| <= 1e-5):
  |e| <= 1e-5 + 1e-2 (x - 10.25), |e'| <= 1e-2, |e''| <= 1e-2/x; reported per dyadic x-range, with the scale-free products
  n4, x n3, x^2 n2 (T8 Lemma N3 needs a0, a1 from these).
Base jets per cell over the slab (tc-balls): kappa^(0..3) (exact polys), p, q, n, t1 (T6 bivariate), w = w* + A, th_t = Phi + rho_theta,
V = V* + r (T6 identities), c, b.
With --se the far cells (x >= 0.6 ell_lo) use the base jets widened to the (SE)-corrected base (se7.py).
usage: python3 nl7.py --slabs a-b --out DIR [--se]"""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import sympy as sp
from flint import arb, arb_poly, ctx
ctx.prec = 128
E = sp.symbols('e0:5'); K = sp.symbols('k0:4'); p, q, n, t1, w, tht, V, c, b = sp.symbols('p q n t1 w tht V c b')


def derive():
    rules = {p: -K[0] * q, q: 1 + K[0] * p, n: -K[0] * t1, t1: K[0] * n}
    def D(f):
        out = sum(sp.diff(f, E[i]) * E[i + 1] for i in range(4)) + sum(sp.diff(f, K[i]) * K[i + 1] for i in range(3))
        out += sum(sp.diff(f, v) * r for v, r in rules.items())
        return out
    al = 1 - K[0] * E[0]; alp = D(al); J = sp.sqrt(al ** 2 + E[1] ** 2)
    kap = (al * E[2] - E[1] * alp + K[0] * (al ** 2 + E[1] ** 2)) / J ** 3
    ks = D(kap) / J; kss = D(ks) / J
    F = (J / al) * (-kss + c * (-q * E[1] + (p + E[0]) * al) / J + b * (-t1 * E[1] + n * al) / J) - V + E[1] * (w - E[0] * tht) / al
    A4 = sp.diff(F, E[4]); A3 = sp.diff(F.subs(E[4], 0), E[3]); A2 = F.subs({E[4]: 0, E[3]: 0})
    chk = sp.simplify(sp.diff(F, E[4], 2)) == 0 and sp.simplify(sp.diff(F.subs(E[4], 0), E[3], 2)) == 0 and \
        all(sp.diff(A4, v) == 0 for v in (E[2], E[3])) and sp.diff(A3, E[3]) == 0
    z = (E[0], E[1], E[2])
    G4 = [sp.diff(A4, v) for v in z[:2]]; G3 = [sp.diff(A3, v) for v in z]; H2 = [sp.diff(A2, u, v) for i, u in enumerate(z) for v in z[i:]]
    # polynomial structure in e'' (no e'' box needed): A3 affine, A2 cubic in e''
    poly_chk = sp.simplify(sp.diff(A3, E[2], 2)) == 0 and sp.simplify(sp.diff(A2, E[2], 4)) == 0
    A30 = A3.subs(E[2], 0); A31 = sp.diff(A3, E[2])
    A20 = A2.subs(E[2], 0); A21 = sp.diff(A2, E[2]).subs(E[2], 0)
    A22 = (sp.diff(A2, E[2], 2) / 2).subs(E[2], 0); A23 = sp.diff(A2, E[2], 3) / 6
    zz = (E[0], E[1])
    P = dict(G30=[sp.diff(A30, v) for v in zz], A31=[A31], H20=[sp.diff(A20, u, v) for i, u in enumerate(zz) for v in zz[i:]],
             G21=[sp.diff(A21, v) for v in zz], A22=[A22], A23=[A23])
    return dict(A4=A4, A3=A3, A2=A2, G4=G4, G3=G3, H2=H2, P=P, linear_check=bool(chk), poly_check=bool(poly_chk))


class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        bb, ex = e.as_base_exp(); s = f"({self._print(bb)})"
        if ex.is_Integer:
            m = int(ex); prod = "*".join([s] * abs(m)); return f"({prod})" if m > 0 else f"(1/({prod}))"
        if ex == sp.Rational(1, 2): return f"({s}.sqrt())"
        if ex == -sp.Rational(1, 2): return f"(1/{s}.sqrt())"
        return f"({s}**(arb({ex.p})/arb({ex.q})))"


ETA_E_B = 1e-3                     # sharp far box for Lemma N3: |e| <= 1e-3 on x >= 10.25 (|k~ e| <~ 1e-4; T8 closes it by A_j Ybar)
ARGS = list(E[:3]) + list(K) + [p, q, n, t1, w, tht, V, c, b]
_F = None
def funcs():
    global _F
    if _F is None:
        d = derive()
        mk = lambda ex: eval("lambda " + ",".join(str(a) for a in ARGS) + ": " + _Pr().doprint(sp.expand(ex) if False else ex), {"arb": arb})
        _F = dict(G4=[mk(g) for g in d['G4']], G3=[mk(g) for g in d['G3']], H2=[mk(g) for g in d['H2']], check=d['linear_check'],
                  A4=mk(d['A4']), A3=mk(d['A3']), P={k: [mk(g) for g in v] for k, v in d['P'].items()}, poly_check=d['poly_check'])
    return _F


def run_slab(s, nsub=2, verbose=True, se=False, prev=None):
    """prev: the output of an earlier run on the polished base (same nsub).  With se and prev, the far cells with
    34.25 < x and x_hi < x_c (where the corrected base equals the polished one) are not recomputed: their per-range maxima
    are taken from prev (max with the new cells, valid since those cells are unchanged)."""
    import base7, c2r3 as R3
    t0 = time.time(); S = base7.Slab(s); FF = funcs()
    SEW = None; keep = None
    if se:
        import se7; SEW = se7.widths(s)
        if prev is not None: keep = SEW['xc']
    TIPK = ('n4', 'n3_0', 'm3', 'n2_0', 'n21', 'm2', 'm23')
    SB = R3.SlabBP(S, arb(10 ** 7))                         # all elements (T6 chaining)
    tip = {k: 0.0 for k in TIPK}; far = {}
    eta1, eta2, etaF = 1e-5, 1e-2, 1e-2
    for j in range(nsub):
        t = arb((2 * j + 1) / nsub - 1, 1 / nsub)
        keys = ('k0', 'k1', 'k2', 'p', 'q', 'n', 't1', 'r', 'A', 'Q', 'thx', 'Vs')
        els = []
        for el in SB.els:
            d = {k: base7.bp_at(el[k], t) for k in keys}
            d['k3'] = base7.bp_at(el['k2'].ds().tmul(arb(2) / el['h'], SB.B.D), t) if False else None
            els.append(d)
        sc = {k: v(t) for k, v in SB.sc.items()}
        for e, (el, d) in enumerate(zip(SB.els, els)):
            xe, he = SB.B.x[e](t), SB.B.h[e](t)
            xa = float(xe.lower())
            if keep is not None and xa > 34.25 and float((xe + he).upper()) < keep: continue   # unchanged far cells: from --prev
            m = 8 if xa < 40 else 4
            for qq in range(m):
                sb = arb(-1 + (2 * qq + 1) / m, 1 / m)            # s in [-1, 1]
                x = xe + he * (sb + 1) / 2
                k0, k1, k2 = d['k0'](sb), d['k1'](sb), d['k2'](sb)
                k3 = (el['k2'].ds().tmul(SB.B.ih[e] * 2, SB.B.D)) if False else None
                kk3 = base7.bp_at(el['k2'].ds(), t)(sb) * 2 / he
                pp, qv, nn, tt = d['p'](sb), d['q'](sb), d['n'](sb), d['t1'](sb)
                cc, bb = sc['c'], sc['b']
                wv = -k0 * k1 + d['Q'](sb) + cc * (qv - x) + bb * tt + d['A'](sb)
                tht_ = d['thx'](sb); Vv = d['Vs'](sb) + d['r'](sb)
                xl = float(x.lower()); xh = float(x.upper())
                if SEW is not None and xh >= SEW['xc']:              # (SE)-corrected base: widen every jet on the far cells
                    k0, k1, k2, kk3 = (v + arb(0, SEW['k'][m]) for m, v in enumerate((k0, k1, k2, kk3)))
                    pp, qv = pp + arb(0, SEW['p']), qv + arb(0, SEW['q']); nn, tt = nn + arb(0, SEW['n']), tt + arb(0, SEW['t1'])
                    wv = wv + arb(0, SEW['varpi']); tht_ = tht_ + arb(0, SEW['tht']); Vv = Vv + arb(0, SEW['V'])
                regs = []
                if xl <= 34.25:                                      # tip: box on (e, e') only, e'' free (polynomial structure)
                    a = (arb(0, eta1), arb(0, eta1), arb(0), k0, k1, k2, kk3, pp, qv, nn, tt, wv, tht_, Vv, cc, bb)
                    mx = lambda key: max(float(abs(g(*a)).upper()) for g in FF['P'][key])
                    vals = dict(n4=max(float(abs(g(*a)).upper()) for g in FF['G4']), n3_0=mx('G30'), m3=mx('A31'), n2_0=0.5 * mx('H20'),
                                n21=mx('G21'), m2=mx('A22'), m23=mx('A23'))
                    for k_, v_ in vals.items(): tip[k_] = max(tip[k_], v_)
                if xh >= 10.25:
                    regs.append(('far', arb(0, eta1 + etaF * max(xh - 10.25, 0)), arb(0, etaF), arb(0, etaF / max(xl, 10.25))))
                    regs.append(('farB', arb(0, ETA_E_B), arb(0, etaF), arb(0, etaF / max(xl, 10.25))))
                for reg, e0, e1, e2 in regs:
                    a = (e0, e1, e2, k0, k1, k2, kk3, pp, qv, nn, tt, wv, tht_, Vv, cc, bb)
                    n4 = max(float(abs(g(*a)).upper()) for g in FF['G4'])
                    n3 = max(float(abs(g(*a)).upper()) for g in FF['G3'])
                    n2 = 0.5 * max(float(abs(g(*a)).upper()) for g in FF['H2'])
                    key = (reg, int(math.floor(math.log2(max(xl, 10.25) / 10.25))))
                    f = far.setdefault(key, dict(box=reg, x=[10.25 * 2 ** key[1], 10.25 * 2 ** (key[1] + 1)], n4=0.0, n3=0.0, n2=0.0, xn3=0.0, x2n2=0.0))
                    f['n4'] = max(f['n4'], n4); f['n3'] = max(f['n3'], n3); f['n2'] = max(f['n2'], n2)
                    f['xn3'] = max(f['xn3'], n3 * xh); f['x2n2'] = max(f['x2n2'], n2 * xh * xh)
    fmax = lambda b_: dict(n4=max(f['n4'] for f in far.values() if f['box'] == b_), n3=max(f['n3'] for f in far.values() if f['box'] == b_),
                           n2=max(f['n2'] for f in far.values() if f['box'] == b_), xn3=max(f['xn3'] for f in far.values() if f['box'] == b_),
                           x2n2=max(f['x2n2'] for f in far.values() if f['box'] == b_))
    tip['n3'] = max(tip['n3_0'], tip['m3']); tip['n2'] = max(tip['n2_0'], tip['n21'] / 2, tip['m2'])
    if keep is not None:                                            # merge the unchanged far cells of the earlier run
        assert prev['slab'] == s
        for f0 in prev['far']:
            key = (f0['box'], int(round(math.log2(f0['x'][0] / 10.25))))
            f = far.setdefault(key, dict(box=f0['box'], x=list(f0['x']), n4=0.0, n3=0.0, n2=0.0, xn3=0.0, x2n2=0.0))
            for k_ in ('n4', 'n3', 'n2', 'xn3', 'x2n2'): f[k_] = max(f[k_], f0[k_])
    out = dict(slab=s, tau=[S.tau0, S.tau0 + S.dtf], k=S.k0, tip=tip, far=[far[k] for k in sorted(far, key=lambda k_: (k_[0], k_[1]))],
               far_max=fmax('far'), farB_max=fmax('farB'),
               box=dict(tip='|e|,|e1| <= 1e-5 on [0, 34.25], no condition on e2: |N| <= n4 w|e4| + n3 (w+|e2|)|e3| + n2 (w+|e2|)^2 + m23 |e2|^3',
                        far='T5 box: |e| <= 1e-5 + 1e-2 (x - 10.25), |e1| <= 1e-2, |e2| <= 1e-2/x',
                        farB=f'sharp box: |e| <= {ETA_E_B:g} (x >= 10.25; |k e| <~ 1e-4), |e1| <= 1e-2, |e2| <= 1e-2/x'),
               linear_check=FF['check'], poly_check=FF['poly_check'], se=bool(se), far_cells_from_prev=keep is not None,
               t=time.time() - t0)
    if verbose:
        print(f"slab {s} (k {S.k0:.1f}): tip n4 {tip['n4']:.3f} n3 {tip['n3']:.3f} n2 {tip['n2']:.3f} m23 {tip['m23']:.3f} | far (T5 box) n4 {out['far_max']['n4']:.3f} "
              f"x n3 {out['far_max']['xn3']:.3f} x^2 n2 {out['far_max']['x2n2']:.3g} | far (sharp box) n4 {out['farB_max']['n4']:.3g} "
              f"x n3 {out['farB_max']['xn3']:.3g} x^2 n2 {out['farB_max']['x2n2']:.3g} ({out['t']:.0f}s)", flush=True)
    return out



def parse_slabs(spec):
    """'a-b' or comma-separated slabs / ranges, e.g. '5-20' or '5,13,21-23'"""
    out = []
    for part in spec.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; out += list(range(a_, b_ + 1))
        elif part.strip(): out.append(int(part))
    return out

if __name__ == '__main__':
    import base7
    ap = argparse.ArgumentParser(); ap.add_argument('--slabs', default='5-5'); ap.add_argument('--out', default=os.path.join(HERE, 'out_nl'))
    ap.add_argument('--nsub', type=int, default=1); ap.add_argument('--se', action='store_true', help='far cells on the (SE)-corrected base')
    ap.add_argument('--prev', default=None, help='with --se: directory of an earlier run on the polished base (same nsub); '
                    'its far cells below x_c = 0.6 ell_lo are reused')
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    for s in parse_slabs(A.slabs):
        fn = os.path.join(A.out, 'nl_%03d.json' % s)
        if os.path.exists(fn): continue
        pf = os.path.join(A.prev, 'nl_%03d.json' % s) if A.prev else None
        prev = json.load(open(pf)) if pf and os.path.exists(pf) else None
        if prev is not None and prev.get('nsub', 1) != A.nsub: prev = None
        o = run_slab(s, A.nsub, se=A.se, prev=prev); o['nsub'] = A.nsub
        json.dump(o, open(fn + '.tmp', 'w'), indent=1); os.replace(fn + '.tmp', fn)
