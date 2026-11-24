"""Review A1: Lemma O (H4), int (h0^(5))^2 dX <= M^2 per cylinder (cylinder units), from the good-time level-2 dissipation.

Localized level-2 dissipation on a cylinder support S = [x_lo, x_hi] (y_0 := x_lo - X0 - 1 >= 22, so T_2 = 28 Y_l^2 >= 28 y_0^2
on S and S lies in the core of every level):
    int_S w_4^2 <= D_S / (28 y_0^2),   D_S := int_S T_2 w_4^2 <= D_2,
and at the good time tau* of the window W (LEMMA_F C5, camp handoff):  D_2(tau*) <= [F_2(t_a) + int_W (2 g_2 sqrt F_2 + h_2)
+ inflow] / (theta_2 |W|).  A localized energy E_S = int T_2 psi^2 w_2^2 (psi = 1 on S) has the same leading term: its value at
t_a is bounded only by F_2(t_a), and its inflow (the ramp, or the transport through an inner cut-off) is that of F_2.  So the
window length, not the localization, decides; the localized forcing is <= the global one, which is used here.
In cylinder units (X = x/L, h0 = h/L): int (h0^(5))^2 dX = L^7 int_S (h^(5))^2 dx, and
    error part <= L^7 * 1.001^2 * D_2(tau*) / (28 y_0^2),    base part <= L^7 * |S| * sup_S |h_base^(5)|^2,
with h^(5) of the base graph over the centre tangent (sympy below; kappa^(<=3) from the T6 hulls of the window slabs).
Cylinders as in p2/t9.py and the review: cone r_c in [200, 380] (rows of width 10, L = r_hi/4, support [r_lo/4 - 1, 7 r_hi/4 + 1],
M = 0.12); far centres 376 and 602 (L = 113, half-width 2L, M = 1).
usage: python3 -I h4loc.py AMP_JSON [AMP_JSON ...]     (each a camp output; its handoff window and D2(tau*) are used)"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 't6'))
import sympy as sp
from flint import arb, ctx
import reader as RD
ctx.prec = 128
JD = os.path.join(os.path.normpath(os.path.abspath(os.environ['T6DIR'])) if os.environ.get('T6DIR')
                  else os.path.normpath(os.path.join(HERE, '..', 't6')), 'jets')          # env T6DIR as cconst
X0 = arb('2.2490776515')
U = lambda v: math.nextafter(float(v.upper()), math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))
def absb(a): return arb(U(abs(a)))


def h5_form():
    """h^(5) of a graph h(x) over a line, from k^(<=3) (arclength derivatives) and the angle psi (h' = tan psi); sympy gives
    h^(5) = [-60 k^4 S + 105 k^4 S/C^2 - 86 k^2 k1 C + 105 k^2 k1/C + 15 k k2 S + 10 k1^2 S + k3 C]/C^7  (S = sin psi, C = cos psi),
    checked below. Returned as a termwise majorant in K_m = sup|k^(m)| >= 0 and Pm = sup|psi| < pi/2."""
    s = sp.symbols('s'); k = sp.Function('k')(s); psi = sp.Function('psi')(s)
    D = lambda f: sp.diff(f, s).subs(sp.Derivative(psi, s), k) / sp.cos(psi)
    h = sp.tan(psi)
    for _ in range(4): h = D(h)
    K = sp.symbols('K0:4'); P = sp.symbols('P')
    rep = {sp.Derivative(k, (s, 3)): K[3], sp.Derivative(k, (s, 2)): K[2], sp.Derivative(k, s): K[1]}
    h = h.subs(rep).subs(k, K[0]).subs(psi, P)
    Sn, Cs = sp.sin(P), sp.cos(P)
    ref = (-60 * K[0] ** 4 * Sn + 105 * K[0] ** 4 * Sn / Cs ** 2 - 86 * K[0] ** 2 * K[1] * Cs + 105 * K[0] ** 2 * K[1] / Cs
           + 15 * K[0] * K[2] * Sn + 10 * K[1] ** 2 * Sn + K[3] * Cs) / Cs ** 7
    assert sp.simplify(h - ref) == 0

    def bound(K0, K1, K2, K3, Pm):
        sn, cs = Pm.sin(), Pm.cos()
        return (60 * K0 ** 4 * sn + 105 * K0 ** 4 * sn / cs ** 2 + 86 * K0 ** 2 * K1 + 105 * K0 ** 2 * K1 / cs + 15 * K0 * K2 * sn
                + 10 * K1 ** 2 * sn + K3) / cs ** 7
    return bound


def cells_of(slabs):
    out = []
    for s in slabs:
        J = RD.jets(s, d=JD)
        for e in range(J['x'].shape[0]):
            for c in range(J['x'].shape[1]):
                g = lambda nm: iv(*J[nm][e, c])
                out.append(dict(k=[g('k0'), g('k1'), g('k2'), g('k3')], p=g('p'), q=g('q'), x=(float(J['x'][e, c, 0]), float(J['x'][e, c, 1]))))
    for c in out: c['r'] = (c['p'] ** 2 + c['q'] ** 2).sqrt()
    return out


def cylinders():
    cyl = [dict(name=f'cone r_c {r}-{r + 10}', rlo=r, rhi=r + 10, L=arb(r + 10) / 4, X=3, M=arb('0.12')) for r in range(200, 380, 10)]
    cyl += [dict(name=f'far centre {c}', rlo=c, rhi=c, L=arb(113), X=2, M=arb(os.environ.get('FAR_M', '1'))) for c in (376, 602)]   # env FAR_M: far M (v3.2: 0.15)
    return cyl


def main():
    amps = [json.load(open(p)) for p in sys.argv[1:]]
    h5 = h5_form()
    rows = {}
    for path, A in zip(sys.argv[1:], amps):
        h = A['handoff']; D2 = arb(repr(h['at_Ybar']['D2_goodtime'])); win = h['window_slabs']; parts = h.get('D2_parts', {})
        cells = cells_of(range(win[0], win[1] + 1))
        res = []
        for cy in cylinders():
            L, X = cy['L'], cy['X']
            if cy['name'].startswith('cone'):                             # union over the row of [r_c/4, 7 r_c/4] (as p2/t9.py)
                xlo, xhi = 0.25 * cy['rlo'] - 1, 1.75 * cy['rhi'] + 1
            else:
                xlo = float((arb(cy['rlo']) - X * L).lower()) - 1; xhi = float((arb(cy['rhi']) + X * L).upper()) + 1
            y0 = arb(xlo) - X0 - 1
            sel = [c for c in cells if float(c['r'].upper()) >= xlo - 1 and float(c['r'].lower()) <= xhi + 1]
            ks = [max((absb(c['k'][m]) for c in sel), key=U) for m in range(4)]
            P = ks[0] * (X * L + 2) * arb('1.01')                         # |psi| <= sup|k| (half-width + 2) * 1.01, as p2/t9.py
            b5 = h5(ks[0], ks[1], ks[2], ks[3], arb(U(P)))
            base = L ** 7 * (arb(xhi) - arb(xlo)) * b5 * b5
            err = L ** 7 * arb('1.001') ** 2 * D2 / (28 * y0 * y0)
            tot = (err.sqrt() + base.sqrt()) ** 2
            res.append(dict(name=cy['name'], L=U(L), y0=float(y0.lower()), error=U(err), base=U(base), total=U(tot), M2=U(cy['M'] ** 2),
                            ratio=U(tot / cy['M'] ** 2)))
        worst = max(res, key=lambda r: r['ratio'])
        rows[path] = dict(window=win, window_length=parts.get('window_length'), D2=U(D2), D2_parts=parts, cylinders=res, worst=worst)
        print(f"{os.path.basename(path)}: window slabs {win}, |W| = {parts.get('window_length', float('nan')):.3e}, D2(tau*) <= {U(D2):.3e} "
              f"(Ybar-free {parts.get('zero', float('nan')):.2e}, linear {parts.get('linear', float('nan')):.2e}, "
              f"quadratic {parts.get('quadratic', float('nan')):.2e})")
        for r in res:
            if r['name'].startswith('far') or r['name'].endswith(('200-210', '290-300', '370-380')):
                print(f"   {r['name']:>18s}: L {r['L']:6.1f}  y0 {r['y0']:6.1f}  error {r['error']:.3e}  base {r['base']:.1e}  "
                      f"total/M^2 = {r['ratio']:.3g}")
        print(f"   worst: {worst['name']} total/M^2 = {worst['ratio']:.3g}")
    env = {k: os.environ[k] for k in ('FAR_M', 'SUBW', 'T9JSON', 'T9WIN', 'T6DIR', 'H4OUT', 'TOPLOSS', 'NLFAC', 'BASE_T1', 'YBAR', 'T6RES', 'SLAB0', 'YOUT') if k in os.environ}
    out = dict(args=dict(argv=sys.argv[1:], amp=sys.argv[1:], jets=JD, chart_factor=1.001, env=env), windows=rows)
    fo = os.environ.get('H4OUT') or os.path.join(HERE, 'constants', 'h4loc.json')          # env H4OUT: output file
    out['args']['out'] = os.path.abspath(fo)
    json.dump(out, open(fo, 'w'), indent=1)


if __name__ == '__main__':
    main()
