"""v3.3 review M9, item 3: Lemma O (H4), int (h0^(5))^2 dX <= M^2, for EVERY member of the far family that v3.3 uses (L = 113, X = 2,
M = 3/20 from lemmaO_v32_consts.json far_v32; centres at hand-off radius >= 490, Z5 / zone3_interface.json), on the arm and, through
the far end, on the doubled curve (the mirror arm is the reflection: same |kappa^(m)|, |e|).  v3.2 (p2/t8se/h4loc.py) certified only the
two centres 376 and 602.  Per-cylinder bound: h4loc.py's (imported read only: err_c4, chart_h5, h5_form, cells_of; same formula)
    sqrt(int (h0^(5))^2 dX) <= L^3.5 sqrt(M_hi) [A5 sqrt(D_2(tau*) / (28 y_e^2)) + sqrt(|S|) R + sqrt(|S|) sup|h_base^(5)|],
applied to ROWS of centres: a row is an interval [a, b] of base arclength of the centre, its support S = [a - X L - s, b + X L + s]
(s = SLACK = 3: h4loc's excess term X L (1/M_lo - 1) + (1 - lam_lo) x_hi + |x_c| is asserted <= s; h4loc used 1, which fails beyond
x ~ 2e5 because of the frame scale), cells = slab-180 T6 jet cells meeting S by arclength or by radius (as h4loc's far cylinders);
if S passes the far end (x_hi > ell_lo) it continues on the mirror arm: cells with x >= 2 ell_lo - x_hi are added.  The sups over the
row's cells hold for every centre of the row, so one evaluation certifies the whole row.
Changes against h4loc.py, all conservative:
  - y_e := min(y0, y_l/2), y0 = (lowest arclength of S on the arm) - X0 - 1, y_l = ell_lo - X0 (constants_b1.json, slab 180): the
    level-2 weight is T_2 = 28 Y_l^2 with the far-end flattening Y_l = y_l Phi(y/y_l) >= min(y, y_l/2) (LEMMA_F 1.3), so
    int_S w_4^2 <= D_2 / (28 y_e^2); the C4 pointwise bounds e2, e3, e4 are taken at y_e (they decrease in the weight; h4loc: at y0);
  - slack s = 3 instead of 1 (asserted as above).
First centre: hand-off radius rho = |Gamma~ + e N - x_c| / lam >= 490 implies |Gamma~| >= 490 lam_lo - |x_c| - |e|; XSTART := the least
x_lo over the cells whose radius hull reaches that bound (cell-wise |e| as h4loc's E0); every far centre lies at base arclength >= XSTART.
Rows: from floor(XSTART) with widths max(113, a/5), the last row ending at or beyond the largest cell x_hi (the far end; centres beyond
the far end are mirror images).  Also, as fidelity checks with h4loc's exact cylinder geometry (s = 1, y_e = y0): the v3.2 centres
376 and 602, compared with the pinned h4loc_v32.json totals.
Thresholds: total <= M^2 = 0.0225 per row; margin = M^2 / total (rounded down).
usage (cwd p2/t8se, as h4loc): python3 -I -B ../v33/h4far_v33.py constants/amp_A1_w180.json   (env T6DIR=../t6c, T9JSON=../out/t9_cone_kch.json)
writes p2/out_v33/h4far_v33.json"""
import os, sys, json, math, hashlib
from fractions import Fraction
HERE = '/home/claude/s1spec/step4/p2/t8se'
sys.path.insert(0, HERE)
import numpy as np
import h4loc as H
from flint import arb, ctx
ctx.prec = 128
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U, Lw, X0, TC = H.U, H.Lw, H.X0, H.TC
P2 = '/home/claude/s1spec/step4/p2'
LOC = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json'; CB1 = f'{P2}/t8se/constants/constants_b1.json'
H4V32 = f'{P2}/t8se/constants/h4loc_v32.json'; Z3 = '/home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json'
far = json.load(open(LOC))['far_v32']; case = far['case']
assert far['L'] == 113 and case['X'] == '2'
L, X = arb(113), 2; MF = Fraction(case['M']); M = arb(MF.numerator) / MF.denominator
z3 = json.load(open(Z3))['cone_family']['covering']; assert 'centres >= 490' in z3['far_members']
RHO0 = 490; SLACK = 3.0


def main():
    amp_path = sys.argv[1]; A = json.load(open(amp_path)); h = A['handoff']; assert h['window_slabs'] == [180, 180]
    D2 = arb(repr(h['at_Ybar']['D2_goodtime'])); F = [arb(repr(v)) ** 2 for v in h['at_Ybar']['sqrtF']]
    h5 = H.h5_form()
    T9 = json.load(open(H.T9J)); fr = T9['frame']; bx = T9['b5']['boxes']
    LAM = H.iv(*fr['lam']); XCF = arb(max(abs(v) for v in fr['xc']))
    assert LAM.upper() <= 1
    YBOX = arb(bx['Y']); B0, B1, B2 = arb(repr(bx['b0'])), arb(repr(bx['b1'])), arb(repr(bx['b2'])); CE = arb(repr(bx['c_e'])); SF0 = arb(repr(bx['sqrtF0_at_Y']))
    ell = [p for p in json.load(open(CB1))['per_slab'] if p['slab'] == 180][0]['scalars']['ell']
    ELL_LO = arb(repr(ell[0])); YL = ELL_LO - X0
    cells = H.cells_of([180])
    xl = np.array([c['x'][0] for c in cells]); xh = np.array([c['x'][1] for c in cells])
    ELLHI = float(xh.max()); assert ELLHI >= float(ELL_LO.lower())
    Ecell = [CE * YBOX + SF0 * ((arb(c['x'][1]) - X0) ** 3 / 84).sqrt() for c in cells]
    # first centre: radius >= 490 in hand-off units
    reach = [i for i, c in enumerate(cells) if U(c['r']) >= Lw(arb(RHO0) * LAM.lower() - XCF - Ecell[i])]
    XSTART = min(float(xl[i]) for i in reach)
    assert all(float(xl[i]) >= XSTART for i in reach)

    def cyl(name, xlo, xhi, s, flat=True):
        """h4loc's far-cylinder bound on the support [xlo, xhi] (base arclength on the doubled curve), slack s"""
        mir = xhi + 1 > float(ELL_LO.lower())
        xm = min(xlo, float((2 * ELL_LO - arb(xhi)).lower())) if mir else xlo           # lowest arm arclength of S
        sel = [c for c in cells if (float(c['r'].upper()) >= xlo - 1 and float(c['r'].lower()) <= xhi + 1) or (c['x'][1] >= xm - 1 and c['x'][0] <= xhi + 1)]
        ks = [max((H.absb(c['k'][m]) for c in sel), key=U) for m in range(5)]
        E0 = max((CE * YBOX + SF0 * ((arb(c['x'][1]) - X0) ** 3 / 84).sqrt() for c in sel), key=U)
        y0 = arb(xm) - X0 - 1
        assert y0 >= 28
        ye = min(y0, YL / 2, key=U) if flat else y0
        assert ye.lower() > 28
        e2, e3, e4 = H.err_c4(F, D2, ye)
        Bx = [E0 / y0, B0, min(B1, y0 * e2, key=U), min(B2, y0 ** 2 * e3, key=U), y0 ** 3 * e4]
        beta = B0 / (1 - ks[0] * E0)
        P = ks[0] * (X * L + 2) * arb('1.01') + beta
        q, Mq = H.chart_h5(y0, arb(U(P)), [U(v) for v in ks], Bx[:4])
        A5 = arb(U(abs(q.l.get(5, arb(0))) + q.r.get(5, arb(0))))
        R = arb(U(sum(((abs(q.l.get(i, arb(0))) + q.r.get(i, arb(0))) * Bx[i] for i in range(5)), arb(0)) / y0 ** 4))
        Mt = TC.enclose(Mq); mlo, mhi = arb(Lw(Mt)), arb(U(Mt)); assert mlo > 0
        ext = X * L * (1 / mlo - 1) + (1 - LAM.lower()) * arb(xhi) + XCF
        assert ext <= s and X * L / mlo <= (X * L + 2) * arb('1.01'), (name, U(ext))
        b5 = h5(ks[0], ks[1], ks[2], ks[3], arb(U(P)))
        S_ = arb(xhi) - arb(xlo)
        base = L ** 7 * mhi * S_ * b5 * b5
        err = L ** 7 * mhi * A5 * A5 * D2 / (28 * ye * ye)
        errR = L ** 7 * mhi * S_ * R * R
        tot = (err.sqrt() + errR.sqrt() + base.sqrt()) ** 2
        return dict(name=name, support=[xlo, xhi], mirror=mir, cells=len(sel), y0=Lw(y0), y_e=Lw(ye), E0=U(E0), k=[U(v) for v in ks], P=U(P), ext=U(ext),
                    chart_factor=U(A5 * mhi.sqrt()), error=U(err), error_lower=U(errR), base=U(base), total=U(tot), M2=U(M * M),
                    ratio=U(tot / (M * M)), margin=Lw((M * M) / tot), ok=bool(tot < M * M))

    # fidelity: h4loc's two v3.2 far cylinders (s = 1, y_e = y0, its support formula)
    ref = {r['name']: r for r in json.load(open(H4V32))['windows'][amp_path]['cylinders']}
    fid = []
    for c in (376, 602):
        xlo = float((arb(c) - X * L).lower()) - 1; xhi = float((arb(c) + X * L).upper()) + 1
        r = cyl(f'far centre {c}', xlo, xhi, 1.0, flat=False); p = ref[f'far centre {c}']
        fid.append(dict(name=r['name'], total=r['total'], pinned_total=p['total'], identical=r['total'] == p['total'], rel=r['total'] / p['total'] - 1))
    assert all(abs(f['rel']) < 1e-12 for f in fid), fid
    # rows covering every far centre
    rows = []; a = float(math.floor(XSTART))
    while True:
        b = a + max(113.0, a / 5)
        if b >= ELLHI: b = math.ceil(ELLHI)
        rows.append(cyl(f'far centres x in [{a:.0f}, {b:.0f}]', a - 226 - SLACK, b + 226 + SLACK, SLACK) | dict(centres=[a, b]))
        if b >= ELLHI: break
        a = b
    assert rows[0]['centres'][0] <= XSTART and rows[-1]['centres'][1] >= ELLHI
    assert all(rows[i]['centres'][1] == rows[i + 1]['centres'][0] for i in range(len(rows) - 1))
    worst = max(rows, key=lambda r: r['ratio'])
    out = dict(item='Lemma O (H4) for the whole far family of v3.3 (L = 113, X = 2, M = 3/20, centres at hand-off radius >= 490, to the far end)',
               label='certified', M=str(MF), M2=U(M * M), rho0=RHO0, XSTART=XSTART, ell_lo=Lw(ELL_LO), ell_cells_max=ELLHI, y_l_half=Lw(YL / 2), slack=SLACK,
               rows=rows, n_rows=len(rows), all_ok=all(r['ok'] for r in rows), min_margin=min(r['margin'] for r in rows), worst=worst['name'],
               margin_first_row=rows[0]['margin'], fidelity_vs_h4loc_v32=fid,
               args=dict(argv=sys.argv[1:], env={k: os.environ[k] for k in ('T6DIR', 'T9JSON') if k in os.environ}, jets=H.JD, t9json=H.T9J,
                         inputs={p: sha(p) for p in (os.path.abspath(__file__), f'{HERE}/h4loc.py', f'{P2}/t9chart.py', os.path.abspath(amp_path), H.T9J,
                                                     os.path.join(H.JD, 'jets_180.npz'), LOC, CB1, H4V32, Z3)}))
    fo = '/home/claude/s1spec/step4/p2/out_v33/h4far_v33.json'
    json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
    for f in fid: print(f"fidelity {f['name']}: total {f['total']:.6e} vs pinned {f['pinned_total']:.6e} (identical {f['identical']})")
    print(f"XSTART = {XSTART:.3f} (radius {RHO0}); ell_lo = {float(ELL_LO.mid()):.1f}; {len(rows)} rows")
    for r in rows: print(f"  {r['name']:>36s}: cells {r['cells']:3d} y_e {r['y_e']:9.1f}  err {r['error']:.3e} lower {r['error_lower']:.1e} base {r['base']:.2e}  total/M^2 {r['ratio']:.4g}  x{r['margin']:.2f}")
    print(f"min margin x{out['min_margin']:.3f} ({out['worst']})")
    assert out['all_ok'], 'far H4 fails'


if __name__ == '__main__':
    main()
