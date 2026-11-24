"""v3.3 task 6, (H2f) stage 2 (one t_c ball): base part of
    G_j = int_{ln 60}^{ln 700} (28/y) (D^j (f - U_1))^2 ds,  D = y d/dy, s = ln y,  U_1 = 1.3056 alpha y^1.0246,   j = 0..3,
on the t_c ball m of slab 180, with step 4's window machinery p2/t9h/win.py (imported read only; T9TC, T9IFT as in win_cmds.txt, but
T9IFT = the widened ball p2/out_v33/h2f/ift_m.json of h2f_ift.py).  The geometry is built to x = 730 (win.py: 400) and the cells cover
x in [45, 725] (win.py: [45, 390]).  `cells_far` below is win.energies(g, 'win') with the x-range as a parameter; asserted identical to
win.energies on win's own range [45, 390] in this run (same cells, same floats), so the per-cell bound is win.py's.
Per cell (base arclength piece): the cell's y-hull [ylo, yhi] contains the image of the piece; the integrand bound is
U(28 z_j^2 / y^2) on the cell, so int over the piece <= U(28 z^2/y^2) * (yhi - ylo) (width rounded up here; win.py subtracts floats).
Coverage (asserted): cos(phi) > 0 on every cell (y strictly increasing in x, so the cell images tile [y(45), y(725)]), some cell lies
below y = 60 and some above y = 700.  G_j <= the arb sum (rounded up) over the cells meeting [60, 700].
Informational: sup|f|, sup|y f_y|, sup|z_j| on y in [379.5, 700]; the parts on [60, 379.5] and [379.5, 700].
usage: python3 -I h2f_win.py M       (writes p2/out_v33/h2f/win_M.json)"""
import sys, os, json, math, time, hashlib
m = int(sys.argv[1])
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33/h2f'
tc = -1 + (2 * m + 1) / 32
os.environ['T9TC'] = f'{tc},0.03125'; os.environ['T9IFT'] = f'{OV}/ift_{m}.json'      # read by ift.py / win.py at import
sys.path.insert(0, f'{P2}/t9h')
import win
from flint import arb
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U_ = win.U_
XMAX = 730.0; XR = (45.0, XMAX - 5.0); YA, YB, YMID = 60.0, 700.0, 379.5


def cells_far(g, xlo, xhi):
    """win.energies(g, 'win') with x-range (xlo, xhi); returns win's dict per cell plus diagnostics"""
    out = []
    for i, sc, sr, h in win.cells(g, xlo, xhi):
        qc = g.at(i, arb(sc)); qb = g.at(i, arb(sc, sr))
        yc, rc, dc = win.quant(qc); yb, rb, db = win.quant(qb)
        hx = arb(0, sr * h / 2 * 1.0000001)
        Uc = arb('1.3056') * win.AL * yc ** arb('1.0246'); Ub = arb('1.3056') * win.AL * yb ** arb('1.0246')
        zc = [dc[j] - arb('1.0246') ** j * Uc for j in range(5)]
        zb = [db[j] - arb('1.0246') ** j * Ub for j in range(5)]
        rp = 1 + win.LAM * yb * (qb['th'] - win.TH).tan() / (qb['th'] - win.TH).cos() * qb['k0']
        zcell = []
        for j in range(4):
            if j + 2 <= 4:
                z2 = (zb[j + 2] - zb[j + 1] * rp) / (rb * rb)
                zcell.append(zc[j] + zc[j + 1] / rc * hx + arb(0, U_(z2)) * hx * hx / 2)
            else:
                zcell.append(zc[j] + arb(0, U_(zb[j + 1] / rb)) * hx)
        yr = yb
        c = [U_(arb(28) * zcell[j] * zcell[j] / (yr * yr)) * (float(yr.upper()) - float(yr.lower())) for j in range(4)]   # win.py's value
        wd = math.nextafter(float((arb(float(yr.upper())) - arb(float(yr.lower()))).upper()), math.inf)                       # width, rounded up
        cout = [math.nextafter(float((arb(U_(arb(28) * zcell[j] * zcell[j] / (yr * yr))) * arb(wd)).upper()), math.inf) for j in range(4)]
        cphi = (qb['th'] - win.TH).cos()
        out.append(dict(ylo=float(yr.lower()), yhi=float(yr.upper()), c=c, cout=cout, cos_lo=float(cphi.lower()),
                        f=U_(db[0]), yfy=U_(db[1]), z=[U_(zcell[j]) for j in range(4)]))
    return out


if __name__ == '__main__':
    t0 = time.time()
    g = win.basegeo.Geo(180, XMAX, win.ift.TC)
    W = cells_far(g, *XR)
    # fidelity: identical to win.energies on win's own range
    ref = win.energies(g, 'win'); mine = cells_far(g, 45.0, 390.0)
    same = len(ref) == len(mine) and all(a['ylo'] == b['ylo'] and a['yhi'] == b['yhi'] and a['c'] == b['c'] for a, b in zip(ref, mine))
    assert same, 'cells_far differs from win.energies on [45, 390]'
    assert all(c['cos_lo'] > 0 for c in W), 'cos(phi) <= 0 on a cell'
    assert min(c['yhi'] for c in W) < YA and max(c['ylo'] for c in W) > YB, 'cells do not cover [60, 700]'
    def tot(lo, hi):
        s = [arb(0)] * 4
        for c in W:
            if c['yhi'] >= lo and c['ylo'] <= hi:
                for j in range(4): s[j] += arb(c['cout'][j])
        return [math.nextafter(float(v.upper()), math.inf) for v in s]
    far = [c for c in W if c['yhi'] >= YMID and c['ylo'] <= YB]
    out = dict(ball=m, tc_ball=[tc - 1 / 32, tc + 1 / 32], G_60_700=tot(YA, YB), G_60_379p5=tot(YA, YMID), G_379p5_700=tot(YMID, YB),
               ncells=len(W), y_range=[min(c['ylo'] for c in W), max(c['yhi'] for c in W)], min_cos_phi=min(c['cos_lo'] for c in W),
               fidelity_vs_win_energies_45_390=dict(identical=same, cells=len(ref)),
               info_379p5_700=dict(sup_f=max(c['f'] for c in far), sup_yfy=max(c['yfy'] for c in far), sup_z=[max(c['z'][j] for c in far) for j in range(4)]),
               alpha=[float(win.AL.lower()), float(win.AL.upper())], cpu=round(time.time() - t0, 1),
               args=dict(argv=sys.argv[1:], env={k: os.environ[k] for k in ('T9TC', 'T9IFT') if k in os.environ}, xmax=XMAX, x_range=list(XR),
                         y_range=[YA, YB], inputs={p: sha(p) for p in (os.path.abspath(__file__), os.environ['T9IFT'], f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py',
                                                                       f'{P2}/t9h/s1chart.py', '/home/claude/s1spec/step4/t1/polished_v2/slab_180.pkl')}))
    tmp = f'{OV}/win_{m}.json.tmp'
    json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, f'{OV}/win_{m}.json')
    print(f"ball {m}: G_60_700 = {['%.4e' % v for v in out['G_60_700'][:3]]}  ({out['cpu']}s)")
