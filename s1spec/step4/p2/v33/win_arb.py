"""v3.3 review M9, item 1: arb (outward) replacement of the cell bounds of p2/t9h/win.py (closed v3.2; imported read only for the
unchanged parts: parameter balls AL, XC, LAM, TH from T9IFT, quant() = D^j f by sympy/lambdify, STIR, C4 profiles, s1chart.SIG_TIP).
One t_c ball and one parameter-ball set per process; geometry built to x = 730 (as h2f_win.py), so one run gives every value v3.3 reads:
  (H3) windows  [y, 1.5 y], y in [50, 253]: E_j = int (28/y) (D^j (f - U))^2 ds, j <= 3, sliding cover step 0.5 (as win.sliding);
  (H2) shells   [14, 60]: E_j for w' = f - alpha psi3^' - alpha^2 h2^' (+ the level-0 chart ball 2e-7), j <= 2 (as win 'shell' mode);
  (H2f)         G_j on [60, 700] for f - U_1 (as h2f_win.py), parts [60, 379.5], [379.5, 700].
Changes against win.py / h2f_win.py (everything else is win.py's formula, cell by cell):
  (a) cell boundaries in the element parameter s: float breakpoints b_m = -1 + m (2/n), with b_0 = -1 and b_n = 1 exactly; the cell
      ball is arb(b_m) u arb(b_{m+1}) (contains [b_m, b_{m+1}]; neighbours share the exact endpoint, so the cells tile [-1, 1] with no
      rounding gaps); the centre point s_c is the ball's (exact) midpoint;
  (b) the Taylor half-width |x - x_c| <= |hh_i| rad(s-ball) / 2 with hh_i the arb element length on the t_c ball (win.py: float(mid hh)
      times 1.0000001, not asserted);
  (c) per-cell integral <= upper(28 z_j^2 / y^2 on the cell) * (upper(y) - lower(y)), all in arb; the window, shell and [60, 700] sums
      are arb sums, rounded up once at the end (win.py: float products and float sums);
  (d) the cell y-hull [ylo, yhi] is recorded rounded outward (win.py: float(), nearest), so the selection 'cell meets [A, B]' is safe;
  (e) coverage asserted per range: the cells used are consecutive in curve order, cos(phi) > 0 on every cell (y strictly increasing
      along the arc, so the arc is a graph over y there and the cell images tile the y-range), the first cell lies below the range and
      the last above it: windows [50, 380.25], shells [14, 60], (H2f) [60, 700].
usage: python3 -I -B win_arb.py SET M     SET = v32 (p2/out/h1slab/ift_M.json, v3.2's balls) or wide (p2/out_v33/h2f/ift_M.json, widened)
writes p2/out_v33/m9/win_arb_SET_M.json"""
import sys, os, json, math, time, hashlib
SET, m = sys.argv[1], int(sys.argv[2]); assert SET in ('v32', 'wide') and 0 <= m < 32
P2 = '/home/claude/s1spec/step4/p2'; OM = f'{P2}/out_v33/m9'
IFTP = f'{P2}/out/h1slab/ift_{m}.json' if SET == 'v32' else f'{P2}/out_v33/h2f/ift_{m}.json'
tc = -1 + (2 * m + 1) / 32
os.environ['T9TC'] = f'{tc},0.03125'; os.environ['T9IFT'] = IFTP          # read by ift.py / win.py at import (as win_cmds.txt)
sys.path.insert(0, f'{P2}/t9h')
import win
from flint import arb
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
UP = lambda v: math.nextafter(float(v.upper()), math.inf); LO = lambda v: math.nextafter(float(v.lower()), -math.inf)
XMAX = 730.0; LC = 0.5


def cells(g, xlo, xhi):
    """win.cells with exact shared breakpoints; yields (global index, element i, s-ball, centre s_c, |x - x_c| radius as arb)"""
    k = 0
    for i, el in enumerate(g.els):
        a, h = g.xe_mid[i], float(g.hh[i].mid())                      # floats only to choose n and to skip far-away cells
        n = max(1, math.ceil(h / LC)); ds = 2.0 / n
        bp = [-1.0] + [-1 + j * ds for j in range(1, n)] + [1.0]
        assert all(bp[j] < bp[j + 1] for j in range(n))
        if a + h < xlo or a > xhi: k += n; continue
        for j in range(n):
            kk = k + j
            s0, s1 = bp[j], bp[j + 1]
            if a + h * (s1 + 1) / 2 < xlo or a + h * (s0 + 1) / 2 > xhi: continue
            sb = arb(s0).union(arb(s1)); sc = sb.mid()
            assert sb.lower() <= s0 and sb.upper() >= s1
            hx = arb(0, (abs(g.hh[i]) * sb.rad() / 2).upper())          # sup over the t_c ball of hh |s - s_c| / 2
            yield kk, i, sb, sc, hx
        k += n


def energies(g, mode, xr):
    """win.energies(g, mode) cell by cell with (a)-(d); returns per-cell dicts with arb contributions"""
    out = []
    for kk, i, sb, sc, hx in cells(g, *xr):
        qc = g.at(i, sc); qb = g.at(i, sb)
        yc, rc, dc = win.quant(qc); yb, rb, db = win.quant(qb)
        if mode == 'win':
            Uc = arb('1.3056') * win.AL * yc ** arb('1.0246'); Ub = arb('1.3056') * win.AL * yb ** arb('1.0246')
            zc = [dc[j] - arb('1.0246') ** j * Uc for j in range(5)]
            zb = [db[j] - arb('1.0246') ** j * Ub for j in range(5)]
            J = 4
        else:
            xs = yc - win.STIP; xsb = yb - win.STIP
            pc = win.C4.PSI.jets_any(xs, float(xs.mid()), 5); hc = win.C4.H2.jets_any(xs, float(xs.mid()), 5)
            pb_ = win.C4.PSI.jets_any(xsb, float(xsb.mid()), 5); hb = win.C4.H2.jets_any(xsb, float(xsb.mid()), 5)
            def Wj(j, y, P, H):
                if j == 0: return win.AL * P[1] + win.AL * win.AL * H[1]
                return sum((c * y ** mm * (win.AL * P[mm + 1] + win.AL * win.AL * H[mm + 1]) for mm, c in win.STIR[j].items() if mm + 1 < len(P)), arb(0))
            Bs = arb(0, 2e-7)
            zc = [dc[j] - Wj(j, yc, pc, hc) + (Bs if j == 0 else arb(0)) for j in range(4)]
            zb = [db[j] - Wj(j, yb, pb_, hb) + (Bs if j == 0 else arb(0)) for j in range(4)]
            J = 3
        rp = 1 + win.LAM * yb * (qb['th'] - win.TH).tan() / (qb['th'] - win.TH).cos() * qb['k0']
        zcell = []
        for j in range(J):
            if mode == 'win' and j + 2 <= 4:
                z2 = (zb[j + 2] - zb[j + 1] * rp) / (rb * rb)
                zcell.append(zc[j] + zc[j + 1] / rc * hx + arb(0, abs(z2).upper()) * hx * hx / 2)
            else:
                zcell.append(zc[j] + arb(0, abs(zb[j + 1] / rb).upper()) * hx)
        wy = yb.upper() - yb.lower()                                      # exact
        c = [(arb(28) * zcell[j] * zcell[j] / (yb * yb)).upper() * wy for j in range(J)]
        out.append(dict(k=kk, ylo=LO(yb), yhi=UP(yb), c=c, cos_lo=LO((qb['th'] - win.TH).cos()),
                        hrel=UP(g.hh[i].rad() / abs(g.hh[i]).lower())))
    return out


def covered(cl, A, B, name):
    ks = [c['k'] for c in cl]
    assert ks == list(range(ks[0], ks[0] + len(ks))), (name, 'cells not consecutive')
    assert all(c['cos_lo'] > 0 for c in cl), (name, 'cos(phi) <= 0')
    assert cl[0]['yhi'] < A and cl[-1]['ylo'] > B, (name, 'range not covered', cl[0]['yhi'], cl[-1]['ylo'])
    return dict(cells=len(cl), first_yhi=cl[0]['yhi'], last_ylo=cl[-1]['ylo'], min_cos_phi=min(c['cos_lo'] for c in cl))


def tot(cl, lo, hi, J):
    s = [arb(0)] * J
    for c in cl:
        if c['yhi'] >= lo and c['ylo'] <= hi:
            for j in range(J): s[j] += c['c'][j]
    return s


def sliding(cl, ya0=50.0, ya1=253.0, step=0.5):
    best = [arb(0)] * 4; where = [None] * 4; ya = ya0; n = 0
    while ya <= ya1 + 1e-9:
        t = tot(cl, ya, 1.5 * (ya + step), 4); n += 1
        for j in range(4):
            if t[j].upper() > best[j].upper(): best[j], where[j] = t[j], ya
        ya += step
    return best, where, n


if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(OM, exist_ok=True)
    g = win.basegeo.Geo(180, XMAX, win.ift.TC)
    W = energies(g, 'win', (45.0, XMAX - 5.0))
    Sh = energies(g, 'shell', (13.5, 64.0))
    cov = dict(windows=covered(W, 50.0, 380.25, 'windows'),
               shells=covered(Sh, 14.0, 60.0, 'shells'), h2f=covered(W, 60.0, 700.0, 'h2f'))
    bw, ww, nwin = sliding(W)
    shell = tot(Sh, 14.0, 60.0, 3)
    G = tot(W, 60.0, 700.0, 4); Ga = tot(W, 60.0, 379.5, 4); Gb = tot(W, 379.5, 700.0, 4)
    out = dict(set=SET, ball=m, tc_ball=[tc - 1 / 32, tc + 1 / 32], label='certified (arb, outward; M9 cells)',
               windows_max=[UP(v) for v in bw], windows_at=ww, n_windows=nwin, shells=[UP(v) for v in shell],
               G_60_700=[UP(v) for v in G], G_60_379p5=[UP(v) for v in Ga], G_379p5_700=[UP(v) for v in Gb],
               coverage=cov, ncells=[len(W), len(Sh)], max_rel_rad_hh=max(c['hrel'] for c in W + Sh),
               alpha=[LO(win.AL), UP(win.AL)],
               args=dict(argv=sys.argv[1:], env={k: os.environ[k] for k in ('T9TC', 'T9IFT')}, xmax=XMAX,
                         inputs={p: sha(p) for p in (os.path.abspath(__file__), IFTP, f'{P2}/t9h/win.py', f'{P2}/t9h/basegeo.py', f'{P2}/t9h/ift.py',
                                                     f'{P2}/t9h/s1chart.py', '/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py',
                                                     '/home/claude/s1spec/step4/t1/polished_v2/slab_180.pkl')}))
    fo = f'{OM}/win_arb_{SET}_{m}.json'
    json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
    print(f"{SET} ball {m}: windows {['%.4e' % v for v in out['windows_max']]} shells {['%.4e' % v for v in out['shells']]} "
          f"G {['%.4e' % v for v in out['G_60_700'][:3]]} ({time.time() - t0:.1f}s)")
