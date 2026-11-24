"""T9 items 2-3 (base values, arb) at the t_c-ball of ift.py (slab 180, t_c in [1 - 1/512, 1]):
windows [y, 1.5y], y in [50, 253]: E_j = int (28/y) (D^j v)^2 ds, D = y d/dy, s = ln y, v = f - U, U = 1.3056 alpha y^1.0246
(D^j U = 1.0246^j U); shells [14, 60]: v = w' = f - alpha psi3^'(x) - alpha^2 h2^'(x), x = y - sigma_tip (S1 = ray on sigma >= 14 up to
|theta1 - Theta| <= B(14), included), j <= 2.  f = tan(phi), phi = theta~ - Theta, y = Re(((Gamma~ - xc)/lam) e^{-i Theta});
D^j f from sympy (D = r d/dx, r = lam y / cos phi, dy/dx = cos phi / lam).  Inputs: alpha, xc, lam, Theta from out/t9_ift.json.
Per cell (length <= 0.5): centre values (point s, t_c-ball) + Taylor forms in x with hull bounds of higher D^j."""
import sys, os, json, math, time
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
import sympy as sp
from flint import arb, acb, ctx
ctx.prec = 128
import basegeo, s1chart, ift
import c4fun as C4
U_ = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
I = acb(0, 1)
R = json.load(open(os.environ.get('T9IFT', '/home/claude/s1spec/step4/p2/out/t9_ift.json')))
iv = lambda a: arb(a[0]).union(arb(a[1]))
AL, XC, LAM, TH = iv(R['alpha']), iv(R['xc']), iv(R['lam']), iv(R['Theta'])
STIP = s1chart.SIG_TIP

# ---- D^j f symbolically
yv, ph, lm = sp.symbols('y phi lam'); kk = sp.symbols('k0:5')
def ddx(e):
    out = sp.diff(e, yv) * sp.cos(ph) / lm + sp.diff(e, ph) * kk[0]
    return out + sum(sp.diff(e, kk[i]) * kk[i + 1] for i in range(4))
r_ = lm * yv / sp.cos(ph)
Df = [sp.tan(ph)]
for j in range(4): Df.append(sp.simplify(r_ * ddx(Df[-1])))
mods = [{'cos': lambda a: a.cos(), 'sin': lambda a: a.sin(), 'tan': lambda a: a.tan(), 'sec': lambda a: 1 / a.cos()}, 'math']
fDf = [sp.lambdify((yv, ph, lm, kk[0], kk[1], kk[2], kk[3]), e, modules=mods) for e in Df]
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}, 4: {1: 1, 2: 7, 3: 6, 4: 1}}

def quant(q):
    phi = q['th'] - TH
    y = (((q['Gam'] - XC) / LAM) * (-I * TH).exp()).real
    r = LAM * y / phi.cos()
    d = [fDf[j](y, phi, LAM, q['k0'], q['k1'], q['k2'], q['k3']) for j in range(5)]
    return y, r, d

def cells(g, xlo, xhi, Lc=0.5):
    for i, el in enumerate(g.els):
        a, h = g.xe_mid[i], float(g.hh[i].mid())
        if a + h < xlo or a > xhi: continue
        n = max(1, math.ceil(h / Lc)); ds = 2.0 / n
        for m in range(n):
            s0, s1 = -1 + m * ds, -1 + (m + 1) * ds
            if a + h * (s1 + 1) / 2 < xlo or a + h * (s0 + 1) / 2 > xhi: continue
            yield i, (s0 + s1) / 2, ds / 2, h

def energies(g, mode):
    """mode 'win': z_j = D^j f - 1.0246^j U (j <= 3); 'shell': z_j = D^j f - sum S(j,m) y^m (a psi^(m+1) + a^2 h2^(m+1)) (j <= 2)"""
    out = []
    xr = (45.0, 390.0) if mode == 'win' else (13.5, 64.0)
    for i, sc, sr, h in cells(g, *xr):
        qc = g.at(i, arb(sc)); qb = g.at(i, arb(sc, sr))
        yc, rc, dc = quant(qc); yb, rb, db = quant(qb)
        hx = arb(0, sr * h / 2 * 1.0000001)                       # x - x_c
        if mode == 'win':
            Uc = arb('1.3056') * AL * yc ** arb('1.0246'); Ub = arb('1.3056') * AL * yb ** arb('1.0246')
            zc = [dc[j] - arb('1.0246') ** j * Uc for j in range(5)]
            zb = [db[j] - arb('1.0246') ** j * Ub for j in range(5)]
            J = 4
        else:
            xs = yc - STIP; xsb = yb - STIP
            pc = C4.PSI.jets_any(xs, float(xs.mid()), 5); hc = C4.H2.jets_any(xs, float(xs.mid()), 5)
            pb_ = C4.PSI.jets_any(xsb, float(xsb.mid()), 5); hb = C4.H2.jets_any(xsb, float(xsb.mid()), 5)
            def Wj(j, y, P, H):
                if j == 0: return AL * P[1] + AL * AL * H[1]
                return sum((c * y ** m * (AL * P[m + 1] + AL * AL * H[m + 1]) for m, c in STIR[j].items() if m + 1 < len(P)), arb(0))
            Bs = arb(0, 2e-7)                                      # |theta1 - Theta| <= B(14) and |k1 u| f terms (sigma >= 14)
            zc = [dc[j] - Wj(j, yc, pc, hc) + (Bs if j == 0 else arb(0)) for j in range(4)]
            zb = [db[j] - Wj(j, yb, pb_, hb) + (Bs if j == 0 else arb(0)) for j in range(4)]
            J = 3
        rp = 1 + LAM * yb * (qb['th'] - TH).tan() / (qb['th'] - TH).cos() * qb['k0']      # r' on the cell
        zcell = []
        for j in range(J):
            if mode == 'win' and j + 2 <= 4:                       # second order: z_j' = z_{j+1}/r, z_j'' = (z_{j+2} - z_{j+1} r')/r^2
                z2 = (zb[j + 2] - zb[j + 1] * rp) / (rb * rb)
                zcell.append(zc[j] + zc[j + 1] / rc * hx + arb(0, U_(z2)) * hx * hx / 2)
            else:
                zcell.append(zc[j] + arb(0, U_(zb[j + 1] / rb)) * hx)
        yr = yb                                                    # y-range of the cell
        out.append(dict(ylo=float(yr.lower()), yhi=float(yr.upper()), c=[U_(arb(28) * zcell[j] * zcell[j] / (yr * yr)) * (float(yr.upper()) - float(yr.lower())) for j in range(J)]))
    return out

def sliding(cl, ya0=50.0, ya1=253.0, step=0.5):
    best = [0.0] * 4; where = [None] * 4
    ya = ya0
    while ya <= ya1 + 1e-9:
        lo, hi = ya, 1.5 * (ya + step)
        tot = [0.0] * 4
        for c in cl:
            if c['yhi'] >= lo and c['ylo'] <= hi:
                for j in range(len(c['c'])): tot[j] += c['c'][j]
        for j in range(4):
            if tot[j] > best[j]: best[j], where[j] = tot[j], ya
        ya += step
    return best, where

if __name__ == '__main__':
    t0 = time.time()
    g = basegeo.Geo(180, 400.0, ift.TC)
    print(f'geometry to x = 400 ({time.time() - t0:.0f}s); alpha {AL.str(5)}', flush=True)
    W = energies(g, 'win'); bw, ww = sliding(W)
    Sh = energies(g, 'shell')
    shell = [sum(c['c'][j] for c in Sh if c['yhi'] >= 14.0 and c['ylo'] <= 60.0) for j in range(3)]
    out = dict(args=dict(argv=sys.argv[1:], tc_ball=[float(ift.TC.lower()), float(ift.TC.upper())], env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'H1TAG', 'H1NSUB', 'H1REF', 'TREF', 'PDOT', 'T9TC', 'T9IFT', 'WINOUT', 'SUBW', 'T9JSON', 'TOPLOSS', 'N3JSON') if k in os.environ}), windows_max=bw, windows_at=ww, shells=shell, ncells=[len(W), len(Sh)], alpha=R['alpha'], tc_ball=R['tc_ball'], cpu=time.time() - t0)
    json.dump(out, open(os.environ.get('WINOUT', '/home/claude/s1spec/step4/p2/out/t9_win.json'), 'w'), indent=1)
    print('windows E0..E3 max:', ['%.3e' % v for v in bw], 'at y =', ww, '| shells E0..E2:', ['%.3e' % v for v in shell], f'({time.time() - t0:.0f}s)')
