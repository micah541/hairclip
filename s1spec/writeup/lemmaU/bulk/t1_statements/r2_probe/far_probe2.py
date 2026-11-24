"""Feasibility probe (read-only use of step-4 scripts; output only in scratchpad).
Base-curve energy F_j = int_{379.5}^{700} (28/y)(D^j(f-U))^2 ds on one t_c ball, same cell method as t9h/win.py.
Not a pinned stage; base curve only (no error term e)."""
import sys, os, json, math, time
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
import win
from flint import arb
t0=time.time()
XMAX=float(os.environ.get('XMAX','730'))
g = win.basegeo.Geo(180, XMAX, win.ift.TC)
# monkeypatch range: reuse energies() with mode 'win' but wider x-range
import types
src_xr=(45.0, XMAX-5)
def energies_far(g):
    out=[]
    for i, sc, sr, h in win.cells(g, 45.0, XMAX-5):
        qc = g.at(i, arb(sc)); qb = g.at(i, arb(sc, sr))
        yc, rc, dc = win.quant(qc); yb, rb, db = win.quant(qb)
        hx = arb(0, sr * h / 2 * 1.0000001)
        Uc = arb('1.3056') * win.AL * yc ** arb('1.0246'); Ub = arb('1.3056') * win.AL * yb ** arb('1.0246')
        zc = [dc[j] - arb('1.0246') ** j * Uc for j in range(5)]
        zb = [db[j] - arb('1.0246') ** j * Ub for j in range(5)]
        rp = 1 + win.LAM * yb * (qb['th'] - win.TH).tan() / (qb['th'] - win.TH).cos() * qb['k0']
        zcell=[]
        for j in range(4):
            if j+2<=4:
                z2 = (zb[j + 2] - zb[j + 1] * rp) / (rb * rb)
                zcell.append(zc[j] + zc[j + 1] / rc * hx + arb(0, win.U_(z2)) * hx * hx / 2)
            else:
                zcell.append(zc[j] + arb(0, win.U_(zb[j + 1] / rb)) * hx)
        yr=yb
        out.append(dict(ylo=float(yr.lower()), yhi=float(yr.upper()), z=[win.U_(zcell[j]) for j in range(4)], f=win.U_(db[0]), f1=win.U_(db[1]),
                        c=[win.U_(arb(28) * zcell[j] * zcell[j] / (yr * yr)) * (float(yr.upper()) - float(yr.lower())) for j in range(4)]))
    return out
W=energies_far(g)
def tot(lo,hi): return [sum(c['c'][j] for c in W if c['yhi']>=lo and c['ylo']<=hi) for j in range(4)]
res=dict(tc=os.environ.get('T9TC'), ymax=max(c['yhi'] for c in W), ymin=min(c['ylo'] for c in W), ncells=len(W),
  G_60_700=tot(60.0,700.0), G_60_379=tot(60.0,379.5), G_60_350=tot(60.0,350.0), F_379_700=tot(379.5,700.0), F_350_700=tot(350.0,700.0),
  wins=[tot(50*1.5**k, 50*1.5**(k+1)) for k in range(5)],
  sup_f_60_380=max(c['f'] for c in W if c['yhi']>=60 and c['ylo']<=379.5), sup_f=max(c['f'] for c in W if c['yhi']>=379.5 and c['ylo']<=700), sup_Df=max(c['f1'] for c in W if c['yhi']>=379.5 and c['ylo']<=700),
  sup_z=[max(c['z'][j] for c in W if c['yhi']>=379.5 and c['ylo']<=700) for j in range(4)], cpu=time.time()-t0)
print(json.dumps(res))
