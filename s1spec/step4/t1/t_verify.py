"""after the nrm fix: C4 jumps in arb, and exact direct residual vs kernel on a polished chain (P11 slabs 8..12)"""
import numpy as np
from flint import arb, acb, arb_poly, ctx
from base2 import Base; from arbops import A; import resid
from polish_chain import polish_chain
from rep import DS1
B = Base('P11'); SPs = polish_chain(B, 8, 5, verbose=False); SP = SPs[-1]
for th_ in (0.0, 0.37):
    Y, Yd = SP.Y_at(th_); els, gl = SP.op.elements(Y, Yd)
    jk = max(float(abs(els[e]['kp'].derivative().derivative().derivative()(arb(1)) / els[e]['h']**3 - els[e + 1]['kp'].derivative().derivative().derivative()(arb(0)) / els[e + 1]['h']**3).mid()) for e in range(len(els) - 1))
    jP = max(float(abs(els[e]['Phi'](arb(1)) - els[e + 1]['Phi'](arb(0))).mid()) for e in range(len(els) - 1))
    a = gl['a']; b = gl['k2t'] + (arb(1) / 4 - a) * A(DS1)
    NG = 50; gr = [arb.legendre_p_root(NG, k, weight=True) for k in range(NG)]; gu = [(g + 1) / 2 for g, w in gr]; gw = [w / 2 for g, w in gr]
    th = arb(-1) * arb.pi() / 2; tht = arb(0); Gm = acb(-A(DS1)); Gt = acb(0)
    sm = resid.samples(SP, th_); o, f = resid.norms(sm); X = sm['X']
    targets = [2.0, 5.0, 10.0, 20.0, 60.0]; ti = 0; out = []
    for el in els:
        he = el['h']; xe = el['x']
        thp = arb_poly([th]) + el['kp'].integral() * he; ttp = arb_poly([tht]) + el['kt'].integral() * he
        while ti < len(targets) and float((xe + he).mid()) >= targets[ti]:
            us = (A(targets[ti]) - xe) / he; Gx, Gtx = Gm, Gt
            for uq, wq in zip(gu, gw):
                uu = us * uq; T = acb(thp(uu).cos(), thp(uu).sin()); Gx += he * us * wq * T; Gtx += he * us * wq * acb(0, 1) * ttp(uu) * T
            tx = thp(us); nu = acb(0, 1) * acb(tx.cos(), tx.sin()); dot = lambda z, w: z.real * w.real + z.imag * w.imag
            r = dot(Gtx, nu) + el['ku'].derivative()(us) / he**2 - (arb(1) / 4 - a) * dot(Gx, nu) - b * nu.real
            jj = np.unravel_index(np.argmin(np.abs(X - targets[ti])), X.shape)
            out.append(f"x={targets[ti]}: direct {float(r.mid()):+.3e} kernel {f['r'][jj]:+.3e}"); ti += 1
        for uq, wq in zip(gu, gw):
            T = acb(thp(uq).cos(), thp(uq).sin()); Gm += he * wq * T; Gt += he * wq * acb(0, 1) * ttp(uq) * T
        th = thp(arb(1)); tht = ttp(arb(1))
        if ti >= len(targets): break
    print(f'theta={th_}: max C^3 jump {jk:.1e}, max Phi jump {jP:.1e} | ' + ' | '.join(out), flush=True)
    print('   norms:', {k: f'{o[k]:.2e}' for k in ('R2', 'R1', 'Rs', 'RIs', 'RF3', 'RF4')})
